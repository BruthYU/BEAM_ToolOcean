"""Offline SpatialFormer inference using explicitly provisioned official model bundles."""

import hashlib
import importlib
import json
from pathlib import Path
from threading import Lock

from .exceptions import ToolUnavailableError, ToolValidationError
from .spatial_transcriptomics_tool import (
    _SpatialTool, _artifact, _expression, _labels, _read_h5ad, _result_directory,
)
from .tool_registry import register_tool

_INFERENCE_LOCK = Lock()


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bundle(directory, mode):
    root = Path(directory).expanduser().resolve()
    try:
        manifest = json.loads((root / "manifest.json").read_text())
    except (OSError, ValueError) as exc:
        raise ToolValidationError(f"Cannot read model bundle manifest: {exc}") from exc
    if manifest.get("mode") != mode:
        raise ToolValidationError(f"This operation requires a {mode!r} model bundle.")
    for key in ("model_id", "source_commit", "package_version", "files"):
        if not manifest.get(key):
            raise ToolValidationError(f"Model manifest is missing {key!r}.")
    files = {}
    for name in ("checkpoint", "config", "tokenizer", "gene_medians"):
        entry = manifest["files"].get(name, {})
        if not entry.get("path") or not entry.get("sha256"):
            raise ToolValidationError(f"Model manifest requires path and sha256 for {name}.")
        path = (root / entry["path"]).resolve()
        if not path.is_file() or _sha256(path) != entry["sha256"]:
            raise ToolValidationError(f"Missing or checksum-mismatched model resource: {name}.")
        files[name] = path
    return manifest, files


def _prepare(request, mode):
    import numpy as np

    if request.get("species", "human") != "human":
        raise ToolValidationError("Pinned SpatialFormer inference supports the Human context token only.")
    manifest, files = _bundle(request["model_bundle_path"], mode)
    vocabulary = json.loads(files["tokenizer"].read_text())
    allowed_ids = {"tissue": set(range(8, 25)), "condition": {25, 26}, "assay": {5, 6, 7}}
    for field, token_ids in allowed_ids.items():
        if vocabulary.get(request[field]) not in token_ids:
            raise ToolValidationError(f"{field}={request[field]!r} is not a valid context token in this model vocabulary.")
    adata = _read_h5ad(request["input_path"])
    matrix = _expression(adata, request.get("count_layer"), counts=True)
    gene_key = request.get("gene_name_key", "gene_name")
    if gene_key == "var_names":
        genes = adata.var_names.astype(str)
    elif gene_key in adata.var and not adata.var[gene_key].isna().any():
        genes = adata.var[gene_key].astype(str)
    else:
        raise ToolValidationError(f"Missing or incomplete gene-name column: {gene_key!r}.")
    if len(set(genes)) != len(genes) or any(not name.strip() for name in genes):
        raise ToolValidationError("Gene names must be unique and nonempty; harmonize them before inference.")
    # Official 5k vocabulary reserves IDs 0..26 for biological context and special tokens.
    known = np.array([name in vocabulary and vocabulary[name] >= 27 and not name.startswith("<") for name in genes])
    if not known.any():
        raise ToolValidationError("No input genes match this model's gene vocabulary.")
    totals = np.asarray(matrix[:, known].sum(axis=1)).ravel()
    if (totals <= 0).any():
        raise ToolValidationError(f"{int((totals <= 0).sum())} cells have no counts in matched genes; filter them explicitly.")
    # Work on a copy: never replace original counts with ranks or truncate the output data.
    work = adata[:, known].copy()
    work.X = matrix[:, known].copy()
    work.var["gene_name"] = np.asarray(genes)[known]
    warnings = []
    if request["assay"] != "Xenium":
        warnings.append("Cross-platform inference: the bundled checkpoint was pretrained on Xenium; validate on this assay.")
    report = {
        "model_id": manifest["model_id"], "source_commit": manifest["source_commit"],
        "training_mode": manifest.get("training_mode", "unspecified"),
        "mode": mode, "species": "human", "input_cells": adata.n_obs, "input_genes": adata.n_vars,
        "matched_genes": int(known.sum()), "gene_match_fraction": float(known.mean()),
        "assay": request["assay"], "tissue": request["tissue"], "condition": request["condition"],
        "max_genes_per_cell": request.get("max_genes_per_cell", 500),
        "warnings": warnings,
        "runtime_adaptations": ["single_collate_structural_dispatch", "verified_checkpoint_legacy_torch_load"],
    }
    return adata, work, manifest, files, report



def _single_collate(batch):
    """Official single-input padding, with structural rather than gene-count dispatch."""
    import torch

    length = max(len(tokens) + len(prefix) + len(end) for tokens, prefix, end in batch)
    indices, masks, types, lengths = [], [], [], []
    for tokens, prefix, end in batch:
        sequence = prefix + tokens + end
        size = len(sequence)
        indices.append(sequence + [0] * (length - size))
        masks.append([True] * size + [False] * (length - size))
        types.append([1] * length)
        lengths.append([size, 0])
    return {"indices": torch.tensor(indices), "attention_mask": torch.tensor(masks),
            "token_type_ids": torch.tensor(types), "sequence_length": torch.tensor(lengths)}


class _CheckpointTorchProxy:
    """Scope PyTorch's legacy-checkpoint compatibility to this verified checkpoint."""
    def __init__(self, torch_module, checkpoint):
        self._torch = torch_module
        self._checkpoint = checkpoint.resolve()

    def __getattr__(self, name):
        return getattr(self._torch, name)

    def load(self, path, *args, **kwargs):
        if Path(path).resolve() != self._checkpoint:
            raise ToolValidationError("Unexpected checkpoint load during inference.")
        # The model bundle was SHA256-checked before entering upstream code.
        kwargs.setdefault("weights_only", False)
        return self._torch.load(path, *args, **kwargs)

def _infer(work, request, manifest, files, mode, **extra):
    import torch
    import spatialformer

    runtime = importlib.import_module("spatialformer.tools.get_embeddings")
    if getattr(spatialformer, "__version__", None) != manifest["package_version"]:
        raise ToolUnavailableError("Installed SpatialFormer version does not match the model bundle.")
    package_root = Path(spatialformer.__file__).parent
    for relative, expected_hash in manifest.get("runtime_source_sha256", {}).items():
        code_path = (package_root / relative).resolve()
        if not code_path.is_relative_to(package_root.resolve()) or not code_path.is_file() or _sha256(code_path) != expected_hash:
            raise ToolUnavailableError("Installed SpatialFormer source differs from the pinned model bundle.")
    if not manifest.get("runtime_source_sha256"):
        raise ToolUnavailableError("Model bundle is missing pinned runtime source checksums.")
    spatial_resource = package_root / "spatial_embeddings/xenium_5k_pandavid_gene_embeddings.pkl"
    expected = manifest.get("spatial_embeddings_sha256")
    if not expected or not spatial_resource.is_file() or _sha256(spatial_resource) != expected:
        raise ToolUnavailableError("SpatialFormer package spatial embeddings are missing or do not match the bundle.")
    requested_device = request.get("device", "auto")
    device = "cuda" if requested_device == "auto" and torch.cuda.is_available() else requested_device
    if device == "auto":
        device = "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        raise ToolUnavailableError("CUDA was requested but is unavailable.")
    if device not in ("cpu", "cuda"):
        raise ToolValidationError("device must be auto, cpu, or cuda.")
    if request.get("use_flash_attention", False) and device != "cuda":
        raise ToolValidationError("FlashAttention requires CUDA in this integration.")
    # Upstream uses a module-global device; isolate explicit device selection across calls.
    with _INFERENCE_LOCK:
        previous = runtime.device
        previous_torch = runtime.torch
        previous_collate = runtime.collate_fn
        runtime.device = torch.device(device)
        runtime.torch = _CheckpointTorchProxy(torch, files["checkpoint"])
        if mode == "single":
            runtime.collate_fn = _single_collate
        try:
            result = runtime.embed_data(
                adata=work, tissue=request["tissue"], condition=request["condition"],
                assay=request["assay"], method=request.get("method", "cls"), mode=mode,
                model_ckp_path=str(files["checkpoint"]), config_path=str(files["config"]),
                token_path=str(files["tokenizer"]), gene_median_path=str(files["gene_medians"]),
                batch_size=request.get("batch_size", 4), num_workers=request.get("num_workers", 0),
                max_len=request.get("max_genes_per_cell", 500), resume_before_5k=False,
                reverse_check=False, use_flash_attn=request.get("use_flash_attention", False),
                **extra,
            )
        finally:
            runtime.device = previous
            runtime.torch = previous_torch
            runtime.collate_fn = previous_collate
    return result, device


@register_tool("SpatialFormerTool")
class SpatialFormerTool(_SpatialTool):
    """Prediction only: no network access, training, or checkpoint mutation."""

    operations = ("extract_embeddings", "predict_pairs")

    def extract_embeddings(self, request):
        import numpy as np

        adata, work, manifest, files, report = _prepare(request, "single")
        if work.n_obs * work.n_vars * 8 > request.get("max_dense_bytes", 2_000_000_000):
            raise ToolValidationError("Official inference densifies input; split the file or increase max_dense_bytes explicitly.")
        output, device = _infer(work, request, manifest, files, "single")
        embeddings = np.asarray(output.obsm["X_SpaF"])
        if embeddings.ndim != 2 or embeddings.shape[0] != adata.n_obs or not np.isfinite(embeddings).all():
            raise ToolUnavailableError("SpatialFormer returned invalid cell embeddings.")
        if list(output.obs_names) != list(adata.obs_names):
            raise ToolUnavailableError("SpatialFormer changed the cell order.")
        report.update(device=device, embedding_dim=int(embeddings.shape[1]))
        with _result_directory(request["output_dir"], "spatialformer_embeddings") as out:
            adata.obsm["X_SpaF"] = embeddings
            adata.uns["spatialformer"] = json.dumps(report)
            path = out / "embeddings.h5ad"
            adata.write_h5ad(path)
            summary = out / "summary.json"
            summary.write_text(json.dumps(report, indent=2))
            return {"status": "success", "output_path": str(path), "summary_path": str(summary),
                    "embedding": _artifact(path, embeddings.shape, ["cell", "feature"]), **report}

    def predict_pairs(self, request):
        import numpy as np
        import pandas as pd

        adata, work, manifest, files, report = _prepare(request, "pair")
        path = Path(request["pairs_path"]).expanduser()
        pairs = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, dtype=str)
        required = ["sample_id", "cell_id_a", "cell_id_b"]
        if not set(required).issubset(pairs) or pairs.empty:
            raise ToolValidationError("Pair table must contain nonempty sample_id, cell_id_a, cell_id_b columns.")
        if pairs[required].isna().any().any() or pairs[required].astype(str).apply(lambda col: col.str.strip().eq("")).any().any():
            raise ToolValidationError("Pair IDs must not be missing or empty.")
        pairs = pairs[required].astype(str).copy()
        if len(pairs) > request.get("max_pairs", 100_000):
            raise ToolValidationError("Pair count exceeds max_pairs; provide an explicit smaller pair table.")
        if pairs.duplicated(required).any() or (pairs.cell_id_a == pairs.cell_id_b).any():
            raise ToolValidationError("Duplicate ordered pairs and self-pairs are not allowed.")
        samples = _labels(adata, request.get("sample_key", "sample_id"))
        for side in ("cell_id_a", "cell_id_b"):
            if not pairs[side].isin(adata.obs_names).all():
                raise ToolValidationError(f"Unknown cell IDs in {side}; match exact standardized observation IDs.")
            if not np.array_equal(samples.loc[pairs[side]].to_numpy(), pairs.sample_id.to_numpy()):
                raise ToolValidationError("Pair sample_id must match both cells; cross-sample pairs are not allowed.")
        if len(pairs) * work.n_vars * 16 > request.get("max_dense_bytes", 2_000_000_000):
            raise ToolValidationError("Official pair inference densifies both sides; reduce pair count or raise max_dense_bytes.")
        (embeddings, batches), device = _infer(
            work, request, manifest, files, "pair", left_cell=pairs.cell_id_a.tolist(),
            right_cell=pairs.cell_id_b.tolist(),
        )
        probabilities = np.concatenate([np.asarray(batch) for batch in batches], axis=0)
        if probabilities.shape != (len(pairs), 2) or not np.isfinite(probabilities).all():
            raise ToolUnavailableError("SpatialFormer returned an invalid binary probability matrix.")
        if (probabilities < 0).any() or (probabilities > 1).any() or not np.allclose(probabilities.sum(axis=1), 1, atol=1e-4):
            raise ToolUnavailableError("SpatialFormer pair scores are not normalized probabilities.")
        pairs["proximity_score"] = probabilities[:, 1]
        pairs["model_id"] = manifest["model_id"]
        pairs["status"] = "success"
        report.update(device=device, n_pairs=len(pairs), score_definition="P(positive), class index 1; not calibrated to population prevalence", reverse_check=False)
        with _result_directory(request["output_dir"], "spatialformer_pairs") as out:
            prediction = out / "pair_predictions.csv"
            pairs.to_csv(prediction, index=False)
            summary = out / "summary.json"
            summary.write_text(json.dumps(report, indent=2))
            result = {"status": "success", "predictions_path": str(prediction), "summary_path": str(summary), **report}
            if request.get("save_embeddings", False):
                embedding_path = out / "pair_embeddings.npy"
                np.save(embedding_path, np.asarray(embeddings))
                result["embeddings_path"] = str(embedding_path)
            return result
