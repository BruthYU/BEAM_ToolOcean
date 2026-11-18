"""Standalone IMPA worker with native output size and strict weights."""

import importlib.metadata
import importlib.util
import json
import sys
import time
from pathlib import Path

if __package__:
    from .image_evaluation_io import sha256
    from .image_evaluation_worker import read_cell
else:
    from image_evaluation_io import sha256
    from image_evaluation_worker import read_cell


def load_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def match_fingerprints(matrix, fingerprints):
    """Recover IMPA training labels by exact stored fingerprint equality."""
    import numpy as np

    labels = []
    values = fingerprints.to_numpy(dtype=np.float32)
    for row in matrix:
        hits = np.flatnonzero(np.all(values == row, axis=1))
        if len(hits) != 1:
            raise ValueError(
                "Each frozen embedding row must match exactly one named fingerprint."
            )
        labels.append(str(fingerprints.index[hits[0]]))
    if len(set(labels)) != len(labels):
        raise ValueError("Duplicate compounds in frozen IMPA embeddings.")
    return labels


def prepare_impa(job, device):
    import numpy as np
    import pandas as pd
    import torch

    source = load_source(
        "_impa_inference_model", Path(job["model_dir"]) / "IMPA/model.py"
    )
    state = torch.load(job["checkpoint_path"], map_location="cpu", weights_only=True)
    embeddings = pd.read_csv(job["embeddings_path"], index_col=0)
    if (
        not embeddings.index.is_unique
        or embeddings.shape[1] != 1024
        or not np.isfinite(embeddings.to_numpy(dtype=np.float32)).all()
    ):
        raise ValueError(
            "IMPA requires unique named, finite 1024-feature fingerprints."
        )
    frozen = torch.load(
        job["embedding_checkpoint"], map_location="cpu", weights_only=True
    )["embedding_matrix"]["weight"].numpy()
    training_labels = match_fingerprints(frozen, embeddings)
    missing = set(job["compounds"]) - set(embeddings.index)
    unseen = sorted(set(job["compounds"]) - set(training_labels))
    if missing:
        raise ValueError(f"Compounds missing from fingerprint CSV: {sorted(missing)}")
    if unseen and not job.get("allow_unseen_compounds", False):
        raise ValueError(
            f"Compounds absent from frozen training vocabulary: {unseen}; enable allow_unseen_compounds explicitly for OOD experiments."
        )
    mapping_state = state["mapping_network"]
    shape = mapping_state["mapping_network.layers.0.weight"].shape
    noise_dim = shape[1] - 1024
    if not 1 <= noise_dim <= 128 or len(mapping_state) != 2:
        raise ValueError(
            "Unsupported IMPA mapping architecture; expected one fingerprint+noise projection."
        )
    generator = source.Generator(
        img_size=96, style_dim=shape[0], in_channels=3, dim_in=64
    )
    mapping = source.MappingNetwork(
        latent_dim=shape[1],
        style_dim=shape[0],
        num_domains=len(training_labels),
        num_layers=1,
        single_style=True,
        multimodal=False,
    )
    generator.load_state_dict(state["generator"], strict=True)
    mapping.load_state_dict(mapping_state, strict=True)
    del state
    generator.eval().requires_grad_(False).to(device)
    mapping.eval().requires_grad_(False).to(device)

    def predict(images, compound):
        pixels = torch.tensor(images, dtype=torch.float32, device=device).permute(
            0, 3, 1, 2
        )
        pixels = ((pixels + torch.rand_like(pixels)) / 255.0 - 0.5) / 0.5
        condition = torch.tensor(
            embeddings.loc[compound].to_numpy(dtype=np.float32), device=device
        ).expand(len(images), -1)
        noise = torch.randn(len(images), noise_dim, device=device)
        style = mapping(torch.cat([condition, noise], dim=1))
        _, generated = generator(pixels, style)
        if not torch.isfinite(generated).all():
            raise ValueError("IMPA generated non-finite pixels.")
        return (
            ((generated * 0.5 + 0.5) * 255 + 0.5)
            .clamp(0, 255)
            .permute(0, 2, 3, 1)
            .to("cpu", torch.uint8)
            .numpy()
        )

    return predict, {
        "model": "IMPA",
        "output_size": 96,
        "noise_dimension": int(noise_dim),
        "training_condition_labels": training_labels,
        "unseen_compounds": unseen,
        "condition_mapping_verified": True,
        "preprocessing": "native 96x96 uint8 controls; official random dequantization and [-1,1] normalization",
        "method": "fingerprint_plus_gaussian_style_mapping",
    }


def main(request_path):
    import numpy as np
    import torch
    from PIL import Image

    job = json.loads(Path(request_path).read_text())
    for source in job["source_files"]:
        if sha256(source["path"]) != source["sha256"]:
            raise ValueError("An asset changed after validation.")
    for control in job["controls"]:
        read_cell(control)
    torch.set_num_threads(job.get("cpu_threads", 4))
    seed = job.get("seed", 42)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    rng = np.random.default_rng(seed)
    device = job.get("device", "auto")
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA is unavailable on this host.")
    if job["backend"] != "impa":
        raise ValueError("Unsupported image generation backend.")
    predict, model_info = prepare_impa(job, device)
    print(f"Loaded {model_info['model']} on {device}", flush=True)
    directory = Path(job["run_dir"])
    image_dir = directory / "images"
    image_dir.mkdir()
    records = []
    start = time.monotonic()
    with torch.inference_mode():
        for batch in sorted({row["batch"] for row in job["controls"]}):
            controls = [row for row in job["controls"] if row["batch"] == batch]
            count = job.get("samples_per_condition", 8)
            for compound in job["compounds"]:
                chosen = rng.choice(len(controls), count, replace=count > len(controls))
                for offset in range(0, count, job.get("batch_size", 4)):
                    selected = [
                        controls[i]
                        for i in chosen[offset : offset + job.get("batch_size", 4)]
                    ]
                    generated = predict(
                        np.stack([read_cell(row) for row in selected]), compound
                    )
                    if generated.shape != (
                        len(selected),
                        model_info["output_size"],
                        model_info["output_size"],
                        3,
                    ):
                        raise ValueError("Unexpected generated image dimensions.")
                    for row, pixels in zip(selected, generated):
                        sample_id = f"cell_{len(records):06d}"
                        path = image_dir / f"{sample_id}.png"
                        Image.fromarray(pixels).save(path)
                        records.append(
                            {
                                "sample_id": sample_id,
                                "image_path": str(path),
                                "condition": compound,
                                "batch": batch,
                                "control_sample_id": row["sample_id"],
                                "control_image_path": row["image_path"],
                                "control_sha256": row["sha256"],
                                "shape": list(pixels.shape),
                                "axes": "YXC",
                                "dtype": "uint8",
                                "condition_mapping_verified": model_info[
                                    "condition_mapping_verified"
                                ],
                                "sha256": sha256(path),
                            }
                        )
                    print(
                        f"{batch}: {compound}: {min(offset + len(selected), count)}/{count}",
                        flush=True,
                    )
    for source in job["source_files"]:
        if sha256(source["path"]) != source["sha256"]:
            raise ValueError("An asset changed during inference.")
    index = directory / "generated_images.jsonl"
    index.write_text("".join(json.dumps(row) + "\n" for row in records))
    source_record = Path(job["model_dir"]) / "SOURCE.json"
    packages = ["torch", "numpy", "Pillow", "munch"]
    result = {
        "schema_version": "1.0",
        "run_id": job["run_id"],
        **model_info,
        "n_images": len(records),
        "image_dir": str(image_dir),
        "image_index": str(index),
        "run_manifest": str(directory / "run.json"),
        "compounds": job["compounds"],
        "seed": seed,
        "device": device,
        "inference_seconds": time.monotonic() - start,
        "input_fingerprints": {
            row["name"]: row["sha256"] for row in job["source_files"]
        },
        "source": json.loads(source_record.read_text())
        if source_record.is_file()
        else None,
        "versions": {name: importlib.metadata.version(name) for name in packages},
        "benchmark_verified": False,
        "interpretation": "Unpaired condition-level population predictions, not observed single-cell trajectories.",
        "preview": records[:5],
    }
    Path(result["run_manifest"]).write_text(
        json.dumps(result, indent=2, allow_nan=False)
    )
    (directory / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False)
    )


if __name__ == "__main__":
    main(sys.argv[1])
