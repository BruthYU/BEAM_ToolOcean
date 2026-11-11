"""Local spatial transcriptomics tools for AnnData files."""

from contextlib import contextmanager
from pathlib import Path
import shutil
import tempfile
from threading import Lock

from .base_tool import BaseTool
from .exceptions import ToolUnavailableError, ToolValidationError
from .tool_registry import register_tool

_PYRO_LOCK = Lock()


def _read_h5ad(file_path):
    from anndata import read_h5ad

    try:
        adata = read_h5ad(Path(file_path).expanduser())
    except (OSError, ValueError) as exc:
        raise ToolValidationError(f"Cannot read AnnData file {file_path!r}: {exc}") from exc
    if not adata.n_obs or not adata.n_vars:
        raise ToolValidationError(f"{file_path!r} must contain observations and genes.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ToolValidationError(f"{file_path!r} must have unique observation and gene IDs.")
    return adata


def _expression(adata, layer=None, *, counts=False):
    import numpy as np
    from scipy.sparse import issparse

    if layer is not None and layer not in adata.layers:
        raise ToolValidationError(f"Count layer {layer!r} was not found.")
    matrix = adata.X if layer is None else adata.layers[layer]
    if matrix is None:
        raise ToolValidationError("The selected expression matrix is empty.")
    values = matrix.data if issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ToolValidationError("Expression values must be finite and nonnegative.")
    if counts and (values != np.floor(values)).any():
        raise ToolValidationError("The selected matrix must contain raw integer counts.")
    return matrix


def _labels(adata, key):
    if key not in adata.obs:
        raise ToolValidationError(f"Observation column {key!r} was not found.")
    labels = adata.obs[key]
    if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
        raise ToolValidationError(f"Observation column {key!r} contains missing labels.")
    return labels.astype(str)


def _common_genes(left, right):
    genes = left.intersection(right, sort=False)
    if not len(genes):
        raise ToolValidationError("The inputs must have at least one shared gene ID.")
    return genes


@contextmanager
def _result_directory(output_dir, operation):
    # Isolate each call so retries and concurrent calls cannot overwrite results.
    parent = Path(output_dir).expanduser().resolve()
    try:
        parent.mkdir(parents=True, exist_ok=True)
        result_dir = Path(tempfile.mkdtemp(prefix=f"{operation}_", dir=parent))
    except OSError as exc:
        raise ToolValidationError(f"Cannot create an output directory at {parent}: {exc}") from exc
    try:
        yield result_dir
    except BaseException:
        shutil.rmtree(result_dir, ignore_errors=True)
        raise


def _artifact(path, shape, axes):
    return {"path": str(path), "shape": list(shape), "axes": list(axes)}


class _SpatialTool(BaseTool):
    operations = ()

    def validate_parameters(self, arguments):
        error = super().validate_parameters(arguments)
        request = arguments.get("request")
        if error is None or not isinstance(request, dict):
            return error

        from jsonschema import Draft202012Validator
        from jsonschema.exceptions import best_match

        operation = request.get("operation")
        if operation not in self.operations:
            return ToolValidationError(
                f"Choose a request.operation from {', '.join(self.operations)}."
            )
        branches = self.tool_config["parameter"]["properties"]["request"]["oneOf"]
        for branch in branches:
            properties = branch["properties"]
            if operation not in properties["operation"]["enum"]:
                continue
            if "mode" in properties and request.get("mode", "cells") not in properties["mode"]["enum"]:
                continue
            problem = best_match(Draft202012Validator(branch).iter_errors(request))
            if problem is not None:
                path = self._format_instance_path(["request", *problem.path])
                return ToolValidationError(f"Invalid {operation} parameter '{path}'. {problem.message}")
        return error

    def run(self, arguments):
        request = arguments["request"]
        operation = request["operation"]
        if operation not in self.operations:
            raise ToolValidationError(f"Unsupported operation: {operation!r}.")
        try:
            return {"operation": operation, **getattr(self, operation)(request)}
        except ImportError as exc:
            raise ToolUnavailableError(
                f"{operation} could not import its local dependencies: {exc}",
                retriable=False,
                next_steps=["Install the dependencies in the Spatial Transcriptomics guide."],
            ) from exc


@register_tool("Cell2locationTool")
class Cell2locationTool(_SpatialTool):
    operations = ("fit_reference", "map_abundance")

    def fit_reference(self, request):
        import pyro
        from cell2location.models import RegressionModel

        adata = _read_h5ad(request["reference_path"])
        layer = request.get("count_layer")
        _expression(adata, layer, counts=True)
        label_key = request["cell_type_key"]
        adata.obs[label_key] = _labels(adata, label_key)
        batch_key = request.get("batch_key")
        if batch_key is not None:
            _labels(adata, batch_key)

        with _PYRO_LOCK, pyro.get_param_store().scope():
            RegressionModel.setup_anndata(
                adata, layer=layer, labels_key=label_key, batch_key=batch_key
            )
            model = RegressionModel(adata)
            model.train(
                max_epochs=request.get("max_epochs", 250),
                accelerator="cpu",
                enable_progress_bar=False,
            )
            adata = model.export_posterior(adata, sample_kwargs={"accelerator": "cpu"})
            types = list(adata.uns["mod"]["factor_names"])
            key = "means_per_cluster_mu_fg"
            columns = [f"{key}_{cell_type}" for cell_type in types]
            source = adata.varm[key] if key in adata.varm else adata.var
            signatures = source.loc[:, columns].copy()
            signatures.columns = types

        with _result_directory(request["output_dir"], "fit_reference") as output:
            path = output / "reference_signatures.csv"
            signatures.to_csv(path, index_label="gene_id")
        return {"reference_signatures": _artifact(path, signatures.shape, ("gene", "cell_type"))}

    def map_abundance(self, request):
        import csv

        import pandas as pd
        import pyro
        from cell2location.models import Cell2location

        adata = _read_h5ad(request["spatial_path"])
        layer = request.get("count_layer")
        _expression(adata, layer, counts=True)
        batch_key = request.get("batch_key")
        if batch_key is not None:
            _labels(adata, batch_key)

        signature_path = Path(request["reference_signatures_path"]).expanduser()
        try:
            with signature_path.open(newline="", encoding="utf-8-sig") as stream:
                header = next(csv.reader(stream), [])
                if len(header) < 2 or len(set(header[1:])) != len(header) - 1:
                    raise ValueError("The CSV must have a gene ID column and unique cell type columns.")
                if any(not label.strip() for label in header[1:]):
                    raise ValueError("Cell type names must not be empty.")
                stream.seek(0)
                table = pd.read_csv(stream, dtype=str, keep_default_na=False)
            signatures = table.iloc[:, 1:].astype(float)
            signatures.columns = header[1:]
            signatures.index = table.iloc[:, 0]
        except (OSError, ValueError) as exc:
            raise ToolValidationError(f"Cannot read reference signatures: {exc}") from exc
        import numpy as np

        if not signatures.index.is_unique or (signatures.index.str.strip() == "").any():
            raise ToolValidationError("Reference signature gene IDs must be unique and nonempty.")
        values = signatures.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise ToolValidationError("Reference signatures must be finite and nonnegative.")
        genes = _common_genes(adata.var_names, signatures.index)
        aligned = adata[:, genes].copy()

        with _PYRO_LOCK, pyro.get_param_store().scope():
            Cell2location.setup_anndata(aligned, layer=layer, batch_key=batch_key)
            model = Cell2location(
                aligned,
                cell_state_df=signatures.loc[genes],
                N_cells_per_location=request["n_cells_per_location"],
                detection_alpha=request.get("detection_alpha", 20),
            )
            model.train(
                max_epochs=request.get("max_epochs", 30000),
                accelerator="cpu",
                enable_progress_bar=False,
            )
            aligned = model.export_posterior(
                aligned,
                add_to_obsm=["means", "q05"],
                sample_kwargs={"accelerator": "cpu"},
            )

        keys = ["means_cell_abundance_w_sf", "q05_cell_abundance_w_sf"]
        types = list(aligned.uns["mod"]["factor_names"])
        for key in keys:
            abundance = aligned.obsm[key].copy()
            abundance.columns = types
            adata.obsm[key] = abundance
        with _result_directory(request["output_dir"], "map_abundance") as output:
            path = output / "spatial_abundance.h5ad"
            adata.write_h5ad(path)
        return {
            "spatial_abundance": _artifact(path, adata.obsm[keys[0]].shape, ("spot", "cell_type")),
            "abundance_keys": [f"obsm/{key}" for key in keys],
            "n_shared_genes": len(genes),
        }


@register_tool("TangramTool")
class TangramTool(_SpatialTool):
    operations = ("map_to_space",)

    def map_to_space(self, request):
        import tangram as tg

        reference = _read_h5ad(request["reference_path"])
        spatial = _read_h5ad(request["spatial_path"])
        for adata in (reference, spatial):
            _expression(adata)
            if not adata.var_names.str.lower().is_unique:
                raise ToolValidationError("Tangram requires gene IDs unique after lowercasing.")
        mode = request.get("mode", "cells")
        label_key = request.get("cell_type_key")
        project_types = request.get("project_cell_types", False)
        if mode == "clusters" or project_types:
            reference.obs[label_key] = _labels(reference, label_key)

        # Basic mapping uses expression alone. Keep the original spatial metadata
        # for the output without computing an unused Squidpy neighbor graph.
        mapping_spatial = spatial.copy()
        if "spatial" in mapping_spatial.obsm:
            del mapping_spatial.obsm["spatial"]
        tg.pp_adatas(reference, mapping_spatial, genes=request.get("training_genes"))
        if not mapping_spatial.uns["training_genes"]:
            raise ToolValidationError("Tangram requires nonzero shared training genes.")
        cluster_label = label_key if mode == "clusters" else None
        mapping = tg.map_cells_to_space(
            reference,
            mapping_spatial,
            mode=mode,
            cluster_label=cluster_label,
            device="cpu",
            num_epochs=request.get("max_epochs", 1000),
            verbose=False,
        )
        if project_types:
            tg.project_cell_annotations(mapping, spatial, annotation=label_key)
        projected = None
        if request.get("project_genes", False):
            projected = tg.project_genes(mapping, reference, cluster_label=cluster_label)
            if "spatial" in spatial.obsm:
                projected.obsm["spatial"] = spatial.obsm["spatial"].copy()

        with _result_directory(request["output_dir"], "map_to_space") as output:
            path = output / "mapping.h5ad"
            mapping.write_h5ad(path)
            axis = "cell" if mode == "cells" else "cell_group"
            result = {
                "mapping": _artifact(path, mapping.shape, (axis, "spot")),
                "n_training_genes": len(mapping_spatial.uns["training_genes"]),
            }
            if project_types:
                path = output / "cell_type_projection.h5ad"
                spatial.write_h5ad(path)
                projection = spatial.obsm["tangram_ct_pred"]
                result["cell_type_projection"] = _artifact(
                    path, projection.shape, ("spot", "cell_type")
                )
                result["cell_type_key"] = "obsm/tangram_ct_pred"
            if projected is not None:
                path = output / "gene_projection.h5ad"
                projected.write_h5ad(path)
                result["gene_projection"] = _artifact(path, projected.shape, ("spot", "gene"))
        return result


def _spatial_input(file_path):
    import numpy as np

    adata = _read_h5ad(file_path)
    coordinates = adata.obsm.get("spatial")
    if coordinates is None:
        raise ToolValidationError(f"{file_path!r} must contain coordinates in obsm['spatial'].")
    try:
        coordinates = np.asarray(coordinates, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ToolValidationError(f"{file_path!r} must contain numeric spatial coordinates.") from exc
    if coordinates.shape != (adata.n_obs, 2) or not np.isfinite(coordinates).all():
        raise ToolValidationError(f"{file_path!r} must have two finite coordinates per spot.")
    adata.obsm["spatial"] = coordinates
    _expression(adata, counts=True)
    return adata


def _graph_input(adata):
    from anndata import AnnData

    # Train on fresh data so a prior analysis's graph, scaled features, or HVG
    # flags cannot change the model's preprocessing.
    spatial = {"spatial": adata.obsm["spatial"]} if "spatial" in adata.obsm else {}
    return AnnData(
        X=adata.X.copy(),
        obs=adata.obs.iloc[:, :0].copy(),
        var=adata.var.iloc[:, :0].copy(),
        obsm=spatial,
    )


def _graph_domains(adata, request):
    import scanpy as sc
    from GraphST.GraphST import GraphST

    if adata.n_obs < 4:
        raise ToolValidationError("GraphST requires at least four spots to build its spatial graph.")
    model_input = _graph_input(adata)
    trained = GraphST(model_input, epochs=request.get("max_epochs", 600)).train()
    sc.pp.neighbors(trained, use_rep="emb", n_neighbors=min(15, trained.n_obs - 1))
    sc.tl.leiden(
        trained, resolution=request.get("resolution", 1.0), key_added="domain",
        flavor="igraph", n_iterations=2, directed=False, random_state=0,
    )
    adata.obsm["emb"] = trained.obsm["emb"]
    adata.obs["domain"] = trained.obs["domain"]
    return adata


def _write_domains(adata, request):
    with _result_directory(request["output_dir"], request["operation"]) as output:
        path = output / "spatial_domains.h5ad"
        adata.write_h5ad(path)
    return {
        "domains": _artifact(path, adata.obsm["emb"].shape, ("spot", "embedding_feature")),
        "embedding_key": "obsm/emb",
        "domain_key": "obs/domain",
        "n_domains": int(adata.obs["domain"].nunique()),
    }


@register_tool("GraphSTTool")
class GraphSTTool(_SpatialTool):
    operations = ("cluster_domains", "integrate_slices", "deconvolve")

    def cluster_domains(self, request):
        spatial = _spatial_input(request["spatial_path"])
        return _write_domains(_graph_domains(spatial, request), request)

    def integrate_slices(self, request):
        from anndata import concat

        sections = {}
        shared = None
        for entry in request["slices"]:
            slice_id = entry["slice_id"]
            if slice_id in sections:
                raise ToolValidationError(f"Duplicate slice_id: {slice_id!r}.")
            section = _spatial_input(entry["spatial_path"])
            section.obs["original_spot_id"] = section.obs_names
            sections[slice_id] = section
            shared = section.var_names if shared is None else _common_genes(shared, section.var_names)

        combined = concat(sections, join="inner", label="slice_id", index_unique=":")
        combined.obs_names_make_unique()
        result = _write_domains(_graph_domains(combined, request), request)
        result.update({"slice_key": "obs/slice_id", "n_slices": len(sections), "n_shared_genes": len(shared)})
        return result

    def deconvolve(self, request):
        import numpy as np
        import pandas as pd
        from GraphST.GraphST import GraphST
        from GraphST.preprocess import filter_with_overlap_gene, get_feature, preprocess

        spatial = _spatial_input(request["spatial_path"])
        if spatial.n_obs < 4:
            raise ToolValidationError("GraphST requires at least four spots to build its spatial graph.")
        reference = _read_h5ad(request["reference_path"])
        _expression(reference, counts=True)
        labels = _labels(reference, request["cell_type_key"])
        _common_genes(spatial.var_names, reference.var_names)
        spatial_input, reference_input = _graph_input(spatial), _graph_input(reference)
        preprocess(spatial_input)
        preprocess(reference_input)
        _common_genes(
            spatial_input.var_names[spatial_input.var["highly_variable"]],
            reference_input.var_names[reference_input.var["highly_variable"]],
        )
        spatial_input, reference_input = filter_with_overlap_gene(spatial_input, reference_input)
        spatial_input, reference_input = spatial_input.copy(), reference_input.copy()
        get_feature(spatial_input)
        model = GraphST(
            spatial_input, reference_input, deconvolution=True,
            epochs=request.get("max_epochs", 1200),
        )
        mapped, mapped_reference = model.train_map()
        weights = mapped.obsm["map_matrix"]
        ordered_labels = labels.loc[mapped_reference.obs_names]
        # Upstream keeps floor(10% of cells) per spot. Keep at least one
        # for very small reference sets so each spot has a composition.
        n_keep = max(1, weights.shape[1] // 10)
        top = np.argsort(weights, axis=1, kind="stable")[:, -n_keep:]
        row = np.arange(weights.shape[0])[:, None]
        retained = np.zeros_like(weights)
        retained[row, top] = weights[row, top]
        indicators = pd.get_dummies(ordered_labels, dtype=float)
        fractions = pd.DataFrame(
            retained @ indicators.to_numpy(), index=spatial.obs_names, columns=indicators.columns
        )
        fractions = fractions.div(fractions.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        spatial.obsm["map_matrix"] = weights
        spatial.obsm["cell_type_proportions"] = fractions
        spatial.uns["map_reference_cell_ids"] = mapped_reference.obs_names.to_numpy(dtype=str)

        with _result_directory(request["output_dir"], "deconvolve") as output:
            path = output / "spatial_deconvolution.h5ad"
            spatial.write_h5ad(path)
        return {
            "mapping": _artifact(path, weights.shape, ("spot", "reference_cell")),
            "cell_type_proportions": _artifact(path, fractions.shape, ("spot", "cell_type")),
            "mapping_key": "obsm/map_matrix",
            "proportions_key": "obsm/cell_type_proportions",
            "reference_cell_ids_key": "uns/map_reference_cell_ids",
            "n_shared_genes": int(spatial_input.n_vars),
        }
