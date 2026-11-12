"""Validate a spatial gene projection against measured expression."""

from .exceptions import ToolValidationError
from .spatial_transcriptomics_tool import _expression, _read_h5ad


def validate_gene_projection(request):
    import numpy as np
    from scipy.sparse import issparse
    from scipy.stats import spearmanr

    measured = _read_h5ad(request["measured_path"])
    predicted = _read_h5ad(request["projection_path"])
    for adata in (measured, predicted):
        if request.get("gene_id_case", "uppercase") == "uppercase":
            adata.var_names = adata.var_names.str.upper()
        if not adata.var_names.is_unique:
            raise ToolValidationError("Gene IDs must remain unique after case normalization.")
        _expression(adata)

    if len(measured.obs_names.difference(predicted.obs_names)) or len(
        predicted.obs_names.difference(measured.obs_names)
    ):
        raise ToolValidationError("Measured and projected files must contain the same spot IDs.")
    predicted = predicted[measured.obs_names]
    uppercase = request.get("gene_id_case", "uppercase") == "uppercase"

    def gene_ids(values):
        return [value.upper() if uppercase else value for value in values]

    programs = {name: gene_ids(genes) for name, genes in request["gene_sets"].items()}
    targets = {gene for genes in programs.values() for gene in genes}
    shared = measured.var_names.intersection(predicted.var_names)
    panel = gene_ids(request["normalization_genes"]) if "normalization_genes" in request else sorted(set(shared) - targets)
    if not panel or len(set(panel)) != len(panel) or set(panel) & targets:
        raise ToolValidationError("The normalization panel must contain distinct genes and exclude all evaluation genes.")
    if not set(panel).issubset(shared):
        raise ToolValidationError("All normalization genes must be present in both files.")
    if any(len(set(genes)) != len(genes) for genes in programs.values()):
        raise ToolValidationError("Each gene set must remain unique after case normalization.")

    measured_depth = np.asarray(measured[:, panel].X.sum(axis=1)).ravel()
    predicted_depth = np.asarray(predicted[:, panel].X.sum(axis=1)).ravel()
    if not np.isfinite(measured_depth).all() or not np.isfinite(predicted_depth).all():
        raise ToolValidationError("Normalization panel totals must be finite.")
    valid = (measured_depth > 0) & (predicted_depth > 0)
    valid_spots = int(valid.sum())
    training = predicted.uns.get("training_genes")
    training_ids = (
        set(gene_ids(list(training))) if isinstance(training, (list, tuple, np.ndarray))
        and len(training) and all(isinstance(gene, str) for gene in training) else None
    )

    def program_scores(matrix, depth, floor):
        matrix = matrix.toarray() if issparse(matrix) else matrix
        return np.log1p(matrix * 1e4 / np.maximum(depth[:, None], floor)).mean(axis=1)

    results = {}
    for name, genes in programs.items():
        missing_measured = sorted(set(genes) - set(measured.var_names))
        missing_projected = sorted(set(genes) - set(predicted.var_names))
        overlap = sorted(set(genes) & training_ids) if training_ids is not None else None
        score = None
        if missing_measured or missing_projected:
            status = "missing_genes"
        elif valid_spots < 2:
            status = "insufficient_spots"
        else:
            observed = program_scores(measured[valid, genes].X, measured_depth[valid], 1)
            projected = program_scores(predicted[valid, genes].X, predicted_depth[valid], 1e-12)
            if np.ptp(observed) == 0 or np.ptp(projected) == 0:
                status = "constant_scores"
            else:
                score = float(spearmanr(observed, projected).statistic)
                status = "ok"
        results[name] = {
            "spearman": score, "status": status,
            "missing_measured_genes": missing_measured, "missing_projected_genes": missing_projected,
            "training_gene_overlap": overlap,
            "holdout_status": "unknown" if overlap is None else "training_overlap" if overlap else "held_out",
        }

    return {
        "method": "Spearman correlation of mean log1p expression after scaling each spot to 10000 using its shared-panel total.",
        "n_spots": measured.n_obs, "valid_spots": valid_spots,
        "normalization_gene_count": len(panel), "gene_set_results": results,
    }
