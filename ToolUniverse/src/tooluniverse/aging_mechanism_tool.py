"""Donor-level count modelling and local gene-set analysis for Aging."""
from pathlib import Path
from .aging_common import AgingTool, result_directory, write_table, read_table, read_expression
from .exceptions import ToolValidationError
from .tool_registry import register_tool


def _labels(frame, columns):
    for col in columns:
        if col not in frame:
            raise ToolValidationError(f"Missing metadata column: {col}")
        if frame[col].isna().any() or frame[col].astype(str).str.strip().eq("").any():
            raise ToolValidationError(f"Missing values in {col}")


def _expression(request):
    """Read observations x genes; counts validation belongs to caller."""
    import pandas as pd
    from scipy.sparse import issparse
    adata, matrix = read_expression(request["input_path"], layer=request.get("layer"))
    if adata.n_obs * adata.n_vars > request.get("max_dense_elements", 25000000):
        raise ToolValidationError("Expression exceeds max_dense_elements; subset/aggregate first")
    counts = pd.DataFrame(matrix.toarray() if issparse(matrix) else matrix,
                          index=adata.obs_names, columns=adata.var_names)
    metadata = adata.obs.copy()
    sample_key = request.get("sample_key", "sample_id")
    if sample_key not in metadata:
        metadata[sample_key] = metadata.index.astype(str)
    if request.get("metadata_path"):
        metadata = read_table(request["metadata_path"])
        key = request.get("sample_key", "sample_id")
        _labels(metadata, [key])
        metadata[key] = metadata[key].astype(str)
        if metadata[key].duplicated().any():
            raise ToolValidationError("Metadata must contain one row per sample")
        metadata = metadata.set_index(key, drop=False)
        if not counts.index.isin(metadata.index).all():
            raise ToolValidationError("Metadata is missing expression samples")
        metadata = metadata.loc[counts.index]
    return counts, metadata


def _design(request, metadata):
    """Explicit full-rank design, without evaluation of user formulas."""
    import numpy as np
    import pandas as pd
    op = request["operation"]
    donor = request.get("donor_key", "donor_id")
    covariates = request.get("covariates", [])
    categorical = request.get("categorical_covariates", [])
    _labels(metadata, [donor, *covariates, *categorical])
    if len(set(covariates + categorical)) != len(covariates + categorical):
        raise ToolValidationError("Covariates must be distinct")
    paired = op == "analyze_paired_change"
    if not paired and metadata[donor].astype(str).duplicated().any():
        raise ToolValidationError("Use one sample per donor (and one cell type) or paired analysis")
    design = pd.DataFrame({"intercept": np.ones(len(metadata))}, index=metadata.index)
    if op == "associate_with_score":
        key = request["score_key"]
        _labels(metadata, [key])
        try:
            target = pd.to_numeric(metadata[key], errors="raise").astype(float)
        except (TypeError, ValueError) as exc:
            raise ToolValidationError("Score must be numeric") from exc
        direction = f"per unit increase in {key}"
    else:
        key = request.get("group_key", "group")
        _labels(metadata, [key])
        groups = request["compare_groups"]
        if len(groups) != 2 or len(set(groups)) != 2:
            raise ToolValidationError("compare_groups requires distinct reference and test levels")
        labels = metadata[key].astype(str)
        if set(labels) != set(groups):
            raise ToolValidationError("Input must contain exactly the two requested groups")
        if labels.value_counts().min() < 2:
            raise ToolValidationError("At least two independent donors per group are required")
        target = (labels == groups[1]).astype(float)
        direction = f"{groups[1]} versus {groups[0]}"
        if paired:
            pairs = metadata.assign(_group=labels).groupby(metadata[donor], observed=True)["_group"]
            if (pairs.size() != 2).any() or (pairs.nunique() != 2).any():
                raise ToolValidationError("Each paired donor must have exactly one sample in each group")
            if metadata[donor].nunique() < 3:
                raise ToolValidationError("At least three paired donors are required")
            dummies = pd.get_dummies(metadata[donor].astype(str), prefix="donor", drop_first=True, dtype=float)
            design = pd.concat([design, dummies], axis=1)
    if key in covariates + categorical or donor in covariates + categorical:
        raise ToolValidationError("Target and donor cannot also be supplied as covariates")
    for col in covariates:
        try:
            design[f"numeric_{col}"] = pd.to_numeric(metadata[col], errors="raise").astype(float)
        except (TypeError, ValueError) as exc:
            raise ToolValidationError(f"Covariate {col} must be numeric") from exc
    for col in categorical:
        dummies = pd.get_dummies(metadata[col].astype(str), prefix=f"category_{col}", drop_first=True, dtype=float)
        design = pd.concat([design, dummies], axis=1)
    design["target"] = target
    if not np.isfinite(design.to_numpy()).all():
        raise ToolValidationError("Design values must be finite")
    if np.linalg.matrix_rank(design.to_numpy()) != design.shape[1]:
        raise ToolValidationError("Design is rank deficient: remove constant or confounded covariates")
    if len(design) - design.shape[1] < 2:
        raise ToolValidationError("At least two residual degrees of freedom are required")
    contrast = np.zeros(design.shape[1])
    contrast[-1] = 1.0
    return design, contrast, direction


@register_tool("AgingDifferentialExpressionTool")
class AgingDifferentialExpressionTool(AgingTool):
    """Negative-binomial differential expression at biological replicate level."""
    operations = ("compare_groups", "associate_with_score", "analyze_paired_change")

    def _fit(self, request):
        import numpy as np
        import pandas as pd
        from importlib.metadata import version
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        counts, metadata = _expression(request)
        if "cell_type" in request:
            cell_key = request.get("cell_type_key", "cell_type")
            _labels(metadata, [cell_key])
            selected = metadata[cell_key].astype(str) == request["cell_type"]
            if not selected.any():
                raise ToolValidationError("Requested cell type has no samples")
            counts, metadata = counts.loc[selected], metadata.loc[selected]
        if request.get("expression_scale") != "counts":
            raise ToolValidationError("Differential expression requires expression_scale=counts")
        values = counts.to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
            raise ToolValidationError("Finite nonnegative raw integer counts are required")
        if (values.sum(axis=1) <= 0).any():
            raise ToolValidationError("Every sample must have a positive library size")
        design, contrast, direction = _design(request, metadata)
        minimum = request.get("min_total_count", 10)
        keep = (values.sum(axis=0) >= minimum) & (values.sum(axis=0) > 0)
        if keep.sum() < 2:
            raise ToolValidationError("At least two genes must remain after count filtering")
        fit_counts = counts.loc[:, keep].astype("int64")
        dds = DeseqDataSet(counts=fit_counts, metadata=metadata, design=design,
                           refit_cooks=True, size_factors_fit_type="poscounts",
                           n_cpus=request.get("n_cpus", 1), quiet=True)
        dds.deseq2()
        alpha = request.get("alpha", 0.05)
        stats = DeseqStats(dds, contrast=contrast, alpha=alpha,
                           n_cpus=request.get("n_cpus", 1), quiet=True,
                           cooks_filter=True, independent_filter=True)
        stats.summary()
        result = stats.results_df.copy()
        result.index.name = "gene_id"
        result = result.rename(columns={"log2FoldChange": "log2_fold_change", "lfcSE": "lfc_se"})
        from scipy.stats import norm
        critical = norm.ppf(1 - alpha / 2)
        result["ci_low"] = result["log2_fold_change"] - critical * result["lfc_se"]
        result["ci_high"] = result["log2_fold_change"] + critical * result["lfc_se"]
        result["status"] = np.where(result["pvalue"].notna(), "tested", "not_testable_or_cooks_filtered")
        result.loc[result.pvalue.notna() & result.padj.isna(), "status"] = "independent_filtered"
        missing = pd.DataFrame(index=counts.columns[~keep])
        missing.index.name = "gene_id"
        missing["status"] = "low_count_filtered"
        result = pd.concat([result, missing]).reindex(counts.columns)
        result.index.name = "gene_id"
        result = result.replace([np.inf, -np.inf], np.nan).reset_index()
        with result_directory(request["output_dir"], request["operation"]) as out:
            return {"results": write_table(result, out, "differential_expression.csv"),
                    "design": write_table(design.reset_index(names="sample_id"), out, "design.csv"),
                    "method": "PyDESeq2 negative-binomial Wald test",
                    "software_version": version("pydeseq2"), "contrast": direction,
                    "n_samples": len(counts), "n_donors": int(metadata[request.get("donor_key", "donor_id")].nunique()),
                    "n_genes_tested": int(result.pvalue.notna().sum()),
                    "n_significant": int((result.padj < alpha).sum()),
                    "interpretation": "Exploratory association; model-derived scores may reuse tested expression features."}

    compare_groups = _fit
    associate_with_score = _fit
    analyze_paired_change = _fit


def _gene_sets(request):
    sets = request.get("gene_sets")
    if bool(sets) == bool(request.get("gmt_path")):
        raise ToolValidationError("Provide exactly one of gene_sets or gmt_path")
    if request.get("gmt_path"):
        sets = {}
        with Path(request["gmt_path"]).expanduser().open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or not fields[0] or fields[0] in sets:
                    raise ToolValidationError("GMT needs unique set names, descriptions and genes")
                sets[fields[0]] = fields[2:]
    cleaned = {}
    for name, genes in sets.items():
        if not str(name).strip() or not genes or any(not isinstance(g, str) or not g.strip() for g in genes):
            raise ToolValidationError("Gene sets need nonempty names and gene IDs")
        cleaned[name] = sorted(set(genes))
    if not cleaned:
        raise ToolValidationError("No gene sets found")
    return cleaned


def _bh(pvalues):
    import numpy as np
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(1, adjusted)
    return result


@register_tool("AgingPathwayTool")
class AgingPathwayTool(AgingTool):
    """Local-background enrichment and explicitly defined expression scores."""
    operations = ("run_preranked_gsea", "run_overrepresentation", "score_gene_sets")

    def run_overrepresentation(self, request):
        import pandas as pd
        from scipy.stats import hypergeom
        sets = _gene_sets(request)
        background = set(request["background_genes"])
        selected = set(request["selected_genes"])
        if not background or not selected or not selected.issubset(background):
            raise ToolValidationError("Selected genes must be a nonempty subset of the explicit tested background")
        rows = []
        for name, members in sets.items():
            tested = set(members) & background
            overlap = tested & selected
            if not tested:
                continue
            rows.append({"gene_set": name, "set_size_in_background": len(tested),
                         "selected_size": len(selected), "background_size": len(background),
                         "overlap_size": len(overlap), "overlap_genes": ";".join(sorted(overlap)),
                         "pvalue": float(hypergeom.sf(len(overlap) - 1, len(background), len(tested), len(selected)))})
        if not rows:
            raise ToolValidationError("No gene set overlaps the tested background")
        table = pd.DataFrame(rows)
        table["padj"] = _bh(table.pvalue)
        with result_directory(request["output_dir"], request["operation"]) as out:
            return {"results": write_table(table.sort_values("padj"), out, "overrepresentation.csv"),
                    "method": "hypergeometric upper-tail test with BH correction",
                    "n_gene_sets": len(table), "gene_set_source": request["gene_set_source"]}

    def run_preranked_gsea(self, request):
        import numpy as np
        import pandas as pd
        import gseapy
        sets = _gene_sets(request)
        table = read_table(request["input_path"])
        gene_key, rank_key = request.get("gene_key", "gene_id"), request.get("rank_key", "stat")
        _labels(table, [gene_key])
        if rank_key not in table:
            raise ToolValidationError(f"Missing rank column {rank_key}")
        if table[gene_key].duplicated().any():
            raise ToolValidationError("GSEA requires unique gene IDs")
        ranks = pd.to_numeric(table[rank_key], errors="coerce")
        valid = np.isfinite(ranks.to_numpy())
        rnk = pd.DataFrame({"gene": table[gene_key].astype(str), "rank": ranks}).loc[valid]
        rnk = rnk.sort_values(["rank", "gene"], ascending=[False, True], kind="stable")
        if len(rnk) < 3 or rnk["rank"].nunique() < 2:
            raise ToolValidationError("GSEA needs at least three finite ranked genes and variable ranks")
        minimum, maximum = request.get("min_size", 5), request.get("max_size", 500)
        if minimum > maximum:
            raise ToolValidationError("min_size must not exceed max_size")
        universe = set(rnk.gene)
        usable = {name: [g for g in genes if g in universe] for name, genes in sets.items()}
        usable = {name: genes for name, genes in usable.items() if minimum <= len(genes) <= maximum and len(genes) < len(universe)}
        if not usable:
            raise ToolValidationError("No gene sets pass overlap size limits")
        with result_directory(request["output_dir"], request["operation"]) as out:
            result = gseapy.prerank(rnk=rnk, gene_sets=usable, min_size=minimum, max_size=maximum,
                                   permutation_num=request.get("n_permutations", 1000),
                                   seed=request.get("random_seed", 0), threads=request.get("n_cpus", 1),
                                   outdir=None, no_plot=True, verbose=False)
            enrichment = result.res2d.rename(columns={"Term": "gene_set", "ES": "es", "NES": "nes",
                                                      "NOM p-val": "pvalue", "FDR q-val": "fdr_qvalue",
                                                      "Lead_genes": "leading_edge_genes"})
            enrichment = enrichment.replace([np.inf, -np.inf], np.nan)
            return {"results": write_table(enrichment, out, "gsea.csv"),
                    "method": "GSEApy preranked gene-set permutation GSEA", "software_version": gseapy.__version__,
                    "n_gene_sets": len(enrichment), "n_ranked_genes": len(rnk),
                    "n_excluded_nonfinite": int((~valid).sum()), "n_tied_ranks": int(rnk["rank"].duplicated().sum()),
                    "gene_set_source": request["gene_set_source"]}

    def score_gene_sets(self, request):
        import numpy as np
        import pandas as pd
        sets = _gene_sets(request)
        expression, _ = _expression(request)
        if request.get("expression_scale") != "log_normalized":
            raise ToolValidationError("Mean expression scoring requires expression_scale=log_normalized")
        scores = pd.DataFrame({"observation_id": expression.index})
        coverage = []
        for name, genes in sets.items():
            if name == "observation_id":
                raise ToolValidationError("Gene-set name observation_id is reserved")
            matched = expression.columns.intersection(genes, sort=False)
            fraction = len(matched) / len(genes)
            usable = len(matched) >= request.get("min_genes", 1) and fraction >= request.get("min_coverage", 0.5)
            scores[name] = expression.loc[:, matched].mean(axis=1).to_numpy() if usable else np.nan
            coverage.append({"gene_set": name, "n_genes": len(genes), "n_matched": len(matched),
                             "coverage": fraction, "status": "scored" if usable else "insufficient_coverage"})
        with result_directory(request["output_dir"], request["operation"]) as out:
            return {"results": write_table(scores, out, "gene_set_scores.csv"),
                    "coverage": write_table(pd.DataFrame(coverage), out, "gene_set_coverage.csv"),
                    "method": "mean of supplied log-normalized expression; not a calibrated senescence probability",
                    "n_gene_sets": len(sets), "gene_set_source": request["gene_set_source"]}
