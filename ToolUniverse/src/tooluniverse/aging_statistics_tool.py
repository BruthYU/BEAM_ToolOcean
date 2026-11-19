"""Donor-level aging summaries, associations, composition and validation."""
from .aging_common import AgingTool, result_directory, read_table, write_table
from .exceptions import ToolValidationError
from .tool_registry import register_tool


def _columns(df, keys):
    missing = [k for k in keys if k not in df]
    if missing:
        raise ToolValidationError(f"Missing columns: {missing}")
    if df.empty:
        raise ToolValidationError("Input table is empty")


def _labels(df, keys):
    _columns(df, keys)
    for key in keys:
        if df[key].isna().any() or df[key].astype(str).str.strip().eq("").any():
            raise ToolValidationError(f"Missing or empty labels in {key}")


def _numeric(df, keys):
    import numpy as np
    import pandas as pd
    _columns(df, keys)
    for key in keys:
        df[key] = pd.to_numeric(df[key], errors="coerce")
        if not np.isfinite(df[key].to_numpy(dtype=float)).all():
            raise ToolValidationError(f"{key} must contain finite numeric values")


def _unique(df, keys):
    _labels(df, keys)
    if df.duplicated(keys).any():
        raise ToolValidationError(f"Expected one independent observation per {keys}; aggregate cells or specify strata/repeated-measure design")


def _strata(df, keys):
    _labels(df, keys)
    if not keys:
        return [({}, df.copy())]
    return [(dict(zip(keys, k if isinstance(k, tuple) else (k,))), v.copy())
            for k, v in df.groupby(keys, observed=True, sort=False, dropna=False)]


def _bh(values):
    import numpy as np
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    if len(order):
        out[order] = np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order)+1))[::-1])[::-1].clip(0, 1)
    return out


def _output(r, frame, **extra):
    with result_directory(r["output_dir"], r["operation"]) as out:
        return {"table": write_table(frame, out, "results.csv"), **extra}


def _fit(df, outcome, predictor, covariates=(), categorical=()):
    """Full-rank ordinary least squares with HC3 standard errors, t intervals."""
    import numpy as np
    import pandas as pd
    from scipy.stats import t
    covariates, categorical = list(covariates), list(categorical)
    if len(set([predictor, *covariates])) != 1 + len(covariates):
        raise ToolValidationError("Predictor and covariates must be distinct")
    if outcome in [predictor, *covariates]:
        raise ToolValidationError("Outcome cannot also be a predictor or covariate")
    if not set(categorical).issubset(covariates):
        raise ToolValidationError("categorical_covariates must be a subset of covariates")
    _numeric(df, [outcome, predictor, *[k for k in covariates if k not in categorical]])
    _labels(df, categorical)
    parts = [np.ones((len(df), 1)), df[[predictor]].to_numpy(dtype=float)]
    for key in covariates:
        if key in categorical:
            part = pd.get_dummies(df[key].astype(str), drop_first=True, dtype=float)
            if not part.shape[1]:
                raise ToolValidationError(f"Categorical covariate {key} has only one level")
            parts.append(part.to_numpy())
        else:
            parts.append(df[[key]].to_numpy(dtype=float))
    x, y = np.column_stack(parts), df[outcome].to_numpy(dtype=float)
    if len(y) <= x.shape[1] or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ToolValidationError("Design is rank deficient or has no residual degrees of freedom")
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    bread = np.linalg.inv(x.T @ x)
    leverage = np.sum((x @ bread) * x, axis=1)
    if np.any(leverage >= 1 - 1e-10):
        raise ToolValidationError("Design has unit leverage; HC3 inference is not estimable")
    variance = bread @ (x.T @ (x * ((residual / (1-leverage))**2)[:, None])) @ bread
    se, estimate, df_resid = float(np.sqrt(max(0, variance[1, 1]))), float(beta[1]), len(y)-x.shape[1]
    if se <= 1e-14:
        raise ToolValidationError("Residual uncertainty is zero; statistical inference is not estimable")
    ci = float(t.ppf(.975, df_resid) * se)
    return {"estimate": estimate, "standard_error": se, "ci_lower": estimate-ci,
            "ci_upper": estimate+ci, "p_value": float(2*t.sf(abs(estimate/se), df_resid)),
            "n_donors": len(y), "df_residual": df_resid, "method": "OLS_HC3_t"}


def _association(r, frame=None):
    import pandas as pd
    df = read_table(r["input_path"]) if frame is None else frame.copy()
    donor, strata, outcomes = r.get("donor_key", "donor_id"), r.get("strata", []), r["value_columns"]
    _unique(df, [donor, *strata])
    _numeric(df, outcomes)
    group = r.get("group_key")
    predictor = group if group else r["predictor_key"]
    _labels(df, [predictor])
    rows = []
    for labels, part in _strata(df, strata):
        if group:
            levels = r["contrast"]
            if len(levels) != 2 or len(set(levels)) != 2:
                raise ToolValidationError("contrast must be [reference, treatment]")
            part = part[part[group].astype(str).isin(levels)].copy()
            counts = part.groupby(group, observed=True).size()
            if len(counts) != 2 or min(counts) < 2:
                raise ToolValidationError("Each contrast group needs at least two donors per stratum")
            part["__aging_predictor"] = (part[group].astype(str) == levels[1]).astype(float)
            fit_predictor = "__aging_predictor"
        else:
            fit_predictor = predictor
        for outcome in outcomes:
            result = _fit(part.copy(), outcome, fit_predictor, r.get("covariates", []), r.get("categorical_covariates", []))
            rows.append({**labels, "readout": outcome, "predictor": predictor,
                         "contrast": " - ".join(reversed(r["contrast"])) if group else "per unit", **result})
    result = pd.DataFrame(rows)
    result["q_value"] = _bh(result.p_value)
    return result


@register_tool("AgingJointSummaryTool")
class AgingJointSummaryTool(AgingTool):
    operations = ("join_model_results", "summarize_by_donor", "summarize_by_cell_type")

    def join_model_results(self, r):
        keys, merged, coverage, names = r["join_keys"], None, [], set()
        all_readouts, all_metadata = set(), set()
        for spec in r["inputs"]:
            name = spec["name"]
            if name in names:
                raise ToolValidationError("Input names must be unique")
            names.add(name)
            df = read_table(spec["path"])
            _unique(df, keys)
            columns = spec["value_columns"]
            metadata = [k for k in spec.get("metadata_columns", []) if k not in keys]
            if set(columns) & set([*keys, *metadata]):
                raise ToolValidationError("Readouts must be distinct from join keys and metadata columns")
            new_readouts = {f"{name}__{k}" for k in columns}
            if set(metadata) & all_readouts or new_readouts & (all_metadata | set(keys)):
                raise ToolValidationError("Metadata names collide with prefixed model readouts")
            all_readouts.update(new_readouts)
            all_metadata.update(metadata)
            _numeric(df, columns)
            _labels(df, metadata)
            df = df[keys + metadata + columns].rename(columns={k: f"{name}__{k}" for k in columns})
            if not df.columns.is_unique:
                raise ToolValidationError("Metadata names collide with prefixed model readouts")
            coverage.append({"name": name, "n_rows": len(df)})
            if merged is None:
                merged = df
            else:
                overlaps = [k for k in metadata if k in merged]
                collisions = (set(df.columns) & set(merged.columns)) - set(keys) - set(overlaps)
                if collisions:
                    raise ToolValidationError(f"Conflicting model readout/metadata names: {sorted(collisions)}")
                suffix = "__aging_incoming_metadata"
                if any(k.endswith(suffix) for k in [*df.columns, *merged.columns]):
                    raise ToolValidationError("Reserved metadata merge suffix used in a column name")
                merged = merged.merge(df, on=keys, how=r.get("how", "outer"), validate="one_to_one", suffixes=("", suffix))
                for key in overlaps:
                    incoming = key + suffix
                    available = merged[key].notna() & merged[incoming].notna()
                    if not merged.loc[available, key].astype(str).eq(merged.loc[available, incoming].astype(str)).all():
                        raise ToolValidationError(f"Conflicting metadata for joined observations: {key}")
                    merged[key] = merged[key].combine_first(merged[incoming])
                    merged = merged.drop(columns=incoming)
        if merged is None or merged.empty:
            raise ToolValidationError("No matched results; verify common data grain and explicit IDs")
        for entry in coverage:
            cols = [k for k in merged if k.startswith(entry["name"] + "__")]
            entry["missing_joined_rows"] = int(merged[cols].isna().all(axis=1).sum())
        return _output(r, merged, coverage=coverage, joined_rows=len(merged))

    def summarize_by_donor(self, r):
        return self._summarize(r, False)

    def summarize_by_cell_type(self, r):
        return self._summarize(r, True)

    def _summarize(self, r, cell_type):
        df = read_table(r["input_path"])
        keys = [r.get("donor_key", "donor_id"), *r.get("strata", [])]
        if cell_type and r.get("cell_type_key", "cell_type") not in keys:
            keys.append(r.get("cell_type_key", "cell_type"))
        _labels(df, keys)
        if r.get("observation_key"):
            _unique(df, [r["observation_key"]])
        values, metadata = r["value_columns"], r.get("metadata_columns", [])
        if set(values) & set([*keys, *metadata]) or "n_observations" in [*values, *keys, *metadata]:
            raise ToolValidationError("Readouts, summary keys and metadata must not overlap or use reserved n_observations")
        _numeric(df, values)
        _labels(df, metadata)
        grouped = df.groupby(keys, observed=True, sort=False, dropna=False)
        for key in metadata:
            if grouped[key].nunique().gt(1).any():
                raise ToolValidationError(f"Metadata {key} varies within a summary unit; add a sample/timepoint stratum")
        agg = grouped[values].agg(r.get("statistic", "mean"))
        agg["n_observations"] = grouped.size()
        for key in metadata:
            if key not in keys:
                agg[key] = grouped[key].first()
        return _output(r, agg.reset_index(), unit=keys, statistic=r.get("statistic", "mean"))


@register_tool("AgingAssociationTool")
class AgingAssociationTool(AgingTool):
    operations = ("compare_groups", "analyze_age_association", "analyze_paired_change", "compare_model_readouts")

    def compare_groups(self, r):
        return _output(r, _association(r), independent_unit=r.get("donor_key", "donor_id"), correction="BH across all reported tests")

    def analyze_age_association(self, r):
        return self.compare_groups(r)

    def compare_model_readouts(self, r):
        import pandas as pd
        from scipy.stats import spearmanr
        df, strata = read_table(r["input_path"]), r.get("strata", [])
        _unique(df, [r.get("donor_key", "donor_id"), *strata])
        a, b = r["readouts"]
        _numeric(df, [a, b])
        rows = []
        for labels, part in _strata(df, strata):
            if len(part) < 3 or part[a].nunique() < 2 or part[b].nunique() < 2:
                raise ToolValidationError("Correlation needs at least three donors and nonconstant readouts")
            rho, p = spearmanr(part[a], part[b])
            rows.append({**labels, "readout_a": a, "readout_b": b, "spearman_rho": float(rho), "p_value": float(p), "n_donors": len(part)})
        result = pd.DataFrame(rows)
        result["q_value"] = _bh(result.p_value)
        return _output(r, result, interpretation="Association of related expression-derived readouts; not independent validation")

    def analyze_paired_change(self, r):
        import numpy as np
        import pandas as pd
        from scipy.stats import t
        df = read_table(r["input_path"])
        donor, time, strata = r.get("donor_key", "donor_id"), r["time_key"], r.get("strata", [])
        _unique(df, [donor, time, *strata])
        _numeric(df, r["value_columns"])
        levels, rows, dropped = r["contrast"], [], []
        for labels, part in _strata(df, strata):
            part = part[part[time].astype(str).isin(levels)].copy()
            part[time] = part[time].astype(str)
            for value in r["value_columns"]:
                wide = part.pivot(index=donor, columns=time, values=value).reindex(columns=levels)
                complete = wide.dropna()
                dropped.append({**labels, "readout": value, "excluded_donors": wide.index[wide.isna().any(axis=1)].astype(str).tolist()})
                if len(complete) < 3:
                    raise ToolValidationError("Paired inference requires at least three complete donor pairs")
                delta = complete[levels[1]]-complete[levels[0]]
                se = float(delta.std(ddof=1)/np.sqrt(len(delta)))
                if se <= 1e-14:
                    raise ToolValidationError("Paired changes have zero variance; inference is not estimable")
                estimate, ci = float(delta.mean()), float(t.ppf(.975, len(delta)-1)*se)
                rows.append({**labels, "readout": value, "estimate": estimate, "ci_lower": estimate-ci,
                             "ci_upper": estimate+ci, "standard_error": se, "p_value": float(2*t.sf(abs(estimate/se), len(delta)-1)),
                             "n_donors": len(delta), "contrast": f"{levels[1]} - {levels[0]}", "method": "paired_t"})
        result = pd.DataFrame(rows)
        result["q_value"] = _bh(result.p_value)
        return _output(r, result, excluded_pairs=dropped)


@register_tool("AgingCompositionTool")
class AgingCompositionTool(AgingTool):
    operations = ("summarize_proportions", "test_composition_changes", "assess_composition_sensitivity")

    def summarize_proportions(self, r):
        import numpy as np
        df = read_table(r["input_path"])
        donor, cell = r.get("donor_key", "donor_id"), r.get("cell_type_key", "cell_type")
        keys = [donor, *r.get("strata", [])]
        _labels(df, [*keys, cell])
        if r.get("count_key"):
            count = r["count_key"]
            _numeric(df, [count])
            if (df[count] < 0).any() or (df[count] != np.floor(df[count])).any():
                raise ToolValidationError("Cell counts must be nonnegative integers")
        else:
            _unique(df, [r.get("observation_key", "cell_id")])
            count = "__count"
            df[count] = 1
        metadata = [k for k in r.get("metadata_columns", []) if k not in keys]
        if set(metadata) & {cell, "cell_count", "total_cells", "proportion"}:
            raise ToolValidationError("Metadata collides with composition output columns")
        _labels(df, metadata)
        groups = df.groupby(keys, observed=True, dropna=False)
        for key in metadata:
            if groups[key].nunique().gt(1).any():
                raise ToolValidationError(f"Metadata {key} varies within composition unit")
        counts = df.groupby([*keys, cell], observed=True)[count].sum().unstack(cell, fill_value=0)
        if counts.sum(axis=1).le(0).any():
            raise ToolValidationError("Composition unit has zero cells")
        long = counts.rename_axis(columns=cell).stack().rename("cell_count").reset_index()
        long["total_cells"] = long.groupby(keys, observed=True)["cell_count"].transform("sum")
        long["proportion"] = long.cell_count / long.total_cells
        if metadata:
            meta = groups[metadata].first().reset_index()
            long = long.merge(meta, on=keys, validate="many_to_one")
        return _output(r, long, denominator="All supplied cell types per donor/stratum; absent types are zero-filled")

    def test_composition_changes(self, r):
        df = read_table(r["input_path"])
        _numeric(df, [r.get("proportion_key", "proportion")])
        if not df[r.get("proportion_key", "proportion")].between(0, 1).all():
            raise ToolValidationError("Proportions must lie in [0,1]")
        request = dict(r, value_columns=[r.get("proportion_key", "proportion")])
        request["strata"] = list(dict.fromkeys([*r.get("strata", []), r.get("cell_type_key", "cell_type")]))
        return _output(r, _association(request, df), interpretation="Exploratory donor-level proportion differences; compositional dependence remains")

    def assess_composition_sensitivity(self, r):
        import pandas as pd
        df, composition = read_table(r["input_path"]), r["composition_columns"]
        _numeric(df, composition)
        if not df[composition].ge(0).all().all() or not df[composition].le(1).all().all():
            raise ToolValidationError("Composition covariates must be proportions in [0,1]")
        baseline = _association(r, df)
        baseline["adjustment"] = "baseline"
        request = dict(r, covariates=[*r.get("covariates", []), *composition])
        adjusted = _association(request, df)
        adjusted["adjustment"] = "composition_adjusted"
        return _output(r, pd.concat([baseline, adjusted], ignore_index=True), interpretation="Sensitivity analysis, not causal mediation; omit a reference proportion to avoid collinearity")


@register_tool("AgingValidationTool")
class AgingValidationTool(AgingTool):
    operations = ("evaluate_age_prediction", "validate_external_cohort", "assess_robustness", "evaluate_perturbation")

    def evaluate_age_prediction(self, r):
        import numpy as np
        import pandas as pd
        from scipy.stats import spearmanr
        df, strata = read_table(r["input_path"]), r.get("strata", [])
        _unique(df, [r.get("donor_key", "donor_id"), *strata])
        actual, predicted = r["age_key"], r["prediction_key"]
        _numeric(df, [actual, predicted])
        rows = []
        for labels, part in _strata(df, strata):
            error = part[predicted]-part[actual]
            rho = float(spearmanr(part[actual], part[predicted]).statistic) if len(part) >= 3 and part[actual].nunique() > 1 and part[predicted].nunique() > 1 else None
            rows.append({**labels, "n_donors": len(part), "mae": float(error.abs().mean()), "rmse": float(np.sqrt((error**2).mean())),
                         "mean_error": float(error.mean()), "spearman_rho": rho})
        return _output(r, pd.DataFrame(rows), interpretation="Chronological age prediction accuracy; not a validation of biological age or senescence")

    def validate_external_cohort(self, r):
        discovery, external = read_table(r["discovery_results_path"]), read_table(r["external_results_path"])
        keys = r.get("match_keys", ["readout"])
        for df in (discovery, external):
            _unique(df, keys)
            _numeric(df, ["estimate", "p_value"])
            if not df.p_value.between(0, 1).all():
                raise ToolValidationError("p_value must lie in [0,1]")
        for field in ("contrast", "predictor"):
            if field in discovery and field in external and field not in keys:
                check = discovery[keys+[field]].merge(external[keys+[field]], on=keys, suffixes=("_a", "_b"), validate="one_to_one")
                if not check[field+"_a"].eq(check[field+"_b"]).all():
                    raise ToolValidationError(f"External cohort {field} differs from discovery")
        result = discovery[keys+["estimate", "p_value"]].merge(external[keys+["estimate", "p_value"]], on=keys, suffixes=("_discovery", "_external"), how="inner", validate="one_to_one")
        if result.empty:
            raise ToolValidationError("No matching hypotheses between cohorts")
        result["q_value_external"] = _bh(result.p_value_external)
        result["same_direction"] = (result.estimate_discovery * result.estimate_external) > 0
        result["direction_and_fdr_supported"] = result.same_direction & (result.q_value_external <= r.get("alpha", .05))
        return _output(r, result, interpretation="Requires independent cohorts, identical outcome scales and contrasts; supplied results are not refitted", unmatched_discovery=len(discovery)-len(result), unmatched_external=len(external)-len(result))

    def assess_robustness(self, r):
        import pandas as pd
        df, donor = read_table(r["input_path"]), r.get("donor_key", "donor_id")
        base = _association(r, df)
        base["excluded_donor"] = "none"
        results, failures = [base], []
        for identity in df[donor].unique():
            try:
                one = _association(r, df[df[donor] != identity].copy())
                one["excluded_donor"] = str(identity)
                results.append(one)
            except ToolValidationError as exc:
                failures.append({"excluded_donor": str(identity), "reason": str(exc)})
        return _output(r, pd.concat(results, ignore_index=True), failed_refits=failures, interpretation="Leave-one-donor-out sensitivity; p values are descriptive and not independent replication")

    def evaluate_perturbation(self, r):
        return _output(r, _association(r), interpretation="Observed treatment-group association; causal interpretation requires an appropriate experimental design")
