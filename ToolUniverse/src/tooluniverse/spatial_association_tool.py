"""Patient-level spatial association of two abundance profiles."""

from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path

from .exceptions import ToolValidationError
from .spatial_transcriptomics_data_tool import _coordinates, _open_h5ad
from .spatial_transcriptomics_tool import _SpatialTool, _result_directory
from .tool_registry import register_tool


def _input_table(request):
    import numpy as np
    import pandas as pd

    path = Path(request["input_path"]).expanduser()
    fields = {
        "patient": request["patient_key"], "group": request["group_key"],
        "domain": request["domain_key"],
    }
    kinds = request["cell_types"]
    columns = ["a", "b"] if len(kinds) == 2 else [f"feature_{i}" for i in range(len(kinds))]
    slice_key = request.get("slice_key")
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        if "x_key" in request or "y_key" in request:
            raise ToolValidationError("x_key and y_key are only supported for CSV inputs.")
        with _open_h5ad(path) as adata:
            if not adata.obs_names.is_unique:
                raise ToolValidationError("The h5ad must have unique observation IDs.")
            for key in fields.values():
                if key not in adata.obs:
                    raise ToolValidationError(f"Metadata column {key!r} was not found in obs.")
            if slice_key is None:
                slice_key = next((key for key in ("slice_id", "library_id") if key in adata.obs), None)
            if slice_key is not None and slice_key not in adata.obs:
                raise ToolValidationError(f"Slice column {slice_key!r} was not found in obs.")
            table = pd.DataFrame({name: adata.obs[key] for name, key in fields.items()})
            table["slice"] = "1" if slice_key is None else adata.obs[slice_key]
            table[["x", "y"]] = _coordinates(adata, request.get("coordinate_key", "spatial"))
            abundance_key = request.get("abundance_key")
            abundance = adata.obs if abundance_key is None else adata.obsm.get(abundance_key)
            if not isinstance(abundance, pd.DataFrame) or not abundance.columns.is_unique:
                raise ToolValidationError("Abundances must be named columns in obs or in a DataFrame in obsm.")
            for column, kind in zip(columns, kinds):
                if kind not in abundance:
                    raise ToolValidationError(f"Abundance column {kind!r} was not found.")
                table[column] = abundance[kind]
    elif suffix == ".csv":
        if "abundance_key" in request or "coordinate_key" in request:
            raise ToolValidationError("abundance_key and coordinate_key are only supported for h5ad inputs.")
        try:
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        except (OSError, ValueError) as exc:
            raise ToolValidationError(f"Cannot read the input CSV: {exc}") from exc
        if slice_key is None:
            slice_key = next((key for key in ("slice_id", "library_id") if key in frame), None)
        if slice_key is not None:
            fields["slice"] = slice_key
        fields.update(
            x=request.get("x_key", "x"), y=request.get("y_key", "y"),
        )
        fields.update(zip(columns, kinds))
        missing = [key for key in fields.values() if key not in frame]
        if missing:
            raise ToolValidationError(f"CSV columns were not found: {missing}.")
        table = pd.DataFrame({name: frame[key] for name, key in fields.items()})
        if slice_key is None:
            table["slice"] = "1"
    else:
        raise ToolValidationError("The input file must end in .h5ad or .csv.")

    if table.empty:
        raise ToolValidationError("The input contains no spots.")
    for key in ("patient", "group", "slice", "domain"):
        values = table[key].astype("string").str.strip()
        if values.isna().any() or values.eq("").any():
            raise ToolValidationError(f"The {key} field contains missing or empty values.")
        table[key] = values
    if table.groupby("patient")["group"].nunique().gt(1).any():
        raise ToolValidationError("Each patient must belong to exactly one group.")
    for key in ["x", "y", *columns]:
        table[key] = pd.to_numeric(table[key], errors="coerce")
    if not np.isfinite(table[["x", "y", *columns]]).to_numpy().all():
        raise ToolValidationError("Coordinates and abundances must contain finite numbers.")
    if table[columns].lt(0).to_numpy().any():
        raise ToolValidationError("Cell abundances must be nonnegative.")
    return table


def _neighborhood(xy, n_neighbors, max_distance):
    import numpy as np
    from scipy.sparse import csr_matrix, diags
    from scipy.spatial import cKDTree

    if len(np.unique(xy, axis=0)) != len(xy):
        raise ToolValidationError("Coordinates must be unique within each patient and slice for neighborhood analysis.")
    n = len(xy)
    k = min(n_neighbors, n - 1)
    distances, neighbors = cKDTree(xy).query(xy, k=k + 1)
    sources = np.broadcast_to(np.arange(n)[:, None], neighbors.shape)
    mask = sources != neighbors
    if max_distance is not None:
        mask &= distances <= max_distance
    graph = csr_matrix((np.ones(int(mask.sum())), (sources[mask], neighbors[mask])), shape=(n, n))
    graph = graph.maximum(graph.T)
    degree = np.asarray(graph.sum(axis=1)).ravel()
    inverse = np.divide(1.0, np.sqrt(degree), out=np.zeros_like(degree), where=degree > 0)
    diagonal = diags(inverse)
    return diagonal @ graph @ diagonal


def _slice_score(frame, request, rng):
    import numpy as np

    n = len(frame)
    if n < 3:
        return None, None, "fewer_than_three_spots"
    vectors = frame[["a", "b"]].to_numpy(dtype=float, copy=True)
    # Rescale first to avoid overflow when centering high abundance values.
    maximum = vectors.max(axis=0)
    if (maximum == 0).any():
        return None, None, "constant_abundance"
    vectors /= maximum
    vectors -= vectors.mean(axis=0)
    norms = np.linalg.norm(vectors, axis=0)
    if (norms <= np.finfo(float).eps).any():
        return None, None, "constant_abundance"
    a, b = (vectors / norms).T
    mode = request.get("method", "neighborhood")
    if mode == "neighborhood":
        graph = _neighborhood(
            frame[["x", "y"]].to_numpy(dtype=float), request.get("n_neighbors", 6),
            request.get("max_neighbor_distance"),
        )
        if not graph.nnz:
            return None, None, "no_neighbor_edges"
        a = graph @ a

    observed = float(np.clip(a @ b, -1, 1))
    blocks = list(frame.groupby("domain", sort=True).indices.values())
    if not any(len(block) > 1 and np.ptp(b[block]) > 0 for block in blocks):
        return observed, None, "no_permutable_variation"
    # The graph and the A profile stay fixed. Shuffle B only within each domain.
    permutations = request.get("n_permutations", 999)
    simulated = np.empty(permutations)
    for index in range(permutations):
        shuffled = b.copy()
        for block in blocks:
            shuffled[block] = rng.permutation(b[block])
        simulated[index] = np.clip(a @ shuffled, -1, 1)
    return observed, simulated, "ok"


def _patient_comparison(patients, groups, n_permutations, rng):
    import numpy as np

    base = {"unit": "patient", "groups": groups, "difference": None, "p_value": None, "n_permutations": 0}
    counts = {group: sum(row["group"] == group and row["score"] is not None for row in patients)
              for group in groups}
    base["n_patients_by_group"] = counts
    if len(groups) < 2 or min(counts.values()) < 2:
        return {**base, "status": "insufficient_patients", "permutation_method": None}
    left, right = (
        np.array([row["score"] for row in patients if row["group"] == group and row["score"] is not None])
        for group in groups
    )
    scores = np.r_[left, right]
    total = scores.sum()
    observed = float(left.mean() - right.mean())
    exact = comb(len(scores), len(left)) <= n_permutations
    assignments = (
        combinations(range(len(scores)), len(left)) if exact else
        (rng.choice(len(scores), len(left), replace=False) for _ in range(n_permutations))
    )
    extreme = tested = 0
    for indices in assignments:
        subtotal = scores[list(indices)].sum()
        delta = subtotal / len(left) - (total - subtotal) / len(right)
        extreme += abs(delta) >= abs(observed) - 1e-12
        tested += 1
    p_value = extreme / tested if exact else (extreme + 1) / (tested + 1)
    return {
        **base, "status": "ok", "difference": observed, "p_value": float(p_value),
        "n_permutations": tested, "permutation_method": "exact" if exact else "monte_carlo",
    }


def _plot_patients(patients, groups, path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(6.4, 3.8), layout="constrained")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    for index, group in enumerate(groups):
        values = sorted(row["score"] for row in patients if row["group"] == group and row["score"] is not None)
        if values:
            offsets = [(point - (len(values) - 1) / 2) * 0.3 / max(1, len(values)) for point in range(len(values))]
            axis.scatter([index + offset for offset in offsets], values, color="#207887")
            axis.hlines(sum(values) / len(values), index - 0.22, index + 0.22, color="#202020")
    axis.axhline(0, color="#aaaaaa", linewidth=0.8, linestyle="--")
    axis.set(xticks=range(len(groups)), xticklabels=groups, ylabel="Patient association score",
             title="Spatial association by group")
    figure.savefig(path, dpi=150)


@register_tool("SpatialAssociationTool")
class SpatialAssociationTool(_SpatialTool):
    operations = ("analyze_spatial_association", "screen_spatial_associations")

    def analyze_spatial_association(self, request):
        import numpy as np
        import pandas as pd

        max_distance = request.get("max_neighbor_distance")
        if max_distance is not None and not isfinite(max_distance):
            raise ToolValidationError("max_neighbor_distance must be finite.")
        table = _input_table(request)
        groups = request.get("compare_groups")
        available = sorted(table["group"].unique().tolist())
        if groups is None:
            if len(available) > 2:
                raise ToolValidationError("Specify compare_groups when the input has more than two groups.")
            groups = available
        elif not set(groups).issubset(available):
            raise ToolValidationError(f"Unknown comparison groups. Available groups: {available}.")
        table = table[table["group"].isin(groups)]
        if "domains" in request:
            missing = sorted(set(request["domains"]) - set(table["domain"]))
            if missing:
                raise ToolValidationError(f"Requested domains were not found in the selected groups: {missing}.")
        rng = np.random.default_rng(request.get("random_seed", 0))
        n_permutations = request.get("n_permutations", 999)
        patient_rows, slice_rows = [], []
        for (patient, group), patient_table in table.groupby(["patient", "group"], sort=True):
            scores, nulls = [], []
            n_spots = 0
            has_permutations = False
            for slice_id, frame in patient_table.groupby("slice", sort=True):
                if "domains" in request:
                    frame = frame[frame["domain"].isin(request["domains"])]
                score, null, status = _slice_score(frame, request, rng)
                slice_rows.append({
                    "patient_id": patient, "group": group, "slice_id": slice_id,
                    "n_spots": len(frame), "score": score, "status": status,
                })
                if score is not None:
                    scores.append(score)
                    nulls.append(np.full(n_permutations, score) if null is None else null)
                    has_permutations |= null is not None
                    n_spots += len(frame)
            observed = float(np.mean(scores)) if scores else None
            p_value = None
            if has_permutations:
                combined = np.mean(nulls, axis=0)
                lower = int((combined <= observed + 1e-12).sum())
                upper = int((combined >= observed - 1e-12).sum())
                p_value = min(1.0, 2 * (min(lower, upper) + 1) / (n_permutations + 1))
            patient_rows.append({
                "patient_id": patient, "group": group, "score": observed,
                "spatial_p_value": p_value, "n_slices_used": len(scores), "n_spots_used": n_spots,
            })
        comparison = _patient_comparison(patient_rows, groups, n_permutations, rng)
        with _result_directory(request["output_dir"], "spatial_association") as output:
            patient_path, slice_path = output / "patients.csv", output / "slices.csv"
            comparison_path, figure_path = output / "group_comparison.json", output / "patients.png"
            pd.DataFrame(patient_rows).to_csv(patient_path, index=False)
            pd.DataFrame(slice_rows).to_csv(slice_path, index=False)
            comparison_path.write_text(json.dumps(comparison, indent=2) + "\n", encoding="utf-8")
            _plot_patients(patient_rows, groups, figure_path)
        return {
            "method": request.get("method", "neighborhood"), "cell_types": request["cell_types"],
            "patient_summary_path": str(patient_path), "slice_summary_path": str(slice_path),
            "group_comparison_path": str(comparison_path), "figure_path": str(figure_path),
            "group_comparison": comparison,
        }

    def screen_spatial_associations(self, request):
        """Screen a declared family of pairs with patient-level inference."""
        import numpy as np
        import pandas as pd
        from scipy.stats import false_discovery_control

        datasets = request["datasets"]
        identifiers = [item["dataset_id"] for item in datasets]
        if len(set(identifiers)) != len(identifiers):
            raise ToolValidationError("dataset_id values must be unique.")
        kinds = request["cell_types"]
        columns = ["a", "b"] if len(kinds) == 2 else [f"feature_{i}" for i in range(len(kinds))]
        groups = request["compare_groups"]
        rng = np.random.default_rng(request.get("random_seed", 0))
        rows, patients, slices = [], [], []
        known_groups = {}
        for dataset in datasets:
            table = _input_table({**request, **dataset})
            if not set(groups).issubset(set(table["group"])):
                raise ToolValidationError("Both compare_groups must be present in every dataset.")
            for patient, group in table.groupby("patient")["group"].first().items():
                if patient in known_groups and known_groups[patient] != group:
                    raise ToolValidationError("Patient group labels must agree across datasets.")
                known_groups[patient] = group
            table = table[table["group"].isin(groups)]
            if "domains" in request and not set(request["domains"]).issubset(set(table["domain"])):
                raise ToolValidationError("Requested domains were not found in the selected groups.")
            for i, j in combinations(range(len(kinds)), 2):
                identity = {"dataset_id": dataset["dataset_id"], "cell_a": kinds[i], "cell_b": kinds[j]}
                pair_patients = []
                for (patient, group), frame in table.groupby(["patient", "group"], sort=True):
                    scores = []
                    for slice_id, section in frame.groupby("slice", sort=True):
                        if "domains" in request:
                            section = section[section["domain"].isin(request["domains"])]
                        section = section[["domain", columns[i], columns[j]]].rename(
                            columns={columns[i]: "a", columns[j]: "b"})
                        # Zero spatial shuffles: inference here permutes patient labels only.
                        score, _, status = _slice_score(section, {"method": "colocalization", "n_permutations": 0}, rng)
                        slices.append({**identity, "patient_id": patient, "group": group,
                                       "slice_id": slice_id, "n_spots": len(section), "score": score, "status": status})
                        if score is not None:
                            scores.append(score)
                    row = {**identity, "patient_id": patient, "group": group,
                           "score": float(np.mean(scores)) if scores else None, "n_slices_used": len(scores)}
                    pair_patients.append(row)
                    patients.append(row)
                comparison = _patient_comparison(pair_patients, groups, request.get("n_permutations", 9999), rng)
                valid = [row for row in pair_patients if row["score"] is not None]
                loo = []
                for omitted in range(len(valid)):
                    remaining = valid[:omitted] + valid[omitted + 1:]
                    values = [[row["score"] for row in remaining if row["group"] == group] for group in groups]
                    if all(values):
                        loo.append(float(np.mean(values[0]) - np.mean(values[1])))
                rows.append({**identity, **comparison, "n_patients_excluded": len(pair_patients) - len(valid),
                             "loo_min": min(loo) if loo else None, "loo_max": max(loo) if loo else None})
        # Untestable pairs remain in the declared family, conservatively assigned p=1 for adjustment.
        adjusted = false_discovery_control(
            [row["p_value"] if row["p_value"] is not None else 1.0 for row in rows],
            method=request.get("fdr_method", "bh"))
        for row, q in zip(rows, adjusted):
            row["q_value"] = float(q) if row["p_value"] is not None else None
        ranked = sorted(rows, key=lambda row: (row["p_value"] if row["p_value"] is not None else 2,
                                               -abs(row["difference"] or 0), row["dataset_id"], row["cell_a"], row["cell_b"]))
        summary = {
            "method": "colocalization", "unit": "patient", "compare_groups": groups,
            "number_tests": len(rows), "number_testable": sum(row["p_value"] is not None for row in rows),
            "fdr_method": request.get("fdr_method", "bh"), "fdr_threshold": request.get("fdr_threshold", 0.05),
            "number_significant": sum(row["q_value"] is not None and row["q_value"] <= request.get("fdr_threshold", 0.05) for row in rows),
            "top_pairs": ranked[:request.get("top_n", 20)],
            "interpretation": "Exploratory same-spot Pearson association, not cell contact or causality. Equal slice means per patient; difference is first group minus second. FDR covers all pairs across all supplied datasets. Leave-one-patient-out ranges are descriptive, not confidence intervals.",
        }
        with _result_directory(request["output_dir"], "spatial_screen") as output:
            paths = {"screen_table_path": output / "screen.csv", "patient_summary_path": output / "patients.csv",
                     "slice_summary_path": output / "slices.csv", "summary_path": output / "summary.json"}
            pd.DataFrame(ranked).to_csv(paths["screen_table_path"], index=False)
            pd.DataFrame(patients).to_csv(paths["patient_summary_path"], index=False)
            pd.DataFrame(slices).to_csv(paths["slice_summary_path"], index=False)
            paths["summary_path"].write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        return {**summary, **{key: str(value) for key, value in paths.items()}}
