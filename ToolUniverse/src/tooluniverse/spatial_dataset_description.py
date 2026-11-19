"""Describe existing spatial datasets without inferring biological identities."""

import csv
from pathlib import Path

from .exceptions import ToolValidationError

_ALIASES = {
    "patient_id": ("patient_id", "patient", "donor_id", "donor"),
    "sample_id": ("sample_id", "sample", "slice_id", "library_id"),
    "group": ("group", "response", "condition"),
    "cell_type": ("cell_type", "cell_type_label", "celltype"),
    "domain": ("domain", "author_region", "region"),
    "species": ("species", "organism"),
    "disease": ("disease",),
    "tissue": ("tissue", "organ"),
    "treatment": ("treatment",),
    "sampling_time": ("sampling_time", "timepoint"),
    "source": ("source", "accession"),
    "dataset_path": ("dataset_path",),
}
_BIOLOGY = ("species", "disease", "tissue", "treatment", "sampling_time", "source")


def _columns(frame, mapping, source):
    selected = {}
    for meaning, column in mapping.items():
        if column not in frame:
            raise ToolValidationError(
                f"Column {column!r} was not found in {source}. Update metadata_columns to use an existing column."
            )
        selected[meaning] = column
    for meaning, aliases in _ALIASES.items():
        if meaning in selected:
            continue
        # Prefer the canonical spelling, then the documented alias order.
        selected[meaning] = next((key for key in aliases if key in frame), None)
    return {meaning: column for meaning, column in selected.items() if column is not None}


def _values(series):
    values = series.astype("string").str.strip().dropna()
    return values[values.ne("")]


def _describe_values(series, meaning, source, limit):
    values = _values(series)
    unique = list(dict.fromkeys(values.tolist()))
    if not unique:
        return f"{meaning.replace('_', ' ').capitalize()}: Not recorded in {source}."
    shown = ", ".join(repr(value) for value in unique[:limit])
    remainder = f" Only the first {limit} distinct values are shown." if len(unique) > limit else ""
    missing = len(series) - len(values)
    coverage = f" {missing} rows have no value." if missing else ""
    noun = "value" if len(unique) == 1 else "values"
    return (
        f"{meaning.replace('_', ' ').capitalize()} recorded in {source}: {shown} "
        f"({len(unique)} distinct {noun}).{remainder}{coverage}"
    )


def _background(request, obs, obs_columns, input_path):
    import pandas as pd

    path = Path(request["metadata_path"]).expanduser().resolve()
    try:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    except (OSError, ValueError) as exc:
        raise ToolValidationError(f"Cannot read the background CSV: {exc}") from exc
    if frame.empty:
        raise ToolValidationError("The background CSV is empty. Supply at least one annotated sample.")
    columns = _columns(frame, request.get("metadata_columns", {}).get("background", {}), "the background CSV")
    if "dataset_path" in columns:
        key = "dataset_path"
        raw = frame[columns[key]]
        if raw.str.strip().eq("").any():
            raise ToolValidationError("Background dataset paths must not be empty.")
        identities = raw.map(lambda value: str((path.parent / Path(value).expanduser()).resolve()))
        wanted = {str(input_path)}
    else:
        key = next((key for key in ("sample_id", "patient_id") if key in columns and key in obs_columns), None)
        if key is None:
            raise ToolValidationError(
                "The background cannot be linked to this dataset. Provide a matching sample_id or patient_id column "
                "in both files, or a dataset_path column in the background CSV."
            )
        identities = frame[columns[key]].astype("string").str.strip()
        if identities.eq("").any():
            raise ToolValidationError("Background sample or patient identifiers must not be empty.")
        wanted = set(_values(obs[obs_columns[key]]))
        if not wanted:
            raise ToolValidationError("The dataset has no recorded identifiers for linking the background CSV.")
    if identities.duplicated().any():
        raise ToolValidationError(
            f"The background CSV contains repeated {key} values. Use one row per identifier or link by a more specific sample ID."
        )
    matched = identities.isin(wanted)
    if not matched.any():
        raise ToolValidationError("No background rows match this dataset. Check the identifiers or dataset paths.")
    missing = len(wanted - set(identities[matched]))
    message = (
        f"Background descriptions use {int(matched.sum())} matching rows from {str(path)!r}, linked by {key}. "
        f"{int((~matched).sum())} unrelated rows were excluded."
    )
    if missing:
        message += f" Background information is missing for {missing} dataset identifiers."
    return frame.loc[matched], columns, message


def describe_spatial_dataset(request):
    import numpy as np
    import pandas as pd

    from .spatial_transcriptomics_data_tool import _open_h5ad

    path = Path(request["input_path"]).expanduser().resolve()
    limit = request.get("max_values", 20)
    mapping = request.get("metadata_columns", {})
    if mapping.get("background") and "metadata_path" not in request:
        raise ToolValidationError("Provide metadata_path when mapping background CSV columns.")
    obs = pd.DataFrame()
    uns = {}
    contents = []
    cells = []
    if path.suffix.lower() == ".h5ad":
        with _open_h5ad(path) as adata:
            obs = adata.obs.copy()
            if not obs.index.is_unique:
                raise ToolValidationError("Observation identifiers must be unique before describing sample annotations.")
            obs_columns = _columns(obs, mapping.get("obs", {}), "obs")
            uns = {key: value for key, value in adata.uns.items() if key in _ALIASES and isinstance(value, (str, int, float))}
            summary = f"This AnnData file contains {adata.n_obs} observations and {adata.n_vars} gene or feature identifiers."
            contents.append("An expression matrix is present in X." if adata.X is not None else "X does not contain an expression matrix.")
            contents.append("Expression layers: " + (", ".join(adata.layers.keys()) or "none recorded") + ".")
            contents.append("Observation annotation columns: " + (", ".join(obs.columns) or "none recorded") + ".")
            for key in adata.obsm:
                contents.append(f"obsm[{key!r}] contains an observation-associated matrix with shape {adata.obsm[key].shape}.")
            if "spatial" in adata.obsm:
                contents.append("Spatial coordinates are recorded in obsm['spatial'].")
            if "cell_type" in obs_columns:
                cells.append(_describe_values(obs[obs_columns["cell_type"]], "cell_type",
                                               f"obs[{obs_columns['cell_type']!r}]", limit))
            contents.append("The names of expression layers describe stored fields; raw-count or normalization status has not been inferred.")
    elif path.suffix.lower() == ".csv":
        if mapping.get("obs"):
            raise ToolValidationError("obs column mappings apply to h5ad files. For a reference CSV, provide dataset-linked background metadata.")
        try:
            with path.open(encoding="utf-8-sig", newline="") as stream:
                header = next(csv.reader(stream), [])
            if len(header) < 2 or len(set(header[1:])) != len(header) - 1 or any(not value.strip() for value in header[1:]):
                raise ValueError("Expected a gene ID column followed by distinct, nonempty cell type columns.")
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            values = frame.iloc[:, 1:].astype(float).to_numpy()
            ids = frame.iloc[:, 0]
            if frame.empty or ids.duplicated().any() or ids.str.strip().eq("").any():
                raise ValueError("Gene identifiers must be nonempty and unique.")
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError("Reference expressions must be finite and nonnegative.")
        except (OSError, ValueError) as exc:
            raise ToolValidationError(f"Cannot describe this reference-signature CSV: {exc}") from exc
        obs_columns = {}
        summary = f"This reference-signature table contains {len(frame)} genes and {len(header) - 1} cell type columns."
        cells.append(_describe_values(pd.Series(header[1:]), "cell_type", "the reference CSV column names", limit))
        contents.append("Rows are gene identifiers and columns contain reference expression values for the named cell types.")
        contents.append("Training history and the original cells are not recorded in this CSV.")
    else:
        raise ToolValidationError("Provide an h5ad dataset or a reference-signature CSV.")

    samples = []
    for meaning in ("patient_id", "sample_id", "group", "domain"):
        if meaning in obs_columns:
            column = obs_columns[meaning]
            samples.append(_describe_values(obs[column], meaning, f"obs[{column!r}]", limit))
    if "patient_id" in obs_columns and "group" in obs_columns:
        patients = _values(obs[obs_columns["patient_id"]])
        groups = _values(obs[obs_columns["group"]])
        pairs = pd.DataFrame({"patient": patients, "group": groups}).dropna().drop_duplicates()
        if pairs.groupby("patient")["group"].nunique().gt(1).any():
            samples.append("Some patients have multiple recorded group labels; review these labels before a patient-group comparison.")
        else:
            counts = pairs.groupby("group")["patient"].nunique()
            group_description = ", ".join(
                f"{label!r}: {count}" for label, count in counts.iloc[:limit].items()
            )
            samples.append(f"Patients with recorded group labels: {group_description}.")
            if len(counts) > limit:
                samples.append(f"Only the first {limit} patient groups are shown.")

    background = None
    background_columns = {}
    if "metadata_path" in request:
        background, background_columns, message = _background(request, obs, obs_columns, path)
        samples.append(message)
        for meaning in ("patient_id", "sample_id", "group", "cell_type", "domain"):
            if meaning in background_columns:
                description = _describe_values(background[background_columns[meaning]], meaning,
                                               f"matched background CSV column {background_columns[meaning]!r}", limit)
                (cells if meaning == "cell_type" else samples).append(description)
    biology = []
    for meaning in _BIOLOGY:
        found = False
        if meaning in obs_columns:
            column = obs_columns[meaning]
            biology.append(_describe_values(obs[column], meaning, f"obs[{column!r}]", limit))
            found = True
        if meaning in uns:
            biology.append(f"{meaning.replace('_', ' ').capitalize()} recorded in uns[{meaning!r}]: {uns[meaning]!r}.")
            found = True
        if meaning in background_columns:
            column = background_columns[meaning]
            biology.append(_describe_values(background[column], meaning, f"matched background CSV column {column!r}", limit))
            found = True
        if not found:
            biology.append(f"{meaning.replace('_', ' ').capitalize()}: Not recorded in the supplied data.")

    return {
        "dataset_summary": summary,
        "biological_metadata": biology,
        "sample_annotations": samples or ["Sample and patient annotations: Not recorded in the supplied data."],
        "cell_annotations": cells or ["Cell type annotations: Not recorded in the supplied data."],
        "data_contents": contents,
    }
