"""Tabular metadata readers without experiment-specific interpretation."""

import csv
import gzip
from pathlib import Path

from .perturbation_dataset_tool import DatasetLoadError


def read_metadata_table(source, fmt="auto", sheet_name=None):
    import pandas as pd

    path = Path(source)
    if not path.is_file():
        raise DatasetLoadError("SOURCE_NOT_FOUND", f"Metadata source not found: {path}")
    if fmt == "auto":
        name = path.name.lower()
        fmt = next(
            (
                kind
                for extension, kind in [
                    (".csv.gz", "csv"),
                    (".tsv.gz", "tsv"),
                    (".csv", "csv"),
                    (".tsv", "tsv"),
                    (".xlsx", "xlsx"),
                    (".parquet", "parquet"),
                    (".h5ad", "h5ad_obs"),
                ]
                if name.endswith(extension)
            ),
            None,
        )
    if fmt != "xlsx" and sheet_name is not None:
        raise DatasetLoadError("INVALID_ARGUMENT", "sheet_name applies only to XLSX.")
    if fmt in ("csv", "tsv"):
        separator = "," if fmt == "csv" else "\t"
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt", newline="", encoding="utf-8-sig") as stream:
            header = next(csv.reader(stream, delimiter=separator), [])
        if (
            not header
            or any(not name for name in header)
            or len(set(header)) != len(header)
        ):
            raise DatasetLoadError(
                "INVALID_DATA", "Metadata headers must be unique and non-empty."
            )
        table = pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False)
    elif fmt == "xlsx":
        sheet = sheet_name if sheet_name is not None else 0
        header_frame = pd.read_excel(
            path, sheet_name=sheet, header=None, nrows=1, keep_default_na=False
        )
        if header_frame.empty:
            raise DatasetLoadError("INVALID_DATA", "Metadata spreadsheet is empty.")
        header = header_frame.iloc[0].tolist()
        if any(not str(name) for name in header) or len(set(header)) != len(header):
            raise DatasetLoadError(
                "INVALID_DATA", "Metadata headers must be unique and non-empty."
            )
        table = pd.read_excel(path, sheet_name=sheet, dtype=str, keep_default_na=False)
    elif fmt == "parquet":
        table = pd.read_parquet(path)
    elif fmt == "h5ad_obs":
        import h5py
        from anndata.io import read_elem

        with h5py.File(path, "r") as handle:
            if "obs" not in handle:
                raise DatasetLoadError(
                    "INVALID_DATA", "H5AD has no observation metadata."
                )
            table = read_elem(handle["obs"])
        if "__obs_index__" in table:
            raise DatasetLoadError(
                "INVALID_DATA", "Reserved __obs_index__ column already exists."
            )
        table = table.copy()
        table["__obs_index__"] = table.index
    else:
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT",
            "Supported metadata: CSV/TSV (gzip), XLSX, Parquet, H5AD obs.",
        )
    if not isinstance(table, pd.DataFrame) or table.empty:
        raise DatasetLoadError("INVALID_DATA", "Metadata table is empty.")
    if not table.columns.is_unique or any(
        not isinstance(name, str) or not name for name in table.columns
    ):
        raise DatasetLoadError(
            "INVALID_DATA", "Metadata column names must be unique non-empty strings."
        )
    return table, fmt
