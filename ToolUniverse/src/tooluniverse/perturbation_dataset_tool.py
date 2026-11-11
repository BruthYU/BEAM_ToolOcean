"""Local expression data loading and persistent references for perturbation workflows."""

import hashlib
import json
import re
import uuid
from pathlib import Path

from .base_tool import BaseTool
from .tool_registry import register_tool


class DatasetLoadError(ValueError):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


def _fingerprint(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise DatasetLoadError("INVALID_DATA", "Source changed while being read.")
    return {
        "path": str(path),
        "bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": digest.hexdigest(),
    }


def _align_metadata(obs, metadata, left, right):
    """Align unique metadata keys to observation order, supporting composite keys."""
    import pandas as pd

    lookup = metadata.set_index(right)
    index = (
        pd.MultiIndex.from_frame(obs[left]) if len(left) > 1 else pd.Index(obs[left[0]])
    )
    missing = int((~index.isin(lookup.index)).sum())
    return lookup.reindex(index), missing


def _mtx_members(source):
    """Require one unambiguous matrix, feature list and barcode list."""
    choices = [
        ["matrix.mtx", "matrix.mtx.gz"],
        ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"],
        ["barcodes.tsv", "barcodes.tsv.gz"],
    ]
    members = []
    for names in choices:
        found = [source / name for name in names if (source / name).is_file()]
        if len(found) != 1:
            raise DatasetLoadError(
                "INVALID_DATA", f"Expected exactly one of {names} in {source}."
            )
        members.append(found[0])
    return members


def _read_expression(arguments):
    """Return AnnData, detected format and fingerprints without dense conversion."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    source = Path(arguments["source"]).expanduser().resolve()
    if not source.exists():
        raise DatasetLoadError("SOURCE_NOT_FOUND", f"Source not found: {source}")
    fmt = arguments.get("format", "auto")
    if fmt == "auto":
        fmt = (
            "10x_mtx"
            if source.is_dir()
            else {".h5ad": "h5ad", ".h5": "10x_h5", ".csv": "csv", ".tsv": "tsv"}.get(
                source.suffix.lower()
            )
        )
        if source.name.lower().endswith((".csv.gz", ".tsv.gz")):
            fmt = "csv" if source.name.lower().endswith(".csv.gz") else "tsv"
    if fmt not in ("h5ad", "10x_h5", "10x_mtx", "csv", "tsv"):
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT", "Supported formats: h5ad, 10x_h5, 10x_mtx, csv, tsv."
        )
    if (fmt == "10x_mtx") != source.is_dir():
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "10x_mtx requires a directory; other formats require a file.",
        )
    if fmt != "h5ad" and (
        "backed" in arguments or arguments.get("matrix_key", "X") != "X"
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT", "backed and non-X matrix selection apply only to H5AD."
        )
    if fmt not in ("csv", "tsv") and any(
        k in arguments for k in ("matrix_orientation", "index_column")
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "matrix_orientation and index_column apply only to CSV/TSV.",
        )
    paths = _mtx_members(source) if fmt == "10x_mtx" else [source]
    fingerprints = [_fingerprint(path) for path in paths]
    if fmt == "h5ad":
        return (
            ad.read_h5ad(source, backed="r" if arguments.get("backed", True) else None),
            fmt,
            fingerprints,
        )
    if fmt == "10x_h5":
        import h5py
        from scipy.sparse import csc_matrix

        with h5py.File(source, "r") as handle:
            if "matrix" not in handle:
                raise DatasetLoadError(
                    "UNSUPPORTED_FORMAT",
                    "Expected modern 10x H5 /matrix layout; legacy or generic HDF5 is unsupported.",
                )
            group = handle["matrix"]
            required = [
                "data",
                "indices",
                "indptr",
                "shape",
                "barcodes",
                "features/id",
                "features/name",
            ]
            if any(key not in group for key in required):
                raise DatasetLoadError(
                    "INVALID_DATA", "Incomplete 10x H5 matrix layout."
                )
            shape = tuple(int(v) for v in group["shape"][:])
            if len(shape) != 2 or min(shape) <= 0:
                raise DatasetLoadError("INVALID_DATA", "Invalid 10x matrix dimensions.")
            matrix = csc_matrix(
                (group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape
            )
            matrix.check_format(full_check=True)
            barcodes = group["barcodes"][:].astype(str)
            genes = group["features/id"][:].astype(str)
            symbols = group["features/name"][:].astype(str)
            variables = pd.DataFrame({"gene_symbols": symbols}, index=genes)
            if "feature_type" in group["features"]:
                variables["feature_type"] = group["features/feature_type"][:].astype(
                    str
                )
    elif fmt == "10x_mtx":
        from scipy.io import mmread
        from scipy.sparse import issparse

        matrix = mmread(paths[0])
        if not issparse(matrix):
            raise DatasetLoadError(
                "INVALID_DATA", "10x MTX must be a sparse coordinate matrix."
            )
        matrix = matrix.tocsc()
        matrix.check_format(full_check=True)
        features = pd.read_csv(
            paths[1], sep="\t", header=None, dtype=str, keep_default_na=False
        )
        if features.shape[1] not in (2, 3):
            raise DatasetLoadError(
                "INVALID_DATA", "10x features file must have two or three columns."
            )
        barcode_table = pd.read_csv(
            paths[2], sep="\t", header=None, dtype=str, keep_default_na=False
        )
        if barcode_table.shape[1] != 1:
            raise DatasetLoadError(
                "INVALID_DATA", "10x barcode file must have one column."
            )
        barcodes, genes = barcode_table[0].to_numpy(), features[0].to_numpy()
        variables = pd.DataFrame({"gene_symbols": features[1].to_numpy()}, index=genes)
        if features.shape[1] == 3:
            variables["feature_type"] = features[2].to_numpy()
    else:
        import csv

        if "matrix_orientation" not in arguments or "index_column" not in arguments:
            raise DatasetLoadError(
                "INVALID_ARGUMENT",
                "CSV/TSV requires matrix_orientation and index_column.",
            )
        separator = "\t" if fmt == "tsv" else ","
        import gzip

        opener = gzip.open if source.suffix.lower() == ".gz" else open
        with opener(source, "rt", newline="", encoding="utf-8-sig") as stream:
            header = next(csv.reader(stream, delimiter=separator), [])
        if (
            not header
            or any(not value for value in header)
            or len(header) != len(set(header))
        ):
            raise DatasetLoadError("INVALID_DATA", "Table column names must be unique.")
        index_column = arguments["index_column"]
        if index_column not in header:
            raise DatasetLoadError(
                "INVALID_ARGUMENT", "index_column is not in the table header."
            )
        frame = pd.read_csv(
            source,
            sep=separator,
            index_col=index_column,
            dtype={index_column: str},
            keep_default_na=False,
        )
        values = frame.apply(pd.to_numeric, errors="raise").to_numpy()
        if not np.isfinite(values).all():
            raise DatasetLoadError(
                "INVALID_DATA", "Expression table must contain finite numeric values."
            )
        orientation = arguments["matrix_orientation"]
        if orientation == "features_by_observations":
            frame = frame.T
            values = values.T
        elif orientation != "observations_by_features":
            raise DatasetLoadError("INVALID_ARGUMENT", "Invalid matrix_orientation.")
        data = ad.AnnData(
            X=values,
            obs=pd.DataFrame(index=frame.index),
            var=pd.DataFrame(index=frame.columns),
        )
        return data, fmt, fingerprints
    if matrix.shape != (len(genes), len(barcodes)):
        raise DatasetLoadError(
            "INVALID_DATA", "Matrix dimensions do not match feature and barcode files."
        )
    if matrix.data.size and (
        not np.isfinite(matrix.data).all() or (matrix.data < 0).any()
    ):
        raise DatasetLoadError(
            "INVALID_DATA", "10x counts must be finite and non-negative."
        )
    if (
        "feature_type" in variables
        and (variables["feature_type"] != "Gene Expression").any()
    ):
        raise DatasetLoadError(
            "UNSUPPORTED_FORMAT",
            "Mixed feature types require explicit modality selection, not yet supported.",
        )
    return (
        ad.AnnData(X=matrix.T.tocsr(), obs=pd.DataFrame(index=barcodes), var=variables),
        fmt,
        fingerprints,
    )


def _inspect(arguments):
    adata, fmt, inputs = _read_expression(arguments)
    matrix_key = arguments.get("matrix_key", "X")
    try:
        if matrix_key == "X":
            matrix, variables = adata.X, adata.var
        elif matrix_key == "raw/X" and adata.raw is not None:
            matrix, variables = adata.raw.X, adata.raw.var
        elif matrix_key.startswith("layers/") and matrix_key[7:] in adata.layers:
            matrix, variables = adata.layers[matrix_key[7:]], adata.var
        else:
            raise DatasetLoadError(
                "INVALID_ARGUMENT", f"Matrix does not exist: {matrix_key}"
            )
        if matrix is None or min(matrix.shape) == 0:
            raise DatasetLoadError(
                "INVALID_DATA", "Selected expression matrix is empty."
            )
        gene_column = arguments.get("gene_id_column")
        if gene_column and gene_column not in variables:
            raise DatasetLoadError(
                "METADATA_MISMATCH", f"Missing gene ID column: {gene_column}"
            )
        genes = variables[gene_column] if gene_column else variables.index
        if (
            not adata.obs_names.is_unique
            or not genes.is_unique
            or bool(genes.isna().any())
            or adata.obs_names.isna().any()
            or (genes.astype(str) == "").any()
            or (adata.obs_names == "").any()
        ):
            raise DatasetLoadError(
                "INVALID_DATA", "Observation and gene IDs must be unique and non-null."
            )
        obs = adata.obs.copy()
        # An explicit reserved key exposes the observation index for external joins.
        if "__obs_index__" in obs:
            raise DatasetLoadError(
                "INVALID_DATA", "Reserved column __obs_index__ already exists."
            )
        obs["__obs_index__"] = adata.obs_names
        warnings = []
        metadata_path = arguments.get("metadata_path")
        if metadata_path:
            import pandas as pd

            path = Path(metadata_path).expanduser().resolve()
            if not path.is_file():
                raise DatasetLoadError(
                    "SOURCE_NOT_FOUND", f"Metadata file not found: {path}"
                )
            inputs.append(_fingerprint(path))
            if path.suffix.lower() == ".parquet":
                metadata = pd.read_parquet(path)
            elif path.suffix.lower() in (".csv", ".tsv"):
                metadata = pd.read_csv(
                    path, sep="\t" if path.suffix.lower() == ".tsv" else ","
                )
            else:
                raise DatasetLoadError(
                    "UNSUPPORTED_FORMAT", "Metadata must be CSV, TSV or Parquet."
                )
            keys = arguments.get("metadata_join_keys", {})
            left, right = keys.get("data", []), keys.get("metadata", [])
            if not left or len(left) != len(right):
                raise DatasetLoadError(
                    "INVALID_ARGUMENT",
                    "Provide equal non-empty metadata join key lists.",
                )
            if any(key not in obs for key in left) or any(
                key not in metadata for key in right
            ):
                raise DatasetLoadError(
                    "METADATA_MISMATCH", "A metadata join key is missing."
                )
            if (
                metadata.duplicated(right).any()
                or metadata[right].isna().any().any()
                or obs[left].isna().any().any()
            ):
                raise DatasetLoadError(
                    "METADATA_MISMATCH",
                    "Metadata join keys must be non-null and metadata keys unique.",
                )
            payload = [col for col in metadata if col not in right]
            if set(payload) & set(obs.columns):
                raise DatasetLoadError(
                    "METADATA_MISMATCH",
                    "Metadata payload columns overlap observation columns.",
                )
            aligned, missing = _align_metadata(obs, metadata, left, right)
            if missing and arguments.get("metadata_match", "strict") == "strict":
                raise DatasetLoadError(
                    "METADATA_MISMATCH",
                    f"{missing} observations have no metadata match.",
                )
            for column in payload:
                obs[column] = aligned[column].to_numpy()
            if missing:
                warnings.append({"code": "UNMATCHED_METADATA", "count": missing})
        elif arguments.get("metadata_join_keys"):
            raise DatasetLoadError(
                "INVALID_ARGUMENT", "metadata_join_keys requires metadata_path."
            )
        mapping = arguments.get("obs_column_map", {})
        for column in mapping.values():
            if column not in obs:
                raise DatasetLoadError(
                    "METADATA_MISMATCH", f"Missing mapped column: {column}"
                )
        if "control_values" in arguments and "perturbation" not in mapping:
            raise DatasetLoadError(
                "INVALID_ARGUMENT",
                "control_values requires an explicit perturbation column mapping.",
            )
        if "control_values" in arguments:
            values = arguments["control_values"]
            count = int(obs[mapping["perturbation"]].isin(values).sum())
            if count == 0:
                raise DatasetLoadError(
                    "METADATA_MISMATCH", "No observations match control_values."
                )
        summary = {
            "modality": "expression",
            "format": fmt,
            "shape": [int(value) for value in matrix.shape],
            "dtype": str(matrix.dtype),
            "matrix_key": matrix_key,
            "orientation": "observations_by_features",
            "expression_scale": arguments.get("expression_scale", "unknown"),
            "metadata_fields": list(mapping),
        }
        if "control_values" in arguments:
            summary["n_controls"] = count
        if matrix_key.startswith("layers/") and arguments.get("backed", True):
            warnings.append(
                {
                    "code": "LAYER_MEMORY",
                    "message": "AnnData may materialize layers despite backed=True.",
                }
            )
        # Verify the source still matches the data that was inspected.
        for entry in inputs:
            if _fingerprint(Path(entry["path"])) != entry:
                raise DatasetLoadError(
                    "INVALID_DATA", "Input changed during inspection."
                )
        return summary, inputs, warnings
    finally:
        if adata.isbacked:
            adata.file.close()


def resolve_dataset_ref(dataset_ref):
    """Validate a persistent reference, including contents of all original inputs.

    Returns a descriptor with reproducible read options; does not open matrix handles.
    """
    descriptor = json.loads(Path(dataset_ref).read_text())
    if descriptor.get("schema_version") != 1:
        raise DatasetLoadError("INVALID_DATA", "Unsupported dataset reference version.")
    if descriptor.get("modality") != "expression" or not descriptor.get("inputs"):
        raise DatasetLoadError("INVALID_DATA", "Invalid expression dataset reference.")
    options = descriptor.get("read_options", {})
    if descriptor.get("format") == "10x_mtx":
        current = {str(path) for path in _mtx_members(Path(options["source"]))}
        recorded = {entry["path"] for entry in descriptor["inputs"][:3]}
        if current != recorded:
            raise DatasetLoadError(
                "INVALID_DATA", "MTX directory members changed; load again."
            )
    for entry in descriptor["inputs"]:
        path = Path(entry["path"])
        if not path.is_file():
            raise DatasetLoadError(
                "SOURCE_NOT_FOUND", f"Referenced input is missing: {path}"
            )
        if _fingerprint(path) != entry:
            raise DatasetLoadError(
                "INVALID_DATA", "Referenced input changed; load the dataset again."
            )
    return descriptor


@register_tool("PerturbationTool")
class PerturbationTool(BaseTool):
    def run(self, arguments):
        directory = None
        try:
            dataset_id = arguments.get("dataset_id", uuid.uuid4().hex)
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,127}", dataset_id):
                raise DatasetLoadError(
                    "INVALID_ARGUMENT",
                    "dataset_id must contain 1-128 letters, digits, underscores or hyphens, starting with a letter or digit.",
                )
            summary, inputs, warnings = _inspect(arguments)
            destination = (
                Path(arguments["output_dir"]).expanduser().resolve() / dataset_id
            )
            destination.parent.mkdir(parents=True, exist_ok=True)
            try:
                destination.mkdir()
            except FileExistsError:
                raise DatasetLoadError(
                    "DATASET_EXISTS", f"Dataset already exists: {destination}"
                ) from None
            directory = destination
            options = dict(arguments)
            options["source"] = str(Path(arguments["source"]).expanduser().resolve())
            options["format"] = summary["format"]
            if options.get("metadata_path"):
                options["metadata_path"] = str(
                    Path(arguments["metadata_path"]).expanduser().resolve()
                )
            descriptor = {
                "schema_version": 1,
                "dataset_id": dataset_id,
                **summary,
                "inputs": inputs,
                "read_options": options,
                "warnings": warnings,
            }
            temporary = directory / "manifest.json.part"
            temporary.write_text(json.dumps(descriptor, indent=2) + "\n")
            manifest = directory / "manifest.json"
            temporary.replace(manifest)
            return {
                "status": "success",
                "data": {
                    "dataset_ref": str(manifest),
                    "dataset_id": dataset_id,
                    **summary,
                },
                "warnings": warnings,
            }
        except DatasetLoadError as exc:
            return {"status": "error", "error_code": exc.code, "error": str(exc)}
        except ImportError as exc:
            return {
                "status": "error",
                "error_code": "DEPENDENCY_MISSING",
                "error": str(exc),
            }
        except OSError as exc:
            return {"status": "error", "error_code": "IO_ERROR", "error": str(exc)}
        except (ValueError, KeyError, TypeError) as exc:
            return {"status": "error", "error_code": "INVALID_DATA", "error": str(exc)}
        finally:
            if directory is not None and not (directory / "manifest.json").exists():
                temporary = directory / "manifest.json.part"
                temporary.unlink(missing_ok=True)
                directory.rmdir()
