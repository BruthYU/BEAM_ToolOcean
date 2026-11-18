"""Standardize experimental perturbation annotations independently of a dataset."""

import json
import re
import uuid
from pathlib import Path

from .base_tool import BaseTool
from .metadata_io import read_metadata_table
from .perturbation_dataset_tool import DatasetLoadError, _fingerprint
from .tool_registry import register_tool

LEVEL_KEYS = {
    "cell": ["cell_id"],
    "sample": ["sample_id"],
    "well": ["plate_id", "well_id"],
    "field": ["plate_id", "well_id", "site_id"],
}
FIELDS = [
    "perturbation",
    "perturbation_type",
    "cell_id",
    "sample_id",
    "plate_id",
    "well_id",
    "site_id",
    "cell_type",
    "dose",
    "dose_unit",
    "time",
    "time_unit",
    "batch",
    "replicate",
    "is_control",
]


def _standardize(table, arguments):
    import numpy as np
    import pandas as pd

    level = arguments["experiment_level"]
    if level not in LEVEL_KEYS:
        raise DatasetLoadError("INVALID_ARGUMENT", "Unknown experiment_level.")
    mapping = arguments["field_map"]
    constants = arguments.get("constants", {})
    keys = arguments["association_keys"]
    supplied = set(mapping) | set(constants)
    if set(mapping) & set(constants) or supplied - set(FIELDS):
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            "Invalid or overlapping standard fields in field_map/constants.",
        )
    required = set(LEVEL_KEYS[level]) | {"perturbation"}
    if (
        not required <= supplied
        or not required - {"perturbation"} <= set(keys)
        or not set(keys) <= supplied
    ):
        raise DatasetLoadError(
            "INVALID_ARGUMENT",
            f"{level} requires fields {sorted(required)} and level identifiers in association_keys.",
        )
    if len(keys) != len(set(keys)):
        raise DatasetLoadError("INVALID_ARGUMENT", "association_keys must be unique.")
    if set(mapping.values()) - set(table.columns):
        raise DatasetLoadError(
            "METADATA_MISMATCH", "A field_map source column does not exist."
        )
    standardized = pd.DataFrame(index=range(len(table)))
    for field, column in mapping.items():
        standardized[field] = table[column].to_numpy()
    for field, value in constants.items():
        standardized[field] = value
    warnings = []
    labels = standardized["perturbation"].astype("string")
    missing_labels = labels.isna() | labels.str.fullmatch(r"\s*", na=False)
    if (
        missing_labels.any()
        and arguments.get("missing_perturbation", "error") == "drop_rows"
    ):
        count = int(missing_labels.sum())
        standardized = standardized.loc[~missing_labels].reset_index(drop=True)
        warnings.append({"code": "DROPPED_UNANNOTATED_ROWS", "count": count})
        if standardized.empty:
            raise DatasetLoadError(
                "INVALID_DATA", "No annotated experiment records remain."
            )
    for field in standardized:
        values = standardized[field]
        if values.map(lambda v: not pd.api.types.is_scalar(v)).any():
            raise DatasetLoadError("INVALID_DATA", f"Non-scalar values in {field}.")
        text = values.astype("string")
        missing = values.isna() | text.str.fullmatch(r"\s*", na=False)
        if (field in keys or field in ("perturbation", "is_control")) and missing.any():
            raise DatasetLoadError(
                "INVALID_DATA", f"Required field {field} contains missing values."
            )
        if field in ("dose", "time"):
            if values.map(lambda v: isinstance(v, (bool, np.bool_))).any():
                raise DatasetLoadError(
                    "INVALID_DATA", f"{field} cannot contain booleans."
                )
            numeric = pd.to_numeric(values.mask(missing), errors="coerce").astype(
                "Float64"
            )
            bad = (~missing & numeric.isna()) | (
                numeric.notna()
                & (
                    ~np.isfinite(numeric.to_numpy(dtype=float, na_value=np.nan))
                    | (numeric < 0)
                )
            )
            if bad.fillna(False).any():
                raise DatasetLoadError(
                    "INVALID_DATA",
                    f"{field} must contain finite non-negative numbers or explicit missing values.",
                )
            standardized[field] = numeric
        elif field == "is_control":
            normalized = text.str.strip().str.lower()
            parsed = normalized.map(
                {
                    "true": True,
                    "false": False,
                    "1": True,
                    "0": False,
                    "yes": True,
                    "no": False,
                }
            )
            if parsed.isna().any():
                raise DatasetLoadError(
                    "INVALID_DATA",
                    "is_control must contain booleans, true/false, 1/0 or yes/no.",
                )
            standardized[field] = parsed.astype(bool)
        else:
            standardized[field] = text.mask(missing)
        count = int(missing.sum())
        if count and arguments.get("missing_optional", "preserve") == "error":
            raise DatasetLoadError(
                "INVALID_DATA", f"Missing optional field values: {field}"
            )
    if "control_values" in arguments:
        derived = standardized["perturbation"].isin(arguments["control_values"])
        if not derived.any():
            raise DatasetLoadError(
                "METADATA_MISMATCH", "No perturbation labels match control_values."
            )
        if "is_control" in standardized and not standardized["is_control"].equals(
            derived
        ):
            raise DatasetLoadError(
                "METADATA_MISMATCH",
                "Explicit is_control conflicts with control_values.",
            )
        standardized["is_control"] = derived
    if "is_control" not in standardized:
        warnings.append(
            {
                "code": "CONTROL_UNSPECIFIED",
                "message": "No control labels were declared; controls were not inferred.",
            }
        )
    before = len(standardized)
    if standardized.duplicated(keys).any():
        if arguments.get("duplicate_keys", "error") == "error":
            raise DatasetLoadError(
                "METADATA_MISMATCH",
                "Duplicate association keys. Choose the correct level or explicitly collapse identical annotations.",
            )
        unique = standardized.drop_duplicates()
        if unique.duplicated(keys).any():
            raise DatasetLoadError(
                "METADATA_MISMATCH",
                "Duplicate association keys have conflicting experiment annotations.",
            )
        standardized = unique.reset_index(drop=True)
        warnings.append(
            {"code": "COLLAPSED_IDENTICAL_ROWS", "count": before - len(standardized)}
        )
    # Report optional-value counts for the rows actually persisted, rather than
    # repeating counts for image fields collapsed to a single well annotation.
    warnings = [
        w
        for w in warnings
        if w["code"] not in ("MISSING_OPTIONAL_VALUES", "UNIT_UNSPECIFIED")
    ]
    for field in standardized:
        count = int(standardized[field].isna().sum())
        if count:
            warnings.append(
                {"code": "MISSING_OPTIONAL_VALUES", "field": field, "count": count}
            )
    for numeric, unit in [("dose", "dose_unit"), ("time", "time_unit")]:
        if numeric in standardized:
            unknown = standardized[numeric].notna() & (
                standardized[unit].isna() if unit in standardized else True
            )
            if unknown.any():
                warnings.append(
                    {
                        "code": "UNIT_UNSPECIFIED",
                        "field": numeric,
                        "count": int(unknown.sum()),
                    }
                )
    return standardized, warnings


def resolve_perturbation_metadata_ref(metadata_ref):
    descriptor = json.loads(Path(metadata_ref).read_text())
    if (
        descriptor.get("schema_version") != 1
        or descriptor.get("modality") != "perturbation_metadata"
    ):
        raise DatasetLoadError(
            "INVALID_DATA", "Unsupported perturbation metadata reference."
        )
    for stamp in [descriptor["source"], descriptor["table"]]:
        path = Path(stamp["path"])
        if not path.is_file():
            raise DatasetLoadError(
                "SOURCE_NOT_FOUND", f"Referenced file missing: {path}"
            )
        if _fingerprint(path) != stamp:
            raise DatasetLoadError(
                "INVALID_DATA",
                "Metadata source or standardized table changed; load again.",
            )
    return descriptor


@register_tool("PerturbationMetadataTool")
class PerturbationMetadataTool(BaseTool):
    def run(self, arguments):
        directory = None
        try:
            preview_rows = arguments.get("preview_rows", 5)
            if type(preview_rows) is not int or not 0 <= preview_rows <= 20:
                raise DatasetLoadError(
                    "INVALID_ARGUMENT", "preview_rows must be an integer from 0 to 20."
                )
            for key, default, allowed in [
                ("missing_optional", "preserve", ("preserve", "error")),
                ("missing_perturbation", "error", ("error", "drop_rows")),
                ("duplicate_keys", "error", ("error", "collapse_identical")),
            ]:
                if arguments.get(key, default) not in allowed:
                    raise DatasetLoadError("INVALID_ARGUMENT", f"Invalid {key} policy.")
            identifier = arguments.get("metadata_id", uuid.uuid4().hex)
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,127}", identifier):
                raise DatasetLoadError("INVALID_ARGUMENT", "Invalid metadata_id.")
            source = Path(arguments["source"]).expanduser().resolve()
            if not source.is_file():
                raise DatasetLoadError(
                    "SOURCE_NOT_FOUND", f"Metadata source not found: {source}"
                )
            stamp = _fingerprint(source)
            original, fmt = read_metadata_table(
                source, arguments.get("format", "auto"), arguments.get("sheet_name")
            )
            table, warnings = _standardize(original, arguments)
            if _fingerprint(source) != stamp:
                raise DatasetLoadError(
                    "INVALID_DATA", "Metadata source changed while loading."
                )
            summary = {
                "experiment_level": arguments["experiment_level"],
                "n_source_rows": len(original),
                "n_rows": len(table),
                "fields": list(table.columns),
                "association_keys": arguments["association_keys"],
                "n_perturbations": int(table["perturbation"].nunique()),
                "n_controls": int(table["is_control"].sum())
                if "is_control" in table
                else None,
            }
            destination = (
                Path(arguments["output_dir"]).expanduser().resolve() / identifier
            )
            destination.parent.mkdir(parents=True, exist_ok=True)
            try:
                destination.mkdir()
            except FileExistsError:
                raise DatasetLoadError(
                    "DATASET_EXISTS", f"Metadata already exists: {destination}"
                ) from None
            directory = destination
            normalized = directory / "metadata.jsonl"
            table.to_json(
                normalized,
                orient="records",
                lines=True,
                force_ascii=False,
                double_precision=15,
            )
            options = {
                **arguments,
                "source": str(source),
                "format": fmt,
                "output_dir": str(destination.parent),
            }
            descriptor = {
                "schema_version": 1,
                "modality": "perturbation_metadata",
                "metadata_id": identifier,
                "source": stamp,
                "table": _fingerprint(normalized),
                "format": fmt,
                "read_options": options,
                "summary": summary,
                "warnings": warnings,
            }
            preview = json.loads(
                table.head(arguments.get("preview_rows", 5)).to_json(
                    orient="records", force_ascii=False, double_precision=15
                )
            )
            temporary = directory / "manifest.json.part"
            temporary.write_text(
                json.dumps(descriptor, indent=2, ensure_ascii=False) + "\n"
            )
            manifest = directory / "manifest.json"
            temporary.replace(manifest)
            return {
                "status": "success",
                "data": {
                    "metadata_ref": str(manifest),
                    "table_path": str(normalized),
                    "metadata_id": identifier,
                    **summary,
                    "preview": preview,
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
                for name in ["metadata.jsonl", "manifest.json.part"]:
                    (directory / name).unlink(missing_ok=True)
                directory.rmdir()
