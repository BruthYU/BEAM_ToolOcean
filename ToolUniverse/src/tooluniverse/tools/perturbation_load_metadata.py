"""
perturbation_load_metadata

Load and standardize perturbation experiment annotations for expression or image experiments: per...
"""

from collections.abc import Callable
from typing import Any, Optional

from ._shared_client import get_shared_client


def perturbation_load_metadata(
    source: str,
    output_dir: str,
    experiment_level: str,
    field_map: dict[str, Any],
    association_keys: list[str],
    format: str | None = None,
    sheet_name: str | None = None,
    constants: dict[str, Any] | None = None,
    control_values: list[str] | None = None,
    duplicate_keys: str | None = None,
    missing_optional: str | None = None,
    preview_rows: int | None = None,
    metadata_id: str | None = None,
    missing_perturbation: str | None = None,
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> Any:
    """
    Load and standardize perturbation experiment annotations for expression or image experiments: per...

    Parameters
    ----------
    source : str
        Local perturbation experiment metadata table or H5AD whose obs is read withou...
    output_dir : str
        Directory for standardized JSONL and persistent metadata reference.
    experiment_level : str
        One row per cell, sample, well or image field; level IDs are required and mus...
    field_map : dict[str, Any]
        Standard field to original column. perturbation and level IDs must be mapped ...
    association_keys : list[str]
        Standardized columns uniquely identifying records. Include cell_id; sample_id...
    format : str
        Default auto; CSV/TSV also support gzip. H5AD reads obs only.
    sheet_name : str
        XLSX sheet name; default first sheet, invalid for other formats.
    constants : dict[str, Any]
        Explicit standard field constants, e.g. dose_unit=uM. Cannot overlap field_ma...
    control_values : list[str]
        Exact perturbation labels denoting controls. Must match at least one row; mus...
    duplicate_keys : str
        Default error. collapse_identical only collapses identical standardized rows;...
    missing_optional : str
        Default preserve: optional blanks become null with counts. Keys, perturbation...
    preview_rows : int
        Number of actual standardized records returned for preview; default 5.
    metadata_id : str
        Unique metadata artifact ID; auto-generated when omitted.
    missing_perturbation : str
        Default error. Explicit drop_rows removes records whose perturbation label is...
    stream_callback : Callable, optional
        Callback for streaming output
    use_cache : bool, default False
        Enable caching
    validate : bool, default True
        Validate parameters

    Returns
    -------
    Any
    """
    # Handle mutable defaults to avoid B006 linting error

    # Strip None values so optional parameters don't trigger schema validation errors
    _args = {
        k: v
        for k, v in {
            "source": source,
            "output_dir": output_dir,
            "experiment_level": experiment_level,
            "field_map": field_map,
            "association_keys": association_keys,
            "format": format,
            "sheet_name": sheet_name,
            "constants": constants,
            "control_values": control_values,
            "duplicate_keys": duplicate_keys,
            "missing_optional": missing_optional,
            "preview_rows": preview_rows,
            "metadata_id": metadata_id,
            "missing_perturbation": missing_perturbation,
        }.items()
        if v is not None
    }
    return get_shared_client().run_one_function(
        {
            "name": "perturbation_load_metadata",
            "arguments": _args,
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["perturbation_load_metadata"]
