"""
perturbation_load_dataset

Load local expression data from H5AD, modern 10x H5, 10x MTX directories, or numeric CSV/TSV (inc...
"""

from collections.abc import Callable
from typing import Any, Optional

from ._shared_client import get_shared_client


def perturbation_load_dataset(
    source: str,
    output_dir: str,
    dataset_id: str | None = None,
    format: str | None = "auto",
    matrix_key: str | None = None,
    gene_id_column: str | None = None,
    obs_column_map: dict[str, Any] | None = None,
    control_values: list[str] | None = None,
    expression_scale: str | None = "unknown",
    backed: bool | None = None,
    metadata_path: str | None = None,
    metadata_join_keys: dict[str, Any] | None = None,
    metadata_match: str | None = "strict",
    matrix_orientation: str | None = None,
    index_column: str | None = None,
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> Any:
    """
    Load local expression data from H5AD, modern 10x H5, 10x MTX directories, or numeric CSV/TSV (inc...

    Parameters
    ----------
    source : str
        Local expression file or 10x MTX directory. MTX directory must directly conta...
    output_dir : str
        Directory for persistent dataset references.
    dataset_id : str
        Optional unique dataset name; generated when omitted.
    format : str
        auto detects H5AD, modern 10x H5, 10x MTX directory or CSV/TSV (optionally gz...
    matrix_key : str
        Expression matrix: X (default), raw/X or layers/<name>.
    gene_id_column : str
        Variable annotation column containing unique gene IDs; default variable index.
    obs_column_map : dict[str, Any]
        Standard observation fields mapped to existing original columns.
    control_values : list[str]
        Control labels; requires obs_column_map.perturbation.
    expression_scale : str
        Caller-declared expression scale; no normalization is performed.
    backed : bool
        Read H5AD using backed mode. Layers can still be materialized by AnnData.
    metadata_path : str
        Optional CSV, TSV or Parquet metadata file.
    metadata_join_keys : dict[str, Any]
        Equal-length join key lists; data key __obs_index__ refers to observation IDs.
    metadata_match : str
        Whether missing metadata matches fail or produce a counted warning.
    matrix_orientation : str
        Required for CSV/TSV. Orientation of the input table; output always uses obse...
    index_column : str
        Required for CSV/TSV. Exact header name of the column holding row IDs. Remain...
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
            "dataset_id": dataset_id,
            "format": format,
            "matrix_key": matrix_key,
            "gene_id_column": gene_id_column,
            "obs_column_map": obs_column_map,
            "control_values": control_values,
            "expression_scale": expression_scale,
            "backed": backed,
            "metadata_path": metadata_path,
            "metadata_join_keys": metadata_join_keys,
            "metadata_match": metadata_match,
            "matrix_orientation": matrix_orientation,
            "index_column": index_column,
        }.items()
        if v is not None
    }
    return get_shared_client().run_one_function(
        {
            "name": "perturbation.load_dataset",
            "arguments": _args,
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["perturbation_load_dataset"]
