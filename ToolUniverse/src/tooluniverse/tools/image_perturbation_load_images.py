"""
image_perturbation_load_images

Index microscopy TIFF/OME-TIFF filesets, BigTIFF, CZI, ND2, LIF, PNG/JPEG/BMP and local or anonym...
"""

from collections.abc import Callable
from typing import Any, Optional

from ._shared_client import get_shared_client


def image_perturbation_load_images(
    source: str,
    source_type: str,
    output_dir: str,
    dataset_id: str | None = None,
    format: str | None = None,
    file_pattern: str | None = None,
    image_root: str | None = None,
    image_path_columns: dict[str, Any] | None = None,
    metadata_column_map: dict[str, Any] | None = None,
    field_keys: list[str] | None = None,
    missing_images: str | None = None,
    axes: str | None = None,
    channel_map: dict[str, Any] | None = None,
    lazy: bool | None = None,
    series_index: int | None = None,
    resolution_level: int | None = None,
    group_path: str | None = None,
    position_index: int | None = None,
    remote_timeout: float | None = None,
    remote_endpoint: str | None = None,
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> Any:
    """
    Index microscopy TIFF/OME-TIFF filesets, BigTIFF, CZI, ND2, LIF, PNG/JPEG/BMP and local or anonym...

    Parameters
    ----------
    source : str
        Local image, vendor file, OME fileset entry, local store, manifest, or anonym...
    source_type : str
        file, directory or manifest for local sources; remote for one anonymous HTTP(...
    output_dir : str
        Directory for image index and dataset reference.
    dataset_id : str
        Unique dataset ID; generated if omitted.
    format : str
        Default auto; applies to image files, not manifest table format.
    file_pattern : str
        Directory glob, default top-level supported raster files and OME-Zarr stores;...
    image_root : str
        Root for manifest relative image paths; default manifest directory.
    image_path_columns : dict[str, Any]
        Manifest wide channel name to path-column mapping. Each row represents one fi...
    metadata_column_map : dict[str, Any]
        Standard metadata field to manifest column name.
    field_keys : list[str]
        Mapped metadata fields that uniquely identify each manifest row; otherwise ro...
    missing_images : str
        Default error; explicit skip_field skips entire incomplete fields and reports...
    axes : str
        Explicit axis order from TCZYX; required if TIFF axes are ambiguous.
    channel_map : dict[str, Any]
        Internal zero-based channel index to name. Must cover every channel; incompat...
    lazy : bool
        Default true indexes metadata without pixel decoding. False decodes one selec...
    series_index : int
        Zero-based TIFF series or OME-Zarr multiscales index. Required when more than...
    resolution_level : int
        Zero-based pyramid resolution, default zero (full resolution).
    group_path : str
        OME-Zarr image group inside source store (e.g. A/1/0); default root. No path ...
    position_index : int
        Zero-based mosaic position for LIF; required when the selected image has mult...
    remote_timeout : float
        Remote request timeout in seconds, default 30; remote mode only.
    remote_endpoint : str
        Optional anonymous S3-compatible endpoint URL, e.g. https://uk1s3.embassy.ebi...
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
            "source_type": source_type,
            "output_dir": output_dir,
            "dataset_id": dataset_id,
            "format": format,
            "file_pattern": file_pattern,
            "image_root": image_root,
            "image_path_columns": image_path_columns,
            "metadata_column_map": metadata_column_map,
            "field_keys": field_keys,
            "missing_images": missing_images,
            "axes": axes,
            "channel_map": channel_map,
            "lazy": lazy,
            "series_index": series_index,
            "resolution_level": resolution_level,
            "group_path": group_path,
            "position_index": position_index,
            "remote_timeout": remote_timeout,
            "remote_endpoint": remote_endpoint,
        }.items()
        if v is not None
    }
    return get_shared_client().run_one_function(
        {
            "name": "image_perturbation_load_images",
            "arguments": _args,
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["image_perturbation_load_images"]
