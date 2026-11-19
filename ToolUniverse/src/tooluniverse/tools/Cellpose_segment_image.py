# ruff: noqa: N999
"""Local Cellpose single/batch segmentation and persistent complete results."""

from collections.abc import Callable
from typing import Any

from ._shared_client import get_shared_client


def Cellpose_segment_image(
    image_path: str | None = None,
    model_type: str | None = None,
    diameter: float | None = None,
    channels: list[Any] | None = None,
    save_mask: bool | None = None,
    mask_output_path: str | None = None,
    image_paths: list[Any] | None = None,
    image_ref: str | None = None,
    image_indices: list[Any] | None = None,
    output_dir: str | None = None,
    run_id: str | None = None,
    resume: bool | None = None,
    device: str | None = None,
    pretrained_model: str | None = None,
    batch_size: int | None = None,
    flow_threshold: float | None = None,
    cellprob_threshold: float | None = None,
    channel_indices: list[Any] | None = None,
    z_index: int | None = None,
    t_index: int | None = None,
    roi: dict[str, Any] | None = None,
    max_pixels: int | None = None,
    mask_format: str | None = None,
    save_input: bool | None = None,
    save_preview: bool | None = None,
    save_flows: bool | None = None,
    save_styles: bool | None = None,
    preview_rows: int | None = None,
    on_error: str | None = None,
    channel_names: list[str] | None = None,
    max_run_seconds: float | None = None,
    cancel_file: str | None = None,
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> dict[str, Any]:
    """Provide one input mode; batch/reference or complete results require output_dir.

    Without output_dir, image_path preserves the legacy statistics/optional PNG
    response. Complete mode returns run_ref/results_path and bounded previews;
    masks, object tables, selected inputs, overlays and requested flows/styles
    live in files. See examples/cellpose for exact options and resume semantics.
    """
    arguments = {
        "max_run_seconds": max_run_seconds,
        "cancel_file": cancel_file,
        "channel_names": channel_names,
        "image_path": image_path,
        "model_type": model_type,
        "diameter": diameter,
        "channels": channels,
        "save_mask": save_mask,
        "mask_output_path": mask_output_path,
        "image_paths": image_paths,
        "image_ref": image_ref,
        "image_indices": image_indices,
        "output_dir": output_dir,
        "run_id": run_id,
        "resume": resume,
        "device": device,
        "pretrained_model": pretrained_model,
        "batch_size": batch_size,
        "flow_threshold": flow_threshold,
        "cellprob_threshold": cellprob_threshold,
        "channel_indices": channel_indices,
        "z_index": z_index,
        "t_index": t_index,
        "roi": roi,
        "max_pixels": max_pixels,
        "mask_format": mask_format,
        "save_input": save_input,
        "save_preview": save_preview,
        "save_flows": save_flows,
        "save_styles": save_styles,
        "preview_rows": preview_rows,
        "on_error": on_error,
    }
    return get_shared_client().run_one_function(
        {
            "name": "Cellpose_segment_image",
            "arguments": {k: v for k, v in arguments.items() if v is not None},
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["Cellpose_segment_image"]
