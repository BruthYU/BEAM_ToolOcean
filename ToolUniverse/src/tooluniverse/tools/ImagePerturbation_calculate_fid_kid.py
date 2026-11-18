# ruff: noqa: N999
"""
ImagePerturbation_calculate_fid_kid

Compute model-independent FID and KID for generated versus real treated BBBC021 images using one ...
"""

from collections.abc import Callable
from typing import Any

from ._shared_client import get_shared_client


def ImagePerturbation_calculate_fid_kid(
    predictions: str,
    model_name: str,
    evaluation_profile: str,
    output_dir: str,
    reference_images: str,
    feature_extractor_weights: str,
    split: str | None = "test",
    python_executable: str | None = None,
    device: str | None = "auto",
    batch_size: int | None = 16,
    seed: int | None = 42,
    cpu_threads: int | None = 4,
    timeout_seconds: int | None = 1800,
    match_keys: list[str] | None = None,
    kid_subset_size: int | None = 100,
    kid_subsets: int | None = 100,
    per_condition: bool | None = True,
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> dict[str, Any]:
    """
    Compute model-independent FID and KID for generated versus real treated BBBC021 images using one ...

    Parameters
    ----------
    predictions : str
        CSV/JSONL image index: unique sample_id,image_path,condition; batch when matc...
    model_name : str
        Report label such as CellFlux, IMPA, or a custom model. The label does not select evaluator logic or weights.
    evaluation_profile : str
        BBBC021 protocol intent. Standard/unseen profiles require native 96x96 images...
    output_dir : str
        Parent directory for a unique evaluation run, metrics reports, input fingerpr...
    split : str
        Exact reference split to select. Reference tables must contain split. No impl...
    python_executable : str
        Evaluation environment Python; defaults to the current interpreter. Requires ...
    device : str
        Device for frozen Inception features. Cached metrics run on CPU; auto uses CU...
    batch_size : int
        Feature extraction minibatch size.
    seed : int
        Random seed, including reproducible KID subsets.
    cpu_threads : int
        PyTorch CPU thread count in the worker.
    timeout_seconds : int
        Worker timeout. Failed or incomplete reports are returned as errors.
    reference_images : str
        Real treated-image CSV/JSONL with sample_id,image_path,condition,split and re...
    feature_extractor_weights : str
        Local torch-fidelity inception-v3-compat state_dict. Must use the same frozen...
    match_keys : list[str]
        Reference matching fields; must include condition. Default also enforces batc...
    kid_subset_size : int
        KID subset size. Must not exceed either cohort; small tests must explicitly l...
    kid_subsets : int
        Number of reproducible KID subset draws. Unbiased KID may be negative.
    per_condition : bool
        Also compute metrics per condition where sample counts meet the configured mi...
    stream_callback : Callable, optional
        Callback for streaming output
    use_cache : bool, default False
        Enable caching
    validate : bool, default True
        Validate parameters

    Returns
    -------
    dict[str, Any]
    """
    # Handle mutable defaults to avoid B006 linting error
    if match_keys is None:
        match_keys = ["condition", "batch"]
    # Strip None values so optional parameters don't trigger schema validation errors
    _args = {
        k: v
        for k, v in {
            "predictions": predictions,
            "model_name": model_name,
            "evaluation_profile": evaluation_profile,
            "output_dir": output_dir,
            "split": split,
            "python_executable": python_executable,
            "device": device,
            "batch_size": batch_size,
            "seed": seed,
            "cpu_threads": cpu_threads,
            "timeout_seconds": timeout_seconds,
            "reference_images": reference_images,
            "feature_extractor_weights": feature_extractor_weights,
            "match_keys": match_keys,
            "kid_subset_size": kid_subset_size,
            "kid_subsets": kid_subsets,
            "per_condition": per_condition,
        }.items()
        if v is not None
    }
    return get_shared_client().run_one_function(
        {
            "name": "ImagePerturbation_calculate_fid_kid",
            "arguments": _args,
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["ImagePerturbation_calculate_fid_kid"]
