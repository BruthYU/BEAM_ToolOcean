"""
Compatibility SDK wrapper for ImagePerturbation_evaluate_moa.

Prefer ImagePerturbation_evaluate_moa for new code; the legacy wire name remains an alias.
image_perturbation_evaluate_moa

Evaluate images from any BBBC021-compatible generator with an independently supplied frozen MoA c...
"""

from collections.abc import Callable
from typing import Any

from ._shared_client import get_shared_client


def image_perturbation_evaluate_moa(
    predictions: str,
    model_name: str,
    evaluation_profile: str,
    output_dir: str,
    reference_metadata: str,
    evaluator_checkpoint: str,
    label_map: str,
    split: str | None = "test",
    python_executable: str | None = None,
    device: str | None = "auto",
    batch_size: int | None = 16,
    seed: int | None = 42,
    cpu_threads: int | None = 4,
    timeout_seconds: int | None = 1800,
    join_on: str | None = "condition",
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> dict[str, Any]:
    """
    Evaluate images from any BBBC021-compatible generator with an independently supplied frozen MoA c...

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
    reference_metadata : str
        Ground-truth CSV/JSONL with condition,moa,split; sample_id additionally requi...
    evaluator_checkpoint : str
        Local frozen CellFlux-compatible MoA evaluator checkpoint containing model_st...
    label_map : str
        JSON ordered list of 13 unique MoA class names, or {labels:[...]}, in this ev...
    join_on : str
        Join labels by condition only when every condition has one unambiguous MoA. F...
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
            "reference_metadata": reference_metadata,
            "evaluator_checkpoint": evaluator_checkpoint,
            "label_map": label_map,
            "join_on": join_on,
        }.items()
        if v is not None
    }
    return get_shared_client().run_one_function(
        {
            "name": "image_perturbation_evaluate_moa",
            "arguments": _args,
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["image_perturbation_evaluate_moa"]
