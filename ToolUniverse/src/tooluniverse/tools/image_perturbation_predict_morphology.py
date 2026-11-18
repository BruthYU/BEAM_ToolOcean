"""
Compatibility SDK wrapper for CellFlux_predict_morphology.

Prefer CellFlux_predict_morphology for new code; the legacy wire name remains an alias.
image_perturbation_predict_morphology

Generate unpaired treated-cell population images from preprocessed BBBC021 DMSO controls and offi...
"""

from collections.abc import Callable
from typing import Any

from ._shared_client import get_shared_client


def image_perturbation_predict_morphology(
    control_manifest: str,
    compounds: list[str],
    model_dir: str,
    checkpoint_path: str,
    embeddings_path: str,
    output_dir: str,
    batch: str | None = None,
    python_executable: str | None = None,
    device: str | None = "auto",
    samples_per_condition: int | None = 8,
    batch_size: int | None = 4,
    sampling_steps: int | None = 50,
    cfg_scale: float | None = 0.2,
    noise_level: float | None = 1.0,
    seed: int | None = 42,
    cpu_threads: int | None = 4,
    timeout_seconds: int | None = 1800,
    *,
    stream_callback: Callable[[str], None] | None = None,
    use_cache: bool = False,
    validate: bool = True,
) -> dict[str, Any]:
    """
    Generate unpaired treated-cell population images from preprocessed BBBC021 DMSO controls and offi...

    Parameters
    ----------
    control_manifest : str
        Control-only CSV with sample_id,image_path,batch,condition. condition must be...
    compounds : list[str]
        Exact compound identifiers in the supplied official 1024-feature embedding CS...
    model_dir : str
        Local trusted official CellFlux checkout (models/ and training/).
    checkpoint_path : str
        Official BBBC021 EMA checkpoint.pth, not an LFS pointer.
    embeddings_path : str
        Official emb_fp.csv, compound identifiers as first/index column and 1024 orde...
    output_dir : str
        Parent directory for a unique run containing PNG predictions, JSONL index, ma...
    batch : str
        Optional exact control batch filter. Otherwise generate separately for every ...
    python_executable : str
        Python executable in the model environment; defaults to the current interpreter.
    device : str
        auto selects CUDA when available, otherwise CPU.
    samples_per_condition : int
        Images per compound per batch. At most 4096 images per call. Controls resampl...
    batch_size : int
        Inference minibatch size; reduce for limited memory.
    sampling_steps : int
        Official EDM schedule parameter; heun2 solver. Low values are for smoke tests...
    cfg_scale : float
        Official classifier-free guidance scale.
    noise_level : float
        Gaussian noise added to normalized controls (official use_initial=2).
    seed : int
        Seed for control selection, dequantization and sampling noise.
    cpu_threads : int
        PyTorch intra-op CPU threads in the isolated worker.
    timeout_seconds : int
        Worker timeout; a timed-out or failed run is never returned as successful.
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
            "control_manifest": control_manifest,
            "compounds": compounds,
            "model_dir": model_dir,
            "checkpoint_path": checkpoint_path,
            "embeddings_path": embeddings_path,
            "output_dir": output_dir,
            "batch": batch,
            "python_executable": python_executable,
            "device": device,
            "samples_per_condition": samples_per_condition,
            "batch_size": batch_size,
            "sampling_steps": sampling_steps,
            "cfg_scale": cfg_scale,
            "noise_level": noise_level,
            "seed": seed,
            "cpu_threads": cpu_threads,
            "timeout_seconds": timeout_seconds,
        }.items()
        if v is not None
    }
    return get_shared_client().run_one_function(
        {
            "name": "image_perturbation_predict_morphology",
            "arguments": _args,
        },
        stream_callback=stream_callback,
        use_cache=use_cache,
        validate=validate,
    )


__all__ = ["image_perturbation_predict_morphology"]
