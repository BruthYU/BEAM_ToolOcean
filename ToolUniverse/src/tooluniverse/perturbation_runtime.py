"""Small shared guards for local perturbation workers and result publication."""

import json
import math
import os
import signal
import subprocess
from pathlib import Path


def finite_arguments(value, path="arguments"):
    """JSON Schema permits non-finite Python floats; reject them before any I/O."""
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{path} must be finite.")
    if isinstance(value, dict):
        for key, item in value.items():
            finite_arguments(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            finite_arguments(item, f"{path}[{index}]")


def wait_for_worker(process, timeout):
    """Reap the complete process group on timeout OR interruption of its caller.

    Workers must be launched with start_new_session=True. Preserve the original
    exception if the process already exited while cleanup was being requested.
    """
    try:
        return process.wait(timeout=timeout)
    except BaseException:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=10)
        except (ProcessLookupError, subprocess.TimeoutExpired):
            pass
        raise


def worker_result(directory):
    """Never publish a missing, malformed, foreign, or non-finite worker result."""
    directory = Path(directory).resolve()
    path = directory / "result.json"
    if not path.is_file() or path.stat().st_size > 8 * 1024**2:
        raise ValueError("Missing or oversized worker result.json.")
    result = json.loads(path.read_text())
    if not isinstance(result, dict):
        raise TypeError("Worker result must be a JSON object.")
    finite_arguments(result, "result")
    if result.get("schema_version") != "1.0" or result.get("run_id") != directory.name:
        raise ValueError("Worker result has an unexpected schema_version or run_id.")
    if result.get("status", "success") != "success":
        raise ValueError("Worker did not publish a successful result.")
    if not isinstance(result.get("device"), str) or not result["device"]:
        raise ValueError("Worker result is missing its actual device.")
    return result


def published_file(directory, value):
    if not isinstance(value, str):
        raise TypeError("Worker result is missing an output path.")
    path = Path(value).resolve()
    if not path.is_relative_to(Path(directory).resolve()) or not path.is_file():
        raise ValueError("Worker output is missing or outside the run directory.")
    return path


def generation_result(directory, expected_count):
    from .image_evaluation_io import image_rows

    result = worker_result(directory)
    if type(result.get("n_images")) is not int or result["n_images"] != expected_count:
        raise ValueError("Worker returned an unexpected image count.")
    index = published_file(directory, result.get("image_index"))
    published_file(directory, result.get("run_manifest"))
    image_dir = result.get("image_dir")
    if (
        not isinstance(image_dir, str)
        or not Path(image_dir).resolve().is_relative_to(Path(directory).resolve())
        or not Path(image_dir).is_dir()
    ):
        raise ValueError("Worker result is missing its local image directory.")
    if not isinstance(result.get("model"), str) or not result["model"]:
        raise ValueError("Worker result is missing the model identity.")
    rows = image_rows(index)
    if len(rows) != expected_count:
        raise ValueError("Worker image index count differs from its result.")
    for row in rows:
        published_file(directory, row["image_path"])
    return result


def evaluation_result(directory, operation, expected_count=None):
    result = worker_result(directory)
    if result.get("operation") != operation:
        raise ValueError("Worker returned a different evaluation operation.")
    count = result.get("n_predictions")
    if (
        type(count) is not int
        or count < 1
        or (expected_count is not None and count != expected_count)
    ):
        raise ValueError("Worker prediction count differs from the requested cohort.")
    if operation == "fid_kid" and (
        type(result.get("n_reference_images")) is not int
        or result["n_reference_images"] < 2
    ):
        raise ValueError("Worker returned an invalid reference image count.")
    published_file(directory, result.get("report_path"))
    published_file(directory, result.get("per_condition_path"))
    metrics = result.get("metrics")
    required = (
        ("accuracy", "macro_f1", "weighted_f1")
        if operation == "moa"
        else ("fid", "kid_mean", "kid_std")
    )
    if not isinstance(metrics, dict) or any(
        type(metrics.get(key)) not in (int, float) for key in required
    ):
        raise ValueError("Worker returned incomplete/non-numeric evaluation metrics.")
    if operation == "moa":
        if any(not 0 <= metrics[key] <= 1 for key in required):
            raise ValueError("Classification metrics must be in [0,1].")
        published_file(directory, result.get("details_path"))
        published_file(directory, result.get("predictions_path"))
    return result
