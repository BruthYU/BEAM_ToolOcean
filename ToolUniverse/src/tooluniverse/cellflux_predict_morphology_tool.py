"""Run the official CellFlux BBBC021 model in a separate Python process."""

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import jsonschema

from .base_tool import BaseTool
from .perturbation_runtime import (
    configured_worker_arguments,
    generation_result,
    wait_for_worker,
)
from .tool_registry import register_tool


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"File does not exist: {path}")
    return path


def _controls(manifest, batch=None):
    """Resolve an explicit control-only manifest, without guessing labels."""
    if manifest.stat().st_size > 128 * 1024**2:
        raise ValueError("Control CSV exceeds the 128 MiB manifest limit.")
    with manifest.open(newline="") as stream:
        reader = csv.DictReader(stream)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or any(
            not key or not key.strip() for key in header
        ):
            raise ValueError("Control CSV headers must be unique and nonempty.")
        if not {"sample_id", "image_path", "batch", "condition"} <= set(
            reader.fieldnames or []
        ):
            raise ValueError(
                "Control CSV requires sample_id,image_path,batch,condition columns."
            )
        rows = []
        for row in reader:
            if len(rows) >= 100000:
                raise ValueError("Control CSV must contain 1 to 100000 rows.")
            rows.append(row)
    if not rows or len(rows) > 100000:
        raise ValueError("Control CSV must contain 1 to 100000 rows.")
    seen = set()
    selected = []
    for row in rows:
        if not all(
            isinstance(row.get(key), str) and row[key].strip()
            for key in ("sample_id", "image_path", "batch")
        ):
            raise ValueError("Missing control identity, image path or batch.")
        if row["condition"] != "DMSO":
            raise ValueError("Control CSV must contain DMSO cells only.")
        if row["sample_id"] in seen:
            raise ValueError(f"Duplicate sample_id: {row['sample_id']}")
        seen.add(row["sample_id"])
        if batch is not None and row["batch"] != batch:
            continue
        image = Path(row["image_path"]).expanduser()
        image = _file(image if image.is_absolute() else manifest.parent / image)
        if image.suffix.lower() not in (".npy", ".png"):
            raise ValueError("CellFlux v1 accepts preprocessed NPY or PNG cells only.")
        if image.stat().st_size > 1024 * 1024:
            raise ValueError("Control image exceeds the single-cell size limit.")
        selected.append({**row, "image_path": str(image), "sha256": _sha256(image)})
    if not selected:
        raise ValueError("No controls remain for the requested batch.")
    return selected


@register_tool("CellFluxTool")
class CellFluxTool(BaseTool):
    """Predict cell populations with fixed BBBC021 preprocessing and EMA weights."""

    def run(self, arguments):
        directory = None
        try:
            args = configured_worker_arguments(
                self.tool_config, arguments, ("model_dir", "checkpoint_path", "embeddings_path")
            )
            manifest = _file(args["control_manifest"])
            checkpoint = _file(args["checkpoint_path"])
            embeddings = _file(args["embeddings_path"])
            repository = Path(args["model_dir"]).expanduser().resolve()
            for relative in (
                "models/model_configs.py",
                "models/ema.py",
                "training/edm_time_discretization.py",
            ):
                _file(repository / relative)
            source_files = []
            for relative in (
                "models/model_configs.py",
                "models/ema.py",
                "models/unet.py",
                "models/nn.py",
                "training/edm_time_discretization.py",
            ):
                path = repository / relative
                if path.is_file():
                    source_files.append({"path": str(path), "sha256": _sha256(path)})
            executable = args.get("python_executable", sys.executable)
            executable = shutil.which(executable)
            if not executable:
                raise ValueError("python_executable was not found.")
            manifest_hash = _sha256(manifest)
            controls = _controls(manifest, args.get("batch"))
            if _sha256(manifest) != manifest_hash:
                raise ValueError("Control manifest changed while being read.")
            batches = sorted({row["batch"] for row in controls})
            count = (
                len(batches)
                * len(args["compounds"])
                * args.get("samples_per_condition", 8)
            )
            if count > 4096:
                raise ValueError("A call may generate at most 4096 images.")
            output = Path(args["output_dir"]).expanduser().resolve()
            run_id = "cellflux_" + uuid.uuid4().hex
            directory = output / run_id
            directory.mkdir(parents=True)
            job = {
                **args,
                "model_dir": str(repository),
                "checkpoint_path": str(checkpoint),
                "embeddings_path": str(embeddings),
                "controls": controls,
                "run_id": run_id,
                "run_dir": str(directory),
                "input_fingerprints": {
                    "control_manifest": manifest_hash,
                    "checkpoint": _sha256(checkpoint),
                    "embeddings": _sha256(embeddings),
                },
            }
            source_files.extend(
                {"path": str(path), "sha256": job["input_fingerprints"][name]}
                for name, path in (
                    ("control_manifest", manifest),
                    ("checkpoint", checkpoint),
                    ("embeddings", embeddings),
                )
            )
            job["source_files"] = source_files
            request = directory / "request.json"
            request.write_text(json.dumps(job, indent=2))
            worker = Path(__file__).with_name("cellflux_inference_worker.py")
            log = directory / "worker.log"
            with log.open("w") as stream:
                process = subprocess.Popen(
                    [executable, str(worker), str(request)],
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                    env={**os.environ, "PYTHONUNBUFFERED": "1"},
                )
                try:
                    code = wait_for_worker(process, args.get("timeout_seconds", 1800))
                except subprocess.TimeoutExpired:
                    return {
                        "status": "error",
                        "error": "CellFlux inference timed out; output is incomplete.",
                        "error_code": "TIMEOUT",
                        "log_path": str(log),
                    }
            if code:
                return {
                    "status": "error",
                    "error_code": "INFERENCE_FAILED",
                    "error": log.read_text()[-4000:],
                    "log_path": str(log),
                }
            try:
                result = generation_result(directory, count)
            except (ValueError, OSError, TypeError, KeyError) as exc:
                return {
                    "status": "error",
                    "error_code": "INVALID_OUTPUT",
                    "error": str(exc),
                    "run_dir": str(directory),
                    "log_path": str(log),
                }
            return {"status": "success", **result, "log_path": str(log)}
        except (
            ValueError,
            OSError,
            TypeError,
            KeyError,
            jsonschema.ValidationError,
        ) as exc:
            result = {
                "status": "error",
                "error_code": "INVALID_INPUT",
                "error": str(exc),
            }
            if directory is not None:
                result["run_dir"] = str(directory)
            return result
