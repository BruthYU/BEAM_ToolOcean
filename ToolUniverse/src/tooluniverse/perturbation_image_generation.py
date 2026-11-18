"""Shared process boundary for separately registered image generation backends."""

import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import jsonschema

from .base_tool import BaseTool
from .image_evaluation_io import file_path, image_rows, sha256
from .perturbation_runtime import finite_arguments, generation_result, wait_for_worker


class BaselineImageGenerationTool(BaseTool):
    """Validate explicit local assets and publish complete runs only."""

    backend = None

    def run(self, arguments):
        directory = None
        try:
            finite_arguments(arguments)
            jsonschema.validate(arguments, self.tool_config["parameter"])
            args = dict(arguments)
            sources = []

            def track(key):
                path = file_path(args[key])
                sources.append({"name": key, "path": str(path), "sha256": sha256(path)})
                args[key] = str(path)
                return path

            manifest = track("control_manifest")
            rows = image_rows(manifest)
            if any(row["condition"] != "DMSO" or not row.get("batch") for row in rows):
                raise ValueError(
                    "Controls must all be DMSO with explicit nonempty batch identifiers."
                )
            controls = [
                row
                for row in rows
                if not args.get("batch") or row["batch"] == args["batch"]
            ]
            if not controls:
                raise ValueError("No DMSO controls for the selected batch.")
            total = (
                len({row["batch"] for row in controls})
                * len(args["compounds"])
                * args.get("samples_per_condition", 8)
            )
            if total > 4096:
                raise ValueError("At most 4096 predictions are allowed per call.")
            repository = Path(args["model_dir"]).expanduser().resolve()
            args["model_dir"] = str(repository)
            if self.backend == "impa":
                relative_sources = ["IMPA/model.py"]
                for key in (
                    "checkpoint_path",
                    "embeddings_path",
                    "embedding_checkpoint",
                ):
                    track(key)
            else:
                raise ValueError("Unsupported image generation backend.")
            for name in relative_sources:
                path = file_path(repository / name)
                sources.append(
                    {"name": name, "path": str(path), "sha256": sha256(path)}
                )
            if any(sha256(row["path"]) != row["sha256"] for row in sources):
                raise ValueError("An input changed during preparation.")
            executable = shutil.which(args.get("python_executable", sys.executable))
            if not executable:
                raise ValueError("Model Python executable was not found.")
            run_id = self.backend + "_" + uuid.uuid4().hex
            directory = Path(args["output_dir"]).expanduser().resolve() / run_id
            directory.mkdir(parents=True)
            job = {
                **args,
                "backend": self.backend,
                "controls": controls,
                "source_files": sources,
                "run_id": run_id,
                "run_dir": str(directory),
            }
            request = directory / "request.json"
            request.write_text(json.dumps(job, indent=2))
            log = directory / "worker.log"
            with log.open("w") as stream:
                process = subprocess.Popen(
                    [
                        executable,
                        str(
                            Path(__file__).with_name("perturbation_baseline_worker.py")
                        ),
                        str(request),
                    ],
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
                        "error_code": "TIMEOUT",
                        "error": "Generation timed out; no complete run is published.",
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
                result = generation_result(directory, total)
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
            if directory:
                result["run_dir"] = str(directory)
            return result
