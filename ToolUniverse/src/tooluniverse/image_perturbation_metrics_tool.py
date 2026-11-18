"""Independent, model-agnostic MoA and FID/KID evaluation tools."""

import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import ClassVar

import jsonschema

from .base_tool import BaseTool
from .image_evaluation_io import (
    align_labels,
    evaluation_cohorts,
    file_path,
    image_rows,
    label_names,
    match_references,
    read_table,
    sha256,
)
from .perturbation_runtime import evaluation_result, finite_arguments, wait_for_worker
from .tool_registry import register_tool


@register_tool("ImagePerturbationEvaluationTool")
class ImagePerturbationEvaluationTool(BaseTool):
    """Shared evaluation workflow and compatibility for existing tool configs."""

    operation: ClassVar[str | None] = None

    def run(self, arguments):
        directory = None
        try:
            finite_arguments(arguments)
            jsonschema.validate(arguments, self.tool_config["parameter"])
            args = dict(arguments)
            operation = (
                self.operation
                or {
                    "ImagePerturbation_evaluate_moa": "moa",
                    "ImagePerturbation_calculate_fid_kid": "fid_kid",
                    "image_perturbation_evaluate_moa": "moa",
                    "image_perturbation_evaluate_fid_kid": "fid_kid",
                }[self.tool_config["name"]]
            )
            split = args.get("split", "test")
            if args["evaluation_profile"] != "bbbc021_smoke" and split != "test":
                raise ValueError(
                    "Standard/held-out profiles require the explicit test split."
                )
            fingerprints = {}
            sources = []

            def track(key):
                path = file_path(args[key])
                fingerprint = sha256(path)
                fingerprints[key] = fingerprint
                sources.append({"path": str(path), "sha256": fingerprint})
                return path

            predictions_path = track("predictions")
            predictions = image_rows(predictions_path)
            if args["evaluation_profile"] != "bbbc021_smoke" and any(
                row.get("condition_mapping_verified") is False for row in predictions
            ):
                raise ValueError(
                    "Predictions with an unverified condition map require the smoke evaluation profile."
                )
            if any(row.get("split", split) != split for row in predictions):
                raise ValueError(
                    "Prediction split conflicts with the selected reference split."
                )
            job = {
                **args,
                "operation": operation,
                "predictions": predictions,
                "split": split,
            }
            if operation == "moa":
                reference = track("reference_metadata")
                label_map = track("label_map")
                checkpoint = track("evaluator_checkpoint")
                labels = label_names(label_map)
                metadata = read_table(reference, {"condition", "moa"}, split)
                job["predictions"] = align_labels(
                    predictions, metadata, labels, args.get("join_on", "condition")
                )
                job["labels"] = labels
                job["evaluator_checkpoint"] = str(checkpoint)
            else:
                reference = track("reference_images")
                references = image_rows(reference, split)
                keys = args.get("match_keys", ["condition", "batch"])
                matched = match_references(predictions, references, keys)
                predictions, selected, composition = evaluation_cohorts(
                    predictions,
                    matched,
                    keys,
                    args.get("cohort_balance", "as_provided"),
                    args.get("seed", 42),
                )
                job.update(
                    predictions=predictions,
                    reference_images=selected,
                    reference_rows_excluded=len(references) - len(matched),
                    **composition,
                )
                count = min(len(predictions), len(selected))
                if count < 2:
                    raise ValueError(
                        "FID/KID requires at least two images in each cohort."
                    )
                if args.get("kid_subset_size", 100) > count:
                    raise ValueError(
                        "kid_subset_size exceeds a cohort size; set it explicitly for a small smoke test."
                    )
                weights = track("feature_extractor_weights")
                job["feature_extractor_weights"] = str(weights)
            if any(sha256(source["path"]) != source["sha256"] for source in sources):
                raise ValueError("An input changed during evaluation preparation.")
            executable = shutil.which(args.get("python_executable", sys.executable))
            if not executable:
                raise ValueError("Evaluation Python executable was not found.")
            run_id = "image_eval_" + uuid.uuid4().hex
            directory = Path(args["output_dir"]).expanduser().resolve() / run_id
            directory.mkdir(parents=True)
            job.update(
                {
                    "run_id": run_id,
                    "run_dir": str(directory),
                    "input_fingerprints": fingerprints,
                    "source_files": sources,
                }
            )
            request = directory / "request.json"
            request.write_text(json.dumps(job, indent=2))
            log = directory / "worker.log"
            with log.open("w") as stream:
                process = subprocess.Popen(
                    [
                        executable,
                        str(Path(__file__).with_name("image_evaluation_worker.py")),
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
                        "error": "Evaluation timed out; no complete report was published.",
                        "log_path": str(log),
                    }
            if code:
                return {
                    "status": "error",
                    "error_code": "EVALUATION_FAILED",
                    "error": log.read_text()[-4000:],
                    "log_path": str(log),
                }
            try:
                result = evaluation_result(directory, operation, len(predictions))
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
            KeyError,
            TypeError,
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
