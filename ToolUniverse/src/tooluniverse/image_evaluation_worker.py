"""Standalone image evaluator. No generator model is imported or selected here."""

import csv
import hashlib
import importlib.metadata
import io
import json
import math
import sys
import time
from pathlib import Path

if __package__:
    from .image_evaluation_io import sha256
else:
    from image_evaluation_io import sha256

    # Executing a file inside tooluniverse otherwise places its package directory
    # first on sys.path. That would shadow stdlib profile during torchvision imports.
    _package_directory = Path(__file__).resolve().parent
    sys.path = [
        entry for entry in sys.path if Path(entry).resolve() != _package_directory
    ]


def read_cell(row, allowed_sizes=(96,)):
    import numpy as np
    from PIL import Image

    path = Path(row["image_path"])
    with path.open("rb") as stream:
        data = stream.read(1024 * 1024 + 1)
    if len(data) > 1024 * 1024:
        raise ValueError("Image exceeds the preprocessed single-cell size limit.")
    if hashlib.sha256(data).hexdigest() != row["sha256"]:
        raise ValueError("Image changed after manifest validation.")
    if path.suffix.lower() == ".npy":
        image = np.load(io.BytesIO(data), allow_pickle=False)
    else:
        with Image.open(io.BytesIO(data)) as source:
            if (
                source.format != "PNG"
                or source.mode != "RGB"
                or source.size not in {(size, size) for size in allowed_sizes}
            ):
                raise ValueError(
                    f"Expected a native RGB PNG with size in {allowed_sizes}; no implicit conversion is applied."
                )
            image = np.array(source)
    if (
        image.shape not in {(size, size, 3) for size in allowed_sizes}
        or image.dtype != np.uint8
    ):
        raise ValueError(
            f"Expected an HWC uint8 preprocessed cell with size in {allowed_sizes}."
        )
    return image


def extract_features(rows, extractor, device, batch_size, allowed_sizes=(96,)):
    import numpy as np
    import torch

    blocks = []
    cohort_shape = None
    with torch.inference_mode():
        for start in range(0, len(rows), batch_size):
            cells = [
                read_cell(row, allowed_sizes)
                for row in rows[start : start + batch_size]
            ]
            if cohort_shape is None:
                cohort_shape = cells[0].shape
            if any(cell.shape != cohort_shape for cell in cells):
                raise ValueError(
                    "Each image cohort must have one uniform native shape."
                )
            images = np.stack(cells)
            pixels = torch.from_numpy(images).permute(0, 3, 1, 2).to(device)
            # The official MoA evaluator converts uint8 -> float /255 -> uint8.
            # Direct uint8 preserves these stored pixels without extra dithering.
            features = extractor(pixels)[0].float().cpu()
            if (
                features.ndim != 2
                or features.shape[1] != 2048
                or not torch.isfinite(features).all()
            ):
                raise ValueError("Unexpected/non-finite Inception features.")
            blocks.append(features)
    return torch.cat(blocks)


def classification_metrics(true, predicted, labels):
    from sklearn.metrics import (
        accuracy_score,
        classification_report,
        confusion_matrix,
        f1_score,
    )

    per_class = classification_report(
        true,
        predicted,
        labels=list(range(len(labels))),
        target_names=labels,
        zero_division=0,
        output_dict=True,
    )
    return {
        "accuracy": float(accuracy_score(true, predicted)),
        "macro_f1": float(f1_score(true, predicted, average="macro", zero_division=0)),
        "weighted_f1": float(
            f1_score(true, predicted, average="weighted", zero_division=0)
        ),
        "per_class": {label: per_class[label] for label in labels},
        "confusion_matrix": confusion_matrix(
            true, predicted, labels=list(range(len(labels)))
        ).tolist(),
        "confusion_matrix_labels": labels,
        "macro_f1_label_policy": "union_of_true_and_predicted_labels; zero_division=0",
    }


def distribution_metrics(real, predicted, subset_size, subsets, seed):
    import torch
    from torch import nn
    from torchmetrics.image.fid import FrechetInceptionDistance
    from torchmetrics.image.kid import KernelInceptionDistance

    if min(len(real), len(predicted)) < max(2, subset_size):
        raise ValueError(
            "Insufficient samples for the requested FID/KID configuration."
        )

    class PrecomputedFeatures(nn.Module):
        def __init__(self, dimensions):
            super().__init__()
            self.num_features = dimensions

        def forward(self, features):
            return features

    # Public custom-extractor interfaces let both metrics consume the same frozen
    # features. FID covariance arithmetic is float64; KID follows float32 features.
    fid = FrechetInceptionDistance(
        feature=PrecomputedFeatures(real.shape[1])
    ).set_dtype(torch.float64)
    fid.update(real.double(), real=True)
    fid.update(predicted.double(), real=False)
    kid = KernelInceptionDistance(
        feature=PrecomputedFeatures(real.shape[1]),
        subset_size=subset_size,
        subsets=subsets,
    )
    kid.update(real, real=True)
    kid.update(predicted, real=False)
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(seed)
        fid_value = float(fid.compute())
        kid_mean, kid_std = kid.compute()
    result = {"fid": fid_value, "kid_mean": float(kid_mean), "kid_std": float(kid_std)}
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError("Non-finite FID/KID; no numeric score is published.")
    return result


def main(request_path):
    import torch
    from torch import nn
    from torch_fidelity.feature_extractor_inceptionv3 import FeatureExtractorInceptionV3

    job = json.loads(Path(request_path).read_text())
    for source in job["source_files"]:
        if sha256(source["path"]) != source["sha256"]:
            raise ValueError("An evaluation input changed after validation.")
    torch.set_num_threads(job.get("cpu_threads", 4))
    seed = job.get("seed", 42)
    torch.manual_seed(seed)
    device = job.get("device", "auto")
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA is unavailable on this host.")
    directory = Path(job["run_dir"])
    temporary_weights = None
    if job["operation"] == "moa":
        checkpoint = torch.load(
            job["evaluator_checkpoint"], map_location="cpu", weights_only=True
        )
        state = checkpoint["model_state"]
        prefix = "feature_extractor.inception."
        inception = {
            key[len(prefix) :]: value
            for key, value in state.items()
            if key.startswith(prefix)
        }
        if not inception:
            raise ValueError("Evaluator checkpoint has no frozen Inception weights.")
        temporary_weights = directory / "temporary_inception.pth"
        torch.save(inception, temporary_weights)
        weights_path = str(temporary_weights)
        head = nn.Sequential(
            nn.Linear(2048, 512),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(512, len(job["labels"])),
        )
        head.load_state_dict(
            {
                key[len("classifier.") :]: value
                for key, value in state.items()
                if key.startswith("classifier.")
            },
            strict=True,
        )
        head.eval()
        del checkpoint, state, inception
    else:
        weights_path = job["feature_extractor_weights"]
    # Passing a local file prevents implicit weight downloads at evaluation time.
    extractor = FeatureExtractorInceptionV3(
        "inception-v3-compat", ["2048"], feature_extractor_weights_path=weights_path
    )
    extractor.eval().requires_grad_(False).to(device)
    if temporary_weights:
        temporary_weights.unlink()
    start = time.monotonic()
    print(f"Evaluating {len(job['predictions'])} predictions on {device}", flush=True)
    allowed_sizes = (
        (96, 128)
        if job["evaluation_profile"] in ("bbbc021_native_inception", "bbbc021_smoke")
        else (96,)
    )
    predicted_features = extract_features(
        job["predictions"], extractor, device, job.get("batch_size", 16), allowed_sizes
    )
    result = {
        "schema_version": "1.0",
        "run_id": job["run_id"],
        "operation": job["operation"],
        "model_name": job["model_name"],
        "evaluation_profile": job["evaluation_profile"],
        "split": job["split"],
        "n_predictions": len(job["predictions"]),
        "device": device,
        "seed": seed,
        "input_fingerprints": job["input_fingerprints"],
        "preprocessing": "stored uint8 HWC; original channel order; no external resizing, augmentation or contrast normalization; frozen Inception performs its own internal resize",
        "prediction_shape": list(read_cell(job["predictions"][0], allowed_sizes).shape),
        "feature_extractor": "torch-fidelity inception-v3-compat 2048",
        "benchmark_verified": False,
        "prediction_condition_mapping_verified": all(
            row.get("condition_mapping_verified") is True for row in job["predictions"]
        )
        if any("condition_mapping_verified" in row for row in job["predictions"])
        else None,
        "scope": "integration_smoke"
        if job["evaluation_profile"] == "bbbc021_smoke"
        else "user_supplied_evaluation_cohort",
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("torch", "torchmetrics", "torch-fidelity", "scikit-learn")
        },
    }
    if job["operation"] == "moa":
        with torch.inference_mode():
            scores = head(predicted_features).softmax(dim=1)
        if not torch.isfinite(scores).all():
            raise ValueError(
                "Non-finite classifier scores; no predictions are published."
            )
        predicted = scores.argmax(dim=1).tolist()
        true = [row["true_label_id"] for row in job["predictions"]]
        details = classification_metrics(true, predicted, job["labels"])
        result["metrics"] = {
            key: details[key] for key in ("accuracy", "macro_f1", "weighted_f1")
        }
        result["label_map"] = job["labels"]
        result["join_on"] = job.get("join_on", "condition")
        details_path = directory / "moa_metrics.json"
        details_path.write_text(json.dumps(details, indent=2, allow_nan=False))
        result["details_path"] = str(details_path)
        predictions_path = directory / "moa_predictions.csv"
        with predictions_path.open("w", newline="") as stream:
            fields = [
                "sample_id",
                "image_path",
                "condition",
                "true_moa",
                "predicted_moa",
                "confidence",
                "correct",
            ]
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for i, row in enumerate(job["predictions"]):
                writer.writerow(
                    {
                        **{key: row[key] for key in fields[:4]},
                        "predicted_moa": job["labels"][predicted[i]],
                        "confidence": float(scores[i].max()),
                        "correct": true[i] == predicted[i],
                    }
                )
        result["predictions_path"] = str(predictions_path)
        per_condition = {}
        for condition in sorted({row["condition"] for row in job["predictions"]}):
            positions = [
                i
                for i, row in enumerate(job["predictions"])
                if row["condition"] == condition
            ]
            metrics = classification_metrics(
                [true[i] for i in positions],
                [predicted[i] for i in positions],
                job["labels"],
            )
            per_condition[condition] = {
                "n": len(positions),
                **{
                    key: metrics[key] for key in ("accuracy", "macro_f1", "weighted_f1")
                },
            }
    else:
        real_features = extract_features(
            job["reference_images"],
            extractor,
            device,
            job.get("batch_size", 16),
            allowed_sizes,
        )
        subset = job.get("kid_subset_size", 100)
        subsets = job.get("kid_subsets", 100)
        result["metrics"] = distribution_metrics(
            real_features, predicted_features, subset, subsets, seed
        )
        result.update(
            {
                "n_reference_images": len(job["reference_images"]),
                "reference_shape": list(
                    read_cell(job["reference_images"][0], allowed_sizes).shape
                ),
                "reference_rows_excluded": job["reference_rows_excluded"],
                "match_keys": job.get("match_keys", ["condition", "batch"]),
                "kid_subset_size": subset,
                "kid_subsets": subsets,
            }
        )
        per_condition = {}
        if job.get("per_condition", True):
            for condition in sorted({row["condition"] for row in job["predictions"]}):
                pred_idx = [
                    i
                    for i, row in enumerate(job["predictions"])
                    if row["condition"] == condition
                ]
                real_idx = [
                    i
                    for i, row in enumerate(job["reference_images"])
                    if row["condition"] == condition
                ]
                item = {
                    "n_predictions": len(pred_idx),
                    "n_reference_images": len(real_idx),
                }
                if min(len(pred_idx), len(real_idx)) < max(2, subset):
                    item.update({"status": "insufficient_samples", "metrics": None})
                else:
                    item.update(
                        {
                            "status": "success",
                            "metrics": distribution_metrics(
                                real_features[real_idx],
                                predicted_features[pred_idx],
                                subset,
                                subsets,
                                seed,
                            ),
                        }
                    )
                per_condition[condition] = item
    condition_path = directory / "per_condition.json"
    condition_path.write_text(json.dumps(per_condition, indent=2, allow_nan=False))
    result["per_condition_path"] = str(condition_path)
    result["evaluation_seconds"] = time.monotonic() - start
    for source in job["source_files"]:
        if sha256(source["path"]) != source["sha256"]:
            raise ValueError("An evaluation input changed during execution.")
    result["report_path"] = str(directory / "report.json")
    Path(result["report_path"]).write_text(
        json.dumps(result, indent=2, allow_nan=False)
    )
    (directory / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False)
    )


if __name__ == "__main__":
    main(sys.argv[1])
