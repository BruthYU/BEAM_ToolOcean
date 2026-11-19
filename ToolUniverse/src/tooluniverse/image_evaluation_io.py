"""Model-independent image evaluation manifests and label alignment."""

import csv
import hashlib
import json
import random
from collections import defaultdict
from pathlib import Path


def file_path(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"Missing local file: {path}")
    return path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path, required, split=None):
    """Read a bounded canonical CSV/JSONL table, filtering explicit splits only."""
    path = file_path(path)
    if path.stat().st_size > 128 * 1024 * 1024:
        raise ValueError(
            "Manifest exceeds 128 MiB; supply a smaller evaluation cohort."
        )
    required = set(required) | ({"split"} if split is not None else set())
    with path.open(newline="") as stream:
        if path.suffix.lower() == ".csv":
            reader = csv.DictReader(stream)
            header = reader.fieldnames or []
            if len(header) != len(set(header)) or any(
                not key or not key.strip() for key in header
            ):
                raise ValueError("Manifest headers must be unique and nonempty.")
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f"Missing columns in {path}: {sorted(required)}")
        elif path.suffix.lower() == ".jsonl":
            reader = (json.loads(line) for line in stream if line.strip())
        else:
            raise ValueError("Evaluation manifests must be CSV or JSONL.")
        rows = []
        for row in reader:
            if not isinstance(row, dict) or not required <= row.keys():
                raise ValueError(f"Record missing required fields: {sorted(required)}")
            if split is not None and row["split"] != split:
                continue
            for key in required:
                if not isinstance(row[key], str) or not row[key].strip():
                    raise ValueError(f"Field {key} must be a nonempty string.")
            rows.append(row)
            if len(rows) > 50000:
                raise ValueError("Evaluation cohort exceeds 50000 selected rows.")
    if not rows:
        raise ValueError("No rows remain in the requested evaluation split.")
    return rows


def image_rows(path, split=None):
    path = file_path(path)
    rows = read_table(path, {"sample_id", "image_path", "condition"}, split)
    ids, paths = set(), set()
    result = []
    for row in rows:
        row = dict(row)
        for key in ("condition_mapping_verified", "benchmark_verified"):
            if key in row:
                value = row[key]
                if isinstance(value, str) and value.strip().lower() in (
                    "true",
                    "false",
                ):
                    value = value.strip().lower() == "true"
                if not isinstance(value, bool):
                    raise ValueError(f"{key} must be a boolean (CSV: true/false).")
                row[key] = value
        relative = Path(row["image_path"]).expanduser()
        image = file_path(
            relative if relative.is_absolute() else path.parent / relative
        )
        if image.suffix.lower() not in (".png", ".npy"):
            raise ValueError("Images must be preprocessed uint8 PNG or NPY cells.")
        if image.stat().st_size > 1024 * 1024:
            raise ValueError("Image exceeds the preprocessed single-cell size limit.")
        if row["sample_id"] in ids or str(image) in paths:
            raise ValueError(
                "Duplicate sample_id or image path in the selected cohort."
            )
        ids.add(row["sample_id"])
        paths.add(str(image))
        digest = sha256(image)
        if row.get("sha256") and row["sha256"] != digest:
            raise ValueError(f"Image fingerprint mismatch: {image}")
        result.append({**row, "image_path": str(image), "sha256": digest})
    return result


def label_names(path):
    labels = json.loads(file_path(path).read_text())
    if isinstance(labels, dict):
        labels = labels.get("labels")
    if (
        not isinstance(labels, list)
        or len(labels) != 13
        or any(not isinstance(x, str) or not x.strip() for x in labels)
        or len(set(labels)) != len(labels)
    ):
        raise ValueError(
            "BBBC021 label_map must contain 13 unique names in classifier-output order."
        )
    return labels


def align_labels(predictions, metadata, labels, join_on):
    """Ground truth comes exclusively from references, never generated annotations."""
    label_to_id = {label: i for i, label in enumerate(labels)}
    lookup = {}
    for row in metadata:
        key = row.get(join_on)
        if not isinstance(key, str) or not key.strip():
            raise ValueError(
                f"Reference metadata requires {join_on} for label alignment."
            )
        if row["moa"] not in label_to_id:
            raise ValueError(f"Unknown reference MoA label: {row['moa']}")
        if key in lookup and (
            join_on == "sample_id" or lookup[key]["moa"] != row["moa"]
        ):
            raise ValueError(f"Ambiguous reference label for {key}.")
        lookup[key] = row
    aligned = []
    for row in predictions:
        key = (
            row.get("reference_sample_id", row["sample_id"])
            if join_on == "sample_id"
            else row["condition"]
        )
        if key not in lookup:
            raise ValueError(f"No reference label for {join_on}={key}.")
        reference = lookup[key]
        if row["condition"] != reference["condition"]:
            raise ValueError("Prediction and reference conditions disagree.")
        aligned.append(
            {
                **row,
                "true_moa": reference["moa"],
                "true_label_id": label_to_id[reference["moa"]],
            }
        )
    return aligned


def match_references(predictions, references, keys):
    if "condition" not in keys:
        raise ValueError("match_keys must include condition.")

    def group(row):
        if any(
            not isinstance(row.get(key), str) or not row[key].strip() for key in keys
        ):
            raise ValueError(f"Missing match fields: {keys}")
        return tuple(row[key] for key in keys)

    requested = {group(row) for row in predictions}
    available = {group(row) for row in references}
    if not requested <= available:
        raise ValueError(
            f"No real images for prediction groups: {sorted(requested - available)}"
        )
    selected = [row for row in references if group(row) in requested]
    if {row["image_path"] for row in predictions} & {
        row["image_path"] for row in selected
    }:
        raise ValueError("Predictions and reference cohorts share an image path.")
    if {row.get("sha256") for row in predictions if row.get("sha256")} & {
        row.get("sha256") for row in selected if row.get("sha256")
    }:
        raise ValueError(
            "Predictions and reference cohorts share identical image content."
        )
    return selected


def evaluation_cohorts(predictions, references, keys, policy, seed):
    """Report the group composition and optionally subsample each matched pair."""
    prediction_groups, reference_groups = defaultdict(list), defaultdict(list)
    for index, row in enumerate(predictions):
        prediction_groups[tuple(row[key] for key in keys)].append(index)
    for index, row in enumerate(references):
        reference_groups[tuple(row[key] for key in keys)].append(index)

    randomizer = random.Random(seed)
    prediction_indices, reference_indices, groups = [], [], []
    for group in sorted(prediction_groups):
        predicted, real = prediction_groups[group], reference_groups[group]
        selected_count = min(len(predicted), len(real))
        if policy == "match_group_counts":
            prediction_indices.extend(randomizer.sample(predicted, selected_count))
            reference_indices.extend(randomizer.sample(real, selected_count))
        groups.append(
            {
                "match": dict(zip(keys, group, strict=True)),
                "n_predictions_available": len(predicted),
                "n_references_available": len(real),
                "n_predictions_used": selected_count
                if policy == "match_group_counts"
                else len(predicted),
                "n_references_used": selected_count
                if policy == "match_group_counts"
                else len(real),
            }
        )

    if policy == "match_group_counts":
        predictions = [predictions[index] for index in sorted(prediction_indices)]
        references = [references[index] for index in sorted(reference_indices)]
    proportions_match = all(
        group["n_predictions_used"] * len(references)
        == group["n_references_used"] * len(predictions)
        for group in groups
    )
    return (
        predictions,
        references,
        {
            "cohort_balance": policy,
            "group_proportions_match": proportions_match,
            "cohort_groups": groups,
        },
    )
