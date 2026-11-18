"""Prepare an official frozen MoA head and its independent Inception weights.

Run with the evaluation Python. The ordered example label vocabulary follows
the upstream np.unique ordering of the original BBBC021 MoA annotations.
For benchmark claims, verify it against the training metadata ANNOT vocabulary.
"""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import requests

COMMIT = "041c4c071056896316338872145cef164b07a5fc"
CHECKPOINTS = {
    "standard": (
        "checkpoint.pth",
        108438340,
        "807e7f4c082b37f083d579e0ac2473d459350d7f027be66369ecf55c8d5351f2",
    ),
    "ood": (
        "checkpoint_ood.pth",
        108403820,
        "8e5623c77c3bbca1624d25f8b5fb1492b858ac072b5c8a1058d4bb894422cdd2",
    ),
}
ANNOTATIONS = "https://data.broadinstitute.org/bbbc/BBBC021/BBBC021_v1_moa.csv"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    import torch

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", required=True, type=Path)
    parser.add_argument("--variant", choices=list(CHECKPOINTS), default="standard")
    args = parser.parse_args()
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    filename, size, expected = CHECKPOINTS[args.variant]
    checkpoint = args.asset_dir / f"moa_{args.variant}.pth"
    url = f"https://media.githubusercontent.com/media/yuhui-zh15/CellFlux/{COMMIT}/moa/{filename}"
    if not checkpoint.exists():
        with requests.get(url, stream=True, timeout=(30, 120)) as response:
            response.raise_for_status()
            with checkpoint.with_suffix(".part").open("wb") as stream:
                for block in response.iter_content(4 * 1024 * 1024):
                    stream.write(block)
        if (
            checkpoint.with_suffix(".part").stat().st_size != size
            or sha256(checkpoint.with_suffix(".part")) != expected
        ):
            raise ValueError(
                "Downloaded evaluator does not match its official fingerprint."
            )
        checkpoint.with_suffix(".part").rename(checkpoint)
    if sha256(checkpoint) != expected:
        raise ValueError("Existing evaluator checkpoint SHA256 mismatch.")
    state = torch.load(checkpoint, map_location="cpu", weights_only=True)["model_state"]
    prefix = "feature_extractor.inception."
    inception = {
        key[len(prefix) :]: value
        for key, value in state.items()
        if key.startswith(prefix)
    }
    if not inception:
        raise ValueError("No frozen Inception weights in the evaluator.")
    feature_weights = args.asset_dir / f"inception_from_moa_{args.variant}.pth"
    torch.save(inception, feature_weights)
    response = requests.get(ANNOTATIONS, timeout=60)
    response.raise_for_status()
    annotations = response.content
    (args.asset_dir / "BBBC021_v1_moa.csv").write_bytes(annotations)
    labels = sorted(
        {row["moa"] for row in csv.DictReader(io.StringIO(annotations.decode()))}
    )
    if len(labels) != state["classifier.3.weight"].shape[0] or len(labels) != 13:
        raise ValueError("Annotation vocabulary and classifier dimensions disagree.")
    (args.asset_dir / "bbbc021_labels.json").write_text(
        json.dumps({"labels": labels}, indent=2)
    )
    provenance = {
        "checkpoint_url": url,
        "checkpoint_sha256": expected,
        "source_commit": COMMIT,
        "variant": args.variant,
        "feature_extractor_weights": str(feature_weights.resolve()),
        "feature_weights_sha256": sha256(feature_weights),
        "annotation_source": ANNOTATIONS,
        "annotation_sha256": hashlib.sha256(annotations).hexdigest(),
        "label_order": "sorted original BBBC021 moa vocabulary; verify against original training ANNOT before benchmark claims",
    }
    (args.asset_dir / f"sources_{args.variant}.json").write_text(
        json.dumps(provenance, indent=2)
    )
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
