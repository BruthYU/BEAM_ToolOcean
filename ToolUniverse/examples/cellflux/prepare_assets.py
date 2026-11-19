"""Download official CellFlux assets and a small real BBBC021 input subset.

The subset preserves original pixels and annotations but is NOT an official
benchmark split. HTTP range access avoids downloading the full 12.5 GB bundle.
"""

import argparse
import csv
import hashlib
import io
import json
import shutil
import struct
import tarfile
import zipfile
import zlib
from pathlib import Path

import requests

CELLFLUX_COMMIT = "041c4c071056896316338872145cef164b07a5fc"
ZENODO = "https://zenodo.org/records/8307629/files/IMPA_reproducibility.zip?download=1"
WEIGHTS = "https://huggingface.co/suyc21/CellFlux/resolve/ac33c8d5d74b06062bb3de571ff14d47314eb27e/checkpoints/cellflux/bbbc021/checkpoint.pth"
WEIGHTS_SHA256 = "a3bf1e8471a6c01565a866faa1f5cf59899165d37dd7d46aeb3ea8552ec59445"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def get(url, byte_range=None):
    headers = {"Range": f"bytes={byte_range[0]}-{byte_range[1]}"} if byte_range else {}
    with requests.get(url, headers=headers, timeout=(30, 120), stream=True) as response:
        response.raise_for_status()
        if byte_range and (
            response.status_code != 206
            or response.headers.get("Content-Range", "").split("/")[0]
            != f"bytes {byte_range[0]}-{byte_range[1]}"
        ):
            raise RuntimeError("Server did not return the requested HTTP range.")
        return response.content


def prepare_subset(destination):
    import numpy as np
    import pandas as pd

    destination.mkdir(parents=True, exist_ok=True)
    # Fixed offsets belong only to the immutable Zenodo record 8307629.
    offset = 2816169106
    compressed = get(ZENODO, (offset, offset + 2 * 1024 * 1024 - 1))
    header = struct.unpack_from("<4s5H3L2H", compressed)
    begin = 30 + header[-2] + header[-1]
    member_name = compressed[30 : 30 + header[-2]].decode()
    if (
        header[0] != b"PK\x03\x04"
        or member_name != "IMPA_reproducibility/datasets/bbbc021_all.tar"
    ):
        raise RuntimeError("Unexpected immutable archive layout.")
    gzip_prefix = zlib.decompressobj(-15).decompress(compressed[begin:])
    tar_prefix = zlib.decompressobj(31).decompress(gzip_prefix)
    annotation_url = "https://data.broadinstitute.org/bbbc/BBBC021/BBBC021_v1_image.csv"
    annotation = get(annotation_url)
    (destination / "BBBC021_v1_image.csv").write_bytes(annotation)
    table = pd.read_csv(io.BytesIO(annotation)).set_index(
        ["TableNumber", "ImageNumber"]
    )
    records = []
    with tarfile.open(fileobj=io.BytesIO(tar_prefix), mode="r|") as archive:
        try:
            for member in archive:
                if not member.isfile() or not member.name.endswith(".npy"):
                    continue
                try:
                    pixels = archive.extractfile(member).read()
                except tarfile.ReadError:
                    break  # The last member in a prefix may be incomplete.
                if len(pixels) != member.size:
                    break
                image = np.load(io.BytesIO(pixels), allow_pickle=False)
                if image.shape != (96, 96, 3) or image.dtype != np.uint8:
                    raise ValueError("Unexpected official image shape or dtype.")
                relative = Path(member.name)
                numbers = relative.stem.split("_")
                key = (int(numbers[0]), int(numbers[1]))
                row = table.loc[key]
                if isinstance(row, pd.DataFrame):
                    raise TypeError("Original image annotations are ambiguous.")
                plate = str(row["Image_Metadata_Plate_DAPI"])
                if plate != "_".join(relative.parts[1:3]):
                    raise ValueError(
                        "Cell crop and original plate annotation disagree."
                    )
                path = destination / "images" / relative.name
                path.parent.mkdir(exist_ok=True)
                path.write_bytes(pixels)
                records.append(
                    {
                        "sample_id": plate + "_" + relative.stem,
                        "image_path": str(path.relative_to(destination)),
                        "batch": plate,
                        "condition": str(row["Image_Metadata_Compound"]),
                        "dose": float(row["Image_Metadata_Concentration"]),
                        "split": "unassigned_smoke_subset",
                        "source_member": member.name,
                        "sha256": digest(pixels),
                    }
                )
        except tarfile.ReadError:
            pass
    if not records or not any(row["condition"] == "DMSO" for row in records):
        raise ValueError("Downloaded prefix contains no usable real DMSO controls.")
    for name, rows in (
        ("samples.csv", records),
        ("controls.csv", [row for row in records if row["condition"] == "DMSO"]),
    ):
        with (destination / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(rows)
    provenance = {
        "source": ZENODO,
        "range": [offset, offset + len(compressed) - 1],
        "range_sha256": digest(compressed),
        "annotation_source": annotation_url,
        "annotation_sha256": digest(annotation),
        "n_cells": len(records),
        "n_controls": sum(row["condition"] == "DMSO" for row in records),
        "benchmark_eligible": False,
        "note": "Archive-prefix smoke subset; no official train/test split assigned. Original stored channel order is preserved.",
    }
    (destination / "sources.json").write_text(json.dumps(provenance, indent=2))
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--skip-weights", action="store_true")
    args = parser.parse_args()
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    repository = args.asset_dir / "CellFlux"
    if not repository.exists():
        url = f"https://codeload.github.com/yuhui-zh15/CellFlux/zip/{CELLFLUX_COMMIT}"
        code = get(url)
        with zipfile.ZipFile(io.BytesIO(code)) as archive:
            for member in archive.infolist():
                relative = Path(*Path(member.filename).parts[1:])
                if member.is_dir() or relative == Path("."):
                    continue
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("Unsafe archive member.")
                target = repository / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.read(member))
        (repository / "SOURCE.json").write_text(
            json.dumps(
                {"url": url, "commit": CELLFLUX_COMMIT, "sha256": digest(code)},
                indent=2,
            )
        )
    embeddings_url = (
        "https://raw.githubusercontent.com/theislab/IMPA/main/embeddings/csv/emb_fp.csv"
    )
    embeddings = get(embeddings_url)
    (args.asset_dir / "emb_fp.csv").write_bytes(embeddings)
    if not args.skip_weights:
        path = args.asset_dir / "checkpoint.pth"
        if not path.exists():
            with requests.get(WEIGHTS, stream=True, timeout=(30, 120)) as response:
                response.raise_for_status()
                with path.with_suffix(".part").open("wb") as stream:
                    shutil.copyfileobj(response.raw, stream)
            if path.with_suffix(".part").stat().st_size != 899852463:
                raise ValueError("Unexpected checkpoint size.")
            path.with_suffix(".part").rename(path)
        checksum = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                checksum.update(block)
        if checksum.hexdigest() != WEIGHTS_SHA256:
            raise ValueError("Checkpoint SHA256 does not match the official asset.")
    print(json.dumps(prepare_subset(args.data_dir), indent=2))
    (args.asset_dir / "sources.json").write_text(
        json.dumps(
            {
                "cellflux_commit": CELLFLUX_COMMIT,
                "checkpoint_url": WEIGHTS,
                "checkpoint_sha256": WEIGHTS_SHA256,
                "embeddings_url": embeddings_url,
                "embeddings_sha256": digest(embeddings),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
