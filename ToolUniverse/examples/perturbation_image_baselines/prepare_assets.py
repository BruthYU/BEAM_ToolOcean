"""Download pinned pretrained IMPA inference assets."""

import argparse
import concurrent.futures
import hashlib
import io
import json
import shutil
import zipfile
from pathlib import Path

import requests

HF_REVISION = "ac33c8d5d74b06062bb3de571ff14d47314eb27e"
REPOSITORIES = {
    "IMPA": ("theislab/IMPA", "13eb82f582268b43ffe5173e75c3ba8b3f42c8d8"),
}
ASSETS = {
    "impa/bbbc021/nets.ckpt": (
        212153202,
        "9e66e7ee61f5a78afc7433224f66b70ea785f3f1a3b8a143ec6e1a4a0b120bb5",
    ),
    "impa/bbbc021/embeddings.ckpt": (
        107854,
        "32564991acd6f19fd3d39f41f7a7ff8c5c0e951e6ee3d590ac1a9195127ad292",
    ),
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(path, url, size, expected):
    if (
        path.exists()
        and path.stat().st_size == size
        and (not expected or sha256(path) == expected)
    ):
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if size < 1024 * 1024:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        path.write_bytes(response.content)
    else:
        response = requests.head(url, allow_redirects=True, timeout=60)
        response.raise_for_status()
        resolved = response.url
        block_size = 8 * 1024 * 1024

        def part(index):
            start, end = index * block_size, min(size, (index + 1) * block_size) - 1
            target = path.with_name(path.name + f".part{index}")
            if target.exists() and target.stat().st_size == end - start + 1:
                return target
            with requests.get(
                resolved,
                headers={"Range": f"bytes={start}-{end}"},
                stream=True,
                timeout=(30, 120),
            ) as reply:
                reply.raise_for_status()
                if (
                    reply.status_code != 206
                    or reply.headers.get("Content-Range", "").split("/")[0]
                    != f"bytes {start}-{end}"
                ):
                    raise ValueError("Server did not return the requested byte range.")
                with target.open("wb") as stream:
                    for data in reply.iter_content(1024 * 1024):
                        stream.write(data)
            if target.stat().st_size != end - start + 1:
                raise ValueError("Incomplete download chunk.")
            return target

        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            parts = list(pool.map(part, range((size + block_size - 1) // block_size)))
        temporary = path.with_name(path.name + ".assembling")
        with temporary.open("wb") as stream:
            for item in parts:
                with item.open("rb") as source:
                    shutil.copyfileobj(source, stream)
        temporary.rename(path)
        for item in parts:
            item.unlink()
    if path.stat().st_size != size or (expected and sha256(path) != expected):
        raise ValueError(f"Downloaded asset size/SHA256 mismatch: {path}")


def main():
    import numpy as np
    import pandas as pd
    import torch

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.asset_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    for name, (repo, commit) in REPOSITORIES.items():
        target = root / name
        if target.exists():
            metadata = json.loads((target / "SOURCE.json").read_text())
            if metadata.get("commit") != commit:
                raise ValueError(
                    f"Existing {name} source does not match the pinned revision."
                )
            continue
        url = f"https://codeload.github.com/{repo}/zip/{commit}"
        response = requests.get(url, timeout=(30, 120))
        response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            for member in archive.infolist():
                relative = Path(*Path(member.filename).parts[1:])
                if member.is_dir() or relative == Path("."):
                    continue
                if ".." in relative.parts or relative.is_absolute():
                    raise ValueError("Unsafe source archive path.")
                path = target / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(archive.read(member))
        (target / "SOURCE.json").write_text(
            json.dumps(
                {
                    "url": url,
                    "commit": commit,
                    "sha256": hashlib.sha256(response.content).hexdigest(),
                },
                indent=2,
            )
        )
    checksums = {}
    for name, (size, expected) in ASSETS.items():
        url = f"https://huggingface.co/suyc21/CellFlux/resolve/{HF_REVISION}/checkpoints/{name}?download=true"
        path = root / "weights" / name
        download(path, url, size, expected)
        checksums[name] = sha256(path)
    fingerprints = pd.read_csv(root / "IMPA/embeddings/csv/emb_fp.csv", index_col=0)
    frozen = torch.load(
        root / "weights/impa/bbbc021/embeddings.ckpt",
        map_location="cpu",
        weights_only=True,
    )["embedding_matrix"]["weight"].numpy()
    labels = []
    values = fingerprints.to_numpy(dtype=np.float32)
    for row in frozen:
        hits = np.flatnonzero(np.all(values == row, axis=1))
        if len(hits) != 1:
            raise ValueError("Frozen IMPA rows do not match unique named fingerprints.")
        labels.append(str(fingerprints.index[hits[0]]))
    if labels != sorted(labels) or len(set(labels)) != 26:
        raise ValueError("Unexpected released training vocabulary.")
    shutil.copyfile(root / "IMPA/embeddings/csv/emb_fp.csv", root / "emb_fp.csv")
    (root / "impa_training_conditions.json").write_text(
        json.dumps(
            {
                "labels": labels,
                "verified": True,
                "provenance": "Exact equality between frozen IMPA embedding rows and named Morgan fingerprints.",
            },
            indent=2,
        )
    )
    summary = {
        "hf_revision": HF_REVISION,
        "sha256": checksums,
        "impa_condition_mapping_verified": True,
    }
    (root / "assets.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
