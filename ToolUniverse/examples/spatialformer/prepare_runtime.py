#!/usr/bin/env python3
"""Download pinned official SpatialFormer code and checkpoints for offline inference.

Run on the networked development host, never as a tool operation on a GPU worker.
Uses only the Python standard library. No training or scientific computation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError
from urllib.parse import urlsplit

COMMIT = "b3f57ea858acc99983c92415f8e56ac521352406"
REPO = "TerminatorJ/Spatialformer"
MODELS = {
    "single": (31146238, 61331527, "fb2c33044244bf21a962e106951769f5"),
    "pair": (31146247, 61331557, "e02503e10a45dd1df577c1853703890d"),
}


def download(url, path):
    """Resume partial transfers and validate byte ranges before appending."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    for attempt in range(6):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "ToolUniverse-SpatialFormer", "Cache-Control": "no-cache"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        # Figshare issues short-lived redirect URLs; do not reuse cached redirects.
        target = url
        if urlsplit(url).hostname == "ndownloader.figshare.com":
            target += ("&" if "?" in target else "?") + f"transfer={time.time_ns()}"
        try:
            try:
                response = urlopen(Request(target, headers=headers), timeout=120)
            except HTTPError as exc:
                if exc.code == 416 and offset:
                    total = exc.headers.get("Content-Range", "").removeprefix("bytes */")
                    if total.isdigit() and int(total) == offset:
                        partial.replace(path)
                        return  # Caller still verifies the published checksum.
                raise
            with response:
                total = None
                if response.status == 206:
                    value = response.headers.get("Content-Range", "")
                    interval, total_text = value.removeprefix("bytes ").split("/")
                    start, end = map(int, interval.split("-"))
                    total = int(total_text)
                    if start != offset or end < start or end >= total:
                        raise ValueError("Server returned an incompatible byte range")
                    mode = "ab" if offset else "wb"
                elif response.status == 200:
                    # Server ignored Range: restart instead of corrupting the file.
                    mode = "wb"
                    if response.headers.get("Content-Length"):
                        total = int(response.headers["Content-Length"])
                else:
                    raise ValueError(f"Unexpected HTTP status: {response.status}")
                with partial.open(mode) as out:
                    shutil.copyfileobj(response, out, 64 * 1024)
                if total is not None and partial.stat().st_size != total:
                    raise ValueError("Incomplete transfer; retained partial file for resume")
            partial.replace(path)
            return
        except Exception:
            if attempt == 5:
                raise
            time.sleep(2 ** attempt)


def digest(path, algorithm="sha256"):
    value = hashlib.new(algorithm)
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--source-only", action="store_true")
    parser.add_argument("--include-single-checkpoint", action="store_true", help="Optional separate single-trained checkpoint; default uses the paired checkpoint for both inference modes.")
    args = parser.parse_args()
    root = Path(args.output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    tree_path = root / "source_tree.json"
    if not tree_path.exists():
        download(f"https://api.github.com/repos/{REPO}/git/trees/{COMMIT}?recursive=1", tree_path)
    tree = json.loads(tree_path.read_text())["tree"]
    entries = [entry for entry in tree if entry["type"] == "blob" and
               (entry["path"].startswith("spatialformer/") or entry["path"] in
                {"LICENSE", "README.md", "setup.py", "requirements.txt", "MANIFEST.in"}) and
               not any(part in entry["path"] for part in ("__pycache__", ".ipynb_checkpoints"))]

    def source_file(entry):
        path = root / "source" / entry["path"]
        def matches():
            if not path.exists() or path.stat().st_size != entry["size"]:
                return False
            h = hashlib.sha1(f"blob {entry['size']}\0".encode())
            h.update(path.read_bytes())
            return h.hexdigest() == entry["sha"]
        if not matches():
            download(f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{entry['path']}", path)
        if not matches():
            raise ValueError(f"Official Git blob checksum mismatch: {entry['path']}")
        return entry["path"]

    with ThreadPoolExecutor(max_workers=6) as executor:
        for name in executor.map(source_file, entries):
            print(f"Verified source: {name}", flush=True)
    print(f"Official source: {root / 'source'}", flush=True)
    if args.source_only:
        return
    package = root / "source/spatialformer"

    def model_bundle(item):
        mode, (article_id, file_id, expected_md5) = item
        bundle = root / mode
        bundle.mkdir(exist_ok=True)
        metadata_path = bundle / "figshare_metadata.json"
        download(f"https://api.figshare.com/v2/articles/{article_id}", metadata_path)
        metadata = json.loads(metadata_path.read_text())
        file = next(f for f in metadata["files"] if f["id"] == file_id)
        if file["computed_md5"] != expected_md5:
            raise ValueError("Published checkpoint checksum changed; review before updating.")
        checkpoint = bundle / "model.ckpt"
        if not checkpoint.exists() or digest(checkpoint, "md5") != expected_md5:
            download(file["download_url"], checkpoint)
        if checkpoint.stat().st_size != file["size"] or digest(checkpoint, "md5") != expected_md5:
            raise ValueError(f"Checkpoint verification failed: {mode}")
        resources = {
            "checkpoint": checkpoint,
            # Both upstream configs have identical architecture; embed_data sets mode at inference.
            "config": package / "config/_config_train_large_pair.json",
            "tokenizer": package / "tokenizer/tokenv5.json",
            "gene_medians": package / "data/gene_median.pkl",
        }
        manifest = {"model_id": f"spatialformer-{mode}-71slides-5k", "mode": mode, "training_mode": mode,
                    "source_commit": COMMIT, "package_version": "0.1.8",
                    "source_url": metadata["url_public_html"], "checkpoint_license": metadata["license"],
                    "runtime_source_sha256": {str(p.relative_to(package)): digest(p) for p in package.rglob("*.py")},
                    "spatial_embeddings_sha256": digest(package / "spatial_embeddings/xenium_5k_pandavid_gene_embeddings.pkl"),
                    "files": {key: {"path": str(path.relative_to(bundle)) if path.is_relative_to(bundle)
                                     else str(Path("../source/spatialformer") / path.relative_to(package)),
                                    "sha256": digest(path)} for key, path in resources.items()}}
        (bundle / "manifest.json").write_text(json.dumps(manifest, indent=2))
        print(f"Verified model bundle: {bundle}", flush=True)
    with ThreadPoolExecutor(max_workers=2) as executor:
        selected = MODELS.items() if args.include_single_checkpoint else [("pair", MODELS["pair"])]
        list(executor.map(model_bundle, selected))
    if not args.include_single_checkpoint:
        # Official embed_data accepts the paired pretrained model in single-input mode.
        # Reuse the same checkpoint without copying or changing any trained weights.
        alias = root / "single"
        alias.mkdir(exist_ok=True)
        manifest = json.loads((root / "pair/manifest.json").read_text())
        manifest["mode"] = "single"
        manifest["training_mode"] = "pair"
        manifest["files"]["checkpoint"]["path"] = "../pair/model.ckpt"
        manifest["note"] = "Single-cell feature extraction with the paired pretrained checkpoint; not a single-trained checkpoint."
        (alias / "manifest.json").write_text(json.dumps(manifest, indent=2))
        print(f"Single-input alias of paired checkpoint: {alias}", flush=True)


if __name__ == "__main__":
    main()
