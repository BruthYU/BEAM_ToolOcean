"""Download original expression fixtures; no observations are generated."""

import hashlib
import json
import shutil
import tarfile
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent / "test_data" / "real"
SOURCES = [
    (
        "molecular/adamson_2016.h5ad",
        "https://zenodo.org/api/records/13350497/files/AdamsonWeissman2016_GSM2406675_10X001.h5ad/content",
        "232f7e3756d41602bbe434b50662a76f",
    ),
    (
        "format_reference/pbmc3k_mtx.tar.gz",
        "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
        None,
    ),
    (
        "format_reference/pbmc_1k.h5",
        "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5",
        None,
    ),
]


def main():
    records = []
    for relative, url, md5 in SOURCES:
        path = ROOT / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            temporary = path.with_name(path.name + ".part")
            with requests.get(url, stream=True, timeout=(20, 120)) as response:
                response.raise_for_status()
                with temporary.open("wb") as output:
                    for block in response.iter_content(1024 * 1024):
                        output.write(block)
            temporary.replace(path)
        data = path.read_bytes()
        if md5 and hashlib.md5(data).hexdigest() != md5:
            raise ValueError(f"Published MD5 mismatch: {path}")
        records.append(
            {
                "path": relative,
                "url": url,
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
        )
        print(f"Ready: {relative}", flush=True)
    destination = ROOT / "format_reference/pbmc3k_mtx"
    with tarfile.open(ROOT / SOURCES[1][0]) as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            target = (destination / member.name).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise ValueError("Unsafe archive member")
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(member) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output)
    (ROOT / "expression_sources.json").write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
