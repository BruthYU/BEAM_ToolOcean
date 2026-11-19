"""Standalone CellFlux worker; optional model dependencies stay out of SDK imports."""

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_cell(path, expected_sha256=None):
    import numpy as np
    from PIL import Image

    with Path(path).open("rb") as stream:
        data = stream.read(1024 * 1024 + 1)
    if len(data) > 1024 * 1024:
        raise ValueError("Control image exceeds the single-cell size limit.")
    if (
        expected_sha256 is not None
        and hashlib.sha256(data).hexdigest() != expected_sha256
    ):
        raise ValueError("A control image changed after validation.")
    if Path(path).suffix.lower() == ".npy":
        if Path(path).stat().st_size > 131072:
            raise ValueError("Input NPY exceeds the size of a 96x96x3 uint8 cell.")
        image = np.load(io.BytesIO(data), allow_pickle=False)
    else:
        with Image.open(io.BytesIO(data)) as source:
            if (
                source.mode != "RGB"
                or source.format != "PNG"
                or source.size != (96, 96)
            ):
                raise ValueError("PNG input must already have three ordered channels.")
            image = np.asarray(source)
    if image.shape != (96, 96, 3) or image.dtype != np.uint8:
        raise ValueError(
            f"Expected HWC uint8 (96,96,3), got {image.shape} {image.dtype}"
        )
    return np.array(image, copy=True)


def main(request_path):
    import numpy as np
    import pandas as pd
    import torch
    from flow_matching.solver.ode_solver import ODESolver
    from flow_matching.utils import ModelWrapper
    from PIL import Image

    job = json.loads(Path(request_path).read_text())
    repository = Path(job["model_dir"])
    for source in job.get("source_files", []):
        if sha256(source["path"]) != source["sha256"]:
            raise ValueError("A model/input source changed after validation.")
    sys.path.insert(0, str(repository))
    from models.model_configs import instantiate_model
    from training.edm_time_discretization import get_time_discretization

    device_name = job.get("device", "auto")
    if device_name == "auto":
        device_name = "cuda" if torch.cuda.is_available() else "cpu"
    if device_name == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA was requested but no CUDA device is available.")
    device = torch.device(device_name)
    torch.set_num_threads(job.get("cpu_threads", 4))
    seed = job.get("seed", 42)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    generator = np.random.default_rng(seed)
    embeddings = pd.read_csv(job["embeddings_path"], index_col=0)
    if not embeddings.index.is_unique or embeddings.shape[1] != 1024:
        raise ValueError(
            "Expected unique compound rows and 1024 ordered fingerprint features."
        )
    missing = set(job["compounds"]) - set(embeddings.index)
    if missing:
        raise ValueError(
            f"Compounds missing from the supplied encoding: {sorted(missing)}"
        )
    values = embeddings.loc[job["compounds"]].to_numpy(dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError("Condition encoding contains non-finite values.")
    # Validate all selected pixels and fingerprints before loading the large model.
    for row in job["controls"]:
        load_cell(row["image_path"], row["sha256"])
    for name, key in (
        ("checkpoint", "checkpoint_path"),
        ("embeddings", "embeddings_path"),
    ):
        if sha256(job[key]) != job["input_fingerprints"][name]:
            raise ValueError(f"{name} changed after validation.")
    print(f"Loading official BBBC021 EMA model on {device}", flush=True)
    # Official checkpoints store their training arguments in a namespace object.
    # Keep weights_only loading; do not fall back to unrestricted pickle loading.
    with torch.serialization.safe_globals([argparse.Namespace, SimpleNamespace]):
        checkpoint = torch.load(
            job["checkpoint_path"], map_location="cpu", weights_only=True
        )
    model = instantiate_model("bbbc021", is_discrete=False, use_ema=True)
    model.load_state_dict(checkpoint["model"], strict=True)
    epoch = checkpoint.get("epoch")
    del checkpoint
    model.to(device)
    model.train(
        False
    )  # Official EMA wrapper copies shadow parameters on this transition.

    class GuidedVelocity(ModelWrapper):
        def __init__(self, network, embedding, scale):
            super().__init__(network)
            self.embedding = embedding
            self.scale = scale
            self.calls = 0

        def forward(self, x, t, **kwargs):
            self.calls += 1
            times = torch.zeros(x.shape[0], device=x.device) + t
            mixed_precision = (
                torch.autocast("cuda", dtype=torch.float16)
                if device.type == "cuda"
                else contextlib.nullcontext()
            )
            with mixed_precision:
                conditional = self.model(
                    x, times, extra={"concat_conditioning": self.embedding}
                )
                if self.scale:
                    unconditional = self.model(x, times, extra={})
                    conditional = (
                        1 + self.scale
                    ) * conditional - self.scale * unconditional
            return conditional.float()

    run_dir = Path(job["run_dir"])
    image_dir = run_dir / "images"
    image_dir.mkdir()
    records = []
    count = job.get("samples_per_condition", 8)
    batch_size = job.get("batch_size", 4)
    steps = job.get("sampling_steps", 50)
    grid = get_time_discretization(nfes=steps).to(device)
    start = time.monotonic()
    with torch.inference_mode():
        for batch in sorted({row["batch"] for row in job["controls"]}):
            controls = [row for row in job["controls"] if row["batch"] == batch]
            for compound in job["compounds"]:
                positions = generator.choice(
                    len(controls), size=count, replace=count > len(controls)
                )
                condition = torch.tensor(
                    embeddings.loc[compound].to_numpy(dtype=np.float32), device=device
                )
                for offset in range(0, count, batch_size):
                    selected = [
                        controls[i] for i in positions[offset : offset + batch_size]
                    ]
                    images = np.stack(
                        [
                            load_cell(row["image_path"], row["sha256"])
                            for row in selected
                        ]
                    )
                    x = torch.tensor(
                        images, dtype=torch.float32, device=device
                    ).permute(0, 3, 1, 2)
                    # Match CustomTransform's uint8 dequantization and [-1,1] scaling.
                    x = ((x + torch.rand_like(x)) / 255.0 - 0.5) / 0.5
                    initial = x + torch.randn_like(x) * job.get("noise_level", 1.0)
                    velocity = GuidedVelocity(
                        model,
                        condition.expand(len(selected), -1),
                        job.get("cfg_scale", 0.2),
                    )
                    prediction = ODESolver(velocity_model=velocity).sample(
                        x_init=initial,
                        time_grid=grid,
                        method="heun2",
                        step_size=None,
                        return_intermediates=False,
                        atol=1e-5,
                        rtol=1e-5,
                    )
                    if not torch.isfinite(prediction).all():
                        raise ValueError("Model produced non-finite pixels.")
                    # Match torchvision.save_image's clipping and rounding.
                    pixels = ((prediction * 0.5 + 0.5) * 255 + 0.5).clamp(0, 255)
                    pixels = pixels.permute(0, 2, 3, 1).to("cpu", torch.uint8).numpy()
                    for row, image in zip(selected, pixels):
                        image_id = f"cell_{len(records):06d}"
                        path = image_dir / f"{image_id}.png"
                        Image.fromarray(image).save(path)
                        records.append(
                            {
                                "sample_id": image_id,
                                "image_path": str(path),
                                "condition": compound,
                                "batch": batch,
                                "control_sample_id": row["sample_id"],
                                "control_image_path": row["image_path"],
                                "control_sha256": row["sha256"],
                                "sha256": sha256(path),
                                "shape": [96, 96, 3],
                                "axes": "YXC",
                                "dtype": "uint8",
                                "solver_calls": velocity.calls,
                            }
                        )
                    print(
                        f"{batch}: {compound}: {min(offset + batch_size, count)}/{count}",
                        flush=True,
                    )
    for source in job.get("source_files", []):
        if sha256(source["path"]) != source["sha256"]:
            raise ValueError("A model/input source changed during inference.")
    index = run_dir / "generated_images.jsonl"
    index.write_text("".join(json.dumps(row) + "\n" for row in records))
    commit = (
        (
            subprocess.run(
                ["git", "-C", str(repository), "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                check=False,
            ).stdout.strip()
            or None
        )
        if (repository / ".git").exists()
        else None
    )
    source_record = repository / "SOURCE.json"
    if not commit and source_record.is_file():
        commit = json.loads(source_record.read_text()).get("commit")
    source_hashes = {
        name: sha256(repository / name)
        for name in (
            "models/model_configs.py",
            "models/unet.py",
            "models/nn.py",
            "models/ema.py",
            "training/edm_time_discretization.py",
        )
    }
    result = {
        "schema_version": "1.0",
        "run_id": job["run_id"],
        "model": "CellFlux",
        "dataset": "BBBC021",
        "n_images": len(records),
        "image_index": str(index),
        "image_dir": str(image_dir),
        "run_manifest": str(run_dir / "run.json"),
        "device": str(device),
        "compounds": job["compounds"],
        "sampling_steps": steps,
        "ema": True,
        "solver": "heun2",
        "time_schedule": "official_edm",
        "seed": seed,
        "cfg_scale": job.get("cfg_scale", 0.2),
        "noise_level": job.get("noise_level", 1.0),
        "checkpoint_epoch": epoch,
        "input_fingerprints": job["input_fingerprints"],
        "source_commit": commit,
        "source_sha256": source_hashes,
        "inference_seconds": time.monotonic() - start,
        "versions": {
            key: importlib.metadata.version(key)
            for key in ("torch", "numpy", "flow-matching", "torchdiffeq")
        },
        "interpretation": "Unpaired condition-level population predictions, not observed trajectories.",
        "preview": records[:5],
    }
    Path(result["run_manifest"]).write_text(json.dumps(result, indent=2))
    (run_dir / "result.json").write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
