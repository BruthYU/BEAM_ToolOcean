# CellFlux BBBC021 morphology prediction

`CellFlux_predict_morphology` (type `CellFluxTool`) runs the official
pretrained BBBC021 EMA model. Model dependencies run in a separate process and
are not imported during ToolUniverse discovery. `_lazy_registry_static.py` does
not need modification: source installations discover the decorated class by AST.

This first backend predicts unpaired treated-cell populations from DMSO cells
and the official 1024-dimensional compound fingerprints. It does not infer a
physical time course, dose response, or a particular cell's observed future.
An embedding row being present does not establish reliable generalization to
that compound. Benchmark and held-out claims require the original splits and
the matching training configuration, checkpoint and evaluator.

## Prepare the model and real input data

From the directory containing `ToolUniverse/` and `perturb_tools/`:

```bash
python -m venv --system-site-packages perturb_tools/models/cellflux-env
perturb_tools/models/cellflux-env/bin/python -m pip install \
  -r ToolUniverse/examples/cellflux/requirements-inference.txt

python ToolUniverse/examples/cellflux/prepare_assets.py \
  --asset-dir perturb_tools/models/cellflux_assets \
  --data-dir perturb_tools/test_data/real/cellflux_bbbc021
```

Install the PyTorch build appropriate for the target CPU/CUDA host. A clean
environment can be used instead of inheriting installed packages. The supplied
setup inherits the existing PyTorch installation to avoid downloading it again.

The preparation script downloads a pinned official CellFlux source archive,
the approximately 900 MB BBBC021 checkpoint, IMPA's `emb_fp.csv`, and a small
real cell subset from the immutable IMPA Zenodo record. It joins original
BBBC021 image annotations using TableNumber and ImageNumber and verifies the
plate identifiers. Original NPY pixels and channel order are preserved. The
subset contains DMSO controls and treated cells; `controls.csv` contains only
DMSO. Source URLs and hashes are saved alongside assets/data.

The subset uses an HTTP range of the archive, **not the official benchmark
train/test split**. It is suitable for integration checks, not paper score
comparison. `samples.csv` labels its split `unassigned_smoke_subset`. Obtain the
complete prepared dataset and official metadata separately for benchmark work.

## Invoke the registered tool

```bash
PYTHONPATH=ToolUniverse/src python ToolUniverse/examples/cellflux/run_prediction.py \
  --asset-dir perturb_tools/models/cellflux_assets \
  --data-dir perturb_tools/test_data/real/cellflux_bbbc021 \
  --python "$PWD/perturb_tools/models/cellflux-env/bin/python" \
  --output-dir perturb_tools/test_runs/cellflux \
  --compound floxuridine --samples 1 --sampling-steps 50 --device auto
```

For a quick plumbing check use `--sampling-steps 2`; do not compare its images
or scores with the official 50-step setup. CPU is supported but slower. CUDA is
selected automatically when available; explicitly requesting unavailable CUDA
returns an error. A trained model is always required; there is no synthetic
fallback. Calls perform no asset downloads.

The SDK also exports:

```python
from tooluniverse.tools import CellFlux_predict_morphology

result = CellFlux_predict_morphology(
    control_manifest="/data/controls.csv",
    compounds=["floxuridine"],
    model_dir="/models/CellFlux",
    checkpoint_path="/models/checkpoint.pth",
    embeddings_path="/models/emb_fp.csv",
    python_executable="/envs/cellflux/bin/python",
    output_dir="/results/cellflux",
    samples_per_condition=8,
    sampling_steps=50,
)
```

## Input and output contracts

The control CSV requires `sample_id,image_path,batch,condition`. Paths may be
absolute or relative to the CSV directory. Every row must be DMSO, identifiers
must be unique, and images must already be **96×96×3 HWC uint8** NPY or lossless
RGB PNG in official stored channel order. Full microscopy fields, uint16 arrays,
resized previews and independently contrast-normalized images are rejected or
outside this contract. No automatic cropping or channel reassignment occurs.

Controls are selected separately within each supplied batch. Optional `batch`
filters one batch; otherwise each requested compound is generated for each
batch. `samples_per_condition` counts images per compound per batch; controls
are resampled when needed. The model receives compound fingerprints, not dose
or MoA labels. Selected controls and their hashes are recorded per prediction.

A successful response provides `image_dir`, `image_index`, `run_manifest`,
`n_images`, `device`, and a bounded preview. Each unique run contains:

- `images/cell_*.png`: generated images in original channel order.
- `generated_images.jsonl`: condition, batch, source control, output paths/hashes.
- `run.json`: model/data fingerprints, versions, seed, solver, schedule and EMA.
- `request.json`, `worker.log`: exact request and execution diagnostics.

The image directory can subsequently be indexed with
`image_perturbation_load_images(source=..., source_type="directory", ...)`.
Timeouts and worker failures return `status="error"`; partial files are retained
for diagnosis and are not published as a successful run. Cache is disabled by
default because calls create new run artifacts.

## Sources

- [CellFlux official implementation](https://github.com/yuhui-zh15/CellFlux),
  source commit `041c4c071056896316338872145cef164b07a5fc`.
- [Official model weights](https://huggingface.co/suyc21/CellFlux).
- [IMPA prepared data](https://zenodo.org/records/8307629) and
  [compound fingerprints](https://github.com/theislab/IMPA/blob/main/embeddings/csv/emb_fp.csv).
- [BBBC021 original data and annotations](https://bbbc.broadinstitute.org/BBBC021).

Keep upstream license notices with downloaded code. Some upstream model source
files carry their own notices in addition to the repository-level license.

## Worker lifecycle and publication checks

Numeric arguments must be finite. Worker timeouts and caller interruptions kill
and reap the worker process group. A zero exit code alone is insufficient:
result identity, required output files and generated image counts/hashes (or
required evaluation metrics) are checked before returning success. Invalid or
incomplete outputs return `INVALID_OUTPUT` with the run/log paths.

CSV verification flags use explicit `true`/`false` values and are parsed as
booleans, preserving the same validation rules as JSONL. Duplicate CSV headers
are rejected. FID/KID cohorts cannot share an image path or identical file hash.

## Tool name compatibility

The canonical tool/SDK name is `CellFlux_predict_morphology`.
`image_perturbation_predict_morphology` and
`image_perturbation.predict_morphology` remain wire aliases; the old Python SDK
import is also retained. New code should use the model-prefixed name. Arguments,
outputs and model behavior are unchanged.
