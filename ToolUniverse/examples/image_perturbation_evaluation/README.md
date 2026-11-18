# Independent perturbation image evaluation

Each tool registers a backend named for its operation. The backends share the evaluation workflow in `image_perturbation_metrics_tool.py`.

| Tool | Backend module | Output |
| --- | --- | --- |
| `ImagePerturbation_evaluate_moa` | `image_perturbation_evaluate_moa_tool.py` | Accuracy, Macro-F1, Weighted-F1; per-image/class/condition results |
| `ImagePerturbation_calculate_fid_kid` | `image_perturbation_calculate_fid_kid_tool.py` | FID, KID mean/std; sample counts and per-condition results |

`model_name` is a report label. It never selects a generator, checkpoint,
preprocessing branch or metric. The tools neither import CellFlux nor call any
generator. CellFlux, IMPA, and other models can supply
the same canonical prediction index. Current profiles implement the **BBBC021
96×96 three-channel uint8 contract**; other image modalities need explicit
profiles and compatible evaluators.

An additional explicit `bbbc021_native_inception` profile accepts uniform native
96×96 or 128×128 three-channel cohorts. It preserves
pixels and uses Inception's internal resize, with actual prediction/reference
shapes recorded. Standard/unseen profiles remain strict 96×96. The smoke
profile also accepts these native shapes for integration tests. Predictions
flagged with an unverified condition map are allowed only in the smoke profile.

## Environment and evaluation assets

Use an evaluation environment with a suitable CPU/CUDA PyTorch and torchvision
build; keep these optional dependencies outside ToolUniverse's base environment:

```bash
perturb_tools/models/cellflux-env/bin/python -m pip install \
  -r ToolUniverse/examples/image_perturbation_evaluation/requirements-evaluation.txt

perturb_tools/models/cellflux-env/bin/python \
  ToolUniverse/examples/image_perturbation_evaluation/prepare_assets.py \
  --asset-dir perturb_tools/models/image_evaluation --variant standard
```

The helper verifies the official frozen MoA checkpoint's SHA256 and extracts its
Inception state dict for independent FID/KID calls. `--variant ood` prepares the
separate released OOD classifier. Evaluations themselves download nothing and
never train or replace a classifier when weights are missing.

`bbbc021_labels.json` is an ordered vocabulary: entry *i* is output *i* of the
classifier head. The example helper sorts the original BBBC021 MoA annotations,
following the upstream `np.unique(ANNOT)` convention. Before benchmark claims,
verify this vocabulary against the actual classifier training metadata. The
tools require an explicit map; they never infer it from the evaluated labels.

## MoA evaluation

```python
from tooluniverse.tools import ImagePerturbation_evaluate_moa

result = ImagePerturbation_evaluate_moa(
    predictions="/results/model/generated_images.jsonl",
    model_name="CellFlux",  # arbitrary report label
    evaluation_profile="bbbc021_standard",
    reference_metadata="/data/test_metadata.csv",
    evaluator_checkpoint="/models/moa_standard.pth",
    label_map="/models/bbbc021_labels.json",
    python_executable="/envs/image-evaluation/bin/python",
    output_dir="/results/evaluation",
    split="test",
)
```

Prediction CSV/JSONL requires `sample_id,image_path,condition`. Relative image
paths resolve beside the index. An optional `sha256` is checked; sample IDs and
paths must be unique. Both 96×96×3 HWC uint8 NPY cells and lossless RGB PNG are
accepted. Channel order is preserved. There is no implicit resizing, channel
conversion, random flip, dequantization or contrast normalization.

Ground-truth CSV/JSONL requires:

```csv
condition,moa,split
floxuridine,DNA replication,test
mitoxantrone,DNA replication,test
```

The default `join_on="condition"` requires one unambiguous MoA per condition.
Ground truth comes exclusively from this reference file. For a sample-specific
cohort use `join_on="sample_id"`: reference metadata must additionally contain
`sample_id`; predictions use `reference_sample_id` when provided, otherwise
their own `sample_id`. Row order is never used to align labels. Different cells
may represent the reference and generated populations; there is no pixel-pair
assumption. Dose-dependent labels require distinct condition identifiers or
explicit sample alignment.

The supplied MoA checkpoint must contain frozen Inception weights plus a
2048→512→13 MLP in the released evaluator state format. This is an evaluation
backend contract, independent of the generator architecture.

## FID and KID evaluation

```python
from tooluniverse.tools import ImagePerturbation_calculate_fid_kid

result = ImagePerturbation_calculate_fid_kid(
    predictions="/results/model/generated_images.jsonl",
    model_name="IMPA",
    evaluation_profile="bbbc021_standard",
    reference_images="/data/real_test_images.csv",
    feature_extractor_weights="/models/inception_from_moa_standard.pth",
    python_executable="/envs/image-evaluation/bin/python",
    output_dir="/results/evaluation",
    kid_subset_size=100,
    kid_subsets=100,
)
```

Real-image CSV/JSONL requires `sample_id,image_path,condition,split` and the
selected matching fields. Example:

```csv
sample_id,image_path,condition,batch,split
real_001,images/cell_001.npy,floxuridine,Week9_39282,test
real_002,images/cell_002.npy,floxuridine,Week9_39282,test
```

Default `match_keys=["condition", "batch"]` requires real references for every
generated condition/batch group. References outside those groups are excluded
and counted. `match_keys=["condition"]` explicitly permits cross-batch pooled
references; the chosen policy is written to the report. Cohort sizes may differ
because these are distribution comparisons. Shared image paths are rejected.

The default `cohort_balance="as_provided"` retains both matched cohorts and
reports their available and used counts per match group. When group proportions
differ, the pooled scores reflect both image differences and cohort composition.
Set `cohort_balance="match_group_counts"` to use the smaller available count
on both sides of each group. Selection is seeded and without replacement. The
global and per-condition metrics use the selected cohorts. Use the same group
counts and seed across model comparisons, and set `kid_subset_size` for the
available cohort sizes.

One frozen torch-fidelity Inception-v3 extracts 2048-dimensional features once
per image. TorchMetrics computes FID with float64 covariance arithmetic and
KID with its default polynomial kernel (degree 3, gamma 1/2048, coefficient 1).
KID is unbiased and may be negative. The seed controls subset draws. At least
two images are needed in each cohort and `kid_subset_size` cannot exceed either
cohort size. Small tests must explicitly lower it. Per-condition scores are
marked `insufficient_samples` when counts do not meet this configuration.

## Reports and profiles

Each call creates a unique run directory and returns `metrics`, `report_path`,
`per_condition_path`, counts, asset fingerprints, profile and versions. MoA
additionally returns `predictions_path` (per-image CSV) and `details_path`
(per-class metrics and labeled confusion matrix). Accuracy/F1 are fractions
in [0,1]. Macro-F1 uses the union of true and predicted classes, with undefined
terms set to zero; unused classifier outputs are not added to its denominator.

`bbbc021_standard` and `bbbc021_unseen_compounds` require reference `split=test`.
`bbbc021_smoke` permits explicit integration subsets, such as
`split="unassigned_smoke_subset"`. Profiles declare preprocessing/protocol
intent. They do not independently prove training exclusions, label-map
correctness, official cohort identity or absence of data leakage. Reports keep
`benchmark_verified=false`; compare paper results only after those assets and
cohort details are matched. Never use the smoke subset as a benchmark result.

Use identical reference cohorts, ordered label maps, evaluator weights,
preprocessing, KID parameters and seeds across models. Change the prediction
index and `model_name` when evaluating another generator.

Failed, changed-input and timed-out runs return errors with log paths. No
partial result is returned as successful. Inputs are limited to 50,000 selected
rows per manifest; large KID cohorts retain their feature matrices in CPU RAM.

For a CLI invocation see `run_evaluation.py --help`. Tools load with
`ToolUniverse().load_tools(categories=["image_perturbation_evaluation"])`; SDK
functions and dotted aliases `image_perturbation.evaluate_moa` /
`image_perturbation.evaluate_fid_kid` are registered. Static registry files are
not manually edited; source discovery reads the class decorator.

Sources: [released MoA evaluator](https://github.com/yuhui-zh15/CellFlux/tree/main/moa),
[TorchMetrics FID](https://lightning.ai/docs/torchmetrics/stable/image/frechet_inception_distance),
[TorchMetrics KID](https://lightning.ai/docs/torchmetrics/stable/image/kernel_inception_distance),
[torch-fidelity](https://github.com/toshas/torch-fidelity).

## Worker lifecycle and publication checks

Numeric arguments must be finite. Worker timeouts and caller interruptions kill
and reap the worker process group. A zero exit code alone is insufficient:
result identity, required output files and generated image counts/hashes (or
required evaluation metrics) are checked before returning success. Invalid or
incomplete outputs return `INVALID_OUTPUT` with the run/log paths.

CSV verification flags use explicit `true`/`false` values and are parsed as
booleans, preserving the same validation rules as JSONL. Duplicate CSV headers
are rejected. FID/KID cohorts cannot share an image path or identical file hash.

## Canonical evaluation names and compatibility

Use `ImagePerturbation_evaluate_moa` and
`ImagePerturbation_calculate_fid_kid` as the canonical tools and Python SDK
functions. The old Python imports `image_perturbation_evaluate_moa` and
`image_perturbation_evaluate_fid_kid` remain available. Both old snake-case wire
names and dotted aliases (`image_perturbation.evaluate_moa`,
`image_perturbation.evaluate_fid_kid`) continue to resolve to the canonical tools.
The internal worker operations remain `moa` and `fid_kid`; parameters, frozen
models and metric definitions are unchanged. No generator backend is selected
by either evaluation name or the report-only model_name argument.
