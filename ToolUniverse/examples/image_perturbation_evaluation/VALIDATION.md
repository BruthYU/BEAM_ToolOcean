# Evaluation tool validation — 2026-11-12

Both independent tools were invoked through ToolUniverse and the actual Python
SDK with real images and a released frozen evaluator. No generator is loaded
by either evaluator; only the prediction manifests and evaluator assets are used.

| Check | Result |
| --- | --- |
| Registration | Both tool names resolve to `ImagePerturbationEvaluationTool` |
| MoA | One existing 50-step CellFlux prediction evaluated successfully |
| FID/KID | Two genuine 2-step CellFlux predictions versus seven same-plate floxuridine cells |
| Model-name independence | Repeating both evaluations with `alternate_report_label` produced identical metrics |
| Ground truth | Original BBBC021 compound/MoA annotations; explicit ordered 13-class map |
| Split | `bbbc021_smoke`, `unassigned_smoke_subset` |
| KID | Explicit subset size 2, 10 subsets, seed 42 |
| Device | CPU; GPU execution not tested |
| Tests | 40 passed across evaluator, CellFlux and existing loader SDK tests |
| Static checks | Ruff and format checks passed |

The single MoA prediction was classified as Kinase inhibitors; its reference
label is DNA replication. Consequently this one-image smoke run has Accuracy,
Macro-F1 and Weighted-F1 of 0.0. This observation does not estimate overall model
accuracy and has not been compared with paper results.

The tiny FID/KID run returned FID 332.21593296003647, KID mean
0.2847291827201843 and KID std 0.04752727970480919. These values only verify the
calculation path: the generated cohort has two images, uses low-step sampling,
and is not the official benchmark cohort. They must not be reported as
benchmark performance.

The example label order follows sorted original BBBC021 MoA names, consistent
with upstream `np.unique` ordering. Original training `ANNOT` vocabulary and
complete official test-set identity still need verification for benchmark use.
Every report therefore records `benchmark_verified=false`.

Official standard evaluator checkpoint:

```text
bytes: 108438340
sha256: 807e7f4c082b37f083d579e0ac2473d459350d7f027be66369ecf55c8d5351f2
source: yuhui-zh15/CellFlux@041c4c071056896316338872145cef164b07a5fc
```

Environment: PyTorch 2.8.0, TorchMetrics 1.7.4, torch-fidelity 0.3.0,
scikit-learn 1.9.1. The worker path setup was checked to avoid shadowing
Python's standard-library `profile` module during torchvision imports.

Local artifacts (intentionally not versioned):

- `perturb_tools/models/image_evaluation/`: checkpoint, independent Inception weights, labels and provenance.
- `perturb_tools/test_data/real/cellflux_bbbc021/evaluation_reference.csv`: canonical real references/labels.
- `perturb_tools/test_runs/image_evaluation/sdk_moa_model_name_check.json`: successful MoA SDK result.
- `perturb_tools/test_runs/image_evaluation/sdk_fid_kid_model_name_check.json`: successful FID/KID SDK result.
- `perturb_tools/test_runs/image_evaluation/validation_summary.json`: verification summary and report links.

Incomplete diagnostic runs from the initial worker import-path issue are
retained with error logs; the successful report paths are listed in the summary.

## Follow-up verification — 2026-11-18

After the rebase to origin/main `bf20982`, both real evaluation workers passed
on two new IMPA predictions and original BBBC021 references. Changing only the
report model_name preserved the MoA metrics. Controlled corruption of a copy
of the released classifier produced non-finite logits and was rejected before
publishing any metrics report. This check used real images and the actual
frozen network, with one deliberately corrupted classifier bias.

CSV true/false provenance flags now receive the same interpretation as JSONL;
invalid flags/duplicate headers are rejected. Prediction/reference files cannot
share a content hash under different paths. Worker outputs are checked before
success, including finite numeric metrics and prediction counts. These remain
small, unassigned integration cohorts, not a verified benchmark.

Local evidence: `perturb_tools/test_runs/tool_batch_quality/18f343fec00d44b7b5a8353b7ead2e68/summary.json`
and `perturb_tools/test_runs/tool_batch_invalid_classifier/`.

## Canonical evaluation name migration

Primary names are now `ImagePerturbation_evaluate_moa` and
`ImagePerturbation_calculate_fid_kid`. Old snake-case/dotted wire names and old
SDK imports remain compatible. Internal operation dispatch explicitly supports
the canonical names and pre-migration configs; worker operations are unchanged.

75 targeted regression tests passed, covering SDK generation/exports/signatures,
wire aliases, legacy configs, request operation and CLI selection. Two optional
metric-library unit checks were skipped in the SDK environment; the real metric
workers ran in their dedicated evaluation environment.

Canonical and old SDK names both evaluated the same two actual pretrained IMPA
predictions against original BBBC021 annotations/images. MoA metrics and seeded
FID/KID values matched exactly under both names. This is a small CPU integration
check and does not certify benchmark scores.

Local evidence: `perturb_tools/test_runs/evaluation_name_validation/740eb653185c4b1eaa9b926b0f22a249/summary.json`.
