# Integration validation — 2026-11-12

Verified locally with the actual official checkpoint and real preprocessed
BBBC021 cells. No generated fixture or substitute model was used for inference.

| Check | Result |
| --- | --- |
| Checkpoint | 899,852,463 bytes; SHA256 matches the official asset |
| Source | CellFlux `041c4c071056896316338872145cef164b07a5fc` |
| Input data | 86 original HWC uint8 96×96×3 crops; 7 DMSO cells on Week9_39282 |
| Compound encoding | IMPA emb_fp.csv, 35 compounds × 1,024 features |
| ToolUniverse dispatch | Real 2-step CPU inference succeeded |
| Python SDK | Real 50-step CPU inference succeeded for floxuridine |
| Output | One 96×96×3 uint8 PNG, JSONL index, run manifest and log |
| Sampling | Official EMA, EDM schedule, heun2, guidance 0.2, noise 1.0, seed 42 |
| Sampling time | 9.00 seconds at 2 steps; 222.59 seconds at 50 steps, excluding model loading |
| Solver calls | 4 and 100 respectively; schedule parameter is not actual function-call count |
| Unit/SDK checks | 16 passed, including existing loader SDK checks |
| Static checks | Ruff check and format check passed |

Runtime: PyTorch 2.8.0, NumPy 2.2.6, flow-matching 1.0.10,
torchdiffeq 0.2.5; CPU with 4 PyTorch threads. CUDA is unavailable on this host
and has not been tested. Checkpoint epoch metadata is 99.

Checkpoint SHA256:
`a3bf1e8471a6c01565a866faa1f5cf59899165d37dd7d46aeb3ea8552ec59445`.

Local artifacts, intentionally outside version control:

- `perturb_tools/models/cellflux_assets/`: source, checkpoint, embeddings, provenance.
- `perturb_tools/test_data/real/cellflux_bbbc021/`: cells, control manifest, original annotations.
- `perturb_tools/test_runs/cellflux/sdk_default50_result.json`: SDK result.
- `perturb_tools/test_runs/cellflux/validation_summary.json`: validation summary.
- `perturb_tools/test_runs/cellflux/comparison.png`: control/generated/unpaired-real visualization.

These checks establish model loading, inference and callable-tool integration.
They do not establish MoA accuracy, FID/KID, held-out performance or equivalence
to paper scores. The downloaded subset has no official train/test assignment.
One observed treated cell is included in the visualization for context only.

## Follow-up verification — 2026-11-18

After rebasing onto origin/main `bf20982`, real SDK inference with the released
checkpoint and original BBBC021 controls passed at 50 steps (one image,
333.68 seconds including preparation/loading). A second real 2-step run passed
with source-file checks before/after inference and byte-bound selected-pixel
hash validation enabled. A real one-second worker timeout returned TIMEOUT.
No model or data download occurred during these calls.

The shared runtime now rejects non-finite arguments and incomplete worker
outputs, cleans process groups on caller interruption, and validates the image
index/files before publication. The configuration disables framework result
caching even when a caller requests it. These are CPU integration checks;
GPU and benchmark validation remain outside the verified scope.

Local evidence: `perturb_tools/test_runs/tool_batch_quality/18f343fec00d44b7b5a8353b7ead2e68/summary.json`
and `perturb_tools/test_runs/tool_batch_source_integrity/`.

## Canonical name migration verification

The primary CellFlux name is now `CellFlux_predict_morphology`; old wire names
and SDK imports remain compatible. Canonical and legacy SDK calls both ran
the actual pretrained network on the same original BBBC021 controls. With the
same seed, generated PNG SHA256 values matched exactly (CellFlux used two
sampling steps for this migration check). 39 focused regression tests passed.

Local evidence: `perturb_tools/test_runs/morphology_name_validation/9b8f9f84ea4649f9a3466cad1e683ea8/summary.json`.
