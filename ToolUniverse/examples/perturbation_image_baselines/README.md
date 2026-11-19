# IMPA image generation

## Assets and environment

Use a separate Python environment with a matching CPU or CUDA PyTorch build.

```bash
python -m venv --system-site-packages perturb_tools/models/baselines-env
perturb_tools/models/baselines-env/bin/python -m pip install \
  -r ToolUniverse/examples/perturbation_image_baselines/requirements-inference.txt

perturb_tools/models/baselines-env/bin/python \
  ToolUniverse/examples/perturbation_image_baselines/prepare_assets.py \
  --asset-dir perturb_tools/models/baseline_assets
```

A clean environment is preferable for deployment. The development setup above
reuses installed PyTorch. The preparation script verifies the released IMPA weight
hashes and recovers the training vocabulary from the frozen embeddings.
Model calls use local assets.

## IMPA

```bash
PYTHONPATH=ToolUniverse/src python \
  ToolUniverse/examples/perturbation_image_baselines/run_prediction.py \
  --asset-dir perturb_tools/models/baseline_assets \
  --control-manifest perturb_tools/test_data/real/cellflux_bbbc021/controls.csv \
  --python "$PWD/perturb_tools/models/baselines-env/bin/python" \
  --output-dir perturb_tools/test_runs/image_baselines \
  --compound floxuridine --samples 2 --device cpu
```

The implementation loads the official generator and mapping network strictly.
The released weights use a 1024-dimensional fingerprint plus **16** random style
dimensions, determined from the checkpoint's 64×1040 mapping matrix. Current
repository example YAMLs use different parameters, so they are not substituted
for the actual checkpoint architecture. Training-only networks are not run.

The matching frozen embedding table contains 26 training compounds. Their names
are recovered by exact equality against the named Morgan fingerprint CSV,
preserving stored row order. Unknown/ambiguous rows fail rather than receiving
invented names. By default only these training compounds can be requested.
`allow_unseen_compounds=True` explicitly enables fingerprint-based OOD
experiments; it does not establish predictive validity for those compounds.

Native 96×96 cells use upstream dequantization and [-1,1] normalization. The
generator architecture supports that native input shape; no resize is applied.

## Python SDK and output

```python
from tooluniverse.tools import IMPA_predict_morphology
```

Provide `control_manifest`, `compounds`, `model_dir`, `output_dir`,
`checkpoint_path`, `embeddings_path`, and `embedding_checkpoint`.
The control CSV uses `sample_id,image_path,batch,condition`. Use preprocessed
96×96×3 uint8 DMSO cells in the original channel order. Relative image paths
resolve beside the CSV.

Each successful run writes PNGs and a `generated_images.jsonl` index, along with
its manifest and worker log. Pass the returned `image_index` to
`ImagePerturbation_evaluate_moa` or `ImagePerturbation_calculate_fid_kid`.
Use `bbbc021_smoke` with the unassigned integration subset.

The old SDK import and wire names `image_perturbation_predict_morphology_impa`
and `image_perturbation.predict_morphology_impa` remain compatible.
