#!/usr/bin/env bash
# Execute inside a GPU worker; no downloads or dependency installation here.
set -euo pipefail
P=${1:?Pass the absolute BEAM_ToolOcean project directory}
OUT=${2:?Pass a new evidence directory}
V=/opt/aging-venv/bin/python
mkdir -p "$OUT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 WANDB_MODE=disabled
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
unset PYTHONPATH
nvidia-smi -L
"$V" - "$P" <<'PY'
import hashlib,sys,pathlib,torch,tooluniverse
root=pathlib.Path(sys.argv[1]);installed=pathlib.Path(tooluniverse.__file__).parent
assert str(installed).startswith('/opt/'),installed
assert torch.cuda.is_available()
print('Installed package',installed,'torch',torch.__version__,'CUDA',torch.cuda.is_available())
for name in ['aging_common.py','aging_data_tool.py','aging_models_tool.py','aging_statistics_tool.py','aging_mechanism_tool.py','default_config.py','data/aging_tools.json']:
    source=root/'ToolUniverse/src/tooluniverse'/name
    assert hashlib.sha256(source.read_bytes()).digest()==hashlib.sha256((installed/name).read_bytes()).digest(),name
print('Installed Aging source hashes match reviewed workspace')
PY
"$V" "$P/ToolUniverse/examples/aging/model_smoke.py" --bundles "$P/cluster-runtime/aging/models" --upstream "$P/cluster-runtime/aging/upstream" --output "$OUT/model-smoke"
cd "$P/ToolUniverse"
"$V" -m pytest -q -o addopts= tests/unit/test_aging_data.py tests/unit/test_aging_models.py tests/unit/test_aging_statistics.py tests/unit/test_aging_mechanism.py tests/integration/test_aging_category.py --basetemp="$OUT/pytest" --junitxml="$OUT/pytest.xml"
