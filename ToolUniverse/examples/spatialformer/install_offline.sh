#!/usr/bin/env bash
# Run inside slime:v0 on the GPU worker; downloads happen on the development host.
set -euo pipefail
runtime=${1:?Usage: install_offline.sh RUNTIME_DIR WHEELHOUSE}
wheelhouse=${2:?Usage: install_offline.sh RUNTIME_DIR WHEELHOUSE}
python -m venv --system-site-packages /opt/spatialformer-venv
venv_python=/opt/spatialformer-venv/bin/python
mkdir -p /opt/spatialformer-vendor
cp -a "$runtime/source/spatialformer" /opt/spatialformer-vendor/
printf '%s\n' /opt/spatialformer-vendor /opt/toolocean-venv/lib/python3.12/site-packages /opt/toolocean-vendor > /opt/spatialformer-venv/lib/python3.12/site-packages/toolocean.pth
"$venv_python" -m pip install --no-index --no-deps --no-build-isolation --find-links="$wheelhouse" communities==3.0.0 transformers==4.57.3 huggingface-hub==0.36.0
"$venv_python" - <<'PY'
import spatialformer, torch, transformers
print('SpatialFormer', spatialformer.__version__, 'transformers', transformers.__version__, 'CUDA', torch.cuda.is_available())
PY
