#!/usr/bin/env bash
# Submit from the development machine. Pin IMAGE to a verified digest for releases.
set -euo pipefail
P=${1:?Pass the absolute BEAM_ToolOcean project directory}
OUT=${2:?Pass a new evidence directory on the shared filesystem}
IMAGE=${IMAGE:-registry.h.pjlab.org.cn/ailab-deepdivegzy-deepdivegzy_gpu/slime:v0}
env -u https_proxy -u HTTPS_PROXY -u HTTP_PROXY -u http_proxy -u ALL_PROXY -u all_proxy \
 rlaunch --namespace=ailab-llmagent --charged-group=llmagent_gpu_pool \
 --private-machine=no --gpu=1 --cpu=4 --memory=16000 --forward-agent=false \
 --image="$IMAGE" --mount=gpfs://gpfs2/gpfs2-agent-math:/mnt/shared-storage-gpfs2/gpfs2-agent-math \
 --max-wait-duration=8m -- bash "$P/ToolUniverse/examples/aging/validate_worker.sh" "$P" "$OUT"
