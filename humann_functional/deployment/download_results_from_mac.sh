#!/usr/bin/env bash
set -euo pipefail
HOST="zyjiang@login.hoffman2.oarc.ucla.edu"
REMOTE="/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline/humann_functional"
LOCAL="${1:-./HUMAnN_cluster_runs_download}"
mkdir -p "$LOCAL"
scp -r "$HOST:$REMOTE/cluster_runs/" "$LOCAL/"
echo "Downloaded cluster runs into $LOCAL"
