#!/usr/bin/env bash
# Run from your Mac in the directory containing the extracted folder.
# This script uses your *local* OpenSSH interactive authentication; it does not store passwords.
set -euo pipefail
SRC="${1:-HUMAnN_Hoffman2_Integration}"
DEST="/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline/humann_functional"
HOST="zyjiang@login.hoffman2.oarc.ucla.edu"
[[ -d "$SRC" ]] || { echo "Missing local directory: $SRC" >&2; exit 1; }
echo "You will authenticate as $HOST to create $DEST and upload the analysis."
ssh "$HOST" "mkdir -p '$DEST'"
scp -r "$SRC"/. "$HOST:$DEST/"
ssh "$HOST" "cd '$DEST' && sha256sum -c SHA256SUMS && echo 'Upload verified.' && echo 'To run: sbatch submit_humann.slurm'"
