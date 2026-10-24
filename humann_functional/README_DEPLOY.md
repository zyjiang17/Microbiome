# Hoffman2 deployment / HUMAnN pathway age analysis

**Prepared, not yet validated as a Hoffman2 rerun.** The scripts and exploratory
reference outputs are versioned here so the cluster run can be compared to the
local snapshot.

## Package layout

- `../humann_inputs/`: original two Madeline Shono TSVs (large, exact original bytes).
- `code/`: age/pathway/statistical plots, reference comparison, input QC and optional official MetaCyc export.
- `reference_results/`: CSV/JSON from preliminary local analysis; these are **not cluster reruns**.
- `reference_figures/`: original 12 PNG figures.
- `reports/`: detailed methodology/results, original report, figure captions, and PDF summary.
- `submit_humann.slurm`: job for cluster Slurm environment.
- `deployment/`: scripts for interactive upload/download from a user-controlled Mac.

## Destination

`/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline/humann_functional`

The existing lab pipeline remains intact. This is an additional folder, not a change to `workflow/run_all.py` or earlier manuscripts.

## Copy from your Mac (Slurm-preview)

1. Pull the current repository on the Slurm-preview login node.
2. Submit from the `humann_functional` directory.

```bash
git pull --ff-only
```

Or manually:

```bash
ssh zyjiang@login.hoffman2.oarc.ucla.edu 'mkdir -p /u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline/humann_functional'
scp -r HUMAnN_Hoffman2_Integration/. zyjiang@login.hoffman2.oarc.ucla.edu:/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline/humann_functional/
```

## Submit and check

```bash
ssh zyjiang@login.hoffman2.oarc.ucla.edu
cd /u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline/humann_functional
sbatch submit_humann.slurm
squeue -u zyjiang
# After job completes:
ls cluster_runs/slurm_*/reference_comparison.json
```

On the Slurm preview cluster the host is `login.hoffman2.oarc.ucla.edu`, consistent with the existing `submit_all.slurm` workflow. If the existing `.venv_preview_20261002` is missing, create another virtual environment with `requirements_humann.txt` and submit with `HUMANN_VENV=/path/to/venv sbatch submit_humann.slurm`. On production `hoffman2.idre.ucla.edu`, use Grid Engine `qsub` rather than this Slurm submission script.

**Check SSH host keys against UCLA documentation** before accepting a first-time connection. Do not paste passwords or private SSH keys into ChatGPT.

## Retrieve outputs

```bash
bash HUMAnN_Hoffman2_Integration/deployment/download_results_from_mac.sh
```

## Not-yet-resolved study QC

- 695 TSV samples vs Madeline's emailed n=625. This package **does not silently discard 70 dogs**.
- Raw sequencing depth/plate/PD status are missing; do not declare periodontal disease findings.
- Original MetaPhlAn manuscript uses separate QC1/QC2/QC3/QC4 cohorts (1125/1043/675 and other figure-specific subsets). Join to HUMAnN by IDs before comparisons.
