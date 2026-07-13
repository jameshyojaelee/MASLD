#!/usr/bin/env bash
#SBATCH --job-name=qtlmerge
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_merge_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_merge_%j.err

# Phase-5 src/79 merge: assemble direction_features_{caqtl,sqtl}.tsv from the panel
# intermediates. Re-expresses every SIGNED feature into the label ref->alt frame
# (orient_mult) and reports descriptive sign-concordance (corroboration only, NEVER a
# scored convergence channel). Run AFTER chrombpnet/alphagenome/borzoi/spliceai finish.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="${ROOT}/GWAS/finemapping/src"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
set +u; micromamba activate rnaseq; set -u
python "${SRC}/79_score_panels.py" --panel merge
