#!/usr/bin/env bash
#SBATCH --job-name=ws4a_t2_reldig
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/Hypotheses/REL1-20261110T210214Z/logs/ws4a_t2_%j.out
#SBATCH --error=Analysis/Hypotheses/REL1-20261110T210214Z/logs/ws4a_t2_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

REL1_OUT="$PWD/Analysis/Hypotheses/REL1-20261110T210214Z/out_20261118"
DIG19_OUT="$PWD/Analysis/Hypotheses/DIG1-20261024T113951Z/out_19_20261118"
DIG35_OUT="$PWD/Analysis/Hypotheses/DIG1-20261024T113951Z/out_35_20261118"

"$R_BIN" --vanilla Analysis/Hypotheses/REL1-20261110T210214Z/scripts/rel1_20261118.R "$REL1_OUT"
"$R_BIN" --vanilla Analysis/Hypotheses/DIG1-20261024T113951Z/scripts/19_protein_rna_discordance_20261118.R "$DIG19_OUT"
"$R_BIN" --vanilla Analysis/Hypotheses/DIG1-20261024T113951Z/scripts/35_ci_composition_20261118.R "$DIG35_OUT"
