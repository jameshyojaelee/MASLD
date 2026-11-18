#!/bin/bash
#SBATCH --job-name=trail4_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/TRAIL4-20261111T025501Z/logs/trail4_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/TRAIL4-20261111T025501Z/logs/trail4_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/TRAIL4-20261111T025501Z
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
module load PLINK/1.9
OUT=$RUN/out_20261118
if [[ -e "$OUT" ]]; then echo "$OUT exists; not overwriting" >&2; exit 1; fi
mkdir -p "$OUT"
echo "host $(hostname) | $(date) | $(plink --version | head -1)"
sha256sum $RUN/PRESPEC_TRAIL4.md $RUN/scripts/01_trail4_coloc_20261118.R $RUN/scripts/run_TRAIL4_20261118.sbatch $RUN/downloads/*/*.tsv.gz > $OUT/script_and_input_sha256.txt
cd "$OUT"
Rscript $RUN/scripts/01_trail4_coloc_20261118.R "$OUT"
echo "done $(date)"
