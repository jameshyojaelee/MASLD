#!/bin/bash -l
#SBATCH --job-name=SuSiE
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/SuSiE_%A_%a.out
#SBATCH --error=logs/SuSiE_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$ROOT/GWAS/finemapping"
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping
export MASLD_PROJECT_ROOT="$ROOT" LD_PANEL=topld SKIP_CARMA=1
MANIFEST=results/seqfunc/functional_prior_sensitivity/pilot_manifest.tsv
IFS=$'\t' read -r study ld_pop locus N_tot N_cases window_mb ancestry < <(awk -F'\t' -v OFS='\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1 {print $1,$2,$3,$4,$5,$6,$7}' "$MANIFEST")
test -n "$study"
Rscript src/108_functional_prior_susie.R "$study" "$ld_pop" "$locus" "$N_tot" "$N_cases" "$window_mb" "$ancestry"
