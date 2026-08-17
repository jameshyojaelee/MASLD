#!/usr/bin/env bash
#SBATCH --job-name=fig2h_xanc_audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fig2h_xanc_audit_%j.out
#SBATCH --error=logs/fig2h_xanc_audit_%j.err

set -euo pipefail
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
stamp="$(date +%Y%m%d_%H%M%S)"
export FIG2H_AUDIT_DIR="${project_root}/figures/candidates/fig2h_cross_ancestry_${stamp}_${SLURM_JOB_ID}"
mkdir -p "${FIG2H_AUDIT_DIR}"
cd "${project_root}"
python3 scripts/figures/audit_fig2h_cross_ancestry_candidates.py
sha256sum GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv \
  > "${FIG2H_AUDIT_DIR}/input_sha256.txt"
printf '%s\n' "${FIG2H_AUDIT_DIR}"
