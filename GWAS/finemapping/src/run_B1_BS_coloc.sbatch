#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --array=1-22%8
#SBATCH --output=GWAS/finemapping/logs/B1_BS_%x_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/B1_BS_%x_%A_%a.err

# B1 BOTH-SEXES baseline COLOC: Neale R2 both_sexes ALT/AST/GGT
# Set GWAS_NAME env via sbatch --export=GWAS_NAME=PanUKBB_BS_{ALT,AST,GGT}
# Uses cpu @ 96G (verified safe by F/M B1 jobs; MaxRSS ~60GB peak)

set -e

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate finemapping

export CHR_FILTER="${SLURM_ARRAY_TASK_ID}"
export LD_PANEL="${LD_PANEL:-polyfun}"
export N_WORKERS="${N_WORKERS:-1}"

if [ -z "${GWAS_NAME:-}" ]; then
  echo "ERROR: GWAS_NAME env var must be set (PanUKBB_BS_ALT|PanUKBB_BS_AST|PanUKBB_BS_GGT)"
  exit 2
fi

echo "=== B1 BS baseline SuSiE-COLOC: ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "  LD panel: ${LD_PANEL}"
echo "  Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Start: $(date)"

Rscript src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}"

rm -f "results/susie_coloc/${GWAS_NAME}/checkpoint_chr${CHR_FILTER}.csv" 2>/dev/null || true
echo "End: $(date)"
