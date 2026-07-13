#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/perm_gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/perm_gate_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
export MASLD_PROJECT_ROOT=${ROOT}
micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/seqfunc_v2/splicing/06_permutation_qc.R"
if test -f "${OUT}/permutations/expand_to_20.flag"; then
  jid=$(sbatch --wait --parsable --array=10-19%10 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/splicing/05_permutation_ds.sbatch")
  printf '%s\n' "${jid}" > "${OUT}/permutations/conditional_expansion_job.txt"
  micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/seqfunc_v2/splicing/06_permutation_qc.R"
  test -f "${OUT}/permutations/PASS"
fi
