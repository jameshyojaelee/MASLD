#!/bin/bash
#SBATCH --job-name=Rprof
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/prof_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/prof_%j.err
#
# Settles "is the COLOC pipeline slow because it is R?" with measurement rather
# than arithmetic. Profiles a genuine 06_susie_coloc.R run on chr1 -- the worst
# case, where the LD blocks are ~1 GB each and the median loaded matrix is
# 23,644 x 23,644 while SuSiE fits on a median of 3,567 variants.
#
# 20 minutes at a 20 ms sampling interval is ~60,000 stack samples over roughly
# a dozen genes: plenty to separate I/O and matrix reshaping from the statistics.
#
# Identical environment to the production rerun, so the profile describes the
# real thing and not a lookalike.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
FM="${PROJECT_ROOT}/GWAS/finemapping"
LDREF="${FM}/data/ld_ref"
export MASLD_PROJECT_ROOT

export AFR_LD_DIR="${LDREF}/1kg_afr_gram"
export AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram"
export SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Isolated output. Same guard as the identity harness: an empty suffix would
# resolve to results/susie_coloc/, the canonical directory behind Fig2.
export COLOC_OUT_SUFFIX="_profile"
case "${COLOC_OUT_SUFFIX}" in
  _profile) : ;;
  *) echo "FATAL: refusing to run with suffix '${COLOC_OUT_SUFFIX}'" >&2; exit 1 ;;
esac

export GWAS_NAME="2019_31311600_NAFLD_EUR"
export CHR_FILTER=1
export PROF_SECONDS=1200
export PROF_OUT="${FM}/logs/coloc_perf/coloc_prof_${SLURM_JOB_ID}.out"

echo "[$(date -u +%FT%TZ)] profiling ${GWAS_NAME} chr${CHR_FILTER} for ${PROF_SECONDS}s"
cd "${FM}/src"
Rscript "${FM}/src/perf/prof_wrapper.R"
echo "[$(date -u +%FT%TZ)] profile complete"
