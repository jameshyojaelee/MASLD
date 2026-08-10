#!/bin/bash
#SBATCH --job-name=plink2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/panel_af/panelaf_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/panel_af/panelaf_%j.err

# All five populations in ONE task, sequentially.  Not an array: every task
# would otherwise race on the same panel_af_manifest.tsv / panel_af_contract.json
# and the last writer would win with only its own rows.  The whole build is a
# couple of hours, so the parallelism is not worth the race.
#
# The output root is shared data (data/ld_ref/panel_af), not run-scoped, so the
# fixed log path here is safe -- unlike the run-scoped finemap/*.sbatch scripts.
set -euo pipefail

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SRC="${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap"

# The rnaseq env's binutils activation hook references unbound variables, which
# trips `set -u`.  Relax it just for activation.
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "[$(date -u +%FT%TZ)] building panel AF sidecars for all populations"
Rscript "${SRC}/00a_build_panel_af.R"
echo "[$(date -u +%FT%TZ)] done"
