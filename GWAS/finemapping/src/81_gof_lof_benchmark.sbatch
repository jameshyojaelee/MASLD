#!/bin/bash
#SBATCH --job-name=logofunc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/81_gof_lof_benchmark_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/81_gof_lof_benchmark_%j.log
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
# logofunc env supplies pandas1.5/lightgbm3.2.1/sklearn1.1.2 (version-locked to the pickled models);
# rnaseq/bin only added for tabix (htslib) so the fast PNPLA3 lookup path is used.
export PATH="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH"
LOGOFUNC_PY=$ROOT/.envs/logofunc/bin/python
echo "host=$(hostname) tabix=$(which tabix) py=$LOGOFUNC_PY"
cd $ROOT/GWAS/finemapping/src
$LOGOFUNC_PY 81_gof_lof_benchmark.py
echo "EXIT=$?"
