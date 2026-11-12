#!/bin/bash
#SBATCH --job-name=DIR6
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIR6-20261111T024808Z/logs/dir6_%j.out
set -eo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/DIR6-20261111T024808Z
OUT=$RUN/out_${SLURM_JOB_ID}
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
source /etc/profile.d/modules.sh 2>/dev/null || true
module load BCFtools/1.19-GCC-13.2.0
bcftools --version | head -1
sha256sum $RUN/PRESPEC_DIR6.md $RUN/scripts/01_extract.py $RUN/scripts/02_classify.R $RUN/scripts/run_dir6.sbatch > $OUT/script_and_prespec_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -m pip freeze > $OUT/pip_freeze_spatial.txt 2>/dev/null || true
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I $RUN/scripts/01_extract.py $ROOT $OUT
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/02_classify.R $ROOT $OUT
rm -rf $OUT/work
echo "DIR6 DONE"
