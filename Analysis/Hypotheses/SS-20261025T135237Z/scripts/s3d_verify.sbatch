#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3d/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3d/%x_%j.err
# SS S3d verification on existing cohorts (no SteatoSITE data). 1) mocks; 2) s3d_programs.R on the four Fig 4E
# cohorts with the stored release factors and with within-cohort RLE; 3) Fig 4E four-cohort arm refit (target 16
# fibrosis calls); 4) H54 set A91 with the H54 design (targets in stage2/set_statistics.tsv, stage1/NC1_summary.tsv);
# 5) smoke run of the SteatoSITE design (3 fake batches) on GSE135251.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
V=$SS/verify/s3d-$(date -u +%Y%m%dT%H%M%SZ); [ -e $V ] && exit 1; mkdir -p $V
$R $SS/scripts/s3d_mock_inputs.R $V/mock
for c in GSE130970 GSE135251 GSE162694 GSE174478; do
  $R $SS/scripts/s3d_programs.R --input $V/mock/fig4e_$c --stage $V/mock/fig4e_$c/stage.tsv --out $V/stored/fig4e_$c --norm stored --h54 skip
  $R $SS/scripts/s3d_programs.R --input $V/mock/fig4e_$c --stage $V/mock/fig4e_$c/stage.tsv --out $V/rle/fig4e_$c --norm rle --h54 skip
done
$R $SS/scripts/s3d_verify_fig4e.R $V/mock $V/stored $V/fig4e_stored
$R $SS/scripts/s3d_verify_fig4e.R $V/mock $V/rle $V/fig4e_rle
$R $SS/scripts/s3d_programs.R --input $V/mock/h54_A91 --stage $V/mock/h54_A91/stage.tsv --out $V/h54_A91 --h54_design h54
cp -r $V/mock/fig4e_GSE135251 $V/mock/smoke_batch
$R -e 'u <- data.table::fread("'$V'/mock/smoke_batch/units.tsv"); u[, batch := paste0("b", seq_len(.N) %% 3)]; data.table::fwrite(u, "'$V'/mock/smoke_batch/units.tsv", sep = "\t", na = "NA")'
$R $SS/scripts/s3d_programs.R --input $V/mock/smoke_batch --stage $V/mock/smoke_batch/stage.tsv --out $V/smoke_batch
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I -m pip freeze > $V/pip_freeze.txt 2>/dev/null || true
echo EXIT_OK $V
