#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3d/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3d/%x_%j.err
# SS S3d verification, part 2 (after the --h54 argument fix): H54 set A91 with the H54 design on the mocks of
# verify/s3d-20261025T141917Z, and the composition (ILR) path on GSE130970 with a placeholder ILR table (random
# normal values, seed 1; checks the code path only).
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
M=$SS/verify/s3d-20261025T141917Z/mock
V=$SS/verify/s3d2-$(date -u +%Y%m%dT%H%M%SZ); [ -e $V ] && exit 1; mkdir -p $V
$R $SS/scripts/s3d_programs.R --input $M/h54_A91 --stage $M/h54_A91/stage.tsv --out $V/h54_A91 --h54_design h54
$R -e 'set.seed(1); u <- data.table::fread("'$M'/fig4e_GSE130970/units.tsv"); x <- data.table::data.table(sample_alias = u$sample_alias, matrix(rnorm(6 * nrow(u)), ncol = 6, dimnames = list(NULL, paste0("ilr", 1:6)))); data.table::fwrite(x, "'$V'/placeholder_ilr.tsv", sep = "\t")'
$R $SS/scripts/s3d_programs.R --input $M/fig4e_GSE130970 --stage $M/fig4e_GSE130970/stage.tsv --out $V/ilr_path_GSE130970 --composition $V/placeholder_ilr.tsv --h54 skip
echo EXIT_OK $V
