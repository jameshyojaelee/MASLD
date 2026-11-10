#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3/%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
# SS jobs S3f, S3g, S3h (PRESPEC_SS_addendum_S3f/g/h.md): allele counts from the S2 VCFs, genotype table, then the
# three model scripts. Needs S2DIR, STAGE, OUT (new dir; sub-dirs counts, geno, s3f, s3g, s3h).
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
S=$SS/scripts; T=$SS/out_s0d_targets; mkdir "$OUT"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript; PY="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I"
sha256sum $S/s3_allele_counts.py $S/s3f_genotype_table.py $S/s3f_genotype_models.R $S/s3g_expression.R $S/s3h_stage_ase.R > "$OUT/script_sha256.txt"
$PY $S/s3_allele_counts.py --input "$S2DIR" --targets "$T" --out "$OUT/counts"
$PY $S/s3f_genotype_table.py --input "$S2DIR" --counts "$OUT/counts" --targets "$T" --out "$OUT/geno"
$RS $S/s3f_genotype_models.R --input "$S2DIR" --stage "$STAGE" --geno "$OUT/geno" --counts "$OUT/counts" --targets "$T" --out "$OUT/s3f"
$RS $S/s3g_expression.R --input "$S2DIR" --stage "$STAGE" --geno "$OUT/geno" --targets "$T" --out "$OUT/s3g"
$RS $S/s3h_stage_ase.R --input "$S2DIR" --stage "$STAGE" --counts "$OUT/counts" --targets "$T" --out "$OUT/s3h"
$PY -m pip freeze > "$OUT/pip_freeze.txt"
echo EXIT_OK
