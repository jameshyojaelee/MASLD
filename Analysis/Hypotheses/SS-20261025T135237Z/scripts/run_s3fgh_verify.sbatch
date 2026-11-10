#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS S3f/S3g/S3h verification on the 774-unit release mock (no SteatoSITE data): S3f with --dataset-term must give
# DIG1 job 63 PNPLA3 +0.176 (SE 0.041, n 757); S3g with --dataset-term on the job 41a stage set must give TNFSF10
# fibrosis rho -0.335 (n 629); S3h on simulated counts must recover planted effects.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
S=$SS/scripts; V=$SS/out_verify_s3fgh_$(date -u +%Y%m%dT%H%M%SZ); mkdir "$V"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript; PY="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I"
T=$SS/out_s0d_targets
sha256sum $S/s3_mock_release.R $S/s3_allele_counts.py $S/s3f_genotype_table.py $S/s3f_genotype_models.R $S/s3g_expression.R \
  $S/s3h_stage_ase.R $S/s3h_synthetic_test.R $S/run_s3fgh_verify.sbatch > "$V/script_sha256.txt"
$RS $S/s3_mock_release.R "$V/mock"
$PY $S/s3_allele_counts.py --input "$V/mock" --targets "$T" --out "$V/counts"
$PY $S/s3f_genotype_table.py --input "$V/mock" --counts "$V/counts" --targets "$T" --out "$V/geno" --dataset-term
echo "=== S3f (dataset term) ==="
$RS $S/s3f_genotype_models.R --input "$V/mock" --stage "$V/mock/stage_all.tsv" --geno "$V/geno" --counts "$V/counts" --targets "$T" --out "$V/s3f" --dataset-term
echo "=== S3g (dataset term) ==="
$RS $S/s3g_expression.R --input "$V/mock" --stage "$V/mock/stage_disease4.tsv" --geno "$V/geno" --targets "$T" --out "$V/s3g" --dataset-term
echo "=== S3h synthetic ==="
$RS $S/s3h_synthetic_test.R "$V/s3h_synthetic"
echo "verify dir: $V"
