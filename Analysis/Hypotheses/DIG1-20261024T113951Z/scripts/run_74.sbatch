#!/bin/bash
#SBATCH --job-name=DIG1cat74
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --prefer=v6
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum $RUN/scripts/74_catalog_enzyme_locus_sidecar.py > $RUN/script_sha256_job74.txt
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I $RUN/scripts/74_catalog_enzyme_locus_sidecar.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/evidence_passports/candidates/gene-catalog-r8-anchor-disease-sidecar-20261025T131731Z
