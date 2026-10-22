#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/logs/%j.out
# AN3 task 1: GSE126848 filterByExpr refit and Fig 3C union (seed 20260923 set in the R script).
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/an3_gse126848_filtered_union.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/* > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/output_sha256.tsv
