#!/bin/bash
#SBATCH --job-name=reads_indep
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=Analysis/Hypotheses/READS-20261009T141017Z/tests_v2/indep_check/logs/indep-%j.out
#SBATCH --error=Analysis/Hypotheses/READS-20261009T141017Z/tests_v2/indep_check/logs/indep-%j.err
# Independent pileup recount of the named variants for one library. Usage: sbatch run_indep.sbatch BAM
set -euo pipefail
BAM="$1"
IC=Analysis/Hypotheses/READS-20261009T141017Z/tests_v2/indep_check
VAR=Analysis/Hypotheses/READS-20261009T141017Z/resources/build_v2_nsamples50/variants_hg38.tsv
LIB=$(basename "$BAM" | cut -d. -f1)
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I "$IC/indep_variant_check.py" "$BAM" "$IC/idx" "$VAR" "$IC/$LIB.indep_variants.tsv"
