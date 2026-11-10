#!/bin/bash
#SBATCH --job-name=DIG1t63
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/%j.err
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/scripts/63_genotype_state.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/gdir /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/out63 > /dev/null
ls /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/out63; head -3 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/out63/genotype_models.tsv | cut -c1-150; cat /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/out63/variants_used.txt; cat /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_job63_dryrun_20261025T001806Z/out63/hsd17b13_junction_qc.tsv
