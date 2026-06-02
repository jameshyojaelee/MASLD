#!/bin/bash
#SBATCH --job-name=liver_pipeline_GSE240729
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/logs/smk_%j.out
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729

mkdir -p logs/cluster

snakemake \
    -s ../../shared/Snakefile \
    --configfile workflow/config.yaml \
    --cluster-config workflow/cluster_config.yaml \
    --jobs 15 \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second 10 \
    --cluster "sbatch \
                --partition={cluster.partition} \
                --cpus-per-task={cluster.cpusPerTask} \
                --mem={cluster.mem} \
                --time={cluster.time} \
                --job-name=liver_{cluster.jobName} \
                --output={cluster.logOut} \
                --error={cluster.logErr} \
                --parsable" \
    --latency-wait 60 \
    --keep-going
