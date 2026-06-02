#!/bin/bash
#SBATCH --job-name=master_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --partition=cpu
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# Activate Environment
source ~/.bashrc
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl

# Snakemake Execution
# We use --cluster to submit rule jobs to Slurm
# We limit concurrent jobs to 50 to avoid queue flooding

snakemake --snakefile workflow/Snakefile \
          --configfile workflow/config.yaml \
          --jobs 50 \
          --latency-wait 60 \
          --default-resources "partition=cpu" "mem_mb=4000" \
          --cluster "sbatch --partition={resources.partition} --mem={resources.mem_mb}M --cpus-per-task={threads} --output=logs/slurm/%x_%j.out --error=logs/slurm/%x_%j.err --parsable" \
          --keep-going \
          --rerun-incomplete \
          --printshellcmds

