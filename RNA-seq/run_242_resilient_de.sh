#!/bin/bash
#SBATCH --job-name=resilient_de
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/audit_sensitivity/healthy_control_audit/logs/242_%j.out
#SBATCH --error=RNA-seq/results/audit_sensitivity/healthy_control_audit/logs/242_%j.err

mkdir -p RNA-seq/results/audit_sensitivity/healthy_control_audit/logs

source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript RNA-seq/242_resilient_de.R
