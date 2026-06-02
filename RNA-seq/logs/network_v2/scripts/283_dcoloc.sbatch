#!/bin/bash
#SBATCH --job-name=netv2_283_dcoloc
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "283 | Start: $(date)"
Rscript 283_coloc_locus_pairs.R
echo "283 | Exit: $? | End: $(date)"
