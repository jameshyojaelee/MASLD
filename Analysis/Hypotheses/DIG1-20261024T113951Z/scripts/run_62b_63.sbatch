#!/bin/bash
#SBATCH --job-name=DIG1gt63
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum "$RUN/scripts/62b_genotype_table.py" "$RUN/scripts/63_genotype_state.R" > "$RUN/script_sha256_job62b_63.txt"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I "$RUN/scripts/62b_genotype_table.py" "$RUN/out_rna_genotypes" "$RUN/out_genotype_table"
echo "=== job 63 ==="
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$RUN/scripts/63_genotype_state.R" "$RUN/out_genotype_table" "$RUN/out_genotype_state"
