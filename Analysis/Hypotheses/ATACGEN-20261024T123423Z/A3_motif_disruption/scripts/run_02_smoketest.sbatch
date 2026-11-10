#!/bin/bash
#SBATCH --job-name=A3_annot_smoke
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption/logs/%x_%j.err
set -euo pipefail
A3=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption
E=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.mamba/atac-context-v3-chromvar-r46
A3_OUT=$A3/work/smoketest $E/bin/Rscript $A3/scripts/02_annotate_hg38.R
