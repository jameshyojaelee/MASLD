#!/bin/bash
#SBATCH --job-name=A3_smoke
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption/logs/%x_%j.err
# End-to-end smoke test of steps 2-5 on the subset universe in work/smoketest (outputs stay there).
set -euo pipefail
A3=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption
E=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.mamba/atac-context-v3-chromvar-r46
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export A3_OUT=$A3/work/smoketest PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
echo "step2 $(date -u +%T)"; $E/bin/Rscript $A3/scripts/02_annotate_hg38.R
echo "step3 $(date -u +%T)"; $PY -s $A3/scripts/03_score_motifs.py
echo "step4 $(date -u +%T)"; $PY -s $A3/scripts/04_cluster_motifs.py
echo "step5 $(date -u +%T)"; $PY -s $A3/scripts/05_tests.py
echo "end $(date -u +%T)"
