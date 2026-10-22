#!/bin/bash
#SBATCH --job-name=9844
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/job_9844_%j.log
# Calibration follow-up: 4D axis-title check (rotated text read from chars) and 4G coordinate-tie resolution on the calib copies.
export OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
C=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/calib/fig4
$PY g4d_geometry.py $C/fig4d_same_atlas_lineage_specificity.pdf calib_4d | tail -4
$PY g4g_resolve.py calib_4g
$PY g4h_geometry.py $C/fig4h_molecular_systems_stage_vs_continuum.pdf calib_4h | grep -E "points checked|printed rho|x title"
