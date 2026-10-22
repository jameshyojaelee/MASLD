#!/bin/bash
#SBATCH --job-name=9843
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/job_9843_%j.log
# Calibration: the same label/geometry checks as on the main panels, run on each calib copy; compared to my recomputed values.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
C=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/calib/fig4
for pair in "a:g4a_check.py:fig4a_atlas_program_schematic" "b:g4b_geometry.py:fig4b_scrna_umap_embeddable" "c:g4c_geometry.py:fig4c_hotspot_stage_heatmap" \
            "d:g4d_geometry.py:fig4d_same_atlas_lineage_specificity" "e:g4e_geometry.py:fig4e_bulk_tissue_state_transport" \
            "f:g4f_geometry.py:fig4f_continuum_program_trajectories" "g:g4g_geometry.py:fig4g_molecular_systems_continuum_effect_map" \
            "h:g4h_geometry.py:fig4h_molecular_systems_stage_vs_continuum"; do
  IFS=: read p s f <<< "$pair"
  echo "================ calib 4${p^^} ($f.pdf)"
  $PY $s $C/$f.pdf calib_4$p 2>&1 | grep -v -i "warning" || echo "CHECK FAILED for 4$p"
done
echo "================ fonts"
$PY fonts_check.py $C/*.pdf
for f in $C/*.pdf; do echo "== pdftotext $(basename $f)"; pdftotext -layout $f - | sed '/^\s*$/d' | head -40; done
