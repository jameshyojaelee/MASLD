#!/bin/bash
# =============================================================================
# run_ncrna_pipeline.sh
# SLURM Orchestrator: ncRNA Analysis Pipeline
#
# Dependency graph:
#   M1 (landscape, 2h) ──┬──> M2 (ceRNA, 6h) ──────────┐
#                         ├──> M4 (conservation, 4h) ────┤
#                         └──> M5 (epigenomic, 2h) ──────┼─> M6 (integration, 1h) ─> M7 (figure, 2h)
#                                                        │
#   M3 (ELATUS/scRNA, 4h, GPU) ─────────────────────────┘
#       [independent]
#
# Critical path: M1 (2h) -> M2 (6h) -> M6 (1h) -> M7 (2h) = ~11h
# Total with parallelism: ~10h
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOGDIR="${BASE}/RNA-seq/logs/ncrna_pipeline"
SCRIPTDIR="${BASE}/RNA-seq/logs/ncrna_pipeline/scripts"
mkdir -p "${LOGDIR}" "${SCRIPTDIR}"

echo "=== ncRNA Analysis Pipeline Orchestrator ==="
echo "Start: $(date)"
echo "Log directory: ${LOGDIR}"
echo ""

# ── Module 1: ncRNA Landscape Characterization ──
cat > "${SCRIPTDIR}/m1_landscape.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M1_landscape
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 1: ncRNA Landscape | Start: $(date)"
Rscript 53_ncrna_landscape.R
echo "Module 1: Exit code: $? | End: $(date)"
SBEOF

M1=$(sbatch --parsable \
  --output="${LOGDIR}/M1_landscape_%j.out" \
  --error="${LOGDIR}/M1_landscape_%j.err" \
  "${SCRIPTDIR}/m1_landscape.sbatch")
echo "M1 (landscape):     Job ${M1}"

# ── Module 3: ELATUS/scRNA Profiling (INDEPENDENT — GPU) ──
cat > "${SCRIPTDIR}/m3_elatus.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M3_elatus
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell
export LD_LIBRARY_PATH=$(python -c 'import nvidia.cusparselt.lib; print(nvidia.cusparselt.lib.__path__[0])' 2>/dev/null):${LD_LIBRARY_PATH:-}
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "Module 3: ELATUS/scRNA | Start: $(date)"
python Analysis/SingleCell/scripts/elatus_lncrna_profiling.py
echo "Module 3: Exit code: $? | End: $(date)"
SBEOF

M3=$(sbatch --parsable \
  --output="${LOGDIR}/M3_elatus_%j.out" \
  --error="${LOGDIR}/M3_elatus_%j.err" \
  "${SCRIPTDIR}/m3_elatus.sbatch")
echo "M3 (ELATUS/scRNA):  Job ${M3} (independent, GPU)"

# ── Module 2: ceRNA Network (depends on M1) ──
cat > "${SCRIPTDIR}/m2_cerna.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M2_cerna
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 2: ceRNA Network | Start: $(date)"
Rscript 54_cerna_network.R
echo "Module 2: Exit code: $? | End: $(date)"
SBEOF

M2=$(sbatch --parsable \
  --dependency=afterok:${M1} \
  --output="${LOGDIR}/M2_cerna_%j.out" \
  --error="${LOGDIR}/M2_cerna_%j.err" \
  "${SCRIPTDIR}/m2_cerna.sbatch")
echo "M2 (ceRNA):         Job ${M2} (after M1:${M1})"

# ── Module 4: Synteny Conservation (depends on M1) ──
cat > "${SCRIPTDIR}/m4_conservation.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M4_conservation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 4: Synteny Conservation | Start: $(date)"
Rscript 55_ncrna_conservation.R
echo "Module 4: Exit code: $? | End: $(date)"
SBEOF

M4=$(sbatch --parsable \
  --dependency=afterok:${M1} \
  --output="${LOGDIR}/M4_conservation_%j.out" \
  --error="${LOGDIR}/M4_conservation_%j.err" \
  "${SCRIPTDIR}/m4_conservation.sbatch")
echo "M4 (conservation):  Job ${M4} (after M1:${M1})"

# ── Module 5: Epigenomic Regulation (depends on M1) ──
cat > "${SCRIPTDIR}/m5_epigenomic.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M5_epigenomic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 5: Epigenomic Regulation | Start: $(date)"
Rscript 56_ncrna_epigenomic.R
echo "Module 5: Exit code: $? | End: $(date)"
SBEOF

M5=$(sbatch --parsable \
  --dependency=afterok:${M1} \
  --output="${LOGDIR}/M5_epigenomic_%j.out" \
  --error="${LOGDIR}/M5_epigenomic_%j.err" \
  "${SCRIPTDIR}/m5_epigenomic.sbatch")
echo "M5 (epigenomic):    Job ${M5} (after M1:${M1})"

# ── Module 6: Atlas Integration (depends on M2, M3, M4, M5) ──
cat > "${SCRIPTDIR}/m6_integration.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M6_integration
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "Module 6: Atlas Integration | Start: $(date)"
Rscript 57_ncrna_atlas_integration.R
echo "Module 6: Exit code: $? | End: $(date)"
SBEOF

M6=$(sbatch --parsable \
  --dependency=afterok:${M2}:${M3}:${M4}:${M5} \
  --output="${LOGDIR}/M6_integration_%j.out" \
  --error="${LOGDIR}/M6_integration_%j.err" \
  "${SCRIPTDIR}/m6_integration.sbatch")
echo "M6 (integration):   Job ${M6} (after M2:${M2}, M3:${M3}, M4:${M4}, M5:${M5})"

# ── Module 7: Supplementary Figure (depends on M6) ──
cat > "${SCRIPTDIR}/m7_figure.sbatch" << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=ncRNA_M7_figure
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "Module 7: Supplementary Figure | Start: $(date)"
Rscript scripts/figures/figS_ncrna.R
echo "Module 7: Exit code: $? | End: $(date)"
SBEOF

M7=$(sbatch --parsable \
  --dependency=afterok:${M6} \
  --output="${LOGDIR}/M7_figure_%j.out" \
  --error="${LOGDIR}/M7_figure_%j.err" \
  "${SCRIPTDIR}/m7_figure.sbatch")
echo "M7 (figure):        Job ${M7} (after M6:${M6})"

echo ""
echo "=== Pipeline Submitted ==="
echo "  M1 landscape:     ${M1}"
echo "  M2 ceRNA:         ${M2} (after M1)"
echo "  M3 ELATUS/scRNA:  ${M3} (independent, GPU)"
echo "  M4 conservation:  ${M4} (after M1)"
echo "  M5 epigenomic:    ${M5} (after M1)"
echo "  M6 integration:   ${M6} (after M2,M3,M4,M5)"
echo "  M7 figure:        ${M7} (after M6)"
echo ""
echo "Monitor: squeue -u \$USER --name='ncRNA_*'"
echo "Logs:    ${LOGDIR}/"
echo "End:     $(date)"
