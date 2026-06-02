#!/bin/bash
# =============================================================================
# 10-cohort PARALLEL integration pipeline
#
# Submits SLURM jobs with dependency chains to maximize parallelism.
# Run from login node (lightweight — just submits jobs).
#
# Dependency graph:
#   Phase 1 (00→01→02→03): sequential metadata + QC + counts
#   Phase 2a: 04 (varpart), 05 (dream), 26 (sex), 09 (UMAP),
#             13 (NAFL/NASH), 14 (fibrosis), 15 (NAS) — all parallel
#   Phase 2b: 25 (deconv attribution) — after 05
#   Phase 3:  07 (consensus), 08 (pathway), 10 (volcano) — after 05
#   Phase 4:  12 (library intersection) — after 07
#             16 (disease signatures) — after 07+13+14+15
#             17 (annotate) — after 16
#
# Usage:
#   bash run_10cohort_parallel.sh
#   # Or with a prior dependency:
#   WAIT_FOR=12345678 bash run_10cohort_parallel.sh
# =============================================================================

set -euo pipefail

S="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
L="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs"
R="${S}/run_rscript.sh"
mkdir -p "${L}"

dep_flag() {
    # Build --dependency flag from one or more job IDs
    local ids="$1"
    echo "--dependency=afterok:${ids}"
}

submit() {
    local name="$1" cpus="$2" mem="$3" time="$4" partition="$5" dep="$6"
    shift 6
    local scripts=("$@")

    local dep_arg=""
    [ -n "$dep" ] && dep_arg=$(dep_flag "$dep")

    sbatch --parsable \
        --job-name="int_${name}" \
        --partition="${partition}" \
        --cpus-per-task="${cpus}" \
        --mem="${mem}" \
        --time="${time}" \
        --output="${L}/${name}_%j.out" \
        --error="${L}/${name}_%j.err" \
        ${dep_arg} \
        "${R}" "${scripts[@]}"
}

echo "=== 10-Cohort Parallel Integration ==="
echo "Submitting jobs: $(date)"
echo ""

# Optional: wait for a prior job (e.g., featureCounts merge)
PRIOR="${WAIT_FOR:-}"

# ── Phase 1: Sequential metadata → QC → per-study DE → merge counts ──
P1=$(submit "phase1_meta" 16 64G 3:00:00 cpu "$PRIOR" \
    "${S}/00_harmonize_metadata.R" \
    "${S}/01_sample_qc.R" \
    "${S}/02_per_study_de.R" \
    "${S}/03_integrate_counts.R")
echo "Phase 1 (00-03): $P1"

# ── Phase 2a: Parallel after Phase 1 (all need merged_dge.rds only) ──
J04=$(submit "varpart" 32 128G 4:00:00 cpu "$P1" \
    "${S}/04_variance_partition.R")
echo "  04 variance partition: $J04"

J05=$(submit "dream" 32 200G 6:00:00 cpu "$P1" \
    "${S}/05_dream_mega_analysis.R")
echo "  05 dream mega-analysis: $J05"

J26=$(submit "sex_strat" 32 200G 8:00:00 cpu "$P1" \
    "${S}/26_sex_stratified_analysis.R")
echo "  26 sex stratification: $J26"

J09=$(submit "umap" 4 32G 1:00:00 cpu "$P1" \
    "${S}/09_batch_correction_umap.R")
echo "  09 UMAP: $J09"

J13=$(submit "nafl_nash" 8 64G 2:00:00 cpu "$P1" \
    "${S}/13_nafl_vs_nash_de.R")
echo "  13 NAFL vs NASH: $J13"

J14=$(submit "fibrosis" 8 64G 2:00:00 cpu "$P1" \
    "${S}/14_fibrosis_progression_de.R")
echo "  14 fibrosis progression: $J14"

J15=$(submit "nas_comp" 8 64G 2:00:00 cpu "$P1" \
    "${S}/15_nas_component_de.R")
echo "  15 NAS components: $J15"

# ── Phase 2b: After dream (needs dream_results.csv) ──
J25=$(submit "deconv" 32 200G 6:00:00 cpu "$J05" \
    "${S}/25_deconv_attribution.R")
echo "  25 deconv attribution: $J25 (after dream)"

# ── Phase 3: After dream ──
J06=$(submit "meta_analysis" 4 32G 1:00:00 cpu "$J05" \
    "${S}/06_meta_analysis.R")
echo "  06 meta-analysis (reference): $J06 (after dream)"

J07=$(submit "consensus" 4 32G 1:00:00 cpu "$J05" \
    "${S}/07_consensus_degs.R")
echo "  07 consensus DEGs: $J07"

J08=$(submit "pathway" 4 32G 1:00:00 cpu "$J05" \
    "${S}/08_pathway_analysis.R")
echo "  08 pathway analysis: $J08"

J10=$(submit "volcano" 4 32G 1:00:00 cpu "$J05" \
    "${S}/10_volcano_plots.R")
echo "  10 volcano plots: $J10"

# ── Phase 4: After consensus + disease subtypes ──
J12=$(submit "library_int" 4 32G 1:00:00 cpu "$J07" \
    "${S}/12_library_intersection.R")
echo "  12 library intersection: $J12"

J16=$(submit "disease_sig" 4 32G 2:00:00 cpu "${J07}:${J13}:${J14}:${J15}" \
    "${S}/16_disease_signatures_consensus.R")
echo "  16 disease signatures: $J16 (after 07+13+14+15)"

J17=$(submit "annotate" 4 32G 1:00:00 cpu "$J16" \
    "${S}/17_annotate_gene_symbols.R")
echo "  17 annotate symbols: $J17"

echo ""
echo "=== All jobs submitted ==="
echo ""
echo "Dependency chain:"
echo "  Phase 1: $P1 (00→01→02→03)"
echo "  Phase 2: $J04 $J05 $J26 $J09 $J13 $J14 $J15 (parallel)"
echo "  Phase 2b: $J25 (after $J05)"
echo "  Phase 3: $J06 $J07 $J08 $J10 (after $J05)"
echo "  Phase 4: $J12 $J16 $J17 (after consensus+subtypes)"
echo ""
echo "Monitor: squeue -u \$USER --name='int_%'"
echo "Estimated wall time: ~6-8h (vs ~12h sequential)"
