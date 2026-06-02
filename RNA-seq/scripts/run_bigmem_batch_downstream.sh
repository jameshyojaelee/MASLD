#!/bin/bash
#SBATCH --job-name=downstream-batch
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/bigmem_batch_%j.out
#SBATCH --error=RNA-seq/logs/bigmem_batch_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export SLURM_CPUS_PER_TASK=8
MAX_PARALLEL=8
FAIL_COUNT=0

run_script() {
    local script="$1"
    local name="$2"
    local logdir="RNA-seq/logs/bigmem_batch"
    mkdir -p "$logdir"
    echo "[$(date)] Starting: $name"
    if Rscript "$script" > "$logdir/${name}.out" 2> "$logdir/${name}.err"; then
        echo "[$(date)] COMPLETED: $name"
    else
        echo "[$(date)] FAILED: $name"
        ((FAIL_COUNT++)) || true
    fi
}

# Wave 1: 8 parallel (MUST-RUN + high priority)
run_script "RNA-seq/20_drug_repurposing_v3.R" "LINCS" &
run_script "RNA-seq/35_broadaway_coloc.R" "broadaway-COLOC" &
run_script "RNA-seq/37_integrate_zenodo_coloc.R" "zenodo-COLOC" &
run_script "RNA-seq/38_ieqtl_analysis.R" "ieQTL" &
run_script "RNA-seq/39_pleiotropy_analysis.R" "pleiotropy" &
run_script "RNA-seq/41_ferroptosis_lipotoxicity.R" "ferroptosis" &
run_script "RNA-seq/20_drug_repurposing_cmap.R" "CMap" &
run_script "RNA-seq/30_gwas_overlay_v2.R" "gwas-overlay" &
wait
echo "[$(date)] === Wave 1 done (8 scripts) ==="

# Wave 2: 8 parallel
run_script "RNA-seq/32_network_proximity.R" "network-proximity" &
run_script "RNA-seq/34_combat_seq_sensitivity.R" "ComBat-seq" &
run_script "RNA-seq/36_pseudobulk_de.R" "pseudobulk-DE" &
run_script "RNA-seq/51_decoupler_functional_activity.R" "decoupleR" &
run_script "RNA-seq/58_diamond_kda.R" "DIAMOND" &
run_script "RNA-seq/80_celltype_intrinsic_attribution.R" "celltype-attribution" &
run_script "RNA-seq/85_secretome_plasma_chain.R" "secretome" &
run_script "RNA-seq/86_cellstate_signatures_bulk.R" "cellstate-fgsea" &
wait
echo "[$(date)] === Wave 2 done (8 scripts) ==="

# Wave 3: remaining
run_script "RNA-seq/87_celltype_drug_target_enrichment.R" "celltype-drug" &
run_script "RNA-seq/201_celltype_heritability.R" "celltype-heritability" &
run_script "RNA-seq/202_geneset_enrichment.R" "geneset-enrichment" &
run_script "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/206_published_vs_perstudy.R" "published-benchmark" &
run_script "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/206c_combined_figure.R" "benchmark-combined" &
run_script "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/206d_dream_vs_published.R" "dream-vs-published" &
run_script "RNA-seq/215_pharmacogenomic_targets.R" "pharmacogenomic" &
run_script "RNA-seq/218_ashr_lfc_threshold.R" "ashr-threshold" &
wait
echo "[$(date)] === Wave 3 done (8 scripts) ==="

# Wave 4: stragglers
run_script "RNA-seq/218_switch_gene_classification.R" "switch-genes" &
run_script "RNA-seq/53_ncrna_landscape.R" "ncRNA-landscape" &
run_script "RNA-seq/76b_govaere2026_crossmodal.R" "govaere-crossmodal" &
run_script "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/36b_pseudobulk_bulk_concordance.R" "pseudobulk-concordance" &
wait
echo "[$(date)] === Wave 4 done (4 scripts) ==="

echo "==============================================="
echo "[$(date)] ALL WAVES COMPLETE. Failures: $FAIL_COUNT"
echo "==============================================="
