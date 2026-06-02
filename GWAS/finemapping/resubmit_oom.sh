#!/bin/bash -l
# Resubmit OOM chromosomes at 128G
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

submit_retry() {
    local gwas="$1"
    local chrs="$2"
    local JOB
    JOB=$(GWAS_NAME="$gwas" sbatch --parsable \
        --array="$chrs" --mem=128G --partition=cpu \
        --cpus-per-task=1 --time=48:00:00 \
        --job-name="coloc_retry_${gwas}" \
        --output="logs/coloc_retry_${gwas}_%A_%a.out" \
        --error="logs/coloc_retry_${gwas}_%A_%a.err" \
        src/06_susie_coloc.sh)
    echo "  ${gwas} [chr ${chrs}] -> Job ${JOB} (128G)"
}

echo "=== Resubmitting OOM chromosomes at 128G ==="
submit_retry "2019_31311600_NAFLD_EUR"           "5,6,7,8,9,10,11,12,13,17,18"
submit_retry "2020_32298765_NAFLD_EUR"           "2,5,6,7,8,10,11,12,17,18"
submit_retry "2021_34128465_PDFF_EUR"            "5,6,7,8,10,11,12,13,17,18"
submit_retry "2021_34841290_NAFLD_EUR"           "5,6,8,11,12,13,17,18"
submit_retry "2021_34957434_PDFF_EUR"            "5,6,7,8,17,18"
submit_retry "2022_36402844_PDFF_EUR"            "5,7,8,10,11,17,18"
submit_retry "2023_36280732_NAFLD_deCode_EUR"    "8"
echo "Done."
