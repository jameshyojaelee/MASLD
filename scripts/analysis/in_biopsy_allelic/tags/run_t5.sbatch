#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
# Model A analysis steps after the PRESPEC is hashed. One step per submission:
#   t5a  participant table, main effects, G2, guard, stage MDE
#   t5b  stage + lineage tests, G3            (after t5a)
#   t6   predictive L0/L1 cross-fitting      (after t5a)
#   a9   coding variants                     (after t5a)
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_t5.sbatch <step> <exec_dir> <restricted_dir> <identity_dir> <seal_dir> <t2_json> <composition_tsv> <coding_ad_dir>
set -euo pipefail
step=$1; ex=$2; rs=$3; idt=$4; seal=$5; t2=$6; comp=$7; cad=$8
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
here=$repo/scripts/analysis/in_biopsy_allelic
py="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
meta=$repo/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv
c130=$repo/RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv
p193=$repo/Analysis/MASLD_Model_Benchmark/executions/excluded_cohort_axis_placement-20260829T232000Z/gse193066_donor_keyed_placement.tsv
sha256sum --check "$here/PRESPEC.sha256"
$py -m pip freeze > "$ex/pip_freeze_$step.txt"
case $step in
  t5a)
    $py "$here/tags/t5a_allelic_table.py" --tags "$ex/t3/tags_final.tsv" --ad "$rs/tag_allele_depths.tsv.gz" --t2 "$t2" \
      --identity-dir "$idt" --seal-dir "$seal" --metadata "$meta" --gse130970-controls "$c130" --gse193066-placement "$p193" \
      --composition "$comp" \
      --gtex-egenes "$repo/data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz" \
      --gtf /gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf \
      --rediportal-hits "$ex/t3/tags_at_rediportal_sites.tsv" \
      --restricted-out "$rs/t5a" --out "$ex/t5a" ;;
  t5b)
    $py "$here/tags/t5b_stage_lineage.py" --participant-gene "$rs/t5a/participant_gene_allelic.tsv.gz" \
      --t5a-dir "$ex/t5a" --identity-dir "$idt" --out "$ex/t5b" ;;
  t6)
    $py "$here/tags/t6_predict.py" --participant-gene "$rs/t5a/participant_gene_allelic.tsv.gz" --t5a-dir "$ex/t5a" \
      --identity-dir "$idt" \
      --pseudobulk "$repo/Analysis/MASLD_Model_Benchmark/executions/deconv_estimand_harmonization-20260828T141918Z/intermediate/reference_pseudobulk_gene_by_celltype.tsv" \
      --out "$ex/t6" ;;
  a9)
    $py "$here/a9_coding/a9_analysis.py" --ad-dir "$cad" --variants "$here/a9_coding/coding_variants.tsv" --t2 "$t2" \
      --identity-dir "$idt" --seal-dir "$seal" --metadata "$meta" --gse130970-controls "$c130" --gse193066-placement "$p193" \
      --programs "$repo/Analysis/MASLD_Model_Benchmark/executions/aspect-latent-b2-inputs-20260924T003835Z/b_inputs_corr.tsv.gz" \
      --out "$ex/a9" ;;
  *) echo "unknown step $step"; exit 2 ;;
esac
