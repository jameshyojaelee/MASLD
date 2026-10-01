#!/bin/bash
#SBATCH --job-name=portal_inputs
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/logs/gene_catalog_inputs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/logs/gene_catalog_inputs_%j.err

# Gene Catalog inputs on the corrected bulk fit and corrected COLOC release.
# Bulk and COLOC reach the Catalog only through GEN_RAW (frozen_evidence_classes.tsv)
# and its identity-adjudicated copy GEN_CLASSES. This job rebuilds those two and
# writes an UNSIGNED passport_input_selection.tsv in a new passport root:
#   1. 11_write_corrected_input_contract.py  v1-shaped contract + expected counts from the inputs
#   2. 01_freeze_and_rederive.py             primary frozen classes: canonical state rule
#                                            (padj < 0.05 and |logFC| > 0.5), SuSiE PP4 > 0.5,
#                                            genetic_untestable_shared_posterior as its own class
#   3. the same with the TREAT rule          named sensitivity, sensitivity_treat/
#   4. 12_compare_state_rules.py             class differences between the two rules
#   5. adjudicate_gen_ensembl_identity.py    GENCODE v49 identity adjudication (primary)
#   6. prepare_passport_selection.py         unsigned selection; GEN_CLOSURE and MYOJIN_CLASS
#                                            limitations are marked stale
# There is no unsigned build mode: every passport build requires the signed
# attestation, so the chain stops at the unsigned selection.
# Required (absolute paths):
#   GEN_BULK_DEG                corrected F_five deg_results.csv
#   GEN_COLOC_TIER12            susie_coloc_r2_<ts>/tier12/gene_level_coloc_tier12.csv
#   GEN_COLOC_GENE_LEVEL        susie_coloc_r2_<ts>/gene_level_coloc.csv
#   GEN_COLOC_ALL_GWAS          susie_coloc_r2_<ts>/susie_coloc_all_gwas.csv
#   GEN_CANDIDATE_ROOT          new Analysis/Multimodal_Program_Projection/candidates/<name>
#   PASSPORT_CANDIDATE_ROOT     new RNA-seq/results/evidence_passports/candidates/<name>
#   PASSPORT_ANALYSIS_RELEASE_ID
# Signing (prepare_passport_attestations.py selection) is the user's decision at
# adoption; run_gene_catalog.sbatch then runs with the same PASSPORT_CANDIDATE_ROOT
# and PASSPORT_ANALYSIS_RELEASE_ID exported.
set -euo pipefail
umask 027

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
for name in GEN_BULK_DEG GEN_COLOC_TIER12 GEN_COLOC_GENE_LEVEL GEN_COLOC_ALL_GWAS \
            GEN_CANDIDATE_ROOT PASSPORT_CANDIDATE_ROOT PASSPORT_ANALYSIS_RELEASE_ID; do
  : "${!name:?${name} is required}"
done
for name in GEN_BULK_DEG GEN_COLOC_TIER12 GEN_COLOC_GENE_LEVEL GEN_COLOC_ALL_GWAS \
            GEN_CANDIDATE_ROOT PASSPORT_CANDIDATE_ROOT; do
  [[ "${!name}" == /* ]] || { echo "${name} must be absolute" >&2; exit 2; }
done
for input in "$GEN_BULK_DEG" "$GEN_COLOC_TIER12" "$GEN_COLOC_GENE_LEVEL" "$GEN_COLOC_ALL_GWAS"; do
  test -s "$input"
done
test ! -e "$GEN_CANDIDATE_ROOT"
test ! -e "$PASSPORT_CANDIDATE_ROOT"

gen_scripts="$project_root/Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2"
gen_root="$GEN_CANDIDATE_ROOT/genetics_context"
treat_root="$GEN_CANDIDATE_ROOT/sensitivity_treat/genetics_context"
export MASLD_PROJECT_ROOT="$project_root"
export PASSPORT_CANDIDATE_ROOT PASSPORT_ANALYSIS_RELEASE_ID
export PASSPORT_GEN_RAW_ROOT="$gen_root"
export PASSPORT_GEN_RAW_RELEASE="$(basename "$GEN_CANDIDATE_ROOT"):GEN"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

# The corrected COLOC release must carry the method-specific ABF driver.
tier12_header="$(head -n 1 "$GEN_COLOC_TIER12")"
[[ ",${tier12_header}," == *",abf_driving_gwas,"* ]] || { echo "tier-1/2 table lacks abf_driving_gwas" >&2; exit 2; }

cd "$gen_scripts"
for rule in canonical treat; do
  inputs="$GEN_CANDIDATE_ROOT/gen_inputs/$rule"
  root="$gen_root"; [[ "$rule" == treat ]] && root="$treat_root"
  micromamba run -n spatial python -B 11_write_corrected_input_contract.py \
    --bulk "$GEN_BULK_DEG" --tier12 "$GEN_COLOC_TIER12" --full "$GEN_COLOC_GENE_LEVEL" \
    --all-gwas "$GEN_COLOC_ALL_GWAS" --state-rule "$rule" --untestable-state \
    --contract-out "$inputs/input_contract.tsv" --expected-out "$inputs/expected_counts.json"
  GEN_CONTEXT_ROOT_REL="$root" micromamba run -n spatial python -B 01_freeze_and_rederive.py \
    --candidate-root "$root" --contract "$inputs/input_contract.tsv" \
    --expected "$inputs/expected_counts.json" --state-rule "$rule" --untestable-state
done
micromamba run -n spatial python -B 12_compare_state_rules.py \
  --canonical "$gen_root/frozen_evidence_classes.tsv" \
  --treat "$treat_root/frozen_evidence_classes.tsv" \
  --out-dir "$GEN_CANDIDATE_ROOT/state_rule_comparison"

cd "$project_root/scripts/portal"
micromamba run -n spatial python -B adjudicate_gen_ensembl_identity.py \
  --source "$gen_root/frozen_evidence_classes.tsv" \
  --identity "$project_root/data/gencode_v49_gene_metadata.tsv.gz" \
  --output-root "$PASSPORT_CANDIDATE_ROOT/preflight/gen_identity_adjudication_v1" \
  --release-id "$PASSPORT_ANALYSIS_RELEASE_ID"
micromamba run -n spatial python -B prepare_passport_selection.py \
  --output "$PASSPORT_CANDIDATE_ROOT/passport_input_selection.tsv"
printf 'item\tstatus\tdetail\n%s\n%s\n%s\n' \
  $'GEN_CLOSURE\tstale_text\tallowed_downstream_use all_joint_34_of_447 counts the v1 frozen classes' \
  $'MYOJIN_CLASS\tstale_not_rerun\tclass test used v1 classes (228/463/15); Myojin is not rerun' \
  $'SIGNATURE\tuser_decision\tunsigned selection; no unsigned build mode exists' \
  > "$GEN_CANDIDATE_ROOT/gen_inputs/open_items.tsv"
micromamba run -n spatial python -m pip freeze > "$GEN_CANDIDATE_ROOT/gen_inputs/pip_freeze.txt"
echo "UNSIGNED_SELECTION $PASSPORT_CANDIDATE_ROOT/passport_input_selection.tsv"
