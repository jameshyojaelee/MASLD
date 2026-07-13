#!/bin/bash -l
#SBATCH --job-name=PLINK
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-55
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/leadsnp_twotier_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/leadsnp_twotier_%A_%a.err

# 02_lead_snps_twotier_all.sh
# Two-tier lead-SNP identification applied CONSISTENTLY to EVERY GWAS in the
# finemapping registry (config/gwas_registry.tsv), generalizing the MVP-only
# 02_mvp_lead_snps_twotier.sh template to all ancestry strata.
#
# Why this script exists / what it supersedes:
#   The canonical GWAS/finemapping/src/01_identify_lead_snps.sh only special-cases
#   BBJ_* (distance-based); EVERY other GWAS — including PanUKBB_AFR_* and
#   PanUKBB_CSA_* — fell into the "EUR" branch and was clumped against the EUR
#   1000G panel at a flat 5e-8. That is wrong on BOTH axes for non-EUR strata
#   (European LD + European threshold). This script fixes both: each GWAS is
#   clumped against its ANCESTRY-MATCHED 1000G panel and each index SNP is
#   labeled by a per-ancestry genome-wide-significance tier.
#
#   Two tiers (per index-SNP p-value):
#     Tier 1 "significant": index p < min(5e-8, Kanai_ancestry) = 5e-8 for
#                           EUR/EAS/AMR/SAS, 3.24e-8 for AFR (the only stratum where
#                           the panel-matched Kanai value is stricter than 5e-8).
#     Tier 2 "suggestive" : index p < 1e-5 (Lander-Kruglyak) — NON-EUR strata ONLY;
#                           EUR clumps at Tier-1 (no suggestive expansion).
#   One PLINK clump pass at the per-ancestry LOOSEST p1 (5e-8 for EUR, 1e-5 for
#   non-EUR), then each index SNP is labeled by tier from its own p-value.
#
#   This SUPERSEDES 01 for the existing GWAS: the plain output path
#   data/lead_snps/<study>_leadSNPs.tsv is SHARED with 01 and is intentionally
#   overwritten. EUR is UNCHANGED (still a 5e-8 clump — main-figure locus set
#   preserved); non-EUR now use ANCESTRY-MATCHED LD + the 1e-5 suggestive union,
#   and AFR the stricter 3.24e-8 Tier-1. Canonical 01 and the registry untouched.
#
# Array mapping:
#   SLURM_ARRAY_TASK_ID = k  ->  data row k of config/gwas_registry.tsv
#                              (= file line k+1; line 1 is the header).
#   study_name = registry col 1, sumstats_path = col 2, ancestry = col 4.
#   Array sized 1-55 (55 data rows: 23 original + 32 MVP). The sumstats input
#   is taken DIRECTLY from the registry sumstats_path column (resolved under
#   <FM>/), which already points to the right file per GWAS — usually
#   data/sumstats/<study>_reformatted_hg19.tsv, but data/sumstats/<study>_preprocessed.tsv
#   for the 6 EUR rows (2019/2020/2021 NAFLD/PDFF) — giving full 55/55 coverage
#   without special-casing. Studies whose sumstats file is absent are skipped.
#
# Outputs (per GWAS) under data/lead_snps/:
#   <study>_leadSNPs_tiered.tsv  (CHR BP locus index_p tier  — reporting deliverable)
#   <study>_leadSNPs.tsv         (CHR BP locus  — registry/fine-map compat, union @1e-6)

set -eo pipefail   # NOT set -u (breaks `module load` / lmod shell functions)
module load PLINK/1.9

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FM="${ROOT}/GWAS/finemapping"
REG="${FM}/config/gwas_registry.tsv"
cd "$FM"
OUTDIR="data/lead_snps"
TMP="data/lead_snps/tmp_twotier_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$OUTDIR" "$TMP"

# Tier-1 genome-wide significance = min(5e-8, Kanai_ancestry): field-standard 5e-8
# everywhere, tightened to the panel-matched Kanai 2016 AFR value (3.24e-8) — the
# only stratum where Kanai is STRICTER than 5e-8. Kanai's non-AFR values are more
# lenient than 5e-8 (small 1000G sub-panels under-count independent tests) so are
# floored at 5e-8. Full rationale: GWAS/MR_Data/MVP/THRESHOLD_FRAMEWORK.md.
declare -A TIER1=( [EUR]=5e-8 [AFR]=3.24e-8 [EAS]=5e-8 [AMR]=5e-8 [SAS]=5e-8 )
# Tier-2 suggestive (1e-5; Lander-Kruglyak) is applied to NON-EUR strata ONLY. EUR is
# the well-powered discovery backbone (clumps at Tier-1, no suggestive expansion) and
# COLOC tests all eGenes regardless, so sub-threshold EUR colocalizations are still
# captured. CLUMP_P1 = the loosest p clumped per ancestry (= Tier-1 for EUR, 1e-5 else).
declare -A CLUMP_P1=( [EUR]=5e-8 [AFR]=1e-5 [EAS]=1e-5 [AMR]=1e-5 [SAS]=1e-5 )

# --- pick this array task's registry row (skip header) ---
ROW=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${ROW}p" "$REG"); [ -z "$LINE" ] && { echo "no registry row $ROW"; exit 0; }
NAME=$(echo "$LINE" | cut -f1); SSPATH=$(echo "$LINE" | cut -f2); ANC=$(echo "$LINE" | cut -f4)

# Input sumstats = registry sumstats_path column (col 2), taken verbatim. The
# column is FM-relative (e.g. data/sumstats/<study>_{reformatted_hg19,preprocessed}.tsv)
# and resolves directly since we cd'd to $FM; an absolute path would also work as-is.
SS="$SSPATH"

# robust skip: sumstats not generated yet (e.g. MVP files still mid-build)
[ -f "$SS" ] || { echo "SKIP ${NAME}: missing sumstats $SS"; exit 0; }

T1=${TIER1[$ANC]}; CP1=${CLUMP_P1[$ANC]}
[ -z "$T1" ] && { echo "SKIP ${NAME}: unknown ancestry '$ANC' (no Tier-1 threshold)"; exit 0; }
echo "=== ${NAME}  ancestry=${ANC}  tier1=${T1}  clump_p1(loosest)=${CP1} ==="

# ancestry-matched PLINK bed prefix template (uses {CHR}); ABSOLUTE paths:
#   EUR panel lives at project-root data/1kg_eur (NOT under finemapping/data),
#   the others under finemapping/data/ld_ref/1kg_<anc>.
case "$ANC" in
  EUR) BEDT="${ROOT}/data/1kg_eur/chr{CHR}_eur" ;;
  AFR) BEDT="${FM}/data/ld_ref/1kg_afr/chr{CHR}_afr" ;;
  AMR) BEDT="${FM}/data/ld_ref/1kg_amr/chr{CHR}_amr" ;;
  EAS) BEDT="${FM}/data/ld_ref/1kg_eas/chr{CHR}_eas" ;;
  SAS) BEDT="${FM}/data/ld_ref/1kg_sas/chr{CHR}_sas" ;;
  *)   echo "SKIP ${NAME}: unmapped ancestry $ANC"; exit 0 ;;
esac

# PLINK clump input: both allele orderings of the chr:bp:a1:a2 ID so we match the
# bed's variant IDs regardless of A1/A2 strand orientation. Sumstats cols:
#   1 chromosome  2 position  3 allele1  4 allele2  5 beta  6 se  7 pval
awk -F'\t' 'NR>1 { print $1":"$2":"$3":"$4, $7; print $1":"$2":"$4":"$3, $7 }' "$SS" \
  | sort -k1,1 -u | sed '1i SNP P' > "$TMP/clump_in.txt"

> "$TMP/leads.txt"   # accumulates: CHR  BP  P(index)
for CHR in $(seq 1 22); do
  # only bother with chromosomes that carry at least one suggestive (<1e-6) hit
  CHR_SUGG=$(awk -F'\t' -v c="$CHR" -v p="$CP1" 'NR>1 && $1==c && $7<p' "$SS" | wc -l)
  [ "$CHR_SUGG" -eq 0 ] && continue

  BED=${BEDT/\{CHR\}/$CHR}
  if [ ! -f "${BED}.bed" ]; then
    # fallback: no ancestry-matched bed for this chr -> distance-based single best
    echo "  WARNING chr${CHR}: no bed ${BED}.bed; distance-based single best-per-chr"
    awk -F'\t' -v c="$CHR" -v p="$CP1" 'NR>1 && $1==c && $7<p {print $2"\t"$7}' "$SS" \
      | sort -k2,2g | head -1 | awk -v c="$CHR" '{print c"\t"$1"\t"$2}' >> "$TMP/leads.txt"
    continue
  fi

  # Re-key the bed's variant IDs to chr:pos:a1:a2 so they match our clump-input IDs.
  # The ancestry LD beds (ld_ref/1kg_<anc>) carry rsIDs, which do NOT match the
  # chr:bp:a1:a2 IDs we build from the sumstats -> without this, clumping matches
  # zero variants and every non-EUR stratum returns 0 loci. (EUR data/1kg_eur
  # already uses chr:pos:a1:a2, so this re-derivation is a no-op there.)
  RB="$TMP/reid_chr${CHR}"
  awk 'BEGIN{OFS="\t"} {print $1, $1":"$4":"$5":"$6, $3, $4, $5, $6}' "${BED}.bim" > "${RB}.bim"
  ln -sf "$(readlink -f "${BED}.bed")" "${RB}.bed"
  ln -sf "$(readlink -f "${BED}.fam")" "${RB}.fam"

  plink --bfile "$RB" --clump "$TMP/clump_in.txt" \
        --clump-p1 "$CP1" --clump-r2 0.1 --clump-kb 500 --chr "$CHR" \
        --out "$TMP/c${CHR}" --allow-extra-chr 2>/dev/null || true
  if [ -f "$TMP/c${CHR}.clumped" ]; then
    # PLINK 1.9 .clumped cols: CHR(1) F(2) SNP(3=chr:bp:a1:a2) BP(4) P(5) ...
    awk 'NR>1 && $1!="" {print $1"\t"$4"\t"$5}' "$TMP/c${CHR}.clumped" >> "$TMP/leads.txt"
  fi
done

# --- label each index SNP by tier and write the two outputs ---
TIERED="$OUTDIR/${NAME}_leadSNPs_tiered.tsv"
PLAIN="$OUTDIR/${NAME}_leadSNPs.tsv"
echo -e "CHR\tBP\tlocus\tindex_p\ttier" > "$TIERED"
echo -e "CHR\tBP\tlocus" > "$PLAIN"
sort -k1,1n -k2,2n "$TMP/leads.txt" | awk -v k="$T1" -F'\t' '
  $1!="" { tier = ($3+0 < k+0) ? "significant" : "suggestive";
           print $1"\t"$2"\t"$1"."$2"\t"$3"\t"tier >> "'"$TIERED"'";
           print $1"\t"$2"\t"$1"."$2                >> "'"$PLAIN"'" }'

NSIG=$(awk -F'\t' 'NR>1 && $5=="significant"' "$TIERED" | wc -l)
NSUG=$(awk -F'\t' 'NR>1 && $5=="suggestive"'  "$TIERED" | wc -l)
echo "  ${NAME}: significant=${NSIG}  suggestive=${NSUG}  (total loci $((NSIG+NSUG)))"
echo "  wrote: $TIERED"
echo "  wrote: $PLAIN"
rm -rf "$TMP"
