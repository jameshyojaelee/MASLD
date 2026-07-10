#!/bin/bash -l
#SBATCH --job-name=PLINK
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-32
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/leadsnp_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/leadsnp_%A_%a.err

# 02_mvp_lead_snps_twotier.sh
#
# ============================ DEPRECATED — DO NOT RUN (2026-07-02) ============================
# Superseded by the production scripts/mvp/02_lead_snps_twotier_all.sh, which implements the
# LOCKED per-ancestry scheme: Tier-1 = min(5e-8, Kanai) [EUR/EAS/AMR/SAS = 5e-8, AFR = 3.24e-8];
# Tier-2 "suggestive" = 1e-5 NON-EUR ONLY. This template uses WRONG thresholds (raw Kanai
# EUR 9.26e-8, and SUGG=1e-6 applied to ALL ancestries incl. EUR — both too lenient / wrong tier).
# Retained for provenance only. See docs/reference/causal_threshold_logic.md (Layer 1).
# =============================================================================================
#
# Two-tier lead-SNP identification for MVP strata (standalone; does NOT touch the
# canonical 01_identify_lead_snps.sh used by the existing 23 GWAS).
#   Tier 1 "significant": index-SNP p < per-ancestry Kanai 2016 genome-wide threshold
#                         (AFR 3.24e-8 / EUR 9.26e-8 / EAS 1.61e-7 / AMR 1.83e-7)
#   Tier 2 "suggestive" : index-SNP p < 1e-6 (exploratory; for underpowered AFR/AMR/EAS)
# One PLINK clump pass at p1=1e-6 (loosest), then label each locus by tier.
# LD clumping uses the ancestry-matched 1000G PLINK panel.
# Outputs:
#   data/lead_snps/MVP_<name>_leadSNPs.tsv         (CHR BP locus  — registry/fine-map compat, union @1e-6)
#   data/lead_snps/MVP_<name>_leadSNPs_tiered.tsv  (CHR BP locus index_p tier — reporting deliverable)

set -eo pipefail
module load PLINK/1.9

FM="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
MAN="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/mvp_manifest.tsv"
cd "$FM"
OUTDIR="data/lead_snps"; TMP="data/lead_snps/tmp_mvp_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$OUTDIR" "$TMP"

# Kanai 2016 (J Hum Genet) per-ancestry genome-wide significance thresholds
declare -A KANAI=( [EUR]=9.26e-8 [AFR]=3.24e-8 [AMR]=1.83e-7 [EAS]=1.61e-7 )
SUGG=1e-6

ROW=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${ROW}p" "$MAN"); [ -z "$LINE" ] && { echo "no row $ROW"; exit 0; }
NAME=$(echo "$LINE" | cut -f1); ANC=$(echo "$LINE" | cut -f5)
SS="data/sumstats/${NAME}_reformatted_hg19.tsv"
[ -f "$SS" ] || { echo "missing sumstats $SS"; exit 1; }
KTHR=${KANAI[$ANC]}
echo "=== $NAME  ancestry=$ANC  Kanai=$KTHR  suggestive=$SUGG ==="

# ancestry-matched PLINK bed prefix template (uses {CHR})
case "$ANC" in
  EUR) BEDT="data/1kg_eur/chr{CHR}_eur" ;;
  AFR) BEDT="data/ld_ref/1kg_afr/chr{CHR}_afr" ;;
  AMR) BEDT="data/ld_ref/1kg_amr/chr{CHR}_amr" ;;
  EAS) BEDT="data/ld_ref/1kg_eas/chr{CHR}_eas" ;;
  *) echo "unknown ancestry $ANC"; exit 1 ;;
esac

# PLINK clump input: both allele orderings so we match regardless of strand
awk -F'\t' 'NR>1 { print $1":"$2":"$3":"$4, $7; print $1":"$2":"$4":"$3, $7 }' "$SS" \
  | sort -k1,1 -u | sed '1i SNP P' > "$TMP/clump_in.txt"

> "$TMP/leads.txt"   # CHR BP P(index)
for CHR in $(seq 1 22); do
  CHR_SUGG=$(awk -F'\t' -v c="$CHR" -v p="$SUGG" 'NR>1 && $1==c && $7<p' "$SS" | wc -l)
  [ "$CHR_SUGG" -eq 0 ] && continue
  BED=${BEDT/\{CHR\}/$CHR}
  if [ ! -f "${BED}.bed" ]; then
    # fallback: distance-based single best per chr
    awk -F'\t' -v c="$CHR" -v p="$SUGG" 'NR>1 && $1==c && $7<p {print $2"\t"$7}' "$SS" \
      | sort -k2,2g | head -1 | awk -v c="$CHR" '{print c"\t"$1"\t"$2}' >> "$TMP/leads.txt"
    continue
  fi
  plink --bfile "$BED" --clump "$TMP/clump_in.txt" \
        --clump-p1 "$SUGG" --clump-r2 0.1 --clump-kb 500 --chr "$CHR" \
        --out "$TMP/c${CHR}" --allow-extra-chr 2>/dev/null || true
  if [ -f "$TMP/c${CHR}.clumped" ]; then
    # .clumped cols: CHR(1) F(2) SNP(3=chr:bp:a1:a2) BP(4) P(5) ...
    awk 'NR>1 && $1!="" {print $1"\t"$4"\t"$5}' "$TMP/c${CHR}.clumped" >> "$TMP/leads.txt"
  fi
done

# Label tier by index-SNP p, write outputs
TIERED="$OUTDIR/${NAME}_leadSNPs_tiered.tsv"
PLAIN="$OUTDIR/${NAME}_leadSNPs.tsv"
echo -e "CHR\tBP\tlocus\tindex_p\ttier" > "$TIERED"
echo -e "CHR\tBP\tlocus" > "$PLAIN"
sort -k1,1n -k2,2n "$TMP/leads.txt" | awk -v k="$KTHR" -F'\t' '
  $1!="" { tier = ($3+0 < k+0) ? "significant" : "suggestive";
           print $1"\t"$2"\t"$1"."$2"\t"$3"\t"tier >> "'"$TIERED"'";
           print $1"\t"$2"\t"$1"."$2 >> "'"$PLAIN"'" }'
NSIG=$(awk -F'\t' 'NR>1 && $5=="significant"' "$TIERED" | wc -l)
NSUG=$(awk -F'\t' 'NR>1 && $5=="suggestive"'  "$TIERED" | wc -l)
echo "  $NAME: significant=$NSIG  suggestive=$NSUG  (total loci $((NSIG+NSUG)))"
rm -rf "$TMP"
