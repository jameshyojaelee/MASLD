#!/bin/bash -l
# run_fm_nonEUR_expanded.sh
# ---------------------------------------------------------------------------
# SuSiE-only fine-mapping of the EXPANDED (significant + suggestive) lead-SNP set
# for all NON-EUR strata (MVP AFR/AMR/EAS + BBJ EAS + PanUKBB AFR/CSA), each with
# its ancestry-matched 1kg LD panel (1kg_afr / 1kg_amr / 1kg_eas / 1kg_csa).
#
# Per study it submits:
#   1. a prep job (02_prep_locus_ss.sh) -> per-locus summary-stat windows, then
#   2. a dependent 03 array (03_run_fm_per_locus.sh) running SuSiE with the
#      convergence + purity + lambda_s (LD-consistency) guards. CARMA is SKIPPED
#      (SKIP_CARMA=1) -- SuSiE + guards suffice for this exploratory discovery pass.
#
# The suggestive-tier / non-EUR credible sets are EXPLORATORY (weak signal + small
# out-of-sample 1kg LD panels: AMR n~347, EAS 504, AFR ~660) and are labeled as
# such downstream via lambda_s + the coloc provenance framework. EUR / MVP_EUR are
# deliberately EXCLUDED (redundant with the already-fine-mapped UKBB/FinnGen EUR).
#
# Usage:
#   bash run_fm_nonEUR_expanded.sh [STUDY_FILTER]
#   QOS_OVERRIDE=interactive bash run_fm_nonEUR_expanded.sh PanUKBB_CSA_ALT   # validate one study now
#
# NOTE: default QOS = nslab; while the MVP COLOC saturates the 7 TB nslab cap these
# jobs PEND until COLOC frees memory (expected). QOS_OVERRIDE=interactive (max 4
# concurrent) is only for validating a small study immediately.
# ---------------------------------------------------------------------------
set -o pipefail
FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "$FM_DIR"
REGISTRY="config/gwas_registry.tsv"
WIN_MB=0.5
STUDY_FILTER="${1:-}"
mkdir -p src/logs

QOS_ARG=""
[ -n "${QOS_OVERRIDE:-}" ] && QOS_ARG="--qos=${QOS_OVERRIDE}"

echo "=== non-EUR expanded SuSiE-only fine-mapping ${QOS_OVERRIDE:+(QOS=$QOS_OVERRIDE)} ==="

tail -n +2 "$REGISTRY" | while IFS=$'\t' read -r study ss_path lead_path anc trait N_tot N_cases ld_panel win; do
  case "$anc" in EAS|AFR|AMR|SAS) ;; *) continue ;; esac        # non-EUR only
  [ -n "$STUDY_FILTER" ] && [ "$study" != "$STUDY_FILTER" ] && continue
  win="${win:-0.5}"
  lead="data/lead_snps/${study}_leadSNPs.tsv"
  [ -f "$lead" ] || { echo "  skip $study: missing $lead"; continue; }
  N_LOCI=$(( $(wc -l < "$lead") - 1 ))
  [ "$N_LOCI" -lt 1 ] && { echo "  skip $study: 0 loci"; continue; }
  anc_lower=$(echo "$anc" | tr '[:upper:]' '[:lower:]')
  LD_POP="1kg_${anc_lower}"
  [ -z "$N_cases" ] && N_cases="NA"

  echo "--- $study | anc=$anc | ${N_LOCI} loci | LD=$LD_POP ---"

  # 1. prep per-locus ss windows. HYBRID: run prep NOW on bigmem + interactive
  #    (bigmem nodes are open + interactive bypasses the COLOC-saturated nslab 7TB
  #    cap). bigmem enforces a 500G per-job MIN under interactive, so size at 500G.
  #    Skip if the windows already exist (idempotent; e.g. the AMR validation run).
  SSDIR="output/${study}/${LD_POP}_${win}Mb/ss"
  N_SS=$(ls "$SSDIR" 2>/dev/null | wc -l)
  DEP=""
  if [ "$N_SS" -ge "$N_LOCI" ]; then
    echo "    prep: SKIP (${N_SS} windows already present)"
  else
    PREP=$(sbatch --parsable --job-name=ssprep --qos=interactive --partition=bigmem --mem=500G \
      src/02_prep_locus_ss.sh "$study" "$ss_path" "$lead" "$LD_POP" "$WIN_MB")
    DEP="--dependency=afterok:${PREP}"
    echo "    prep (bigmem/interactive): $PREP"
  fi

  # 2. 03 array (SuSiE-only) on NSLAB -> pends behind COLOC, then runs at high
  #    concurrency (~200 tasks) once COLOC frees the 7TB cap. LD_PANEL=1kg routes
  #    each ancestry to 1kg_<anc>; the *_LD_DIR overrides are cleared so dispatch
  #    (not a pinned dir) selects the panel.
  ARR=$(sbatch --parsable --job-name=susie --mem=32G $DEP \
    --export=ALL,LD_PANEL=1kg,SKIP_CARMA=1,EAS_LD_DIR=,AFR_LD_DIR=,AMR_LD_DIR=,SAS_LD_DIR= \
    --array=1-${N_LOCI} \
    src/03_run_fm_per_locus.sh "$study" "$LD_POP" "$lead" "$N_tot" "$N_cases" "$WIN_MB" "$anc")
  echo "    03 array (nslab, pend->fast): $ARR (1-$N_LOCI) ${DEP:+[dep]}"
done
echo "=== all matching non-EUR studies submitted ==="
