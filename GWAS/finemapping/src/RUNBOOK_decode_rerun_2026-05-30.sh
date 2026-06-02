#!/bin/bash
# RUNBOOK — deCODE/Sveinbjornsson-2022 NAFLD re-download + proper re-run (2026-05-30/31)
# =============================================================================================
# WHY (corrected 2026-05-31 after inspecting the authoritative files): the old deCode + Intermountain
# *_preprocessed.tsv sumstats were unusable because (1) the preprocessing put `Effect` (an ODDS RATIO)
# straight into a `beta` field with NO log() transform, and (2) deCODE's own release fills unestimable
# ultra-rare variants with a placeholder Effect (deCODE/Intermountain: OR=0.018 on ~31%/21% of rows at
# MAF~1e-5; UKBB: OR=1.000, p=1) that was never filtered. So the "beta=0.018 sentinel / |beta|~102"
# was placeholder ODDS RATIOS, not a download artifact. They anchored coloc_best_gwas for ~1,903/18,975
# genes (~10% of gene-level COLOC). Source: Sveinbjornsson et al. 2022 Nat Genet, PMID 36280732.
# FIX (validated against the data): reformat_decode() does beta=log(Effect), back-computes SE from
# beta+pval (case-control files have no SE), af=MAF_PC/100; the QC gate then drops the placeholder via
# sentinel detection (log(0.018) repeated >1%) AND MAF>=1% (confirmed: 0 placeholders among MAF>=1%
# common variants in deCODE; common-variant ORs are sane, range [0.225,1.998]). UKBB OR=1.000 -> beta=0
# -> dropped by the gate's beta!=0.
#
# SLURM note (2026-05-31, per PI): default to the nslab QOS. Use --qos=interactive ONLY if the 7TB
# nslab cap blocks a job (reason QOSMaxMemoryPerUser). Do NOT run heavy work on the login node.
# =============================================================================================
set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$BASE/GWAS/finemapping/src           # (this RUNBOOK + 00 live on branch gwas-review-fixes-2026-05-30)
DECODE_DIR=$BASE/GWAS/MR_Data/Sveinbjornsson2022   # NOTE: GWAS/MR_Data -> GWAS/data rename pending (post-runs)

echo "This is a runbook. Read it; run the steps deliberately. It will not auto-execute heavy jobs."
exit 0   # safety stop — remove only if you intend to source individual blocks

# ---------------------------------------------------------------------------------------------
# STEP 0 — download the authoritative bundle.  [DONE 2026-05-31 — array job 16938648]
# ---------------------------------------------------------------------------------------------
#   decode.com/summarydata -> "Multiomics study of NAFLD" data-use form -> signed token.
#   Parallel array download: $DECODE_DIR/_download_bundle_array.sh (one file per task, io partition).
#   12 files landed (~63GB): NAFL_{deCODE,INTERMOUNTAIN,UKBB}_sumstat.txt + NAFL_FINNGEN_sumstats.txt,
#   Cirrhosis_{deCODE,INTERMOUNTAIN,UKBB}_sumstat.txt + Cirrhosis_FINNGEN_sumstats.txt,
#   HCC_{deCODE,INTERMOUNTAIN,UKBB}_sumstat.txt, Proton_density_fat_fraction_UKBB_sumstats.txt.
#   No md5 sidecars exist at the token path (.tmd5sum + NAFLD2022.zmd5sum both 404); verified via
#   header sanity + row counts + the QC sentinel detector. readme: $DECODE_DIR/nafld2022readme.txt.

# ---------------------------------------------------------------------------------------------
# STEP 1 — reformat_decode() in 00_reformat_gwas.R.  [DONE 2026-05-31]
# ---------------------------------------------------------------------------------------------
#   Format confirmed from readme: case-control NAFL = Pval, Effect(OR), MAF_PC(%), Info, Chrom, Pos(hg38),
#   Amin(minor/effect), Amaj(major/other); NO SE. reformat_decode handles it (beta=log(OR), SE back-comp,
#   af=MAF_PC/100, hg38->hg19, QC gate). Targeted re-run via REFORMAT_SET=decode (only the 3 deCODE NAFL).

# ---------------------------------------------------------------------------------------------
# STEP 2 — fix the registry N/case counts + repoint paths.  [PENDING — do before STEP 4 COLOC]
#   The OLD registry N are wrong. Intermountain = 2,301 cases / 24,715 controls (CMDKP), NOT 111,000/3,757.
#   CONFIRM deCode + UKBB case/control N from the Sveinbjornsson-2022 paper (readme gives columns, not N).
# ---------------------------------------------------------------------------------------------
#   $SRC/../config/gwas_registry.tsv rows: 2023_36280732_NAFLD_{deCode,Intermountain,UKBB}_EUR
#     N_tot / N_cases  -> authoritative values ;  sumstats_path -> the *_reformatted_hg19.tsv from STEP 3.

# ---------------------------------------------------------------------------------------------
# STEP 3 — reformat (liftover hg38->hg19 + QC gate).  [RUNNING 2026-05-31 — job 16938684]
# ---------------------------------------------------------------------------------------------
sbatch --partition=cpu --mem=96G --cpus-per-task=4 --time=48:00:00 \
  --job-name=reformat_decode --output=$BASE/GWAS/finemapping/logs/reformat_decode_%j.out \
  --wrap="REFORMAT_SET=decode micromamba run -n rnaseq Rscript $SRC/00_reformat_gwas.R"
#   VERIFY the QC log: sentinel warning FIRES for deCode/Intermountain (drops log(0.018)); sane common-
#   variant count (~8-10M) for the 3 deCODE cohorts; UKBB OR=1.000 placeholders gone via beta!=0.

# ---------------------------------------------------------------------------------------------
# STEP 4 — re-run SuSiE/ABF COLOC for the 3 sub-cohorts (LD_PANEL=polyfun, EUR). 22-chr array each.
#          AFTER STEP 2 (registry N) + STEP 3 (reformatted files exist). bigmem, nslab QOS.
# ---------------------------------------------------------------------------------------------
for G in 2023_36280732_NAFLD_deCode_EUR 2023_36280732_NAFLD_Intermountain_EUR 2023_36280732_NAFLD_UKBB_EUR; do
  GWAS_NAME=$G LD_PANEL=polyfun sbatch --partition=bigmem --mem=128G \
    --time=48:00:00 $SRC/06_susie_coloc.sh
done
#   If a job pends with reason QOSMaxMemoryPerUser (7TB nslab cap hit), add --qos=interactive (max 4
#   concurrent). Otherwise leave it on nslab.

# ---------------------------------------------------------------------------------------------
# STEP 5 — re-aggregate gene-level COLOC (picks up the corrected deCODE/Intermountain). io, nslab QOS.
# ---------------------------------------------------------------------------------------------
sbatch --partition=io --mem=32G --cpus-per-task=4 --time=8:00:00 \
  --job-name=07_refresh --wrap="cd $BASE/GWAS/finemapping && micromamba run -n rnaseq Rscript src/07_combine_susie_coloc.R"
#   VERIFY: the ~1,903 genes that had coloc_best_gwas = deCode/Intermountain now draw from CLEAN
#   sumstats; SuSiE/ABF threshold counts in gene_level_coloc.csv shift accordingly.

# ---------------------------------------------------------------------------------------------
# STEP 6 — atlas rebuild (27a->...->75->217->27b) is OWNED BY THE scATAC/atlas SESSION. Notify them
#          to rebuild once gene_level_coloc.csv is refreshed. Then re-run 46d -> convergence_evidence.
#          (Do NOT run RNA-seq/run_atlas_rebuild.sh from this gwas branch — scATAC owns it.)
# ---------------------------------------------------------------------------------------------
echo "RUNBOOK complete — coordinate the atlas rebuild with the scATAC session."

# ---------------------------------------------------------------------------------------------
# ASIDE — overlap/difference analysis of the full bundle vs our existing GWAS (job 16938685):
#   $SRC/00c_compare_overlapping_gwas.R -> results/gwas_overlap/overlap_comparison.csv.
#   Pairs: Sveinb FinnGen NAFL vs our R12 NAFLD (release diff); Sveinb PDFF-UKBB vs Pazoki PDFF;
#   Sveinb deCODE/UKBB NAFL vs Ghodsian meta (Ghodsian INCLUDES deCODE+UKBB -> sample overlap, so they
#   are NOT independent COLOC evidence — null-Z corr quantifies it). Decides what to add vs what double-counts.
# ---------------------------------------------------------------------------------------------
