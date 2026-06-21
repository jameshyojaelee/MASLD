#!/usr/bin/env Rscript
# ============================================================================
# q2_continuous_concordance_vs_vacca_prog.R
#
# Q2: The 172/1323 gene overlap + the conserved-core are bottlenecked by arbitrary
# human+mouse LFC/padj cutoffs. Here we go CUTOFF-FREE: correlate OUR continuous
# per-gene cross-species concordance metrics against Vacca's continuous human
# progression logFC (MOESM4, mean Severe-vs-Mild) over ALL shared genes.
#
# OUR continuous concordance metrics tested:
#   (A) translatability_score  -- pre-computed (concordance_atlas_unified.csv); unsigned
#   (B) hm_agreement (computed here) -- signed human x mouse significance agreement:
#         human_ss  = sign(t_human) * -log10(P_human)
#         mouse_ss  = mean_over_4_diets( sign(t_diet) * -log10(P_diet) )
#         hm_agreement = human_ss * mouse_ss   (>0 same-direction, magnitude = joint signif)
#       This is a continuous concordance: high+ = both species strongly agree.
#
# Vacca target (continuous, NO cutoffs):
#   vacca_prog_lfc = rowMeans(L2FC_UCAM/VCU:Severe vs Mild, L2FC_EPoS:Severe vs Mild)
#
# We report Spearman + Pearson + n vs BOTH the signed Vacca logFC and its abs value.
#   - translatability (unsigned) is naturally compared to |vacca_prog_lfc|.
#   - hm_agreement (signed concordance) is compared to both.
#   - we also build hm_signed_human = human_ss alone and a directional concordance to
#     sanity-check that the signed axis tracks Vacca's signed progression.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CONC  <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
ANNOT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
MPD   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
strip <- function(x) gsub("\\..*", "", x)
sp <- function(a,b) suppressWarnings(cor(a,b,method="spearman",use="complete.obs"))
pe <- function(a,b) suppressWarnings(cor(a,b,method="pearson", use="complete.obs"))
nfin <- function(a,b) sum(is.finite(a) & is.finite(b))

# ── Vacca human progression logFC (continuous, no cutoff) ────────────────────
s4 <- as.data.table(read_excel(file.path(VACCA,"42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet="Table S4"))
prog_cols <- c("L2FC_UCAM/VCU: Severe vs Mild","L2FC_EPoS: Severe vs Mild")
for (cc in prog_cols) suppressWarnings(s4[, (cc) := as.numeric(get(cc))])
s4[, vacca_prog_lfc := rowMeans(as.matrix(.SD), na.rm=TRUE), .SDcols=prog_cols]
vac <- s4[, .(human_symbol = GeneSymbol, vacca_prog_lfc)]
vac <- vac[!is.na(human_symbol) & human_symbol != "" & is.finite(vacca_prog_lfc)]
vac <- vac[!duplicated(human_symbol)]
cat(sprintf("Vacca progression logFC: %d genes with finite mean Severe-vs-Mild\n", nrow(vac)))

# ── (A) translatability_score from unified atlas ─────────────────────────────
atl <- fread(file.path(CONC,"concordance_atlas_unified.csv"))
atl <- atl[!duplicated(human_symbol)]
cat(sprintf("Atlas: %d genes; translatability_score finite=%d\n",
            nrow(atl), sum(is.finite(atl$translatability_score))))

# ── (B) compute hm_agreement from human t + per-diet mouse t ─────────────────
hum <- fread(file.path(INT_I,"canonical_deg_results.csv"))
hum[, hbase := strip(gene)]
hum[, human_ss := sign(t) * -log10(pmax(P.Value, .Machine$double.xmin))]
ortho <- fread(file.path(ANNOT,"ortholog_mapping.tsv"))
ortho[, `:=`(hbase = strip(human_gene_id), mbase = strip(mouse_gene_id))]
ortho1 <- ortho[!duplicated(mbase) & !duplicated(hbase)]   # 1:1 backbone for the join

diets <- c("MCD","HFD","CDAHFD","FPC")
mouse_ss <- NULL
for (d in diets) {
  md <- fread(file.path(MPD, paste0(d,"_de_results.csv")))
  md[, mbase := strip(gene)]
  md[, ss := sign(t) * -log10(pmax(P.Value, .Machine$double.xmin))]
  sub <- md[is.finite(ss) & !duplicated(mbase), .(mbase, ss)]
  setnames(sub, "ss", paste0("ss_", d))
  mouse_ss <- if (is.null(mouse_ss)) sub else merge(mouse_ss, sub, by="mbase", all=TRUE)
}
ss_cols <- paste0("ss_", diets)
mouse_ss[, mouse_ss := rowMeans(as.matrix(.SD), na.rm=TRUE), .SDcols=ss_cols]

# bridge mouse->human via 1:1 ortholog, attach human symbol + human_ss
mh <- merge(mouse_ss[, .(mbase, mouse_ss)], ortho1[, .(mbase, hbase, human_symbol)], by="mbase")
mh <- merge(mh, hum[, .(hbase, human_ss, human_symbol_h = symbol)], by="hbase")
mh[, hm_agreement := human_ss * mouse_ss]
mh <- mh[!duplicated(human_symbol)]
cat(sprintf("hm_agreement computed for %d genes (1:1 ortholog, finite both arms=%d)\n",
            nrow(mh), nfin(mh$human_ss, mh$mouse_ss)))

# ── Merge each metric with Vacca over SHARED genes (no cutoffs) ───────────────
mA <- merge(vac, atl[, .(human_symbol, translatability_score, mean_h_lfc,
                         n_concordant, n_diets_sig)], by="human_symbol")
mB <- merge(vac, mh[, .(human_symbol, hm_agreement, human_ss, mouse_ss)], by="human_symbol")

cat(sprintf("\nSHARED gene counts: translatability∩Vacca=%d ; hm_agreement∩Vacca=%d\n",
            nrow(mA), nrow(mB)))

# Helper: emit a result row
mk <- function(metric, target_lab, x, y) data.table(
  metric=metric, vacca_target=target_lab,
  spearman=round(sp(x,y),4), pearson=round(pe(x,y),4), n=nfin(x,y))

res <- rbindlist(list(
  # (A) translatability_score (unsigned magnitude of concordance)
  mk("translatability_score","vacca_prog_lfc_signed", mA$translatability_score, mA$vacca_prog_lfc),
  mk("translatability_score","abs_vacca_prog_lfc",    mA$translatability_score, abs(mA$vacca_prog_lfc)),
  # (B) hm_agreement (signed human x mouse significance product)
  mk("hm_agreement","vacca_prog_lfc_signed", mB$hm_agreement, mB$vacca_prog_lfc),
  mk("hm_agreement","abs_vacca_prog_lfc",    mB$hm_agreement, abs(mB$vacca_prog_lfc)),
  # sanity: our human signed significance alone vs Vacca signed logFC (should track)
  mk("human_ss_only","vacca_prog_lfc_signed", mB$human_ss, mB$vacca_prog_lfc),
  # sanity: our mouse mean signed significance vs Vacca signed (cross-species directional)
  mk("mouse_ss_only","vacca_prog_lfc_signed", mB$mouse_ss, mB$vacca_prog_lfc),
  # extra continuous concordance proxies from atlas (signed n_concordant scaling)
  mk("atlas_mean_h_lfc","vacca_prog_lfc_signed", mA$mean_h_lfc, mA$vacca_prog_lfc),
  mk("atlas_n_concordant","abs_vacca_prog_lfc",  mA$n_concordant, abs(mA$vacca_prog_lfc))
))

cat("\n=== Q2: continuous per-gene concordance vs Vacca progression logFC (NO cutoffs) ===\n")
print(res)
fwrite(res, file.path(OUT,"q2_continuous_concordance_vs_vacca_prog.csv"))

# Also persist the merged per-gene table for the two PRIMARY metrics
prim <- merge(
  mA[, .(human_symbol, vacca_prog_lfc, translatability_score, mean_h_lfc)],
  mB[, .(human_symbol, hm_agreement, human_ss, mouse_ss)],
  by="human_symbol", all=TRUE)
fwrite(prim, file.path(OUT,"q2_pergene_merged.csv"))
cat(sprintf("\nWrote %s\nWrote %s (per-gene, %d rows)\n",
            file.path(OUT,"q2_continuous_concordance_vs_vacca_prog.csv"),
            file.path(OUT,"q2_pergene_merged.csv"), nrow(prim)))
cat(strrep("=",78),"\n")
