#!/usr/bin/env Rscript
# ============================================================================
# q2adv_govaere_holdout.R   (Q2 ADVERSARIAL — shared-cohort circularity)
#
# CHARGE: Vacca's EPoS human arm == our Govaere cohort (GSE135251). GSE135251 is
# also the LARGEST cohort in our 5-cohort canonical mega DE (25.5% of 847 samples,
# 359/847). So the Vacca-951 vs our-human concordance could be partly CIRCULAR —
# the same human transcriptomes feeding both sides.
#
# TEST: rebuild our human disease-vs-control signature WITHOUT Govaere and re-run
# the IDENTICAL Vacca-951 overlap (Fisher OR) + fgsea (NES) as q2_contrast_match.
# Govaere-excluded signature = the canonical-method (limma-voom-qw C2) leave-one-
# cohort-out refit already on disk: loo_cv_C2/lvqw_loo_C2_GSE135251.csv
# (4 training cohorts: GSE126848/GSE130970/GSE162694/GSE213621, n=631; Govaere held out).
# Same method/design as canonical_deg_results.csv → an apples-to-apples drop-Govaere.
#
# Head-to-head: FULL (canonical, Govaere IN) vs NOGOV (Govaere OUT) on:
#   (a) Fisher OR + overlap, common universe per (sig × Vacca target)
#   (b) fgsea NES + padj, our t-ranked genes, Vacca target as gene set
# Vacca targets: V_all_951, V_prog_526, V_early_271 (S4 stage partition).
# VERDICT = does the OR/NES survive (ratio NOGOV/FULL) after removing the shared cohort.
# Env: rnaseq
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
set.seed(42)
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO   <- file.path(INT_I, "loo_cv_C2")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
RES   <- file.path(BASE, "results/vacca_benchmark"); dir.create(RES, recursive=TRUE, showWarnings=FALSE)
strip <- function(x) gsub("\\..*", "", x)

# ── Vacca S4: 951 signature partitioned by stage (human symbols) ─────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"), sheet="Table S4"))
setnames(s4, c("GeneSymbol","Gene_used_in_DSEA(0:No/1:Yes)",
               "Early/All/Late disease stage (if used in DSEA otherwise NA)"),
         c("sym","in_dsea","stage"), skip_absent=TRUE)
s4[, sym := toupper(trimws(sym))]
s4_uni  <- unique(s4[!is.na(sym) & sym != "", sym])
v_all   <- unique(s4[in_dsea==1, sym])
v_prog  <- unique(s4[in_dsea==1 & stage=="Disease progression", sym])
v_early <- unique(s4[in_dsea==1 & stage=="Early disease development", sym])
cat(sprintf("Vacca S4 universe=%d | 951=%d | progression=%d | early=%d\n",
            length(s4_uni), length(v_all), length(v_prog), length(v_early)))
vacca_targets <- list(V_all_951=v_all, V_prog_526=v_prog, V_early_271=v_early)

# ── Symbol map (ENSEMBL->symbol) from canonical ──────────────────────────────
map <- fread(file.path(INT_I, "canonical_deg_results.csv"))[, .(gene=strip(gene), symbol)]
map <- map[!duplicated(gene)]

# normalize a DE table -> (sym, logFC, t, padj), dedup by symbol keeping strongest |t|
load_sig <- function(path, has_symbol) {
  d <- fread(path)
  d[, gene := strip(gene)]
  if (!has_symbol) d <- merge(d, map, by="gene", all.x=TRUE)
  out <- d[, .(sym=toupper(trimws(symbol)),
               logFC=as.numeric(logFC), t=as.numeric(t), padj=as.numeric(padj))]
  out <- out[!is.na(sym) & sym!="" & is.finite(t)]
  out[, at := abs(t)]; setorder(out, -at); out[, at := NULL]   # keep strongest |t| per symbol
  out[!duplicated(sym)]
}
SIGS <- list(
  FULL  = load_sig(file.path(INT_I, "canonical_deg_results.csv"),        has_symbol=TRUE),
  NOGOV = load_sig(file.path(LOO,   "lvqw_loo_C2_GSE135251.csv"),         has_symbol=FALSE)
)
for (nm in names(SIGS)) cat(sprintf("%-6s n_tested=%d  n_sig(padj<.05,|lfc|>.5)=%d\n",
   nm, nrow(SIGS[[nm]]), SIGS[[nm]][padj<0.05 & abs(logFC)>0.5, .N]))

# ── (a) Fisher OR + overlap, common universe per (sig × vacca target) ────────
fisher_block <- rbindlist(lapply(names(SIGS), function(sn) {
  hsig <- SIGS[[sn]]
  rbindlist(lapply(names(vacca_targets), function(vn) {
    vt  <- vacca_targets[[vn]]
    uni <- intersect(hsig$sym, s4_uni)
    sig <- intersect(hsig[padj<0.05 & abs(logFC)>0.5, sym], uni)
    tgt <- intersect(vt, uni)
    a <- length(intersect(sig, tgt)); b <- length(setdiff(sig, tgt))
    c <- length(setdiff(tgt, sig));  d <- length(uni) - a - b - c
    ft <- fisher.test(matrix(c(a,b,c,d), 2), alternative="greater")
    data.table(sig_set=sn, vacca_target=vn, universe=length(uni),
               n_our_sig=length(sig), n_vacca=length(tgt), overlap=a,
               OR=unname(ft$estimate), fisher_p=ft$p.value, jaccard=a/length(union(sig,tgt)))
  }))
}))

# ── (b) fgsea: our t-ranked genes, Vacca target = gene set ───────────────────
gsea_block <- rbindlist(lapply(names(SIGS), function(sn) {
  hsig <- SIGS[[sn]]; st <- setNames(hsig$t, hsig$sym); st <- st[is.finite(st)]
  sets <- lapply(vacca_targets, function(vt) intersect(vt, names(st)))
  fg <- as.data.table(suppressWarnings(fgsea(sets, st, scoreType="std", nPermSimple=10000)))
  fg[, sig_set := sn][, .(sig_set, vacca_target=pathway, set_size=size,
                          NES=NES, gsea_padj=padj, leadingEdge_n=lengths(leadingEdge))]
}))

res <- merge(fisher_block, gsea_block, by=c("sig_set","vacca_target"))
setorder(res, vacca_target, sig_set)

cat("\n=== (a)+(b) FULL (Govaere IN) vs NOGOV (Govaere OUT) ===\n")
print(res[, .(vacca_target, sig_set, universe, overlap, OR=round(OR,2),
              fisher_p=signif(fisher_p,3), NES=round(NES,3), gsea_padj=signif(gsea_padj,3))])

# ── survival ratios NOGOV/FULL per Vacca target ──────────────────────────────
wide <- dcast(res, vacca_target ~ sig_set, value.var=c("overlap","OR","NES","gsea_padj","n_our_sig"))
wide[, OR_survival   := round(OR_NOGOV   / OR_FULL,   3)]
wide[, NES_survival  := round(NES_NOGOV  / NES_FULL,  3)]
wide[, overlap_ratio := round(overlap_NOGOV / overlap_FULL, 3)]
cat("\n=== SURVIVAL: NOGOV / FULL ratios per Vacca target ===\n")
print(wide[, .(vacca_target,
               OR_FULL=round(OR_FULL,2),  OR_NOGOV=round(OR_NOGOV,2),  OR_survival,
               NES_FULL=round(NES_FULL,2),NES_NOGOV=round(NES_NOGOV,2),NES_survival,
               overlap_FULL=overlap_FULL, overlap_NOGOV=overlap_NOGOV, overlap_ratio)])

# ── gene-level robustness: does dropping Govaere move the human ranking? ──────
m <- merge(SIGS$FULL[, .(sym, t_full=t, lfc_full=logFC)],
           SIGS$NOGOV[, .(sym, t_nogov=t, lfc_nogov=logFC)], by="sym")
rho_t   <- cor(m$t_full,   m$t_nogov,   method="spearman", use="complete.obs")
rho_lfc <- cor(m$lfc_full, m$lfc_nogov, method="spearman", use="complete.obs")
dir_conc<- mean(sign(m$t_full)==sign(m$t_nogov))
cat(sprintf("\n=== Human-signature stability (drop-Govaere vs full) over %d shared genes ===\n", nrow(m)))
cat(sprintf("Spearman t = %.4f | Spearman logFC = %.4f | sign(t) concordance = %.3f\n",
            rho_t, rho_lfc, dir_conc))

# ── write ────────────────────────────────────────────────────────────────────
fwrite(res,  file.path(RES, "q2adv_govaere_holdout.csv"))
fwrite(res,  file.path(OUT, "q2adv_govaere_holdout.csv"))
fwrite(wide, file.path(RES, "q2adv_govaere_holdout_survival.csv"))
fwrite(data.table(metric=c("spearman_t_full_vs_nogov","spearman_lfc_full_vs_nogov",
                           "sign_t_concordance","n_shared_genes"),
                  value=c(rho_t, rho_lfc, dir_conc, nrow(m))),
       file.path(RES, "q2adv_govaere_holdout_stability.csv"))

cat("\nWrote q2adv_govaere_holdout.csv / _survival.csv / _stability.csv\n")
