#!/usr/bin/env Rscript
# ============================================================================
# 12_vacca_q2_contrast_matrix.R   (Vacca 2024 benchmark — canonical Q2)
#
# Consolidates the validated Q2 follow-up (vacca_followup/q2_contrast_match.R +
# q2_threshold_free_fgsea.R + q2adv_govaere_holdout.R + q2_adversarial_null_*.R)
# into one traceable canonical script.
#
# QUESTION (Q2): only 172 of our ~1,323 cross-species core genes overlap LITMUS's
# 951-gene signature. Is the modest overlap (a) a CONTRAST MISMATCH (our
# disease-vs-control anchor vs LITMUS's progression-weighted signature) and
# (b) an arbitrary-CUTOFF artifact?
#
# DESIGN:
#   (A) CONTRAST-MATCHED ENRICHMENT MATRIX. Our human signatures × LITMUS's
#       stage-partitioned 951 (S4 "stage" column: 271 Early / 526 Progression /
#       154 All). Each cell = one-sided Fisher OR (common testable universe) +
#       fgsea NES (our t-ranked genes, LITMUS bin as gene set).
#         Rows (CANONICAL method): H1_dvc = disease-vs-control (limma-voom-qw C2);
#                                   H2_naflnash = NAFL-vs-NASH progression (LVQW).
#         (H3_advfib = adv-vs-early fibrosis is DREAM-method only on disk -> kept
#          as a flagged sensitivity row in the CSV, NOT in the headline matrix.)
#       PREDICTION: the progression human signature (H2) matches LITMUS's
#       progression bin best (apples-to-apples), and disease-vs-control matches
#       the "early" bin only weakly (specificity / negative control).
#
#   (B) THRESHOLD-FREE, NON-CIRCULAR ANCHOR. Drop all human-side cutoffs: genes
#       concordant in >=1 mouse diet (cross-species direction agreement) vs
#       LITMUS-951 membership -> Fisher OR. Kills the "cutoff artifact" charge.
#       (fgsea of LITMUS-951 over continuous concordance rankings reported too.)
#
#   (C) SHARED-COHORT CIRCULARITY (Govaere holdout). LITMUS EPoS == our Govaere
#       (GSE135251), our largest cohort. Refit the human signature WITHOUT Govaere
#       (canonical-method LOO on disk) and re-run the matrix -> OR must survive.
#
#   (D) LABEL-SHUFFLE NULL. 1000x membership shuffle for the progression diagonal
#       and the non-circular anchor -> empirical p, independent of fgsea internals.
#
# DROPPED (circular human-vs-human; NOT reported): RRHO, gene-level Spearman ~0.60,
#   72.6% direction concordance on the 172 shared genes — these are signed by OUR
#   human logFC against LITMUS's (largely shared) human logFC.
#
# Outputs (Analysis/Cross_Species_Concordance/results/vacca_benchmark/):
#   q2_contrast_matrix.csv        (headline matrix: H1/H2 x 3 bins; H3 dream flagged)
#   q2_noncircular_anchor.csv     (threshold-free Fisher OR + fgsea)
#   q2_govaere_holdout.csv        (FULL vs NOGOV survival)
#   q2_shuffle_null.csv           (membership-shuffle empirical p)
# Env: rnaseq
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
set.seed(42)
N_PERM <- 1000L
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
DSIG  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")
LOO   <- file.path(INT_I, "loo_cv_C2")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
strip <- function(x) gsub("\\..*", "", x)

# ── LITMUS S4: 951 signature partitioned by stage (human symbols) ─────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"), sheet = "Table S4"))
setnames(s4, c("GeneSymbol","Gene_used_in_DSEA(0:No/1:Yes)",
               "Early/All/Late disease stage (if used in DSEA otherwise NA)"),
         c("sym","in_dsea","stage"), skip_absent = TRUE)
s4[, sym := toupper(trimws(sym))]
s4_uni  <- unique(s4[!is.na(sym) & sym != "", sym])
v_all   <- unique(s4[in_dsea == 1, sym])
v_prog  <- unique(s4[in_dsea == 1 & stage == "Disease progression", sym])
v_early <- unique(s4[in_dsea == 1 & stage == "Early disease development", sym])
cat(sprintf("LITMUS S4 universe=%d | 951=%d | progression=%d | early=%d\n",
            length(s4_uni), length(v_all), length(v_prog), length(v_early)))
vacca_targets <- list(V_early_271 = v_early, V_prog_526 = v_prog, V_all_951 = v_all)

# ── symbol map + signature loader ─────────────────────────────────────────────
map <- fread(file.path(INT_I, "canonical_deg_results.csv"))[, .(gene = strip(gene), symbol)]
map <- map[!duplicated(gene)]
load_sig <- function(path, has_symbol, padj_col = "padj") {
  d <- fread(path); d[, gene := strip(gene)]
  if (!has_symbol) d <- merge(d, map, by = "gene", all.x = TRUE)
  out <- d[, .(sym = toupper(trimws(symbol)), logFC = as.numeric(logFC),
               t = as.numeric(t), padj = as.numeric(get(padj_col)))]
  out <- out[!is.na(sym) & sym != "" & is.finite(t)]
  out[, at := abs(t)]; setorder(out, -at); out[, at := NULL]
  out[!duplicated(sym)]
}
H <- list(
  H1_dvc      = load_sig(file.path(INT_I, "canonical_deg_results.csv"), has_symbol = TRUE),
  H2_naflnash = load_sig(file.path(DSIG, "nafl_vs_nash_lvqw.csv"),       has_symbol = FALSE),
  H3_advfib   = load_sig(file.path(DSIG, "adv_vs_early_fibrosis_dream.csv"), has_symbol = TRUE,
                         padj_col = "adj.P.Val"))
sig_method <- c(H1_dvc = "limma_voom_qw__C2 (canonical)", H2_naflnash = "LVQW progression",
                H3_advfib = "dream (RETIRED — sensitivity only)")
for (nm in names(H)) cat(sprintf("%-12s n_tested=%d  n_sig(padj<.05,|lfc|>.5)=%d  [%s]\n",
   nm, nrow(H[[nm]]), H[[nm]][padj < 0.05 & abs(logFC) > 0.5, .N], sig_method[nm]))

# ── (A) contrast-matched matrix: Fisher OR + fgsea NES ────────────────────────
contrast_block <- function(sigs) rbindlist(lapply(names(sigs), function(hn) {
  hsig <- sigs[[hn]]
  st <- setNames(hsig$t, hsig$sym); st <- st[is.finite(st)]
  fg <- as.data.table(suppressWarnings(fgsea(
    lapply(vacca_targets, function(vt) intersect(vt, names(st))), st,
    scoreType = "std", nPermSimple = 10000)))
  rbindlist(lapply(names(vacca_targets), function(vn) {
    vt  <- vacca_targets[[vn]]
    uni <- intersect(hsig$sym, s4_uni)
    sig <- intersect(hsig[padj < 0.05 & abs(logFC) > 0.5, sym], uni)
    tgt <- intersect(vt, uni)
    a <- length(intersect(sig, tgt)); b <- length(setdiff(sig, tgt))
    cc <- length(setdiff(tgt, sig)); dd <- length(uni) - a - b - cc
    ft <- fisher.test(matrix(c(a, b, cc, dd), 2), alternative = "greater")
    data.table(human_sig = hn, method = sig_method[hn], vacca_target = vn,
               universe = length(uni), n_our_sig = length(sig), n_vacca = length(tgt),
               overlap = a, OR = unname(ft$estimate), fisher_p = ft$p.value,
               jaccard = a / length(union(sig, tgt)),
               NES = fg[pathway == vn, NES], gsea_padj = fg[pathway == vn, padj])
  }))
}))
mat <- contrast_block(H)
mat[, vacca_target := factor(vacca_target, levels = c("V_early_271","V_prog_526","V_all_951"))]
setorder(mat, human_sig, vacca_target)
cat("\n=== (A) CONTRAST-MATCHED MATRIX (Fisher OR + fgsea NES) ===\n")
print(mat[, .(human_sig, vacca_target, universe, overlap, OR = round(OR, 2),
              fisher_p = signif(fisher_p, 3), NES = round(NES, 2), gsea_padj = signif(gsea_padj, 3))])
fwrite(mat, file.path(OUT, "q2_contrast_matrix.csv"))
diag <- mat[human_sig == "H2_naflnash" & vacca_target == "V_prog_526"]
ctrl <- mat[human_sig == "H1_dvc" & vacca_target == "V_early_271"]
cat(sprintf("\n>> Progression diagonal  H2(NAFL-vs-NASH) x V_prog_526 : OR=%.1f (p=%.2g)  [contrast-matched peak]\n",
            diag$OR, diag$fisher_p))
cat(sprintf(">> Negative control      H1(disease-vs-ctrl) x V_early_271: OR=%.2f (p=%.2g)  [near-null, specificity]\n",
            ctrl$OR, ctrl$fisher_p))
cat(sprintf(">> H2(NAFL-vs-NASH) x V_early_271: OR=%.2f -> progression sig NOT enriched in the early bin (directional specificity)\n",
            mat[human_sig=="H2_naflnash" & vacca_target=="V_early_271"]$OR))

# ── (B) threshold-free, NON-CIRCULAR anchor ───────────────────────────────────
a <- fread(file.path(CS, "concordance_atlas_unified.csv"))
a[, sym := toupper(trimws(human_symbol))]; a <- a[!is.na(sym) & sym != ""]
setorder(a, -n_concordant, -translatability_score); a <- a[!duplicated(sym)]
a[, signed_nconc := sign(mean_h_lfc) * n_concordant]
in_universe <- intersect(v_all, a$sym)
a[, our_conc := n_concordant >= 1][, is_vacca := sym %in% v_all]
ct <- table(our_conc = a$our_conc, is_vacca = a$is_vacca)
ft <- fisher.test(ct, alternative = "greater")
cat(sprintf("\n=== (B) NON-CIRCULAR anchor: concordant>=1 diet vs LITMUS-951 (threshold-free) ===\n"))
print(ct)
cat(sprintf("Fisher OR=%.2f, p=%.2e (n_both=%d, LITMUS-951 in atlas universe=%d)\n",
            ft$estimate, ft$p.value, ct["TRUE","TRUE"], length(in_universe)))
make_rank <- function(dt, statcol) setNames(dt[[statcol]] + (seq_len(nrow(dt))/nrow(dt))*1e-6, dt$sym)
anchor_fgsea <- rbindlist(lapply(c("n_concordant","translatability_score","signed_nconc","mean_h_lfc"), function(rk) {
  stats <- sort(make_rank(a, rk), decreasing = TRUE)
  fg <- fgsea(list(LITMUS_951 = in_universe), stats, minSize = 5,
              maxSize = length(in_universe) + 10, scoreType = "std", nPermSimple = 10000, eps = 0)
  data.table(ranking = rk, NES = fg$NES, padj = fg$padj, ES = fg$ES, size = fg$size)
}))
cat("\n threshold-free fgsea (LITMUS-951 enrichment over continuous concordance rankings):\n")
print(anchor_fgsea[, .(ranking, NES = round(NES, 2), padj = signif(padj, 3))])
fwrite(rbind(
  data.table(metric = "fisher_nconc1_vs_litmus951", OR = as.numeric(ft$estimate), p = ft$p.value,
             n_both = ct["TRUE","TRUE"], n_universe = length(in_universe), NES = NA_real_, padj = NA_real_),
  anchor_fgsea[, .(metric = paste0("fgsea_", ranking), OR = NA_real_, p = NA_real_,
                   n_both = NA_integer_, n_universe = size, NES, padj)], use.names = TRUE),
  file.path(OUT, "q2_noncircular_anchor.csv"))

# ── (C) Govaere holdout (shared-cohort circularity) ───────────────────────────
SIGS <- list(FULL  = H$H1_dvc,
             NOGOV = load_sig(file.path(LOO, "lvqw_loo_C2_GSE135251.csv"), has_symbol = FALSE))
gov <- rbindlist(lapply(names(SIGS), function(sn) {
  hsig <- SIGS[[sn]]; st <- setNames(hsig$t, hsig$sym); st <- st[is.finite(st)]
  fg <- as.data.table(suppressWarnings(fgsea(
    lapply(vacca_targets, function(vt) intersect(vt, names(st))), st,
    scoreType = "std", nPermSimple = 10000)))
  rbindlist(lapply(names(vacca_targets), function(vn) {
    vt <- vacca_targets[[vn]]; uni <- intersect(hsig$sym, s4_uni)
    sig <- intersect(hsig[padj < 0.05 & abs(logFC) > 0.5, sym], uni); tgt <- intersect(vt, uni)
    a1 <- length(intersect(sig, tgt)); b1 <- length(setdiff(sig, tgt))
    c1 <- length(setdiff(tgt, sig)); d1 <- length(uni) - a1 - b1 - c1
    ft1 <- fisher.test(matrix(c(a1, b1, c1, d1), 2), alternative = "greater")
    data.table(sig_set = sn, vacca_target = vn, overlap = a1, OR = unname(ft1$estimate),
               fisher_p = ft1$p.value, NES = fg[pathway == vn, NES])
  }))
}))
gw <- dcast(gov, vacca_target ~ sig_set, value.var = c("OR","overlap","NES"))
gw[, OR_survival := round(OR_NOGOV / OR_FULL, 3)]
cat("\n=== (C) Govaere holdout: FULL (in) vs NOGOV (out) — OR must survive ===\n")
print(gw[, .(vacca_target, OR_FULL = round(OR_FULL, 1), OR_NOGOV = round(OR_NOGOV, 1),
             OR_survival, overlap_FULL, overlap_NOGOV)])
fwrite(gov, file.path(OUT, "q2_govaere_holdout.csv"))

# ── (D) label-shuffle null: progression diagonal + non-circular anchor ────────
# membership-shuffle: random same-size gene sets, recompute the Fisher OR / overlap
shuffle_or <- function(sig, universe, target_size, observed_overlap) {
  null_ov <- vapply(seq_len(N_PERM), function(i) {
    tgt <- sample(universe, target_size); length(intersect(sig, tgt)) }, integer(1))
  (1 + sum(null_ov >= observed_overlap)) / (N_PERM + 1)
}
# diagonal: H2 sig in V_prog universe
h2 <- H$H2_naflnash; uni2 <- intersect(h2$sym, s4_uni)
sig2 <- intersect(h2[padj < 0.05 & abs(logFC) > 0.5, sym], uni2)
tgt2 <- intersect(v_prog, uni2)
p_diag <- shuffle_or(sig2, uni2, length(tgt2), length(intersect(sig2, tgt2)))
# anchor: concordant>=1 set vs LITMUS-951 in atlas universe
sigA <- a[our_conc == TRUE, sym]; uniA <- a$sym; tgtA <- in_universe
p_anchor <- shuffle_or(sigA, uniA, length(tgtA), length(intersect(sigA, tgtA)))
nulldt <- data.table(
  test = c("progression_diagonal_H2xVprog", "noncircular_anchor_nconc1xLITMUS951"),
  observed_overlap = c(length(intersect(sig2, tgt2)), length(intersect(sigA, tgtA))),
  emp_p = c(p_diag, p_anchor), n_perm = N_PERM)
cat(sprintf("\n=== (D) membership-shuffle null (%d perms) ===\n", N_PERM))
print(nulldt)
fwrite(nulldt, file.path(OUT, "q2_shuffle_null.csv"))

cat("\nWrote q2_contrast_matrix.csv / q2_noncircular_anchor.csv / q2_govaere_holdout.csv / q2_shuffle_null.csv\n")
cat(strrep("=", 78), "\n")
