#!/usr/bin/env Rscript
# ============================================================================
# q2_contrast_match.R  (Q2 follow-up to Vacca 2024 benchmark)
#
# QUESTION: Vacca's 951-gene signature is itself stage-stratified — its largest
# bin (526 genes) is a PROGRESSION (Moderate/Severe-vs-Mild) signature, while our
# Conserved_Core's prior overlap (172) was computed against our DISEASE-vs-CONTROL
# anchor. Is the modest 172 partly a CONTRAST MISMATCH (disease-vs-ctrl vs
# progression), not a true biological ceiling?
#
# DESIGN: 3 HUMAN signatures × 3 VACCA targets, each evaluated two ways:
#   Human signatures (ranked by their own t/logFC, padj<0.05 & |logFC|>0.5 = sig set):
#     H1  disease-vs-control     canonical_deg_results.csv
#     H2  nafl_vs_nash (LVQW)    nafl_vs_nash_lvqw.csv          [progression contrast]
#     H3  adv-vs-early fibrosis  adv_vs_early_fibrosis_dream.csv [progression contrast]
#   Vacca targets (human symbols, the 951 sig partitioned by their stage column):
#     V_all   951  (full DSEA signature)
#     V_prog  526  ("Disease progression")          <- apples-to-apples w/ H2/H3
#     V_early 271  ("Early disease development")     <- apples-to-apples w/ H1
#     (the remaining 154 = "All disease stages")
#
#   (a) GENE-LEVEL OVERLAP + Fisher OR: universe = genes testable in our human
#       signature AND present in Vacca S4 (so each Vacca target is compared on a
#       fair common universe). a=overlap of OUR-sig ∩ VACCA-target.
#   (b) GSEA: rank OUR human genes by t-stat; Vacca target (as human-symbol set)
#       as the gene set; report NES + padj (directionless membership enrichment).
#
# EXPECTATION (to test): the progression human signatures (H2/H3) should match
# Vacca's progression bin (V_prog) better than disease-vs-ctrl matches it, and
# better than H1 matches V_prog — i.e. apples-to-apples lifts the overlap/OR/NES.
# Env: rnaseq
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
set.seed(42)
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
DSIG  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
RES   <- file.path(BASE, "results/vacca_benchmark"); dir.create(RES, recursive=TRUE, showWarnings=FALSE)
strip <- function(x) gsub("\\..*", "", x)

# ── Vacca S4: 951 signature partitioned by stage (human symbols) ─────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"), sheet="Table S4"))
setnames(s4, c("GeneSymbol","Gene_used_in_DSEA(0:No/1:Yes)",
               "Early/All/Late disease stage (if used in DSEA otherwise NA)"),
         c("sym","in_dsea","stage"), skip_absent=TRUE)
s4[, sym := toupper(trimws(sym))]
s4_uni  <- unique(s4[!is.na(sym) & sym != "", sym])                 # Vacca-testable universe
v_all   <- unique(s4[in_dsea==1, sym])
v_prog  <- unique(s4[in_dsea==1 & stage=="Disease progression", sym])
v_early <- unique(s4[in_dsea==1 & stage=="Early disease development", sym])
cat(sprintf("Vacca S4 universe=%d | 951=%d | progression=%d | early=%d\n",
            length(s4_uni), length(v_all), length(v_prog), length(v_early)))
vacca_targets <- list(V_all_951=v_all, V_prog_526=v_prog, V_early_271=v_early)

# ── Load 3 human signatures, normalize to (symbol, logFC, t, padj) ───────────
load_sig <- function(path, lfc_col, t_col, padj_col, sym_col=NULL) {
  d <- fread(path)
  if (is.null(sym_col)) {                      # need to add symbol from canonical map
    sym_col <- "symbol"
  }
  if (!"symbol" %in% names(d)) {
    # map ENSEMBL -> symbol via canonical_deg_results.csv
    map <- fread(file.path(INT_I, "canonical_deg_results.csv"))[, .(gene=strip(gene), symbol)]
    d[, gene := strip(gene)]
    d <- merge(d, map, by="gene", all.x=TRUE)
  }
  out <- d[, .(sym=toupper(trimws(get(sym_col))),
               logFC=as.numeric(get(lfc_col)),
               t=as.numeric(get(t_col)),
               padj=as.numeric(get(padj_col)))]
  out <- out[!is.na(sym) & sym!="" & is.finite(t)]
  out[!duplicated(sym)]
}
H <- list(
  H1_dvc       = load_sig(file.path(INT_I,"canonical_deg_results.csv"),"logFC","t","padj","symbol"),
  H2_naflnash  = load_sig(file.path(DSIG,"nafl_vs_nash_lvqw.csv"),"logFC","t","padj"),  # no symbol col -> map
  H3_advfib    = load_sig(file.path(DSIG,"adv_vs_early_fibrosis_dream.csv"),"logFC","t","adj.P.Val","symbol")
)
for (nm in names(H)) cat(sprintf("%-12s n_tested=%d  n_sig(padj<.05,|lfc|>.5)=%d\n",
   nm, nrow(H[[nm]]), H[[nm]][padj<0.05 & abs(logFC)>0.5, .N]))

# ── (a) overlap + Fisher OR, common universe per (human sig × vacca target) ──
fisher_block <- rbindlist(lapply(names(H), function(hn) {
  hsig <- H[[hn]]
  rbindlist(lapply(names(vacca_targets), function(vn) {
    vt <- vacca_targets[[vn]]
    uni <- intersect(hsig$sym, s4_uni)                 # genes testable in OUR sig AND in Vacca S4
    sig <- intersect(hsig[padj<0.05 & abs(logFC)>0.5, sym], uni)
    tgt <- intersect(vt, uni)
    a <- length(intersect(sig, tgt)); b <- length(setdiff(sig, tgt))
    c <- length(setdiff(tgt, sig));  d <- length(uni) - a - b - c
    ft <- fisher.test(matrix(c(a,b,c,d), 2), alternative="greater")
    data.table(human_sig=hn, vacca_target=vn, universe=length(uni),
               n_our_sig=length(sig), n_vacca=length(tgt), overlap=a,
               OR=unname(ft$estimate), fisher_p=ft$p.value,
               jaccard=a/length(union(sig,tgt)))
  }))
}))
cat("\n=== (a) GENE-LEVEL OVERLAP + Fisher OR (one-sided greater) ===\n")
print(fisher_block[, .(human_sig, vacca_target, universe, n_our_sig, n_vacca, overlap,
                       OR=round(OR,2), fisher_p=signif(fisher_p,3), jaccard=round(jaccard,3))])

# ── (b) fgsea: rank OUR genes by t; Vacca target = gene set ──────────────────
gsea_block <- rbindlist(lapply(names(H), function(hn) {
  hsig <- H[[hn]]; st <- setNames(hsig$t, hsig$sym); st <- st[is.finite(st)]
  sets <- lapply(vacca_targets, function(vt) intersect(vt, names(st)))
  fg <- as.data.table(suppressWarnings(fgsea(sets, st, scoreType="std", nPermSimple=10000)))
  fg[, human_sig := hn][, .(human_sig, vacca_target=pathway, set_size=size,
                            NES=round(NES,3), gsea_padj=signif(padj,3),
                            leadingEdge_n=lengths(leadingEdge))]
}))
cat("\n=== (b) GSEA: our t-ranked genes, Vacca target as gene set ===\n")
print(gsea_block)

# ── merge + verdict ──────────────────────────────────────────────────────────
res <- merge(fisher_block, gsea_block, by=c("human_sig","vacca_target"))
fwrite(res, file.path(RES, "q2_contrast_match.csv"))
fwrite(res, file.path(OUT, "q2_contrast_match.csv"))

cat("\n=== VERDICT: which human signature best matches Vacca progression (V_prog_526)? ===\n")
print(res[vacca_target=="V_prog_526"][order(-OR),
      .(human_sig, overlap, OR=round(OR,2), fisher_p=signif(fisher_p,2),
        NES=round(NES,2), gsea_padj=signif(gsea_padj,2))])
cat("\n=== best human match per Vacca target (by OR) ===\n")
print(res[, .SD[which.max(OR)], by=vacca_target,
      .SDcols=c("human_sig","overlap","OR","fisher_p","NES","gsea_padj")][
      , .(vacca_target, human_sig, overlap, OR=round(OR,2),
          fisher_p=signif(fisher_p,2), NES=round(NES,2), gsea_padj=signif(gsea_padj,2))])
cat("\nWrote:", file.path(RES,"q2_contrast_match.csv"), "\n")
