#!/usr/bin/env Rscript
# 206_published_vs_perstudy.R
# Apple-to-apple comparison of published DEG results from Govaere 2020
# (GSE135251) and Hoang 2019 (GSE130970) against our pipeline results.
#
# For each published analysis, we reproduce the SAME contrast on the SAME
# data using limma-voom, then also compare to dream mega-analysis.
#
# Hoang: NAS ordinal + fibrosis ordinal regression → our NAS/fib ordinal +
# dream fibrosis ordinal
# Govaere: NASH F2-F4 vs NAFL, NASH F3-F4 vs F0/1, NAS>=4, SAF>=2 →
# our reproduced pairwise contrasts + 25-gene signature

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
INTB <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PUB  <- file.path(BASE, "data/published_degs")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 206: Published vs Per-Study DEG Comparison ===\n")

# ── Step 1: Gene mapping ────────────────────────────────────────────────────
cat("Loading gene annotation...\n")
annot <- fread(file.path(INT, "gene_annotation/human_ensg_to_symbol.tsv"))
annot_pc <- annot[gene_type == "protein_coding"]
annot_other <- annot[!gene_base %in% annot_pc$gene_base]
sym2ens <- setNames(annot_pc$gene_base, annot_pc$symbol)
extra <- annot_other[!symbol %in% names(sym2ens)]
extra_dedup <- extra[!duplicated(symbol)]
sym2ens <- c(sym2ens, setNames(extra_dedup$gene_base, extra_dedup$symbol))
sym2ens_upper <- setNames(sym2ens, toupper(names(sym2ens)))
cat(sprintf("  %d symbol -> Ensembl mappings\n", length(sym2ens)))

map_symbols <- function(symbols) {
  symbols <- as.character(symbols)
  mapped <- sym2ens[symbols]
  missing <- is.na(mapped)
  if (any(missing)) mapped[missing] <- sym2ens_upper[toupper(symbols[missing])]
  mapped
}

# ── Step 2: Load raw data for custom contrasts ──────────────────────────────
cat("Loading raw counts and metadata...\n")
counts <- readRDS(file.path(INTB, "results/integration/merged_counts_raw.rds"))
sample_meta <- readRDS(file.path(INTB, "results/integration/meta_matched.rds"))
qc <- fread(file.path(INTB, "qc/sample_qc_report.csv"))
sample_meta <- sample_meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Helper: run limma-voom for a subset of samples with a given design/contrast
run_limma <- function(sample_ids, design_mat, coef_name) {
  idx <- colnames(counts) %in% sample_ids
  dge <- DGEList(counts = counts[, idx])
  dge <- calcNormFactors(dge, method = "TMM")
  keep <- filterByExpr(dge, design = design_mat)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  v <- voom(dge, design_mat, plot = FALSE)
  fit <- lmFit(v, design_mat)
  if (is.character(coef_name) && !coef_name %in% colnames(fit$coefficients)) {
    # Might be a contrast
    fit <- contrasts.fit(fit, coef_name)
    fit <- eBayes(fit)
    tt <- topTable(fit, number = Inf, sort.by = "none")
  } else {
    fit <- eBayes(fit)
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  }
  tt$gene <- rownames(tt)
  dt <- as.data.table(tt)
  dt[, gene_base := sub("\\.\\d+$", "", gene)]
  dt
}

# ── Step 3: Load published DEGs ─────────────────────────────────────────────
cat("\nLoading Hoang et al. 2019 published DEGs...\n")
hoang_xlsx <- file.path(PUB, "GSE130970/MOESM2.xlsx")

std_hoang <- function(dt) {
  nms <- tolower(gsub("\\s+", "_", names(dt)))
  setnames(dt, names(dt), nms)
  if ("gene_symbol" %in% nms) setnames(dt, "gene_symbol", "symbol", skip_absent = TRUE)
  if ("adj_p" %in% nms) setnames(dt, "adj_p", "adj_p_val", skip_absent = TRUE)
  dt[, gene_base := map_symbols(symbol)]
  dt[!is.na(gene_base)]
}

hoang_nas <- std_hoang(as.data.table(read_excel(hoang_xlsx, sheet = "NAS ordinal regression")))
hoang_fib <- std_hoang(as.data.table(read_excel(hoang_xlsx, sheet = "fibrosis ordinal regression")))
cat(sprintf("  NAS: %d mapped, %d DEGs | Fib: %d mapped, %d DEGs\n",
            nrow(hoang_nas), sum(hoang_nas$adj_p_val < 0.01),
            nrow(hoang_fib), sum(hoang_fib$adj_p_val < 0.01)))

cat("\nLoading Govaere et al. 2020 supplementary tables...\n")
govaere_25 <- data.table(
  symbol = c("AKR1B10", "ANKRD29", "CCL20", "CFAP221", "CLIC6", "COL1A1",
             "COL1A2", "DTNA", "DUSP8", "EPB41L4A", "FERMT1", "GDF15",
             "HECW1", "IL32", "ITGBL1", "LTBP2", "PDGFA", "PPAPDC1A",
             "RGS4", "SCTR", "STMN2", "THY1", "TNFRSF12A", "TYMS",
             "HSD17B14"),
  published_direction = c(rep("up", 24), "down")
)
govaere_25[, gene_base := map_symbols(symbol)]

govaere_supp <- list()
govaere_xlsx <- file.path(PUB, "GSE135251/aba4448_supplementary_tables.xlsx")
if (file.exists(govaere_xlsx)) {
  supp_info <- data.table(
    sheet = paste("Table", c("S1", "S3", "S4", "S5", "S6", "S7", "S8", "S9")),
    contrast = c("Cluster A vs B", "NASH_F2_vs_NAFL", "NASH_F3_vs_NAFL",
                 "NASH_F4_vs_NAFL", "NASH_F3_vs_NASH_F01",
                 "NASH_F4_vs_NASH_F01", "NAS_ge4", "SAF_ge2"),
    expected_n = c(250, 50, 907, 1369, 434, 1194, 369, 320)
  )
  avail_sheets <- excel_sheets(govaere_xlsx)
  for (i in seq_len(nrow(supp_info))) {
    sn <- supp_info$sheet[i]
    if (!sn %in% avail_sheets) next
    tryCatch({
      raw <- as.data.table(read_excel(govaere_xlsx, sheet = sn, skip = 1))
      orig_names <- names(raw)
      setnames(raw, orig_names, tolower(gsub("\\s+", "_", orig_names)))
      nms <- names(raw)
      sym_col <- grep("gene_name|external_gene_name", nms, value = TRUE)[1]
      ens_col <- grep("ensembl|gene_id", nms, value = TRUE)[1]
      lfc_col <- grep("log2?fc|logfc", nms, value = TRUE)[1]
      fdr_col <- grep("qvalue|q_value|fdr|adj", nms, value = TRUE)[1]
      if (!is.na(ens_col)) raw[, gene_base := sub("\\.\\d+$", "", as.character(raw[[ens_col]]))]
      if (!is.na(sym_col)) raw[, symbol := as.character(raw[[sym_col]])]
      if (is.na(ens_col) && !is.na(sym_col)) raw[, gene_base := map_symbols(symbol)]
      if (!is.na(lfc_col)) raw[, pub_logFC := as.numeric(raw[[lfc_col]])]
      if (!is.na(fdr_col)) raw[, pub_fdr := as.numeric(raw[[fdr_col]])]
      govaere_supp[[supp_info$contrast[i]]] <- list(
        data = raw[!is.na(gene_base)], contrast = supp_info$contrast[i],
        expected_n = supp_info$expected_n[i],
        has_lfc = !is.na(lfc_col), has_fdr = !is.na(fdr_col))
      cat(sprintf("  %s: %d genes (expected %d)\n", supp_info$contrast[i],
                  sum(!is.na(raw$gene_base)), supp_info$expected_n[i]))
    }, error = function(e) cat(sprintf("  %s: FAILED (%s)\n", sn, conditionMessage(e))))
  }
}

# ── Step 4: Compute apple-to-apple contrasts ────────────────────────────────
cat("\n=== Computing apple-to-apple contrasts ===\n")

# --- Hoang: Load CLM results (same method as Hoang, from 206b) ---
cat("Hoang: Loading CLM results (same method as published)...\n")
clm_nas <- fread(file.path(OUTDIR, "clm_nas_results.csv"))
clm_fib <- fread(file.path(OUTDIR, "clm_fib_results.csv"))

# --- Dream fibrosis ordinal ---
dream_fib_ord <- fread(file.path(INT,
  "../results/progression/c17_fibrosis_ordinal_dream.csv"))
dream_fib_ord[, gene_base := sub("\\.\\d+$", "", gene)]

# --- Dream disease vs control ---
dream <- fread(file.path(INT, "integration/dream_results.csv"))
dream[, gene_base := sub("\\.\\d+$", "", gene)]

# --- Govaere: Reproduce each specific contrast for GSE135251 ---
cat("\nGovaere: Computing reproduced contrasts for GSE135251...\n")
m135 <- sample_meta[dataset == "GSE135251" & !is.na(fibrosis_stage)]

# Build fine-grained groups matching Govaere's paper
m135[, govaere_group := fcase(
  condition == "Control", "Control",
  condition == "NAFL", "NAFL",
  condition %in% c("NASH") & fibrosis_stage %in% 0:1, "NASH_F01",
  condition %in% c("NASH") & fibrosis_stage == 2, "NASH_F2",
  condition %in% c("NASH_Fibrosis") & fibrosis_stage == 3, "NASH_F3",
  condition %in% c("NASH_Fibrosis") & fibrosis_stage == 4, "NASH_F4",
  default = NA_character_
)]
cat("  Group sizes:\n")
print(m135[, .N, by = govaere_group][order(govaere_group)])

# Define contrasts matching Govaere's supplementary tables
govaere_contrasts <- list(
  NASH_F2_vs_NAFL = list(group1 = "NASH_F2", group2 = "NAFL"),
  NASH_F3_vs_NAFL = list(group1 = "NASH_F3", group2 = "NAFL"),
  NASH_F4_vs_NAFL = list(group1 = "NASH_F4", group2 = "NAFL"),
  NASH_F3_vs_NASH_F01 = list(group1 = "NASH_F3", group2 = "NASH_F01"),
  NASH_F4_vs_NASH_F01 = list(group1 = "NASH_F4", group2 = "NASH_F01"),
  NAS_ge4 = "NAS_binary"  # handled separately
)

our_govaere <- list()
for (cname in names(govaere_contrasts)) {
  info <- govaere_contrasts[[cname]]

  if (identical(info, "NAS_binary")) {
    # NAS >= 4 vs NAS < 4
    m_sub <- m135[!is.na(nas_score)]
    m_sub[, nas_ge4 := factor(ifelse(nas_score >= 4, "high", "low"), levels = c("low", "high"))]
    design_c <- model.matrix(~ nas_ge4 + inferred_sex, data = m_sub)
    our_govaere[[cname]] <- run_limma(m_sub$sample_id, design_c, "nas_ge4high")
  } else {
    m_sub <- m135[govaere_group %in% c(info$group1, info$group2)]
    m_sub[, grp := factor(govaere_group, levels = c(info$group2, info$group1))]
    if (nrow(m_sub) < 5) { cat(sprintf("  %s: too few samples, skipping\n", cname)); next }
    design_c <- model.matrix(~ grp + inferred_sex, data = m_sub)
    coef_nm <- paste0("grp", info$group1)
    our_govaere[[cname]] <- run_limma(m_sub$sample_id, design_c, coef_nm)
  }

  n_sig <- sum(our_govaere[[cname]]$adj.P.Val < 0.05)
  cat(sprintf("  %s: %d genes, %d DEGs (padj<0.05)\n", cname,
              nrow(our_govaere[[cname]]), n_sig))
}

# ── Step 5: Comparison functions ─────────────────────────────────────────────

lfc_correlation <- function(dt1, dt2, lfc1 = "logFC", lfc2 = "logFC",
                            id_col = "gene_base") {
  merged <- merge(dt1[, c(id_col, lfc1), with = FALSE],
                  dt2[, c(id_col, lfc2), with = FALSE], by = id_col)
  setnames(merged, c(lfc1, lfc2), c("lfc1", "lfc2"))
  merged <- merged[is.finite(lfc1) & is.finite(lfc2)]
  n <- nrow(merged)
  if (n < 10) return(list(n = n, pearson = NA, spearman = NA, dir_conc = NA, data = merged))
  list(n = n,
       pearson = cor(merged$lfc1, merged$lfc2, method = "pearson"),
       spearman = cor(merged$lfc1, merged$lfc2, method = "spearman"),
       dir_conc = mean(sign(merged$lfc1) == sign(merged$lfc2)) * 100,
       data = merged)
}

deg_overlap <- function(set1, set2, universe_n) {
  overlap <- length(intersect(set1, set2))
  union_n <- length(union(set1, set2))
  jaccard <- if (union_n > 0) overlap / union_n else 0
  a <- overlap; b <- length(set1) - overlap
  c <- length(set2) - overlap; d <- universe_n - length(union(set1, set2))
  ft <- fisher.test(matrix(c(a, b, c, max(d, 0)), nrow = 2))
  list(n1 = length(set1), n2 = length(set2), overlap = overlap,
       jaccard = jaccard, fisher_OR = ft$estimate, fisher_p = ft$p.value)
}

recovery_rate <- function(ref_set, query_list) {
  sapply(query_list, function(q) mean(ref_set %in% q) * 100)
}

# Scatter panel helper
scatter_panel <- function(merged, xlab, ylab, title, anno_text) {
  ggplot(merged, aes(x = lfc1, y = lfc2)) +
    geom_point(alpha = 0.15, size = 0.3, color = masld_colors$ns) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = masld_colors$up) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted",
                color = "grey40") +
    annotate("text", x = -Inf, y = Inf, label = anno_text,
             hjust = -0.1, vjust = 1.3, size = 2, color = "grey30") +
    labs(x = xlab, y = ylab, title = title) +
    theme_masld()
}

make_anno <- function(res) {
  sprintf("r = %.3f\nrho = %.3f\ndir = %.1f%%\nn = %s",
          res$pearson, res$spearman, res$dir_conc, formatC(res$n, big.mark = ","))
}

# ── Step 6: Hoang comparisons ───────────────────────────────────────────────
cat("\n=== Hoang Comparisons ===\n")
metrics <- list()
PADJ <- 0.05; DREAM_PADJ <- 0.1
hoang_nas_degs <- hoang_nas[adj_p_val < 0.01, gene_base]
hoang_fib_degs <- hoang_fib[adj_p_val < 0.01, gene_base]
clm_nas_degs <- clm_nas[!is.na(adj_P) & adj_P < 0.01, gene_base]
clm_fib_degs <- clm_fib[!is.na(adj_P) & adj_P < 0.01, gene_base]
dream_degs <- dream[padj < DREAM_PADJ, gene_base]
dream_fib_degs <- dream_fib_ord[padj < DREAM_PADJ, gene_base]
universe_n <- length(union(clm_nas$gene_base, hoang_nas$gene_base))

# H1: NAS ordinal — same method (CLM), same metric (range_log2FC)
cat("H1: Hoang NAS ordinal vs Our CLM NAS ordinal...\n")
h1 <- lfc_correlation(hoang_nas, clm_nas, "range_log2fc", "range_log2FC", "gene_base")
h1_ov <- deg_overlap(hoang_nas_degs, clm_nas_degs, universe_n)
cat(sprintf("  r=%.3f, rho=%.3f, dir=%.1f%%, overlap=%d, J=%.3f\n",
            h1$pearson, h1$spearman, h1$dir_conc, h1_ov$overlap, h1_ov$jaccard))
metrics[["H1"]] <- list(study="Hoang", comparison="NAS_ordinal_CLM_per_study",
  pub_method="NAS ordinal regression (CLM)", our_method="CLM NAS ordinal per-study",
  pub_n=h1_ov$n1, our_n=h1_ov$n2, overlap=h1_ov$overlap, jaccard=h1_ov$jaccard,
  fisher_OR=h1_ov$fisher_OR, fisher_p=h1_ov$fisher_p,
  pearson=h1$pearson, spearman=h1$spearman, dir_conc=h1$dir_conc)

# H2: Fibrosis ordinal — same method (CLM), same metric (range_log2FC)
cat("H2: Hoang Fib ordinal vs Our CLM Fib ordinal...\n")
h2 <- lfc_correlation(hoang_fib, clm_fib, "range_log2fc", "range_log2FC", "gene_base")
h2_ov <- deg_overlap(hoang_fib_degs, clm_fib_degs, universe_n)
cat(sprintf("  r=%.3f, rho=%.3f, dir=%.1f%%, overlap=%d, J=%.3f\n",
            h2$pearson, h2$spearman, h2$dir_conc, h2_ov$overlap, h2_ov$jaccard))
metrics[["H2"]] <- list(study="Hoang", comparison="Fib_ordinal_CLM_per_study",
  pub_method="Fibrosis ordinal regression (CLM)", our_method="CLM fibrosis ordinal per-study",
  pub_n=h2_ov$n1, our_n=h2_ov$n2, overlap=h2_ov$overlap, jaccard=h2_ov$jaccard,
  fisher_OR=h2_ov$fisher_OR, fisher_p=h2_ov$fisher_p,
  pearson=h2$pearson, spearman=h2$spearman, dir_conc=h2$dir_conc)

# H3: NAS ordinal vs dream disease-vs-control
cat("H3: Hoang NAS ordinal vs Dream (disease vs control)...\n")
h3 <- lfc_correlation(hoang_nas, dream, "range_log2fc", "logFC", "gene_base")
h3_ov <- deg_overlap(hoang_nas_degs, dream_degs, universe_n)
cat(sprintf("  r=%.3f, rho=%.3f, dir=%.1f%%\n", h3$pearson, h3$spearman, h3$dir_conc))
metrics[["H3"]] <- list(study="Hoang", comparison="NAS_ordinal_vs_dream",
  pub_method="NAS ordinal regression", our_method="Dream mega-analysis",
  pub_n=h3_ov$n1, our_n=h3_ov$n2, overlap=h3_ov$overlap, jaccard=h3_ov$jaccard,
  fisher_OR=h3_ov$fisher_OR, fisher_p=h3_ov$fisher_p,
  pearson=h3$pearson, spearman=h3$spearman, dir_conc=h3$dir_conc)

# H4: Fib ordinal vs dream fibrosis ordinal
cat("H4: Hoang Fib ordinal vs Dream fibrosis ordinal...\n")
h4 <- lfc_correlation(hoang_fib, dream_fib_ord, "range_log2fc", "logFC", "gene_base")
h4_ov <- deg_overlap(hoang_fib_degs, dream_fib_degs, universe_n)
cat(sprintf("  r=%.3f, rho=%.3f, dir=%.1f%%\n", h4$pearson, h4$spearman, h4$dir_conc))
metrics[["H4"]] <- list(study="Hoang", comparison="Fib_ordinal_vs_dream_fib_ordinal",
  pub_method="Fibrosis ordinal regression", our_method="Dream fibrosis ordinal (8 cohorts)",
  pub_n=h4_ov$n1, our_n=h4_ov$n2, overlap=h4_ov$overlap, jaccard=h4_ov$jaccard,
  fisher_OR=h4_ov$fisher_OR, fisher_p=h4_ov$fisher_p,
  pearson=h4$pearson, spearman=h4$spearman, dir_conc=h4$dir_conc)

# H5: Recovery rates
cat("H5: Recovery rates...\n")
rec_nas <- recovery_rate(hoang_nas_degs, list(
  "Per-study" = clm_nas_degs,
  "Dream" = dream_degs
))
rec_fib <- recovery_rate(hoang_fib_degs, list(
  "Per-study" = clm_fib_degs,
  "Dream" = dream_fib_degs
))
cat(sprintf("  NAS: per-study=%.1f%%, dream=%.1f%%\n", rec_nas[1], rec_nas[2]))
cat(sprintf("  Fib: per-study=%.1f%%, dream=%.1f%%\n", rec_fib[1], rec_fib[2]))

# ── Step 7: Govaere comparisons ─────────────────────────────────────────────
cat("\n=== Govaere Comparisons ===\n")

# G1: 25-gene signature
cat("G1: 25-gene signature...\n")
ps_135251 <- fread(file.path(INT, "per_study/GSE135251_de_results.csv"))
ps_135251[, gene_base := sub("\\.\\d+$", "", gene)]

govaere_25_mapped <- govaere_25[!is.na(gene_base)]
g1_merge <- merge(govaere_25_mapped,
                  ps_135251[, .(gene_base, logFC, adj.P.Val)],
                  by = "gene_base", all.x = TRUE)
g1_merge <- merge(g1_merge,
                  dream[, .(gene_base, dream_logFC = logFC, dream_padj = padj)],  # C2-OK-sensitivity: dream_results.csv retired sensitivity arm, labeled "Dream" comparator
                  by = "gene_base", all.x = TRUE)
g1_merge[, our_direction := ifelse(logFC > 0, "up", "down")]
g1_merge[, direction_match := published_direction == our_direction]
cat(sprintf("  Direction: %d/%d (%.0f%%), sig per-study: %d, sig dream: %d\n",
            sum(g1_merge$direction_match, na.rm = TRUE), nrow(g1_merge),
            100 * mean(g1_merge$direction_match, na.rm = TRUE),
            sum(g1_merge$adj.P.Val < PADJ, na.rm = TRUE),
            sum(g1_merge$dream_padj < DREAM_PADJ, na.rm = TRUE)))  # C2-OK-sensitivity: dream arm of published comparison

# G2-G7: Pairwise contrast comparisons (published vs our reproduced)
govaere_compare_names <- intersect(names(govaere_supp), names(our_govaere))
for (cname in govaere_compare_names) {
  pub <- govaere_supp[[cname]]
  our <- our_govaere[[cname]]
  cat(sprintf("G-%s: published vs our reproduced...\n", cname))

  pub_degs <- pub$data$gene_base
  our_degs <- our[adj.P.Val < PADJ, gene_base]
  ov <- deg_overlap(pub_degs, our_degs, length(unique(our$gene_base)))

  lfc_res <- list(pearson = NA, spearman = NA, dir_conc = NA, data = NULL)
  if (pub$has_lfc) {
    lfc_res <- lfc_correlation(pub$data, our, "pub_logFC", "logFC", "gene_base")
  }
  cat(sprintf("  pub=%d, ours=%d, overlap=%d, J=%.3f",
              length(pub_degs), length(our_degs), ov$overlap, ov$jaccard))
  if (!is.na(lfc_res$pearson)) {
    cat(sprintf(", r=%.3f, rho=%.3f, dir=%.1f%%",
                lfc_res$pearson, lfc_res$spearman, lfc_res$dir_conc))
  }
  cat("\n")

  metrics[[paste0("G_", cname)]] <- list(
    study = "Govaere", comparison = paste0(cname, "_reproduced"),
    pub_method = "limma-voom (q<0.05, FC>1.5)", our_method = "limma-voom reproduced",
    pub_n = length(pub_degs), our_n = length(our_degs),
    overlap = ov$overlap, jaccard = ov$jaccard,
    fisher_OR = ov$fisher_OR, fisher_p = ov$fisher_p,
    pearson = lfc_res$pearson, spearman = lfc_res$spearman, dir_conc = lfc_res$dir_conc)
}

# ── Step 8: Figure assembly ─────────────────────────────────────────────────
cat("\nGenerating figures...\n")

# --- Page 1: Hoang (5 panels) ---
pa <- scatter_panel(h1$data,
  "Hoang NAS ordinal (range log2FC)", "Our NAS ordinal CLM (range log2FC)",
  "NAS ordinal per-study", make_anno(h1))

pb <- scatter_panel(h2$data,
  "Hoang Fib ordinal (range log2FC)", "Our Fib ordinal CLM (range log2FC)",
  "Fib ordinal per-study", make_anno(h2))

pc <- scatter_panel(h3$data,
  "Hoang NAS ordinal (range log2FC)", "Dream disease vs ctrl (logFC)",
  "NAS ordinal vs Dream", make_anno(h3))

pd <- scatter_panel(h4$data,
  "Hoang Fib ordinal (range log2FC)", "Dream fib ordinal (logFC)",
  "Fib ordinal vs Dream fib ordinal", make_anno(h4))

# Recovery bar chart
rec_dt <- data.table(
  published = rep(c("Hoang NAS\n(2,970)", "Hoang Fib\n(1,656)"), each = 2),
  method = rep(c("Per-study", "Dream"), 2),
  recovery = c(rec_nas, rec_fib)
)
rec_dt[, method := factor(method, levels = c("Per-study", "Dream"))]
pe <- ggplot(rec_dt, aes(x = published, y = recovery, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", recovery)),
            position = position_dodge(0.8), vjust = -0.3, size = 2) +
  scale_fill_manual(values = c("Per-study" = masld_colors$down,
                               "Dream" = masld_colors$up), name = NULL) +
  labs(x = NULL, y = "Recovery (%)", title = "Published DEG recovery") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_masld()

hoang_fig <- (pa | pb | pc) / (pd | pe | plot_spacer()) +
  plot_annotation(tag_levels = "a", title = "Hoang et al. 2019 (GSE130970)") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(hoang_fig, file.path(OUTDIR, "hoang_comparison.pdf"),
         width = fig_full_width, height = 5.5)
cat("  Saved hoang_comparison.pdf\n")

# --- Page 2: Govaere pairwise (up to 6 panels + 25-gene + summary) ---
govaere_panels <- list()
panel_contrasts <- c("NASH_F2_vs_NAFL", "NASH_F3_vs_NAFL", "NASH_F4_vs_NAFL",
                     "NASH_F3_vs_NASH_F01", "NASH_F4_vs_NASH_F01", "NAS_ge4")
contrast_labels <- c("NASH F2 vs NAFL", "NASH F3 vs NAFL", "NASH F4 vs NAFL",
                     "NASH F3 vs F0/1", "NASH F4 vs F0/1", "NAS >= 4")

for (i in seq_along(panel_contrasts)) {
  cname <- panel_contrasts[i]
  if (!cname %in% names(govaere_supp) || !cname %in% names(our_govaere)) {
    govaere_panels[[i]] <- placeholder(paste(contrast_labels[i], "\n(not available)"))
    next
  }
  pub <- govaere_supp[[cname]]
  our <- our_govaere[[cname]]
  if (!pub$has_lfc) {
    govaere_panels[[i]] <- placeholder(paste(contrast_labels[i], "\n(no LFC)"))
    next
  }
  lfc_res <- lfc_correlation(pub$data, our, "pub_logFC", "logFC", "gene_base")
  govaere_panels[[i]] <- scatter_panel(lfc_res$data,
    paste("Govaere", contrast_labels[i], "(logFC)"),
    "Our reproduced (logFC)",
    contrast_labels[i], make_anno(lfc_res))
}

# 25-gene dot plot
g1_plot <- melt(g1_merge[, .(symbol, published_direction, logFC, adj.P.Val,
                              dream_logFC, dream_padj)],  # C2-OK-sensitivity: dream arm of published comparison figure
                id.vars = c("symbol", "published_direction"),
                measure.vars = list(lfc = c("logFC", "dream_logFC"),  # C2-OK-sensitivity
                                    padj = c("adj.P.Val", "dream_padj")))  # C2-OK-sensitivity
g1_plot[, source := fifelse(variable == 1, "Per-study", "Dream")]
g1_plot[, sig := padj < 0.05]
g1_plot[, symbol := factor(symbol, levels = g1_merge[order(logFC)]$symbol)]

p_25gene <- ggplot(g1_plot, aes(x = source, y = symbol)) +
  geom_point(aes(color = ifelse(lfc > 0, "up", "down"),
                 size = pmin(-log10(padj), 10), shape = sig)) +
  scale_color_manual(values = c(up = masld_colors$up, down = masld_colors$down),
                     guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c("NS", "Sig"), name = NULL) +
  scale_size_continuous(range = c(0.5, 3), name = expression(-log[10](padj))) +
  labs(x = NULL, y = NULL, title = "25-gene signature") +
  theme_masld() + theme(axis.text.y = element_text(size = 5))

# Correlation summary across all Govaere contrasts
summary_dt <- rbindlist(lapply(metrics[grep("^G_", names(metrics))], function(x) {
  data.table(contrast = x$comparison, r = x$pearson, rho = x$spearman,
             dir = x$dir_conc, pub_n = x$pub_n, our_n = x$our_n,
             overlap = x$overlap, jaccard = x$jaccard)
}))
if (nrow(summary_dt) > 0) {
  summary_dt[, contrast := gsub("_reproduced$", "", contrast)]
  summary_dt[, contrast := factor(contrast, levels = rev(contrast))]
  p_summary <- ggplot(summary_dt, aes(x = r, y = contrast)) +
    geom_point(aes(size = pub_n), color = masld_colors$up) +
    geom_text(aes(label = sprintf("%.2f", r)), hjust = -0.3, size = 2) +
    scale_size_continuous(range = c(1, 4), name = "Published\nDEGs") +
    labs(x = "Pearson r (published vs reproduced LFC)",
         y = NULL, title = "Correlation summary") +
    xlim(0, 1) +
    theme_masld()
} else {
  p_summary <- placeholder("No Govaere contrasts\navailable")
}

# Assemble Govaere figure
govaere_fig <- (govaere_panels[[1]] | govaere_panels[[2]] | govaere_panels[[3]]) /
               (govaere_panels[[4]] | govaere_panels[[5]] | govaere_panels[[6]]) /
               (p_25gene | p_summary | plot_spacer()) +
  plot_annotation(tag_levels = "a", title = "Govaere et al. 2020 (GSE135251)") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(govaere_fig, file.path(OUTDIR, "govaere_comparison.pdf"),
         width = fig_full_width, height = 8)
cat("  Saved govaere_comparison.pdf\n")

# ── Step 9: Output tables ───────────────────────────────────────────────────
metrics_dt <- rbindlist(lapply(metrics, as.data.table), fill = TRUE)
fwrite(metrics_dt, file.path(OUTDIR, "published_vs_perstudy_metrics.csv"))
fwrite(g1_merge, file.path(OUTDIR, "govaere_25gene_detail.csv"))
cat(sprintf("  Saved metrics (%d rows) and 25-gene detail\n", nrow(metrics_dt)))

cat("\n=== Done ===\n")
