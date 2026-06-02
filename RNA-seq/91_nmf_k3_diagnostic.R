#!/usr/bin/env Rscript
# 91_nmf_k3_diagnostic.R
# NMF k=3 diagnostic: determine whether cached k=3 fit cleanly decomposes
# MASLD samples into Metabolic / Inflammatory / Fibrotic programs.
#
# Reads nmf_results_cache.rds (already contains k=2..6), extracts k=3 W/H,
# runs per-program pathway enrichment + clinical correlations + F-stage boundary test,
# cross-tabs k=3 dominant program vs current S1/S2, and writes a decision
# recommendation (k=3 refactor vs hybrid Topic-Program-3 path).
#
# NOTE: "F2 switch" framing retired 2026-05-09. F-stage boundary tests retained
# as diagnostics but are NOT central to the narrative (see paper_outline.md).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(NMF)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

cache_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache.rds")
meta_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
s12_path   <- file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
out_dir    <- file.path(BASE, "RNA-seq/results/subtypes")
report_md  <- file.path(out_dir, "nmf_k3_diagnostic.md")
report_csv <- file.path(out_dir, "nmf_k3_program_summary.csv")
pdf_out    <- file.path(FIG_MISC, "nmf_k3_diagnostic.pdf")

K <- 3L
FIBROTIC_CANON <- c("COL1A1","COL1A2","COL3A1","ACTA2","TIMP1","LUM",
                    "DCN","FBN1","LOX","SPARC","COL15A1","COL6A3",
                    "CTSK","CDH11","LTBP2","PDGFRB","THY1","VIM")
METABOLIC_CANON <- c("CYP3A4","CYP2E1","ALB","APOB","APOA1","G6PC",
                     "PCK1","HMGCS2","CPT1A","ACOX1","SLC10A1","ABCB11")
INFLAM_CANON    <- c("CCL2","CXCL10","IL6","TNF","NFKB1","STAT1",
                     "IRF1","CD14","CD68","IL1B","NLRP3","S100A8")

# ============================================================
# 1. Load cache + metadata + existing S1/S2 labels
# ============================================================
cat("=== Step 1: Load cache + metadata ===\n")
cached <- readRDS(cache_path)
stopifnot(all(c("nmf_results","mat","mat_nn") %in% names(cached)))
if (!as.character(K) %in% names(cached$nmf_results)) {
  stop(sprintf("k=%d not present in cache (keys: %s)", K,
               paste(names(cached$nmf_results), collapse=",")))
}

res_k <- cached$nmf_results[[as.character(K)]]
W <- basis(res_k)   # genes x K loadings
H <- coef(res_k)    # K x samples coefficients
sample_ids <- colnames(cached$mat_nn)
colnames(H) <- sample_ids
gene_ids   <- rownames(cached$mat_nn)
rownames(W) <- gene_ids
cat(sprintf("  W: %d genes x %d programs\n", nrow(W), ncol(W)))
cat(sprintf("  H: %d programs x %d samples\n", nrow(H), ncol(H)))

meta <- fread(meta_path)
meta <- meta[sample_id %in% sample_ids]
setkey(meta, sample_id)
s12 <- fread(s12_path)
setkey(s12, sample_id)

# Ensembl -> symbol map (strip version)
atlas <- fread(atlas_path, select = c("ensembl_id","human_symbol"))
id2sym <- setNames(atlas$human_symbol, atlas$ensembl_id)
strip_ver <- function(x) sub("\\.\\d+$", "", x)
sym_W <- id2sym[strip_ver(gene_ids)]
names(sym_W) <- gene_ids
n_mapped <- sum(!is.na(sym_W))
cat(sprintf("  Mapped %d/%d genes to symbols\n", n_mapped, length(gene_ids)))

# ============================================================
# 2. Per-program top genes + Hallmark fgsea
# ============================================================
cat("\n=== Step 2: Per-program top genes + fgsea ===\n")
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hallmark$gene_symbol, hallmark$gs_name)

top_genes_per_prog <- list()
fgsea_per_prog     <- list()

# Contrastive loading: W_centered[, j] = W[, j] - mean(W[, -j])
# This gives genes UNIQUELY high in program j positive values and genes
# shared/low in program j negative values — required for meaningful fgsea
# (raw W is all positive, which collapses "std" scoreType).
W_centered <- sapply(seq_len(K), function(j) W[, j] - rowMeans(W[, -j, drop = FALSE]))
rownames(W_centered) <- rownames(W)

for (j in seq_len(K)) {
  w_j <- W_centered[, j]
  ord <- order(w_j, decreasing = TRUE)
  top50_ids  <- gene_ids[ord][1:50]
  top50_sym_raw <- sym_W[top50_ids]
  ok <- !is.na(top50_sym_raw) & top50_sym_raw != ""
  top_genes_per_prog[[j]] <- data.table(
    program = paste0("P", j),
    rank    = unname(which(ok)),
    symbol  = unname(top50_sym_raw[ok]),
    loading = unname(w_j[ord][ok])
  )

  # Rank full contrastive vector by symbol for fgsea
  has_sym <- !is.na(sym_W) & sym_W != ""
  ranked  <- w_j[has_sym]
  names(ranked) <- sym_W[has_sym]
  ranked  <- ranked[!duplicated(names(ranked))]
  ranked  <- sort(ranked, decreasing = TRUE)

  set.seed(42)
  fg <- fgsea(pathways = hall_list, stats = ranked,
              minSize = 10, maxSize = 500, nPermSimple = 10000,
              scoreType = "std")
  fg$program <- paste0("P", j)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, FUN.VALUE=character(1), collapse=";")
  fgsea_per_prog[[j]] <- as.data.table(fg)
  cat(sprintf("  P%d top-10 genes: %s\n", j,
              paste(head(top_genes_per_prog[[j]]$symbol, 10), collapse=", ")))
}

top_genes_dt <- rbindlist(top_genes_per_prog)
fgsea_dt     <- rbindlist(fgsea_per_prog)

# Label each program by best Hallmark match
label_program <- function(fg_sub) {
  sig <- fg_sub[padj < 0.05 & NES > 0][order(-NES)]
  fib_hits <- sig[grepl("EPITHELIAL_MESENCHYMAL|TGF_BETA|MATRISOME|ANGIOGENESIS|COAGULATION|COMPLEMENT", pathway)]
  inf_hits <- sig[grepl("TNFA|INFLAMMATORY|INTERFERON|IL6|IL2|ALLOGRAFT|APOPTOSIS", pathway)]
  met_hits <- sig[grepl("OXIDATIVE_PHOSPH|FATTY_ACID|BILE_ACID|XENOBIOTIC|ADIPOGENESIS|GLYCOLYSIS|PEROXISOME", pathway)]
  scores <- c(
    Fibrotic    = if (nrow(fib_hits)) max(fib_hits$NES) else 0,
    Inflammatory = if (nrow(inf_hits)) max(inf_hits$NES) else 0,
    Metabolic   = if (nrow(met_hits)) max(met_hits$NES) else 0
  )
  if (max(scores) <= 0) return("Unlabeled")
  names(scores)[which.max(scores)]
}
prog_labels <- vapply(seq_len(K), function(j) label_program(fgsea_dt[program == paste0("P", j)]),
                      FUN.VALUE = character(1))
names(prog_labels) <- paste0("P", seq_len(K))
cat("  Program labels:\n"); print(prog_labels)

# Canonical gene hit counts (dedup'd; cap at canon size just in case)
count_canon <- function(j, canon) {
  syms <- unique(as.character(top_genes_per_prog[[j]]$symbol))
  sum(syms %in% canon)
}
canon_tbl <- data.table(
  program     = paste0("P", seq_len(K)),
  label       = prog_labels,
  n_fibrotic  = vapply(seq_len(K), count_canon, canon=FIBROTIC_CANON,  FUN.VALUE=integer(1)),
  n_metabolic = vapply(seq_len(K), count_canon, canon=METABOLIC_CANON, FUN.VALUE=integer(1)),
  n_inflam    = vapply(seq_len(K), count_canon, canon=INFLAM_CANON,    FUN.VALUE=integer(1))
)
cat("\nCanonical gene hits in top-50:\n"); print(canon_tbl)

# ============================================================
# 3. Clinical associations per program
# ============================================================
cat("\n=== Step 3: Clinical associations ===\n")
clin_long <- meta[, .(sample_id, dataset, sex, age, fibrosis_stage, nas_score)]
H_df <- data.frame(sample_id = sample_ids,
                   P1 = H[1, ], P2 = H[2, ], P3 = H[3, ],
                   check.names = FALSE, stringsAsFactors = FALSE)
H_dt <- as.data.table(H_df)
setkey(H_dt, sample_id)
comb <- merge(clin_long, H_dt, by = "sample_id")

clin_assoc <- data.table()
for (j in seq_len(K)) {
  p <- paste0("P", j)
  # Fibrosis: spearman
  sf <- tryCatch(cor.test(comb[[p]], as.numeric(comb$fibrosis_stage), method = "spearman",
                          exact = FALSE), error = function(e) list(estimate=NA, p.value=NA))
  # NAS
  sn <- tryCatch(cor.test(comb[[p]], as.numeric(comb$nas_score), method = "spearman",
                          exact = FALSE), error = function(e) list(estimate=NA, p.value=NA))
  # Age
  sa <- tryCatch(cor.test(comb[[p]], as.numeric(comb$age), method = "spearman",
                          exact = FALSE), error = function(e) list(estimate=NA, p.value=NA))
  # Sex: point-biserial via spearman on factor
  sx <- comb[!is.na(sex) & sex %in% c("M","F")]
  ss <- if (nrow(sx) > 10) {
    tryCatch(cor.test(sx[[p]], as.integer(sx$sex == "F"), method = "spearman",
                      exact = FALSE), error = function(e) list(estimate=NA, p.value=NA))
  } else list(estimate=NA, p.value=NA)
  # Dataset eta^2 from KW
  kw <- tryCatch(kruskal.test(comb[[p]] ~ factor(comb$dataset)),
                 error = function(e) list(statistic=NA, p.value=NA, parameter=NA))
  # eta^2_H = (H - k + 1) / (n - k)  where k = n_groups
  n_ds <- length(unique(comb$dataset))
  eta2 <- if (!is.na(kw$statistic)) max(0, (as.numeric(kw$statistic) - n_ds + 1) / (nrow(comb) - n_ds)) else NA_real_

  clin_assoc <- rbind(clin_assoc, data.table(
    program = p, label = prog_labels[p],
    rho_fibrosis = unname(sf$estimate), p_fibrosis = sf$p.value,
    rho_nas      = unname(sn$estimate), p_nas      = sn$p.value,
    rho_age      = unname(sa$estimate), p_age      = sa$p.value,
    rho_sex_F    = unname(ss$estimate), p_sex      = ss$p.value,
    kw_dataset_stat = as.numeric(kw$statistic),
    kw_dataset_p    = as.numeric(kw$p.value),
    eta2_dataset    = eta2
  ))
}
cat("  Clinical associations:\n"); print(clin_assoc[, .(program, label, rho_fibrosis, rho_nas, rho_sex_F, eta2_dataset)])

# ============================================================
# 4. F-stage boundary test (F0-2 vs F3-4)
# NOTE: "F2 switch" framing retired 2026-05-09; see paper_outline.md
# ============================================================
cat("\n=== Step 4: F-stage boundary test (F0-2 vs F3-4) ===\n")
f_comb <- comb[!is.na(fibrosis_stage)]
f_comb[, f_bin := ifelse(fibrosis_stage <= 2, "F0-F2", "F3-F4")]
switch_tbl <- data.table()
for (j in seq_len(K)) {
  p <- paste0("P", j)
  wt <- tryCatch(wilcox.test(f_comb[[p]] ~ f_comb$f_bin),
                 error = function(e) list(p.value=NA))
  m_lo <- mean(f_comb[f_bin == "F0-F2"][[p]], na.rm = TRUE)
  m_hi <- mean(f_comb[f_bin == "F3-F4"][[p]], na.rm = TRUE)
  switch_tbl <- rbind(switch_tbl, data.table(
    program = p, label = prog_labels[p],
    mean_F02 = m_lo, mean_F34 = m_hi,
    switch_ratio = m_hi / (m_lo + 1e-10),
    wilcox_p = wt$p.value
  ))
}
cat("  F-stage boundary (F0-2 vs F3-4):\n"); print(switch_tbl)

# ============================================================
# 5. Cross-tab dominant k=3 program vs S1/S2
# ============================================================
cat("\n=== Step 5: Cross-tab vs S1/S2 ===\n")
dom_prog <- apply(H, 2, function(col) paste0("P", which.max(col)))
xt <- merge(data.table(sample_id = sample_ids, dom_k3 = dom_prog),
            s12[, .(sample_id, nmf_subtype)], by = "sample_id")
cross_tab <- xt[, .N, by = .(dom_k3, nmf_subtype)]
cross_tab_wide <- dcast(cross_tab, dom_k3 ~ nmf_subtype, value.var = "N", fill = 0L)
cat("  3x2 cross-tab:\n"); print(cross_tab_wide)

# ============================================================
# 6. Orthogonality
# ============================================================
cat("\n=== Step 6: Orthogonality ===\n")
H_cor <- cor(t(H))
rownames(H_cor) <- colnames(H_cor) <- paste0("P", seq_len(K))
cat("  Pearson correlation of H rows:\n"); print(round(H_cor, 3))

# Top-500 gene overlap (Jaccard)
top500 <- lapply(seq_len(K), function(j) {
  sym <- sym_W[order(W[, j], decreasing = TRUE)]
  sym <- sym[!is.na(sym) & sym != ""]
  head(sym, 500)
})
jacc_mat <- matrix(NA_real_, K, K, dimnames = list(paste0("P", seq_len(K)), paste0("P", seq_len(K))))
for (i in seq_len(K)) for (j in seq_len(K)) {
  a <- top500[[i]]; b <- top500[[j]]
  jacc_mat[i, j] <- length(intersect(a, b)) / length(union(a, b))
}
cat("  Top-500 Jaccard:\n"); print(round(jacc_mat, 3))

# ============================================================
# 7. Decision
# ============================================================
cat("\n=== Step 7: Decision ===\n")
fib_prog_idx <- which(prog_labels == "Fibrotic")
decision <- "HYBRID"
decision_reasons <- character()

if (length(fib_prog_idx) != 1L) {
  decision_reasons <- c(decision_reasons,
                        sprintf("FAIL cond 1: expected exactly 1 Fibrotic program, found %d", length(fib_prog_idx)))
} else {
  fp <- paste0("P", fib_prog_idx)
  fg_fp <- fgsea_dt[program == fp & padj < 0.05 & NES > 0]
  fib_nes <- fg_fp[grepl("EPITHELIAL_MESENCHYMAL|TGF_BETA|MATRISOME", pathway)]
  has_fib_sig <- nrow(fib_nes) > 0 && max(fib_nes$NES) >= 2
  canon_hits <- canon_tbl[program == fp, n_fibrotic]
  cond1 <- has_fib_sig && canon_hits >= 3
  if (!cond1) decision_reasons <- c(decision_reasons,
                                    sprintf("FAIL cond 1: Fibrotic program %s Hallmark NES>=2 for EMT/TGFB/Matrisome = %s; canonical fibrotic genes in top-50 = %d",
                                            fp, has_fib_sig, canon_hits))

  fp_clin <- clin_assoc[program == fp]
  cond2 <- abs(fp_clin$rho_fibrosis) >= 0.3 &&
           (is.na(fp_clin$rho_sex_F) || abs(fp_clin$rho_sex_F) < 0.2) &&
           fp_clin$eta2_dataset < 0.2
  if (!cond2) decision_reasons <- c(decision_reasons,
                                    sprintf("FAIL cond 2: Fibrotic %s rho_fibrosis=%.3f, rho_sex=%.3f, eta2_dataset=%.3f",
                                            fp, fp_clin$rho_fibrosis, fp_clin$rho_sex_F, fp_clin$eta2_dataset))

  offdiag <- H_cor[lower.tri(H_cor)]
  cond3 <- all(abs(offdiag) < 0.5)
  if (!cond3) decision_reasons <- c(decision_reasons,
                                    sprintf("FAIL cond 3: max |H correlation| = %.3f (>=0.5)", max(abs(offdiag))))

  fp_sw <- switch_tbl[program == fp]
  cond4 <- !is.na(fp_sw$switch_ratio) && fp_sw$switch_ratio > 1.5
  if (!cond4) decision_reasons <- c(decision_reasons,
                                    sprintf("FAIL cond 4: switch ratio = %.3f (<=1.5)", fp_sw$switch_ratio))

  if (cond1 && cond2 && cond3 && cond4) decision <- "K3_REFACTOR"
}

cat(sprintf("\nDECISION: %s\n", decision))
if (length(decision_reasons)) cat(paste(decision_reasons, collapse = "\n"), "\n")

# ============================================================
# 8. Write CSV + Markdown report
# ============================================================
cat("\n=== Step 8: Write report ===\n")
prog_summary <- merge(merge(canon_tbl, clin_assoc[, .(program, rho_fibrosis, rho_nas, rho_sex_F, eta2_dataset)],
                            by = "program"),
                      switch_tbl[, .(program, switch_ratio, wilcox_p)],
                      by = "program")
fwrite(prog_summary, report_csv)

top_hall <- fgsea_dt[, .SD[order(padj)][1:5], by = program][, .(program, pathway, NES, padj)]

md <- c(
  "# NMF k=3 Diagnostic",
  sprintf("**Generated:** %s  ", format(Sys.time())),
  sprintf("**Input:** `%s`  ", basename(cache_path)),
  sprintf("**k=3 stability:** cophenetic=0.964, silhouette=0.841 (from `nmf_metrics.csv`)"),
  "",
  sprintf("## DECISION: **%s**", decision),
  "",
  if (decision == "K3_REFACTOR") {
    "All four conditions met — proceed with Script 44 k=3 refactor + downstream (208/209/217/figures)."
  } else {
    paste0("At least one condition failed — fall back to **hybrid path** (keep Script 44 at k=2, add Topic Program 3 from Script 230 as continuous `fibrotic_score`).",
           "\n\n### Failed conditions\n",
           paste0("- ", decision_reasons, collapse = "\n"))
  },
  "",
  "## Per-program summary",
  "",
  paste(capture.output(print(prog_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Top Hallmark pathways per program (top 5 by padj)",
  "",
  paste(capture.output(print(top_hall, row.names = FALSE, nrows = 999)), collapse = "\n"),
  "",
  "## Canonical gene hits in top-50 loadings",
  "",
  paste(capture.output(print(canon_tbl, row.names = FALSE)), collapse = "\n"),
  "",
  "## Top-10 loading genes per program",
  ""
)
for (j in seq_len(K)) {
  md <- c(md, sprintf("**P%d (%s):**  %s  ", j, prog_labels[j],
                       paste(top_genes_per_prog[[j]]$symbol[1:10], collapse = ", ")), "")
}
md <- c(md,
  "## Cross-tab dominant k=3 program vs current S1/S2",
  "",
  paste(capture.output(print(cross_tab_wide, row.names = FALSE)), collapse = "\n"),
  "",
  "## Orthogonality",
  "",
  "### Pearson of H rows",
  paste(capture.output(print(round(H_cor, 3))), collapse = "\n"),
  "",
  "### Top-500 loading-gene Jaccard",
  paste(capture.output(print(round(jacc_mat, 3))), collapse = "\n"),
  "",
  "## Decision criteria",
  "1. Exactly one Fibrotic-labeled program with EMT/TGFB/Matrisome Hallmark NES >= 2 AND >=3 canonical fibrotic genes in top-50",
  "2. Fibrotic program: |rho_fibrosis| >= 0.3 AND |rho_sex| < 0.2 AND eta2_dataset < 0.2",
  "3. All pairwise |H correlations| < 0.5",
  "4. F-stage boundary ratio (F3-4 / F0-2) > 1.5 for Fibrotic program"
)
writeLines(md, report_md)
cat(sprintf("  Wrote %s\n", report_md))
cat(sprintf("  Wrote %s\n", report_csv))

# ============================================================
# 9. Diagnostic PDF (3 panels A/B/C)
# ============================================================
cat("\n=== Step 9: Diagnostic PDF ===\n")

# Panel A — top-30 loading heatmap (scaled W per program)
topA <- lapply(seq_len(K), function(j) {
  ord <- order(W[, j], decreasing = TRUE)
  ids <- gene_ids[ord][1:30]
  syms <- sym_W[ids]
  data.table(gene_id = ids, symbol = syms, program = paste0("P", j), rank = seq_along(ids))
}) |> rbindlist()
heat_ids  <- unique(topA$gene_id)
heat_syms <- sym_W[heat_ids]
heat_syms[is.na(heat_syms) | heat_syms == ""] <- heat_ids[is.na(heat_syms) | heat_syms == ""]
heat_mat  <- W[heat_ids, , drop = FALSE]
rownames(heat_mat) <- heat_syms
heat_mat_z <- t(scale(t(heat_mat)))
colnames(heat_mat_z) <- paste0("P", seq_len(K), " (", prog_labels, ")")

# Convert pheatmap to grob via gridExtra — use base grid ggplot-compatible alternative
heat_df <- as.data.table(heat_mat_z, keep.rownames = "symbol")
heat_long <- melt(heat_df, id.vars = "symbol", variable.name = "program", value.name = "z")
# Order rows by program assignment of max z
row_order <- heat_long[, .(max_prog = as.character(program[which.max(z)])), by = symbol]
row_order[, max_prog := factor(max_prog, levels = colnames(heat_mat_z))]
setorder(row_order, max_prog, symbol)
heat_long[, symbol := factor(symbol, levels = row_order$symbol)]
heat_long[, program := factor(program, levels = colnames(heat_mat_z))]

pA <- ggplot(heat_long, aes(program, symbol, fill = z)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#1565C0", mid = "grey95", high = "#C2185B",
                       midpoint = 0, name = "z(W)") +
  labs(x = NULL, y = NULL, title = "A — Top-30 loading genes per program (scaled)") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        legend.position = "right")

# Panel B — H violin by fibrosis stage, faceted by program
long_H <- melt(H_dt, id.vars = "sample_id", variable.name = "program", value.name = "H_val")
long_H <- merge(long_H, meta[, .(sample_id, fibrosis_stage)], by = "sample_id")
long_H <- long_H[!is.na(fibrosis_stage)]
long_H[, fstage := factor(paste0("F", fibrosis_stage), levels = paste0("F", 0:4))]
long_H[, prog_lab := paste0(program, " (", prog_labels[as.character(program)], ")")]

pB <- ggplot(long_H, aes(fstage, H_val, fill = fstage)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, color = "grey30") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, outlier.alpha = 0.4, fill = "white") +
  facet_wrap(~ prog_lab, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("F0"="#90CAF9","F1"="#64B5F6","F2"="#F48FB1",
                               "F3"="#C2185B","F4"="#880E4F"), guide = "none") +
  labs(x = "Fibrosis stage", y = "H coefficient",
       title = "B — Program H coefficient by fibrosis stage") +
  theme_masld() +
  theme(strip.text = element_text(size = 8))

# Panel C — 3x2 contingency (stacked bar)
ct_dt <- cross_tab[, .(dom_k3, nmf_subtype, N)]
ct_dt[, dom_lab := paste0(dom_k3, " (", prog_labels[dom_k3], ")")]
pC <- ggplot(ct_dt, aes(nmf_subtype, N, fill = dom_lab)) +
  geom_col(position = "fill", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = N), position = position_fill(vjust = 0.5), size = 2.8, color = "white") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2", name = "Dominant k=3 program") +
  labs(x = "Current subtype (k=2)", y = "Proportion",
       title = "C — k=3 dominant program composition within each current subtype") +
  theme_masld()

fig <- (pA | pB) / pC +
  plot_annotation(title = sprintf("NMF k=3 diagnostic — decision: %s", decision),
                  theme = theme(plot.title = element_text(size = 11, face = "bold"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(pdf_out, fig, width = 16, height = 11, device = cairo_pdf)
cat(sprintf("  Wrote %s\n", pdf_out))
cat("\nDone.\n")
