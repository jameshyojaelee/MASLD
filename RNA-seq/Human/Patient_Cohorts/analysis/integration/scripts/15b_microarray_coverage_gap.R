#!/usr/bin/env Rscript
# ============================================================
# 15b: Microarray Coverage Gap Analysis
#   Quantifies the fraction of RNA-seq DEGs invisible to the
#   two most common MASLD microarray platforms:
#     - Affymetrix Human Genome U133 Plus 2.0 (GPL570)
#     - Illumina HumanHT-12 v4 BeadChip (GPL10558)
#   Shows that mixed microarray+RNA-seq meta-analyses structurally
#   miss lncRNAs and low-abundance transcripts.
#
# Approach:
#   - Downloads GPL annotation tables from NCBI GEO FTP (stable)
#   - Uses biomaRt as supplementary source (with fallback if offline)
#   - Maps probe gene symbols to Ensembl IDs via atlas + org.Hs.eg.db
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/microarray_gap")
FIGDIR <- file.path(BASE, "figures/supplementary/figS_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGDIR, "panels"), recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(OUTDIR, "probe_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== 15b: Microarray Coverage Gap Analysis ===\n")

# ============================================================
# 1. Load dream results + gene biotype annotations
# ============================================================
cat("Loading dream results...\n")
dream <- fread(file.path(INT, "results/integration/dream_results.csv"))
cat(sprintf("  Total genes in dream: %d\n", nrow(dream)))

# Strip version suffix from Ensembl IDs for matching
dream[, ensembl_clean := sub("\\.\\d+$", "", gene)]

# Load biotype from multi-evidence atlas
cat("Loading gene biotype annotations from atlas...\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("human_symbol", "ensembl_id", "gene_biotype"))
atlas[, ensembl_clean := sub("\\.\\d+$", "", ensembl_id)]

# Merge biotype into dream
dream <- merge(dream, atlas[, .(ensembl_clean, human_symbol, gene_biotype)],
               by = "ensembl_clean", all.x = TRUE)

# Fill missing biotypes
n_missing <- sum(is.na(dream$gene_biotype))
if (n_missing > 0) {
  cat(sprintf("  %d genes missing biotype — marking as unknown\n", n_missing))
  dream[is.na(gene_biotype), gene_biotype := "unknown"]
}

# Define DEGs
dream[, is_deg := padj < 0.1]
n_degs <- sum(dream$is_deg, na.rm = TRUE)
cat(sprintf("  DEGs (padj < 0.1): %d\n", n_degs))

# Biotype summary for DEGs
cat("  DEG biotype breakdown:\n")
print(dream[is_deg == TRUE, .N, by = gene_biotype][order(-N)])

# ============================================================
# 2. Build symbol <-> Ensembl lookup from atlas + org.Hs.eg.db
# ============================================================
cat("\nBuilding comprehensive symbol-to-Ensembl mapping...\n")

# Source 1: atlas
sym_map <- atlas[human_symbol != "" & !is.na(human_symbol),
                 .(symbol = human_symbol, ensembl_clean)]
sym_map <- unique(sym_map)

# Source 2: org.Hs.eg.db (fills gaps for genes not in atlas)
org_map <- tryCatch({
  res <- AnnotationDbi::select(org.Hs.eg.db,
    keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
    keytype = "ENSEMBL", columns = c("ENSEMBL", "SYMBOL"))
  setDT(res)
  res <- res[!is.na(SYMBOL) & SYMBOL != "" & !is.na(ENSEMBL)]
  res[, .(symbol = SYMBOL, ensembl_clean = ENSEMBL)]
}, error = function(e) {
  cat("  Warning: org.Hs.eg.db lookup failed, using atlas only\n")
  data.table(symbol = character(), ensembl_clean = character())
})

sym_map <- unique(rbind(sym_map, org_map))
sym_map <- sym_map[!duplicated(paste(symbol, ensembl_clean))]
cat(sprintf("  Symbol-to-Ensembl map: %d symbol-gene pairs (%d unique symbols, %d unique Ensembl IDs)\n",
            nrow(sym_map), uniqueN(sym_map$symbol), uniqueN(sym_map$ensembl_clean)))

# ============================================================
# 3. Fetch microarray probe annotations
# ============================================================

# --- Helper: parse NCBI GEO GPL annotation file ---
parse_gpl_annot <- function(filepath) {
  raw <- readLines(gzfile(filepath))
  # Skip header comment lines (start with # or !)
  data_start <- which(!grepl("^[#!^]", raw))[1]
  dt <- fread(text = raw[data_start:length(raw)], fill = TRUE, sep = "\t")
  dt
}

# --- 3a. GPL570: Affymetrix HG-U133 Plus 2.0 ---
affy_cache <- file.path(cache_dir, "GPL570_genes.csv.gz")
if (file.exists(affy_cache)) {
  cat("Loading cached GPL570 (Affy U133 Plus 2.0) gene list...\n")
  affy_symbols <- fread(affy_cache)$symbol
} else {
  cat("Downloading GPL570 (Affy HG-U133 Plus 2.0) annotation from NCBI GEO...\n")
  gpl570_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
  tmp <- tempfile(fileext = ".annot.gz")
  download.file(gpl570_url, tmp, quiet = TRUE, method = "libcurl")

  gpl570 <- parse_gpl_annot(tmp)
  unlink(tmp)

  # Identify the gene symbol column
  sym_col <- intersect(c("Gene symbol", "Gene Symbol", "Symbol", "GENE_SYMBOL"), names(gpl570))
  if (length(sym_col) == 0) {
    cat("  GPL570 columns: ", paste(names(gpl570), collapse = ", "), "\n")
    stop("Cannot find gene symbol column in GPL570 annotation")
  }
  sym_col <- sym_col[1]
  cat(sprintf("  Using column '%s' for gene symbols\n", sym_col))

  symbols_raw <- gpl570[[sym_col]]
  symbols_raw <- symbols_raw[!is.na(symbols_raw) & symbols_raw != "" & symbols_raw != "---"]

  # Expand multi-symbol entries (separated by " /// ")
  affy_symbols <- unique(unlist(strsplit(symbols_raw, " /// ")))
  affy_symbols <- affy_symbols[affy_symbols != "" & affy_symbols != "---"]

  fwrite(data.table(symbol = affy_symbols), affy_cache)
  cat(sprintf("  GPL570: %d unique gene symbols\n", length(affy_symbols)))
}
cat(sprintf("  Affy HG-U133 Plus 2.0: %d unique gene symbols\n", length(affy_symbols)))

# --- 3b. GPL10558: Illumina HumanHT-12 v4 ---
ill_cache <- file.path(cache_dir, "GPL10558_genes.csv.gz")
if (file.exists(ill_cache)) {
  cat("Loading cached GPL10558 (Illumina HT-12 v4) gene list...\n")
  ill_symbols <- fread(ill_cache)$symbol
} else {
  cat("Downloading GPL10558 (Illumina HumanHT-12 v4) annotation from NCBI GEO...\n")
  gpl10558_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
  tmp <- tempfile(fileext = ".annot.gz")
  download.file(gpl10558_url, tmp, quiet = TRUE, method = "libcurl")

  gpl10558 <- parse_gpl_annot(tmp)
  unlink(tmp)

  sym_col <- intersect(c("Gene symbol", "Gene Symbol", "Symbol", "GENE_SYMBOL"), names(gpl10558))
  if (length(sym_col) == 0) {
    cat("  GPL10558 columns: ", paste(names(gpl10558), collapse = ", "), "\n")
    stop("Cannot find gene symbol column in GPL10558 annotation")
  }
  sym_col <- sym_col[1]
  cat(sprintf("  Using column '%s' for gene symbols\n", sym_col))

  symbols_raw <- gpl10558[[sym_col]]
  symbols_raw <- symbols_raw[!is.na(symbols_raw) & symbols_raw != "" & symbols_raw != "---"]

  ill_symbols <- unique(unlist(strsplit(symbols_raw, " /// ")))
  ill_symbols <- ill_symbols[ill_symbols != "" & ill_symbols != "---"]

  fwrite(data.table(symbol = ill_symbols), ill_cache)
  cat(sprintf("  GPL10558: %d unique gene symbols\n", length(ill_symbols)))
}
cat(sprintf("  Illumina HumanHT-12 v4: %d unique gene symbols\n", length(ill_symbols)))

# --- 3c. Supplement with biomaRt (if available) ---
# biomaRt gives Ensembl-based probe mappings that may catch genes whose
# HGNC symbols differ between the GPL annotation and our atlas.
biomart_affy_cache <- file.path(cache_dir, "biomart_affy_u133plus2.csv.gz")
biomart_affy_ensembl <- character(0)
if (file.exists(biomart_affy_cache)) {
  cat("Loading cached biomaRt Affy probe mappings...\n")
  biomart_affy_ensembl <- fread(biomart_affy_cache)$ensembl_clean
} else {
  cat("Attempting biomaRt Affy HG-U133 Plus 2.0 probe lookup (supplementary)...\n")
  tryCatch({
    library(biomaRt)
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                          mirror = "useast")
    bm_res <- setDT(getBM(
      attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"),
      filters    = "with_affy_hg_u133_plus_2",
      values     = TRUE,
      mart       = ensembl
    ))
    biomart_affy_ensembl <- unique(bm_res$ensembl_gene_id[bm_res$ensembl_gene_id != ""])
    fwrite(data.table(ensembl_clean = biomart_affy_ensembl), biomart_affy_cache)
    cat(sprintf("  biomaRt Affy: %d Ensembl IDs\n", length(biomart_affy_ensembl)))
  }, error = function(e) {
    cat(sprintf("  biomaRt unavailable (%s) — using GPL annotation only\n",
                conditionMessage(e)))
  })
}

biomart_ill_cache <- file.path(cache_dir, "biomart_illumina_wg6v3.csv.gz")
biomart_ill_ensembl <- character(0)
if (file.exists(biomart_ill_cache)) {
  cat("Loading cached biomaRt Illumina probe mappings...\n")
  biomart_ill_ensembl <- fread(biomart_ill_cache)$ensembl_clean
} else {
  cat("Attempting biomaRt Illumina WG-6 v3 probe lookup (supplementary)...\n")
  tryCatch({
    if (!exists("ensembl")) {
      library(biomaRt)
      ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                            mirror = "useast")
    }
    bm_res <- setDT(getBM(
      attributes = c("ensembl_gene_id", "illumina_humanwg_6_v3"),
      filters    = "with_illumina_humanwg_6_v3",
      values     = TRUE,
      mart       = ensembl
    ))
    biomart_ill_ensembl <- unique(bm_res$ensembl_gene_id[bm_res$ensembl_gene_id != ""])
    fwrite(data.table(ensembl_clean = biomart_ill_ensembl), biomart_ill_cache)
    cat(sprintf("  biomaRt Illumina: %d Ensembl IDs\n", length(biomart_ill_ensembl)))
  }, error = function(e) {
    cat(sprintf("  biomaRt unavailable (%s) — using GPL annotation only\n",
                conditionMessage(e)))
  })
}

# ============================================================
# 4. Map symbols to Ensembl IDs and build platform gene sets
# ============================================================
cat("\nMapping platform symbols to Ensembl IDs...\n")

# Affy: GPL symbol-based + biomaRt Ensembl-based
affy_ensembl_from_sym <- sym_map[symbol %in% affy_symbols, unique(ensembl_clean)]
affy_genes <- unique(c(affy_ensembl_from_sym, biomart_affy_ensembl))
cat(sprintf("  Affy gene set: %d from symbols + %d from biomaRt = %d unique Ensembl IDs\n",
            length(affy_ensembl_from_sym), length(biomart_affy_ensembl), length(affy_genes)))

# Illumina: GPL symbol-based + biomaRt Ensembl-based
ill_ensembl_from_sym <- sym_map[symbol %in% ill_symbols, unique(ensembl_clean)]
ill_genes <- unique(c(ill_ensembl_from_sym, biomart_ill_ensembl))
cat(sprintf("  Illumina gene set: %d from symbols + %d from biomaRt = %d unique Ensembl IDs\n",
            length(ill_ensembl_from_sym), length(biomart_ill_ensembl), length(ill_genes)))

# Union
union_genes <- unique(c(affy_genes, ill_genes))
cat(sprintf("  Union (either platform): %d unique Ensembl IDs\n", length(union_genes)))

# ============================================================
# 5. Cross-reference with dream DEGs
# ============================================================
cat("\n--- Cross-referencing with dream DEGs ---\n")

dream[, detectable_affy := ensembl_clean %in% affy_genes]
dream[, detectable_illumina := ensembl_clean %in% ill_genes]
dream[, detectable_either := ensembl_clean %in% union_genes]

# Simplify biotype for display
dream[, biotype_simple := fifelse(
  gene_biotype == "protein_coding", "Protein-coding",
  fifelse(gene_biotype == "lncRNA", "lncRNA",
    fifelse(gene_biotype %in% c("miRNA", "snoRNA", "snRNA", "misc_RNA",
                                 "rRNA", "Mt_rRNA", "Mt_tRNA", "scaRNA",
                                 "vault_RNA", "ribozyme", "sRNA"),
            "Other ncRNA",
      fifelse(grepl("pseudogene", gene_biotype), "Pseudogene", "Other"))))]

# --- Summary statistics ---
degs <- dream[is_deg == TRUE]

summary_dt <- degs[, .(
  total           = .N,
  detectable_affy = sum(detectable_affy),
  detectable_ill  = sum(detectable_illumina),
  detectable_any  = sum(detectable_either),
  invisible       = sum(!detectable_either)
), by = biotype_simple]

summary_dt[, pct_invisible := round(100 * invisible / total, 1)]
setorder(summary_dt, -total)

cat("\nDEG detectability by biotype:\n")
print(summary_dt)

# Overall stats
overall <- degs[, .(
  total_degs      = .N,
  detectable_affy = sum(detectable_affy),
  detectable_ill  = sum(detectable_illumina),
  detectable_any  = sum(detectable_either),
  invisible       = sum(!detectable_either),
  pct_invisible   = round(100 * sum(!detectable_either) / .N, 1)
)]
cat(sprintf("\n=== OVERALL: %d/%d DEGs (%.1f%%) invisible to microarray ===\n",
            overall$invisible, overall$total_degs, overall$pct_invisible))

# lncRNA-specific
lncrna_degs <- degs[biotype_simple == "lncRNA"]
if (nrow(lncrna_degs) > 0) {
  cat(sprintf("  lncRNA DEGs: %d total, %d detectable (%.1f%%), %d invisible (%.1f%%)\n",
              nrow(lncrna_degs),
              sum(lncrna_degs$detectable_either),
              100 * mean(lncrna_degs$detectable_either),
              sum(!lncrna_degs$detectable_either),
              100 * mean(!lncrna_degs$detectable_either)))
}

# ============================================================
# 6. Expression bias analysis
# ============================================================
cat("\n--- Expression bias: microarray-detectable vs invisible ---\n")
expr_comparison <- degs[, .(
  mean_AveExpr   = mean(AveExpr, na.rm = TRUE),
  median_AveExpr = median(AveExpr, na.rm = TRUE),
  n              = .N
), by = detectable_either]
setnames(expr_comparison, "detectable_either", "microarray_detectable")
print(expr_comparison)

# Wilcoxon test
wtest <- wilcox.test(
  degs[detectable_either == TRUE, AveExpr],
  degs[detectable_either == FALSE, AveExpr]
)
cat(sprintf("  Wilcoxon p-value (AveExpr detectable vs invisible): %.2e\n", wtest$p.value))

# Effect size (Cohen's d)
d_detectable <- degs[detectable_either == TRUE, AveExpr]
d_invisible  <- degs[detectable_either == FALSE, AveExpr]
pooled_sd <- sqrt(((length(d_detectable) - 1) * var(d_detectable) +
                    (length(d_invisible) - 1) * var(d_invisible)) /
                   (length(d_detectable) + length(d_invisible) - 2))
cohens_d <- (mean(d_detectable) - mean(d_invisible)) / pooled_sd
cat(sprintf("  Cohen's d: %.2f\n", cohens_d))

# ============================================================
# 7. Save results CSV
# ============================================================
cat("\nSaving results...\n")

# Per-gene results
out_per_gene <- dream[, .(
  ensembl_id       = gene,
  ensembl_clean    = ensembl_clean,
  human_symbol     = human_symbol,
  gene_biotype     = gene_biotype,
  biotype_simple   = biotype_simple,
  logFC            = logFC,
  AveExpr          = AveExpr,
  padj             = padj,
  is_deg           = is_deg,
  detectable_affy  = detectable_affy,
  detectable_illumina = detectable_illumina,
  detectable_either = detectable_either
)]
fwrite(out_per_gene, file.path(OUTDIR, "microarray_coverage_gap.csv"))
cat(sprintf("  Per-gene results: %s\n", file.path(OUTDIR, "microarray_coverage_gap.csv")))

# Summary table
fwrite(summary_dt, file.path(OUTDIR, "microarray_gap_summary_by_biotype.csv"))

# Platform stats
platform_stats <- data.table(
  metric = c("Affy_U133Plus2_probed_genes", "Illumina_HT12v4_probed_genes",
             "Union_probed_genes", "Total_dream_genes", "Total_DEGs",
             "DEGs_detectable_affy", "DEGs_detectable_illumina",
             "DEGs_detectable_either", "DEGs_invisible",
             "Pct_DEGs_invisible", "lncRNA_DEGs_total",
             "lncRNA_DEGs_detectable", "lncRNA_DEGs_invisible",
             "Pct_lncRNA_invisible",
             "Mean_AveExpr_detectable", "Mean_AveExpr_invisible",
             "Wilcoxon_pvalue", "Cohens_d"),
  value = c(length(affy_genes), length(ill_genes),
            length(union_genes), nrow(dream), n_degs,
            overall$detectable_affy, overall$detectable_ill,
            overall$detectable_any, overall$invisible,
            overall$pct_invisible, nrow(lncrna_degs),
            sum(lncrna_degs$detectable_either),
            sum(!lncrna_degs$detectable_either),
            round(100 * mean(!lncrna_degs$detectable_either), 1),
            round(mean(d_detectable), 3),
            round(mean(d_invisible), 3),
            signif(wtest$p.value, 3),
            round(cohens_d, 3))
)
fwrite(platform_stats, file.path(OUTDIR, "microarray_platform_stats.csv"))
cat(sprintf("  Platform stats: %s\n", file.path(OUTDIR, "microarray_platform_stats.csv")))

# ============================================================
# 8. Generate figure
# ============================================================
cat("\nGenerating figures...\n")

# --- Panel (a): Stacked bar — DEGs detectable vs not, by biotype ---
bar_dt <- degs[, .(n = .N), by = .(biotype_simple, detectable_either)]
bar_dt[, detectability := fifelse(detectable_either,
                                   "Microarray-detectable", "Microarray-invisible")]

# Order biotypes by total DEG count
biotype_order <- bar_dt[, .(total = sum(n)), by = biotype_simple][order(-total), biotype_simple]
bar_dt[, biotype_simple := factor(biotype_simple, levels = biotype_order)]

pa <- ggplot(bar_dt, aes(x = biotype_simple, y = n, fill = detectability)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c(
    "Microarray-detectable" = masld_colors$deg,
    "Microarray-invisible"  = masld_colors$up
  ), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of DEGs",
       title = "DEG detectability by gene biotype") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = "top")

# Add percentage invisible labels above each bar
pct_labels <- degs[, .(
  total = .N,
  n_invisible = sum(!detectable_either)
), by = biotype_simple]
pct_labels[, pct_raw := 100 * n_invisible / total]
pct_labels[, pct_text := fifelse(pct_raw < 1,
                                  "<1% invisible",
                                  sprintf("%.0f%% invisible", pct_raw))]
pct_labels[, biotype_simple := factor(biotype_simple, levels = biotype_order)]

pa <- pa +
  geom_text(data = pct_labels,
            aes(x = biotype_simple, y = total, label = pct_text),
            inherit.aes = FALSE, vjust = -0.3, size = 1.8, fontface = "bold")

# --- Panel (b): Density — AveExpr by detectability ---
dens_dt <- degs[, .(AveExpr = AveExpr,
                     detectability = fifelse(detectable_either,
                                             "Microarray-detectable",
                                             "Microarray-invisible"))]

pb <- ggplot(dens_dt, aes(x = AveExpr, fill = detectability, color = detectability)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Microarray-detectable" = masld_colors$deg,
    "Microarray-invisible"  = masld_colors$up
  ), name = NULL) +
  scale_color_manual(values = c(
    "Microarray-detectable" = masld_colors$deg,
    "Microarray-invisible"  = masld_colors$up
  ), name = NULL) +
  geom_vline(xintercept = median(d_detectable), linetype = "dashed",
             color = masld_colors$deg, linewidth = 0.3) +
  geom_vline(xintercept = median(d_invisible), linetype = "dashed",
             color = masld_colors$up, linewidth = 0.3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Wilcoxon %s\nCohen's d = %.2f",
             ifelse(wtest$p.value < .Machine$double.xmin,
                    "p < 2.2e-308", sprintf("p = %.1e", wtest$p.value)),
             cohens_d),
           size = 2, fontface = "italic") +
  labs(x = "Average expression (log-CPM)", y = "Density",
       title = "Expression bias: detectable vs invisible DEGs") +
  theme_masld() +
  theme(legend.position = "top")

# --- Panel (c): Volcano colored by detectability ---
volc_dt <- copy(dream)
volc_dt[, neg_log10_p := -log10(pmax(padj, 1e-300))]
volc_dt[, category := fifelse(
  !is_deg, "Not significant",
  fifelse(detectable_either, "Detectable DEG", "Invisible DEG")
)]
volc_dt[, category := factor(category,
  levels = c("Not significant", "Detectable DEG", "Invisible DEG"))]

# Downsample NS points for faster rendering
set.seed(42)
ns_idx <- which(volc_dt$category == "Not significant")
keep_ns <- sample(ns_idx, min(5000, length(ns_idx)))
plot_dt <- rbind(
  volc_dt[category != "Not significant"],
  volc_dt[keep_ns]
)
setorder(plot_dt, category)  # NS plotted first (behind)

pc <- ggplot(plot_dt, aes(x = logFC, y = neg_log10_p, color = category)) +
  rasterize_layer(
    geom_point(size = 0.15, alpha = 0.5, shape = 16),
    dpi = 300
  ) +
  scale_color_manual(values = c(
    "Not significant" = masld_colors$ns,
    "Detectable DEG"  = masld_colors$deg,
    "Invisible DEG"   = masld_colors$up
  ), name = NULL) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed",
             linewidth = 0.2, color = "gray50") +
  labs(x = expression(log[2]~"fold change"), y = expression(-log[10]~"adjusted p"),
       title = "Volcano: microarray detectability") +
  theme_masld() +
  theme(legend.position = "top")

# --- Combine panels ---
fig <- (pa | pb) / pc +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

fig_path <- file.path(FIGDIR, "figS_microarray_gap.pdf")
save_fig(fig, fig_path, width = fig_full_width, height = 5.5)
cat(sprintf("  Figure saved: %s\n", fig_path))

# ============================================================
# 9. Print final summary for manuscript
# ============================================================
cat("\n============================================================\n")
cat("MANUSCRIPT SUMMARY\n")
cat("============================================================\n")
cat(sprintf("Microarray platform coverage:\n"))
cat(sprintf("  Affy HG-U133 Plus 2.0: %s genes (Ensembl)\n",
            format(length(affy_genes), big.mark = ",")))
cat(sprintf("  Illumina HumanHT-12 v4: %s genes (Ensembl)\n",
            format(length(ill_genes), big.mark = ",")))
cat(sprintf("  Union (either): %s genes (Ensembl)\n",
            format(length(union_genes), big.mark = ",")))
cat(sprintf("\nOf %s dream DEGs (padj < 0.1):\n",
            format(n_degs, big.mark = ",")))
cat(sprintf("  %s (%.1f%%) detectable on either platform\n",
            format(overall$detectable_any, big.mark = ","),
            100 * overall$detectable_any / overall$total_degs))
cat(sprintf("  %s (%.1f%%) invisible to microarray\n",
            format(overall$invisible, big.mark = ","),
            overall$pct_invisible))
cat(sprintf("\nlncRNA DEGs: %d total, %d (%.1f%%) invisible to microarray\n",
            nrow(lncrna_degs), sum(!lncrna_degs$detectable_either),
            100 * mean(!lncrna_degs$detectable_either)))
cat(sprintf("\nExpression bias: detectable DEGs mean AveExpr = %.2f; invisible = %.2f\n",
            mean(d_detectable), mean(d_invisible)))
cat(sprintf("  Wilcoxon p = %.1e, Cohen's d = %.2f\n", wtest$p.value, cohens_d))
cat("\nDone.\n")
