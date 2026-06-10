#!/usr/bin/env Rscript
# 206d_dream_vs_published.R
# Compare our integrated dream mega-analysis DEGs directly against
# published DEGs from Hoang 2019 and Govaere 2020.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
PUB  <- file.path(BASE, "data/published_degs")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/sensitivity")

cat("=== 206d: Dream vs Published DEG Comparison ===\n")

# ── Gene mapping ────────────────────────────────────────────────────────────
annot <- fread(file.path(INT, "gene_annotation/human_ensg_to_symbol.tsv"))
annot_pc <- annot[gene_type == "protein_coding"]
sym2ens <- setNames(annot_pc$gene_base, annot_pc$symbol)
extra <- annot[!gene_base %in% annot_pc$gene_base & !symbol %in% names(sym2ens)]
sym2ens <- c(sym2ens, setNames(extra[!duplicated(symbol)]$gene_base,
                                extra[!duplicated(symbol)]$symbol))
sym2ens_upper <- setNames(sym2ens, toupper(names(sym2ens)))
map_sym <- function(s) {
  s <- as.character(s); m <- sym2ens[s]
  miss <- is.na(m); if (any(miss)) m[miss] <- sym2ens_upper[toupper(s[miss])]; m
}

# ── Load dream results ──────────────────────────────────────────────────────
cat("Loading dream results...\n")
dream <- fread(file.path(INT, "integration/dream_results.csv"))
dream[, gene_base := sub("\\.\\d+$", "", gene)]

dream_fib <- fread(file.path(INT, "../results/progression/c17_fibrosis_ordinal_dream.csv"))
dream_fib[, gene_base := sub("\\.\\d+$", "", gene)]

# Also load NAFL vs NASH dream and adv vs early fibrosis dream
nafl_nash_dream <- fread(file.path(INT, "disease_signatures/nafl_vs_nash_dream.csv"))
nafl_nash_dream[, gene_base := sub("\\.\\d+$", "", gene)]

adv_fib_dream <- fread(file.path(INT, "disease_signatures/adv_vs_early_fibrosis_dream.csv"))
adv_fib_dream[, gene_base := sub("\\.\\d+$", "", gene)]

DPADJ <- 0.1
dream_degs <- dream[padj < DPADJ, gene_base]
dream_fib_degs <- dream_fib[padj < DPADJ, gene_base]

# ── Load Hoang published ────────────────────────────────────────────────────
cat("Loading Hoang published DEGs...\n")
std_h <- function(dt) {
  setnames(dt, names(dt), tolower(gsub("\\s+", "_", names(dt))))
  if ("gene_symbol" %in% names(dt)) setnames(dt, "gene_symbol", "symbol", skip_absent = TRUE)
  if ("adj_p" %in% names(dt)) setnames(dt, "adj_p", "adj_p_val", skip_absent = TRUE)
  dt[, gene_base := map_sym(symbol)]; dt[!is.na(gene_base)]
}
hoang_nas <- std_h(as.data.table(read_excel(
  file.path(PUB, "GSE130970/MOESM2.xlsx"), sheet = "NAS ordinal regression")))
hoang_fib <- std_h(as.data.table(read_excel(
  file.path(PUB, "GSE130970/MOESM2.xlsx"), sheet = "fibrosis ordinal regression")))

hoang_nas_degs <- hoang_nas[adj_p_val < 0.01, gene_base]
hoang_fib_degs <- hoang_fib[adj_p_val < 0.01, gene_base]

# ── Load Govaere published ──────────────────────────────────────────────────
cat("Loading Govaere supplementary tables...\n")
govaere_xlsx <- file.path(PUB, "GSE135251/aba4448_supplementary_tables.xlsx")
supp_info <- data.table(
  sheet = paste("Table", c("S1", "S3", "S4", "S5", "S6", "S7", "S8", "S9")),
  contrast = c("Cluster_A_vs_B", "NASH_F2_vs_NAFL", "NASH_F3_vs_NAFL",
               "NASH_F4_vs_NAFL", "NASH_F3_vs_F01", "NASH_F4_vs_F01",
               "NAS_ge4", "SAF_ge2"),
  label = c("Cluster A vs B", "NASH F2 vs NAFL", "NASH F3 vs NAFL",
            "NASH F4 vs NAFL", "NASH F3 vs F0/1", "NASH F4 vs F0/1",
            "NAS >= 4", "SAF >= 2"),
  expected_n = c(250, 50, 907, 1369, 434, 1194, 369, 320)
)
avail <- excel_sheets(govaere_xlsx)
gov_pub <- list()
for (i in seq_len(nrow(supp_info))) {
  sn <- supp_info$sheet[i]
  if (!sn %in% avail) next
  raw <- as.data.table(read_excel(govaere_xlsx, sheet = sn, skip = 1))
  setnames(raw, names(raw), tolower(gsub("\\s+", "_", names(raw))))
  nms <- names(raw)
  ens_col <- grep("ensembl|gene_id", nms, value = TRUE)[1]
  lfc_col <- grep("log2?fc|logfc", nms, value = TRUE)[1]
  fdr_col <- grep("qvalue|q_value|fdr", nms, value = TRUE)[1]
  if (!is.na(ens_col)) raw[, gene_base := sub("\\.\\d+$", "", as.character(raw[[ens_col]]))]
  if (!is.na(lfc_col)) raw[, pub_logFC := as.numeric(raw[[lfc_col]])]
  if (!is.na(fdr_col)) raw[, pub_fdr := as.numeric(raw[[fdr_col]])]
  gov_pub[[supp_info$contrast[i]]] <- list(
    data = raw[!is.na(gene_base)], label = supp_info$label[i],
    has_lfc = !is.na(lfc_col))
}

# ── Helpers ─────────────────────────────────────────────────────────────────
scatter <- function(data, xlab, ylab, title) {
  r <- cor(data$lfc1, data$lfc2)
  rho <- cor(data$lfc1, data$lfc2, method = "spearman")
  dir <- mean(sign(data$lfc1) == sign(data$lfc2)) * 100
  anno <- sprintf("r = %.3f\nrho = %.3f\ndir = %.1f%%\nn = %s",
                  r, rho, dir, formatC(nrow(data), big.mark = ","))
  ggplot(data, aes(x = lfc1, y = lfc2)) +
    geom_point(alpha = 0.12, size = 0.2, color = masld_colors$ns) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = masld_colors$up) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dotted", color = "grey40") +
    annotate("text", x = -Inf, y = Inf, label = anno, hjust = -0.05, vjust = 1.2,
             size = 1.6, color = "grey30") +
    labs(x = xlab, y = ylab, title = title) +
    theme_masld(base_size = 6)
}

lfc_merge <- function(dt1, dt2, l1, l2, id = "gene_base") {
  m <- merge(dt1[, c(id, l1), with = FALSE], dt2[, c(id, l2), with = FALSE], by = id)
  setnames(m, c(l1, l2), c("lfc1", "lfc2"))
  m[is.finite(lfc1) & is.finite(lfc2)]
}

# ══════════════════════════════════════════════════════════════════════════
# BUILD PANELS
# ══════════════════════════════════════════════════════════════════════════
cat("\nBuilding panels...\n")
panels <- list()
metrics <- list()

# --- Hoang vs Dream ---
# H1: NAS ordinal vs Dream disease-vs-ctrl
m_h1 <- lfc_merge(hoang_nas, dream, "range_log2fc", "logFC")
panels[["h_nas_dvc"]] <- scatter(m_h1, "Hoang NAS ordinal (range log2FC)",
                                  "Dream disease vs ctrl (logFC)",
                                  "Hoang NAS vs Dream (disease vs ctrl)")

# H2: Fib ordinal vs Dream fibrosis ordinal
m_h2 <- lfc_merge(hoang_fib, dream_fib, "range_log2fc", "logFC")
panels[["h_fib_dfo"]] <- scatter(m_h2, "Hoang Fib ordinal (range log2FC)",
                                  "Dream fib ordinal (logFC)",
                                  "Hoang Fib vs Dream (fib ordinal)")

# H3: NAS ordinal vs Dream fib ordinal (cross-phenotype)
m_h3 <- lfc_merge(hoang_nas, dream_fib, "range_log2fc", "logFC")
panels[["h_nas_dfo"]] <- scatter(m_h3, "Hoang NAS ordinal (range log2FC)",
                                  "Dream fib ordinal (logFC)",
                                  "Hoang NAS vs Dream (fib ordinal)")

for (nm in c("h_nas_dvc", "h_fib_dfo", "h_nas_dfo")) {
  d <- panels[[nm]]$data
  if (is.null(d)) d <- switch(nm, h_nas_dvc=m_h1, h_fib_dfo=m_h2, h_nas_dfo=m_h3)
  metrics[[nm]] <- data.table(
    study = "Hoang", contrast = nm,
    r = cor(d$lfc1, d$lfc2), rho = cor(d$lfc1, d$lfc2, method = "spearman"),
    dir = mean(sign(d$lfc1) == sign(d$lfc2)) * 100, n = nrow(d))
}

# --- Govaere vs Dream ---
# For each Govaere contrast, scatter their published logFC vs dream logFC
for (cn in names(gov_pub)) {
  info <- gov_pub[[cn]]
  if (!info$has_lfc) next
  m <- lfc_merge(info$data, dream, "pub_logFC", "logFC")
  panels[[paste0("g_", cn, "_dvc")]] <- scatter(
    m, paste("Govaere", info$label, "(logFC)"),
    "Dream disease vs ctrl (logFC)",
    paste(info$label, "vs Dream"))
  metrics[[paste0("g_", cn)]] <- data.table(
    study = "Govaere", contrast = info$label,
    r = cor(m$lfc1, m$lfc2), rho = cor(m$lfc1, m$lfc2, method = "spearman"),
    dir = mean(sign(m$lfc1) == sign(m$lfc2)) * 100, n = nrow(m))
  cat(sprintf("  %s vs Dream: r=%.3f, dir=%.1f%%, n=%d\n",
              info$label, cor(m$lfc1, m$lfc2),
              mean(sign(m$lfc1) == sign(m$lfc2)) * 100, nrow(m)))
}

# --- Recovery bar chart ---
cat("\nComputing recovery rates...\n")

# For Hoang: what fraction of their DEGs appear in dream DEGs
rec_data <- data.table(
  published = character(), recovery = numeric(), dream_type = character()
)
rec_data <- rbind(rec_data, data.table(
  published = "Hoang NAS\n(2,970)", dream_type = "Disease vs ctrl",
  recovery = mean(hoang_nas_degs %in% dream_degs) * 100))
rec_data <- rbind(rec_data, data.table(
  published = "Hoang NAS\n(2,970)", dream_type = "Fib ordinal",
  recovery = mean(hoang_nas_degs %in% dream_fib_degs) * 100))
rec_data <- rbind(rec_data, data.table(
  published = "Hoang Fib\n(1,656)", dream_type = "Disease vs ctrl",
  recovery = mean(hoang_fib_degs %in% dream_degs) * 100))
rec_data <- rbind(rec_data, data.table(
  published = "Hoang Fib\n(1,656)", dream_type = "Fib ordinal",
  recovery = mean(hoang_fib_degs %in% dream_fib_degs) * 100))

# For Govaere: each contrast's DEGs recovered by dream
for (cn in names(gov_pub)) {
  info <- gov_pub[[cn]]
  pub_genes <- info$data$gene_base
  rec <- mean(pub_genes %in% dream_degs) * 100
  lab <- paste0("Govaere\n", info$label, "\n(", length(pub_genes), ")")
  rec_data <- rbind(rec_data, data.table(
    published = lab, dream_type = "Disease vs ctrl", recovery = rec))
}
cat("  Recovery rates computed\n")

# Recovery plot - Hoang
rec_hoang <- rec_data[grepl("Hoang", published)]
rec_hoang[, dream_type := factor(dream_type, levels = c("Disease vs ctrl", "Fib ordinal"))]
p_rec_h <- ggplot(rec_hoang, aes(x = published, y = recovery, fill = dream_type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", recovery)),
            position = position_dodge(0.7), vjust = -0.3, size = 1.6) +
  scale_fill_manual(values = c("Disease vs ctrl" = masld_colors$up,
                               "Fib ordinal" = masld_colors$down), name = "Dream contrast") +
  labs(x = NULL, y = "Recovery (%)", title = "Hoang DEGs in Dream") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_masld(base_size = 6) + theme(legend.key.size = unit(0.3, "cm"))

# Recovery plot - Govaere
rec_gov <- rec_data[grepl("Govaere", published)]
rec_gov[, published := gsub("Govaere\n", "", published)]
p_rec_g <- ggplot(rec_gov, aes(x = reorder(published, -recovery), y = recovery)) +
  geom_col(fill = masld_colors$up, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", recovery)), vjust = -0.3, size = 1.6) +
  labs(x = NULL, y = "Recovery (%)", title = "Govaere DEGs in Dream") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4))

# --- Correlation summary dot plot ---
met_dt <- rbindlist(metrics)
met_dt[, contrast := factor(contrast, levels = rev(contrast))]
p_corsum <- ggplot(met_dt, aes(x = r, y = contrast, color = study)) +
  geom_point(aes(size = n)) +
  geom_text(aes(label = sprintf("%.3f", r)), hjust = -0.3, size = 1.8) +
  scale_color_manual(values = c(Hoang = masld_colors$down, Govaere = masld_colors$up), name = NULL) +
  scale_size_continuous(range = c(1, 3), name = "Genes\npaired") +
  labs(x = "Pearson r (published vs Dream logFC)", y = NULL,
       title = "Published vs Dream correlation") +
  xlim(0, 1.1) +
  theme_masld(base_size = 6) + theme(legend.key.size = unit(0.3, "cm"))

# ══════════════════════════════════════════════════════════════════════════
# ASSEMBLE FIGURE
# ══════════════════════════════════════════════════════════════════════════
cat("\nAssembling figure...\n")

# Row 1: Hoang vs Dream (3 scatter + recovery)
row1 <- panels[["h_nas_dvc"]] | panels[["h_fib_dfo"]] | panels[["h_nas_dfo"]] | p_rec_h

# Row 2: Govaere vs Dream - NAFL contrasts (3 scatter)
row2 <- panels[["g_NASH_F2_vs_NAFL_dvc"]] |
        panels[["g_NASH_F3_vs_NAFL_dvc"]] |
        panels[["g_NASH_F4_vs_NAFL_dvc"]]

# Row 3: Govaere vs Dream - F0/1 contrasts + NAS (3 scatter)
row3 <- panels[["g_NASH_F3_vs_F01_dvc"]] |
        panels[["g_NASH_F4_vs_F01_dvc"]] |
        panels[["g_NAS_ge4_dvc"]]

# Row 4: Summary panels
row4 <- p_rec_g | p_corsum | plot_spacer()

fig <- (row1 / row2 / row3 / row4) +
  plot_annotation(tag_levels = "a",
                  title = "Published DEGs vs Dream mega-analysis") &
  theme(plot.tag = element_text(size = 7, face = "bold"))

save_fig(fig, file.path(OUTDIR, "dream_vs_published.pdf"),
         width = fig_full_width, height = 10)
cat(sprintf("  Saved %s\n", file.path(OUTDIR, "dream_vs_published.pdf")))

# Save metrics
fwrite(met_dt, file.path(OUTDIR, "dream_vs_published_metrics.csv"))
fwrite(rec_data, file.path(OUTDIR, "dream_vs_published_recovery.csv"))
cat("  Saved metrics and recovery tables\n")

cat("\n=== Done ===\n")
