#!/usr/bin/env Rscript
# 57_gwas_atac_figures.R
# Enrichment tests, figure panels, and atlas integration for GWAS-ATAC overlap
# Requires: 55/56 outputs
# Usage: Rscript 57_gwas_atac_figures.R

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(patchwork)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")
FIG_DIR  <- file.path(BASE_DIR, "figures/supplementary/figS05_epigenomic_spatial")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
PANELS_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANELS_DIR, recursive = TRUE, showWarnings = FALSE)

# Load publication theme
theme_file <- file.path(BASE_DIR, "scripts/figures/publication_theme.R")
if (file.exists(theme_file)) {
  source(theme_file)
} else {
  theme_masld <- function(base_size = 7) theme_classic(base_size = base_size)
  save_fig <- function(p, f, w = 7, h = 5, ...) ggsave(f, p, width = w, height = h, device = cairo_pdf)
  fig_full_width <- 7.09
  fig_half_width <- 3.46
}

cat("============================================================\n")
cat("57_gwas_atac_figures.R\n")
cat("Enrichment tests, figures, and atlas integration\n")
cat("============================================================\n\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
ann_file <- file.path(OUT_DIR, "gwas_atac_variant_annotation.csv")
sum_file <- file.path(OUT_DIR, "variant_overlap_summary.csv")
mot_file <- file.path(OUT_DIR, "motif_disruption_scores.csv")

if (!file.exists(ann_file) || !file.exists(sum_file)) {
  stop("Script 55 outputs not found. Run 55_gwas_atac_variant_overlap.R first.")
}

ann     <- fread(ann_file)
# FIX (review F096/F097): fread reads empty cells as "" (not NA), so `!is.na("")` is TRUE and the
# gene-assignment case_when/ifelse blocks treat "" as a valid gene — collapsing 3,555/3,723 variants
# into one "" pseudo-gene (gene_level_gwas_atac.csv had 7 rows). Coerce "" -> NA on the gene columns.
for (.c in intersect(c("scenic_target_gene", "linked_gene", "nearest_gene"), names(ann)))
  ann[get(.c) == "", (.c) := NA_character_]
var_sum <- fread(sum_file)
cat("Loaded annotation:", nrow(ann), "variant-peak overlaps\n")
cat("Loaded variant summary:", nrow(var_sum), "unique variants\n")

if (file.exists(mot_file)) {
  motif <- fread(mot_file)
  for (.c in intersect(c("scenic_target_gene", "linked_gene", "nearest_gene"), names(motif)))
    motif[get(.c) == "", (.c) := NA_character_]  # FIX (review F096/F097): "" -> NA (same bug, motif joins)
  cat("Loaded motif disruption:", nrow(motif), "variant-motif pairs\n")
  has_motif <- TRUE
} else {
  cat("WARNING: Motif disruption file not found\n")
  motif <- data.table()
  has_motif <- FALSE
}

# ── 2. Enrichment tests ─────────────────────────────────────────────────────
cat("\n--- Enrichment analysis ---\n")

peak_dir <- file.path(ATAC_DIR, "results/label_transfer/cell_type_peak_sets_v2")
peak_files <- list.files(peak_dir, pattern = "_peaks\\.bed$", full.names = TRUE)
cell_types <- gsub("_peaks\\.bed$", "", basename(peak_files))
GENOME_SIZE <- 3.088e9

enrich_results <- list()
for (i in seq_along(peak_files)) {
  ct <- cell_types[i]
  bed <- fread(peak_files[i], header = FALSE)
  if (ncol(bed) == 3) colnames(bed) <- c("chr", "start", "end")
  else { bed <- bed[, 2:4]; colnames(bed) <- c("chr", "start", "end") }

  peak_coverage <- sum(bed$end - bed$start)
  bg_fraction   <- peak_coverage / GENOME_SIZE
  n_total   <- nrow(var_sum)
  n_overlap <- sum(grepl(ct, var_sum$cell_types_overlapping))
  obs_fraction <- n_overlap / n_total
  fold_enrich  <- obs_fraction / bg_fraction

  fisher_res <- fisher.test(matrix(c(
    n_overlap, n_total - n_overlap,
    round(bg_fraction * 1e6), round((1 - bg_fraction) * 1e6)
  ), nrow = 2, byrow = TRUE))

  enrich_results[[ct]] <- data.table(
    cell_type = ct, n_peaks = nrow(bed),
    peak_coverage_mb = round(peak_coverage / 1e6, 1),
    bg_fraction = bg_fraction, n_variants = n_total,
    n_overlapping = n_overlap, obs_fraction = obs_fraction,
    fold_enrichment = fold_enrich,
    fisher_p = fisher_res$p.value, fisher_or = fisher_res$estimate)
}

enrich_dt <- rbindlist(enrich_results)

# FIX (review A11#6 / L12): replace the anticonservative genome-coverage Fisher (1e6 pseudocount,
# which treats LD-clustered credible-set variants as independent draws) with the LD-block-matched
# permutation null from Script 56b (enrichment_ld_null.csv). The Fisher columns are kept as a
# deprecated diagnostic; the canonical Panel-A p-value and fold-enrichment now come from the null.
.ld_null_file <- file.path(OUT_DIR, "enrichment_ld_null.csv")
if (file.exists(.ld_null_file)) {
  .ldn <- fread(.ld_null_file)
  .norm <- function(x) tolower(gsub("[ _]+", "", x))
  enrich_dt[, .ctn := .norm(cell_type)]
  .ldn[, .ctn := .norm(cell_type)]
  enrich_dt <- merge(enrich_dt,
                     .ldn[, .(.ctn, ld_p_emp = p_emp, ld_p_emp_bh = p_emp_bh,
                              ld_log2FE = log2FE, ld_ci_lo = ci_lo, ld_ci_hi = ci_hi)],
                     by = ".ctn", all.x = TRUE)
  enrich_dt[, `:=`(fisher_p_deprecated = fisher_p, fisher_or_deprecated = fisher_or)]
  enrich_dt[!is.na(ld_p_emp), `:=`(fisher_p = ld_p_emp, fold_enrichment = 2^ld_log2FE)]
  enrich_dt[, .ctn := NULL]
  cat("  Panel A: enrichment p + fold-enrichment now sourced from the 56b LD-permutation null.\n")
} else {
  warning("enrichment_ld_null.csv not found — Panel A falls back to the anticonservative Fisher; run Script 56b first.", call. = FALSE)
}

enrich_dt$fisher_padj <- p.adjust(enrich_dt$fisher_p, method = "BH")
enrich_dt <- enrich_dt %>% arrange(fisher_p)
fwrite(enrich_dt, file.path(OUT_DIR, "enrichment_statistics.csv"))
cat("Enrichment statistics written\n")
print(as.data.frame(enrich_dt[, .(cell_type, n_overlapping, fold_enrichment, fisher_padj)]))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A: Cell-type chromatin landscape
# Cleveland dot plot — fold enrichment, significance, counts
# ══════════════════════════════════════════════════════════════════════════════
# FIX (pipeline bug 2026-05-31): wrap ALL figure panels so a cosmetic plotting failure
# (e.g. a duplicate ggplot factor level) can NEVER block the downstream atlas-integration data
# product. gene_level_gwas_atac.csv was stale at 7 rows for weeks because Panel B crashed before
# the atlas write at the bottom of this script. Figures are non-essential; the atlas CSV is not.
tryCatch({
cat("\n--- Panel A ---\n")

ep <- as.data.frame(enrich_dt)
ep <- ep[ep$n_overlapping > 0, ]
ep$sig <- ep$fisher_padj < 0.05
ep$cell_type <- reorder(factor(ep$cell_type), ep$fold_enrichment)
ep$padj_label <- ifelse(ep$fisher_padj < 1e-4,
  formatC(ep$fisher_padj, format = "e", digits = 1),
  ifelse(ep$fisher_padj < 0.05, paste0("p=", round(ep$fisher_padj, 3)), "ns"))

panel_a <- ggplot(ep, aes(x = fold_enrichment, y = cell_type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_segment(aes(xend = 1, yend = cell_type), color = "grey80", linewidth = 0.3) +
  geom_point(aes(size = n_overlapping, color = sig), shape = 16) +
  geom_text(aes(label = padj_label), hjust = -0.15, size = 1.9, color = "grey30") +
  geom_text(aes(label = paste0("n=", n_overlapping)), x = 0.97, hjust = 1,
            size = 1.9, color = "grey50") +
  scale_size_continuous(range = c(1.5, 7), name = "Variants",
                        breaks = c(25, 100, 200, 400)) +
  scale_color_manual(values = c("TRUE" = "#C2185B", "FALSE" = "#BDBDBD"), guide = "none") +
  scale_x_continuous(limits = c(0.9, max(ep$fold_enrichment) * 1.3),
                     breaks = seq(1, 3, 0.5)) +
  labs(x = "Fold enrichment vs genome background", y = NULL,
       title = "GWAS variant enrichment in cell-type ATAC peaks") +
  theme_masld(base_size = 7)

save_fig(panel_a, file.path(PANELS_DIR, "figS_gwas_atac_panel_a.pdf"),
         width = fig_half_width, height = 2.8)
cat("Panel A saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B: Top regulatory variants at MASLD loci
# Gene-centric — PIP on x-axis, shapes for disease regulon motif disruption
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Panel B ---\n")

hep_df <- as.data.frame(ann[ann$cell_type == "Hepatocytes", ])
hep_uniq <- hep_df[!duplicated(hep_df$variant_id), ]

# Motif disruption per variant
if (has_motif && nrow(motif) > 0) {
  mot_df <- as.data.frame(motif)
  mot_agg <- do.call(rbind, lapply(split(mot_df, mot_df$SNP_id), function(x) {
    reg <- x[x$motif_in_disease_regulon == TRUE, ]
    data.frame(SNP_id = x$SNP_id[1],
               n_tfs = length(unique(x$tf_name)),
               has_regulon = nrow(reg) > 0,
               top_regulon = if (nrow(reg) > 0) paste(unique(reg$tf_name), collapse = "/") else NA_character_,
               stringsAsFactors = FALSE)
  }))
  rownames(mot_agg) <- NULL
  hep_uniq <- merge(hep_uniq, mot_agg, by.x = "variant_id", by.y = "SNP_id", all.x = TRUE)
  hep_uniq$n_tfs[is.na(hep_uniq$n_tfs)] <- 0
  hep_uniq$has_regulon[is.na(hep_uniq$has_regulon)] <- FALSE
} else {
  hep_uniq$n_tfs <- 0; hep_uniq$has_regulon <- FALSE; hep_uniq$top_regulon <- NA
}

# Cell-type sharing
ct_per_var <- aggregate(cell_type ~ variant_id, data = as.data.frame(ann),
                         FUN = function(x) length(unique(x)))
colnames(ct_per_var)[2] <- "n_ct"
hep_uniq <- merge(hep_uniq, ct_per_var, by = "variant_id", all.x = TRUE)

# Clean gene labels
hep_uniq$gene <- hep_uniq$nearest_gene
hep_uniq$gene[grepl("^ENSG", hep_uniq$gene) | hep_uniq$gene == ""] <- NA
na_idx <- is.na(hep_uniq$gene)
hep_uniq$gene[na_idx] <- hep_uniq$linked_gene[na_idx]
hep_uniq$gene[grepl("^ENSG", hep_uniq$gene) | hep_uniq$gene == ""] <- NA

# Top 20 named variants by PIP
hep_named <- hep_uniq[!is.na(hep_uniq$gene), ]
hep_named <- hep_named[order(-hep_named$max_pip), ]
top <- head(hep_named, 20)
top <- top[order(top$max_pip), ]
# Make labels unique: gene name + short chr:pos for duplicates
dup_genes <- top$gene[duplicated(top$gene) | duplicated(top$gene, fromLast = TRUE)]
top$label <- ifelse(top$gene %in% dup_genes,
  paste0(top$gene, " (chr", top$chromosome, ":", round(top$position / 1e6, 2), "M)"),
  top$gene)
# FIX (pipeline bug 2026-05-31): two variants of the same gene rounding to the same 0.01-Mb
# position still collide → identical label → "factor level is duplicated" crash. make.unique()
# appends .1/.2 so every row gets a distinct level (and a distinct y-axis position).
top$label <- make.unique(as.character(top$label))
top$label <- factor(top$label, levels = top$label)
top$shape_cat <- ifelse(top$has_regulon, "Disrupts disease\nregulon motif", "No regulon hit")

# Right-margin TF annotation
top$tf_label <- ifelse(!is.na(top$top_regulon), top$top_regulon, "")

panel_b <- ggplot(top, aes(x = max_pip, y = label)) +
  geom_point(aes(color = n_ct, shape = shape_cat), size = 3, stroke = 0.4) +
  geom_text(aes(label = tf_label), x = max(top$max_pip) * 1.08, hjust = 0,
            size = 1.6, color = "#C2185B", fontface = "italic") +
  scale_color_gradient(low = "#BBDEFB", high = "#0D47A1",
                       name = "Cell types\nwith peak", breaks = c(1, 3, 5, 7)) +
  scale_shape_manual(values = c("Disrupts disease\nregulon motif" = 17, "No regulon hit" = 16),
                     name = NULL) +
  scale_x_continuous(limits = c(0, max(top$max_pip) * 1.35),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "Posterior inclusion probability (PIP)", y = NULL,
       title = "Top fine-mapped variants in hepatocyte open chromatin") +
  theme_masld(base_size = 7) +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(5, 35, 5, 5))

save_fig(panel_b, file.path(PANELS_DIR, "figS_gwas_atac_panel_b.pdf"),
         width = fig_half_width + 0.5, height = 4.5)
cat("Panel B saved\n")

# Save table
top_out <- top[, c("variant_id", "gene", "max_pip", "n_ct", "hep_da_logFC",
                     "hep_da_padj", "coloc_best_pp4", "n_tfs", "has_regulon", "top_regulon")]
top_out <- top_out[order(-top_out$max_pip), ]
fwrite(top_out, file.path(OUT_DIR, "top_regulatory_variants.csv"))
cat("Top regulatory variants:", nrow(top_out), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C: Disease regulon motif disruption
# Focused on the 12 SCENIC+ disease regulons disrupted by GWAS variants
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Panel C ---\n")

# Cross-modality disease master-regulator set (hepatocyte SCENIC+ regulon TF that is a
# well-powered bulk MASLD DEG or COLOC hit). Replaces the FDR-gated disease_regulons.csv,
# which is empty at the donor level (n=18 underpowered; the old non-zero set was
# cell-level pseudoreplication). Env-overridable to stay in lockstep with Script 56.
regulon_file <- Sys.getenv("REGULON_FILE",
  unset = file.path(ATAC_DIR, "scenic_plus/disease_master_regulators.csv"))
regulons <- if (file.exists(regulon_file)) fread(regulon_file) else data.table()

if (has_motif && nrow(motif) > 0) {
  reg_hits <- as.data.frame(motif[motif$motif_in_disease_regulon == TRUE, ])

  if (nrow(reg_hits) > 0) {
    # Variants per disease regulon TF
    rs <- aggregate(SNP_id ~ tf_name, data = reg_hits, FUN = function(x) length(unique(x)))
    colnames(rs)[2] <- "n_variants"
    rs <- rs[order(-rs$n_variants), ]

    # Join regulon activity
    if (nrow(regulons) > 0) {
      ra <- as.data.frame(regulons)[, c("tf_name", "regulon_activity_diff",
                                          "activity_padj", "n_target_genes")]
      rs <- merge(rs, ra, by = "tf_name", all.x = TRUE)
    }

    # Drug target annotation
    drug_tfs <- c("THRB", "NR1H4", "PPARA", "PPARG")
    rs$is_drug <- rs$tf_name %in% drug_tfs

    # Significance
    rs$stars <- ifelse(rs$activity_padj < 1e-10, "***",
                 ifelse(rs$activity_padj < 1e-5, "**",
                 ifelse(rs$activity_padj < 0.05, "*", "")))

    rs <- rs[order(rs$n_variants), ]
    rs$tf_name <- factor(rs$tf_name, levels = rs$tf_name)

    panel_c <- ggplot(rs, aes(x = n_variants, y = tf_name)) +
      geom_col(aes(fill = regulon_activity_diff), width = 0.65) +
      geom_text(aes(label = stars), hjust = -0.3, size = 2.5, color = "grey30") +
      geom_point(data = rs[rs$is_drug, ], aes(x = -0.4), shape = 18,
                 size = 2.5, color = "#00695C") +
      scale_fill_gradient2(low = "#0D47A1", mid = "#90CAF9", high = "#C2185B",
                           midpoint = 0, name = "Regulon\nactivity\n(disease)") +
      scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      labs(x = "GWAS variants disrupting TF motif", y = NULL,
           title = "Disease regulon motifs disrupted by GWAS variants",
           subtitle = expression(paste(diamond, " = FDA drug target TF; * = regulon padj"))) +
      theme_masld(base_size = 7) +
      theme(
        axis.text.y = element_text(
          face = ifelse(levels(rs$tf_name) %in% drug_tfs, "bold.italic", "plain")),
        plot.subtitle = element_text(size = 5, color = "grey40"))

    save_fig(panel_c, file.path(PANELS_DIR, "figS_gwas_atac_panel_c.pdf"),
             width = fig_half_width, height = 3.8)
    cat("Panel C saved:", nrow(rs), "regulon TFs\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# PANEL D: Regulatory chain examples
# Variant → ATAC peak → motif disruption → dysregulated regulon
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Panel D ---\n")

if (has_motif && nrow(motif) > 0 && nrow(regulons) > 0) {
  # Build chains: variant in hep peak + disease regulon motif disruption
  hep_v <- as.data.frame(ann[ann$cell_type == "Hepatocytes", ])
  hep_v <- hep_v[!duplicated(hep_v$variant_id), ]
  reg_v <- as.data.frame(motif[motif$motif_in_disease_regulon == TRUE, ])
  reg_r <- as.data.frame(regulons)

  chains <- merge(
    hep_v[, c("variant_id", "max_pip", "nearest_gene", "linked_gene",
              "peak_chr", "peak_start", "peak_end")],
    reg_v[, c("SNP_id", "tf_name", "alleleDiff")],
    by.x = "variant_id", by.y = "SNP_id")

  chains <- merge(chains,
    reg_r[, c("tf_name", "regulon_activity_diff", "activity_padj",
              "n_target_genes", "target_genes")],
    by = "tf_name", all.x = TRUE)

  # Clean gene name
  chains$gene <- chains$nearest_gene
  chains$gene[grepl("^ENSG", chains$gene)] <- chains$linked_gene[grepl("^ENSG", chains$gene)]
  chains$gene[chains$gene == ""] <- NA

  # Deduplicate: keep best variant per gene-TF pair
  chains <- chains[order(-chains$max_pip), ]
  chains$key <- paste(chains$gene, chains$tf_name)
  chains <- chains[!duplicated(chains$key), ]

  # Pick top 5
  top_chains <- head(chains[!is.na(chains$gene), ], 5)

  if (nrow(top_chains) > 0) {
    # Build strip data — one row per chain, columns as steps
    rows <- lapply(seq_len(nrow(top_chains)), function(i) {
      r <- top_chains[i, ]
      # Truncate target genes
      tg <- strsplit(r$target_genes, ";")[[1]]
      tg_str <- if (length(tg) <= 3) paste(tg, collapse = ", ")
                else paste(c(head(tg, 3), paste0("+", length(tg) - 3, " more")), collapse = ", ")
      data.frame(
        row = paste0(r$gene, " / ", r$tf_name),
        step = factor(c("Variant", "ATAC peak", "Motif disrupted", "Regulon"),
                       levels = c("Variant", "ATAC peak", "Motif disrupted", "Regulon")),
        label = c(
          paste0(r$variant_id, "\nPIP = ", round(r$max_pip, 3)),
          paste0(r$peak_chr, ":", formatC(r$peak_start, big.mark = ","), "-",
                 formatC(r$peak_end, big.mark = ",")),
          paste0(r$tf_name, " motif\n\u0394 = ", round(r$alleleDiff, 2)),
          paste0(r$tf_name, " regulon ", ifelse(r$regulon_activity_diff < 0, "\u2193", "\u2191"),
                 "\npadj = ", formatC(r$activity_padj, format = "e", digits = 0),
                 "\n", r$n_target_genes, " targets")),
        pip = r$max_pip,
        stringsAsFactors = FALSE)
    })
    strip_df <- do.call(rbind, rows)
    strip_df$row <- factor(strip_df$row, levels = rev(unique(strip_df$row)))

    panel_d <- ggplot(strip_df, aes(x = step, y = row)) +
      geom_tile(aes(fill = step), color = "white", linewidth = 1) +
      geom_text(aes(label = label), size = 1.7, lineheight = 0.9) +
      scale_fill_manual(values = c(
        "Variant"         = "#E3F2FD",
        "ATAC peak"       = "#BBDEFB",
        "Motif disrupted" = "#F8BBD0",
        "Regulon"         = "#FCE4EC"), guide = "none") +
      labs(x = NULL, y = NULL,
           title = "Regulatory variant-to-function chains",
           subtitle = "Variant \u2192 open chromatin \u2192 TF motif disruption \u2192 dysregulated regulon") +
      theme_masld(base_size = 7) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(face = "bold", size = 6),
            plot.subtitle = element_text(size = 5.5, color = "grey40", face = "italic"))

    save_fig(panel_d, file.path(PANELS_DIR, "figS_gwas_atac_panel_d.pdf"),
             width = fig_full_width, height = 3)
    cat("Panel D saved:", nrow(top_chains), "chains\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED FIGURE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Combining panels ---\n")

has_all <- exists("panel_a") && exists("panel_b") && exists("panel_c") && exists("panel_d")

if (has_all) {
  top_row    <- panel_a + panel_b + plot_layout(widths = c(1, 1.2))
  bottom_row <- panel_c + panel_d + plot_layout(widths = c(0.8, 1.2))
  combined   <- top_row / bottom_row +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "A")

  save_fig(combined, file.path(FIG_DIR, "figS_gwas_atac_combined.pdf"),
           width = fig_full_width, height = 7.5)
  cat("Combined figure saved\n")
} else {
  panels <- list()
  for (p in c("panel_a", "panel_b", "panel_c", "panel_d")) {
    if (exists(p)) panels[[length(panels) + 1]] <- get(p)
  }
  if (length(panels) >= 2) {
    combined <- wrap_plots(panels, ncol = 2) + plot_annotation(tag_levels = "A")
    save_fig(combined, file.path(FIG_DIR, "figS_gwas_atac_combined.pdf"),
             width = fig_full_width, height = 4 * ceiling(length(panels) / 2))
  }
  cat("Combined figure saved (fallback)\n")
}
}, error = function(e)
  cat("WARNING: figure generation failed (non-fatal) — atlas integration still runs:",
      conditionMessage(e), "\n"))

# ── 8. Atlas integration ────────────────────────────────────────────────────
cat("\n--- Atlas integration ---\n")

atlas_file <- file.path(BASE_DIR, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (!file.exists(atlas_file)) {
  cat("WARNING: Multi-evidence atlas not found — skipping atlas integration\n")
} else {
  atlas <- fread(atlas_file)
  cat("Loaded atlas:", nrow(atlas), "genes x", ncol(atlas), "columns\n")

  gene_col <- if ("gene_symbol" %in% colnames(atlas)) "gene_symbol" else
              if ("human_symbol" %in% colnames(atlas)) "human_symbol" else "gene"

  # FIX (audit bug 2): assigned_gene/nearest_gene carry GENCODE gene_name, which is an
  # ENSG ID when no HGNC symbol exists — so those genes never join to the atlas
  # human_symbol. Build an ENSG(base, version-stripped) -> HGNC symbol map from the
  # GENCODE v49 metadata and resolve any "^ENSG" gene labels that DO have a real symbol
  # before the join. Values that are already symbols pass through unchanged, and ENSGs
  # with no GENCODE symbol are left as-is (the metadata's gene_name == ENSG for them).
  # NOTE (verified 2026-06-07 on live data): all current ENSG-labelled assigned genes are
  # GENCODE symbol-less lncRNAs/pseudogenes, AND the 63 that don't already match the atlas
  # are absent from the atlas gene universe entirely — so on today's data this resolves 0
  # to a new symbol. The map is retained because it is the correct, robust behaviour and
  # will recover any future ENSG that does carry a real GENCODE symbol.
  meta_file <- file.path(BASE_DIR, "data/gencode_v49_gene_metadata.tsv.gz")
  resolve_ensg <- function(x) x  # identity fallback if metadata missing
  if (file.exists(meta_file)) {
    gmeta <- fread(meta_file)
    # ensembl_base is the version-stripped ENSG; gene_name is the HGNC symbol (or ENSG)
    ensg_base   <- if ("ensembl_base" %in% names(gmeta)) gmeta$ensembl_base
                   else sub("\\.\\d+$", "", gmeta$gene_id)
    sym_for_base <- gmeta$gene_name
    # Only keep entries whose gene_name is a real symbol (not itself an ENSG)
    keep_sym <- !grepl("^ENSG", sym_for_base) & !is.na(sym_for_base) & sym_for_base != ""
    ensg2sym <- setNames(sym_for_base[keep_sym], ensg_base[keep_sym])
    resolve_ensg <- function(x) {
      out <- x
      is_ensg <- !is.na(x) & grepl("^ENSG", x)
      base    <- sub("\\.\\d+$", "", x[is_ensg])
      mapped  <- ensg2sym[base]
      # only overwrite where a symbol was found; leave unmapped ENSG as-is
      hit <- !is.na(mapped)
      out[is_ensg][hit] <- mapped[hit]
      out
    }
    cat("Loaded GENCODE metadata:", length(ensg2sym), "ENSG->symbol mappings\n")
  } else {
    cat("WARNING: GENCODE metadata not found at", meta_file,
        "— ENSG gene labels will not be resolved\n")
  }

  ann_df <- as.data.frame(ann)
  gene_variant_map <- ann_df %>%
    mutate(
      assigned_gene = case_when(
        !is.na(scenic_target_gene) ~ scenic_target_gene,
        !is.na(linked_gene)        ~ linked_gene,
        !is.na(nearest_gene) & distance_to_tss <= 500000 ~ nearest_gene,
        TRUE ~ NA_character_
      ),
      assigned_gene = resolve_ensg(assigned_gene)
    ) %>%
    filter(!is.na(assigned_gene))

  # Report ENSG->symbol recovery against the atlas symbol column
  .atlas_syms <- unique(atlas[[gene_col]])
  .raw_assigned <- ann_df %>%
    mutate(raw_gene = case_when(
      !is.na(scenic_target_gene) ~ scenic_target_gene,
      !is.na(linked_gene)        ~ linked_gene,
      !is.na(nearest_gene) & distance_to_tss <= 500000 ~ nearest_gene,
      TRUE ~ NA_character_)) %>%
    filter(!is.na(raw_gene), grepl("^ENSG", raw_gene)) %>%
    distinct(raw_gene) %>% pull(raw_gene)
  .recovered <- sum(resolve_ensg(.raw_assigned) %in% .atlas_syms)
  cat("ENSG-labelled assigned genes:", length(.raw_assigned),
      "; recovered to atlas symbol after mapping:", .recovered, "\n")

  cat("Gene-to-variant mapping:", nrow(gene_variant_map), "entries\n")

  gvm <- gene_variant_map
  gene_gwas_atac <- do.call(rbind, lapply(split(gvm, gvm$assigned_gene), function(x) {
    data.frame(
      assigned_gene         = x$assigned_gene[1],
      gwas_variant_in_peak  = TRUE,
      gwas_variant_cell_types = paste(sort(unique(x$cell_type)), collapse = ","),
      gwas_variant_da_peak  = any(x$hep_da_padj < 0.05, na.rm = TRUE),
      gwas_n_variants       = length(unique(x$variant_id)),
      gwas_max_pip_in_peak  = max(x$max_pip, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  rownames(gene_gwas_atac) <- NULL

  if (has_motif && nrow(motif) > 0) {
    motif_gene_map <- motif %>%
      left_join({
          tmp <- ann_df
          tmp$assigned_gene <- ifelse(!is.na(tmp$scenic_target_gene), tmp$scenic_target_gene,
                               ifelse(!is.na(tmp$linked_gene), tmp$linked_gene,
                               ifelse(!is.na(tmp$nearest_gene) & tmp$distance_to_tss <= 500000,
                                      tmp$nearest_gene, NA_character_)))
          tmp$assigned_gene <- resolve_ensg(tmp$assigned_gene)  # FIX (audit bug 2): ENSG -> symbol
          tmp <- tmp[, c("variant_id", "assigned_gene")]
          tmp <- tmp[!duplicated(paste(tmp$variant_id, tmp$assigned_gene)), ]
          colnames(tmp)[1] <- "SNP_id"
          tmp
        }, by = "SNP_id") %>%
      filter(!is.na(assigned_gene))

    if (nrow(motif_gene_map) > 0) {
      motif_per_gene <- motif_gene_map %>%
        group_by(assigned_gene) %>%
        summarise(gwas_motif_disrupted = paste(sort(unique(tf_name)), collapse = ";"),
                  .groups = "drop")
      gene_gwas_atac <- gene_gwas_atac %>%
        left_join(motif_per_gene, by = "assigned_gene")
    } else {
      gene_gwas_atac$gwas_motif_disrupted <- NA_character_
    }
  } else {
    gene_gwas_atac$gwas_motif_disrupted <- NA_character_
  }

  da_per_gene <- gene_variant_map %>%
    filter(!is.na(hep_da_logFC)) %>%
    group_by(assigned_gene) %>%
    summarise(max_da_abs = max(abs(hep_da_logFC), na.rm = TRUE), .groups = "drop")

  motif_per_gene_score <- if (has_motif && nrow(motif) > 0) {
    motif %>%
      left_join({
          tmp <- ann_df
          tmp$assigned_gene <- ifelse(!is.na(tmp$scenic_target_gene), tmp$scenic_target_gene,
                               ifelse(!is.na(tmp$linked_gene), tmp$linked_gene,
                               ifelse(!is.na(tmp$nearest_gene) & tmp$distance_to_tss <= 500000,
                                      tmp$nearest_gene, NA_character_)))
          tmp$assigned_gene <- resolve_ensg(tmp$assigned_gene)  # FIX (audit bug 2): ENSG -> symbol
          tmp <- tmp[, c("variant_id", "assigned_gene")]
          tmp <- tmp[!duplicated(paste(tmp$variant_id, tmp$assigned_gene)), ]
          colnames(tmp)[1] <- "SNP_id"
          tmp
        }, by = "SNP_id") %>%
      filter(!is.na(assigned_gene)) %>%
      group_by(assigned_gene) %>%
      summarise(max_motif_diff = max(abs(alleleDiff), na.rm = TRUE), .groups = "drop")
  } else {
    data.table(assigned_gene = character(0), max_motif_diff = numeric(0))
  }

  gene_gwas_atac <- gene_gwas_atac %>%
    left_join(da_per_gene, by = "assigned_gene") %>%
    left_join(motif_per_gene_score, by = "assigned_gene") %>%
    mutate(
      pip_norm   = pmin(gwas_max_pip_in_peak, 1),
      da_norm    = pmin(replace_na(max_da_abs, 0) / 2, 1),
      motif_norm = pmin(replace_na(max_motif_diff, 0) / 1, 1),
      gwas_atac_regulatory_score = pip_norm * (0.4 + 0.3 * da_norm + 0.3 * motif_norm)
    ) %>%
    select(-pip_norm, -da_norm, -motif_norm, -max_da_abs, -max_motif_diff)

  cat("Gene-level GWAS-ATAC annotations for", nrow(gene_gwas_atac), "genes\n")

  drop_cols <- grep("^gwas_(variant_|motif_|atac_|n_variants|max_pip_in)", colnames(atlas), value = TRUE)
  if (length(drop_cols) > 0) {
    cat("Removing", length(drop_cols), "existing GWAS-ATAC columns (re-run)\n")
    atlas <- atlas[, !..drop_cols]
  }

  atlas_updated <- atlas %>%
    left_join(gene_gwas_atac, by = setNames("assigned_gene", gene_col))
  atlas_updated <- atlas_updated %>%
    mutate(
      gwas_variant_in_peak = replace_na(gwas_variant_in_peak, FALSE),
      gwas_variant_da_peak = replace_na(gwas_variant_da_peak, FALSE),
      gwas_n_variants      = replace_na(gwas_n_variants, 0L))

  fwrite(atlas_updated, atlas_file)
  cat("Updated atlas:", nrow(atlas_updated), "x", ncol(atlas_updated), "columns\n")

  n_with <- sum(atlas_updated$gwas_variant_in_peak)
  cat("  Genes with GWAS variant in peak:", n_with, "\n")
  fwrite(gene_gwas_atac, file.path(OUT_DIR, "gene_level_gwas_atac.csv"))
}

cat("\n============================================================\n")
cat("Script 57 complete\n")
cat("Figures in:", FIG_DIR, "\n")
cat("Results in:", OUT_DIR, "\n")
cat("============================================================\n")
