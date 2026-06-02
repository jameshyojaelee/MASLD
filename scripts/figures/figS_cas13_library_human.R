#!/usr/bin/env Rscript
# =============================================================================
# Figure S_lib_5 -- Human evidence integration for Cas13 library
# KEY MESSAGE: Mouse DEG backbone is strongly enriched for human transcriptomic
#   and genetic evidence, with concordant effect directions and COLOC validation
#   across multiple diet models.
#
# Panels:
#   A  Alluvial: diet replication -> human evidence -> library tier
#   B  Scatter: mouse logFC vs human dream logFC (colored by n_diets, shaped by
#      COLOC, alpha by library membership)
#   C  Dot plot: per-diet enrichment for HumanDE / COLOC / Finemap
#
# Output: figures/supplementary/figS_cas13_library/S_lib_5_human_evidence.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(ggrastr)
  library(ggrepel)
  library(viridis)
  library(patchwork)
})

# -- Project paths & theme ----------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -- Constants ----------------------------------------------------------------
DIETS <- c("MCD", "CDAHFD", "Western", "HFD")  # NASH+Western merged 2026-05-29 (4 groups)
DIET_LABELS <- c(MCD = "MCD", CDAHFD = "CDAHFD", Western = "Western", HFD = "HFD")
PERDIET_DIR <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
ATLAS_PATH  <- file.path(ME, "multi_evidence_atlas.csv")
ORTHO_TABLE <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
COLOC_PP4_THR <- 0.5
PADJ_THR <- 0.05
LFC_THR  <- 0.5   # per-diet DE threshold (mouse)
# Mouse UP-DEG definition: ashr-shrunk effect size (M02c), power-fair across
# diets of different N. lfsr<0.05 & shrunk_logFC>0.5 (matches human Tier-1).
LFSR_THR       <- 0.05
SHRUNK_LFC_THR <- 0.5

# Tier colors
tier_pal <- c(Core       = "#880E4F",
              Extended   = "#E91E63",
              Exploratory = "#42A5F5")

# Human evidence state colors
evidence_pal <- c("HumanDE+COLOC" = "#00695C",
                  "HumanDE only"  = "#C9265E",
                  "COLOC only"    = "#7B1FA2",
                  "Neither"       = "#9E9E9E")

# =============================================================================
# 1. Load per-diet mouse DE, compute per-gene diet replication stats
# =============================================================================
message("Loading per-diet mouse DE ...")
perdiet_list <- lapply(DIETS, function(d) {
  f <- file.path(PERDIET_DIR, paste0(d, "_de_results.csv"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  dt <- fread(f)
  dt[, gene_id_base := sub("\\..*", "", gene)]
  dt[, diet := d]
  dt
})
perdiet <- rbindlist(perdiet_list[!sapply(perdiet_list, is.null)], fill = TRUE)

# Significant UP genes per diet: ashr-shrunk effect size (lfsr<0.05 &
# shrunk_logFC>0.5), power-fair across diets of different sample size.
perdiet[, sig_up := lfsr < LFSR_THR & shrunk_logFC > SHRUNK_LFC_THR]

# Per-gene: n_diets_up and max (shrunk) logFC across the 5 diets
gene_stats <- perdiet[sig_up == TRUE, .(
  n_diets_up = uniqueN(diet),
  max_logFC  = max(shrunk_logFC, na.rm = TRUE),
  diet_list  = paste(sort(unique(diet)), collapse = ",")
), by = gene_id_base]

# Also get max (shrunk) logFC for ALL genes (not just sig) for Panel B scatter
gene_max_lfc <- perdiet[, .(
  max_logFC_any = max(shrunk_logFC, na.rm = TRUE)
), by = gene_id_base]

message(sprintf("  %d genes UP in >= 1 diet (6-diet basis)", nrow(gene_stats)))

# =============================================================================
# 2. Load ortholog bridge (master table, min_tier M)
# =============================================================================
message("Loading ortholog bridge ...")
if (file.exists(ORTHO_TABLE)) {
  ortho_raw <- fread(ORTHO_TABLE, select = c("mouse_ensembl", "human_ensembl",
                                              "mouse_symbol", "human_symbol",
                                              "confidence_tier"))
  ortho_raw[, mouse_ensembl := sub("\\..*", "", mouse_ensembl)]
  ortho_raw[, human_ensembl := sub("\\..*", "", human_ensembl)]
  # Keep M and H tier
  ortho <- ortho_raw[confidence_tier %in% c("H", "M")]
  # Deduplicate: keep one (mouse, human) pair, prefer H
  ortho <- ortho[order(confidence_tier)]  # H before M
  ortho <- unique(ortho, by = c("mouse_ensembl", "human_ensembl"))
  message(sprintf("  Ortholog bridge: %d pairs (H+M tier)", nrow(ortho)))
} else {
  stop("Master ortholog table not found: ", ORTHO_TABLE)
}

# =============================================================================
# 3. Load human atlas
# =============================================================================
message("Loading human atlas ...")
atlas <- fread(ATLAS_PATH, select = c("human_symbol", "ensembl_id",
                                       "dream_logFC", "dream_padj",
                                       "coloc_best_susie_pp4_polyfun",
                                       "susiex_max_pip", "ctwas_pip"))
atlas[, human_ensembl := sub("\\..*", "", ensembl_id)]

# ashr human effect sizes -- library human-arm basis (lfsr<0.05 & shrunk_logFC>0.2).
# The library human arm uses ashr (not raw dream logFC), so this figure's
# human-evidence classification must use the same instrument.
ASHR_PATH <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                       "results/integration/dream_results_ashr.csv")
ashr_h <- fread(ASHR_PATH, select = c("gene", "shrunk_logFC", "lfsr"))
ashr_h[, human_ensembl := sub("\\..*", "", gene)]
setnames(ashr_h, c("shrunk_logFC", "lfsr"), c("human_shrunk_logFC", "human_lfsr"))
ashr_h <- unique(ashr_h, by = "human_ensembl")
atlas <- merge(atlas, ashr_h[, .(human_ensembl, human_shrunk_logFC, human_lfsr)],
               by = "human_ensembl", all.x = TRUE)

# Human DE flag -- ashr UP, matching the library human arm (shrunk_logFC>0.2)
atlas[, is_human_de := !is.na(human_lfsr) & human_lfsr < 0.05 &
        !is.na(human_shrunk_logFC) & human_shrunk_logFC > 0.2]

# COLOC flag
atlas[, is_coloc := !is.na(coloc_best_susie_pp4_polyfun) &
        coloc_best_susie_pp4_polyfun > COLOC_PP4_THR]

# Finemap flag (any PIP > 0.5 from SuSiEx or cTWAS)
atlas[, is_finemap := (!is.na(susiex_max_pip) & susiex_max_pip > 0.5) |
        (!is.na(ctwas_pip) & ctwas_pip > 0.5)]

message(sprintf("  Human DE: %d genes; COLOC: %d; Finemap: %d",
                sum(atlas$is_human_de, na.rm = TRUE),
                sum(atlas$is_coloc, na.rm = TRUE),
                sum(atlas$is_finemap, na.rm = TRUE)))

# =============================================================================
# 4. Bridge mouse genes to human
# =============================================================================
message("Bridging mouse -> human ...")
# Join gene_stats (mouse sig UP) with ortho bridge
bridged <- merge(gene_stats, ortho, by.x = "gene_id_base", by.y = "mouse_ensembl",
                 allow.cartesian = TRUE)
# Join with atlas
bridged <- merge(bridged, atlas, by = "human_ensembl", all.x = TRUE,
                 suffixes = c("_ortho", "_atlas"))
# Coalesce human_symbol from ortho (preferred) and atlas
if ("human_symbol_ortho" %in% names(bridged)) {
  bridged[, human_symbol := fifelse(!is.na(human_symbol_ortho) & human_symbol_ortho != "",
                                     human_symbol_ortho, human_symbol_atlas)]
} else if ("human_symbol" %in% names(bridged)) {
  # no conflict — do nothing
} else if ("human_symbol_atlas" %in% names(bridged)) {
  bridged[, human_symbol := human_symbol_atlas]
}

# De-dup: keep one row per mouse gene (take best evidence)
bridged[, has_any_human := !is.na(dream_logFC)]
bridged <- bridged[order(-is_coloc, -is_human_de, -abs(human_shrunk_logFC))]
bridged <- bridged[!duplicated(gene_id_base)]

message(sprintf("  Bridged: %d mouse genes -> human", nrow(bridged)))

# Classify human evidence state
bridged[, evidence_state := fifelse(
  is_human_de & is_coloc, "HumanDE+COLOC",
  fifelse(is_human_de & !is_coloc, "HumanDE only",
          fifelse(!is_human_de & is_coloc, "COLOC only", "Neither"))
)]

# Library tier based on diet replication
bridged[, library_tier := fifelse(
  n_diets_up >= 3, "Core",
  fifelse(n_diets_up == 2, "Extended", "Exploratory")
)]

# Diet replication bin (for alluvial left axis)
bridged[, diet_bin := fifelse(
  n_diets_up >= 5, "5+ diets",
  fifelse(n_diets_up == 4, "4 diets",
          fifelse(n_diets_up == 3, "3 diets",
                  fifelse(n_diets_up == 2, "2 diets", "1 diet")))
)]

# Include un-bridged mouse genes (those without human ortholog or human data)
# that were sig UP in at least 1 diet
unbridged_ids <- setdiff(gene_stats$gene_id_base, bridged$gene_id_base)
if (length(unbridged_ids) > 0) {
  unbridged <- gene_stats[gene_id_base %in% unbridged_ids]
  unbridged[, evidence_state := "Neither"]
  unbridged[, library_tier := fifelse(n_diets_up >= 3, "Core",
                                       fifelse(n_diets_up == 2, "Extended",
                                               "Exploratory"))]
  unbridged[, diet_bin := fifelse(
    n_diets_up >= 5, "5+ diets",
    fifelse(n_diets_up == 4, "4 diets",
            fifelse(n_diets_up == 3, "3 diets",
                    fifelse(n_diets_up == 2, "2 diets", "1 diet")))
  )]
  # Fill missing columns so rbind works
  for (col in setdiff(names(bridged), names(unbridged))) {
    unbridged[, (col) := NA]
  }
  bridged <- rbind(bridged, unbridged[, names(bridged), with = FALSE], fill = TRUE)
}

message(sprintf("  Total genes for panels: %d", nrow(bridged)))

# =============================================================================
# Panel A -- Alluvial: diet replication -> human evidence -> library tier
# =============================================================================
message("Generating Panel A (alluvial) ...")

# Build frequency table for ggalluvial
alluv_data <- bridged[, .N, by = .(diet_bin, evidence_state, library_tier)]

# Order factors
alluv_data[, diet_bin := factor(diet_bin,
  levels = c("4 diets", "3 diets", "2 diets", "1 diet"))]
alluv_data[, evidence_state := factor(evidence_state,
  levels = c("HumanDE+COLOC", "HumanDE only", "COLOC only", "Neither"))]
alluv_data[, library_tier := factor(library_tier,
  levels = c("Core", "Extended", "Exploratory"))]

pA <- ggplot(alluv_data,
       aes(y = N,
           axis1 = diet_bin,
           axis2 = evidence_state,
           axis3 = library_tier)) +
  geom_alluvium(aes(fill = library_tier), width = 1/6, alpha = 0.75,
                decreasing = FALSE) +
  geom_stratum(width = 1/4, fill = "grey92", color = "grey40", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 1.8, lineheight = 0.9) +
  scale_x_discrete(limits = c("Diet replication", "Human evidence", "Library tier"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = tier_pal, name = "Library tier") +
  labs(y = "Number of genes",
       title = "A") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 6, face = "bold"),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 6))

# =============================================================================
# Panel B -- Scatter: mouse max logFC vs human dream logFC
# =============================================================================
message("Generating Panel B (scatter) ...")

# Scatter data: bridged genes that have both mouse and human (ashr) shrunk logFC
scatter_dt <- bridged[!is.na(human_shrunk_logFC) & !is.na(max_logFC)]
scatter_dt[, is_library := library_tier %in% c("Core", "Extended", "Exploratory")]
scatter_dt[, shape_coloc := fifelse(is_coloc, "COLOC PP4>0.5", "No COLOC")]

# Identify label candidates: COLOC genes with concordant direction, sorted by
# combined evidence strength
scatter_dt[, concordant := sign(max_logFC) == sign(human_shrunk_logFC)]
label_cands <- scatter_dt[is_coloc == TRUE & concordant == TRUE][
  order(-abs(human_shrunk_logFC) * abs(max_logFC))
]

# Pick top 8 for labeling, prefer those with human_symbol
label_cands <- label_cands[!is.na(human_symbol) & human_symbol != ""]
label_genes <- head(label_cands, 8)

pB <- ggplot(scatter_dt, aes(x = max_logFC, y = human_shrunk_logFC)) +
  rasterize_layer(
    geom_point(
      aes(color = n_diets_up,
          shape = shape_coloc,
          alpha = fifelse(is_library, 0.85, 0.12)),
      size = 0.8, stroke = 0.15
    ), dpi = 300
  ) +
  scale_alpha_identity() +
  scale_color_viridis_c(name = "N diets UP", option = "viridis", direction = 1) +
  scale_shape_manual(name = "COLOC status",
                     values = c("COLOC PP4>0.5" = 16, "No COLOC" = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey50") +
  geom_text_repel(
    data = label_genes,
    aes(label = human_symbol),
    size = 1.8, max.overlaps = 20,
    segment.size = 0.2, segment.color = "grey40",
    fontface = "italic", color = "black",
    box.padding = 0.3, point.padding = 0.15, min.segment.length = 0.1
  ) +
  labs(x = "Mouse max shrunk logFC (across 4 diets)",
       y = "Human ashr shrunk logFC",
       title = "B") +
  theme_masld() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(face = "bold", size = 9),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        legend.title = element_text(size = 5.5),
        legend.text = element_text(size = 5)) +
  guides(color = guide_colorbar(barwidth = 0.4, barheight = 2.5),
         shape = guide_legend(override.aes = list(size = 1.5)))

# Correlation annotation
rho_val <- cor(scatter_dt$max_logFC, scatter_dt$human_shrunk_logFC,
               use = "complete.obs", method = "spearman")
n_concordant <- sum(scatter_dt$concordant, na.rm = TRUE)
n_total <- nrow(scatter_dt)

pB <- pB + annotate("text",
  x = max(scatter_dt$max_logFC, na.rm = TRUE) * 0.65,
  y = min(scatter_dt$human_shrunk_logFC, na.rm = TRUE) * 0.85,
  label = sprintf("rho = %.3f\n%d/%d concordant", rho_val, n_concordant, n_total),
  size = 1.8, hjust = 0, color = "grey30")

# =============================================================================
# Panel C -- Dot plot: per-diet enrichment for HumanDE / COLOC / Finemap
# =============================================================================
message("Generating Panel C (dotplot enrichment) ...")

# For each diet, test enrichment of that diet's UP DEGs among:
# (a) human DE genes, (b) COLOC genes, (c) finemapped genes
# Universe: all testable genes (present in at least one diet DE + bridged)

# Build universe: all genes with mouse DE output across any diet, bridged to human
all_mouse_genes <- unique(perdiet$gene_id_base)
universe <- merge(
  data.table(gene_id_base = all_mouse_genes),
  ortho[, .(gene_id_base = mouse_ensembl, human_ensembl)],
  by = "gene_id_base", allow.cartesian = TRUE
)
universe <- merge(universe, atlas[, .(human_ensembl, is_human_de, is_coloc, is_finemap)],
                  by = "human_ensembl", all.x = FALSE)
# Dedup mouse gene
universe <- universe[!duplicated(gene_id_base)]

enrich_results <- list()
for (d_i in DIETS) {
  # sig UP genes for this diet
  sig_genes <- perdiet[diet == d_i & sig_up == TRUE, unique(gene_id_base)]
  sig_in_univ <- intersect(sig_genes, universe$gene_id_base)
  nonsig_in_univ <- setdiff(universe$gene_id_base, sig_in_univ)

  for (evidence_type in c("HumanDE", "COLOC", "Finemap")) {
    ev_col <- switch(evidence_type,
                     HumanDE = "is_human_de",
                     COLOC = "is_coloc",
                     Finemap = "is_finemap")

    # 2x2 table
    a_val <- sum(universe[gene_id_base %in% sig_in_univ, get(ev_col)], na.rm = TRUE)
    b_val <- length(sig_in_univ) - a_val
    c_val <- sum(universe[gene_id_base %in% nonsig_in_univ, get(ev_col)], na.rm = TRUE)
    d_val <- length(nonsig_in_univ) - c_val

    ft <- fisher.test(matrix(c(a_val, b_val, c_val, d_val), nrow = 2))
    enrich_results[[length(enrich_results) + 1]] <- data.table(
      diet = d_i,
      evidence = evidence_type,
      overlap = a_val,
      diet_sig = length(sig_in_univ),
      evidence_total = a_val + c_val,
      universe_size = nrow(universe),
      odds_ratio = ft$estimate,
      p_value = ft$p.value
    )
  }
}
enrich_dt <- rbindlist(enrich_results)

# BH correction
enrich_dt[, padj := p.adjust(p_value, method = "BH")]
enrich_dt[, neg_log10p := -log10(pmax(p_value, 1e-300))]

# Cap odds ratio for display
enrich_dt[, or_display := pmin(odds_ratio, 10)]

# Significance marker
enrich_dt[, sig_label := fifelse(padj < 0.001, "***",
                                  fifelse(padj < 0.01, "**",
                                          fifelse(padj < 0.05, "*", "")))]

# Factor ordering
enrich_dt[, diet := factor(diet, levels = DIETS,
                           labels = DIET_LABELS[DIETS])]
enrich_dt[, evidence := factor(evidence, levels = c("HumanDE", "COLOC", "Finemap"))]

pC <- ggplot(enrich_dt, aes(x = diet, y = evidence)) +
  geom_point(aes(size = neg_log10p, fill = or_display),
             shape = 21, color = "grey30", stroke = 0.3) +
  geom_text(aes(label = sig_label), vjust = -0.6, size = 2, color = "black") +
  scale_size_continuous(name = expression(-log[10](p)),
                        range = c(1, 7), breaks = c(1, 3, 5, 10)) +
  scale_fill_gradient2(name = "Odds ratio",
                       low = "#1565C0", mid = "white", high = "#C9265E",
                       midpoint = 1, limits = c(0, 10),
                       oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "C") +
  theme_masld() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 5.5, angle = 35, hjust = 1),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 5.5),
        legend.text = element_text(size = 5)) +
  guides(size = guide_legend(override.aes = list(fill = "grey80")),
         fill = guide_colorbar(barwidth = 0.4, barheight = 2))

# =============================================================================
# Compose and save
# =============================================================================
message("Composing final figure ...")

layout <- "
AABB
AABB
CCCC
"

composite <- pA + pB + pC +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c("", "", "")))  # tags already in titles

out_path <- file.path(OUT_DIR, "S_lib_5_human_evidence.pdf")
save_fig(composite, out_path,
         width = 180 / 25.4, height = 100 / 25.4)

message("Saved: ", out_path)

# Print summary stats
message("\n--- Summary statistics ---")
message(sprintf("Total mouse UP genes (6-diet): %d", nrow(gene_stats)))
message(sprintf("Bridged to human: %d", sum(!is.na(bridged$human_ensembl))))
message(sprintf("Human evidence breakdown:"))
message(sprintf("  HumanDE+COLOC: %d", sum(bridged$evidence_state == "HumanDE+COLOC", na.rm = TRUE)))
message(sprintf("  HumanDE only: %d", sum(bridged$evidence_state == "HumanDE only", na.rm = TRUE)))
message(sprintf("  COLOC only: %d", sum(bridged$evidence_state == "COLOC only", na.rm = TRUE)))
message(sprintf("  Neither: %d", sum(bridged$evidence_state == "Neither", na.rm = TRUE)))
message(sprintf("Library tiers: Core=%d, Extended=%d, Exploratory=%d",
                sum(bridged$library_tier == "Core"),
                sum(bridged$library_tier == "Extended"),
                sum(bridged$library_tier == "Exploratory")))
message(sprintf("Scatter rho (mouse vs human logFC): %.3f", rho_val))
message(sprintf("Enrichment results:"))
for (i in seq_len(nrow(enrich_dt))) {
  r <- enrich_dt[i]
  message(sprintf("  %s x %s: OR=%.2f, p=%.2e, padj=%.2e %s",
                  r$diet, r$evidence, r$odds_ratio, r$p_value, r$padj, r$sig_label))
}
message("Done.")
