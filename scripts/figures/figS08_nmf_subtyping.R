#!/usr/bin/env Rscript
# =============================================================================
# figS08_nmf_subtyping.R — Supplementary figure: k-program molecular subtyping
#
# Reads chosen k from program_labels.csv and renders five panels:
#   (a) k-selection quality metrics (cophenetic correlation and silhouette score) across k=3–10
#   (b) Program × Hallmark-pathway NES heatmap (per-program biological profile)
#   (c) Program composition across fibrosis stages (stacked bar, dominant program per donor)
#   (d) Per-program F1–F3 usage inflection (ratio of F3 to F1 mean usage)
#   (e) Bulk NMF program × cell-type validation: program z-score Δ-vs-Healthy
#       trajectory per cell-type umbrella, scored on the integrated scRNA atlas
#       (1.23M cells; bulk top-10 program genes z-scored). Sanity check that
#       bulk NMF programs are dominated by their constituent cell types and
#       that those cell-type signals shift with disease stage.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(NMF)
  library(RColorBrewer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

subtype_dir  <- file.path(BASE, "RNA-seq/results/subtypes")
cache_path   <- Sys.getenv("NMF_CACHE_PATH",
  unset = file.path(subtype_dir, "nmf_results_cache_clean.rds"))
assign_path  <- file.path(subtype_dir, "nmf_assignments.csv")
labels_path  <- file.path(subtype_dir, "program_labels.csv")
pathways_path <- file.path(subtype_dir, "subtype_pathways.csv")
markers_path <- file.path(subtype_dir, "subtype_markers.csv")
meta_path    <- file.path(INTEGRATION, "metadata/unified_metadata.csv")

out_dir <- FIGS08_DIR
dir.create(file.path(out_dir, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== figS08_nmf_subtyping (k-program refactor) ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# Load
nmf_cache <- readRDS(cache_path)
metrics   <- as.data.table(nmf_cache$metrics)
assignments <- fread(assign_path)
labels <- fread(labels_path)
meta <- fread(meta_path)
chosen_k <- nrow(labels)
program_codes <- labels$program_code
program_bios  <- labels$biological_label
cat(sprintf("  chosen_k = %d\n", chosen_k))
cat(sprintf("  Programs: %s\n", paste(program_bios, collapse=", ")))

# Color palette for k programs
prog_palette <- setNames(
  colorRampPalette(brewer.pal(min(max(chosen_k,3),8), "Set2"))(chosen_k),
  program_bios
)

# ============================================================
# Panel (a): rank survey with chosen_k highlight
# ============================================================
cat("\n--- Panel (a): rank survey ---\n")
ml <- melt(metrics[, .(k, cophenetic, silhouette)], id.vars = "k",
           variable.name = "metric", value.name = "value")
ml[, metric_label := ifelse(metric == "cophenetic",
                            "Cophenetic coefficient", "Mean silhouette")]

pa <- ggplot(ml, aes(k, value, color = metric_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_k, linetype = "dashed",
             color = "#C2185B", linewidth = 0.4) +
  annotate("text", x = chosen_k + 0.15, y = 0.97,
           label = sprintf("chosen\nk = %d", chosen_k),
           hjust = 0, size = 2.2, color = "#C2185B", fontface = "bold") +
  scale_x_continuous(breaks = metrics$k) +
  scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
  scale_color_manual(values = c("Cophenetic coefficient" = masld_colors$up,
                                "Mean silhouette" = masld_colors$down)) +
  labs(x = "Number of programs (k)", y = "Score", color = NULL,
       title = "NMF rank survey (k sweep)") +
  theme_masld() +
  theme(legend.position = c(0.7, 0.3), legend.background = element_blank())
save_fig(pa, file.path(out_dir, "panels", "panel_a.pdf"),
         width = fig_half_width, height = 2.8)

# ============================================================
# Panel (b): per-program Hallmark NES heatmap (top 15 pathways)
# ============================================================
cat("\n--- Panel (b): program × Hallmark NES ---\n")
if (file.exists(pathways_path)) {
  paths <- fread(pathways_path)
  # paths has: pathway, padj, NES, program (P1..Pk), program_bio
  # Pick top 15 pathways by |NES| across all programs (union)
  top_paths <- paths[padj < 0.05][order(-abs(NES))][1:min(.N, 30)]
  # Ensure we have at least one hit per program
  top_by_prog <- paths[padj < 0.05, .SD[order(-abs(NES))][1:5], by = program]
  top_pw <- unique(c(top_paths$pathway, top_by_prog$pathway))
  hm <- paths[pathway %in% top_pw, .(program, program_bio, pathway, NES, padj)]
  hm[, pathway_short := gsub("^HALLMARK_", "", pathway)]
  hm[, pathway_short := gsub("_", " ", pathway_short)]
  hm[, pathway_short := tools::toTitleCase(tolower(pathway_short))]
  # Order programs per labels
  hm[, program_bio := factor(program_bio, levels = program_bios)]

  pb <- ggplot(hm, aes(program_bio, pathway_short, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = ifelse(padj < 0.05, sprintf("%.1f", NES), "")),
              size = 1.8) +
    scale_fill_gradient2(low = masld_colors$down, mid = "grey95",
                         high = masld_colors$up, midpoint = 0, name = "NES") +
    labs(x = NULL, y = NULL,
         title = sprintf("NMF programs × Hallmark enrichment (k=%d)", chosen_k)) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          axis.text.y = element_text(size = 6.5))
} else {
  pb <- ggplot() + labs(title = "(b) pathways file missing") + theme_masld()
}
save_fig(pb, file.path(out_dir, "panels", "panel_b.pdf"),
         width = fig_half_width, height = 5)

# ============================================================
# Panel (c): dominant program composition by fibrosis stage
# ============================================================
cat("\n--- Panel (c): dominant program × fibrosis stage ---\n")
comp <- merge(assignments[, .(sample_id, dominant_program)],
              meta[, .(sample_id, fibrosis_stage)], by = "sample_id")
comp <- comp[!is.na(fibrosis_stage)]
comp[, fibrosis_label := paste0("F", fibrosis_stage)]
prop_dt <- comp[, .N, by = .(fibrosis_label, dominant_program)]
prop_dt[, total := sum(N), by = fibrosis_label]
prop_dt[, pct := N / total * 100]
prop_dt[, fibrosis_label := factor(fibrosis_label, levels = paste0("F", 0:4))]
prop_dt[, dominant_program := factor(dominant_program, levels = program_bios)]

pc <- ggplot(prop_dt, aes(fibrosis_label, pct, fill = dominant_program)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = prog_palette, name = "Program") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Fibrosis stage", y = "Proportion",
       title = sprintf("NMF program composition by fibrosis stage (k=%d)", chosen_k)) +
  theme_masld() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))
save_fig(pc, file.path(out_dir, "panels", "panel_c.pdf"),
         width = fig_half_width, height = 3.2)

# ============================================================
# Panel (d): mid-stage (F1-F3) inflection ratio per program
# ============================================================
cat("\n--- Panel (d): mid-stage (F1-F3) inflection per program ---\n")
comp[, f2_bin := ifelse(fibrosis_stage <= 2, "F0-F2", "F3-F4")]
sw_dt <- comp[, .(.N), by = .(dominant_program, f2_bin)][
  CJ(dominant_program = program_bios, f2_bin = c("F0-F2","F3-F4"), unique = TRUE),
  on = c("dominant_program","f2_bin")
]
sw_dt[is.na(N), N := 0]
sw_dt[, pct := N / sum(N) * 100, by = f2_bin]
sw_wide <- dcast(sw_dt, dominant_program ~ f2_bin, value.var = "pct", fill = 0)
sw_wide[, switch_ratio := `F3-F4` / pmax(`F0-F2`, 0.001)]
sw_wide[, dominant_program := factor(dominant_program, levels = program_bios)]
setorder(sw_wide, -switch_ratio)

pd <- ggplot(sw_wide, aes(switch_ratio, reorder(dominant_program, switch_ratio),
                          fill = dominant_program)) +
  geom_col(width = 0.7, color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("%.2fx", switch_ratio)),
            hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = prog_palette, guide = "none") +
  expand_limits(x = max(sw_wide$switch_ratio, na.rm = TRUE) * 1.2) +
  labs(x = "F3-F4 vs F0-F2 dominance ratio", y = NULL,
       title = "NMF program mid-stage (F1-F3) inflection") +
  theme_masld()
save_fig(pd, file.path(out_dir, "panels", "panel_d.pdf"),
         width = fig_half_width, height = 3)

# ============================================================
# Panel e: bulk NMF program x cell-type validation trajectory
#   Aggregated by Analysis/SingleCell/scripts/aggregate_nmf_celltype_stage.py
#   from `bulk_nmf_celltype_stage_mean.csv` (top-10 program genes z-scored on
#   the integrated atlas, per (cell_type x disease_stage_coarse)).  Pooled
#   into 6 lineage umbrellas (Hepatocytes / Cholangiocytes / Endothelial /
#   Stromal-Fibroblasts / Phagocytes / Lymphocytes) by n_cells-weighted mean,
#   then expressed as Delta-vs-Healthy.
# ============================================================
bulk_nmf_cts_path <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/disease_signatures/bulk_nmf_celltype_stage_mean.csv")

pe <- tryCatch({
  if (!file.exists(bulk_nmf_cts_path))
    stop("aggregator output missing: ", bulk_nmf_cts_path)

  bulk_nmf_palette <- c("P1" = "#C2185B", "P2" = "#7B1FA2", "P3" = "#1565C0",
                        "P4" = "#00897B", "P5" = "#558B2F", "P6" = "#E65100")
  bulk_nmf_legend <- c(
    "P1" = "P1 Pro-inflammatory", "P2" = "P2 Innate-immune",
    "P3" = "P3 Parenchymal",      "P4" = "P4 lncRNA",
    "P5" = "P5 Hepatic-metabolic","P6" = "P6 Stellate-myofibroblast"
  )
  umbrella_map_e <- list(
    "Hepatocytes"         = "Hepatocytes",
    "Cholangiocytes"      = "Cholangiocytes",
    "Endothelial"         = "Endothelial cells",
    "Stromal/Fibroblasts" = "Fibroblasts",
    "Phagocytes"          = c("Macrophages", "Mono+mono derived cells",
                              "cDC1s", "cDC2s", "pDCs",
                              "Neutrophils", "Basophils"),
    "Lymphocytes"         = c("T cells", "Resident NK", "Circulating NK/NKT",
                              "B cells", "Plasma cells")
  )
  ct_to_umb <- setNames(rep(names(umbrella_map_e), lengths(umbrella_map_e)),
                        unlist(umbrella_map_e))
  stage_lvls <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
  prog_lvls <- paste0("P", 1:6)

  d <- fread(bulk_nmf_cts_path)
  d[, umbrella := ct_to_umb[cell_type]]
  d <- d[!is.na(umbrella)]

  long_raw <- melt(d, id.vars = c("umbrella", "stage", "n_cells"),
                   measure.vars = prog_lvls,
                   variable.name = "program", value.name = "score")
  agg <- long_raw[, .(score   = sum(score * n_cells) / sum(n_cells),
                       n_cells = sum(n_cells)),
                  by = .(umbrella, stage, program)]
  agg <- agg[n_cells >= 400]
  agg[, score_healthy := score[stage == "Healthy"][1L],
      by = .(umbrella, program)]
  agg[, delta := score - score_healthy]
  agg <- agg[!is.na(delta)]
  agg[, umbrella := factor(umbrella, levels = names(umbrella_map_e))]
  agg[, stage    := factor(stage, levels = stage_lvls)]
  agg[, program  := factor(program, levels = prog_lvls)]
  setorder(agg, umbrella, program, stage)

  ggplot(agg, aes(x = stage, y = delta, color = program, group = program)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray70", linewidth = 0.25) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    geom_point(size = 0.9, alpha = 0.9) +
    facet_wrap(~ umbrella, ncol = 3) +
    scale_color_manual(values = bulk_nmf_palette,
                       labels = bulk_nmf_legend, name = NULL,
                       guide = guide_legend(ncol = 1, keyheight = unit(8, "pt"))) +
    labs(x = NULL,
         y = expression(Delta * " mean program z-score (vs Healthy)"),
         title = "Bulk NMF program shift vs Healthy, by cell-type umbrella") +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.background = element_rect(fill = "gray95", color = NA),
          panel.spacing = unit(4, "pt"),
          legend.position = "right")
}, error = function(e) {
  cat("  Panel e skipped: ", e$message, "\n")
  ggplot() + theme_void()
})

save_fig(pe, file.path(out_dir, "panels", "panel_e.pdf"),
         width = fig_full_width, height = 4)

# ============================================================
# Composite
# ============================================================
tryCatch({
  composite <- ((pa | pb) / (pc | pd) / pe) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  save_fig(composite, file.path(out_dir, "figS08_nmf_subtyping.pdf"),
           width = fig_full_width, height = 12)
  cat("Saved figS08_nmf_subtyping.pdf\n")
}, error = function(e) cat("Composite assembly failed:", e$message, "\n"))
cat("Done.\n")
