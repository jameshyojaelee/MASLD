#!/usr/bin/env Rscript
# figS02_progression_characterization.R
# Outputs each panel as an individual PDF to figures/supplementary/figS02_progression/.
# No composite layouts — assemble in Illustrator.
#
# Outputs:
#   cohort_diagnosis_composition.pdf     — stacked bar by cohort and disease category
#   nas_fibrosis_sample_heatmap.pdf      — NAS × Fibrosis sample count heatmap
#   deg_counts_nas_vs_baseline.pdf       — DEG counts each NAS vs NAS 0
#   deg_counts_fibrosis_vs_baseline.pdf  — DEG counts each fibrosis stage vs F0
#   deg_counts_fibrosis_transitions.pdf  — incremental DEGs F0→F1, F1→F2, F2→F3, F3→F4
#   deg_counts_nas_transitions.pdf       — incremental DEGs NAS0→1, 1→2, 2→3, …
#
# All DEG panels: padj < 0.05, NO LFC filter (Tier-2 progression convention; LFC filter dropped 2026-06-09 harmonization)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PADJ <- 0.05
# Tier-2 progression convention: padj < 0.05, NO LFC floor (binary stage
# contrasts dilute per-gene fold changes). No |logFC| threshold is applied.
OUT  <- FIGS02_DIR

save_panel <- function(p, name, w, h) {
  path <- file.path(OUT, paste0(name, ".pdf"))
  ggsave(path, p, width = w, height = h, device = "pdf")
  cat("Saved:", basename(path), "\n")
}

# ── helper: bidirectional bar chart ──────────────────────────────────────────
make_bar <- function(dt, x_col, xlab, flip_angle = 0) {
  dt[, direction := fifelse(logFC > 0, "Up", "Down")]
  up   <- dt[direction == "Up",   .(n =  .N, direction = "Up"),   by = c(x_col)]
  down <- dt[direction == "Down", .(n = -.N, direction = "Down"), by = c(x_col)]
  bar  <- rbindlist(list(up, down))

  all_lvls <- unique(dt[[x_col]])
  complete_grid <- CJ(lvl = all_lvls, direction = c("Up", "Down"))
  setnames(complete_grid, "lvl", x_col)
  bar <- merge(complete_grid, bar, by = c(x_col, "direction"), all.x = TRUE)
  bar[is.na(n), n := 0]
  bar[, (x_col) := factor(get(x_col))]
  bar[, lbl := fifelse(n == 0L, "", as.character(abs(n)))]

  ggplot(bar, aes(x = .data[[x_col]], y = n, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = lbl, vjust = ifelse(n >= 0, -0.3, 1.3)),
              size = GEOM_TEXT_6PT, color = "black") +
    scale_y_continuous(labels = function(x) comma(abs(x)),
                       expand = expansion(mult = c(0.12, 0.15))) +
    scale_fill_manual(values = c(Up = masld_colors$up, Down = masld_colors$down), name = NULL) +
    labs(x = xlab,
         y = paste0("DEGs  (padj < ", PADJ, ")")) +
    theme_masld() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"),
          legend.text     = element_text(size = 6),
          axis.text.x     = element_text(size = 6, angle = flip_angle,
                                         hjust = ifelse(flip_angle > 0, 1, 0.5)),
          axis.title      = element_text(size = 6))
}

# ─────────────────────────────────────────────────────────────────────────────
# cohort_diagnosis_composition
# ─────────────────────────────────────────────────────────────────────────────
cat("=== cohort_diagnosis_composition ===\n")
meta     <- load_metadata()
qc       <- fread(file.path(INTEGRATION, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
meta_qc  <- meta[sample_id %in% pass_ids]

meta_qc[, diag_cat := fcase(
  group_binary == "Control",                          "Healthy",
  diagnosis_harmonized == "NAFL",                     "MASL",
  diagnosis_harmonized %in% c("NASH", "Borderline"),  "MASH",
  grepl("NASH",   condition),                         "MASH",
  grepl("F[34]",  condition),                         "MASH",
  grepl("F3F4",   condition),                         "MASH",
  grepl("F[012]", condition),                         "MASL",
  grepl("F0F1",   condition),                         "MASL",
  default = "MASH"
)]
meta_qc[, diag_cat := factor(diag_cat, levels = c("Healthy", "MASL", "MASH"))]

ds_order   <- meta_qc[, .N, by = dataset][order(-N), dataset]
meta_qc[, dataset := factor(dataset, levels = ds_order)]
n_ds       <- length(ds_order)
ds_palette <- setNames(
  colorRampPalette(c("#0D47A1", "#1565C0", "#42A5F5", "#7B1FA2",
                     "#AD1457", "#C2185B", "#E91E63", "#F48FB1",
                     "#00695C", "#F57F17"))(n_ds),
  ds_order
)

p_cohort <- ggplot(meta_qc, aes(x = diag_cat, fill = dataset)) +
  geom_bar(width = 0.7) +
  scale_fill_manual(values = ds_palette, name = "Cohort") +
  labs(x = NULL, y = "Samples") +
  theme_masld() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title    = element_text(size = 6),
        legend.text     = element_text(size = 6),
        axis.text.x     = element_text(size = 6))

save_panel(p_cohort, "cohort_diagnosis_composition", w = 5, h = 4)

# ─────────────────────────────────────────────────────────────────────────────
# nas_fibrosis_sample_heatmap
# ─────────────────────────────────────────────────────────────────────────────
cat("=== nas_fibrosis_sample_heatmap ===\n")
sd <- load_sample_distribution()

if (!is.null(sd) && nrow(sd) > 0) {
  sd[, nas_group := factor(nas_group)]
  sd[, fib_stage := factor(fib_stage)]

  p_heatmap <- ggplot(sd, aes(x = nas_group, y = fib_stage, fill = N)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = N), size = GEOM_TEXT_6PT, color = "black") +
    scale_fill_gradient(low = "#E3F2FD", high = "#0D47A1",
                        name = "n", na.value = "gray95") +
    labs(x = "NAS Score", y = "Fibrosis Stage") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          panel.grid      = element_blank(),
          axis.line       = element_blank())

  save_panel(p_heatmap, "nas_fibrosis_sample_heatmap", w = 5, h = 4)
} else {
  cat("WARNING: sample distribution data not found — skipping heatmap\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# deg_counts_nas_vs_baseline
# ─────────────────────────────────────────────────────────────────────────────
cat("=== deg_counts_nas_vs_baseline ===\n")
nas_dream <- load_nas_score_progression()
nas_sizes <- load_nas_score_sample_sizes()

if (!is.null(nas_dream) && nrow(nas_dream) > 0) {
  nas_sig <- nas_dream[padj < PADJ]
  all_lvls <- sort(unique(nas_dream$nas_level))
  nas_sig[, stage_label := factor(paste0("NAS", nas_level),
                                  levels = paste0("NAS", all_lvls))]

  p_nas_base <- make_bar(nas_sig[, .(logFC, stage_label)],
                         x_col = "stage_label",
                         xlab  = "NAS score (vs NAS 0)")

  if (!is.null(nas_sizes)) {
    ns_lbl <- nas_sizes[nas_group %in% all_lvls,
                        .(stage_label = factor(paste0("NAS", nas_group),
                                               levels = paste0("NAS", all_lvls)),
                          N)]
    max_dn <- nas_sig[logFC < 0, .N, by = nas_level][, if (.N) max(N) else 0]
    p_nas_base <- p_nas_base +
      geom_text(data = ns_lbl,
                aes(x = stage_label, y = -max_dn * 1.18, label = paste0("n=", N)),
                inherit.aes = FALSE, size = GEOM_TEXT_6PT, color = "black")
  }

  save_panel(p_nas_base, "deg_counts_nas_vs_baseline", w = 5, h = 4)
} else {
  cat("WARNING: NAS score dream data not found\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# deg_counts_fibrosis_vs_baseline
# ─────────────────────────────────────────────────────────────────────────────
cat("=== deg_counts_fibrosis_vs_baseline ===\n")
fib_dream <- load_fibrosis_stage_dream()
fib_sizes <- load_fibrosis_stage_sample_sizes()

if (!is.null(fib_dream) && nrow(fib_dream) > 0) {
  fib_sig  <- fib_dream[padj < PADJ]
  all_fib  <- sort(unique(fib_dream$fib_stage))
  fib_sig[, stage_label := factor(paste0("F", fib_stage),
                                  levels = paste0("F", all_fib))]

  p_fib_base <- make_bar(fib_sig[, .(logFC, stage_label)],
                         x_col = "stage_label",
                         xlab  = "Fibrosis stage (vs F0)")

  if (!is.null(fib_sizes)) {
    fs_lbl <- fib_sizes[fib_stage %in% all_fib,
                        .(stage_label = factor(paste0("F", fib_stage),
                                               levels = paste0("F", all_fib)),
                          N)]
    max_dn <- fib_sig[logFC < 0, .N, by = fib_stage][, if (.N) max(N) else 0]
    p_fib_base <- p_fib_base +
      geom_text(data = fs_lbl,
                aes(x = stage_label, y = -max_dn * 1.18, label = paste0("n=", N)),
                inherit.aes = FALSE, size = GEOM_TEXT_6PT, color = "black")
  }

  save_panel(p_fib_base, "deg_counts_fibrosis_vs_baseline", w = 4, h = 4)
} else {
  cat("WARNING: fibrosis stage dream data not found\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# deg_counts_fibrosis_transitions  (F0→F1, F1→F2, F2→F3, F3→F4)
# ─────────────────────────────────────────────────────────────────────────────
cat("=== deg_counts_fibrosis_transitions ===\n")
fib_consec <- load_fibrosis_consecutive()

if (!is.null(fib_consec) && nrow(fib_consec) > 0) {
  fib_inc <- fib_consec[padj < PADJ]
  fib_inc[, trans_label := gsub("F(\\d+)_vs_F(\\d+)", "F\\2→F\\1", contrast)]
  fib_inc[, order_idx   := as.integer(gsub(".*F(\\d+)$", "\\1", trans_label))]
  trans_lvls <- fib_inc[, .(order_idx = first(order_idx)), by = trans_label
                        ][order(order_idx), trans_label]
  fib_inc[, trans_label := factor(trans_label, levels = trans_lvls)]

  p_fib_trans <- make_bar(fib_inc[, .(logFC, stage_label = trans_label)],
                          x_col = "stage_label",
                          xlab  = NULL)

  # sample size annotation
  ts     <- fread(file.path(INTEGRATION, "results/progression/transition_summary.csv"))
  ts_fib <- ts[stage_type == "fibrosis",
               .(trans_label = factor(gsub("_to_", "→", transition),
                                      levels = trans_lvls),
                 n_samples)]
  max_dn <- fib_inc[logFC < 0, .N, by = trans_label][, if (.N) max(N) else 0]
  p_fib_trans <- p_fib_trans +
    geom_text(data = ts_fib[!is.na(trans_label)],
              aes(x = trans_label, y = -max_dn * 1.18, label = paste0("n=", n_samples)),
              inherit.aes = FALSE, size = GEOM_TEXT_6PT, color = "black")

  save_panel(p_fib_trans, "deg_counts_fibrosis_transitions", w = 4.5, h = 4)
} else {
  cat("WARNING: fibrosis consecutive dream data not found\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# deg_counts_nas_transitions  (NAS0→1, 1→2, 2→3, …)
# ─────────────────────────────────────────────────────────────────────────────
cat("=== deg_counts_nas_transitions ===\n")
nas_consec <- load_nas_consecutive()

if (!is.null(nas_consec) && nrow(nas_consec) > 0) {
  nas_inc <- nas_consec[padj < PADJ]
  nas_inc[, trans_label := gsub("NAS(\\d+)_vs_NAS(\\d+)", "NAS\\2→\\1", contrast)]
  nas_inc[, order_idx   := as.integer(gsub(".*→(\\d+)$", "\\1", trans_label))]
  nas_lvls <- nas_inc[, .(order_idx = first(order_idx)), by = trans_label
                      ][order(order_idx), trans_label]
  nas_inc[, trans_label := factor(trans_label, levels = nas_lvls)]

  p_nas_trans <- make_bar(nas_inc[, .(logFC, stage_label = trans_label)],
                          x_col = "stage_label",
                          xlab  = NULL)

  save_panel(p_nas_trans, "deg_counts_nas_transitions", w = 6, h = 4)
} else {
  cat("WARNING: NAS consecutive dream data not found\n")
}

cat("\nAll panels saved to:", OUT, "\n")
