#!/usr/bin/env Rscript
# ============================================================================
# 343u_fstage_final_figures.R
#
# Manuscript-ready figures for the F-stage inference investigation:
#   Panel A: SH F-stage distribution per method (median + bootstrap 95% CI)
#   Panel B: Bootstrap gate-pass heatmap (4 methods × 5 gates)
#   Panel C: Head-to-head method comparison (LOOCV QWK / external rho / Hlt->F4)
#   Panel D: Anchor jackknife influence (ranked per-anchor)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
SUPP_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

up <- masld_colors$up        # "#C2185B" disease/up magenta
dn <- masld_colors$down      # "#1565C0" control/down blue
ns <- masld_colors$ns        # "#BDBDBD" gray
mash <- masld_colors$mash    # "#C2185B"
masl <- masld_colors$masl    # "#F48FB1"

# F-stage gradient: F0 gray (control rule); F1->F4 pink->dark magenta
fstage_palette <- c(
  "F0" = ns,
  "F1" = "#F8BBD0",
  "F2" = masl,
  "F3" = mash,
  "F4" = "#880E4F"
)
method_palette <- c(
  "scvi_ord_augmented" = mash,     # SOTA
  "scvi_ord_vanilla"   = masl,
  "scanvi_vanilla"     = dn,
  "scanvi_augmented"   = ns
)
method_label <- c(
  "scvi_ord_augmented" = "scVI ord. augmented (SOTA)",
  "scvi_ord_vanilla"   = "scVI ord. vanilla",
  "scanvi_vanilla"     = "scANVI vanilla",
  "scanvi_augmented"   = "scANVI augmented"
)

# ---------------------------------------------------------------------------
# Panel A: SH F-stage distribution per method (bootstrap median + 95% CI)
# ---------------------------------------------------------------------------
boot <- fread(file.path(OUT_DIR, "bootstrap_sh_distribution_stability.tsv"))
boot[, fstage := paste0("F", fstage)]
boot[, method := factor(method, levels = names(method_label),
                                labels = method_label)]
boot[, fstage := factor(fstage, levels = paste0("F", 0:4))]

boot_summary <- boot[, .(median = median(fraction),
                          q025   = quantile(fraction, 0.025),
                          q975   = quantile(fraction, 0.975)),
                     by = .(method, fstage)]

p_dist <- ggplot(boot_summary, aes(x = fstage, y = median, fill = fstage)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.25,
                color = "grey25", linewidth = 0.3) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(values = fstage_palette, guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0),
                     limits = c(0, 0.85)) +
  labs(x = NULL, y = "fraction of SH donor cells",
       title = "Steatohepatitis F-stage distribution",
       subtitle = "bootstrap median +/- 95% CI (B=500); SH donors n=101") +
  theme_masld(base_size = 7) +
  theme(strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6),
        plot.title.position = "plot",
        plot.margin = margin(3, 3, 3, 3))

# ---------------------------------------------------------------------------
# Panel B: Bootstrap gate-pass heatmap
# ---------------------------------------------------------------------------
gates <- fread(file.path(OUT_DIR, "bootstrap_gates_redesign.tsv"))
gate_cols <- c("gateA_F2plus_50pct", "gateB_F3plus_majority",
               "gateC_F3plus_mode", "gateD_no_bimodal",
               "gateE_F3_unique_mode")
gates_long <- melt(gates, id.vars = "method", measure.vars = gate_cols,
                   variable.name = "gate", value.name = "pass_pct")
gates_long[, gate := factor(gate, levels = gate_cols, labels = c(
  "A: F2+ ≥ 50%",
  "B: F3+F4 majority",
  "C: F3+F4 > each early",
  "D: no F1/F4 bimodal",
  "E: F3 unique mode (strict)"
))]
gates_long[, method := factor(method, levels = names(method_label),
                                       labels = method_label)]
gates_long[, pass_pct100 := pass_pct * 100]
gates_long[, label := sprintf("%.0f%%", pass_pct100)]

p_gate <- ggplot(gates_long, aes(x = gate, y = method, fill = pass_pct100)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label, color = pass_pct100 > 50), size = 2.4,
            fontface = "bold") +
  scale_fill_gradient2(low = dn, mid = "white", high = up,
                       midpoint = 50, limits = c(0, 100),
                       name = "% of B=500\nbootstraps passing") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey25"),
                     guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Bootstrap stability gates (biology-meaningful)",
       subtitle = "pass = >= 95% of bootstraps satisfy the gate") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.key.width = unit(0.25, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        plot.title.position = "plot",
        plot.margin = margin(3, 3, 3, 3))

# ---------------------------------------------------------------------------
# Panel C: Head-to-head metric comparison
# ---------------------------------------------------------------------------
metrics <- data.table(
  method   = factor(names(method_label), levels = names(method_label),
                                          labels = method_label),
  loocv_qwk = c(0.720, 0.742, 0.782, 0.635),
  external_rho = c(0.639, 0.391, 0.412, 0.543),
  healthy_F4 = c(2, 27, 13, 10)
)
# methods order: scvi_ord_augmented, scvi_ord_vanilla, scanvi_vanilla, scanvi_augmented
metrics_long <- melt(metrics, id.vars = "method",
                     variable.name = "metric", value.name = "value")
metrics_long[, metric := factor(metric, levels = c(
  "loocv_qwk", "external_rho", "healthy_F4"),
  labels = c("Andrews LOOCV QWK", "External Spearman rho",
             "Healthy -> F4 errors"))]

p_metric <- ggplot(metrics_long,
                   aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(metric == "Healthy -> F4 errors",
                               sprintf("%.0f", value),
                               sprintf("%.3f", value))),
            size = 2.3, vjust = -0.4, fontface = "bold",
            color = "grey20") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = setNames(method_palette[names(method_label)],
                                       method_label), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = NULL,
       title = "Head-to-head: 4 F-stage methods",
       subtitle = "Andrews LOOCV intra-cohort; external rho cross-cohort (n_ext=260)") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
        strip.text = element_text(face = "bold", size = 7),
        plot.title.position = "plot",
        plot.margin = margin(3, 3, 3, 3))

# ---------------------------------------------------------------------------
# Panel D: Anchor jackknife (top 15 most-influential)
# ---------------------------------------------------------------------------
jk_path <- file.path(OUT_DIR, "anchor_jackknife_influence.tsv")
if (file.exists(jk_path)) {
  jk <- fread(jk_path)
  # Schema (confirmed): rank, sample, dataset, origin, anchor_F,
  # disease_stage_coarse, influence_argmax_mean_abs, influence_proba_l1_mean,
  # n_sh_flipped, n_sh_flipped_higher, n_sh_flipped_lower, max_abs_delta_F
  jk[, influence    := influence_argmax_mean_abs]
  jk[, anchor_class := origin]
  jk[, anchor_id    := sample]
  jk <- jk[order(-influence)]
  jk_top <- jk[1:min(15, nrow(jk))]
  jk_top[, anchor_id := factor(anchor_id, levels = rev(anchor_id))]

  anchor_class_pal <- c(
    "documented"            = mash,
    "clean_healthy_anchor"  = ns,
    "cirrhosis_anchor"      = "#880E4F"
  )

  p_jk <- ggplot(jk_top, aes(x = influence, y = anchor_id,
                              fill = anchor_class)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = anchor_class_pal, name = "anchor class") +
    labs(x = "leave-one-out influence (mean dF-stage on 101 SH donors)",
         y = NULL,
         title = "Top-15 most influential anchor donors",
         subtitle = "Top-5 share 16.8% of total influence (< 30% threshold)") +
    theme_masld(base_size = 7) +
    theme(axis.text.y = element_text(size = 5.5),
          legend.position = "top",
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          plot.title.position = "plot",
          plot.margin = margin(3, 3, 3, 3))
} else {
  p_jk <- ggplot() + theme_void() +
    labs(title = "(anchor_jackknife_influence.tsv missing)")
}

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------
top_row    <- p_dist
middle_row <- p_metric
bot_row    <- p_gate | p_jk
bot_row    <- bot_row + plot_layout(widths = c(1, 1))

panel_full <- (top_row / middle_row / bot_row) +
  plot_layout(heights = c(1, 1, 1.2)) +
  plot_annotation(
    title = "F-stage inference: final method comparison + stability verdict",
    subtitle = "SOTA = scVI ord. augmented. All biology-meaningful gates pass; F3 strict-mode gate fails (F3+F4 close).",
    theme = theme(plot.title = element_text(face = "bold", size = 10),
                  plot.subtitle = element_text(size = 8, color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 10))

out_pdf <- file.path(SUPP_DIR, "figS_fstage_final_methods.pdf")
ggsave(out_pdf, panel_full, width = 11, height = 11, device = cairo_pdf)
cat(sprintf("[output] %s\n", out_pdf))

# ---------------------------------------------------------------------------
# Smaller second figure: SH distribution alone, manuscript headline
# ---------------------------------------------------------------------------
sh_one <- p_dist +
  labs(title = "Steatohepatitis F-stage distribution: scVI ordinal augmented",
       subtitle = "bootstrap median +/- 95% CI (B=500; per-method panels for comparison)") +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, color = "grey30"))

out_pdf2 <- file.path(SUPP_DIR, "figS_sh_fstage_distribution_final.pdf")
ggsave(out_pdf2, sh_one, width = 9, height = 4, device = cairo_pdf)
cat(sprintf("[output] %s\n", out_pdf2))

cat("\n[done] final F-stage figures written\n")
