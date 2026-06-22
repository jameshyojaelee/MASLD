#!/usr/bin/env Rscript
# figS_robustness_atlas.R
# Composite summary of the per-gene robustness atlas across all 4 pillars.
# Panels:
#   E1: 4-pillar pass-summary stacked bar (canonical DEGs only)
#   E2: 4,724 DEGs heat-strip — gene × pillar membership (top 200 by sum_pillars)
#   E3: Pillar-pass joint distribution (UpSet-style, simplified)
#   E4: Headline-numbers panel (text)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(grid); library(gridExtra)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

AUDIT <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

atlas <- fread(file.path(AUDIT, "dream_robustness_atlas.csv"))
deg <- atlas[is_canonical_DEG == TRUE]

# E1 — n_pillars_passed stacked bar
e1_dt <- deg[, .N, by = n_pillars_passed]
e1_dt[, pct := round(100 * N / sum(N), 1)]
p1 <- ggplot(e1_dt, aes(x = factor(n_pillars_passed), y = N, fill = factor(n_pillars_passed))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", comma(N), pct)), vjust = -0.3) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "E1 — Pillar-pass count among canonical DEGs",
       subtitle = sprintf("Out of %s canonical DEGs (Pillar B is signature-level, not counted here)", comma(nrow(deg))),
       x = "# Pillars passed (A, C, D)", y = "DEGs", fill = "") +
  theme_pub() + theme(legend.position = "none")

# E2 — heat strip: genes × pillar pass
deg[, sum_pass := as.integer(A_pass) + as.integer(C_pass) + as.integer(D_pass)]
top200 <- deg[order(-sum_pass, -abs(dream_logFC))][1:min(200, nrow(deg))]  # C2-OK-sensitivity
top200[, gene_idx := seq_len(.N)]
heat_long <- melt(top200[, .(gene_idx, A = A_pass, C = C_pass, D = D_pass)],
                  id.vars = "gene_idx", variable.name = "pillar", value.name = "pass")
p2 <- ggplot(heat_long, aes(x = gene_idx, y = pillar, fill = pass)) +
  geom_tile() +
  scale_fill_manual(values = c("FALSE" = "grey92", "TRUE" = "#4B0082"), na.value = "grey95") +
  labs(title = "E2 — Pillar membership: top 200 DEGs by joint pass",
       x = "Gene rank (top 200)", y = "Pillar", fill = "Pass") +
  theme_pub() + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())

# E3 — UpSet-style: counts of pillar-pass intersections
pat <- deg[, .(A = A_pass, C = C_pass, D = D_pass)]
pat[, pat_str := paste0(ifelse(A, "A", "."), ifelse(C, "C", "."), ifelse(D, "D", "."))]
ix <- pat[, .N, by = pat_str][order(-N)]
p3 <- ggplot(ix, aes(x = reorder(pat_str, -N), y = N, fill = pat_str)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(N)), vjust = -0.3) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "E3 — Pillar pass-pattern distribution",
       x = "Pattern (A/C/D pass; '.' = no pass)", y = "DEGs", fill = "") +
  theme_pub() + theme(legend.position = "none",
                       axis.text.x = element_text(family = "mono"))

# E4 — Headline numbers
hd_lines <- c(
  sprintf("Canonical DEGs: %s", comma(nrow(deg))),
  "",
  sprintf("Pillar A (sample stability):"),
  sprintf("  Mean CPSS π̂ = %.3f", mean(deg$cpss_pi_hat, na.rm = TRUE)),
  sprintf("  π̂ ≥ 0.7: %.1f%%", 100 * mean(deg$cpss_pi_hat >= 0.7, na.rm = TRUE)),
  sprintf("  Mean bootstrap freq = %.3f", mean(deg$bootstrap_freq, na.rm = TRUE)),
  sprintf("  Bootstrap freq ≥ 0.9: %.1f%%", 100 * mean(deg$bootstrap_freq >= 0.9, na.rm = TRUE)),
  sprintf("  Mean kfold-10 recur: %.2f / 10", mean(deg$kfold10_recur, na.rm = TRUE)),
  sprintf("  Mean kfold-50 recur: %.2f / 50", mean(deg$kfold50_recur, na.rm = TRUE)),
  "",
  sprintf("Pillar C (empirical null):"),
  sprintf("  Canonical DEGs with |perm z| ≥ 3: %.1f%%",
          100 * mean(abs(deg$perm_z) >= 3, na.rm = TRUE)),
  "",
  sprintf("Pillar D (batch RE):"),
  sprintf("  Canonical DEGs SVA-concordant: %.1f%%",
          100 * mean(deg$sva_concordant, na.rm = TRUE)),
  sprintf("  Median residual cohort var: %.4f",
          median(deg$residual_cohort_var, na.rm = TRUE)),
  "",
  sprintf("Joint pass (all 3 of A/C/D): %s (%.1f%%)",
          comma(sum(deg$n_pillars_passed == 3, na.rm = TRUE)),
          100 * mean(deg$n_pillars_passed == 3, na.rm = TRUE))
)
p4 <- ggplot() +
  annotate("text", x = 0, y = seq_along(hd_lines),
           label = rev(hd_lines), hjust = 0, family = "mono", size = 3.4) +
  xlim(-0.05, 1) + ylim(0, length(hd_lines) + 1) +
  labs(title = "E4 — Headline numbers", x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

combined <- (p1 | p3) / (p2) / (p4) + plot_layout(heights = c(1, 0.8, 1.2))
out_pdf <- file.path(OUT_DIR, "figS_robustness_atlas.pdf")
ggsave(out_pdf, combined, width = 12, height = 12, device = cairo_pdf)
cat("Saved:", out_pdf, "\n")
