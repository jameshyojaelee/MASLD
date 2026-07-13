#!/usr/bin/env Rscript
# ============================================================================
# fibrosis_stage_directional_asymmetry.R
# Supp Fig S02 — up:down DEG ratio collapses as fibrosis advances.
#
# Per fibrosis-stage contrast (F1..F4 vs F0) count significant up- and
# down-regulated genes (C2 canonical: analytical TREAT, fdr_treat < 0.05 at lfc=0.25)
# and plot the up:down ratio as a bar chart. Early fibrosis is induction-dominated
# (F1 = induction only, 0 down; F2 ~8.8:1); advanced fibrosis approaches parity as
# downregulation (parenchymal collapse) overtakes induction (F4 ~2.1:1).
#
# Source: fibrosis_stage_dream.csv (limma_voom_qw per-stage-vs-F0 contrasts)
#
# Output: figures/supplementary/figS02_progression/fibrosis_stage_directional_asymmetry.pdf
#   sized 70 x 65 mm
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS02_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DATA_DIR <- file.path(OUT_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF  <- file.path(OUT_DIR, "fibrosis_stage_directional_asymmetry.pdf")

# ---------------------------------------------------------------------------
# Load DE results (one row per gene per stage contrast)
# ---------------------------------------------------------------------------
de <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv"))

# C2 canonical stage-vs-F0 DEG gate: analytical TREAT, fdr_treat < 0.05 at lfc=0.25.
# treat() tests H0:|true logFC|<=lfc, so the effect-size floor is folded INTO the test
# (no separate |logFC| filter). The per-stage CSV carries logFC + SE + t + P.Value but
# no treat_fdr, so we reconstruct the moderated-t TREAT statistic here (mirrors
# rebuild_cas13_library.R::add_treat_fdr, proven identical to limma::treat), inferring
# df.total per contrast from (t, P.Value). Direction from logFC sign.
TREAT_LFC <- 0.25
infer_df_total <- function(dt) {
  pr <- dt[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
  idx <- unique(round(seq(1, nrow(pr), length.out = min(nrow(pr), 12))))
  median(vapply(idx, function(i)
    uniroot(function(df) 2 * pt(-abs(pr$t[i]), df = df) - pr$P.Value[i], c(0.1, 1e6))$root,
    numeric(1)))
}
add_treat_fdr <- function(dt, lfc, se_col = "SE") {
  out <- copy(dt)
  se <- if (se_col %in% names(out)) out[[se_col]] else abs(out$logFC / out$t)
  se[!is.finite(se) | se <= 0] <- NA_real_
  df_use <- infer_df_total(out)
  out[, p_treat := pt((abs(logFC) - lfc) / se, df = df_use, lower.tail = FALSE) +
                   pt((abs(logFC) + lfc) / se, df = df_use, lower.tail = FALSE)]
  out[, fdr_treat := p.adjust(p_treat, method = "BH")]
  out
}
de <- de[!is.na(logFC) & !is.na(SE)]
de <- rbindlist(lapply(split(de, by = "contrast"), add_treat_fdr, lfc = TREAT_LFC))
sig <- de[fdr_treat < 0.05]
sig[, dir := ifelse(logFC > 0, "up", "down")]

# ---------------------------------------------------------------------------
# Per-stage up/down counts + ratio
# ---------------------------------------------------------------------------
tab <- sig[, .(
  n_up   = sum(dir == "up"),
  n_down = sum(dir == "down")
), by = contrast]
tab[, ratio := n_up / n_down]
tab[, stage := factor(contrast, levels = c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"),
                      labels = c("F1", "F2", "F3", "F4"))]
setorder(tab, stage)

# Under TREAT, F1 has 0 down DEGs -> ratio is undefined (Inf). Cap the BAR HEIGHT just
# above the max finite ratio so the bar still renders (induction-only = off the chart),
# and label it with the TRUE counts ("N:0") rather than fabricating a finite ratio.
fin_max <- max(tab[is.finite(ratio), ratio], na.rm = TRUE)
tab[, ratio_plot := ifelse(is.finite(ratio), ratio, fin_max * 1.12)]
tab[, ratio_lab  := ifelse(is.finite(ratio), sprintf("%.1f:1", ratio),
                           sprintf("%d:0", n_up))]
if (any(!is.finite(tab$ratio)))
  message(sprintf("CAPTION: %s has 0 down DEGs (pure induction); ratio undefined, bar capped at panel max, labelled with raw up:down counts.",
                  paste(tab[!is.finite(ratio), as.character(stage)], collapse = ", ")))

fwrite(tab, file.path(DATA_DIR, "fibrosis_stage_directional_asymmetry.csv"))

cat("[hero] up:down DEG ratio per stage (TREAT fdr<0.05, lfc=0.25):\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("    %s: %d up / %d down = %.2f:1\n",
              tab$stage[i], tab$n_up[i], tab$n_down[i], tab$ratio[i]))
}

# ---------------------------------------------------------------------------
# Bar chart (up:down ratio per stage)
# ---------------------------------------------------------------------------
message("[caption] Down-regulation overtakes induction as fibrosis advances")
p <- ggplot(tab, aes(x = stage, y = ratio_plot)) +
  geom_col(aes(fill = ratio_plot), width = 0.68) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#9E9E9E", linewidth = 0.3) +
  geom_text(aes(label = ratio_lab),
            vjust = -0.6, size = 6 / .pt, fontface = "plain") +
  scale_fill_gradient(low = "#C9265E", high = "#1565C0", guide = "none") +
  scale_y_continuous(name = "Up : down DEG ratio",
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Fibrosis stage (vs F0)") +
  theme_masld(base_size = 7) +
  theme(
    axis.text.x = element_text(size = 6, face = "plain")
  )

ggsave(OUT_PDF, p,
       width  = 70 / 25.4,
       height = 65 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
