#!/usr/bin/env Rscript
# figS_cas13_library_lfc_vs_tpm.R
# ---------------------------------------------------------------------------
# Final figure (10): how to reach a ~1,900-gene Cas13 KD library.
#   LEFT  : library size vs human ashr LFC cutoff, with vs without the
#           biotype-aware TPM gate (PC>=1.0 / lncRNA>=0.5) -> the "TPM cutoff
#           option". Marks the two chosen libraries A and B.
#   RIGHT : head-to-head of the two ~1,900-gene libraries (median disease-liver
#           TPM and median human shrunk LFC).
#     A = lower LFC (ashr>0.20) + split TPM gate
#     B = higher LFC (ashr>0.40), no TPM gate
# Sources: Cas13_Library_Design/data/lfc_vs_tpm_size_sweep.csv
#          Cas13_Library_Design/data/lfc_vs_tpm_headtohead.csv
#          (built by Cas13_Library_Design/scripts/compare_lfc_vs_tpm_strategy.R)
#
# Output: 10_lfc_vs_tpm_strategy.pdf -> Cas13_Library_Design/figures/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT_DIR    <- FIGS_CAS13LIB_DIR
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

DATD <- file.path(BASE, "Cas13_Library_Design/data")
TPM_REF <- Sys.getenv("TPM_REF", "all")
SUF     <- if (TPM_REF == "HFD") "_hfd" else ""
REF_LAB <- if (TPM_REF == "HFD") "HFD diet only" else "all 4 DEG diet groups"
sweep <- fread(file.path(DATD, paste0("lfc_vs_tpm_size_sweep",  SUF, ".csv")))
cmp   <- fread(file.path(DATD, paste0("lfc_vs_tpm_headtohead", SUF, ".csv")))
cmp[, strat := ifelse(grepl("no TPM", strategy), "B", "A")]

# chosen thresholds parsed from the strategy labels
thr_of <- function(lbl) as.numeric(sub(".*ashr>([0-9.]+).*", "\\1", lbl))
thr_A <- thr_of(cmp[strat == "A", strategy]); thr_B <- thr_of(cmp[strat == "B", strategy])
nA <- cmp[strat == "A", n_genes];             nB <- cmp[strat == "B", n_genes]

pal       <- c(A = "#1565C0", B = "#9E9E9E")
gate_pal  <- c("Split TPM gate (PC≥1.0 / lncRNA≥0.5)" = "#1565C0", "No TPM gate" = "#9E9E9E")
strat_lab <- c(A = sub("^A: ", "A  ", cmp[strat == "A", strategy]),
               B = sub("^B: ", "B  ", cmp[strat == "B", strategy]))

# ---- LEFT: size sweep (the TPM-cutoff option) ------------------------------
sw <- melt(sweep[, .(thr, `No TPM gate` = n_lib,
                     `Split TPM gate (PC≥1.0 / lncRNA≥0.5)` = n_lib_tpm1)],
           id.vars = "thr", variable.name = "gate", value.name = "n")
sw[, gate := factor(gate, levels = names(gate_pal))]

# the two chosen libraries as labelled points
picks <- data.table(
  thr   = c(thr_A, thr_B),
  n     = c(nA, nB),
  gate  = factor(c("Split TPM gate (PC≥1.0 / lncRNA≥0.5)", "No TPM gate"), levels = names(gate_pal)),
  label = c(sprintf("A (%s)", format(nA, big.mark = ",")),
            sprintf("B (%s)", format(nB, big.mark = ","))))

p_sweep <- ggplot(sw, aes(thr, n, color = gate)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  geom_point(data = picks, aes(thr, n), color = "black", size = 2.4, shape = 21,
             fill = c("#1565C0", "#9E9E9E"), stroke = 0.7, inherit.aes = FALSE) +
  geom_text(data = picks, aes(thr, n, label = label), color = "black",
            vjust = -1.1, size = 2.6, fontface = "bold", inherit.aes = FALSE) +
  scale_color_manual(values = gate_pal, name = NULL) +
  scale_x_continuous(name = "Human ashr shrunk-LFC cutoff", breaks = sweep$thr,
                     expand = expansion(mult = c(0.08, 0.03))) +
  scale_y_continuous(name = "Library size (genes)",
                     expand = expansion(mult = c(0.04, 0.15)),
                     labels = scales::comma) +
  labs(title = "TPM gate across LFC cutoffs") +
  theme_masld() + theme_pub() +
  theme(legend.position = "top", legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.32, "cm"),
        plot.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 6),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.22))

# ---- RIGHT: A vs B head-to-head bars ---------------------------------------
specs <- data.table(
  metric = c("median_TPM", "median_human_LFC"),
  label  = c("Median disease\nliver TPM", "Median human\nshrunk LFC"),
  fmt    = c("%.2f", "%.3f"))
specs[, label := factor(label, levels = label)]
long <- melt(cmp[, .(strat, median_TPM, median_human_LFC)],
             id.vars = "strat", variable.name = "metric", value.name = "value")
long <- merge(long, specs, by = "metric"); long[, label := factor(label, levels = specs$label)]
long[, txt := sprintf(fmt, value)]

p_bars <- ggplot(long, aes(strat, value, fill = strat)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = txt), vjust = -0.35, size = 2.5, fontface = "bold") +
  facet_wrap(~ label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal, labels = strat_lab, name = NULL) +
  scale_x_discrete(labels = c(A = "A", B = "B")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = NULL, title = "A vs B") +
  theme_masld() + theme_pub() +
  theme(legend.position = "top", legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 7, face = "bold", lineheight = 0.95),
        strip.background = element_rect(fill = "grey94", color = NA),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.22),
        plot.title = element_text(size = 8, face = "bold"))

fig <- p_sweep + p_bars + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    title = "LFC cutoff + TPM gate + GWAS tier",
    theme = theme(plot.title = element_text(size = 9, face = "bold")))

# caption detail (was in the on-figure title) -> stdout for the manuscript caption
cat(sprintf("TPM reference: %s\n", REF_LAB))

ggsave(file.path(OUT_DIR, paste0("lfc_vs_tpm_strategy", SUF, ".pdf")), fig,
       width = 8.4, height = 3.7, device = pdf_device)
cat(sprintf("Saved: 10_lfc_vs_tpm_strategy%s.pdf\n", SUF))
