#!/usr/bin/env Rscript
# DIRECTIONAL forward-validation panel (TEST 1b).
#   directional_reversal_g = -sign(atlas_disease_logFC) * signed_disease_axis_effect
#   POSITIVE = CRISPRi knockdown reversed g's disease direction (expected for a true driver)
# Panels:
#   A: directional_reversal vs convergence_score (genes with a clear atlas direction)
#   B: high- vs low-convergence directional boxplot
#   C: disease-UP subset one-sample distribution (is knockdown reversing disease?)
# House style: theme_masld, cairo PDF.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(patchwork)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(ROOT, "scripts/figures/publication_theme.R"))

OUTDIR <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/perturb_forward_validation")
d  <- read.csv(file.path(OUTDIR, "saunders_forward_validation_per_gene.csv"),
               stringsAsFactors = FALSE)
st <- fromJSON(file.path(OUTDIR, "forward_validation_stats.json"))
dr <- st$directional

tier_pal <- c("1_Genetic_validated" = "#880E4F", "2_Strong" = "#C2185B",
              "3_Suggestive" = "#E91E63", "4_Weak" = "#90A4AE", "Excluded" = "#CFD8DC")
d$tier[is.na(d$tier) | d$tier == ""] <- "Excluded"
d$tier <- factor(d$tier, levels = names(tier_pal))

clear <- d[!is.na(d$directional_reversal) & d$has_clear_dir == "True", ]
if (nrow(clear) == 0) clear <- d[!is.na(d$directional_reversal) & d$has_clear_dir == TRUE, ]

# ---- Panel A: directional_reversal vs convergence_score ----
bA <- clear[!is.na(clear$convergence_score), ]
lab <- head(bA[order(-bA$convergence_score), "human_ortholog"], 10)
lab <- unique(c(lab, head(bA[order(-abs(bA$directional_reversal)), "human_ortholog"], 5)))
bA$label <- ifelse(bA$human_ortholog %in% lab, bA$human_ortholog, NA)
annA <- sprintf("raw rho = %+.2f (p = %.3f)\npartial rho = %+.2f (p = %.3f)\nn = %d clear-direction genes",
                dr$spearman_raw_rho, dr$spearman_raw_perm_p,
                dr$spearman_partial_rho, dr$spearman_partial_perm_p, dr$spearman_raw_n)

pA <- ggplot(bA, aes(convergence_score, directional_reversal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#B0BEC5", linewidth = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "#37474F", fill = "#ECEFF1",
              linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(fill = tier, size = n_cells), shape = 21, color = "white",
             stroke = 0.25, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = label), size = 1.9, color = "#212121",
                           segment.size = 0.2, segment.color = "#9E9E9E",
                           min.segment.length = 0, max.overlaps = 30, na.rm = TRUE) +
  scale_fill_manual(values = tier_pal, name = "Atlas tier", drop = FALSE) +
  scale_size_continuous(range = c(0.6, 3.0), name = "n cells") +
  annotate("text", x = min(bA$convergence_score), y = max(bA$directional_reversal),
           label = annA, hjust = 0, vjust = 1, size = 1.9, color = "#263238") +
  labs(x = "Atlas convergence score",
       y = "Directional reversal\n(+ = knockdown reverses disease)",
       title = "A  Directional reversal vs convergence") +
  theme_masld() + theme(plot.title = element_text(size = 7, face = "bold"))

# ---- Panel B: high vs low convergence directional boxplot ----
bB <- clear[clear$conv_group %in% c("high", "low"), ]
bB$conv_group <- factor(bB$conv_group, levels = c("low", "high"))
annB <- sprintf("MWU high>low\np = %.3f", dr$mwu_p_greater)
pB <- ggplot(bB, aes(conv_group, directional_reversal, fill = conv_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#B0BEC5", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.55, linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 0.5, alpha = 0.5, color = "#37474F") +
  scale_fill_manual(values = c(low = "#90A4AE", high = "#C2185B"), guide = "none") +
  annotate("text", x = 1.5, y = max(bB$directional_reversal, na.rm = TRUE),
           label = annB, size = 2.0, color = "#263238", vjust = 1) +
  labs(x = "Convergence group", y = "Directional reversal",
       title = "B  High vs low convergence") +
  theme_masld() + theme(plot.title = element_text(size = 7, face = "bold"))

# ---- Panel C: disease-UP subset one-sample distribution ----
up <- clear[clear$atlas_dir > 0, ]
medUP <- median(up$directional_reversal, na.rm = TRUE)
annC <- sprintf("disease-UP genes (cleanest)\nn=%d  median=%+.3f\nsign p(>0)=%.2f  Wilcoxon p(>0)=%.2f",
                dr$onesample_disease_up$n, dr$onesample_disease_up$median,
                dr$onesample_disease_up$sign_p_greater, dr$onesample_disease_up$wilcoxon_p_greater)
pC <- ggplot(up, aes(directional_reversal)) +
  geom_histogram(bins = 22, fill = "#90A4AE", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#37474F", linewidth = 0.4) +
  geom_vline(xintercept = medUP, color = "#C2185B", linewidth = 0.5) +
  annotate("text", x = min(up$directional_reversal), y = Inf, label = annC,
           hjust = 0, vjust = 1.3, size = 1.9, color = "#263238") +
  labs(x = "Directional reversal (disease-UP genes)", y = "Perturbed genes",
       title = "C  Knockdown of disease-UP genes") +
  theme_masld() + theme(plot.title = element_text(size = 7, face = "bold"))

p <- (pA | pB | pC) + plot_layout(widths = c(1.5, 0.8, 1.1)) +
  plot_annotation(
    title = "Directional forward validation: atlas convergence vs Saunders 2025 in-vivo CRISPRi",
    subtitle = "Positive reversal = CRISPRi knockdown moved hepatocytes OPPOSITE to the gene's atlas disease direction",
    theme = theme(plot.title = element_text(size = 7.5, face = "bold"),
                  plot.subtitle = element_text(size = 6, color = "#546E7A")))

out_pdf <- file.path(ROOT, "figures/supplementary/figS_perturb_forward_validation_directional.pdf")
save_fig(p, out_pdf, width = fig_full_width, height = 3.4)
cat("Wrote figure:", out_pdf, "\n")
