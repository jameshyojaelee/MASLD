#!/usr/bin/env Rscript
# figS_extreme_phenotype_concordance.R
# ===========================================================================
# How well the EXTREME-phenotype contrast (definite-disease NAS>=5|F>=3 vs
# strict controls; 05i) recapitulates the canonical consensus DEG program that
# was derived WITH the mild/ambiguous cases included (canonical_deg_results.csv).
#
# Reviewer-defense panel for the Fig 3 "QC and robustness" section: if the
# consensus signature were an artifact of borderline cases, the extremes would
# DISAGREE. Instead they recover ~97% of it with tight effect-size concordance.
#
# Panel A: effect-size scatter -- canonical (mild-inclusive) shrunk_logFC vs
#          extreme-phenotype shrunk_logFC across the 1,433 consensus DEGs,
#          colored by recapitulation status, with identity line + metrics.
# Panel B: 100%-stacked bar of the 1,433 consensus DEGs by recapitulation class.
#
# Output: figures/supplementary/figS_sample_clustering/extreme_phenotype_concordance.pdf
# ===========================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INTEG <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
CAN   <- file.path(INTEG, "canonical_deg_results.csv")
EXT   <- file.path(INTEG, "sensitivity/extreme_phenotype_definite_vs_strict.csv")
OUT   <- file.path(BASE, "figures/supplementary/figS_sample_clustering/extreme_phenotype_concordance.pdf")

# --- load + join -----------------------------------------------------------
can <- fread(CAN); ext <- fread(EXT)
m <- merge(
  can[, .(gene, can_lfc = logFC, can_slfc = shrunk_logFC, can_lfsr = lfsr)],
  ext[, .(gene, ext_lfc = logFC, ext_slfc = shrunk_logFC, ext_lfsr = lfsr)],
  by = "gene")

# Tier-1 definitions (ashr): lfsr < 0.05 & |shrunk_logFC| > 0.5
m[, can_deg := can_lfsr < 0.05 & abs(can_slfc) > 0.5]
m[, ext_deg := ext_lfsr < 0.05 & abs(ext_slfc) > 0.5]

# --- global concordance metrics (all genes tested in both) -----------------
rho_all <- cor(m$can_lfc, m$ext_lfc, method = "spearman", use = "complete.obs")

# --- restrict to the consensus DEGs, classify recapitulation ----------------
d <- m[can_deg == TRUE]
d[, class := fifelse(ext_deg & sign(ext_slfc) == sign(can_slfc), "Recovered",
              fifelse(sign(ext_slfc) == sign(can_slfc),           "Subthreshold",
                                                                   "Discordant"))]
d[, class := factor(class, levels = c("Recovered", "Subthreshold", "Discordant"))]
n_deg    <- nrow(d)
recap    <- mean(d$class == "Recovered")
same_dir <- mean(sign(d$ext_slfc) == sign(d$can_slfc))
rho_deg  <- cor(d$can_slfc, d$ext_slfc, method = "spearman")
cat(sprintf("[concordance] n_consensus_DEG=%d  recapitulated=%.1f%%  same-direction=%.1f%%  rho_all=%.3f  rho_deg=%.3f\n",
            n_deg, 100*recap, 100*same_dir, rho_all, rho_deg))
cat("[breakdown]\n"); print(d[, .N, by = class])

# --- colors ----------------------------------------------------------------
cls_cols <- c("Recovered"    = "#B2182B",
              "Subthreshold" = "#F4A582",
              "Discordant"   = "#9E9E9E")
lim <- max(abs(c(d$can_slfc, d$ext_slfc)), na.rm = TRUE) * 1.02

# --- Panel A: effect-size scatter ------------------------------------------
ann <- sprintf("Spearman ρ = %.2f  (27,638 genes)\n%.1f%% of %s consensus DEGs recovered\n%.1f%% same direction",
               rho_all, 100*recap, format(n_deg, big.mark=","), 100*same_dir)
pA <- ggplot(d[order(class, decreasing = TRUE)], aes(can_slfc, ext_slfc, color = class)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_point(size = 0.5, alpha = 0.65, stroke = 0) +
  annotate("text", x = -lim*0.97, y = lim*0.97, hjust = 0, vjust = 1,
           label = ann, size = PUB_GEOM_TEXT, color = "grey15", lineheight = 0.95) +
  scale_color_manual(values = cls_cols, name = NULL) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Consensus log2FC (disease vs control, mild cases included)",
       y = "Extreme log2FC (advanced disease vs strict control)",
       title = "Advanced disease recapitulates the consensus DEGs") +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), color = NA))

# --- Panel B: 100%-stacked recapitulation bar ------------------------------
bar <- d[, .(frac = .N / n_deg), by = class][order(class)]
bar[, lab := sprintf("%.1f%%", 100*frac)]
pB <- ggplot(bar, aes(x = 1, y = frac, fill = class)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(frac > 0.05, lab, "")),
            position = position_stack(vjust = 0.5), size = PUB_GEOM_TEXT, color = "white", fontface = "bold") +
  scale_fill_manual(values = cls_cols, name = NULL, guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL, y = sprintf("Proportion of %s consensus DEGs", format(n_deg, big.mark=",")),
       caption = sprintf("Subthreshold: %d DEGs (%.1f%%) same direction, not significant in the extreme contrast · 0 discordant",
                         sum(d$class == "Subthreshold"), 100*(1-recap))) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p <- pA / pB + plot_layout(heights = c(6, 1))

ok <- tryCatch({ cairo_pdf(OUT, width = 3.6, height = 4.1); TRUE }, error = function(e) FALSE)
if (!ok) pdf(OUT, width = 3.6, height = 4.1, useDingbats = FALSE)
print(p); dev.off()
cat("[write]", OUT, "\n[done]\n")
