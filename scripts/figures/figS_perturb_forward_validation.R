#!/usr/bin/env Rscript
# Forward-validation TEST 1 scatter:
#   x = atlas convergence_score
#   y = hepatocyte disease-axis perturbation magnitude (|effect|)
#   points = Saunders 2025 panel genes, colored by atlas tier
#   labels = top convergent perturbed genes
#   annotation = Spearman rho (raw), permutation p, partial rho
# House style: theme_masld, PDF, cairo, useDingbats=FALSE-equivalent.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(ROOT, "scripts/figures/publication_theme.R"))

OUTDIR <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/perturb_forward_validation")
csv  <- file.path(OUTDIR, "saunders_forward_validation_per_gene.csv")
sj   <- file.path(OUTDIR, "forward_validation_stats.json")
stopifnot(file.exists(csv), file.exists(sj))

d <- read.csv(csv, stringsAsFactors = FALSE)
st <- fromJSON(sj)

# keep genes with both axes
d <- d[!is.na(d$convergence_score) & !is.na(d$abs_disease_axis_effect), ]

# tier color mapping (atlas 4-tier scheme present in overlap CSV)
tier_pal <- c(
  "1_Genetic_validated" = "#880E4F",
  "2_Strong"            = "#C2185B",
  "3_Suggestive"        = "#E91E63",
  "4_Weak"              = "#90A4AE",
  "Excluded"            = "#CFD8DC"
)
d$tier[is.na(d$tier) | d$tier == ""] <- "Excluded"
d$tier <- factor(d$tier, levels = names(tier_pal))

# label the most convergent perturbed genes (top by convergence_score) + a few
# with the strongest disease-axis effect, so the labels carry signal.
n_lab_conv <- 12
lab_conv <- head(d[order(-d$convergence_score), "human_ortholog"], n_lab_conv)
lab_eff  <- head(d[order(-d$abs_disease_axis_effect), "human_ortholog"], 6)
d$label <- ifelse(d$human_ortholog %in% unique(c(lab_conv, lab_eff)),
                  d$human_ortholog, NA)

ann <- sprintf(
  "Spearman rho = %.2f  (perm p = %.3f)\nPartial rho (ess + log n) = %.2f  (p = %.3f)\nn = %d perturbed genes",
  st$rho_raw, st$perm_p_raw, st$rho_partial_ess_ncells, st$perm_p_partial, st$n_in_core)

x_rng <- range(d$convergence_score, na.rm = TRUE)
y_rng <- range(d$abs_disease_axis_effect, na.rm = TRUE)

p <- ggplot(d, aes(convergence_score, abs_disease_axis_effect)) +
  geom_smooth(method = "lm", se = TRUE, color = "#37474F",
              fill = "#ECEFF1", linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(fill = tier, size = n_cells), shape = 21,
             color = "white", stroke = 0.25, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = label), size = GEOM_TEXT_6PT, color = "#212121",
                           segment.size = 0.2, segment.color = "#9E9E9E",
                           min.segment.length = 0, max.overlaps = 30,
                           box.padding = 0.3, na.rm = TRUE) +
  scale_fill_manual(values = tier_pal, name = "Atlas tier", drop = FALSE) +
  scale_size_continuous(range = c(0.6, 3.2), name = "n cells") +
  annotate("text", x = x_rng[1], y = y_rng[2],
           label = ann, hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, color = "#263238") +
  labs(
    x = "Atlas convergence score",
    y = "Hepatocyte disease-axis perturbation magnitude  |effect| (NC-standardized)") +
  theme_masld() +
  theme(legend.position = "right")

message("[caption] Forward validation: atlas convergence vs. Saunders 2025 in-vivo CRISPRi effect (independent mouse hepatocyte Perturb-seq, not designed around the atlas)")
out_pdf <- file.path(ROOT, "figures/supplementary/figS_perturb_forward_validation.pdf")
save_fig(p, out_pdf, width = fig_full_width, height = 4.6)
cat("Wrote figure:", out_pdf, "\n")
