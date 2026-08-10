#!/usr/bin/env Rscript

# Figure 4E: native-PCC gsMap spatial GWAS-risk, evidence-class comparison
# (boundary/annotation panel). For each direct-disease / PDFF trait in two Visium
# cohorts, the matched (expression + GSS quintile) convergent-minus-single-map spatial
# GWAS-risk PCC difference. Differences are small and straddle zero; NO trait shows a
# significant convergent advantage replicated in both cohorts (release headline
# gsmap_class_validation headline_pass = FALSE). Inherited spatial risk is NOT
# exclusive to the convergent class -- shown WITHOUT an exclusivity claim.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RELEASE_ID  <- Sys.getenv("MANUSCRIPT_RELEASE_ID", "2026-07-15-r2")
RELEASE_DIR <- file.path(BASE, "RNA-seq/results/manuscript_release", RELEASE_ID)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

v <- fread(file.path(RELEASE_DIR, "gsmap_class_validation.tsv"))
v <- v[trait_scope == "direct_disease"]
v[, comp_label := factor(
  fifelse(comparison == "convergent_vs_disease_state_only",
          "vs disease-state-only", "vs genetic-only"),
  levels = c("vs genetic-only", "vs disease-state-only"))]
v[, cohort_label := fifelse(cohort == "vu", "Vu", "GSE192741")]

# order traits by mean difference for a readable forest
ord <- v[, .(m = mean(mean_pcc_difference)), by = trait][order(m), trait]
v[, trait := factor(trait, levels = ord)]

p4e <- ggplot(v, aes(x = mean_pcc_difference, y = trait, color = cohort_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_point(size = 1.4, alpha = 0.9) +
  facet_wrap(~ comp_label) +
  scale_color_manual(values = c("GSE192741" = "#5E35B1", "Vu" = "#EF6C00"), name = NULL) +
  labs(x = "Matched convergent − single-map spatial-risk PCC difference", y = NULL) +
  theme_masld(base_size = 6) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(FIG4_DIR, "fig4e_gsmap_class_pcc.pdf"), p4e,
       width = 5.6, height = 3.4, device = cairo_pdf, bg = "white")

n_sig <- v[positive_significant == TRUE, .N]
cat(sprintf("[release %s] Figure 4E gsMap native-PCC written; %d/%d direct-disease trait-cohort-comparisons significant\n",
            RELEASE_ID, n_sig, nrow(v)))
