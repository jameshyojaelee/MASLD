#!/usr/bin/env Rscript

# Figure 4D: representative genes illustrating distinct evidence-class behavior.
# A binary evidence-signature heatmap (gene x modality): a cell is filled when the
# gene carries that evidence -- genetic colocalization (SuSiE or ABF PP4 > 0.5),
# canonical bulk DE (effect-size-aware interval-null FDR gate), significant liver
# protein change (padj < 0.05), significant single-cell pseudobulk change
# (padj < 0.05), and spatial variable expression (SVG). Genes are chosen DISTINCT
# from the Fig 5 vignettes (THRB/HKDC1/GLP1R/GCGR/RORA/MTARC1/PNPLA3) so this panel
# characterizes the classes rather than duplicating Fig 5's gene-level stories.
# Descriptive only; no ranking or target claim (that is Fig 5).

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RELEASE_ID  <- Sys.getenv("MANUSCRIPT_RELEASE_ID", "2026-07-15-r2")
RELEASE_DIR <- file.path(BASE, "RNA-seq/results/manuscript_release", RELEASE_ID)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

class_colors <- c(genetic_only = "#1565C0", disease_state_only = "#C9265E",
                  convergent = "#00695C")
class_labels <- c(genetic_only = "Genetic only", disease_state_only = "Disease-state only",
                  convergent = "Convergent")

ect <- fread(file.path(RELEASE_DIR, "evidence_class_table.tsv"))
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("human_symbol", "best_protein_padj", "sc_best_padj", "spatial_is_svg"))
setnames(atlas, "human_symbol", "symbol")
d <- merge(ect, atlas, by = "symbol", all.x = TRUE)

# per-axis binary evidence
d[, genetic_strength := pmax(max_susie_pp4_all, max_abf_pp4_all, na.rm = TRUE)]
d[!is.finite(genetic_strength), genetic_strength := NA_real_]
d[, ev_Genetic    := !is.na(genetic_strength) & genetic_strength > 0.5]
d[, ev_Bulk       := in_treat_deg == TRUE]
d[, ev_Proteomics := !is.na(best_protein_padj) & best_protein_padj < 0.05]
d[, `ev_Single-cell` := !is.na(sc_best_padj) & sc_best_padj < 0.05]
d[, `ev_Spatial SVG` := spatial_is_svg %in% TRUE]
d[, n_tissue := ev_Bulk + ev_Proteomics + `ev_Single-cell` + `ev_Spatial SVG`]

FIG5 <- c("THRB", "HKDC1", "GLP1R", "GCGR", "RORA", "MTARC1", "PNPLA3", "TM6SF2", "GPAM")
d <- d[!(symbol %in% FIG5) & joint_testable == TRUE & !is.na(symbol) & symbol != ""]

# representative gene selection (4 per class), by class-appropriate strength
g_gen <- d[primary_evidence_class == "genetic_only"][order(-genetic_strength)][1:4, symbol]
g_dis <- d[primary_evidence_class == "disease_state_only"][order(-n_tissue, bulk_treat_fdr)][1:4, symbol]
g_con <- d[primary_evidence_class == "convergent"][order(-genetic_strength, -n_tissue)][1:4, symbol]
sel <- d[symbol %in% c(g_gen, g_dis, g_con)]

ax <- c("Genetic", "Bulk", "Proteomics", "Single-cell", "Spatial SVG")
long <- melt(sel,
             id.vars = c("symbol", "primary_evidence_class"),
             measure.vars = paste0("ev_", ax),
             variable.name = "modality", value.name = "present")
long[, modality := factor(sub("^ev_", "", modality), levels = ax)]
long[, class_label := factor(class_labels[primary_evidence_class],
                             levels = class_labels[c("genetic_only", "disease_state_only", "convergent")])]
# order genes within class by evidence signature
gene_order <- c(g_gen, g_dis, g_con)
long[, symbol := factor(symbol, levels = rev(gene_order))]
long[, fill_key := fifelse(present, primary_evidence_class, "absent")]

fill_vals <- c(class_colors, absent = "#ECEFF1")

p4d <- ggplot(long, aes(x = modality, y = symbol, fill = fill_key)) +
  geom_tile(color = "white", linewidth = 0.6) +
  facet_grid(class_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = fill_vals, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 30, hjust = 0),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(FIG4_DIR, "fig4d_representative_genes.pdf"), p4d,
       width = 4.6, height = 3.2, device = cairo_pdf, bg = "white")

cat(sprintf("[release %s] Figure 4D representative genes written\n", RELEASE_ID))
cat("  genetic-only:", paste(g_gen, collapse = ", "), "\n")
cat("  disease-state-only:", paste(g_dis, collapse = ", "), "\n")
cat("  convergent:", paste(g_con, collapse = ", "), "\n")
