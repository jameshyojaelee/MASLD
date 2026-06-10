##############################################################################
# Fig 4 Panel E: GWAS-ATAC drug-target convergence
#
# Two sub-panels (Liang convention — paired data types):
#
#   Sub-panel 1 (LEFT) -- 3-way Venn (HNF4A / THRB / RORA) of fine-mapped
#     GWAS variants that disrupt each TF's motif. Highlights convergent
#     variants where drug-target TFs (THRB) share disrupted peaks with
#     master hepatocyte TFs (HNF4A).
#
#   Sub-panel 2 (RIGHT) -- Forest plot of GWAS-eQTL COLOC PP.H4 for the
#     four FDA-approved / clinical-trial MASLD drug-target genes (THRB,
#     NR1H4 (FXR), PPARA, PPARG). THRB is Strong (PP4>0.99); the others
#     are expression-supported only. Drug-target highlight in gold.
#
# Data sources:
#   GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#
# Output:
#   figures/main/fig4_validation/panels/fig4_panel_E_drug_regulon_disruption.pdf
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

OUT_DIR <- file.path(BASE, "figures/main/fig4_validation/panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(OUT_DIR, "fig4_panel_E_drug_regulon_disruption.pdf")

# Liang aesthetic palette
COL_HNF4A   <- "#C9265E"   # warm magenta (discovery / master TF)
COL_THRB    <- "#FFB300"   # gold (drug target highlight)
COL_RORA    <- "#7B1FA2"   # warm purple
COL_STRONG  <- "#C9265E"   # magenta for Strong-tier
COL_EXPSUP  <- "#1565C0"   # cool blue for expression-supported
COL_REFLINE <- "#9E9E9E"   # neutral gray

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
GWAS_ATAC_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
SUSIE_DIR     <- file.path(BASE, "GWAS/finemapping/results/susie_coloc")

motif_file <- file.path(GWAS_ATAC_DIR, "motif_disruption_scores.csv")
coloc_file <- file.path(SUSIE_DIR, "gene_level_coloc.csv")

stopifnot(file.exists(motif_file), file.exists(coloc_file))

motif <- fread(motif_file)
coloc <- fread(coloc_file)

# Unique disrupting variants per TF (collapse multi-motif duplicates)
tf_keep <- c("HNF4A", "THRB", "RORA")
m_keep  <- unique(motif[tf_name %in% tf_keep, .(SNP_id, tf_name)])

hnf4a_set <- m_keep[tf_name == "HNF4A", unique(SNP_id)]
thrb_set  <- m_keep[tf_name == "THRB",  unique(SNP_id)]
rora_set  <- m_keep[tf_name == "RORA",  unique(SNP_id)]

# Pairwise + 3-way intersections
n_HNF4A_only        <- length(setdiff(hnf4a_set, union(thrb_set, rora_set)))
n_THRB_only         <- length(setdiff(thrb_set,  union(hnf4a_set, rora_set)))
n_RORA_only         <- length(setdiff(rora_set,  union(hnf4a_set, thrb_set)))
n_HNF4A_THRB        <- length(setdiff(intersect(hnf4a_set, thrb_set), rora_set))
n_HNF4A_RORA        <- length(setdiff(intersect(hnf4a_set, rora_set), thrb_set))
n_THRB_RORA         <- length(setdiff(intersect(thrb_set,  rora_set), hnf4a_set))
n_all3              <- length(Reduce(intersect, list(hnf4a_set, thrb_set, rora_set)))

n_HNF4A <- length(hnf4a_set)
n_THRB  <- length(thrb_set)
n_RORA  <- length(rora_set)

convergent_HT <- setdiff(intersect(hnf4a_set, thrb_set), rora_set)
convergent_TR <- setdiff(intersect(thrb_set, rora_set), hnf4a_set)

cat("Unique variants per TF (motif disruption, disease regulons):\n")
cat(sprintf("  HNF4A=%d  THRB=%d  RORA=%d\n", n_HNF4A, n_THRB, n_RORA))
cat(sprintf("  HNF4A&THRB=%d  HNF4A&RORA=%d  THRB&RORA=%d  all3=%d\n",
            n_HNF4A_THRB, n_HNF4A_RORA, n_THRB_RORA, n_all3))
cat("Convergent HNF4A&THRB variants:\n  ", paste(convergent_HT, collapse=", "), "\n")
cat("Convergent THRB&RORA variants:\n  ", paste(convergent_TR, collapse=", "), "\n")

# ---------------------------------------------------------------------------
# Sub-panel 1: Venn diagram (manual ggplot2 circles)
# ---------------------------------------------------------------------------
# 3-circle Venn — pre-computed circle centres for visual balance
venn_circles <- data.table(
  TF     = c("HNF4A", "THRB", "RORA"),
  x      = c(-0.55, 0.55, 0.00),
  y      = c( 0.32, 0.32, -0.55),
  r      = 0.85,
  fill   = c(COL_HNF4A, COL_THRB, COL_RORA)
)

# Build polygons for each circle (transparent overlapping fills)
make_circle <- function(cx, cy, r, n = 100) {
  th <- seq(0, 2 * pi, length.out = n)
  data.table(x = cx + r * cos(th), y = cy + r * sin(th))
}

circle_polys <- rbindlist(lapply(seq_len(nrow(venn_circles)), function(i) {
  p <- make_circle(venn_circles$x[i], venn_circles$y[i], venn_circles$r[i])
  p[, TF := venn_circles$TF[i]]
  p
}))

# Region labels (counts in each Venn region)
region_labels <- data.table(
  x = c(-0.95,  0.95,  0.00,
         0.00, -0.55,  0.55,
         0.00),
  y = c( 0.55,  0.55, -0.90,
         0.55, -0.40, -0.40,
         0.05),
  label = c(as.character(n_HNF4A_only),
            as.character(n_THRB_only),
            as.character(n_RORA_only),
            as.character(n_HNF4A_THRB),
            as.character(n_HNF4A_RORA),
            as.character(n_THRB_RORA),
            as.character(n_all3))
)

# TF labels outside circles
tf_labels <- data.table(
  x = c(-1.25,  1.25,  0.00),
  y = c( 1.05,  1.05, -1.45),
  label = c(sprintf("HNF4A\n(n=%d)", n_HNF4A),
            sprintf("THRB\n(n=%d)", n_THRB),
            sprintf("RORA\n(n=%d)", n_RORA)),
  colour = c(COL_HNF4A, COL_THRB, COL_RORA)
)

# Optional callout of convergent variants (HNF4A∩THRB) — drug-target rationale
callout_text <- if (length(convergent_HT) > 0 && length(convergent_HT) <= 3) {
  paste("HNF4A ∩ THRB:",
        paste(sub(":(.):.:.$", ":\\1>...", convergent_HT), collapse = ", "))
} else if (length(convergent_HT) > 0) {
  sprintf("HNF4A ∩ THRB: %d variants", length(convergent_HT))
} else {
  ""
}

p_venn <- ggplot() +
  geom_polygon(data = circle_polys,
               aes(x = x, y = y, group = TF, fill = TF),
               colour = "black", linewidth = 0.4, alpha = 0.32) +
  scale_fill_manual(values = c(HNF4A = COL_HNF4A, THRB = COL_THRB, RORA = COL_RORA),
                    guide = "none") +
  geom_text(data = region_labels, aes(x = x, y = y, label = label),
            size = 3.2, fontface = "bold", colour = "black") +
  geom_text(data = tf_labels, aes(x = x, y = y, label = label),
            colour = tf_labels$colour, fontface = "bold", size = 2.6,
            lineheight = 0.9) +
  annotate("text", x = 0, y = -1.85, label = callout_text,
           size = 2.1, fontface = "italic", colour = "gray20") +
  coord_equal(xlim = c(-1.7, 1.7), ylim = c(-2.0, 1.4)) +
  labs(title = "Fine-mapped variants disrupt\ndisease-regulon TF motifs",
       subtitle = "Unique variants per TF (disease-regulon-active motifs)") +
  theme_void(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold", hjust = 0.5,
                                     margin = margin(b = 2)),
        plot.subtitle = element_text(size = 6, hjust = 0.5, colour = "gray30",
                                     margin = margin(b = 6)),
        plot.margin   = margin(4, 4, 4, 4))

# ---------------------------------------------------------------------------
# Sub-panel 2: Forest plot of COLOC PP.H4 for MASLD drug-target genes
# ---------------------------------------------------------------------------
DRUGS <- data.table(
  gene      = c("THRB", "NR1H4", "PPARA", "PPARG"),
  drug      = c("Resmetirom (FDA 2024)",
                "Obeticholic acid",
                "Elafibranor",
                "Lanifibranor"),
  receptor  = c("THRβ", "FXR", "PPARα", "PPARγ"),
  tier      = c("Strong", "Expression-supported",
                "Expression-supported", "Expression-supported")
)

# Pull canonical coloc_best_pp4 (ABF) + susie when present.
coloc_sub <- coloc[gene %in% DRUGS$gene,
                   .(gene, coloc_best_pp4, coloc_best_gwas,
                     coloc_best_susie_pp4, coloc_best_susie_gwas)]
setkey(coloc_sub, gene)
setkey(DRUGS,    gene)
forest <- DRUGS[coloc_sub]

# Choose best PP.H4 across SuSiE + ABF (largest)
forest[, pp4_best := pmax(coloc_best_pp4, coloc_best_susie_pp4, na.rm = TRUE)]
forest[, pp4_gwas := ifelse(
          !is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 >= coloc_best_pp4,
          coloc_best_susie_gwas, coloc_best_gwas)]

# Order top -> bottom: Strong first (THRB), then others by PP4 desc
forest[, tier := factor(tier, levels = c("Strong", "Expression-supported"))]
forest <- forest[order(tier, -pp4_best)]
forest[, y := seq.int(.N, 1)]   # top=THRB

forest[, label_gene := sprintf("%s\n(%s)", gene, receptor)]

cat("\nDrug-target COLOC PP.H4 values:\n")
print(forest[, .(gene, drug, tier, pp4_best, pp4_gwas)])

p_forest <- ggplot(forest, aes(y = y)) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             colour = COL_REFLINE, linewidth = 0.4) +
  # hollow circle markers, thin black outline (Liang convention)
  geom_point(aes(x = pp4_best, fill = tier),
             shape = 21, colour = "black",
             size = 4, stroke = 0.5) +
  # gold ring on THRB to flag Strong drug-target
  geom_point(data = forest[gene == "THRB"],
             aes(x = pp4_best),
             shape = 1, colour = COL_THRB, size = 6, stroke = 0.9) +
  scale_fill_manual(values = c(Strong = COL_STRONG,
                               `Expression-supported` = COL_EXPSUP),
                    name = "Validation tier") +
  geom_text(aes(x = pp4_best, label = sprintf("%.3f", pp4_best)),
            hjust = -0.35, size = 2.4, fontface = "bold", colour = "black") +
  geom_text(aes(x = -0.04, label = drug),
            hjust = 1, size = 2.3, colour = "gray25") +
  scale_y_continuous(breaks = forest$y, labels = forest$label_gene,
                     limits = c(0.4, nrow(forest) + 0.6),
                     expand = expansion(add = 0)) +
  scale_x_continuous(limits = c(-0.55, 1.18),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                     expand = c(0, 0)) +
  annotate("text", x = 0.5, y = nrow(forest) + 0.45,
           label = "PP.H4 > 0.5",
           size = 1.9, colour = "gray35", fontface = "italic") +
  labs(title = "GWAS-eQTL colocalization at\nMASLD drug-target loci",
       subtitle = "Liver eQTL × GWAS (best across SuSiE + ABF, 23 GWAS portfolio)",
       x = "Colocalization PP.H4 (best GWAS)",
       y = NULL) +
  theme_masld(base_size = 7) +
  theme(plot.title       = element_text(size = 8, face = "bold", hjust = 0,
                                         margin = margin(b = 2)),
        plot.subtitle    = element_text(size = 6, hjust = 0, colour = "gray30",
                                         margin = margin(b = 6)),
        axis.title.x     = element_text(face = "bold", size = 7,
                                         margin = margin(t = 4)),
        axis.text.y      = element_text(face = "bold", size = 6.5,
                                         colour = "black", lineheight = 0.85),
        axis.text.x      = element_text(size = 6, colour = "black"),
        axis.line.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        legend.position  = c(0.78, 0.18),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.title     = element_text(size = 6, face = "bold"),
        legend.text      = element_text(size = 5.5),
        legend.key.size  = unit(0.25, "cm"),
        plot.margin      = margin(4, 6, 4, 6))

# ---------------------------------------------------------------------------
# Compose: side-by-side with shared header
# ---------------------------------------------------------------------------
header_title <- ggplot() +
  annotate("text", x = 0, y = 0,
           label = "GWAS-ATAC convergence on FDA-approved MASLD drug targets",
           hjust = 0.5, vjust = 0.5, fontface = "bold", size = 3.4) +
  theme_void() +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme(plot.margin = margin(2, 2, 2, 2))

composite <- header_title /
  (p_venn | p_forest) +
  plot_layout(heights = c(0.08, 1)) +
  plot_annotation(tag_levels = list(c("", "i", "ii"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(OUT_PDF, composite, width = 8, height = 4.2, device = pdf_device)

cat(sprintf("\nWrote: %s\n", OUT_PDF))
cat(sprintf("Size:  %.1f KB\n", file.info(OUT_PDF)$size / 1024))
