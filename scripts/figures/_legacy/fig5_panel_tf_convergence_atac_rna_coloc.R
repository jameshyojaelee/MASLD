# =============================================================================
# Fig 5 (Convergence) — TF convergence scatter
#
# 2D scatter anchoring 4 evidence layers on the SAME TF axis:
#   x    = SCENIC+ hepatocyte regulon activity Δ (disease − healthy)
#   y    = C2-canonical logFC of the TF gene itself (disease vs control;
#          canonical_deg_results.csv)
#   size = # variants disrupting that TF's motif at 2-of-3 concordance
#          (allele_concordance_summary.csv :: n_concord_2of3)
#   color = max COLOC PP.H4 across 23-GWAS portfolio
#           (gene_level_coloc.csv :: max(coloc_best_pp4, coloc_best_susie_pp4))
#
# Output: figures/main/fig5_convergence/panels/fig5_tf_convergence_scatter.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE_DIR, "scripts/figures/load_figure_data.R"))
source(file.path(BASE_DIR, "scripts/figures/publication_theme.R"))

OUT_DIR <- file.path(FIG5_DIR, "panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. SCENIC+ hepatocyte regulons (primary x-axis)
# -----------------------------------------------------------------------------
hep_reg_path <- file.path(ATAC_DIR, "scenic_plus", "hepatocyte_regulons.csv")
hep <- fread(hep_reg_path)
# one row per (tf, target); collapse to unique TFs
hep_u <- unique(hep[, .(tf_name, regulon_activity_diff, activity_padj,
                        n_target_genes, n_enhancers)])
setnames(hep_u, "regulon_activity_diff", "activity_diff")
cat(sprintf("SCENIC+ hepatocyte regulons: %d unique TFs\n", nrow(hep_u)))

# -----------------------------------------------------------------------------
# 2. C2-canonical TF logFC (y-axis) — canonical_deg_results.csv
# -----------------------------------------------------------------------------
deg_path <- file.path(INT_RESULTS, "canonical_deg_results.csv")
deg <- fread(deg_path, select = c("symbol", "logFC", "padj"))
# Some symbols missing — drop those
deg <- deg[!is.na(symbol) & symbol != ""]
# Average duplicate symbols (rare)
deg <- deg[, .(logFC = mean(logFC, na.rm = TRUE),
                   padj  = min(padj,  na.rm = TRUE)), by = symbol]
setnames(deg, "symbol", "tf_name")

# -----------------------------------------------------------------------------
# 3. Motif disruption (size) — allele_concordance_summary.csv
# -----------------------------------------------------------------------------
mot_path <- file.path(BASE_DIR,
  "GWAS/finemapping/results/gwas_atac/allele_concordance_summary.csv")
mot <- fread(mot_path,
  select = c("tf_name", "n_concord_2of3", "n_concord_3of3",
             "is_activator", "in_disease_regulon"))
# allele_concordance may carry compound motifs (e.g. "PPARA::RXRA");
# split to individual TFs and take the MAX hit count per TF
expand_motif <- function(dt) {
  rows <- dt[, .(tf = unlist(strsplit(tf_name, "::"))), by = .(seq_len(nrow(dt)))]
  setnames(rows, "tf", "tf_single")
  rows <- merge(rows, dt[, .SD, .SDcols = c("tf_name", setdiff(names(dt), "tf_name"))],
                by.x = "seq_len", by.y = NULL, allow.cartesian = TRUE)
  rows
}
mot_split <- mot[, .(tf_single = unlist(strsplit(tf_name, "::"))),
                 by = .(tf_name, n_concord_2of3, n_concord_3of3,
                        is_activator, in_disease_regulon)]
mot_tf <- mot_split[, .(n_concord_2of3 = max(n_concord_2of3, na.rm = TRUE),
                        n_concord_3of3 = max(n_concord_3of3, na.rm = TRUE),
                        in_disease_regulon = any(in_disease_regulon)),
                    by = tf_single]
setnames(mot_tf, "tf_single", "tf_name")

# -----------------------------------------------------------------------------
# 4. COLOC PP4 (color) — gene_level_coloc.csv
# -----------------------------------------------------------------------------
coloc_path <- file.path(BASE_DIR,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
coloc <- fread(coloc_path,
  select = c("gene", "coloc_best_pp4", "coloc_best_susie_pp4",
             "coloc_n_gwas_h4_05", "coloc_n_gwas_susie_h4_05"))
coloc <- coloc[!is.na(gene) & gene != ""]
coloc[, max_pp4 := pmax(coloc_best_pp4, coloc_best_susie_pp4, na.rm = TRUE)]
coloc[is.na(max_pp4) | !is.finite(max_pp4), max_pp4 := 0]
coloc[, n_gwas_h4 := pmax(coloc_n_gwas_h4_05,
                          coloc_n_gwas_susie_h4_05, na.rm = TRUE)]
coloc[is.na(n_gwas_h4), n_gwas_h4 := 0]
setnames(coloc, "gene", "tf_name")

# -----------------------------------------------------------------------------
# 5. Merge — keep ALL hepatocyte SCENIC+ TFs as the universe
# -----------------------------------------------------------------------------
df <- merge(hep_u, deg[, .(tf_name, logFC, padj)], by = "tf_name", all.x = TRUE)
df <- merge(df, mot_tf, by = "tf_name", all.x = TRUE)
df <- merge(df, coloc[, .(tf_name, max_pp4, n_gwas_h4)],
            by = "tf_name", all.x = TRUE)

# Fill missing
df[is.na(n_concord_2of3), n_concord_2of3 := 0]
df[is.na(n_concord_3of3), n_concord_3of3 := 0]
df[is.na(max_pp4),        max_pp4        := 0]
df[is.na(n_gwas_h4),      n_gwas_h4      := 0]
df[is.na(logFC),          logFC          := 0]
df[is.na(padj),           padj           := 1]

# -----------------------------------------------------------------------------
# 6. Annotation flags
# -----------------------------------------------------------------------------
HEADLINE_TFS <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")
DRUG_TARGETS <- c("THRB", "NR1H4", "PPARA", "PPARG")  # FDA-approved/clinical
LABEL_TFS    <- unique(c(HEADLINE_TFS, DRUG_TARGETS,
                         "FOXA1", "FOXA2", "CEBPB", "ESR1", "KLF15", "RXRA",
                         "HNF1A", "AHR", "XBP1", "SMAD4"))

df[, is_headline    := tf_name %in% HEADLINE_TFS]
df[, is_drug_target := tf_name %in% DRUG_TARGETS]
df[, coloc_sig      := max_pp4 >= 0.5]
df[, deg_sig        := !is.na(padj) & padj < 0.05]

# "Combined evidence score" for ranking labels
df[, evidence_score := abs(activity_diff) * (1 + log1p(n_concord_2of3)) *
                       (1 + 2 * max_pp4)]

# Always label headline TFs + drug targets + top-evidence
df[, do_label := tf_name %in% LABEL_TFS]
top_extra <- df[do_label == FALSE][order(-evidence_score)][1:5, tf_name]
df[tf_name %in% top_extra, do_label := TRUE]

# -----------------------------------------------------------------------------
# 7. Quadrant + region summary (for caption)
# -----------------------------------------------------------------------------
cat("\n=== Quadrant occupancy ===\n")
df[, quadrant := fcase(
  activity_diff >  0 & logFC >  0, "Q1 activated + DEG-up",
  activity_diff <  0 & logFC >  0, "Q2 suppressed regulon, DEG-up",
  activity_diff <  0 & logFC <  0, "Q3 suppressed + DEG-down",
  activity_diff >  0 & logFC <  0, "Q4 activated regulon, DEG-down",
  default = "axis"
)]
print(df[, .N, by = quadrant])
cat("\n=== Headline TFs ===\n")
print(df[is_headline == TRUE,
         .(tf_name, activity_diff, logFC, padj,
           n_concord_2of3, max_pp4, n_gwas_h4)])

# -----------------------------------------------------------------------------
# 8. Plot — Liang aesthetic (2026-05-18)
#   - Hollow circles (no fill, thin black outline) for ALL points
#   - Color of OUTLINE encodes COLOC PP4 (gray <0.5 → warm magenta)
#   - Thick gold outer RING for FDA drug targets (under main points)
#   - Italic gray quadrant labels in plot corners
# -----------------------------------------------------------------------------
breaks_pp4 <- c(0, 0.25, 0.5, 0.75, 1.0)

# Map color of point outline to COLOC PP4 (gray under 0.5, magenta above)
df[, outline_color := {
  cr <- colorRamp(c("#D9D9D9", "#BDBDBD", "#E8B6CA", "#C9265E", "#7A1140"))
  rgb_mat <- cr(scales::rescale(pmin(pmax(max_pp4, 0), 1),
                                from = c(0, 1)))
  rgb(rgb_mat[, 1], rgb_mat[, 2], rgb_mat[, 3], maxColorValue = 255)
}]

# Drop CYP26A1 from labels (not a TF; gene-level, lives in Panel K)
df[tf_name == "CYP26A1", do_label := FALSE]

# Only label highest-confidence drug + headline TFs (top-15 by evidence score),
# with THRB / HNF4A / RORA in bold
df_lab <- copy(df[do_label == TRUE])
setorder(df_lab, -evidence_score)
df_lab <- df_lab[seq_len(min(15, .N))]
df_lab[, label_face := ifelse(tf_name %in% c("THRB", "HNF4A", "RORA"),
                              "bold", "plain")]

p <- ggplot(df, aes(x = activity_diff, y = logFC)) +
  # quadrant guides
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  # thick gold OUTER RING for FDA drug targets (drawn under main points)
  geom_point(data = df[is_drug_target == TRUE, ],
             aes(size = log1p(n_concord_2of3)),
             shape = 21, fill = NA, color = "#FFB300", stroke = 2.2,
             show.legend = FALSE) +
  # main hollow points — color of outline = max_pp4 (no fill)
  geom_point(aes(size = log1p(n_concord_2of3), color = max_pp4),
             shape = 1, stroke = 0.6) +
  # labels
  geom_text_repel(
    data = df_lab,
    aes(label = tf_name, fontface = label_face),
    size = 2.8, min.segment.length = 0,
    box.padding = 0.45, point.padding = 0.25,
    seed = 42, max.overlaps = Inf,
    segment.size = 0.25, segment.color = "grey45",
    family = "Helvetica"
  ) +
  # color scale: gray below 0.5 → warm magenta above (Liang convention)
  scale_color_gradientn(
    colors = c("#D9D9D9", "#BDBDBD", "#E8B6CA", "#C9265E", "#7A1140"),
    values = scales::rescale(c(0, 0.49, 0.5, 0.8, 1.0)),
    limits = c(0, 1), breaks = breaks_pp4,
    name = "Max COLOC PP4\n(23 GWAS)",
    guide = guide_colorbar(barwidth = 0.6, barheight = 4, ticks = FALSE)
  ) +
  scale_size_continuous(
    range = c(1.5, 7),
    breaks = log1p(c(0, 2, 4, 6, 8)),
    labels = c("0", "2", "4", "6", "8"),
    name = "Motif-disrupting\nvariants (2-of-3)"
  ) +
  # PI directive (2026-06-11): short single title line only. The legend keys
  # (hepatocyte regulons n = 67 TFs; bold = 4-way validated TF; gold ring =
  # FDA drug target) move to the figure caption. y-axis is C2-canonical
  # logFC (stale method token removed — data = canonical_deg_results.csv).
  labs(
    x = "SCENIC+ regulon activity Δ (disease − healthy)",
    y = "TF logFC (disease vs control)",
    title = "TF convergence"
  ) +
  theme_pub() +
  theme(
    plot.title       = element_text(size = 9, face = "bold", color = "black",
                                    family = "Helvetica"),
    plot.subtitle    = element_text(size = 7.5, color = "grey30"),
    axis.title       = element_text(size = 8, face = "bold", color = "black",
                                    family = "Helvetica"),
    axis.title.y     = element_text(angle = 90),
    axis.text        = element_text(size = 6.5, color = "black"),
    legend.position  = "right",
    legend.title     = element_text(size = 7, face = "bold"),
    legend.text      = element_text(size = 6.5),
    legend.key.size  = unit(0.35, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# Quadrant annotations — italic gray, Liang style
xlims <- range(df$activity_diff)
ylims <- range(df$logFC)
xpad  <- diff(xlims) * 0.03
ypad  <- diff(ylims) * 0.05

ann_df <- data.frame(
  x = c(xlims[2] - xpad, xlims[1] + xpad, xlims[2] - xpad, xlims[1] + xpad),
  y = c(ylims[2] - ypad, ylims[2] - ypad, ylims[1] + ypad, ylims[1] + ypad),
  hjust = c(1, 0, 1, 0),
  vjust = c(1, 1, 0, 0),
  label = c("Active, GWAS-anchored",
            "Active, no GWAS",
            "Suppressed + COLOC",
            "Suppressed, no GWAS"),
  color = c("#C9265E", "grey55", "#C9265E", "grey55")
)
p <- p + geom_text(data = ann_df,
                   aes(x = x, y = y, label = label,
                       hjust = hjust, vjust = vjust, color = NULL),
                   inherit.aes = FALSE,
                   size = 2.4, color = ann_df$color, fontface = "italic",
                   lineheight = 0.85)

out_pdf <- file.path(OUT_DIR, "fig5_tf_convergence_scatter.pdf")
ggsave(out_pdf, p, width = 6.2, height = 4.8, device = cairo_pdf)
cat(sprintf("\nWrote %s\n", out_pdf))

# Source data
src_csv <- file.path(OUT_DIR, "fig5_tf_convergence_source.csv")
fwrite(df[, .(tf_name, activity_diff, activity_padj, logFC, padj,
              n_concord_2of3, n_concord_3of3, in_disease_regulon,
              max_pp4, n_gwas_h4,
              is_headline, is_drug_target, coloc_sig, deg_sig, quadrant,
              evidence_score, do_label)],
       src_csv)
cat(sprintf("Wrote source %s\n", src_csv))

# Print top-right (suppressed regulon + DEG-up + COLOC) — biological interest
cat("\n=== Suppressed regulon + DEG-up + COLOC>=0.5 (canonical loss-of-function with genetic causality) ===\n")
print(df[activity_diff < 0 & logFC > 0 & max_pp4 >= 0.5][order(-max_pp4),
        .(tf_name, activity_diff, logFC, n_concord_2of3, max_pp4)])

cat("\n=== Multi-line evidence (>=3 layers: regulon padj<0.05, DEG padj<0.05, motif>=2, COLOC>=0.5) ===\n")
n_layers <- with(df, (activity_padj < 0.05) + (padj < 0.05) +
                     (n_concord_2of3 >= 2) + (max_pp4 >= 0.5))
df[, n_evidence_layers := n_layers]
print(df[n_evidence_layers >= 3][order(-n_evidence_layers, -max_pp4),
        .(tf_name, activity_diff, logFC, padj,
          n_concord_2of3, max_pp4, n_evidence_layers)])

cat("\nDone.\n")
