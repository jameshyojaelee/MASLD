# =============================================================================
# Fig 5 (Convergence) — 4-way TF survival lollipop
#
# Per-TF horizontal "stacked" lollipop showing how many evidence layers each
# TF carries:
#   - motif disruption    (segment 0 -> 1)
#   - eQTL                (segment 1 -> 2)
#   - DA peak             (segment 2 -> 3)
#   - Currin caQTL        (segment 3 -> 4)
#
# Lollipop terminator (circle): size encodes per-TF caQTL motif-agreement %.
# Headline 4-way TFs (THRB, HNF4A, RORA, MLXIPL, MAX) shown in bold.
# Drug-target TFs (THRB, NR1H4, PPARA, PPARG) get a gold ring on the terminator.
#
# Output: figures/main/fig5_convergence/panels/fig5_tf_4way_survival_lollipop.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE_DIR, "scripts/figures/load_figure_data.R"))
source(file.path(BASE_DIR, "scripts/figures/publication_theme.R"))

OUT_DIR <- file.path(FIG5_DIR, "panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_ATAC <- file.path(BASE_DIR, "GWAS/finemapping/results/gwas_atac")

# -----------------------------------------------------------------------------
# 1. Per-TF evidence counts
# -----------------------------------------------------------------------------
allele <- fread(file.path(GWAS_ATAC, "allele_concordance_summary.csv"))
# Split compound TF names "X::Y" -> single TFs
allele_split <- allele[, .(tf_single = unlist(strsplit(tf_name, "::"))),
                        by = .(tf_name, n_total, n_with_eqtl, n_with_da,
                               n_concord_3of3, n_concord_2of3,
                               is_activator, in_disease_regulon)]
allele_tf <- allele_split[, .(
  n_total      = sum(n_total),
  n_with_eqtl  = sum(n_with_eqtl),
  n_with_da    = sum(n_with_da),
  n_concord_2of3 = sum(n_concord_2of3),
  n_concord_3of3 = sum(n_concord_3of3),
  in_disease_regulon = any(in_disease_regulon)
), by = .(tf_name = tf_single)]

# Motif disruption count per TF (number of variant rows where TF motif is hit)
motif <- fread(file.path(GWAS_ATAC, "motif_disruption_scores.csv"),
               select = c("SNP_id", "tf_name", "alleleDiff", "max_pip"))
motif_split <- motif[, .(tf_single = unlist(strsplit(tf_name, "::"))),
                      by = .(tf_name, SNP_id)]
motif_per_tf <- motif_split[, .(n_motif_hits = uniqueN(SNP_id)),
                             by = .(tf_name = tf_single)]

# Currin caQTL agreement
caqtl <- fread(file.path(GWAS_ATAC, "currin_caqtl_per_tf.csv"))
caqtl_per_tf <- caqtl[, .(tf_name,
                          n_variants_caqtl    = n_variants,
                          n_caqtl_tested      = n_caqtl_tested,
                          n_caqtl_p1e3        = n_caqtl_p1e3,
                          n_motif_caqtl_agree = n_motif_caqtl_agree,
                          pct_motif_caqtl_agree = pct_motif_caqtl_agree,
                          any_disease_regulon = any_disease_regulon)]

# -----------------------------------------------------------------------------
# 2. Merge into a single TF table, derive evidence flags
# -----------------------------------------------------------------------------
df <- merge(motif_per_tf, allele_tf, by = "tf_name", all = TRUE)
df <- merge(df, caqtl_per_tf, by = "tf_name", all = TRUE)
for (col in c("n_motif_hits", "n_total", "n_with_eqtl", "n_with_da",
              "n_concord_2of3", "n_concord_3of3", "n_caqtl_tested",
              "n_motif_caqtl_agree", "pct_motif_caqtl_agree")) {
  if (col %in% names(df)) df[is.na(get(col)), (col) := 0]
}
df[is.na(in_disease_regulon), in_disease_regulon := FALSE]

# Evidence layer presence
df[, has_motif := n_motif_hits >= 1]
df[, has_eqtl  := n_with_eqtl  >= 1]
df[, has_da    := n_with_da    >= 1]
df[, has_caqtl := n_caqtl_tested >= 1 & n_motif_caqtl_agree >= 1]

df[, n_layers := as.integer(has_motif) + as.integer(has_eqtl) +
                 as.integer(has_da)    + as.integer(has_caqtl)]

# Combined evidence score for ranking (weight 4-way layers + motif hits)
df[, evidence_score := n_layers * 10 +
                       log1p(n_motif_hits) +
                       log1p(n_concord_2of3) +
                       0.05 * pct_motif_caqtl_agree]

# -----------------------------------------------------------------------------
# 3. Pick top 20 by combined evidence
# -----------------------------------------------------------------------------
HEADLINE_TFS <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")
DRUG_TARGETS <- c("THRB", "NR1H4", "PPARA", "PPARG")

# Headline + drug-target TFs are always included; pad to 20 by evidence_score
must_include <- unique(c(HEADLINE_TFS, DRUG_TARGETS))
must_include <- must_include[must_include %in% df$tf_name]

remaining <- df[!(tf_name %in% must_include)][order(-evidence_score)]
top_extra <- remaining[1:max(0, 20 - length(must_include)), tf_name]

top20 <- df[tf_name %in% c(must_include, top_extra)]
top20 <- top20[order(-n_layers, -evidence_score)]
top20[, tf_name := factor(tf_name, levels = rev(tf_name))]   # top of plot = best

cat("Top 20 TFs by 4-way evidence:\n")
print(top20[, .(tf_name, n_layers, has_motif, has_eqtl, has_da, has_caqtl,
                n_motif_hits, n_concord_2of3, pct_motif_caqtl_agree)])

# -----------------------------------------------------------------------------
# 4. Build long segment table — one row per (TF, layer) where layer is present
# -----------------------------------------------------------------------------
LAYER_ORDER <- c("motif", "eqtl", "da", "caqtl")
LAYER_COLOR <- c(motif = "#C9265E",  # warm magenta
                 eqtl  = "#9C27B0",  # violet
                 da    = "#00695C",  # teal
                 caqtl = "#FFB300")  # gold

build_segments <- function(dt) {
  out <- list()
  for (i in seq_len(nrow(dt))) {
    row <- dt[i]
    layers_present <- c(
      motif = row$has_motif,
      eqtl  = row$has_eqtl,
      da    = row$has_da,
      caqtl = row$has_caqtl
    )
    # Cumulative survival: x extends as long as each consecutive layer is present
    x_left <- 0
    for (k in seq_along(LAYER_ORDER)) {
      lyr <- LAYER_ORDER[k]
      if (isTRUE(layers_present[lyr])) {
        out[[length(out) + 1]] <- data.frame(
          tf_name = row$tf_name,
          layer   = lyr,
          x_from  = x_left,
          x_to    = x_left + 1L,
          y       = as.numeric(row$tf_name)
        )
        x_left <- x_left + 1L
      } else {
        break   # cascade: stop at first missing layer
      }
    }
  }
  rbindlist(out)
}

segs <- build_segments(top20)
segs[, layer_color := LAYER_COLOR[layer]]
segs[, layer := factor(layer, levels = LAYER_ORDER)]

# -----------------------------------------------------------------------------
# 5. Lollipop terminator: at x = n_layers, size = pct_motif_caqtl_agree
# -----------------------------------------------------------------------------
top20[, x_end := n_layers]
top20[, y_pos := as.numeric(tf_name)]
top20[, is_headline := tf_name %in% HEADLINE_TFS]
top20[, is_drug     := tf_name %in% DRUG_TARGETS]

# Terminator fill: deeper magenta if 4-way validated, lighter if partial
top20[, term_fill := fifelse(n_layers >= 4, "#8E1B43",
                       fifelse(n_layers >= 3, "#C9265E",
                         fifelse(n_layers >= 2, "#E091AE", "#BDBDBD")))]

# Size range: pct 0..100 -> point radius
size_range <- c(1.6, 5.0)

# -----------------------------------------------------------------------------
# 6. Plot
# -----------------------------------------------------------------------------
n_tf <- nrow(top20)
y_breaks <- seq_len(n_tf)
y_labels <- as.character(levels(top20$tf_name))

# Bold font for headline TFs (axis text)
y_face <- ifelse(y_labels %in% HEADLINE_TFS, "bold.italic", "italic")
y_color <- ifelse(y_labels %in% HEADLINE_TFS, "#8E1B43", "black")

# Vertical gridlines at 1..4 (evidence count axis)
xgrid_df <- data.frame(x = 0:4)

p <- ggplot() +
  # vertical guides
  geom_vline(data = xgrid_df, aes(xintercept = x),
             color = "grey92", linewidth = 0.25) +
  # horizontal stem
  geom_segment(data = top20,
               aes(x = 0, xend = x_end, y = y_pos, yend = y_pos),
               color = "grey85", linewidth = 0.25) +
  # evidence segments (colored stripes per layer)
  geom_segment(data = segs,
               aes(x = x_from, xend = x_to, y = y, yend = y,
                   color = layer),
               linewidth = 2.5, lineend = "butt") +
  scale_color_manual(
    values = LAYER_COLOR,
    breaks = LAYER_ORDER,
    labels = c(motif = "Motif disruption",
               eqtl  = "eQTL",
               da    = "DA peak",
               caqtl = "Currin caQTL"),
    name = "Evidence layer (cumulative)"
  ) +
  # Gold ring under FDA drug-target terminators
  geom_point(data = top20[is_drug == TRUE],
             aes(x = x_end, y = y_pos),
             shape = 21, fill = NA, color = "#FFB300",
             stroke = 1.1, size = 6.0) +
  # Lollipop terminators (size = caQTL motif-agreement %)
  geom_point(data = top20,
             aes(x = x_end, y = y_pos, size = pct_motif_caqtl_agree),
             shape = 21, fill = top20$term_fill,
             color = "black", stroke = 0.4) +
  scale_size_continuous(
    range = size_range, breaks = c(0, 50, 75, 100),
    name = "Motif↔caQTL\nagreement (%)",
    limits = c(0, 100)
  ) +
  # Per-TF agreement % label, on the right side
  geom_text(data = top20[pct_motif_caqtl_agree > 0],
            aes(x = x_end + 0.18, y = y_pos,
                label = sprintf("%.0f%%", pct_motif_caqtl_agree)),
            size = 1.9, hjust = 0, color = "grey25") +
  scale_x_continuous(
    breaks = 0:4, limits = c(-0.05, 4.9),
    expand = c(0, 0),
    name = "Cumulative evidence layers"
  ) +
  scale_y_continuous(
    breaks = y_breaks, labels = y_labels, expand = c(0.02, 0.02),
    name = NULL
  ) +
  # PI directive (2026-06-11): short single title line only. The cumulative-
  # layer order (Motif -> eQTL -> DA peak -> Currin caQTL) is encoded by the
  # color legend; the "top 20 TFs" detail moves to the figure caption.
  labs(
    title = "4-way TF survival"
  ) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(
    axis.text.y     = element_text(face = y_face, color = y_color,
                                   size = PUB_AXIS_TEXT + 0.2),
    axis.title.x    = element_text(face = "bold"),
    plot.title      = element_text(size = PUB_TITLE, face = "bold"),
    plot.subtitle   = element_text(size = PUB_SUBTITLE, color = "grey30"),
    legend.position = "right",
    legend.box      = "vertical",
    legend.spacing.y = unit(0.1, "cm"),
    panel.grid      = element_blank(),
    plot.margin     = margin(8, 8, 8, 8)
  ) +
  guides(
    color = guide_legend(order = 1,
                         override.aes = list(linewidth = 4)),
    size  = guide_legend(order = 2,
                         override.aes = list(shape = 21,
                                             fill = "#C9265E",
                                             color = "black"))
  )

out_pdf <- file.path(OUT_DIR, "fig5_tf_4way_survival_lollipop.pdf")
ggsave(out_pdf, p, width = 6.0, height = 7.0, device = cairo_pdf)
cat(sprintf("\nWrote %s\n", out_pdf))

# Source data
src_csv <- file.path(OUT_DIR, "fig5_tf_4way_survival_lollipop_source.csv")
fwrite(top20[, .(tf_name, n_layers, has_motif, has_eqtl, has_da, has_caqtl,
                 n_motif_hits, n_with_eqtl, n_with_da, n_concord_2of3,
                 n_caqtl_tested, n_motif_caqtl_agree, pct_motif_caqtl_agree,
                 in_disease_regulon, is_headline, is_drug, evidence_score)],
       src_csv)
cat(sprintf("Wrote source %s\n", src_csv))
cat("\nDone.\n")
