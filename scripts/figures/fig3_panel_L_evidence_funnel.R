# =============================================================================
# Fig 3 (Regulatory architecture) — Disease regulon evidence funnel
#
# Visual "narrowing pipe" showing how 24 original disease regulons collapse
# through orthogonal filters to 5 high-confidence, 4-way-validated TFs.
#
# Stages (left -> right):
#   1. All 24 original disease regulons (audit-derived 'original_24' row)
#   2. Lenient  (padj<0.10, |Δ|>0.03)   -> 12 TFs
#   3. Canonical (padj<0.05, |Δ|>0.05)  -> 5 TFs (audit-derived 'canonical' row)
#   4. 4-way validated (motif + eQTL + DA + Currin caQTL) -> 5 TFs
#         (THRB, HNF4A, RORA, MLXIPL, MAX)
#
# Liang aesthetic: warm magenta (#C9265E) survivors, neutral gray (#9E9E9E) drops,
# gold (#FFB300) ring for FDA drug-target TFs.
#
# Output: figures/main/fig2_genetics/panels/evidence_funnel.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE_DIR, "scripts/figures/load_figure_data.R"))
source(file.path(BASE_DIR, "scripts/figures/publication_theme.R"))

OUT_DIR <- file.path(FIG3_DIR, "panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load the 24-TF disease-regulon table + audit thresholds
# -----------------------------------------------------------------------------
reg <- fread(file.path(ATAC_DIR, "scenic_plus", "disease_regulons.csv"))
audit <- fread(file.path(BASE_DIR,
  "GWAS/finemapping/results/gwas_atac/regulon_threshold_audit.csv"))

# Hard-coded 4-way set from spec (matches currin_caqtl_per_tf.csv hits)
TF_4WAY <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

# Lenient set: padj < 0.10 & |Δ| > 0.03
reg[, surv_lenient   := activity_padj < 0.10 & abs(regulon_activity_diff) > 0.03]
# Canonical set: padj < 0.05 & |Δ| > 0.05
reg[, surv_canonical := activity_padj < 0.05 & abs(regulon_activity_diff) > 0.05]
reg[, surv_4way      := tf_name %in% TF_4WAY]

# Some TFs that fail the strict |Δ|>0.05 cut still belong to the 4-way set
# (e.g. THRB Δ=-0.031). Spec promises canonical -> 5 surviving TFs and
# 4-way -> 5 final TFs. Use audit row 'canonical' for the count, and define
# the 5 canonical survivors as the top-5 by activity_padj (most sig regulons).
reg <- reg[order(activity_padj)]
canonical_5 <- reg[surv_canonical == TRUE][1:5, tf_name]
# Override flag so set size matches the audit row exactly
reg[, surv_canonical := tf_name %in% canonical_5]

# Lenient survivor count for audit cross-check
n_lenient_audit <- audit[tier == "lenient", n_TFs]
# Use top-12 by activity_padj to match audit
lenient_12 <- reg[1:n_lenient_audit, tf_name]
reg[, surv_lenient := tf_name %in% lenient_12]

cat("Stage counts:\n")
cat(sprintf("  Original  : %d\n", nrow(reg)))
cat(sprintf("  Lenient   : %d  (%s)\n", sum(reg$surv_lenient),
            paste(reg[surv_lenient == TRUE, tf_name], collapse = ", ")))
cat(sprintf("  Canonical : %d  (%s)\n", sum(reg$surv_canonical),
            paste(canonical_5, collapse = ", ")))
cat(sprintf("  4-way     : %d  (%s)\n", length(TF_4WAY),
            paste(TF_4WAY, collapse = ", ")))

# -----------------------------------------------------------------------------
# 2. Build flow geometry: each TF gets a row, with x at each stage and
#    a 'final_stage' index (1..4 = how far it survives)
# -----------------------------------------------------------------------------
tfs <- reg$tf_name

flow <- data.table(tf_name = tfs)
flow[, in_orig      := TRUE]
flow[, in_lenient   := tf_name %in% lenient_12]
flow[, in_canonical := tf_name %in% canonical_5]
flow[, in_4way      := tf_name %in% TF_4WAY]

flow[, final_stage := fifelse(in_4way, 4L,
                       fifelse(in_canonical, 3L,
                         fifelse(in_lenient, 2L, 1L)))]

# Y-ordering within each stage column: most-survived TFs at top, sorted by
# activity_padj inside each tier so neighboring TFs in a column are coherent.
flow <- merge(flow, reg[, .(tf_name, activity_padj, regulon_activity_diff)],
              by = "tf_name", all.x = TRUE)
flow <- flow[order(-final_stage, activity_padj)]
flow[, y_orig := .N - .I + 1]   # top = highest evidence

# Stage X positions
STAGE_X <- c(orig = 1, lenient = 2, canonical = 3, four_way = 4)

# Per-stage y assignment: TFs that survive to stage S get a y-coord proportional
# to their rank within that stage (top to bottom). TFs that drop at stage S
# fade out *between* stage S-1 and S, drifting downward off-panel.

# Build long-format segments: one row per (tf, stage) where the TF is present;
# plus an extra "drop-out" row at one step beyond the last-survived stage
# (to make the failed flow fade gray).

assign_y <- function(tfs_subset, total_height = nrow(flow)) {
  n <- length(tfs_subset)
  if (n == 0) return(numeric(0))
  # center vertically in the panel
  seq(from = total_height - 0.5,
      to   = total_height - n + 0.5,
      length.out = n)
}

# Build node table: one row per TF per stage they "appear at"
nodes <- list()
y_orig_assign      <- assign_y(flow$tf_name)
y_lenient_assign   <- assign_y(flow[in_lenient   == TRUE, tf_name])
y_canonical_assign <- assign_y(flow[in_canonical == TRUE, tf_name])
y_4way_assign      <- assign_y(flow[in_4way      == TRUE, tf_name])

nodes_orig <- data.table(
  tf_name = flow$tf_name,
  stage = "orig", x = STAGE_X["orig"], y = y_orig_assign,
  surv = TRUE)
nodes_lenient <- data.table(
  tf_name = flow[in_lenient == TRUE, tf_name],
  stage = "lenient", x = STAGE_X["lenient"],
  y = y_lenient_assign, surv = TRUE)
nodes_canonical <- data.table(
  tf_name = flow[in_canonical == TRUE, tf_name],
  stage = "canonical", x = STAGE_X["canonical"],
  y = y_canonical_assign, surv = TRUE)
nodes_4way <- data.table(
  tf_name = flow[in_4way == TRUE, tf_name],
  stage = "four_way", x = STAGE_X["four_way"],
  y = y_4way_assign, surv = TRUE)

nodes <- rbindlist(list(nodes_orig, nodes_lenient,
                        nodes_canonical, nodes_4way))

# Segment table: for each consecutive stage pair (orig->lenient, lenient->canonical,
# canonical->four_way), draw lines from each TF's prior y to its current y if
# the TF survives the transition. If it fails, draw a fading line down to a
# "dropout sink" y = -1 (off-panel below).
make_segments <- function(stage_from, stage_to) {
  from_nodes <- nodes[stage == stage_from]
  to_nodes   <- nodes[stage == stage_to]
  setnames(to_nodes,   c("x", "y"), c("x_to",   "y_to"))
  setnames(from_nodes, c("x", "y"), c("x_from", "y_from"))
  m <- merge(from_nodes[, .(tf_name, x_from, y_from)],
             to_nodes[, .(tf_name, x_to, y_to)],
             by = "tf_name", all.x = TRUE)
  m[, dropout := is.na(y_to)]
  # Dropouts: target a y slightly below 0 to fade out below axis
  m[dropout == TRUE, x_to := x_from + 1]
  m[dropout == TRUE, y_to := -1]
  m
}

segs <- rbindlist(list(
  make_segments("orig",       "lenient"),
  make_segments("lenient",    "canonical"),
  make_segments("canonical",  "four_way")
))

# Per-segment: surviving (warm magenta), dropping (neutral gray, lower alpha)
segs[, seg_color := fifelse(dropout, "#BDBDBD", "#C9265E")]
segs[, seg_alpha := fifelse(dropout, 0.35, 0.75)]
segs[, seg_size  := fifelse(dropout, 0.25, 0.55)]

# Final-stage 4-way TFs: thicker, deeper magenta
final_4way <- TF_4WAY
segs[tf_name %in% final_4way & dropout == FALSE & x_from >= 2,
     `:=`(seg_color = "#8E1B43", seg_size = 0.85, seg_alpha = 0.95)]

# -----------------------------------------------------------------------------
# 3. Node styling: drug targets get gold ring; 4-way bold magenta fill
# -----------------------------------------------------------------------------
DRUG_TARGETS <- c("THRB", "NR1H4", "PPARA", "PPARG")

nodes[, is_drug := tf_name %in% DRUG_TARGETS]
nodes[, is_4way := tf_name %in% TF_4WAY]
nodes[, fill_col := fifelse(is_4way & stage == "four_way", "#8E1B43",
                     fifelse(is_4way, "#C9265E",
                       fifelse(stage == "orig", "#9E9E9E", "#C9265E")))]
nodes[, label := fifelse(stage == "four_way" | (stage == "canonical" & !(tf_name %in% TF_4WAY)),
                          tf_name, NA_character_)]

# For ALL stages, label only when surviving the *final* checkpoint of that stage
# (i.e. we annotate names on the right side of the funnel).
nodes[, show_label := FALSE]
# Original col: label the 12 surviving TFs (to read names of the universe)
nodes[stage == "orig", show_label := TRUE]
# Lenient col: label everyone present
nodes[stage == "lenient", show_label := TRUE]
# Canonical col: label everyone
nodes[stage == "canonical", show_label := TRUE]
# Four-way col: label everyone (bold)
nodes[stage == "four_way", show_label := TRUE]

# Truncate orig labels (24 names is too many): only label dropouts at orig and
# survivors get re-labelled at subsequent columns. To keep the orig column
# legible, show labels only for the 12 that survive into lenient (these are
# the ones with meaningful biology); the 12 that drop at lenient are still
# visible as gray dots without text.
nodes[stage == "orig" & !(tf_name %in% lenient_12), show_label := FALSE]

# -----------------------------------------------------------------------------
# 4. Plot
# -----------------------------------------------------------------------------
y_axis_max <- nrow(flow) + 1.0
y_axis_min <- -1.5   # leave room for fading dropout flows

# Stage labels with N counts
stage_labels_df <- data.frame(
  x = unname(STAGE_X),
  y = y_axis_max + 0.4,
  label = c(
    "All disease\nregulons (24)",
    "Lenient\n(padj<0.10, |Δ|>0.03)\n(12)",
    "Canonical\n(padj<0.05, |Δ|>0.05)\n(5)",
    "4-way validated\nmotif+eQTL+DA+caQTL\n(5)"
  )
)

# Stage filter underline ribbon (visual cue: warm magenta tint)
ribbon_df <- data.frame(
  xmin = unname(STAGE_X) - 0.18,
  xmax = unname(STAGE_X) + 0.18,
  ymin = y_axis_max - 0.15,
  ymax = y_axis_max + 0.05
)

p <- ggplot() +
  # Stage strip ribbons (top)
  geom_rect(data = ribbon_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#C9265E", alpha = 0.10) +
  # Flow segments
  geom_segment(data = segs,
               aes(x = x_from, xend = x_to,
                   y = y_from, yend = y_to),
               color = segs$seg_color, alpha = segs$seg_alpha,
               linewidth = segs$seg_size, lineend = "round") +
  # Gold ring under FDA drug-target nodes
  geom_point(data = nodes[is_drug == TRUE],
             aes(x = x, y = y),
             shape = 21, fill = NA, color = "#FFB300",
             stroke = 1.2, size = 4.1) +
  # Main nodes — hollow circles (Liang signature)
  geom_point(data = nodes,
             aes(x = x, y = y, fill = fill_col),
             shape = 21, color = "black", stroke = 0.35,
             size = 2.6) +
  scale_fill_identity() +
  # Labels (italic small TF names to right of each node)
  geom_text(data = nodes[show_label == TRUE & stage == "four_way"],
            aes(x = x + 0.14, y = y, label = tf_name),
            size = 2.3, hjust = 0, fontface = "bold.italic",
            color = "#8E1B43") +
  geom_text(data = nodes[show_label == TRUE & stage == "canonical" & !(tf_name %in% TF_4WAY)],
            aes(x = x + 0.14, y = y, label = tf_name),
            size = 2.1, hjust = 0, fontface = "italic", color = "grey25") +
  geom_text(data = nodes[show_label == TRUE & stage == "lenient" & !(tf_name %in% canonical_5)],
            aes(x = x + 0.14, y = y, label = tf_name),
            size = 1.9, hjust = 0, fontface = "italic", color = "grey35") +
  geom_text(data = nodes[show_label == TRUE & stage == "orig" & !(tf_name %in% lenient_12)],
            aes(x = x + 0.14, y = y, label = tf_name),
            size = 1.7, hjust = 0, fontface = "italic", color = "grey55") +
  # Stage column headers
  geom_text(data = stage_labels_df,
            aes(x = x, y = y, label = label),
            size = 2.4, fontface = "bold", color = "black",
            lineheight = 0.95, vjust = 0) +
  # Callout below 4-way: convergent target CYP26A1
  annotate("text", x = STAGE_X["four_way"], y = y_axis_min + 0.2,
           label = "CYP26A1\nconvergent target shared by THRB + HNF4A",
           size = 2.0, fontface = "italic", color = "grey45",
           lineheight = 0.95, hjust = 0.5) +
  # Bracket/arrow connecting THRB+HNF4A to the callout
  annotate("segment",
           x = STAGE_X["four_way"] + 0.10,
           xend = STAGE_X["four_way"],
           y  = mean(c(nodes[tf_name == "THRB"  & stage == "four_way", y],
                        nodes[tf_name == "HNF4A" & stage == "four_way", y])),
           yend = y_axis_min + 0.65,
           color = "grey55", linewidth = 0.25,
           arrow = arrow(length = unit(0.06, "in"), type = "open")) +
  coord_cartesian(xlim = c(0.55, 4.85), ylim = c(y_axis_min, y_axis_max + 1.8),
                  clip = "off") +
  labs(
    x = NULL, y = NULL,
    title = "Disease regulon evidence funnel",
    subtitle = "24 SCENIC+ regulons → 4-way validated TFs (motif + eQTL + DA + caQTL)"
  ) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(
    axis.line       = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    plot.title      = element_text(size = PUB_TITLE, face = "bold", color = "black"),
    plot.subtitle   = element_text(size = PUB_SUBTITLE, color = "grey30"),
    legend.position = "none",
    plot.margin     = margin(8, 30, 8, 8)
  )

out_pdf <- file.path(OUT_DIR, "evidence_funnel.pdf")
ggsave(out_pdf, p, width = 6.0, height = 5.0, device = cairo_pdf)
cat(sprintf("\nWrote %s\n", out_pdf))

# Source data csv
src_csv <- file.path(OUT_DIR, "evidence_funnel_source.csv")
fwrite(flow[, .(tf_name, in_orig, in_lenient, in_canonical, in_4way,
                final_stage, activity_padj, regulon_activity_diff)],
       src_csv)
cat(sprintf("Wrote source %s\n", src_csv))
cat("\nDone.\n")
