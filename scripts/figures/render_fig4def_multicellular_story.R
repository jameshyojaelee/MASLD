#!/usr/bin/env Rscript
# KEY MESSAGE: Two frozen MASLD programs remain disease-associated after explicit
# ambient qualification, localize to different non-parenchymal lineages in
# donor-paired contrasts, and transport into cross-sectional bulk fibrosis state.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
OUT <- Sys.getenv(
  "FIG4DEF_CANDIDATE_ROOT",
  file.path(BASE, "figures/candidates/fig4def-multicellular-story-2026-08-17-v2")
)
if (dir.exists(OUT)) {
  stop("Candidate output already exists: ", OUT, call. = FALSE)
}

source(file.path(BASE, "scripts/figures/publication_theme.R"))
dir.create(file.path(OUT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "source_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "captions"), recursive = TRUE, showWarnings = FALSE)

V8 <- file.path(
  BASE,
  "figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8"
)
inputs <- data.table(
  role = c(
    "style_anchor", "style_anchor", "input", "input", "input", "theme"
  ),
  path = c(
    file.path(BASE, "figures/main/fig4_singlecell_programs/panels/fig4b_scrna_umap_embeddable.pdf"),
    file.path(BASE, "figures/main/fig4_singlecell_programs/panels/fig4c_hotspot_stage_heatmap.pdf"),
    file.path(V8, "source_tables/fig4d_ambient_effect_transport.tsv"),
    file.path(V8, "source_tables/fig4e_same_atlas_lineage_specificity.tsv"),
    file.path(V8, "source_tables/fig4f_bulk_stage_transport.tsv"),
    file.path(BASE, "scripts/figures/publication_theme.R")
  ),
  expected_sha256 = c(
    "8b22e80ae170a8668e022d779b52caacbbad2329500a38423d1ef5502b2cb7df",
    "5d3e91b3cb23f37189ca6e642c86768c0ad964a0a0c4b3796781b4a9bc8a0526",
    "ad338a3adc155e658f876c3b086f9051543ff7ac5c87e143cc06223e3e4d2730",
    "4e03b9b9fb4f8ccbade689b4a10d1afe307dc59ac092751faeacb5afedcd3879",
    "1a1065bd37739f133bfd6286733d150f9753cec691a2b1c09a739c81a64bd936",
    "0305135a2d68e37b1fde19008ea659e5aa39b533373f5e0d701f702b915d2874"
  )
)

sha256 <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(value)) stop("sha256sum failed for ", path, call. = FALSE)
  strsplit(value[[1]], "[[:space:]]+")[[1]][[1]]
}
if (any(!file.exists(inputs$path))) {
  stop("One or more candidate inputs are missing", call. = FALSE)
}
inputs[, observed_sha256 := vapply(path, sha256, character(1))]
inputs[, bytes := file.info(path)$size]
if (any(inputs$observed_sha256 != inputs$expected_sha256)) {
  print(inputs[observed_sha256 != expected_sha256])
  stop("Input drift detected", call. = FALSE)
}
fwrite(inputs, file.path(OUT, "input_manifest.tsv"), sep = "\t")

d <- fread(inputs[role == "input"]$path[[1]], na.strings = c("", "NA"))
e <- fread(inputs[role == "input"]$path[[2]], na.strings = c("", "NA"))
f <- fread(inputs[role == "input"]$path[[3]], na.strings = c("", "NA"))

hero_uid <- c(
  ECM = "hotspot_hepatocytes_f05c535ae5bbc0b9",
  Ductular = "hotspot_hepatocytes_48f39dd4d817a10e"
)
program_colors <- c(
  "ECM/IGFBP7" = "#FF6F00",
  "Ductular-injury/BICC1" = "#00BCD4"
)
program_shapes <- c(
  "ECM/IGFBP7" = 21L,
  "Ductular-injury/BICC1" = 22L
)

numeric_d <- c(
  "raw_beta", "corrected_beta", "delta_beta", "raw_qvalue",
  "corrected_qvalue", "delta_qvalue", "raw_hc3_qvalue",
  "corrected_hc3_qvalue", "delta_hc3_qvalue"
)
d[, (numeric_d) := lapply(.SD, as.numeric), .SDcols = numeric_d]
d[, plotted := is.finite(raw_beta) & is.finite(corrected_beta)]
d[, biological_unit := "biological donor"]
d[, multiple_testing_family := paste0(
  "117 frozen programs, BH separately for raw, corrected, and paired-delta models"
)]
d[, state_reason := fifelse(
  evidence_state == "untestable",
  "complete-atlas integer transport does not reproduce the native non-integer T-cell scoring substrate",
  fifelse(
    evidence_state == "attenuation_indeterminate",
    "paired correction-induced change does not pass the 117-program BH family under ordinary or HC3 inference",
    "raw and ambient-corrected donor-level effects are estimable"
  )
)]
d[, unresolved_alternative := paste0(
  "correction-sensitive signal may reflect ambient RNA, score transport, or shared multicellular disease state; cellular origin remains unresolved"
)]
d[, claim_boundary := paste0(
  "ambient sensitivity of donor-level cross-sectional stage association; not lineage origin or cell autonomy"
)]
d[, panel_role := fifelse(
  program_uid == hero_uid[["ECM"]], "hero_ECM_IGFBP7",
  fifelse(program_uid == hero_uid[["Ductular"]], "hero_ductular_BICC1",
          fifelse(plotted, "context_program", "untestable_program"))
)]
fwrite(d, file.path(OUT, "source_tables/fig4d_ambient_effect_transport.tsv"), sep = "\t")

e[, biological_unit := "paired biological donor"]
e[, multiple_testing_family := paste0(
  "10 prespecified program-by-lineage contrasts, BH across the complete family"
)]
e[, evidence_state := fifelse(
  grepl("reference", comparison, fixed = TRUE), "reference",
  fifelse(is.finite(hc3_qvalue) & hc3_qvalue < 0.05, "hc3_family_supported",
          "tested_nonsignificant")
)]
e[, state_reason := fifelse(
  evidence_state == "reference", "hepatocyte reference score fixed at zero",
  fifelse(
    evidence_state == "hc3_family_supported",
    "lineage-minus-hepatocyte contrast passes the ten-comparison BH family under ordinary and HC3 inference",
    "lineage-minus-hepatocyte contrast does not pass the ten-comparison HC3 BH family"
  )
)]
e[, unresolved_alternative := paste0(
  "ambient carryover, within-lineage state composition, and a shared multicellular response remain possible"
)]
e[, claim_boundary := paste0(
  "same-atlas descriptive localization; not transcript origin, secretion, or cell autonomy"
)]
fwrite(e, file.path(OUT, "source_tables/fig4e_same_atlas_lineage_specificity.tsv"), sep = "\t")

f[, biological_unit := "bulk RNA-seq participant"]
f[, multiple_testing_family := paste0(
  "23,370 genes within each fibrosis-stage-versus-F0 contrast; no program-level significance test"
)]
f[, evidence_state := "descriptive_tissue_transport"]
f[, state_reason := paste0(
  "fixed L1-weighted member-gene effect with ", n_bh, "/", n_genes,
  " available members passing the gene-level BH family"
)]
f[, unresolved_alternative := paste0(
  "bulk tissue composition and shared fibrosis state cannot establish cellular origin"
)]
f[, claim_boundary := paste0(
  "cross-sectional fixed-member tissue-state transport; not longitudinal progression or a program-level inferential test"
)]
fwrite(f, file.path(OUT, "source_tables/fig4f_bulk_stage_transport.tsv"), sep = "\t")

validation <- list()
check <- function(name, condition, observed, expected) {
  validation[[length(validation) + 1L]] <<- data.table(
    check = name,
    status = if (isTRUE(condition)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected)
  )
  if (!isTRUE(condition)) {
    stop(name, ": observed ", observed, "; expected ", expected, call. = FALSE)
  }
}
check("fixed_program_rows", nrow(d) == 117L, nrow(d), 117L)
check("ambient_evaluable", sum(d$plotted) == 106L, sum(d$plotted), 106L)
check("ambient_untestable", sum(!d$plotted) == 11L, sum(!d$plotted), 11L)
check("hero_programs", sum(d$program_uid %in% hero_uid) == 2L,
      sum(d$program_uid %in% hero_uid), 2L)
check("hero_delta_hc3_not_supported",
      all(d[program_uid %in% hero_uid]$delta_hc3_qvalue >= 0.05),
      paste(signif(d[program_uid %in% hero_uid]$delta_hc3_qvalue, 3), collapse = ";"),
      "both >= 0.05")
check("lineage_planned_family", sum(!grepl("reference", e$comparison, fixed = TRUE)) == 10L,
      sum(!grepl("reference", e$comparison, fixed = TRUE)), 10L)
check("lineage_references", sum(grepl("reference", e$comparison, fixed = TRUE)) == 2L,
      sum(grepl("reference", e$comparison, fixed = TRUE)), 2L)
check("bulk_transport_rows", nrow(f) == 8L, nrow(f), 8L)
check("bulk_transport_stages", uniqueN(f$stage) == 4L, uniqueN(f$stage), 4L)

theme_panel <- function() {
  theme_masld(base_size = 6) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 6, face = "plain", margin = margin(b = 2)),
      legend.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      plot.margin = margin(3, 4, 3, 3)
    )
}

# 4D: identity-line transport shows whether ambient correction changes each
# frozen program's donor-level disease association.
d_plot <- d[plotted == TRUE]
d_hero <- d_plot[program_uid %in% hero_uid]
d_hero[, program_name := fifelse(
  program_uid == hero_uid[["ECM"]], "ECM/IGFBP7", "Ductular-injury/BICC1"
)]
d_hero[, `:=`(
  color = unname(program_colors[program_name]),
  shape = unname(program_shapes[program_name]),
  parsed_label = fifelse(program_name == "ECM/IGFBP7", "italic(IGFBP7)", "italic(BICC1)"),
  label_x = raw_beta + 0.055,
  label_y = corrected_beta + fifelse(program_name == "ECM/IGFBP7", -0.055, 0.055)
)]
xy <- range(c(d_plot$raw_beta, d_plot$corrected_beta), finite = TRUE)
xy <- c(floor(min(xy[[1]], -0.72) * 10) / 10, ceiling(max(xy[[2]], 0.82) * 10) / 10)
p_d <- ggplot(d_plot, aes(x = raw_beta, y = corrected_beta)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_point(shape = 21, size = 1.55, stroke = 0.25, fill = "white", color = "#BDBDBD") +
  geom_segment(
    data = d_hero,
    aes(x = raw_beta, xend = raw_beta, y = raw_beta, yend = corrected_beta),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = d_hero$color,
    arrow = arrow(length = unit(0.045, "inches"), type = "closed")
  ) +
  geom_point(
    data = d_hero,
    aes(shape = program_name, fill = program_name),
    size = 2.25, stroke = 0.3, color = "black"
  ) +
  geom_text(
    data = d_hero,
    aes(x = label_x, y = label_y, label = parsed_label),
    inherit.aes = FALSE,
    parse = TRUE, hjust = 0, size = 6 / ggplot2::.pt, family = "Helvetica"
  ) +
  annotate(
    "text", x = xy[[1]] + 0.04, y = xy[[2]] - 0.05,
    label = "11 T-cell programs: untestable", hjust = 0, vjust = 1,
    size = 6 / ggplot2::.pt, family = "Helvetica", color = "#555555"
  ) +
  scale_shape_manual(values = program_shapes, guide = "none") +
  scale_fill_manual(values = program_colors, guide = "none") +
  scale_x_continuous(limits = xy, breaks = pretty(xy, n = 4), expand = c(0, 0)) +
  scale_y_continuous(limits = xy, breaks = pretty(xy, n = 4), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  labs(x = "Raw stage effect", y = "Ambient-corrected stage effect") +
  theme_panel()

# 4E: two compact, donor-paired forests retain the full prespecified lineage
# family and use the atlas colors of the highest-scoring lineages.
lineage_order <- c(
  "Hepatocytes", "Fibroblasts", "Cholangiocytes",
  "Endothelial cells", "Macrophages", "T cells"
)
e_plot <- copy(e)
e_plot[, lineage_label := factor(comparison_lineage, levels = rev(lineage_order))]
e_plot[, program_label := fifelse(
  program_name == "ECM/IGFBP7",
  "'ECM'~'/'~italic(IGFBP7)",
  "'Ductular injury'~'/'~italic(BICC1)"
)]
e_plot[, target := comparison_lineage == target_lineage]
e_plot[, target_color := fifelse(
  program_name == "ECM/IGFBP7", program_colors[["ECM/IGFBP7"]],
  program_colors[["Ductular-injury/BICC1"]]
)]
e_plot[, target_shape := fifelse(program_name == "ECM/IGFBP7", 21L, 22L)]
p_e <- ggplot(e_plot, aes(y = lineage_label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_segment(
    data = e_plot[evidence_state != "reference"],
    aes(x = hc3_ci_low, xend = hc3_ci_high, yend = lineage_label),
    linewidth = 0.35, color = "#707070"
  ) +
  geom_point(
    data = e_plot[target == FALSE], aes(x = beta),
    shape = 21, size = 1.65, stroke = 0.25, fill = "white", color = "#707070"
  ) +
  geom_point(
    data = e_plot[target == TRUE], aes(x = beta, shape = program_name, fill = program_name),
    size = 2.25, stroke = 0.3, color = "black"
  ) +
  facet_grid(
    . ~ program_label, scales = "free_x",
    labeller = labeller(program_label = label_parsed)
  ) +
  scale_shape_manual(values = c("ECM/IGFBP7" = 21L, "Ductular-injury/BICC1" = 22L), guide = "none") +
  scale_fill_manual(values = program_colors, guide = "none") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0.04, 0.08))) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Lineage − hepatocyte score", y = NULL) +
  theme_panel() +
  theme(
    panel.spacing.x = unit(0.18, "inches"),
    axis.ticks.y = element_blank()
  )

# 4F: fixed-member bulk transport is shown in its native log2FC unit. Point
# area reports the fraction of fixed program weight supported by the complete
# per-contrast gene-level BH family.
f_plot <- copy(f)
f_plot[, stage_index := match(stage, c("F1", "F2", "F3", "F4"))]
f_plot[, program_key := fifelse(
  program_name == "ECM/IGFBP7", "ECM/IGFBP7", "Ductular-injury/BICC1"
)]
f_plot[, gene_label := fifelse(
  program_key == "ECM/IGFBP7", "italic(IGFBP7)", "italic(BICC1)"
)]
end_labels <- f_plot[stage == "F4"]
end_labels[, label_y := effect + fifelse(program_key == "ECM/IGFBP7", -0.015, 0.015)]
p_f <- ggplot(
  f_plot,
  aes(
    x = stage_index, y = effect, group = program_key,
    color = program_key, shape = program_key
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#BDBDBD") +
  geom_line(linewidth = 0.65) +
  geom_point(
    aes(size = bh_weight_fraction, fill = program_key),
    stroke = 0.3, color = "black"
  ) +
  geom_text(
    data = end_labels,
    aes(x = 4.08, y = label_y, label = gene_label),
    inherit.aes = FALSE,
    parse = TRUE, hjust = 0, size = 6 / ggplot2::.pt, family = "Helvetica"
  ) +
  scale_color_manual(values = program_colors, guide = "none") +
  scale_fill_manual(values = program_colors, guide = "none") +
  scale_shape_manual(values = program_shapes, guide = "none") +
  scale_size_area(
    limits = c(0, 1), max_size = 3.0,
    breaks = c(0.25, 0.50, 0.75),
    labels = percent_format(accuracy = 1),
    name = "BH-supported\nprogram weight"
  ) +
  scale_x_continuous(
    breaks = 1:4, labels = c("F1", "F2", "F3", "F4"),
    limits = c(0.85, 4.65), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.2, 0.3), limits = c(-0.03, 1.20),
    expand = c(0, 0)
  ) +
  labs(x = "Fibrosis stage versus F0", y = expression("Bulk member-gene " * log[2] * "FC")) +
  theme_panel() +
  theme(
    legend.position = c(0.24, 0.74),
    legend.justification = c(0, 0.5),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.height = unit(0.10, "inches"),
    legend.key.width = unit(0.12, "inches")
  ) +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "white", color = "#707070")))

panel_files <- c(
  fig4d = file.path(OUT, "panels/fig4d_ambient_effect_transport.pdf"),
  fig4e = file.path(OUT, "panels/fig4e_same_atlas_lineage_specificity.pdf"),
  fig4f = file.path(OUT, "panels/fig4f_bulk_tissue_state_transport.pdf")
)
save_fig(p_d, panel_files[["fig4d"]], width = 2.48, height = 2.43, dpi = 600)
save_fig(p_e, panel_files[["fig4e"]], width = 3.30, height = 2.43, dpi = 600)
save_fig(p_f, panel_files[["fig4f"]], width = 2.72, height = 2.43, dpi = 600)
check("panel_files_written", all(file.exists(panel_files)),
      sum(file.exists(panel_files)), length(panel_files))

caption <- c(
  "Figure 4D–F candidate. D, Raw versus ambient-corrected donor-level cross-sectional stage effects for the frozen program registry. Of 117 programs, 106 are transport-evaluable and 11 T-cell programs are untestable because integer-count transport does not reproduce their native non-integer substrate. Arrows mark the two recurring programs; neither paired correction-induced change passes ordinary or HC3 BH correction across 117 programs. E, Same-atlas donor-paired lineage-minus-hepatocyte contrasts for the ECM/IGFBP7 and ductular-injury/BICC1 programs. Error bars are HC3 95% confidence intervals. BH correction spans ten prespecified program-by-lineage contrasts. ECM/IGFBP7 is highest in fibroblasts (53 paired donors across four datasets), whereas ductular-injury/BICC1 is highest in cholangiocytes (51 paired donors across five datasets). These are descriptive localizations, not transcript-origin or cell-autonomy tests. F, Fixed program-member genes transported into bulk fibrosis-stage contrasts. Lines show the L1-weighted member-gene log2 fold-change for F1–F4 versus F0; point area is the fraction of available fixed program weight passing the complete 23,370-gene BH family within each contrast. This is cross-sectional tissue-state transport, not longitudinal progression or a program-level significance test.",
  "",
  "Color continuity: orange follows the fibroblast-localized ECM/IGFBP7 program and cyan follows the cholangiocyte-localized ductular-injury/BICC1 program from panel E into panel F."
)
writeLines(caption, file.path(OUT, "captions/figure4def.md"))

writeLines(c(
  "# Figure 4D–F multicellular-story candidate",
  "",
  "Candidate-only render. The current main PDFs are unchanged.",
  "",
  "The renderer uses the validated v8 ambient, cross-lineage, and bulk-transport tables without changing the frozen 117-program registry or any membership or weight. Panels use 6-point Helvetica, cairo PDF, useDingbats=FALSE, and individual-panel output for Illustrator assembly."
), file.path(OUT, "README.md"))

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
validation_table <- rbindlist(validation, fill = TRUE)
fwrite(validation_table, file.path(OUT, "validation_report.tsv"), sep = "\t")

metadata <- c(
  "{",
  paste0('  "candidate_id": "', basename(OUT), '",'),
  '  "release_state": "candidate_pending_visual_review",',
  '  "source_candidate": "aging-inspired-resource-rebuild-candidate-2026-08-15-v8",',
  '  "frozen_program_count": 117,',
  '  "seed": 20260817,',
  '  "font_size_pt": 6,',
  '  "pdf_device": "cairo_pdf",',
  '  "current_main_pdfs_overwritten": false',
  "}"
)
writeLines(metadata, file.path(OUT, "build_metadata.json"))

message("Rendered candidate panels in ", OUT)
