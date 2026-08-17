#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(42043)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(prefix) {
  hit <- args[startsWith(args, paste0("--", prefix, "="))]
  if (length(hit) != 1L) stop("Require exactly one --", prefix, "= argument")
  sub(paste0("^--", prefix, "="), "", hit)
}

INPUT <- normalizePath(value_arg("input"), mustWork = TRUE)
OUTPUT <- value_arg("output")
SPECIFICITY <- normalizePath(value_arg("specificity"), mustWork = TRUE)
if (file.exists(OUTPUT)) stop("Output already exists: ", OUTPUT)

MEMBERS_FILE <- file.path(INPUT, "credible_set_members_annotated.tsv")
ARCH_FILE <- file.path(INPUT, "credible_set_pip_architecture.tsv")
CONTEXT_FILE <- file.path(INPUT, "credible_set_regulatory_context.tsv")
VERDICT_FILE <- file.path(INPUT, "genetics_noncoding_verdict.tsv")
stopifnot(all(file.exists(c(MEMBERS_FILE, ARCH_FILE, CONTEXT_FILE, VERDICT_FILE))))

members <- fread(MEMBERS_FILE)
architecture <- fread(ARCH_FILE)
context <- fread(CONTEXT_FILE)
specificity <- fread(SPECIFICITY)

scope_labels <- c(
  tier1_direct_masld_pdff = "Direct MASLD / liver fat",
  tier2_liver_enzyme = "Liver enzymes"
)
scope_order <- unname(scope_labels)

members <- members[pip_sum_gate_passed == TRUE]
architecture <- architecture[pip_sum_gate_passed == TRUE]
context <- context[pip_sum_gate_passed == TRUE]
stopifnot(nrow(architecture) == uniqueN(members$credible_set_uid))
stopifnot(nrow(context) == nrow(architecture))
stopifnot(all(architecture$annotated_pip_mass >= 0.95))

members[, trait_scope_label := scope_labels[trait_scope]]
architecture[, trait_scope_label := scope_labels[trait_scope]]
context[, trait_scope_label := scope_labels[trait_scope]]
stopifnot(!anyNA(members$trait_scope_label))

members[, is_utr := grepl("(^|,)[35]_prime_UTR_variant($|,)", consequence)]
members[, is_synonymous := grepl("(^|,)synonymous_variant($|,)", consequence)]
members[, is_noncoding_dna := consequence_category == "other_noncoding" |
  (consequence_category == "synonymous_or_utr" & is_utr & !is_synonymous)]
members[, accessible_lineages := fifelse(
  lineage_accessibility_status == "resolved" & lineage_accessible == TRUE,
  lineage_cell_types,
  ""
)]

lineage_map <- list(
  Hepatocytes = c("Hepatocytes"),
  Fibroblasts = c("Fibroblasts"),
  Macrophages = c("Macrophages"),
  Cholangiocytes = c("Cholangiocytes"),
  Endothelial = c("Endothelial_cells", "Endothelial"),
  Lymphoid = c("T_cells", "Resident_NK", "Plasma_cells", "Lymphoid")
)
for (lineage in names(lineage_map)) {
  patterns <- paste(lineage_map[[lineage]], collapse = "|")
  members[, (paste0("accessible_", lineage)) :=
    lineage_accessibility_status == "resolved" & grepl(patterns, accessible_lineages)]
}

category_labels <- c(
  protein_altering_pip_mass = "Protein-altering",
  canonical_splice_pip_mass = "Canonical splice",
  synonymous_or_utr_pip_mass = "Synonymous / UTR",
  other_noncoding_pip_mass = "Other noncoding",
  unresolved_pip_mass = "Unresolved"
)
category_colors <- c(
  "Protein-altering" = "#C9265E",
  "Canonical splice" = "#8E3B76",
  "Synonymous / UTR" = "#56B4A9",
  "Other noncoding" = "#1565C0",
  "Unresolved" = "#9E9E9E"
)
for (column in names(category_labels)) {
  architecture[, (column) := as.numeric(get(column))]
}

arch_long <- melt(
  architecture,
  id.vars = c("credible_set_uid", "trait_scope", "trait_scope_label"),
  measure.vars = names(category_labels),
  variable.name = "consequence_category",
  value.name = "pip_mass"
)
arch_summary <- arch_long[, .(
  mean_pip_mass = mean(pip_mass),
  n_credible_sets = uniqueN(credible_set_uid)
), by = .(trait_scope, trait_scope_label, consequence_category)]
arch_summary[, consequence_label := category_labels[consequence_category]]
arch_summary[, trait_scope_label := factor(trait_scope_label, levels = rev(scope_order))]
arch_summary[, consequence_label := factor(
  consequence_label,
  levels = names(category_colors)
)]
arch_summary[, trait_scope_display := sprintf(
  "%s  (n=%s sets)", as.character(trait_scope_label),
  format(n_credible_sets, big.mark = ",", trim = TRUE)
)]
arch_summary[, trait_scope_display := factor(
  trait_scope_display,
  levels = c(
    "Liver enzymes  (n=6,208 sets)",
    "Direct MASLD / liver fat  (n=207 sets)"
  )
)]

arch_plot <- ggplot(
  arch_summary,
  aes(x = mean_pip_mass, y = trait_scope_display, fill = consequence_label)
) +
  geom_col(width = 0.66, color = "white", linewidth = 0.2) +
  geom_text(
    data = arch_summary[mean_pip_mass >= 0.025],
    aes(label = sprintf("%.0f%%", 100 * mean_pip_mass)),
    position = position_stack(vjust = 0.5),
    size = 6 / .pt,
    color = "white"
  ) +
  scale_fill_manual(values = category_colors, name = NULL, drop = FALSE) +
  scale_x_continuous(
    labels = function(x) sprintf("%.0f%%", 100 * x),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Mean SuSiE posterior mass per reliable credible set", y = NULL) +
  theme_masld() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

cs_noncoding <- members[is_noncoding_dna == TRUE, .(
  noncoding_pip_mass = sum(normalized_susie_pip),
  promoter_mass = sum(normalized_susie_pip[promoter_context_status == "resolved" & promoter_proximal == TRUE]),
  abc_mass = sum(normalized_susie_pip[abc_enhancer_overlap == TRUE]),
  accessible_mass = sum(normalized_susie_pip[lineage_accessibility_status == "resolved" & lineage_accessible == TRUE]),
  resolved_no_context_mass = sum(normalized_susie_pip[
    promoter_context_status == "resolved" &
      lineage_accessibility_status == "resolved" &
      promoter_proximal == FALSE &
      abc_enhancer_overlap == FALSE &
      lineage_accessible == FALSE
  ]),
  unresolved_context_mass = sum(normalized_susie_pip[
    promoter_context_status != "resolved" | lineage_accessibility_status != "resolved"
  ])
), by = .(credible_set_uid, trait_scope, trait_scope_label)]

all_cs <- unique(architecture[, .(credible_set_uid, trait_scope, trait_scope_label)])
cs_noncoding <- merge(all_cs, cs_noncoding, by = c(
  "credible_set_uid", "trait_scope", "trait_scope_label"
), all.x = TRUE)
mass_cols <- c(
  "noncoding_pip_mass", "promoter_mass", "abc_mass", "accessible_mass",
  "resolved_no_context_mass", "unresolved_context_mass"
)
for (column in mass_cols) set(cs_noncoding, which(is.na(cs_noncoding[[column]])), column, 0)

context_labels <- c(
  promoter_mass = "Within 2 kb of a GENCODE TSS",
  abc_mass = "Liver ABC enhancer",
  accessible_mass = "Open in at least one liver lineage",
  resolved_no_context_mass = "No deposited regulatory context",
  unresolved_context_mass = "Context unresolved"
)
context_summary <- rbindlist(lapply(names(context_labels), function(metric) {
  cs_noncoding[, .(
    value = sum(get(metric)) / sum(noncoding_pip_mass),
    n_credible_sets = .N,
    total_noncoding_pip_mass = sum(noncoding_pip_mass)
  ), by = .(trait_scope, trait_scope_label)][, metric := metric]
}))
context_summary[, context_label := context_labels[metric]]
context_summary[, trait_scope_label := factor(trait_scope_label, levels = scope_order)]
context_summary[, context_label := factor(context_label, levels = rev(unname(context_labels)))]

context_plot <- ggplot(context_summary, aes(
  x = trait_scope_label, y = context_label, fill = value
)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * value)), size = 6 / .pt) +
  scale_fill_gradient(
    low = "white", high = "#1565C0",
    limits = c(0, max(context_summary$value)),
    labels = function(x) sprintf("%.0f%%", 100 * x),
    name = "Noncoding\nPIP mass"
  ) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right"
  )

lineage_cols <- paste0("accessible_", names(lineage_map))
lineage_cs <- rbindlist(lapply(seq_along(lineage_cols), function(index) {
  column <- lineage_cols[index]
  lineage <- names(lineage_map)[index]
  values <- members[is_noncoding_dna == TRUE, .(
    lineage_mass = sum(normalized_susie_pip[get(column) == TRUE]),
    noncoding_pip_mass = sum(normalized_susie_pip)
  ), by = .(credible_set_uid, trait_scope, trait_scope_label)]
  values[, lineage := lineage]
  values
}))
lineage_summary <- lineage_cs[, .(
  value = sum(lineage_mass) / sum(noncoding_pip_mass),
  n_credible_sets = uniqueN(credible_set_uid),
  total_noncoding_pip_mass = sum(noncoding_pip_mass)
), by = .(trait_scope, trait_scope_label, lineage)]
lineage_summary[, trait_scope_label := factor(trait_scope_label, levels = scope_order)]
lineage_summary[, lineage := factor(lineage, levels = rev(names(lineage_map)))]

lineage_plot <- ggplot(lineage_summary, aes(
  x = trait_scope_label, y = lineage, fill = value
)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * value)), size = 6 / .pt) +
  scale_fill_gradient(
    low = "white", high = "#C9265E",
    limits = c(0, max(lineage_summary$value)),
    labels = function(x) sprintf("%.0f%%", 100 * x),
    name = "Noncoding PIP mass\nin open peaks"
  ) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right"
  )

member_pip <- members[, .(
  max_susie_pip = max(susie_pip),
  max_normalized_pip = max(normalized_susie_pip),
  direct_masld_member = any(trait_scope == "tier1_direct_masld_pdff")
), by = variant_id]
specificity <- merge(specificity, member_pip, by = "variant_id", all = FALSE)
specificity[, source_group := fifelse(
  direct_masld_member,
  "Direct MASLD / liver fat",
  "Liver-enzyme only"
)]
setorder(specificity, source_group, -max_susie_pip, variant_id)
top_variants <- specificity[, head(.SD, 15L), by = source_group]
specificity_cols <- paste0("norm_", names(lineage_map))
specificity_long <- melt(
  top_variants,
  id.vars = c("variant_id", "max_susie_pip", "direct_masld_member"),
  measure.vars = specificity_cols,
  variable.name = "lineage",
  value.name = "normalized_accessibility"
)
specificity_long[, lineage := sub("^norm_", "", lineage)]
specificity_long[, variant_label := sprintf(
  "%s  (PIP %.2f; %s)", variant_id, max_susie_pip,
  fifelse(direct_masld_member, "direct", "enzyme")
)]
variant_levels <- unique(specificity_long$variant_label)
specificity_long[, variant_label := factor(variant_label, levels = rev(variant_levels))]
specificity_long[, lineage := factor(lineage, levels = names(lineage_map))]

specificity_plot <- ggplot(specificity_long, aes(
  x = lineage, y = variant_label, fill = log1p(normalized_accessibility)
)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white", high = "#C9265E",
    labels = function(x) sprintf("%.1f", exp(x) - 1),
    name = "Accessibility /\nlineage background"
  ) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

distribution <- cs_noncoding[, .(
  credible_set_uid,
  trait_scope,
  trait_scope_label,
  noncoding_pip_mass
)]
distribution[, trait_scope_label := factor(trait_scope_label, levels = scope_order)]
distribution_plot <- ggplot(distribution, aes(
  x = noncoding_pip_mass, color = trait_scope_label
)) +
  stat_ecdf(linewidth = 0.6) +
  scale_color_manual(values = c(
    "Direct MASLD / liver fat" = "#C9265E",
    "Liver enzymes" = "#1565C0"
  ), name = NULL) +
  scale_x_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Noncoding posterior mass per credible set",
    y = "Cumulative fraction of credible sets"
  ) +
  theme_masld() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

write_tsv_once <- function(table, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  fwrite(table, path, sep = "\t", quote = FALSE, na = "")
}
md5_file <- function(path) {
  unname(tools::md5sum(path))
}

output_parent <- dirname(OUTPUT)
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
tmp <- file.path(
  output_parent,
  paste0(".", basename(OUTPUT), ".tmp.", Sys.getenv("SLURM_JOB_ID", Sys.getpid()))
)
if (file.exists(tmp)) stop("Temporary output exists: ", tmp)
dir.create(tmp, recursive = FALSE)

plots <- list(
  fig2_credible_set_pip_architecture = list(
    plot = arch_plot, width = 5.2, height = 2.45, source = arch_summary
  ),
  fig2_noncoding_regulatory_context = list(
    plot = context_plot, width = 4.3, height = 2.45, source = context_summary
  ),
  fig4_noncoding_variant_accessibility = list(
    plot = lineage_plot, width = 4.0, height = 2.55, source = lineage_summary
  ),
  figS_high_pip_variant_lineage_specificity = list(
    plot = specificity_plot, width = 5.8, height = 5.8, source = specificity_long
  ),
  figS_credible_set_noncoding_distribution = list(
    plot = distribution_plot, width = 3.65, height = 2.55, source = distribution
  )
)

for (name in names(plots)) {
  item <- plots[[name]]
  ggsave(
    file.path(tmp, paste0(name, ".pdf")), item$plot,
    width = item$width, height = item$height, units = "in",
    device = cairo_pdf
  )
  write_tsv_once(item$source, file.path(tmp, paste0(name, "_source.tsv")))
}

input_paths <- c(MEMBERS_FILE, ARCH_FILE, CONTEXT_FILE, VERDICT_FILE, SPECIFICITY)
write_tsv_once(data.table(
  path = normalizePath(input_paths),
  size_bytes = file.info(input_paths)$size,
  md5 = vapply(input_paths, md5_file, character(1))
), file.path(tmp, "input_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))
artifacts <- setdiff(list.files(tmp), "output_manifest.tsv")
write_tsv_once(data.table(
  relative_path = artifacts,
  size_bytes = file.info(file.path(tmp, artifacts))$size,
  md5 = vapply(file.path(tmp, artifacts), md5_file, character(1))
), file.path(tmp, "output_manifest.tsv"))
if (!file.rename(tmp, OUTPUT)) stop("Atomic candidate publication failed")
cat("NONCODING_DNA_PLOTS_COMPLETE\t", OUTPUT, "\n", sep = "")
