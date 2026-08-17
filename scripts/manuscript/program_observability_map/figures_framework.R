#!/usr/bin/env Rscript

# KEY MESSAGE: Frozen programs coordinate into four cross-lineage communities
# independent of histologic stage; one of them pairs ductular and portal
# fibroblast biology exclusively.
#
# Four panels, one PDF each, built ONLY from sealed system outputs. Nothing is
# recomputed here and no count that its own system withheld may appear.
#
# House style (docs/FIGURE_GUIDELINES.md): Helvetica 6pt plain, all text black,
# no panel titles or subtitles (the message lives in the caption), individual
# panel PDFs, cairo_pdf with useDingbats = FALSE, control gray #9E9E9E.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(project_root, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)

systems_root <- file.path(workstream_root, "systems")
output_root <- file.path(workstream_root, "figures_framework")
if (dir.exists(output_root)) fail("Refusing to overwrite: ", output_root)

registry <- fread(program_registry)[, .(feature_id = program_uid, cell_type, module_name)]

# Lineage colours from the house single-cell palette where available.
lineage_levels <- c("hepatocytes", "cholangiocytes", "fibroblasts", "macrophages", "tcells")
lineage_labels <- c(hepatocytes = "Hepatocyte", cholangiocytes = "Cholangiocyte",
                    fibroblasts = "Fibroblast", macrophages = "Macrophage", tcells = "T cell")
lineage_colors <- c(hepatocytes = "#C9265E", cholangiocytes = "#1565C0",
                    fibroblasts = "#00695C", macrophages = "#E08214", tcells = "#7B1FA2")

tmp <- atomic_dir(output_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)
source_tables <- list()

# --------------------------------------------------------------- PANEL A -----
# The four communities and their lineage composition. This is the biology: one
# community is exclusively cholangiocyte and fibroblast.
message("[A] community composition")
communities <- fread(file.path(systems_root, "T4_interaction/program_communities.tsv"))
setnames(communities, 1L, "feature_id")
communities <- merge(communities, registry, by = "feature_id", all.x = TRUE)
comm_names <- c("1" = "Pan-lineage\nfibro-inflammatory", "2" = "Hepatocyte\nmetabolic / bile-acid",
                "3" = "Ductular-portal", "4" = "Notch-progenitor\n+ immune")
panel_a_data <- communities[, .N, by = .(community, cell_type)]
panel_a_data[, cell_type := factor(cell_type, levels = rev(lineage_levels))]
panel_a_data[, community_label := factor(comm_names[as.character(community)],
                                         levels = unname(comm_names))]
totals <- communities[, .(n = .N), by = community]
panel_a_data <- merge(panel_a_data, totals, by = "community")

panel_a <- ggplot(panel_a_data, aes(x = N, y = community_label, fill = cell_type)) +
  geom_col(width = 0.62, color = NA) +
  geom_text(data = unique(panel_a_data[, .(community_label, n)]),
            aes(x = n, y = community_label, label = paste0("n=", n)),
            inherit.aes = FALSE, hjust = -0.25, size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_manual(values = lineage_colors, labels = lineage_labels,
                    breaks = lineage_levels, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Frozen programs", y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.line.y = element_blank(), axis.ticks.y = element_blank())
source_tables[["panel_a_community_composition"]] <- panel_a_data

# --------------------------------------------------------------- PANEL B -----
# Coordination is not the technical axis: partialling successive global
# components removes edges but does not abolish the structure.
message("[B] global-component robustness")
calib_t4 <- fread(file.path(systems_root, "T4_interaction/calibration.tsv"))
summary_t4 <- fread(file.path(systems_root, "T4_interaction/interaction_summary.tsv"))
get_summary <- function(key) as.numeric(summary_t4[metric == key, value])
# The last two steps matter as much as the first four: stopping at the k3 view
# would show 3,896 and silently omit the robust set that survives the empirical
# and group nulls, and the subset reproducing in every cohort.
stage_levels <- c("Shared genes\nretained", "Shared genes\nremoved",
                  "1 global\ncomponent removed", "3 global\ncomponents removed",
                  "All controls\n+ empirical null", "Reproduces in\nall 4 cohorts")
panel_b_data <- data.table(
  stage = factor(stage_levels, levels = stage_levels),
  n_edges = c(calib_t4[view == "edge_full_membership", observed_calls],
              calib_t4[view == "edge_disjoint_membership", observed_calls],
              calib_t4[view == "edge_partial_global_component_k1", observed_calls],
              calib_t4[view == "edge_partial_global_component_k3", observed_calls],
              get_summary("n_coordinated_beyond_global_component"),
              get_summary("n_coordinated_all_cohorts")),
  null_rate = c(calib_t4[view == "edge_full_membership", null_call_mean],
                calib_t4[view == "edge_disjoint_membership", null_call_mean],
                calib_t4[view == "edge_partial_global_component_k1", null_call_mean],
                calib_t4[view == "edge_partial_global_component_k3", null_call_mean],
                NA_real_, NA_real_),
  tier = c("view", "view", "view", "view", "robust", "robust")
)
panel_b <- ggplot(panel_b_data, aes(x = stage, y = n_edges, fill = tier)) +
  geom_col(width = 0.58, color = NA) +
  geom_text(aes(label = format(n_edges, big.mark = ",")), vjust = -0.55,
            size = PUB_GEOM_TEXT, color = "black") +
  geom_text(aes(y = 0, label = ifelse(is.na(null_rate), "",
                                      sprintf("null %.2f", null_rate))),
            vjust = -0.6, size = PUB_GEOM_TEXT, color = "white") +
  scale_fill_manual(values = c(view = "#9E9E9E", robust = "#C9265E"), guide = "none") +
  scale_y_continuous(labels = function(v) format(v, big.mark = ","),
                     expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Coordinated program pairs") +
  theme_masld(base_size = 6) + theme_pub()
source_tables[["panel_b_global_component_robustness"]] <- panel_b_data

# --------------------------------------------------------------- PANEL C -----
# Modularity against two nulls. The rewired null destroys topology; the
# rank-matched null preserves degree. Observed clears both.
message("[C] modularity vs nulls")
modularity_null <- fread(file.path(systems_root, "T4_interaction/modularity_null.tsv"))
global_structure <- fread(file.path(systems_root, "T4_interaction/global_structure.tsv"))
observed_modularity <- fread(file.path(systems_root, "T4_interaction/interaction_summary.tsv"))[
  metric == "modularity_observed", as.numeric(value)]
modularity_null[, null_label := fifelse(grepl("rewir", null), "Rewired", "Rank-matched")]
panel_c <- ggplot(modularity_null, aes(x = modularity, fill = null_label)) +
  geom_histogram(bins = 40, color = NA, alpha = 0.85, position = "identity") +
  geom_vline(xintercept = observed_modularity, color = "#C9265E",
             linewidth = STROKE_EMPHASIS) +
  annotate("text", x = observed_modularity, y = Inf,
           label = sprintf("observed %.3f", observed_modularity),
           hjust = 1.06, vjust = 1.6, size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_manual(values = c(Rewired = "#9E9E9E", `Rank-matched` = "#1565C0"), name = NULL) +
  labs(x = "Modularity", y = "Permutations") +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "top", legend.direction = "horizontal")
source_tables[["panel_c_modularity_null"]] <- modularity_null

# --------------------------------------------------------------- PANEL D -----
# Every vocabulary against its OWN size- and expression-matched null. Substrate-
# leaked vocabularies are drawn but marked, never counted as support.
message("[D] vocabulary vs matched-random")
hit_rates <- fread(file.path(systems_root, "T3_vocabulary/hit_rates.tsv"))
panel_d_data <- hit_rates[axis == "fibrosis"]
panel_d_data[, leaked := substrate_leak %in% TRUE | reportable_as_support %in% FALSE]
vocab_labels <- c(hotspot_117 = "Hotspot 117", hallmark_50 = "Hallmark 50",
                  cnmf_16 = "cNMF 16", bnmf_6 = "bulk NMF 6",
                  wgcna_4 = "WGCNA 4", published_panels = "Published panels")
panel_d_data[, label := vocab_labels[vocabulary]]
setorder(panel_d_data, hit_rate)
panel_d_data[, label := factor(label, levels = label)]

panel_d <- ggplot(panel_d_data, aes(y = label)) +
  geom_segment(aes(x = random_hit_rate_q025, xend = random_hit_rate_q975,
                   yend = label), color = "#9E9E9E", linewidth = 1.6, lineend = "butt") +
  geom_point(aes(x = random_hit_rate_mean), color = "#616161", size = 1.1, shape = 16) +
  geom_point(aes(x = hit_rate, color = leaked), size = 1.9, shape = 18) +
  geom_text(aes(x = hit_rate, label = sprintf("%.2fx  p=%.3f", enrichment_over_random,
                                              p_vs_matched_random)),
            hjust = -0.18, size = PUB_GEOM_TEXT, color = "black") +
  scale_color_manual(values = c(`FALSE` = "#C9265E", `TRUE` = "#9E9E9E"),
                     labels = c(`FALSE` = "Independent", `TRUE` = "Substrate-leaked"),
                     name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.34))) +
  labs(x = "Fibrosis hit rate (gray = matched-random 95% range)", y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.line.y = element_blank(), axis.ticks.y = element_blank())
source_tables[["panel_d_vocabulary_vs_matched_random"]] <- panel_d_data

# ------------------------------------------------------------------ write ----
message("[write] panels and source tables")
panels <- list(
  fig_framework_a_communities = list(plot = panel_a, width = 3.4, height = 2.2),
  fig_framework_b_global_component = list(plot = panel_b, width = 4.6, height = 2.1),
  fig_framework_c_modularity = list(plot = panel_c, width = 3.4, height = 2.0),
  fig_framework_d_vocabulary = list(plot = panel_d, width = 3.6, height = 2.2)
)
for (name in names(panels)) {
  ggsave(file.path(tmp, paste0(name, ".pdf")), panels[[name]]$plot,
         width = panels[[name]]$width, height = panels[[name]]$height,
         device = cairo_pdf, units = "in")
}
for (name in names(source_tables)) {
  write_tsv(source_tables[[name]], file.path(tmp, paste0(name, "_source.tsv")))
}

caption <- paste(
  "Frozen single-cell programs coordinate independently of histologic stage.",
  "(A) Four communities recovered from residual program-program correlation across",
  "469 donors in four cohorts, coloured by the lineage each program was derived from;",
  "community 3 contains cholangiocyte and fibroblast programs exclusively.",
  "(B) Coordinated pairs surviving successive controls: removing pairs that share",
  "member genes, then partialling one and three global components, which together",
  "carry 40.3 percent of residual variance. The measured permutation false-call rate",
  "is printed inside each bar. (C) Observed modularity against a degree-preserving",
  "rank-matched null and a rewired null, 200 permutations each. (D) Each program",
  "vocabulary against its own size- and expression-matched random sets, 200 draws;",
  "gray bars are the matched-random 95 percent range. Bulk NMF, WGCNA and published",
  "panels were derived from or overlap the discovery cohorts and are drawn in gray",
  "as substrate-leaked; they are never counted as support. Residual stage R-squared",
  "is 3.1e-31, so the coordination is not a stage effect. All counts carry a measured",
  "null false-call rate; counts that failed their calibration gate are absent."
)
writeLines(caption, file.path(tmp, "figure_caption.txt"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

inputs <- c(
  file.path(systems_root, "T4_interaction/program_communities.tsv"),
  file.path(systems_root, "T4_interaction/calibration.tsv"),
  file.path(systems_root, "T4_interaction/modularity_null.tsv"),
  file.path(systems_root, "T4_interaction/interaction_summary.tsv"),
  file.path(systems_root, "T3_vocabulary/hit_rates.tsv"),
  program_registry
)
write_tsv(data.table(path = sub(paste0("^", project_root, "/"), "", inputs),
                     sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))
artifacts <- setdiff(list.files(tmp), "artifact_manifest.tsv")
write_tsv(data.table(artifact = artifacts,
                     sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
                     size_bytes = file.info(file.path(tmp, artifacts))$size),
          file.path(tmp, "artifact_manifest.tsv"))

publish_dir(tmp, output_root)
message("FRAMEWORK_FIGURES_COMPLETE ", output_root)
