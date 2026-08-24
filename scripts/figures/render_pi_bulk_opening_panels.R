# KEY MESSAGE: The five-cohort fragment-count release separates cohort coverage, phenotype structure, and reproducible disease-state expression.
# Candidate Figure 3A-C renderer. Sourced by nas_fib_grid.R so the existing
# stage-distribution entrypoint remains the implementation route.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(patchwork)
})

base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
if (!nzchar(candidate_root) || !dir.exists(candidate_root)) stop("Active candidate root required", call. = FALSE)
source(file.path(base, "scripts/figures/publication_theme.R"))
stage_extension_root <- Sys.getenv(
  "STAGE_EXTENSION_ROOT",
  file.path(candidate_root, "analysis", "stage_extensions")
)
analysis_root <- normalizePath(stage_extension_root, mustWork = TRUE)
out_dir <- file.path(candidate_root, "figure3", "panels")
src_dir <- file.path(candidate_root, "source_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

samples <- fread(file.path(analysis_root, "five_cohort_sample_manifest.tsv"))
summary <- fread(file.path(analysis_root, "five_cohort_summary.tsv"))
cohorts <- summary[order(-n_participants), dataset]

# 3A: synchronized five-cohort composition plus metadata observability.
availability <- rbindlist(lapply(cohorts, function(cohort) {
  d <- samples[dataset == cohort]
  data.table(
    dataset = cohort,
    field = c("Control", "Sex", "NAS", "Fibrosis"),
    fraction = c(mean(d$group_binary == "Control"), mean(!is.na(d$inferred_sex)),
                 mean(!is.na(d$nas_score)), mean(!is.na(d$fibrosis_stage)))
  )
}))
availability[, state := cut(fraction, breaks = c(-Inf, 0, 0.99, Inf),
                            labels = c("Absent", "Partial", "Available"))]
availability[, dataset := factor(dataset, levels = rev(cohorts))]
availability[, field := factor(field, levels = c("Control", "Sex", "NAS", "Fibrosis"))]
p_meta <- ggplot(availability, aes(field, dataset, fill = state)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(Absent = "#ECECEC", Partial = cat_palette[[5]], Available = cat_palette[[1]]),
                    name = NULL) +
  labs(x = NULL, y = NULL) + theme_masld_compact() +
  theme(axis.text.x = element_text(size = 6, face = "plain", angle = 35, hjust = 1),
        axis.text.y = element_text(size = 6, face = "plain"), legend.position = "bottom")
composition <- summary[, .(dataset, Control = n_control, Disease = n_disease)]
composition <- melt(composition, id.vars = "dataset", variable.name = "group", value.name = "n")
composition[, dataset := factor(dataset, levels = rev(cohorts))]
p_bar <- ggplot(composition, aes(dataset, n, fill = group)) +
  geom_col(width = 0.72) + coord_flip() +
  scale_fill_manual(values = c(Control = "#9E9E9E", Disease = masld_colors$mash), name = NULL) +
  labs(x = NULL, y = "Participants") + theme_masld_compact() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
panel_a <- (p_meta | p_bar) + plot_layout(widths = c(1.35, 0.8), guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig3a_cohort_metadata_matrix.pdf"), panel_a,
       width = 3.0, height = 2.15, device = cairo_pdf)

# 3B: exact-stage x NAS table from the same five-cohort manifest. The coarse
# GSE213621 fibrosis bins are excluded from this exact Kleiner display.
grid_samples <- samples[dataset != "GSE213621" & fibrosis_stage %in% 0:4 & nas_score %in% 0:8]
grid <- merge(CJ(fibrosis_stage = 0:4, nas_score = 0:8),
              grid_samples[, .N, by = .(fibrosis_stage, nas_score)],
              by = c("fibrosis_stage", "nas_score"), all.x = TRUE)
grid[is.na(N), N := 0L]
p_grid <- ggplot(grid, aes(factor(fibrosis_stage), factor(nas_score), fill = N)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(data = grid[N > 0], aes(label = N, color = N > quantile(grid$N, 0.8)), size = 6 / .pt) +
  scale_color_manual(values = c(`FALSE` = "grey20", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient(low = "#EAF2FB", high = "#1565C0", name = "Participants") +
  labs(x = "Fibrosis stage", y = "NAS") + theme_masld_compact() +
  theme(axis.text = element_text(size = 6, face = "plain"), legend.position = "bottom")
ggsave(file.path(out_dir, "fig3b_nas_fib_grid.pdf"), p_grid,
       width = 2.2, height = 2.15, device = cairo_pdf)

# 3C: per-cohort DEG membership. The five complete BH families are kept
# separate; a gene enters a cohort stream at FDR<0.05 and |log2FC|>0.5.
cohort_de <- fread(file.path(analysis_root, "cohort_disease_all_gene_results.tsv"))
if (any(cohort_de[, uniqueN(gene_id_versioned), by = dataset]$V1 != unique(cohort_de$bh_family_size))) {
  stop("Incomplete per-cohort BH family", call. = FALSE)
}
hits <- cohort_de[FDR < 0.05 & abs(logFC) > 0.5, .(dataset, gene_id_versioned)]
hit_counts <- hits[, .(n_cohorts = uniqueN(dataset)), by = gene_id_versioned]
alluvial <- merge(hits, hit_counts, by = "gene_id_versioned")[, .N, by = .(dataset, n_cohorts)]
alluvial[, tier := factor(
  fifelse(n_cohorts == 1L, "Unique to study", paste0(n_cohorts, " cohorts")),
  levels = c("Unique to study", "2 cohorts", "3 cohorts", "4 cohorts", "5 cohorts")
)]
cohort_order_top_down <- alluvial[, .(N = sum(N)), by = dataset][order(-N), dataset]
tier_order_top_down <- c("Unique to study", "2 cohorts", "3 cohorts", "4 cohorts", "5 cohorts")
alluvial[, dataset := factor(dataset, levels = cohort_order_top_down)]
alluvial[, tier := factor(tier, levels = tier_order_top_down)]

# Restore the original Fig. 3C alluvial visual language at its original
# 364 x 285 pt footprint. Labels remain outside the streams and only the
# unique tier carries a count.
cohort_pal_pastel <- c(
  GSE126848 = "#a8c7e5", GSE130970 = "#c9b8e8", GSE135251 = "#e8bde9",
  GSE162694 = "#f3bec1", GSE213621 = "#fbdebd"
)
tier_pal_pastel <- c(
  `Unique to study` = "#D5D5D5", `2 cohorts` = "#e3bee8",
  `3 cohorts` = "#beace0", `4 cohorts` = "#a3acd3", `5 cohorts` = "#7f87b6"
)
left_dt <- alluvial[, .(N = sum(N)), by = dataset][order(dataset)]
total_y <- sum(left_dt$N)
left_dt[, y_top := total_y - cumsum(c(0, head(N, -1L)))]
left_dt[, y_bot := y_top - N]
left_dt[, y_mid := (y_top + y_bot) / 2]
right_dt <- alluvial[, .(N = sum(N)), by = tier][order(tier)]
right_dt[, y_top := total_y - cumsum(c(0, head(N, -1L)))]
right_dt[, y_bot := y_top - N]
right_dt[, y_mid := (y_top + y_bot) / 2]
unique_n <- right_dt[as.character(tier) == "Unique to study", N]

p_alluvial <- ggplot(alluvial, aes(y = N, axis1 = dataset, axis2 = tier)) +
  geom_alluvium(aes(fill = dataset), width = 1 / 10, alpha = 0.55,
                knot.pos = 0.42, curve_type = "sigmoid", linewidth = 0) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1 / 10,
               color = "white", linewidth = 0.45) +
  annotate("text", x = 1 - 0.085, y = left_dt$y_mid,
           label = as.character(left_dt$dataset), hjust = 1,
           size = 6 / .pt, fontface = "plain", color = "black") +
  annotate("text", x = 2 + 0.085, y = right_dt$y_mid,
           label = ifelse(as.character(right_dt$tier) == "Unique to study",
                          sprintf("Unique to study\n(%s)", format(unique_n, big.mark = ",")),
                          as.character(right_dt$tier)),
           hjust = 0, size = 6 / .pt, fontface = "plain", color = "black") +
  scale_x_discrete(limits = c("Study", "Overlap"), expand = c(0.30, 0.30), position = "top") +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  scale_fill_manual(values = c(cohort_pal_pastel, tier_pal_pastel), guide = "none") +
  labs(x = NULL, y = NULL) + theme_void() +
  theme(axis.text.x.top = element_text(size = 6, face = "plain", color = "black",
                                       margin = margin(b = 4)),
        plot.margin = margin(8, 14, 6, 14))
ggsave(file.path(out_dir, "fig3c_cohort_alluvial.pdf"), p_alluvial,
       width = 364 / 72, height = 285 / 72, device = cairo_pdf)

fwrite(summary, file.path(src_dir, "fig3a_five_cohort_summary.tsv"), sep = "\t")
fwrite(grid, file.path(src_dir, "fig3b_fibrosis_nas_grid.tsv"), sep = "\t")
fwrite(alluvial, file.path(src_dir, "fig3c_cohort_deg_replication.tsv"), sep = "\t")
cat("[saved] Figure 3A-B and main cohort-membership alluvial\n")
