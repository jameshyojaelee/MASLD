#!/usr/bin/env Rscript
# KEY MESSAGE: k-means on XIST/DDX3Y log1p-CPM cleanly separates F/M samples
# and concords with annotated sex in the 8 cohorts that ship sex labels.
# GSE135251 (Govaere) and GSE213621 (Chen) USE this k-means as their sex
# call, so concordance = 100% by construction (flagged).
#
# Outputs:
#   figS01/panels/figS01_sex_inference.pdf  (2-panel scatter + concordance)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIGS01_DIR, "panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dge <- load_merged_dge()
stopifnot(!is.null(dge))

# ----- log1p-CPM for XIST + DDX3Y -----
XIST_ID  <- grep("^ENSG00000229807", rownames(dge$counts), value = TRUE)[1]
DDX3Y_ID <- grep("^ENSG00000067048", rownames(dge$counts), value = TRUE)[1]
stopifnot(!is.na(XIST_ID), !is.na(DDX3Y_ID))

cpm_mat <- edgeR::cpm(dge, log = FALSE)
sex_dt <- data.table(
  sample_id  = rownames(dge$samples),
  dataset    = dge$samples$dataset,
  sex_stored = dge$samples$sex,
  XIST_logcpm  = log1p(cpm_mat[XIST_ID,  ]),
  DDX3Y_logcpm = log1p(cpm_mat[DDX3Y_ID, ])
)
sex_dt[sex_stored == "" | is.na(sex_stored), sex_stored := NA_character_]

# k-means on (XIST, DDX3Y) log1p-CPM, 2 centers
km <- kmeans(sex_dt[, .(XIST_logcpm, DDX3Y_logcpm)], centers = 2, nstart = 25)
sex_dt[, km_cluster := km$cluster]
# Female cluster = higher mean XIST
female_cluster <- which.max(tapply(sex_dt$XIST_logcpm, sex_dt$km_cluster, mean))
sex_dt[, sex_inferred := fifelse(km_cluster == female_cluster, "F", "M")]
sex_dt[, sex_inferred := factor(sex_inferred, levels = c("F", "M"))]
sex_dt[, sex_stored   := factor(sex_stored,   levels = c("F", "M"))]

# Source: "Inferred" if this cohort has no stored sex label for any sample
# (covers both the cohorts explicitly inferred upstream — Govaere/Chen — and
# Verschuren which has no SRA sex annotation at all).
inferred_datasets <- sex_dt[, .(any_stored = any(!is.na(sex_stored))),
                            by = dataset][any_stored == FALSE, dataset]
sex_dt[, sex_source := fifelse(dataset %in% inferred_datasets,
                               "Inferred (k-means)", "Annotated (SRA)")]

# ----- Panel A: scatter -----
COLS <- c(F = masld_colors$female, M = masld_colors$male)
p_a <- ggplot(sex_dt, aes(DDX3Y_logcpm, XIST_logcpm, color = sex_inferred,
                          shape = sex_stored)) +
  geom_point(size = 0.7, alpha = 0.75) +
  scale_color_manual(values = COLS, name = "k-means", drop = FALSE,
                     na.value = "grey60",
                     labels = c(F = "Female", M = "Male")) +
  scale_shape_manual(values = c(F = 16, M = 4), name = "Stored label",
                     drop = FALSE, na.value = 1,
                     labels = c(F = "Female", M = "Male"),
                     na.translate = TRUE) +
  labs(x = "DDX3Y log1p-CPM", y = "XIST log1p-CPM") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_masld()

# ----- Panel B: per-cohort concordance bar -----
sex_dt[, concordant := !is.na(sex_stored) & sex_inferred == sex_stored]
concord <- sex_dt[, .(
  n             = .N,
  n_stored      = sum(!is.na(sex_stored)),
  n_concordant  = sum(concordant, na.rm = TRUE),
  n_discordant  = sum(!is.na(sex_stored) & !concordant),
  pct_concord   = round(100 * sum(concordant, na.rm = TRUE) /
                          pmax(sum(!is.na(sex_stored)), 1L), 1)
), by = .(dataset, sex_source)]

cat("=== Sex inference concordance per cohort ===\n")
print(concord[order(-pct_concord)])

# Reshape for stacked bar
bar_dt <- melt(concord,
               id.vars        = c("dataset", "sex_source", "n_stored",
                                  "pct_concord"),
               measure.vars   = c("n_concordant", "n_discordant"),
               variable.name  = "agreement", value.name = "N")
bar_dt[, agreement := factor(
  fifelse(agreement == "n_concordant", "Concordant", "Discordant"),
  levels = c("Concordant", "Discordant"))]
bar_dt[, dataset := factor(dataset,
                           levels = concord[order(pct_concord)]$dataset)]

p_b <- ggplot(bar_dt, aes(x = dataset, y = N, fill = agreement)) +
  geom_col(width = 0.75) +
  geom_text(data = concord,
            aes(x = dataset, y = n_stored + 3,
                label = ifelse(sex_source == "Inferred (k-means)",
                               "k-means used",
                               sprintf("%.0f%%", pct_concord))),
            inherit.aes = FALSE, size = GEOM_TEXT_6PT, hjust = 0) +
  coord_flip() +
  scale_fill_manual(values = c("Concordant" = "#26A69A",
                                "Discordant" = "#E53935")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = NULL, y = "Samples with stored sex label",
       fill = NULL) +
  theme_masld() +
  theme(legend.position = "bottom")

# ----- Compose + save -----
fig <- (p_a | p_b) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin(2, 2, 2, 2)))

out_pdf <- file.path(OUT_DIR, "figS01_sex_inference.pdf")
out_csv <- file.path(FIGS01_DIR, "figS01_sex_inference_concordance.csv")
ggsave(out_pdf, fig, width = fig_full_width, height = 4.0, device = cairo_pdf)
fwrite(concord, out_csv)
message(sprintf("[caption] a: XIST/DDX3Y k-means clustering (n=%s). b: Annotated vs k-means concordance per cohort.",
                format(nrow(sex_dt), big.mark = ",")))
cat("\nWrote:\n  ", out_pdf, "\n  ", out_csv, "\n", sep = "")
