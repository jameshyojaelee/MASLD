#!/usr/bin/env Rscript
# =============================================================================
# 08_finemapping_figures.R
# Publication-quality figures summarizing MASLD fine-mapping results
# (SuSiE + CARMA across 29 GWAS studies)
#
# Panels:
#   A — Fine-mapping overview heatmap (study x locus, colored by max PIP)
#   B — PNPLA3 cross-study concordance (PIP tracks per study)
#   C — SuSiE vs CARMA concordance scatter + per-locus Pearson r inset
#   D — Credible set size distribution (violin by convergence x ancestry)
#   E — Multi-signal loci barplot (# CS per locus, stacked by phenotype)
#   F — Cross-ancestry comparison at PNPLA3 (EUR vs EAS side-by-side)
#
# Output: figures/finemapping/panel_{A-F}.pdf + combined_finemapping.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_BASE <- file.path(BASE, "GWAS/finemapping")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS05_epigenomic_spatial")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PANELS_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANELS_DIR, recursive = TRUE, showWarnings = FALSE)

# Source publication theme
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# Illustrator-safe PDFs
grDevices::pdf.options(useDingbats = FALSE)

# ---------------------------------------------------------------------------
# 0. Read GWAS registry and build phenotype map
# ---------------------------------------------------------------------------
registry <- fread(file.path(FM_BASE, "config/gwas_registry.tsv"))

# Phenotype extraction from study_name
phenotype_map <- function(study) {
  case_when(
    grepl("NAFLD|NAFL", study, ignore.case = TRUE) ~ "NAFLD",
    grepl("NASH", study, ignore.case = TRUE) ~ "NASH",
    grepl("Cirrhosis", study, ignore.case = TRUE) ~ "Cirrhosis",
    grepl("HCC", study, ignore.case = TRUE) ~ "HCC",
    grepl("PDFF", study, ignore.case = TRUE) ~ "PDFF",
    grepl("_ALT$|_ALT_", study) ~ "ALT",
    grepl("_AST$|_AST_", study) ~ "AST",
    grepl("_GGT$|_GGT_", study) ~ "GGT",
    grepl("Obesity", study, ignore.case = TRUE) ~ "Obesity",
    TRUE ~ "Other"
  )
}

# Ordered phenotype factor levels
pheno_levels <- c("NAFLD", "NASH", "Cirrhosis", "HCC", "PDFF",
                  "ALT", "AST", "GGT", "Obesity")

# Phenotype colors (colorblind-friendly)
pheno_colors <- c(
  NAFLD     = "#0D47A1",
  NASH      = "#C2185B",
  Cirrhosis = "#880E4F",
  HCC       = "#E91E63",
  PDFF      = "#7B1FA2",
  ALT       = "#00695C",
  AST       = "#2E7D32",
  GGT       = "#F57F17",
  Obesity   = "#78909C"
)

# ---------------------------------------------------------------------------
# 1. Read all SuSiE results
# ---------------------------------------------------------------------------
cat("Reading SuSiE results...\n")

read_susie_locus <- function(tsv_path, study_name) {
  dt <- fread(tsv_path, sep = "\t")
  # Detect convergence from filename
  converged <- grepl("adjustedvar", basename(tsv_path))
  # Extract locus from the locus column
  dt$study <- study_name
  dt$susie_converged <- converged
  # Rename columns for consistency
  setnames(dt, old = c("PIP", "CS"), new = c("susie_PIP", "susie_CS"),
           skip_absent = TRUE)
  dt
}

all_susie <- list()
studies <- list.dirs(file.path(FM_BASE, "output"), recursive = FALSE, full.names = FALSE)
studies <- studies[studies != "gpfs"]  # exclude spurious dir

for (study in studies) {
  for (anc in c("EUR", "EAS")) {
    susie_dir <- file.path(FM_BASE, "output", study, paste0(anc, "_0.5Mb"), "susie")
    if (!dir.exists(susie_dir)) next
    tsv_files <- list.files(susie_dir, pattern = "\\.tsv$", full.names = TRUE)
    if (length(tsv_files) == 0) next
    for (f in tsv_files) {
      tryCatch({
        dt <- read_susie_locus(f, study)
        dt$ancestry <- anc
        all_susie[[length(all_susie) + 1]] <- dt
      }, error = function(e) {
        message("  WARN: Could not read SuSiE file: ", basename(f), " — ", e$message)
      })
    }
  }
}

susie_dt <- rbindlist(all_susie, fill = TRUE)
cat("  Total SuSiE variants:", nrow(susie_dt), "across",
    uniqueN(susie_dt$study), "studies\n")

# ---------------------------------------------------------------------------
# 2. Read all CARMA results
# ---------------------------------------------------------------------------
cat("Reading CARMA results...\n")

read_carma_locus <- function(path, study_name) {
  if (grepl("\\.txt\\.gz$", path)) {
    dt <- tryCatch(fread(cmd = paste("zcat", shQuote(path)), sep = "\t"),
                   error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt$study <- study_name
    setnames(dt, old = c("PIP", "CS"), new = c("carma_PIP", "carma_CS"),
             skip_absent = TRUE)
    return(dt)
  }
  NULL
}

read_carma_rds <- function(rds_path, ss_path, study_name) {
  carma_obj <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(carma_obj)) return(NULL)
  # Need to read corresponding SS file for variant identifiers
  if (!file.exists(ss_path)) return(NULL)
  ss <- fread(ss_path, sep = "\t")
  pips <- carma_obj[[1]]$PIPs
  if (length(pips) != nrow(ss)) return(NULL)
  ss$carma_PIP <- pips
  # Assign CS from CARMA credible sets
  ss$carma_CS <- 0L
  cs_list <- carma_obj[[1]][["Credible set"]]
  if (length(cs_list) >= 2 && length(cs_list[[2]]) > 0) {
    for (j in seq_along(cs_list[[2]])) {
      ss$carma_CS[cs_list[[2]][[j]]] <- j
    }
  }
  ss$study <- study_name
  ss
}

all_carma <- list()
for (study in studies) {
  for (anc in c("EUR", "EAS")) {
    carma_dir <- file.path(FM_BASE, "output", study, paste0(anc, "_0.5Mb"), "CARMA")
    ss_dir    <- file.path(FM_BASE, "output", study, paste0(anc, "_0.5Mb"), "ss")
    if (!dir.exists(carma_dir)) next

    gz_files <- list.files(carma_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
    rds_files <- list.files(carma_dir, pattern = "\\.rds$", full.names = TRUE)

    # Try txt.gz first; fall back to rds + ss if empty
    for (gz in gz_files) {
      tryCatch({
        dt <- read_carma_locus(gz, study)
        if (!is.null(dt) && nrow(dt) > 0) {
          dt$ancestry <- anc
          all_carma[[length(all_carma) + 1]] <- dt
        } else {
          # Fall back to rds
          # Extract locus from filename: {study}_{LDpanel}_{window}_{chr}.{pos}.txt.gz
          locus_str <- sub(".*_(\\d+\\.\\d+)\\.txt\\.gz$", "\\1", basename(gz))
          rds_match <- grep(locus_str, rds_files, value = TRUE)
          if (length(rds_match) == 1) {
            ss_match <- list.files(ss_dir, pattern = paste0("_", locus_str, "\\.txt$"),
                                   full.names = TRUE)
            if (length(ss_match) == 1) {
              dt2 <- read_carma_rds(rds_match, ss_match, study)
              if (!is.null(dt2) && nrow(dt2) > 0) {
                dt2$ancestry <- anc
                all_carma[[length(all_carma) + 1]] <- dt2
              }
            }
          }
        }
      }, error = function(e) {
        message("  WARN: Could not read CARMA: ", basename(gz), " — ", e$message)
      })
    }

    # Handle rds-only loci (no corresponding gz)
    for (rds in rds_files) {
      locus_str <- sub(".*_(\\d+\\.\\d+)\\.rds$", "\\1", basename(rds))
      gz_exists <- any(grepl(locus_str, gz_files))
      if (!gz_exists) {
        tryCatch({
          ss_match <- list.files(ss_dir, pattern = paste0("_", locus_str, "\\.txt$"),
                                 full.names = TRUE)
          if (length(ss_match) == 1) {
            dt2 <- read_carma_rds(rds, ss_match, study)
            if (!is.null(dt2) && nrow(dt2) > 0) {
              dt2$ancestry <- anc
              all_carma[[length(all_carma) + 1]] <- dt2
            }
          }
        }, error = function(e) {
          message("  WARN: Could not read CARMA rds: ", basename(rds), " — ", e$message)
        })
      }
    }
  }
}

carma_dt <- rbindlist(all_carma, fill = TRUE)
cat("  Total CARMA variants:", nrow(carma_dt), "across",
    uniqueN(carma_dt$study), "studies\n")

# ---------------------------------------------------------------------------
# 3. Build per-locus summary tables
# ---------------------------------------------------------------------------
cat("Building per-locus summaries...\n")

# Ensure SNP column exists in SuSiE data
if (!"SNP" %in% names(susie_dt) && "chromosome" %in% names(susie_dt)) {
  susie_dt[, SNP := paste(chromosome, position, allele1, allele2, sep = ":")]
}

# SuSiE per-locus summary
susie_locus <- susie_dt[, .(
  max_susie_PIP     = max(susie_PIP, na.rm = TRUE),
  n_variants        = .N,
  n_cs              = max(susie_CS, na.rm = TRUE),
  converged         = any(susie_converged, na.rm = TRUE),
  chr               = chromosome[1]
), by = .(study, locus, ancestry)]

# Add phenotype
susie_locus[, phenotype := phenotype_map(study)]
susie_locus[, phenotype := factor(phenotype, levels = pheno_levels)]

# CARMA per-locus summary
carma_locus <- carma_dt[, .(
  max_carma_PIP = max(carma_PIP, na.rm = TRUE),
  n_carma_cs    = max(carma_CS, na.rm = TRUE)
), by = .(study, locus, ancestry)]

# Merge
locus_summary <- merge(susie_locus, carma_locus,
                       by = c("study", "locus", "ancestry"), all.x = TRUE)

# Mark skipped loci
skip_files <- list.files(file.path(FM_BASE, "output"), pattern = "SKIPPED_",
                         recursive = TRUE, full.names = TRUE)
skipped_loci <- data.table(
  locus = sub("SKIPPED_", "", sub("\\.txt$", "", basename(skip_files)))
)
# Extract study from path
skipped_loci[, study := sapply(skip_files, function(p) {
  parts <- strsplit(p, "/")[[1]]
  idx <- which(parts == "output") + 1
  if (idx <= length(parts)) parts[idx] else NA_character_
})]

cat("  Locus summary:", nrow(locus_summary), "loci across",
    uniqueN(locus_summary$study), "studies\n")

# ---------------------------------------------------------------------------
# PANEL A: Fine-mapping overview heatmap
# ---------------------------------------------------------------------------
cat("Generating Panel A: Overview heatmap...\n")

# Add phenotype to locus summary for ordering
locus_summary[, phenotype := factor(phenotype_map(study), levels = pheno_levels)]

# Create study short labels
study_labels <- function(s) {
  s <- gsub("2019_31311600_", "", s)
  s <- gsub("2020_32298765_", "", s)
  s <- gsub("2020_32514122_", "", s)
  s <- gsub("2021_34128465_", "", s)
  s <- gsub("2021_34841290_", "", s)
  s <- gsub("2021_34957434_", "", s)
  s <- gsub("2022_36402844_", "", s)
  s <- gsub("2023_36280732_", "", s)
  s <- gsub("2023_36653562_", "", s)
  s <- gsub("_EUR$|_EAS$", "", s)
  s
}

# Order studies by phenotype then name
study_order <- locus_summary[, .(phenotype = phenotype[1]),
                              by = study][order(phenotype, study)]$study

locus_summary[, study_label := study_labels(study)]
locus_summary[, study_label := factor(study_label,
               levels = unique(study_labels(study_order)))]

# Create status column: converged (dot shape 16), not converged (X shape 4)
locus_summary[, status := fifelse(converged, "Converged", "Not converged")]

# For the heatmap, use chr as column grouping
locus_summary[, chr := as.integer(chr)]

# Aggregate: max PIP per study x chromosome
heatmap_dt <- locus_summary[, .(
  max_PIP  = max(max_susie_PIP, na.rm = TRUE),
  status   = fifelse(any(converged), "Converged", "Not converged"),
  n_loci   = .N
), by = .(study_label, chr, phenotype)]

# Cap infinities
heatmap_dt[is.infinite(max_PIP), max_PIP := 0]

p_A <- ggplot(heatmap_dt, aes(x = factor(chr), y = study_label)) +
  geom_tile(aes(fill = max_PIP), color = "white", linewidth = 0.3) +
  geom_point(aes(shape = status), size = 0.8, color = "black") +
  scale_fill_gradientn(
    colors = c("white", "#E3F2FD", "#42A5F5", "#1565C0", "#C2185B", "#880E4F"),
    values = rescale(c(0, 0.1, 0.3, 0.5, 0.8, 1.0)),
    limits = c(0, 1),
    name = "Max PIP"
  ) +
  scale_shape_manual(
    values = c("Converged" = 16, "Not converged" = 4),
    name = "SuSiE status"
  ) +
  labs(x = "Chromosome", y = NULL,
       title = "Fine-mapping overview across 29 MASLD GWAS") +
  theme_masld(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5),
    axis.text.y = element_text(size = 5),
    legend.position = "right",
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(0.5, "cm"),
    panel.grid = element_blank()
  )

save_fig(p_A, file.path(PANELS_DIR, "panel_A_overview_heatmap.pdf"),
         width = fig_full_width, height = 5.5)
cat("  Saved panels/panel_A_overview_heatmap.pdf\n")

# ---------------------------------------------------------------------------
# PANEL B: PNPLA3 cross-study concordance
# ---------------------------------------------------------------------------
cat("Generating Panel B: PNPLA3 cross-study concordance...\n")

# PNPLA3 region: chr22, positions ~43.8M - 44.9M (with 0.5Mb window)
# rs738409 is at position 44324727 or 44324855 (depending on build/study)
RS738409_POS <- 44324727  # canonical hg19 position

# Filter SuSiE for chr22 loci near PNPLA3
pnpla3_susie <- susie_dt[chromosome == 22 &
                           position >= 43800000 & position <= 44900000]

# Select top 8 studies by max PIP at PNPLA3
pnpla3_top_studies <- pnpla3_susie[, .(max_pip = max(susie_PIP, na.rm = TRUE)),
                                     by = study][order(-max_pip)][1:min(8, .N)]$study

pnpla3_plot <- pnpla3_susie[study %in% pnpla3_top_studies]
pnpla3_plot[, study_short := study_labels(study)]

# Also get CARMA for these studies
pnpla3_carma <- carma_dt[chromosome == 22 &
                           position >= 43800000 & position <= 44900000 &
                           study %in% pnpla3_top_studies]
pnpla3_carma[, study_short := study_labels(study)]

# Combine for plotting
pnpla3_combined <- rbind(
  pnpla3_plot[, .(position, PIP = susie_PIP, study_short, method = "SuSiE")],
  pnpla3_carma[, .(position, PIP = carma_PIP, study_short, method = "CARMA")]
)

# Convert position to Mb for cleaner axis
pnpla3_combined[, pos_mb := position / 1e6]

p_B <- ggplot(pnpla3_combined, aes(x = pos_mb, y = PIP, color = study_short)) +
  geom_point(aes(shape = method), size = 0.3, alpha = 0.6) +
  geom_vline(xintercept = RS738409_POS / 1e6, linetype = "dashed",
             color = "gray30", linewidth = 0.4) +
  annotate("text", x = RS738409_POS / 1e6, y = 1.02, label = "rs738409",
           size = 2, fontface = "italic", hjust = 0.5) +
  scale_color_viridis_d(option = "turbo", name = "Study") +
  scale_shape_manual(values = c("SuSiE" = 16, "CARMA" = 1), name = "Method") +
  labs(x = "Position (Mb, chr22)", y = "PIP",
       title = "PNPLA3 locus: cross-study fine-mapping concordance") +
  theme_masld(base_size = 7) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm")
  ) +
  coord_cartesian(ylim = c(0, 1.05))

save_fig(p_B, file.path(PANELS_DIR, "panel_B_pnpla3_concordance.pdf"),
         width = fig_full_width, height = 3.5)
cat("  Saved panels/panel_B_pnpla3_concordance.pdf\n")

# ---------------------------------------------------------------------------
# PANEL C: SuSiE vs CARMA concordance scatter
# ---------------------------------------------------------------------------
cat("Generating Panel C: SuSiE vs CARMA concordance...\n")

# Build SNP key in both datasets for merging
if (!"SNP" %in% names(susie_dt)) {
  susie_dt[, SNP := paste(chromosome, position, allele1, allele2, sep = ":")]
}
if (!"SNP" %in% names(carma_dt)) {
  carma_dt[, SNP := paste(chromosome, position, allele1, allele2, sep = ":")]
}

# Merge on study + locus + SNP
concordance_dt <- merge(
  susie_dt[, .(SNP, study, locus, ancestry, susie_PIP)],
  carma_dt[, .(SNP, study, locus, ancestry, carma_PIP)],
  by = c("SNP", "study", "locus", "ancestry"),
  all = FALSE
)

cat("  Concordance pairs:", nrow(concordance_dt), "\n")

# Per-locus Pearson r
locus_cor <- concordance_dt[, .(
  r = tryCatch(cor(susie_PIP, carma_PIP, use = "complete.obs"),
               error = function(e) NA_real_),
  n = .N
), by = .(study, locus, ancestry)]
locus_cor <- locus_cor[!is.na(r)]
cat("  Per-locus correlations:", nrow(locus_cor),
    "median r =", round(median(locus_cor$r, na.rm = TRUE), 3), "\n")

# Downsample for scatter (keep all PIP > 0.01, random 5% of rest)
set.seed(42)
high_pip <- concordance_dt[susie_PIP > 0.01 | carma_PIP > 0.01]
low_pip  <- concordance_dt[susie_PIP <= 0.01 & carma_PIP <= 0.01]
low_sample <- low_pip[sample(.N, min(.N, 20000))]
scatter_dt <- rbind(high_pip, low_sample)

p_C_main <- ggplot(scatter_dt, aes(x = susie_PIP, y = carma_PIP, color = ancestry)) +
  rasterize_layer(
    geom_point(size = 0.2, alpha = 0.3)
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray40", linewidth = 0.4) +
  scale_color_manual(
    values = c("EUR" = "#1565C0", "EAS" = "#F57F17"),
    name = "Ancestry"
  ) +
  labs(x = "SuSiE PIP", y = "CARMA PIP",
       title = "SuSiE vs CARMA concordance") +
  theme_masld(base_size = 7) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1))

# Inset: histogram of per-locus r
p_C_inset <- ggplot(locus_cor, aes(x = r)) +
  geom_histogram(fill = "#7B1FA2", color = "white",
                 bins = 30, linewidth = 0.2) +
  geom_vline(xintercept = median(locus_cor$r, na.rm = TRUE),
             linetype = "dashed", color = "#C2185B", linewidth = 0.4) +
  annotate("text",
           x = median(locus_cor$r, na.rm = TRUE) - 0.05,
           y = Inf, vjust = 1.5,
           label = paste0("median r = ", round(median(locus_cor$r, na.rm = TRUE), 3)),
           size = 1.8, color = "#C2185B") +
  labs(x = "Pearson r", y = "Count") +
  theme_masld(base_size = 5) +
  theme(plot.background = element_rect(fill = "white", color = "gray70", linewidth = 0.3),
        plot.margin = margin(2, 2, 2, 2))

p_C <- p_C_main +
  inset_element(p_C_inset,
                left = 0.55, bottom = 0.02, right = 0.98, top = 0.42)

save_fig(p_C, file.path(PANELS_DIR, "panel_C_concordance_scatter.pdf"),
         width = fig_half_width, height = fig_half_width)
cat("  Saved panels/panel_C_concordance_scatter.pdf\n")

# ---------------------------------------------------------------------------
# PANEL D: Credible set size distribution
# ---------------------------------------------------------------------------
cat("Generating Panel D: Credible set size distribution...\n")

# Count variants per credible set
# SuSiE CS: each variant has CS = integer (0 = not in CS, 1..N = CS index)
cs_sizes_susie <- susie_dt[susie_CS > 0, .(
  cs_size = .N
), by = .(study, locus, ancestry, susie_CS, converged = susie_converged)]

cs_sizes_susie[, convergence := fifelse(converged, "Converged", "Not converged")]

# CARMA CS
cs_sizes_carma <- carma_dt[carma_CS > 0, .(
  cs_size = .N
), by = .(study, locus, ancestry, carma_CS)]
cs_sizes_carma[, convergence := "CARMA"]

# Combine SuSiE only for the violin (cleaner story)
cs_violin_dt <- cs_sizes_susie[, .(cs_size, ancestry, convergence)]

# Cap at 200 for visualization
cs_violin_dt[cs_size > 200, cs_size := 200]

p_D <- ggplot(cs_violin_dt, aes(x = interaction(convergence, ancestry, sep = "\n"),
                                  y = cs_size, fill = ancestry)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3) +
  stat_summary(fun = median, geom = "point", size = 1.2, color = "white") +
  scale_fill_manual(values = c("EUR" = "#1565C0", "EAS" = "#F57F17"),
                    name = "Ancestry") +
  scale_y_log10(breaks = c(1, 5, 10, 50, 100, 200),
                labels = c("1", "5", "10", "50", "100", "200+")) +
  labs(x = NULL, y = "Credible set size (variants)",
       title = "SuSiE 95% credible set sizes") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(size = 6))

save_fig(p_D, file.path(PANELS_DIR, "panel_D_cs_sizes.pdf"),
         width = fig_half_width, height = 3)
cat("  Saved panels/panel_D_cs_sizes.pdf\n")

# ---------------------------------------------------------------------------
# PANEL E: Multi-signal loci barplot
# ---------------------------------------------------------------------------
cat("Generating Panel E: Multi-signal loci barplot...\n")

# Number of CS per locus (SuSiE, converged only)
multi_signal <- locus_summary[converged == TRUE, .(
  n_cs = pmax(n_cs, 0),
  phenotype
)]
multi_signal[, n_cs_bin := fifelse(n_cs >= 5, "5+", as.character(n_cs))]
multi_signal[, n_cs_bin := factor(n_cs_bin, levels = c("0", "1", "2", "3", "4", "5+"))]

bar_dt <- multi_signal[, .N, by = .(n_cs_bin, phenotype)]

p_E <- ggplot(bar_dt, aes(x = n_cs_bin, y = N, fill = phenotype)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pheno_colors, name = "Phenotype") +
  labs(x = "Number of credible sets per locus",
       y = "Number of loci",
       title = "Multi-signal fine-mapping loci") +
  theme_masld(base_size = 7) +
  theme(legend.position = "right")

save_fig(p_E, file.path(PANELS_DIR, "panel_E_multisignal.pdf"),
         width = fig_half_width, height = 3)
cat("  Saved panels/panel_E_multisignal.pdf\n")

# ---------------------------------------------------------------------------
# PANEL F: Cross-ancestry comparison at PNPLA3
# ---------------------------------------------------------------------------
cat("Generating Panel F: Cross-ancestry PNPLA3 comparison...\n")

# EUR representative: Ghodsian NAFLD (2021_34841290_NAFLD_EUR)
# EAS representative: BBJ_ALT
eur_study <- "2021_34841290_NAFLD_EUR"
eas_study <- "BBJ_ALT"

pnpla3_eur_susie <- susie_dt[study == eur_study & chromosome == 22 &
                               position >= 43800000 & position <= 44900000]
pnpla3_eas_susie <- susie_dt[study == eas_study & chromosome == 22 &
                               position >= 43800000 & position <= 44900000]
pnpla3_eur_carma <- carma_dt[study == eur_study & chromosome == 22 &
                               position >= 43800000 & position <= 44900000]
pnpla3_eas_carma <- carma_dt[study == eas_study & chromosome == 22 &
                               position >= 43800000 & position <= 44900000]

# Combine into long format
panel_f_dt <- rbind(
  pnpla3_eur_susie[, .(position, PIP = susie_PIP,
                         panel = "EUR (Ghodsian NAFLD)", method = "SuSiE")],
  pnpla3_eur_carma[, .(position, PIP = carma_PIP,
                         panel = "EUR (Ghodsian NAFLD)", method = "CARMA")],
  pnpla3_eas_susie[, .(position, PIP = susie_PIP,
                         panel = "EAS (BBJ ALT)", method = "SuSiE")],
  pnpla3_eas_carma[, .(position, PIP = carma_PIP,
                         panel = "EAS (BBJ ALT)", method = "CARMA")]
)

panel_f_dt[, pos_mb := position / 1e6]

# Find disagreement: where SuSiE and CARMA top variants differ substantially in EAS
eas_top_susie <- pnpla3_eas_susie[which.max(susie_PIP), .(position, susie_PIP)]
eas_top_carma <- pnpla3_eas_carma[which.max(carma_PIP), .(position, carma_PIP)]

p_F <- ggplot(panel_f_dt, aes(x = pos_mb, y = PIP, color = method)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_vline(xintercept = RS738409_POS / 1e6, linetype = "dashed",
             color = "gray40", linewidth = 0.3) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c("SuSiE" = "#1565C0", "CARMA" = "#C2185B"),
    name = "Method"
  ) +
  # Mark disagreement in EAS panel
  {
    if (nrow(eas_top_susie) > 0 && nrow(eas_top_carma) > 0 &&
        abs(eas_top_susie$position - eas_top_carma$position) > 10000) {
      list(
        annotate("segment",
                 x = eas_top_susie$position / 1e6, xend = eas_top_carma$position / 1e6,
                 y = 0.95, yend = 0.95,
                 arrow = arrow(ends = "both", length = unit(0.05, "cm")),
                 color = "#E91E63", linewidth = 0.3),
        annotate("text",
                 x = mean(c(eas_top_susie$position, eas_top_carma$position)) / 1e6,
                 y = 1.0,
                 label = "Method\ndiscordance",
                 size = 1.8, color = "#E91E63", fontface = "italic")
      )
    }
  } +
  labs(x = "Position (Mb, chr22)", y = "PIP",
       title = "Cross-ancestry fine-mapping: PNPLA3 locus") +
  annotate("text", x = RS738409_POS / 1e6 + 0.01, y = 0.02,
           label = "rs738409", size = 1.8, color = "gray40",
           fontface = "italic", hjust = 0) +
  theme_masld(base_size = 7) +
  coord_cartesian(ylim = c(0, 1.05))

save_fig(p_F, file.path(PANELS_DIR, "panel_F_cross_ancestry_pnpla3.pdf"),
         width = fig_full_width, height = 4)
cat("  Saved panels/panel_F_cross_ancestry_pnpla3.pdf\n")

# ---------------------------------------------------------------------------
# Combined figure
# ---------------------------------------------------------------------------
cat("Assembling combined figure...\n")

# Layout: top row = A (full width)
#         middle row = B (full width)
#         bottom left = C + D, bottom right = E + F
combined <- (p_A / p_B / ((p_C | p_D) / (p_E | p_F))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(combined, file.path(OUT_DIR, "combined_finemapping.pdf"),
              width = fig_full_width, height = 16)
cat("  Saved combined_finemapping.pdf\n")

# ---------------------------------------------------------------------------
# Summary statistics table
# ---------------------------------------------------------------------------
cat("\n=== Fine-mapping Figure Summary ===\n")
cat("Total studies:", uniqueN(locus_summary$study), "\n")
cat("Total loci (SuSiE):", nrow(locus_summary), "\n")
cat("  Converged:", sum(locus_summary$converged), "\n")
cat("  Not converged:", sum(!locus_summary$converged), "\n")
cat("Multi-signal loci (>1 CS):", sum(locus_summary$n_cs > 1, na.rm = TRUE), "\n")
cat("SuSiE-CARMA concordance pairs:", nrow(concordance_dt), "\n")
cat("Per-locus median r:", round(median(locus_cor$r, na.rm = TRUE), 3), "\n")
cat("Median CS size (SuSiE converged):",
    median(cs_sizes_susie[converged == TRUE]$cs_size), "\n")
cat("\nAll panels saved to:", OUT_DIR, "\n")
cat("Done.\n")
