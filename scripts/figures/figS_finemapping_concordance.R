#!/usr/bin/env Rscript
# figS_finemapping_concordance.R
# Cross-study fine-mapping concordance beyond PNPLA3
# 3 panels: (A) multi-locus top-variant heatmap, (B) PIP concordance at shared loci,
#           (C) credible set overlap rates
# Output: figures/supplementary/figS05_epigenomic_spatial/figS_finemapping_concordance.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
FM <- file.path(BASE, "GWAS/finemapping")
outdir <- FIGS05_DIR
dir.create(FIGS05_DIR, showWarnings = FALSE, recursive = TRUE)
grDevices::pdf.options(useDingbats = FALSE)

# ---------------------------------------------------------------------------
# Load all SuSiE + CARMA results
# ---------------------------------------------------------------------------
cat("Loading fine-mapping results...\n")
susie_files <- c(
  Sys.glob(file.path(FM, "output/*/EUR_0.5Mb/susie/*_cov0.95_*.tsv")),
  Sys.glob(file.path(FM, "output/*/EUR_0.5Mb/susie/*_notconverged.tsv")),
  Sys.glob(file.path(FM, "output/*/EAS_0.5Mb/susie/*_cov0.95_*.tsv")),
  Sys.glob(file.path(FM, "output/*/EAS_0.5Mb/susie/*_notconverged.tsv"))
)

all_susie <- rbindlist(lapply(susie_files, function(f) {
  dt <- fread(f, select = c("chromosome", "position", "allele1", "allele2",
                             "beta", "se", "pval", "locus", "converged", "PIP", "CS"))
  parts <- strsplit(f, "/")[[1]]
  dt$study <- parts[which(parts == "output") + 1]
  dt$ancestry <- ifelse(grepl("EAS", f), "EAS", "EUR")
  dt
}), fill = TRUE)

# ARCHIVED 2026-04-09: Whitfield 36653562 removed (provenance unverified); Anstee removed (duplicate)
# RESTORED 2026-04-09: FinnGen_NAFLD, FinnGen_NASH restored (verified R12 provenance)
# DROPPED 2026-06-06: FinnGen_HCC and all cirrhosis/HCC + PanUKBB sex-strat duplicates removed
#                     (23-GWAS MASLD portfolio; etiology-mixed / duplicate, Broadaway-excluded)
archived_studies <- c("Anstee_NAFLD",
                      "2023_36653562_NAFLD_EUR", "2023_36653562_NASH_EUR",
                      "2023_36653562_Cirrhosis_EUR", "2023_36653562_HCC_EUR",
                      "2023_36653562_Obesity_EUR",
                      # cirrhosis/HCC dropped 2026-06-06 from the 23-GWAS MASLD portfolio
                      # (etiology-mixed, not MASLD-specific; Broadaway-excluded)
                      "FinnGen_HCC", "Ghouse_Cirrhosis", "Ghouse_HCC",
                      "2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS",
                      "Sveinbjornsson2022_Cirrhosis_meta_EUR", "Sveinbjornsson2022_Hcc_meta_EUR",
                      # PanUKBB sex-stratified / both-sex EUR duplicates dropped 2026-06-06
                      "PanUKBB_F_ALT", "PanUKBB_F_AST", "PanUKBB_F_GGT",
                      "PanUKBB_M_ALT", "PanUKBB_M_AST", "PanUKBB_M_GGT",
                      "PanUKBB_BS_ALT", "PanUKBB_BS_AST", "PanUKBB_BS_GGT")
all_susie <- all_susie[!study %in% archived_studies]
cat("Loaded", nrow(all_susie), "variants from", uniqueN(all_susie$study), "studies (after removing archived)\n")

# Also load CARMA for concordance
carma_files <- c(
  Sys.glob(file.path(FM, "output/*/EUR_0.5Mb/CARMA/*.txt.gz")),
  Sys.glob(file.path(FM, "output/*/EAS_0.5Mb/CARMA/*.txt.gz"))
)

all_carma <- rbindlist(lapply(carma_files, function(f) {
  dt <- fread(f, select = c("chromosome", "position", "PIP", "CS"))
  parts <- strsplit(f, "/")[[1]]
  dt$study <- parts[which(parts == "output") + 1]
  dt
}), fill = TRUE)
setnames(all_carma, "PIP", "CARMA_PIP")
setnames(all_carma, "CS", "CARMA_CS")
all_carma <- all_carma[!study %in% archived_studies]
cat("Loaded", nrow(all_carma), "CARMA variants (after removing archived)\n")

# Phenotype map
phenotype_map <- function(s) {
  if (grepl("ALT", s)) return("ALT")
  if (grepl("AST", s)) return("AST")
  if (grepl("GGT", s)) return("GGT")
  if (grepl("PDFF", s)) return("PDFF")
  if (grepl("NAFLD", s)) return("NAFLD")
  if (grepl("NASH", s)) return("NASH")
  if (grepl("Cirrhosis", s)) return("Cirrhosis")
  if (grepl("HCC", s)) return("HCC")
  if (grepl("Obesity", s)) return("Obesity")
  return("Other")
}
all_susie[, phenotype := sapply(study, phenotype_map)]

# Short study labels
study_short <- function(s) {
  # Biobank/method-based naming convention (2026-04-09)
  s <- gsub("2019_31311600_NAFLD_EUR", "eMERGE NAFLD", s)
  # ARCHIVED 2026-04-08: Anstee NAFLD (2020_32298765_NAFLD_EUR) removed — filtered above
  s <- gsub("2020_32298765_NAFLD_EUR", "Biopsy NAFLD", s)
  s <- gsub("2021_34841290_NAFLD_EUR", "EHR NAFLD", s)
  s <- gsub("2021_34128465_PDFF_EUR", "Abdominal MRI PDFF", s)
  s <- gsub("2021_34957434_PDFF_EUR", "ML-derived PDFF", s)
  s <- gsub("2022_36402844_PDFF_EUR", "Whole-body MRI PDFF", s)
  # ARCHIVED 2026-04-09: Whitfield 36653562 removed (provenance unverified)
  s <- gsub("2023_36280732_NAFLD_deCode_EUR", "deCODE NAFLD", s)
  s <- gsub("2023_36280732_NAFLD_Intermountain_EUR", "Intermountain NAFLD", s)
  s <- gsub("2023_36280732_NAFLD_UKBB_EUR", "UKBB NAFLD", s)
  # DROPPED 2026-06-06: 2020_32514122_*_EAS (BBJ Cirr/HCC) and Ghouse relabels removed with the cirrhosis/HCC GWAS
  # RESTORED 2026-04-09: FinnGen label mapping
  s <- gsub("FinnGen_", "FinnGen ", s)
  # ARCHIVED 2026-04-09: Pazoki_PDFF removed — duplicate of 2022_36402844_PDFF_EUR (same GCST90267352)
  s <- gsub("UKBB_", "UKBB ", s)
  s <- gsub("BBJ_", "BBJ ", s)
  return(s)
}

# ---------------------------------------------------------------------------
# Identify top shared loci with gene annotations
# ---------------------------------------------------------------------------
# Known gene annotations for major MASLD loci
locus_genes <- data.table(
  region = c("chr22:44.3", "chr8:126.5", "chr4:88.2", "chr19:19.5",
             "chr19:45.4", "chr1:221", "chr2:27.8", "chr12:113.3",
             "chr19:19.4", "chr4:100.2", "chr6:32.7", "chr1:150.3",
             "chr9:136.1", "chr1:155.2", "chr19:41.3", "chr11:93.9",
             "chr10:113.9"),
  gene = c("PNPLA3", "TRIB1/LPL", "TM6SF2/NCAN", "TM4SF4/GAMT",
           "ERLIN1/CHUK", "LYPLAL1", "GCKR", "HNF1A",
           "SUGP1/NCAN", "MTTP", "HLA", "IL6R",
           "ABO", "FMO3", "IFNL4", "APOA5/BUD13",
           "TCF7L2")
)

# Get top-PIP variant per locus per study
top_variants <- all_susie[, .(
  top_pos = position[which.max(PIP)],
  top_pip = max(PIP),
  top_allele1 = allele1[which.max(PIP)],
  top_allele2 = allele2[which.max(PIP)],
  n_in_cs = sum(CS > 0),
  converged = converged[1]
), by = .(study, locus)]

top_variants[, chr := as.integer(sub("[.].*", "", locus))]
top_variants[, region := paste0("chr", chr, ":", round(top_pos / 1e6, 1))]
top_variants[, study_short := sapply(study, study_short)]
top_variants[, phenotype := sapply(study, phenotype_map)]

# Find regions with 3+ studies
region_study_count <- top_variants[, .(n_studies = uniqueN(study)), by = region][order(-n_studies)]
shared_regions <- region_study_count[n_studies >= 3]$region

cat("Shared regions (3+ studies):", length(shared_regions), "\n")

# Merge gene annotations
top_shared <- top_variants[region %in% shared_regions]
top_shared <- merge(top_shared, locus_genes, by = "region", all.x = TRUE)
top_shared[is.na(gene), gene := region]
top_shared[, locus_label := paste0(gene, " (", region, ")")]

# Order loci by number of studies
locus_order <- top_shared[, .(n = uniqueN(study)), by = locus_label][order(-n)]$locus_label
top_shared[, locus_label := factor(locus_label, levels = rev(locus_order))]

# Order studies by phenotype
study_order <- top_shared[, .N, by = .(study_short, phenotype)][order(phenotype, study_short)]$study_short
top_shared[, study_short := factor(study_short, levels = study_order)]

# ---------------------------------------------------------------------------
# Panel A: Top-variant PIP heatmap across studies and loci
# ---------------------------------------------------------------------------
cat("Building Panel A...\n")

pA <- ggplot(top_shared, aes(x = study_short, y = locus_label, fill = top_pip)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(aes(shape = converged), size = 0.8, color = "black") +
  scale_fill_gradientn(
    colors = c("gray95", "#BBDEFB", "#1565C0", "#880E4F"),
    values = scales::rescale(c(0, 0.3, 0.7, 1.0)),
    limits = c(0, 1),
    name = "Top variant\nPIP",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("TRUE" = "Converged", "FALSE" = "Not conv."),
                     name = "SuSiE") +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right")
message("[caption] Cross-study fine-mapping: top variant PIP per locus")

# ---------------------------------------------------------------------------
# Panel B: Per-variant PIP concordance across study pairs at shared loci
# ---------------------------------------------------------------------------
cat("Building Panel B...\n")

# For each shared region, get variant-level PIPs across all studies that have it
# Then compute pairwise study correlations
shared_loci_list <- top_shared[, unique(locus)]

# Get all variants in shared loci
shared_variants <- all_susie[locus %in% shared_loci_list,
                             .(study, chromosome, position, PIP, locus)]
shared_variants[, chr := as.integer(sub("[.].*", "", locus))]
shared_variants[, region := paste0("chr", chr, ":", round(position / 1e6, 1))]

# Compute pairwise correlations per locus
locus_cors <- list()
for (loc in shared_loci_list) {
  loc_data <- shared_variants[locus == loc]
  studies_here <- unique(loc_data$study)
  if (length(studies_here) < 2) next

  # Wide format: position x study
  wide <- dcast(loc_data, position ~ study, value.var = "PIP", fun.aggregate = max)
  study_cols <- setdiff(names(wide), "position")
  if (length(study_cols) < 2) next

  # Pairwise correlations
  for (i in 1:(length(study_cols) - 1)) {
    for (j in (i + 1):length(study_cols)) {
      x <- wide[[study_cols[i]]]
      y <- wide[[study_cols[j]]]
      both_present <- !is.na(x) & !is.na(y)
      if (sum(both_present) < 10) next
      r <- cor(x[both_present], y[both_present], method = "pearson")
      locus_cors[[length(locus_cors) + 1]] <- data.table(
        locus = loc,
        study1 = study_cols[i], study2 = study_cols[j],
        r = r, n = sum(both_present)
      )
    }
  }
}
pair_cors <- rbindlist(locus_cors)
# locus is like "22.44324855", region is like "chr22:44.3" — merge via region
pair_cors[, chr_num := as.integer(sub("[.].*", "", locus))]
pair_cors[, locus_pos := as.numeric(sub("^[0-9]+[.]", "", locus))]
pair_cors[, region := paste0("chr", chr_num, ":", round(locus_pos / 1e6, 1))]
pair_cors <- merge(pair_cors, locus_genes, by = "region", all.x = TRUE)

cat("Pairwise study correlations:", nrow(pair_cors), "pairs across",
    uniqueN(pair_cors$locus), "loci\n")
cat("Median r:", round(median(pair_cors$r, na.rm = TRUE), 3), "\n")

# Classify pairs by phenotype concordance
pair_cors[, pheno1 := sapply(study1, phenotype_map)]
pair_cors[, pheno2 := sapply(study2, phenotype_map)]
pair_cors[, pair_type := ifelse(pheno1 == pheno2, "Same phenotype", "Cross-phenotype")]

# Use density + rug instead of histogram for small N
med_r <- median(pair_cors$r, na.rm = TRUE)

pB <- ggplot(pair_cors, aes(x = r)) +
  geom_histogram(aes(fill = pair_type), binwidth = 0.1, color = "white",
                 linewidth = 0.2, position = "stack") +
  geom_rug(aes(color = pair_type), alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = med_r, linetype = "dashed", color = "red", linewidth = 0.4) +
  annotate("text", x = med_r - 0.02, y = Inf, vjust = 1.5, hjust = 1,
           label = paste0("median r = ", round(med_r, 3)),
           size = GEOM_TEXT_6PT, color = "red") +
  annotate("text", x = -0.15, y = Inf, vjust = 1.5, hjust = 0,
           label = paste0("n = ", nrow(pair_cors), " pairs"),
           size = GEOM_TEXT_6PT, color = "gray40") +
  scale_fill_manual(values = c("Same phenotype" = "#1565C0",
                                "Cross-phenotype" = "#F57F17"),
                    name = NULL) +
  scale_color_manual(values = c("Same phenotype" = "#1565C0",
                                 "Cross-phenotype" = "#F57F17"),
                     guide = "none") +
  scale_x_continuous(limits = c(-0.25, 1.05), breaks = seq(-0.2, 1, 0.2)) +
  labs(x = "Pearson r (variant-level PIP)", y = "Number of study pairs") +
  theme_masld(base_size = 6) +
  theme(legend.position = c(0.25, 0.85),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())
message("[caption] PIP concordance across study pairs at shared loci")

# ---------------------------------------------------------------------------
# Panel C: Credible set overlap at shared loci
# ---------------------------------------------------------------------------
cat("Building Panel C...\n")

# For each shared locus, compute CS overlap between study pairs
cs_overlaps <- list()
for (loc in shared_loci_list) {
  loc_data <- all_susie[locus == loc & CS > 0]
  studies_here <- unique(loc_data$study)
  if (length(studies_here) < 2) next

  for (i in 1:(length(studies_here) - 1)) {
    for (j in (i + 1):length(studies_here)) {
      cs1 <- loc_data[study == studies_here[i]]$position
      cs2 <- loc_data[study == studies_here[j]]$position
      if (length(cs1) == 0 || length(cs2) == 0) next
      overlap <- length(intersect(cs1, cs2))
      union_size <- length(union(cs1, cs2))
      jaccard <- overlap / union_size
      cs_overlaps[[length(cs_overlaps) + 1]] <- data.table(
        locus = loc,
        study1 = studies_here[i], study2 = studies_here[j],
        cs1_size = length(cs1), cs2_size = length(cs2),
        overlap = overlap, jaccard = jaccard
      )
    }
  }
}
cs_overlap_dt <- rbindlist(cs_overlaps)

# Annotate with gene names and region
cs_overlap_dt[, chr := as.integer(sub("[.].*", "", locus))]
loci_info <- all_susie[locus %in% shared_loci_list,
                       .(med_pos = median(position)), by = locus]
cs_overlap_dt <- merge(cs_overlap_dt, loci_info, by = "locus", all.x = TRUE)
cs_overlap_dt[, region := paste0("chr", chr, ":", round(med_pos / 1e6, 1))]
cs_overlap_dt <- merge(cs_overlap_dt, locus_genes, by = "region", all.x = TRUE)
cs_overlap_dt[is.na(gene), gene := region]

cat("CS overlap pairs:", nrow(cs_overlap_dt), "\n")
cat("Median Jaccard:", round(median(cs_overlap_dt$jaccard, na.rm = TRUE), 3), "\n")

# Per-locus summary
cs_per_locus <- cs_overlap_dt[, .(
  median_jaccard = median(jaccard, na.rm = TRUE),
  mean_jaccard = mean(jaccard, na.rm = TRUE),
  n_pairs = .N,
  median_cs_size = median(c(cs1_size, cs2_size))
), by = .(gene)][order(-median_jaccard)]

# Filter to loci with >=3 pairs for meaningful statistics
cs_plot <- cs_per_locus[n_pairs >= 3]

pC <- ggplot(cs_plot, aes(x = reorder(gene, median_jaccard), y = median_jaccard)) +
  geom_col(aes(fill = median_jaccard), width = 0.7) +
  geom_text(aes(label = paste0("n=", n_pairs, " pairs")), hjust = -0.05, size = GEOM_TEXT_6PT) +
  scale_fill_gradientn(
    colors = c("#FFCDD2", "#E53935", "#880E4F"),
    limits = c(0, 1), name = "Jaccard",
    guide = "none"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
  coord_flip() +
  labs(x = NULL, y = "Median Jaccard index (credible set overlap)") +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(size = 6))
message("[caption] Credible set concordance across studies (loci with >=3 study pairs; ",
        nrow(cs_plot), " loci shown)")

# ---------------------------------------------------------------------------
# Save individual panels
# ---------------------------------------------------------------------------
panels_dir <- outdir
save_fig(pB, file.path(panels_dir, "panel_G_pip_concordance.pdf"),
         width = fig_half_width, height = 3.5)
cat("Saved: panels/panel_G_pip_concordance.pdf\n")

save_fig(pC, file.path(panels_dir, "panel_H_cs_overlap.pdf"),
         width = fig_half_width, height = 3.5)
cat("Saved: panels/panel_H_cs_overlap.pdf\n")

# Summary stats
cat("\n=== Summary ===\n")
cat("Shared loci (3+ studies):", length(shared_regions), "\n")
cat("Total study pairs tested:", nrow(pair_cors), "\n")
cat("Median PIP correlation:", round(median(pair_cors$r, na.rm = TRUE), 3), "\n")
cat("Median CS Jaccard:", round(median(cs_overlap_dt$jaccard, na.rm = TRUE), 3), "\n")
cat("Loci with Jaccard > 0.5:", sum(cs_per_locus$median_jaccard > 0.5), "/",
    nrow(cs_per_locus), "\n")
