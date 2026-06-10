#!/usr/bin/env Rscript
# =============================================================================
# 43_zonation_classification.R
# Module B1: Hepatocyte Zonation Classification
#
# Classifies DEGs by periportal (PP) vs. pericentral (PC) hepatocyte program
# membership using Halpern 2017 zonation profiles + spatial validation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(biomaRt)
})

# Ensure dplyr select takes priority over data.table
select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

atlas_path    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
ortholog_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/ortholog_mapping.tsv")
dis_reg_path  <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
hep_reg_path  <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")

# Spatial zonation data (from GSE192741 Visium)
spatial_zonation_path <- file.path(BASE, "Analysis/Spatial/results/zonation/deg_zonation_classification.csv")

# Disease signatures
nafl_nash_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")

out_dir <- file.path(BASE, "RNA-seq/results/zonation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module B1: Hepatocyte Zonation Classification ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Build Zonation Reference from Halpern 2017
# =============================================================================
cat("--- Building zonation reference ---\n")

# Halpern 2017 canonical zonation markers (curated from Table S1, Cell 2017)
# Pericentral (Zone 3, PC) markers — drug metabolism, lipogenesis
pc_markers <- c(
  # Drug metabolism (CYP450)
  "Cyp2e1", "Cyp1a2", "Cyp2f2", "Cyp2c29", "Cyp2c37", "Cyp2c50",
  "Cyp2c54", "Cyp2c70", "Cyp3a11", "Cyp3a25", "Cyp3a41a",
  # Glutamine synthesis
  "Glul", "Rhbg", "Oat",
  # Wnt targets / PC markers
  "Axin2", "Lgr5", "Tbx3", "Rnf43",
  # Bile acid synthesis
  "Cyp7a1", "Cyp8b1", "Cyp27a1",
  # Lipogenesis
  "Fasn", "Acly", "Scd1", "Acaca",
  # Other PC
  "Gstm1", "Gstm2", "Gstm3", "Gsta3",
  "Aldh1a1", "Aldh3a2", "Hsd17b6"
)

# Periportal (Zone 1, PP) markers — urea cycle, gluconeogenesis, oxidative
pp_markers <- c(
  # Urea cycle
  "Cps1", "Otc", "Ass1", "Asl", "Arg1",
  # Gluconeogenesis
  "Pck1", "G6pc", "Fbp1", "Pcx",
  # Albumin/secretory
  "Alb", "Ttr", "Serpina1e", "Serpina1a", "Serpina1b",
  # Complement
  "C2", "C3", "C4b", "C5", "C8a", "C8b", "C8g", "C9",
  # Oxidative metabolism
  "Hal", "Hgd", "Fah", "Sds",
  # Amino acid catabolism
  "Aox3", "Aox1", "Aldob",
  # Other PP
  "Hnf4a", "Sult1a1", "Igfbp1",
  # Bile acid uptake
  "Slco1a1", "Slco1b2", "Slc10a1"
)

# Map mouse to human orthologs
cat("Loading ortholog map...\n")
ortho_map <- fread(ortholog_path)
cat("  Orthologs: ", nrow(ortho_map), " pairs\n")

# Also use Ensembl for comprehensive mapping
# The ortholog map has mouse_symbol and human_symbol columns
ortho_cols <- names(ortho_map)
cat("  Ortholog columns: ", paste(ortho_cols, collapse = ", "), "\n")

# Standardize column names (ortholog_mapping.tsv has mouse_symbol, human_symbol already)
if ("mouse_gene_name" %in% ortho_cols) {
  ortho_map <- ortho_map %>% rename(mouse_symbol = mouse_gene_name, human_symbol = human_gene_name)
} else if ("human_symbol" %in% ortho_cols && "mouse_symbol" %in% ortho_cols) {
  # Already correct
} else if ("Mouse.gene.name" %in% ortho_cols) {
  ortho_map <- ortho_map %>% rename(mouse_symbol = Mouse.gene.name, human_symbol = Gene.name)
}

# Create mouse→human mapping (case-insensitive for mouse symbols)
mouse2human <- ortho_map %>%
  select(mouse_symbol, human_symbol) %>%
  distinct() %>%
  filter(!is.na(mouse_symbol) & !is.na(human_symbol) &
         mouse_symbol != "" & human_symbol != "")

cat("  Unique mouse→human pairs: ", nrow(mouse2human), "\n")

# Map PC markers to human
map_to_human <- function(mouse_genes, m2h) {
  # Try exact match first
  mapped <- m2h %>%
    filter(mouse_symbol %in% mouse_genes) %>%
    pull(human_symbol) %>%
    unique()

  # Try case-insensitive
  if (length(mapped) < length(mouse_genes) / 2) {
    mapped_ci <- m2h %>%
      filter(tolower(mouse_symbol) %in% tolower(mouse_genes)) %>%
      pull(human_symbol) %>%
      unique()
    mapped <- unique(c(mapped, mapped_ci))
  }

  mapped
}

pc_human <- map_to_human(pc_markers, mouse2human)
pp_human <- map_to_human(pp_markers, mouse2human)

cat(sprintf("  PC markers mapped to human: %d / %d\n", length(pc_human), length(pc_markers)))
cat(sprintf("  PP markers mapped to human: %d / %d\n", length(pp_human), length(pp_markers)))

# Extended zonation reference: add well-known human markers
# (to compensate for incomplete ortholog mapping)
pc_human_extra <- c("CYP2E1", "CYP1A2", "GLUL", "AKR1B10",
                     "CYP3A4", "CYP2C9", "CYP2C19", "CYP2D6",
                     "FASN", "ACLY", "SCD", "ACACA", "ACACB",
                     "TBX3", "AXIN2", "CYP7A1", "CYP8B1",
                     "ALDH1A1", "GSTA1", "GSTM1",
                     "ABCC2", "AKR1C4", "ADH1B", "ADH4")

pp_human_extra <- c("CPS1", "OTC", "ASS1", "ASL", "ARG1",
                     "PCK1", "G6PC", "FBP1", "PCX",
                     "ALB", "TTR", "SERPINA1",
                     "C2", "C3", "C4A", "C4B", "C5", "C8A", "C9",
                     "HAL", "HGD", "FAH", "SDS",
                     "HNF4A", "SLC10A1", "SLCO1B1", "SLCO1B3",
                     "ALDOB", "IGFBP1", "SULT1A1")

pc_human <- unique(c(pc_human, pc_human_extra))
pp_human <- unique(c(pp_human, pp_human_extra))

# Remove any overlap (mark as ambiguous)
overlap <- intersect(pc_human, pp_human)
if (length(overlap) > 0) {
  cat(sprintf("  Removing %d ambiguous genes in both sets: %s\n",
              length(overlap), paste(overlap, collapse = ", ")))
  pc_human <- setdiff(pc_human, overlap)
  pp_human <- setdiff(pp_human, overlap)
}

cat(sprintf("  Final PC markers (human): %d\n", length(pc_human)))
cat(sprintf("  Final PP markers (human): %d\n", length(pp_human)))

# =============================================================================
# 2. Build Continuous Zonation Score
# =============================================================================
cat("\n--- Computing zonation scores ---\n")

# Use extended Halpern-derived markers to assign a continuous score
# Score = 1 (fully PC) to -1 (fully PP), 0 = non-zoned/mid-zonal

# For genes in our atlas, classify based on marker membership
atlas <- fread(atlas_path)
cat("  Atlas: ", nrow(atlas), " genes\n")

atlas <- atlas %>%
  mutate(
    is_deg = !is.na(bulk_padj) & bulk_padj < 0.1,  # Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.3 (Script 05b)
    is_deg_up = is_deg & !is.na(bulk_logFC) & bulk_logFC > 0,
    is_deg_down = is_deg & !is.na(bulk_logFC) & bulk_logFC < 0,
    is_pc = human_symbol %in% pc_human,
    is_pp = human_symbol %in% pp_human,
    zonation_class = case_when(
      is_pc ~ "Pericentral",
      is_pp ~ "Periportal",
      TRUE ~ "Non-zoned"
    )
  )

cat("\nZonation classification:\n")
cat(sprintf("  Pericentral: %d genes\n", sum(atlas$is_pc)))
cat(sprintf("  Periportal: %d genes\n", sum(atlas$is_pp)))
cat(sprintf("  Non-zoned: %d genes\n", sum(atlas$zonation_class == "Non-zoned")))

# =============================================================================
# 3. Integrate Spatial Zonation Data
# =============================================================================
cat("\n--- Integrating spatial zonation ---\n")

if (file.exists(spatial_zonation_path)) {
  spatial_zon <- fread(spatial_zonation_path)
  cat("  Spatial zonation data: ", nrow(spatial_zon), " genes\n")
  cat("  Columns: ", paste(names(spatial_zon), collapse = ", "), "\n")

  # Merge spatial zonation scores into atlas
  spatial_cols <- intersect(names(spatial_zon),
    c("gene", "zonation_class", "kruskal_stat", "kruskal_pval",
      "spearman_rho", "spearman_pval",
      "mean_PP1", "mean_PP2", "mean_Mid", "mean_PC2", "mean_PC1"))

  spatial_slim <- spatial_zon %>%
    select(any_of(spatial_cols)) %>%
    rename(spatial_zonation_class = zonation_class)

  if ("gene" %in% names(spatial_slim)) {
    atlas <- atlas %>%
      left_join(spatial_slim, by = c("human_symbol" = "gene"))
    cat("  Merged spatial zonation for ",
        sum(!is.na(atlas$spatial_zonation_class)), " genes\n")
  }
} else {
  cat("  Spatial zonation file not found\n")
}

# =============================================================================
# 4. Enrichment Testing
# =============================================================================
cat("\n--- Enrichment testing ---\n")

run_fisher <- function(label, zone_genes, deg_genes, universe) {
  in_zone_deg <- sum(universe %in% zone_genes & universe %in% deg_genes)
  in_zone_nodeg <- sum(universe %in% zone_genes & !(universe %in% deg_genes))
  nozone_deg <- sum(!(universe %in% zone_genes) & universe %in% deg_genes)
  nozone_nodeg <- sum(!(universe %in% zone_genes) & !(universe %in% deg_genes))

  mat <- matrix(c(in_zone_deg, in_zone_nodeg, nozone_deg, nozone_nodeg), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")

  data.frame(
    test = label,
    n_zone = sum(universe %in% zone_genes),
    n_deg_in_zone = in_zone_deg,
    odds_ratio = ft$estimate,
    pvalue = ft$p.value,
    ci_lower = ft$conf.int[1],
    stringsAsFactors = FALSE
  )
}

universe <- atlas$human_symbol
deg_all <- atlas$human_symbol[atlas$is_deg]
deg_up <- atlas$human_symbol[atlas$is_deg_up]
deg_down <- atlas$human_symbol[atlas$is_deg_down]

# Hepatocyte-intrinsic
hep_intrinsic <- atlas$human_symbol[!is.na(atlas$attribution_class) &
  atlas$attribution_class == "Hepatocyte_Intrinsic"]

# Conserved
conserved <- atlas$human_symbol[!is.na(atlas$is_conserved) &
  atlas$is_conserved == TRUE]

enrichment <- bind_rows(
  # PC enrichment
  run_fisher("PC_vs_AllDEG", pc_human, deg_all, universe),
  run_fisher("PC_vs_UpDEG", pc_human, deg_up, universe),
  run_fisher("PC_vs_DownDEG", pc_human, deg_down, universe),
  run_fisher("PC_vs_HepIntrinsic", pc_human, hep_intrinsic, universe),
  run_fisher("PC_vs_ConservedCore", pc_human, conserved, universe),

  # PP enrichment
  run_fisher("PP_vs_AllDEG", pp_human, deg_all, universe),
  run_fisher("PP_vs_UpDEG", pp_human, deg_up, universe),
  run_fisher("PP_vs_DownDEG", pp_human, deg_down, universe),
  run_fisher("PP_vs_HepIntrinsic", pp_human, hep_intrinsic, universe),
  run_fisher("PP_vs_ConservedCore", pp_human, conserved, universe)
)

enrichment$padj <- p.adjust(enrichment$pvalue, method = "BH")

cat("\nEnrichment results:\n")
print(enrichment %>%
  select(test, n_zone, n_deg_in_zone, odds_ratio, pvalue, padj) %>%
  arrange(pvalue), row.names = FALSE)

# =============================================================================
# 5. Disease Regulon × Zonation
# =============================================================================
cat("\n--- Disease regulon × zonation ---\n")

dis_regs <- fread(dis_reg_path)
hep_regs <- fread(hep_reg_path)

disease_tfs <- dis_regs$tf_name
disease_targets <- hep_regs %>%
  filter(tf_name %in% disease_tfs) %>%
  pull(target_gene) %>%
  unique()

reg_zon <- data.frame(
  gene = disease_targets,
  stringsAsFactors = FALSE
) %>%
  mutate(
    zonation_class = case_when(
      gene %in% pc_human ~ "Pericentral",
      gene %in% pp_human ~ "Periportal",
      TRUE ~ "Non-zoned"
    )
  )

cat("  Disease regulon targets by zone:\n")
print(table(reg_zon$zonation_class))

# Also classify the TFs themselves
tf_zon <- data.frame(
  tf = disease_tfs,
  zonation_class = case_when(
    disease_tfs %in% pc_human ~ "Pericentral",
    disease_tfs %in% pp_human ~ "Periportal",
    TRUE ~ "Non-zoned"
  )
)
cat("\n  Disease TFs by zone:\n")
print(table(tf_zon$zonation_class))

# =============================================================================
# 6. Cross-Species Zonation Concordance
# =============================================================================
cat("\n--- Cross-species zonation ---\n")

# Check if Conserved genes are zonation-biased
core_genes <- atlas %>% filter(is_conserved == TRUE)
core_pc <- sum(core_genes$is_pc, na.rm = TRUE)
core_pp <- sum(core_genes$is_pp, na.rm = TRUE)
core_non <- nrow(core_genes) - core_pc - core_pp

cat(sprintf("  Conserved zonation: %d PC, %d PP, %d Non-zoned (total %d)\n",
            core_pc, core_pp, core_non, nrow(core_genes)))

# Fisher test: are Conserved genes enriched for PC?
if (core_pc > 0) {
  core_pc_enrichment <- run_fisher("ConservedCore_PC", pc_human,
    core_genes$human_symbol, universe)
  cat(sprintf("  Core × PC enrichment: OR=%.2f, p=%.4f\n",
              core_pc_enrichment$odds_ratio, core_pc_enrichment$pvalue))
}

# =============================================================================
# 7. Disease Progression × Zonation
# =============================================================================
cat("\n--- Disease progression × zonation ---\n")

# Check if NAFL vs NASH signatures differ by zone
nafl_nash_file <- file.path(nafl_nash_path, "nafl_vs_nash_dream.csv")
if (file.exists(nafl_nash_file)) {
  nafl_nash <- fread(nafl_nash_file)
  cat("  NAFL vs NASH DEGs: ", nrow(nafl_nash), "\n")

  # 2026-05-29 (STAR -s2): dream output now uses 'gene' (versioned ENSEMBL), not 'symbol'.
  # Map ENSEMBL -> human_symbol via the atlas (loaded above) so the zone comparison
  # (against human-symbol PC/PP markers) stays valid. Falls back to legacy 'symbol' if present.
  if ("symbol" %in% names(nafl_nash)) {
    nafl_nash <- nafl_nash %>% rename(human_symbol_nn = symbol)
  } else {
    nafl_nash[, .ens := sub("\\.[0-9]+$", "", gene)]
    .amap <- unique(data.table(.ens = sub("\\.[0-9]+$", "", atlas$ensembl_id),
                               human_symbol_nn = atlas$human_symbol))
    nafl_nash <- merge(nafl_nash, .amap, by = ".ens", all.x = TRUE)
  }

  # Rename adj.P.Val to avoid R parsing issues
  if ("adj.P.Val" %in% names(nafl_nash)) {
    setnames(nafl_nash, "adj.P.Val", "padj_nn")
  } else if ("padj" %in% names(nafl_nash)) {
    setnames(nafl_nash, "padj", "padj_nn")
  }

  nafl_nash_deg <- nafl_nash %>%
    filter(!is.na(padj_nn) & padj_nn < 0.1 & !is.na(human_symbol_nn)) %>%
    rename(human_symbol = human_symbol_nn)

  nn_pc <- sum(nafl_nash_deg$human_symbol %in% pc_human)
  nn_pp <- sum(nafl_nash_deg$human_symbol %in% pp_human)

  cat(sprintf("  NAFL→NASH DEGs: %d PC, %d PP (of %d total)\n",
              nn_pc, nn_pp, nrow(nafl_nash_deg)))

  # Direction of change per zone
  nn_pc_up <- sum(nafl_nash_deg$human_symbol %in% pc_human & nafl_nash_deg$logFC > 0)
  nn_pc_down <- sum(nafl_nash_deg$human_symbol %in% pc_human & nafl_nash_deg$logFC < 0)
  nn_pp_up <- sum(nafl_nash_deg$human_symbol %in% pp_human & nafl_nash_deg$logFC > 0)
  nn_pp_down <- sum(nafl_nash_deg$human_symbol %in% pp_human & nafl_nash_deg$logFC < 0)

  cat(sprintf("  PC DEGs: %d up, %d down in NASH vs NAFL\n", nn_pc_up, nn_pc_down))
  cat(sprintf("  PP DEGs: %d up, %d down in NASH vs NAFL\n", nn_pp_up, nn_pp_down))

  progression_zon <- data.frame(
    zone = c("Pericentral", "Pericentral", "Periportal", "Periportal"),
    direction = c("Up_in_NASH", "Down_in_NASH", "Up_in_NASH", "Down_in_NASH"),
    n_degs = c(nn_pc_up, nn_pc_down, nn_pp_up, nn_pp_down)
  )
} else {
  cat("  NAFL vs NASH file not found\n")
  progression_zon <- data.frame()
}

# =============================================================================
# 8. Canonical Marker Validation
# =============================================================================
cat("\n--- Canonical marker validation ---\n")

# Validate that known markers are correctly classified
canonical_pc <- c("CYP2E1", "CYP1A2", "GLUL", "CYP3A4", "FASN", "SCD", "TBX3")
canonical_pp <- c("ALB", "CPS1", "ASS1", "PCK1", "G6PC", "HNF4A", "SLC10A1")

cat("PC markers:\n")
for (g in canonical_pc) {
  ar <- atlas[atlas$human_symbol == g, ]
  if (nrow(ar) > 0) {
    ar <- ar[1, ]
    deg_str <- if (ar$is_deg) sprintf("DEG (LFC=%.2f, p=%.2e)", ar$bulk_logFC, ar$bulk_padj) else "Not DEG"
    cat(sprintf("  %s: %s, zone=%s\n", g, deg_str, ar$zonation_class))
  } else {
    cat(sprintf("  %s: Not in atlas\n", g))
  }
}

cat("\nPP markers:\n")
for (g in canonical_pp) {
  ar <- atlas[atlas$human_symbol == g, ]
  if (nrow(ar) > 0) {
    ar <- ar[1, ]
    deg_str <- if (ar$is_deg) sprintf("DEG (LFC=%.2f, p=%.2e)", ar$bulk_logFC, ar$bulk_padj) else "Not DEG"
    cat(sprintf("  %s: %s, zone=%s\n", g, deg_str, ar$zonation_class))
  } else {
    cat(sprintf("  %s: Not in atlas\n", g))
  }
}

# =============================================================================
# 9. Write Outputs
# =============================================================================
cat("\n--- Writing outputs ---\n")

# Zonation reference
zon_ref <- data.frame(
  human_symbol = c(pc_human, pp_human),
  zonation_class = c(rep("Pericentral", length(pc_human)),
                     rep("Periportal", length(pp_human))),
  source = "Halpern2017_extended",
  stringsAsFactors = FALSE
)
fwrite(zon_ref, file.path(out_dir, "zonation_reference.csv"))
cat("  Written: zonation_reference.csv (", nrow(zon_ref), " genes)\n")

# DEG zonation classification
deg_zon <- atlas %>%
  select(human_symbol, ensembl_id, zonation_class, is_pc, is_pp,
         bulk_logFC, bulk_padj, bulk_tstat, is_deg, is_deg_up, is_deg_down,
         mouse_meta_logFC, is_conserved, attribution_class, sex_class) %>%
  arrange(zonation_class, bulk_padj)

# Merge spatial if available
if ("spatial_zonation_class" %in% names(atlas)) {
  deg_zon <- deg_zon %>%
    left_join(atlas %>% select(human_symbol, spatial_zonation_class), by = "human_symbol")
}

fwrite(deg_zon, file.path(out_dir, "deg_zonation_classification.csv"))
cat("  Written: deg_zonation_classification.csv (", nrow(deg_zon), " genes)\n")

# Enrichment results
fwrite(enrichment, file.path(out_dir, "zonation_enrichment.csv"))
cat("  Written: zonation_enrichment.csv\n")

# Disease progression × zonation
if (nrow(progression_zon) > 0) {
  fwrite(progression_zon, file.path(out_dir, "zonation_progression.csv"))
  cat("  Written: zonation_progression.csv\n")
}

# =============================================================================
# 10. Summary
# =============================================================================
cat("\n--- Final Summary ---\n")

deg_pc <- sum(atlas$is_deg & atlas$is_pc, na.rm = TRUE)
deg_pp <- sum(atlas$is_deg & atlas$is_pp, na.rm = TRUE)
deg_up_pc <- sum(atlas$is_deg_up & atlas$is_pc, na.rm = TRUE)
deg_up_pp <- sum(atlas$is_deg_up & atlas$is_pp, na.rm = TRUE)
deg_down_pc <- sum(atlas$is_deg_down & atlas$is_pc, na.rm = TRUE)
deg_down_pp <- sum(atlas$is_deg_down & atlas$is_pp, na.rm = TRUE)

cat(sprintf("DEGs by zonation: %d PC, %d PP (of %d total DEGs)\n",
            deg_pc, deg_pp, sum(atlas$is_deg, na.rm = TRUE)))
cat(sprintf("  Up DEGs: %d PC, %d PP\n", deg_up_pc, deg_up_pp))
cat(sprintf("  Down DEGs: %d PC, %d PP\n", deg_down_pc, deg_down_pp))

cat("\n=== Module B1 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
