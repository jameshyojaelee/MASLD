#!/usr/bin/env Rscript
# _build_reference_caches.R
# Build cached gene-set tables for Hotspot novelty comparison.
#
# Output files: data/reference_signatures/{hallmark_v2025_1,scenic_hep_regulons,liver_curated}.tsv
# Each TSV has 3 cols: panel, program, gene
#
# Sources:
#   hallmark     — msigdbr v10, Hallmark 2025.1 (50 gene sets, ~7,500 rows)
#   scenic_hep   — Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
#                  (139 rows; tf_name = TF/regulon, target_gene = gene)
#   curated_liver:
#     govaere25  — 25-gene fibrosis progression panel (Govaere 2020, Sci Transl Med)
#     govaere_inflammatory — inflammatory subtype markers (Govaere 2023, Nat Metab)
#     tzouanas_LI / _LD / _SU / _SD — 4 hepatocyte stress programs (Tzouanas 2026, Cell)
#     feng       — 27-gene core progression signature (Feng 2026, bioRxiv)
#
# Usage: micromamba run -n rnaseq Rscript _build_reference_caches.R
# Login node OK — lightweight R only.

suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(tibble)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(ROOT, "data/reference_signatures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Hallmark v2025.1 -------------------------------------------------------
# msigdbr v10 uses `collection`/`subcollection` (not category/subcategory)
message("\n=== Building Hallmark ===")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |>
  transmute(panel   = "hallmark",
            program = gs_name,
            gene    = gene_symbol) |>
  distinct()
write_tsv(hallmark, file.path(OUTDIR, "hallmark_v2025_1.tsv"))
message(sprintf("Hallmark: %d rows, %d programs",
                nrow(hallmark), n_distinct(hallmark$program)))

# ---- SCENIC+ hepatocyte regulons --------------------------------------------
# Source: Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
# Schema: tf_name, target_gene, regulon_id, n_target_genes, target_genes,
#         n_enhancers, mean_activity_masld, mean_activity_normal, ...
message("\n=== Building SCENIC+ hepatocyte regulons ===")
SCENIC_REGULON_PATH <- file.path(
  ROOT, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv"
)
scenic_written <- FALSE
if (file.exists(SCENIC_REGULON_PATH)) {
  scenic_raw <- read_csv(SCENIC_REGULON_PATH, show_col_types = FALSE)
  # Use tf_name as program, target_gene as gene (both are present in this file)
  tf_col   <- intersect(c("tf_name", "TF", "tf", "regulator", "source"),
                        colnames(scenic_raw))[1]
  gene_col <- intersect(c("target_gene", "target", "gene", "Gene"),
                        colnames(scenic_raw))[1]
  if (!is.na(tf_col) && !is.na(gene_col)) {
    scenic <- scenic_raw |>
      # Drop rows where target gene is an Ensembl ID (non-symbol rows)
      filter(!grepl("^ENSG", .data[[gene_col]])) |>
      transmute(panel   = "scenic_hep",
                program = .data[[tf_col]],
                gene    = .data[[gene_col]]) |>
      distinct()
    write_tsv(scenic, file.path(OUTDIR, "scenic_hep_regulons.tsv"))
    message(sprintf("SCENIC+: %d rows, %d regulons",
                    nrow(scenic), n_distinct(scenic$program)))
    scenic_written <- TRUE
  } else {
    warning(sprintf(
      "SCENIC+ file found but expected columns not detected. Found: %s",
      paste(colnames(scenic_raw), collapse = ", ")
    ))
  }
}
if (!scenic_written) {
  warning("SCENIC+ regulon file not found or schema not recognized; writing empty placeholder.")
  write_tsv(
    tibble(panel = character(), program = character(), gene = character()),
    file.path(OUTDIR, "scenic_hep_regulons.tsv")
  )
}

# ---- Curated liver signatures -----------------------------------------------
# Sources: all gene lists are hardcoded from published papers (no file dependency).
# Govaere 25-gene panel from Sci Transl Med 2020 (PMID 32546671)
govaere_25 <- c(
  "AKR1B10", "DUSP6", "GDF15", "THBS2", "A2M", "CDH2",
  "COL1A1", "COL3A1", "COL4A1", "COL4A2", "COL6A3", "DCN",
  "FBN1", "FSTL1", "IGFBP7", "LUM", "MFAP4", "MMP2",
  "POSTN", "SPARC", "SPP1", "TAGLN", "THY1", "TIMP1", "VCAN"
)

# Govaere inflammatory subtype markers (Nat Metab 2023) — from 312a benchmark
govaere_inflammatory <- c(
  "ADAMTSL2", "AKR1B10", "CFHR4", "TREM2", "COL1A1", "LUM",
  "VCAN", "THBS2", "CCL2", "CCL20", "CXCL8", "IL1B", "SPP1"
)

# Tzouanas et al. (Cell 2026) — 4 hepatocyte stress programs from 312a benchmark
tzouanas_LI <- c(
  "TNF", "IL6", "CXCL8", "CCL2", "COL1A1", "TGFB1", "ACTA2", "TIMP1",
  "IL1B", "CXCL1", "CCL20", "SERPINE1", "MMP9", "ICAM1", "VCAM1"
)
tzouanas_LD <- c(
  "HMGCS2", "CYP2E1", "CYP3A4", "ALB", "APOB", "PCK1", "G6PC",
  "ALDOB", "HAL", "ASS1", "TAT", "HPD", "SDS", "AGXT", "OTC",
  "ARG1", "CPS1", "ASGR1", "HGD", "ABCB11"
)
tzouanas_SU <- c(
  "LGALS3", "SPP1", "CD9", "TREM2", "GPNMB", "FABP5", "CTSB",
  "CTSD", "LAMP1", "LIPA", "NPC2", "GRN", "APOE", "LPL"
)
tzouanas_SD <- c(
  "GLUL", "AXIN2", "WNT2", "RSPO3", "LGR5", "TBX3", "RNF43",
  "ODAM", "CYP1A2", "CYP2A6"
)

# Feng et al. (bioRxiv 2026) — 27-gene core progression signature
feng_top27 <- c(
  "EFHD1", "MLIP", "TREM2", "SPP1", "GPNMB", "CCL2",
  "CCL20", "CXCL1", "CXCL6", "IL1B", "IL32",
  "COL1A1", "COL1A2", "COL3A1", "FN1", "LOX", "LOXL2",
  "ACTA2", "PDGFRB", "TGFB1", "SERPINE1", "MMP9",
  "AKR1B10", "GDF15", "THY1", "THBS2", "LUM"
)

message("\n=== Building curated liver signatures ===")
curated <- bind_rows(
  tibble(panel = "curated_liver", program = "govaere25",
         gene = govaere_25),
  tibble(panel = "curated_liver", program = "govaere_inflammatory",
         gene = govaere_inflammatory),
  tibble(panel = "curated_liver", program = "tzouanas_LI",
         gene = tzouanas_LI),
  tibble(panel = "curated_liver", program = "tzouanas_LD",
         gene = tzouanas_LD),
  tibble(panel = "curated_liver", program = "tzouanas_SU",
         gene = tzouanas_SU),
  tibble(panel = "curated_liver", program = "tzouanas_SD",
         gene = tzouanas_SD),
  tibble(panel = "curated_liver", program = "feng_top27",
         gene = feng_top27)
) |> distinct()

write_tsv(curated, file.path(OUTDIR, "liver_curated.tsv"))
message(sprintf("Curated: %d rows, %d programs",
                nrow(curated), n_distinct(curated$program)))

message("\n=== Done ===")
message(sprintf("Outputs written to: %s", OUTDIR))
message(sprintf("  hallmark_v2025_1.tsv      (%d rows)", nrow(hallmark)))
message(sprintf("  scenic_hep_regulons.tsv   (%s)", if (scenic_written) sprintf("%d rows", nrow(scenic)) else "EMPTY PLACEHOLDER"))
message(sprintf("  liver_curated.tsv         (%d rows)", nrow(curated)))
