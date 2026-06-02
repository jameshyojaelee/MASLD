#!/usr/bin/env Rscript
# ============================================================================
# 77_hpa_annotation.R
# Tier-1 atlas layer: Human Protein Atlas target characterization.
#
# Adds tissue-validated protein context (classification/ranking layer; IHC is
# semi-quantitative, so treated as annotation, not quantitative evidence):
#   - hpa_rna_tissue_specificity   (Tissue enriched / enhanced / Low specificity ...)
#   - hpa_liver_elevated           (bool: liver among elevated tissues)
#   - hpa_secretome_location       (e.g., "Secreted to blood" -> plasma-biomarker tractable)
#   - hpa_subcellular_main         (drug-accessibility context)
#   - hpa_protein_class            (e.g., "Plasma proteins, Enzymes, Transporters, ...")
#   - hpa_liver_hcc_prognostic     (TCGA Liver HCC prognostic direction)
#
# Output: RNA-seq/results/multi_evidence/hpa_atlas_columns.tsv
# Merged by 27a's post-assembly block (drop-then-merge on human_symbol).
# ============================================================================

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HPA  <- file.path(BASE, "data/external/hpa/proteinatlas.tsv")
OUT  <- file.path(BASE, "RNA-seq/results/multi_evidence/hpa_atlas_columns.tsv")

cat("== 77_hpa_annotation ==\n")

# NOTE: the TCGA Liver-HCC prognostic column is deliberately NOT ingested — HCC is
# a cancer endpoint, not MASLD, and would be indefensible as MASLD evidence (v2
# relevance audit 2026-06-01). Only disease-agnostic functional characterization
# (T3 annotation tier) is kept: tissue specificity, secretome, subcellular, class.
cols <- c("Gene", "Ensembl", "Protein class", "RNA tissue specificity",
          "RNA tissue specific nTPM", "Secretome location",
          "Subcellular main location")
hpa <- fread(HPA, select = cols, quote = "\"")

out <- hpa[, .(
  human_symbol            = Gene,
  hpa_rna_tissue_specificity = `RNA tissue specificity`,
  # liver among the tissue-specific/elevated nTPM list
  hpa_liver_elevated      = grepl("Liver", `RNA tissue specific nTPM`, ignore.case = TRUE),
  hpa_secretome_location  = `Secretome location`,
  hpa_subcellular_main    = `Subcellular main location`,
  hpa_protein_class       = `Protein class`
)]
out <- out[human_symbol != "" & !is.na(human_symbol)]
out <- out[!duplicated(human_symbol)]
# blank -> NA for cleaner downstream handling
for (c in setdiff(names(out), c("human_symbol", "hpa_liver_elevated")))
  out[get(c) == "", (c) := NA_character_]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT, sep = "\t")
cat(sprintf("WROTE %s : %d genes x %d cols\n", OUT, nrow(out), ncol(out)))
cat(sprintf("  liver-elevated: %d  secreted-to-blood: %d\n",
            sum(out$hpa_liver_elevated, na.rm = TRUE),
            sum(grepl("blood", out$hpa_secretome_location, ignore.case = TRUE), na.rm = TRUE)))
for (g in c("THRB", "PNPLA3", "ALB", "SERPINE1")) {
  r <- out[human_symbol == g]
  if (nrow(r)) cat(sprintf("  [%s] spec=%s liver_elev=%s secretome=%s subcell=%s\n",
                           g, r$hpa_rna_tissue_specificity[1], r$hpa_liver_elevated[1],
                           substr(r$hpa_secretome_location[1], 1, 25),
                           substr(r$hpa_subcellular_main[1], 1, 25)))
}
