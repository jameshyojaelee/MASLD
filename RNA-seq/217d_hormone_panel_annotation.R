#!/usr/bin/env Rscript
# 217d_hormone_panel_annotation.R — Expanded curated sex-hormone panel annotation
#
# Agent: B5 (Sex-Hormone TF x COLOC) -- HEADLINE A7 (C2 critique fix)
#
# C2 critique: original A2 panel missed hepatic sex-dimorphic genes.
# Add CYP3A4 (F>M), CYP1A2 (M>F), PRLR, SULT2A1.
# Keep all 26 original A2 panel genes.
#
# Inputs:
#   - sex_deg_classification.csv (per-sex LFC + padj + interaction_padj)
#   - gene_level_coloc.csv (PP4)
#   - multi_evidence_atlas.csv (symbol <-> ensembl)
#
# Output:
#   - RNA-seq/results/stratified_causal/sex_hormone_panel_expanded.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Expanded panel
# ---------------------------------------------------------------------------
# Sex steroid receptors + biosynthesis (original A2)
ORIGINAL_A2 <- c(
  # Receptors
  "AR", "ESR1", "ESR2", "PGR", "GPER1",
  # Biosynthesis / metabolism
  "CYP19A1", "SRD5A2", "CYP17A1", "STAR",
  "HSD17B1", "HSD17B2", "HSD17B3", "HSD17B6", "HSD17B7", "HSD17B12",
  "HSD3B1",
  # Transport / binding
  "SHBG", "GHR",
  # Sex-chromosome regulators
  "KDM6A", "DDX3X", "XIST",
  "DDX3Y", "EIF1AY", "RPS4Y1",
  # GH/pulsatile axis TFs
  "BCL6", "STAT5A", "STAT5B", "HNF4A", "CUX2"
)

# C2 additions: canonical hepatic sex-dimorphic genes
C2_ADDITIONS <- c(
  "CYP3A4",   # F>M, canonical female-biased hepatic CYP
  "CYP1A2",   # M>F
  "PRLR",     # prolactin receptor
  "SULT2A1"   # DHEA sulfotransferase
)

PANEL <- unique(c(ORIGINAL_A2, C2_ADDITIONS))
cat("=== 217d: Expanded sex-hormone panel ===\n")
cat("  N panel genes:", length(PANEL), "\n")

# ---------------------------------------------------------------------------
# Load sex DEG classification
# ---------------------------------------------------------------------------
sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_degs <- fread(sex_file)
sex_degs[, ensembl_id := sub("\\..*", "", gene)]

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
sex_degs <- merge(sex_degs, atlas, by = "ensembl_id", all.x = TRUE)

# ---------------------------------------------------------------------------
# Load COLOC
# ---------------------------------------------------------------------------
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
}
coloc <- coloc[gene != "" & !is.na(gene)]
coloc_lookup <- unique(coloc[, .(human_symbol = gene, coloc_best_pp4, coloc_best_gwas)])
coloc_lookup <- coloc_lookup[, .(coloc_best_pp4 = max(coloc_best_pp4, na.rm = TRUE),
                                   coloc_best_gwas = coloc_best_gwas[which.max(coloc_best_pp4)]),
                              by = human_symbol]

# ---------------------------------------------------------------------------
# Build panel table
# ---------------------------------------------------------------------------
panel_dt <- data.table(human_symbol = PANEL)
panel_dt[, in_A2_original := human_symbol %in% ORIGINAL_A2]
panel_dt[, in_C2_addition := human_symbol %in% C2_ADDITIONS]

# Join sex stats
sex_sub <- sex_degs[!is.na(human_symbol) & human_symbol %in% PANEL,
                     .(human_symbol, ensembl_id,
                       logFC_M, logFC_F, padj_M, padj_F,
                       interaction_logFC, interaction_padj,
                       sex_class)]
# Dedup (keep most significant per symbol)
sex_sub <- sex_sub[order(human_symbol, interaction_padj, na.last = TRUE)]
sex_sub <- sex_sub[!duplicated(human_symbol)]
panel_dt <- merge(panel_dt, sex_sub, by = "human_symbol", all.x = TRUE)

# Join COLOC
panel_dt <- merge(panel_dt, coloc_lookup, by = "human_symbol", all.x = TRUE)

# Annotate hormone category
panel_dt[, category := fcase(
  human_symbol %in% c("AR", "ESR1", "ESR2", "PGR", "GPER1", "PRLR"), "receptor",
  human_symbol %in% c("CYP19A1", "SRD5A2", "CYP17A1", "STAR",
                       "HSD17B1", "HSD17B2", "HSD17B3", "HSD17B6", "HSD17B7",
                       "HSD17B12", "HSD3B1"), "biosynthesis",
  human_symbol %in% c("CYP3A4", "CYP1A2", "SULT2A1"), "hepatic_metabolism",
  human_symbol %in% c("SHBG", "GHR"), "transport",
  human_symbol %in% c("KDM6A", "DDX3X", "XIST",
                       "DDX3Y", "EIF1AY", "RPS4Y1"), "sex_chromosome",
  human_symbol %in% c("BCL6", "STAT5A", "STAT5B", "HNF4A", "CUX2"), "TF_GH_axis",
  default = "other"
)]

# Reorder
setcolorder(panel_dt,
            c("human_symbol", "ensembl_id", "category",
              "in_A2_original", "in_C2_addition",
              "logFC_M", "padj_M", "logFC_F", "padj_F",
              "interaction_logFC", "interaction_padj", "sex_class",
              "coloc_best_pp4", "coloc_best_gwas"))
panel_dt <- panel_dt[order(category, human_symbol)]

fwrite(panel_dt, file.path(outdir, "sex_hormone_panel_expanded.csv"))
cat("  Wrote sex_hormone_panel_expanded.csv (", nrow(panel_dt), "rows)\n")
print(panel_dt[, .(human_symbol, category, sex_class, interaction_padj,
                    coloc_best_pp4)])

cat("\n=== 217d complete ===\n")
