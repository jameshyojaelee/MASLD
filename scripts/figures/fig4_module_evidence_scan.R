#!/usr/bin/env Rscript
# ============================================================================
# fig4_module_evidence_scan.R
#
# Computes a 7-feature cross-modal evidence matrix for all 30 hepatocyte
# Hotspot modules. Used for Panel 4a (evidence heatmap) and for determining
# which modules are highlighted in the manuscript text.
#
# Features (all computed from existing result files; no heavy compute):
#   F1 — Severity association:  -log10(F_stage_q)
#   F2 — Bulk replication score: |bulk_module_score|
#   F3 — LOO stability:          stability_score
#   F4 — Bulk replicated (binary)
#   F5 — COLOC fraction:         fraction of module genes with susie PP4 > 0.5
#   F6 — SVG fraction:           fraction of module genes that are disease-emergent SVGs
#   F7 — Plasma protein count:   number of module genes detected in PXD052937
#
# Tier A = ≥3 features above threshold (main figure candidates)
# Tier B = 2 features (supplementary)
# Tier C = ≤1 feature (context only)
#
# Output: figures/main/fig4_validation/panels/data/module_evidence_matrix.csv
# Runtime: ~1 min on login node
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

HS_RES   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
OUT_DIR  <- file.path(FIG4_DIR, "panels", "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Module-level features ────────────────────────────────────────────────
mods <- fread(file.path(HS_RES, "all_modules.tsv"))
hep  <- mods[cell_type == "hepatocytes",
             .(module     = as.integer(module),
               F_stage_q  = suppressWarnings(as.numeric(F_stage_q)),
               bulk_score = suppressWarnings(as.numeric(bulk_module_score)),
               stability  = suppressWarnings(as.numeric(stability_score)),
               bulk_rep   = as.integer(bulk_replicated %in% c("TRUE", TRUE, 1L)))]

# ── 2. Module gene lists ─────────────────────────────────────────────────────
mod_genes <- fread(file.path(HS_RES, "hepatocytes", "module_genes.tsv"))
setnames(mod_genes, c("gene", "module", "weight"))
mod_genes[, module := as.integer(module)]

# ── 3. Atlas join for COLOC ──────────────────────────────────────────────────
atlas_cols <- c("human_symbol", "coloc_susie_best_pp4")
atlas_f    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas      <- fread(atlas_f, select = atlas_cols)
setnames(atlas, "human_symbol", "gene")
atlas[, coloc_susie_best_pp4 := suppressWarnings(as.numeric(coloc_susie_best_pp4))]

mg_atlas <- atlas[mod_genes, on = "gene"]
coloc_feat <- mg_atlas[, .(
  coloc_frac = mean(!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5, na.rm = TRUE),
  coloc_n    = sum(!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5, na.rm = TRUE),
  n_genes    = .N
), by = module]

# ── 4. SVG feature: fraction of module genes that are disease-emergent ───────
svg_f <- file.path(SPATIAL_DIR, "svg", "differential_svgs.csv")
if (file.exists(svg_f)) {
  diff_svg <- fread(svg_f)
  setnames(diff_svg, 1, "gene")
  emergent_genes <- diff_svg[grepl("emergent", category, ignore.case = TRUE), gene]
} else {
  # Fallback: known disease-emergent SVGs
  emergent_genes <- c("JUN","ATF3","KLF6","SERPINE1","HSPA1A","JUND","LDLR","COL1A1",
                      "COL3A1","TIMP1","SPP1","ACTA2","LUM","PDGFRB","THBS2")
}
svg_feat <- mod_genes[, .(
  svg_frac = mean(gene %in% emergent_genes, na.rm = TRUE),
  svg_n    = sum(gene %in% emergent_genes, na.rm = TRUE)
), by = module]

# ── 5. Plasma feature: PXD052937 detected protein genes ─────────────────────
# First try the pre-mapped HGNC file written by fig4_validation.R
conc_hgnc_f <- file.path(PROTEOMICS_DIR, "pxd052937_mrna_protein_concordance.csv")
conc_v3_f   <- file.path(PROTEOMICS_DIR, "protein_transcript_concordance_v3.csv")

if (file.exists(conc_hgnc_f)) {
  plasma_dt    <- fread(conc_hgnc_f)
  plasma_genes <- unique(plasma_dt[!is.na(gene), gene])
} else if (file.exists(conc_v3_f)) {
  # Use raw concordance; 'gene' col may be UniProt in PXD052937 rows
  prot         <- fread(conc_v3_f)
  plasma_genes <- unique(prot[dataset == "PXD052937" & !is.na(gene), gene])
} else {
  # Hard fallback: known plasma proteins
  plasma_genes <- c("CHI3L1","LGALS3","TIMP1","A2M","SERPINE1","HGF","APOB",
                    "APOA1","THBS2","IGFBP7","IGFBP1","SOD2","LEPR","GDF15",
                    "GLS","BICC1","HKDC1","TXNRD1","SQSTM1","AKR1B10")
}
plasma_feat <- mod_genes[, .(
  plasma_n = sum(gene %in% plasma_genes, na.rm = TRUE)
), by = module]

# ── 6. Assemble feature table ─────────────────────────────────────────────────
feat <- hep[coloc_feat[, .(module, coloc_frac, coloc_n, n_genes)], on = "module"]
feat <- feat[svg_feat, on = "module"]
feat <- feat[plasma_feat, on = "module"]

# ── 7. Normalize 0-1 within the 30-module set ─────────────────────────────────
norm01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(ifelse(is.na(x), NA_real_, 0.5))
  (x - r[1]) / diff(r)
}
feat[, fstage_z  := norm01(-log10(pmax(F_stage_q, 1e-25, na.rm = TRUE)))]
feat[, bulk_z    := norm01(abs(bulk_score))]
feat[, stab_z    := norm01(stability)]
feat[, brep_z    := as.numeric(bulk_rep)]
feat[, coloc_z   := norm01(coloc_frac)]
feat[, svg_z     := norm01(svg_frac)]
feat[, plasma_z  := norm01(plasma_n)]

# ── 8. Per-feature pass thresholds ────────────────────────────────────────────
# Severity is MANDATORY (must pass for Tier A or B).
# Validation passes are counted on top: need ≥2 independent layers for Tier A.
feat[, fstage_pass := !is.na(F_stage_q) & F_stage_q < 0.01]
feat[, bulk_pass   := !is.na(bulk_score) & abs(bulk_score) >
                        quantile(abs(na.omit(feat$bulk_score)), 0.6)]  # top 40%
feat[, stab_pass   := !is.na(stability) & stability > 0.55]
feat[, brep_pass   := bulk_rep == 1L]
feat[, coloc_pass  := !is.na(coloc_frac) & coloc_frac >= 0.05]  # ≥5% genes PP4>0.5
feat[, svg_pass    := !is.na(svg_frac)   & svg_frac   >= 0.05]  # ≥5% disease-emergent
feat[, plasma_pass := !is.na(plasma_n)   & plasma_n   >= 3L]    # ≥3 genes detected

# Validation-layer n_pass excludes fstage (mandatory) and brep (low info)
val_cols <- c("bulk_pass","stab_pass","coloc_pass","svg_pass","plasma_pass")
feat[, n_val_pass := rowSums(.SD == TRUE, na.rm = TRUE), .SDcols = val_cols]

# All-in pass for reporting
all_pass_cols <- c("fstage_pass","bulk_pass","stab_pass","brep_pass",
                   "coloc_pass","svg_pass","plasma_pass")
feat[, n_pass := rowSums(.SD == TRUE, na.rm = TRUE), .SDcols = all_pass_cols]

# ── 9. Composite score and tier assignment ────────────────────────────────────
feat[, composite := rowMeans(cbind(fstage_z, bulk_z, stab_z, brep_z,
                                   coloc_z, svg_z, plasma_z), na.rm = TRUE)]
# Tier A: severity passes + ≥2 independent validation layers
# Tier B: severity passes + 1 validation layer
# Tier C: everything else
feat[, tier := fcase(
  fstage_pass == TRUE & n_val_pass >= 2L, "A_strong",
  fstage_pass == TRUE & n_val_pass == 1L, "B_moderate",
  default                                = "C_weak"
)]

setorder(feat, -composite)

fwrite(feat, file.path(OUT_DIR, "module_evidence_matrix.csv"))
message(sprintf(
  "Wrote module_evidence_matrix.csv | %d modules | Tier A: %d | Tier B: %d",
  nrow(feat),
  sum(feat$tier == "A_strong"),
  sum(feat$tier == "B_moderate")
))
message("Top 5 by composite score:")
print(feat[1:5, .(module, composite, n_pass, tier, F_stage_q, bulk_rep, coloc_frac, svg_frac, plasma_n)])
