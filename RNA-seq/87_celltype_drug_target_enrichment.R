#!/usr/bin/env Rscript
# 87_celltype_drug_target_enrichment.R
#
# Analysis F3 (v1) — Cell-type-specific drug target enrichment.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load drug target gene lists:
#      - OpenTargets known drugs (RNA-seq/results/drug_repurposing/)
#      - LINCS CGP reversal compounds
#      - Curated MASLD drug targets (THRB, NR1H4/FXR, GLP1R, PPARA/G/D)
#   2. Load A1 primary cell-type attribution.
#   3. For each drug (or drug class), enrich its targets in cell-type-specific
#      attribution: which cell types does drug X primarily act on?
#   4. Overlay with GWAS / LINCS reversal to prioritize new targets.
#
# Env: rnaseq
# Outputs: RNA-seq/results/celltype_drug/

suppressPackageStartupMessages({
  library(data.table)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DRUG    <- file.path(BASE, "RNA-seq/results/drug_repurposing")
ATT     <- file.path(BASE, "RNA-seq/results/celltype_attribution")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR  <- file.path(BASE, "RNA-seq/results/celltype_drug")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading drug target data...")
ot <- fread(file.path(DRUG, "opentargets_known_drugs.csv"))
message(sprintf("  OpenTargets rows: %d", nrow(ot)))
# Standardize target symbol column
sym_col <- intersect(c("target_symbol","symbol","gene_symbol","gene"), names(ot))[1]
drug_col <- intersect(c("drug_name","drug","compound_name","drug_id"), names(ot))[1]
moa_col  <- intersect(c("mechanism_of_action","moa","action_type"), names(ot))[1]
setnames(ot, c(sym_col, drug_col), c("target_symbol","drug_name"))

# Curated MASLD therapeutic targets
curated_targets <- data.table(
  target_symbol = c("THRB","NR1H4","GLP1R","PPARA","PPARG","PPARD","FGF21","FGFR1",
                    "FGFR4","DGAT2","SCD","ACC1","ACC2","ACACA","ACACB","FASN",
                    "CPT1A","CPT1B","HMGCR","SREBF1","SREBF2","NR5A2","NR1I2","NR1I3",
                    "RORA","CEBPB","HNF4A","LPL","APOC3","APOB","MTTP","PCSK9"),
  drug_class = c("Thyroid receptor beta","FXR","GLP-1","PPARalpha","PPARgamma","PPARdelta",
                 "FGF21","FGFR1","FGFR4","DGAT2","SCD1","ACC","ACC","ACC","ACC","FASN",
                 "CPT1A","CPT1B","HMGCR (statin)","SREBP1","SREBP2","LRH-1","PXR","CAR",
                 "RORalpha","CEBPB","HNF4A","LPL","APOC3","APOB","MTTP","PCSK9")
)

message("[2] Loading A1 primary attribution...")
att <- fread(file.path(ATT, "celltype_primary_attribution.csv"))

message("[3] Joining drug targets with cell-type attribution (plus bulk LFC fallback)...")
# Load canonical bulk DE for fallback: if target isn't NS-labelled but not in A1
# primary, we can still use bulk logFC direction.
bulk <- fread(file.path(INT_RES, "canonical_deg_results.csv"),
              select = c("symbol","logFC","padj"))
bulk <- bulk[!is.na(symbol) & symbol != ""]
setnames(bulk, c("logFC","padj"), c("bulk_lfc_full","bulk_padj_full"))
bulk <- unique(bulk[order(-abs(bulk_lfc_full))][!duplicated(symbol)])

# OpenTargets drugs
ot_sub <- unique(ot[, .(target_symbol, drug_name,
                         moa = if (!is.na(moa_col)) get(moa_col) else NA_character_)])
ot_joined <- merge(ot_sub, att, by.x = "target_symbol", by.y = "symbol", all.x = TRUE)
ot_joined <- merge(ot_joined, bulk, by.x = "target_symbol", by.y = "symbol", all.x = TRUE)
fwrite(ot_joined, file.path(OUTDIR, "opentargets_drugs_celltype_annotated.csv"))

# Curated
cur_joined <- merge(curated_targets, att, by.x = "target_symbol", by.y = "symbol", all.x = TRUE)
cur_joined <- merge(cur_joined, bulk, by.x = "target_symbol", by.y = "symbol", all.x = TRUE)
fwrite(cur_joined, file.path(OUTDIR, "curated_masld_targets_celltype.csv"))

message("[4] Cell-type enrichment per drug class (OpenTargets)...")
# For each drug, count primary_celltype frequency
drug_ct <- ot_joined[!is.na(primary_celltype),
                     .N, by = .(drug_name, primary_celltype)]
# Top cell type per drug
drug_top_ct <- drug_ct[order(-N)][!duplicated(drug_name)]
fwrite(drug_top_ct, file.path(OUTDIR, "drug_top_celltype.csv"))

# Aggregate per cell type: total drug-target hits + DEG count
ct_drug_summary <- drug_ct[, .(n_drug_targets = sum(N),
                                n_distinct_drugs = uniqueN(drug_name)),
                           by = primary_celltype][order(-n_drug_targets)]

# LINCS top50 reversal compounds
if (file.exists(file.path(DRUG, "top50_reversal_compounds.csv"))) {
  lincs <- fread(file.path(DRUG, "top50_reversal_compounds.csv"))
  message(sprintf("  LINCS top 50 columns: %s", paste(names(lincs), collapse = ", ")))
  fwrite(lincs, file.path(OUTDIR, "lincs_top50_reference.csv"))
}

message("[5] Summary...")
summary_lines <- c(
  sprintf("OpenTargets drugs: %d unique drug names, %d gene targets",
          uniqueN(ot_joined$drug_name), uniqueN(ot_joined$target_symbol)),
  sprintf("Curated MASLD targets: %d", nrow(curated_targets)),
  "",
  "=== Curated MASLD target celltype attribution ===",
  capture.output(print(cur_joined[!is.na(primary_celltype),
                                   .(target_symbol, drug_class, primary_celltype,
                                     attribution_class, bulk_lfc, bulk_padj)],
                       nrows = 40)),
  "",
  "=== Top 20 drugs by primary celltype ===",
  capture.output(print(drug_top_ct[1:20], nrows = 20)),
  "",
  "=== Cell-type-level drug target load ===",
  capture.output(print(ct_drug_summary, nrows = 20)))
writeLines(summary_lines, file.path(OUTDIR, "celltype_drug_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
