#!/usr/bin/env Rscript
# 326_hep_metabolic_macrophage_ccc.R
#
# Analysis E1 (v1) — Hepatocyte metabolic state × macrophage CCC.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load hepatocyte subtype markers + meta-subtype mapping.
#   2. Score each subtype by metabolic modules: FAO (CPT1A,ACOX1,PPARA,HMGCS2),
#      Lipogenic (FASN,SCD,ACACA), ER-stressed (HSPA5,ATF4,DDIT3),
#      Senescent (CDKN2A,CDKN1A,TP53), Lipotoxic (SAA1,SAA2,HP,CRP).
#   3. Identify which Hep subtypes have highest expression of each
#      Hep->Mac LIANA ligand (from D2 output).
#   4. Cross-reference with meta-subtype (Healthy/Disease-Progressor/etc.) to
#      determine if specific metabolic states preferentially drive
#      macrophage-directed CCC.
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/metabolic_ccc/

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HEP  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
CCC  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
OUTDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/metabolic_ccc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading hepatocyte subtype markers + meta-subtype mapping...")
markers <- fread(file.path(HEP, "subtype_markers.csv"))
setnames(markers, "names", "gene")
meta_map <- fread(file.path(HEP, "meta_subtype_mapping.csv"))

# Metabolic modules
modules <- list(
  FAO         = c("CPT1A","CPT1B","CPT2","ACOX1","ACOX2","HMGCS2","ACADVL","ACAT1","PPARA","PPARGC1A","HADHA","HADHB"),
  Lipogenic   = c("FASN","SCD","ACACA","ACACB","SREBF1","INSIG1","MLXIPL","DGAT1","DGAT2"),
  ER_Stress   = c("HSPA5","ATF4","DDIT3","XBP1","EIF2AK3","ERN1","ATF6"),
  Senescent   = c("CDKN2A","CDKN1A","TP53","GDF15","SERPINE1","IGFBP3","LMNB1"),
  Lipotoxic_DAMP = c("SAA1","SAA2","SAA4","HP","CRP","ORM1","ORM2","C3","LBP"),
  Ferroptotic = c("GPX4","SLC7A11","ALOX15","PTGS2","HMOX1","FTH1","FTL","NCOA4")
)

message("[2] Aggregating subtype markers to meta-subtype level...")
markers <- merge(markers, meta_map[, .(subtype, meta_subtype)], by = "subtype", all.x = TRUE)

# For each meta-subtype, average logfoldchanges across its constituent subtypes
meta_markers <- markers[, .(meta_lfc = mean(logfoldchanges, na.rm = TRUE),
                             n_subtypes_expressing = sum(logfoldchanges > 0)),
                         by = .(meta_subtype, gene)]

message("[3] Scoring meta-subtypes by metabolic modules...")
module_scores <- rbindlist(lapply(names(modules), function(mod) {
  module_genes <- modules[[mod]]
  meta_markers[gene %in% module_genes,
               .(module = mod,
                 mean_lfc = mean(meta_lfc, na.rm = TRUE),
                 n_genes = uniqueN(gene)),
               by = meta_subtype]
}), fill = TRUE)
fwrite(module_scores, file.path(OUTDIR, "hep_metabolic_module_scores_by_meta.csv"))

# Pivot to wide
mod_wide <- dcast(module_scores, meta_subtype ~ module, value.var = "mean_lfc")
fwrite(mod_wide, file.path(OUTDIR, "hep_metabolic_module_scores_wide.csv"))

message("[4] Loading LIANA Hep->Mac ligands (from D2)...")
hep_mac <- fread(file.path(CCC, "D2_hep_Mac_bidirectional.csv"))
hep_to_mac <- hep_mac[direction == "Hep_to_Mac" & score_diff > 0.1]
hep_to_mac_ligands <- unique(hep_to_mac$ligand)
message(sprintf("  MASLD-enriched Hep->Mac ligands: %d", length(hep_to_mac_ligands)))

message("[5] For each meta-subtype, which Hep->Mac ligands is it enriched for?")
lig_by_meta <- meta_markers[gene %in% hep_to_mac_ligands &
                            meta_lfc > 0.1,
                            .(meta_subtype, gene, meta_lfc, n_subtypes_expressing)][
                            order(meta_subtype, -meta_lfc)]
fwrite(lig_by_meta, file.path(OUTDIR, "hep_to_mac_ligands_by_meta.csv"))

# Summary per meta-subtype
meta_summary <- lig_by_meta[, .(n_ligands_enriched = uniqueN(gene),
                                 top_3_ligands = paste(gene[1:min(3, .N)], collapse = ",")),
                             by = meta_subtype]

message("[6] Identify metabolic-specific CCC: for each Hep->Mac ligand, which module it belongs to?")
# Is each ligand a metabolic module gene?
lig_module_map <- rbindlist(lapply(hep_to_mac_ligands, function(g) {
  mods_containing <- names(modules)[sapply(modules, function(m) g %in% m)]
  if (length(mods_containing) == 0) return(NULL)
  data.table(gene = g, module = mods_containing)
}), fill = TRUE)

# Cross: LIANA score_diff for each metabolic-module ligand
hep_to_mac_metabolic <- merge(hep_to_mac, lig_module_map,
                              by.x = "ligand", by.y = "gene", all.y = TRUE)
fwrite(hep_to_mac_metabolic, file.path(OUTDIR, "hep_to_mac_metabolic_ligand_axes.csv"))

summary_lines <- c(
  sprintf("Hep subtypes: %d, meta-subtypes: %d",
          length(unique(meta_map$subtype)), length(unique(meta_map$meta_subtype))),
  "",
  "=== Hep metabolic module scores by meta-subtype (mean_lfc) ===",
  capture.output(print(mod_wide)),
  "",
  sprintf("MASLD-enriched Hep->Mac ligands: %d", length(hep_to_mac_ligands)),
  "",
  "=== Per meta-subtype: Hep->Mac ligand enrichment ===",
  capture.output(print(meta_summary)),
  "",
  "=== Metabolic module ligands in Hep->Mac LIANA CCC ===",
  capture.output(print(hep_to_mac_metabolic[!is.na(module),
                        .(ligand, receptor, module, score_diff, both_concordant)][order(-score_diff)],
                       nrows = 40))
)
writeLines(summary_lines, file.path(OUTDIR, "hep_metabolic_macrophage_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
