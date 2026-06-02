#!/usr/bin/env Rscript
# compute_translational_priority.R
#
# Fig 2 panel C3/C5 intermediate: cross-species x variant-disruption intersection.
#
# Inputs:
#   - Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/conserved_masld_ccc_priority.csv
#       (LIANA differential L-R pairs annotated with cross-species concordance tier;
#        score_diff = MASLD - control LIANA score on the human side; mouse signal is
#        embedded in axis_conservation via the dvc_category atlas, NOT a separate
#        per-axis score.)
#   - Analysis/SingleCell/results_gpu_v2/chromatin_ccc/priority_variant_disrupted_CCC_axes.csv
#       (variant-disrupted ligand/receptor flags from the GWAS-ATAC pipeline.)
#   - GWAS/finemapping/results/susie_coloc_1kg/gene_level_coloc_1kg.csv
#       (canonical 1KG-LD SuSiE-COLOC gene-level table; used to attach a gwas_locus
#        tag = the gene symbol itself when PP.H4 > 0.5 in any GWAS.)
#
# Algorithm:
#   1. Inner join crossspecies x variant-disrupted on
#        (source, target, ligand_complex, receptor_complex). Heteromer receptors
#        like B2M_FCGRT are kept as a single string.
#   2. Classify priority_class:
#        top_priority : Fully_Conserved AND (ligand or receptor motif-disrupted)
#                       AND score_diff (human) > 0.5
#        strong       : Fully_Conserved AND (ligand or receptor ATAC-regulated)
#        moderate     : Partially_Conserved AND any disruption
#        weak         : everything else
#   3. Attach gwas_locus from gene-level COLOC PP.H4 > 0.5 (ligand or receptor).
#
# Outputs:
#   - Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/translational_priority_axes.tsv
#   - Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/translational_priority_summary.tsv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

XSP_PRIORITY_F <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/conserved_masld_ccc_priority.csv")
# Full annotated LIANA table (all conservation tiers, not pre-filtered to
# Fully_Conserved); needed so that Partially_Conserved -> "moderate" and other
# tiers -> "weak" classes are reachable in the priority_class output.
XSP_FULL_F <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/liana_crossspecies_annotated.csv")
VAR_F <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/chromatin_ccc/priority_variant_disrupted_CCC_axes.csv")
COLOC_F <- file.path(BASE,
  "GWAS/finemapping/results/susie_coloc_1kg/gene_level_coloc_1kg.csv")

OUTDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/crossspecies_ccc")
OUT_AXES    <- file.path(OUTDIR, "translational_priority_axes.tsv")
OUT_SUMMARY <- file.path(OUTDIR, "translational_priority_summary.tsv")

# ---------- 1. Load --------------------------------------------------------
message("[1] Loading inputs...")
xsp <- fread(XSP_FULL_F)        # full annotated LIANA (all conservation tiers)
xsp_pri <- fread(XSP_PRIORITY_F) # Fully_Conserved priority subset (sanity check)
var <- fread(VAR_F)
message(sprintf("  crossspecies rows (full annotated): %d", nrow(xsp)))
message(sprintf("  crossspecies rows (Fully_Conserved priority subset): %d",
                nrow(xsp_pri)))
message(sprintf("  variant-disrupted rows: %d", nrow(var)))

# Canonical join keys: source, target, ligand_complex, receptor_complex.
# All four columns are character, heteromers preserved as single string (B2M_FCGRT).
join_keys <- c("source", "target", "ligand_complex", "receptor_complex")

# Some xsp columns we need; also keep score_diff and conservation tier.
xsp_keep <- xsp[, .(source, target, ligand_complex, receptor_complex,
                    ligand_first, receptor_first,
                    score_diff_human = score_diff,
                    ligand_h_lfc, receptor_h_lfc,
                    ligand_category, receptor_category,
                    axis_conservation)]

var_keep <- var[, .(source, target, ligand_complex, receptor_complex,
                    var_score_diff = score_diff,
                    ligand_atac_regulated, receptor_atac_regulated,
                    ligand_motif_disrupted, receptor_motif_disrupted)]

# ---------- 2. Inner join --------------------------------------------------
message("[2] Inner join on (source, target, ligand_complex, receptor_complex)...")
joined <- merge(xsp_keep, var_keep, by = join_keys)
message(sprintf("  joined rows: %d", nrow(joined)))

# Coerce TRUE/FALSE columns to logical (fread reads as "TRUE"/"FALSE" strings or
# logical depending on the file; force logical).
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  out <- toupper(as.character(x))
  out <- out == "TRUE"
  out
}
joined[, ligand_atac_regulated   := to_logical(ligand_atac_regulated)]
joined[, receptor_atac_regulated := to_logical(receptor_atac_regulated)]
joined[, ligand_motif_disrupted  := to_logical(ligand_motif_disrupted)]
joined[, receptor_motif_disrupted := to_logical(receptor_motif_disrupted)]

# ---------- 3. GWAS locus annotation --------------------------------------
message("[3] Attaching GWAS locus tags from gene-level COLOC...")
gwas_locus <- character(nrow(joined))
gwas_locus[] <- NA_character_

if (file.exists(COLOC_F)) {
  cl <- fread(COLOC_F, select = c("gene", "coloc_best_pp4", "coloc_best_susie_pp4",
                                  "coloc_best_gwas", "coloc_best_susie_gwas"))
  # PP.H4 > 0.5 in either ABF or SuSiE branch
  cl[, pp4_max := pmax(coloc_best_pp4, coloc_best_susie_pp4, na.rm = TRUE)]
  cl <- cl[!is.na(pp4_max) & pp4_max > 0.5]
  cl[, gwas := fifelse(!is.na(coloc_best_susie_pp4) &
                         coloc_best_susie_pp4 >= replace(coloc_best_pp4, is.na(coloc_best_pp4), -1),
                       coloc_best_susie_gwas, coloc_best_gwas)]
  coloc_genes <- setNames(cl$gwas, cl$gene)

  # For each row, check if ligand_first OR receptor_first (split heteromers) hit COLOC.
  split_unit <- function(s) unique(unlist(strsplit(s, "_", fixed = TRUE)))

  joined[, gwas_locus := vapply(seq_len(.N), function(i) {
    units <- unique(c(split_unit(ligand_complex[i]), split_unit(receptor_complex[i])))
    hits <- units[units %in% names(coloc_genes)]
    if (length(hits) == 0L) return(NA_character_)
    paste0(hits, "(", coloc_genes[hits], ")", collapse = ";")
  }, character(1))]
  message(sprintf("  COLOC hits attached: %d / %d rows",
                  sum(!is.na(joined$gwas_locus)), nrow(joined)))
} else {
  joined[, gwas_locus := NA_character_]
  message("  COLOC table missing; gwas_locus = NA for all rows.")
}

# ---------- 4. Priority class ---------------------------------------------
message("[4] Computing priority_class...")
joined[, any_motif := ligand_motif_disrupted | receptor_motif_disrupted]
joined[, any_atac  := ligand_atac_regulated  | receptor_atac_regulated]
joined[, any_disruption := any_motif | any_atac]

joined[, priority_class := fifelse(
  axis_conservation == "Fully_Conserved" & any_motif & score_diff_human > 0.5,
  "top_priority",
  fifelse(axis_conservation == "Fully_Conserved" & any_atac,
    "strong",
    fifelse(axis_conservation == "Partially_Conserved" & any_disruption,
      "moderate",
      "weak")))]

# ---------- 5. Tidy columns + write ---------------------------------------
out_cols <- c(
  "source", "target", "ligand_complex", "receptor_complex",
  "ligand_first", "receptor_first",
  "score_diff_human",
  "ligand_h_lfc", "receptor_h_lfc",
  "axis_conservation",
  "ligand_motif_disrupted", "receptor_motif_disrupted",
  "ligand_atac_regulated",  "receptor_atac_regulated",
  "gwas_locus",
  "priority_class"
)
setnames(joined, "source", "source_ct", skip_absent = TRUE)
setnames(joined, "target", "target_ct", skip_absent = TRUE)
out_cols[1:2] <- c("source_ct", "target_ct")

out <- joined[, ..out_cols]

# Order: top_priority > strong > moderate > weak, then by score_diff_human.
class_rank <- c("top_priority" = 1L, "strong" = 2L, "moderate" = 3L, "weak" = 4L)
out[, .rank := class_rank[priority_class]]
setorder(out, .rank, -score_diff_human)
out[, .rank := NULL]

fwrite(out, OUT_AXES, sep = "\t")
message(sprintf("[5] Wrote %s (%d rows)", OUT_AXES, nrow(out)))

# Summary
summary_dt <- out[, .(n_axes = .N), by = priority_class]
summary_dt[, rank := class_rank[priority_class]]
setorder(summary_dt, rank)
summary_dt[, rank := NULL]
fwrite(summary_dt, OUT_SUMMARY, sep = "\t")

message("[6] Counts per priority_class:")
print(summary_dt)

message("\n[7] Top 10 top_priority axes:")
top10 <- out[priority_class == "top_priority"][1:min(10L, .N)]
print(top10[, .(source_ct, target_ct, ligand_complex, receptor_complex,
                score_diff_human, gwas_locus)])

message("\nDone.")
