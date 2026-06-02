#!/usr/bin/env Rscript
#SBATCH --job-name=net_284_dlr
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_284_dlr_%j.out
#SBATCH --error=logs/net_284_dlr_%j.err
#
# Script 284: Build D-LR edges from LIANA/CCC differential ligand-receptor data.
# Uses RNA-seq/.../progression/ccc_top_rewired.csv when available, else falls
# back to Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv.
#
# Observed schema for ccc_top_rewired.csv:
#   transition, ligand, receptor, axis, lr_pair,
#   mean_from, mean_to, lfc, pval, n_from, n_to, nonzero_frac, padj, pathway
# transition values include "F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4" etc.

suppressPackageStartupMessages({
  library(data.table)
})

PROJ    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CCC_P   <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/ccc_top_rewired.csv")
LIANA_P <- file.path(PROJ, "Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv")
NODES   <- file.path(PROJ, "RNA-seq/results/network/network_nodes.csv")
OUT     <- file.path(PROJ, "RNA-seq/results/network/edges_d_lr.csv")

# ---- Node set V -----------------------------------------------------------
nodes <- fread(NODES)
V <- unique(nodes[[1]])
cat(sprintf("[284] Node set V: %d genes\n", length(V)))

# ---- Load source ----------------------------------------------------------
if (file.exists(CCC_P)) {
  cat(sprintf("[284] Using ccc_top_rewired.csv\n"))
  src <- fread(CCC_P)
  src_type <- "ccc_top_rewired"
} else {
  cat(sprintf("[284] ccc_top_rewired.csv not found; using LIANA fallback\n"))
  src <- fread(LIANA_P)
  src_type <- "liana"
}
cat(sprintf("[284] Source rows: %d   cols: %s\n", nrow(src),
            paste(colnames(src), collapse = ",")))

# ---- Parse per-transition scores into F01 / F2 / F34 stages ---------------
split_complex <- function(x) {
  # ligand/receptor may be a complex like "ITGA1_ITGB1" -> split on "_"
  strsplit(as.character(x), "_", fixed = TRUE)
}

if (src_type == "ccc_top_rewired") {
  stopifnot(all(c("transition", "ligand", "receptor", "axis", "lfc", "padj")
                %in% colnames(src)))

  # Map transitions -> stage comparison classes
  # F0_to_F1 -> within F01
  # F1_to_F2 -> F01 -> F2 (onset)
  # F2_to_F3 -> F2 -> F34 (post-switch)
  # F3_to_F4 -> within F34
  src[, stage_emerging := fcase(
    transition %in% c("F1_to_F2"), "F2_emerging_onset",
    transition %in% c("F2_to_F3"), "F2_resolving_or_post",
    transition == "F0_to_F1",     "pre_F2",
    transition == "F3_to_F4",     "post_F2",
    default = NA_character_
  )]

  # F2-specific if it is onset (F1->F2) OR if |lfc| at F1_to_F2 is dominant
  # over adjacent transitions. We compute a per-LR "F2 specificity" score.
  # Build wide table: lfc per transition per lr_pair+axis.
  wide <- dcast(src, lr_pair + ligand + receptor + axis ~ transition,
                value.var = "lfc", fun.aggregate = function(x) x[which.max(abs(x))],
                fill = NA_real_)

  # Stage-bin scores (F01, F2, F34) from available transitions
  # score_F01 = lfc(F0_to_F1), score_F2 = lfc(F1_to_F2), score_F34 = lfc(F3_to_F4)
  for (col in c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"))
    if (!col %in% colnames(wide)) wide[, (col) := NA_real_]

  wide[, score_F01 := get("F0_to_F1")]
  wide[, score_F2  := get("F1_to_F2")]
  wide[, score_F34 := get("F3_to_F4")]

  # score_diff: F2 transition magnitude vs mean of adjacent (absolute)
  wide[, adj_mean := rowMeans(cbind(abs(score_F01), abs(score_F34)), na.rm = TRUE)]
  wide[, score_diff := abs(score_F2) - adj_mean]

  # F2-specific: either explicit F2 transition with padj support, or top 10%
  # by score_diff.
  thr <- quantile(wide$score_diff, 0.90, na.rm = TRUE)
  cat(sprintf("[284] score_diff top-10%% threshold: %.3f\n", thr))

  # bring in padj for F1_to_F2 (onset) for labelling
  padj_dt <- src[transition == "F1_to_F2",
                 .(lr_pair, axis, padj_onset = padj, lfc_onset = lfc)]
  wide <- merge(wide, padj_dt, by = c("lr_pair", "axis"), all.x = TRUE)

  wide[, is_f2_specific := (
    (!is.na(padj_onset) & padj_onset < 0.05 & abs(lfc_onset) > 0) |
    (!is.na(score_diff) & score_diff >= thr)
  )]

  f2 <- wide[is_f2_specific == TRUE]
  cat(sprintf("[284] F2-specific LR pairs (raw rows): %d\n", nrow(f2)))

  # Determine emergence_stage label
  f2[, emergence_stage := fcase(
    !is.na(padj_onset) & padj_onset < 0.05 & lfc_onset > 0, "F2_emerging",
    !is.na(padj_onset) & padj_onset < 0.05 & lfc_onset < 0, "F2_declining",
    default = "F2_specific"
  )]

  # Expand complexes; create canonical gene_a (ligand) / gene_b (receptor)
  expanded <- f2[, {
    lig_parts <- split_complex(ligand)[[1]]
    rec_parts <- split_complex(receptor)[[1]]
    grid <- CJ(ga = lig_parts, gb = rec_parts)
    list(gene_a = grid$ga, gene_b = grid$gb,
         score_F01 = score_F01, score_F2 = score_F2, score_F34 = score_F34,
         score_diff = score_diff,
         cell_type_pair = axis,
         emergence_stage = emergence_stage)
  }, by = .(lr_pair, axis, score_F01, score_F2, score_F34,
            score_diff, emergence_stage, ligand, receptor)]

  edges <- expanded[, .(gene_a, gene_b, score_diff, score_F01, score_F2,
                        score_F34, cell_type_pair, emergence_stage)]

} else {
  # Fallback: LIANA differential file - inspect columns defensively
  cn <- colnames(src)
  lig_col <- intersect(c("ligand", "ligand_complex", "source_genesymbol", "ligand_symbol"), cn)[1]
  rec_col <- intersect(c("receptor", "receptor_complex", "target_genesymbol", "receptor_symbol"), cn)[1]
  stopifnot(!is.na(lig_col), !is.na(rec_col))

  diff_col <- intersect(c("score_diff", "lfc", "mean_diff", "delta", "diff"), cn)[1]
  ct_col   <- intersect(c("cell_type_pair", "axis", "source_target", "source"), cn)[1]
  sF01 <- intersect(c("score_F01", "score_control"), cn)[1]
  sF2  <- intersect(c("score_F2"), cn)[1]
  sF34 <- intersect(c("score_F34", "score_masld"), cn)[1]

  if (is.na(diff_col)) {
    # Build difference from MASLD vs control if present
    if (!is.na(sF34) && !is.na(sF01)) {
      src[, score_diff := get(sF34) - get(sF01)]
      diff_col <- "score_diff"
    } else stop("Cannot find differential score column in LIANA file")
  }

  thr <- quantile(abs(src[[diff_col]]), 0.90, na.rm = TRUE)
  f2 <- src[abs(get(diff_col)) >= thr]

  expanded <- f2[, {
    lig_parts <- split_complex(get(lig_col))[[1]]
    rec_parts <- split_complex(get(rec_col))[[1]]
    grid <- CJ(ga = lig_parts, gb = rec_parts)
    list(gene_a = grid$ga, gene_b = grid$gb,
         score_diff = get(diff_col),
         score_F01 = if (!is.na(sF01)) get(sF01) else NA_real_,
         score_F2  = if (!is.na(sF2))  get(sF2)  else NA_real_,
         score_F34 = if (!is.na(sF34)) get(sF34) else NA_real_,
         cell_type_pair = if (!is.na(ct_col)) as.character(get(ct_col)) else NA_character_,
         emergence_stage = "F2_specific")
  }, by = seq_len(nrow(f2))]

  edges <- expanded[, .(gene_a, gene_b, score_diff, score_F01, score_F2,
                        score_F34, cell_type_pair, emergence_stage)]
}

cat(sprintf("[284] Gene-level edges before filter: %d\n", nrow(edges)))

# ---- Filter to V ----------------------------------------------------------
edges <- edges[gene_a %in% V & gene_b %in% V & gene_a != "" & gene_b != ""]

# Canonical ordering (ligand stays as gene_a by biological convention, but
# enforce alphabetical to match other edge tables and enable deduplication)
swap <- edges$gene_a > edges$gene_b
if (any(swap)) {
  tmp_a <- edges$gene_a[swap]
  edges[swap, gene_a := gene_b[swap]]
  edges[swap, gene_b := tmp_a]
}

# Deduplicate (keep max |score_diff|); aggregate cell-type pairs
edges[, .abs_score_diff := abs(score_diff)]
setorder(edges, gene_a, gene_b, -.abs_score_diff)
edges[, .abs_score_diff := NULL]
dedup <- edges[, {
  keep <- 1L
  list(
    score_diff      = score_diff[keep],
    score_F01       = score_F01[keep],
    score_F2        = score_F2[keep],
    score_F34       = score_F34[keep],
    cell_type_pairs = paste(unique(na.omit(cell_type_pair)), collapse = ";"),
    emergence_stage = emergence_stage[keep]
  )
}, by = .(gene_a, gene_b)]

cat(sprintf("[284] Final D-LR edges: %d\n", nrow(dedup)))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(dedup, OUT)
cat(sprintf("[284] Wrote %s\n", OUT))
cat(sprintf("  median |score_diff|: %.3f\n",
            median(abs(dedup$score_diff), na.rm = TRUE)))
cat(sprintf("  emergence_stage distribution:\n"))
print(table(dedup$emergence_stage))
