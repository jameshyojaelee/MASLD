#!/usr/bin/env Rscript
# ============================================================================
# WS3_02_liana_stagedeg_enrichment.R
#
# LIANA cell-cell-communication axis of the stage-DEG routing analysis.
#
# QUESTION (DATA-DECIDES): do the stage-specific bulk DEGs route as
# ligands/receptors in stage-ADVANCING signaling circuits, and which
# source->target cell-type circuits do they participate in as stage progresses?
#
# We REUSE the per-donor LIANA mixed-model slopes (Script 346) — we do NOT
# re-run the GPU per-donor LIANA. A "stage-strengthening" LR pair is one whose
# per-donor magnitude signal (-log10 magnitude_rank) increases with stage:
#   - augmented (fine) axis : F_stage_numeric slope > 0 AND padj_within_ct < 0.05
#   - coarse axis           : Cirrhosis-vs-Healthy contrast Estimate > 0 AND
#                             padj_within_ct < 0.05 (latest-stage advancement)
#
# For each stage_deg_set we ask, SEPARATELY for the ligand role and the
# receptor role, whether the set's UP-DEGs are over-represented among the
# genes that sit on the ligand (resp. receptor) side of stage-strengthening
# LR pairs. Background = all genes that appear in that role in the consensus
# resource (i.e. across all LMM-tested LR pairs on that axis). Direction is
# matched: a stage-UP bulk DEG should sit in a strengthening (up-with-stage)
# circuit.
#
# Outputs (Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing/):
#   liana_stagedeg_ligand_receptor_enrichment.csv
#       stage_set x {ligand,receptor} x axis x OR/p/q/n
#   stagedeg_circuit_participation.csv
#       gene x role x source->target x stage-slope (for stage-DEGs that ARE a
#       ligand/receptor in a significant stage-strengthening LR pair)
#   circuit_celltype_flow_by_stage.csv
#       source->target x stage_set x axis x n stage-DEG-bearing strengthening
#       LR pairs (how the dominant circuit shifts across the stage axis)
#
# Env: rnaseq.  CPU only (reuses parquets + LMM TSVs; no GPU).
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CCC_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
OUT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_THRESH <- 0.05   # strengthening LR pair significance (padj_within_ct)

# ----------------------------------------------------------------------------
# 1. Inputs
# ----------------------------------------------------------------------------
deg_path <- file.path(OUT_DIR, "stage_deg_sets.csv")
cat(sprintf("[input] stage DEG sets: %s\n", deg_path))
deg <- fread(deg_path)
# Keep only the gene-level fields we need; symbol is the LR vocabulary key.
deg <- deg[, .(stage_set, axis_set = axis, gene, symbol, direction, is_deg)]
cat(sprintf("[input] %d stage_set x gene rows, %d stage_sets\n",
            nrow(deg), uniqueN(deg$stage_set)))

# LMM slope tables, one per CCC axis we test.
lmm_files <- list(
  coarse    = file.path(CCC_DIR, "stage_lr_lmm_coarse.tsv"),
  augmented = file.path(CCC_DIR, "stage_lr_lmm_fstage_augmented.tsv")
)
# Fall back to the inferred-axis canonical fstage if the augmented re-fit is
# absent (so the script still runs end-to-end), but warn loudly.
if (!file.exists(lmm_files$augmented)) {
  fallback <- file.path(CCC_DIR, "stage_lr_lmm_fstage.tsv")
  warning("[axis] augmented LMM not found at ", lmm_files$augmented,
          "; falling back to canonical (inferred-axis) ", fallback,
          ". Re-run WS3_01_refit_lmm_augmented.sh for the true augmented axis.")
  lmm_files$augmented <- fallback
}

# ----------------------------------------------------------------------------
# 2. Complex parsing — split LIANA complexes on '_' (as Script 346 does)
# ----------------------------------------------------------------------------
split_complex <- function(x) strsplit(x, "_", fixed = TRUE)

# ----------------------------------------------------------------------------
# 3. Helper: from one LMM table, derive the per-(role) strengthening universe.
#    Returns a list with:
#      - lig_bg / rec_bg  : all genes that appear as ligand / receptor in any
#                           LMM-tested LR pair on this axis (the consensus
#                           resource background, restricted to what was tested)
#      - lig_str / rec_str: genes appearing as ligand / receptor in a
#                           STRENGTHENING LR pair (positive slope, padj<thr)
#      - str_pairs        : data.table of strengthening LR pairs with the
#                           role->gene expansion + slope + ct_pair
# ----------------------------------------------------------------------------
build_axis_universe <- function(lmm_path, axis) {
  cat(sprintf("\n[axis %s] reading %s\n", axis, lmm_path))
  L <- fread(lmm_path)

  # Identify the stage-advancement term per axis.
  if (axis == "coarse") {
    # Latest-stage contrast = Cirrhosis vs Healthy (monotone stage advancement).
    L <- L[term == "disease_stage_coarseCirrhosis"]
    slope_label <- "Cirrhosis_vs_Healthy"
  } else {
    L <- L[term == "F_stage_numeric"]
    slope_label <- "F_stage_slope"
  }
  cat(sprintf("[axis %s] %d LR-pair rows on the stage-advancement term (%s)\n",
              axis, nrow(L), slope_label))

  # One row per (ct_pair, lr_pair) on this term.
  L[, slope := Estimate]
  L[, padj  := padj_within_ct]
  L[, strengthening := (slope > 0) & is.finite(padj) & (padj < PADJ_THRESH)]
  n_str <- sum(L$strengthening, na.rm = TRUE)
  cat(sprintf("[axis %s] %d strengthening LR-pair circuits (slope>0, padj<%.2f) of %d tested\n",
              axis, n_str, PADJ_THRESH, nrow(L)))

  # Expand ligand and receptor complexes to constituent genes, keeping the
  # ct_pair so we can do circuit attribution.
  expand_role <- function(DT, complex_col, role) {
    g <- DT[, .(gene = unlist(split_complex(get(complex_col)))),
            by = .(ct_pair, source, target, lr_pair,
                   ligand_complex, receptor_complex,
                   slope, padj, strengthening)]
    g[, role := role]
    g
  }
  lig <- expand_role(L, "ligand_complex",   "ligand")
  rec <- expand_role(L, "receptor_complex", "receptor")

  # Background = genes tested in this role (consensus resource ∩ tested LR pairs)
  lig_bg <- unique(lig$gene)
  rec_bg <- unique(rec$gene)
  # Strengthening-role genes
  lig_str <- unique(lig[strengthening == TRUE]$gene)
  rec_str <- unique(rec[strengthening == TRUE]$gene)

  cat(sprintf("[axis %s] ligand background %d genes (%d in strengthening circuits)\n",
              axis, length(lig_bg), length(lig_str)))
  cat(sprintf("[axis %s] receptor background %d genes (%d in strengthening circuits)\n",
              axis, length(rec_bg), length(rec_str)))

  list(
    axis = axis, slope_label = slope_label,
    lig_bg = lig_bg, rec_bg = rec_bg,
    lig_str = lig_str, rec_str = rec_str,
    lig_expand = lig, rec_expand = rec,
    str_pairs = L[strengthening == TRUE]
  )
}

# ----------------------------------------------------------------------------
# 4. Fisher enrichment: are a stage_set's UP-DEGs over-represented among the
#    genes on the {ligand|receptor} side of strengthening LR pairs?
#    Universe restricted to genes that appear in that role in the resource.
# ----------------------------------------------------------------------------
fisher_role <- function(set_up_genes, role_bg, role_str, stage_set, role, axis) {
  # Restrict the DEG set to the role universe (only genes that can act in role).
  up_in_univ <- intersect(set_up_genes, role_bg)
  str_set    <- role_str                       # strengthening-role genes (subset of bg)
  N  <- length(role_bg)                         # universe
  K  <- length(str_set)                         # strengthening-role genes
  n  <- length(up_in_univ)                       # set up-DEGs in universe
  k  <- length(intersect(up_in_univ, str_set))   # set up-DEGs that are strengthening
  if (N == 0 || K == 0 || n == 0) {
    return(data.table(stage_set = stage_set, role = role, axis = axis,
                      n_universe = N, n_role_strengthening = K,
                      n_set_in_universe = n, n_overlap = k,
                      OR = NA_real_, p = NA_real_,
                      expected = NA_real_, log2_enrichment = NA_real_))
  }
  # 2x2: rows = is strengthening (yes/no); cols = is set up-DEG (yes/no)
  a <- k                       # strengthening & set
  b <- K - k                   # strengthening & not-set
  c <- n - k                   # not-strengthening & set
  d <- N - K - (n - k)         # not-strengthening & not-set
  mat <- matrix(c(a, c, b, d), nrow = 2,
                dimnames = list(c("set", "notset"),
                                c("strength", "notstrength")))
  ft <- fisher.test(mat, alternative = "greater")  # over-representation
  expected <- n * K / N
  data.table(stage_set = stage_set, role = role, axis = axis,
             n_universe = N, n_role_strengthening = K,
             n_set_in_universe = n, n_overlap = k,
             OR = unname(ft$estimate), p = ft$p.value,
             expected = expected,
             log2_enrichment = log2((k + 0.5) / (expected + 0.5)))
}

# ----------------------------------------------------------------------------
# 5. Drive both axes
# ----------------------------------------------------------------------------
universes <- list(
  coarse    = build_axis_universe(lmm_files$coarse,    "coarse"),
  augmented = build_axis_universe(lmm_files$augmented, "augmented")
)

stage_sets <- unique(deg$stage_set)

enr_rows  <- list()
circ_rows <- list()
flow_rows <- list()

for (axis in names(universes)) {
  U <- universes[[axis]]
  for (ss in stage_sets) {
    sub <- deg[stage_set == ss]
    # UP-DEGs of this set (direction-matched to strengthening circuits)
    up_genes <- unique(sub[is_deg == TRUE & direction == "up"]$symbol)
    up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]

    enr_rows[[length(enr_rows) + 1]] <-
      fisher_role(up_genes, U$lig_bg, U$lig_str, ss, "ligand",   axis)
    enr_rows[[length(enr_rows) + 1]] <-
      fisher_role(up_genes, U$rec_bg, U$rec_str, ss, "receptor", axis)

    # ----- Circuit attribution: which strengthening circuits do this set's
    # up-DEGs participate in (as ligand or as receptor)? --------------------
    lig_hits <- U$lig_expand[strengthening == TRUE & gene %in% up_genes]
    rec_hits <- U$rec_expand[strengthening == TRUE & gene %in% up_genes]
    if (nrow(lig_hits)) {
      circ_rows[[length(circ_rows) + 1]] <- lig_hits[, .(
        stage_set = ss, axis = axis, gene, role = "ligand",
        source, target, ct_pair, lr_pair,
        ligand_complex, receptor_complex, slope, padj)]
    }
    if (nrow(rec_hits)) {
      circ_rows[[length(circ_rows) + 1]] <- rec_hits[, .(
        stage_set = ss, axis = axis, gene, role = "receptor",
        source, target, ct_pair, lr_pair,
        ligand_complex, receptor_complex, slope, padj)]
    }

    # ----- Cell-type flow: per source->target, how many DISTINCT
    # strengthening LR pairs carry one of this set's up-DEGs (as either role)?
    all_hits <- rbindlist(list(
      if (nrow(lig_hits)) lig_hits[, .(ct_pair, source, target, lr_pair)] else NULL,
      if (nrow(rec_hits)) rec_hits[, .(ct_pair, source, target, lr_pair)] else NULL
    ), use.names = TRUE)
    if (nrow(all_hits)) {
      fl <- unique(all_hits)[, .(
        n_stage_deg_LR_pairs = uniqueN(lr_pair),
        n_genes              = NA_integer_),
        by = .(ct_pair, source, target)]
      # n distinct set up-DEGs participating in this ct_pair
      genes_by_ct <- unique(rbindlist(list(
        if (nrow(lig_hits)) lig_hits[, .(ct_pair, gene)] else NULL,
        if (nrow(rec_hits)) rec_hits[, .(ct_pair, gene)] else NULL
      )))[, .(n_genes = uniqueN(gene)), by = ct_pair]
      fl[, n_genes := genes_by_ct[.SD, on = "ct_pair", x.n_genes]]
      fl[, `:=`(stage_set = ss, axis = axis)]
      flow_rows[[length(flow_rows) + 1]] <- fl
    }
  }
}

# ----------------------------------------------------------------------------
# 6. Assemble + multiple-testing correction
# ----------------------------------------------------------------------------
enr <- rbindlist(enr_rows, use.names = TRUE, fill = TRUE)
# BH within each axis x role family (across stage_sets)
enr[, q := p.adjust(p, method = "BH"), by = .(axis, role)]
setorder(enr, axis, role, p)

circ <- if (length(circ_rows))
  rbindlist(circ_rows, use.names = TRUE, fill = TRUE) else data.table()
flow <- if (length(flow_rows))
  rbindlist(flow_rows, use.names = TRUE, fill = TRUE) else data.table()
if (nrow(flow))
  setorder(flow, axis, stage_set, -n_stage_deg_LR_pairs)

# ----------------------------------------------------------------------------
# 7. Write outputs
# ----------------------------------------------------------------------------
f_enr  <- file.path(OUT_DIR, "liana_stagedeg_ligand_receptor_enrichment.csv")
f_circ <- file.path(OUT_DIR, "stagedeg_circuit_participation.csv")
f_flow <- file.path(OUT_DIR, "circuit_celltype_flow_by_stage.csv")
fwrite(enr,  f_enr)
fwrite(circ, f_circ)
fwrite(flow, f_flow)
cat(sprintf("\n[output] %s (%d rows)\n", f_enr,  nrow(enr)))
cat(sprintf("[output] %s (%d rows)\n",   f_circ, nrow(circ)))
cat(sprintf("[output] %s (%d rows)\n",   f_flow, nrow(flow)))

# ----------------------------------------------------------------------------
# 8. VERIFY: known liver ligand recovery + biological sanity of strengthening
# ----------------------------------------------------------------------------
cat("\n================ VERIFICATION ================\n")
known_ligands <- c("TGFB1", "PDGFB", "PDGFA", "PDGFD", "COL1A1",
                   "SPP1", "CCL2", "TNF", "IL6", "CTGF")
for (axis in names(universes)) {
  U <- universes[[axis]]
  found <- intersect(known_ligands, U$lig_bg)
  in_str <- intersect(known_ligands, U$lig_str)
  cat(sprintf("[verify %s] known ligands in resource: %s\n",
              axis, paste(found, collapse = ", ")))
  cat(sprintf("[verify %s] known ligands in STRENGTHENING circuits: %s\n",
              axis, paste(in_str, collapse = ", ")))
}

cat("\n[verify] significant ligand/receptor enrichments (q<0.05):\n")
sig <- enr[is.finite(q) & q < 0.05]
if (nrow(sig)) {
  print(sig[order(axis, role, q),
            .(axis, stage_set, role, n_overlap, n_set_in_universe,
              n_role_strengthening, OR = round(OR, 2),
              p = signif(p, 3), q = signif(q, 3))])
} else {
  cat("  (none at q<0.05; showing top nominal p per axis x role)\n")
  print(enr[order(axis, role, p),
            head(.SD, 3), by = .(axis, role),
            .SDcols = c("stage_set", "OR", "p", "q", "n_overlap")])
}

cat("\n[verify] top stage-strengthening source->target circuits by stage-DEG load (augmented axis):\n")
if (nrow(flow)) {
  fa <- flow[axis == "augmented"]
  if (nrow(fa)) {
    top_flow <- fa[, .(total_LR = sum(n_stage_deg_LR_pairs)),
                   by = .(ct_pair)][order(-total_LR)][1:min(15, .N)]
    print(top_flow)
  }
}

cat("\n[done] WS3 LIANA stage-DEG routing enrichment complete.\n")
