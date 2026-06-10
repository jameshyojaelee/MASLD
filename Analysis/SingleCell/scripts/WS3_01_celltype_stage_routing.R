#!/usr/bin/env Rscript
# ============================================================================
# WS3_01_celltype_stage_routing.R
#
# EXPRESSION axis of the cell-type routing analysis.
#
# Question: for each STAGE-SPECIFIC bulk DEG, which cell type(s) express and
# carry it, STAGE-RESOLVED, on the co-primary stage axis? DATA-DECIDES:
# report the honest hepatocyte vs non-hepatocyte share, whichever way it falls.
#
# Reuses (NO new pseudobulk, NO GPU, NO 348b re-run):
#   * Stage-DEG sets (bulk, Tier-2 stage contrasts):
#       disease_signatures/stagedeg_routing/stage_deg_sets.csv
#   * Per-cell-type COARSE stage DE (disease_stage_coarse bins):
#       pseudobulk_de/{CT}_{Steatosis_vs_Healthy|Steatohepatitis_vs_Steatosis|
#                          Cirrhosis_vs_Steatohepatitis}_de.csv
#   * Per-cell-type FINE/augmented stage DE (F_stage_augmented transitions,
#     348b dream output):
#       disease_signatures/celltype_fstage_dream/
#         per_celltype_per_transition_logFC_augmented.csv
#
# Tiering logic reused verbatim from RNA-seq/80_celltype_intrinsic_attribution.R
# (signed_score = sign(bulk_stage_lfc) * scRNA_celltype_t; sig_concordant =
# celltype_padj<0.05 & same sign; primary_celltype = max positive signed_score;
# tiers <CT>_intrinsic_strong / _intrinsic_likely / multi_celltype / bulk_only).
# The ONLY change vs Script 80: the bulk side is the TRANSITION-specific bulk
# logFC (from stage_deg_sets.csv) and the scRNA side is the TRANSITION-specific
# per-cell-type t — not disease-vs-control.
#
# Expression-presence gate: a cell type can only "carry" a gene it expresses.
# The per-cell-type DE engines (both coarse pseudobulk DE and 348b dream) apply
# the SAME edgeR CPM>1-in->=3-donors filter before fitting, so a gene's presence
# in a {CT, transition} DE table IS the expression-presence gate (mirrors 348b
# line 144: rowSums(cpm(dge) > 1) >= 3). `expressed_<CT>` = gene appears in that
# CT's DE table for that transition.
#
# Outputs (disease_signatures/stagedeg_routing/):
#   stage_deg_celltype_routing_matrix.csv   (stage_set x gene x axis; per-CT
#       signed_score + sig_concordant + expressed flags; primary_celltype;
#       attribution_class)
#   routing_summary_by_transition.csv       (cell_type x stage_set x axis: n and
#       fraction of that set's DEGs carried as intrinsic_strong / _likely /
#       multi_celltype share)
#
# Env: rnaseq  (CPU, no GPU)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

DS_DIR   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
ROUT_DIR <- file.path(DS_DIR, "stagedeg_routing")
PB_DE    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
FINE_DE  <- file.path(DS_DIR, "celltype_fstage_dream",
                      "per_celltype_per_transition_logFC_augmented.csv")
stopifnot(dir.exists(ROUT_DIR), dir.exists(PB_DE), file.exists(FINE_DE))

# Thresholds (mirror Script 80) ------------------------------------------------
PVAL_CT <- 0.05  # nominal per-cell-type padj for concordance (limited CT power)

# 9-lineage LIANA common denominator. Maps the scRNA cell_type labels (as they
# appear in the DE tables) to the canonical lineage name. Cell types outside
# this set (Mono+mono derived cells, Circulating NK/NKT, Basophils, DCs) are NOT
# part of the LIANA-9 denominator and are dropped so the routing share is
# comparable across axes.
LIANA9 <- c(
  "Hepatocytes"        = "Hepatocytes",
  "Macrophages"        = "Macrophages",
  "Fibroblasts"        = "Fibroblasts",
  "Endothelial_cells"  = "Endothelial",
  "Endothelial cells"  = "Endothelial",
  "Cholangiocytes"     = "Cholangiocytes",
  "T_cells"            = "T_cells",
  "T cells"            = "T_cells",
  "B_cells"            = "B_cells",
  "B cells"            = "B_cells",
  "Resident_NK"        = "Resident_NK",
  "Resident NK"        = "Resident_NK",
  "Plasma_cells"       = "Plasma",
  "Plasma cells"       = "Plasma"
)
CELLTYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial",
               "Cholangiocytes", "T_cells", "B_cells", "Resident_NK", "Plasma")

harmonize_ct <- function(x) {
  out <- LIANA9[x]
  unname(out)  # NA for non-LIANA-9 labels
}

# ----------------------------------------------------------------------------
# 1. Load bulk stage-DEG sets
# ----------------------------------------------------------------------------
message("[1] Loading bulk stage-DEG sets ...")
sds <- fread(file.path(ROUT_DIR, "stage_deg_sets.csv"))
# gene = unversioned ENSG; symbol; logFC (bulk transition LFC); is_deg (Tier-2)
sds[, ensg := gene]
message(sprintf("    %d rows across %d stage_sets", nrow(sds), uniqueN(sds$stage_set)))

# ----------------------------------------------------------------------------
# 2. Define axis -> (stage_set -> scRNA transition) routing
# ----------------------------------------------------------------------------
# COARSE axis: bulk coarse_* stage_sets <-> scRNA disease_stage_coarse contrasts
coarse_map <- data.table(
  stage_set    = c("coarse_Steatosis", "coarse_SH", "coarse_Cirrhosis"),
  sc_contrast  = c("Steatosis_vs_Healthy", "Steatohepatitis_vs_Steatosis",
                   "Cirrhosis_vs_Steatohepatitis"),
  axis         = "coarse"
)
# FINE/augmented axis: bulk fine_F* stage_sets <-> 348b F-stage transitions
fine_map <- data.table(
  stage_set    = c("fine_F0_F1", "fine_F1_F2", "fine_F2_F3", "fine_F3_F4"),
  sc_contrast  = c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3"),
  axis         = "augmented"
)

# ----------------------------------------------------------------------------
# 3. Load per-cell-type scRNA stage DE for both axes into a long table:
#    (axis, sc_contrast, lineage, ensg/symbol, sc_logFC, sc_t, sc_padj)
# ----------------------------------------------------------------------------
message("[2] Loading per-cell-type COARSE stage DE ...")
load_coarse <- function() {
  files <- list.files(PB_DE, pattern = "_de\\.csv$", full.names = TRUE)
  contrasts <- coarse_map$sc_contrast
  out <- list()
  for (f in files) {
    bn <- basename(f)
    # match {CT}_{contrast}_de.csv for the 3 coarse contrasts
    hit <- contrasts[vapply(contrasts,
                            function(cn) grepl(paste0("_", cn, "_de\\.csv$"), bn),
                            logical(1))]
    if (length(hit) != 1) next
    ct_raw <- sub(paste0("_", hit, "_de\\.csv$"), "", bn)
    lin <- harmonize_ct(ct_raw)
    if (is.na(lin)) next  # not a LIANA-9 lineage
    dt <- fread(f)
    # coarse DE cols: logFC, AveExpr, t_stat, pvalue, padj, B, gene, ...
    if ("t_stat" %in% names(dt)) setnames(dt, "t_stat", "t", skip_absent = TRUE)
    dt <- dt[, .(gene_id = gene, sc_logFC = logFC, sc_t = t, sc_padj = padj)]
    dt[, `:=`(axis = "coarse", sc_contrast = hit, lineage = lin)]
    out[[bn]] <- dt
  }
  rbindlist(out, use.names = TRUE)
}
coarse_de <- load_coarse()
message(sprintf("    coarse: %d rows, %d lineages x %d contrasts",
                nrow(coarse_de), uniqueN(coarse_de$lineage),
                uniqueN(coarse_de$sc_contrast)))

message("[3] Loading per-cell-type FINE/augmented stage DE (348b dream) ...")
fine_raw <- fread(FINE_DE)  # cols: logFC, AveExpr, t, P.Value, adj.P.Val, B, gene, cell_type, transition
fine_raw[, lineage := harmonize_ct(cell_type)]
fine_de <- fine_raw[!is.na(lineage) & transition %in% fine_map$sc_contrast,
                    .(gene_id = gene, sc_logFC = logFC, sc_t = t,
                      sc_padj = adj.P.Val, axis = "augmented",
                      sc_contrast = transition, lineage)]
message(sprintf("    augmented: %d rows, %d lineages x %d contrasts",
                nrow(fine_de), uniqueN(fine_de$lineage),
                uniqueN(fine_de$sc_contrast)))

# ----------------------------------------------------------------------------
# 4. Resolve the scRNA gene id (mix of symbol + ENSG in BOTH coarse and fine)
#    against the bulk stage-DEG ensg/symbol, per axis.
# ----------------------------------------------------------------------------
# Build a long bulk frame with both keys. BOTH the coarse pseudobulk DE and the
# 348b fine dream output store `gene` as a MIX of HGNC symbols (most genes) and
# bare ENSG ids (genes without a symbol), so both axes need the same dual-key
# resolver: symbol-keyed for non-ENSG ids, ensg-keyed for ENSG ids.
bulk_keys <- unique(sds[, .(ensg, symbol)])

attach_bulk_key <- function(sc, axis_name) {
  sc <- copy(sc)
  sc[, is_ensg := grepl("^ENSG", gene_id)]
  # ENSG-keyed
  m_ensg <- merge(sc[is_ensg == TRUE],
                  bulk_keys[, .(ensg)], by.x = "gene_id", by.y = "ensg")
  m_ensg[, ensg := gene_id]
  # symbol-keyed
  m_sym  <- merge(sc[is_ensg == FALSE],
                  bulk_keys, by.x = "gene_id", by.y = "symbol")
  # m_sym now has ensg from bulk_keys
  cols <- c("ensg", "sc_logFC", "sc_t", "sc_padj", "axis", "sc_contrast", "lineage")
  rbind(m_ensg[, ..cols], m_sym[, ..cols], use.names = TRUE)
}

coarse_keyed <- attach_bulk_key(coarse_de, "coarse")
fine_keyed   <- attach_bulk_key(fine_de,   "augmented")

sc_long <- rbind(coarse_keyed, fine_keyed, use.names = TRUE)
# Collapse duplicate (axis, sc_contrast, lineage, ensg): keep most-significant
setorder(sc_long, axis, sc_contrast, lineage, ensg, sc_padj)
sc_long <- sc_long[, .SD[1], by = .(axis, sc_contrast, lineage, ensg)]
message(sprintf("[4] sc_long after key-join + dedup: %d rows", nrow(sc_long)))

# ----------------------------------------------------------------------------
# 5. Per (axis, stage_set, gene): assemble per-cell-type signed_score,
#    sig_concordant, expressed; classify primary + attribution_class.
# ----------------------------------------------------------------------------
axis_map <- rbind(coarse_map, fine_map)

route_one_axis <- function(this_axis) {
  message(sprintf("[5] Routing axis = %s ...", this_axis))
  amap <- axis_map[axis == this_axis]
  per_set <- list()
  for (i in seq_len(nrow(amap))) {
    ss  <- amap$stage_set[i]
    scc <- amap$sc_contrast[i]

    # Bulk side: stage-DEGs for this stage_set (Tier-2 flag is_deg == TRUE)
    bulk_ss <- sds[stage_set == ss & is_deg == TRUE,
                   .(stage_set, ensg, symbol, bulk_lfc = logFC, direction)]
    if (nrow(bulk_ss) == 0) next

    # scRNA side for this contrast (long over lineage). Drop rows where the DE
    # was non-estimable (sc_t == NA) so they do not masquerade as "expressed".
    sc_ss <- sc_long[axis == this_axis & sc_contrast == scc & !is.na(sc_t),
                     .(ensg, lineage, sc_logFC, sc_t, sc_padj)]

    # NON-ESTIMABLE GUARD: if the scRNA contrast carries no usable per-cell-type
    # statistics for ANY lineage (e.g. coarse Cirrhosis_vs_Steatohepatitis,
    # where 3-4 cirrhosis donors are perfectly confounded with `dataset` and the
    # stage coefficient is not estimable), routing is impossible. Emit the bulk
    # DEGs flagged `non_estimable` so they are NOT mislabelled bulk_only.
    estimable <- nrow(sc_ss) > 0

    base <- unique(bulk_ss[, .(stage_set, ensg, symbol, bulk_lfc, direction)])
    setkey(base, ensg)

    if (!estimable) {
      base[, primary_celltype    := NA_character_]
      base[, primary_score       := NA_real_]
      base[, n_sig_concordant_ct := NA_integer_]
      base[, attribution_class   := "non_estimable"]
      for (ct in CELLTYPES) {
        set(base, j = paste0("score_", ct), value = NA_real_)
        set(base, j = paste0("sig_", ct),   value = NA)
        set(base, j = paste0("expr_", ct),  value = NA)
      }
      set(base, j = "axis", value = this_axis)
      set(base, j = "sc_contrast", value = scc)
      per_set[[ss]] <- base
      message(sprintf("    [%s/%s] NON-ESTIMABLE scRNA contrast (no per-CT stats) -- %d DEGs flagged",
                      this_axis, ss, nrow(base)))
      next
    }

    # Per-row metrics on the long join
    m <- merge(base[, .(ensg, bulk_lfc)], sc_ss, by = "ensg",
               all.x = FALSE, allow.cartesian = TRUE)
    m[, signed_score   := sign(bulk_lfc) * sc_t]
    m[, expressed      := TRUE]  # row present => CT expressed gene (passed CPM gate)
    m[, sig_concordant := !is.na(sc_padj) & sc_padj < PVAL_CT &
                          !is.na(sc_logFC) & sign(sc_logFC) == sign(bulk_lfc)]

    # Pivot lineage -> wide per-CT columns, per gene
    score_w <- dcast(m, ensg ~ lineage, value.var = "signed_score")
    sig_w   <- dcast(m, ensg ~ lineage, value.var = "sig_concordant")
    exp_w   <- dcast(m, ensg ~ lineage, value.var = "expressed")

    ensure_cols <- function(dt) {
      for (ct in CELLTYPES) if (!ct %in% names(dt)) set(dt, j = ct, value = NA)
      dt
    }
    score_w <- ensure_cols(score_w)
    sig_w   <- ensure_cols(sig_w)
    exp_w   <- ensure_cols(exp_w)
    setkey(score_w, ensg); setkey(sig_w, ensg); setkey(exp_w, ensg)

    # Align matrices to ALL bulk DEGs (genes with no expressed CT -> all-NA row,
    # classified not_expressed; never silently dropped).
    sc_mat  <- as.matrix(score_w[base$ensg, ..CELLTYPES])
    sig_mat <- as.matrix(sig_w[base$ensg, ..CELLTYPES]) == TRUE
    sig_mat[is.na(sig_mat)] <- FALSE
    exp_mat <- as.matrix(exp_w[base$ensg, ..CELLTYPES]) == TRUE
    exp_mat[is.na(exp_mat)] <- FALSE
    rownames(sc_mat) <- rownames(sig_mat) <- rownames(exp_mat) <- base$ensg

    # classify (Script 80 logic; + not_expressed for genes absent in every CT)
    classify_row <- function(scores, sigs, anyexpr) {
      if (!anyexpr) return(c(NA_character_, NA_character_, "0", "not_expressed"))
      scores[is.na(scores)] <- -Inf
      best <- which.max(scores)
      best_score <- scores[best]
      best_ct <- CELLTYPES[best]
      n_sig <- sum(sigs, na.rm = TRUE)
      class <- if (best_score <= 0 || is.infinite(best_score)) {
        "bulk_only"
      } else if (n_sig >= 2) {
        "multi_celltype"
      } else if (isTRUE(sigs[best])) {
        paste0(best_ct, "_intrinsic_strong")
      } else {
        paste0(best_ct, "_intrinsic_likely")
      }
      c(best_ct, as.character(best_score), as.character(n_sig), class)
    }
    any_expr <- rowSums(exp_mat) > 0
    cl <- vapply(seq_len(nrow(base)),
                 function(j) classify_row(sc_mat[j, ], sig_mat[j, ], any_expr[j]),
                 character(4))

    base[, primary_celltype    := cl[1, ]]
    base[, primary_score       := as.numeric(cl[2, ])]
    base[, n_sig_concordant_ct := as.integer(cl[3, ])]
    base[, attribution_class   := cl[4, ]]
    base[is.na(primary_score) | primary_score <= 0 | is.infinite(primary_score),
         primary_celltype := NA]

    # attach per-CT signed_score / sig / expressed columns
    for (ct in CELLTYPES) {
      set(base, j = paste0("score_", ct), value = sc_mat[, ct])
      set(base, j = paste0("sig_", ct),   value = sig_mat[, ct])
      set(base, j = paste0("expr_", ct),  value = exp_mat[, ct])
    }
    set(base, j = "axis", value = this_axis)
    set(base, j = "sc_contrast", value = scc)
    per_set[[ss]] <- base
  }
  rbindlist(per_set, use.names = TRUE, fill = TRUE)
}

routing_coarse <- route_one_axis("coarse")
routing_fine   <- route_one_axis("augmented")
routing <- rbind(routing_coarse, routing_fine, use.names = TRUE, fill = TRUE)

# column ordering: identity, classification, then per-CT blocks
id_cols  <- c("axis", "stage_set", "sc_contrast", "ensg", "symbol",
              "bulk_lfc", "direction",
              "primary_celltype", "primary_score", "n_sig_concordant_ct",
              "attribution_class")
ct_cols  <- as.vector(rbind(paste0("score_", CELLTYPES),
                            paste0("sig_",   CELLTYPES),
                            paste0("expr_",  CELLTYPES)))
setcolorder(routing, c(id_cols, ct_cols))

out_matrix <- file.path(ROUT_DIR, "stage_deg_celltype_routing_matrix.csv")
fwrite(routing, out_matrix)
message(sprintf("[out] wrote %s (%d rows)", out_matrix, nrow(routing)))

# ----------------------------------------------------------------------------
# 6. routing_summary_by_transition.csv
#    cell_type x stage_set x axis: n + fraction of that set's DEGs it carries as
#    intrinsic_strong / intrinsic_likely; plus multi_celltype share (set-level).
# ----------------------------------------------------------------------------
message("[6] Building routing_summary_by_transition ...")
summ_rows <- list()
for (ax in c("coarse", "augmented")) {
  for (ss in unique(routing[axis == ax, stage_set])) {
    sub <- routing[axis == ax & stage_set == ss]
    n_total <- nrow(sub)
    sc_contrast_lbl <- sub$sc_contrast[1]
    estimable_set <- !all(sub$attribution_class == "non_estimable")
    n_multi <- sum(sub$attribution_class == "multi_celltype")
    n_bulkonly <- sum(sub$attribution_class == "bulk_only")
    n_notexpr <- sum(sub$attribution_class == "not_expressed")
    n_nonest  <- sum(sub$attribution_class == "non_estimable")
    for (ct in CELLTYPES) {
      # carried as STRONG = primary_celltype == ct & intrinsic_strong
      n_strong <- sum(sub$primary_celltype == ct &
                      sub$attribution_class == paste0(ct, "_intrinsic_strong"),
                      na.rm = TRUE)
      n_likely <- sum(sub$primary_celltype == ct &
                      sub$attribution_class == paste0(ct, "_intrinsic_likely"),
                      na.rm = TRUE)
      # ANY sig_concordant in this CT (regardless of primary) — broader carrier
      n_sig_any <- sum(sub[[paste0("sig_", ct)]], na.rm = TRUE)
      n_expr    <- sum(sub[[paste0("expr_", ct)]], na.rm = TRUE)
      summ_rows[[paste(ax, ss, ct, sep = "|")]] <- data.table(
        axis = ax, stage_set = ss, sc_contrast = sc_contrast_lbl,
        estimable = estimable_set,
        cell_type = ct, n_set_degs = n_total,
        n_intrinsic_strong = n_strong,
        n_intrinsic_likely = n_likely,
        frac_intrinsic_strong = ifelse(n_total > 0, n_strong / n_total, NA_real_),
        frac_intrinsic_likely = ifelse(n_total > 0, n_likely / n_total, NA_real_),
        n_sig_concordant_any  = n_sig_any,
        frac_sig_concordant_any = ifelse(n_total > 0, n_sig_any / n_total, NA_real_),
        n_expressed = n_expr,
        frac_expressed = ifelse(n_total > 0, n_expr / n_total, NA_real_),
        n_multi_celltype = n_multi,
        frac_multi_celltype = ifelse(n_total > 0, n_multi / n_total, NA_real_),
        n_bulk_only = n_bulkonly,
        frac_bulk_only = ifelse(n_total > 0, n_bulkonly / n_total, NA_real_),
        n_not_expressed = n_notexpr,
        n_non_estimable = n_nonest
      )
    }
  }
}
routing_summary <- rbindlist(summ_rows, use.names = TRUE)
out_summary <- file.path(ROUT_DIR, "routing_summary_by_transition.csv")
fwrite(routing_summary, out_summary)
message(sprintf("[out] wrote %s (%d rows)", out_summary, nrow(routing_summary)))

# ----------------------------------------------------------------------------
# 7. Hep-vs-non-hep honest share per transition (console + CSV companion)
# ----------------------------------------------------------------------------
message("[7] Hep vs summed-non-hep carrier share per transition ...")
hep_share <- routing[, {
  n_hep   <- sum(primary_celltype == "Hepatocytes" &
                 grepl("^Hepatocytes_intrinsic", attribution_class), na.rm = TRUE)
  n_nonhep <- sum(!is.na(primary_celltype) &
                  primary_celltype != "Hepatocytes" &
                  grepl("_intrinsic_(strong|likely)$", attribution_class),
                  na.rm = TRUE)
  .(estimable = !all(attribution_class == "non_estimable"),
    n_set_degs = .N,
    n_multi_celltype = sum(attribution_class == "multi_celltype"),
    n_bulk_only = sum(attribution_class == "bulk_only"),
    n_not_expressed = sum(attribution_class == "not_expressed"),
    n_non_estimable = sum(attribution_class == "non_estimable"),
    n_hep_intrinsic = n_hep,
    n_nonhep_intrinsic = n_nonhep,
    hep_share_of_attributed = ifelse((n_hep + n_nonhep) > 0,
                                     n_hep / (n_hep + n_nonhep), NA_real_))
}, by = .(axis, stage_set, sc_contrast)]
setorder(hep_share, axis, stage_set)
out_hep <- file.path(ROUT_DIR, "routing_hep_vs_nonhep_share.csv")
fwrite(hep_share, out_hep)
message(sprintf("[out] wrote %s", out_hep))

cat("\n================ HEP vs NON-HEP INTRINSIC SHARE ================\n")
print(hep_share)

cat("\n================ CARRIER FRACTIONS (intrinsic_strong) ================\n")
strong_wide <- dcast(routing_summary, axis + stage_set + sc_contrast ~ cell_type,
                     value.var = "frac_intrinsic_strong")
print(strong_wide)

message("\nDone.")
