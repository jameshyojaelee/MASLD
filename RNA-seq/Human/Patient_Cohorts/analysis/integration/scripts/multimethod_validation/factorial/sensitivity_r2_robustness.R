#!/usr/bin/env Rscript
# ===========================================================================
# SENSITIVITY robustness of the interim winner — does the top method survive
# (A) DIRECTIONAL mouse precision and (B) TOP-N size-matched DEG sets?
# Read-only: recomputes ONLY the R2 axis under refined definitions, holds the
# frozen LOCO/R1 (loco_repro_scalar) fixed, re-ranks with the SAME composite
# formula. Does NOT touch the frozen primary outputs. Writes
# sensitivity_r2_robustness.csv.
# ===========================================================================
suppressMessages({library(data.table)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DEGX <- file.path(BASE, "RNA-seq/results/degx_factorial")
GENEMETA <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
OT_FILE  <- file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv")
ATLAS    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
TOPN <- 1000L
strip_v <- function(x) sub("[.][0-9]+$", "", x)

# ---- gene symbol -> ensembl-base map ----
gm <- fread(GENEMETA); gm[, eb := strip_v(gene_id)]
sym2eb <- gm[!is.na(gene_name) & gene_name != "", .(gene_name, eb)]
map_sym <- function(s) sym2eb[match(s, gene_name), eb]

# ---- OpenTargets panel ----
ot <- fread(OT_FILE, skip = "gene_symbol")
ot_eb <- unique(na.omit(map_sym(unique(ot$gene_symbol))))

# ---- atlas: mouse (padj + logFC) + coloc, universe ----
atlas <- fread(ATLAS, select = c("ensembl_id","mouse_meta_padj","mouse_meta_logFC","coloc_best_susie_pp4_polyfun"))
atlas[, gb := strip_v(ensembl_id)]; atlas <- atlas[!is.na(gb) & gb != ""][!duplicated(gb)]
atlas[, mouse_tested := !is.na(mouse_meta_padj)]
atlas[, is_mouse_deg := mouse_tested & mouse_meta_padj < 0.05]
atlas[, is_coloc := !is.na(coloc_best_susie_pp4_polyfun) & coloc_best_susie_pp4_polyfun > 0.5]
universe   <- atlas$gb
mouse_univ <- atlas[mouse_tested == TRUE, gb]
mouse_up   <- atlas[is_mouse_deg & mouse_meta_logFC > 0, gb]
mouse_dn   <- atlas[is_mouse_deg & mouse_meta_logFC < 0, gb]
mouse_any  <- atlas[is_mouse_deg == TRUE, gb]
coloc_truth<- atlas[is_coloc == TRUE, gb]
ot_truth   <- intersect(ot_eb, universe)

fisher_or <- function(hit, truth, univ) {
  hit <- intersect(hit, univ); truth <- intersect(truth, univ)
  a <- length(intersect(hit, truth)); b <- length(setdiff(hit, truth))
  c <- length(setdiff(truth, hit)); d <- length(univ) - a - b - c
  ft <- tryCatch(fisher.test(matrix(c(a,b,c,max(d,0)),2)), error=function(e) NULL)
  if (is.null(ft)) NA_real_ else unname(ft$estimate)
}

# ---- eligible + gate-pass survivors from the primary winners table ----
win <- fread(file.path(DEGX, "degmethod_ranked_winners.csv"))
surv <- win[gate_pass == TRUE & eligible == TRUE,
            .(cell_id, loco_repro_scalar,
              ot_up_prim = ot_enrichment_or_up, mouse_prim = mouse_precision,
              coloc_prim = coloc_or, n_deg_prim = n_deg)]
cat(sprintf("[sens] %d gate-pass eligible survivors\n", nrow(surv)))

rows <- list()
for (i in seq_len(nrow(surv))) {
  cid <- surv$cell_id[i]
  f <- file.path(DEGX, "cells", paste0("cell_", cid, ".csv"))
  if (!file.exists(f)) next
  d <- fread(f); d[, gb := strip_v(gene)]; d <- d[!is.na(gb) & gb != ""]
  d <- d[!is.na(pval)]
  # primary DEG set (replicate): padj<0.05 & |logFC|>0.5
  deg    <- d[padj < 0.05 & abs(logFC) > 0.5, gb]
  deg_up <- d[padj < 0.05 & logFC >  0.5, gb]
  # top-N by p-value (size-matched), directional split
  setorder(d, pval)
  topN   <- d[seq_len(min(TOPN, .N))]
  tn     <- topN$gb
  tn_up  <- topN[logFC > 0, gb]
  # ---- (A) directional mouse precision on the PRIMARY DEG set ----
  dm <- merge(data.table(gb = deg), d[, .(gb, hlfc = logFC)], by = "gb")
  dm <- merge(dm, atlas[, .(gb, is_mouse_deg, mlfc = mouse_meta_logFC)], by = "gb")
  dm <- dm[!is.na(is_mouse_deg)]                       # restrict to mouse-tested
  mouse_dir <- if (nrow(dm)) mean(dm$is_mouse_deg & sign(dm$hlfc) == sign(dm$mlfc)) else NA_real_
  mouse_undir_chk <- if (nrow(dm)) mean(dm$is_mouse_deg) else NA_real_
  # ---- (B) top-N size-matched R2 ----
  ot_up_tn  <- fisher_or(tn_up, ot_truth, universe)
  coloc_tn  <- fisher_or(tn,    coloc_truth, universe)
  dmn <- merge(data.table(gb = tn), d[, .(gb, hlfc = logFC)], by = "gb")
  dmn <- merge(dmn, atlas[, .(gb, is_mouse_deg, mlfc = mouse_meta_logFC)], by = "gb")[!is.na(is_mouse_deg)]
  mouse_tn_dir <- if (nrow(dmn)) mean(dmn$is_mouse_deg & sign(dmn$hlfc) == sign(dmn$mlfc)) else NA_real_
  rows[[length(rows)+1]] <- data.table(
    cell_id = cid, loco = surv$loco_repro_scalar[i],
    ot_up_prim = surv$ot_up_prim[i], mouse_prim = surv$mouse_prim[i], coloc_prim = surv$coloc_prim[i],
    mouse_dir = mouse_dir, mouse_undir_chk = mouse_undir_chk,
    ot_up_tn = ot_up_tn, mouse_tn_dir = mouse_tn_dir, coloc_tn = coloc_tn,
    n_deg = surv$n_deg_prim[i])
}
res <- rbindlist(rows)

# ---- re-rank under three R2 definitions (R1 = primary loco, held fixed) ----
rk <- function(x, desc = TRUE) frank(if (desc) -x else x, ties.method = "average", na.last = TRUE)
r1 <- rk(res$loco)
composite <- function(ot, mou, col) 0.5*r1 + 0.5*((rk(ot)+rk(mou)+rk(col))/3)
res[, comp_primary    := composite(ot_up_prim, mouse_prim, coloc_prim)]   # replicate of frozen
res[, comp_dir_mouse  := composite(ot_up_prim, mouse_dir,  coloc_prim)]   # (A) directional mouse swapped in
res[, comp_topN       := composite(ot_up_tn,   mouse_tn_dir, coloc_tn)]   # (B) size-matched + directional
for (cc in c("comp_primary","comp_dir_mouse","comp_topN")) res[, paste0("rank_", sub("comp_","",cc)) := frank(get(cc), ties.method="min")]

setorder(res, comp_primary)
fwrite(res, file.path(DEGX, "sensitivity_r2_robustness.csv"))
cat("\n=== TOP 5 by each composite (lower=better) ===\n")
show <- function(by) { setorderv(res, by); cat(sprintf("\n-- ranked by %s --\n", by))
  print(res[1:5, .(cell_id, loco=round(loco,3), ot_up_prim=round(ot_up_prim,2),
    mouse_prim=round(mouse_prim,3), mouse_dir=round(mouse_dir,3),
    ot_up_tn=round(ot_up_tn,2), mouse_tn_dir=round(mouse_tn_dir,3), coloc_tn=round(coloc_tn,2))], row.names=FALSE) }
show("comp_primary"); show("comp_dir_mouse"); show("comp_topN")
cat("\n=== where does limma_trend__C2__kna land under each? ===\n")
print(res[cell_id=="limma_trend__C2__kna", .(cell_id, rank_primary, rank_dir_mouse, rank_topN)], row.names=FALSE)
cat("\n=== and limma_voom_qw__C2__kna / metafor_re__C2__kna ===\n")
print(res[cell_id %in% c("limma_voom_qw__C2__kna","metafor_re__C2__kna"), .(cell_id, rank_primary, rank_dir_mouse, rank_topN)], row.names=FALSE)
