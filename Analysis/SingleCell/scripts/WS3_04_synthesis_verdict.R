#!/usr/bin/env Rscript
# WS3_04_synthesis_verdict.R
# ---------------------------------------------------------------------------
# Synthesize the three cell-type ROUTING AXES for stage-DEGs into a single
# cell-by-cell verdict, on BOTH stage axes (coarse + augmented/F-stage).
#
#   Axis 1  Expression  — does a cell type carry the stage-DEG as an intrinsic
#                         strong/likely signal, or share it multi-cellularly?
#                         (routing_summary_by_transition.csv)
#   Axis 2  LIANA       — does a stage-DEG route as a ligand/receptor in a
#                         stage-significant CCC circuit sourced/targeted by
#                         that cell type? (stagedeg_circuit_participation.csv)
#   Axis 3  Hotspot     — is the stage-DEG captured by >=1 stage-associated
#                         autocorrelation module of that cell type?
#                         (hotspot_celltype_stagedeg_capture.csv)
#
# DATA-DECIDES framing. The verdict is computed purely from the data: a
# transition is hepatocyte_dominant only if hepatocytes hold >50% of the
# attributed signal on >=2 of the 3 axes; otherwise multi_cellular.
# Coarse-vs-augmented agreement is reported as the robustness result.
#
# Env: rnaseq.  Login-node light (small CSV I/O + tabulation only).
# Outputs:
#   stagedeg_celltype_verdict.csv      (cell_type x axis x transition: 3 fracs + verdict)
#   stagedeg_verdict_by_transition.csv (per transition: hep share each axis,
#                                       carriers, verdict, coarse/aug agreement)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing")
stopifnot(dir.exists(DIR))

# ---------------------------------------------------------------------------
# Conventions: cell-type harmonization + the 5-step transition cascade
# ---------------------------------------------------------------------------
# Canonical cell-type labels (used in expression axis); LIANA uses long forms,
# hotspot uses lowercase. Normalize everything to these 9 canonical labels.
canon_ct <- function(x) {
  x <- as.character(x)
  m <- c(
    "Hepatocytes"="Hepatocytes", "hepatocytes"="Hepatocytes",
    "Macrophages"="Macrophages", "macrophages"="Macrophages",
    "Fibroblasts"="Fibroblasts", "fibroblasts"="Fibroblasts",
    "Endothelial"="Endothelial", "Endothelial cells"="Endothelial", "endothelial"="Endothelial",
    "Cholangiocytes"="Cholangiocytes", "cholangiocytes"="Cholangiocytes",
    "T_cells"="T cells", "T cells"="T cells", "tcells"="T cells", "Tcells"="T cells",
    "B_cells"="B cells", "B cells"="B cells", "bcells"="B cells",
    "Resident_NK"="Resident NK", "Resident NK"="Resident NK", "NK"="Resident NK",
    "Plasma"="Plasma", "plasma"="Plasma"
  )
  out <- m[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

# The two stage axes and their stage_set members, in cascade order.
# coarse axis ("coarse") -> 3 estimable disease-stage transitions
# augmented axis ("augmented") -> 4 fine F-stage transitions (Kleiner)
COARSE_STAGES <- c("coarse_Steatosis", "coarse_SH", "coarse_Cirrhosis")
AUG_STAGES    <- c("fine_F0_F1", "fine_F1_F2", "fine_F2_F3", "fine_F3_F4")

# Human-readable transition label per stage_set
trans_label <- c(
  coarse_Steatosis = "Healthy->Steatosis",
  coarse_SH        = "Steatosis->Steatohepatitis",
  coarse_Cirrhosis = "Steatohepatitis->Cirrhosis",
  fine_F0_F1       = "F0->F1",
  fine_F1_F2       = "F1->F2",
  fine_F2_F3       = "F2->F3",
  fine_F3_F4       = "F3->F4"
)

# Cross-axis transition pairing for coarse-vs-augmented agreement.
# Map each coarse transition to its nearest augmented (F-stage) analog.
#   Healthy->Steatosis        ~ F0->F1   (disease onset / early)
#   Steatosis->Steatohep      ~ F1->F2 + F2->F3 (mid-stage inflammatory expansion)
#   Steatohep->Cirrhosis      ~ F3->F4   (advanced fibrosis -> cirrhosis)
# Pairing used only for the agreement diagnostic; the verdict itself is per
# (axis, transition).
agree_pairs <- list(
  coarse_Steatosis = "fine_F0_F1",
  coarse_SH        = "fine_F2_F3",   # mid-stage inflammatory analog
  coarse_Cirrhosis = "fine_F3_F4"
)

# Non-estimable / low-power flags (documented in routing_hep_vs_nonhep_share.csv)
NON_ESTIMABLE <- c("coarse_Cirrhosis")          # estimable == FALSE on expression axis
LOW_POWER     <- c("fine_F3_F4")                # n_set_degs=1472, smallest set

NONHEP <- function() setdiff(
  c("Hepatocytes","Macrophages","Fibroblasts","Endothelial","Cholangiocytes",
    "T cells","B cells","Resident NK","Plasma"), "Hepatocytes")

# ===========================================================================
# AXIS 1 — Expression routing
#   expr_frac_carried = (n_intrinsic_strong + n_intrinsic_likely + n_multi_celltype)
#                        / n_set_degs
#   i.e. the fraction of this transition's stage-DEGs that the cell type carries
#   as an intrinsic strong/likely signal OR shares as a multi-cell-type signal.
# ===========================================================================
expr <- fread(file.path(DIR, "routing_summary_by_transition.csv"))
expr[, cell_type := canon_ct(cell_type)]
# n_multi_celltype is a transition-level constant (genes shared across >1 CT);
# attribute it to each cell type that participates. To stay honest we treat the
# multi-cellular pool as a shared component every cell type can carry, so it
# counts toward "carried multi-cellularly" for each CT. The intrinsic terms are
# CT-specific.
expr[, expr_frac_carried := (n_intrinsic_strong + n_intrinsic_likely + n_multi_celltype) / n_set_degs]
expr[, stage_axis := fifelse(axis == "coarse", "coarse", "augmented")]
expr_long <- expr[, .(stage_axis, transition = stage_set, cell_type,
                      expr_frac_carried, n_set_degs, estimable)]

# ===========================================================================
# AXIS 2 — LIANA circuit participation
#   liana_frac_as_LR = (# distinct stage-DEGs that route as a ligand on a
#                        circuit SOURCED by that CT, or as a receptor on a
#                        circuit TARGETED by that CT, in a stage-significant
#                        circuit) / n_set_degs
# ===========================================================================
part <- fread(file.path(DIR, "stagedeg_circuit_participation.csv"))
part[, stage_axis := fifelse(axis == "coarse", "coarse", "augmented")]
# Attribute each gene to the CT that emits it (ligand->source) or receives it
# (receptor->target); count distinct genes per CT.
part[, ct_attr := canon_ct(fifelse(role == "ligand", source, target))]
liana_ct <- part[, .(n_lr_genes = uniqueN(gene)),
                 by = .(stage_axis, transition = stage_set, cell_type = ct_attr)]
# Denominator = stage-DEG set size for that transition (from expression summary,
# which equals stage_deg_sets_summary n_DEG).
set_sizes <- unique(expr_long[, .(stage_axis, transition, n_set_degs)])
liana_ct <- merge(liana_ct, set_sizes, by = c("stage_axis","transition"), all.x = TRUE)
liana_ct[, liana_frac_as_LR := n_lr_genes / n_set_degs]

# ===========================================================================
# AXIS 3 — Hotspot module capture
#   hotspot_frac_modularized = frac_deg_captured_prog
#     (fraction of the transition's stage-DEGs captured by >=1 STAGE-ASSOCIATED
#      / progression module of that cell type). Use frac_deg_captured_prog as
#      the primary; frac_deg_captured_any retained for reference.
# ===========================================================================
hot <- fread(file.path(DIR, "hotspot_celltype_stagedeg_capture.csv"))
hot[, cell_type := canon_ct(cell_type)]
hot[, stage_axis := fifelse(axis == "disease_stage_coarse", "coarse", "augmented")]
hot_long <- hot[, .(stage_axis, transition = stage_set, cell_type,
                    hotspot_frac_modularized = frac_deg_captured_prog,
                    hotspot_frac_any = frac_deg_captured_any)]

# ===========================================================================
# COMBINE — long table: cell_type x stage_axis x transition with all 3 fracs
# ===========================================================================
# Restrict to the cascade transitions only (drop ordinal_/prog_ pseudo-sets).
keep_trans <- c(COARSE_STAGES, AUG_STAGES)

m <- merge(expr_long[transition %in% keep_trans],
           liana_ct[, .(stage_axis, transition, cell_type, n_lr_genes, liana_frac_as_LR)],
           by = c("stage_axis","transition","cell_type"), all.x = TRUE)
m <- merge(m, hot_long[transition %in% keep_trans],
           by = c("stage_axis","transition","cell_type"), all.x = TRUE)

# Missing axis values (e.g. a CT with no LIANA circuit, or no hotspot run for
# B/NK/Plasma) -> honest zero participation on that axis.
m[is.na(n_lr_genes),        n_lr_genes := 0L]
m[is.na(liana_frac_as_LR),  liana_frac_as_LR := 0]
m[is.na(hotspot_frac_modularized), hotspot_frac_modularized := NA_real_]  # NA = axis not run for this CT
m[is.na(hotspot_frac_any),  hotspot_frac_any := NA_real_]

m[, transition_label := trans_label[transition]]
m[, non_estimable := transition %in% NON_ESTIMABLE]
m[, low_power      := transition %in% LOW_POWER]

setcolorder(m, c("stage_axis","transition","transition_label","cell_type",
                 "n_set_degs","estimable","non_estimable","low_power",
                 "expr_frac_carried","n_lr_genes","liana_frac_as_LR",
                 "hotspot_frac_modularized","hotspot_frac_any"))
setorder(m, stage_axis, transition, -expr_frac_carried)

fwrite(m, file.path(DIR, "stagedeg_celltype_verdict.csv"))

# ===========================================================================
# PER-TRANSITION VERDICT
#   Hepatocyte vs SUMMED non-hepatocyte share PER AXIS:
#     hep_share = hep_value / sum(all CT values on that axis)
#   For LIANA & hotspot we sum CT-level participation (genes / captured).
#   For expression, n_multi_celltype is a shared pool double-counted across
#   CTs; to compute an HONEST hep-vs-nonhep share we instead use the
#   *intrinsic* attribution (n_hep_intrinsic vs n_nonhep_intrinsic) which is
#   exactly what routing_hep_vs_nonhep_share.csv reports. We read that file
#   directly for the expression-axis hep share (authoritative, non-double-counted).
# ===========================================================================
share_file <- fread(file.path(DIR, "routing_hep_vs_nonhep_share.csv"))
# expression-axis hep share = hep_share_of_attributed (intrinsic attribution)
share_file[, stage_axis := axis]
expr_hepshare <- share_file[stage_set %in% keep_trans,
                            .(stage_axis, transition = stage_set,
                              estimable,
                              expr_hep_share = as.numeric(hep_share_of_attributed),
                              n_hep_intrinsic, n_nonhep_intrinsic)]

build_axis_share <- function(dt, value_col) {
  # per (stage_axis, transition): hep value, summed-nonhep value, hep share
  v <- dt[, .(val = sum(get(value_col), na.rm = TRUE)),
          by = .(stage_axis, transition, cell_type)]
  hep <- v[cell_type == "Hepatocytes", .(stage_axis, transition, hep = val)]
  tot <- v[, .(tot = sum(val, na.rm = TRUE)), by = .(stage_axis, transition)]
  out <- merge(tot, hep, by = c("stage_axis","transition"), all.x = TRUE)
  out[is.na(hep), hep := 0]
  out[, share := fifelse(tot > 0, hep / tot, NA_real_)]
  out
}

liana_share <- build_axis_share(m, "n_lr_genes")
setnames(liana_share, c("hep","tot","share"),
         c("liana_hep_genes","liana_tot_genes","liana_hep_share"))

# hotspot capture: use captured DEG COUNTS (n_deg_captured_prog) for share, so
# CTs with more captured DEGs weigh more. Recompute from hot table directly.
hot_counts <- hot[stage_set %in% keep_trans,
                  .(stage_axis, transition = stage_set, cell_type,
                    n_cap = n_deg_captured_prog)]
hs <- build_axis_share(hot_counts, "n_cap")
setnames(hs, c("hep","tot","share"),
         c("hotspot_hep_captured","hotspot_tot_captured","hotspot_hep_share"))

# Assemble per-transition verdict table
trans_tbl <- unique(m[, .(stage_axis, transition, transition_label,
                          n_set_degs, non_estimable, low_power)])
trans_tbl <- merge(trans_tbl, expr_hepshare[, .(stage_axis, transition, expr_hep_share,
                                                 n_hep_intrinsic, n_nonhep_intrinsic)],
                   by = c("stage_axis","transition"), all.x = TRUE)
trans_tbl <- merge(trans_tbl, liana_share[, .(stage_axis, transition, liana_hep_share,
                                              liana_hep_genes, liana_tot_genes)],
                   by = c("stage_axis","transition"), all.x = TRUE)
trans_tbl <- merge(trans_tbl, hs[, .(stage_axis, transition, hotspot_hep_share,
                                     hotspot_hep_captured, hotspot_tot_captured)],
                   by = c("stage_axis","transition"), all.x = TRUE)

# Dominant carriers per transition (top non-hep + overall), from expression
# intrinsic+multi attribution (axis 1) — the carrier-shift narrative axis.
carriers <- m[, .(stage_axis, transition, cell_type, expr_frac_carried)]
top_carrier <- carriers[order(stage_axis, transition, -expr_frac_carried),
                        .(dominant_carriers = paste(head(cell_type, 3), collapse = "; "),
                          top_nonhep = head(cell_type[cell_type != "Hepatocytes"], 1)),
                        by = .(stage_axis, transition)]
trans_tbl <- merge(trans_tbl, top_carrier, by = c("stage_axis","transition"), all.x = TRUE)

# ---- VERDICT: hepatocyte_dominant if hep>50% on >=2 of 3 axes, else multi_cellular ----
trans_tbl[, n_axes_hepdom := (expr_hep_share > 0.5 & !is.na(expr_hep_share)) +
                             (liana_hep_share > 0.5 & !is.na(liana_hep_share)) +
                             (hotspot_hep_share > 0.5 & !is.na(hotspot_hep_share))]
trans_tbl[, n_axes_estimable := (!is.na(expr_hep_share)) +
                                (!is.na(liana_hep_share)) +
                                (!is.na(hotspot_hep_share))]
trans_tbl[, verdict := fifelse(n_axes_hepdom >= 2, "hepatocyte_dominant", "multi_cellular")]
# Honest override: a transition with the expression axis non-estimable (coarse
# Cirrhosis) carries reduced confidence — flag it explicitly.
trans_tbl[non_estimable == TRUE, verdict := paste0(verdict, " (expr_nonestimable)")]
trans_tbl[low_power == TRUE,     verdict := paste0(verdict, " (low_power)")]

# ---- coarse-vs-augmented AGREEMENT (robustness) ----
# Compare each coarse transition's base verdict to its paired augmented analog.
base_verdict <- function(v) sub(" \\(.*$", "", v)
v_coarse <- trans_tbl[stage_axis == "coarse",
                      .(transition, cverdict = base_verdict(verdict))]
v_aug    <- trans_tbl[stage_axis == "augmented",
                      .(transition, averdict = base_verdict(verdict))]
agree_dt <- rbindlist(lapply(names(agree_pairs), function(ct) {
  at <- agree_pairs[[ct]]
  cv <- v_coarse[transition == ct, cverdict]
  av <- v_aug[transition == at, averdict]
  data.table(coarse_transition = ct,
             augmented_analog  = at,
             coarse_verdict    = if (length(cv)) cv else NA_character_,
             augmented_verdict = if (length(av)) av else NA_character_,
             agree             = if (length(cv) && length(av)) cv == av else NA)
}))

# Attach agreement back onto coarse rows (per-transition convenience), and the
# augmented analog onto augmented rows.
trans_tbl[, coarse_vs_augmented_agreement := NA_character_]
for (i in seq_len(nrow(agree_dt))) {
  ct <- agree_dt$coarse_transition[i]; at <- agree_dt$augmented_analog[i]
  ag <- agree_dt$agree[i]
  lab <- if (is.na(ag)) "NA" else if (ag) "AGREE" else "DISAGREE"
  trans_tbl[stage_axis == "coarse"    & transition == ct,
            coarse_vs_augmented_agreement := sprintf("%s (vs %s)", lab, at)]
  trans_tbl[stage_axis == "augmented" & transition == at,
            coarse_vs_augmented_agreement := sprintf("%s (vs %s)", lab, ct)]
}

setcolorder(trans_tbl, c(
  "stage_axis","transition","transition_label","n_set_degs",
  "non_estimable","low_power",
  "expr_hep_share","liana_hep_share","hotspot_hep_share",
  "n_axes_hepdom","n_axes_estimable","verdict",
  "coarse_vs_augmented_agreement",
  "dominant_carriers","top_nonhep",
  "n_hep_intrinsic","n_nonhep_intrinsic",
  "liana_hep_genes","liana_tot_genes",
  "hotspot_hep_captured","hotspot_tot_captured"))
setorder(trans_tbl, stage_axis, transition)

fwrite(trans_tbl, file.path(DIR, "stagedeg_verdict_by_transition.csv"))

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" STAGE-DEG CELL-TYPE ROUTING VERDICT (3-axis synthesis)\n")
cat("=============================================================\n\n")

pf <- function(x) ifelse(is.na(x), "  NA ", sprintf("%.2f", x))
cat(sprintf("%-10s %-26s %6s %6s %6s  %-3s  %s\n",
            "axis","transition","exprHep","lianH","hotsH","#hd","verdict"))
cat(strrep("-", 92), "\n")
for (i in seq_len(nrow(trans_tbl))) {
  r <- trans_tbl[i]
  cat(sprintf("%-10s %-26s %6s %6s %6s  %-3s  %s\n",
              r$stage_axis, r$transition_label,
              pf(r$expr_hep_share), pf(r$liana_hep_share), pf(r$hotspot_hep_share),
              r$n_axes_hepdom, r$verdict))
}

cat("\n--- Coarse-vs-augmented agreement (robustness) ---\n")
print(agree_dt)

n_multi <- trans_tbl[!grepl("hepatocyte_dominant", verdict), .N]
n_total <- trans_tbl[, .N]
cat(sprintf("\nVerdict tally: %d / %d transition-rows = multi_cellular; %d = hepatocyte_dominant.\n",
            n_multi, n_total, n_total - n_multi))
cat(sprintf("Coarse/augmented agreement: %d / %d paired analogs agree.\n",
            sum(agree_dt$agree, na.rm = TRUE), sum(!is.na(agree_dt$agree))))

cat("\nWrote:\n  ", file.path(DIR, "stagedeg_celltype_verdict.csv"),
    "\n  ", file.path(DIR, "stagedeg_verdict_by_transition.csv"), "\n")
