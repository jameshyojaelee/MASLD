#!/usr/bin/env Rscript
# ============================================================================
# ccc_v3_dc_rebuild.R  — DONOR-COLLAPSED rebuild of the Fig3H CCC L-R heatmap
#
# Pseudoreplication fix (2026-07-12). The live Fig3H panel
# (fig3h_ccc_lr_heatmap.pdf) is built from a RUN-LEVEL v2/v3 CCC chain in which
# each sequencing run (SRR/GSM) was treated as an independent donor. This
# script rebuilds the panel + trajectory table on TRUE-donor-collapsed LIANA
# scores (06 donor-collapse fix) and re-fits the coarse-axis stage LMM so the
# per-donor statistics use biological donors, not runs. Everything is written
# to NEW *_dc paths; the canonical run-level artifacts are left untouched for
# before/after review (team lead promotes after).
#
# BEFORE/AFTER = a clean A/B: the coarse-axis stage LMM is REFIT with ONE model
# spec (07b's dataset_pooled random effect + log10 cell-count covariates; the
# two-level BH-within-ct x n_families family-q) on (i) the run-level scores and
# (ii) the donor-collapsed scores. Only the donor grouping changes, so the
# p/Estimate deltas isolate the pseudoreplication effect. age/sex covariates are
# omitted (scRNA metadata is NA-dominated for them); the published headline
# before-numbers are reported alongside as reference.
#
# Inputs:
#   run-level : all_donor_lr_scores_v2.tsv.gz     + donor_metadata_v2.tsv
#   collapsed : all_donor_lr_scores_v2_dc.tsv.gz  + donor_metadata_v2_dc.tsv
#   headline  : stage_lr_headline_v3.tsv  (READ-ONLY; 13-pair selection FIXED)
#
# Outputs:
#   - stage_trajectory_v2/stage_lr_lmm_coarse_v2_dc.tsv        (DC family-wide LMM)
#   - stage_trajectory_v2/stage_lr_lmm_coarse_v2_runrefit.tsv  (run-level refit)
#   - figures/main/fig3_RNAseq/panels/data/ccc_trajectories_data_dc.csv
#   - figures/main/fig3_RNAseq/panels/data/fig3h_liana_matrix_dc.csv
#   - figures/main/fig3_RNAseq/panels/fig3h_ccc_lr_heatmap_dc.pdf
#   - stage_trajectory_v2/ccc_dc_before_after_report.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
options(mc.cores = N_CORES)
cat(sprintf("[parallel] %d cores\n", N_CORES))

V2_STAGE <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
MCP <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs")

LR_DC    <- file.path(V2_STAGE, "all_donor_lr_scores_v2_dc.tsv.gz")
META_DC  <- file.path(MCP, "donor_metadata_v2_dc.tsv")
LR_RUN   <- file.path(V2_STAGE, "all_donor_lr_scores_v2.tsv.gz")
META_RUN <- file.path(MCP, "donor_metadata_v2.tsv")

PANEL_DIR <- file.path(BASE, "figures/main/fig3_RNAseq/panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_TRAJ_DC <- file.path(DATA_DIR, "ccc_trajectories_data_dc.csv")
OUT_MAT_DC  <- file.path(DATA_DIR, "fig3h_liana_matrix_dc.csv")
OUT_FIG_DC  <- file.path(PANEL_DIR, "fig3h_ccc_lr_heatmap_dc.pdf")
OUT_LMM_DC  <- file.path(V2_STAGE, "stage_lr_lmm_coarse_v2_dc.tsv")
OUT_LMM_RUN <- file.path(V2_STAGE, "stage_lr_lmm_coarse_v2_runrefit.tsv")
OUT_REPORT  <- file.path(V2_STAGE, "ccc_dc_before_after_report.tsv")

STAGES       <- c("Healthy", "Steatosis", "Steatohepatitis")           # Fig3H axis
STAGE_LEVELS <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
MIN_DONORS_PER_LR      <- 30
MIN_DONORS_PER_STRATUM <- 10
TARGET_TERM  <- "disease_stage_coarseSteatohepatitis"
SMALL_COHORTS <- c("GSE174748", "GSE189600")   # 07b micro-cohort pooling

# ============================================================================
# helpers: build lr_long, fit coarse-axis LMM (07b model spec, minus age/sex)
# ============================================================================
build_lr_long <- function(lr_path, meta_path, lr_has_dataset) {
  lr <- fread(lr_path)
  meta <- fread(meta_path)
  if (lr_has_dataset) {
    lr <- merge(lr, meta[, .(sample, disease_stage_coarse)], by = "sample", all.x = FALSE)
  } else {
    lr <- merge(lr, meta[, .(sample, dataset, disease_stage_coarse)],
                by = "sample", all.x = FALSE)
  }
  lr[, dataset_pooled := ifelse(dataset %in% SMALL_COHORTS, "Other_small", as.character(dataset))]
  lr[, score := -log10(pmax(magnitude_rank, 1e-4))]
  lr[, disease_stage_coarse := factor(disease_stage_coarse, levels = STAGE_LEVELS)]
  lr <- lr[!is.na(disease_stage_coarse)]
  lr[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
  lr[, ct_pair := paste(source, target, sep = "->")]
  lr[]
}

fit_one <- function(d, axis_term) {
  if (nrow(d) < MIN_DONORS_PER_LR) return(NULL)
  n_pools <- length(unique(d$dataset_pooled))
  use_random <- n_pools >= 3
  fixed <- c(axis_term)
  if (!use_random && n_pools > 1) fixed <- c(fixed, "dataset_pooled")
  fixed <- c(fixed, "log10(pmax(n_source_cells,1))", "log10(pmax(n_target_cells,1))")
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    rhs <- paste(rhs, "+ (1 | dataset_pooled)")
    fit <- try(lmerTest::lmer(stats::as.formula(paste("score ~", rhs)),
                              data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)), data = d), silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co); hits <- rn[startsWith(rn, axis_term)]
  if (!length(hits)) return(NULL)
  out <- as.data.table(co[hits, , drop = FALSE], keep.rownames = "term")
  ncs <- ncol(out) - 1L
  if (ncs == 5L) setnames(out, c("term","Estimate","StdErr","df","tval","pval"))
  else if (ncs == 4L) { setnames(out, c("term","Estimate","StdErr","tval","pval")); out[, df := NA_real_] }
  else return(NULL)
  out[, n_donors := nrow(d)]
  out
}

fit_coarse <- function(lr_long, tag) {
  splits <- split(lr_long, by = c("ct_pair", "lr_pair"), drop = TRUE)
  cat(sprintf("[lmm:%s] %d (ct_pair,lr_pair) groups\n", tag, length(splits)))
  wrapper <- function(d) {
    tab <- table(d$disease_stage_coarse)
    if (any(tab < MIN_DONORS_PER_STRATUM)) return(NULL)
    fr <- fit_one(d, "disease_stage_coarse")
    if (is.null(fr)) return(NULL)
    fr[, ct_pair := d$ct_pair[1]]; fr[, lr_pair := d$lr_pair[1]]
    fr[, ligand_complex := d$ligand_complex[1]]; fr[, receptor_complex := d$receptor_complex[1]]
    fr
  }
  t0 <- Sys.time()
  rl <- if (N_CORES > 1) parallel::mclapply(splits, wrapper, mc.cores = N_CORES) else lapply(splits, wrapper)
  cat(sprintf("[lmm:%s] took %.1f min\n", tag, as.numeric(Sys.time() - t0, units = "mins")))
  out <- rbindlist(rl[!sapply(rl, is.null)], fill = TRUE)
  if (nrow(out) == 0) return(out)
  # 07b two-level family-q: BH within (ct_pair, term) x n_families (distinct ct_pairs)
  out[, q_within_ct := stats::p.adjust(pval, method = "BH"), by = .(ct_pair, term)]
  nf <- length(unique(out$ct_pair))
  out[, n_families := nf]
  out[, q_bonferroni_family := pmin(q_within_ct * nf, 1)]
  out
}

# ============================================================================
# 1. Fit coarse LMM on RUN-LEVEL and DONOR-COLLAPSED data (identical spec)
# ============================================================================
cat("\n========== RUN-LEVEL refit (before) ==========\n")
lr_run <- build_lr_long(LR_RUN, META_RUN, lr_has_dataset = FALSE)
cat(sprintf("[run] LR rows=%d donors=%d ; per-stage donors:\n", nrow(lr_run), uniqueN(lr_run$sample)))
print(lr_run[, .(n = uniqueN(sample)), by = disease_stage_coarse])
res_run <- fit_coarse(lr_run, "run")
fwrite(res_run, OUT_LMM_RUN, sep = "\t")
sh_run <- res_run[term == TARGET_TERM]

cat("\n========== DONOR-COLLAPSED fit (after) ==========\n")
lr_dc <- build_lr_long(LR_DC, META_DC, lr_has_dataset = TRUE)
cat(sprintf("[dc] LR rows=%d donors=%d ; per-stage donors:\n", nrow(lr_dc), uniqueN(lr_dc$sample)))
print(lr_dc[, .(n = uniqueN(sample)), by = disease_stage_coarse])
res_dc <- fit_coarse(lr_dc, "dc")
fwrite(res_dc, OUT_LMM_DC, sep = "\t")
sh_dc <- res_dc[term == TARGET_TERM]
cat(sprintf("[lmm] SH-term testable pairs: run=%d  dc=%d\n", nrow(sh_run), nrow(sh_dc)))

# ============================================================================
# 2. HEADLINE SELECTION (held fixed from run-level v3 headline; ccc_v3_panels.R)
# ============================================================================
hl <- fread(file.path(V3_DIR, "stage_lr_headline_v3.tsv"))
hl[, sender_ct := tstrsplit(ct_pair, "->", fixed = TRUE)[[1]]]
hl[, receiver_ct := tstrsplit(ct_pair, "->", fixed = TRUE)[[2]]]
para <- hl[sender_ct != receiver_ct & gate_significance == TRUE & gate_bootstrap == TRUE &
           gate_permutation == TRUE & (gate_leverage == TRUE | is.na(gate_leverage))]
para[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
ct_short_map <- c("Hepatocytes"="Hep","Endothelial cells"="Endo","Fibroblasts"="Fib",
                  "Macrophages"="Mac","Cholangiocytes"="Chol")
para <- para[sender_ct != "T cells" & receiver_ct != "T cells"]
para[, c("sender","receiver") := tstrsplit(ct_pair, "->", fixed = TRUE)]
para[, sender_short := ct_short_map[sender]]
para[, receiver_short := ct_short_map[receiver]]
para[, lr_label := sub("__", "→", lr_pair)]
para[, eff_estimate := fcoalesce(coarse_Estimate_SH, cont_Estimate, aug_Estimate, doc_Estimate)]
para[, q_chord := fcoalesce(q_min_tippett, q_min_bonferroni)]
para[, ct_pair_label := paste(sender_short, receiver_short, sep = "→")]
para[, headline_label := paste0(lr_label, "  (", ct_pair_label, ")")]
para[, q_label := sprintf("q=%.1e", q_chord)]
setorder(para, q_chord)
cat(sprintf("\n[headline] %d paracrine headline pairs held fixed\n", nrow(para)))

# ============================================================================
# 3. TRAJECTORY recompute on DC data (ccc_v3_panels.R B_summary logic)
# ============================================================================
pair_keys <- para[, .(ct_pair, lr_pair, headline_label, q_label)]
B_data <- merge(lr_dc[, .(sample, ct_pair, lr_pair, disease_stage_coarse, score)],
                pair_keys, by = c("ct_pair", "lr_pair"))
B_data <- B_data[!is.na(disease_stage_coarse)]
B_summary <- B_data[, .(mean_score = mean(score, na.rm = TRUE),
                        se = sd(score, na.rm = TRUE) / sqrt(.N),
                        n = .N), by = .(headline_label, q_label, disease_stage_coarse)]
B_summary[, lo := mean_score - 1.96 * se]
B_summary[, hi := mean_score + 1.96 * se]
B_summary[, is_nampt := grepl("NAMPT", headline_label)]
setorder(B_summary, headline_label, disease_stage_coarse)
fwrite(B_summary, OUT_TRAJ_DC)
cat(sprintf("[output] %s\n", OUT_TRAJ_DC))
cat("[traj] donor N per stage among headline pairs (post-collapse):\n")
print(B_summary[disease_stage_coarse %in% STAGES,
                .(median_n = as.integer(median(n)), min_n = min(n), max_n = max(n)),
                by = disease_stage_coarse])

# ============================================================================
# 4. Fig3H heatmap on DC data (replicate singlecell_module_heatmap.R block)
# ============================================================================
G6 <- function(...) gpar(fontsize = 6, fontfamily = "Helvetica", ...)
lia <- B_summary[disease_stage_coarse %in% STAGES]
LW  <- dcast(lia, headline_label ~ factor(disease_stage_coarse, levels = STAGES),
             value.var = "mean_score")
lr_lab <- LW$headline_label
lr_lab <- gsub("ITGA(\\w+)_ITGB(\\w+)", "ITGA\\1/B\\2", lr_lab)
lr_lab <- gsub("_", "/", lr_lab)
lr_lab <- gsub(" +\\(", " (", lr_lab)
lmat <- t(scale(t(as.matrix(LW[, ..STAGES])))); rownames(lmat) <- lr_lab
fwrite(data.table(headline_label = lr_lab, lmat), OUT_MAT_DC)
na_cells <- sum(is.na(lmat))
if (na_cells > 0) cat(sprintf("[fig3h] WARNING %d NA z-cells (stage with no donors post-collapse)\n", na_cells))
lmat <- lmat[order(-(lmat[, "Steatohepatitis"] - lmat[, "Healthy"])), , drop = FALSE]
lmat <- t(lmat)
ZCOL <- colorRamp2(c(-1.5, 0, 1.5), c("#1565C0", "white", "#C9265E"))
ht_l <- Heatmap(lmat, name = "L-R z", col = ZCOL, na_col = "grey92",
  cluster_rows = FALSE, cluster_columns = FALSE, row_title = NULL,
  row_names_side = "left", row_names_gp = G6(),
  column_names_gp = G6(), column_names_rot = 45, column_names_side = "bottom",
  column_title = NULL,
  height = unit(3 * 0.50, "cm"), width = unit(ncol(lmat) * 0.42, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(title_gp = G6(), labels_gp = G6(), grid_width = unit(2.5, "mm")))
cairo_pdf(OUT_FIG_DC, width = 3.54, height = 1.95)
draw(ht_l, heatmap_legend_side = "right", padding = unit(c(2, 1, 2, 1), "mm"))
dev.off()
cat(sprintf("[output] %s  (%d L-R pairs)\n", OUT_FIG_DC, ncol(lmat)))

# ============================================================================
# 5. BEFORE / AFTER report
# ============================================================================
# (a) Fig3H z-matrix concordance before vs after (aligned on tidy label)
before_mat <- fread(file.path(DATA_DIR, "fig3h_liana_matrix.csv"))
after_mat  <- fread(OUT_MAT_DC)
mm <- merge(before_mat, after_mat, by = "headline_label", suffixes = c("_before", "_after"))
vb <- as.matrix(mm[, .(Healthy_before, Steatosis_before, Steatohepatitis_before)])
va <- as.matrix(mm[, .(Healthy_after, Steatosis_after, Steatohepatitis_after)])
ok <- is.finite(vb) & is.finite(va)
r_cell <- cor(as.vector(vb)[as.vector(ok)], as.vector(va)[as.vector(ok)])
cat(sprintf("\n[report] Fig3H z-matrix Pearson r(before,after)=%.3f over %d pairs x 3 stages (n_cells=%d, NA_after=%d)\n",
            r_cell, nrow(mm), sum(ok), sum(!is.finite(va))))

# (b) per-pair coarse SH stat: run-refit vs DC-refit (clean A/B) + published ref
pub <- para[, .(ct_pair, lr_pair, headline_label,
                pub_Estimate = coarse_Estimate_SH, pub_pval = coarse_pval_SH,
                pub_q = coarse_q_bonferroni_family_SH, pub_n = coarse_n_donors_SH)]
runst <- sh_run[, .(ct_pair, lr_pair, run_Estimate = Estimate, run_pval = pval,
                    run_q = q_bonferroni_family, run_n = n_donors)]
dcst  <- sh_dc[, .(ct_pair, lr_pair, dc_Estimate = Estimate, dc_pval = pval,
                   dc_q = q_bonferroni_family, dc_n = n_donors)]
rep_dt <- Reduce(function(a,b) merge(a,b, by=c("ct_pair","lr_pair"), all.x=TRUE),
                 list(pub, runst, dcst))
rep_dt[, testable_dc := !is.na(dc_pval)]
rep_dt[, sign_preserved := !is.na(dc_Estimate) & !is.na(run_Estimate) &
                            sign(dc_Estimate) == sign(run_Estimate)]
rep_dt[, dc_p05 := !is.na(dc_pval) & dc_pval < 0.05]
rep_dt[, dc_q05 := !is.na(dc_q) & dc_q < 0.05]
rep_dt[, run_p05 := !is.na(run_pval) & run_pval < 0.05]
rep_dt[, run_q05 := !is.na(run_q) & run_q < 0.05]
setorder(rep_dt, pub_pval)
fwrite(rep_dt, OUT_REPORT, sep = "\t")

cat("\n================= BEFORE/AFTER (13 headline pairs) =================\n")
cat("cols: pub=published headline(07b full model) | run=refit run-level | dc=refit donor-collapsed\n")
print(rep_dt[, .(headline_label,
                 runE=round(run_Estimate,2), runp=signif(run_pval,2), runq=signif(run_q,2), runN=run_n,
                 dcE=round(dc_Estimate,2), dcp=signif(dc_pval,2), dcq=signif(dc_q,2), dcN=dc_n,
                 sgn=sign_preserved, dc_p05, dc_q05)])
cat(sprintf("\n[verdict] 13 headline pairs:\n"))
cat(sprintf("  run-refit  coarse p<0.05: %d/13 ; family-q<0.05: %d/13\n",
            sum(rep_dt$run_p05, na.rm=TRUE), sum(rep_dt$run_q05, na.rm=TRUE)))
cat(sprintf("  DC-refit   testable: %d/13 ; sign_preserved: %d/13 ; p<0.05: %d/13 ; family-q<0.05: %d/13\n",
            sum(rep_dt$testable_dc), sum(rep_dt$sign_preserved, na.rm=TRUE),
            sum(rep_dt$dc_p05, na.rm=TRUE), sum(rep_dt$dc_q05, na.rm=TRUE)))
cat(sprintf("  published  family-q<0.05: %d/13 (reference; 07b full model incl sex where available)\n",
            sum(pub$pub_q < 0.05, na.rm=TRUE)))
cat("\n[done] ccc_v3_dc_rebuild\n")
