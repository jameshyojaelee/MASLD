#!/usr/bin/env Rscript
# 08: the single-nucleus atlas, at donor level.
#
# All-cells pseudobulk is the primary column, for comparability with the other
# assays. Per-lineage pseudobulk is a within-lineage donor-level arm: it removes
# confounding by broad lineage abundance but can retain subtype, zonation,
# state and technical composition differences, so it is not called
# "composition-free". Its tests form one family over modules x lineages and are
# corrected as such.
#
# ENDPOINT. The coarse disease stage (Healthy < Steatosis < Steatohepatitis <
# Cirrhosis) is a mixed severity category, not a fibrosis grade; it is mapped
# to the fibrosis family for the panel and labelled as mixed in every row.
#
# BIOLOGICAL UNIT. Sequencing runs are not donors. Runs are collapsed to donors
# through the recorded pairing tables before any model is fitted.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages(library(edgeR))
source(file.path(CAM_PROJECT_ROOT, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))

contract <- cam_contract()
out <- cam_dir("assays")
cam_assert(file.exists(file.path(cam_out_root(), "discovery", "READY")), "Run 05 first")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)
n_family <- length(module_ids)

sn_meta <- fread(cam_input("snrna_metadata", contract))
pb_dir <- cam_input("snrna_pseudobulk_dir", contract)
srr2donor <- build_srr_to_donor_map(CAM_PROJECT_ROOT)
cam_say("donor pairing map covers ", length(srr2donor), " runs")

stage_levels <- c(Healthy = 0, Steatosis = 1, Steatohepatitis = 2, Cirrhosis = 3)

collapse_runs <- function(mat) {
  donor <- ifelse(colnames(mat) %in% names(srr2donor),
                  srr2donor[colnames(mat)], colnames(mat))
  t(rowsum(t(mat), group = donor))
}

score_block <- function(mat, tag) {
  if (ncol(mat) < 12L) return(NULL)
  logcpm <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  z <- cam_zscore_rows(logcpm)
  measured <- rownames(z)
  run_meta <- unique(sn_meta[, .(sample, dataset, disease_stage_coarse)])
  run_meta[, donor := ifelse(sample %in% names(srr2donor), srr2donor[sample], sample)]
  donor_meta <- run_meta[, .(
    dataset = names(sort(table(dataset), decreasing = TRUE))[1],
    stage = names(sort(table(disease_stage_coarse), decreasing = TRUE))[1]
  ), by = donor]
  dm <- donor_meta[match(colnames(z), donor)]
  y <- cam_endpoint_z(unname(stage_levels[dm$stage]), data.frame(dataset = dm$dataset))
  rows <- list()
  for (i in seq_along(module_ids)) {
    tb <- cam_testability(modules[[i]], measured, contract)
    idx <- intersect(modules[[i]], measured)
    score <- if (length(idx)) standardize_vector(colMeans(z[idx, , drop = FALSE], na.rm = TRUE))
             else rep(NA_real_, ncol(z))
    d <- data.table(score = score, y = y, dataset = dm$dataset)
    fit <- if (tb$testable) cam_fit_score(d, "score", "y", "dataset") else cam_fit_empty()
    rows[[length(rows) + 1L]] <- data.table(
      assay = tag, module_id = module_ids[i],
      endpoint = "disease_stage_coarse_ordinal", endpoint_family = "fibrosis",
      endpoint_note = "coarse diagnosis category Healthy<Steatosis<Steatohepatitis<Cirrhosis, mixed severity, z-scored",
      unit = "biological_donor", n_units = fit$n, df = fit$df,
      testable = tb$testable, n_members = tb$n_members,
      n_measured = tb$n_measured, fraction_measured = tb$fraction_measured,
      beta = fit$beta, se = fit$se, p_two_sided = fit$p_two_sided)
  }
  list(rows = rbindlist(rows), n_donors = ncol(z), n_with_stage = sum(is.finite(y)))
}

results <- list(); census <- list()

all_cells <- readRDS(file.path(cam_out_root(), "universe", "snrna_all_cells_pseudobulk.rds"))
donor_all <- collapse_runs(all_cells)
cam_say("all-cells: ", ncol(all_cells), " runs -> ", ncol(donor_all), " donors")
cam_assert(ncol(donor_all) < ncol(all_cells),
           "Donor collapse did not reduce the run count; the pairing tables were not applied")
blk <- score_block(donor_all, "snrna_all_cells")
results[[length(results) + 1L]] <- blk$rows
census[[length(census) + 1L]] <- data.table(
  block = "all_cells", n_runs = ncol(all_cells), n_donors = blk$n_donors,
  n_donors_with_stage = blk$n_with_stage)

for (f in list.files(pb_dir, pattern = "^pseudobulk_.*\\.tsv\\.gz$", full.names = TRUE)) {
  lineage <- sub("^pseudobulk_", "", sub("\\.tsv\\.gz$", "", basename(f)))
  x <- fread(f); setnames(x, 1L, "gene")
  m <- as.matrix(x[, -1L]); rownames(m) <- x$gene
  dm <- collapse_runs(m)
  b <- score_block(dm, paste0("snrna_lineage_", lineage))
  if (is.null(b)) { cam_say("  ", lineage, ": ", ncol(dm), " donors, too few to model"); next }
  results[[length(results) + 1L]] <- b$rows
  census[[length(census) + 1L]] <- data.table(
    block = lineage, n_runs = ncol(m), n_donors = b$n_donors,
    n_donors_with_stage = b$n_with_stage)
  cam_say("  ", lineage, ": ", ncol(m), " runs -> ", b$n_donors, " donors")
}

res <- rbindlist(results)
# The per-lineage arm is ONE family: modules x lineages, corrected together.
lin <- res$assay != "snrna_all_cells"
res[, lineage_family_q := NA_real_]
res[lin, lineage_family_q := ml_complete_bh(p_two_sided, sum(lin))]
cam_assert_no_prohibited_columns(res, contract)
cam_write_tsv(res, file.path(out, "snrna_results.tsv"))
cam_write_tsv(rbindlist(census), file.path(out, "snrna_donor_census.tsv"))
cam_write_json(list(
  n_blocks = length(unique(res$assay)),
  all_cells_donors = census[[1]]$n_donors,
  modules_testable_all_cells = sum(res[assay == "snrna_all_cells", testable]),
  lineage_blocks = setdiff(unique(res$assay), "snrna_all_cells"),
  lineage_family_size = sum(lin),
  lineage_family_q05 = sum(res$lineage_family_q < 0.05, na.rm = TRUE),
  direction_applied_at_scoring = FALSE,
  note = "per-lineage arm is within-lineage donor-level evidence, corrected over modules x lineages"
), file.path(out, "snrna_summary.json"))
writeLines("snrna scored", file.path(out, "READY_snrna"))
cam_say("08 complete")
