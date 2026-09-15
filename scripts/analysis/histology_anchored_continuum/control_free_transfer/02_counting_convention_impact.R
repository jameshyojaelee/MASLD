#!/usr/bin/env Rscript
# CFT-v1b: does the read-vs-fragment counting convention move the anchor?
# merged_dge (the substrate the frozen axis was scored on) is FRAGMENT-counted.
# The current gene_counts.txt files are READ-counted for every -p cohort (ratio ~0.50).
# Arm 1: re-score the released paired-end cohorts from READ counts, compare anchors.
# Arm 2: score the new control-free cohorts from floor(read/2), a fragment surrogate.
suppressPackageStartupMessages({ library(data.table); library(limma) })

OUT <- Sys.getenv("CFT_OUT_DIR"); stopifnot(nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ROOT     <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RELEASE  <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z")
DGE      <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/arms/F_five/results/integration/merged_dge.rds")
MANIFEST <- file.path(ROOT, "figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/stage_extensions/five_cohort_sample_manifest.tsv")
COHORT_RES <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results")
META_ROOT  <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/pipelines/custom")

abort <- function(...) stop(paste0("CFT1B-ABORT: ", sprintf(...)), call. = FALSE)
chk <- function(cond, ...) if (!isTRUE(cond)) abort(...)
base_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))
sp_rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }

load_tab <- fread(file.path(RELEASE, "projection", "fixed_projection_loadings.tsv"))
load_tab[, gene_id_base := base_id(gsub('"', "", gene_id_base))]
sig139 <- load_tab$gene_id_base[load_tab$used_for_fixed_projection %in% c(TRUE, "TRUE")]
chk(length(sig139) == 139L, "frozen signature is not 139 genes")
model <- readRDS(file.path(RELEASE, "reproduction", "discovery_signature_model.rds"))
loading <- model$discovery_loading_oriented; names(loading) <- base_id(names(loading))
centre  <- model$discovery_center;            names(centre)  <- base_id(names(centre))

dge <- readRDS(DGE)
universe <- rownames(dge$counts)

project_set <- function(mat, g, ld, ct)
  as.numeric(crossprod(ld[g], sweep(mat[g, , drop = FALSE], 1L, ct[g], FUN = "-")))
normalise_cohort <- function(counts, id) {
  keep <- rowSums(counts) > ncol(counts)
  chk(sum(keep) >= 1000L, "%s: fewer than 1000 genes survive the row-sum filter", id)
  le <- log2(counts[keep, , drop = FALSE] + 1)
  nz <- normalizeQuantiles(le); dimnames(nz) <- dimnames(le); nz
}
read_fc <- function(cohort) {
  f <- file.path(COHORT_RES, cohort, "counts/featurecounts/gene_counts.txt")
  tb <- fread(f, skip = 1L, sep = "\t", header = TRUE)
  m <- as.matrix(tb[, 7:ncol(tb), with = FALSE])
  colnames(m) <- sub("\\.Aligned.*$", "", basename(colnames(m)))
  rownames(m) <- tb[[1]]
  storage.mode(m) <- "integer"; m
}
score_from <- function(cm, id) {
  rownames(cm) <- base_id(rownames(cm))
  cm <- cm[!duplicated(rownames(cm)), , drop = FALSE]
  X <- normalise_cohort(cm, id)
  g <- sig139[sig139 %in% rownames(X)]
  list(score = setNames(project_set(X, g, loading, centre), colnames(X)), n_sig = length(g))
}

## ---- strata for the released cohorts, exactly as the sealed gate defined them
man <- fread(MANIFEST)
man[, fibrosis_stage := suppressWarnings(as.numeric(fibrosis_stage))]
proj <- fread(file.path(RELEASE, "projection", "participant_scores.tsv"))
for (cc in c("sample_id", "dataset", "axis_id")) proj[[cc]] <- gsub('"', "", proj[[cc]])
rel <- proj[axis_id == "fixed_projection", .(sample_id, dataset, released = as.numeric(axis_raw))]

rows <- list(); persample <- list()
for (cid in c("GSE213621", "GSE130970", "GSE135251", "GSE162694")) {
  say("[arm1] %s", cid)
  fc <- read_fc(cid)
  smp <- rownames(dge$samples)[dge$samples$dataset == cid]
  g <- intersect(universe, rownames(fc))
  R <- score_from(fc[g, smp, drop = FALSE], paste0(cid, "_read"))
  H <- score_from(matrix(as.integer(floor(fc[g, smp, drop = FALSE] / 2L)),
                         nrow = length(g), dimnames = list(g, smp)), paste0(cid, "_half"))
  d <- merge(man[dataset == cid], rel[dataset == cid], by = c("sample_id", "dataset"))
  chk(nrow(d) > 0, "%s: zero-row manifest join -- RAISING", cid)
  d[, score_read_counted := R$score[sample_id]]
  d[, score_half_counted := H$score[sample_id]]
  persample[[cid]] <- d
  for (st in c("S2_within_disease", "S3_drop_stage_zero")) {
    dd <- d[is.finite(fibrosis_stage)]
    dd <- if (st == "S2_within_disease")
      dd[is.na(diagnosis_harmonized) | diagnosis_harmonized != "Control"] else dd[fibrosis_stage != 0]
    if (nrow(dd) < 10L) next
    rows[[length(rows) + 1L]] <- data.table(
      cohort = cid, stratum = st, n = nrow(dd),
      rho_released_fragment_counted = sp_rho(dd$released, dd$fibrosis_stage),
      rho_recomputed_read_counted   = sp_rho(dd$score_read_counted, dd$fibrosis_stage),
      rho_recomputed_half_counted   = sp_rho(dd$score_half_counted, dd$fibrosis_stage),
      score_spearman_released_vs_read = sp_rho(dd$released, dd$score_read_counted),
      score_spearman_released_vs_half = sp_rho(dd$released, dd$score_half_counted))
  }
  rm(fc, R, H); gc()
}
arm1 <- rbindlist(rows)
arm1[, delta_rho_read_minus_released := rho_recomputed_read_counted - rho_released_fragment_counted]
arm1[, delta_rho_half_minus_released := rho_recomputed_half_counted - rho_released_fragment_counted]
fwrite(arm1, file.path(OUT, "counting_convention_anchor_impact_released_cohorts.tsv"), sep = "\t")
fwrite(rbindlist(persample, fill = TRUE),
       file.path(OUT, "counting_convention_per_sample_released_cohorts.tsv"), sep = "\t")

## ---- arm 2: new control-free cohorts, read counts vs fragment surrogate ----
meta <- list(
  GSE174478 = { d <- fread(file.path(META_ROOT, "GSE174478/metadata/SraRunTable.csv"))
                data.table(sample_id = d$Run, fibrosis_stage = as.numeric(d$fibrosis_stage),
                           nas_score = as.numeric(d$nas_score)) },
  GSE240729 = { d <- fread(file.path(META_ROOT, "GSE240729/metadata/SraRunTable.csv"))
                data.table(sample_id = d$Run,
                           fibrosis_stage = as.numeric(sub("^F", "", d$fibrosisscore)),
                           nas_score = NA_real_) })
rows2 <- list(); ps2 <- list()
for (cid in names(meta)) {
  say("[arm2] %s", cid)
  fc <- read_fc(cid)
  g <- intersect(universe, rownames(fc))
  chk(length(g) > 10000L, "%s: only %d genes shared with the merged_dge substrate", cid, length(g))
  R <- score_from(fc[g, , drop = FALSE], paste0(cid, "_read"))
  H <- score_from(matrix(as.integer(floor(fc[g, , drop = FALSE] / 2L)),
                         nrow = length(g), dimnames = list(g, colnames(fc))), paste0(cid, "_half"))
  d <- meta[[cid]]
  common <- intersect(d$sample_id, names(R$score))
  chk(length(common) > 0, "%s: zero-row join between metadata and counts -- RAISING", cid)
  d <- d[match(common, sample_id)]
  d[, score_read_counted := R$score[common]]
  d[, score_half_counted := H$score[common]]
  d[, dataset := cid]
  ps2[[cid]] <- d
  for (st in c("S1_all_sample", "S3_drop_stage_zero")) {
    dd <- if (st == "S1_all_sample") d[is.finite(fibrosis_stage)] else d[is.finite(fibrosis_stage) & fibrosis_stage != 0]
    rows2[[length(rows2) + 1L]] <- data.table(
      cohort = cid, stratum = st, outcome = "fibrosis_stage", n = nrow(dd),
      rho_read_counted = sp_rho(dd$score_read_counted, dd$fibrosis_stage),
      rho_half_counted = sp_rho(dd$score_half_counted, dd$fibrosis_stage),
      score_agreement_spearman = sp_rho(dd$score_read_counted, dd$score_half_counted))
  }
  if (cid == "GSE174478") {
    dd <- d[is.finite(nas_score)]
    rows2[[length(rows2) + 1L]] <- data.table(
      cohort = cid, stratum = "S1_all_sample", outcome = "nas_score", n = nrow(dd),
      rho_read_counted = sp_rho(dd$score_read_counted, dd$nas_score),
      rho_half_counted = sp_rho(dd$score_half_counted, dd$nas_score),
      score_agreement_spearman = sp_rho(dd$score_read_counted, dd$score_half_counted))
  }
  rm(fc, R, H); gc()
}
arm2 <- rbindlist(rows2)
arm2[, delta_rho_half_minus_read := rho_half_counted - rho_read_counted]
fwrite(arm2, file.path(OUT, "counting_convention_anchor_impact_new_cohorts.tsv"), sep = "\t")
fwrite(rbindlist(ps2, fill = TRUE),
       file.path(OUT, "counting_convention_per_sample_new_cohorts.tsv"), sep = "\t")

writeLines(c("ARM 1 released cohorts:", capture.output(print(arm1)),
             "", "ARM 2 new control-free cohorts:", capture.output(print(arm2))),
           file.path(OUT, "COUNTING_CONVENTION_IMPACT.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_counting_convention.txt"))
say("DONE")
