#!/usr/bin/env Rscript
# 01: build the three expression matrices every ordering arm consumes, audit what
# ComBat does to real disease signal, and prove the rebuilt pipeline reproduces
# known values from the validated substrate.
#
# The three matrices differ only in preprocessing:
#   E_raw  voom log-CPM, no cross-cohort transform
#   E_qn   quantile-normalised across samples          (Kamzolas step 1)
#   E_cb   ComBat(batch = cohort, mod = ~ sex) on E_qn (Kamzolas step 2)
#
# Arm A1 uses E_cb because that is what the paper did. Arm A2 uses E_qn because
# cohort is confounded with disease in this substrate in a way it was not in
# theirs, so ComBat can remove real signal. The audit below measures that rather
# than assuming it either way.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(sva)
})
# NOT preprocessCore. Its normalize.quantiles() is compiled against pthreads and
# probes the node's core count rather than the cgroup allocation, so under SLURM
# it dies with "return code from pthread_create() is 22" (EINVAL). limma's pure-R
# normalizeQuantiles() computes the same thing, is already a dependency, and has
# no threading. Verified equivalent in the smoke test.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
assert_true(!dir.exists(file.path(OUT, "arms")), paste0("Output already exists: ", OUT))
dir.create(file.path(OUT, "arms"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
set.seed(pre$seeds$master)
stamp_exemption(OUT)

log_step("reading substrate")
dge <- readRDS(substrate_path("dge"))
meta <- load_manifest()

assert_true(identical(dim(dge), c(23370L, 844L)),
            sprintf("Expected a 23370 x 844 DGEList; found %d x %d", nrow(dge), ncol(dge)))
assert_true(all(meta$sample_id %in% colnames(dge)),
            "A manifest participant is absent from the DGEList")
counts <- dge$counts
assert_true(all(counts == floor(counts)), "Counts are not integer; this is not the raw fragment matrix")

dge <- dge[, meta$sample_id]
assert_true(identical(colnames(dge), meta$sample_id), "Column order drifted from the manifest")

write_input_manifest(list(
  merged_dge = substrate_path("dge"),
  five_cohort_manifest = substrate_path("manifest"),
  stage_extension_results = substrate_path("stage_results")
), OUT)

# --------------------------------------------------------------- expression --
# No gene filtering: keeping all 23,370 makes every BH family in this workstream
# identical to stage_extension_all_gene_results.tsv, so gene counts are directly
# comparable to the eleven existing contrasts.
design0 <- model.matrix(~ dataset + inferred_sex, data = meta)

# voomWithQualityWeights on 23,370 x 844 costs ~27 minutes, so it is cached and
# reusable. CAB_REUSE_ROOT points at an earlier run's arms/ directory; the matrix
# is re-validated against this run's manifest before it is trusted, and the new
# output root still gets its own copy so each root stays self-contained.
reuse <- Sys.getenv("CAB_REUSE_ROOT", "")
reuse_path <- if (nzchar(reuse)) file.path(reuse, "arms", "E_raw.rds") else ""
if (nzchar(reuse_path) && file.exists(reuse_path)) {
  log_step("reusing cached voom matrix from ", reuse_path)
  E_raw <- readRDS(reuse_path)
  assert_true(identical(dim(E_raw), c(23370L, 844L)), "Cached voom matrix has the wrong shape")
  assert_true(identical(colnames(E_raw), meta$sample_id),
              "Cached voom matrix columns do not match this run's manifest")
  assert_true(identical(rownames(E_raw), rownames(dge)),
              "Cached voom matrix genes do not match this run's DGEList")
} else {
  log_step("voomWithQualityWeights on the full 23,370-gene family (~27 min)")
  v0 <- voomWithQualityWeights(dge, design0, plot = FALSE)
  E_raw <- v0$E
  assert_true(identical(dim(E_raw), c(23370L, 844L)), "voom output changed shape")
}
# Persist immediately. The previous run lost 27 minutes of voom because nothing
# was written until all three matrices existed and the second one threw.
saveRDS(E_raw, file.path(OUT, "arms", "E_raw.rds"))
log_step("E_raw persisted")

log_step("quantile normalisation (limma, pure R)")
E_qn <- normalizeQuantiles(as.matrix(E_raw))
dimnames(E_qn) <- dimnames(E_raw)
saveRDS(E_qn, file.path(OUT, "arms", "E_qn.rds"))

log_step("ComBat: batch = cohort, mod = ~ inferred_sex (as published)")
E_cb <- ComBat(dat = E_qn, batch = as.character(meta$dataset),
               mod = model.matrix(~ inferred_sex, data = meta), par.prior = TRUE)
dimnames(E_cb) <- dimnames(E_raw)
saveRDS(E_cb, file.path(OUT, "arms", "E_cb.rds"))
log_step("all three matrices persisted")

# PCA/Slingshot input only. filterByExpr here does NOT change any BH family; it
# only decides which genes define the axis.
keep <- filterByExpr(dge, design = design0)
pca_genes <- rownames(dge)[keep]
write_tsv_once(data.table(gene_id = pca_genes), file.path(OUT, "arms", "pca_input_genes.tsv"))
log_step("PCA input genes: ", length(pca_genes), " of ", nrow(dge))

# ------------------------------------------------------ ComBat damage audit --
# The published pipeline was applied to cohorts that all spanned the disease
# spectrum. Here GSE126848 (53) is entirely unstaged control/obese and GSE213621
# (361) is entirely staged disease, so cohort partly IS disease. Removing cohort
# without disease in mod must therefore remove some disease signal. Measure it.
log_step("ComBat damage audit on F4 vs F0")
popC <- population(meta, "POP-C")
fib_meta <- popC[fibrosis_stage %in% c(0L, 4L)]
fib_meta[, grp := factor(ifelse(fibrosis_stage == 4L, "F4", "F0"), levels = c("F0", "F4"))]

fit_on <- function(mat, info) {
  d <- model.matrix(~ dataset + inferred_sex + grp, data = info)
  d <- d[, qr(d)$pivot[seq_len(qr(d)$rank)], drop = FALSE]
  fit <- eBayes(lmFit(mat[, info$sample_id], d))
  cf <- grep("^grpF4$", colnames(d), value = TRUE)
  assert_true(length(cf) == 1L, "F4 coefficient is not estimable")
  topTable(fit, coef = cf, number = Inf, sort.by = "none")$logFC
}
lfc_raw <- fit_on(E_raw, fib_meta)
lfc_cb <- fit_on(E_cb, fib_meta)
slope <- coef(lm(lfc_cb ~ 0 + lfc_raw))[[1]]

audit <- data.table(
  quantity = c("n_F0", "n_F4", "attenuation_slope_cb_on_raw",
               "pearson_r_cb_vs_raw", "median_abs_lfc_raw", "median_abs_lfc_cb"),
  value = c(sum(fib_meta$grp == "F0"), sum(fib_meta$grp == "F4"), slope,
            cor(lfc_cb, lfc_raw), median(abs(lfc_raw)), median(abs(lfc_cb)))
)

# Axis diagnostics comparable to their reported 9% PC1 variance.
pc_diag <- rbindlist(lapply(c("E_raw", "E_qn", "E_cb"), function(nm) {
  mat <- get(nm)[pca_genes, ]
  p <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pc1 <- p$x[, 1]
  stage <- meta$fibrosis_stage
  ok <- !is.na(stage)
  data.table(matrix_id = nm,
             pc1_pct_variance = 100 * p$sdev[1]^2 / sum(p$sdev^2),
             kruskal_eta2_pc1_on_cohort = {
               kw <- kruskal.test(pc1 ~ meta$dataset)
               (unname(kw$statistic) - length(unique(meta$dataset)) + 1) /
                 (length(pc1) - length(unique(meta$dataset)))
             },
             spearman_pc1_fibrosis = cor(pc1[ok], stage[ok], method = "spearman"))
}))
write_tsv_once(rbind(audit, data.table(quantity = "kamzolas_reported_pc1_pct_variance", value = 9)),
               file.path(OUT, "arms", "combat_damage_audit.tsv"))
write_tsv_once(pc_diag, file.path(OUT, "arms", "pc1_diagnostics.tsv"))

# ------------------------------------------------------- pipeline controls --
# If the rebuilt pipeline cannot reproduce known values from the validated
# substrate, nothing downstream of it is trustworthy.
log_step("pipeline controls against stage_extension_all_gene_results.tsv")
ref <- fread(substrate_path("stage_results"),
             select = c("gene_id_versioned", "FDR", "contrast"))
controls <- rbindlist(lapply(pre$pipeline_controls$must_reproduce, function(ctl) {
  q <- ctl$quantity
  if (grepl("BH family cardinality", q)) {
    obs <- ref[, .N, by = contrast][, unique(N)]
    obs <- if (length(obs) == 1L) obs else NA_integer_
  } else if (grepl("^NAS 1-2", q)) {
    obs <- ref[contrast == "NAS3_4_vs_NAS1_2" & FDR < 0.05, .N]
  } else {
    obs <- ref[contrast == "F4_vs_F0" & FDR < 0.05, .N]
  }
  data.table(quantity = q, expected = as.integer(ctl$expected),
             observed = as.integer(obs), pass = identical(as.integer(obs), as.integer(ctl$expected)))
}))
write_tsv_once(controls, file.path(OUT, "arms", "pipeline_controls.tsv"))
print(controls)
assert_true(all(controls$pass), "A pipeline control failed; the substrate or reference has drifted")

write_run_parameters(OUT, list(n_pca_genes = length(pca_genes),
                               combat_attenuation_slope = slope))
log_step("SUBSTRATE_COMPLETE: ", OUT)
