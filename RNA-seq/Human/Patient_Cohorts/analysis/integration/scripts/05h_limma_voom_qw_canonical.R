#!/usr/bin/env Rscript
# ===========================================================================
# 05h — production limma_voom_qw__C2 canonical DEG call (benchmark winner).
# Mirrors 05_dream's input loading EXACTLY (merged_dge.rds + yaml include_in_mega
# + meta_matched.rds), but fits the benchmark-winning method:
#   voomWithQualityWeights -> lmFit -> eBayes,  design ~ dataset + sex + group.
# Emits BOTH raw (logFC, padj) and ashr-shrunk (shrunk_logFC, lfsr) columns,
# matching the dream_results.csv two-tier schema.
#
# Write order (BG-004, 2026-08-09 — this script DOES write the canonical file):
#   1. always  -> results/integration/limma_voom_qw_C2_results.csv   (method-named)
#   2. always  -> results/integration/limma_voom_qw_C2_sanity.csv    (dream concordance)
#   3. always  -> results/integration/limma_voom_qw_C2_promotion_gate.csv
#   4. ONLY if every promotion-gate check passes, by atomic rename:
#              -> results/integration/canonical_deg_results.csv
# A blocked run exits 1 leaving the canonical file untouched. Override with
# CANONICAL_PROMOTE_FORCE='<reason>' (the reason is recorded in the gate file).
# The header previously claimed this script never overwrote a canonical file; it did,
# ~80 lines before any diagnostic ran.
# ===========================================================================
suppressMessages({library(edgeR); library(limma); library(ashr); library(data.table)})
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BASE <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts")
SOURCE_RDIR <- file.path(BASE, "analysis/integration/results/integration")
RUN_ROOT <- Sys.getenv("MASLD_RUN_ROOT", "")
REMEDIATION_MODE <- nzchar(RUN_ROOT)
if (REMEDIATION_MODE) {
  RUN_ROOT <- normalizePath(RUN_ROOT, mustWork = TRUE)
  if (!file.exists(file.path(RUN_ROOT, ".bg001_candidate_root"))) stop("Missing BG-001 candidate sentinel")
  RDIR <- file.path(RUN_ROOT, "results/integration")
  dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)
} else {
  RDIR <- SOURCE_RDIR
}
DGE_PATH <- Sys.getenv("MASLD_DGE_INPUT", file.path(RDIR, "merged_dge.rds"))
META_PATH <- Sys.getenv("MASLD_META_INPUT", file.path(RDIR, "meta_matched.rds"))
DEG_OUTPUT <- Sys.getenv("MASLD_DEG_OUTPUT", if (REMEDIATION_MODE) file.path(RDIR, "deg_results.csv") else file.path(RDIR, "canonical_deg_results.csv"))
if (REMEDIATION_MODE) {
  output_parent <- normalizePath(dirname(DEG_OUTPUT), mustWork = TRUE)
  if (!startsWith(paste0(output_parent, .Platform$file.sep), paste0(RUN_ROOT, .Platform$file.sep))) {
    stop("MASLD_DEG_OUTPUT escapes MASLD_RUN_ROOT")
  }
  if (basename(DEG_OUTPUT) == "canonical_deg_results.csv") stop("Candidate mode may not write canonical_deg_results.csv")
}

candidate_temp <- function(destination) {
  if (!REMEDIATION_MODE) stop("candidate_temp is remediation-only")
  # Sys.readlink() yields NA for an absent path and nzchar(NA) is TRUE by
  # default, so a bare nzchar() refused every normal (absent) destination.
  # keepNA = TRUE + isTRUE() keeps the symlink test while tolerating absence.
  if (file.exists(destination) ||
      isTRUE(nzchar(Sys.readlink(destination), keepNA = TRUE))) {
    stop("Refusing existing or symlinked candidate output: ", destination)
  }
  uuid <- system2(
    "python3",
    c("-c", shQuote("import uuid; print(uuid.uuid4())")),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(uuid, "status")
  if (!is.null(status) && status != 0L) stop("Failed to allocate candidate output UUID")
  file.path(dirname(destination), paste0(".", basename(destination), ".", uuid[[1L]], ".tmp"))
}

publish_candidate <- function(temporary, destination) {
  publisher <- file.path(
    ROOT,
    "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/publish_candidate_file.py"
  )
  output <- system2(
    "python3",
    c(shQuote(publisher), "--run-root", shQuote(RUN_ROOT),
      "--temporary", shQuote(temporary), "--destination", shQuote(destination)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Fail-closed candidate publication failed for ", destination, ": ", paste(output, collapse = " "))
  }
  invisible(destination)
}

candidate_write_table <- function(object, destination, sep = "\t", compress = FALSE) {
  temporary <- candidate_temp(destination)
  # R parses "wbx" as a TEXT connection (the trailing x defeats binary
  # detection), which breaks saveRDS and gzcon. "wxb" is binary AND
  # exclusive-create. gzcon() additionally refuses an exclusive-mode
  # connection ("can only use read- or write- binary connections"), so the
  # compressed path creates the temp exclusively, closes it, then reopens with
  # gzfile(): exclusive creation is still enforced and the output is real gzip.
  if (compress) {
    close(file(temporary, open = "wxb"))
    connection <- gzfile(temporary, open = "wb")
  } else {
    connection <- file(temporary, open = "wxb")
  }
  old_options <- options(digits = 17, scipen = 999)
  on.exit({
    options(old_options)
    try(close(connection), silent = TRUE)
    if (file.exists(temporary)) unlink(temporary)
  }, add = TRUE)
  write.table(
    as.data.frame(object), connection, sep = sep, row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA", eol = "\n"
  )
  close(connection)
  publish_candidate(temporary, destination)
}

candidate_save_rds <- function(object, destination) {
  temporary <- candidate_temp(destination)
  # R parses "wbx" as a TEXT connection (the trailing x defeats binary
  # detection), which breaks saveRDS ("binary-mode connection required for
  # ascii=FALSE") and gzcon. "wxb" is binary AND still exclusive-create.
  connection <- file(temporary, open = "wxb")
  on.exit({
    try(close(connection), silent = TRUE)
    if (file.exists(temporary)) unlink(temporary)
  }, add = TRUE)
  saveRDS(object, connection)
  close(connection)
  publish_candidate(temporary, destination)
}

candidate_write_lines <- function(lines, destination) {
  temporary <- candidate_temp(destination)
  # R parses "wbx" as a TEXT connection (the trailing x defeats binary
  # detection), which breaks saveRDS ("binary-mode connection required for
  # ascii=FALSE") and gzcon. "wxb" is binary AND still exclusive-create.
  connection <- file(temporary, open = "wxb")
  on.exit({
    try(close(connection), silent = TRUE)
    if (file.exists(temporary)) unlink(temporary)
  }, add = TRUE)
  writeLines(lines, connection, sep = "\n", useBytes = TRUE)
  close(connection)
  publish_candidate(temporary, destination)
}

# ---- input loading (identical to 05_dream lines 44-99) ----
dge <- readRDS(DGE_PATH)
config_path <- Sys.getenv("MASLD_CONFIG_PATH", file.path(ROOT, "config/human_datasets.yaml"))
ycfg <- yaml::read_yaml(config_path)$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep <- dge$samples$dataset %in% mega
dge_mega <- dge[, keep]
meta_new <- readRDS(META_PATH)
sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]
if (anyNA(sex)) stop("Missing inferred_sex for model samples")
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control","Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(sex))
rownames(info) <- colnames(dge_mega)
cat(sprintf("[in] %d genes x %d samples; %d cohorts; %s\n", nrow(dge_mega), ncol(dge_mega),
            nlevels(info$dataset), paste(names(table(info$group_binary)), table(info$group_binary), collapse=" ")))

# ---- C2 design: ~ dataset + sex + group  (cohort + sex fixed effects) ----
design <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
coef_name <- "group_binaryDisease"
if (!coef_name %in% colnames(design)) stop("Disease coefficient is not estimable")
if (qr(design)$rank != ncol(design)) stop("C2 design matrix is rank deficient")

# ---- limma_voom_qw (benchmark-winning engine) ----
v    <- voomWithQualityWeights(dge_mega, design, save.plot = REMEDIATION_MODE)
if (REMEDIATION_MODE) {
  sw <- v$targets$sample.weights
  target_ids <- rownames(v$targets)
  if (is.null(sw) || length(sw) != ncol(dge_mega) ||
      !identical(target_ids, colnames(dge_mega)) ||
      any(!is.finite(sw)) || any(sw <= 0)) {
    stop("voomWithQualityWeights did not return exact ordered positive finite sample weights")
  }
  names(sw) <- target_ids
  validate_voom_curve <- function(curve, label) {
    if (is.null(curve) || is.null(curve$x) || is.null(curve$y) ||
        !is.numeric(curve$x) || !is.numeric(curve$y) ||
        length(curve$x) == 0L || length(curve$x) != length(curve$y) ||
        any(!is.finite(curve$x)) || any(!is.finite(curve$y))) {
      stop("voomWithQualityWeights returned an invalid ", label, " curve")
    }
  }
  validate_voom_curve(v$voom.xy, "observed mean-variance")
  validate_voom_curve(v$voom.line, "fitted mean-variance")
}
fit0 <- lmFit(v, design)
fit  <- eBayes(fit0)
res  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)
dt <- as.data.table(res); setnames(dt, "adj.P.Val", "padj")
coef_index <- match(coef_name, colnames(fit$coefficients))
moderated_se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
names(moderated_se) <- rownames(fit$coefficients)
dt[, SE := moderated_se[match(gene, names(moderated_se))]]

# ---- TREAT: canonical DEG significance gate (2026-06-29) ----
# treat() tests H0: |true logFC| <= TREAT_LFC (McCarthy & Smyth 2009). FDR<0.05 on this
# test IS the canonical DEG call -- the effect-size floor is folded INTO the test, so NO
# separate |logFC| filter is applied downstream. Supersedes the ashr lfsr+|shrunk|>0.3 gate.
TREAT_LFC <- as.numeric(Sys.getenv("CANONICAL_TREAT_LFC", "0.25"))
if (!is.finite(TREAT_LFC) || length(TREAT_LFC) != 1L) {
  stop("CANONICAL_TREAT_LFC must be a single finite number; got: ",
       Sys.getenv("CANONICAL_TREAT_LFC"), call. = FALSE)
}
# BG-004 (2026-08-09): the canonical contract is lfc = 0.25. CANONICAL_TREAT_LFC used to
# be able to silently change the canonical threshold while keeping the canonical FILENAME,
# so a differently-gated table could sit at the headline path indistinguishably. The value
# is still tunable for exploration, but a non-contract value now blocks promotion (see the
# promotion gate at the end of this script).
CANONICAL_CONTRACT_LFC <- 0.25
TREAT_LFC_IS_CANONICAL <- isTRUE(all.equal(TREAT_LFC, CANONICAL_CONTRACT_LFC))
if (REMEDIATION_MODE && !TREAT_LFC_IS_CANONICAL) {
  stop("BG-001 remediation requires CANONICAL_TREAT_LFC=0.25")
}
if (!TREAT_LFC_IS_CANONICAL) {
  cat(sprintf("\n[WARN] CANONICAL_TREAT_LFC=%.4f differs from the canonical contract %.2f.\n",
              TREAT_LFC, CANONICAL_CONTRACT_LFC),
      "       This run will NOT be promoted to the canonical filename.\n", sep = "")
}
ttm <- topTreat(treat(fit0, lfc = TREAT_LFC), coef = coef_name, number = Inf, sort.by = "none")
dt[, treat_lfc := TREAT_LFC]
dt[, treat_p   := ttm[gene, "P.Value"]]
dt[, treat_fdr := ttm[gene, "adj.P.Val"]]

# ---- ashr shrinkage (kept SIDE-BY-SIDE with raw) ----
ok <- is.finite(dt$logFC) & is.finite(dt$SE) & dt$SE > 0
ash <- ashr::ash(dt$logFC[ok], dt$SE[ok], mixcompdist = "normal")
dt[, shrunk_logFC := NA_real_][ok, shrunk_logFC := ash$result$PosteriorMean]
dt[, lfsr := NA_real_][ok, lfsr := ash$result$lfsr]

out <- dt[, .(gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, treat_lfc, treat_p, treat_fdr, AveExpr)]

# ---- symbol column (drop-in parity with dream_results_ashr.csv `symbol`) ----
# Class-B downstream readers (e.g. 80_celltype_intrinsic_attribution.R, 310a) select
# `symbol`; map via GENCODE v49 metadata so canonical_deg_results.csv is a true drop-in.
gene_metadata_path <- Sys.getenv("MASLD_GENE_METADATA_PATH", file.path(ROOT, "data/gencode_v49_gene_metadata.tsv.gz"))
gm0 <- fread(gene_metadata_path)
gm0[, eb := sub("[.][0-9]+$", "", gene_id)]
out[, symbol := gm0[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]

if (REMEDIATION_MODE) {
  core_numeric <- c("logFC", "SE", "t", "P.Value", "padj", "treat_lfc", "treat_p", "treat_fdr", "AveExpr")
  probability_fields <- c("P.Value", "padj", "treat_p", "treat_fdr")
  if (!identical(out$gene, rownames(dge_mega)) ||
      any(vapply(core_numeric, function(field) any(!is.finite(out[[field]])), logical(1))) ||
      any(out$SE <= 0) || any(out$treat_lfc != 0.25) ||
      any(vapply(probability_fields, function(field) {
        any(out[[field]] < 0 | out[[field]] > 1)
      }, logical(1)))) {
    stop("Candidate coefficient/TREAT output failed the exact per-arm numerical contract")
  }
  # NA-safe symlink test (see candidate_temp above): a bare nzchar() here
  # refused every fresh arm, because the candidate DEG output does not exist yet.
  if (file.exists(DEG_OUTPUT) ||
      isTRUE(nzchar(Sys.readlink(DEG_OUTPUT), keepNA = TRUE))) {
    stop("Refusing an existing candidate DEG output")
  }
  sample_weights <- data.table(sample_id = colnames(dge_mega))
  sample_weights[, sample_weight := as.numeric(sw)]
  candidate_write_table(sample_weights, file.path(RDIR, "lvqw_sample_weights.tsv"))
  design_table <- as.data.table(design)
  design_table[, sample_id := colnames(dge_mega)]
  setcolorder(design_table, c("sample_id", setdiff(names(design_table), "sample_id")))
  candidate_write_table(design_table, file.path(RDIR, "model_design.tsv"))
  candidate_save_rds(
    list(voom_xy = v$voom.xy, voom_line = v$voom.line),
    file.path(RDIR, "voom_curve.rds")
  )
  candidate_write_table(
    data.table(mean_log_count = v$voom.xy$x, sqrt_residual_sd = v$voom.xy$y),
    file.path(RDIR, "voom_mean_variance.tsv")
  )
  candidate_write_table(
    out[, .(gene, symbol, logFC, SE, t, P.Value, padj, AveExpr)],
    file.path(RDIR, "coefficient_table.tsv.gz"), compress = TRUE
  )
  candidate_write_table(
    out[, .(gene, symbol, logFC, SE, AveExpr, treat_lfc, treat_p, treat_fdr)],
    file.path(RDIR, "treat_table.tsv.gz"), compress = TRUE
  )
  manifest <- data.table(
    dge_path = normalizePath(DGE_PATH, mustWork = TRUE),
    meta_path = normalizePath(META_PATH, mustWork = TRUE),
    config_path = normalizePath(config_path, mustWork = TRUE),
    n_genes = nrow(dge_mega),
    n_samples = ncol(dge_mega),
    n_cohorts = nlevels(info$dataset),
    n_control = sum(info$group_binary == "Control"),
    n_disease = sum(info$group_binary == "Disease"),
    treat_lfc = TREAT_LFC,
    gene_metadata_path = normalizePath(gene_metadata_path, mustWork = TRUE)
  )
  candidate_write_table(manifest, file.path(RDIR, "model_input_manifest.tsv"))
} else {
  # BG-004 (2026-08-09): write ONLY the method-named artifact here. Until this change the
  # canonical filename was also written at this point -- before a single diagnostic had run
  # -- so a run that later died, or one launched with a nonstandard threshold, replaced the
  # headline DEG table anyway. Promotion to DEG_OUTPUT is now the LAST action in the script
  # and is conditional on the promotion gate.
  fwrite(out, file.path(RDIR, "limma_voom_qw_C2_results.csv"))
}

# ---- DEG counts: canonical TREAT + legacy reference definitions ----
treat1 <- out[treat_fdr < 0.05]
raw1 <- out[padj < 0.05 & abs(logFC) > 0.5]
ash1 <- out[lfsr < 0.05 & abs(shrunk_logFC) > 0.5]
cat(sprintf("\n[C2] genes=%d\n  CANONICAL TREAT (treat_fdr<.05, lfc=%.2f)  : %d  (up %d / dn %d)\n  RAW  ref (padj<.05 & |logFC|>.5)           : %d  (up %d / dn %d)\n  ashr ref (lfsr<.05 & |shrunk_logFC|>.5)    : %d  (up %d / dn %d)\n",
  nrow(out), TREAT_LFC, nrow(treat1), sum(treat1$logFC>0), sum(treat1$logFC<0),
  nrow(raw1), sum(raw1$logFC>0), sum(raw1$logFC<0),
  nrow(ash1), sum(ash1$shrunk_logFC>0), sum(ash1$shrunk_logFC<0)))

# ======================= SANITY GATE =======================
if (REMEDIATION_MODE) {
  reference_path <- Sys.getenv("MASLD_SANITY_REFERENCE", file.path(SOURCE_RDIR, "canonical_deg_results.csv"))
  reference <- fread(reference_path)
  common <- merge(
    out[, .(gene, candidate_logFC = logFC, candidate_SE = SE, candidate_treat_fdr = treat_fdr)],
    reference[, .(gene, reference_logFC = logFC, reference_SE = SE, reference_treat_fdr = treat_fdr)],
    by = "gene"
  )
  sanity <- data.table(
    metric = c("candidate_genes", "reference_genes", "common_genes", "logFC_spearman", "logFC_max_abs_delta", "treat_jaccard"),
    value = c(
      nrow(out), nrow(reference), nrow(common),
      cor(common$candidate_logFC, common$reference_logFC, method = "spearman", use = "complete.obs"),
      max(abs(common$candidate_logFC - common$reference_logFC), na.rm = TRUE),
      {
        a <- out[is.finite(treat_fdr) & treat_fdr < 0.05, gene]
        b <- reference[is.finite(treat_fdr) & treat_fdr < 0.05, gene]
        length(intersect(a, b)) / length(union(a, b))
      }
    )
  )
  candidate_write_table(sanity, file.path(RDIR, "candidate_vs_canonical_sanity.tsv"))
  candidate_write_lines(capture.output(sessionInfo()), file.path(RDIR, "session_info.txt"))
  required_columns <- c("gene", "logFC", "SE", "P.Value", "padj", "treat_lfc", "treat_p", "treat_fdr", "AveExpr", "symbol")
  if (!all(required_columns %in% names(out)) || anyDuplicated(out$gene) || nrow(out) != nrow(dge_mega) ||
      any(!is.finite(out$treat_lfc)) || any(out$treat_lfc != 0.25)) {
    stop("Candidate DEG schema or gene cardinality sanity check failed")
  }
  candidate_write_table(out, DEG_OUTPUT, sep = ",")
  cat("\n[done] wrote candidate DEG and provenance under ", RDIR, "\n", sep = "")
  quit(save = "no", status = 0L)
}

dr <- fread(file.path(RDIR, "dream_results.csv"))
m  <- merge(out[, .(gene, c2_lfc=logFC, c2_padj=padj, c2_slfc=shrunk_logFC, c2_lfsr=lfsr)],
            dr[, .(gene, dr_lfc=logFC, dr_padj=padj,
                   dr_slfc=if ("shrunk_logFC" %in% names(dr)) shrunk_logFC else NA_real_,
                   dr_lfsr=if ("lfsr" %in% names(dr)) lfsr else NA_real_)], by="gene")
rho  <- cor(m$c2_lfc, m$dr_lfc, method="spearman", use="complete.obs")
pear <- cor(m$c2_lfc, m$dr_lfc, method="pearson",  use="complete.obs")
dirA <- mean(sign(m$c2_lfc)==sign(m$dr_lfc), na.rm=TRUE)
# dream canonical Tier-1 (ashr if present, else raw)
dr1 <- if (all(c("shrunk_logFC","lfsr") %in% names(dr))) dr[lfsr<0.05 & abs(shrunk_logFC)>0.5, gene] else dr[padj<0.05 & abs(logFC)>0.5, gene]
jac <- function(a,b) length(intersect(a,b))/length(union(a,b))
# positive-control recovery
pc <- fread(file.path(ROOT, "RNA-seq/results/validation/positive_control_validation.csv"))
gm <- fread(gene_metadata_path)
gm[, eb := sub("[.][0-9]+$","",gene_id)]
pc_eb <- unique(na.omit(gm[match(pc$gene, gene_name), eb]))
out[, eb := sub("[.][0-9]+$","",gene)]
pc_in_raw  <- length(intersect(pc_eb, out[padj<0.05 & abs(logFC)>0.5, eb]))
pc_in_ash  <- length(intersect(pc_eb, out[lfsr<0.05 & abs(shrunk_logFC)>0.5, eb]))
pc_tested  <- length(intersect(pc_eb, out$eb))

san <- data.table(
  metric = c("genes_tested","C2_raw_Tier1","C2_ashr_Tier1","dream_Tier1",
             "logFC_spearman_vs_dream","logFC_pearson_vs_dream","direction_agree_vs_dream",
             "jaccard_C2ashr_vs_dreamTier1","PC_tested","PC_in_C2raw","PC_in_C2ashr"),
  value  = c(nrow(out), nrow(raw1), nrow(ash1), length(dr1),
             round(rho,4), round(pear,4), round(dirA,4),
             round(jac(ash1$gene, dr1),4), pc_tested, pc_in_raw, pc_in_ash))
fwrite(san, file.path(RDIR, "limma_voom_qw_C2_sanity.csv"))
cat("\n=== SANITY GATE vs current dream canonical ===\n"); print(san, row.names=FALSE)

# ======================= BG-004 PROMOTION GATE =======================
# These are catastrophe detectors, not scientific adjudication: they are set loose enough
# that a legitimate methodological change still passes, and tight enough that a scrambled,
# truncated, misaligned or wrongly-thresholded table cannot reach the canonical filename.
gate <- data.table(check = character(), observed = character(),
                   required = character(), pass = logical())
add_gate <- function(check, observed, required, pass) {
  gate <<- rbind(gate, data.table(check = check, observed = as.character(observed),
                                  required = as.character(required), pass = isTRUE(pass)))
}
required_columns <- c("gene", "logFC", "SE", "P.Value", "padj", "shrunk_logFC", "lfsr",
                      "treat_lfc", "treat_p", "treat_fdr", "AveExpr", "symbol")
add_gate("treat_lfc_contract", TREAT_LFC, CANONICAL_CONTRACT_LFC, TREAT_LFC_IS_CANONICAL)
add_gate("schema_complete", sum(required_columns %in% names(out)), length(required_columns),
         all(required_columns %in% names(out)))
add_gate("gene_cardinality", nrow(out), nrow(dge_mega), nrow(out) == nrow(dge_mega))
add_gate("no_duplicate_genes", anyDuplicated(out$gene), 0, anyDuplicated(out$gene) == 0L)
add_gate("logFC_all_finite", sum(!is.finite(out$logFC)), 0, all(is.finite(out$logFC)))
add_gate("SE_positive_finite", sum(!(is.finite(out$SE) & out$SE > 0)), 0,
         all(is.finite(out$SE) & out$SE > 0))
add_gate("treat_lfc_column_locked", length(unique(out$treat_lfc)), 1,
         length(unique(out$treat_lfc)) == 1L &&
           isTRUE(all.equal(unique(out$treat_lfc), TREAT_LFC)))
add_gate("tier1_nonempty", nrow(treat1), ">0", nrow(treat1) > 0L)
add_gate("positive_controls_recovered", pc_in_raw, ">0", pc_in_raw > 0L)

# Continuity against whatever currently occupies the canonical path.
if (file.exists(DEG_OUTPUT)) {
  inc <- fread(DEG_OUTPUT)
  cm  <- merge(out[, .(gene, new_lfc = logFC)], inc[, .(gene, old_lfc = logFC)], by = "gene")
  rho_inc <- suppressWarnings(cor(cm$new_lfc, cm$old_lfc, method = "spearman", use = "complete.obs"))
  dir_inc <- mean(sign(cm$new_lfc) == sign(cm$old_lfc), na.rm = TRUE)
  add_gate("vs_incumbent_common_genes", nrow(cm), ">=10000", nrow(cm) >= 10000L)
  add_gate("vs_incumbent_logFC_spearman", round(rho_inc, 4), ">=0.95",
           is.finite(rho_inc) && rho_inc >= 0.95)
  add_gate("vs_incumbent_direction_agree", round(dir_inc, 4), ">=0.90",
           is.finite(dir_inc) && dir_inc >= 0.90)
} else {
  add_gate("vs_incumbent", "no incumbent file", "skipped", TRUE)
}

FORCE <- Sys.getenv("CANONICAL_PROMOTE_FORCE", "")
gate <- rbind(gate, data.table(check = "force_override",
                               observed = if (nzchar(FORCE)) FORCE else "(unset)",
                               required = "(unset)", pass = TRUE))
fwrite(gate, file.path(RDIR, "limma_voom_qw_C2_promotion_gate.csv"))
cat("\n=== PROMOTION GATE ===\n"); print(gate, row.names = FALSE)

failed <- gate[check != "force_override" & !pass, check]
if (length(failed) && !nzchar(FORCE)) {
  cat("\n[BLOCKED] canonical NOT promoted; failing checks: ",
      paste(failed, collapse = ", "), "\n",
      "  Candidate retained at ", file.path(RDIR, "limma_voom_qw_C2_results.csv"), "\n",
      "  ", basename(DEG_OUTPUT), " is UNCHANGED.\n",
      "  To promote anyway, set CANONICAL_PROMOTE_FORCE='<reason>' (recorded in the gate file).\n",
      sep = "")
  quit(save = "no", status = 1L)
}
if (length(failed)) {
  cat("\n[FORCED] promoting despite failing checks (", paste(failed, collapse = ", "),
      ") on the recorded reason: ", FORCE, "\n", sep = "")
}

# Archive the incumbent BEFORE replacing it, so promotion is reversible. The stamp uses the
# incumbent's own mtime, not wall-clock, so the archived name identifies WHICH canonical it
# was rather than when it happened to be archived.
if (file.exists(DEG_OUTPUT)) {
  archive_dir <- file.path(RDIR, "archive")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(file.info(DEG_OUTPUT)$mtime, "%Y%m%dT%H%M%S")
  archived <- file.path(archive_dir,
                        sub("\\.csv$", paste0(".superseded_", stamp, ".csv"), basename(DEG_OUTPUT)))
  if (!file.exists(archived) && !file.copy(DEG_OUTPUT, archived)) {
    stop("Refusing to promote: could not archive the incumbent to ", archived, call. = FALSE)
  }
  cat("\n[archive] incumbent preserved at ", archived, "\n", sep = "")
}

# Atomic promotion: same-directory rename, so a reader never observes a partial file.
promote_tmp <- file.path(dirname(DEG_OUTPUT), paste0(".", basename(DEG_OUTPUT), ".promote.tmp"))
if (file.exists(promote_tmp)) unlink(promote_tmp)
fwrite(out, promote_tmp)
if (!file.rename(promote_tmp, DEG_OUTPUT)) {
  unlink(promote_tmp)
  stop("Atomic promotion failed; ", basename(DEG_OUTPUT), " left unchanged.", call. = FALSE)
}
cat("\n[done] wrote limma_voom_qw_C2_results.csv + limma_voom_qw_C2_sanity.csv",
    "\n[done] PROMOTED ", basename(DEG_OUTPUT), " (all promotion-gate checks passed)\n", sep = "")
