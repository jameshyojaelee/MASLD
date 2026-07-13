#!/usr/bin/env Rscript

# Experimental calibration of the seqfunc substrate against Hu et al. 2025
# MASLD MPRA DAVs (HepG2/LX-2, control and PAOA/TGF-beta conditions).
#
# This is an apply-only benchmark. It never writes the atlas or convergence score.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(jsonlite)
})

ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
FM <- file.path(ROOT, "GWAS/finemapping")
SF <- file.path(FM, "results/seqfunc")
SRC_DATA <- Sys.getenv(
  "HU2025_MPRA_DIR",
  file.path(FM, "data/seqfunc_external/hu2025_mpra")
)
RAW <- file.path(SRC_DATA, "raw")
OUT <- file.path(SF, "mpra_benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NBOOT <- as.integer(Sys.getenv("MPRA_NBOOT", "500"))

required <- c(
  file.path(SRC_DATA, "tested_oligo_intervals.tsv"),
  file.path(RAW, paste0("TableS", 2:6, ".xlsx")),
  file.path(SF, "variant_substrate_hg38.tsv"),
  file.path(SF, "chrombpnet_accessibility.tsv")
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required MPRA inputs:\n", paste(missing, collapse = "\n"),
       "\nRun 92_download_hu2025_mpra.sh first.")
}

num <- function(x) suppressWarnings(as.numeric(x))
bool <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")

parse_interval <- function(x) {
  out <- data.table(
    coord = x,
    chr = sub(":.*", "", x),
    start = suppressWarnings(as.integer(sub(".*:([0-9]+)-.*", "\\1", x))),
    end = suppressWarnings(as.integer(sub(".*-", "", x)))
  )
  out <- out[!is.na(start) & !is.na(end) & start <= end & grepl("^chr", chr)]
  out[, pos_hg38 := as.integer((start + end) / 2)]
  out
}

read_dav <- function(path, context, cell_model, stimulus) {
  x <- as.data.table(read_excel(path, skip = 2))
  if (ncol(x) < 5) stop("Unexpected DAV table shape: ", path)
  setnames(x, names(x)[1:5], c("coord", "mpra_log2fc", "mpra_p", "mpra_fdr", "rsid"))
  x <- x[!is.na(coord)]
  iv <- parse_interval(as.character(x$coord))
  x <- merge(x, iv[, .(coord, chr, pos_hg38)], by = "coord", all = FALSE)
  x[, `:=`(
    mpra_log2fc = num(mpra_log2fc),
    mpra_p = num(mpra_p),
    mpra_fdr = num(mpra_fdr),
    context = context,
    cell_model = cell_model,
    stimulus = stimulus,
    dav_source_call = TRUE
  )]
  unique(x, by = c("chr", "pos_hg38", "context"))
}

contexts <- data.table(
  table = 2:5,
  context = c("HepG2_ctrl", "LX2_ctrl", "HepG2_PAOA", "LX2_TGFb"),
  cell_model = c("HepG2", "LX2", "HepG2", "LX2"),
  stimulus = c("control", "control", "PAOA", "TGFb")
)
davs <- rbindlist(lapply(seq_len(nrow(contexts)), function(i) {
  z <- contexts[i]
  read_dav(
    file.path(RAW, paste0("TableS", z$table, ".xlsx")),
    z$context, z$cell_model, z$stimulus
  )
}))

tested_raw <- fread(
  file.path(SRC_DATA, "tested_oligo_intervals.tsv"),
  header = FALSE, col.names = "coord"
)
tested <- unique(parse_interval(tested_raw$coord), by = c("chr", "pos_hg38"))
tested[, mpra_barcode_mapped := TRUE]

substrate <- fread(file.path(SF, "variant_substrate_hg38.tsv"))
substrate[, max_pip := num(max_pip)]
substrate[, pos_hg38 := as.integer(pos_hg38)]
substrate[, coloc_bool := bool(colocalizes)]
setorder(substrate, chr, pos_hg38, -max_pip, variant_id_hg19)
substrate_one <- substrate[, .SD[1], by = .(chr, pos_hg38)]

universe <- merge(tested, substrate_one, by = c("chr", "pos_hg38"), all = FALSE)
universe[, block_1mb := paste(chr, floor(pos_hg38 / 1e6), sep = ":")]

cbp <- fread(file.path(SF, "chrombpnet_accessibility.tsv"))
cbp[, `:=`(
  cbp_abs_logfc = num(cbp_abs_logfc),
  cbp_jsd = num(cbp_jsd),
  cbp_abs_logfc_x_jsd = num(cbp_abs_logfc_x_jsd),
  cbp_logfc = num(cbp_logfc),
  max_pip_cbp = num(max_pip),
  in_hepg2_peak = bool(in_hepg2_peak)
)]
setorder(cbp, variant_id_hg19, -cbp_abs_logfc)
cbp <- cbp[, .SD[1], by = variant_id_hg19]
universe <- merge(
  universe,
  cbp[, .(variant_id_hg19, cbp_abs_logfc, cbp_jsd, cbp_abs_logfc_x_jsd,
          cbp_logfc, in_hepg2_peak)],
  by = "variant_id_hg19", all.x = TRUE
)

borzoi_path <- file.path(SF, "borzoi_magnitude_candidates.tsv")
if (file.exists(borzoi_path)) {
  bz <- fread(borzoi_path)
  bz[, `:=`(
    borzoi_abs_logsed = num(borzoi_abs_logsed),
    borzoi_logsed_signed = num(borzoi_logsed_signed),
    magnitude_percentile = num(magnitude_percentile)
  )]
  setorder(bz, variant_id_hg19, -borzoi_abs_logsed)
  bz <- bz[, .SD[1], by = variant_id_hg19]
  universe <- merge(
    universe,
    bz[, .(variant_id_hg19, borzoi_abs_logsed, borzoi_logsed_signed,
           magnitude_percentile)],
    by = "variant_id_hg19", all.x = TRUE
  )
}

nom_path <- file.path(SF, "eqtl_absent_nominations.tsv")
if (file.exists(nom_path)) {
  nom <- fread(nom_path)
  nom_ids <- unique(as.character(nom$credible_variant))
  universe[, is_final_nomination_lead := variant_id_hg19 %in% nom_ids]
} else {
  universe[, is_final_nomination_lead := FALSE]
}

# Expand the experimentally tested substrate universe across the four contexts.
truth <- universe[, .SD[rep(seq_len(.N), each = nrow(contexts))]]
truth[, `:=`(
  context = rep(contexts$context, times = nrow(universe)),
  cell_model = rep(contexts$cell_model, times = nrow(universe)),
  stimulus = rep(contexts$stimulus, times = nrow(universe))
)]
truth <- merge(
  truth,
  davs[, .(chr, pos_hg38, context, rsid, mpra_log2fc, mpra_p, mpra_fdr,
           dav_source_call)],
  by = c("chr", "pos_hg38", "context"), all.x = TRUE
)
truth[, dav := as.integer(!is.na(dav_source_call) & dav_source_call)]
truth[, dav_source_call := NULL]

# Supplementary Table S6: CROP-seq truth. The table supplies selected candidate
# target genes, not a complete target-gene universe.
crop_raw <- as.data.table(read_excel(file.path(RAW, "TableS6.xlsx"),
                                     skip = 4, col_names = FALSE))
if (ncol(crop_raw) >= 10) {
  setnames(crop_raw, names(crop_raw)[1:10], c(
    "rsid", "ctrl_log2fc", "ctrl_fdr", "paoa_log2fc", "paoa_fdr",
    "unused", "hepatocyte_peak", "stellate_peak", "kupffer_peak",
    "candidate_target_genes"
  ))
  crop <- crop_raw[!is.na(rsid), .(
    rsid = as.character(rsid),
    ctrl_log2fc = num(ctrl_log2fc), ctrl_fdr = num(ctrl_fdr),
    paoa_log2fc = num(paoa_log2fc), paoa_fdr = num(paoa_fdr),
    hepatocyte_peak = as.character(hepatocyte_peak),
    stellate_peak = as.character(stellate_peak),
    kupffer_peak = as.character(kupffer_peak),
    candidate_target_genes = as.character(candidate_target_genes),
    truth_status = "selected_CROPseq_candidate_set_not_complete_gene_truth"
  )]
  fwrite(crop, file.path(OUT, "cropseq_candidate_truth.tsv"), sep = "\t", na = "")
}

auc_rank <- function(y, s) {
  ok <- is.finite(s) & !is.na(y)
  y <- y[ok]; s <- s[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

average_precision <- function(y, s) {
  ok <- is.finite(s) & !is.na(y)
  y <- y[ok]; s <- s[ok]
  if (!sum(y == 1) || !sum(y == 0)) return(NA_real_)
  o <- order(s, decreasing = TRUE)
  y <- y[o]
  precision <- cumsum(y == 1) / seq_along(y)
  mean(precision[y == 1])
}

block_boot_metric <- function(d, score_col, metric_fun, nboot = NBOOT, seed = 42L) {
  d <- d[is.finite(get(score_col))]
  blocks <- unique(d$block_1mb)
  obs <- metric_fun(d$dav, d[[score_col]])
  if (length(blocks) < 3 || !is.finite(obs)) {
    return(c(estimate = obs, lo = NA_real_, hi = NA_real_))
  }
  set.seed(seed)
  vals <- replicate(nboot, {
    samp <- sample(blocks, length(blocks), replace = TRUE)
    z <- rbindlist(lapply(seq_along(samp), function(i) {
      q <- copy(d[block_1mb == samp[i]])
      q[, boot_copy := i]
      q
    }))
    metric_fun(z$dav, z[[score_col]])
  })
  c(
    estimate = obs,
    lo = unname(quantile(vals, 0.025, na.rm = TRUE)),
    hi = unname(quantile(vals, 0.975, na.rm = TRUE))
  )
}

features <- intersect(c(
  "cbp_abs_logfc", "cbp_jsd", "cbp_abs_logfc_x_jsd", "in_hepg2_peak",
  "borzoi_abs_logsed", "magnitude_percentile", "max_pip"
), names(truth))

metrics <- rbindlist(lapply(unique(truth$context), function(ctx) {
  d <- truth[context == ctx]
  rbindlist(lapply(features, function(f) {
    d[, score_tmp := num(get(f))]
    au <- block_boot_metric(d, "score_tmp", auc_rank)
    ap <- block_boot_metric(d, "score_tmp", average_precision, seed = 84L)
    data.table(
      context = ctx,
      cell_model = unique(d$cell_model),
      stimulus = unique(d$stimulus),
      feature = f,
      n_complete = sum(is.finite(d$score_tmp)),
      n_positive = sum(d$dav == 1 & is.finite(d$score_tmp)),
      n_negative = sum(d$dav == 0 & is.finite(d$score_tmp)),
      auroc = au["estimate"], auroc_lo = au["lo"], auroc_hi = au["hi"],
      average_precision = ap["estimate"], ap_lo = ap["lo"], ap_hi = ap["hi"],
      positive_rate = mean(d$dav[is.finite(d$score_tmp)] == 1),
      bootstrap_unit = "1Mb genomic block",
      n_bootstrap = NBOOT,
      source_unit = "barcode-mapped MPRA interval"
    )
  }))
}))

paired_diff <- function(ctx_a, ctx_b, feature, nboot = NBOOT, seed = 126L) {
  a <- truth[context == ctx_a, .(chr, pos_hg38, block_1mb, y_a = dav,
                                score = num(get(feature)))]
  b <- truth[context == ctx_b, .(chr, pos_hg38, y_b = dav)]
  d <- merge(a, b, by = c("chr", "pos_hg38"))
  d <- d[is.finite(score)]
  blocks <- unique(d$block_1mb)
  obs <- auc_rank(d$y_a, d$score) - auc_rank(d$y_b, d$score)
  set.seed(seed)
  vals <- replicate(nboot, {
    samp <- sample(blocks, length(blocks), replace = TRUE)
    z <- rbindlist(lapply(seq_along(samp), function(i) {
      q <- copy(d[block_1mb == samp[i]])
      q[, boot_copy := i]
      q
    }))
    auc_rank(z$y_a, z$score) - auc_rank(z$y_b, z$score)
  })
  data.table(
    context_a = ctx_a, context_b = ctx_b, feature = feature,
    n_complete = nrow(d), difference_auroc = obs,
    difference_lo = unname(quantile(vals, .025, na.rm = TRUE)),
    difference_hi = unname(quantile(vals, .975, na.rm = TRUE)),
    n_bootstrap = NBOOT,
    interpretation = paste0("AUROC(", ctx_a, ") - AUROC(", ctx_b, ")")
  )
}

comparisons <- data.table(
  a = c("HepG2_ctrl", "HepG2_PAOA", "HepG2_ctrl", "LX2_ctrl"),
  b = c("LX2_ctrl", "LX2_TGFb", "HepG2_PAOA", "LX2_TGFb")
)
paired <- rbindlist(lapply(seq_len(nrow(comparisons)), function(i) {
  rbindlist(lapply(features, function(f) paired_diff(comparisons$a[i], comparisons$b[i], f)))
}))

# All published DAVs that overlap the full substrate, whether or not their oligo
# survived the barcode-map file used for the negative universe.
all_dav_overlap <- merge(
  unique(davs[, .(chr, pos_hg38, context, cell_model, stimulus, rsid,
                  mpra_log2fc, mpra_p, mpra_fdr)]),
  substrate_one,
  by = c("chr", "pos_hg38"), all = FALSE
)
nom_audit <- all_dav_overlap[, .(
  n_contexts = uniqueN(context),
  contexts = paste(sort(unique(context)), collapse = ","),
  rsid = paste(sort(unique(na.omit(rsid))), collapse = ","),
  mpra_log2fc_min = min(mpra_log2fc, na.rm = TRUE),
  mpra_log2fc_max = max(mpra_log2fc, na.rm = TRUE),
  primary_gene = primary_gene[1],
  max_pip = max(max_pip, na.rm = TRUE),
  colocalizes = coloc_bool[1],
  is_final_nomination_lead = any(variant_id_hg19 %in% universe[is_final_nomination_lead == TRUE,
                                                                variant_id_hg19]),
  variant_id_hg19 = variant_id_hg19[1]
), by = .(chr, pos_hg38)]

overlap_summary <- data.table(
  metric = c(
    "source_MPRA_library_variants_reported",
    "barcode_mapped_unique_intervals",
    "barcode_mapped_intervals_overlapping_seqfunc_substrate",
    "published_DAV_positions_overlapping_full_seqfunc_substrate",
    "published_DAV_variant_by_context_overlaps",
    "published_DAV_overlaps_no_observed_bulk_liver_eQTL_coloc",
    "published_DAV_overlaps_maxPIP_ge_0.1",
    "published_DAV_overlaps_maxPIP_ge_0.5",
    "published_DAV_overlaps_maxPIP_ge_0.9",
    "published_DAV_overlaps_retained_as_final_nomination_lead"
  ),
  value = c(
    5369,
    nrow(tested),
    uniqueN(universe[, .(chr, pos_hg38)]),
    uniqueN(all_dav_overlap[, .(chr, pos_hg38)]),
    nrow(unique(all_dav_overlap[, .(chr, pos_hg38, context)])),
    uniqueN(all_dav_overlap[coloc_bool == FALSE, .(chr, pos_hg38)]),
    uniqueN(all_dav_overlap[max_pip >= 0.1, .(chr, pos_hg38)]),
    uniqueN(all_dav_overlap[max_pip >= 0.5, .(chr, pos_hg38)]),
    uniqueN(all_dav_overlap[max_pip >= 0.9, .(chr, pos_hg38)]),
    uniqueN(nom_audit[is_final_nomination_lead == TRUE, .(chr, pos_hg38)])
  ),
  note = c(
    "Source-paper library size; includes variants not retained in the barcode-map file",
    "Reference/Mut interval names collapsed after barcode filtering",
    "Primary benchmark universe before model complete-case filtering",
    "Significant DAV tables S2-S5; exact hg38 interval-center match",
    "Same variant can be a DAV in multiple cell/treatment contexts",
    "No observed bulk-liver eQTL colocalization under the current analysis",
    "Maximum PIP in current substrate",
    "Maximum PIP in current substrate",
    "Maximum PIP in current substrate",
    "Retrospective audit only; selection and MPRA sampling universes differ"
  )
)

fwrite(truth, file.path(OUT, "mpra_substrate_truth.tsv"), sep = "\t", na = "")
fwrite(metrics, file.path(OUT, "mpra_model_benchmark.tsv"), sep = "\t", na = "")
fwrite(paired, file.path(OUT, "mpra_paired_context_differences.tsv"), sep = "\t", na = "")
fwrite(nom_audit, file.path(OUT, "mpra_nomination_audit.tsv"), sep = "\t", na = "")
fwrite(overlap_summary, file.path(OUT, "mpra_overlap_summary.tsv"), sep = "\t", na = "")

best <- metrics[feature == "cbp_abs_logfc"][order(context)]
verdict <- list(
  status = "complete_experimental_calibration_v1",
  apply_only_firewall = TRUE,
  source = list(
    citation = "Hu et al. 2025 preprint, PMCID PMC12633503",
    geo = c("GSE281364", "GSE281367", "GSE281160"),
    source_unit = "barcode-mapped MPRA interval",
    caveat = paste(
      "The negative universe is limited to intervals retained in the GEO barcode map.",
      "DAV labels are official significant calls from Supplementary Tables S2-S5."
    )
  ),
  exact_overlap = as.list(setNames(overlap_summary$value, overlap_summary$metric)),
  chrombpnet_hepg2 = split(best, seq_len(nrow(best))),
  inference_limits = c(
    "HepG2 and LX-2 are cell lines, not adult primary liver cells.",
    "MPRA is episomal and does not preserve endogenous chromatin or 3D contact.",
    "Effect direction requires explicit verification of the source logFC allele convention.",
    "No AlphaGenome output is used to fit or calibrate a model."
  )
)
write_json(verdict, file.path(OUT, "mpra_benchmark_verdict.json"),
           pretty = TRUE, auto_unbox = TRUE, na = "null")

readme <- c(
  "# Hu 2025 MASLD MPRA calibration",
  "",
  "Experimental calibration of the existing seqfunc substrate against HepG2 and LX-2",
  "MPRA DAVs in control and PAOA/TGF-beta conditions. This output is apply-only and",
  "is never read by the atlas or convergence-score builders.",
  "",
  "## Primary units",
  "",
  "- Positive: official significant DAV call in Supplementary Tables S2-S5.",
  "- Negative: barcode-mapped interval without a significant DAV call in that context.",
  "- Resampling: 1-Mb genomic blocks.",
  "- The full source library contains 5,369 variants; the GEO barcode map retains a",
  paste0("  subset of ", nrow(tested), " unique intervals, of which ", nrow(universe),
         " overlap the current seqfunc substrate before model complete-case filtering."),
  "",
  "## ChromBPNet HepG2 zero-shot results",
  "",
  "| Context | N | DAV | AUROC [block-bootstrap 95% CI] | AUPRC |",
  "|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(best)), function(i) {
    r <- best[i]
    sprintf("| %s | %d | %d | %.3f [%.3f, %.3f] | %.3f |",
            r$context, r$n_complete, r$n_positive,
            r$auroc, r$auroc_lo, r$auroc_hi, r$average_precision)
  }, character(1)),
  "",
  "## Interpretation",
  "",
  "- This is a retrospective external calibration, not prospective validation.",
  "- Cell-type matching is evaluated directly: a HepG2 accessibility model should not",
  "  be assumed to transfer to LX-2 stellate-cell DAVs.",
  "- The published DAV-only tables cannot define the full source-paper library universe;",
  "  the GEO barcode map supplies the experimentally retained negative universe.",
  "- MPRA activity does not by itself validate an endogenous adult-liver mechanism.",
  "",
  "## Files",
  "",
  "- `mpra_substrate_truth.tsv`: tested substrate variant x context table.",
  "- `mpra_model_benchmark.tsv`: per-model endpoint metrics.",
  "- `mpra_paired_context_differences.tsv`: matched context comparisons.",
  "- `mpra_nomination_audit.tsv`: exact DAV overlap with current nominations.",
  "- `cropseq_candidate_truth.tsv`: source CROP-seq candidate target-gene annotations.",
  "- `mpra_benchmark_verdict.json`: machine-readable summary."
)
writeLines(readme, file.path(OUT, "README.md"))

message("[92] wrote MPRA benchmark to ", OUT)
message("[92] exact published DAV positions in full substrate: ",
        uniqueN(all_dav_overlap[, .(chr, pos_hg38)]))
message("[92] barcode-mapped substrate benchmark universe: ", nrow(universe))
