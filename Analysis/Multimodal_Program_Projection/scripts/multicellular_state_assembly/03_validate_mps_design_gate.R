#!/usr/bin/env Rscript

# Validate the Plan 44 cue-by-lineage design from GEO metadata only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: 03_validate_mps_design_gate.R <GSE168285 SOFT.gz> <output.tsv> <preflight_manifest_sha256>")
}

soft_path <- args[[1L]]
output_path <- args[[2L]]
manifest_sha256 <- args[[3L]]

stopifnot(file.exists(soft_path), grepl("^[0-9a-f]{64}$", manifest_sha256))
lines <- readLines(gzfile(soft_path), warn = FALSE)
titles <- sub("^!Sample_title = ", "", grep("^!Sample_title = ", lines, value = TRUE))
treatments <- sub(
  "^!Sample_characteristics_ch1 = treatment: ",
  "",
  grep("^!Sample_characteristics_ch1 = treatment: ", lines, value = TRUE)
)
stopifnot(length(titles) == 179L, length(treatments) == 179L)

sample_design <- data.frame(
  experiment = sub("_.*$", "", titles),
  treatment = treatments,
  stringsAsFactors = FALSE
)
condition_design <- unique(sample_design)
stopifnot(nrow(condition_design) == 60L)

condition_design$high_NPC <- as.integer(!grepl("^Reduced", condition_design$treatment))
condition_design$fat <- as.integer(grepl("(^|_)Fat($|_)", condition_design$treatment))
condition_design$fructose <- as.integer(grepl("Fructose", condition_design$treatment))
condition_design$cholesterol <- as.integer(grepl("Cholesterol", condition_design$treatment))
condition_design$LPS <- as.integer(grepl("LPS", condition_design$treatment))
condition_design$TGF_beta <- as.integer(grepl("TGF-B", condition_design$treatment))

long_design <- do.call(
  rbind,
  lapply(
    c("hepatocyte", "fibroblast", "macrophage"),
    function(lineage_name) transform(condition_design, lineage = lineage_name)
  )
)
long_design$experiment <- factor(long_design$experiment)
long_design$lineage <- factor(
  long_design$lineage,
  levels = c("hepatocyte", "fibroblast", "macrophage")
)

design <- model.matrix(
  ~ experiment + lineage + high_NPC + fat + fructose + cholesterol + LPS +
    TGF_beta + lineage:(high_NPC + fat + LPS + TGF_beta),
  data = long_design
)
design_rank <- qr(design)$rank
condition_number <- kappa(design)

expanded <- design[, colnames(design) != "(Intercept)", drop = FALSE]
vifs <- vapply(seq_len(ncol(expanded)), function(column_index) {
  response <- expanded[, column_index]
  predictors <- cbind(1, expanded[, -column_index, drop = FALSE])
  fit <- lm.fit(predictors, response)
  residual_ss <- sum(fit$residuals^2)
  total_ss <- sum((response - mean(response))^2)
  if (total_ss <= 0 || residual_ss <= 0) Inf else total_ss / residual_ss
}, numeric(1L))
names(vifs) <- colnames(expanded)
primary_columns <- grep(
  "lineage.*:(high_NPC|fat|TGF_beta)",
  names(vifs),
  value = TRUE
)
primary_vif_max <- max(vifs[primary_columns])

gate_pass <- design_rank == ncol(design) && is.finite(primary_vif_max) && primary_vif_max < 10
result <- data.frame(
  gate = "mps_primary_cue_by_lineage_design",
  status = if (gate_pass) "pass_full_rank_estimable" else "fail_rank_or_vif",
  n_source_samples = nrow(sample_design),
  n_condition_means = nrow(condition_design),
  n_long_rows = nrow(design),
  n_model_columns = ncol(design),
  design_rank = design_rank,
  condition_number = format(condition_number, digits = 12),
  primary_vif_max = format(primary_vif_max, digits = 12),
  primary_columns = paste(primary_columns, collapse = ";"),
  preflight_manifest_sha256 = manifest_sha256,
  new_expression_outcomes_accessed = "false",
  stringsAsFactors = FALSE
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(result, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
if (!gate_pass) stop("Plan 44 MPS design gate failed")
cat(
  "PLAN44_MPS_DESIGN_GATE_PASS",
  paste0("rows=", nrow(design)),
  paste0("columns=", ncol(design)),
  paste0("rank=", design_rank),
  paste0("primary_vif_max=", format(primary_vif_max, digits = 6)),
  sep = "\t"
)
cat("\n")
