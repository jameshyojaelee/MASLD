#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

dge <- readRDS(precontract_nonholdout_dge)
meta <- as.data.table(readRDS(precontract_nonholdout_meta))
meta <- meta[match(colnames(dge), sample_id)]
assert_true(!holdout_cohort %in% as.character(dge$samples$dataset),
            "HOLDOUT LEAK: precontract audit DGE contains GSE193066")
assert_true(!holdout_cohort %in% meta$dataset,
            "HOLDOUT LEAK: precontract audit metadata contains GSE193066")
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbols <- collapse_symbols(
  edgeR::cpm(dge, log = TRUE, prior.count = 1), rownames(dge), annotation
)
membership <- fread(program_membership)[, .(
  feature_id = program_uid,
  gene_symbol = mapped_symbol,
  weight = as.numeric(original_l1_weight)
)]
total <- membership[, .(total_weight = sum(weight)), by = feature_id]
rows <- list()
for (cohort in unique(meta$dataset)) {
  j <- which(meta$dataset == cohort)
  variable <- rownames(symbols)[
    apply(symbols[, j, drop = FALSE], 1L, function(x) {
      value <- stats::sd(x, na.rm = TRUE)
      is.finite(value) && value > 0
    })
  ]
  observed <- membership[
    !is.na(gene_symbol) & gene_symbol %in% variable,
    .(observed_weight = sum(weight)), by = feature_id
  ]
  audit <- merge(total, observed, by = "feature_id", all.x = TRUE)
  audit[is.na(observed_weight), observed_weight := 0]
  audit[, `:=`(
    cohort = cohort,
    weight_coverage = observed_weight / total_weight
  )]
  rows[[cohort]] <- audit
}
audit <- rbindlist(rows)
registry <- fread(program_registry)[, .(feature_id = program_uid, cell_type, module_name)]
audit <- merge(audit, registry, by = "feature_id", all.x = TRUE)
setorder(audit, feature_id, cohort)
audit[, holdout_expression_accessed := FALSE]

discovery_audit <- audit[cohort %in% discovery_cohorts]
program_summary <- discovery_audit[, .(
  minimum_discovery_coverage = min(weight_coverage),
  all_discovery_cohorts_testable = all(weight_coverage >= program_weight_coverage)
), by = .(feature_id, cell_type, module_name)]
setorder(program_summary, minimum_discovery_coverage)

dir.create(dirname(precontract_testability_audit), recursive = TRUE, showWarnings = FALSE)
if (file.exists(precontract_testability_audit)) {
  fail("Refusing to overwrite precontract testability audit: ", precontract_testability_audit)
}
fwrite(
  audit, precontract_testability_audit,
  sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
)
Sys.chmod(precontract_testability_audit, mode = "0440")

cat("PROGRAM_SUMMARY\n")
print(program_summary[1:10])
cat("PASS_ALL_DISCOVERY", sum(program_summary$all_discovery_cohorts_testable), "\n")
cat("FAIL_ALL_DISCOVERY", sum(!program_summary$all_discovery_cohorts_testable), "\n")
cat("COHORT_COUNTS\n")
print(discovery_audit[, .(
  testable = sum(weight_coverage >= program_weight_coverage),
  untestable = sum(weight_coverage < program_weight_coverage)
), by = cohort])
