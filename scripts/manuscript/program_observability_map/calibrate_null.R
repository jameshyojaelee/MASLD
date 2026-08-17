#!/usr/bin/env Rscript

# Null calibration for the three axis views.
#
# The first v9 discovery returned counts that cannot all be right: the joint
# 2-df test supported 67 fibrosis and 71 NAS programs while the quadratic term
# alone supported none, and the categorical omnibus supported 53 fibrosis
# programs against 27 for the linear term even though fibrosis is close to
# monotone, where a categorical test should lose power rather than gain it.
#
# Two candidate causes. First, the joint meta test referenced chi-square on 2 df
# while the linear arm referenced a Knapp-Hartung t on 3 df, so the two were
# never comparable. Second, an HC3 Wald test on a categorical block with sparse
# levels is known to be anticonservative, and NAS carries up to nine levels per
# cohort.
#
# This script answers the question empirically. It permutes fibrosis and NAS
# jointly within cohort, which preserves their correlation and every design
# feature except the association being tested, and reports how many programs
# each view calls supported under the null. A calibrated view should return
# close to zero after BH.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

n_replicates <- as.integer(Sys.getenv("MASLD_NULL_REPLICATES", "10"))
output_root <- file.path(workstream_root, "diagnostics", "null_calibration")
if (dir.exists(output_root)) fail("Refusing to overwrite: ", output_root)

set.seed(seed)
message("[1/4] Loading sealed non-holdout input")
dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
assert_true(!holdout_cohort %in% meta$dataset, "HOLDOUT LEAK")

logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()

registry <- fread(program_registry)
membership <- fread(program_membership)
program_scores <- score_programs(symbol_matrix, meta, membership, registry,
                                 program_weight_coverage)
testable <- program_scores$coverage[testable == TRUE, feature_id]
S <- program_scores$primary[testable, , drop = FALSE]
rm(symbol_matrix); gc()

discovery_meta <- meta[dataset %in% discovery_cohorts &
                         !is.na(fibrosis_stage) & !is.na(nas_score)]
assert_true(nrow(discovery_meta) == 469L, "Discovery census drift")

count_supported <- function(md, label) {
  linear <- fit_feature_models(S, md, discovery_cohorts, "primary")
  shape <- fit_axis_shape_models(S, md, discovery_cohorts)
  omni <- categorical_omnibus(S, md, discovery_cohorts)

  lm_meta <- meta_analyze(linear)
  lm_meta[, q := p.adjust(p_value_meta, method = "BH"), by = axis]
  sh_meta <- meta_analyze_shape(shape)[estimable == TRUE]
  # Both references, so the comparison is explicit rather than assumed.
  sh_meta[, p_joint_chisq := p_joint_meta]
  sh_meta[, p_joint_F := stats::pf(joint_chisq / 2, df1 = 2, df2 = n_cohorts - 1,
                                   lower.tail = FALSE)]
  sh_meta[, `:=`(q_joint_chisq = p.adjust(p_joint_chisq, method = "BH"),
                 q_joint_F = p.adjust(p_joint_F, method = "BH")), by = axis]
  om <- combine_omnibus(omni)
  om[, q := p.adjust(p_combined, method = "BH"), by = axis]

  rbindlist(lapply(c("fibrosis", "nas"), function(a) data.table(
    replicate = label, axis = a,
    linear = lm_meta[axis == a & q < 0.05, .N],
    joint_chisq = sh_meta[axis == a & q_joint_chisq < 0.05, .N],
    joint_F = sh_meta[axis == a & q_joint_F < 0.05, .N],
    omnibus = om[axis == a & q < 0.05, .N],
    omnibus_nominal = omni[axis == a & estimable == TRUE & p_value < 0.05, .N]
  )))
}

message("[2/4] Observed data")
rows <- list(count_supported(discovery_meta, "observed"))

message("[3/4] ", n_replicates, " joint within-cohort permutations")
for (b in seq_len(n_replicates)) {
  md <- copy(discovery_meta)
  md[, .perm := sample.int(.N), by = dataset]
  md[, `:=`(fibrosis_stage = fibrosis_stage[.perm], nas_score = nas_score[.perm]),
     by = dataset]
  md[, .perm := NULL]
  rows[[length(rows) + 1L]] <- count_supported(md, paste0("permuted_", b))
  message("      replicate ", b, " done")
}

result <- rbindlist(rows)
tmp <- atomic_dir(output_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)
write_tsv(result, file.path(tmp, "null_calibration.tsv"))

summary_table <- result[replicate != "observed", .(
  n_replicates = .N,
  mean_linear = mean(linear), max_linear = max(linear),
  mean_joint_chisq = mean(joint_chisq), max_joint_chisq = max(joint_chisq),
  mean_joint_F = mean(joint_F), max_joint_F = max(joint_F),
  mean_omnibus = mean(omnibus), max_omnibus = max(omnibus),
  mean_omnibus_nominal = mean(omnibus_nominal)
), by = axis]
write_tsv(summary_table, file.path(tmp, "null_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))
publish_dir(tmp, output_root)

cat("\n=== OBSERVED ===\n")
print(result[replicate == "observed"])
cat("\n=== NULL (", n_replicates, " joint within-cohort permutations) ===\n", sep = "")
print(summary_table)
cat("\nA calibrated view returns close to zero under the null.\n")
