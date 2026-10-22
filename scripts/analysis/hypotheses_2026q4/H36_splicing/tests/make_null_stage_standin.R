# Builds a NULL stand-in for stage_test/ from the first completed chunks: "observed" := permutation 1,
# "permuted" := permutations 2..10. Code-path test of h36_05 only; no real stage outcome is used.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE); real <- args[1]; fake <- args[2]; nch <- as.integer(args[3])
dir.create(fake, recursive = TRUE, showWarnings = FALSE)
R <- rbindlist(lapply(sprintf("%s/chunk_%03d.rds", real, seq_len(nch)), readRDS), fill = TRUE)
R <- R[status == "ok"]
R[, `:=`(lr = lr_perm1, p = pchisq(lr_perm1, df, lower.tail = FALSE), max_abs_dpsi = dpsi_perm1, conv1 = conv_perm1)]
R[, c("lr_perm1", "dpsi_perm1", "conv_perm1") := NULL]
for (k in 2:10) setnames(R, c(paste0("lr_perm", k), paste0("dpsi_perm", k), paste0("conv_perm", k)),
                         c(paste0("lr_perm", k - 1), paste0("dpsi_perm", k - 1), paste0("conv_perm", k - 1)))
cl_tab <- fread(file.path(args[4], "clusters.tsv.gz"))
R <- merge(cl_tab[, .(cl, chr, str, start_min, end_max, J_tested_in_build = J)], R, by = "cl")
allp <- grep("_perm[0-9]+$", names(R), value = TRUE)
fwrite(R[, setdiff(names(R), allp), with = FALSE], file.path(fake, "stage_test_observed.tsv.gz"), sep = "\t")
fwrite(R[, c("cl", "df", allp), with = FALSE], file.path(fake, "stage_test_perm.tsv.gz"), sep = "\t")
fwrite(data.table(cl = 1:5, n_persons = 700L, n_cohorts = 4L, lr = 1, df = 3L, p = 0.8), file.path(fake, "cohort_effect_pipeline_check.tsv"), sep = "\t")
cat("fake stage rows:", nrow(R), "\n")
