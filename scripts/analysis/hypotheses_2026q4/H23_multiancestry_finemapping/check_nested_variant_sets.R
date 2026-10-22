#!/usr/bin/env Rscript
# H23 check: do the reference fit (study A on the variants shared with study B) and the nested EUR fits use the same variants?
# Usage: Rscript check_nested_variant_sets.R <fits_dir> <new_output_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: check_nested_variant_sets.R <fits_dir> <new_output_dir>")
fits_dir <- args[[1]]; out_dir <- args[[2]]
if (file.exists(out_dir)) stop("Refusing to overwrite existing output dir: ", out_dir)
dir.create(out_dir, recursive = TRUE)
obj <- readRDS(file.path(fits_dir, "unit_fits.rds"))
nv <- function(f, nm) if (is.null(f[[nm]])) NA_integer_ else as.integer(f[[nm]]$n_var)
tab <- rbindlist(lapply(obj$fits, function(f) data.table(uid = f$uid, n_VMA = as.integer(f$n_VMA), n_A_aligned = as.integer(f$n_A_aligned),
                                                          n_VB = as.integer(f$n_VB), n_ma = nv(f, "ma"), n_A_ma = nv(f, "A_ma"),
                                                          n_Ab = nv(f, "A_b"), n_B = nv(f, "B"), n_Mf = nv(f, "Mf"), n_Mm = nv(f, "Mm"))), fill = TRUE)
tab <- merge(obj$units[, .(uid, set, trait)], tab, by = "uid")
tab[, `:=`(Ab_minus_Mf = n_Ab - n_Mf, Ab_minus_Mm = n_Ab - n_Mm, B_minus_Ab = n_B - n_Ab, VMA_minus_VB = n_VMA - n_VB)]
fwrite(tab, file.path(out_dir, "variant_set_sizes.tsv"), sep = "\t")
s <- tab[!is.na(n_Mm), .(n_units = .N, n_Ab_ne_Mf = sum(Ab_minus_Mf != 0), n_Ab_ne_Mm = sum(Ab_minus_Mm != 0), max_abs_Ab_minus_Mm = max(abs(Ab_minus_Mm)),
                          n_B_ne_Ab = sum(B_minus_Ab != 0), median_VMA_minus_VB = median(VMA_minus_VB), max_VMA_minus_VB = max(VMA_minus_VB)), by = .(set, trait)]
fwrite(s, file.path(out_dir, "variant_set_summary.tsv"), sep = "\t")
print(s)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
