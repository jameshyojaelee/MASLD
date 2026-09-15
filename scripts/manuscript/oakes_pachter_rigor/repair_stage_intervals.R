#!/usr/bin/env Rscript
# Reproduce all 11 existing stage fits and repair only their missing intervals.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
setDTthreads(1)
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
if (dir.exists(out)) stop("Output must be new: ", out)
dir.create(out, recursive = TRUE)
stage <- file.path(root, "figures/candidates/pi-figure-redesign-2026-08-13-v8/analysis/stage_extensions")
pool <- file.path(root, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-POOLED-REPRO")
mm <- fread(file.path(pool, "model_input_manifest.tsv"))
dge_all <- readRDS(mm$dge_path)
stopifnot(identical(dim(dge_all), c(23370L, 844L)))
samples <- fread(file.path(stage, "stage_extension_sample_manifest.tsv"))
original <- fread(file.path(stage, "stage_extension_all_gene_results.tsv"))
stopifnot(!anyDuplicated(original, by = c("contrast", "gene_id_versioned")))
checks <- list(); repaired <- list()
for (cc in unique(original$contrast)) {
  ref <- unique(original[contrast == cc, reference]); cmp <- unique(original[contrast == cc, comparison])
  stopifnot(length(ref) == 1L, length(cmp) == 1L)
  sel <- copy(samples[contrast == cc]); setorder(sel, sample_id)
  stopifnot(!anyDuplicated(sel$sample_id), all(sel$sample_id %in% colnames(dge_all)))
  info <- data.frame(dataset = droplevels(factor(sel$dataset)), inferred_sex = droplevels(factor(sel$inferred_sex)),
    group = factor(sel$group, levels = c(ref, cmp)), row.names = sel$sample_id)
  design <- model.matrix(~ dataset + inferred_sex + group, info)
  stopifnot(qr(design)$rank == ncol(design))
  dge <- dge_all[, match(sel$sample_id, colnames(dge_all))]
  fit <- eBayes(lmFit(voomWithQualityWeights(dge, design, plot = FALSE), design))
  ci <- match(paste0("group", cmp), colnames(design))
  tt <- topTable(fit, coef = ci, number = Inf, sort.by = "none")
  se <- fit$stdev.unscaled[, ci] * sqrt(fit$s2.post)
  critical <- qt(.975, df = fit$df.total)
  names(critical) <- rownames(fit)
  old <- copy(original[contrast == cc]); g <- old$gene_id_versioned
  fields <- list(logFC = tt[g, "logFC"], SE = se[g], t = tt[g, "t"],
    P.Value = tt[g, "P.Value"], FDR = p.adjust(tt[g, "P.Value"], "BH"))
  err <- vapply(names(fields), function(nm) max(abs(old[[nm]] - fields[[nm]])), numeric(1))
  stopifnot(all(is.finite(err)), max(err) < 1e-8)
  # Keep the recorded effects/P/FDR bytes as values; fill only the confidence limits.
  old[, CI_low := logFC - critical[g] * SE]
  old[, CI_high := logFC + critical[g] * SE]
  stopifnot(all(is.finite(old$CI_low)), all(is.finite(old$CI_high)),
    all(old$CI_low <= old$logFC), all(old$CI_high >= old$logFC))
  checks[[length(checks) + 1L]] <- data.table(contrast = cc, n_samples = nrow(sel), n_genes = nrow(old),
    design_rank = qr(design)$rank, max_abs_reproduction_error = max(err), missing_intervals = 0L,
    max_abs_BH_error = max(abs(old$FDR - p.adjust(old$P.Value, "BH"))))
  repaired[[length(repaired) + 1L]] <- old
  cat(cc, ": ", nrow(old), " intervals; max effect/SE/t/P/FDR error ", max(err), "\n", sep = "")
}
result <- rbindlist(repaired)
stopifnot(nrow(result) == nrow(original), length(repaired) == 11L)
fwrite(result, file.path(out, "stage_extension_all_gene_results.tsv.gz"), sep = "\t", na = "NA")
fwrite(rbindlist(checks), file.path(out, "reproduction_checks.tsv"), sep = "\t")
paths <- c(mm$dge_path, file.path(pool, "model_input_manifest.tsv"),
  file.path(stage, "stage_extension_sample_manifest.tsv"), file.path(stage, "stage_extension_all_gene_results.tsv"))
hashes <- vapply(paths, function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), " ")[[1]][1], character(1))
fwrite(data.table(source = paths, sha256 = hashes), file.path(out, "input_hashes.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
