#!/usr/bin/env Rscript
# R2: build two BG001-style run roots from the F_five arm inputs.
#   orig/ : the F_five merged_counts_raw.rds unchanged (must reproduce 1,347 DEGs)
#   corr/ : the same matrix with counts replaced ONLY for primary genes that have a
#           same-name copy on an ALT/patch/scaffold contig: multimapper-inclusive
#           (-M --fraction) count of the primary gene plus its ALT copies, rounded.
# All other genes and all other cohorts' columns are left unchanged.
# Usage: Rscript r2_build_roots.R <multi_counts_dir> <gtf> <out_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
multi_dir <- args[[1]]; gtf <- args[[2]]; out <- args[[3]]
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z/arms/F_five"
cohorts <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
primary_chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

m <- readRDS(file.path(arm, "results/integration/merged_counts_raw.rds"))
stopifnot(is.matrix(m), is.integer(m[1, 1]))

g <- fread(cmd = paste("zcat", shQuote(gtf), "| awk -F'\\t' '$3 == \"gene\"'"), header = FALSE, sep = "\t", quote = "")
g[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", V9)]
g[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", V9)]
g[, primary := V1 %in% primary_chr]
alt_names <- unique(g[primary == FALSE, gene_name])
prim <- g[primary == TRUE & gene_name %in% alt_names]
one_primary <- prim[, .N, by = gene_name][N == 1, gene_name]
fold <- g[gene_name %in% one_primary, .(gene_id, gene_name, primary)]
fold_primary <- fold[primary == TRUE]
cat("ALT-copy names:", length(alt_names), "| foldable (one primary gene):", length(one_primary), "\n")

corr <- m
report <- list()
for (c in cohorts) {
  f <- fread(file.path(multi_dir, paste0(c, ".multi.txt")), skip = "Geneid")
  runs <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(names(f)[-(1:6)]))
  setnames(f, c(names(f)[1:6], runs))
  keep_runs <- intersect(runs, colnames(m))
  mm <- as.matrix(f[, ..keep_runs]); rownames(mm) <- f$Geneid
  # same inputs and flags: genes without any ALT copy must match the original counts closely
  other <- setdiff(intersect(rownames(mm), rownames(m)), fold$gene_id)
  tot_orig <- rowSums(m[other, keep_runs, drop = FALSE]); tot_multi <- rowSums(mm[other, , drop = FALSE])
  ok <- tot_orig >= 100
  ratio_other <- median(tot_multi[ok] / tot_orig[ok])
  # fold ALT copies into their primary gene
  ids <- intersect(fold$gene_id, rownames(mm))
  by_name <- rowsum(mm[ids, , drop = FALSE], fold[match(ids, gene_id), gene_name])
  tgt <- fold_primary[gene_id %in% rownames(m) & gene_name %in% rownames(by_name)]
  new_vals <- round(by_name[tgt$gene_name, keep_runs, drop = FALSE])
  before <- rowSums(m[tgt$gene_id, keep_runs, drop = FALSE])
  corr[tgt$gene_id, keep_runs] <- matrix(as.integer(new_vals), nrow = nrow(new_vals))
  after <- rowSums(corr[tgt$gene_id, keep_runs, drop = FALSE])
  report[[c]] <- data.table(cohort = c, gene_id = tgt$gene_id, gene_name = tgt$gene_name,
                            count_before = before, count_after = after)
  cat(c, ": runs", length(keep_runs), "| genes replaced", nrow(tgt),
      "| median multi/orig for non-ALT genes (>=100 reads):", round(ratio_other, 4), "\n")
}
rep <- rbindlist(report)
unchanged_other <- all(corr[setdiff(rownames(m), fold_primary$gene_id), ] == m[setdiff(rownames(m), fold_primary$gene_id), ])
cat("non-ALT rows identical to original:", unchanged_other, "\n")
stopifnot(unchanged_other)

meta <- file.path(arm, "results/integration/meta_matched.rds")
qc <- file.path(arm, "qc/sample_qc_report.csv")
for (label in c("orig", "corr")) {
  root <- file.path(out, label)
  dir.create(file.path(root, "results/integration"), recursive = TRUE)
  dir.create(file.path(root, "qc"))
  writeLines(paste(basename(out), label, sep = "\t"), file.path(root, ".bg001_candidate_root"))
  file.copy(meta, file.path(root, "results/integration/meta_matched.rds"))
  file.copy(qc, file.path(root, "qc/sample_qc_report.csv"))
  saveRDS(if (label == "orig") m else corr, file.path(root, "results/integration/merged_counts_raw.rds"))
}
fwrite(rep, file.path(out, "replaced_gene_counts.tsv.gz"), sep = "\t")
cat("roots written under", out, "\n")
