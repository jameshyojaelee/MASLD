#!/usr/bin/env Rscript
# Model A composition: InstaPrism initial posterior (fraction of bulk mRNA) for every
# library of the 9 cohorts in the kallisto matrix, reproducing the code path of
# deconv_estimand_harmonization-20260828T141918Z (same GTF symbol map, same
# aggregation and rounding, same ScaleSC reference, Post.ini.ct@theta).
# Guard: the 521 libraries already in that table must reproduce to 1e-6.
suppressPackageStartupMessages({
  library(data.table); library(SingleCellExperiment); library(SummarizedExperiment)
  library(Matrix); library(InstaPrism)
})
set.seed(20260828)
OUT <- Sys.getenv("COMP_OUT"); stopifnot(nzchar(OUT))
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
KALL <- file.path(BASE, "RNA-seq/results/kallisto_sensitivity/all_cohorts_gene_counts.tsv.gz")
GTF  <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
REF  <- file.path(BASE, "Analysis/Deconvolution/reference/reference_human_scalesc_sce.rds")
PRIOR <- file.path(BASE, "Analysis/MASLD_Model_Benchmark/executions/deconv_estimand_harmonization-20260828T141918Z/intermediate")
BAI  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/a1-full-20260924T000606Z/bai")
cohorts <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523",
             "GSE174478", "GSE193066", "GSE213621", "GSE240729")

gtf_lines <- fread(cmd = paste("zcat", GTF, "| grep -P '^chr.*\\tgene\\t'"), header = FALSE, sep = "\t")
ens2sym <- unique(data.table(ensembl = sub("\\..*$", "", gsub('.*gene_id "([^"]+)".*', "\\1", gtf_lines$V9)),
                             symbol = gsub('.*gene_name "([^"]+)".*', "\\1", gtf_lines$V9)), by = "ensembl")
lib <- rbindlist(lapply(cohorts, function(c)
  data.table(cohort = c, run = sub("\\.bam$", "", basename(Sys.glob(file.path(BAI, c, "*.bam")))))))
kall_hdr <- names(fread(cmd = paste("zcat", KALL, "| head -1"), nrows = 0))
miss <- setdiff(lib$run, kall_hdr)
if (length(miss)) stop(sprintf("%d libraries absent from kallisto matrix", length(miss)))
kc <- fread(cmd = paste("zcat", KALL), select = c("gene_id", lib$run))
sym <- ens2sym$symbol[match(sub("\\..*$", "", kc$gene_id), ens2sym$ensembl)]
kc[, gene_id := NULL]
kc <- data.table(SYMBOL = sym, kc)[!is.na(SYMBOL)]
kc <- kc[, lapply(.SD, sum), by = SYMBOL, .SDcols = lib$run]
for (col in lib$run) set(kc, j = col, value = as.integer(round(kc[[col]])))
cat(sprintf("bulk: %d symbols x %d libraries\n", nrow(kc), length(lib$run)))

sce <- readRDS(REF)
ref_counts <- assay(sce, "counts"); if (is.null(ref_counts)) ref_counts <- assay(sce, 1)
ct <- as.character(colData(sce)[["celltype"]])
ref_obj <- refPrepare(sc_Expr = as.matrix(ref_counts), cell.type.labels = ct, cell.state.labels = ct)
rm(sce, ref_counts); gc()

props <- list()
for (c in cohorts) {
  runs <- lib[cohort == c, run]
  bulk <- as.matrix(kc[, c("SYMBOL", runs), with = FALSE], rownames = "SYMBOL")
  res <- InstaPrism(bulk_Expr = bulk, refPhi_cs = ref_obj)
  p <- t(res@Post.ini.ct@theta)
  stopifnot(nrow(p) == length(runs))
  props[[c]] <- data.table(cohort = c, run = rownames(p), p)
  cat(c, nrow(p), "libraries\n")
}
P <- rbindlist(props, use.names = TRUE)
fwrite(P, file.path(OUT, "instaprism_initial_mrna_fraction_all_libraries.tsv"), sep = "\t")

# reproduction guard against the deposited 5-cohort initial posteriors
dev <- c()
for (f in Sys.glob(file.path(PRIOR, "*_kallisto_instaprism_INITIAL_proportions.tsv"))) {
  q <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  m <- as.matrix(P[match(rownames(q), P$run), colnames(q), with = FALSE])
  stopifnot(!anyNA(m))
  dev[basename(f)] <- max(abs(m - q))
}
print(dev)
writeLines(sprintf("%s\t%.3e", names(dev), dev), file.path(OUT, "reproduction_vs_deposited_max_abs_diff.tsv"))
if (any(dev > 1e-6)) stop("deposited initial posteriors not reproduced")
sink(file.path(OUT, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("done\n")
