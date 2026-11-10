# SS S3h check on simulated counts (no real data): 400 individuals, 9 batches, stage codes 1-4; 30 "catalog" sites
# (5 with ALT logit fraction +0.5 per stage step, 25 without) and 300 null sites, all with a reference bias of
# logit -0.08 and beta-binomial overdispersion (rho 0.02). Writes S2-like inputs, then s3h_stage_ase.R must recover
# the 5 planted genes at BH < 0.05, call at most 2 of the 25 others, and give lambda_null near 1.
# Usage: Rscript s3h_synthetic_test.R <new_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1])); w <- args[1]
dir.create(file.path(w, "in"), recursive = TRUE); dir.create(file.path(w, "counts")); dir.create(file.path(w, "targets"))
set.seed(20261025)
n <- 400
U <- data.table(sample_alias = sprintf("S%04d", 1:n), batch = sample(paste0("b", 1:9), n, TRUE), inferred_sex = sample(c("F", "M"), n, TRUE), include = TRUE)
fwrite(U, file.path(w, "in/units.tsv"), sep = "\t")
fwrite(data.table(sample_accession = U$sample_alias, sample_alias = U$sample_alias, nash_crn_kleiner_fibrosis_stage = sample(1:4, n, TRUE, c(0.45, 0.15, 0.2, 0.2))),
       file.path(w, "in/stage.tsv"), sep = "\t")
stg <- as.numeric(fread(file.path(w, "in/stage.tsv"))$nash_crn_kleiner_fibrosis_stage)
TG <- data.table(chrom = "chr1", pos = seq(1e6, by = 1000, length.out = 330), ref = "A", alt = "G",
                 role = rep(c("s3h", "s3h_null"), c(30, 300)), gene = c(paste0("CAT", 1:30), paste0("NUL", 1:300)),
                 id = c(sprintf("ENSG%011d", 1:30), sprintf("ENSG9%010d", 1:300)), af_eur = runif(330, 0.08, 0.5), note = "")
fwrite(TG, file.path(w, "targets/s1b_targets_annotated.tsv"), sep = "\t")
writeLines(system2("sha256sum", normalizePath(file.path(w, "targets/s1b_targets_annotated.tsv")), stdout = TRUE), file.path(w, "targets/SHA256SUMS"))
rows <- list()
for (k in seq_len(nrow(TG))) {
  het <- runif(n) < 2 * TG$af_eur[k] * (1 - TG$af_eur[k]); dp <- rnbinom(n, mu = 40, size = 3)
  eff <- if (k <= 5) 0.5 * (stg - 2.5) else 0
  pr <- plogis(-0.08 + eff); pi <- rbeta(n, pr * (1 / 0.02 - 1), (1 - pr) * (1 / 0.02 - 1))
  nalt <- ifelse(het, rbinom(n, dp, pi), rbinom(n, dp, 0.003))
  rows[[k]] <- data.table(sample_alias = U$sample_alias, chrom = "chr1", pos = TG$pos[k], nref = dp - nalt, nalt = nalt)[dp > 0]
}
fwrite(rbindlist(rows), file.path(w, "counts/ase_counts.tsv.gz"), sep = "\t")
RS <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
src <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3h_stage_ase.R"
stopifnot(system2(RS, c(src, "--input", file.path(w, "in"), "--stage", file.path(w, "in/stage.tsv"), "--counts", file.path(w, "counts"),
                        "--targets", file.path(w, "targets"), "--out", file.path(w, "out"))) == 0)
F <- fread(file.path(w, "out/f_ss8_catalog_genes.tsv"))
hit <- F[bh < 0.05, unit]; planted <- sprintf("ENSG%011d", 1:5)
cat("planted recovered:", sum(planted %in% hit), "of 5; false calls:", sum(!hit %in% planted), "of 25\n")
cat(readLines(file.path(w, "out/summary.txt")), sep = "\n")
stopifnot(all(planted %in% hit), sum(!hit %in% planted) <= 2)
cat("S3h synthetic test PASSED\n")
