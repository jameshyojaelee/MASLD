# SS S3a verification mock (not a result). Builds an S2-layout input from the 774-unit release so that s3a_replication.R
# can reproduce H51 G0 (`~ inferred_sex + bin` on the 258 GSE213621 disease units, strata accession quartile x sex).
# batch = one level ("mock"); acc_q = H51 accession quartile (strata column); stage code from GSE213621 `condition`
# (F0F1 -> 1, F2 -> 2, F3F4 -> 3; no code 4). Normalization factors are the stored release factors (--norm stored).
# Quality columns: log10_assigned from stored library size; unique_pct and intron_exon_ratio are synthetic (seed 1),
# only to exercise the quality-covariate code path. ILR: H51 7-part lineage coding of the corrected-774 BayesPrism
# fractions (floor 1e-6, Helmert basis), only to exercise the composition path.
# Usage: Rscript s3a_build_mock.R --out <new dir>
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
A <- parse_args(commandArgs(trailingOnly = TRUE), "out")
OUT <- new_outdir(A$out)
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
DEC <- file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/results")
stopifnot(sha256(file.path(REL, "merged_dge.rds")) == "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b",
          sha256(file.path(REL, "meta_matched.rds")) == "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3")
dge_all <- readRDS(file.path(REL, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
g <- meta[dataset == "GSE213621" & group_binary == "Disease"]
stopifnot(nrow(g) == 258L, !anyNA(g$inferred_sex))
setorder(g, sample_id)
code_of <- c(Fibrosis_F0F1 = 1L, Fibrosis_F2 = 2L, Fibrosis_F3F4 = 3L)
g[, stage_code := unname(code_of[condition])]
stopifnot(!anyNA(g$stage_code))
g[, acc_rank := frank(as.numeric(sub("SRR", "", sample_id)))]
g[, acc_q := as.character(cut(acc_rank, quantile(acc_rank, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4)))]

d <- dge_all[, g$sample_id]
set.seed(1L)
units <- data.table(sample_alias = g$sample_id, analysis_unit_id = g$sample_id, batch = "mock", inferred_sex = g$inferred_sex,
  assigned_fragments = d$samples$lib.size, log10_assigned = log10(d$samples$lib.size),
  unique_pct = round(80 + rnorm(nrow(g), 0, 4), 2), intron_exon_ratio = round(0.3 + abs(rnorm(nrow(g), 0, 0.1)), 4),
  q1_pass = TRUE, q2_pass = TRUE, q3_pass = TRUE, include = TRUE, exclusion_reason = "", acc_q = g$acc_q)
d$samples <- cbind(d$samples[, c("group", "lib.size", "norm.factors")], as.data.frame(units))
saveRDS(d, file.path(OUT, "dge_included.rds"))
fwrite(units, file.path(OUT, "units.tsv"), sep = "\t")
st <- c("sample_accession\tsample_alias\tnash_crn_kleiner_fibrosis_stage", paste0("MOCK", seq_len(nrow(g)), "\t", g$sample_id, "\t", g$stage_code))
writeBin(charToRaw(paste0(paste(st, collapse = "\r\n"), "\r\n")), file.path(OUT, "stage_mock.tsv"))

LINEAGE <- list(hep_chol = c("Hepatocytes", "Cholangiocytes"), endothelial = "Endothelial cells", fibroblasts = "Fibroblasts",
                mac_mono = c("Macrophages", "Mono+mono derived cells"), other_myeloid = c("cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"),
                t_nk = c("T cells", "Circulating NK/NKT", "Resident NK"), b_plasma = c("B cells", "Plasma cells"))
f <- read.delim(file.path(DEC, "GSE213621/GSE213621_bayesprism_proportions.tsv"), check.names = FALSE)
stopifnot(all(g$sample_id %in% rownames(f)))
f <- as.matrix(f[g$sample_id, , drop = FALSE])
P <- sapply(LINEAGE, function(cols) rowSums(f[, cols, drop = FALSE])); P <- P / rowSums(P)
P <- pmax(P, 1e-6); P <- P / rowSums(P)
L <- log(P); clr <- L - rowMeans(L)
V <- matrix(0, 7, 6); for (k in 1:6) { V[1:k, k] <- -1 / sqrt(k * (k + 1)); V[k + 1, k] <- k / sqrt(k * (k + 1)) }
ilr <- clr %*% V; colnames(ilr) <- paste0("ilr", 1:6)
fwrite(data.table(sample_alias = g$sample_id, ilr), file.path(OUT, "ilr_mock.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("Mock written:", OUT, "units", nrow(units), "; stage", paste(names(table(g$stage_code)), table(g$stage_code), collapse = " "))
