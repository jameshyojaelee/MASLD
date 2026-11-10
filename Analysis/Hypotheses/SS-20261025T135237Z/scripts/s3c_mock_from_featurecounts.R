# SS S3c reproduction mock: S2-layout input built from a control-free cohort of the 2026-08-28 transport check
# (GSE174478 or GSE240729), so that s3c_continuum_transport.R can be run on it and compared with
# cft-v1-20260828T141625Z (GSE174478 rho 0.510, n 93; GSE240729 rho 0.710, n 66).
# Writes <out>/mock_s2/{units.tsv, counts_all_genes.rds}, <out>/stage.tsv (S2 stage-file layout, codes 0-4 as
# recorded) and <out>/fractions.tsv (the 2026-08 BayesPrism proportions used by that run, S3e column layout).
# batch = one level; inferred_sex = recorded sex (GSE240729 has none: NA). Not a SteatoSITE input.
# Usage: Rscript s3c_mock_from_featurecounts.R <GSE174478|GSE240729> <new_out_dir>
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 2)
cid <- args[1]; out <- args[2]; stopifnot(cid %in% c("GSE174478", "GSE240729"), !dir.exists(out))
dir.create(file.path(out, "mock_s2"), recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
f <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results", cid, "counts/featurecounts/gene_counts.txt")
tb <- fread(f, skip = 1L, sep = "\t", header = TRUE)
mat <- as.matrix(tb[, 7:ncol(tb), with = FALSE]); storage.mode(mat) <- "integer"
colnames(mat) <- sub("\\.Aligned.*$", "", basename(colnames(mat))); rownames(mat) <- tb[[1]]
md <- fread(file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/pipelines/custom", cid, "metadata/SraRunTable.csv"))
if (cid == "GSE174478") {
  md <- data.table(sample_alias = md$Run, stage = as.numeric(md$fibrosis_stage),
                   inferred_sex = c(male = "M", female = "F")[tolower(md$sex)])
} else {
  md <- data.table(sample_alias = md$Run, stage = as.numeric(sub("^F", "", md$fibrosisscore)), inferred_sex = NA_character_)
}
smp <- intersect(colnames(mat), md$sample_alias); mat <- mat[, smp]; md <- md[match(smp, sample_alias)]
cat(cid, ": ", length(smp), " samples with counts and metadata, ", nrow(mat), " genes\n", sep = "")
lib <- colSums(mat)
units <- data.table(sample_alias = smp, analysis_unit_id = paste0("mock::", smp), batch = "B1",
                    inferred_sex = md$inferred_sex, assigned_fragments = lib, log10_assigned = log10(lib),
                    unique_pct = NA_real_, intron_exon_ratio = NA_real_, q1_pass = TRUE, q2_pass = TRUE, q3_pass = TRUE,
                    include = TRUE, exclusion_reason = NA_character_)
fwrite(units, file.path(out, "mock_s2/units.tsv"), sep = "\t")
saveRDS(mat, file.path(out, "mock_s2/counts_all_genes.rds"))
fwrite(data.table(sample_accession = smp, sample_alias = smp, nash_crn_kleiner_fibrosis_stage = md$stage),
       file.path(out, "stage.tsv"), sep = "\t", na = "NA")
pf <- file.path(ROOT, "Analysis/Deconvolution/results", cid, paste0(cid, "_bayesprism_proportions.tsv"))
hdr <- strsplit(readLines(pf, n = 1L), "\t", fixed = TRUE)[[1]]
pr <- fread(pf, skip = 1L, header = FALSE, sep = "\t"); stopifnot(ncol(pr) == length(hdr) + 1L)
setnames(pr, c("sample_alias", hdr))
pr <- cbind(pr[, .(sample_alias)], analysis_unit_id = paste0("mock::", pr$sample_alias), batch = "B1", pr[, -1])
fwrite(pr[sample_alias %in% smp], file.path(out, "fractions.tsv"), sep = "\t")
cat("DONE ->", out, "\n")
