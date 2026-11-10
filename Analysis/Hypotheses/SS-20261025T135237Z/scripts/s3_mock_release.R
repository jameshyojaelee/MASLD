# SS reproduction mock (S2_CONTRACT layout) built from the 774-unit release, for checking the S3f/S3g/S3h scripts
# against DIG1. No SteatoSITE data. units.tsv: sample_alias = release sample_id, batch = one level ("B1"),
# inferred_sex = sex_final, include = TRUE for all 774, dataset kept (used only with --dataset-term).
# dge_included.rds = release merged_dge.rds (24,115 genes). genotypes/<id>.vcf.gz = links to the DIG1 job 62 VCFs
# (job 62 sites only). sj_paths.tsv = release STAR SJ.out.tab. Two stage files (columns as the ENA file):
#   stage_all.tsv      fibrosis_stage of every unit with one (job 63 call-validity set);
#   stage_disease4.tsv fibrosis_stage of disease units in GSE130970, GSE135251, GSE162694, GSE213621 (job 41a set).
# Usage: Rscript s3_mock_release.R <new_out_dir>
suppressPackageStartupMessages({library(edgeR); library(data.table)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1])); out <- args[1]
dir.create(file.path(out, "genotypes"), recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
D62 <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_rna_genotypes/vcf")
B <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results")
d <- readRDS(file.path(REL, "merged_dge.rds")); m <- readRDS(file.path(REL, "meta_matched.rds"))
m <- m[match(colnames(d$counts), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(d$counts)))
U <- m[, .(sample_alias = sample_id, analysis_unit_id, batch = "B1", inferred_sex = sex_final, include = TRUE, dataset)]
fwrite(U, file.path(out, "units.tsv"), sep = "\t")
d$samples <- cbind(d$samples[, intersect(c("group", "lib.size", "norm.factors"), names(d$samples))], as.data.frame(U))
saveRDS(d, file.path(out, "dge_included.rds"))
for (s in U$sample_alias) { f <- file.path(D62, paste0(s, ".vcf.gz")); if (file.exists(f)) file.symlink(f, file.path(out, "genotypes", paste0(s, ".vcf.gz"))) }
fwrite(U[, .(sample_alias, path = file.path(B, dataset, "alignments/star", sample_alias, paste0(sample_alias, ".SJ.out.tab")))],
       file.path(out, "sj_paths.tsv"), sep = "\t")
st <- function(sel) data.table(sample_accession = m$sample_id, sample_alias = m$sample_id,
                               nash_crn_kleiner_fibrosis_stage = ifelse(sel & !is.na(m$fibrosis_stage), m$fibrosis_stage, NA))
fwrite(st(rep(TRUE, nrow(m))), file.path(out, "stage_all.tsv"), sep = "\t", na = "NA")
fwrite(st(m$group_binary == "Disease" & m$dataset %in% c("GSE130970", "GSE135251", "GSE162694", "GSE213621")),
       file.path(out, "stage_disease4.tsv"), sep = "\t", na = "NA")
cat("units", nrow(U), "vcf links", length(list.files(file.path(out, "genotypes"))), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
