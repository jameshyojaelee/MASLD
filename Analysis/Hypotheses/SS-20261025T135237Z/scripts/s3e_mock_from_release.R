# SS S3e reproduction mock: one 774-release cohort written in the S2 layout (S2_CONTRACT.md), so s3e_composition.R
# can be run on data whose fractions, continuum position and IgA counts are already stored. batch = one level
# ("B1"); sample_alias = sample_id; stage file = release fibrosis_stage (0-4, NA for controls), written with CRLF
# like the SteatoSITE file. dge_included.rds = release merged_dge columns; counts_all_genes.rds = the 79,095-gene
# corrected counts behind the stored fractions. No SteatoSITE file is read.
# Usage: Rscript s3e_mock_from_release.R <dataset> <new_out_dir>
suppressPackageStartupMessages({library(edgeR); library(data.table)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 2)
ds <- args[1]; out <- args[2]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
RAW <- file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")
d <- readRDS(file.path(REL, "merged_dge.rds")); m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
m <- m[dataset == ds]; ids <- m$sample_id; stopifnot(length(ids) > 0, all(ids %in% colnames(d)))
raw <- readRDS(RAW); stopifnot(nrow(raw) == 79095L, all(ids %in% colnames(raw)))
u <- data.table(sample_alias = ids, analysis_unit_id = m$analysis_unit_id, batch = "B1", inferred_sex = m$inferred_sex,
                assigned_fragments = colSums(raw[, ids]), unique_pct = NA_real_, intron_exon_ratio = NA_real_,
                q1_pass = TRUE, q2_pass = TRUE, q3_pass = TRUE, include = TRUE, exclusion_reason = NA_character_)
u[, log10_assigned := log10(assigned_fragments)]
fwrite(u, file.path(out, "units.tsv"), sep = "\t")
dg <- DGEList(d$counts[, ids]); dg$samples <- cbind(dg$samples[, c("group", "lib.size", "norm.factors")], as.data.frame(u))
saveRDS(dg, file.path(out, "dge_included.rds"))
saveRDS(raw[, ids], file.path(out, "counts_all_genes.rds"))
st <- data.table(sample_accession = ids, sample_alias = ids, nash_crn_kleiner_fibrosis_stage = m$fibrosis_stage)
writeLines(paste0(c(paste(names(st), collapse = "\t"),
                    st[, paste(sample_accession, sample_alias, ifelse(is.na(nash_crn_kleiner_fibrosis_stage), "NA", nash_crn_kleiner_fibrosis_stage), sep = "\t")]), "\r"),
           file.path(out, "mock_stage.tsv"))
cat(ds, length(ids), "units written to", out, "\n")
