#!/usr/bin/env Rscript
# F1: program scores for Model F (PRESPEC_F sections 2-3). Scoring only; no outcome
# statistic is computed here. Human: four staged cohorts of the five-cohort bulk set.
# Mouse: the 17 prespecified diet-versus-lean-control contrasts, mouse genes renamed to
# the human one-to-one ortholog symbol. Same pinned scorer on both sides
# (collapse_symbols -> score_programs, coverage 0.80), z-scored within cohort or contrast.
# Usage: Rscript f1_score_programs.R <out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
out <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
here <- file.path(repo, "scripts/analysis/preclinical_fidelity")
owd <- setwd(here); stopifnot(system2("sha256sum", c("-c", "--quiet", "PRESPEC_F.sha256")) == 0); setwd(owd)
source(file.path(repo, "scripts/manuscript/program_observability_map/analysis_lib.R"))
rel <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs")
hot <- file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
registry <- fread(file.path(hot, "program_registry_v2.tsv")); membership <- fread(file.path(hot, "program_membership_v2.tsv"))

score <- function(lib, annotation, id) {
  sm <- collapse_symbols(lib, rownames(lib), annotation)
  s <- score_programs(sm, data.table(sample_id = colnames(sm), dataset = id), membership, registry, 0.80)
  list(scores = t(s$primary), coverage = cbind(id = id, s$coverage))
}
logcpm <- function(counts) {  # genes with CPM >= 1 in at least 10% of samples
  d <- calcNormFactors(DGEList(round(counts)), method = "TMM")
  cpm(d, log = TRUE, prior.count = 1)[rowMeans(cpm(d) >= 1) >= 0.10, , drop = FALSE]
}
logrpkm <- function(rpkm) log2(rpkm[rowMeans(rpkm >= 1) >= 0.10, , drop = FALSE] + 1)

# ---- human ----
h_ann <- fread(file.path(rel, "BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"), select = c("gene_id", "gene_name"))
dge <- readRDS(file.path(rel, "BG001-DECISION/arms/F_five/results/integration/merged_dge.rds"))
man <- fread(file.path(repo, "figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/stage_extensions/five_cohort_sample_manifest.tsv"))
man <- man[dataset %in% c("GSE130970", "GSE135251", "GSE162694", "GSE213621") & !is.na(fibrosis_stage)]
stopifnot(all(man$sample_id %in% colnames(dge$counts)))
cov_rows <- list(); human <- list()
for (ds in unique(man$dataset)) {
  m <- man[dataset == ds]
  r <- score(logcpm(dge$counts[, m$sample_id]), h_ann, ds)
  human[[ds]] <- cbind(m[, .(sample_id, dataset, fibrosis_stage, inferred_sex)], r$scores[m$sample_id, ])
  cov_rows[[ds]] <- r$coverage
}
fwrite(rbindlist(human), file.path(out, "human_program_scores.tsv"), sep = "\t")

# ---- mouse ----
orth <- fread(file.path(repo, "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"))[ortholog_type == "ortholog_one2one"]
m_ann <- unique(orth[, .(gene_id = mouse_gene_ensembl, gene_name = human_gene_symbol)])
m_ann <- m_ann[!duplicated(gene_id) & !duplicated(gene_name)]
M <- file.path(repo, "RNA-seq/Mouse"); G <- file.path(repo, "RNA-seq/results/glp1ra/mouse_reversal/_inputs")
read_fc <- function(paths) {  # featureCounts tables; sample = BAM file stem
  ans <- NULL
  for (p in paths[file.exists(paths)]) {
    x <- fread(p, skip = "Geneid"); m <- as.matrix(x[, -(1:6)]); rownames(m) <- x$Geneid
    colnames(m) <- sub("\\.Aligned.*$", "", basename(colnames(m)))
    ans <- if (is.null(ans)) m else cbind(ans, m[rownames(ans), setdiff(colnames(m), colnames(ans)), drop = FALSE])
  }
  ans
}
read_mat <- function(p) {  # header read as text: numeric animal IDs are otherwise taken for a data row
  h <- gsub('"', "", strsplit(readLines(p, n = 1), "\t", fixed = TRUE)[[1]])
  x <- fread(p, skip = 1, header = FALSE)
  if (length(h) == ncol(x) - 1) h <- c("gene", h)
  stopifnot(length(h) == ncol(x))
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; colnames(m) <- h[-1]; m
}
pdm <- fread(file.path(M, "Public_Diet_Models/metadata/samples.tsv"))[, .(sample = sample_id, condition, sex, dataset)]
pdm_fc <- file.path(M, "Public_Diet_Models/counts/featurecounts", c("gene_counts.txt", "gene_counts.txt.pe", "gene_counts.txt.se"))
uni <- fread(file.path(M, "Unified_Integration/metadata/unified_mouse_metadata.csv"))
wd <- function(g, idcol) {  # Western_Diet_Datasets sample table with SRR run joined by GSM
  m <- fread(file.path(M, "Western_Diet_Datasets", g, "metadata/sample_metadata.csv"))
  if ("srr" %in% names(m)) return(m[, sample := srr])
  map <- fread(file.path(M, "Western_Diet_Datasets", g, "metadata/gsm_to_srr.tsv"))
  m[, sample := map$srr[match(get(idcol), map$gsm)]]
}
wd_fc <- function(g) file.path(M, "Western_Diet_Datasets", g, "counts/featurecounts/gene_counts.txt")
arms <- function(d, case, control, cov = NULL) {
  d <- d[, arm := fifelse(case, "case", fifelse(control, "control", NA_character_))][!is.na(arm)]
  d[, c("sample", "arm", cov), with = FALSE]
}
u <- function(g) uni[dataset == g][, sample := sample_id]
p <- function(g) pdm[dataset == g]
gl <- function(g) fread(file.path(G, paste0("coldata_", g, ".tsv")))[, sample := as.character(sample_id)]
x4 <- p("GSE274914")[, timepoint := sub("_.*$", "", condition)]
w88 <- wd("GSE246088", "geo_accession"); w84 <- wd("GSE305484", "gsm_accession")
contrasts <- list(
  c01 = list("GSE156918", "MCD", file.path(M, "Unified_Integration/counts/featurecounts/gene_counts_gse156918.txt"), "fc",
             with(u("GSE156918"), arms(u("GSE156918"), group_binary == "Disease", group_binary == "Control"))),
  c02 = list("GSE205974", "MCD", file.path(M, "Unified_Integration/counts/featurecounts/gene_counts_gse205974.txt"), "fc",
             with(u("GSE205974"), arms(u("GSE205974"), group_binary == "Disease", group_binary == "Control"))),
  c03 = list("GSE224069", "HFD", pdm_fc, "fc",
             with(p("GSE224069"), arms(p("GSE224069"), condition == "High_Fat_Diet_Vehicle", condition == "Chow_Diet_Vehicle"))),
  c04 = list("GSE274914", "HFD", pdm_fc, "fc", with(x4, arms(x4, grepl("_HFD$", condition), grepl("_LFD$", condition), c("timepoint", "sex")))),
  c05 = list("GSE246088", "HFD", wd_fc("GSE246088"), "fc", with(w88, arms(copy(w88), condition == "HFD_Control", condition == "Chow_Control"))),
  c06 = list("GSE246088", "WESTERN", wd_fc("GSE246088"), "fc", with(w88, arms(copy(w88), condition == "WD_Control", condition == "Chow_Control"))),
  c07 = list("GSE162876", "CD", pdm_fc, "fc", with(p("GSE162876"), arms(p("GSE162876"), condition == "7wks_CDAHFD", condition == "7wks_LFD_Ctrl"))),
  c08 = list("GSE294630", "CD", file.path(G, "GSE294630_processed_file.tsv.gz"), "mat",
             with(gl("GSE294630"), arms(gl("GSE294630"), group == "CDAA_vehicle_QD_6_18w", group == "Chow_18w"))),
  c09 = list("GSE162876", "WESTERN", pdm_fc, "fc", with(p("GSE162876"), arms(p("GSE162876"), condition == "20wks_FPC", condition == "20wks_LFD_Ctrl"))),
  c10 = list("GSE220575", "WESTERN", wd_fc("GSE220575"), "fc",
             {d <- wd("GSE220575", "gsm"); with(d, arms(d, condition == "MASH", condition == "Control"))}),
  c11 = list("GSE292565", "WESTERN", wd_fc("GSE292565"), "fc",
             {d <- wd("GSE292565", "gsm"); with(d, arms(d, condition == "Disease", condition == "Control"))}),
  # c12: local featureCounts covers 4 of 5 WT chow animals; the submitter raw counts cover all (PRESPEC_F amendment 1)
  c12 = list("GSE305484", "WESTERN", file.path(M, "Western_Diet_Datasets/GSE305484/counts/GSE305484_GEO_CD163KO_rawcounts.txt.gz"), "mat",
             with(w84, arms(copy(w84)[, sample := sample_id], genotype == "CD163WT" & condition == "WD", genotype == "CD163WT" & condition == "Chow"))),
  c13 = list("GSE196908", "WESTERN", file.path(G, "GSE196908_rpkm.txt.gz"), "rpkm",
             with(gl("GSE196908"), arms(gl("GSE196908"), group == "DIONASH_Vehicle" & timepoint %in% "12w", group == "Chow_Vehicle" & timepoint %in% "12w"))),
  c14 = list("GSE226496", "WESTERN", file.path(G, "GSE226496_merged.counts.tsv.gz"), "mat",
             with(gl("GSE226496"), arms(gl("GSE226496"), group == "F37", group == "C37"))),
  c15 = list("GSE243976", "WESTERN", file.path(G, "GSE243976_GAN-DIO-NASH-HCC_counts.txt.gz"), "mat",
             with(gl("GSE243976"), arms(gl("GSE243976"), group == "Parenchyme", group == "Chow"))),
  c16 = list("GSE294629", "WESTERN", file.path(G, "GSE294629_processed_file.tsv.gz"), "mat",
             with(gl("GSE294629"), arms(gl("GSE294629"), group == "Vehicle_24_w", group == "Chow_Vehicle_24_w"))),
  c17 = list("GSE159911", "LIDPAD", pdm_fc, "fc", with(p("GSE159911"), arms(p("GSE159911"), condition == "LIDPAD", condition == "Control")))
)
census <- list(); mouse <- list(); cache <- list()
for (id in names(contrasts)) {
  cc <- contrasts[[id]]; s <- cc[[5]]; key <- paste(cc[[3]], collapse = "|")
  if (is.null(cache[[key]])) cache[[key]] <- if (cc[[4]] == "fc") read_fc(cc[[3]]) else read_mat(cc[[3]])
  counts <- cache[[key]]
  stopifnot(!anyNA(s$sample), all(s$sample %in% colnames(counts)), min(table(s$arm)) >= 3)
  lib <- if (cc[[4]] == "rpkm") logrpkm(counts[, s$sample]) else logcpm(counts[, s$sample])
  r <- score(lib, m_ann, id)
  for (v in c("timepoint", "sex")) if (!v %in% names(s)) s[, (v) := NA_character_]
  mouse[[id]] <- cbind(s[, .(contrast = id, series = cc[[1]], diet_class = cc[[2]], sample, arm, timepoint, sex)], r$scores[s$sample, ])
  cov_rows[[id]] <- r$coverage
  census[[id]] <- data.table(contrast = id, series = cc[[1]], diet_class = cc[[2]], input = cc[[4]],
                             n_case = sum(s$arm == "case"), n_control = sum(s$arm == "control"))
}
fwrite(rbindlist(mouse), file.path(out, "mouse_program_scores.tsv"), sep = "\t")
fwrite(rbindlist(census), file.path(out, "mouse_contrast_census.tsv"), sep = "\t")
fwrite(rbindlist(cov_rows, fill = TRUE), file.path(out, "program_coverage.tsv"), sep = "\t")
print(rbindlist(census))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_f1_score_programs.txt"))
