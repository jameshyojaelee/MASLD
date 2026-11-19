# SS job S3a LOCO (2026-11-18): Composition-adjusted SteatoSITE replication using
# 22-autosome Leave-One-Chromosome-Out (LOCO) separated 8-part ILR coordinates (ilr1..ilr7).
# Reproduces frozen FIT$bin bit-for-bit on the 493 staged units, computes standard 8-part ILR
# (comp_standard_8part) and chr1..chr22 held-out fits, and sets FIT$comp per gene chromosome.

SS_SCRIPTS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts"
source(file.path(SS_SCRIPTS, "s3a_common.R"))
suppressPackageStartupMessages({
  library(parallel)
  library(jsonlite)
})

A <- parse_args(commandArgs(TRUE), required = character(), defaults = list())
OUT <- A$out
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))

IN_DIR <- if (!is.null(A$input)) A$input else file.path(SS, "out_s2_qc")
STAGE_TSV <- if (!is.null(A$stage)) A$stage else file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv")
LOCO_CAND <- file.path(ROOT, "Analysis/Deconvolution/candidates/loco_20261118")
ILR_STD <- if (!is.null(A$ilr_std)) A$ilr_std else file.path(LOCO_CAND, "standard_steatosite_ilr8.tsv")
LOCO_DIR <- if (!is.null(A$loco_dir)) A$loco_dir else file.path(LOCO_CAND, "chr")
REF_S3A <- if (!is.null(A$ref_s3a)) A$ref_s3a else file.path(SS, "out_s3a-20261026T234211Z")

stopifnot(dir.exists(IN_DIR), file.exists(STAGE_TSV), file.exists(ILR_STD), dir.exists(LOCO_DIR), dir.exists(REF_S3A))

ann <- load_ann()
units <- fread(file.path(IN_DIR, "units.tsv"))
dge <- readRDS(file.path(IN_DIR, "dge_included.rds"))
if (!"group" %in% names(dge$samples)) dge$samples <- cbind(group = factor(rep(1L, ncol(dge))), dge$samples)
stopifnot(nrow(dge) == 24115L, setequal(rownames(dge), ann$gene_id_versioned), !anyDuplicated(units$sample_alias))
ann <- ann[match(rownames(dge), gene_id_versioned)]

st <- read_stage(STAGE_TSV)
u <- merge(units[include == TRUE], st, by = "sample_alias", all.x = TRUE)
u <- u[!is.na(stage_code) & !is.na(inferred_sex)]
setorder(u, sample_alias)
u <- add_labels(u, u$stage_code)
fwrite(u, file.path(OUT, "units_used.tsv"), sep = "\t")

dg_all <- dge[, u$sample_alias]
keep_l <- filterByExpr(dg_all, group = u$bin)
KEEP <- rownames(dg_all)[keep_l]
dg <- calcNormFactors(dg_all[keep_l, , keep.lib.sizes = TRUE], method = "RLE")
log_msg("Staged units:", nrow(u), "KEEP genes:", length(KEEP))

# Run BIN_CHECK, STD8, and chr1..chr22 fits in parallel across NCORES (all deterministic, zero RNG)
ref_fit <- readRDS(file.path(REF_S3A, "observed_fits.rds"))
pack <- function(f, cfs) list(coef = f$coefficients[, cfs, drop = FALSE], tab = gene_tab(f, cfs))
cfs <- c("binF2", "binF3F4")
chrs <- paste0("chr", 1:22)
tasks <- c("BIN_CHECK", "STD8", chrs)

run_ilr_fit <- function(path) {
  il <- fread(path)
  icol <- grep("^ilr[0-9]+$", names(il), value = TRUE)
  stopifnot("sample_alias" %in% names(il), length(icol) == 7L, all(u$sample_alias %in% il$sample_alias))
  ilr_dt <- as.data.table(il[match(u$sample_alias, sample_alias), ..icol])
  info <- cbind(u[, .(batch, inferred_sex, bin)], ilr_dt)
  fit_qw(dg, design_of(info, "bin", icol))
}

log_msg(sprintf("Fitting BIN_CHECK + STD8 + 22 LOCO chromosomes (%d tasks) across %d cores...", length(tasks), NCORES))
raw_fits <- mclapply(tasks, function(tk) {
  if (tk == "BIN_CHECK") return(fit_qw(dg, design_of(u, "bin")))
  if (tk == "STD8") return(run_ilr_fit(ILR_STD))
  run_ilr_fit(file.path(LOCO_DIR, paste0(tk, "_steatosite_ilr8.tsv")))
}, mc.cores = min(NCORES, length(tasks)), mc.preschedule = FALSE)
names(raw_fits) <- tasks
stopifnot(!any(vapply(raw_fits, inherits, TRUE, "try-error")))

bin_diff <- max(abs(raw_fits[["BIN_CHECK"]]$coefficients[, cfs] - ref_fit$bin$coef[, cfs]))
log_msg(sprintf("FIT$bin max abs diff vs reference S3a: %.3e", bin_diff))
stopifnot(bin_diff < 1e-11)

comp_std <- pack(raw_fits[["STD8"]], cfs)
chr_fits <- lapply(raw_fits[chrs], pack, cfs = cfs)
names(chr_fits) <- chrs

# Assemble gene-level LOCO comp fit: gene g on chr_k gets fit from chr_fits[[chr_k]]
keep_chr <- ann$chromosome[match(KEEP, ann$gene_id_versioned)]
loco_coef <- comp_std$coef
loco_tab_list <- list()

for (cf in cfs) {
  tab_std <- copy(comp_std$tab[coef == cf])
  stopifnot(identical(tab_std$gene, KEEP))
  tab_out <- copy(tab_std)
  for (ck in chrs) {
    idx <- which(keep_chr == ck)
    if (!length(idx)) next
    tab_k <- chr_fits[[ck]]$tab[coef == cf]
    stopifnot(identical(tab_k$gene, KEEP))
    tab_out[idx] <- tab_k[idx]
    loco_coef[idx, cf] <- chr_fits[[ck]]$coef[idx, cf]
  }
  tab_out[, FDR := p.adjust(P, method = "BH")]
  loco_tab_list[[cf]] <- tab_out
}
comp_loco <- list(coef = loco_coef, tab = rbindlist(loco_tab_list))

# Preserve reference S3a object layout and update comp + comp_standard_8part
FIT <- ref_fit
FIT$comp <- comp_loco
FIT$comp_standard_8part <- comp_std
saveRDS(FIT, file.path(OUT, "observed_fits.rds"))

# Copy reference non-comp outputs and save comparison tables
for (fn in c("gene_bin.tsv.gz", "gene_step.tsv.gz", "gene_num.tsv.gz", "per_gene_states.tsv.gz",
             "perm_stage_codes.tsv.gz", "boot_indices.tsv.gz", "F_SS1_primary.tsv", "F_SS1_beside.tsv",
             "F_SS2_steps.tsv", "F_SS3_hero_igfbp7.tsv", "F_SS3_marks.tsv", "F_SS3_stage_supported.tsv")) {
  src <- file.path(REF_S3A, fn)
  if (file.exists(src)) file.copy(src, file.path(OUT, fn), overwrite = TRUE)
}

cmp <- merge(
  comp_loco$tab[coef == "binF3F4", .(gene, logFC_loco = logFC, SE_loco = SE, t_loco = t, P_loco = P, FDR_loco = FDR)],
  comp_std$tab[coef == "binF3F4", .(gene, logFC_std8 = logFC, SE_std8 = SE, t_std8 = t, P_std8 = P, FDR_std8 = FDR)],
  by = "gene"
)
cmp[, chromosome := ann$chromosome[match(gene, ann$gene_id_versioned)]]
cmp[, symbol := ann$gene_name[match(gene, ann$gene_id_versioned)]]
old_comp <- as.data.table(ref_fit$comp$tab)[coef == "binF3F4", .(gene, logFC_old7 = logFC, SE_old7 = SE, P_old7 = P, FDR_old7 = FDR)]
cmp <- merge(cmp, old_comp, by = "gene", all.x = TRUE)
fwrite(cmp, file.path(OUT, "comp_loco_vs_standard.tsv.gz"), sep = "\t")

summary_dt <- data.table(
  metric = c(
    "n_staged_units", "n_keep_genes", "n_autosomal_loco_genes", "n_non_autosomal_std_genes",
    "bin_max_abs_diff_vs_ref",
    "loco_F3F4_FDR_lt_0.05", "std8_F3F4_FDR_lt_0.05", "old7_F3F4_FDR_lt_0.05",
    "spearman_logFC_loco_vs_std8", "spearman_logFC_loco_vs_old7", "max_abs_diff_logFC_loco_vs_std8"
  ),
  value = c(
    nrow(u), length(KEEP), sum(keep_chr %in% chrs), sum(!(keep_chr %in% chrs)),
    bin_diff,
    sum(cmp$FDR_loco < 0.05), sum(cmp$FDR_std8 < 0.05), sum(cmp$FDR_old7 < 0.05, na.rm = TRUE),
    cor(cmp$logFC_loco, cmp$logFC_std8, method = "spearman"),
    cor(cmp$logFC_loco, cmp$logFC_old7, method = "spearman"),
    max(abs(cmp$logFC_loco - cmp$logFC_std8))
  )
)
fwrite(summary_dt, file.path(OUT, "loco_s3a_summary.tsv"), sep = "\t")
print(summary_dt)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("S3a LOCO complete.")
