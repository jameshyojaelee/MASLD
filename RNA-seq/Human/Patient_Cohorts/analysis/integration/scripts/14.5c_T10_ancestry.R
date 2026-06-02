#!/usr/bin/env Rscript
# 14.5c_T10_ancestry.R
# ---------------------------------------------------------------------------
# A7 T10: Ancestry-proxy interaction term
#
# Annotate cohorts with ancestry proxy (Chen = East Asian, others = European)
# and refit the sex-interaction dream model with an additional fixed effect
# `ancestry` (we cannot use cohort:ancestry interaction because cohort is
# the random-effect grouping; instead we add `ancestry` as a fixed effect,
# which captures ancestry-level mean shifts).
#
# Compare canonical vs ancestry-adjusted interaction LFC.
#
# Outputs:
#   audit_sensitivity/sex_ancestry/ancestry_interaction_dream.csv
#   audit_sensitivity/sex_ancestry/concordance_to_canonical.csv
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml)
})

ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_ancestry")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("CPUs:", ncpus, "\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load (same as Script 26) ---
cat("Loading data...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]

qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) dge_mega <- dge_mega[, !sex_fail]

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# --- Ancestry proxy: Chen = EAS; Suppli/Hoang/Govaere/Bril = EUR ---
ancestry_map <- function(ds) {
  ds_chr <- as.character(ds)
  ifelse(grepl("Chen|GSE164760", ds_chr, ignore.case = TRUE), "EAS", "EUR")
}
ancestry_vec <- ancestry_map(dge_mega$samples$dataset)
cat("Ancestry assignment:\n"); print(table(dge_mega$samples$dataset, ancestry_vec))

info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  ancestry     = factor(ancestry_vec, levels = c("EUR", "EAS")),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) { dge_mega <- dge_mega[, !na_sex]; info <- info[!na_sex, , drop = FALSE] }

cat("group x sex x ancestry table:\n"); print(table(info$group_binary, info$inferred_sex, info$ancestry))

# --- filterByExpr + TMM ---
keep_g <- filterByExpr(dge_mega, group = info$group_binary)
dge_filt <- dge_mega[keep_g, , keep.lib.sizes = FALSE]
dge_filt <- calcNormFactors(dge_filt, method = "TMM")
cat("Genes after filter:", nrow(dge_filt), "\n")

# --- Model: ~ group_binary * inferred_sex + ancestry + (1|dataset) ---
# Note: cohort:ancestry would alias with dataset random effect (ancestry is
# a deterministic function of cohort in this proxy). We use ancestry as
# a fixed effect to test whether adjusting for it changes sex-interaction LFC.
form <- ~ group_binary * inferred_sex + ancestry + (1 | dataset)
cat("Formula:", deparse(form), "\n")

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_filt, form, info, BPPARAM = param))
cat("Running dream...\n")
fit <- suppressWarnings(dream(v, form, info, BPPARAM = param))

all_coefs <- colnames(fit$coefficients)
cat("Available coefficients:", paste(all_coefs, collapse = ", "), "\n")
int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary", all_coefs, value = TRUE)[1]
cat("Using interaction coefficient:", int_coef, "\n")

res <- topTable(fit, coef = int_coef, number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res); setnames(res_dt, "adj.P.Val", "padj")
fwrite(res_dt, file.path(OUT, "ancestry_interaction_dream.csv"))
cat("Saved ancestry_interaction_dream.csv (n=", nrow(res_dt), ")\n")

# --- Compare to canonical sex_interaction_dream.csv ---
canon <- fread(file.path(RDIR, "sex_interaction_dream.csv"))
mrg <- merge(
  res_dt[, .(gene, logFC_anc = logFC, padj_anc = padj)],
  canon[,  .(gene, logFC_canon = logFC, padj_canon = padj)],
  by = "gene"
)
cat("Genes shared with canonical:", nrow(mrg), "\n")

rho <- cor(mrg$logFC_anc, mrg$logFC_canon, method = "spearman", use = "pairwise.complete.obs")
rho_p <- cor(mrg$logFC_anc, mrg$logFC_canon, method = "pearson", use = "pairwise.complete.obs")

# Top-100 sex-dimorphic by canonical |LFC| among padj<0.05
top100_canon <- canon[!is.na(padj) & padj < 0.05][order(-abs(logFC))][1:min(100, .N), gene]
top100_anc   <- res_dt[!is.na(padj) & padj < 0.05][order(-abs(logFC))][1:min(100, .N), gene]
jacc <- length(intersect(top100_canon, top100_anc)) /
        max(length(union(top100_canon, top100_anc)), 1L)

n_dim_canon <- sum(canon$padj < 0.05, na.rm = TRUE)
n_dim_anc   <- sum(res_dt$padj < 0.05, na.rm = TRUE)

concord <- data.table(
  metric = c("spearman_rho_logFC", "pearson_r_logFC",
             "top100_jaccard", "n_sex_dimorphic_canonical_padj05",
             "n_sex_dimorphic_ancestry_padj05", "n_genes_shared"),
  value  = c(rho, rho_p, jacc, n_dim_canon, n_dim_anc, nrow(mrg))
)
fwrite(concord, file.path(OUT, "concordance_to_canonical.csv"))
cat("\nConcordance to canonical:\n"); print(concord)
cat("Done:", as.character(Sys.time()), "\n")
