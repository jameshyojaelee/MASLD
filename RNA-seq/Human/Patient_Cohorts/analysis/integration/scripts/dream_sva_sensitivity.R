#!/usr/bin/env Rscript
# dream_sva_sensitivity.R
# ---------------------------------------------------------------------------
# Pillar D — SVA-augmented dream sensitivity.
# Estimate hidden surrogate variables with `sva::svaseq` using housekeeping
# negative controls (Eisenberg & Levanon 2013 + Risso 2014 RUVg framework).
# Re-fit dream with SVs as additional fixed effects; compare DEG list.
# Pass: dream ∩ dream+SVA Jaccard >= 0.85 → no large hidden batch.
#
# Out: results/audit_sensitivity/dream_results_sva.csv
#      results/audit_sensitivity/pillar_D_sva_concordance.csv
# ---------------------------------------------------------------------------

t0 <- proc.time()
suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml); library(sva)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({ unlockBinding(fn, ns_lme4)
          assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
          lockBinding(fn, ns_lme4) }, silent = TRUE)
  }
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge <- dge[, dge$samples$dataset %in% mega_cohorts]
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge), meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = droplevels(factor(dge$samples$dataset)),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE)
rownames(info) <- colnames(dge)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) { dge <- dge[, !na_sex]; info <- info[!na_sex, , drop = FALSE] }

# --- Negative-control housekeeping list (Eisenberg & Levanon 2013) ---
# Embedded so the script is self-contained; standard reference list of
# stable housekeeping genes used widely in RUVSeq / SVA tutorials.
# DGE rownames are GENCODE v49 Ensembl IDs with version (e.g. "ENSG00000075624.18"),
# so we map symbols → versioned Ensembl IDs via data/gencode_v49_gene_metadata.tsv.gz.
HK_SYMBOLS <- c(
  "ACTB","B2M","GAPDH","HPRT1","RPL13A","RPL27","RPL37","RPS3","RPS27A","RPS29",
  "TBP","TFRC","HMBS","SDHA","UBC","YWHAZ","PPIA","PUM1","GUSB","PSMB2",
  "PSMB4","PSMC4","PSMD2","RAB7A","REEP5","SNRPD3","VCP","VPS29","CHMP2A","C1orf43",
  "EMC7","GPI","UBA52","UBB","RPL19","RPL27A","RPS17","RPS18","EEF1A1")
gene_meta <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
hk_ensembl <- gene_meta[gene_name %in% HK_SYMBOLS, gene_id]
hk_in_dge <- intersect(hk_ensembl, rownames(dge))
cat(sprintf("Housekeeping controls: %d symbols → %d Ensembl IDs → %d present in dge\n",
            length(HK_SYMBOLS), length(hk_ensembl), length(hk_in_dge)))
if (length(hk_in_dge) < 10) stop("Need >=10 housekeeping genes; aborting (broken gene IDs?)")

# --- svaseq estimates n.sv hidden factors ---
# Cap n.sv at 10 — the cohort fixed effect already captures the dominant
# structure, so SVs target residual hidden batch (lib prep, RIN). Risso
# 2014 RUVg paper uses k=5; we use up to 10 as a conservative ceiling.
# Larger n.sv has been observed to trip a "subscript out of bounds" error
# in supervised svaseq when the SVD returns fewer components than requested.
counts <- dge$counts
mod  <- model.matrix(~ group_binary + inferred_sex + dataset, data = info)
mod0 <- model.matrix(~ inferred_sex + dataset, data = info)
ctl  <- rownames(counts) %in% hk_in_dge
n_sv_be <- tryCatch(num.sv(log2(counts + 1), mod, method = "be"),
                    error = function(e) { cat("num.sv failed:", conditionMessage(e), "\n"); NA_integer_ })
n_sv <- min(if (!is.na(n_sv_be)) n_sv_be else 5L, 10L)
n_sv <- max(n_sv, 1L)
cat(sprintf("svaseq num.sv(be) = %s, capped at %d\n", as.character(n_sv_be), n_sv))
sv_obj <- tryCatch(
  svaseq(counts, mod, mod0, n.sv = n_sv, controls = ctl),
  error = function(e) {
    # If supervised svaseq fails, fall back to a smaller n.sv and unsupervised
    fallback_n <- max(1L, n_sv %/% 2L)
    cat(sprintf("svaseq failed at n.sv=%d (%s); retrying unsupervised at n.sv=%d\n",
                n_sv, conditionMessage(e), fallback_n))
    svaseq(counts, mod, mod0, n.sv = fallback_n)
  })
n_sv <- ncol(sv_obj$sv)
cat(sprintf("Final n.sv used = %d\n", n_sv))

# --- Refit dream with SVs as fixed effects ---
sv_df <- as.data.frame(sv_obj$sv); colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
info_sva <- cbind(info, sv_df)
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_sva <- as.formula(paste("~ group_binary + inferred_sex +", sv_terms, "+ (1|dataset)"))
cat("Refit formula:", deparse(form_sva), "\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()
v_sva <- suppressWarnings(voomWithDreamWeights(dge, form_sva, info_sva, BPPARAM = param))
fit_sva <- suppressWarnings(dream(v_sva, form_sva, info_sva, BPPARAM = param))
res_sva <- topTable(fit_sva, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_sva$gene <- rownames(res_sva); res_sva_dt <- as.data.table(res_sva)
setnames(res_sva_dt, "adj.P.Val", "padj")
fwrite(res_sva_dt, file.path(OUT_DIR, "dream_results_sva.csv"))

# --- Concordance vs primary (~4,370 DEGs, kallisto canonical) ---
PRIMARY_PADJ <- 0.05; PRIMARY_LFC <- 0.5
primary <- fread(file.path(RDIR, "dream_results.csv"))
# primary uses padj column already (set in 05_dream_mega_analysis.R)
deg_primary <- primary[padj < PRIMARY_PADJ & abs(logFC) > PRIMARY_LFC, gene]
deg_sva     <- res_sva_dt[padj < PRIMARY_PADJ & abs(logFC) > PRIMARY_LFC, gene]

n_p <- length(deg_primary); n_s <- length(deg_sva)
inter <- length(intersect(deg_primary, deg_sva))
union_n <- length(union(deg_primary, deg_sva))
jaccard <- if (union_n > 0) inter / union_n else NA_real_

# Per-gene merge: which primary DEGs are concordant with dream+SVA?
m <- merge(primary[, .(gene, logFC_primary = logFC, padj_primary = padj)],
           res_sva_dt[, .(gene, logFC_sva = logFC, padj_sva = padj)], by = "gene")
m[, sig_primary := padj_primary < PRIMARY_PADJ & abs(logFC_primary) > PRIMARY_LFC]
m[, sig_sva     := padj_sva     < PRIMARY_PADJ & abs(logFC_sva)     > PRIMARY_LFC]
m[, sva_concordant := sig_primary & sig_sva &
                       sign(logFC_primary) == sign(logFC_sva)]
fwrite(m, file.path(OUT_DIR, "pillar_D_sva_concordance.csv"))

# Summary line
sum_dt <- data.table(
  n_sv = n_sv,
  n_deg_primary = n_p, n_deg_sva = n_s,
  intersect = inter, union = union_n, jaccard = round(jaccard, 4),
  rho_logFC = round(cor(m$logFC_primary, m$logFC_sva, use = "pairwise.complete.obs", method = "spearman"), 4))
fwrite(sum_dt, file.path(OUT_DIR, "pillar_D_sva_summary.csv"))
cat("Summary:\n"); print(sum_dt)
cat(sprintf("Pass (Jaccard >= 0.85): %s\n",
            ifelse(!is.na(jaccard) && jaccard >= 0.85, "PASS", "FAIL")))
cat(sprintf("Elapsed %.1f min  -- SVA sensitivity done\n", (proc.time()-t0)["elapsed"]/60))
