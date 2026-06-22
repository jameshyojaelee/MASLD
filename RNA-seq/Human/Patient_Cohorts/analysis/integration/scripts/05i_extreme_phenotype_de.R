#!/usr/bin/env Rscript
# ===========================================================================
# 05i — EXTREME-PHENOTYPE sensitivity DE: definite-disease vs strict-control.
#
# A confirmatory / sensitivity arm of the canonical disease-vs-control call
# (05h limma_voom_qw__C2). Instead of pooling ALL disease vs ALL control, it
# contrasts the two histologic extremes:
#   definite-disease  := NAS >= 5  OR  fibrosis >= 3        (definite MASH / advanced fibrosis)
#   strict-control    := labeled Control & NAS <= 1 & F0/NA (pristine histology)
# This is the DE analog of the figS_pca_definitive.R "definitive disease" axis.
#
# DESIGN RATIONALE
#   - Same engine as the canonical (voomWithQualityWeights -> lmFit -> eBayes,
#     ashr mixcompdist="normal"), so it is a true method-matched sensitivity arm.
#   - Restricted to mega cohorts that carry BOTH arms (>=5 each) so the cohort
#     fixed effect is estimable and NOT confounded with the phenotype. Cohorts
#     with only one arm (e.g. GSE126848 control-only) are dropped automatically.
#   - Extreme-phenotype sampling INFLATES effect sizes -> the DEG count is NOT
#     directly comparable to the canonical 1,853. This is a confirmatory arm,
#     not a canonical replacement. It (a) tests whether the consensus signal is
#     diluted by ambiguous mild cases, and (b) yields a high-confidence
#     "definite-MASH core" gene subset.
#
# Output: results/integration/sensitivity/
#   extreme_phenotype_definite_vs_strict.csv    (DE table, raw + ashr)
#   extreme_phenotype_definite_vs_strict_meta.csv (per-sample arm/cohort labels
#                                                   consumed by the paired heatmap)
#   extreme_phenotype_sanity.csv                (overlap vs canonical 1,853)
# ===========================================================================
suppressMessages({library(edgeR); library(limma); library(ashr); library(data.table)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
SDIR <- file.path(RDIR, "sensitivity"); dir.create(SDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# ---- inputs (identical loading to 05h) ------------------------------------
dge  <- readRDS(file.path(RDIR, "merged_dge.rds"))
mm   <- readRDS(file.path(RDIR, "meta_matched.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))

i    <- match(colnames(dge), mm$sample_id)
nas  <- suppressWarnings(as.numeric(mm$nas_score[i]))
fib  <- suppressWarnings(as.numeric(mm$fibrosis_stage[i]))
sex  <- mm$inferred_sex[i]
gb   <- as.character(dge$samples$group_binary)
ds   <- as.character(dge$samples$dataset)

# ---- extreme-phenotype arms -----------------------------------------------
definite <- (!is.na(nas) & nas >= 5) | (!is.na(fib) & fib >= 3)
strict   <- gb == "Control" & (is.na(nas) | nas <= 1) & (is.na(fib) | fib <= 0)

# ---- keep mega cohorts carrying BOTH arms (>=5 each) ----------------------
arm_tab  <- data.frame(dataset = mega,
  definite = sapply(mega, function(d) sum(definite & ds == d)),
  strict   = sapply(mega, function(d) sum(strict   & ds == d)))
keep_co  <- arm_tab$dataset[arm_tab$definite >= 5 & arm_tab$strict >= 5]
sel      <- ds %in% keep_co & (definite | strict)
cat("[design] arm counts per mega cohort:\n"); print(arm_tab, row.names = FALSE)
cat(sprintf("[design] cohorts with both arms: %s\n", paste(keep_co, collapse = ", ")))

dge_x <- dge[, sel]
info  <- data.frame(
  phenotype = factor(ifelse(definite[sel], "Definite", "StrictControl"),
                     levels = c("StrictControl", "Definite")),
  dataset   = droplevels(factor(ds[sel])),
  inferred_sex = factor(sex[sel]))
rownames(info) <- colnames(dge_x)
cat(sprintf("[in] %d genes x %d samples; %d cohorts; %s\n",
            nrow(dge_x), ncol(dge_x), nlevels(info$dataset),
            paste(names(table(info$phenotype)), table(info$phenotype), collapse = " ")))

# ---- re-normalize within the extreme-phenotype subset ---------------------
# Recompute TMM factors on the subset so library scaling matches the samples
# actually being compared (the 9-cohort universe genes are retained as-is).
dge_x <- calcNormFactors(dge_x)

# ---- C2 design (phenotype = effect of interest; cohort + sex nuisance) -----
design    <- model.matrix(~ dataset + inferred_sex + phenotype, data = info)
coef_name <- "phenotypeDefinite"
stopifnot(coef_name %in% colnames(design))
stopifnot(limma::is.fullrank(design))

# ---- limma_voom_qw engine (benchmark winner) ------------------------------
v   <- voomWithQualityWeights(dge_x, design)
fit <- eBayes(lmFit(v, design))
res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)
dt  <- as.data.table(res); setnames(dt, "adj.P.Val", "padj")
dt[, SE := abs(logFC / t)]

# ---- ashr shrinkage (side-by-side with raw, mixcompdist="normal") ----------
ok  <- is.finite(dt$logFC) & is.finite(dt$SE) & dt$SE > 0
ash <- ashr::ash(dt$logFC[ok], dt$SE[ok], mixcompdist = "normal")
dt[, shrunk_logFC := NA_real_][ok, shrunk_logFC := ash$result$PosteriorMean]
dt[, lfsr := NA_real_][ok, lfsr := ash$result$lfsr]

# ---- symbol column (GENCODE v49) ------------------------------------------
gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
gm[, eb := sub("[.][0-9]+$", "", gene_id)]
dt[, symbol := gm[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]

out <- dt[, .(gene, symbol, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr)]
out[, method := "limma_voom_qw__C2__extreme_phenotype"]
fwrite(out, file.path(SDIR, "extreme_phenotype_definite_vs_strict.csv"))

# per-sample labels for the paired heatmap
meta_out <- data.table(sample_id = colnames(dge_x),
                       dataset = as.character(info$dataset),
                       inferred_sex = as.character(info$inferred_sex),
                       phenotype = as.character(info$phenotype),
                       nas_score = nas[sel], fibrosis_stage = fib[sel])
fwrite(meta_out, file.path(SDIR, "extreme_phenotype_definite_vs_strict_meta.csv"))

# ---- Tier-1 counts (both definitions) -------------------------------------
raw1 <- out[padj < 0.05 & abs(logFC) > 0.5]
ash1 <- out[lfsr < 0.05 & abs(shrunk_logFC) > 0.5]
cat(sprintf("\n[extreme] genes=%d\n  RAW  Tier-1 (padj<.05 & |logFC|>.5)        : %d  (up %d / dn %d)\n  ashr Tier-1 (lfsr<.05 & |shrunk_logFC|>.5) : %d  (up %d / dn %d)\n",
  nrow(out), nrow(raw1), sum(raw1$logFC > 0), sum(raw1$logFC < 0),
  nrow(ash1), sum(ash1$shrunk_logFC > 0), sum(ash1$shrunk_logFC < 0)))

# ---- sanity: concordance vs canonical 1,853 -------------------------------
can <- fread(file.path(RDIR, "canonical_deg_results.csv"))
m   <- merge(out[, .(gene, x_lfc = logFC, x_slfc = shrunk_logFC, x_lfsr = lfsr)],
             can[, .(gene, c_lfc = logFC, c_slfc = shrunk_logFC, c_lfsr = lfsr)], by = "gene")
rho  <- cor(m$x_lfc, m$c_lfc, method = "spearman", use = "complete.obs")
dirA <- mean(sign(m$x_lfc) == sign(m$c_lfc), na.rm = TRUE)
can1 <- can[lfsr < 0.05 & abs(shrunk_logFC) > 0.5, gene]
ext1 <- ash1$gene
jac  <- length(intersect(ext1, can1)) / length(union(ext1, can1))
recov<- length(intersect(ext1, can1)) / length(can1)

san <- data.table(
  metric = c("genes_tested", "n_definite", "n_strict_control", "cohorts",
             "extreme_raw_Tier1", "extreme_ashr_Tier1", "canonical_ashr_Tier1",
             "logFC_spearman_vs_canonical", "direction_agree_vs_canonical",
             "jaccard_vs_canonical", "canonical_recovered_frac"),
  value  = c(nrow(out), sum(info$phenotype == "Definite"), sum(info$phenotype == "StrictControl"),
             nlevels(info$dataset), nrow(raw1), nrow(ash1), length(can1),
             round(rho, 4), round(dirA, 4), round(jac, 4), round(recov, 4)))
fwrite(san, file.path(SDIR, "extreme_phenotype_sanity.csv"))
cat("\n=== SANITY vs canonical disease-vs-control (1,853) ===\n"); print(san, row.names = FALSE)
cat("\n[done] wrote extreme_phenotype_definite_vs_strict.csv (+ _meta, _sanity)\n")
