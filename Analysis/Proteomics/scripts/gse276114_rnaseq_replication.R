#!/usr/bin/env Rscript
# ============================================================================
# gse276114_rnaseq_replication.R
#
# Compute independent RNA-seq replication of the 10-cohort dream mega-analysis
# using the Zeybel et al. 2025 GSE276114 bulk liver RNA-seq cohort (177 liver
# biopsies/explants; NOT the 10 cohorts). GSE276114 is advertised in the
# Zeybel paper for the NAFL-vs-NASH distinction; we replicate the dream
# NAFL-vs-NASH logFC vector at shared symbols.
#
# T0.6 context: GSE276114 was previously mis-annotated as SomaScan proteomics.
# It is actually bulk RNA-seq (see data/GSE276114/README.md). This script is
# the RNA-seq-replication analog of what used to be labelled "proteomics
# validation" — it does not re-enter the dream integration; it is an
# independent replication cohort at the mRNA layer.
#
# Primary output: Spearman rho + 95% CI between GSE276114 NAFL-vs-NASH
# DEG logFCs and the canonical dream NAFL-vs-NASH logFC vector, at shared
# HGNC gene symbols.
#
# Written: 2026-04-22
# Output:  RNA-seq/results/gse276114_rnaseq_replication/replication_summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

OUT_DIR <- file.path(BASE, "RNA-seq/results/gse276114_rnaseq_replication")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== GSE276114 RNA-seq replication ===\n")
cat("Time:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Load GSE276114 count matrix
# ---------------------------------------------------------------------------
count_f <- file.path(BASE, "data/GSE276114/GSE276114_raw.count.txt.gz")
if (!file.exists(count_f)) {
  stop("GSE276114 count matrix not found at ", count_f)
}
counts <- fread(count_f)
gene_col <- names(counts)[1]
cat("Count matrix: ", nrow(counts), " genes x ", ncol(counts) - 1, " samples\n",
    sep = "")

# Symbol column
genes <- counts[[gene_col]]
mat   <- as.matrix(counts[, -1, with = FALSE])
rownames(mat) <- genes
# Drop zero-sum genes
mat <- mat[rowSums(mat) > 0, ]
cat("After zero filter: ", nrow(mat), " genes\n", sep = "")

# ---------------------------------------------------------------------------
# Sample metadata — use the REAL disease labels deposited on disk
# (gse276114_disease_metadata.csv), NOT a k-means-on-PC1 surrogate.
#
# The metadata keys each count-matrix column via `matrix_column_name`
# (e.g. "Liver sample 1") and carries the real fibrosis-stage label in
# `disease_group` (F0-2 / F3 / F4). The dream replication target is the
# NAFL-vs-NASH disease-severity signature, so we map the real fibrosis-stage
# groups onto the mild-vs-advanced severity axis that the metadata directly
# supports: F0-2 = mild (the "NAFL" arm of the contrast), F3 + F4 = advanced
# (the "NASH" arm). Samples that lack a real label are DROPPED, not guessed.
# ---------------------------------------------------------------------------
meta_f <- file.path(BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv")
if (!file.exists(meta_f)) {
  stop("GSE276114 disease metadata not found at ", meta_f)
}
meta <- fread(meta_f)

# Map real fibrosis-stage groups onto the mild-vs-advanced (NAFL-vs-NASH)
# severity contrast. Only F0-2 / F3 / F4 are recognized; anything else → NA.
meta[, group := fcase(
  disease_group == "F0-2",            "NAFL",
  disease_group %in% c("F3", "F4"),   "NASH",
  default = NA_character_
)]

# Key the count-matrix columns to the metadata via matrix_column_name.
col2group <- meta$group
names(col2group) <- meta$matrix_column_name

group <- unname(col2group[colnames(mat)])

# Drop samples lacking a real label rather than guessing.
keep <- !is.na(group)
n_dropped <- sum(!keep)
if (n_dropped > 0) {
  cat("Dropping ", n_dropped,
      " sample(s) without a real disease label (no metadata match / unmapped stage)\n",
      sep = "")
}
mat   <- mat[, keep, drop = FALSE]
group <- group[keep]

dge <- DGEList(counts = mat)
dge <- calcNormFactors(dge)

cat("Real-label groups: NAFL(F0-2)=", sum(group == "NAFL"),
    " NASH(F3/F4)=", sum(group == "NASH"), "\n", sep = "")

# ---------------------------------------------------------------------------
# limma-voom NAFL vs NASH
# ---------------------------------------------------------------------------
group <- factor(group, levels = c("NAFL", "NASH"))
design <- model.matrix(~ group)
v   <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

tt <- topTable(fit, coef = "groupNASH", number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
setDT(tt)
setnames(tt, "adj.P.Val", "padj", skip_absent = TRUE)
setnames(tt, "P.Value",   "pvalue", skip_absent = TRUE)

fwrite(tt, file.path(OUT_DIR, "gse276114_nafl_vs_nash_limma.csv"))
cat("Saved GSE276114 NAFL-vs-NASH DE: ",
    file.path(OUT_DIR, "gse276114_nafl_vs_nash_limma.csv"), "\n", sep = "")

# ---------------------------------------------------------------------------
# Load canonical dream NAFL-vs-NASH
# ---------------------------------------------------------------------------
dream_f <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration",
  "results/disease_signatures/nafl_vs_nash_dream.csv")
if (!file.exists(dream_f)) {
  stop("Dream NAFL-vs-NASH file not found: ", dream_f)
}
dream <- fread(dream_f)

# The dream file is keyed by Ensembl IDs with optional 'symbol' column. Map.
if (!"symbol" %in% names(dream)) {
  # Try to map via the MASLD Gene Catalog
  atlas_f <- file.path(BASE,
    "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
  if (file.exists(atlas_f)) {
    atlas_lite <- fread(atlas_f, select = c("ensembl_id", "human_symbol"))
    atlas_lite[, ensembl_clean := sub("\\..*", "", ensembl_id)]
    gene_key <- intersect(c("gene", "ensembl_id"), names(dream))[1]
    dream[, ensembl_clean := sub("\\..*", "", get(gene_key))]
    dream <- merge(dream, atlas_lite[, .(ensembl_clean, symbol = human_symbol)],
                   by = "ensembl_clean", all.x = TRUE)
  }
}

# Pick logFC column (may be logFC or similar)
lfc_col <- intersect(c("logFC", "log2FoldChange"), names(dream))[1]
padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), names(dream))[1]

dream_slim <- dream[!is.na(symbol) & symbol != "",
                    .(symbol,
                      bulk_logFC = get(lfc_col),
                      bulk_padj  = get(padj_col))]
dream_slim <- dream_slim[!duplicated(symbol)]

# ---------------------------------------------------------------------------
# Merge and compute replication rho
# ---------------------------------------------------------------------------
rep_tab <- merge(tt[, .(gene, gse_logFC = logFC, gse_padj = padj)],
                 dream_slim,
                 by.x = "gene", by.y = "symbol", all = FALSE)
cat("Shared genes (dream x GSE276114): ", nrow(rep_tab), "\n", sep = "")

fwrite(rep_tab, file.path(OUT_DIR, "gse276114_vs_dream_shared.csv"))

compute_ci <- function(x, y, R = 1000, conf = 0.95) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 20) return(list(rho = NA, ci_lo = NA, ci_hi = NA, n = length(x)))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(42)
  rep_rhos <- replicate(R, {
    idx <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  ci <- quantile(rep_rhos, probs = c((1 - conf) / 2, 1 - (1 - conf) / 2),
                 na.rm = TRUE)
  list(rho = rho, ci_lo = ci[[1]], ci_hi = ci[[2]], n = length(x))
}

# Overall replication
ovr <- compute_ci(rep_tab$gse_logFC, rep_tab$bulk_logFC)

# DEG subset (bulk padj < 0.05, |bulk_logFC| > 0.3)
deg_mask <- rep_tab$bulk_padj < 0.05 & abs(rep_tab$bulk_logFC) > 0.3
deg_subset <- compute_ci(rep_tab$gse_logFC[deg_mask],
                         rep_tab$bulk_logFC[deg_mask])

# Sign concordance
dir_all <- mean(sign(rep_tab$gse_logFC) == sign(rep_tab$bulk_logFC),
                na.rm = TRUE) * 100
dir_deg <- mean(sign(rep_tab$gse_logFC[deg_mask]) ==
                sign(rep_tab$bulk_logFC[deg_mask]),
                na.rm = TRUE) * 100

summary_dt <- data.table(
  stratum     = c("All shared genes", "Dream DEGs (padj<0.05, |LFC|>0.3)"),
  n           = c(ovr$n, deg_subset$n),
  spearman_rho = c(ovr$rho, deg_subset$rho),
  ci_lower    = c(ovr$ci_lo, deg_subset$ci_lo),
  ci_upper    = c(ovr$ci_hi, deg_subset$ci_hi),
  direction_concordance_pct = c(dir_all, dir_deg)
)

fwrite(summary_dt, file.path(OUT_DIR, "replication_summary.csv"))
cat("\n=== Replication summary ===\n")
print(summary_dt)
cat("\nSaved: ", file.path(OUT_DIR, "replication_summary.csv"), "\n", sep = "")

cat("\n=== GSE276114 replication complete ===\n")
cat("Time:", as.character(Sys.time()), "\n")
