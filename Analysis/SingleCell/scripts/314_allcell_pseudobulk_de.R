#!/usr/bin/env Rscript
# 314_allcell_pseudobulk_de.R
#
# Aggregate all cell-type pseudobulk count matrices per donor, run limma-voom
# MASLD vs Healthy DE, and compare logFC to bulk dream results.
#
# Rationale: bulk RNA-seq is a mixture of all liver cell types. Comparing
# bulk to hepatocyte-only sc pseudobulk conflates cell-type composition
# effects with per-cell expression changes. All-cell pseudobulk is the
# apples-to-apples comparison and is expected to show substantially higher ρ.
#
# Donor-collapse fix (2026-07-12): the per-cell-type pseudobulk matrices are
# keyed on obs/sample = SEQUENCING RUN, not biological donor. For GSE244832
# (117 runs -> 18 donors), GSE202379 (67 -> 46) and GSE185477 (21 -> 3) this
# pseudoreplicates. We now collapse RUN -> DONOR (sum raw counts) via
# lib_donor_collapse.R BEFORE summing across cell types, so DE n = donors.
#
# CAVEAT (does NOT fix the scientific confound): the MASLD group is entirely
# SRR-named datasets and the Healthy anchors are partly GSM-named (megareview
# V4-0077). The `~ condition` design has no dataset/batch covariate, so the
# all-cell contrast tracks platform batch and the bulk-vs-allcell rho is
# NEGATIVE. Donor-collapse reduces pseudoreplication but leaves this confound.
# Do not treat allcell_pseudobulk_de.csv as a validated disease axis.
#
# Outputs (in Analysis/SingleCell/results_gpu_v2/pseudobulk_de/donor_collapsed/):
#   allcell_pseudobulk_de.csv        — limma-voom DE (all-cell aggregated, donor-level)
#   allcell_vs_bulk_comparison.csv   — per-gene: bulk_lfc, sc_allcell_lfc, sc_hep_lfc
#   allcell_pseudobulk_summary.txt   — ρ comparison stats
#
# Environment: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
PSEUDO_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
SC_DE_DIR   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
PROP_PATH   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures",
                         "celltype_proportions_per_sample.csv")
DREAM_PATH  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")

# Write to a NEW location so the pre-fix run-level output is preserved.
out_dir <- file.path(SC_DE_DIR, "donor_collapsed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_lines <- character()
add_log <- function(...) {
  msg <- paste0(...)
  message(msg)
  log_lines <<- c(log_lines, msg)
}

# ---------------------------------------------------------------------------
# 1. Discover cell-type pseudobulk matrices (deduplicate by canonical name)
# ---------------------------------------------------------------------------
pseudo_files <- list.files(PSEUDO_DIR, pattern = "_pseudobulk\\.csv$", full.names = TRUE)
cell_types   <- sub("_pseudobulk\\.csv$", "", basename(pseudo_files))

# Canonical names: replace spaces/+ with underscore, lowercase
canonical    <- tolower(gsub("[^A-Za-z0-9]", "_", cell_types))
# Keep only first occurrence of each canonical name
keep_idx     <- !duplicated(canonical)
pseudo_files <- pseudo_files[keep_idx]
cell_types   <- cell_types[keep_idx]
add_log(sprintf("Found %d unique cell-type matrices: %s",
                length(cell_types), paste(cell_types, collapse = ", ")))

# ---------------------------------------------------------------------------
# 2. Load matrices, aggregate using union of donors
# ---------------------------------------------------------------------------
add_log("Loading pseudobulk matrices (collapsing run -> donor per cell type)...")
mat_list <- lapply(pseudo_files, function(f) {
  dt <- fread(f)
  genes <- dt[[1]]
  mat   <- as.matrix(dt[, -1])
  rownames(mat) <- genes
  # Collapse SRR run columns to biological donor (sum raw counts); columns from
  # 1-run-per-donor datasets pass through unchanged. Done BEFORE cross-cell-type
  # summing so a donor's counts unify regardless of which cell type contributed.
  collapse_counts_to_donor(mat, BASE)
})
names(mat_list) <- cell_types

# Union of all donors; missing cell types for a donor contribute zero
all_donors <- unique(unlist(lapply(mat_list, colnames)))
add_log(sprintf("  Union donors across all cell types: %d", length(all_donors)))

# Gene universe: genes present in all matrices
all_genes <- Reduce(intersect, lapply(mat_list, rownames))
add_log(sprintf("  Genes in all matrices: %d", length(all_genes)))

# Build per-donor aggregated matrix (sum across cell types; missing = 0)
agg_mat <- matrix(0, nrow = length(all_genes), ncol = length(all_donors),
                  dimnames = list(all_genes, all_donors))
for (m in mat_list) {
  shared_donors <- intersect(colnames(m), all_donors)
  agg_mat[, shared_donors] <- agg_mat[, shared_donors] + m[all_genes, shared_donors]
}
add_log(sprintf("  Aggregated all-cell pseudobulk: %d genes x %d donors",
                nrow(agg_mat), ncol(agg_mat)))
common_donors <- all_donors

# ---------------------------------------------------------------------------
# 3. Condition metadata
# ---------------------------------------------------------------------------
add_log("Loading condition metadata (collapsing run -> donor)...")
props <- fread(PROP_PATH, select = c("sample", "condition"))
props <- unique(props)
# Collapse run-level condition metadata to donor (majority-vote condition), so
# `sample` here shares the donor namespace with the collapsed count columns.
props <- collapse_sample_meta_to_donor(props, BASE)
props <- props[condition %in% c("MASLD", "Healthy")]

# Map common_donors to conditions
sample_meta <- props[sample %in% common_donors]
# Drop donors without condition
keep_donors <- sample_meta$sample
agg_mat     <- agg_mat[, keep_donors, drop = FALSE]
condition   <- sample_meta[match(keep_donors, sample), condition]

n_masld   <- sum(condition == "MASLD")
n_healthy <- sum(condition == "Healthy")
add_log(sprintf("  MASLD: %d donors, Healthy: %d donors", n_masld, n_healthy))

if (n_masld < 3 || n_healthy < 3) {
  stop("Too few donors per condition for DE — check metadata mapping")
}

# ---------------------------------------------------------------------------
# 4. limma-voom DE
# ---------------------------------------------------------------------------
add_log("Running limma-voom MASLD vs Healthy...")
dge <- DGEList(counts = agg_mat, group = condition)
dge <- dge[filterByExpr(dge, group = condition), ]
add_log(sprintf("  After filterByExpr: %d genes", nrow(dge)))
dge <- calcNormFactors(dge)

design <- model.matrix(~ condition, data = data.frame(condition = factor(condition, levels = c("Healthy", "MASLD"))))
v      <- voom(dge, design, plot = FALSE)
fit    <- lmFit(v, design)
fit    <- eBayes(fit, trend = TRUE)
res    <- topTable(fit, coef = "conditionMASLD", number = Inf, sort.by = "none")
res    <- as.data.table(res, keep.rownames = "gene")
setnames(res, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
              c("lfc",   "ave",     "t_stat", "pvalue", "padj",    "B"))

out_de <- file.path(out_dir, "allcell_pseudobulk_de.csv")
fwrite(res, out_de)
add_log(sprintf("  Wrote %s (%d genes)", out_de, nrow(res)))

n_sig <- sum(res$padj < 0.05 & abs(res$lfc) > 0.5, na.rm = TRUE)
add_log(sprintf("  Significant DEGs (padj<0.05, |LFC|>0.5): %d", n_sig))

# ---------------------------------------------------------------------------
# 5. Compare ρ: all-cell pseudobulk vs bulk vs hepatocyte-only
# ---------------------------------------------------------------------------
add_log("Computing bulk vs sc logFC correlations...")
bulk <- fread(DREAM_PATH, select = c("gene", "logFC", "padj", "symbol"))
bulk[, ensg_base := sub("\\.\\d+$", "", gene)]
setnames(bulk, c("logFC", "padj"), c("bulk_lfc", "bulk_padj"))

res[, ensg_base := sub("\\.\\d+$", "", gene)]

# Load hepatocyte-only for comparison
sc_hep <- fread(file.path(SC_DE_DIR, "Hepatocytes_de.csv"),
                select = c("gene", "logFC", "padj"))
sc_hep[, ensg_base := sub("\\.\\d+$", "", gene)]
setnames(sc_hep, c("logFC", "padj"), c("hep_lfc", "hep_padj"))

# Merge all three on base ENSG
comparison <- Reduce(
  function(a, b) merge(a, b, by = "ensg_base"),
  list(
    bulk[, .(ensg_base, symbol, bulk_lfc, bulk_padj)],
    res[,  .(ensg_base, allcell_lfc = lfc, allcell_padj = padj)],
    sc_hep[, .(ensg_base, hep_lfc, hep_padj)]
  )
)
comparison <- comparison[!is.na(bulk_lfc) & !is.na(allcell_lfc) & !is.na(hep_lfc)]
add_log(sprintf("  Three-way shared genes: %d", nrow(comparison)))

rho_allcell <- cor(comparison$bulk_lfc, comparison$allcell_lfc,
                   method = "spearman", use = "complete.obs")
rho_hep     <- cor(comparison$bulk_lfc, comparison$hep_lfc,
                   method = "spearman", use = "complete.obs")

# Direction concordance among bulk DEGs
bulk_degs <- comparison[bulk_padj < 0.05 & abs(bulk_lfc) > 0.5]
dir_allcell <- mean(sign(bulk_degs$bulk_lfc) == sign(bulk_degs$allcell_lfc), na.rm = TRUE)
dir_hep     <- mean(sign(bulk_degs$bulk_lfc) == sign(bulk_degs$hep_lfc), na.rm = TRUE)

out_comp <- file.path(out_dir, "allcell_vs_bulk_comparison.csv")
fwrite(comparison, out_comp)
add_log(sprintf("  Wrote comparison table: %s", out_comp))

summary_lines <- c(
  "[donor-collapsed run; ~condition design has NO dataset covariate]",
  "[CAVEAT: MASLD=SRR datasets vs Healthy anchors partly GSM-named (megareview",
  " V4-0077) => allcell contrast is batch-confounded; negative rho is expected",
  " and is NOT rescued by donor-collapse. Hepatocyte-only ref below is the",
  " pre-fix run-level Hepatocytes_de.csv (pending its own donor-collapse rerun).]",
  "",
  sprintf("All-cell pseudobulk vs bulk:   ρ = %.3f  (n=%d)", rho_allcell, nrow(comparison)),
  sprintf("Hepatocyte-only  vs bulk:      ρ = %.3f  (n=%d)", rho_hep,     nrow(comparison)),
  sprintf("Direction concordance (bulk DEGs, n=%d):", nrow(bulk_degs)),
  sprintf("  All-cell pseudobulk: %.1f%%", dir_allcell * 100),
  sprintf("  Hepatocyte-only:     %.1f%%", dir_hep * 100)
)
add_log("")
for (l in summary_lines) add_log(l)

out_summary <- file.path(out_dir, "allcell_pseudobulk_summary.txt")
writeLines(log_lines, out_summary)
add_log(sprintf("Wrote summary: %s", out_summary))
message("Done.")
