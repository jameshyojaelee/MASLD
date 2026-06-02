#!/usr/bin/env Rscript
# B1 step 5 - Compare kallisto vs STAR + featureCounts.
# Inputs:
#   - kallisto pooled TPM/counts: RNA-seq/results/kallisto/all_cohorts_gene_{tpm,counts}.tsv.gz
#   - canonical merged_dge.rds (STAR + featureCounts) and dream_results.csv
#   - kallisto dream output: RNA-seq/results/kallisto/dream_results_kallisto.csv
# Output: RNA-seq/results/kallisto_vs_star_comparison/REPORT.md + tables/

suppressPackageStartupMessages({
  library(data.table); library(edgeR)
})

WT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
RDIR_STAR <- file.path(BASE, "analysis/integration/results/integration")
KALL <- file.path(WT, "RNA-seq/results/kallisto")
OUT  <- file.path(WT, "RNA-seq/results/kallisto_vs_star_comparison")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

cat("[", as.character(Sys.time()), "] loading inputs\n")
# 1) Kallisto TPM (gene_id with ENSG version)
ktpm <- fread(file.path(KALL, "all_cohorts_gene_tpm.tsv.gz"))
kgid <- sub("\\..*$", "", ktpm$gene_id)
ktpm[, gene_id := NULL]
ktpm_mat <- as.matrix(ktpm); rownames(ktpm_mat) <- kgid
ktpm_mat[is.na(ktpm_mat)] <- 0
# Collapse duplicate ENSG (post-version-strip) by sum
ktpm_mat <- rowsum(ktpm_mat, group = rownames(ktpm_mat), reorder = FALSE)

# 2) STAR canonical TPM: derive from merged_dge.rds via featureCounts counts + gene length
dge <- readRDS(file.path(RDIR_STAR, "merged_dge.rds"))
star_cnts <- dge$counts
gene_len <- dge$genes$Length
if (is.null(gene_len)) {
  # try alternate column
  if ("length" %in% colnames(dge$genes)) gene_len <- dge$genes$length
}
if (is.null(gene_len)) {
  # Compute from GTF — featureCounts gene length = sum of exon widths per gene.
  GTF_PATH <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
  cat("[", as.character(Sys.time()), "] Length column missing; computing from GTF...\n")
  suppressPackageStartupMessages(library(rtracklayer))
  gtf <- rtracklayer::import(GTF_PATH)
  exons <- gtf[gtf$type == "exon"]
  # Reduce per-gene exon ranges + sum widths
  gene_widths <- tapply(width(exons), exons$gene_id, sum)
  # Match to rownames(dge$counts)
  # First try with version, then without
  star_gids <- rownames(star_cnts)
  gw_match <- gene_widths[star_gids]
  if (sum(is.na(gw_match)) > length(gw_match) * 0.1) {
    star_gids_v <- sub("\\..*$", "", star_gids)
    gw_match2 <- gene_widths[star_gids_v]
    if (sum(is.na(gw_match2)) < sum(is.na(gw_match))) gw_match <- gw_match2
  }
  gene_len <- as.numeric(gw_match)
  names(gene_len) <- star_gids
  cat("  matched", sum(!is.na(gene_len)), "/", length(gene_len), "genes to GTF lengths\n")
  if (sum(!is.na(gene_len)) < 1000) stop("Too few gene lengths matched from GTF")
  # For genes still missing, fill with median (won't materially affect Spearman correlation)
  gene_len[is.na(gene_len)] <- median(gene_len, na.rm = TRUE)
}
rpk <- sweep(star_cnts, 1, gene_len / 1e3, "/")
star_tpm <- sweep(rpk, 2, colSums(rpk, na.rm = TRUE) / 1e6, "/")
rownames(star_tpm) <- sub("\\..*$", "", rownames(star_tpm))

# 3) Match samples + genes between kallisto and STAR
common_samples <- intersect(colnames(ktpm_mat), colnames(star_tpm))
common_genes   <- intersect(rownames(ktpm_mat), rownames(star_tpm))
cat("  common samples:", length(common_samples),
    "| common genes:", length(common_genes), "\n")
ktpm_m <- ktpm_mat[common_genes, common_samples]
stpm_m <- star_tpm[common_genes, common_samples]

# 4) Per-gene Spearman of log1p(TPM) across samples (only genes with non-zero variance)
cat("[", as.character(Sys.time()), "] per-gene Spearman of log1p(TPM)\n")
lkt <- log1p(ktpm_m); lst <- log1p(stpm_m)
gene_var <- apply(lkt, 1, var) + apply(lst, 1, var)
keep_g <- which(gene_var > 0 & complete.cases(lkt) & complete.cases(lst))
spearman_per_gene <- sapply(keep_g, function(i) suppressWarnings(
  cor(lkt[i, ], lst[i, ], method = "spearman")))
gene_spearman_dt <- data.table(gene_id = rownames(lkt)[keep_g],
                               spearman_log_tpm = spearman_per_gene)
fwrite(gene_spearman_dt, file.path(OUT, "tables/per_gene_spearman_log_tpm.csv"))

# 5) Per-sample Spearman of log1p(TPM) across genes
cat("[", as.character(Sys.time()), "] per-sample Spearman of log1p(TPM)\n")
sample_spearman <- sapply(seq_len(ncol(lkt)), function(j)
  suppressWarnings(cor(lkt[, j], lst[, j], method = "spearman")))
sample_spearman_dt <- data.table(sample_id = colnames(lkt),
                                 spearman_log_tpm_across_genes = sample_spearman)
fwrite(sample_spearman_dt, file.path(OUT, "tables/per_sample_spearman_log_tpm.csv"))

# 6) Dream LFC comparison
cat("[", as.character(Sys.time()), "] dream LFC comparison\n")
star_dream_path <- file.path(RDIR_STAR, "dream_results.csv")
kall_dream_path <- file.path(KALL, "dream_results_kallisto.csv")
have_dream <- file.exists(star_dream_path) && file.exists(kall_dream_path)
lfc_pearson <- NA_real_; lfc_spearman <- NA_real_
deg_jaccard <- NA_real_; n_deg_star <- NA_integer_; n_deg_kall <- NA_integer_
n_deg_intersect <- NA_integer_

if (have_dream) {
  star <- fread(star_dream_path)
  kall <- fread(kall_dream_path)
  if (!"gene" %in% names(star)) star$gene <- rownames(star)
  star[, gene := sub("\\..*$", "", gene)]
  kall[, gene := sub("\\..*$", "", gene)]
  common_dream <- intersect(star$gene, kall$gene)
  setkey(star, gene); setkey(kall, gene)
  m <- merge(star[common_dream, .(gene, logFC, padj)],
             kall[common_dream, .(gene, logFC, padj)],
             by = "gene", suffixes = c("_star", "_kall"))
  lfc_pearson  <- cor(m$logFC_star, m$logFC_kall,  use = "complete.obs", method = "pearson")
  lfc_spearman <- cor(m$logFC_star, m$logFC_kall,  use = "complete.obs", method = "spearman")
  star_deg <- m$gene[!is.na(m$padj_star) & m$padj_star < 0.05 & abs(m$logFC_star) >= 0.5]
  kall_deg <- m$gene[!is.na(m$padj_kall) & m$padj_kall < 0.05 & abs(m$logFC_kall) >= 0.5]
  n_deg_star <- length(star_deg); n_deg_kall <- length(kall_deg)
  n_deg_intersect <- length(intersect(star_deg, kall_deg))
  n_deg_union     <- length(union(star_deg, kall_deg))
  deg_jaccard <- if (n_deg_union > 0) n_deg_intersect / n_deg_union else NA_real_
  fwrite(m, file.path(OUT, "tables/dream_lfc_padj_merged.csv"))
} else {
  cat("[WARN] dream results not yet present; LFC comparison skipped\n")
}

# 7) Multimapper recovery: genes detected (TPM > 0) only in kallisto, only in
#    STAR, or in both. Kallisto allocates ambiguous multimappers via EM; STAR
#    + featureCounts (default -M off) discards. Genes with TPM > 0 in kallisto
#    + TPM == 0 in STAR are candidate recovered ambig-multimappers.
cat("[", as.character(Sys.time()), "] multimapper recovery\n")
kall_any <- rowSums(ktpm_m > 0) > 0
star_any <- rowSums(stpm_m > 0) > 0
multi <- data.table(
  metric = c("genes detected in both (TPM > 0 any sample)",
             "kallisto-only (recovered ambig-multimapper candidates)",
             "STAR-only",
             "total common genes considered"),
  count = c(sum(kall_any &  star_any),
            sum(kall_any & !star_any),
            sum(!kall_any &  star_any),
            length(kall_any))
)
fwrite(multi, file.path(OUT, "tables/multimapper_recovery.csv"))

# 8) Write REPORT.md
cat("[", as.character(Sys.time()), "] writing REPORT.md\n")
report <- c(
  "# Kallisto vs STAR + featureCounts comparison (B1)",
  paste0("Generated: ", Sys.time()),
  "",
  "Worktree: rna-validation-2026-05-20",
  "",
  "## 1. TPM concordance",
  sprintf("- Common samples: %d", length(common_samples)),
  sprintf("- Common genes  : %d (ENSG-version-stripped)", length(common_genes)),
  sprintf("- Per-gene Spearman of log1p(TPM) [median]: %.3f",
          median(gene_spearman_dt$spearman_log_tpm, na.rm = TRUE)),
  sprintf("- Per-gene Spearman of log1p(TPM) [Q25-Q75]: %.3f - %.3f",
          quantile(gene_spearman_dt$spearman_log_tpm, 0.25, na.rm = TRUE),
          quantile(gene_spearman_dt$spearman_log_tpm, 0.75, na.rm = TRUE)),
  sprintf("- Per-sample Spearman of log1p(TPM) [median]: %.3f",
          median(sample_spearman_dt$spearman_log_tpm_across_genes, na.rm = TRUE)),
  "",
  "## 2. Dream LFC concordance (Disease vs Control mega-analysis)",
  sprintf("- Pearson(LFC_star, LFC_kallisto)  : %s",
          if (!is.na(lfc_pearson))  sprintf("%.3f", lfc_pearson)  else "PENDING"),
  sprintf("- Spearman(LFC_star, LFC_kallisto) : %s",
          if (!is.na(lfc_spearman)) sprintf("%.3f", lfc_spearman) else "PENDING"),
  sprintf("- DEGs STAR (padj<0.05, |LFC|>=0.5): %s", if (!is.na(n_deg_star)) as.character(n_deg_star) else "PENDING"),
  sprintf("- DEGs kallisto                    : %s", if (!is.na(n_deg_kall)) as.character(n_deg_kall) else "PENDING"),
  sprintf("- Intersection                     : %s", if (!is.na(n_deg_intersect)) as.character(n_deg_intersect) else "PENDING"),
  sprintf("- Jaccard                          : %s", if (!is.na(deg_jaccard))  sprintf("%.3f", deg_jaccard)  else "PENDING"),
  "",
  "## 3. Multimapper recovery (TPM > 0 in any sample)",
  paste(capture.output(print(multi)), collapse = "\n"),
  "",
  "## 4. Recommendation",
  "",
  "Provisional rule:",
  "- If LFC Pearson >= 0.95 AND DEG Jaccard >= 0.80 -> kallisto re-quant CONCORDANT with STAR canonical;",
  "  recommend KEEPING STAR as canonical and citing kallisto as sensitivity arm in the Pachter-defense supplement.",
  "- If LFC Pearson 0.85-0.95 OR Jaccard 0.60-0.80 -> partial concordance; report both in main text;",
  "  trust kallisto recovery for genes with high multimapping (paralogs, immunoglobulins, MHC).",
  "- If LFC Pearson < 0.85 OR Jaccard < 0.60 -> material divergence; investigate gene-length normalization,",
  "  TMM, and ambig-multimapper structure before any swap.",
  "",
  "Final recommendation will be appended after dream_results_kallisto.csv lands.",
  ""
)
writeLines(report, file.path(OUT, "REPORT.md"))
cat("\nWrote", file.path(OUT, "REPORT.md"), "\n")
