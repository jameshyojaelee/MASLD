#!/usr/bin/env Rscript
# =============================================================================
# 53_ncrna_landscape.R
# Module 1: ncRNA Landscape Characterization
#
# Characterizes non-coding RNA (lncRNA, miRNA, snoRNA, snRNA, misc_RNA, scaRNA)
# in the MASLD transcriptomic atlas. Produces per-biotype DEG summaries,
# statistical comparisons vs protein-coding genes, and validates known MASLD
# lncRNAs against dream mega-analysis results.
# =============================================================================

# HISTORICAL ONLY. This mixed-ncRNA, outcome-selected workflow is not part of
# the standalone MASLD Resource production path.
if (Sys.getenv("ALLOW_HISTORICAL_NCRNA_PIPELINE") != "1") {
  stop("Retired Resource workflow. Set ALLOW_HISTORICAL_NCRNA_PIPELINE=1 only for provenance recovery.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(rtracklayer)
  library(GenomicRanges)
})

select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
dream_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
gtf_path   <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
out_dir    <- file.path(BASE, "RNA-seq/results/ncrna")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module 1: ncRNA Landscape Characterization ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================
cat("--- 1. Loading atlas and dream results ---\n")

atlas <- fread(atlas_path)
dream <- fread(dream_path)
# W3/W4 fix: trim whitespace from gene identifiers to prevent merge failures
atlas[, human_symbol := trimws(human_symbol)]
dream[, gene := trimws(gene)]
cat("  Atlas:", nrow(atlas), "genes x", ncol(atlas), "columns\n")
cat("  Dream:", nrow(dream), "genes\n")

# Define ncRNA biotype classes
ncrna_biotypes <- c("lncRNA", "miRNA", "snoRNA", "snRNA", "misc_RNA", "scaRNA")

# Biotype counts
bt_counts <- atlas[, .N, by = gene_biotype][order(-N)]
cat("  Biotype distribution:\n")
for (i in seq_len(min(15, nrow(bt_counts)))) {
  cat(sprintf("    %-35s %6d\n", bt_counts$gene_biotype[i], bt_counts$N[i]))
}

# =============================================================================
# 2. Per-Biotype DEG Summary
# =============================================================================
cat("\n--- 2. Per-biotype DEG summary ---\n")

# Merge dream stats into atlas for genes with dream results
# Atlas already has bulk_logFC, bulk_padj, bulk_tstat
# Define DEG status
# Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.3 (Script 05b)
atlas[, is_deg := !is.na(bulk_padj) & bulk_padj < 0.1]
atlas[, is_deg_strict := !is.na(bulk_treat_fdr) & bulk_treat_fdr < 0.05]
atlas[, deg_direction := fifelse(is_deg & bulk_logFC > 0, "up",
                          fifelse(is_deg & bulk_logFC < 0, "down", "ns"))]

# Include protein_coding for comparison
summary_biotypes <- c("protein_coding", ncrna_biotypes)

landscape_summary <- atlas[gene_biotype %in% summary_biotypes, .(
  n_total       = .N,
  n_tested      = sum(!is.na(bulk_padj)),
  n_deg         = sum(is_deg, na.rm = TRUE),
  n_deg_strict  = sum(is_deg_strict, na.rm = TRUE),
  n_up          = sum(deg_direction == "up", na.rm = TRUE),
  n_down        = sum(deg_direction == "down", na.rm = TRUE),
  frac_deg      = sum(is_deg, na.rm = TRUE) / max(sum(!is.na(bulk_padj)), 1),
  median_abs_lfc = as.double(median(abs(bulk_logFC[is_deg]), na.rm = TRUE)),
  median_aveexpr = as.double(median(bulk_tstat[!is.na(bulk_tstat)], na.rm = TRUE)),
  median_sources = as.double(median(sources_active[!is.na(sources_active)], na.rm = TRUE))
), by = gene_biotype]

# Recompute median AveExpr from dream results directly
# W4 fix: dream uses versioned Ensembl IDs (ENSG*.1), atlas has unversioned ensembl_id
dream[, gene_base := sub("\\..*", "", gene)]
dream_biotype <- merge(dream, atlas[, .(gene_base = ensembl_id, gene_biotype)],
                        by = "gene_base", all.x = FALSE)
cat(sprintf("  Dream-atlas biotype merge: %d / %d dream genes matched\n",
            nrow(dream_biotype), nrow(dream)))
aveexpr_summary <- dream_biotype[, .(median_aveexpr = median(AveExpr, na.rm = TRUE)),
                                  by = gene_biotype]
landscape_summary[, median_aveexpr := NULL]
landscape_summary <- merge(landscape_summary, aveexpr_summary, by = "gene_biotype", all.x = TRUE)

setorder(landscape_summary, -n_total)
cat("  Per-biotype summary:\n")
print(landscape_summary)

fwrite(landscape_summary, file.path(out_dir, "ncrna_landscape_summary.csv"))
cat("  Saved ncrna_landscape_summary.csv\n")

# =============================================================================
# 3. Statistical Comparisons: lncRNA vs Protein-Coding
# =============================================================================
cat("\n--- 3. lncRNA vs protein-coding comparisons ---\n")

pc_genes <- atlas[gene_biotype == "protein_coding" & !is.na(bulk_padj)]
lnc_genes <- atlas[gene_biotype == "lncRNA" & !is.na(bulk_padj)]

comparisons <- list()

# 3a. |LFC| comparison (among DEGs)
pc_degs <- pc_genes[is_deg == TRUE]
lnc_degs <- lnc_genes[is_deg == TRUE]

if (nrow(pc_degs) > 0 & nrow(lnc_degs) > 0) {
  lfc_test <- wilcox.test(abs(lnc_degs$bulk_logFC), abs(pc_degs$bulk_logFC))
  comparisons$lfc <- data.table(
    comparison = "abs_logFC_among_DEGs",
    lncrna_median = median(abs(lnc_degs$bulk_logFC)),
    pc_median = median(abs(pc_degs$bulk_logFC)),
    wilcox_p = lfc_test$p.value,
    direction = ifelse(median(abs(lnc_degs$bulk_logFC)) > median(abs(pc_degs$bulk_logFC)),
                       "lncRNA_higher", "PC_higher")
  )
  cat(sprintf("  |LFC| among DEGs: lncRNA median=%.3f, PC median=%.3f, p=%.2e\n",
              comparisons$lfc$lncrna_median, comparisons$lfc$pc_median, comparisons$lfc$wilcox_p))
}

# 3b. AveExpr comparison (all tested)
# W3/W4 fix: dream uses versioned Ensembl IDs; match via ensembl_id
if (!"gene_base" %in% names(dream)) dream[, gene_base := sub("\\..*", "", gene)]
aveexpr_lnc <- merge(lnc_genes[, .(gene_base = ensembl_id)], dream, by = "gene_base")
aveexpr_pc <- merge(pc_genes[, .(gene_base = ensembl_id)], dream, by = "gene_base")
cat(sprintf("  Dream merge: lncRNA %d/%d matched, PC %d/%d matched\n",
            nrow(aveexpr_lnc), nrow(lnc_genes), nrow(aveexpr_pc), nrow(pc_genes)))

if (nrow(aveexpr_lnc) > 0 & nrow(aveexpr_pc) > 0) {
  expr_test <- wilcox.test(aveexpr_lnc$AveExpr, aveexpr_pc$AveExpr)
  comparisons$expr <- data.table(
    comparison = "AveExpr_all_tested",
    lncrna_median = median(aveexpr_lnc$AveExpr, na.rm = TRUE),
    pc_median = median(aveexpr_pc$AveExpr, na.rm = TRUE),
    wilcox_p = expr_test$p.value,
    direction = ifelse(median(aveexpr_lnc$AveExpr, na.rm = TRUE) > median(aveexpr_pc$AveExpr, na.rm = TRUE),
                       "lncRNA_higher", "PC_higher")
  )
  cat(sprintf("  AveExpr: lncRNA median=%.3f, PC median=%.3f, p=%.2e\n",
              comparisons$expr$lncrna_median, comparisons$expr$pc_median, comparisons$expr$wilcox_p))
}

# 3c. Variance comparison
if (nrow(aveexpr_lnc) > 0 & nrow(aveexpr_pc) > 0) {
  var_lnc <- var(aveexpr_lnc$AveExpr, na.rm = TRUE)
  var_pc <- var(aveexpr_pc$AveExpr, na.rm = TRUE)
  ftest <- var.test(aveexpr_lnc$AveExpr, aveexpr_pc$AveExpr)
  comparisons$variance <- data.table(
    comparison = "AveExpr_variance",
    lncrna_median = var_lnc,
    pc_median = var_pc,
    wilcox_p = ftest$p.value,
    direction = ifelse(var_lnc > var_pc, "lncRNA_higher", "PC_higher")
  )
  cat(sprintf("  Variance: lncRNA=%.3f, PC=%.3f, F-test p=%.2e\n",
              var_lnc, var_pc, ftest$p.value))
}

# 3d. sources_active comparison (among DEGs)
if (nrow(lnc_degs) > 0 & nrow(pc_degs) > 0) {
  src_lnc <- lnc_degs[!is.na(sources_active)]$sources_active
  src_pc <- pc_degs[!is.na(sources_active)]$sources_active
  if (length(src_lnc) > 0 & length(src_pc) > 0) {
    src_test <- wilcox.test(src_lnc, src_pc)
    comparisons$sources <- data.table(
      comparison = "sources_active_among_DEGs",
      lncrna_median = median(src_lnc),
      pc_median = median(src_pc),
      wilcox_p = src_test$p.value,
      direction = ifelse(median(src_lnc) > median(src_pc), "lncRNA_higher", "PC_higher")
    )
    cat(sprintf("  sources_active among DEGs: lncRNA median=%.1f, PC median=%.1f, p=%.2e\n",
                median(src_lnc), median(src_pc), src_test$p.value))
  }
}

# 3e. DEG enrichment (Fisher's test)
deg_table <- matrix(c(
  sum(lnc_genes$is_deg), nrow(lnc_genes) - sum(lnc_genes$is_deg),
  sum(pc_genes$is_deg), nrow(pc_genes) - sum(pc_genes$is_deg)
), nrow = 2, byrow = TRUE)
fisher_deg <- fisher.test(deg_table)
comparisons$deg_enrichment <- data.table(
  comparison = "DEG_enrichment_lncRNA_vs_PC",
  lncrna_median = sum(lnc_genes$is_deg) / nrow(lnc_genes),
  pc_median = sum(pc_genes$is_deg) / nrow(pc_genes),
  wilcox_p = fisher_deg$p.value,
  direction = ifelse(fisher_deg$estimate > 1, "lncRNA_enriched", "PC_enriched")
)
cat(sprintf("  DEG enrichment: lncRNA %.1f%% vs PC %.1f%%, OR=%.2f, p=%.2e\n",
            100 * sum(lnc_genes$is_deg) / nrow(lnc_genes),
            100 * sum(pc_genes$is_deg) / nrow(pc_genes),
            fisher_deg$estimate, fisher_deg$p.value))

comparison_dt <- rbindlist(comparisons, fill = TRUE)
fwrite(comparison_dt, file.path(out_dir, "ncrna_vs_pc_comparison.csv"))
cat("  Saved ncrna_vs_pc_comparison.csv\n")

# =============================================================================
# 4. Annotate Top lncRNA DEGs with Transcript Class
# =============================================================================
cat("\n--- 4. Annotating top lncRNA DEGs ---\n")

# Load GENCODE GTF for transcript class annotation
cat("  Loading GENCODE v49 GTF...\n")
gtf <- import(gtf_path)
gtf_genes <- gtf[gtf$type == "gene"]

# Extract lncRNA gene annotations
lnc_gtf <- gtf_genes[gtf_genes$gene_type == "lncRNA"]
lnc_annot <- data.table(
  gene_name = lnc_gtf$gene_name,
  gene_id = lnc_gtf$gene_id,
  chr = as.character(seqnames(lnc_gtf)),
  start = start(lnc_gtf),
  end = end(lnc_gtf),
  strand = as.character(strand(lnc_gtf))
)

# Get transcript-level info for classification
lnc_tx <- gtf[gtf$type == "transcript" & gtf$gene_type == "lncRNA"]
lnc_tx_dt <- data.table(
  gene_name = lnc_tx$gene_name,
  transcript_type = lnc_tx$transcript_type
)

# Classify by transcript type
# GENCODE transcript types for lncRNA: lincRNA (intergenic), antisense, etc.
tx_class <- lnc_tx_dt[, .(
  transcript_types = paste(unique(transcript_type), collapse = ";")
), by = gene_name]

tx_class[, transcript_class := fcase(
  grepl("antisense", transcript_types), "antisense",
  grepl("lincRNA", transcript_types), "intergenic",
  grepl("bidirectional", transcript_types), "bidirectional",
  grepl("sense_intronic", transcript_types), "intronic",
  grepl("sense_overlapping", transcript_types), "sense_overlapping",
  default = "other_lncRNA"
)]

# Identify host genes for intronic lncRNAs
# Find lncRNAs overlapping protein-coding genes
pc_gtf <- gtf_genes[gtf_genes$gene_type == "protein_coding"]
lnc_gr <- GRanges(lnc_annot$chr, IRanges(lnc_annot$start, lnc_annot$end),
                   gene_name = lnc_annot$gene_name)
pc_gr <- GRanges(seqnames(pc_gtf), ranges(pc_gtf), gene_name = pc_gtf$gene_name)

overlaps <- findOverlaps(lnc_gr, pc_gr, type = "within")
host_genes <- data.table(
  lncrna = lnc_gr$gene_name[queryHits(overlaps)],
  host_gene = pc_gr$gene_name[subjectHits(overlaps)]
)
host_genes <- host_genes[, .(host_gene = paste(unique(host_gene), collapse = ";")), by = lncrna]

# Get top 50 lncRNA DEGs by |t-statistic|
lnc_deg_all <- atlas[gene_biotype == "lncRNA" & is_deg == TRUE]
lnc_deg_all[, abs_tstat := abs(bulk_tstat)]
setorder(lnc_deg_all, -abs_tstat)
top50 <- head(lnc_deg_all, 50)

# Annotate
top50 <- merge(top50, tx_class[, .(gene_name, transcript_class)],
               by.x = "human_symbol", by.y = "gene_name", all.x = TRUE)
top50 <- merge(top50, host_genes, by.x = "human_symbol", by.y = "lncrna", all.x = TRUE)
top50[is.na(transcript_class), transcript_class := "unclassified"]

cat(sprintf("  Top 50 lncRNA DEGs annotated. Transcript class distribution:\n"))
print(top50[, .N, by = transcript_class][order(-N)])

# =============================================================================
# 5. Validate Known MASLD lncRNAs
# =============================================================================
cat("\n--- 5. Known MASLD lncRNA validation ---\n")

known_masld_lncrnas <- data.table(
  gene = c("NEAT1", "MALAT1", "MEG3", "H19", "GAS5", "HULC",
           "HOTAIR", "TUG1", "DANCR", "SNHG16", "CRNDE", "UCA1"),
  literature_direction = c("up", "up", "down", "up", "down", "up",
                           "up", "up", "up", "up", "up", "up"),
  literature_function = c(
    "NF-kB/inflammasome activation; lipid accumulation",
    "Lipid metabolism; hepatic stellate cell activation",
    "Hepatocyte apoptosis; miR-21 sponge",
    "Hepatic lipogenesis; insulin signaling",
    "mTOR inhibition; ferroptosis regulation",
    "Hepatocellular carcinoma; lipid metabolism",
    "Epigenetic silencing via PRC2; fibrosis",
    "Lipid accumulation; oxidative stress",
    "Hepatocyte proliferation; Wnt signaling",
    "miRNA sponge; NF-kB pathway",
    "Wnt/beta-catenin; proliferation",
    "PI3K/AKT; apoptosis resistance"
  )
)

# Match to atlas
known_validated <- merge(known_masld_lncrnas,
                          atlas[, .(human_symbol, bulk_logFC, bulk_padj, bulk_tstat,
                                    is_conserved, sources_active, gene_biotype)],
                          by.x = "gene", by.y = "human_symbol", all.x = TRUE)

known_validated[, in_atlas := !is.na(gene_biotype)]
known_validated[, is_deg := !is.na(bulk_padj) & bulk_padj < 0.1]
known_validated[, dream_direction := fifelse(bulk_logFC > 0, "up", "down")]
known_validated[, direction_concordant := (dream_direction == literature_direction)]

cat("  Known MASLD lncRNA validation:\n")
for (i in seq_len(nrow(known_validated))) {
  row <- known_validated[i]
  status <- ifelse(!row$in_atlas, "NOT_IN_ATLAS",
              ifelse(!row$is_deg, sprintf("NS (LFC=%.2f, p=%.3f)", row$bulk_logFC, row$bulk_padj),
                ifelse(row$direction_concordant,
                  sprintf("CONCORDANT (LFC=%.2f, p=%.2e)", row$bulk_logFC, row$bulk_padj),
                  sprintf("DISCORDANT (LFC=%.2f, p=%.2e)", row$bulk_logFC, row$bulk_padj))))
  cat(sprintf("  %-10s expected=%-5s %s\n", row$gene, row$literature_direction, status))
}

n_validated <- sum(known_validated$is_deg & known_validated$direction_concordant, na.rm = TRUE)
n_in_atlas <- sum(known_validated$in_atlas)
cat(sprintf("\n  Summary: %d/%d in atlas, %d/%d concordant DEGs\n",
            n_in_atlas, nrow(known_validated), n_validated, n_in_atlas))

# =============================================================================
# 6. Full ncRNA DEG Annotation Table
# =============================================================================
cat("\n--- 6. Building full ncRNA DEG annotation table ---\n")

ncrna_all <- atlas[gene_biotype %in% ncrna_biotypes]
cat(sprintf("  Total ncRNAs in atlas: %d\n", nrow(ncrna_all)))
cat(sprintf("  ncRNA DEGs (padj<0.1): %d\n", sum(ncrna_all$is_deg, na.rm = TRUE)))

# Add transcript class for lncRNAs
ncrna_all <- merge(ncrna_all, tx_class[, .(gene_name, transcript_class)],
                    by.x = "human_symbol", by.y = "gene_name", all.x = TRUE)
ncrna_all[gene_biotype != "lncRNA", transcript_class := gene_biotype]
ncrna_all[is.na(transcript_class), transcript_class := "unclassified"]

# Add host gene info
ncrna_all <- merge(ncrna_all, host_genes, by.x = "human_symbol", by.y = "lncrna", all.x = TRUE)

# Flag known MASLD lncRNAs
ncrna_all[, is_known_masld_lncrna := human_symbol %in% known_masld_lncrnas$gene]

# Add genomic coordinates
ncrna_all <- merge(ncrna_all, lnc_annot[, .(gene_name, chr, start, end, strand)],
                    by.x = "human_symbol", by.y = "gene_name", all.x = TRUE)

# Select key columns for output
deg_cols <- c("human_symbol", "ensembl_id", "gene_biotype", "transcript_class",
              "bulk_logFC", "bulk_padj", "bulk_tstat", "is_deg", "deg_direction",
              "host_gene", "is_known_masld_lncrna",
              "sources_active", "is_conserved",
              "chr", "start", "end", "strand",
              "mouse_ortholog", "mouse_meta_logFC", "mouse_meta_padj")

# Keep only columns that exist
deg_cols <- intersect(deg_cols, names(ncrna_all))
ncrna_annotated <- ncrna_all[, ..deg_cols]
setorder(ncrna_annotated, bulk_padj, na.last = TRUE)

fwrite(ncrna_annotated, file.path(out_dir, "ncrna_deg_annotated.csv"))
cat(sprintf("  Saved ncrna_deg_annotated.csv (%d rows)\n", nrow(ncrna_annotated)))

# =============================================================================
# 7. Transcript Class Distribution Among lncRNA DEGs
# =============================================================================
cat("\n--- 7. Transcript class distribution ---\n")

lnc_class_summary <- ncrna_all[gene_biotype == "lncRNA", .(
  n_total = .N,
  n_deg = sum(is_deg, na.rm = TRUE),
  frac_deg = sum(is_deg, na.rm = TRUE) / max(.N, 1),
  median_abs_lfc = median(abs(bulk_logFC[is_deg]), na.rm = TRUE)
), by = transcript_class][order(-n_total)]

cat("  lncRNA transcript class breakdown:\n")
print(lnc_class_summary)

# =============================================================================
# 8. miRNA Host Gene Analysis
# =============================================================================
cat("\n--- 8. miRNA host gene analysis ---\n")

# Find miRNA genes embedded in protein-coding or lncRNA host genes
mirna_gtf <- gtf_genes[gtf_genes$gene_type == "miRNA"]
mirna_gr <- GRanges(seqnames(mirna_gtf), ranges(mirna_gtf),
                     gene_name = mirna_gtf$gene_name)

# Overlap with protein-coding genes
mirna_pc_overlaps <- findOverlaps(mirna_gr, pc_gr, type = "within")
mirna_hosts <- data.table(
  mirna = mirna_gr$gene_name[queryHits(mirna_pc_overlaps)],
  host_gene = pc_gr$gene_name[subjectHits(mirna_pc_overlaps)],
  host_type = "protein_coding"
)

# Overlap with lncRNA genes
mirna_lnc_overlaps <- findOverlaps(mirna_gr, lnc_gr, type = "within")
if (length(mirna_lnc_overlaps) > 0) {
  mirna_lnc_hosts <- data.table(
    mirna = mirna_gr$gene_name[queryHits(mirna_lnc_overlaps)],
    host_gene = lnc_gr$gene_name[subjectHits(mirna_lnc_overlaps)],
    host_type = "lncRNA"
  )
  mirna_hosts <- rbind(mirna_hosts, mirna_lnc_hosts)
}

# Annotate with DEG status
mirna_hosts <- merge(mirna_hosts,
                      atlas[, .(human_symbol, bulk_padj, bulk_logFC, is_deg)],
                      by.x = "mirna", by.y = "human_symbol", all.x = TRUE)
mirna_hosts <- merge(mirna_hosts,
                      atlas[, .(human_symbol, dream_padj_host = bulk_padj,
                                dream_logFC_host = bulk_logFC,
                                is_deg_host = is_deg)],
                      by.x = "host_gene", by.y = "human_symbol", all.x = TRUE)

cat(sprintf("  miRNAs with host genes: %d\n", length(unique(mirna_hosts$mirna))))
cat(sprintf("  miRNA-host co-DEG pairs: %d\n",
            sum(mirna_hosts$is_deg == TRUE & mirna_hosts$is_deg_host == TRUE, na.rm = TRUE)))

# =============================================================================
# 9. Summary Statistics
# =============================================================================
cat("\n--- 9. Summary ---\n")

n_ncrna_total <- nrow(atlas[gene_biotype %in% ncrna_biotypes])
n_ncrna_tested <- nrow(atlas[gene_biotype %in% ncrna_biotypes & !is.na(bulk_padj)])
n_ncrna_deg <- sum(atlas[gene_biotype %in% ncrna_biotypes]$is_deg, na.rm = TRUE)
n_lncrna_deg <- sum(atlas[gene_biotype == "lncRNA"]$is_deg, na.rm = TRUE)
n_mirna_deg <- sum(atlas[gene_biotype == "miRNA"]$is_deg, na.rm = TRUE)

cat(sprintf("  Total ncRNAs in atlas: %d\n", n_ncrna_total))
cat(sprintf("  ncRNAs tested (dream): %d\n", n_ncrna_tested))
cat(sprintf("  ncRNA DEGs (padj<0.1): %d (%.1f%%)\n", n_ncrna_deg, 100 * n_ncrna_deg / max(n_ncrna_tested, 1)))
cat(sprintf("    lncRNA DEGs: %d\n", n_lncrna_deg))
cat(sprintf("    miRNA DEGs: %d\n", n_mirna_deg))
cat(sprintf("  Known MASLD lncRNAs validated: %d/%d\n", n_validated, n_in_atlas))

# Save known validation table
fwrite(known_validated, file.path(out_dir, "known_masld_lncrna_validation.csv"))

cat("\n=== Module 1 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
cat("Outputs:\n")
cat("  ", file.path(out_dir, "ncrna_landscape_summary.csv"), "\n")
cat("  ", file.path(out_dir, "ncrna_vs_pc_comparison.csv"), "\n")
cat("  ", file.path(out_dir, "ncrna_deg_annotated.csv"), "\n")
cat("  ", file.path(out_dir, "known_masld_lncrna_validation.csv"), "\n")
