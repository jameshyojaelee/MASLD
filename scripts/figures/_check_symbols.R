library(data.table)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]

nas <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_score_dream.csv"))
nas[, ensembl_clean := sub("\\..*", "", gene)]
nas <- merge(nas, unique(atlas[, .(ensembl_clean, symbol = human_symbol)], by = "ensembl_clean"),
             by = "ensembl_clean", all.x = TRUE)

cat("Total unique genes:", uniqueN(nas$gene), "\n")
unmapped_genes <- nas[is.na(symbol), unique(gene)]
cat("Unmapped (still Ensembl):", length(unmapped_genes), "\n")

# Check how many NAS sig DEGs are unmapped
sig <- nas[abs(logFC) > 0.5 & padj < 0.1]
sig_unmapped <- sig[is.na(symbol), unique(gene)]
cat("Sig DEGs unmapped:", length(sig_unmapped), "\n")

# Check top-variance unmapped genes (these would get used as padding labels)
all_degs <- unique(sig$gene)
mat_dt <- dcast(nas[gene %in% all_degs], gene ~ as.character(nas_level), value.var = "logFC")
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- mat_dt$gene
mat[is.na(mat)] <- 0
rv <- sort(apply(mat, 1, var), decreasing = TRUE)
cat("\nTop 30 by variance (would be padding labels):\n")
top30 <- names(rv)[1:30]
for (g in top30) {
  sym <- atlas[ensembl_clean == sub("\\..*", "", g), human_symbol]
  sym_str <- if (length(sym) > 0) sym[1] else "UNMAPPED"
  cat(g, "->", sym_str, "\n")
}

# Check target markers
targets <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
             "CXCL10", "CCL2", "TREM2", "TNF", "IL1B",
             "COL1A1", "TIMP1", "TGFB1", "ACTA2", "LUM",
             "FASN", "SCD", "PPARA", "CYP7A1", "AKR1B10")
cat("\nTarget marker presence in atlas:\n")
for (t in targets) {
  n <- nrow(atlas[human_symbol == t])
  cat(t, ":", n, "rows\n")
}
