#!/usr/bin/env Rscript
# WD per-dataset DE from Kallisto tximport counts
# Saves to per_diet/ directory for M03/M04 auto-discovery

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
})

PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
txi <- readRDS(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets/tximport/txi_wd.rds"))
DEDIR <- file.path(PROJECT, "RNA-seq/Mouse/Unified_Integration/results/per_diet")

cat("WD txi:", nrow(txi$counts), "genes x", ncol(txi$counts), "samples\n\n")

run_de <- function(counts, lengths, meta, disease_col, control_val, disease_vals, diet_name) {
  cat("============================================================\n")
  cat("  ", diet_name, "\n")
  cat("============================================================\n")

  meta$group <- ifelse(meta[[disease_col]] %in% disease_vals, "Disease",
                ifelse(meta[[disease_col]] %in% control_val, "Control", NA))
  meta <- meta[!is.na(meta$group), ]
  cat("  Disease:", sum(meta$group == "Disease"), "  Control:", sum(meta$group == "Control"), "\n")

  if (sum(meta$group == "Disease") < 2 | sum(meta$group == "Control") < 2) {
    cat("  SKIPPING: insufficient samples\n\n")
    return(NULL)
  }

  # Subset counts + lengths to matching samples
  shared <- intersect(meta$srr, colnames(counts))
  cat("  Matched samples:", length(shared), "\n")
  if (length(shared) < 4) { cat("  SKIPPING: too few matched\n"); return(NULL) }

  sub_counts <- round(counts[, shared])
  sub_lengths <- lengths[, shared]
  sub_meta <- meta[match(shared, meta$srr), ]

  group <- factor(sub_meta$group, levels = c("Control", "Disease"))
  dge <- DGEList(counts = sub_counts)
  keep <- filterByExpr(dge, group = group, min.count = 5, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("  Genes after filter:", nrow(dge), "\n")

  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("^group", "", colnames(design))

  v <- voom(dge, design, plot = FALSE)
  # Apply tximport length offsets
  len_sub <- sub_lengths[rownames(dge), shared]
  len_sub[is.na(len_sub) | len_sub == 0] <- 1
  v$offset <- log(len_sub)
  cat("  [offset] Applied tximport length offsets\n")

  fit <- lmFit(v, design)
  contrasts <- makeContrasts(Disease - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res <- as.data.table(res)
  setnames(res, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene"))
  setcolorder(res, "gene")

  sig <- res[adj.P.Val < 0.05]
  cat("  DEGs (padj<0.05):", nrow(sig), " (Up:", nrow(sig[logFC > 0]),
      "Down:", nrow(sig[logFC < 0]), ")\n")

  out_file <- file.path(DEDIR, paste0(diet_name, "_de_results.csv"))
  fwrite(res[order(adj.P.Val)], out_file)
  cat("  Saved:", out_file, "\n\n")
  return(res)
}

# --- Load GSM-to-SRR mappings ---
load_gsm_srr <- function(ds) {
  fread(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets", ds, "metadata/gsm_to_srr.tsv"))
}

# ============== GSE220575 (DIAMOND) ==============
meta220 <- fread(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets/GSE220575/metadata/sample_metadata.csv"))
gsm_srr220 <- load_gsm_srr("GSE220575")
meta220 <- merge(meta220, gsm_srr220, by.x = "gsm", by.y = "gsm", all.x = TRUE)
# Keep only samples in txi
meta220 <- meta220[srr %in% colnames(txi$counts)]
run_de(txi$counts, txi$length, meta220, "condition", "Control", c("MASH", "Fatty_Liver", "MASH/FL"),
       "DIAMOND")

# ============== GSE246088 (WD) ==============
meta246 <- fread(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/metadata/sample_metadata.csv"))
gsm_srr246 <- load_gsm_srr("GSE246088")
meta246 <- merge(meta246, gsm_srr246, by.x = "geo_accession", by.y = "gsm", all.x = TRUE)
meta246 <- meta246[srr %in% colnames(txi$counts)]
# WT only, WD vs Chow
meta246_wt <- meta246[genotype == "Plvap_WT" | genotype == "WT" | genotype == "Ctrl"]
run_de(txi$counts, txi$length, meta246_wt, "diet", c("Chow", "chow"), c("WD", "Western", "western_diet"),
       "Western_Diet")

# ============== GSE305484 (WD + fructose) ==============
meta305 <- fread(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets/GSE305484/metadata/sample_metadata.csv"))
gsm_srr305 <- load_gsm_srr("GSE305484")
meta305 <- merge(meta305, gsm_srr305, by.x = "gsm_accession", by.y = "gsm", all.x = TRUE)
meta305 <- meta305[srr %in% colnames(txi$counts)]
# WT only
meta305_wt <- meta305[genotype == "CD163WT" | genotype == "WT"]
run_de(txi$counts, txi$length, meta305_wt, "treatment", c("Chow", "chow"), c("WD", "Western", "WD+fructose"),
       "Western_Diet_Fructose")

# ============== GSE292565 (FFC diet) ==============
meta292 <- fread(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets/GSE292565/metadata/sample_metadata.csv"))
gsm_srr292 <- load_gsm_srr("GSE292565")
meta292 <- merge(meta292, gsm_srr292, by.x = "gsm", by.y = "gsm", all.x = TRUE)
meta292 <- meta292[srr %in% colnames(txi$counts)]
run_de(txi$counts, txi$length, meta292, "condition", "Control", c("Disease"),
       "FFC")

# ============== GSE246328 (GAN diet) ==============
meta246g <- fread(file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets/GSE246328/metadata/sample_metadata.csv"))
gsm_srr246g <- load_gsm_srr("GSE246328")
meta246g <- merge(meta246g, gsm_srr246g, by.x = "gsm", by.y = "gsm", all.x = TRUE)
meta246g <- meta246g[srr %in% colnames(txi$counts)]
# Exclude Chow-reversal arm; keep Control + Disease only
meta246g_de <- meta246g[condition %in% c("Control", "Disease")]
run_de(txi$counts, txi$length, meta246g_de, "condition", "Control", c("Disease"),
       "GAN")

cat("WD per-diet DE complete.\n")
