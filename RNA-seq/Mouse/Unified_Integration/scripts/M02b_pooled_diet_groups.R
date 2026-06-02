#!/usr/bin/env Rscript
# M02b_pooled_diet_groups.R
# ---------------------------------------------------------------------------
# Pooled per-diet-GROUP DE for the 3 multi-dataset diet groups:
#   1. NASH_diet = FPC (GSE162876) + FFC (GSE292565) + GAN (GSE246328)
#   2. Western_diet = DIAMOND (GSE220575) + WD (GSE246088) + WD_Fructose (GSE305484)
#   3. HFD = GSE224069 + GSE274914 + GSE246088 HFD arm  [folded in 2026-05-29 per PI]
#
# HFD note: GSE224069/GSE274914 come from the main integration (M01 merged counts);
# the GSE246088 HFD arm (5 Disease + 6 Chow controls, Plvap_Control) comes from the
# Western_Diet_Datasets featureCounts. This OVERWRITES the M02 HFD_de_results.csv with
# the pooled 3-dataset version, so M02b MUST run after M02. Control scheme is
# dataset-matched (same as M02 and as NASH/Western); the only change vs M02's HFD is
# the added GSE246088 arm. (Its Chow controls are shared with the GSE246088 WD arm,
# which is a separate DE contrast.)
#
# Supports BOTH featureCounts (STAR) and tximport (Kallisto) modes.
# Mode is auto-detected: if txi.rds exists → tximport; otherwise → featureCounts.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEDIR   <- file.path(PROJECT, "RNA-seq/Mouse/Unified_Integration/results/per_diet")

# Guard: M02 must have run first (it produces the HFD file M02b overwrites).
m02_hfd <- file.path(DEDIR, "HFD_de_results.csv")
if (!file.exists(m02_hfd))
  stop("M02b requires M02 output at ", m02_hfd, " -- run M02_mouse_per_diet_de.R first.")
WD_DIR  <- file.path(PROJECT, "RNA-seq/Mouse/Western_Diet_Datasets")

txi_file <- file.path(PROJECT, "RNA-seq/Mouse/Unified_Integration/counts/tximport/txi.rds")
USE_TXIMPORT <- file.exists(txi_file)
cat("Quantification mode:", ifelse(USE_TXIMPORT, "tximport (Kallisto)", "featureCounts (STAR)"), "\n\n")

# ---------------------------------------------------------------------------
# Helper: load featureCounts gene_counts.txt → named matrix (genes × samples)
# ---------------------------------------------------------------------------
load_featurecounts <- function(path) {
  dt <- fread(path, skip = "Geneid")
  genes <- sub("[.][0-9]+$", "", dt$Geneid)
  sample_cols <- setdiff(names(dt), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  mat <- as.matrix(dt[, ..sample_cols])
  rownames(mat) <- genes
  srr_ids <- basename(dirname(sample_cols))
  colnames(mat) <- srr_ids
  mat
}

# ---------------------------------------------------------------------------
# Helper: run pooled limma-voom DE with batch covariate
# ---------------------------------------------------------------------------
run_pooled_de <- function(counts, meta, batch_col, diet_name,
                          lengths = NULL) {
  cat("============================================================\n")
  cat("  ", diet_name, "\n")
  cat("============================================================\n")

  shared <- intersect(meta$srr, colnames(counts))
  meta <- meta[match(shared, meta$srr), ]
  cat("  Samples:", length(shared), "(Disease:", sum(meta$group == "Disease"),
      "Control:", sum(meta$group == "Control"), ")\n")
  cat("  Datasets:", paste(unique(meta[[batch_col]]), collapse = ", "), "\n")

  if (sum(meta$group == "Disease") < 2 | sum(meta$group == "Control") < 2) {
    cat("  SKIPPING: insufficient samples\n\n")
    return(NULL)
  }

  sub_c <- round(counts[, shared, drop = FALSE])
  group <- factor(meta$group, levels = c("Control", "Disease"))
  batch <- factor(meta[[batch_col]])

  dge <- DGEList(counts = sub_c)
  keep <- filterByExpr(dge, group = group, min.count = 5, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("  Genes after filter:", nrow(dge), "\n")

  if (length(unique(batch)) > 1) {
    design <- model.matrix(~ 0 + group + batch)
  } else {
    design <- model.matrix(~ 0 + group)
  }
  colnames(design) <- gsub("^group", "", colnames(design))

  v <- voom(dge, design, plot = FALSE)

  if (!is.null(lengths)) {
    ll <- lengths[rownames(dge), shared, drop = FALSE]
    ll[is.na(ll) | ll == 0] <- 1
    v$offset <- log(ll)
    cat("  [offset] Applied tximport length offsets\n")
  } else {
    cat("  [no offset] featureCounts mode — plain voom\n")
  }

  fit <- lmFit(v, design)
  contrasts <- makeContrasts(Disease - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res <- as.data.table(res)
  setnames(res, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene"))
  setcolorder(res, "gene")

  # Add unmoderated SE + df_total so M02c_ashr_shrinkage.R uses SE_unmoderated
  # consistently across all 4 diet groups (same as M02_mouse_per_diet_de.R:222-223).
  res[, SE_unmoderated := as.numeric(fit2$stdev.unscaled[, 1]) * fit2$sigma]
  res[, df_total       := fit2$df.total]

  sig <- res[adj.P.Val < 0.05]
  cat("  DEGs (padj<0.05):", nrow(sig),
      " (Up:", nrow(sig[logFC > 0]), "Down:", nrow(sig[logFC < 0]), ")\n")

  out_file <- file.path(DEDIR, paste0(diet_name, "_de_results.csv"))
  fwrite(res[order(adj.P.Val)], out_file)
  cat("  Saved:", out_file, "\n\n")
  return(res)
}

# ===========================================================================
# Load count data (branching on USE_TXIMPORT)
# ===========================================================================

meta_main <- readRDS(file.path(PROJECT,
    "RNA-seq/Mouse/Unified_Integration/results/meta_matched.rds"))

if (USE_TXIMPORT) {
  suppressPackageStartupMessages(library(tximport))
  tx2gene <- fread(file.path(PROJECT, "data/reference/kallisto/tx2gene_vM38.tsv"),
                   header = TRUE)
  setnames(tx2gene, c("tx_id", "gene_id", "gene_name"))
  tx2gene$tx_id   <- sub("[.][0-9]+$", "", tx2gene$tx_id)
  tx2gene$gene_id <- sub("[.][0-9]+$", "", tx2gene$gene_id)

  txi_main <- readRDS(txi_file)
  cat("Main txi:", ncol(txi_main$counts), "samples\n")

  txi_wd_file <- file.path(WD_DIR, "tximport/txi_wd.rds")
  if (file.exists(txi_wd_file)) {
    txi_wd <- readRDS(txi_wd_file)
    cat("WD txi:", ncol(txi_wd$counts), "samples\n")
  } else {
    cat("WD txi not found — building from individual quants\n")
    all_files <- c(); all_names <- c()
    for (ds in c("GSE220575", "GSE246088", "GSE305484", "GSE246328", "GSE292565")) {
      qd <- file.path(WD_DIR, ds, "quant/kallisto")
      if (!dir.exists(qd)) next
      for (s in list.dirs(qd, recursive = FALSE, full.names = FALSE)) {
        f <- file.path(qd, s, "abundance.tsv")
        if (file.exists(f)) { all_files <- c(all_files, f); all_names <- c(all_names, s) }
      }
    }
    names(all_files) <- all_names
    cat("Building WD txi from", length(all_files), "quant files\n")
    txi_wd <- tximport(all_files, type = "kallisto",
                       tx2gene = tx2gene[, .(tx_id, gene_id)],
                       countsFromAbundance = "no", ignoreTxVersion = TRUE)
  }

  main_counts  <- txi_main$counts
  main_lengths <- txi_main$length
  wd_counts    <- txi_wd$counts
  wd_lengths   <- txi_wd$length

} else {
  # featureCounts mode: load per-dataset gene_counts.txt from WD datasets
  cat("Loading featureCounts from WD datasets...\n")
  wd_mats <- list()
  for (ds in c("GSE220575", "GSE246088", "GSE246328", "GSE292565", "GSE305484")) {
    fc_file <- file.path(WD_DIR, ds, "counts/featurecounts/gene_counts.txt")
    if (file.exists(fc_file)) {
      wd_mats[[ds]] <- load_featurecounts(fc_file)
      cat("  ", ds, ":", ncol(wd_mats[[ds]]), "samples,", nrow(wd_mats[[ds]]), "genes\n")
    } else {
      cat("  ", ds, ": gene_counts.txt not found — SKIPPING\n")
    }
  }

  # FPC counts come from the M01 merged matrix
  merged_rds <- file.path(PROJECT,
      "RNA-seq/Mouse/Unified_Integration/results/merged_counts_raw.rds")
  main_counts <- readRDS(merged_rds)
  rownames(main_counts) <- sub("[.][0-9]+$", "", rownames(main_counts))
  cat("Main merged counts:", ncol(main_counts), "samples,", nrow(main_counts), "genes\n")

  # Merge all WD featureCounts into a single matrix
  if (length(wd_mats) > 0) {
    common_genes <- Reduce(intersect, lapply(wd_mats, rownames))
    common_genes <- intersect(common_genes, rownames(main_counts))
    wd_counts <- do.call(cbind, lapply(wd_mats, function(m) m[common_genes, , drop = FALSE]))
    main_counts <- main_counts[common_genes, , drop = FALSE]
    cat("Common gene universe:", length(common_genes), "genes\n")
  } else {
    stop("No WD featureCounts files found")
  }
  main_lengths <- NULL
  wd_lengths   <- NULL
}

# ===========================================================================
# Shared gene universe
# ===========================================================================
shared_genes <- intersect(rownames(main_counts), rownames(wd_counts))
main_counts  <- main_counts[shared_genes, , drop = FALSE]
wd_counts    <- wd_counts[shared_genes, , drop = FALSE]
if (!is.null(main_lengths)) {
  main_lengths <- main_lengths[shared_genes, , drop = FALSE]
  wd_lengths   <- wd_lengths[shared_genes, , drop = FALSE]
}
cat("Shared gene universe:", length(shared_genes), "\n\n")

# ===========================================================================
# 1. NASH diet = FPC + FFC + GAN
# ===========================================================================

# FPC arm: select by DATASET to capture its 73 shared controls. FPC controls are
# labeled diet_model == "Control" (shared with CDAHFD), so a diet_model == "FPC"
# filter returns disease-only (bug fixed 2026-05-29). Exclude the CDAHFD disease
# arm of this dual-model dataset (diet_model == "CDAHFD").
fpc_meta <- meta_main[dataset == "GSE162876" &
                      (diet_model == "FPC" | group_binary == "Control"),
                      .(srr = sample_id, group_binary, dataset)]

gan_meta_file <- file.path(WD_DIR, "GSE246328/metadata/sample_metadata.csv")
gan_meta <- if (file.exists(gan_meta_file)) {
  g <- fread(gan_meta_file)
  g <- g[condition %in% c("Control", "Disease")]
  g <- g[srr %in% colnames(wd_counts)]
  g[, dataset := "GSE246328"]
  g[, group_binary := condition]
  g
} else data.table()

ffc_meta_file <- file.path(WD_DIR, "GSE292565/metadata/sample_metadata.csv")
ffc_meta <- if (file.exists(ffc_meta_file)) {
  f <- fread(ffc_meta_file)
  if (!"srr" %in% names(f)) {
    ffc_gsm <- fread(file.path(WD_DIR, "GSE292565/metadata/gsm_to_srr.tsv"))
    f <- merge(f, ffc_gsm, by.x = "gsm", by.y = "gsm", all.x = TRUE)
  }
  f <- f[srr %in% colnames(wd_counts)]
  f[, dataset := "GSE292565"]
  f[, group_binary := condition]
  f
} else data.table()

nash_meta_rows <- list()
for (i in seq_len(nrow(fpc_meta))) {
  s <- fpc_meta$srr[i]
  if (s %in% colnames(main_counts))
    nash_meta_rows[[length(nash_meta_rows) + 1]] <- data.table(
      srr = s, group = fpc_meta$group_binary[i], dataset = "GSE162876_FPC")
}
for (i in seq_len(nrow(gan_meta))) {
  nash_meta_rows[[length(nash_meta_rows) + 1]] <- data.table(
    srr = gan_meta$srr[i], group = gan_meta$group_binary[i], dataset = "GSE246328_GAN")
}
for (i in seq_len(nrow(ffc_meta))) {
  nash_meta_rows[[length(nash_meta_rows) + 1]] <- data.table(
    srr = ffc_meta$srr[i], group = ffc_meta$group_binary[i], dataset = "GSE292565_FFC")
}
nash_meta <- rbindlist(nash_meta_rows)
nash_meta[, group := ifelse(group %in% c("Disease", "disease"), "Disease", "Control")]
cat("NASH diet pool:", nrow(nash_meta), "samples\n")

nash_srrs_main <- intersect(nash_meta$srr, colnames(main_counts))
nash_srrs_wd   <- intersect(nash_meta$srr, colnames(wd_counts))
nash_counts <- cbind(
  main_counts[, nash_srrs_main, drop = FALSE],
  wd_counts[, nash_srrs_wd, drop = FALSE]
)
nash_lengths <- if (!is.null(main_lengths)) {
  cbind(main_lengths[, nash_srrs_main, drop = FALSE],
        wd_lengths[, nash_srrs_wd, drop = FALSE])
} else NULL

run_pooled_de(nash_counts, nash_meta, "dataset", "NASH_diet",
              lengths = nash_lengths)

# ===========================================================================
# 2. Western diet = DIAMOND + WD + WD_Fructose
# ===========================================================================

wd_meta_rows <- list()

dm_meta_file <- file.path(WD_DIR, "GSE220575/metadata/sample_metadata.csv")
if (file.exists(dm_meta_file)) {
  dm <- fread(dm_meta_file)
  if (!"srr" %in% names(dm)) {
    dm_gsm <- fread(file.path(WD_DIR, "GSE220575/metadata/gsm_to_srr.tsv"))
    dm <- merge(dm, dm_gsm, by.x = "gsm", by.y = "gsm", all.x = TRUE)
  }
  dm <- dm[srr %in% colnames(wd_counts) & condition %in% c("Control", "MASH")]
  for (i in seq_len(nrow(dm)))
    wd_meta_rows[[length(wd_meta_rows) + 1]] <- data.table(
      srr = dm$srr[i],
      group = ifelse(dm$condition[i] == "Control", "Control", "Disease"),
      dataset = "GSE220575_DIAMOND")
}

wd246_meta_file <- file.path(WD_DIR, "GSE246088/metadata/sample_metadata.csv")
if (file.exists(wd246_meta_file)) {
  wd246 <- fread(wd246_meta_file)
  if (!"srr" %in% names(wd246)) {
    wd246_gsm <- fread(file.path(WD_DIR, "GSE246088/metadata/gsm_to_srr.tsv"))
    wd246 <- merge(wd246, wd246_gsm, by.x = "geo_accession", by.y = "gsm", all.x = TRUE)
  }
  wd246 <- wd246[genotype == "Plvap_Control" & diet %in% c("Chow", "Western_Diet")]
  wd246 <- wd246[srr %in% colnames(wd_counts)]
  for (i in seq_len(nrow(wd246)))
    wd_meta_rows[[length(wd_meta_rows) + 1]] <- data.table(
      srr = wd246$srr[i],
      group = ifelse(wd246$diet[i] == "Chow", "Control", "Disease"),
      dataset = "GSE246088_WD")
}

wd305_meta_file <- file.path(WD_DIR, "GSE305484/metadata/sample_metadata.csv")
if (file.exists(wd305_meta_file)) {
  wd305 <- fread(wd305_meta_file)
  if (!"srr" %in% names(wd305)) {
    wd305_gsm <- fread(file.path(WD_DIR, "GSE305484/metadata/gsm_to_srr.tsv"))
    wd305 <- merge(wd305, wd305_gsm, by.x = "gsm_accession", by.y = "gsm", all.x = TRUE)
  }
  wd305 <- wd305[genotype == "CD163WT" & srr %in% colnames(wd_counts)]
  for (i in seq_len(nrow(wd305)))
    wd_meta_rows[[length(wd_meta_rows) + 1]] <- data.table(
      srr = wd305$srr[i],
      group = ifelse(wd305$condition[i] == "Chow", "Control", "Disease"),
      dataset = "GSE305484_WD")
}

wd_meta_all <- rbindlist(wd_meta_rows)
cat("\nWestern diet pool:", nrow(wd_meta_all), "samples\n")

wd_diet_counts <- wd_counts[, intersect(wd_meta_all$srr, colnames(wd_counts)), drop = FALSE]
wd_diet_lengths <- if (!is.null(wd_lengths)) {
  wd_lengths[, intersect(wd_meta_all$srr, colnames(wd_lengths)), drop = FALSE]
} else NULL

run_pooled_de(wd_diet_counts, wd_meta_all, "dataset", "Western_diet",
              lengths = wd_diet_lengths)

# ===========================================================================
# 3. HFD = GSE224069 + GSE274914 (main integration) + GSE246088 HFD arm
# ===========================================================================

hfd_meta_rows <- list()

# GSE224069 + GSE274914: select by dataset (HFD controls are labeled
# diet_model == "Control" in meta_main, so a diet_model=="HFD" filter would
# miss them); this reproduces M02's dataset-matched HFD comparison.
hfd_main <- meta_main[dataset %in% c("GSE224069", "GSE274914"),
                      .(srr = sample_id, group_binary, dataset)]
for (i in seq_len(nrow(hfd_main))) {
  s <- hfd_main$srr[i]
  if (s %in% colnames(main_counts))
    hfd_meta_rows[[length(hfd_meta_rows) + 1]] <- data.table(
      srr = s, group = hfd_main$group_binary[i], dataset = hfd_main$dataset[i])
}

# GSE246088 HFD arm: Plvap_Control, High_Fat_Diet = Disease, Chow = Control
hfd246_file <- file.path(WD_DIR, "GSE246088/metadata/sample_metadata.csv")
if (file.exists(hfd246_file)) {
  hfd246 <- fread(hfd246_file)
  if (!"srr" %in% names(hfd246)) {
    hfd246_gsm <- fread(file.path(WD_DIR, "GSE246088/metadata/gsm_to_srr.tsv"))
    hfd246 <- merge(hfd246, hfd246_gsm, by.x = "geo_accession", by.y = "gsm", all.x = TRUE)
  }
  hfd246 <- hfd246[genotype == "Plvap_Control" & diet %in% c("Chow", "High_Fat_Diet")]
  hfd246 <- hfd246[srr %in% colnames(wd_counts)]
  for (i in seq_len(nrow(hfd246)))
    hfd_meta_rows[[length(hfd_meta_rows) + 1]] <- data.table(
      srr = hfd246$srr[i],
      group = ifelse(hfd246$diet[i] == "Chow", "Control", "Disease"),
      dataset = "GSE246088_HFD")
}

hfd_meta_all <- rbindlist(hfd_meta_rows)
hfd_meta_all[, group := ifelse(group %in% c("Disease", "disease"), "Disease", "Control")]
cat("\nHFD pool:", nrow(hfd_meta_all), "samples (datasets:",
    paste(sort(unique(hfd_meta_all$dataset)), collapse = ", "), ")\n")

hfd_srrs_main <- intersect(hfd_meta_all$srr, colnames(main_counts))
hfd_srrs_wd   <- intersect(hfd_meta_all$srr, colnames(wd_counts))
hfd_counts <- cbind(
  main_counts[, hfd_srrs_main, drop = FALSE],
  wd_counts[, hfd_srrs_wd, drop = FALSE]
)
hfd_lengths <- if (!is.null(main_lengths)) {
  cbind(main_lengths[, hfd_srrs_main, drop = FALSE],
        wd_lengths[, hfd_srrs_wd, drop = FALSE])
} else NULL

run_pooled_de(hfd_counts, hfd_meta_all, "dataset", "HFD",
              lengths = hfd_lengths)

# ===========================================================================
# 4. COMBINED Western = NASH_diet + Western_diet pooled (PI decision 2026-05-29)
#    FPC/FFC/GAN + DIAMOND/WD are one Western / metabolic-overload family
#    (verified compositions: all 40-52% fat + fructose ± cholesterol 0.15-2%).
#    This "Western" group REPLACES the separate NASH_diet/Western_diet as the
#    canonical diet axis (those are kept above for provenance/sensitivity).
#    6 datasets: GSE162876_FPC + GSE246328_GAN + GSE292565_FFC +
#    GSE220575_DIAMOND + GSE246088_WD + GSE305484_WD.
# ===========================================================================
western_meta <- rbind(nash_meta, wd_meta_all)
western_meta[, group := ifelse(group %in% c("Disease", "disease"), "Disease", "Control")]
cat("\nCombined Western pool:", nrow(western_meta), "samples (datasets:",
    paste(sort(unique(western_meta$dataset)), collapse = ", "), ")\n")

w_srrs_main <- intersect(western_meta$srr, colnames(main_counts))
w_srrs_wd   <- intersect(western_meta$srr, colnames(wd_counts))
western_counts <- cbind(
  main_counts[, w_srrs_main, drop = FALSE],
  wd_counts[, w_srrs_wd, drop = FALSE]
)
western_lengths <- if (!is.null(main_lengths)) {
  cbind(main_lengths[, w_srrs_main, drop = FALSE],
        wd_lengths[, w_srrs_wd, drop = FALSE])
} else NULL

run_pooled_de(western_counts, western_meta, "dataset", "Western",
              lengths = western_lengths)

cat("\n=== Pooled diet group DE complete ===\n")
cat("Mode:", ifelse(USE_TXIMPORT, "tximport (Kallisto)", "featureCounts (STAR)"), "\n")
