#!/usr/bin/env Rscript
# Pseudobulk DE on consecutive coarse stage contrasts (proxies F-transitions).
#
# Contrasts:
#   1. Steatosis vs Healthy           (proxies F0->F1 / F1->F2 entry)
#   2. Steatohepatitis vs Steatosis   (proxies F2->F3)
#   3. Cirrhosis vs Steatohepatitis   (proxies F3->F4)
#
# Inputs:
#   results_gpu_v2/integrated_atlas.h5ad            (cell metadata for stage mapping)
#   results_gpu_v2/pseudobulk/{CellType}_pseudobulk.csv  (genes x samples raw counts)
#   data/GSE244832/metadata/donor_pairing.csv       (recovers GSE244832 stages)
#
# Outputs:
#   results_gpu_v2/pseudobulk_de/{CellType}_{Contrast}_de.csv
#   Schema matches existing _de.csv files:
#     logFC, AveExpr, t_stat, pvalue, padj, B, gene, cell_type, contrast
#
# Filters: same conventions as pseudobulk_de.R
#   - preparation_method in c("nuclei", "cd45_negative", "unsorted")
#   - >= 50 cells per sample contributing to pseudobulk
#   - >= 2 samples per stage
#   - >= 5 counts in >= 3 samples for gene retention
#   - dataset as covariate when >1 dataset

suppressPackageStartupMessages({
  library(rhdf5)
  library(limma)
  library(edgeR)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
scvi_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
h5ad_file <- file.path(scvi_dir, "integrated_atlas.h5ad")
pb_dir <- file.path(scvi_dir, "pseudobulk")
# Donor-collapsed output (fixes run-level pseudoreplication). Written to a NEW
# directory so the legacy run-level coarse _de.csv files are not clobbered.
out_dir <- file.path(scvi_dir, "pseudobulk_de/coarse_stage_donorcollapsed")
donor_pairing_file <- file.path(BASE, "data/GSE244832/metadata/donor_pairing.csv")
gse202379_pairing_file <- file.path(BASE, "data/GSE202379/metadata/donor_pairing.csv")

# Shared donor-collapse utility (SRR run -> biological donor for GSE244832 /
# GSE185477 / GSE202379). Makes donor-level DE the default going forward.
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Cell types of interest (match existing _de.csv files)
TARGET_CELL_TYPES <- c(
  "Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes",
  "Endothelial_cells", "T_cells", "B_cells", "Mono+mono_derived_cells",
  "Plasma_cells", "Resident_NK", "Circulating_NK_NKT"
)

# Contrast definitions: c(reference, treatment)
CONTRASTS <- list(
  Steatosis_vs_Healthy        = c("Healthy",         "Steatosis"),
  Steatohepatitis_vs_Steatosis = c("Steatosis",      "Steatohepatitis"),
  Cirrhosis_vs_Steatohepatitis = c("Steatohepatitis", "Cirrhosis")
)

# ---------------------------------------------------------------------------
# 1. Build cell metadata from h5ad and map to disease_stage_coarse
# ---------------------------------------------------------------------------
message("Reading h5ad obs from ", h5ad_file)
stopifnot(file.exists(h5ad_file))

items <- h5ls(h5ad_file)

read_categorical <- function(h5file, key) {
  group_path <- paste0("/obs/", key)
  obs_groups <- subset(items, group == "/obs" & otype == "H5I_GROUP")$name
  if (key %in% obs_groups) {
    codes <- h5read(h5file, paste0(group_path, "/codes"))
    cats  <- h5read(h5file, paste0(group_path, "/categories"))
    return(cats[codes + 1L])
  }
  vals <- h5read(h5file, group_path)
  if (is.list(vals) && "codes" %in% names(vals)) {
    return(vals$categories[vals$codes + 1L])
  }
  as.vector(vals)
}

sample_vec  <- read_categorical(h5ad_file, "sample")
dataset_vec <- read_categorical(h5ad_file, "dataset")
cond_vec    <- read_categorical(h5ad_file, "condition")
prep_vec    <- tryCatch(
  read_categorical(h5ad_file, "preparation_method"),
  error = function(e) rep("unsorted", length(sample_vec))
)

cell_meta <- data.table(
  sample = sample_vec,
  dataset = dataset_vec,
  condition = cond_vec,
  preparation_method = prep_vec
)

# Collapse to RUN-level metadata via majority vote on each field.
# NOTE: `sample` here is a SEQUENCING RUN (SRR), not a biological donor, for
# GSE244832 / GSE185477. Donor-collapse + donor-level cell counts happen below,
# after per-SRR stage recovery (which needs the raw SRR ids).
sample_meta <- cell_meta[, .(
  dataset = dataset[1],
  condition = names(sort(table(condition), decreasing = TRUE))[1],
  preparation_method = preparation_method[1]
), by = sample]

# Stage mapping (mirrors extract_atlas_umap_for_fig2.py)
cond_to_stage <- c(
  "Healthy" = "Healthy", "Mixed" = "Healthy",
  "NAFLD" = "Steatosis", "MASL" = "Steatosis",
  "NASH" = "Steatohepatitis", "MASH" = "Steatohepatitis",
  "Cirrhotic" = "Cirrhosis"
)
sample_meta[, disease_stage_coarse := cond_to_stage[condition]]

# Recover GSE202379 from its published diagnosis field. Fibrosis stage is a
# separate axis and must not be used to redefine NAFLD/NASH diagnosis groups.
if (file.exists(gse202379_pairing_file)) {
  dp <- fread(gse202379_pairing_file)
  g202_map <- c(
    "Healthy control" = "Healthy",
    "NAFLD" = "Steatosis",
    "NASH w/o cirrhosis" = "Steatohepatitis",
    "NASH with cirrhosis" = "Cirrhosis",
    "end stage" = "Cirrhosis"
  )
  for (i in seq_len(nrow(dp))) {
    stage_label <- unname(g202_map[dp$condition[i]])
    if (is.na(stage_label)) next
    runs <- trimws(strsplit(dp$rna_srrs[i], ";", fixed = TRUE)[[1]])
    sample_meta[dataset == "GSE202379" & sample %in% runs,
                disease_stage_coarse := stage_label]
  }
}

# Recover GSE244832 stages via SRR -> condition mapping
if (file.exists(donor_pairing_file)) {
  dp <- fread(donor_pairing_file)
  g244_map <- c("NORMAL" = "Healthy", "MASL" = "Steatosis", "MASH" = "Steatohepatitis")
  srr_to_stage <- list()
  for (i in seq_len(nrow(dp))) {
    cond_lbl <- g244_map[dp$condition[i]]
    if (is.na(cond_lbl) || is.na(dp$rna_srrs[i]) || dp$rna_srrs[i] == "") next
    for (srr in strsplit(dp$rna_srrs[i], ";")[[1]]) {
      srr <- trimws(srr)
      if (nzchar(srr)) srr_to_stage[[srr]] <- cond_lbl
    }
  }
  g244_idx <- which(sample_meta$dataset == "GSE244832" &
                    is.na(sample_meta$disease_stage_coarse))
  for (i in g244_idx) {
    s <- srr_to_stage[[sample_meta$sample[i]]]
    if (!is.null(s)) sample_meta$disease_stage_coarse[i] <- s
  }
  message("  GSE244832 stages recovered: ",
          sum(sample_meta$dataset == "GSE244832" &
              !is.na(sample_meta$disease_stage_coarse)),
          " / ",
          sum(sample_meta$dataset == "GSE244832"))
}

# ---------------------------------------------------------------------------
# 1b. Donor-collapse (fix run-level pseudoreplication)
# ---------------------------------------------------------------------------
# `sample` is a sequencing RUN, not a donor, for GSE244832 (117 runs -> 18
# donors) and GSE185477 (21 runs -> 3 donors); all other datasets are 1 run =
# 1 donor and pass through unchanged. The per-run stage recovery above MUST run
# first because it keys on raw SRR ids.
srr_to_donor <- build_srr_to_donor_map(BASE)

# Collapse run-level sample_meta to donor level (majority vote per field). The
# per-run stage recovery above runs FIRST (it keys on raw SRR ids), so the vote
# resolves each donor's stage from its runs; the hardened util returns NA for an
# all-NA field (GSE202379, whose generic "MASLD" has no coarse-stage mapping).
n_runs_pre <- nrow(sample_meta)
sample_meta <- collapse_sample_meta_to_donor(sample_meta, BASE)
message(sprintf("  donor-collapse: %d run-level rows -> %d donor-level rows",
                n_runs_pre, nrow(sample_meta)))

# Donor-level cell counts for the >=50-cell filter: sum each donor's per-run
# cell counts (a run-level n_cells >=50 filter would wrongly drop donors whose
# reads are split across many small runs).
run_cell_counts <- cell_meta[, .N, by = sample]
run_cell_counts[, donor := ifelse(sample %in% names(srr_to_donor),
                                  srr_to_donor[sample], sample)]
donor_cell_counts <- run_cell_counts[, .(n_cells = sum(N)), by = donor]
setnames(donor_cell_counts, "donor", "sample")
sample_meta <- merge(sample_meta, donor_cell_counts, by = "sample", all.x = TRUE)

# Apply prep-method filter consistently (same as pseudobulk_de.R)
de_preps <- c("nuclei", "cd45_negative", "unsorted")
n_before <- nrow(sample_meta)
sample_meta_de <- sample_meta[preparation_method %in% de_preps &
                              !is.na(disease_stage_coarse)]
message(sprintf("  prep-method + stage filter: %d -> %d samples",
                n_before, nrow(sample_meta_de)))
message("  Stage breakdown (DE-eligible samples):")
print(table(sample_meta_de$disease_stage_coarse))
message("  Dataset x stage cross-tab (DE-eligible):")
print(table(sample_meta_de$dataset, sample_meta_de$disease_stage_coarse))

# ---------------------------------------------------------------------------
# 2. Helper: run limma-voom for one cell type x contrast
# ---------------------------------------------------------------------------
run_de <- function(counts_mat_full, ct_name, contrast_name, ref_stage, trt_stage) {
  # counts_mat_full: genes x samples (full pseudobulk matrix for cell type)
  message(sprintf("  [%s] %s vs %s", ct_name, trt_stage, ref_stage))

  col_samples <- colnames(counts_mat_full)
  matched_idx <- match(col_samples, sample_meta_de$sample)
  keep <- !is.na(matched_idx)
  if (sum(keep) < 3) {
    message("    Skipping: <3 DE-eligible samples")
    return(invisible(NULL))
  }
  cm <- counts_mat_full[, keep, drop = FALSE]
  meta <- sample_meta_de[matched_idx[keep]]

  # Drop samples with <50 cells in this cell type's pseudobulk.
  # Cell counts here use sample-level n_cells across ALL cell types as a proxy
  # (matches the existing pipeline's logic in pseudobulk_de.R).
  keep_mc <- meta$n_cells >= 50L
  cm <- cm[, keep_mc, drop = FALSE]
  meta <- meta[keep_mc]
  if (ncol(cm) < 3) {
    message("    Skipping: <3 samples after min-cells filter")
    return(invisible(NULL))
  }

  # Restrict to two stages of this contrast
  keep_stage <- meta$disease_stage_coarse %in% c(ref_stage, trt_stage)
  cm <- cm[, keep_stage, drop = FALSE]
  meta <- meta[keep_stage]

  stage_tab <- table(meta$disease_stage_coarse)
  if (length(stage_tab) < 2 || any(stage_tab < 2)) {
    message("    Skipping: need >=2 samples per stage; got ",
            paste(names(stage_tab), stage_tab, sep = "=", collapse = ", "))
    return(invisible(NULL))
  }
  message(sprintf("    Samples: %d (%s)", ncol(cm),
                  paste(names(stage_tab), stage_tab, sep = "=", collapse = ", ")))

  # Filter low-expression genes
  keep_genes <- rowSums(cm >= 5) >= 3
  cm <- cm[keep_genes, , drop = FALSE]
  if (nrow(cm) < 100) {
    message("    Skipping: too few genes (", nrow(cm), ")")
    return(invisible(NULL))
  }
  message(sprintf("    Genes after filter: %d", nrow(cm)))

  stage <- factor(meta$disease_stage_coarse, levels = c(ref_stage, trt_stage))
  dataset <- factor(meta$dataset)

  if (length(unique(meta$dataset)) > 1) {
    design <- model.matrix(~ dataset + stage)
  } else {
    design <- model.matrix(~ stage)
  }
  coef_name <- paste0("stage", trt_stage)
  if (!coef_name %in% colnames(design)) {
    message("    WARNING: coefficient ", coef_name,
            " not found. Available: ",
            paste(colnames(design), collapse = ", "))
    return(invisible(NULL))
  }

  res <- tryCatch({
    dge <- DGEList(counts = cm)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    coef_idx <- which(colnames(design) == coef_name)
    tt <- topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    setDT(tt)
    setnames(tt,
             old = c("logFC", "adj.P.Val", "P.Value", "t", "B"),
             new = c("logFC", "padj", "pvalue", "t_stat", "B"),
             skip_absent = TRUE)
    tt[, cell_type := ct_name]
    tt[, contrast := contrast_name]
    tt
  }, error = function(e) {
    message("    ERROR in limma-voom: ", e$message)
    NULL
  })

  if (is.null(res)) return(invisible(NULL))

  out_file <- file.path(out_dir, paste0(ct_name, "_", contrast_name, "_de.csv"))
  fwrite(res, out_file)
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  message(sprintf("    Saved %s (%d genes, %d sig padj<0.05)",
                  basename(out_file), nrow(res), n_sig))
  invisible(res)
}

# ---------------------------------------------------------------------------
# 3. Iterate cell types x contrasts
# ---------------------------------------------------------------------------
for (ct_name in TARGET_CELL_TYPES) {
  pb_file <- file.path(pb_dir, paste0(ct_name, "_pseudobulk.csv"))
  if (!file.exists(pb_file)) {
    message(sprintf("\n=== %s: pseudobulk file not found, skipping ===", ct_name))
    next
  }
  message(sprintf("\n=== %s ===", ct_name))
  counts_dt <- fread(pb_file)
  gene_names <- counts_dt[[1]]
  counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
  rownames(counts_mat) <- gene_names
  storage.mode(counts_mat) <- "numeric"
  n_runs_cols <- ncol(counts_mat)
  # Collapse run-level columns (SRR) to donor level (sum raw counts per donor).
  counts_mat <- collapse_counts_to_donor(counts_mat, BASE)
  message(sprintf("  Loaded %d genes x %d samples (donor-collapsed from %d runs)",
                  nrow(counts_mat), ncol(counts_mat), n_runs_cols))

  for (contrast_name in names(CONTRASTS)) {
    pair <- CONTRASTS[[contrast_name]]
    run_de(counts_mat, ct_name, contrast_name,
           ref_stage = pair[1], trt_stage = pair[2])
  }
}

message("\nAll pseudobulk DE coarse-stage contrasts complete.")
message("Outputs in: ", out_dir)
