#!/usr/bin/env Rscript
# 310c: Cross-modal — Plasma proteomics (Olink ssGSEA)
#
# For each hepatocyte subtype, intersect top 50 markers with Olink proteins,
# compute ssGSEA scores per subject, and test association with fibrosis stage.
#
# Input:
#   hepatocyte_subtypes/subtype_markers.csv  (from 309)
#   Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt
#   Analysis/Proteomics/results/gse276114_disease_metadata.csv
#
# Output (to hepatocyte_subtypes/crossmodal/plasma/):
#   subtype_plasma_scores.csv
#   subtype_plasma_fibrosis_association.csv
#   subtype_plasma_overlap.csv
#
# Environment: rnaseq (CPU)

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SUB_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
OLINK_DIR <- file.path(BASE, "Analysis/Proteomics/data/olink_plasma")
PROT_DIR  <- file.path(BASE, "Analysis/Proteomics/results")

OUT_DIR <- file.path(SUB_DIR, "crossmodal", "plasma")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Olink protein matrix ────────────────────────────────────────
message("=== Loading Olink plasma data ===")
olink_file <- file.path(OLINK_DIR, "olink.qc.finished.mendeley.data.txt")
stopifnot(file.exists(olink_file))

olink_raw <- fread(olink_file, header = TRUE, sep = "\t")
olink_proteins <- olink_raw$Assay
olink_mat <- as.matrix(olink_raw[, -1, with = FALSE])
rownames(olink_mat) <- olink_proteins

# NPX matrix: subjects x proteins
npx <- t(olink_mat)
message("  NPX matrix: ", nrow(npx), " subjects x ", ncol(npx), " proteins")

# ── 2. Load fibrosis metadata and match ─────────────────────────────────
message("=== Loading fibrosis metadata ===")
meta_file <- file.path(PROT_DIR, "gse276114_disease_metadata.csv")
meta <- fread(meta_file)
message("  Metadata: ", nrow(meta), " samples")

# Map Olink subjects to metadata
subject_nums <- as.integer(sub("^Subject\\s+", "", rownames(npx)))
label_dt <- data.table(
  subject_idx  = seq_len(nrow(npx)),
  subject_num  = subject_nums,
  subject_name = rownames(npx)
)
label_dt <- merge(label_dt, meta[, .(sample_number, disease_group, disease)],
                  by.x = "subject_num", by.y = "sample_number", all.x = TRUE)

label_dt_matched <- label_dt[!is.na(disease_group)]
message("  Subjects with fibrosis metadata: ", nrow(label_dt_matched))

npx_matched <- npx[label_dt_matched$subject_idx, , drop = FALSE]
rownames(npx_matched) <- label_dt_matched$subject_name

# Extract numeric fibrosis stage
# disease_group values: "F0-2", "F3", "F4" — map explicitly (gsub produces "0-2" → NA)
label_dt_matched[, fibrosis_numeric := fcase(
  disease_group == "F0-2", 1,
  disease_group == "F3",   3,
  disease_group == "F4",   4,
  default = NA_real_
)]

# ── 3. Load subtype markers ─────────────────────────────────────────────
message("=== Loading subtype markers ===")
markers <- fread(file.path(SUB_DIR, "subtype_markers.csv"))
subtypes <- as.character(sort(unique(markers$subtype)))
message("  Subtypes: ", paste(subtypes, collapse = ", "))

olink_protein_names <- colnames(npx_matched)

# Build gene sets (intersected with Olink proteins)
gene_sets <- list()
overlap_rows <- list()

for (st in subtypes) {
  st_genes <- markers[subtype == st, names]
  if (length(st_genes) > 50) st_genes <- st_genes[1:50]

  overlap <- st_genes[st_genes %in% olink_protein_names]
  gene_sets[[paste0("Hep_", st)]] <- overlap

  overlap_rows[[st]] <- data.table(
    subtype = st,
    n_markers = length(st_genes),
    n_in_olink = length(overlap),
    overlap_genes = paste(overlap, collapse = ";")
  )
  message("  ", st, ": ", length(overlap), "/", length(st_genes), " markers in Olink")
}

overlap_dt <- rbindlist(overlap_rows)
fwrite(overlap_dt, file.path(OUT_DIR, "subtype_plasma_overlap.csv"))

# Filter to gene sets with >= 5 genes
gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]
message("Gene sets with >=5 Olink proteins: ", length(gene_sets))

if (length(gene_sets) == 0) {
  message("WARNING: No subtype has >= 5 marker genes in Olink. Exiting.")
  quit(save = "no")
}

# ── 4. ssGSEA on Olink data ─────────────────────────────────────────────
message("=== Running ssGSEA ===")

# GSVA needs genes x samples matrix
npx_t <- t(npx_matched)  # proteins x subjects

# GSVA >= 1.50 (Bioconductor 3.18+) uses ssgseaParam objects
gsva_version <- packageVersion("GSVA")
message("  GSVA version: ", as.character(gsva_version))

if (gsva_version >= "1.50") {
  ssgsea_param <- ssgseaParam(exprData = npx_t, geneSets = gene_sets,
                               normalize = TRUE)
  gsva_res <- gsva(ssgsea_param, verbose = FALSE)
} else {
  gsva_res <- gsva(npx_t, gene_sets, method = "ssgsea", verbose = FALSE)
}
# gsva_res: gene_sets x subjects

# Build output
plasma_scores <- as.data.table(t(gsva_res))
plasma_scores[, subject := rownames(npx_matched)]
plasma_scores[, disease_group := label_dt_matched$disease_group]
plasma_scores[, fibrosis_numeric := label_dt_matched$fibrosis_numeric]

fwrite(plasma_scores, file.path(OUT_DIR, "subtype_plasma_scores.csv"))
message("Saved plasma scores: ", nrow(plasma_scores), " subjects x ", ncol(gsva_res), " subtypes")

# ── 5. Fibrosis association tests ────────────────────────────────────────
message("=== Testing fibrosis association ===")

assoc_rows <- list()
for (nm in rownames(gsva_res)) {
  scores <- gsva_res[nm, ]
  fib <- label_dt_matched$fibrosis_numeric

  # Spearman correlation with fibrosis stage
  rho_test <- cor.test(fib, scores, method = "spearman", exact = FALSE)

  # Kruskal-Wallis across stages
  stage_list <- split(scores, label_dt_matched$disease_group)
  stage_list <- stage_list[sapply(stage_list, length) >= 3]

  if (length(stage_list) >= 3) {
    kw <- kruskal.test(stage_list)
    kw_p <- kw$p.value
  } else {
    kw_p <- NA
  }

  assoc_rows[[nm]] <- data.table(
    subtype = nm,
    spearman_rho = rho_test$estimate,
    spearman_pval = rho_test$p.value,
    kw_pval = kw_p,
    n_subjects = length(scores)
  )
}

assoc_dt <- rbindlist(assoc_rows)
fwrite(assoc_dt, file.path(OUT_DIR, "subtype_plasma_fibrosis_association.csv"))
message("Plasma fibrosis associations:\n")
print(assoc_dt)

message("\n310c complete. Results in: ", OUT_DIR)
