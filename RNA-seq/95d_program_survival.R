#!/usr/bin/env Rscript
# 95d_program_survival.R — cross-cache program correspondence analysis.
#
# Bravo Task 9 second-pass reporting: which programs from the canonical
# clean-cache k=6 survive / merge / split under the PROT_ONLY=TRUE and
# PROT_ONLY=FALSE (extended strip) remediation refits?
#
# Method: Hungarian matching on gene-wise W-column cosine similarity over the
# intersecting gene set. Labels come from 95c-emitted program_labels*.csv.
#
# Writes: RNA-seq/results/subtypes/program_survival_canonical_vs_remediation.tsv

suppressPackageStartupMessages({
  library(data.table); library(NMF); library(clue)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir <- file.path(BASE, "RNA-seq/results/subtypes")
K <- as.integer(Sys.getenv("NMF_CHOSEN_K", "6"))

caches <- list(
  canonical    = file.path(out_dir, "nmf_results_cache_clean.rds"),
  protonly     = file.path(out_dir, "nmf_results_cache_clean_protonly.rds"),
  nonprotonly  = file.path(out_dir, "nmf_results_cache_clean_nonprotonly.rds")
)
labels_files <- list(
  canonical    = file.path(out_dir, "program_labels.csv"),
  protonly     = file.path(out_dir, "program_labels_protonly.csv"),
  nonprotonly  = file.path(out_dir, "program_labels_nonprotonly.csv")
)

load_W <- function(path, k) {
  obj <- readRDS(path)
  fit <- obj$nmf_results[[as.character(k)]]
  W <- basis(fit)
  rownames(W) <- rownames(obj$mat_nn)
  W
}
load_labels <- function(path) {
  # Skip header comment lines (prefixed with '#')
  lines <- readLines(path); lines <- lines[!grepl("^#", lines)]
  fread(text = paste(lines, collapse = "\n"))
}

W_can <- load_W(caches$canonical,   K)
W_pot <- load_W(caches$protonly,    K)
W_non <- load_W(caches$nonprotonly, K)

lab_can <- load_labels(labels_files$canonical)
lab_pot <- load_labels(labels_files$protonly)
lab_non <- load_labels(labels_files$nonprotonly)

# Hungarian match on cosine over intersecting genes
match_pair <- function(W1, W2) {
  common <- intersect(rownames(W1), rownames(W2))
  W1 <- W1[common, ]; W2 <- W2[common, ]
  n1 <- sqrt(colSums(W1^2)); n2 <- sqrt(colSums(W2^2))
  cos_mat <- (t(W1) %*% W2) / (outer(n1, n2))
  rownames(cos_mat) <- paste0("P", seq_len(ncol(W1)))
  colnames(cos_mat) <- paste0("P", seq_len(ncol(W2)))
  assign <- as.integer(clue::solve_LSAP(1 - cos_mat, maximum = FALSE))
  data.table(
    canonical_program = paste0("P", seq_len(nrow(cos_mat))),
    matched_program = paste0("P", assign),
    matched_cosine = vapply(seq_len(nrow(cos_mat)),
                             function(i) cos_mat[i, assign[i]], numeric(1)),
    second_best_cosine = vapply(seq_len(nrow(cos_mat)), function(i) {
      row <- cos_mat[i, ]; row[assign[i]] <- -Inf; max(row)
    }, numeric(1)),
    n_common_genes = length(common)
  )
}

m_pot <- match_pair(W_can, W_pot)
m_non <- match_pair(W_can, W_non)

# Attach labels on both sides
add_labels <- function(match_dt, lab_src, lab_dst, variant) {
  out <- copy(match_dt)
  out[, variant := variant]
  out[, canonical_label := lab_src$biological_label[match(canonical_program, lab_src$program_code)]]
  out[, matched_label   := lab_dst$biological_label[match(matched_program,   lab_dst$program_code)]]
  out[, label_preserved := canonical_label == matched_label]
  out[, margin := round(matched_cosine - second_best_cosine, 3)]
  out[, matched_cosine := round(matched_cosine, 3)]
  out[, second_best_cosine := round(second_best_cosine, 3)]
  setcolorder(out, c("variant","canonical_program","canonical_label",
                     "matched_program","matched_label","label_preserved",
                     "matched_cosine","second_best_cosine","margin","n_common_genes"))
  out
}

res <- rbind(
  add_labels(m_pot, lab_can, lab_pot, "protonly"),
  add_labels(m_non, lab_can, lab_non, "nonprotonly")
)

out_tsv <- file.path(out_dir, "program_survival_canonical_vs_remediation.tsv")
fwrite(res, out_tsv, sep = "\t")
cat("Wrote:", out_tsv, "\n\n")
print(res)

# Summary
cat("\n=== Summary ===\n")
for (v in c("protonly","nonprotonly")) {
  sub <- res[variant == v]
  n_pres <- sum(sub$label_preserved, na.rm = TRUE)
  n_high <- sum(sub$matched_cosine >= 0.7)
  n_mod  <- sum(sub$matched_cosine >= 0.4 & sub$matched_cosine < 0.7)
  n_low  <- sum(sub$matched_cosine < 0.4)
  cat(sprintf("  %s: %d/%d labels preserved; cosine>=0.7: %d, 0.4-0.7: %d, <0.4: %d\n",
              v, n_pres, nrow(sub), n_high, n_mod, n_low))
}
