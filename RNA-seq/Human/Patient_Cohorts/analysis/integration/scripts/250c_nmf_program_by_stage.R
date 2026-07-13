# 250c_nmf_program_by_stage.R
# Patch — compute NMF k=6 program × stage matrix that 250 silently failed on.
# Fix: nmf$nmf_results[["6"]] is the NMFfitX1 directly (no $res wrapper).

suppressPackageStartupMessages({
  library(NMF)
})

PROJECT_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
NMF_CACHE <- file.path(PROJECT_ROOT, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")

cat("== NMF program × stage projection ==\n")

nmf <- readRDS(NMF_CACHE)
k6 <- nmf$nmf_results[["6"]]
H <- NMF::coef(k6)  # 6 programs × n samples
cat(sprintf("H matrix: %d programs × %d samples\n", nrow(H), ncol(H)))

meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))

common <- intersect(colnames(H), meta$sample_id)
cat(sprintf("Samples in H ∩ meta: %d\n", length(common)))

H_match <- H[, common, drop = FALSE]
meta_match <- meta[match(common, meta$sample_id), ]

# Exclude coarse-staged GSE213621 (study reports grouped F0F1/F2/F3F4, mapped F0F1->1,
# F3F4->3 — NOT true Kleiner F0-F4) + dropped PRJNA512027, so the per-stage bins below
# are true individual stages (mirrors the 14b per-stage DE allowlist). 2026-06-24.
keep_clean <- !meta_match$dataset %in% c("GSE213621", "PRJNA512027")
cat(sprintf("Excluding coarse/dropped cohorts: %d -> %d samples\n",
            length(common), sum(keep_clean)))
H_match    <- H_match[, keep_clean, drop = FALSE]
meta_match <- meta_match[keep_clean, ]

prog_names <- paste0("P", seq_len(nrow(H_match)))
rows <- list()
for (s in 0:4) {
  s_idx <- which(meta_match$fibrosis_stage == s)
  if (length(s_idx) < 5) next
  for (p in seq_len(nrow(H_match))) {
    rows[[paste(s, p)]] <- data.frame(
      stage = s, program = prog_names[p],
      n_samples = length(s_idx),
      mean_score = mean(H_match[p, s_idx]),
      median_score = median(H_match[p, s_idx]),
      sd_score = sd(H_match[p, s_idx]),
      stringsAsFactors = FALSE
    )
  }
}
nmf_per_stage <- do.call(rbind, rows)
write.csv(nmf_per_stage, file.path(OUT_DIR, "two_transition_nmf_program_by_stage.csv"),
          row.names = FALSE)

# Wide table
wide <- reshape(nmf_per_stage[, c("stage", "program", "mean_score")],
                idvar = "program", timevar = "stage", direction = "wide")
print(wide)

# Per-transition fold change (mean ratio)
cat("\nPer-transition fold change (mean(B) / mean(A)):\n")
for (transition in list(c(0, 1), c(1, 2), c(2, 3), c(3, 4))) {
  cat(sprintf("\nF%d → F%d:\n", transition[1], transition[2]))
  for (p in seq_len(nrow(H_match))) {
    a <- H_match[p, meta_match$fibrosis_stage == transition[1]]
    b <- H_match[p, meta_match$fibrosis_stage == transition[2]]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 5 || length(b) < 5) next
    fc <- mean(b) / (mean(a) + 1e-9)
    p_val <- wilcox.test(a, b)$p.value
    cat(sprintf("  P%d: mean %s → %s (FC = %.2fx, p = %.2e)\n",
                p, formatC(mean(a), digits = 3),
                formatC(mean(b), digits = 3), fc, p_val))
  }
}

cat("\n== NMF program × stage complete. ==\n")
