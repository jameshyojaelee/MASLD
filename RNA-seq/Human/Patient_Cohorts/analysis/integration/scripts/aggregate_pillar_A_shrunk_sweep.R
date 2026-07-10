#!/usr/bin/env Rscript
# aggregate_pillar_A_shrunk_sweep.R
# ---------------------------------------------------------------------------
# Build the Pillar-A LFC-threshold sweep summaries (the two CSVs consumed by
# scripts/figures/figS_robust_lfc_sweep.R) at BOTH scales:
#   raw    : gate padj<0.05 & |logFC|>cut          (suffix _raw)
#   shrunk : gate lfsr<0.05 & |shrunk_logFC|>cut    (suffix _shrunk)   [PRIMARY]
#
# ashr is applied PER RESAMPLING ITERATION on (logFC, SE) where SE = logFC / t
# (dream moderated t). This needs NO refit — the per-iter dream fits are reused
# as-is. The expensive CPSS/Bootstrap/K-fold dream arrays are NOT re-run.
#
# Inputs (per-iter, complete STAR-backup set; gene,logFC,t,padj[,K,fold]):
#   robustness/<cpss|bootstrap|kfold>_dir/*.csv     (env-overridable)
#   dream_results.csv                                (full ref: logFC,padj,SE)
# Outputs (RNA-seq/results/audit_sensitivity/):
#   pillar_A_lfc_threshold_sweep_{raw,shrunk}.csv
#   pillar_A_lfc_sweep_iter_long_{raw,shrunk}.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ashr); library(parallel)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Per-iter source dirs — default to the complete STAR-backup set (the live
# robustness/{cpss,bootstrap,kfold} dirs hold an interrupted partial re-run).
CPSS_DIR  <- Sys.getenv("CPSS_DIR",  file.path(RDIR, "robustness/cpss_star_backup"))
BOOT_DIR  <- Sys.getenv("BOOT_DIR",  file.path(RDIR, "robustness/bootstrap_star_backup"))
KFOLD_DIR <- Sys.getenv("KFOLD_DIR", file.path(RDIR, "robustness/kfold_star_backup"))
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

PADJ_THR <- 0.05
LFSR_THR <- 0.05
THRS <- c(0, 0.1, 0.3, 0.5, 0.75, 1, 2)
TAGS <- c("nolfc", "lfc1", "lfc3", "lfc5", "lfc075", "lfc100", "lfc200")
LABS <- c("|LFC|>0\n(padj<0.05 only)", "|LFC|>0.1", "|LFC|>0.3",
          "|LFC|>0.5\n(Tier 1)", "|LFC|>0.75", "|LFC|>1.0", "|LFC|>2.0")

# --- Full reference + ashr on the full fit ---------------------------------
full <- fread(file.path(RDIR, "dream_results.csv"))
full <- full[is.finite(logFC) & is.finite(padj)]
se_full <- full$SE
se_full[!is.finite(se_full) | se_full <= 0] <- NA_real_
ok <- !is.na(se_full)
full[, shrunk_logFC := 0.0][, lfsr := 1.0]
a_full <- ash(full$logFC[ok], se_full[ok], mixcompdist = "normal")
full$shrunk_logFC[ok] <- a_full$result$PosteriorMean
full$lfsr[ok]         <- a_full$result$lfsr
cat(sprintf("Full ref: %d genes (%d ashr-eligible)\n", nrow(full), sum(ok)))

# --- Per-iter loader: read, derive SE, ashr-shrink once --------------------
load_iter <- function(f, with_kfold = FALSE) {
  cols <- c("gene", "logFC", "t", "padj")
  if (with_kfold) cols <- c(cols, "K", "fold")
  dt <- fread(f, select = cols)
  dt <- dt[is.finite(logFC) & is.finite(padj)]
  se <- dt$logFC / dt$t
  se[!is.finite(se)] <- NA_real_
  se <- abs(se)
  se[se <= 0] <- NA_real_
  dt[, shrunk_logFC := 0.0][, lfsr := 1.0]
  ok <- !is.na(se)
  if (sum(ok) > 1) {
    a <- ash(dt$logFC[ok], se[ok], mixcompdist = "normal")
    dt$shrunk_logFC[ok] <- a$result$PosteriorMean
    dt$lfsr[ok]         <- a$result$lfsr
  }
  dt
}

cpss_files  <- list.files(CPSS_DIR,  pattern = "^iter_\\d{3}_half[AB]\\.csv$", full.names = TRUE)
boot_files  <- list.files(BOOT_DIR,  pattern = "^iter_\\d{4}\\.csv$",          full.names = TRUE)
kfold_files <- list.files(KFOLD_DIR, pattern = "^iter_\\d{2}\\.csv$",          full.names = TRUE)
cat(sprintf("Per-iter files: cpss=%d boot=%d kfold=%d  (ncores=%d)\n",
            length(cpss_files), length(boot_files), length(kfold_files), NCORES))

cpss_dt  <- mclapply(cpss_files,  load_iter, mc.cores = NCORES)
boot_dt  <- mclapply(boot_files,  load_iter, mc.cores = NCORES)
kfold_dt <- mclapply(kfold_files, load_iter, with_kfold = TRUE, mc.cores = NCORES)
names(cpss_dt) <- basename(cpss_files)
cpss_pair <- sub("^iter_(\\d{3}).*$", "\\1", basename(cpss_files))
all_genes <- full$gene

# --- selection predicate for a scale ---------------------------------------
sel_genes <- function(dt, scale, cut) {
  if (scale == "raw") dt[padj < PADJ_THR & abs(logFC)        > cut, gene]
  else                dt[lfsr < LFSR_THR & abs(shrunk_logFC) > cut, gene]
}
full_deg <- function(scale, cut) {
  if (scale == "raw") full[padj < PADJ_THR & abs(logFC)        > cut, gene]
  else                full[lfsr < LFSR_THR & abs(shrunk_logFC) > cut, gene]
}
eff_col <- function(scale) if (scale == "raw") "logFC" else "shrunk_logFC"

# --- one (scale, threshold) row + iter_long rows ---------------------------
compute_one <- function(scale, cut, tag, lab) {
  deg <- full_deg(scale, cut)
  full_eff <- full[, .(gene, fe = get(eff_col(scale)))]

  # CPSS: per pair union over the two halves, then pi_hat per gene
  pairs <- unique(cpss_pair)
  pmat <- matrix(0L, nrow = length(all_genes), ncol = length(pairs),
                 dimnames = list(all_genes, NULL))
  for (j in seq_along(pairs)) {
    idx <- which(cpss_pair == pairs[j])
    sg  <- unique(unlist(lapply(cpss_dt[idx], sel_genes, scale = scale, cut = cut)))
    pmat[match(sg, all_genes), j] <- 1L
  }
  pi_hat <- rowMeans(pmat); names(pi_hat) <- all_genes
  pi_deg <- pi_hat[deg]; pi_non <- pi_hat[setdiff(all_genes, deg)]

  # Bootstrap: per-gene selection frequency
  bsel <- lapply(boot_dt, sel_genes, scale = scale, cut = cut)
  bcount <- tabulate(match(unlist(bsel), all_genes), nbins = length(all_genes))
  bfreq <- bcount / length(boot_dt); names(bfreq) <- all_genes
  bfreq_deg <- bfreq[deg]

  # K-fold recurrence (K==10)
  k10 <- which(vapply(kfold_dt, function(d) d$K[1] == 10, logical(1)))
  k50 <- which(vapply(kfold_dt, function(d) d$K[1] == 50, logical(1)))
  krec <- function(idx) {
    sg <- unlist(lapply(kfold_dt[idx], sel_genes, scale = scale, cut = cut))
    r <- tabulate(match(sg, all_genes), nbins = length(all_genes)); names(r) <- all_genes; r
  }
  r10 <- krec(k10); r50 <- krec(k50)

  sweep_row <- data.table(
    threshold = cut, tag = tag, n_DEG = length(deg),
    mean_pi_hat   = round(mean(pi_deg, na.rm = TRUE), 3),
    pi_hat_ge_07  = sum(pi_deg >= 0.7, na.rm = TRUE),
    pct_pi_ge_07  = round(100 * mean(pi_deg >= 0.7, na.rm = TRUE), 1),
    mean_boot     = round(mean(bfreq_deg, na.rm = TRUE), 3),
    boot_ge_09    = sum(bfreq_deg >= 0.9, na.rm = TRUE),
    pct_boot_ge_09 = round(100 * mean(bfreq_deg >= 0.9, na.rm = TRUE), 1),
    mean_kfold10  = round(mean(r10[deg], na.rm = TRUE), 2),
    mean_kfold50  = round(mean(r50[deg], na.rm = TRUE), 2),
    mean_pi_nonDEG = round(mean(pi_non, na.rm = TRUE), 3),
    label = lab)

  # iter_long: rho vs full effect + jaccard + recovery, per iteration
  iter_block <- function(dts, src, kvec = NULL) {
    rbindlist(lapply(seq_along(dts), function(i) {
      d <- dts[[i]]
      sg <- sel_genes(d, scale, cut)
      m <- merge(full_eff, d[, .(gene, ie = get(eff_col(scale)))], by = "gene")
      rho <- suppressWarnings(cor(m$fe, m$ie, method = "spearman"))
      inter <- length(intersect(deg, sg)); un <- length(union(deg, sg))
      s <- if (!is.null(kvec)) paste0("kfold_K", kvec[i]) else src
      data.table(threshold = cut, source = s, rho = round(rho, 4),
                 jaccard = round(if (un > 0) inter / un else 0, 4),
                 recovery_pct = round(if (length(deg) > 0) 100 * inter / length(deg) else 0, 2),
                 n_sel = length(sg))
    }))
  }
  iter <- rbind(
    iter_block(cpss_dt, "cpss"),
    iter_block(boot_dt, "bootstrap"),
    iter_block(kfold_dt, NULL, kvec = vapply(kfold_dt, function(d) d$K[1], integer(1))))
  list(sweep = sweep_row, iter = iter)
}

for (scale in c("raw", "shrunk")) {
  cat(sprintf("\n=== scale = %s ===\n", scale))
  res <- lapply(seq_along(THRS),
                function(i) compute_one(scale, THRS[i], TAGS[i], LABS[i]))
  sweep <- rbindlist(lapply(res, `[[`, "sweep"))
  iter  <- rbindlist(lapply(res, `[[`, "iter"))
  fwrite(sweep, file.path(OUT, sprintf("pillar_A_lfc_threshold_sweep_%s.csv", scale)))
  fwrite(iter,  file.path(OUT, sprintf("pillar_A_lfc_sweep_iter_long_%s.csv", scale)))
  cat("Wrote sweep + iter_long for", scale, "\n"); print(sweep)
}
cat("\nDone. (raw recompute is for provenance verification vs the canonical raw table.)\n")
