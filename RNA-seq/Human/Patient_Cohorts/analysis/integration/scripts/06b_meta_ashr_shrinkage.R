#!/usr/bin/env Rscript
# 06b_meta_ashr_shrinkage.R
# ---------------------------------------------------------------------------
# Apply adaptive shrinkage (ashr) to the limma-voom + metafor REML meta-analysis
# effect sizes, mirroring 05b_ashr_shrinkage.R (which shrinks the dream output).
#
# This produces shrunk logFC + local false sign rate (lfsr) on the metafor
# random-effects pooled estimates, so the Cas13 library can consume metafor as
# its canonical human DEG source with the SAME lfsr<0.05 & shrunk_logFC>floor
# thresholding machinery it currently uses on dream_results_ashr.csv.
#
# Processes TWO metafor tables (same recipe, slightly different schemas):
#   A) Disease-vs-Control  : results/integration/meta_analysis_results.csv
#                            (gene version-stripped; SE col = meta_SE)
#        -> results/integration/meta_results_ashr.csv
#   B) MASH-vs-MASL (NASH>NAFL): results/disease_signatures/nafl_vs_nash_meta.csv
#                            (gene versioned; SE col = meta_se)
#        -> results/disease_signatures/nafl_vs_nash_meta_ashr.csv
#
# Output schema (both): gene (version-stripped), logFC (=meta_logFC),
#   shrunk_logFC, shrunk_se, lfsr, svalue, symbol, gene_biotype,
#   meta_padj, meta_I2, n_datasets  -- the first 6 + symbol match the library's
#   select() contract in rebuild_cas13_library.R.
#
# ashr recipe is identical to 05b: ash(betahat, sebetahat,
#   mixcompdist="halfuniform", method="shrink").
#
# Runs AFTER 06_meta_analysis.R and 13_nafl_vs_nash_de.R.
# SLURM: cpu, 4 CPU, 32G, 2h, --job-name=ashr
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ashr)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
INT  <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
DSIG <- file.path(INT, "results/disease_signatures")
GMETA <- file.path(PROJECT_ROOT, "data/gencode_v49_gene_metadata.tsv.gz")

strip_v <- function(x) sub("[.][0-9]+$", "", x)

# ── Gene -> symbol/biotype map (dream-independent; canonical GENCODE v49) ──────
gmeta <- fread(GMETA, select = c("ensembl_base", "gene_name", "gene_biotype"))
setnames(gmeta, c("gene_name", "gene_biotype"), c("symbol", "gene_biotype"))
gmeta <- unique(gmeta, by = "ensembl_base")

# ── Reusable ashr-on-metafor routine ──────────────────────────────────────────
# in_path  : metafor results CSV
# out_path : where to write the ashr-shrunk table
# se_col   : name of the SE column ("meta_SE" or "meta_se")
# i2_col   : name of the I^2 column ("meta_I2" or "I2")
# label    : pretty name for logging
shrink_meta <- function(in_path, out_path, se_col, i2_col, label) {
  cat("\n============================================================\n")
  cat("06b: ashr shrinkage of", label, "\n")
  cat("============================================================\n")
  m <- fread(in_path)
  cat("Loaded:", nrow(m), "genes from", basename(in_path), "\n")
  cat("Columns:", paste(names(m), collapse = ", "), "\n")

  stopifnot("meta_logFC" %in% names(m), se_col %in% names(m))
  m[, gene := strip_v(gene)]
  m[, se := get(se_col)]

  n_bad <- sum(is.na(m$se) | !is.finite(m$se) | m$se <= 0)
  if (n_bad > 0)
    cat("WARNING:", n_bad, "genes with invalid SE — excluded from shrinkage\n")
  valid <- !is.na(m$se) & is.finite(m$se) & m$se > 0

  cat("Running ashr (mixcompdist='halfuniform') on", sum(valid), "genes...\n")
  ash_fit <- ash(
    betahat     = m$meta_logFC[valid],
    sebetahat   = m$se[valid],
    mixcompdist = "halfuniform",
    method      = "shrink"
  )

  m[, `:=`(shrunk_logFC = NA_real_, shrunk_se = NA_real_,
           lfsr = NA_real_, svalue = NA_real_)]
  m[valid, shrunk_logFC := ash_fit$result$PosteriorMean]
  m[valid, shrunk_se    := ash_fit$result$PosteriorSD]
  m[valid, lfsr         := ash_fit$result$lfsr]
  m[valid, svalue       := ash_fit$result$svalue]

  # logFC = raw metafor pooled estimate (library reads a column named `logFC`)
  m[, logFC := meta_logFC]
  # dream_results_ashr.csv-compatible aliases so figure scripts that read
  # `padj` / `P.Value` / `se` are drop-in when their source path is swapped to this file.
  m[, padj    := meta_padj]
  m[, P.Value := meta_pval]
  m[, se      := get(se_col)]

  # Attach symbol + biotype (GENCODE v49, version-stripped join)
  m <- merge(m, gmeta, by.x = "gene", by.y = "ensembl_base", all.x = TRUE)

  # ── Diagnostics ─────────────────────────────────────────────────────────────
  cat("\n--- Shrinkage summary ---\n")
  cat("  lfsr < 0.05:", sum(m$lfsr < 0.05, na.rm = TRUE), "\n")
  cat("  lfsr < 0.01:", sum(m$lfsr < 0.01, na.rm = TRUE), "\n")
  sig <- m[lfsr < 0.05]
  q <- quantile(abs(sig$shrunk_logFC),
                probs = c(0, .05, .10, .25, .50, .75, .90, 1), na.rm = TRUE)
  cat("  |shrunk_logFC| quantiles (lfsr<0.05, n=", nrow(sig), "):\n", sep = "")
  print(round(q, 4))
  rho <- cor(m$meta_logFC[valid], m$shrunk_logFC[valid],
             method = "spearman", use = "complete.obs")
  cat(sprintf("  Spearman rho (raw vs shrunk): %.4f\n", rho))
  for (g in c("COL1A1", "COL1A2", "THRB", "PNPLA3")) {
    r <- m[symbol == g]
    if (nrow(r))
      cat(sprintf("  %-7s logFC=%+.3f shrunk=%+.3f lfsr=%.2g\n",
                  g, r$logFC[1], r$shrunk_logFC[1], r$lfsr[1]))
  }
  # ashr at floors of interest (matches library's HUMAN_SHRUNK_THR)
  for (th in c(0.2, 0.3, 0.5)) {
    n_up <- sum(m$lfsr < 0.05 & m$shrunk_logFC >  th, na.rm = TRUE)
    n_dn <- sum(m$lfsr < 0.05 & m$shrunk_logFC < -th, na.rm = TRUE)
    cat(sprintf("  lfsr<0.05 & |shrunk|>%.1f : up=%d down=%d\n", th, n_up, n_dn))
  }

  # ── Write canonical columns (library-compatible) + passthrough ───────────────
  keep <- c("gene", "logFC", "shrunk_logFC", "shrunk_se", "lfsr", "svalue",
            "symbol", "gene_biotype", "padj", "P.Value", "se",
            "meta_logFC", se_col, "meta_pval", i2_col, "meta_padj", "n_datasets")
  keep <- intersect(keep, names(m))
  setcolorder(m, keep)
  fwrite(m[, ..keep], out_path)
  cat("\nSaved:", out_path, "(", nrow(m), "genes x", length(keep), "cols)\n")
  invisible(m)
}

# ── A. Disease-vs-Control ─────────────────────────────────────────────────────
shrink_meta(
  in_path  = file.path(RDIR, "meta_analysis_results.csv"),
  out_path = file.path(RDIR, "meta_results_ashr.csv"),
  se_col   = "meta_SE",
  i2_col   = "meta_I2",
  label    = "Disease-vs-Control metafor (REML)"
)

# ── B. MASH-vs-MASL (NASH vs NAFL) ────────────────────────────────────────────
shrink_meta(
  in_path  = file.path(DSIG, "nafl_vs_nash_meta.csv"),
  out_path = file.path(DSIG, "nafl_vs_nash_meta_ashr.csv"),
  se_col   = "meta_se",
  i2_col   = "I2",
  label    = "MASH-vs-MASL metafor (REML)"
)

cat("\nScript 06b complete.\n")
