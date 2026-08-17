#!/usr/bin/env Rscript
# 08_run_eqtl_susie.R
# Pre-compute SuSiE fine-mapping fits for Broadaway liver eQTLs (per chromosome).
#
# Invocation:
#   CHR=1 Rscript 08_run_eqtl_susie.R          # explicit chromosome
#   Rscript 08_run_eqtl_susie.R                  # reads SLURM_ARRAY_TASK_ID
#
# Output (per chromosome):
#   results/eqtl_susie/chr{N}/{ENSG}_susie.rds   -- full susie_rss object
#   results/eqtl_susie/chr{N}_summary.tsv        -- per-gene summary (appended incrementally)

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

# ---------------------------------------------------------------------------
# Paths and parameters
# ---------------------------------------------------------------------------

FM_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
EQTL_DIR  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl"
# Env-overridable so a regenerated set can be written side-by-side instead of
# overwriting the RDS that canonical COLOC currently reads.
OUTPUT_DIR <- Sys.getenv("EQTL_SUSIE_DIR", unset = file.path(FM_DIR, "results/eqtl_susie"))

EQTL_N        <- 1183L
SUSIE_L        <- 10L
# REFIT (2026-08-15): iteration cap made env-overridable. 204 eGenes stopped
# converging at max_iter=500 once the LD panel was pinned to polyfun; each one
# costs ~50 gene-GWAS pairs because its fit is reused across all 50 studies.
SUSIE_MAX_ITER <- as.integer(Sys.getenv('SUSIE_MAX_ITER', unset = '500'))
CIS_WINDOW_BP  <- 500000L   # ±500 kb around gene body
MIN_SNPS       <- 30L        # minimum variants after LD alignment
LD_REGULARIZE  <- 1e-3

# ---------------------------------------------------------------------------
# Determine chromosome
# ---------------------------------------------------------------------------

CHR <- as.integer(Sys.getenv("CHR", unset = ""))
if (is.na(CHR) || CHR == 0L) {
  CHR <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""))
}
if (is.na(CHR) || CHR < 1L || CHR > 22L) {
  stop("CHR must be set to an integer 1-22 via env var CHR or SLURM_ARRAY_TASK_ID")
}

cat("=== eQTL SuSiE fine-mapping | chr", CHR, "===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Source helper functions (get_ld_base_dir, find_LD_block, get_ld_per_locus)
# ---------------------------------------------------------------------------

source(file.path(FM_DIR, "src/finemapping_functions.R"))

# ---------------------------------------------------------------------------
# Load eQTL data for this chromosome
# ---------------------------------------------------------------------------

eqtl_file <- file.path(EQTL_DIR, paste0("chr", CHR, "_marginal_summary_results.tsv"))
if (!file.exists(eqtl_file)) {
  stop("eQTL file not found: ", eqtl_file)
}

cat("Loading eQTL data:", eqtl_file, "\n")
dt <- fread(eqtl_file)
# Expected columns: Entrez, Variant, CHR, POS, NEA, EA, EAF, Beta, SE, PVAL, N, Studies, GeneSymbol, ENSG, Gene_Biotype

required_cols <- c("CHR", "POS", "NEA", "EA", "Beta", "SE", "ENSG", "GeneSymbol")
missing_cols  <- setdiff(required_cols, colnames(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Drop rows with missing essentials
dt <- dt[!is.na(POS) & !is.na(Beta) & !is.na(SE) & SE > 0 & ENSG != ""]

cat("Loaded", nrow(dt), "variants across", dt[, uniqueN(ENSG)], "eGenes\n\n")

# ---------------------------------------------------------------------------
# Output directory for this chromosome
# ---------------------------------------------------------------------------

chr_out_dir <- file.path(OUTPUT_DIR, paste0("chr", CHR))
dir.create(chr_out_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(OUTPUT_DIR, paste0("chr", CHR, "_summary.tsv"))

# Write header only if file does not yet exist (supports resume)
if (!file.exists(summary_file)) {
  fwrite(
    data.table(
      gene      = character(),
      ensg      = character(),
      chr       = integer(),
      n_snps    = integer(),
      converged = logical(),
      n_cs      = integer(),
      max_pip   = numeric(),
      status    = character()
    ),
    file = summary_file,
    sep  = "\t"
  )
}

# ---------------------------------------------------------------------------
# Per-eGene SuSiE loop
# ---------------------------------------------------------------------------

ensg_list  <- dt[, unique(ENSG)]
# REFIT: when REFIT_GENES points at a file of ENSG ids, process ONLY those.
# Re-fitting 204 genes instead of 18,889 is the whole point of this copy.
.refit_f <- Sys.getenv('REFIT_GENES', unset = '')
if (nzchar(.refit_f) && file.exists(.refit_f)) {
  .want <- readLines(.refit_f)
  ensg_list <- ensg_list[ensg_list %in% .want]
  cat('*** REFIT MODE: restricted to', length(ensg_list), 'of', length(.want), 'requested genes ***\n')
  if (length(ensg_list) == 0L) { cat('nothing to do on this chromosome\n'); quit(save='no', status=0) }
}
n_genes    <- length(ensg_list)

# Parallel worker support: GENE_START_PCT (0-99) sets where in the list to begin.
# REVERSE_ORDER=1 reverses direction. Checkpoint/skip handles overlap between workers.
start_pct <- as.integer(Sys.getenv("GENE_START_PCT", unset = "0"))
if (Sys.getenv("REVERSE_ORDER", unset = "0") == "1") {
  cat("*** REVERSE ORDER from", start_pct, "% ***\n")
  ensg_list <- rev(ensg_list)
  # For reverse, start_pct is measured from the original end
  start_idx <- max(1L, as.integer(n_genes * (100 - start_pct) / 100))
} else if (start_pct > 0) {
  start_idx <- max(1L, as.integer(n_genes * start_pct / 100))
  cat("*** STARTING FROM gene", start_idx, "/", n_genes, "(", start_pct, "%) ***\n")
} else {
  start_idx <- 1L
}

n_ok       <- 0L
n_cs       <- 0L
n_skipped  <- 0L
n_failed   <- 0L

cat("Processing", n_genes, "eGenes on chr", CHR, "(start_idx =", start_idx, ")\n\n")

for (i in start_idx:n_genes) {
  ensg_id   <- ensg_list[i]
  gene_dt   <- dt[ENSG == ensg_id]
  gene_name <- gene_dt$GeneSymbol[1]

  rds_path  <- file.path(chr_out_dir, paste0(ensg_id, "_susie.rds"))

  # ---- (a) Checkpoint / resume ----
  if (file.exists(rds_path)) {
    cat("[", i, "/", n_genes, "]", gene_name, "(", ensg_id, ") — already done, skipping\n")
    n_skipped <- n_skipped + 1L
    next
  }

  cat("[", i, "/", n_genes, "]", gene_name, "(", ensg_id, ") — ", nrow(gene_dt), "variants\n")

  # ---- (b/c) Define cis-window ----
  START <- min(gene_dt$POS) - CIS_WINDOW_BP
  END   <- max(gene_dt$POS) + CIS_WINDOW_BP
  START <- max(START, 1L)

  # ---- (d) Build sumstats data.frame for get_ld_per_locus ----
  # get_ld_per_locus expects: chromosome, position, allele1, allele2, beta, standard_error
  ss <- data.frame(
    chromosome      = gene_dt$CHR,
    position        = gene_dt$POS,
    allele1         = gene_dt$NEA,   # non-effect allele = ref
    allele2         = gene_dt$EA,    # effect allele = alt
    beta            = gene_dt$Beta,
    standard_error  = gene_dt$SE,
    stringsAsFactors = FALSE
  )

  # ---- (e) LD extraction ----
  ld_result <- tryCatch(
    get_ld_per_locus(ss, LOCUS = ensg_id, CHR = CHR, START = START, END = END, ancestry = "EUR"),
    error = function(e) {
      cat("  LD extraction error:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(ld_result)) {
    cat("  Skipped (LD extraction failed)\n")
    fwrite(data.table(gene = gene_name, ensg = ensg_id, chr = CHR,
                      n_snps = 0L, converged = NA, n_cs = 0L, max_pip = NA_real_,
                      status = "ld_failed"),
           file = summary_file, sep = "\t", append = TRUE, col.names = FALSE)
    n_failed <- n_failed + 1L
    next
  }

  ss_filtered <- ld_result[[1]]
  R_mat       <- ld_result[[2]]
  if (!is.matrix(R_mat)) R_mat <- as.matrix(R_mat)
  rm(ld_result); gc(verbose = FALSE)

  # ---- (f) Minimum SNP check ----
  n_snps <- nrow(ss_filtered)
  if (n_snps < MIN_SNPS) {
    cat("  Skipped (only", n_snps, "variants after LD alignment, min =", MIN_SNPS, ")\n")
    fwrite(data.table(gene = gene_name, ensg = ensg_id, chr = CHR,
                      n_snps = n_snps, converged = NA, n_cs = 0L, max_pip = NA_real_,
                      status = "too_few_snps"),
           file = summary_file, sep = "\t", append = TRUE, col.names = FALSE)
    n_failed <- n_failed + 1L
    next
  }

  # ---- (g/h/i) Z-scores and SNP names ----
  # SNP names in names(z) are transferred to colnames(s$lbf_variable) by susie_rss.
  # coloc.susie() requires these colnames to match between the two datasets.
  z      <- ss_filtered$beta / ss_filtered$standard_error
  snp_ids <- paste0(CHR, ":", ss_filtered$position)
  names(z) <- snp_ids

  # ---- (j) LD regularization (in-place to avoid extra copy) ----
  diag(R_mat) <- diag(R_mat) + LD_REGULARIZE
  colnames(R_mat) <- rownames(R_mat) <- snp_ids

  # ---- (k) susie_rss ----
  # estimate_residual_variance=FALSE + residual_variance=1 is the correct setting
  # for out-of-sample LD reference panels (UKBB LD, not in-sample LD).
  # See: coloc vignette "SuSiE with external LD" and 35s_susie_coloc_broadaway.R lines 716-733.
  s <- tryCatch(
    susie_rss(
      z                          = z,
      R                          = R_mat,
      n                          = EQTL_N,
      L                          = SUSIE_L,
      estimate_residual_variance = FALSE,
      residual_variance          = 1,
      check_R                    = FALSE,
      max_iter                   = SUSIE_MAX_ITER
    ),
    error = function(e) {
      cat("  susie_rss error:", conditionMessage(e), "\n")
      NULL
    }
  )

  # ---- (l) Handle failure ----
  if (is.null(s)) {
    cat("  Skipped (susie_rss returned NULL)\n")
    fwrite(data.table(gene = gene_name, ensg = ensg_id, chr = CHR,
                      n_snps = n_snps, converged = NA, n_cs = 0L, max_pip = NA_real_,
                      status = "susie_failed"),
           file = summary_file, sep = "\t", append = TRUE, col.names = FALSE)
    n_failed <- n_failed + 1L
    next
  }

  # ---- (m) Save RDS ----
  # Stamp the LD provenance into the object.  coloc.susie() assumes the GWAS-side
  # and eQTL-side LBF matrices derive from a COMMON LD reference; the previous
  # fits (2026-04-10/12) predate the production PolyFun layout (2026-04-28) and
  # recorded nothing, so the mismatch was undetectable. Extra list elements are
  # ignored by susieR/coloc.
  s$ld_provenance <- list(
    ld_panel      = Sys.getenv("LD_PANEL", unset = "polyfun"),
    ld_base_dir   = get_ld_base_dir("EUR"),
    ld_regularize = LD_REGULARIZE,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    script        = "08_run_eqtl_susie.R"
  )
  saveRDS(s, rds_path)

  # ---- (n) Append to summary log ----
  gene_n_cs  <- if (!is.null(s$sets$cs)) length(s$sets$cs) else 0L
  gene_maxpip <- if (!is.null(s$pip) && length(s$pip) > 0) max(s$pip, na.rm = TRUE) else NA_real_

  fwrite(data.table(gene = gene_name, ensg = ensg_id, chr = CHR,
                    n_snps = n_snps, converged = isTRUE(s$converged),
                    n_cs = gene_n_cs, max_pip = gene_maxpip,
                    status = "ok"),
         file = summary_file, sep = "\t", append = TRUE, col.names = FALSE)

  n_ok <- n_ok + 1L
  if (gene_n_cs > 0L) n_cs <- n_cs + 1L

  cat("  OK — converged:", isTRUE(s$converged),
      "| CS:", gene_n_cs,
      "| max PIP:", round(gene_maxpip, 3), "\n")

  # Memory cleanup — critical for gene-dense chromosomes (chr6 HLA)
  rm(gene_dt, ss, ss_filtered, R_mat, z, s)
  if (i %% 50 == 0) gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------

cat("\n=== chr", CHR, "complete ===\n")
cat("  Total eGenes   :", n_genes,   "\n")
cat("  Already done   :", n_skipped, "\n")
cat("  SuSiE fit OK   :", n_ok,      "\n")
cat("  With CS (>=1)  :", n_cs,      "\n")
cat("  Failed/skipped :", n_failed,  "\n")
cat("  Summary TSV    :", summary_file, "\n")
cat("End:", format(Sys.time()), "\n")
