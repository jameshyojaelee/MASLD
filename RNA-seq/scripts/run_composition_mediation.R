#!/usr/bin/env Rscript
# ============================================================================
# run_composition_mediation.R
#
# Formal causal mediation decomposition of bulk F-transition DEGs into
# composition-mediated (via cell-type proportion shifts) vs. cell-type-intrinsic
# (residual / direct effect) components.
#
# Methodologically novel for MASLD bulk RNA-seq. Vellame 2021 (PMID 34353365)
# established the precedent for NAFLD methylation. This is the first-time
# application to MASLD bulk RNA-seq.
#
# Algorithm (Imai & Keele 2010; Baron & Kenny 1986 with mediator vector):
#   For each gene g, F-transition T (e.g. F1->F2):
#     1. TE: Y_g ~ stage + dataset + sex          -> beta_TE = stage coef
#     2. DE: Y_g ~ stage + pi_k + dataset + sex   -> beta_DE = stage coef
#                  (mediators pi_k = BayesPrism proportions for major CTs)
#     3. IE = TE - DE   (composition-mediated effect)
#     4. MP = IE / TE clipped to [0,1]
#     5. Per-CT contribution: refit DE_k with each pi_k separately
#        IE_k = TE - DE_k; this attributes mediation to one CT.
#     6. Bootstrap (n=1000): resample sample IDs within (dataset, stage)
#        strata; recompute MP; report 95% percentile CI.
#
# Inputs:
#   - merged_dge.rds                                   (bulk counts)
#   - bayesprism_proportions.csv                       (1444 x 13)
#   - fibrosis_consecutive_dream.csv                   (per-transition DEGs)
#   - dream_results_ashr.csv                           (gene universe + symbols)
#   - unified_metadata.csv + sample_qc_report.csv      (sample filter)
#   - gencode_v49_gene_metadata.tsv.gz                 (protein-coding filter)
#
# Outputs (RNA-seq/.../results/mediation/):
#   - F_transition_mediation_per_gene.csv
#       columns: gene, symbol, transition, beta_TE, beta_DE, IE, MP,
#                MP_CI_lo, MP_CI_hi, top_mediating_CT, top_CT_IE_share,
#                n_donor
#   - mediation_per_celltype.csv
#       columns: gene, symbol, transition, cell_type, IE_k, MP_k
#
# Refs:
#   Vellame 2021 Clin Epigenetics PMID 34353365
#   Imai & Keele 2010 Stat Sci
#   Meng 2023 Brief Bioinformatics PMID 36472568
#   BayesPrism: Chu 2022 Nat Cancer PMID 35469013
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(parallel)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
SIGS   <- file.path(INT, "results/disease_signatures")
META_F <- file.path(INT, "metadata/unified_metadata.csv")
QC_F   <- file.path(INT, "qc/sample_qc_report.csv")
PROP_F <- file.path(INT, "results/progression/cibersortx_celltype_expression",
                    "bayesprism_proportions.csv")
GENC_F <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
OUTDIR <- file.path(INT, "results/mediation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Tunables (override via environment)
N_TOP   <- as.integer(Sys.getenv("MED_N_TOP",  "200"))     # genes per transition
N_BOOT  <- as.integer(Sys.getenv("MED_N_BOOT", "1000"))    # bootstrap iters
N_CPUS  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                 Sys.getenv("MED_N_CPUS", "8")))
PADJ_TH <- as.numeric(Sys.getenv("MED_PADJ", "0.05"))
SEED    <- as.integer(Sys.getenv("MED_SEED", "42"))

set.seed(SEED)
message(sprintf("[cfg] N_TOP=%d  N_BOOT=%d  N_CPUS=%d  PADJ=%g  SEED=%d",
                N_TOP, N_BOOT, N_CPUS, PADJ_TH, SEED))

# Major mediator cell types (proportions). Use 5 lineage-aware umbrellas.
MEDIATORS <- c("Hepatocyte", "Stellate", "Cholangiocyte", "Macrophage", "Endothelial")
TRANSITIONS <- c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3")

# ---------------------------------------------------------------------------
# 1. Metadata + QC
# ---------------------------------------------------------------------------
message("[1] Loading metadata + QC...")
meta <- fread(META_F)
qc   <- fread(QC_F)
meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id",
              all.x = TRUE)
meta_keep <- meta[!is.na(pass_technical) & pass_technical == TRUE &
                    !is.na(fibrosis_stage) & !is.na(sex)]
meta_keep[, fibrosis_stage := as.integer(fibrosis_stage)]
message(sprintf("  QC-pass + staged + sexed: %d samples", nrow(meta_keep)))

# ---------------------------------------------------------------------------
# 2. BayesPrism proportions
# ---------------------------------------------------------------------------
message("[2] Loading BayesPrism proportions...")
prop <- fread(PROP_F)
stopifnot(all(MEDIATORS %in% names(prop)))
prop <- prop[, c("sample_id", MEDIATORS), with = FALSE]
# CLR-style centering (ALR vs Hepatocyte; avoids simplex collinearity in lm)
# Approach: drop one CT to break sum-to-1 collinearity (here: Endothelial),
# keep raw proportions for the remaining 4 mediators (small fractions ~1-30%).
prop_lm <- copy(prop)
DROP_CT <- "Endothelial"  # smallest by mean prop among mediators
prop_lm[, (DROP_CT) := NULL]
mediator_cols <- setdiff(MEDIATORS, DROP_CT)
message(sprintf("  Mediators in lm: %s (dropped %s for identifiability)",
                paste(mediator_cols, collapse = ", "), DROP_CT))

# ---------------------------------------------------------------------------
# 3. Bulk counts -> log2(CPM+1) on protein-coding, detected genes
# ---------------------------------------------------------------------------
message("[3] Loading merged DGE + voom-style logCPM...")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Protein-coding filter
genc <- fread(GENC_F)
pc_ids <- genc[gene_biotype == "protein_coding", gene_id]
keep_pc <- rownames(dge) %in% pc_ids
message(sprintf("  Protein-coding gene retention: %d / %d (%.1f%%)",
                sum(keep_pc), nrow(dge), 100 * mean(keep_pc)))
dge <- dge[keep_pc, , keep.lib.sizes = FALSE]

# Sample subset to QC-pass staged sexed
keep_samp <- colnames(dge) %in% meta_keep$sample_id
dge <- dge[, keep_samp]
message(sprintf("  Bulk samples kept: %d", ncol(dge)))

# Detected: mean log2(CPM+1) > 0.5 (~10 reads/sample)
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
gene_mean <- rowMeans(log_cpm)
keep_det <- gene_mean > 0.5
log_cpm <- log_cpm[keep_det, , drop = FALSE]
message(sprintf("  Detected genes: %d", nrow(log_cpm)))

# Strip ENSG version for downstream matching
ensembl_clean <- sub("\\..*", "", rownames(log_cpm))

# Symbol map from dream_results_ashr.csv (already symbol-mapped) + GENCODE
dream <- fread(file.path(RDIR, "dream_results_ashr.csv"))
sym_map <- unique(dream[, .(gene, symbol)])
sym_map[, ensembl_clean := sub("\\..*", "", gene)]
sym_lookup <- setNames(sym_map$symbol, sym_map$ensembl_clean)

# ---------------------------------------------------------------------------
# 4. Build sample-level data table (sample x [stage, dataset, sex, mediators])
# ---------------------------------------------------------------------------
message("[4] Building sample-level table...")
samp_tbl <- meta_keep[match(colnames(log_cpm), sample_id),
                       .(sample_id, dataset, sex, fibrosis_stage)]
samp_tbl[, dataset := factor(dataset)]
samp_tbl[, sex := factor(sex)]

# Merge mediators
samp_tbl <- merge(samp_tbl, prop_lm, by = "sample_id", sort = FALSE)
# Re-align to log_cpm order
samp_tbl <- samp_tbl[match(colnames(log_cpm), sample_id)]
stopifnot(all(samp_tbl$sample_id == colnames(log_cpm)))

# Cap proportions at small floor to avoid degenerate values
for (c in mediator_cols) {
  samp_tbl[, (c) := pmax(get(c), 1e-6)]
  # Logit-link is over-engineering for small props; raw fraction works.
}

# ---------------------------------------------------------------------------
# 5. Top-N gene sets per transition
# ---------------------------------------------------------------------------
message("[5] Selecting top-N DEGs per F-transition...")
fib <- fread(file.path(SIGS, "fibrosis_consecutive_dream.csv"))
fib[, ensembl_clean := sub("\\..*", "", gene)]
expr_genes <- ensembl_clean
fib <- fib[ensembl_clean %in% expr_genes]
fib <- fib[!is.na(padj) & padj < PADJ_TH]
# Source CSV sometimes carries multiple rows per (gene, contrast) due to
# auxiliary stage P-values stashed in extra columns; collapse to best |t|.
fib <- fib[, .SD[which.max(abs(t))],
           by = .(contrast, ensembl_clean),
           .SDcols = c("gene", "logFC", "padj", "t")]

top_by_trans <- fib[, {
  o <- order(-abs(t))
  head(.SD[o], N_TOP)
}, by = contrast, .SDcols = c("ensembl_clean", "gene", "logFC", "padj", "t")]
message("  Top-N per transition:")
print(top_by_trans[, .N, by = contrast])

# Map ensembl_clean to row indices in log_cpm (use first match per ensembl base)
ens_to_idx <- split(seq_len(nrow(log_cpm)), ensembl_clean)
get_row <- function(ens) {
  idx <- ens_to_idx[[ens]]
  if (is.null(idx)) return(NULL)
  idx[1]
}

# ---------------------------------------------------------------------------
# 6. Mediation core
# ---------------------------------------------------------------------------
# stage_indicator = 1 for "post" stage (e.g. F1 in F1_vs_F0), 0 for "pre"
trans_pre  <- c("F1_vs_F0" = 0L, "F2_vs_F1" = 1L,
                "F3_vs_F2" = 2L, "F4_vs_F3" = 3L)
trans_post <- c("F1_vs_F0" = 1L, "F2_vs_F1" = 2L,
                "F3_vs_F2" = 3L, "F4_vs_F3" = 4L)

# Subset table for given transition
make_trans_data <- function(transition) {
  pre  <- trans_pre[[transition]]
  post <- trans_post[[transition]]
  idx <- which(samp_tbl$fibrosis_stage %in% c(pre, post))
  out <- samp_tbl[idx]
  out[, stage := as.integer(fibrosis_stage == post)]
  list(idx = idx, tbl = out)
}

# Single-gene mediation point estimate
# Returns: beta_TE, beta_DE_full, IE_full, MP_full, IE_per_CT (named vec)
fit_mediation_one <- function(y, tbl) {
  # tbl has cols stage, dataset, sex, mediator_cols
  # If degenerate (e.g. only one dataset level), drop dataset
  use_dataset <- length(unique(tbl$dataset)) > 1
  use_sex     <- length(unique(tbl$sex)) > 1

  base_rhs <- "stage"
  if (use_dataset) base_rhs <- paste(base_rhs, "+ dataset")
  if (use_sex)     base_rhs <- paste(base_rhs, "+ sex")

  fit_te <- lm(as.formula(paste("y ~", base_rhs)), data = tbl)
  beta_te <- coef(fit_te)[["stage"]]

  # Full mediator model
  full_rhs <- paste(base_rhs, "+", paste(mediator_cols, collapse = " + "))
  fit_de <- lm(as.formula(paste("y ~", full_rhs)), data = tbl)
  beta_de <- coef(fit_de)[["stage"]]

  ie_full <- beta_te - beta_de
  mp_full <- ifelse(abs(beta_te) < 1e-8, NA_real_, ie_full / beta_te)

  # Per-CT (each mediator alone)
  ie_per_ct <- vapply(mediator_cols, function(ct) {
    rhs_k <- paste(base_rhs, "+", ct)
    fit_k <- lm(as.formula(paste("y ~", rhs_k)), data = tbl)
    beta_de_k <- coef(fit_k)[["stage"]]
    beta_te - beta_de_k
  }, numeric(1))

  list(beta_te = beta_te, beta_de = beta_de,
       ie = ie_full, mp = mp_full,
       ie_per_ct = ie_per_ct)
}

# Bootstrap CI on MP (resample within strata of (dataset, stage) to preserve
# design balance)
bootstrap_mp <- function(y, tbl, n_boot) {
  strat_key <- paste(tbl$dataset, tbl$stage, sep = "::")
  strat_idx <- split(seq_along(strat_key), strat_key)

  mp_boot <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    rs <- unlist(lapply(strat_idx, function(ix) {
      if (length(ix) == 0) return(integer(0))
      sample(ix, length(ix), replace = TRUE)
    }))
    tbl_b <- tbl[rs]
    y_b   <- y[rs]
    # Quick fits — wrap in tryCatch (singular design from resample possible)
    fit <- tryCatch({
      use_dataset <- length(unique(tbl_b$dataset)) > 1
      use_sex     <- length(unique(tbl_b$sex)) > 1
      base_rhs <- "stage"
      if (use_dataset) base_rhs <- paste(base_rhs, "+ dataset")
      if (use_sex)     base_rhs <- paste(base_rhs, "+ sex")
      fit_te <- lm(as.formula(paste("y_b ~", base_rhs)), data = tbl_b)
      beta_te <- coef(fit_te)[["stage"]]
      full_rhs <- paste(base_rhs, "+", paste(mediator_cols, collapse = " + "))
      fit_de <- lm(as.formula(paste("y_b ~", full_rhs)), data = tbl_b)
      beta_de <- coef(fit_de)[["stage"]]
      if (abs(beta_te) < 1e-8) NA_real_
      else (beta_te - beta_de) / beta_te
    }, error = function(e) NA_real_)
    mp_boot[b] <- fit
  }
  mp_boot <- mp_boot[!is.na(mp_boot)]
  if (length(mp_boot) < 50) return(c(NA_real_, NA_real_))
  unname(quantile(mp_boot, c(0.025, 0.975)))
}

# ---------------------------------------------------------------------------
# 7. Run mediation across (transition, gene) cells, parallelized
# ---------------------------------------------------------------------------
all_jobs <- list()
for (tr in TRANSITIONS) {
  td <- make_trans_data(tr)
  if (nrow(td$tbl) < 30) {
    message(sprintf("  Skipping %s (n=%d)", tr, nrow(td$tbl)))
    next
  }
  genes <- top_by_trans[contrast == tr, ensembl_clean]
  for (g in genes) {
    ri <- get_row(g)
    if (is.null(ri)) next
    all_jobs[[length(all_jobs) + 1L]] <- list(transition = tr, ensembl = g,
                                              row_idx = ri,
                                              tbl_idx = td$idx)
  }
}
message(sprintf("[6] Mediation jobs: %d", length(all_jobs)))

# Pre-make per-transition tbls (avoid repeated subsetting in workers)
trans_tbls <- lapply(TRANSITIONS, function(tr) make_trans_data(tr)$tbl)
names(trans_tbls) <- TRANSITIONS

run_one <- function(job) {
  ri <- job$row_idx
  tr <- job$transition
  tbl <- trans_tbls[[tr]]
  y <- log_cpm[ri, job$tbl_idx]
  pt <- tryCatch(fit_mediation_one(y, tbl), error = function(e) NULL)
  if (is.null(pt)) return(NULL)
  ci <- tryCatch(bootstrap_mp(y, tbl, N_BOOT),
                 error = function(e) c(NA_real_, NA_real_))
  list(transition = tr, ensembl = job$ensembl,
       beta_te = pt$beta_te, beta_de = pt$beta_de,
       ie = pt$ie, mp = pt$mp,
       mp_ci_lo = ci[1], mp_ci_hi = ci[2],
       ie_per_ct = pt$ie_per_ct,
       n_donor = nrow(tbl))
}

t0 <- Sys.time()
message(sprintf("[7] Running mediation in parallel (%d workers, mclapply)...",
                N_CPUS))
results <- mclapply(all_jobs, run_one, mc.cores = N_CPUS,
                    mc.preschedule = TRUE)
results <- Filter(Negate(is.null), results)
message(sprintf("  Completed %d / %d jobs in %.1f min", length(results),
                length(all_jobs),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------------------
# 8. Assemble outputs
# ---------------------------------------------------------------------------
message("[8] Assembling output tables...")
gene_rows <- rbindlist(lapply(results, function(r) {
  ie_k <- r$ie_per_ct
  share <- abs(ie_k) / sum(abs(ie_k))
  top_k <- names(ie_k)[which.max(abs(ie_k))]
  data.table(
    gene        = r$ensembl,
    symbol      = unname(sym_lookup[r$ensembl]),
    transition  = r$transition,
    beta_TE     = r$beta_te,
    beta_DE     = r$beta_de,
    IE          = r$ie,
    MP          = pmin(pmax(r$mp, 0), 1),
    MP_raw      = r$mp,
    MP_CI_lo    = r$mp_ci_lo,
    MP_CI_hi    = r$mp_ci_hi,
    top_mediating_CT = top_k,
    top_CT_IE_share  = share[top_k],
    n_donor     = r$n_donor
  )
}), fill = TRUE)

per_ct_rows <- rbindlist(lapply(results, function(r) {
  data.table(
    gene = r$ensembl,
    symbol = unname(sym_lookup[r$ensembl]),
    transition = r$transition,
    cell_type  = names(r$ie_per_ct),
    IE_k       = unname(r$ie_per_ct),
    MP_k_raw   = ifelse(abs(r$beta_te) < 1e-8, NA_real_,
                        unname(r$ie_per_ct) / r$beta_te),
    MP_k       = pmin(pmax(ifelse(abs(r$beta_te) < 1e-8, NA_real_,
                        unname(r$ie_per_ct) / r$beta_te), 0), 1)
  )
}))

out_main <- file.path(OUTDIR, "F_transition_mediation_per_gene.csv")
out_ct   <- file.path(OUTDIR, "mediation_per_celltype.csv")
fwrite(gene_rows,  out_main)
fwrite(per_ct_rows, out_ct)
message("  Wrote: ", out_main)
message("  Wrote: ", out_ct)

# ---------------------------------------------------------------------------
# 9. Summary
# ---------------------------------------------------------------------------
message("\n== Summary ==")
message(sprintf("Genes-x-transitions analyzed: %d", nrow(gene_rows)))
mp_bins <- gene_rows[!is.na(MP), .(
  pct_high   = mean(MP > 0.7) * 100,
  pct_mid    = mean(MP >= 0.3 & MP <= 0.7) * 100,
  pct_low    = mean(MP < 0.3) * 100
), by = transition]
print(mp_bins)

f2_only <- gene_rows[transition == "F2_vs_F1" & !is.na(MP)]
if (nrow(f2_only) > 0) {
  pct_high <- mean(f2_only$MP > 0.7) * 100
  pct_low  <- mean(f2_only$MP < 0.3) * 100
  pct_mid  <- 100 - pct_high - pct_low
  message(sprintf(
    "F1->F2: %.1f%% composition (>70%%), %.1f%% intrinsic (<30%%), %.1f%% mixed",
    pct_high, pct_low, pct_mid))
  top5_intrinsic <- head(f2_only[order(MP)], 5)
  message("Top-5 most-intrinsic (lowest MP) at F1->F2:")
  print(top5_intrinsic[, .(symbol, MP, beta_TE, top_mediating_CT)])
}

message("Done.")
