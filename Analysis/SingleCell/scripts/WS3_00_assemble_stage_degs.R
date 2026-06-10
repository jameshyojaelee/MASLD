#!/usr/bin/env Rscript
# =============================================================================
# WS3_00_assemble_stage_degs.R
#
# Assemble the master stage-specific-DEG input table for the cell-type routing
# analysis (WS3). Pools the now-canonical (LVQW / limma-voom-qw-C2 harmonized)
# bulk stage / progression DEG sets into ONE tidy long table that the three
# downstream routing axes (expression / LIANA / Hotspot) consume.
#
# READ-ONLY on all DEG source files. WRITES only:
#   - stage_deg_sets.csv         (tidy long master table)
#   - stage_deg_sets_summary.csv (per stage_set DEG counts)
#
# Two axes:
#   COARSE (matches scRNA disease_stage_coarse: Healthy->Steatosis->SH->Cirrhosis)
#       dream_results_stage_{steatosis,sh,cirrhosis}.csv  (versioned ENSG)
#   FINE  (fibrosis F0->F4)
#       disease_signatures/fibrosis_pairwise.csv  (long, `transition`; versioned)
#       progression/c{3,8,9}_*_dream.csv           (UNVERSIONED ENSG)
#       integration/dream_results_ordinal.csv      (long, `transition`+`term`; versioned)
#
# Gene key harmonized to UNVERSIONED ENSG so versioned & unversioned join cleanly.
# Tier-2 DEG definition: padj < 0.05, NO LFC filter.
#
# Env: rnaseq.  Light I/O -> login node, no SLURM.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths ------------------------------------------------------------------
PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RES  <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
META <- file.path(PROJ, "data/gencode_v49_gene_metadata.tsv.gz")
OUT_DIR <- file.path(PROJ, "Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THRESH <- 0.05  # Tier-2: padj only, no LFC filter

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

# ---- Symbol map (unversioned ENSG -> gene_name) -----------------------------
message("Loading GENCODE v49 gene metadata for symbol mapping ...")
meta <- fread(META, sep = "\t")
stopifnot(all(c("gene_id", "gene_name", "ensembl_base") %in% names(meta)))
# ensembl_base is already unversioned; de-dup keeping first symbol per base id
sym_map <- meta[!duplicated(ensembl_base), .(gene = ensembl_base, symbol = gene_name)]
setkey(sym_map, gene)
message(sprintf("  symbol map: %d unique unversioned ENSG", nrow(sym_map)))

# ---- Generic loader: returns standardized long rows for one stage_set -------
# Reads a "wide" per-gene DE table (one row per gene) and tags it stage_set/axis.
load_wide <- function(path, stage_set, axis) {
  if (!file.exists(path)) {
    warning(sprintf("MISSING source: %s (stage_set=%s) -- skipped", path, stage_set))
    return(NULL)
  }
  dt <- fread(path)
  req <- c("gene", "logFC", "padj")
  if (!all(req %in% names(dt))) {
    warning(sprintf("UNEXPECTED SCHEMA in %s -- missing %s -- skipped",
                    path, paste(setdiff(req, names(dt)), collapse = ",")))
    return(NULL)
  }
  dt[, gene := strip_ver(gene)]
  dt[, stage_set := stage_set]
  dt[, axis := axis]
  if (!"shrunk_logFC" %in% names(dt)) dt[, shrunk_logFC := NA_real_]
  if (!"lfsr"        %in% names(dt)) dt[, lfsr := NA_real_]
  dt[, .(stage_set, axis, gene, logFC, padj, shrunk_logFC, lfsr)]
}

# ---- Loader for long tables split by a grouping column ----------------------
# Splits one long file into multiple stage_sets via a label function over rows.
load_long_split <- function(path, axis, label_fun) {
  if (!file.exists(path)) {
    warning(sprintf("MISSING source: %s -- skipped", path))
    return(NULL)
  }
  dt <- fread(path)
  req <- c("gene", "logFC", "padj")
  if (!all(req %in% names(dt))) {
    warning(sprintf("UNEXPECTED SCHEMA in %s -- missing %s -- skipped",
                    path, paste(setdiff(req, names(dt)), collapse = ",")))
    return(NULL)
  }
  dt[, gene := strip_ver(gene)]
  dt[, stage_set := label_fun(dt)]
  dt[, axis := axis]
  if (!"shrunk_logFC" %in% names(dt)) dt[, shrunk_logFC := NA_real_]
  if (!"lfsr"        %in% names(dt)) dt[, lfsr := NA_real_]
  dt <- dt[!is.na(stage_set)]
  dt[, .(stage_set, axis, gene, logFC, padj, shrunk_logFC, lfsr)]
}

# =============================================================================
# COARSE axis : stage-vs-Healthy (one coef per scRNA coarse stage)
# =============================================================================
coarse_pieces <- list(
  load_wide(file.path(RES, "integration/dream_results_stage_steatosis.csv"),
            "coarse_Steatosis", "coarse"),
  load_wide(file.path(RES, "integration/dream_results_stage_sh.csv"),
            "coarse_SH", "coarse"),
  load_wide(file.path(RES, "integration/dream_results_stage_cirrhosis.csv"),
            "coarse_Cirrhosis", "coarse")
)

# =============================================================================
# FINE axis
# =============================================================================
# (a) fibrosis_pairwise: long, `transition` = F0_to_F1 .. F3_to_F4
fine_fib <- load_long_split(
  file.path(RES, "disease_signatures/fibrosis_pairwise.csv"),
  axis = "fine",
  label_fun = function(dt) {
    # F0_to_F1 -> fine_F0_F1
    paste0("fine_", sub("_to_", "_", dt[["transition"]]))
  }
)

# (b) ordinal: long, `term` in {fib_num (linear_trend), fib_grouphigh (F.vF.)}
#     -> ordinal_trend + ordinal_F0_F1 .. ordinal_F3_F4
fine_ord <- load_long_split(
  file.path(RES, "integration/dream_results_ordinal.csv"),
  axis = "fine",
  label_fun = function(dt) {
    tr   <- dt[["transition"]]   # linear_trend | F0vF1 | F1vF2 | F2vF3 | F3vF4
    out  <- character(nrow(dt))
    is_trend <- tr == "linear_trend"
    out[is_trend] <- "ordinal_trend"
    # FxvFy -> ordinal_Fx_Fy
    out[!is_trend] <- paste0("ordinal_", sub("vF", "_F", tr[!is_trend]))
    out
  }
)

# (c) progression key contrasts (UNVERSIONED ENSG already)
prog_pieces <- list(
  load_wide(file.path(RES, "progression/c3_adv_vs_early_fib_dream.csv"),
            "prog_c3_adv_fib", "fine"),
  load_wide(file.path(RES, "progression/c8_cirrhosis_dream.csv"),
            "prog_c8_cirrhosis", "fine"),
  load_wide(file.path(RES, "progression/c9_f2_inflection_dream.csv"),
            "prog_c9_f2_inflection", "fine")
)

# =============================================================================
# Bind everything
# =============================================================================
all_pieces <- c(coarse_pieces, list(fine_fib, fine_ord), prog_pieces)
all_pieces <- all_pieces[!vapply(all_pieces, is.null, logical(1))]
master <- rbindlist(all_pieces, use.names = TRUE)

# ---- Attach symbol ----------------------------------------------------------
master[sym_map, symbol := i.symbol, on = "gene"]
# Fallback: leave symbol NA if not in GENCODE (keeps ENSG-only ncRNA etc.)
master[is.na(symbol) | symbol == "", symbol := NA_character_]

# ---- Derived columns --------------------------------------------------------
master[, direction := fifelse(logFC > 0, "up",
                       fifelse(logFC < 0, "down", NA_character_))]
master[, is_deg := !is.na(padj) & padj < PADJ_THRESH]

# Column order
setcolorder(master, c("stage_set", "axis", "gene", "symbol",
                      "logFC", "padj", "shrunk_logFC", "lfsr",
                      "direction", "is_deg"))

# ---- Write master -----------------------------------------------------------
out_master <- file.path(OUT_DIR, "stage_deg_sets.csv")
fwrite(master, out_master)

# =============================================================================
# Per stage_set summary
# =============================================================================
summ <- master[, .(
  axis          = axis[1],
  n_genes_tested = .N,
  n_DEG         = sum(is_deg, na.rm = TRUE),
  n_up          = sum(is_deg & direction == "up",   na.rm = TRUE),
  n_down        = sum(is_deg & direction == "down", na.rm = TRUE)
), by = stage_set]
setorder(summ, axis, stage_set)

out_summ <- file.path(OUT_DIR, "stage_deg_sets_summary.csv")
fwrite(summ, out_summ)

# =============================================================================
# Verification / report
# =============================================================================
cat("\n================ STAGE-DEG ROUTING TABLE SUMMARY ================\n")
print(summ)

n_with_symbol <- master[!is.na(symbol), uniqueN(gene)]
cat(sprintf("\nUnique genes total            : %d\n", master[, uniqueN(gene)]))
cat(sprintf("Unique genes WITH a symbol    : %d  (target > 20000: %s)\n",
            n_with_symbol, ifelse(n_with_symbol > 20000, "PASS", "FAIL")))

# Confirm coarse and fine join on the unversioned key
coarse_genes <- unique(master[axis == "coarse", gene])
fine_genes   <- unique(master[axis == "fine",   gene])
n_join <- length(intersect(coarse_genes, fine_genes))
cat(sprintf("Coarse/fine shared genes (join): %d  (coarse=%d, fine=%d)\n",
            n_join, length(coarse_genes), length(fine_genes)))

# Confirm no all-NA stage_set (every set has at least some non-NA logFC/padj)
na_check <- master[, .(all_na = all(is.na(logFC)) | all(is.na(padj))), by = stage_set]
bad_sets <- na_check[all_na == TRUE, stage_set]
if (length(bad_sets) > 0) {
  cat(sprintf("\nWARNING: all-NA stage_set(s): %s\n", paste(bad_sets, collapse = ", ")))
} else {
  cat("No all-NA stage_set: PASS\n")
}

cat("\nOutputs:\n")
cat(sprintf("  master : %s  (%d rows)\n", out_master, nrow(master)))
cat(sprintf("  summary: %s  (%d stage_sets)\n", out_summ, nrow(summ)))
cat("================================================================\n")
