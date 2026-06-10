#!/usr/bin/env Rscript
# =============================================================================
# metrics_truth_factorial.R
# -----------------------------------------------------------------------------
# CONSUMER scoring for the factorial DE-method benchmark. PURE consumer: it
# reads per-cell DEG result tables already on disk and computes two metric
# tables. NO model is (re)fit here.
#
#   (A) metrics_external.csv  -- ORTHOGONAL EXTERNAL TRUTH  (role RANK_R2):
#         OpenTargets MASLD recall + fold-enrichment OR, mouse cross-species
#         precision, DEG-COLOC genetic enrichment OR. These are the ONLY
#         metrics permitted to enter the Step-2 R2 rank (preregistration.yaml).
#
#   (B) metrics_report.csv    -- HELD-OUT REPORT          (role REPORT):
#         positive-control precision/recall/F1, genetic-control capture,
#         Govaere + Hoang published-panel concordance. NEVER used for selection.
#
# Plus a tidy long table metrics_truth_long.csv merging A+B with an explicit
# `role` column (RANK_R2 / REPORT) matching the frozen pre-registration so a
# REPORT metric can never be mislabeled as a ranking metric.
#
# Ported / reused from (signatures studied, then ported -- not sourced, to keep
# this a self-contained consumer):
#   * Cas13_Library_Design/scripts/benchmark_cutoffs.R  metrics_of  (L66-72)
#       -> OpenTargets recall/enrich, mouse precision; OT panel + atlas paths.
#   * RNA-seq/27b_benchmark_presets.R  evaluate_selection (L209-237)
#       -> positive-control precision/recall/F1 + genetic-control capture.
#   * RNA-seq/206_published_vs_perstudy.R  lfc_correlation/deg_overlap (L225-249)
#       -> Govaere/Hoang published-panel concordance.
#
# Env: micromamba run -n rnaseq Rscript metrics_truth_factorial.R [args]
# CLI:
#   --cells_dir   dir of per-cell CSVs (cell_<cell_id>.csv | <prefix>*.csv)
#   --manifest    manifest CSV (cell_id, engine, correction_id, k_sv, status);
#                 optional in SMOKE -- if absent, every *.csv in cells_dir is a cell.
#   --out_dir     output dir (default = degx_factorial results dir)
#   --pattern     glob/regex for stand-in cells when no manifest (default cell_*.csv)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# 0. paths + CLI
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEGX_DIR <- file.path(BASE, "RNA-seq/results/degx_factorial")

# Truth-source files (resolved against the live repo; see report in the consumer log).
OT_FILE      <- file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv")
ATLAS_FILE   <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
GENEMETA_FILE<- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
PC_FILE      <- file.path(BASE, "RNA-seq/results/validation/positive_control_validation.csv")
GOVAERE_FILE <- file.path(BASE, "data/published_gene_panels/govaere_2020_panel.tsv")
HOANG_XLSX   <- file.path(BASE, "data/published_degs/GSE130970/MOESM2.xlsx")
ANNOT_FILE   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")

# Canonical atlas truth columns (verified present 2026-06-05):
#  - mouse_meta_padj             : mouse cross-species DEG truth (padj < 0.05)
#  - coloc_best_susie_pp4_polyfun: canonical SUSiE-COLOC posterior (PolyFun EUR LD,
#    production reference since 2026-05-06). The pre-registration names this column
#    `coloc_best_susie_pp4`; that bare column is not materialized in the live atlas,
#    so the production PolyFun equivalent is used and recorded here.
ATLAS_MOUSE_COL <- "mouse_meta_padj"
ATLAS_COLOC_COL <- "coloc_best_susie_pp4_polyfun"
COLOC_THRESH    <- 0.5

# DEG-set definition (canonical Tier-1, prereg comparability mode fixed_threshold).
DEG_PADJ <- 0.05
DEG_LFC  <- 0.5

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(cells_dir = file.path(DEGX_DIR, "cells"),
              manifest  = file.path(DEGX_DIR, "manifest_disease_vs_control.csv"),
              out_dir   = DEGX_DIR,
              pattern   = "cell_*.csv")
  i <- 1
  while (i <= length(a)) {
    key <- sub("^--", "", a[i])
    if (key %in% names(out)) { out[[key]] <- a[i + 1]; i <- i + 2 }
    else { i <- i + 1 }
  }
  out
}
args <- parse_args()
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

strip_v <- function(x) sub("[.][0-9]+$", "", x)

cat("=== metrics_truth_factorial.R (CONSUMER) ===\n")
cat("cells_dir :", args$cells_dir, "\n")
cat("manifest  :", args$manifest, if (file.exists(args$manifest)) "" else "(absent -> SMOKE/glob mode)", "\n")
cat("out_dir   :", args$out_dir, "\n\n")

# ---------------------------------------------------------------------------
# 1. TRUTH SOURCES (load once; shared across all cells)
# ---------------------------------------------------------------------------

# 1a. symbol <-> ensembl base map (GENCODE v49). Used to translate atlas symbols
#     and panel symbols into the ensembl-base id space the cells use.
gm <- fread(GENEMETA_FILE)                      # gene_id, gene_name, ..., ensembl_base
gm[, ensembl_base := strip_v(gene_id)]
sym2ens <- gm[!duplicated(gene_name), setNames(ensembl_base, gene_name)]
sym2ens_up <- setNames(sym2ens, toupper(names(sym2ens)))
map_symbols <- function(symbols) {
  symbols <- as.character(symbols)
  out <- unname(sym2ens[symbols])
  miss <- is.na(out)
  if (any(miss)) out[miss] <- unname(sym2ens_up[toupper(symbols[miss])])
  out
}

# 1b. OpenTargets MASLD gold standard (symbol panel). skip="gene_symbol" drops the
#     comment header (benchmark_cutoffs.R L37). Map to ensembl base for joining.
ot <- fread(OT_FILE, skip = "gene_symbol")
ot_ens <- unique(na.omit(map_symbols(unique(ot$gene_symbol))))
cat(sprintf("OpenTargets: %d symbols -> %d ensembl-base ids\n", length(unique(ot$gene_symbol)), length(ot_ens)))

# 1c. multi-evidence atlas: defines the TESTED-GENE UNIVERSE for the Fisher
#     enrichment tests (mouse + COLOC), plus the mouse / COLOC truth memberships.
atlas <- fread(ATLAS_FILE, select = c("ensembl_id", ATLAS_MOUSE_COL, ATLAS_COLOC_COL))
atlas[, gb := strip_v(ensembl_id)]
atlas <- atlas[!is.na(gb) & gb != ""]
atlas <- atlas[!duplicated(gb)]
# Mouse universe = genes actually TESTED in mouse (non-NA mouse_meta_padj);
# mouse truth = mouse DEG (padj < 0.05).
atlas[, mouse_tested := !is.na(get(ATLAS_MOUSE_COL))]
atlas[, is_mouse_deg := mouse_tested & get(ATLAS_MOUSE_COL) < 0.05]
# COLOC universe = the full atlas gene set (genes carried in the evidence table);
# COLOC truth = SUSiE PP4 > threshold.
atlas[, is_coloc := !is.na(get(ATLAS_COLOC_COL)) & get(ATLAS_COLOC_COL) > COLOC_THRESH]
mouse_universe <- atlas[mouse_tested == TRUE, gb]
mouse_truth    <- atlas[is_mouse_deg == TRUE, gb]
coloc_universe <- atlas[, gb]
coloc_truth    <- atlas[is_coloc == TRUE, gb]
# OpenTargets enrichment is computed over the same atlas-tested universe so the
# enrichment denominator is the tested set, NOT the whole genome.
ot_universe <- atlas[, gb]
ot_truth    <- intersect(ot_ens, ot_universe)
N_OT_UNIV   <- length(ot_truth)
cat(sprintf("Atlas universe: %d genes | mouse-tested %d (DEG %d) | COLOC>%.2f %d | OT-in-univ %d\n",
            nrow(atlas), length(mouse_universe), length(mouse_truth),
            COLOC_THRESH, length(coloc_truth), N_OT_UNIV))

# 1d. positive controls (REPORT). Expression_driven = positive set;
#     GWAS_variant/Genetic_risk = negative controls (expected ~0 captured).
#     control_type splits per the prereg (27b_benchmark_presets.R evaluate_selection).
pc <- fread(PC_FILE)
pos_set_sym  <- unique(pc[control_type == "Expression_driven", gene])
# Negative controls: prereg calls them "Genetic_risk"; this file labels them
# "GWAS_variant". Accept either so the role is robust to the label.
neg_labels   <- c("Genetic_risk", "GWAS_variant")
neg_set_sym  <- unique(pc[control_type %in% neg_labels, gene])
pos_set_ens  <- unique(na.omit(map_symbols(pos_set_sym)))
neg_set_ens  <- unique(na.omit(map_symbols(neg_set_sym)))
cat(sprintf("Positive controls: %d expression-driven (-> %d ens) | %d genetic-risk (-> %d ens)\n",
            length(pos_set_sym), length(pos_set_ens), length(neg_set_sym), length(neg_set_ens)))

# 1e. Govaere 25-gene panel (REPORT). Direction-only (no published per-gene
#     magnitude), so only directional concordance is defined; lfc r/rho = NA.
gov <- fread(GOVAERE_FILE, skip = "gene_symbol")
gov[, gene_base := map_symbols(gene_symbol)]
gov[, dir_sign := ifelse(direction == "up", 1, -1)]
gov <- gov[!is.na(gene_base)]
cat(sprintf("Govaere panel: %d genes mapped\n", nrow(gov)))

# 1f. Hoang panel (REPORT). NAS ordinal-regression sheet carries range_log2FC
#     (magnitude) + coefficient sign + adj_P, so both lfc-correlation AND
#     direction concordance are defined (206_published_vs_perstudy.R).
hoang <- NULL
if (requireNamespace("readxl", quietly = TRUE) && file.exists(HOANG_XLSX)) {
  h <- as.data.table(readxl::read_excel(HOANG_XLSX, sheet = "NAS ordinal regression"))
  setnames(h, tolower(names(h)))
  h[, gene_base := map_symbols(gene_symbol)]
  h <- h[!is.na(gene_base) & !is.na(range_log2fc)]
  # Signed published logFC: range_log2fc is a magnitude; sign comes from coefficient.
  h[, hoang_lfc := range_log2fc * sign(coefficient)]
  hoang <- h[!duplicated(gene_base), .(gene_base, hoang_lfc, hoang_padj = adj_p)]
  cat(sprintf("Hoang NAS panel: %d genes (signed range_log2FC)\n", nrow(hoang)))
} else {
  cat("Hoang panel: readxl/xlsx unavailable -> Hoang concordance = NA\n")
}

# 1g. annotation map for any extra symbol->ensembl gaps (integration build).
#     (Currently only used as a fallback; cells are already ensembl-based.)
if (file.exists(ANNOT_FILE)) {
  an <- fread(ANNOT_FILE)
  extra <- an[!toupper(symbol) %in% toupper(names(sym2ens))]
  if (nrow(extra) > 0) {
    add <- extra[!duplicated(symbol)]
    sym2ens_up <- c(sym2ens_up, setNames(strip_v(add$gene_base), toupper(add$symbol)))
  }
}

# ---------------------------------------------------------------------------
# 2. metric primitives (ported)
# ---------------------------------------------------------------------------

# Fisher 2x2 of DEG-set membership vs truth membership over a SHARED universe.
# (a la 206 deg_overlap.) The universe is the atlas truth-defined gene set
# INTERSECTED with the genes this cell actually tested, so the OR is a clean
# enrichment of DEGs among truth genes WITHIN the genes the method could call --
# never the whole genome, and never atlas genes the method never tested.
fisher_enrich <- function(deg_set, truth_set, universe, cell_universe) {
  universe <- intersect(universe, cell_universe)
  deg_u   <- intersect(deg_set, universe)
  truth_u <- intersect(truth_set, universe)
  a <- length(intersect(deg_u, truth_u))         # DEG & truth
  b <- length(setdiff(deg_u, truth_u))           # DEG & !truth
  c <- length(setdiff(truth_u, deg_u))           # !DEG & truth
  d <- length(universe) - a - b - c              # !DEG & !truth
  m <- matrix(c(a, b, c, max(d, 0)), nrow = 2)
  ft <- tryCatch(fisher.test(m), error = function(e) NULL)
  list(a = a, n_truth_univ = length(truth_u),
       or = if (is.null(ft)) NA_real_ else unname(ft$estimate),
       p  = if (is.null(ft)) NA_real_ else ft$p.value)
}

# ---------------------------------------------------------------------------
# 3. cell discovery (manifest preferred; glob fallback for SMOKE)
# ---------------------------------------------------------------------------
get_cells <- function() {
  # Per-cell fan-out writes one suffixed manifest PER cell
  # (manifest_disease_vs_control__<cell>.csv); the bare manifest is a stale
  # snapshot. Glob ALL shard manifests for this contrast, concatenate, and
  # dedupe by cell_id with the SAME priority as build_ranked_table.R so the
  # RUV-fixed re-runs (most-recent mtime / "ON TOP of cohort" / higher
  # n_covariates) win over stale broken rows.
  man_dir <- dirname(args$manifest)
  mfs <- list.files(man_dir, pattern = "^manifest_disease_vs_control.*\\.csv$",
                    full.names = TRUE)
  if (length(mfs) > 0) {
    parts <- lapply(mfs, function(p) {
      dt <- tryCatch(fread(p), error = function(e) NULL)
      if (is.null(dt) || !("cell_id" %in% names(dt)) || nrow(dt) == 0) return(NULL)
      dt[, `:=`(.__mtime = as.numeric(file.info(p)$mtime))]
      dt
    })
    mf <- data.table::rbindlist(Filter(Negate(is.null), parts), fill = TRUE)
    if (nrow(mf) > 0) {
      if ("status" %in% names(mf)) mf <- mf[status == "ok"]
      if (!("n_covariates" %in% names(mf))) mf[, n_covariates := NA_real_]
      if (!("message" %in% names(mf)))      mf[, message := NA_character_]
      mf[, `:=`(.__ncov   = suppressWarnings(as.numeric(n_covariates)),
                .__ruvfix = as.integer(!is.na(message) &
                              grepl("ON TOP of cohort", message, fixed = TRUE)))]
      data.table::setorder(mf, cell_id, -.__mtime, -.__ruvfix, -.__ncov)
      mf <- mf[!duplicated(cell_id)]
      mf[, file := file.path(args$cells_dir, paste0("cell_", cell_id, ".csv"))]
      mf <- mf[file.exists(file)]
      cat(sprintf("Manifest mode: %d ok cells from %d shard manifests\n",
                  nrow(mf), length(mfs)))
      return(mf[, .(cell_id, file)])
    }
  }
  if (file.exists(args$manifest)) {
    mf <- fread(args$manifest)
    stopifnot("cell_id" %in% names(mf))
    if ("status" %in% names(mf)) mf <- mf[status == "ok"]
    mf[, file := file.path(args$cells_dir, paste0("cell_", cell_id, ".csv"))]
    mf <- mf[file.exists(file)]
    cat(sprintf("Manifest mode (bare): %d ok cells with existing files\n", nrow(mf)))
    return(mf[, .(cell_id, file)])
  }
  # glob fallback -- every matching CSV is treated as a stand-in cell.
  pat <- utils::glob2rx(basename(args$pattern))
  fs <- list.files(args$cells_dir, pattern = pat, full.names = TRUE)
  if (length(fs) == 0) {
    # broaden: any *.csv in the dir
    fs <- list.files(args$cells_dir, pattern = "\\.csv$", full.names = TRUE)
  }
  # never treat a manifest table as a cell
  fs <- fs[!grepl("manifest", basename(fs), ignore.case = TRUE)]
  cid <- sub("\\.csv$", "", basename(fs))
  cid <- sub("^cell_", "", cid)
  cat(sprintf("Glob mode (no manifest): %d stand-in cells from %s\n",
              length(fs), args$cells_dir))
  data.table(cell_id = cid, file = fs)
}
cells <- get_cells()
if (nrow(cells) == 0) stop("No cells found to score.")

# ---------------------------------------------------------------------------
# 4. per-cell scoring
# ---------------------------------------------------------------------------
read_cell <- function(f) {
  dt <- fread(f)
  # tolerate a missing SE column (Rexact stand-ins have no SE)
  req <- c("gene", "logFC", "pval", "padj")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop(sprintf("cell %s missing columns: %s", basename(f), paste(miss, collapse = ",")))
  dt[, gb := strip_v(gene)]
  dt[!is.na(gb) & gb != ""]
}

ext_rows <- list()
rep_rows <- list()

for (k in seq_len(nrow(cells))) {
  cid <- cells$cell_id[k]
  dt  <- read_cell(cells$file[k])

  # canonical Tier-1 DEG set (PRIMARY, bidirectional; prereg comparability mode
  # fixed_threshold). All metrics below are computed on this set unless noted.
  is_deg  <- !is.na(dt$padj) & dt$padj < DEG_PADJ & !is.na(dt$logFC) & abs(dt$logFC) > DEG_LFC
  deg_gb  <- dt[is_deg, gb]
  n_deg   <- length(deg_gb)
  cell_univ <- dt$gb              # genes this cell actually tested
  # UP-directional Tier-1 set: matches the benchmark_cutoffs.R metrics_of code
  # anchor (logFC > x). The library is a knockdown screen selecting up-in-disease
  # genes, and the OpenTargets gold standard (genetics/literature) is best
  # recovered by the up-regulated arm, so the OT enrichment is reported on BOTH
  # the bidirectional primary set and this up-directional set.
  is_deg_up <- !is.na(dt$padj) & dt$padj < DEG_PADJ & !is.na(dt$logFC) & dt$logFC > DEG_LFC
  deg_gb_up <- dt[is_deg_up, gb]

  # ---- (A) EXTERNAL TRUTH (RANK_R2) ----
  # OpenTargets: recall (over OT genes the cell could have called) + fold-enrichment
  # OR (Fisher) over the atlas-tested universe intersected with this cell.
  ot_truth_cell <- intersect(ot_truth, cell_univ)
  ot_in_deg <- length(intersect(deg_gb, ot_truth_cell))
  ot_recall <- if (length(ot_truth_cell) > 0) ot_in_deg / length(ot_truth_cell) else NA_real_
  fe_ot     <- fisher_enrich(deg_gb, ot_truth, ot_universe, cell_univ)
  # up-directional OT enrichment (code-anchor-faithful; benchmark_cutoffs.R)
  fe_ot_up  <- fisher_enrich(deg_gb_up, ot_truth, ot_universe, cell_univ)
  # Mouse cross-species precision: fraction of DEG set (restricted to mouse-tested
  # universe) whose mouse ortholog is a mouse DEG. (benchmark_cutoffs.R prec_mouse.)
  deg_mouse_univ <- intersect(deg_gb, mouse_universe)
  mouse_n_tested <- length(deg_mouse_univ)
  mouse_prec     <- if (mouse_n_tested > 0) length(intersect(deg_mouse_univ, mouse_truth)) / mouse_n_tested else NA_real_
  # COLOC genetic enrichment OR (Fisher) over the atlas universe ∩ this cell.
  fe_coloc <- fisher_enrich(deg_gb, coloc_truth, coloc_universe, cell_univ)

  ext_rows[[k]] <- data.table(
    cell_id        = cid,
    n_deg          = n_deg,
    ot_recall      = ot_recall,
    ot_enrichment_or = fe_ot$or,
    ot_fisher_p    = fe_ot$p,
    ot_enrichment_or_up = fe_ot_up$or,   # up-directional (benchmark_cutoffs.R anchor)
    mouse_precision = mouse_prec,
    mouse_n_tested = mouse_n_tested,
    coloc_or       = fe_coloc$or,
    coloc_fisher_p = fe_coloc$p)

  # ---- (B) HELD-OUT REPORT (REPORT) ----
  # Positive controls (expression-driven) precision/recall/F1 over the DEG set.
  # Restrict to controls present in this cell's tested universe so precision/recall
  # are not penalized for genes the method never tested.
  pos_univ    <- intersect(pos_set_ens, cell_univ)
  tp          <- length(intersect(deg_gb, pos_univ))
  pc_prec     <- if (n_deg > 0) tp / n_deg else 0
  pc_recall   <- if (length(pos_univ) > 0) tp / length(pos_univ) else NA_real_
  pc_f1       <- if (!is.na(pc_recall) && (pc_prec + pc_recall) > 0) 2 * pc_prec * pc_recall / (pc_prec + pc_recall) else 0
  # Genetic-risk controls captured (expectation ~0).
  neg_univ    <- intersect(neg_set_ens, cell_univ)
  genetic_cap <- length(intersect(deg_gb, neg_univ))

  # Govaere concordance: direction only (panel has no magnitude).
  gj <- merge(gov[, .(gb = gene_base, dir_sign)], dt[, .(gb, logFC)], by = "gb")
  gj <- gj[is.finite(logFC)]
  gov_dir_pct <- if (nrow(gj) >= 3) mean(sign(gj$logFC) == gj$dir_sign) * 100 else NA_real_
  gov_lfc_r   <- NA_real_   # no published per-gene magnitude in the Govaere panel
  gov_lfc_rho <- NA_real_

  # Hoang concordance: lfc correlation + (implicit) direction via signed published logFC.
  hoang_r <- NA_real_; hoang_rho <- NA_real_
  if (!is.null(hoang)) {
    hj <- merge(hoang[, .(gb = gene_base, hoang_lfc)], dt[, .(gb, logFC)], by = "gb")
    hj <- hj[is.finite(logFC) & is.finite(hoang_lfc)]
    if (nrow(hj) >= 10) {
      hoang_r   <- cor(hj$logFC, hj$hoang_lfc, method = "pearson")
      hoang_rho <- cor(hj$logFC, hj$hoang_lfc, method = "spearman")
    }
  }

  rep_rows[[k]] <- data.table(
    cell_id              = cid,
    pc_precision         = pc_prec,
    pc_recall            = pc_recall,
    pc_f1                = pc_f1,
    genetic_ctrl_captured = genetic_cap,
    govaere_dir_pct      = gov_dir_pct,
    govaere_lfc_r        = gov_lfc_r,
    govaere_lfc_rho      = gov_lfc_rho,
    hoang_lfc_r          = hoang_r,
    hoang_lfc_rho        = hoang_rho)
}

metrics_external <- rbindlist(ext_rows)
metrics_report   <- rbindlist(rep_rows)
setorder(metrics_external, cell_id)
setorder(metrics_report, cell_id)

# ---------------------------------------------------------------------------
# 5. tidy long output -- one row per (cell_id, family, metric, value, role).
#    Roles are FROZEN per preregistration.yaml truth_assignment:
#      EXTERNAL family (OpenTargets/mouse/COLOC) -> RANK_R2
#      REPORT  family (positive controls/published panels) -> REPORT
# ---------------------------------------------------------------------------
melt_role <- function(dt, id, family, role, value_cols) {
  d2 <- copy(dt[, c(id, value_cols), with = FALSE])
  for (col in value_cols) set(d2, j = col, value = as.numeric(d2[[col]]))  # uniform double
  m <- melt(d2, id.vars = id,
            variable.name = "metric", value.name = "value")
  m[, family := family]
  m[, role := role]
  setcolorder(m, c(id, "family", "metric", "value", "role"))
  m
}
ext_long <- melt_role(metrics_external, "cell_id", "external",
                      "RANK_R2",
                      c("n_deg", "ot_recall", "ot_enrichment_or", "ot_fisher_p",
                        "ot_enrichment_or_up",
                        "mouse_precision", "mouse_n_tested", "coloc_or", "coloc_fisher_p"))
rep_long <- melt_role(metrics_report, "cell_id", "report",
                      "REPORT",
                      c("pc_precision", "pc_recall", "pc_f1", "genetic_ctrl_captured",
                        "govaere_dir_pct", "govaere_lfc_r", "govaere_lfc_rho",
                        "hoang_lfc_r", "hoang_lfc_rho"))
# n_deg / *_fisher_p / mouse_n_tested are descriptive context carried on the
# external table; they are NOT truth metrics, so tag them DESCRIPTIVE so they
# can never be mistaken for a RANK_R2 ranking metric. ot_enrichment_or_up is the
# up-directional sensitivity variant (the authoritative RANK_R2 OT metric is the
# bidirectional ot_enrichment_or on the canonical Tier-1 set).
ext_long[metric %in% c("n_deg", "ot_fisher_p", "mouse_n_tested", "coloc_fisher_p",
                       "ot_enrichment_or_up"),
         role := "DESCRIPTIVE"]
metrics_long <- rbind(ext_long, rep_long)
setorder(metrics_long, cell_id, family, metric)

# ---------------------------------------------------------------------------
# 6. write outputs
# ---------------------------------------------------------------------------
ext_path  <- file.path(args$out_dir, "metrics_external.csv")
rep_path  <- file.path(args$out_dir, "metrics_report.csv")
long_path <- file.path(args$out_dir, "metrics_truth_long.csv")
fwrite(metrics_external, ext_path)
fwrite(metrics_report,   rep_path)
fwrite(metrics_long,     long_path)

cat("\n=== wrote ===\n")
cat(" (A) RANK_R2 external truth :", ext_path, sprintf("(%d cells)\n", nrow(metrics_external)))
cat(" (B) REPORT held-out        :", rep_path, sprintf("(%d cells)\n", nrow(metrics_report)))
cat("     tidy long              :", long_path, sprintf("(%d rows)\n", nrow(metrics_long)))

cat("\n--- head metrics_external.csv ---\n")
print(head(metrics_external, 12))
cat("\n--- head metrics_report.csv ---\n")
print(head(metrics_report, 12))

# sanity: OpenTargets enrichment for real DE methods. The bidirectional Tier-1
# set sits near baseline (the OT gold standard is genetics/literature-anchored
# and many of its genes carry modest fold changes that the |logFC|>0.5 floor
# discards); the up-directional arm -- the library's actual operating direction
# and the benchmark_cutoffs.R code anchor -- is enriched (OR>1).
cat(sprintf("\nSanity (bidirectional Tier-1): OpenTargets OR>1 for %d/%d cells\n",
            sum(metrics_external$ot_enrichment_or > 1, na.rm = TRUE), nrow(metrics_external)))
cat(sprintf("Sanity (up-directional, code-anchor): OpenTargets OR>1 for %d/%d cells\n",
            sum(metrics_external$ot_enrichment_or_up > 1, na.rm = TRUE), nrow(metrics_external)))
de_cells <- metrics_external[grepl("dream|limma", cell_id)]
cat("  dream/limma up-directional OT OR: ",
    paste(sprintf("%s=%.3f", sub("Rexact_disease_vs_control_", "", de_cells$cell_id),
                  de_cells$ot_enrichment_or_up), collapse = ", "), "\n")
cat(sprintf("Sanity: mouse_precision range [%.3f, %.3f]\n",
            min(metrics_external$mouse_precision, na.rm = TRUE),
            max(metrics_external$mouse_precision, na.rm = TRUE)))
cat("\nDone.\n")
