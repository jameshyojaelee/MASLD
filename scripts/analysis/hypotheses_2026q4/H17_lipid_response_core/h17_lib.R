# H17 shared functions: input loading and limma-voom fits for both systems.
# Prespecification: Analysis/Hypotheses/H17-20261009T015231Z/PRESPEC_H17.md

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
})

H17_SEED <- 20260923L
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PROGRAM_RELEASE <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  "program-context-v2-stage-corrected-candidate-2026-08-13-v3")
PATHS <- list(
  hep_pb = file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk/Hepatocytes_pseudobulk.csv"),
  donor_meta = file.path(PROGRAM_RELEASE, "metadata/donor_metadata_stage_corrected_donor.tsv"),
  run_meta = file.path(PROGRAM_RELEASE, "metadata/donor_metadata_stage_corrected.tsv"),
  run_sex = file.path(BASE, "Analysis/SingleCell/metadata/snrna_donor_sex.csv"),
  membership = file.path(PROGRAM_RELEASE, "hotspot/program_membership_v2.tsv"),
  registry = file.path(PROGRAM_RELEASE, "hotspot/program_registry_v2.tsv"),
  bulk = file.path(BASE, "RNA-seq/results/manuscript_release/candidates",
                   "bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"),
  collapse_lib = file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R")
)
POS_CONTROLS <- c("PLIN2", "PLIN4", "CIDEC", "FABP4")
Q_THR <- 0.05

# ---------------------------------------------------------------- iPSC (A)
load_ipsc <- function(hyp_dir) {
  cnt_path <- file.path(hyp_dir, "download_GSE342568_counts/GSE342568_HepFA-mRNA_rawcounts.txt.gz")
  mtx_path <- file.path(hyp_dir, "download_GSE342568_matrix/GSE342568_series_matrix.txt.gz")
  lines <- readLines(gzfile(mtx_path), encoding = "UTF-8")
  grab <- function(tag) {
    l <- lines[startsWith(lines, paste0(tag, "\t"))]
    stopifnot(length(l) == 1L)
    gsub('"', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  }
  title <- grab("!Sample_title")
  lib <- sub("^Library name: ", "", grab("!Sample_description"))
  gsm <- grab("!Sample_geo_accession")
  sheet <- data.table(gsm = gsm, title = title, library = lib,
                      participant = sub("-(HepC|HepFA)$", "", title),
                      treatment = sub("^.*-(HepC|HepFA)$", "\\1", title))
  stopifnot(all(sheet$treatment %in% c("HepC", "HepFA")),
            !anyDuplicated(sheet$library))
  dt <- fread(cmd = sprintf("gzip -dc %s", shQuote(cnt_path)), sep = "\t",
              header = TRUE, check.names = FALSE)
  sym <- as.character(dt[["Gene Symbol"]])
  mat <- as.matrix(dt[, sheet$library, with = FALSE])
  storage.mode(mat) <- "double"
  stopifnot(!anyNA(mat), all(mat >= 0))
  ok <- !is.na(sym) & nzchar(sym)
  dup <- sym %in% sym[duplicated(sym)]
  drop_log <- data.table(n_rows = length(sym), n_missing_symbol = sum(!ok),
                         n_rows_duplicated_symbol = sum(ok & dup))
  keep <- ok & !dup
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- sym[keep]
  colnames(mat) <- sheet$library
  sheet[, treatment := factor(treatment, levels = c("HepC", "HepFA"))]
  sheet[, participant := factor(participant)]
  list(counts = mat, sheet = sheet, drop_log = drop_log)
}

ipsc_design <- function(participant, treatment) {
  model.matrix(~ participant + treatment)
}

# ---------------------------------------------------------------- in vivo (B)
load_invivo <- function(arm = c("steatotic", "steatosis_only")) {
  arm <- match.arg(arm)
  source(PATHS$collapse_lib)
  srr_to_donor <- build_srr_to_donor_map(BASE)
  pb <- fread(PATHS$hep_pb)
  genes <- pb[[1]]
  pbm <- as.matrix(pb[, -1, with = FALSE])
  rownames(pbm) <- genes
  storage.mode(pbm) <- "double"
  run_meta <- fread(PATHS$run_meta)
  donor_meta <- fread(PATHS$donor_meta)
  run_meta[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]
  runs_used <- intersect(colnames(pbm), run_meta$sample)
  runs_dropped <- setdiff(colnames(pbm), run_meta$sample)
  hep_n <- run_meta[sample %in% runs_used, .(n_hep = sum(n_Hepatocytes, na.rm = TRUE),
                                               n_runs = .N), by = donor]
  d <- merge(donor_meta, hep_n, by.x = "sample", by.y = "donor", all.x = TRUE)
  d[is.na(n_hep), n_hep := 0]
  d[, excl := toupper(as.character(exclude_stage_analysis)) == "TRUE"]
  ctrl <- "Healthy"
  case <- if (arm == "steatotic") c("Steatosis", "Steatohepatitis") else "Steatosis"
  d[, group := fifelse(disease_stage_coarse == ctrl, "control",
                       fifelse(disease_stage_coarse %in% case, "steatotic", NA_character_))]
  roster_log <- list()
  roster_log$all <- nrow(d)
  d <- d[excl == FALSE]; roster_log$after_exclude_flag <- nrow(d)
  d <- d[!is.na(group)]; roster_log$after_group <- nrow(d)
  d <- d[n_hep >= 50]; roster_log$after_min50 <- nrow(d)
  both <- d[, .(g = uniqueN(group)), by = dataset][g == 2, dataset]
  d <- d[dataset %in% both]; roster_log$after_both_groups <- nrow(d)
  # sum runs to donors (only runs present in the run-level metadata)
  pbm <- pbm[, runs_used, drop = FALSE]
  run_donor <- run_meta$donor[match(runs_used, run_meta$sample)]
  keep_runs <- run_donor %in% d$sample
  pbm <- pbm[, keep_runs, drop = FALSE]
  pbm <- collapse_counts_to_donor(pbm, BASE)
  stopifnot(setequal(colnames(pbm), d$sample))
  d <- d[match(colnames(pbm), sample)]
  d[, group := factor(group, levels = c("control", "steatotic"))]
  d[, dataset := factor(dataset)]
  # inferred sex per donor (majority of runs); NA if missing or tied
  sx <- fread(PATHS$run_sex)
  sx[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]
  sx_d <- sx[donor %in% d$sample & inferred_sex %in% c("Male", "Female"),
             .(sex = {tb <- sort(table(inferred_sex), decreasing = TRUE)
                      if (length(tb) > 1 && tb[1] == tb[2]) NA_character_ else names(tb)[1]},
               n_sex_values = uniqueN(inferred_sex)), by = donor]
  d <- merge(d, sx_d, by.x = "sample", by.y = "donor", all.x = TRUE, sort = FALSE)
  d <- d[match(colnames(pbm), sample)]
  list(counts = pbm, meta = d, roster_log = as.data.table(roster_log),
       runs_dropped_not_in_meta = runs_dropped)
}

invivo_design <- function(meta, add_sex = FALSE) {
  if (add_sex) model.matrix(~ dataset + sex + group, data = meta)
  else model.matrix(~ dataset + group, data = meta)
}

# ---------------------------------------------------------------- fitting
filter_genes <- function(counts, group) {
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = group)
  rownames(counts)[keep]
}

# Full fit with eBayes, BH within the tested genes.
fit_full <- function(counts, genes, design, coef) {
  dge <- DGEList(counts = counts[genes, , drop = FALSE])
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr,
             t = tt$t, P.Value = tt$P.Value, q = p.adjust(tt$P.Value, "BH"))
}

# logFC only (bootstrap): same normalization and voom weights, no eBayes.
fit_logfc <- function(counts, genes, design, coef) {
  dge <- DGEList(counts = counts[genes, , drop = FALSE])
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  setNames(fit$coefficients[, coef], rownames(fit$coefficients))
}

# ---------------------------------------------------------------- statistics
sign_agreement <- function(sa, sb) mean(sa == sb)

sign_perm_null <- function(sign_A_sel, sign_B_all, sel_idx, n_perm) {
  set.seed(H17_SEED)
  out <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    sb <- sample(sign_B_all)
    out[i] <- mean(sign_A_sel == sb[sel_idx])
  }
  out
}

p_upper <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))
p_lower <- function(obs, null) (1 + sum(null <= obs)) / (1 + length(null))

# Expression-matched random sets: returns a list of index vectors into the
# universe, reproducing the bin composition of `target_idx`.
matched_sets <- function(bin, target_idx, n_sets) {
  set.seed(H17_SEED)
  need <- table(bin[target_idx])
  pools <- split(seq_along(bin), bin)
  lapply(seq_len(n_sets), function(i) {
    unlist(lapply(names(need), function(b) {
      pool <- pools[[b]]
      if (length(pool) == 1L) pool else sample(pool, need[[b]])
    }), use.names = FALSE)
  })
}
