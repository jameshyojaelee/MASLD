#!/usr/bin/env Rscript
# KEY MESSAGE: For the two prespecified focal Hotspot programs, does member-gene
# covariance concentration differ between the bottom and top continuum tercile,
# and does any such difference survive residualizing expression on t?

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(splines)
})

options(digits = 17, scipen = 999)
set.seed(20260827)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 81_focal_program_coherence_mvp.R OUTPUT_DIR", call. = FALSE)
out_dir <- args[[1L]]
if (dir.exists(out_dir)) stop("Refusing to overwrite: ", out_dir, call. = FALSE)
dir.create(out_dir, recursive = TRUE)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
REGISTRY <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates",
                      "hac-continuum-20260818T024923Z/projection/score_registry.tsv")
DGE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates",
                 "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
                 "arms/F_five/results/integration/merged_dge.rds")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates",
                   "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE",
                   "frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
SIGEXC <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/molecular_layers",
                    "hac-molecular-layers-20260818T173348Z/programs/hotspot",
                    "hotspot_signature_excluded_membership.tsv.gz")
MEMV2 <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates",
                   "program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv")

FOCAL <- c(ductular_injury = "hotspot_hepatocytes_48f39dd4d817a10e",
           stromal_ecm     = "hotspot_hepatocytes_f05c535ae5bbc0b9")

stopifnot(all(file.exists(c(REGISTRY, DGE, ANNOT, SIGEXC, MEMV2))))

# ---- axis -------------------------------------------------------------------
reg <- fread(REGISTRY, na.strings = c("", "NA"))
stopifnot(nrow(reg) == 844L, !any(grepl('"', reg$dataset, fixed = TRUE)))
reg <- reg[, .(sample_id, dataset, t = fixed_projection_percentile)]
stopifnot(all(is.finite(reg$t)))

# ---- expression -------------------------------------------------------------
dge <- readRDS(DGE)
stopifnot(identical(dim(dge$counts), c(23370L, 844L)))
ann <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
ann <- unique(ann[!is.na(gene_name) & nzchar(gene_name)])
sym_of <- setNames(toupper(trimws(ann$gene_name)), ann$ensembl_base)

# ---- memberships: two definitions -------------------------------------------
se <- fread(SIGEXC, na.strings = c("", "NA"))
v2 <- fread(MEMV2, na.strings = c("", "NA"))
member_sets <- list()
for (nm in names(FOCAL)) {
  uid <- FOCAL[[nm]]
  g_all <- unique(toupper(trimws(v2[program_uid == uid, mapped_symbol])))
  g_exc <- unique(toupper(trimws(se[program_uid == uid & excluded_signature_gene == FALSE, gene_symbol])))
  member_sets[[paste0(nm, "|membership_v2_full")]] <- list(program = nm, uid = uid,
    definition = "membership_v2_full", genes = g_all[nzchar(g_all) & !is.na(g_all)])
  member_sets[[paste0(nm, "|signature_excluded")]] <- list(program = nm, uid = uid,
    definition = "signature_excluded", genes = g_exc[nzchar(g_exc) & !is.na(g_exc)])
}

# ---- statistics --------------------------------------------------------------
coherence <- function(mat, method) {
  # mat: samples x genes. Returns lambda1 fraction and bias-corrected mean
  # signed off-diagonal Fisher z of the gene correlation matrix.
  keep <- apply(mat, 2L, function(v) is.finite(sd(v)) && sd(v) > 0)
  mat <- mat[, keep, drop = FALSE]
  g <- ncol(mat); n <- nrow(mat)
  if (g < 3L || n < 5L) return(list(g = g, n = n, lev = NA_real_, zbar = NA_real_))
  R <- suppressWarnings(cor(mat, method = method))
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  off <- R[upper.tri(R)]
  off <- pmin(pmax(off, -0.999999), 0.999999)
  list(g = g, n = n,
       lev = max(ev) / sum(ev),
       zbar = mean(atanh(off)) + 1 / (n - 1))  # subtract the -1/(n-1) per-pair bias
}

rows <- list()
for (coh in sort(unique(reg$dataset))) {
  meta <- reg[dataset == coh]
  sub <- dge[, colnames(dge$counts) %in% meta$sample_id]
  meta <- meta[match(colnames(sub$counts), sample_id)]
  stopifnot(identical(meta$sample_id, colnames(sub$counts)))

  logcpm <- edgeR::cpm(sub, log = TRUE, prior.count = 1)          # frozen TMM factors
  syms <- sym_of[sub("\\..*$", "", rownames(logcpm))]
  ok <- !is.na(syms)
  expr <- rowsum(logcpm[ok, , drop = FALSE], syms[ok], reorder = FALSE)  # symbol collapse
  expr <- t(expr)                                                  # samples x genes
  expr <- scale(expr)                                              # z within cohort

  tt <- meta$t
  # residualized arm: remove each gene's smooth trajectory in t, within cohort
  design <- cbind(1, ns(tt, df = 4))
  fitted <- design %*% qr.solve(design, expr)
  resid <- expr - fitted

  cut3 <- quantile(tt, c(1/3, 2/3), names = FALSE, type = 8)
  band <- ifelse(tt <= cut3[1], "bottom_tercile",
          ifelse(tt >= cut3[2], "top_tercile", "middle"))

  for (key in names(member_sets)) {
    ms <- member_sets[[key]]
    gsel <- intersect(ms$genes, colnames(expr))
    for (arm in c("raw", "t_residualized")) {
      X <- if (arm == "raw") expr else resid
      for (method in c("pearson", "spearman")) {
        for (bnd in c("bottom_tercile", "top_tercile")) {
          idx <- which(band == bnd)
          st <- coherence(X[idx, gsel, drop = FALSE], method)
          rows[[length(rows) + 1L]] <- data.table(
            program = ms$program, program_uid = ms$uid,
            membership_definition = ms$definition,
            n_member_genes_declared = length(ms$genes),
            n_member_genes_observed = length(gsel),
            dataset = coh, arm = arm, correlation = method, band = bnd,
            n_samples = st$n, g_used = st$g,
            leading_eigenvalue_fraction = st$lev,
            mean_offdiag_fisher_z = st$zbar)
        }
      }
    }
  }
}

res <- rbindlist(rows)
fwrite(res, file.path(out_dir, "focal_program_coherence_terciles.tsv"), sep = "\t")

wide <- dcast(res, program + membership_definition + dataset + arm + correlation +
                n_member_genes_observed ~ band,
              value.var = c("leading_eigenvalue_fraction", "mean_offdiag_fisher_z",
                            "n_samples", "g_used"))
wide[, `:=`(
  delta_lev = leading_eigenvalue_fraction_top_tercile - leading_eigenvalue_fraction_bottom_tercile,
  delta_zbar = mean_offdiag_fisher_z_top_tercile - mean_offdiag_fisher_z_bottom_tercile)]
setorder(wide, program, membership_definition, correlation, arm, dataset)
fwrite(wide, file.path(out_dir, "focal_program_coherence_contrasts.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\n=== TOP MINUS BOTTOM TERCILE ===\n")
print(wide[correlation == "pearson",
           .(program, membership_definition, dataset, arm,
             g = g_used_top_tercile,
             n_bot = n_samples_bottom_tercile, n_top = n_samples_top_tercile,
             lev_bot = round(leading_eigenvalue_fraction_bottom_tercile, 4),
             lev_top = round(leading_eigenvalue_fraction_top_tercile, 4),
             d_lev = round(delta_lev, 4), d_zbar = round(delta_zbar, 4))],
      nrows = 100)
cat("\nrows written:", nrow(res), "\n")
