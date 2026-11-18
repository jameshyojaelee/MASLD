# SSX-B Joint Inverse-Variance EISA Fix (N16, 2026-11-18)
# Fits joint voom(yE, X) and voom(yI, X) on ~ batch + inferred_sex + bin,
# computes inverse-variance exon-minus-intron weights w_gi = (w_{E,gi}^-1 + w_{I,gi}^-1)^-1,
# sets f$Amean <- pmin(rowMeans(vE$E), rowMeans(vI$E)) before eBayes(f, trend = TRUE),
# and evaluates 200 batch x sex stratified permutations with seed = 20261103.

source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
  library(jsonlite)
  library(parallel)
})

SSX <- file.path(ROOT, "Analysis/Hypotheses/SSX-20261102T133212Z")
REV_EISA <- file.path(ROOT, "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk/eisa_adjusted_v3")
OUT <- file.path(SSX, "out_eisa_20261118")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SSX, "logs"), recursive = TRUE, showWarnings = FALSE)

NAMED <- c("GNMT", "MAT1A", "CYP2C19", "COL1A1", "LUM", "IGFBP7", "EPCAM", "TNFSF10")

U <- fread(file.path(SS, "out_s2_qc/units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
stg <- read_stage(file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"))
U <- merge(U, stg, by = "sample_alias")[!is.na(stage_code)]
U <- add_labels(U, U$stage_code)
setorder(U, sample_alias)
stopifnot(nrow(U) == 493L)

fwrite(U[, .(sample_alias, analysis_unit_id, batch, inferred_sex, stage_code, bin)],
       file.path(OUT, "donor_roster.tsv"), sep = "\t")

d <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds"))
E_full <- as.matrix(d$counts[, U$sample_alias, drop = FALSE])
genes_all <- rownames(E_full)

gmeta <- fread(file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))
sym_map <- gmeta$gene_name[match(genes_all, gmeta$gene_id)]
sym_map[is.na(sym_map) | sym_map == ""] <- genes_all[is.na(sym_map) | sym_map == ""]
names(sym_map) <- genes_all

cat("Loading", nrow(U), "intron count files...\n")
I_list <- lapply(U$sample_alias, function(s) {
  z <- fread(file.path(SSX, "out_b_counts", paste0(s, ".intron.counts.gz")),
             header = FALSE, col.names = c("gene", "n"))
  v <- numeric(length(genes_all))
  idx <- match(z$gene, genes_all)
  ok <- !is.na(idx)
  v[idx[ok]] <- as.numeric(z$n[ok])
  v
})
I_full <- do.call(cbind, I_list)
rownames(I_full) <- genes_all
colnames(I_full) <- U$sample_alias

first_intron_genes <- fread(file.path(SSX, "out_b_counts", paste0(U$sample_alias[1], ".intron.counts.gz")),
                            header = FALSE, col.names = c("gene", "n"))$gene

K <- genes_all[rowMeans(E_full >= 10) >= 0.5 & rowMeans(I_full >= 10) >= 0.5]
cat("Filtered genes K:", length(K), "\n")
stopifnot(length(K) > 10000L)

nf_E <- calcNormFactors(DGEList(E_full))$samples$norm.factors
nf_I <- calcNormFactors(DGEList(I_full))$samples$norm.factors
yE <- DGEList(counts = E_full[K, , drop = FALSE], lib.size = colSums(E_full), norm.factors = nf_E)
yI <- DGEList(counts = I_full[K, , drop = FALSE], lib.size = colSums(I_full), norm.factors = nf_I)

info <- data.frame(
  batch = droplevels(factor(U$batch)),
  inferred_sex = droplevels(factor(U$inferred_sex)),
  bin = factor(U$bin, levels = c("F0F1", "F2", "F3F4")),
  row.names = U$sample_alias
)

fit_eisa_joint <- function(des) {
  vE <- limma::voom(yE, des, plot = FALSE)
  vI <- limma::voom(yI, des, plot = FALSE)
  w <- 1 / (1 / vE$weights + 1 / vI$weights)
  RR <- vE$E - vI$E
  AE <- rowMeans(vE$E)
  AI <- rowMeans(vI$E)
  f <- limma::lmFit(RR, des, weights = w)
  f$Amean <- pmin(AE, AI)
  f <- limma::eBayes(f, trend = TRUE)
  tt <- limma::topTable(f, coef = "binF3F4", number = Inf, sort.by = "none")
  tt$FDR <- p.adjust(tt$P.Value, method = "BH")
  list(tt = tt, vE = vE, vI = vI, w = w, AE = AE, AI = AI)
}

X <- model.matrix(~ batch + inferred_sex + bin, data = info)
stopifnot("binF3F4" %in% colnames(X), qr(X)$rank == ncol(X))

# Separate voomWithQualityWeights exon and intron fits for individual exon/intron classification
vqE <- limma::voomWithQualityWeights(yE, X, plot = FALSE)
fE <- limma::eBayes(limma::lmFit(vqE, X))
tE <- limma::topTable(fE, coef = "binF3F4", number = Inf, sort.by = "none")
tE$FDR <- p.adjust(tE$P.Value, method = "BH")

vqI <- limma::voomWithQualityWeights(yI, X, plot = FALSE)
fI <- limma::eBayes(limma::lmFit(vqI, X))
tI <- limma::topTable(fI, coef = "binF3F4", number = Inf, sort.by = "none")
tI$FDR <- p.adjust(tI$P.Value, method = "BH")

obs <- fit_eisa_joint(X)
t_eisa <- obs$tt

rho_dI_dE <- suppressWarnings(cor(tI[K, "logFC"], tE[K, "logFC"], method = "spearman"))
cat(sprintf("Spearman(dI, dE) = %.4f\n", rho_dI_dE))
stopifnot(rho_dI_dE > 0)

cls <- fifelse(t_eisa[K, "FDR"] < 0.05, "exon_intron_ratio_change_supported",
               fifelse(tI[K, "FDR"] < 0.05, "intron_change_supported_ratio_unresolved", "unresolved"))

R <- data.table(
  gene = K,
  symbol = unname(sym_map[K]),
  dE = tE[K, "logFC"],
  dE_FDR = tE[K, "FDR"],
  dI = tI[K, "logFC"],
  dI_FDR = tI[K, "FDR"],
  dE_joint_voom = rowMeans(obs$vE$E[, info$bin == "F3F4"]) - rowMeans(obs$vE$E[, info$bin == "F0F1"]),
  dI_joint_voom = rowMeans(obs$vI$E[, info$bin == "F3F4"]) - rowMeans(obs$vI$E[, info$bin == "F0F1"]),
  dEmI = t_eisa[K, "logFC"],
  eisa_t = t_eisa[K, "t"],
  eisa_P = t_eisa[K, "P.Value"],
  eisa_FDR = t_eisa[K, "FDR"],
  Amean_min = pmin(obs$AE, obs$AI),
  mean_invvar_weight = rowMeans(obs$w),
  class = cls
)

# Attach comparison columns against previous unweighted adjusted v3 and stage-only if present
prev_v3_path <- file.path(REV_EISA, "eisa_all_filtered_genes.tsv.gz")
if (file.exists(prev_v3_path)) {
  prev <- fread(prev_v3_path)
  prev_sub <- prev[, .(gene, dEmI_unweighted_v3 = dEmI, eisa_FDR_unweighted_v3 = eisa_FDR,
                       dEmI_stage_only = dEmI_stage_only, eisa_FDR_stage_only = eisa_FDR_stage_only)]
  R <- merge(R, prev_sub, by = "gene", all.x = TRUE, sort = FALSE)
}

fwrite(R, file.path(OUT, "eisa_all_filtered_genes.tsv.gz"), sep = "\t")

# Replicated gene classes (matching SteatoSITE catalog sidecar)
SC <- fread(file.path(SS, "out_s3b-20261028T035515Z/steatosite_catalog_sidecar.tsv.gz"))
REP <- merge(
  SC[steatosite_stage_state == "replicated", .(gene = gene_id_versioned, symbol, primary_evidence_class, S3a_F3F4 = steatosite_F3F4_log2FC)],
  R[, !"symbol"],
  by = "gene",
  all.x = TRUE
)
REP[is.na(class), class := fifelse(gene %in% first_intron_genes, "intronic signal too low", "no intronic feature (overlap rule)")]
fwrite(REP, file.path(OUT, "replicated_gene_classes.tsv.gz"), sep = "\t")
fwrite(REP[, .N, by = class][order(-N)], file.path(OUT, "class_counts.tsv"), sep = "\t")

NM <- R[symbol %in% NAMED][order(match(symbol, NAMED))]
fwrite(NM, file.path(OUT, "named_genes.tsv"), sep = "\t")

# 200 permutations of bin within batch x inferred_sex (seed = 20261103)
nperm <- 200L
perm_seed <- 20261103L
set.seed(perm_seed)
strata <- paste(info$batch, info$inferred_sex, sep = "::")
strata_idx <- split(seq_len(nrow(info)), strata)
bin_vec <- info$bin

# Pre-generate all 200 stratified label permutations sequentially so RNG is invariant to core count
perm_list <- vector("list", nperm)
draw_labels <- matrix(NA_character_, nrow = nperm, ncol = nrow(info))
colnames(draw_labels) <- rownames(info)
for (p in seq_len(nperm)) {
  perm_bin <- bin_vec
  for (idx in strata_idx) {
    if (length(idx) > 1L) perm_bin[idx] <- bin_vec[sample(idx)]
  }
  perm_list[[p]] <- perm_bin
  draw_labels[p, ] <- as.character(perm_bin)
}

n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
cat("Running", nperm, "stratified permutations on", n_cores, "cores...\n")
perm_res <- parallel::mclapply(seq_len(nperm), function(p) {
  info_p <- info
  info_p$bin <- perm_list[[p]]
  Xp <- model.matrix(~ batch + inferred_sex + bin, data = info_p)
  res_p <- fit_eisa_joint(Xp)
  calls_p <- sum(res_p$tt$FDR < 0.05)
  max_abs <- max(abs(res_p$tt$logFC))
  if (p %% 25L == 0L) cat(sprintf("Permutation %d/%d: calls=%d\n", p, nperm, calls_p))
  c(calls = calls_p, max_abs = max_abs)
}, mc.cores = n_cores)

stopifnot(all(vapply(perm_res, function(z) is.numeric(z) && length(z) == 2L, logical(1L))))
null_hits <- as.integer(vapply(perm_res, `[[`, numeric(1L), "calls"))
max_abs_dEmI <- as.numeric(vapply(perm_res, `[[`, numeric(1L), "max_abs"))

obs_calls <- sum(R$eisa_FDR < 0.05)
p_global <- (sum(null_hits >= obs_calls) + 1) / (nperm + 1)

fwrite(data.table(permutation = seq_len(nperm), calls = null_hits, max_abs_dEmI = max_abs_dEmI),
       file.path(OUT, "perm_null_hits.tsv"), sep = "\t")
fwrite(data.table(permutation = seq_len(nperm), as.data.table(draw_labels)),
       file.path(OUT, "perm_labels.tsv.gz"), sep = "\t")

summary_obj <- list(
  model = "~ batch + inferred_sex + bin",
  eisa_method = "joint_voom_invvar_weights_with_pmin_Amean_eBayes_trend",
  n_donors = nrow(U),
  n_genes = nrow(R),
  spearman_dI_dE = unname(rho_dI_dE),
  ratio_calls = obs_calls,
  ratio_calls_abs_log2_ge_0_5 = sum(R$eisa_FDR < 0.05 & abs(R$dEmI) >= 0.5),
  median_abs_ratio_change_all_calls = if (obs_calls > 0L) median(abs(R[eisa_FDR < 0.05]$dEmI)) else NA_real_,
  replicated_ratio_calls = sum(REP$eisa_FDR < 0.05, na.rm = TRUE),
  median_abs_ratio_change_replicated_calls = if (sum(REP$eisa_FDR < 0.05, na.rm = TRUE) > 0L) median(abs(REP[eisa_FDR < 0.05]$dEmI)) else NA_real_,
  n_permutations = nperm,
  permutation_seed = perm_seed,
  p_global_call_count = p_global,
  null_median = median(null_hits),
  null_max = max(null_hits),
  first_null_max_abs_dEmI = max_abs_dEmI[1],
  unweighted_v3_calls = if ("eisa_FDR_unweighted_v3" %in% names(R)) sum(R$eisa_FDR_unweighted_v3 < 0.05, na.rm = TRUE) else NA_integer_,
  shared_with_unweighted_v3 = if ("eisa_FDR_unweighted_v3" %in% names(R)) sum(R$eisa_FDR < 0.05 & R$eisa_FDR_unweighted_v3 < 0.05, na.rm = TRUE) else NA_integer_,
  stage_only_calls = if ("eisa_FDR_stage_only" %in% names(R)) sum(R$eisa_FDR_stage_only < 0.05, na.rm = TRUE) else NA_integer_,
  shared_with_stage_only = if ("eisa_FDR_stage_only" %in% names(R)) sum(R$eisa_FDR < 0.05 & R$eisa_FDR_stage_only < 0.05, na.rm = TRUE) else NA_integer_
)
write_json(summary_obj, file.path(OUT, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("DONE EISA joint inverse-variance fix. Calls =", obs_calls, "p_global =", p_global, "\n")
