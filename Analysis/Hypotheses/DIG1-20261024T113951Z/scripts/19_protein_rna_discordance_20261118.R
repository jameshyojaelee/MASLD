#!/usr/bin/env Rscript
# DIG1 job 19 (20261118 candidate, N10/N18):
# Updates over 19_protein_rna_discordance.R:
# 1. Estimate empirical inter-subunit correlation rho_hat_S per gene set directly
#    from covariate-residualized liver protein (P) and bulk RNA (voom$E) matrices,
#    combining them by variance weights w_P^2 = mean(seP^2 / (seP^2 + b^2 seR^2)):
#      rho_hat_S = w_P^2 * rho_{S,P} + (1 - w_P^2) * rho_{S,R},
#    and pass `inter.gene.cor = max(rho_hat_S, 0.01)` into limma::cameraPR() so
#    co-regulated subunit complexes are calibrated with empirical VIF rather than 0.01.
# 2. Add stratified participant bootstrap (B = 500 draws across protein & bulk RNA
#    participants) reporting empirical bootstrap SE, 95% CI, and bootstrap two-sided p
#    on the mean gene-set discordance score `mean_d`.

suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(MASS); library(parallel)})
args <- commandArgs(TRUE)
out <- if (length(args) >= 1L && nzchar(args[1L])) {
  args[1L]
} else {
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/DIG1-20261024T113951Z/out_19_20261118"
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261118)
N_CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
N_BOOT <- 500L

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
GS <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets")
PRIMARY <- c("GOCC_RESPIRATORY_CHAIN_COMPLEX", "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOCC_BLOOD_MICROPARTICLE", "GOCC_LIPID_DROPLET")
SHOW <- c("PLIN1","PLIN2","PLIN3","PLIN4","PLIN5","CIDEB","CIDEC","HSD17B13","PNPLA3","COL1A1","COL1A2","COL3A1","COL6A1","COL6A3","LUM","DCN",
          "MT-CO1","MT-CO2","MT-ND1","MT-ND5","NDUFA9","NDUFS1","SDHA","UQCRC1","COX4I1","ATP5F1A","ALB","APOA1","HP","C3","IGKC","IGHG1",
          "GNMT","MAT1A","CYP2C19","AKR1B10","KRT8","KRT18","KRT7","KRT19")
rd <- function(f) { g <- strsplit(readLines(file.path(GS, f)), "\t"); setNames(lapply(g, function(x) unique(x[-(1:2)])), sapply(g, `[`, 1)) }
cc <- rd("go_cc.gmt"); hm <- rd("hallmark.gmt"); stopifnot(all(PRIMARY %in% names(cc)), length(hm) == 50)

# Protein (PXD051911 single-gene PGs)
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename"); stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & trimws(Genes) != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- trimws(raw$Genes)
x <- x[rowMeans(is.na(x)) <= 0.5, , drop = FALSE]; x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
P <- t(vapply(idx, function(i) if (length(i) == 1L) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(P) <- md$sample_id

pfit <- function(M, mm, coef) {
  fit_obj <- lmFit(M, mm)
  f <- eBayes(fit_obj)
  resid_M <- residuals(fit_obj, M)
  list(
    dt = data.table(gene = rownames(M), logFC = f$coefficients[, coef], se = sqrt(f$s2.post) * f$stdev.unscaled[, coef],
                    t = f$t[, coef], p = f$p.value[, coef]),
    resid = resid_M
  )
}
ok <- complete.cases(md[, .(age, bmi_num, sex)])
P_dis_obj <- pfit(P[, ok], model.matrix(~ masld + batch + age + bmi_num + sex, data = md[ok]), "masld")
okf <- ok & md$masld == 1 & !is.na(md$Fibrosis)
P_fib_obj <- pfit(P[, okf], model.matrix(~ Fibrosis + batch + age + bmi_num + sex, data = md[okf]), "Fibrosis")
Pd <- P_dis_obj$dt; Pf <- P_fib_obj$dt

# RNA (774-unit cohort)
y <- readRDS(file.path(B, "merged_dge.rds")); mb <- readRDS(file.path(B, "meta_matched.rds"))
mb <- mb[match(colnames(y), mb$sample_id)]; stopifnot(identical(mb$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sb <- an$gene_name[match(sub("\\..*$", "", rownames(y)), an$ensembl_base)]
k <- !is.na(sb) & !duplicated(sb); y <- y[k, ]; rownames(y) <- sb[k]
k <- !is.na(mb$sex_final); y <- y[, k]; mb <- mb[k]

rfit <- function(yy, mm, coef) {
  yy <- yy[filterByExpr(yy, design = mm), , keep.lib.sizes = FALSE]; yy <- calcNormFactors(yy)
  v <- voom(yy, mm)
  fit_obj <- lmFit(v, mm)
  f <- eBayes(fit_obj)
  resid_R <- residuals(fit_obj, v$E)
  list(
    dt = data.table(gene = rownames(f), logFC = f$coefficients[, coef], se = sqrt(f$s2.post) * f$stdev.unscaled[, coef],
                    t = f$t[, coef], p = f$p.value[, coef]),
    v = v,
    mm = mm,
    resid = resid_R
  )
}
mb[, disease := as.integer(group_binary == "Disease")]
R_dis_obj <- rfit(y, model.matrix(~ disease + dataset + sex_final, data = mb), "disease")
kf <- mb$disease == 1 & !is.na(mb$fibrosis_stage)
R_fib_obj <- rfit(y[, kf], model.matrix(~ fibrosis_stage + dataset + sex_final, data = mb[kf]), "fibrosis_stage")
Rd <- R_dis_obj$dt; Rf <- R_fib_obj$dt
units <- data.table(contrast = c("disease", "fibrosis_slope"), rna_units = c(ncol(y), sum(kf)), rna_disease_units = c(sum(mb$disease), sum(kf)))

de <- unique(fread(file.path(RUN, "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz")), by = c("lkey", "symbol"))
scF <- de[lkey == "GSE298719:hepatocytes", .(gene = symbol, flex_hep_logFC = logFC)]
scV <- de[lkey == "GSE256398:hepatocytes", .(gene = symbol, v3_hep_logFC = logFC)]

mean_pairwise_cor <- function(mat_g_by_s) {
  if (nrow(mat_g_by_s) < 2L) return(0)
  C <- suppressWarnings(cor(t(mat_g_by_s), use = "pairwise.complete.obs"))
  vals <- C[upper.tri(C)]
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(0)
  mean(vals)
}

run_discordance_vif <- function(P_obj, R_obj, lab, P_sub, meta_P, strata_P, R_sub, meta_R, strata_R) {
  Pt <- P_obj$dt; Rt <- R_obj$dt
  m <- merge(Pt, Rt, by = "gene", suffixes = c("_prot", "_rna"))
  fit <- rlm(logFC_prot ~ logFC_rna, data = m, maxit = 100)
  a <- unname(coef(fit)[1L]); b <- unname(coef(fit)[2L])
  m[, d := (logFC_prot - a - b * logFC_rna) / sqrt(se_prot^2 + b^2 * se_rna^2)]
  m[, w_prot_sq := se_prot^2 / (se_prot^2 + b^2 * se_rna^2)]
  m[, d_p := 2 * pnorm(-abs(d))][, d_bh := p.adjust(d_p, "BH")]
  m <- merge(merge(m, scF, by = "gene", all.x = TRUE), scV, by = "gene", all.x = TRUE)
  d_stats <- setNames(m$d, m$gene)
  resP <- P_obj$resid
  resR <- R_obj$resid

  eval_family <- function(S, fam) {
    S <- lapply(S, intersect, m$gene)
    S <- S[lengths(S) >= 15L & lengths(S) <= 500L]
    res_rows <- rbindlist(lapply(names(S), function(sname) {
      gset <- S[[sname]]
      wp2 <- mean(m[gene %in% gset, w_prot_sq], na.rm = TRUE)
      rho_p <- mean_pairwise_cor(resP[intersect(gset, rownames(resP)), , drop = FALSE])
      rho_r <- mean_pairwise_cor(resR[intersect(gset, rownames(resR)), , drop = FALSE])
      rho_hat <- wp2 * rho_p + (1 - wp2) * rho_r
      rho_used <- max(rho_hat, 0.01)
      idx_list <- ids2indices(list(target = gset), m$gene)
      cam_emp <- as.data.table(cameraPR(d_stats, idx_list, inter.gene.cor = rho_used), keep.rownames = "set")
      cam_01 <- as.data.table(cameraPR(d_stats, idx_list, inter.gene.cor = 0.01), keep.rownames = "set")
      sub_m <- m[gene %in% gset]
      data.table(
        set = sname, family = fam, contrast = lab, NGenes = cam_emp$NGenes[1L],
        Direction = cam_emp$Direction[1L],
        rho_prot = rho_p, rho_rna = rho_r, w_prot_sq = wp2,
        rho_hat_S = rho_hat, inter_gene_cor_used = rho_used,
        vif_empirical = 1 + (cam_emp$NGenes[1L] - 1) * rho_used,
        PValue = cam_emp$PValue[1L],
        PValue_default_0.01 = cam_01$PValue[1L],
        mean_d = mean(sub_m$d, na.rm = TRUE),
        mean_logFC_prot = mean(sub_m$logFC_prot, na.rm = TRUE),
        mean_logFC_rna = mean(sub_m$logFC_rna, na.rm = TRUE),
        mean_flex_hep = mean(sub_m$flex_hep_logFC, na.rm = TRUE),
        mean_v3_hep = mean(sub_m$v3_hep_logFC, na.rm = TRUE),
        frac_d_pos = mean(sub_m$d > 0, na.rm = TRUE)
      )
    }))
    res_rows[, bh_q_within_family := p.adjust(PValue, method = "BH")]
    res_rows
  }
  ccr <- eval_family(cc, "GO_CC")
  hmr <- eval_family(hm, "Hallmark")
  prim <- ccr[set %in% PRIMARY]
  prim[, p_holm_primary := p.adjust(PValue, method = "holm")]

  # Stratified participant bootstrap across protein & RNA participants
  message("[DIG1_19_20261118] Running ", N_BOOT, "-draw stratified participant bootstrap for ", lab, "...")
  p_genes <- intersect(rownames(P_sub), m$gene)
  r_genes <- intersect(rownames(R_obj$v$E), m$gene)
  common_g <- intersect(p_genes, r_genes)
  p_idx <- match(common_g, rownames(P_sub))
  r_idx <- match(common_g, rownames(R_obj$v$E))
  P_common <- P_sub[p_idx, , drop = FALSE]
  E_common <- R_obj$v$E[r_idx, , drop = FALSE]
  W_common <- R_obj$v$weights[r_idx, , drop = FALSE]
  if (lab == "disease") {
    mm_P <- model.matrix(~ masld + batch + age + bmi_num + sex, data = meta_P)
    coef_P <- "masld"
    mm_R <- model.matrix(~ disease + dataset + sex_final, data = meta_R)
    coef_R <- "disease"
  } else {
    mm_P <- model.matrix(~ Fibrosis + batch + age + bmi_num + sex, data = meta_P)
    coef_P <- "Fibrosis"
    mm_R <- model.matrix(~ fibrosis_stage + dataset + sex_final, data = meta_R)
    coef_R <- "fibrosis_stage"
  }
  idx_P_st <- split(seq_len(nrow(meta_P)), strata_P)
  idx_R_st <- split(seq_len(nrow(meta_R)), strata_R)
  all_sets_list <- c(lapply(cc[names(cc) %in% ccr$set], intersect, common_g),
                     lapply(hm[names(hm) %in% hmr$set], intersect, common_g))
  set_names <- names(all_sets_list)
  primary_idx <- match(prim$set, set_names)

  # Weighted OLS per gene using Dirichlet(1,...,1) within strata
  boot_mat <- do.call(rbind, mclapply(seq_len(N_BOOT), function(b_i) {
    set.seed(20261118L + (if (lab == "disease") 0L else 10000L) + b_i)
    wP <- numeric(nrow(meta_P))
    for (ii in idx_P_st) { g <- rexp(length(ii), 1); wP[ii] <- length(ii) * g / sum(g) }
    wR <- numeric(nrow(meta_R))
    for (ii in idx_R_st) { g <- rexp(length(ii), 1); wR[ii] <- length(ii) * g / sum(g) }

    fP <- lmFit(P_common, mm_P, weights = wP)
    # Approximate SE using effective residual variance
    s2P <- rowSums((residuals(fP, P_common)^2) * matrix(wP, nrow = nrow(P_common), ncol = ncol(P_common), byrow = TRUE), na.rm = TRUE) / fP$df.residual
    seP <- sqrt(s2P) * fP$stdev.unscaled[, coef_P]
    lfcP <- fP$coefficients[, coef_P]

    W_boot <- sweep(W_common, 2L, wR, "*")
    fR <- lmFit(E_common, mm_R, weights = W_boot)
    s2R <- rowSums((residuals(fR, E_common)^2) * W_boot, na.rm = TRUE) / fR$df.residual
    seR <- sqrt(s2R) * fR$stdev.unscaled[, coef_R]
    lfcR <- fR$coefficients[, coef_R]

    fit_b <- suppressWarnings(rlm(lfcP ~ lfcR, maxit = 40))
    ab <- coef(fit_b)[1L]; bb <- coef(fit_b)[2L]
    d_b <- (lfcP - ab - bb * lfcR) / sqrt(seP^2 + bb^2 * seR^2)
    names(d_b) <- common_g
    vapply(all_sets_list, function(gs) mean(d_b[gs], na.rm = TRUE), numeric(1))
  }, mc.cores = N_CORES))

  boot_se <- apply(boot_mat, 2L, sd, na.rm = TRUE)
  boot_low <- apply(boot_mat, 2L, quantile, probs = 0.025, na.rm = TRUE)
  boot_high <- apply(boot_mat, 2L, quantile, probs = 0.975, na.rm = TRUE)
  boot_p <- vapply(seq_along(set_names), function(j) {
    v <- boot_mat[, j]
    v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    2 * min(mean(v >= 0), mean(v <= 0), 0.5)
  }, numeric(1))
  boot_dt <- data.table(
    set = set_names, boot_se_mean_d = boot_se,
    boot_ci_low = boot_low, boot_ci_high = boot_high,
    boot_p_two_sided = pmax(boot_p, 1 / (N_BOOT + 1))
  )
  all_sets_dt <- merge(rbind(ccr, hmr), boot_dt, by = "set", all.x = TRUE, sort = FALSE)
  prim_dt <- merge(prim, boot_dt, by = "set", all.x = TRUE, sort = FALSE)
  prim_dt[, p_holm_boot := p.adjust(boot_p_two_sided, method = "holm")]

  list(
    genes = m,
    fit = data.table(contrast = lab, n_genes = nrow(m), intercept = a, slope = b,
                     spearman_logFC = cor(m$logFC_prot, m$logFC_rna, method = "spearman"),
                     n_d_bh05_pos = sum(m$d_bh < 0.05 & m$d > 0), n_d_bh05_neg = sum(m$d_bh < 0.05 & m$d < 0)),
    sets = all_sets_dt,
    prim = prim_dt
  )
}

res <- list(
  disease = run_discordance_vif(
    P_dis_obj, R_dis_obj, "disease",
    P[, ok], md[ok], paste(md$masld[ok], md$batch[ok], sep = "_"),
    y, mb, paste(mb$disease, mb$dataset, sep = "_")
  ),
  fibrosis_slope = run_discordance_vif(
    P_fib_obj, R_fib_obj, "fibrosis_slope",
    P[, okf], md[okf], md$batch[okf],
    y[, kf], mb[kf], mb$dataset[kf]
  )
)

fwrite(merge(rbindlist(lapply(res, `[[`, "fit")), units, by = "contrast"), file.path(out, "fit_summary.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "prim")), file.path(out, "primary_sets_holm.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "sets")), file.path(out, "all_sets_cameraPR.tsv"), sep = "\t")
for (nm in names(res)) fwrite(res[[nm]]$genes[order(-abs(d))], file.path(out, paste0("gene_discordance_", nm, ".tsv.gz")), sep = "\t")
fwrite(rbindlist(lapply(names(res), function(nm) cbind(contrast = nm, res[[nm]]$genes[gene %in% SHOW]))), file.path(out, "named_genes.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(fread(file.path(out, "fit_summary.tsv")))
print(rbindlist(lapply(res, `[[`, "prim"))[, .(contrast, set, NGenes, Direction, rho_hat_S, vif_empirical, PValue, PValue_default_0.01, p_holm_primary, mean_d, boot_ci_low, boot_ci_high, boot_p_two_sided)])
cat("DIG1_19_20261118 DONE\n")
