#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
NF_DIR <- file.path(ROOT, "Analysis/Hypotheses/NUCLEUS-FIX-20261107")
UC_DIR <- file.path(NF_DIR, "unique_counts")
PB_DIR <- file.path(UC_DIR, "pseudobulk")
OUT_DIR <- file.path(NF_DIR, "candidates/sc_20261118")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "de_tables"), recursive = TRUE, showWarnings = FALSE)

donor_meta <- fread(file.path(UC_DIR, "metadata/donor_metadata.tsv"), data.table = TRUE)
lib_meta <- fread(file.path(UC_DIR, "metadata/library_metadata.tsv"), data.table = TRUE)

to_logical <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

donor_meta[, exclude_flag := to_logical(exclude)]
donor_meta[, exclude_stage_flag := to_logical(exclude_stage_analysis)]
donor_meta[, binary_elig_flag := to_logical(binary_eligible)]
lib_meta[, prep_elig_flag := preparation_method %in% c("nuclei", "cd45_negative", "unsorted")]

# Map library_metadata lineage count columns to pseudobulk file prefixes
pb_files <- sort(list.files(PB_DIR, pattern = "_pseudobulk\\.csv$", full.names = TRUE))
stopifnot(length(pb_files) == 16L)

norm_ct_key <- function(s) gsub("[^A-Za-z0-9+]", "_", s)
n_cols <- grep("^n_", names(lib_meta), value = TRUE)
n_cols <- setdiff(n_cols, c("n_cells", "n_source_records", "n_donor_cells"))
col_map <- setNames(n_cols, norm_ct_key(sub("^n_", "", n_cols)))

run_limma_voom <- function(counts_mat, cond_vec, ds_vec = NULL, min_count = 5, min_samples = 3) {
  keep_genes <- rowSums(counts_mat >= min_count) >= min_samples
  counts_sub <- counts_mat[keep_genes, , drop = FALSE]
  if (nrow(counts_sub) < 10L) return(NULL)

  cond_fac <- factor(cond_vec, levels = c("Healthy", "MASLD"))
  if (any(is.na(cond_fac)) || any(table(cond_fac) < 2L)) return(NULL)

  if (!is.null(ds_vec) && length(unique(ds_vec)) > 1L) {
    ds_fac <- factor(ds_vec)
    design <- model.matrix(~ ds_fac + cond_fac)
  } else {
    design <- model.matrix(~ cond_fac)
  }
  coef_name <- "cond_facMASLD"
  if (!(coef_name %in% colnames(design))) return(NULL)

  dge <- DGEList(counts = counts_sub)
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit2 <- eBayes(fit)
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
  dt <- as.data.table(tt)
  dt[, gene := rownames(tt)]
  j <- match(coef_name, colnames(fit2$stdev.unscaled))
  dt[, se := as.numeric(fit2$stdev.unscaled[rownames(tt), j] * sqrt(fit2$s2.post))]
  dt[, z_signed := ifelse(t == 0, 0, sign(t) * qnorm(pmax(P.Value, 1e-300) / 2, lower.tail = FALSE))]
  dt <- dt[, .(gene, logFC, AveExpr, t, se, z_signed, pvalue = P.Value, padj = adj.P.Val, B)]
  dt
}

meta_combine_per_cohort <- function(cohort_dts, cohort_info) {
  # cohort_dts: named list of data.tables by dataset; cohort_info: data.table(dataset, n0, n1)
  if (length(cohort_dts) == 0L) return(NULL)
  all_dt <- rbindlist(lapply(names(cohort_dts), function(ds) {
    d <- copy(cohort_dts[[ds]])
    ci <- cohort_info[dataset == ds]
    d[, `:=`(dataset = ds, n0 = ci$n0[1], n1 = ci$n1[1])]
    d
  }))
  all_dt[, w_df := sqrt((n0 * n1) / (n0 + n1))]
  all_dt[, w_eq := 1.0]
  all_dt[, w_iv := ifelse(is.finite(se) & se > 0, 1.0 / (se^2), 0.0)]

  agg <- all_dt[, {
    k <- .N
    ds_list <- paste(sort(dataset), collapse = ";")
    z_df <- sum(w_df * z_signed) / sqrt(sum(w_df^2))
    p_df <- 2 * pnorm(-abs(z_df))
    z_eq <- sum(w_eq * z_signed) / sqrt(sum(w_eq^2))
    p_eq <- 2 * pnorm(-abs(z_eq))
    if (sum(w_iv) > 0) {
      b_iv <- sum(w_iv * logFC) / sum(w_iv)
      se_iv <- sqrt(1.0 / sum(w_iv))
      z_iv <- b_iv / se_iv
      p_iv <- 2 * pnorm(-abs(z_iv))
      Q_iv <- if (k >= 2L) sum(w_iv * (logFC - b_iv)^2) else NA_real_
      I2_iv <- if (k >= 2L && is.finite(Q_iv) && Q_iv > 0) max(0, (Q_iv - (k - 1)) / Q_iv) else if (k >= 2L) 0 else NA_real_
    } else {
      b_iv <- NA_real_; se_iv <- NA_real_; z_iv <- NA_real_; p_iv <- NA_real_
      Q_iv <- NA_real_; I2_iv <- NA_real_
    }
    same_sign <- all(sign(logFC) == sign(logFC[1]))
    mean_logFC <- mean(logFC)
    .(n_cohorts = k, cohorts = ds_list,
      total_n0 = sum(n0), total_n1 = sum(n1),
      mean_logFC = mean_logFC, same_sign = same_sign,
      z_stouffer_df = z_df, p_stouffer_df = p_df,
      z_stouffer_equal = z_eq, p_stouffer_equal = p_eq,
      beta_ivw = b_iv, se_ivw = se_iv, z_ivw = z_iv, p_ivw = p_iv,
      cochran_Q = Q_iv, I2 = I2_iv)
  }, by = gene]

  # Adjust only across multi-cohort (>=2) when available, else all tested
  agg[, padj_stouffer_df := NA_real_]
  agg[, padj_stouffer_equal := NA_real_]
  agg[, padj_ivw := NA_real_]
  idx2 <- which(agg$n_cohorts >= 2L)
  if (length(idx2) > 0L) {
    agg[idx2, padj_stouffer_df := p.adjust(p_stouffer_df, method = "BH")]
    agg[idx2, padj_stouffer_equal := p.adjust(p_stouffer_equal, method = "BH")]
    agg[idx2, padj_ivw := p.adjust(p_ivw, method = "BH")]
  }
  agg[, padj_stouffer_df_anyk := p.adjust(p_stouffer_df, method = "BH")]
  agg[, padj_stouffer_equal_anyk := p.adjust(p_stouffer_equal, method = "BH")]
  agg[, padj_ivw_anyk := p.adjust(p_ivw, method = "BH")]
  agg[]
}

collapse_lineage_to_donor <- function(pb_path, ct_name) {
  raw <- fread(pb_path, data.table = FALSE)
  genes <- as.character(raw[[1]])
  mat <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  stopifnot(all(abs(mat - round(mat)) < 1e-6))
  mat <- round(mat)
  if (anyDuplicated(genes) > 0L) {
    mat <- rowsum(mat, group = genes, reorder = FALSE)
  } else {
    rownames(mat) <- genes
  }

  caps <- colnames(mat)
  lm <- lib_meta[match(caps, sample)]
  stopifnot(!any(is.na(lm$donor)))

  ncell_col <- col_map[[norm_ct_key(ct_name)]]
  stopifnot(!is.null(ncell_col), ncell_col %in% names(lib_meta))
  lm[, ct_cells := as.integer(get(ncell_col))]

  # Filter capture-level prep_eligible
  keep_cap <- lm$prep_elig_flag
  mat_elig <- mat[, keep_cap, drop = FALSE]
  lm_elig <- lm[keep_cap]

  # Sum counts & cell counts per donor
  donors_u <- unique(lm_elig$donor)
  donor_mat <- sapply(donors_u, function(d) {
    cols <- which(lm_elig$donor == d)
    if (length(cols) == 1L) mat_elig[, cols] else rowSums(mat_elig[, cols, drop = FALSE])
  })
  if (is.null(dim(donor_mat))) {
    donor_mat <- matrix(donor_mat, nrow = nrow(mat_elig), ncol = length(donors_u))
  }
  rownames(donor_mat) <- rownames(mat_elig)
  colnames(donor_mat) <- donors_u

  donor_cells <- lm_elig[, .(n_cells = sum(ct_cells),
                             n_captures = .N,
                             dataset = dataset[1],
                             prep_methods = paste(sort(unique(preparation_method)), collapse = "+")),
                         by = donor]
  dm <- merge(donor_cells, donor_meta[, .(donor, condition_harmonized, disease_stage_coarse,
                                          exclude_flag, exclude_stage_flag, binary_elig_flag)],
              by = "donor", all.x = TRUE, sort = FALSE)
  list(counts = donor_mat, meta = dm)
}

roster_rows <- list()
summary_rows <- list()
arm_cmp_rows <- list()
sig_top100_rows <- list()

# Also store key DE tables for V5 correlation check
de_cache <- list()

for (pb_path in pb_files) {
  ct_name <- sub("_pseudobulk\\.csv$", "", basename(pb_path))
  col_obj <- collapse_lineage_to_donor(pb_path, ct_name)
  dmat <- col_obj$counts
  dmeta <- col_obj$meta

  # Arm definitions at donor n_cells >= 50 (primary)
  arms <- list(
    excl_GSE136103 = dmeta[!exclude_flag & binary_elig_flag & condition_harmonized %in% c("Healthy", "MASLD") & n_cells >= 50L],
    incl_GSE136103 = dmeta[dataset != "Liver_Atlas" & binary_elig_flag & condition_harmonized %in% c("Healthy", "MASLD") & n_cells >= 50L]
  )

  # Also track n_cells >= 30 donor counts for completeness
  dmeta_30 <- dmeta[!exclude_flag & binary_elig_flag & condition_harmonized %in% c("Healthy", "MASLD") & n_cells >= 30L]

  arm_results <- list()

  for (arm_name in names(arms)) {
    subm <- copy(arms[[arm_name]])
    subm[, arm := arm_name]
    subm[, lineage := ct_name]
    roster_rows[[paste0(ct_name, "_", arm_name)]] <- subm

    # Identify per-dataset 2-arm testability
    ds_tab <- subm[, .(n_healthy = sum(condition_harmonized == "Healthy"),
                       n_masld = sum(condition_harmonized == "MASLD"),
                       n_donors = .N), by = dataset]
    two_arm_ds <- ds_tab[n_healthy >= 2L & n_masld >= 2L, dataset]

    # 1) Pooled limma-voom (all arm donors >= 50 cells, ~ dataset + condition)
    pooled_all <- if (sum(subm$condition_harmonized == "Healthy") >= 2L &&
                      sum(subm$condition_harmonized == "MASLD") >= 2L) {
      run_limma_voom(dmat[, subm$donor, drop = FALSE], subm$condition_harmonized, subm$dataset)
    } else NULL

    # 2) Pooled limma-voom restricted to two-arm cohorts
    subm_2arm <- subm[dataset %in% two_arm_ds]
    pooled_2arm <- if (nrow(subm_2arm) >= 4L &&
                       sum(subm_2arm$condition_harmonized == "Healthy") >= 2L &&
                       sum(subm_2arm$condition_harmonized == "MASLD") >= 2L) {
      run_limma_voom(dmat[, subm_2arm$donor, drop = FALSE], subm_2arm$condition_harmonized, subm_2arm$dataset)
    } else NULL

    # 3) Per-cohort limma-voom + Stouffer / IVW meta
    cohort_dts <- list()
    for (ds in two_arm_ds) {
      sd <- subm[dataset == ds]
      fit_ds <- run_limma_voom(dmat[, sd$donor, drop = FALSE], sd$condition_harmonized, NULL)
      if (!is.null(fit_ds)) cohort_dts[[ds]] <- fit_ds
    }
    meta_dt <- meta_combine_per_cohort(cohort_dts, ds_tab[, .(dataset, n0 = n_healthy, n1 = n_masld)])

    # Write DE files
    if (!is.null(pooled_all)) {
      fwrite(pooled_all, file.path(OUT_DIR, "de_tables", sprintf("%s_%s_pooled.tsv.gz", ct_name, arm_name)), sep = "\t")
    }
    if (!is.null(pooled_2arm)) {
      fwrite(pooled_2arm, file.path(OUT_DIR, "de_tables", sprintf("%s_%s_pooled_twoarm.tsv.gz", ct_name, arm_name)), sep = "\t")
    }
    if (!is.null(meta_dt)) {
      fwrite(meta_dt, file.path(OUT_DIR, "de_tables", sprintf("%s_%s_stouffer_meta.tsv.gz", ct_name, arm_name)), sep = "\t")
    }

    arm_results[[arm_name]] <- list(
      subm = subm, ds_tab = ds_tab, two_arm_ds = two_arm_ds,
      pooled_all = pooled_all, pooled_2arm = pooled_2arm, meta = meta_dt, cohort_dts = cohort_dts
    )

    # Record summary row
    n_h <- sum(subm$condition_harmonized == "Healthy")
    n_m <- sum(subm$condition_harmonized == "MASLD")
    n_ds <- uniqueN(subm$dataset)
    cohorts_str <- paste(ds_tab[, sprintf("%s(%dH/%dM)", dataset, n_healthy, n_masld)], collapse = ";")
    summary_rows[[paste0(ct_name, "_", arm_name)]] <- data.table(
      lineage = ct_name,
      arm = arm_name,
      n_donors_ge50 = nrow(subm),
      n_healthy_ge50 = n_h,
      n_masld_ge50 = n_m,
      n_donors_ge30_excl = if (arm_name == "excl_GSE136103") nrow(dmeta_30) else NA_integer_,
      n_datasets = n_ds,
      n_two_arm_datasets = length(two_arm_ds),
      two_arm_datasets = paste(two_arm_ds, collapse = ";"),
      cohorts_breakdown = cohorts_str,
      pooled_testable = !is.null(pooled_all),
      pooled_n_genes = if (!is.null(pooled_all)) nrow(pooled_all) else 0L,
      pooled_n_fdr05 = if (!is.null(pooled_all)) sum(pooled_all$padj < 0.05) else NA_integer_,
      pooled_n_up_fdr05 = if (!is.null(pooled_all)) sum(pooled_all$padj < 0.05 & pooled_all$logFC > 0) else NA_integer_,
      pooled_n_down_fdr05 = if (!is.null(pooled_all)) sum(pooled_all$padj < 0.05 & pooled_all$logFC < 0) else NA_integer_,
      pooled2arm_n_fdr05 = if (!is.null(pooled_2arm)) sum(pooled_2arm$padj < 0.05) else NA_integer_,
      meta_testable_ge2cohort = !is.null(meta_dt) && sum(meta_dt$n_cohorts >= 2L) > 0L,
      meta_n_genes_ge2cohort = if (!is.null(meta_dt)) sum(meta_dt$n_cohorts >= 2L) else 0L,
      stouffer_df_n_fdr05 = if (!is.null(meta_dt)) sum(meta_dt$n_cohorts >= 2L & meta_dt$padj_stouffer_df < 0.05, na.rm = TRUE) else NA_integer_,
      stouffer_equal_n_fdr05 = if (!is.null(meta_dt)) sum(meta_dt$n_cohorts >= 2L & meta_dt$padj_stouffer_equal < 0.05, na.rm = TRUE) else NA_integer_,
      ivw_n_fdr05 = if (!is.null(meta_dt)) sum(meta_dt$n_cohorts >= 2L & meta_dt$padj_ivw < 0.05, na.rm = TRUE) else NA_integer_
    )

    # Extract top-100 UP/DOWN disease signature table for primary excl_GSE136103
    if (arm_name == "excl_GSE136103" && !is.null(pooled_all)) {
      up_dt <- pooled_all[logFC > 0][order(padj, -abs(t))]
      dn_dt <- pooled_all[logFC < 0][order(padj, -abs(t))]
      up_100 <- head(up_dt, 100L)
      dn_100 <- head(dn_dt, 100L)
      if (nrow(up_100) > 0L) {
        up_100[, `:=`(lineage = ct_name, direction = "UP", rank = seq_len(.N), within_fdr20 = padj <= 0.2, within_fdr05 = padj < 0.05)]
        sig_top100_rows[[paste0(ct_name, "_UP")]] <- up_100
      }
      if (nrow(dn_100) > 0L) {
        dn_100[, `:=`(lineage = ct_name, direction = "DOWN", rank = seq_len(.N), within_fdr20 = padj <= 0.2, within_fdr05 = padj < 0.05)]
        sig_top100_rows[[paste0(ct_name, "_DOWN")]] <- dn_100
      }
    }
  }

  de_cache[[ct_name]] <- arm_results

  # Compare excl_GSE136103 vs incl_GSE136103
  ex <- arm_results$excl_GSE136103
  inc <- arm_results$incl_GSE136103

  cmp_pair <- function(dt1, col1, padj1, dt2, col2, padj2, label1, label2) {
    if (is.null(dt1) || is.null(dt2)) return(NULL)
    m <- merge(dt1[, .(gene, v1 = get(col1), q1 = get(padj1))],
               dt2[, .(gene, v2 = get(col2), q2 = get(padj2))], by = "gene")
    m <- m[is.finite(v1) & is.finite(v2)]
    if (nrow(m) < 50L) return(NULL)
    t1_up <- head(m[order(-v1), gene], 100L)
    t2_up <- head(m[order(-v2), gene], 100L)
    t1_dn <- head(m[order(v1), gene], 100L)
    t2_dn <- head(m[order(v2), gene], 100L)
    jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
    s1 <- which(m$q1 < 0.05)
    s2 <- which(m$q2 < 0.05)
    both <- intersect(s1, s2)
    data.table(
      lineage = ct_name,
      comparison = paste0(label1, "_vs_", label2),
      n_common_genes = nrow(m),
      pearson_r = cor(m$v1, m$v2, method = "pearson"),
      spearman_rho = cor(m$v1, m$v2, method = "spearman"),
      sign_concordance_all = mean(sign(m$v1) == sign(m$v2)),
      n_sig1_fdr05 = length(s1),
      n_sig2_fdr05 = length(s2),
      n_sig_overlap = length(both),
      sig_jaccard_fdr05 = if (length(union(s1, s2)) > 0L) length(both) / length(union(s1, s2)) else NA_real_,
      sig_overlap_same_sign = if (length(both) > 0L) mean(sign(m$v1[both]) == sign(m$v2[both])) else NA_real_,
      top100_up_jaccard = jacc(t1_up, t2_up),
      top100_down_jaccard = jacc(t1_dn, t2_dn)
    )
  }

  sub_ge2 <- function(dt) {
    if (is.null(dt) || nrow(dt) == 0L) return(NULL)
    res <- dt[n_cohorts >= 2L]
    if (nrow(res) == 0L) NULL else res
  }

  cmps <- list(
    cmp_pair(ex$pooled_all, "t", "padj", inc$pooled_all, "t", "padj", "excl136103_pooled", "incl136103_pooled"),
    cmp_pair(ex$pooled_2arm, "t", "padj", inc$pooled_2arm, "t", "padj", "excl136103_pooled2arm", "incl136103_pooled2arm"),
    cmp_pair(sub_ge2(ex$meta), "z_stouffer_df", "padj_stouffer_df",
             sub_ge2(inc$meta), "z_stouffer_df", "padj_stouffer_df",
             "excl136103_stouffer_df", "incl136103_stouffer_df"),
    cmp_pair(sub_ge2(ex$meta), "z_stouffer_equal", "padj_stouffer_equal",
             sub_ge2(inc$meta), "z_stouffer_equal", "padj_stouffer_equal",
             "excl136103_stouffer_equal", "incl136103_stouffer_equal"),
    cmp_pair(sub_ge2(ex$meta), "z_ivw", "padj_ivw",
             sub_ge2(inc$meta), "z_ivw", "padj_ivw",
             "excl136103_ivw", "incl136103_ivw"),
    cmp_pair(sub_ge2(ex$meta), "z_stouffer_df", "padj_stouffer_df",
             sub_ge2(ex$meta), "z_stouffer_equal", "padj_stouffer_equal",
             "excl136103_stouffer_df", "excl136103_stouffer_equal"),
    cmp_pair(ex$pooled_all, "t", "padj",
             sub_ge2(ex$meta), "z_stouffer_df", "padj_stouffer_df",
             "excl136103_pooled", "excl136103_stouffer_df")
  )
  cmps <- Filter(Negate(is.null), cmps)
  if (length(cmps) > 0L) arm_cmp_rows[[ct_name]] <- rbindlist(cmps)
}

fwrite(rbindlist(roster_rows, fill = TRUE), file.path(OUT_DIR, "donor_lineage_eligibility_roster.tsv"), sep = "\t")
fwrite(rbindlist(summary_rows, fill = TRUE), file.path(OUT_DIR, "pseudobulk_de_and_stouffer_summary.tsv"), sep = "\t")
fwrite(rbindlist(arm_cmp_rows, fill = TRUE), file.path(OUT_DIR, "gse136103_exclusion_and_weighting_comparisons.tsv"), sep = "\t")
fwrite(rbindlist(sig_top100_rows, fill = TRUE), file.path(OUT_DIR, "disease_signature_top100_up_down.tsv"), sep = "\t")

# ==============================================================================
# Held-out GSE298719 (V5) Lineage Pseudobulk DE Benchmark (excl vs incl GSE136103)
# ==============================================================================
v5_dir <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z")
v5_dq <- fread(file.path(v5_dir, "outputs/donor_qc.tsv"))
v5_dq <- v5_dq[to_logical(included)]
v5_cells <- fread(file.path(v5_dir, "outputs/job3/cells_per_donor_set.tsv"))
v5_feat <- fread(cmd = sprintf("zcat %s", shQuote(file.path(v5_dir, "downloads_untrusted/matrices/GSM9859944_FRPBiopsy_01_Male_Healthy_features.tsv.gz"))),
                 header = FALSE, col.names = c("gene_id_ver", "symbol", "ftype"))
v5_feat[, gene_id := sub("\\..*$", "", gene_id_ver)]
ens2sym <- setNames(as.character(v5_feat$symbol), as.character(v5_feat$gene_id))

v5_map <- list(
  Hepatocytes = "hep",
  Cholangiocytes = "cholangiocytes",
  Fibroblasts = "fibroblasts",
  Macrophages = "macrophages",
  T_cells = "tcells"
)

v5_rows <- list()
for (ct_name in names(v5_map)) {
  v5_set <- v5_map[[ct_name]]
  pb_v5 <- fread(cmd = sprintf("zcat %s", shQuote(file.path(v5_dir, "outputs/job3", sprintf("pseudobulk_%s.tsv.gz", v5_set)))),
                 data.table = FALSE)
  g_ens <- as.character(pb_v5[[1]])
  syms <- ens2sym[g_ens]
  keep_sym <- !is.na(syms) & nzchar(syms)
  mat_v5 <- as.matrix(pb_v5[keep_sym, -1, drop = FALSE])
  storage.mode(mat_v5) <- "numeric"
  syms_keep <- syms[keep_sym]
  # Collapse duplicate symbols by summing (matching V5 01_qc_annotate.py)
  mat_sym <- rowsum(mat_v5, group = syms_keep, reorder = FALSE)

  cset <- v5_cells[set == v5_set & n_cells >= 30L]
  meta_v5 <- merge(v5_dq, cset, by = "gsm")
  meta_v5 <- meta_v5[order(match(gsm, colnames(mat_sym)))]
  cols_use <- meta_v5$gsm
  mat_use <- mat_sym[, cols_use, drop = FALSE]

  keep_g <- rowSums(mat_use >= 5) >= 3L
  mat_use <- mat_use[keep_g, , drop = FALSE]
  dis_f <- factor(ifelse(meta_v5$disease == 1L, "MASLD", "Healthy"), levels = c("Healthy", "MASLD"))
  sex_f <- factor(meta_v5$sex_f)
  des_v5 <- model.matrix(~ dis_f + sex_f)
  dge_v5 <- calcNormFactors(DGEList(counts = mat_use), method = "TMM")
  v_v5 <- voom(dge_v5, des_v5, plot = FALSE)
  fit_v5 <- eBayes(lmFit(v_v5, des_v5))
  tt_v5 <- topTable(fit_v5, coef = "dis_fMASLD", number = Inf, sort.by = "none")
  v5_de <- data.table(gene = rownames(tt_v5), logFC_v5 = tt_v5$logFC, t_v5 = tt_v5$t,
                      pvalue_v5 = tt_v5$P.Value, padj_v5 = tt_v5$adj.P.Val)
  fwrite(v5_de, file.path(OUT_DIR, "de_tables", sprintf("%s_V5_GSE298719_de.tsv.gz", ct_name)), sep = "\t")

  # Evaluate agreement across discovery methods
  ex <- de_cache[[ct_name]]$excl_GSE136103
  inc <- de_cache[[ct_name]]$incl_GSE136103
  sg2 <- function(dt) if (is.null(dt) || nrow(dt) == 0L) NULL else dt[n_cohorts >= 2L]
  cands <- list(
    list(arm = "excl_GSE136103_pooled", dt = ex$pooled_all, col = "t", padj = "padj"),
    list(arm = "incl_GSE136103_pooled", dt = inc$pooled_all, col = "t", padj = "padj"),
    list(arm = "excl_GSE136103_pooled2arm", dt = ex$pooled_2arm, col = "t", padj = "padj"),
    list(arm = "incl_GSE136103_pooled2arm", dt = inc$pooled_2arm, col = "t", padj = "padj"),
    list(arm = "excl_GSE136103_stouffer_df", dt = sg2(ex$meta), col = "z_stouffer_df", padj = "padj_stouffer_df"),
    list(arm = "incl_GSE136103_stouffer_df", dt = sg2(inc$meta), col = "z_stouffer_df", padj = "padj_stouffer_df"),
    list(arm = "excl_GSE136103_stouffer_equal", dt = sg2(ex$meta), col = "z_stouffer_equal", padj = "padj_stouffer_equal"),
    list(arm = "incl_GSE136103_stouffer_equal", dt = sg2(inc$meta), col = "z_stouffer_equal", padj = "padj_stouffer_equal"),
    list(arm = "excl_GSE136103_ivw", dt = sg2(ex$meta), col = "z_ivw", padj = "padj_ivw"),
    list(arm = "incl_GSE136103_ivw", dt = sg2(inc$meta), col = "z_ivw", padj = "padj_ivw")
  )
  for (cc in cands) {
    if (is.null(cc$dt) || nrow(cc$dt) == 0L) next
    m <- merge(cc$dt[, .(gene, stat = get(cc$col), padj_disc = get(cc$padj))], v5_de, by = "gene")
    m <- m[is.finite(stat) & is.finite(t_v5)]
    if (nrow(m) < 50L) next
    sig <- m[padj_disc < 0.05]
    v5_rows[[paste0(ct_name, "_", cc$arm)]] <- data.table(
      lineage = ct_name,
      discovery_arm = cc$arm,
      v5_n_donors = nrow(meta_v5),
      v5_n_healthy = sum(meta_v5$disease == 0L),
      v5_n_masld = sum(meta_v5$disease == 1L),
      n_common_genes = nrow(m),
      pearson_r_vs_v5 = cor(m$stat, m$t_v5, method = "pearson"),
      spearman_rho_vs_v5 = cor(m$stat, m$t_v5, method = "spearman"),
      sign_agree_all = mean(sign(m$stat) == sign(m$t_v5)),
      n_disc_fdr05 = nrow(sig),
      sign_agree_disc_fdr05 = if (nrow(sig) > 0L) mean(sign(sig$stat) == sign(sig$t_v5)) else NA_real_,
      v5_nom_p05_among_disc_fdr05 = if (nrow(sig) > 0L) mean(sig$pvalue_v5 < 0.05 & sign(sig$stat) == sign(sig$t_v5)) else NA_real_,
      top100_up_v5_spearman = {
        u <- head(m[order(-stat)], 100L); if (nrow(u) >= 20L) mean(u$t_v5 > 0) else NA_real_
      },
      top100_dn_v5_frac_neg = {
        d <- head(m[order(stat)], 100L); if (nrow(d) >= 20L) mean(d$t_v5 < 0) else NA_real_
      }
    )
  }
}

fwrite(rbindlist(v5_rows, fill = TRUE), file.path(OUT_DIR, "v5_gse298719_agreement_benchmark.tsv"), sep = "\t")
cat("DONE R DE + Stouffer + V5 benchmark\n")
