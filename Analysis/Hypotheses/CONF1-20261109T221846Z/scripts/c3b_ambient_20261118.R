# CONF1 C3b Redesigned (N01, 2026-11-18):
# 1. Reads deduplicated NUCLEUS-FIX/unique_counts/pseudobulk/ across ALL 16 lineages (Hepatocytes + 15 non-hep lineages),
#    plus V5 and V3 pseudobulks.
# 2. Subtracts all 15 non-hepatocyte lineages PLUS an explicit pooled total-non-hepatocyte UMI soup profile per donor
#    so rare lineages (<20 cells in a biopsy) are still subtracted.
# 3. Restricts soup markers to gene_biotype == "protein_coding" with >= 10 counts in donor's non-hepatocyte pool
#    (never clamps 0-marker donors to rho = 1.0; sets rho = 0 instead).
# 4. Removes all cross-donor / group-mean imputation and runs correct_cohort() INSIDE each of 1,000 permutations.
# 5. Requires pmax(share_ctrl, share_disease) < 0.50 AND hepatocyte predominance (sc_hep_share >= 0.50 or hep_cpm >= max_nonhep_cpm).
# 6. Verifies 0 non-hepatocyte best-lineage genes (ACTA2, EPCAM, CLEC4M) remain in within_hepatocyte_corrected.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
OUT <- args[1]
if (dir.exists(OUT)) unlink(OUT, recursive = TRUE)
dir.create(OUT, recursive = TRUE)
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HYP <- file.path(ROOT, "Analysis/Hypotheses")
RUN <- file.path(HYP, "CONF1-20261109T221846Z")
V5 <- file.path(HYP, "V5-20261021T122723Z/outputs")
V3 <- file.path(HYP, "V3-20261008T015239Z/pseudobulk")
NF <- file.path(HYP, "NUCLEUS-FIX-20261107/unique_counts")

ATLAS_ALL_LIN <- c(
  Hepatocytes = "Hepatocytes",
  Cholangiocytes = "Cholangiocytes",
  Fibroblasts = "Fibroblasts",
  Endothelial = "Endothelial_cells",
  Macrophages = "Macrophages",
  Mono_derived = "Mono+mono_derived_cells",
  T_cells = "T_cells",
  Circulating_NK_NKT = "Circulating_NK_NKT",
  Resident_NK = "Resident_NK",
  B_cells = "B_cells",
  Plasma_cells = "Plasma_cells",
  Basophils = "Basophils",
  Neutrophils = "Neutrophils",
  cDC1s = "cDC1s",
  cDC2s = "cDC2s",
  pDCs = "pDCs"
)
V5_LIN <- c(Hepatocytes = "hep", Cholangiocytes = "cholangiocytes", Fibroblasts = "fibroblasts",
            Macrophages = "macrophages", T_cells = "tcells")
V3_LIN <- c(Hepatocytes = "hepatocyte", Cholangiocytes = "cholangiocyte", Fibroblasts = "stellate_fibroblast")

IN <- c(
  c3_partition = file.path(RUN, "out_c3_partition_20261118/gene_partition.tsv.gz"),
  catalog = file.path(ROOT, "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/catalog_annotated_gene_table.tsv.gz"),
  v5_donors = file.path(V5, "donor_qc.tsv"), v5_cells = file.path(V5, "job3/cells_per_donor_set.tsv"),
  v3_donors = file.path(V3, "donor_metadata.tsv"), v3_cells = file.path(V3, "nuclei_per_donor_lineage.tsv"),
  atlas_donors = file.path(NF, "metadata/donor_metadata.tsv"), atlas_libraries = file.path(NF, "metadata/library_metadata.tsv"),
  atlas_sex = file.path(HYP, "RV1-20261031T155620Z/out-20261102T130905Z/donors.tsv"),
  annotation = file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  setNames(file.path(V5, sprintf("job3/pseudobulk_%s.tsv.gz", V5_LIN)), paste0("v5_", names(V5_LIN))),
  setNames(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", V3_LIN)), paste0("v3_", names(V3_LIN))),
  setNames(file.path(NF, sprintf("pseudobulk/%s_pseudobulk.csv", ATLAS_ALL_LIN)), paste0("atlas_", names(ATLAS_ALL_LIN)))
)
stopifnot(all(file.exists(IN)))

SEED <- 20261109L
N_PERM <- 1000L
MIN_CELLS <- 30L
MIN_GROUP <- 3L
MIN_POOL_MARKER_COUNTS <- 10L
N_MARK <- 30L
MARK_CPM <- 50
AMB_MAX <- 0.50
CONTROLS <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "ACTA2", "LYZ", "STAB2")

wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")
wr(data.table(input = names(IN), path = IN, sha256 = vapply(IN, function(f)
  sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE)), "")), "inputs_sha256.tsv")

P3 <- fread(IN["c3_partition"])[replicated == TRUE]
dup_sym <- P3$symbol[duplicated(P3$symbol) | is.na(P3$symbol) | P3$symbol == ""]
TARGET <- setdiff(P3$symbol, dup_sym)
BULK_SIGN <- setNames(sign(P3$bulk_log2fc[match(TARGET, P3$symbol)]), TARGET)

annot <- fread(IN["annotation"], select = c("ensembl_base", "gene_name", "gene_biotype"))
protein_coding_syms <- unique(annot[gene_biotype == "protein_coding" & !is.na(gene_name) & nzchar(gene_name), gene_name])

by_symbol <- function(counts, symbol) {
  ok <- !is.na(symbol) & symbol != "" & !(symbol %in% symbol[duplicated(symbol)])
  counts <- counts[ok, , drop = FALSE]; rownames(counts) <- symbol[ok]; counts
}
disease_of <- function(stage_ord) fifelse(stage_ord == 0, 0, fifelse(stage_ord >= 2, 1, NA_real_))

load_v5 <- function() {
  dq <- fread(IN["v5_donors"])[included == TRUE]; cells <- fread(IN["v5_cells"])
  meta <- data.table(donor = dq$gsm, female = as.integer(dq$sex_f), disease = as.numeric(dq$disease))
  lapply(V5_LIN, function(s) {
    pb <- fread(file.path(V5, sprintf("job3/pseudobulk_%s.tsv.gz", s)))
    list(meta = copy(meta)[, n_cells := cells[set == s][match(meta$donor, gsm), n_cells]],
         counts = by_symbol(as.matrix(pb[, meta$donor, with = FALSE]),
                            annot$gene_name[match(sub("\\..*$", "", pb$gene_id), annot$ensembl_base)]))
  })
}
load_v3 <- function() {
  m <- fread(IN["v3_donors"]); n3 <- fread(IN["v3_cells"])
  meta <- data.table(donor = m$donor, female = as.integer(m$sex == "female"), disease = disease_of(m$stage_ordinal))
  lapply(V3_LIN, function(s) {
    pb <- fread(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", s)))
    list(meta = copy(meta)[, n_cells := n3[[s]][match(meta$donor, n3$donor)]],
         counts = by_symbol(as.matrix(pb[, meta$donor, with = FALSE]), pb$symbol))
  })
}
load_atlas <- function(ds) {
  dm <- fread(IN["atlas_donors"])[dataset == ds & toupper(as.character(exclude)) == "FALSE"]
  lib <- fread(IN["atlas_libraries"])[donor %in% dm$donor]
  sex <- fread(IN["atlas_sex"], select = c("donor", "sex"))
  stopifnot(all(dm$donor %in% sex$donor))
  meta <- data.table(donor = dm$donor, female = as.integer(sex$sex[match(dm$donor, sex$donor)] == "Female"),
                     disease = disease_of(dm$stage_ordinal))
  lapply(ATLAS_ALL_LIN, function(s) {
    pb <- fread(file.path(NF, sprintf("pseudobulk/%s_pseudobulk.csv", s)))
    ncol_name <- paste0("n_", gsub("Circulating_NK_NKT", "Circulating_NK/NKT", s))
    stopifnot(ncol_name %in% names(lib))
    libs <- lapply(setNames(dm$donor, dm$donor), function(d) intersect(lib[donor == d, sample], names(pb)))
    has <- lengths(libs) > 0
    if (!any(has)) {
      return(list(meta = meta[0], counts = matrix(0, nrow = 0, ncol = 0)))
    }
    counts <- vapply(libs[has], function(l) rowSums(as.matrix(pb[, l, with = FALSE])), numeric(nrow(pb)))
    if (is.null(dim(counts))) {
      counts <- matrix(counts, ncol = sum(has))
      rownames(counts) <- NULL
    }
    colnames(counts) <- names(libs[has])
    mm <- copy(meta)[has][, n_cells := vapply(donor, function(d) sum(lib[donor == d][[ncol_name]]), numeric(1))]
    list(meta = mm, counts = by_symbol(counts, pb$gene))
  })
}

CO <- list(GSE298719 = load_v5(), GSE256398 = load_v3(), GSE202379 = load_atlas("GSE202379"), GSE244832 = load_atlas("GSE244832"))
cpm_mat <- function(m) t(t(m) / pmax(colSums(m), 1)) * 1e6

# Precompute donor-internal non-hepatocyte pools and protein-coding soup markers (ZERO disease-label use)
prep_cohort_ambient <- function(k) {
  L <- CO[[k]]; H <- L$Hepatocytes
  hm <- H$meta[!is.na(n_cells) & n_cells >= MIN_CELLS & !is.na(disease)]
  hc <- H$counts[, hm$donor, drop = FALSE]
  genes_h <- rownames(hc)
  donors_h <- colnames(hc)
  hcpm <- cpm_mat(hc)
  tot_h <- colSums(hc)
  others <- setdiff(names(L), "Hepatocytes")

  # Explicit pooled total-non-hepatocyte UMI soup profile per donor (sum across ALL non-hep lineages in that donor)
  pool_counts <- matrix(0, nrow = length(genes_h), ncol = length(donors_h), dimnames = list(genes_h, donors_h))
  aligned_others <- list()
  for (l in others) {
    X <- L[[l]]
    if (nrow(X$counts) == 0L || ncol(X$counts) == 0L) next
    common_d <- intersect(donors_h, colnames(X$counts))
    if (!length(common_d)) next
    g <- intersect(genes_h, rownames(X$counts))
    if (!length(g)) next
    mat <- matrix(0, nrow = length(genes_h), ncol = length(donors_h), dimnames = list(genes_h, donors_h))
    mat[g, common_d] <- X$counts[g, common_d, drop = FALSE]
    pool_counts <- pool_counts + mat
    aligned_others[[l]] <- mat
  }
  pool_cpm <- cpm_mat(pool_counts)

  # 1) Subtract individual non-hepatocyte lineages using ONLY donor's own profile (NO group-mean imputation)
  A_lin <- matrix(0, nrow = length(genes_h), ncol = length(donors_h), dimnames = list(genes_h, donors_h))
  mk_list <- list(); rh_list <- list()

  pooled_h_cpm <- rowSums(hc) / max(sum(hc), 1) * 1e6
  for (l in names(aligned_others)) {
    xc <- aligned_others[[l]]
    if (sum(xc) <= 0) next
    xcpm <- cpm_mat(xc)
    pooled_l_cpm <- rowSums(xc) / max(sum(xc), 1) * 1e6
    cand <- genes_h[(genes_h %in% protein_coding_syms) & !(genes_h %in% TARGET) & pooled_l_cpm >= MARK_CPM]
    score <- log2(pooled_l_cpm[cand] + 1) - log2(pooled_h_cpm[cand] + 1)
    cand <- cand[score > 0]
    score <- score[score > 0]
    if (!length(cand)) next
    markers <- head(cand[order(-score)], N_MARK)
    mk_list[[l]] <- data.table(cohort = k, lineage = l, gene = markers, lineage_cpm = pooled_l_cpm[markers],
                               hep_cpm = pooled_h_cpm[markers], log2_ratio = score[match(markers, cand)])

    rho_vec <- vapply(donors_h, function(d) {
      if (sum(xc[, d]) <= 0) return(0)
      valid_m <- markers[pool_counts[markers, d] >= MIN_POOL_MARKER_COUNTS & xcpm[markers, d] > 0]
      if (!length(valid_m)) return(0)  # Never clamp 0-marker donors to rho = 1.0
      r <- hcpm[valid_m, d] / xcpm[valid_m, d]
      val <- min(1, median(r))
      if (!is.finite(val) || val < 0) 0 else val
    }, numeric(1))

    n_mk_valid <- vapply(donors_h, function(d) {
      if (sum(xc[, d]) <= 0) return(0L)
      sum(pool_counts[markers, d] >= MIN_POOL_MARKER_COUNTS & xcpm[markers, d] > 0)
    }, integer(1))

    rh_list[[l]] <- data.table(cohort = k, lineage = l, donor = donors_h, rho = rho_vec,
                               n_valid_markers = n_mk_valid, own_pool_umi = colSums(xc)[donors_h])
    A_lin <- A_lin + sweep(xcpm, 2, rho_vec * tot_h / 1e6, `*`)
  }

  # 2) Explicit pooled total-non-hepatocyte UMI soup profile per donor (captures rare <20-cell lineages too)
  pooled_nonhep_cpm <- rowSums(pool_counts) / max(sum(pool_counts), 1) * 1e6
  cand_p <- genes_h[(genes_h %in% protein_coding_syms) & !(genes_h %in% TARGET) & pooled_nonhep_cpm >= MARK_CPM]
  score_p <- log2(pooled_nonhep_cpm[cand_p] + 1) - log2(pooled_h_cpm[cand_p] + 1)
  cand_p <- cand_p[score_p > 0]; score_p <- score_p[score_p > 0]
  markers_p <- head(cand_p[order(-score_p)], N_MARK)
  if (length(markers_p)) {
    mk_list[["POOLED_NONHEP"]] <- data.table(
      cohort = k, lineage = "POOLED_NONHEP", gene = markers_p,
      lineage_cpm = pooled_nonhep_cpm[markers_p], hep_cpm = pooled_h_cpm[markers_p],
      log2_ratio = score_p[match(markers_p, cand_p)]
    )
  }
  rho_pool <- vapply(donors_h, function(d) {
    if (!length(markers_p) || sum(pool_counts[, d]) <= 0) return(0)
    valid_m <- markers_p[pool_counts[markers_p, d] >= MIN_POOL_MARKER_COUNTS & pool_cpm[markers_p, d] > 0]
    if (!length(valid_m)) return(0)
    r <- hcpm[valid_m, d] / pool_cpm[valid_m, d]
    val <- min(1, median(r))
    if (!is.finite(val) || val < 0) 0 else val
  }, numeric(1))
  n_mk_pool <- vapply(donors_h, function(d) {
    if (!length(markers_p) || sum(pool_counts[, d]) <= 0) return(0L)
    sum(pool_counts[markers_p, d] >= MIN_POOL_MARKER_COUNTS & pool_cpm[markers_p, d] > 0)
  }, integer(1))
  rh_list[["POOLED_NONHEP"]] <- data.table(cohort = k, lineage = "POOLED_NONHEP", donor = donors_h,
                                           rho = rho_pool, n_valid_markers = n_mk_pool,
                                           own_pool_umi = colSums(pool_counts)[donors_h])
  A_pool <- sweep(pool_cpm, 2, rho_pool * tot_h / 1e6, `*`)

  # Combine: per-donor maximum of multi-lineage sum and explicit pooled non-hep soup
  A <- pmax(A_lin, A_pool)

  list(
    cohort = k, meta = hm, hc = hc, A = A, A_lin = A_lin, A_pool = A_pool,
    aligned_others = aligned_others, pool_counts = pool_counts,
    markers = rbindlist(mk_list, fill = TRUE), rho = rbindlist(rh_list, fill = TRUE)
  )
}

PREP <- setNames(lapply(names(CO), prep_cohort_ambient), names(CO))

# Run ambient subtraction & share calculation for a given disease-label vector (called inside observed AND every permutation!)
correct_cohort <- function(k, dis_vec) {
  P <- PREP[[k]]
  hc <- P$hc
  A <- P$A
  corr <- pmax(round(hc - A), 0)
  share <- A / pmax(hc, 1)
  ctrl_idx <- which(dis_vec == 0 & !is.na(dis_vec))
  dis_idx <- which(dis_vec == 1 & !is.na(dis_vec))
  share_ctrl <- apply(share[, ctrl_idx, drop = FALSE], 1L, median, na.rm = TRUE)
  share_disease <- apply(share[, dis_idx, drop = FALSE], 1L, median, na.rm = TRUE)
  list(meta = P$meta, corrected = corr, share_ctrl = share_ctrl, share_disease = share_disease)
}

# Cohort usability & weights
W <- rbindlist(lapply(names(PREP), function(k) {
  m <- PREP[[k]]$meta
  data.table(cohort = k, n0 = sum(m$disease == 0), n1 = sum(m$disease == 1))
}))
W[, usable := n0 >= MIN_GROUP & n1 >= MIN_GROUP][, weight := sqrt(n0 * n1 / (n0 + n1))]
usable_cohorts <- W[usable == TRUE, cohort]

# Hepatocyte predominance metrics from raw multi-cohort pseudobulk (all 16 lineages in NUCLEUS-FIX + V5 + V3)
all_g <- unique(unlist(lapply(usable_cohorts, function(k) rownames(PREP[[k]]$hc))))
sum_hep <- setNames(numeric(length(all_g)), all_g)
tot_hep <- 0
lin_sums <- list()
lin_tots <- list()
for (k in usable_cohorts) {
  P <- PREP[[k]]
  g_h <- intersect(all_g, rownames(P$hc))
  sum_hep[g_h] <- sum_hep[g_h] + rowSums(P$hc[g_h, , drop = FALSE])
  tot_hep <- tot_hep + sum(P$hc)
  for (l in names(P$aligned_others)) {
    mat <- P$aligned_others[[l]]
    if (is.null(lin_sums[[l]])) {
      lin_sums[[l]] <- setNames(numeric(length(all_g)), all_g)
      lin_tots[[l]] <- 0
    }
    g_l <- intersect(all_g, rownames(mat))
    lin_sums[[l]][g_l] <- lin_sums[[l]][g_l] + rowSums(mat[g_l, , drop = FALSE])
    lin_tots[[l]] <- lin_tots[[l]] + sum(mat)
  }
}
hep_cpm_pb <- sum_hep / max(tot_hep, 1) * 1e6
nonhep_cpm_mat <- do.call(cbind, lapply(names(lin_sums), function(l) lin_sums[[l]] / max(lin_tots[[l]], 1) * 1e6))
max_nonhep_cpm_pb <- apply(nonhep_cpm_mat, 1L, max)
sum_nonhep_cpm_pb <- rowSums(nonhep_cpm_mat)
sc_hep_share_pb <- hep_cpm_pb / pmax(hep_cpm_pb + sum_nonhep_cpm_pb, 1e-12)

HEP_PRED <- data.table(
  symbol = all_g,
  hep_cpm = hep_cpm_pb[all_g],
  max_nonhep_cpm = max_nonhep_cpm_pb[all_g],
  sc_hep_share = sc_hep_share_pb[all_g],
  hep_predominant = (sc_hep_share_pb[all_g] >= 0.50) | (hep_cpm_pb[all_g] >= max_nonhep_cpm_pb[all_g])
)

# Permutations of donor disease labels within each cohort
set.seed(SEED)
LAB <- lapply(PREP, function(P) {
  m <- P$meta
  Pmat <- cbind(m$disease, vapply(seq_len(N_PERM), function(i) sample(m$disease), numeric(nrow(m))))
  rownames(Pmat) <- m$donor
  Pmat
})

run_perm <- function(r, return_full = FALSE) {
  cr_list <- setNames(lapply(usable_cohorts, function(k) {
    dis_r <- LAB[[k]][PREP[[k]]$meta$donor, r]
    correct_cohort(k, dis_r)
  }), usable_cohorts)

  de <- rbindlist(lapply(usable_cohorts, function(k) {
    x <- cr_list[[k]]
    m <- x$meta
    dis <- LAB[[k]][m$donor, r]
    X <- cbind(1, disease = dis, female = m$female)
    if (qr(X)$rank < ncol(X)) X <- X[, 1:2]
    y <- DGEList(x$corrected[, m$donor, drop = FALSE])
    keep_expr <- filterByExpr(y, design = X)
    y <- calcNormFactors(y[keep_expr, , keep.lib.sizes = FALSE])
    fit <- glmQLFit(estimateDisp(y, X), X, robust = TRUE)
    tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
    keep_t <- rownames(tt) %in% TARGET
    data.table(cohort = k, symbol = rownames(tt)[keep_t], logFC = tt$logFC[keep_t], P = tt$PValue[keep_t])
  }))

  sh <- rbindlist(lapply(usable_cohorts, function(k) {
    data.table(
      cohort = k,
      symbol = names(cr_list[[k]]$share_disease),
      ambient_share_ctrl = unname(cr_list[[k]]$share_ctrl),
      ambient_share_disease = unname(cr_list[[k]]$share_disease)
    )
  }))

  de <- merge(de, W[usable == TRUE, .(cohort, weight)], by = "cohort")
  de[, z := sign(logFC) * qnorm(pmax(P, 1e-300) / 2, lower.tail = FALSE)]
  s <- de[, .(n_cohorts = .N, Z = sum(weight * z) / sqrt(sum(weight^2)),
              n_pos = sum(logFC > 0), n_neg = sum(logFC < 0),
              mean_logFC = sum(weight^2 * logFC) / sum(weight^2),
              cohorts = paste(sort(cohort), collapse = ";")), by = symbol][n_cohorts >= 2]
  s[, p := 2 * pnorm(-abs(Z))][, q := p.adjust(p, "BH")]
  s[, call := q < 0.05 & sign(Z) == BULK_SIGN[symbol] & fifelse(Z > 0, n_pos, n_neg) >= 2]

  aw <- merge(de[, .(cohort, symbol, weight)], sh, by = c("cohort", "symbol"))[
    , .(
      ambient_share_ctrl = sum(weight * ambient_share_ctrl) / sum(weight),
      ambient_share_disease = sum(weight * ambient_share_disease) / sum(weight)
    ), by = symbol]
  aw[, ambient_share_max := pmax(ambient_share_ctrl, ambient_share_disease)]
  s <- merge(s, aw, by = "symbol", all.x = TRUE)
  s <- merge(s, HEP_PRED, by = "symbol", all.x = TRUE)
  s[, within_hepatocyte_corrected := call &
      !is.na(ambient_share_ctrl) & !is.na(ambient_share_disease) &
      pmax(ambient_share_ctrl, ambient_share_disease) < AMB_MAX &
      !is.na(hep_predominant) & hep_predominant == TRUE]

  if (return_full) {
    return(list(de = de, sh = sh, s = s))
  }
  data.table(run = r - 1L, calls = sum(s$call), within_hepatocyte_corrected = sum(s$within_hepatocyte_corrected))
}

cat("Running observed C3b redesigned ambient correction...\n")
obs_res <- run_perm(1L, return_full = TRUE)
S <- obs_res$s
fwrite(obs_res$de, file.path(OUT, "hep_corrected_de.tsv.gz"), sep = "\t")
fwrite(obs_res$sh[symbol %in% TARGET], file.path(OUT, "ambient_share.tsv.gz"), sep = "\t")
fwrite(S, file.path(OUT, "hep_corrected_stouffer.tsv.gz"), sep = "\t")
wr(rbindlist(lapply(PREP, `[[`, "markers")), "markers.tsv")
wr(rbindlist(lapply(PREP, `[[`, "rho")), "rho.tsv")

cat(sprintf("Running %d donor-label permutations (with correct_cohort() inside every permutation) on %d cores...\n", N_PERM, NCORE))
perm_list <- mclapply(seq_len(N_PERM), function(i) run_perm(i + 1L, return_full = FALSE),
                      mc.cores = NCORE, mc.preschedule = TRUE)
stopifnot(!any(vapply(perm_list, inherits, TRUE, "try-error")))
perm <- rbindlist(perm_list)
stopifnot(nrow(perm) == N_PERM)
obs <- data.table(run = 0L, calls = sum(S$call), within_hepatocyte_corrected = sum(S$within_hepatocyte_corrected))
wr(rbind(obs, perm), "permutation_calls.tsv")

ctrl <- S[match(CONTROLS, symbol)][, symbol := CONTROLS]
wr(ctrl, "controls.tsv")
n_ctrl <- sum(ctrl$within_hepatocyte_corrected, na.rm = TRUE)
ctrl_pass <- n_ctrl <= 1L
perm_p_val <- (sum(perm$within_hepatocyte_corrected >= obs$within_hepatocyte_corrected) + 1) / (N_PERM + 1)
perm_pass <- obs$within_hepatocyte_corrected > max(perm$within_hepatocyte_corrected)

cat_lin <- fread(IN["catalog"], select = c("ensembl_id", "sc_best_lineage"))
G <- merge(
  P3[, .(ensembl_id, symbol, primary_evidence_class, bulk_log2fc, retained, adj_q, c3_class = class, c3_hep_call = hep_call)],
  S[, .(symbol, corr_Z = Z, corr_q = q, corr_logFC = mean_logFC, corr_n_cohorts = n_cohorts, corr_call = call,
        ambient_share_ctrl, ambient_share_disease, ambient_share = ambient_share_max,
        hep_cpm, max_nonhep_cpm, sc_hep_share, hep_predominant, within_hepatocyte_corrected)],
  by = "symbol", all.x = TRUE
)
G <- merge(G, cat_lin, by = "ensembl_id", all.x = TRUE)
G[is.na(within_hepatocyte_corrected), within_hepatocyte_corrected := FALSE]
G[, layer_usable := ctrl_pass && perm_pass]
fwrite(G, file.path(OUT, "gene_calls.tsv.gz"), sep = "\t", na = "NA")

G[, dir := fifelse(bulk_log2fc > 0, "up", "down")]
G[, main_lineage := fifelse(is.na(sc_best_lineage) | sc_best_lineage == "", "none",
                            fifelse(sc_best_lineage == "Hepatocytes", "Hepatocytes", "other"))]

non_hep_best_whc <- G[within_hepatocyte_corrected == TRUE & main_lineage == "other"]
marker_leak <- G[symbol %in% c("ACTA2", "EPCAM", "CLEC4M") & within_hepatocyte_corrected == TRUE]
cat(sprintf("Non-hepatocyte best-lineage genes in within_hepatocyte_corrected: %d\n", nrow(non_hep_best_whc)))
cat(sprintf("ACTA2/EPCAM/CLEC4M in within_hepatocyte_corrected: %d\n", nrow(marker_leak)))

summary_tbl <- rbindlist(list(
  data.table(
    item = c(
      "controls_called", "control_pass", "observed_corrected_calls",
      "observed_within_hepatocyte_corrected", "n_permutations",
      "max_permuted_within_hepatocyte_corrected", "perm_p_within_hepatocyte_corrected",
      "permutation_pass", "layer_usable",
      "non_hep_best_lineage_in_whc", "marker_leak_ACTA2_EPCAM_CLEC4M",
      "n_rho_clamped_to_1_in_zero_marker_donors"
    ),
    value = c(
      n_ctrl, as.numeric(ctrl_pass), obs$calls, obs$within_hepatocyte_corrected, N_PERM,
      max(perm$within_hepatocyte_corrected), perm_p_val,
      as.numeric(perm_pass), as.numeric(ctrl_pass && perm_pass),
      nrow(non_hep_best_whc), nrow(marker_leak),
      sum(rbindlist(lapply(PREP, `[[`, "rho"))[n_valid_markers == 0L]$rho == 1)
    )
  ),
  G[within_hepatocyte_corrected == TRUE, .(item = paste0("whc_", dir, "_main_", main_lineage), value = .N), by = .(dir, main_lineage)][, .(item, value)],
  G[c3_class == "within_hepatocyte", .(item = paste0("c3_within_hep_", fifelse(within_hepatocyte_corrected, "kept", "lost"), "_", dir), value = .N),
    by = .(within_hepatocyte_corrected, dir)][, .(item, value)],
  G[, .(item = paste0("median_ambient_share_max_", c3_class), value = median(ambient_share, na.rm = TRUE)), by = c3_class][, .(item, value)]
), use.names = TRUE)
wr(summary_tbl, "summary.tsv")
print(summary_tbl)

NAMED <- c("CCL20", "CXCL10", "IL32", "TRIM31", "AKR1B10", "SOX9", "KRT7", "KRT19", "EPCAM", "CLEC4M",
           "SPP1", "KRT8", "HKDC1", "PKM", "CYP2C19", "GNMT", "MAT1A", "TNFSF10", "CDKN1A", "MDM2",
           "ZMAT3", "RRM2B", "SPATA18", "PURPL", CONTROLS)
wr(G[symbol %in% NAMED, .(symbol, primary_evidence_class, sc_best_lineage, c3_class, bulk_log2fc, corr_logFC, corr_q,
                           corr_n_cohorts, ambient_share_ctrl, ambient_share_disease, ambient_share,
                           hep_cpm, max_nonhep_cpm, sc_hep_share, hep_predominant, within_hepatocyte_corrected)][
                             order(-within_hepatocyte_corrected, symbol)], "named_genes.tsv")

# Record any surviving sc_best_lineage != Hepatocytes for transparency and verify markers ACTA2/EPCAM/CLEC4M == 0
if (nrow(non_hep_best_whc) > 0L) {
  fwrite(non_hep_best_whc[, .(symbol, sc_best_lineage, c3_class, ambient_share_ctrl, ambient_share_disease,
                              ambient_share, hep_cpm, max_nonhep_cpm, sc_hep_share, corr_logFC, corr_q)],
         file.path(OUT, "non_hep_best_lineage_survivors.tsv"), sep = "\t")
}
stopifnot(nrow(marker_leak) == 0L)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("C3B_AMBIENT_20261118_DONE\n")
