# CONF1 C3b (PRESPEC_C3b.md): within-hepatocyte calls after subtracting estimated RNA from other lineages.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HYP <- file.path(ROOT, "Analysis/Hypotheses")
RUN <- file.path(HYP, "CONF1-20261109T221846Z")
V5 <- file.path(HYP, "V5-20261021T122723Z/outputs")
V3 <- file.path(HYP, "V3-20261008T015239Z/pseudobulk")
NF <- file.path(HYP, "NUCLEUS-FIX-20261107/unique_counts")
ATLAS_LIN <- c(Hepatocytes = "Hepatocytes", Cholangiocytes = "Cholangiocytes", Fibroblasts = "Fibroblasts",
               Macrophages = "Macrophages", T_cells = "T_cells", Endothelial = "Endothelial_cells")
V5_LIN <- c(Hepatocytes = "hep", Cholangiocytes = "cholangiocytes", Fibroblasts = "fibroblasts",
            Macrophages = "macrophages", T_cells = "tcells")
V3_LIN <- c(Hepatocytes = "hepatocyte", Cholangiocytes = "cholangiocyte", Fibroblasts = "stellate_fibroblast")
IN <- c(
  c3_partition = file.path(RUN, "out_c3_partition/gene_partition.tsv.gz"),
  catalog = file.path(ROOT, "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/catalog_annotated_gene_table.tsv.gz"),
  v5_donors = file.path(V5, "donor_qc.tsv"), v5_cells = file.path(V5, "job3/cells_per_donor_set.tsv"),
  v3_donors = file.path(V3, "donor_metadata.tsv"), v3_cells = file.path(V3, "nuclei_per_donor_lineage.tsv"),
  atlas_donors = file.path(NF, "metadata/donor_metadata.tsv"), atlas_libraries = file.path(NF, "metadata/library_metadata.tsv"),
  atlas_sex = file.path(HYP, "RV1-20261031T155620Z/out-20261102T130905Z/donors.tsv"),
  annotation = file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  setNames(file.path(V5, sprintf("job3/pseudobulk_%s.tsv.gz", V5_LIN)), paste0("v5_", names(V5_LIN))),
  setNames(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", V3_LIN)), paste0("v3_", names(V3_LIN))),
  setNames(file.path(NF, sprintf("pseudobulk/%s_pseudobulk.csv", ATLAS_LIN)), paste0("atlas_", names(ATLAS_LIN))))
stopifnot(all(file.exists(IN)))
SEED <- 20261109; N_PERM <- 20L; MIN_CELLS <- 30; MIN_GROUP <- 3; MIN_LIN_OWN <- 10
N_MARK <- 30; MARK_CPM <- 50; AMB_MAX <- 0.5
CONTROLS <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "ACTA2", "LYZ", "STAB2")
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")
wr(data.table(input = names(IN), path = IN, sha256 = vapply(IN, function(f)
  sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE)), "")), "inputs_sha256.tsv")

P3 <- fread(IN["c3_partition"])[replicated == TRUE]
stopifnot(nrow(P3) == 7214)
TARGET <- P3[!is.na(symbol) & symbol != "" & !duplicated(symbol) & !(symbol %in% symbol[duplicated(symbol)]), symbol]
BULK_SIGN <- setNames(sign(P3$bulk_log2fc), P3$symbol)

# ---- loaders (as in c3_partition.R) ---------------------------------------------------------
annot <- fread(IN["annotation"], select = c("ensembl_base", "gene_name"))
by_symbol <- function(counts, symbol) {
  ok <- !is.na(symbol) & symbol != "" & !(symbol %in% symbol[duplicated(symbol)])
  counts <- counts[ok, , drop = FALSE]; rownames(counts) <- symbol[ok]; counts
}
disease_of <- function(stage_ord) fifelse(stage_ord == 0, 0, fifelse(stage_ord >= 2, 1, NA_real_))
load_v5 <- function() {
  dq <- fread(IN["v5_donors"])[included == TRUE]; cells <- fread(IN["v5_cells"])
  meta <- data.table(donor = dq$gsm, female = as.integer(dq$sex_f), disease = as.numeric(dq$disease))
  lapply(V5_LIN, function(s) { pb <- fread(file.path(V5, sprintf("job3/pseudobulk_%s.tsv.gz", s)))
    list(meta = copy(meta)[, n_cells := cells[set == s][match(meta$donor, gsm), n_cells]],
         counts = by_symbol(as.matrix(pb[, meta$donor, with = FALSE]),
                            annot$gene_name[match(sub("\\..*$", "", pb$gene_id), annot$ensembl_base)])) })
}
load_v3 <- function() {
  m <- fread(IN["v3_donors"]); n3 <- fread(IN["v3_cells"])
  meta <- data.table(donor = m$donor, female = as.integer(m$sex == "female"), disease = disease_of(m$stage_ordinal))
  lapply(V3_LIN, function(s) { pb <- fread(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", s)))
    list(meta = copy(meta)[, n_cells := n3[[s]][match(meta$donor, n3$donor)]],
         counts = by_symbol(as.matrix(pb[, meta$donor, with = FALSE]), pb$symbol)) })
}
load_atlas <- function(ds) {
  dm <- fread(IN["atlas_donors"])[dataset == ds & toupper(as.character(exclude)) == "FALSE"]
  lib <- fread(IN["atlas_libraries"])[donor %in% dm$donor]
  sex <- fread(IN["atlas_sex"], select = c("donor", "sex"))
  stopifnot(all(dm$donor %in% sex$donor))
  meta <- data.table(donor = dm$donor, female = as.integer(sex$sex[match(dm$donor, sex$donor)] == "Female"),
                     disease = disease_of(dm$stage_ordinal))
  lapply(ATLAS_LIN, function(s) { pb <- fread(file.path(NF, sprintf("pseudobulk/%s_pseudobulk.csv", s)))
    ncol_name <- paste0("n_", s); stopifnot(ncol_name %in% names(lib))
    libs <- lapply(setNames(dm$donor, dm$donor), function(d) intersect(lib[donor == d, sample], names(pb)))
    has <- lengths(libs) > 0
    counts <- vapply(libs[has], function(l) rowSums(as.matrix(pb[, l, with = FALSE])), numeric(nrow(pb)))
    mm <- copy(meta)[has][, n_cells := vapply(donor, function(d) sum(lib[donor == d][[ncol_name]]), numeric(1))]
    list(meta = mm, counts = by_symbol(counts, pb$gene)) })
}
CO <- list(GSE298719 = load_v5(), GSE256398 = load_v3(), GSE202379 = load_atlas("GSE202379"), GSE244832 = load_atlas("GSE244832"))

# ---- ambient estimate and corrected hepatocyte counts -------------------------------------
cpm_mat <- function(m) t(t(m) / pmax(colSums(m), 1)) * 1e6
correct_cohort <- function(k) {
  L <- CO[[k]]; H <- L$Hepatocytes
  hm <- H$meta[!is.na(n_cells) & n_cells >= MIN_CELLS]
  hc <- H$counts[, hm$donor, drop = FALSE]
  hcpm <- cpm_mat(hc)
  others <- setdiff(names(L), "Hepatocytes")
  A <- matrix(0, nrow(hc), ncol(hc), dimnames = dimnames(hc))
  mk <- list(); rh <- list()
  for (l in others) {
    X <- L[[l]]; g <- intersect(rownames(hc), rownames(X$counts))
    xc <- X$counts[g, , drop = FALSE]; xcpm <- cpm_mat(xc)
    own <- X$meta$donor[!is.na(X$meta$n_cells) & X$meta$n_cells >= MIN_LIN_OWN & X$meta$donor %in% colnames(xc)]
    if (length(own) < 2) next
    pooled_l <- rowSums(xc[, own, drop = FALSE]) / sum(xc[, own]) * 1e6
    pooled_h <- rowSums(hc[g, , drop = FALSE]) / sum(hc[g, ]) * 1e6
    cand <- g[!(g %in% TARGET) & pooled_l >= MARK_CPM]
    score <- log2(pooled_l[cand] + 1) - log2(pooled_h[cand] + 1)
    markers <- head(cand[order(-score)], N_MARK)
    mk[[l]] <- data.table(cohort = k, lineage = l, gene = markers, lineage_cpm = pooled_l[markers],
                          hep_cpm = pooled_h[markers], log2_ratio = score[match(markers, cand)])
    # per-donor lineage profile, or the mean profile of donors in the same diagnosis group
    grp <- X$meta[match(own, X$meta$donor), disease]
    prof <- matrix(NA_real_, length(g), ncol(hc), dimnames = list(g, colnames(hc)))
    for (d in colnames(hc)) {
      if (d %in% own) { prof[, d] <- xcpm[, d]; next }
      dg <- hm[donor == d, disease]
      pool <- if (!is.na(dg) && sum(grp == dg, na.rm = TRUE) >= 1) own[!is.na(grp) & grp == dg] else own
      prof[, d] <- rowMeans(xcpm[, pool, drop = FALSE])
    }
    rho <- vapply(colnames(hc), function(d) {
      r <- hcpm[markers, d] / pmax(prof[markers, d], 1e-9); min(1, median(r))
    }, numeric(1))
    rh[[l]] <- data.table(cohort = k, lineage = l, donor = colnames(hc), rho = rho, own_profile = colnames(hc) %in% own)
    A[g, ] <- A[g, ] + sweep(prof, 2, rho * colSums(hc) / 1e6, `*`)
  }
  corr <- pmax(round(hc - A), 0)
  share <- A / pmax(hc, 1)
  dis <- hm$disease == 1 & !is.na(hm$disease)
  list(meta = hm, corrected = corr, share_disease = apply(share[, dis, drop = FALSE], 1, median),
       markers = rbindlist(mk), rho = rbindlist(rh))
}
CR <- setNames(lapply(names(CO), correct_cohort), names(CO))
wr(rbindlist(lapply(CR, `[[`, "markers")), "markers.tsv")
wr(rbindlist(lapply(CR, `[[`, "rho")), "rho.tsv")
SH <- rbindlist(lapply(names(CR), function(k) data.table(cohort = k, symbol = names(CR[[k]]$share_disease),
                                                          ambient_share_disease = CR[[k]]$share_disease)))
fwrite(SH[symbol %in% TARGET], file.path(OUT, "ambient_share.tsv.gz"), sep = "\t")

# ---- hepatocyte test on corrected counts ------------------------------------------------------
W <- rbindlist(lapply(names(CR), function(k) { m <- CR[[k]]$meta[!is.na(disease)]
  data.table(cohort = k, n0 = sum(m$disease == 0), n1 = sum(m$disease == 1)) }))
W[, usable := n0 >= MIN_GROUP & n1 >= MIN_GROUP][, weight := sqrt(n0 * n1 / (n0 + n1))]
set.seed(SEED)
LAB <- lapply(CR, function(x) { m <- x$meta[!is.na(disease)]
  P <- cbind(m$disease, vapply(seq_len(N_PERM), function(i) sample(m$disease), numeric(nrow(m)))); rownames(P) <- m$donor; P })
fit_one <- function(k, r) {
  x <- CR[[k]]; m <- x$meta[!is.na(disease)]
  dis <- LAB[[k]][m$donor, r]
  X <- cbind(1, disease = dis, female = m$female)
  if (qr(X)$rank < ncol(X)) X <- X[, 1:2]
  y <- DGEList(x$corrected[, m$donor, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = X), , keep.lib.sizes = FALSE])
  fit <- glmQLFit(estimateDisp(y, X), X, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  keep <- rownames(tt) %in% TARGET
  data.table(cohort = k, symbol = rownames(tt)[keep], logFC = tt$logFC[keep], P = tt$PValue[keep])
}
combine <- function(de) {
  de <- merge(de, W[usable == TRUE, .(cohort, weight)], by = "cohort")
  de[, z := sign(logFC) * qnorm(pmax(P, 1e-300) / 2, lower.tail = FALSE)]
  s <- de[, .(n_cohorts = .N, Z = sum(weight * z) / sqrt(sum(weight^2)), n_pos = sum(logFC > 0), n_neg = sum(logFC < 0),
              mean_logFC = sum(weight^2 * logFC) / sum(weight^2), cohorts = paste(sort(cohort), collapse = ";")), by = symbol][n_cohorts >= 2]
  s[, p := 2 * pnorm(-abs(Z))][, q := p.adjust(p, "BH")]
  s[, call := q < 0.05 & sign(Z) == BULK_SIGN[symbol] & fifelse(Z > 0, n_pos, n_neg) >= 2]
  amb <- SH[, .(cohort, symbol, ambient_share_disease)]
  aw <- merge(de[, .(cohort, symbol, weight)], amb, by = c("cohort", "symbol"))[, .(ambient_share = sum(weight * ambient_share_disease) / sum(weight)), by = symbol]
  s <- merge(s, aw, by = "symbol", all.x = TRUE)
  s[, within_hepatocyte_corrected := call & !is.na(ambient_share) & ambient_share < AMB_MAX]
  s
}
run <- function(r) combine(rbindlist(lapply(W[usable == TRUE, cohort], fit_one, r = r)))
de_obs <- rbindlist(lapply(W[usable == TRUE, cohort], fit_one, r = 1L))
fwrite(de_obs, file.path(OUT, "hep_corrected_de.tsv.gz"), sep = "\t")
S <- combine(de_obs)
fwrite(S, file.path(OUT, "hep_corrected_stouffer.tsv.gz"), sep = "\t")
perm <- rbindlist(mclapply(seq_len(N_PERM), function(i) { s <- run(i + 1L)
  data.table(run = i, calls = sum(s$call), within_hepatocyte_corrected = sum(s$within_hepatocyte_corrected)) },
  mc.cores = NCORE, mc.preschedule = FALSE))
stopifnot(nrow(perm) == N_PERM)
obs <- data.table(run = 0L, calls = sum(S$call), within_hepatocyte_corrected = sum(S$within_hepatocyte_corrected))
wr(rbind(obs, perm), "permutation_calls.tsv")

ctrl <- S[match(CONTROLS, symbol)][, symbol := CONTROLS]
wr(ctrl, "controls.tsv")
n_ctrl <- sum(ctrl$within_hepatocyte_corrected, na.rm = TRUE)
ctrl_pass <- n_ctrl <= 1
perm_pass <- obs$within_hepatocyte_corrected > max(perm$within_hepatocyte_corrected)

# ---- gene table and summaries ------------------------------------------------------------------
cat_lin <- fread(IN["catalog"], select = c("ensembl_id", "sc_best_lineage"))
G <- merge(P3[, .(ensembl_id, symbol, primary_evidence_class, bulk_log2fc, retained, adj_q, c3_class = class, c3_hep_call = hep_call)],
           S[, .(symbol, corr_Z = Z, corr_q = q, corr_logFC = mean_logFC, corr_n_cohorts = n_cohorts, corr_call = call,
                 ambient_share, within_hepatocyte_corrected)], by = "symbol", all.x = TRUE)
G <- merge(G, cat_lin, by = "ensembl_id", all.x = TRUE)
G[is.na(within_hepatocyte_corrected), within_hepatocyte_corrected := FALSE]
G[, layer_usable := ctrl_pass && perm_pass]
fwrite(G, file.path(OUT, "gene_calls.tsv.gz"), sep = "\t", na = "NA")
G[, dir := fifelse(bulk_log2fc > 0, "up", "down")]
G[, main_lineage := fifelse(is.na(sc_best_lineage) | sc_best_lineage == "", "none", fifelse(sc_best_lineage == "Hepatocytes", "Hepatocytes", "other"))]
wr(rbindlist(list(
  data.table(item = c("controls_called", "control_pass", "observed_corrected_calls", "observed_within_hepatocyte_corrected",
                      "max_permuted_within_hepatocyte_corrected", "permutation_pass", "layer_usable"),
             value = c(n_ctrl, ctrl_pass, obs$calls, obs$within_hepatocyte_corrected, max(perm$within_hepatocyte_corrected), perm_pass, ctrl_pass && perm_pass)),
  G[within_hepatocyte_corrected == TRUE, .(item = paste0("whc_", dir, "_main_", main_lineage), value = .N), by = .(dir, main_lineage)][, .(item, value)],
  G[c3_class == "within_hepatocyte", .(item = paste0("c3_within_hep_", fifelse(within_hepatocyte_corrected, "kept", "lost"), "_", dir), value = .N),
    by = .(within_hepatocyte_corrected, dir)][, .(item, value)],
  G[, .(item = paste0("median_ambient_share_", c3_class), value = median(ambient_share, na.rm = TRUE)), by = c3_class][, .(item, value)]),
  use.names = TRUE), "summary.tsv")
NAMED <- c("CCL20", "CXCL10", "IL32", "TRIM31", "AKR1B10", "SOX9", "KRT7", "KRT19", "EPCAM", "SPP1", "KRT8", "HKDC1", "PKM", "CYP2C19",
           "GNMT", "MAT1A", "TNFSF10", "CDKN1A", "MDM2", "ZMAT3", "RRM2B", "SPATA18", "PURPL", CONTROLS)
wr(G[symbol %in% NAMED, .(symbol, primary_evidence_class, sc_best_lineage, c3_class, bulk_log2fc, corr_logFC, corr_q, corr_n_cohorts,
                           ambient_share, within_hepatocyte_corrected)][order(-within_hepatocyte_corrected, symbol)], "named_genes.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
