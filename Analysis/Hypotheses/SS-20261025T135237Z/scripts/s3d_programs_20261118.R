# SS job S3d LOCO (2026-11-18): Program-level Fig 4D replication using separated 8-part ILR (ilr1..ilr7)
# plus 22-autosome LOCO envelope and dominant-chromosome LOCO checks.
SS_SCRIPTS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts"
source(file.path(SS_SCRIPTS, "ss_common.R"))
A <- ss_args(c("composition", "loco_dir", "norm", "h54_design", "h54"))
NORM <- if (is.null(A[["norm"]])) "rle" else A[["norm"]]; stopifnot(NORM %in% c("rle", "stored"))
H54_DESIGN <- if (is.null(A[["h54_design"]])) "ss" else A[["h54_design"]]; stopifnot(H54_DESIGN %in% c("ss", "h54"))
RUN_H54 <- is.null(A[["h54"]]) || A[["h54"]] == "run"

H54_LIB <- file.path(SS_ROOT, "Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R")
FIG4E <- file.path(SS_ROOT, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis")
SHA <- c(
  "Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R" = "33533ebae942349c14f7441e78f66e43804b16e6e88e7c37f4ec07135f2799a4",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/coordinated_edges.tsv" = "af6a58989639aa382dbdb7a39d3c155fbb933773a5e1200b87eb3c17684716a3",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/program_communities.tsv" = "f4b1b6e51b6e0c3edbdcdcff80654b34ea192957e918c8de4cd5571171bc0480",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/program_interaction_edges.tsv.gz" = "645fce85f9d9af20656de53594857999485c3de9729627e6e1aa1e294b769776",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/sealed_inputs/holdout_dge.rds" = "416ffe9812051308cfbd9b47ec651d0aa0bb3d8ba574f6faf618a6144ab06481",
  "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478/meta_effects.tsv" = "7979a3d16e9d893436a4a25d77698fe76a7928e5e103df7f9be90e9453cd8c9d",
  "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/program_weight_coverage_old_vs_new.tsv" = "a9e4e90f08068ed9be9028127449a3cdc3d5ce08c4ca27ce3f668e5082e989c9")
ss_check_sha(setNames(SHA, file.path(SS_ROOT, names(SHA))))
FZ <- ss_frozen(c("programs/program_registry_v2.tsv", "stage_continuum/hotspot_stage_vs_continuum.tsv"))
source(H54_LIB)
SEED <- 20260923L
STROMAL <- "hotspot_hepatocytes_f05c535ae5bbc0b9"; DUCT <- "hotspot_hepatocytes_48f39dd4d817a10e"
LNC <- c("LINC02331", "MIR4435-2HG")
w <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

U <- ss_units(A$input, A$stage)
dge <- readRDS(file.path(A$input, "dge_included.rds"))
stopifnot(setequal(colnames(dge), U$sample_alias))
U <- U[match(colnames(dge), sample_alias)]
registry <- fread(FZ[1]); membership <- fread(program_membership)
stopifnot(nrow(registry) == 117L, all(c(STROMAL, DUCT) %in% registry$program_uid))
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name", "chromosome"))

qc <- data.table(item = c("n_included", "n_rows_dge", "rows_equal_release_24115", "n_staged", "n_staged_with_sex",
                          "n_batches", "norm"),
                 value = c(nrow(U), nrow(dge), nrow(dge) == 24115L, sum(!is.na(U$stage_num)),
                           sum(!is.na(U$stage_num) & !is.na(U$inferred_sex)), uniqueN(U$batch), NORM))

if (NORM == "rle") {
  dge <- DGEList(counts = dge$counts)
  dge <- calcNormFactors(dge, method = "RLE")
} else stopifnot(any(dge$samples$norm.factors != 1))
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
SM <- collapse_symbols(logcpm, rownames(logcpm), annotation); rm(logcpm)
meta_sc <- data.table(sample_id = U$sample_alias, dataset = "cohort")
sc <- score_programs(SM, meta_sc, membership, registry, program_weight_coverage)
S <- sc$primary
cvg <- merge(sc$coverage, sc$cohort_coverage[, .(feature_id, cohort_weight_coverage, cohort_testable)], by = "feature_id")
cvg[, testable_ss := testable & cohort_testable & rowSums(is.finite(S[feature_id, , drop = FALSE])) == ncol(S)]
TESTABLE <- cvg[testable_ss == TRUE, feature_id]
fwrite(data.table(sample_alias = colnames(S), t(S)), file.path(OUT, "program_scores.tsv.gz"), sep = "\t")
log_msg("scored %d units; %d of 117 programs testable", ncol(S), length(TESTABLE))

D <- U[!is.na(stage_num) & !is.na(inferred_sex)]
D[, sex := factor(inferred_sex)]
rhs <- function(extra = character(0), d = D) {
  t <- c("stage_num", if (uniqueN(d$batch) > 1) "batch", if (uniqueN(d$sex) > 1) "sex", extra)
  as.formula(paste("score ~", paste(t, collapse = " + ")))
}
QUAL <- c("log10_assigned", "unique_pct", "intron_exon_ratio")
has_qual <- all(QUAL %in% names(D)) && all(is.finite(as.matrix(D[, ..QUAL])))

# Dynamically accept ilr1..ilr7 (separated 8-part ILR) or ilr1..ilr6
ILR <- NULL
if (!is.null(A[["composition"]])) {
  cp <- fread(A[["composition"]])
  ic <- grep("^ilr[0-9]+$", names(cp), value = TRUE)
  stopifnot("sample_alias" %in% names(cp), length(ic) >= 6L)
  D <- merge(D, cp[, c("sample_alias", ic), with = FALSE], by = "sample_alias", all.x = TRUE, sort = FALSE)
  ILR <- ic
  log_msg("Using %d ILR coordinates (%s)", length(ILR), paste(ILR, collapse = ","))
}

fit_one <- function(fid, form, d) {
  d <- copy(d); d[, score := as.numeric(S[fid, match(sample_alias, colnames(S))])]
  f <- hc3_fit(form, d, "stage_num")
  if (!isTRUE(f$estimable)) return(data.table(feature_id = fid, estimable = FALSE))
  data.table(feature_id = fid, estimable = TRUE, n = f$n, df = f$df, beta = f$estimate, se_hc3 = f$se,
             t = f$statistic, p_two_sided = f$p_value,
             p_one_sided_increase = pt(f$statistic, df = f$df, lower.tail = FALSE))
}

variants <- list(primary = list(form = rhs(), d = D), without_F4 = list(form = rhs(d = D[stage_num < 4]), d = D[stage_num < 4]))
if (has_qual) variants$quality <- list(form = rhs(QUAL), d = D)
if (!is.null(ILR)) variants$composition_ilr <- list(form = rhs(ILR), d = D[complete.cases(D[, ..ILR])])
fss4 <- rbindlist(lapply(names(variants), function(v) cbind(model = v, rbindlist(lapply(c(STROMAL, DUCT), function(f)
  if (f %in% TESTABLE) fit_one(f, variants[[v]]$form, variants[[v]]$d) else data.table(feature_id = f, estimable = FALSE))))), fill = TRUE)
fss4[, p_holm := NA_real_]
for (v in names(variants)) fss4[model == v & estimable == TRUE, p_holm := p.adjust(p_one_sided_increase, method = "holm")]
fss4[, pass := estimable == TRUE & beta > 0 & p_holm < 0.05]
fss4[, program := fifelse(feature_id == STROMAL, "Stromal ECM (IGFBP7)", "Ductular injury (BICC1)")]
w(fss4, "FSS4_confirmatory_programs.tsv")

fss5_sens <- rbindlist(lapply(setdiff(names(variants), "primary"), function(v) {
  r <- rbindlist(lapply(TESTABLE, fit_one, form = variants[[v]]$form, d = variants[[v]]$d), fill = TRUE)
  r[estimable == TRUE, q_bh := p.adjust(p_two_sided, method = "BH")]
  cbind(model = v, r)
}), fill = TRUE)

# Primary FSS5
fss5 <- rbindlist(lapply(TESTABLE, fit_one, form = rhs(), d = D), fill = TRUE)
fss5[estimable == TRUE, q_bh := p.adjust(p_two_sided, method = "BH")]
fss5 <- merge(registry[, .(feature_id = program_uid, cell_type, module_name)], fss5, by = "feature_id")
stg_map <- fread(FZ[2])
fss5 <- merge(fss5, stg_map[, .(feature_id, frozen_stage_beta = beta, stage_supported,
                                fig4e_4cohort_beta_fibrosis, fig4e_4cohort_q)], by = "feature_id", all.x = TRUE)
w(fss5, "FSS5_all_programs.tsv"); w(fss5_sens, "FSS5_sensitivities.tsv")

# 22-autosome LOCO program sensitivity (envelope across chr1..chr22 + modal-chromosome-held-out LOCO)
if (!is.null(A[["loco_dir"]]) && dir.exists(A[["loco_dir"]])) {
  log_msg("Running 22-autosome LOCO program sensitivity from %s", A[["loco_dir"]])
  sym_chr <- unique(annotation[!is.na(gene_name) & nzchar(gene_name) & chromosome %in% paste0("chr", 1:22),
                               .(symbol = toupper(trimws(gene_name)), chromosome)])
  memb_chr <- merge(
    membership[!is.na(mapped_symbol) & mapped_symbol != "",
               .(weight = sum(as.numeric(original_l1_weight))), by = .(feature_id = program_uid, symbol = toupper(mapped_symbol))],
    sym_chr, by = "symbol", allow.cartesian = TRUE
  )
  prog_chr_wt <- memb_chr[, .(chr_weight = sum(weight)), by = .(feature_id, chromosome)]
  setorder(prog_chr_wt, feature_id, -chr_weight)
  top_chr_map <- prog_chr_wt[, .SD[1L], by = feature_id]

  chrs <- paste0("chr", 1:22)
  loco_all <- rbindlist(lapply(chrs, function(ck) {
    f_ck <- file.path(A[["loco_dir"]], paste0(ck, "_steatosite_ilr8.tsv"))
    stopifnot(file.exists(f_ck))
    cp_k <- fread(f_ck)
    ic_k <- grep("^ilr[0-9]+$", names(cp_k), value = TRUE)
    D_k <- merge(U[!is.na(stage_num) & !is.na(inferred_sex)][, sex := factor(inferred_sex)],
                 cp_k[, c("sample_alias", ic_k), with = FALSE], by = "sample_alias", all.x = TRUE, sort = FALSE)
    D_k <- D_k[complete.cases(D_k[, ..ic_k])]
    form_k <- rhs(ic_k, D_k)
    r <- rbindlist(lapply(TESTABLE, fit_one, form = form_k, d = D_k), fill = TRUE)
    r[estimable == TRUE, q_bh := p.adjust(p_two_sided, method = "BH")]
    r[, held_out_chr := ck]
    r
  }), fill = TRUE)

  fwrite(loco_all, file.path(OUT, "FSS5_loco_22chr.tsv.gz"), sep = "\t")
  w(loco_all[feature_id %in% c(STROMAL, DUCT)], "FSS4_loco_22chr.tsv")

  modal_rows <- merge(loco_all, top_chr_map[, .(feature_id, modal_chr = chromosome, modal_chr_weight = chr_weight)], by = "feature_id")
  modal_only <- modal_rows[held_out_chr == modal_chr]
  env <- loco_all[, .(
    n_chr = .N,
    beta_median = median(beta, na.rm = TRUE),
    beta_min = min(beta, na.rm = TRUE),
    beta_max = max(beta, na.rm = TRUE),
    t_median = median(t, na.rm = TRUE),
    p_median = median(p_two_sided, na.rm = TRUE),
    q_max_across_22chr = max(q_bh, na.rm = TRUE),
    n_chr_q_lt_0.05 = sum(q_bh < 0.05, na.rm = TRUE)
  ), by = feature_id]
  loco_sum <- merge(env, modal_only[, .(feature_id, modal_chr, modal_chr_weight,
                                        beta_modal_loco = beta, se_modal_loco = se_hc3,
                                        t_modal_loco = t, p_modal_loco = p_two_sided, q_modal_loco = q_bh)],
                    by = "feature_id", all.x = TRUE)
  std_comp <- fss5_sens[model == "composition_ilr", .(feature_id, beta_std8 = beta, t_std8 = t, q_std8 = q_bh)]
  loco_sum <- merge(loco_sum, std_comp, by = "feature_id", all.x = TRUE)
  loco_sum <- merge(registry[, .(feature_id = program_uid, cell_type, module_name)], loco_sum, by = "feature_id")
  w(loco_sum, "FSS5_loco_summary.tsv")
  cat(sprintf("FSS5 LOCO summary: programs=%d, median beta cor(modal_loco, std8)=%.4f, focal_q_max: stromal=%.3e, ductular=%.3e\n",
              nrow(loco_sum),
              cor(loco_sum$beta_modal_loco, loco_sum$beta_std8, use = "pairwise.complete.obs"),
              loco_sum[feature_id == STROMAL]$q_max_across_22chr,
              loco_sum[feature_id == DUCT]$q_max_across_22chr))
}

# Permutation null for sign concordance (1,000 permutations within batch x sex, seed 20260923)
supp <- fss5[stage_supported == TRUE & estimable == TRUE]
obs_conc <- mean(sign(supp$beta) == sign(supp$frozen_stage_beta))
set.seed(SEED)
strata <- split(seq_len(nrow(D)), paste(D$batch, D$sex, sep = "::"))
conc_null <- numeric(1000L)
S_supp <- S[supp$feature_id, match(D$sample_alias, colnames(S)), drop = FALSE]
X_cov <- model.matrix(~ batch + sex, data = D)
qr_cov <- qr(X_cov)
res_y <- t(apply(S_supp, 1L, function(y) qr.resid(qr_cov, y)))
for (b in 1:1000L) {
  st_perm <- D$stage_num
  for (ix in strata) if (length(ix) > 1L) st_perm[ix] <- st_perm[sample(ix)]
  rx <- qr.resid(qr_cov, st_perm)
  slope_b <- as.numeric(res_y %*% rx) / sum(rx^2)
  conc_null[b] <- mean(sign(slope_b) == sign(supp$frozen_stage_beta))
}
fss5_sum <- data.table(
  n_testable = length(TESTABLE),
  n_stage_supported = nrow(supp),
  observed_sign_concordance = obs_conc,
  perm_p = (sum(conc_null >= obs_conc) + 1) / 1001,
  std_ilr_n_q_lt_0.05 = sum(fss5_sens[model == "composition_ilr"]$q_bh < 0.05, na.rm = TRUE)
)
w(fss5_sum, "FSS5_summary.tsv")
w(qc, "input_qc.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("S3d LOCO complete.")
