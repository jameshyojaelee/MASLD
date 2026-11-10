# SS job S3d (PRESPEC_SS.md F-SS4, F-SS5; addendum PRESPEC_SS_addendum_S3d.md). Frozen Hotspot programs in one
# held-out bulk cohort (SteatoSITE, FFPE total RNA). Nothing is refit or reselected.
# 1. Scores: release path of Fig 4E (figures/candidates/fig4e-program-axis-774-20261012T010830Z/producers/
#    program_axis_map_774.R): log2 CPM (prior 1) on the 24,115 release rows with RLE factors computed on the included
#    units (the release factors are RLE over its 774 units), collapse_symbols(), score_programs() with frozen
#    original_l1_weight; within-cohort z over all included units. Testable = weight coverage >= 0.80 on the release
#    rows and >= 0.80 of weight on genes that vary in this cohort.
# 2. Observability row per program: release-row coverage, varying-gene coverage, detected-gene coverage (CPM >= 1
#    in >= 50% of included units, full 79,095-gene counts), and the lncRNA members LINC02331, MIR4435-2HG.
# 3. F-SS4: Stromal ECM and Ductular injury, score ~ stage_num + batch + inferred_sex, OLS + HC3 (hc3_fit of
#    analysis_lib.R), one-sided increase, Holm over 2. Beside: + ILR (S3e), + quality covariates, without F4.
# 4. F-SS5: same model for every testable program, two-sided, BH over testable programs; sign agreement with the 64
#    stage-supported programs against 1,000 permutations of stage_num within batch x sex (seed 20260923).
# 5. H54 third replication cohort (PRESPEC_H54.md statistic unchanged; residual design ~ batch + sex +
#    fibrosis factor, NAS is not recorded): S1_k1, S1_k3, label-null E_k, p* = max(p_perm, p_boot), NC1-NC3, gates.
# Unit: one library per genetic individual (S2 include). Cross-sectional stage.
# Usage: Rscript s3d_programs.R --input <S2 dir> --stage <stage tsv> --out <new dir>
#        [--composition <tsv: sample_alias, ilr1..ilr6>] [--norm rle|stored] [--h54_design ss|h54] [--h54 run|skip]
#        (--norm stored and --h54_design h54 exist only for reproduction on mock inputs.)
SS_SCRIPTS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts"
source(file.path(SS_SCRIPTS, "ss_common.R"))
A <- ss_args(c("composition", "norm", "h54_design", "h54"))
NORM <- if (is.null(A[["norm"]])) "rle" else A[["norm"]]; stopifnot(NORM %in% c("rle", "stored"))
H54_DESIGN <- if (is.null(A[["h54_design"]])) "ss" else A[["h54_design"]]; stopifnot(H54_DESIGN %in% c("ss", "h54"))
RUN_H54 <- is.null(A[["h54"]]) || A[["h54"]] == "run"
OUT <- A$out
suppressPackageStartupMessages({library(edgeR); library(Matrix)})
# h54_lib.R sources config.R, analysis_lib.R and t4_interaction_lib.R (all read-only) and defines the H54 statistic.
H54_LIB <- file.path(SS_ROOT, "Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R")
PMAP <- file.path(SS_ROOT, "scripts/manuscript/program_observability_map")
FIG4E <- file.path(SS_ROOT, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis")
SHA <- c(
  "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv" = "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz" = "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189",
  "scripts/manuscript/program_observability_map/analysis_lib.R" = "a6439390506398be9283eb335b3866c45e9a510fa87aec07f1bf7fc21accb71b",
  "scripts/manuscript/program_observability_map/config.R" = "d34dd1a9b1d7146f1734267d4ac0b0d1e69d78a85ead4e784a62cd10de498889",
  "scripts/manuscript/program_observability_map/t4_interaction_lib.R" = "2d64504184a47a13d0c3c9cd2ed05b763e5d5b903a744fb3461dcc6b1d00de24",
  "Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R" = "33533ebae942349c14f7441e78f66e43804b16e6e88e7c37f4ec07135f2799a4",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/coordinated_edges.tsv" = "af6a58989639aa382dbdb7a39d3c155fbb933773a5e1200b87eb3c17684716a3",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/program_communities.tsv" = "f4b1b6e51b6e0c3edbdcdcff80654b34ea192957e918c8de4cd5571171bc0480",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/program_interaction_edges.tsv.gz" = "645fce85f9d9af20656de53594857999485c3de9729627e6e1aa1e294b769776",
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/sealed_inputs/holdout_dge.rds" = "416ffe9812051308cfbd9b47ec651d0aa0bb3d8ba574f6faf618a6144ab06481",
  "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478/meta_effects.tsv" = "7979a3d16e9d893436a4a25d77698fe76a7928e5e103df7f9be90e9453cd8c9d",
  "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/program_weight_coverage_old_vs_new.tsv" = "a9e4e90f08068ed9be9028127449a3cdc3d5ce08c4ca27ce3f668e5082e989c9")
ss_check_sha(setNames(SHA, file.path(SS_ROOT, names(SHA))))
FZ <- ss_frozen(c("programs/program_registry_v2.tsv", "stage_continuum/hotspot_stage_vs_continuum.tsv"))
source(H54_LIB)  # defines fail(), score_programs(), hc3_fit(), t4_*(), H54 helpers; sets MASLD_PROJECT_ROOT
SEED <- 20260923L
STROMAL <- "hotspot_hepatocytes_f05c535ae5bbc0b9"; DUCT <- "hotspot_hepatocytes_48f39dd4d817a10e"
LNC <- c("LINC02331", "MIR4435-2HG")
w <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

# ---- inputs ----------------------------------------------------------------------------------------------------
U <- ss_units(A$input, A$stage)
dge <- readRDS(file.path(A$input, "dge_included.rds"))
stopifnot(setequal(colnames(dge), U$sample_alias))
U <- U[match(colnames(dge), sample_alias)]
registry <- fread(FZ[1]); membership <- fread(program_membership)
stopifnot(nrow(registry) == 117L, all(c(STROMAL, DUCT) %in% registry$program_uid))
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
qc <- data.table(item = c("n_included", "n_rows_dge", "rows_equal_release_24115", "n_staged", "n_staged_with_sex",
                          "n_batches", "norm"),
                 value = c(nrow(U), nrow(dge), nrow(dge) == 24115L, sum(!is.na(U$stage_num)),
                           sum(!is.na(U$stage_num) & !is.na(U$inferred_sex)), uniqueN(U$batch), NORM))

# ---- 1. scores ---------------------------------------------------------------------------------------------------
if (NORM == "rle") {
  dge <- DGEList(counts = dge$counts)  # lib.size = column sums of the release rows, as in the release DGE
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

# ---- 2. observability row ---------------------------------------------------------------------------------------
ca <- readRDS(file.path(A$input, "counts_all_genes.rds"))
ca <- ca[, U$sample_alias, drop = FALSE]
cpm_all <- cpm(ca, lib.size = colSums(ca))
key <- sub("\\..*$", "", rownames(ca))
ann <- unique(annotation[, .(gene_key = sub("\\..*$", "", gene_id), symbol = toupper(trimws(gene_name)))])
sym_all <- ann$symbol[match(key, ann$gene_key)]
det_frac <- rowMeans(cpm_all >= 1)
det_sym <- tapply(det_frac, sym_all, max)  # a symbol is detected if any of its rows is
memb <- membership[!is.na(mapped_symbol) & mapped_symbol != "",
                   .(weight = sum(as.numeric(original_l1_weight))), by = .(feature_id = program_uid, symbol = toupper(mapped_symbol))]
memb[, detected := symbol %in% names(det_sym)[det_sym >= 0.5]]
obs <- merge(cvg, memb[, .(detected_weight_coverage = sum(weight[detected]) / sum(weight)), by = feature_id], by = "feature_id", all.x = TRUE)
rel <- fread(file.path(FIG4E, "program_weight_coverage_old_vs_new.tsv"))
obs <- merge(obs, rel[, .(feature_id, release774_weight_coverage = new_weight_coverage, release774_testable = new_testable)],
             by = "feature_id", all.x = TRUE)
obs <- merge(registry[, .(feature_id = program_uid, cell_type, module_name)], obs, by = "feature_id")
w(obs[, .(feature_id, cell_type, module_name, total_weight, release_row_weight_coverage = weight_coverage,
          varying_weight_coverage = cohort_weight_coverage, detected_weight_coverage, testable = testable_ss,
          release774_weight_coverage, release774_testable)], "observability_by_program.tsv")
lnc <- rbindlist(lapply(LNC, function(g) {
  r <- which(sym_all == g)
  data.table(symbol = g, in_counts = length(r) > 0, in_release_rows = g %in% rownames(SM),
             programs = paste(memb[symbol == g, feature_id], collapse = ";"),
             weight_in_programs = paste(signif(memb[symbol == g, weight], 4), collapse = ";"),
             frac_units_cpm_ge1 = if (length(r)) max(det_frac[r]) else NA_real_,
             median_cpm = if (length(r)) median(colSums(cpm_all[r, , drop = FALSE])) else NA_real_,
             frac_units_ge10_fragments = if (length(r)) mean(colSums(ca[r, , drop = FALSE]) >= 10) else NA_real_)
}))
w(lnc, "observability_lncRNA.tsv"); rm(ca, cpm_all); invisible(gc())

# ---- 3-4. stage models ------------------------------------------------------------------------------------------
D <- U[!is.na(stage_num) & !is.na(inferred_sex)]
D[, sex := factor(inferred_sex)]
nb <- uniqueN(D$batch)
rhs <- function(extra = character(0), d = D) {
  t <- c("stage_num", if (uniqueN(d$batch) > 1) "batch", if (uniqueN(d$sex) > 1) "sex", extra)
  as.formula(paste("score ~", paste(t, collapse = " + ")))
}
QUAL <- c("log10_assigned", "unique_pct", "intron_exon_ratio")
has_qual <- all(QUAL %in% names(D)) && all(is.finite(as.matrix(D[, ..QUAL])))
ILR <- NULL
if (!is.null(A[["composition"]])) {
  cp <- fread(A[["composition"]]); ic <- paste0("ilr", 1:6); stopifnot(all(c("sample_alias", ic) %in% names(cp)))
  D <- merge(D, cp[, c("sample_alias", ic), with = FALSE], by = "sample_alias", all.x = TRUE, sort = FALSE); ILR <- ic
}
fit_one <- function(fid, form, d) {
  d <- copy(d); d[, score := as.numeric(S[fid, match(sample_alias, colnames(S))])]
  f <- hc3_fit(d$score, d, form)
  if (!isTRUE(f$estimable)) return(data.table(feature_id = fid, estimable = FALSE, n = nrow(d)))
  t <- f$statistic[["stage_num"]]
  data.table(feature_id = fid, estimable = TRUE, n = f$n, df = f$df, beta = f$coefficients[["stage_num"]],
             se_hc3 = f$se[["stage_num"]], t = t, p_two_sided = f$p_value[["stage_num"]],
             p_one_sided_increase = stats::pt(t, f$df, lower.tail = FALSE))
}
variants <- list(primary = list(form = rhs(), d = D), without_F4 = list(form = rhs(d = D[stage_num < 4]), d = D[stage_num < 4]))
if (has_qual) variants$quality <- list(form = rhs(QUAL), d = D)
if (!is.null(ILR)) variants$composition_ilr <- list(form = rhs(ILR), d = D[complete.cases(D[, ..ILR])])
fss4 <- rbindlist(lapply(names(variants), function(v) cbind(model = v, rbindlist(lapply(c(STROMAL, DUCT), function(f)
  if (f %in% TESTABLE) fit_one(f, variants[[v]]$form, variants[[v]]$d) else data.table(feature_id = f, estimable = FALSE))))), fill = TRUE)
fss4[model == "primary", p_holm := p.adjust(p_one_sided_increase, "holm")]
fss4[model == "primary", pass := is.finite(p_holm) & p_holm < 0.05]
fss4[, program := fifelse(feature_id == STROMAL, "Stromal ECM (IGFBP7)", "Ductular injury (BICC1)")]
w(fss4, "FSS4_confirmatory_programs.tsv"); print(fss4)

fss5 <- rbindlist(lapply(TESTABLE, fit_one, form = rhs(), d = D))
fss5[, q_bh := p.adjust(p_two_sided, "BH")]
fss5_sens <- rbindlist(lapply(setdiff(names(variants), "primary"), function(v)
  cbind(model = v, rbindlist(lapply(TESTABLE, fit_one, form = variants[[v]]$form, d = variants[[v]]$d)))), fill = TRUE)
# Sign agreement with the frozen stage-supported programs, permutation null of stage_num within batch x sex.
hs <- fread(FZ[2]); stopifnot(nrow(hs) == 117L, sum(hs$stage_supported %in% TRUE) == 64L)
ref <- hs[stage_supported %in% TRUE & program_uid %in% TESTABLE, .(feature_id = program_uid, frozen_beta = beta)]
fss5 <- merge(fss5, hs[, .(feature_id = program_uid, frozen_stage_beta = beta, stage_supported)], by = "feature_id", all.x = TRUE)
fss5 <- merge(registry[, .(feature_id = program_uid, cell_type, module_name)], fss5, by = "feature_id")
fe <- fread(file.path(FIG4E, "new774_plus_GSE174478/meta_effects.tsv"))[model_kind == "primary" & axis == "fibrosis"]
fss5 <- merge(fss5, fe[, .(feature_id, fig4e_4cohort_beta_fibrosis = beta_meta, fig4e_4cohort_q = q_value)], by = "feature_id", all.x = TRUE)
w(fss5, "FSS5_all_programs.tsv"); w(fss5_sens, "FSS5_sensitivities.tsv")
Y <- t(S[ref$feature_id, match(D$sample_alias, colnames(S)), drop = FALSE])
Xc <- model.matrix(delete.response(terms(rhs())), D)
beta_stage <- function(st) { X <- Xc; X[, "stage_num"] <- st; qr.coef(qr(X), Y)["stage_num", ] }
obs_agree <- mean(sign(beta_stage(D$stage_num)) == sign(ref$frozen_beta))
strata <- interaction(D$batch, D$sex, drop = TRUE)
set.seed(SEED)
null_agree <- vapply(seq_len(1000), function(b) {
  st <- D$stage_num; for (g in split(seq_along(st), strata)) st[g] <- st[g][sample.int(length(g))]
  mean(sign(beta_stage(st)) == sign(ref$frozen_beta))
}, numeric(1))
fss5_sum <- data.table(n_units = nrow(D), n_batches = nb, n_testable = length(TESTABLE),
                       n_bh_q05 = sum(fss5$q_bh < 0.05, na.rm = TRUE),
                       n_bh_q05_up = sum(fss5$q_bh < 0.05 & fss5$beta > 0, na.rm = TRUE),
                       n_bh_q05_down = sum(fss5$q_bh < 0.05 & fss5$beta < 0, na.rm = TRUE),
                       n_stage_supported_testable = nrow(ref), sign_agreement_stage_supported = obs_agree,
                       null_mean = mean(null_agree), null_q95 = quantile(null_agree, 0.95, names = FALSE),
                       p_perm_upper = (1 + sum(null_agree >= obs_agree)) / 1001,
                       pearson_vs_frozen_stage_beta = with(fss5[is.finite(beta) & is.finite(frozen_stage_beta)], cor(beta, frozen_stage_beta)),
                       pearson_vs_fig4e_4cohort = with(fss5[is.finite(beta) & is.finite(fig4e_4cohort_beta_fibrosis)],
                                                       cor(beta, fig4e_4cohort_beta_fibrosis)))
w(fss5_sum, "FSS5_summary.tsv"); print(fss5_sum)

# ---- 5. H54 third replication cohort --------------------------------------------------------------------------
if (RUN_H54) {
  frozen <- load_frozen(); pr <- frozen$pairs; pidx <- pair_index_matrix(frozen); features <- frozen$features
  if (H54_DESIGN == "ss") {
    assign("T4_RESIDUAL_FORMULA", if (nb > 1) ~ batch + sex_final + fibrosis_factor else ~ sex_final + fibrosis_factor, envir = globalenv())
    hm <- D[, .(sample_id = sample_alias, dataset = "cohort", fibrosis_stage = as.integer(stage_num), nas_score = NA_real_,
                sex_final = inferred_sex, batch)]
  } else {  # H54 design on a mock with nas_score in units.tsv (reproduction only)
    stopifnot("nas_score" %in% names(D))
    hm <- D[, .(sample_id = sample_alias, dataset = "cohort", fibrosis_stage = as.integer(stage_num), nas_score = as.numeric(nas_score),
                sex_final = inferred_sex)]
  }
  # Expression path of PRESPEC_H54 steps 1-3: frozen 24,196-row universe (unversioned), TMM within the set, log CPM.
  ca <- readRDS(file.path(A$input, "counts_all_genes.rds"))
  ca <- ca[!grepl("_PAR_Y$", rownames(ca)), , drop = FALSE]
  rownames(ca) <- sub("\\..*$", "", rownames(ca)); stopifnot(!anyDuplicated(rownames(ca)))
  build_PM <- function(meta) {
    Mu <- ca[rownames(ca) %in% frozen$universe, meta$sample_id, drop = FALSE]
    lib <- cpm(calcNormFactors(DGEList(Mu), method = "TMM"), log = TRUE, prior.count = 1)
    sm <- collapse_symbols(lib, rownames(lib), frozen$annotation)
    panel <- intersect(sort(unique(frozen$weights$gene_symbol)), rownames(sm))
    list(meta = meta, PM = sm[panel, meta$sample_id, drop = FALSE], panel = panel, n_universe_present = sum(rownames(ca) %in% frozen$universe))
  }
  SETS <- list(all = hm, F0F1 = hm[fibrosis_stage <= 1L], without_F4 = hm[fibrosis_stage < 4L])
  if (H54_DESIGN == "h54") SETS <- SETS[c("all", "F0F1")]
  ROLE <- c(all = "primary", F0F1 = "secondary", without_F4 = "descriptive")
  X <- lapply(SETS, build_PM); rm(ca); invisible(gc())
  # Residual design per set: the formula above, without the fibrosis factor when a set has one stage level.
  form_for <- function(meta) {
    f <- T4_RESIDUAL_FORMULA
    if (uniqueN(meta$fibrosis_stage) == 1L) f <- update(f, ~ . - fibrosis_factor)
    f
  }
  block_for <- function(s, idx) {
    old <- T4_RESIDUAL_FORMULA; assign("T4_RESIDUAL_FORMULA", form_for(X[[s]]$meta), envir = globalenv())
    on.exit(assign("T4_RESIDUAL_FORMULA", old, envir = globalenv()))
    make_block(X[[s]]$PM, X[[s]]$meta, idx, frozen)
  }
  blk0 <- block_for("all", seq_len(nrow(X$all$meta)))
  if (is.null(blk0)) fail("H54 design not estimable on the full set")
  ops0 <- ops_for(blk0, frozen)
  tst <- blk0$coverage[cohort_testable == TRUE, feature_id]
  set.seed(SEED); perms <- draw_lineage_perms(N_PERM, tst, frozen$lineage, features)
  okp <- eligible_pairs(blk0, ops0, frozen)
  lk <- t4_stage_leakage(list(blk0))
  gates <- data.table(n = blk0$n, design_columns = blk0$n_par, residual_df = blk0$df, programs_testable = length(tst),
                      coord_eligible_k1 = sum(pr$is_coord & okp), c3_eligible = sum(pr$is_c3 & okp),
                      G2_pass = blk0$n >= 60 && blk0$df >= 30, G3_pass = length(tst) >= 100 && sum(pr$is_coord & okp) >= 2000,
                      max_raw_stage_r2 = max(lk$raw_stage_r2), stop2_fired = max(lk$raw_stage_r2) > 1e-8,
                      panel_symbols = length(X$all$panel), universe_rows_present = X$all$n_universe_present)
  # NC1 and NC3 on the full set
  set.seed(SEED)
  nc <- rbindlist(lapply(seq_len(N_NC1), function(d) {
    zz <- t4_one_pass(list(blk0), pr, list(ops0), Zlist = list(t4_permute_genes(blk0$Z)))
    z1 <- zz$z_partial_k1; z3 <- zz$z_partial_k3
    o1 <- s1_frac(z1, pr$is_coord)[["s1"]]; o3 <- s1_frac(z3, pr$is_c3)[["s1"]]
    q1 <- p_perm_upper(o1, perm_s1(z1, pr$is_coord, perms, pidx, frozen)); q3 <- p_perm_upper(o3, perm_s1(z3, pr$is_c3, perms, pidx, frozen))
    data.table(draw = d, s1_k1 = o1, s1_k3 = o3, p_perm_k1 = q1, p_perm_k3 = q3, reject = max(q1, q3) < 0.05)
  }))
  w(nc, "H54_NC1_NC3_draws.tsv")
  gates[, `:=`(NC1_mean_s1_k1 = mean(nc$s1_k1), stop3_fired = mean(nc$s1_k1) < 0.44 || mean(nc$s1_k1) > 0.54,
               NC3_rejections = sum(nc$reject), stop8_fired = sum(nc$reject) > 16)]
  run_set <- function(s) {
    b <- block_for(s, seq_len(nrow(X[[s]]$meta)))
    if (is.null(b)) return(list(stat = data.table(set = s, design_ok = FALSE)))
    zz <- t4_one_pass(list(b), pr, list(ops_for(b, frozen)))
    ev <- evaluate_z(zz$z_partial_k1, zz$z_partial_k3, perms, pidx, frozen, N_PERM_E, all_perms = TRUE)
    set.seed(SEED); idx <- draw_stratified_boot(N_BOOT, X[[s]]$meta$fibrosis_stage)
    cache <- new.env()
    bd <- rbindlist(lapply(seq_len(nrow(idx)), function(k) {
      bb <- block_for(s, idx[k, ])
      if (is.null(bb)) return(data.table(b = k, design_ok = FALSE))
      z <- t4_one_pass(list(bb), pr, list(ops_for(bb, frozen, X[[s]]$panel, cache)))
      e <- evaluate_z(z$z_partial_k1, z$z_partial_k3, perms, pidx, frozen, N_PERM_E)
      data.table(b = k, design_ok = TRUE, s1_k1 = e$s1_k1, s1_k3 = e$s1_k3, e_k1 = e$e_k1, e_k3 = e$e_k3)
    }), fill = TRUE)
    ok <- bd[design_ok == TRUE]
    inf_obs <- s1_frac(zz$z_partial_k1, pr$is_infnull); inf_null <- perm_s1(zz$z_partial_k1, pr$is_infnull, perms, pidx, frozen)
    st <- data.table(set = s, role = ROLE[[s]], design_ok = TRUE, n = b$n, residual_df = b$df, n_c1 = ev$n_c1, n_c3 = ev$n_c3,
                     s1_k1 = ev$s1_k1, s1_k3 = ev$s1_k3, null_mean_k1_200 = ev$null_mean_k1, null_mean_k3_200 = ev$null_mean_k3,
                     null_q95_k1 = quantile(ev$null_k1, 0.95, names = FALSE), null_q95_k3 = quantile(ev$null_k3, 0.95, names = FALSE),
                     e_k1 = ev$e_k1, e_k3 = ev$e_k3,
                     p_perm_k1 = p_perm_upper(ev$s1_k1, ev$null_k1), p_perm_k3 = p_perm_upper(ev$s1_k3, ev$null_k3),
                     b_total = nrow(bd), b_ok = nrow(ok),
                     p_boot_k1 = (1 + sum(ok$e_k1 <= 0)) / (nrow(ok) + 1), p_boot_k3 = (1 + sum(ok$e_k3 <= 0)) / (nrow(ok) + 1),
                     e_k1_ci_lo = quantile(ok$e_k1, 0.025, names = FALSE), e_k1_ci_hi = quantile(ok$e_k1, 0.975, names = FALSE),
                     e_k3_ci_lo = quantile(ok$e_k3, 0.025, names = FALSE), e_k3_ci_hi = quantile(ok$e_k3, 0.975, names = FALSE),
                     mde_k1 = MDE_MULT * sd(ok$e_k1), mde_k3 = MDE_MULT * sd(ok$e_k3),
                     base_noncoord = ev$base_noncoord, s1_k1_signmatch = ev$s1_k1_signmatch,
                     NC2_s1_k1 = inf_obs[["s1"]], NC2_null_mean = mean(inf_null), NC2_p = p_perm_upper(inf_obs[["s1"]], inf_null))
    st[, p_star := max(p_perm_k1, p_perm_k3, p_boot_k1, p_boot_k3)]
    list(stat = st, z = data.table(pair_index = pr$pair_index, program_a = pr$program_a, program_b = pr$program_b,
                                   z_k1 = zz$z_partial_k1, z_k3 = zz$z_partial_k3), boot = bd)
  }
  R <- list(all = run_set("all"))
  can_pass <- with(gates, G2_pass && G3_pass && !stop2_fired && !stop3_fired && !stop8_fired) &&
    isTRUE(R$all$stat$b_ok >= 0.95 * N_BOOT)
  R$all$stat[, `:=`(can_pass = can_pass, rejected = can_pass && p_star < 0.05, NC2_passes_test = NC2_p < 0.05)]
  opened_sec <- isTRUE(R$all$stat$rejected) || H54_DESIGN == "h54"
  if (opened_sec) { R$F0F1 <- run_set("F0F1"); R$F0F1$stat[, rejected := p_star < 0.05] }
  if ("without_F4" %in% names(SETS)) R$without_F4 <- run_set("without_F4")
  for (s in names(R)) { if (!is.null(R[[s]]$z)) fwrite(R[[s]]$z, file.path(OUT, paste0("H54_edge_z_", s, ".tsv.gz")), sep = "\t")
                        if (!is.null(R[[s]]$boot)) fwrite(R[[s]]$boot, file.path(OUT, paste0("H54_bootstrap_", s, ".tsv.gz")), sep = "\t") }
  w(gates, "H54_gates.tsv"); H <- rbindlist(lapply(R, `[[`, "stat"), fill = TRUE); w(H, "H54_set_statistics.tsv")
  print(gates); print(H)
}
w(qc, "input_qc.tsv")
ss_finish(OUT)
