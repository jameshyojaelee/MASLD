#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
DIG1 <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
GDIR <- file.path(DIG1, "out_genotype_table_20261118")
OUT <- file.path(DIG1, "out_genotype_state_20261118")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Load inputs
qc <- fread(file.path(GDIR, "site_qc.tsv"))
qc_cohort <- fread(file.path(GDIR, "site_qc_by_cohort.tsv"))
gwide10 <- fread(file.path(GDIR, "genotype_calls_wide.tsv"))
gwide20 <- fread(file.path(GDIR, "genotype_calls_wide_dp20.tsv"))
junc <- fread(file.path(DIG1, "out_hsd17b13_junctions/hsd17b13_skip6.tsv"))
scores <- fread(file.path(DIG1, "out_bulk_state_histology/sample_scores.tsv"))
old_sg <- fread(file.path(DIG1, "out_genotype_state/sample_scores_genotypes.tsv"))
meta_raw_rds <- readRDS(file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds"))
meta_dt <- as.data.table(meta_raw_rds)
anc <- fread(file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv"))

to_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
for (col in c("pass_qc", "pass_qc_dp20", "pass_qc_eur_hwe05", "pass_qc_eur_dp20_hwe05", "pass_qc_hwe_1e4_legacy")) {
  qc[, (col) := to_bool(get(col))]
}

# Variant mapping (dosage of target/risk/protective allele per 62_targets_annotated.tsv & 63)
vmap <- data.table(
  variant = c("PNPLA3_148M", "PNPLA3_434E", "TM6SF2_167K", "GCKR_446L", "MTARC1_165T",
              "MBOAT7_rs641738T", "APOE_e4", "APOE_e2", "HSD17B13_ins", "SERPINA1_Z"),
  rsid = c("rs738409", "rs2294918", "rs58542926", "rs1260326", "rs2642438",
           "rs641738", "rs429358", "rs7412", "rs72613567", "rs28929474"),
  flip_2_minus = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
)
qc_audit <- merge(vmap, qc, by = "rsid", sort = FALSE)
fwrite(qc_audit, file.path(OUT, "strict_hwe05_gate_audit.tsv"), sep = "\t")

# Build variant matrices for dp>=10 and dp>=20
build_var_dt <- function(gw) {
  out <- gw[, .(sample_id, dataset)]
  for (i in seq_len(nrow(vmap))) {
    rs <- vmap$rsid[i]; vn <- vmap$variant[i]
    vals <- as.numeric(gw[[rs]])
    if (vmap$flip_2_minus[i]) vals <- 2.0 - vals
    out[, (vn) := vals]
  }
  out
}
v10 <- build_var_dt(gwide10)
v20 <- build_var_dt(gwide20)
setnames(v20, vmap$variant, paste0(vmap$variant, "_dp20"))

# HSD17B13 junction coding (tot >= 30, skip_ratio <= 0.02 -> 0, 0.04..0.27 -> 1, >= 0.50 -> 2)
junc[, tot := skip6 + canon67]
junc[, HSD17B13_junction := ifelse(tot < 30, NA_real_,
                                   ifelse(skip_ratio <= 0.02, 0.0,
                                          ifelse(skip_ratio >= 0.04 & skip_ratio <= 0.27, 1.0,
                                                 ifelse(skip_ratio >= 0.50, 2.0, NA_real_))))]
# Verify exact concordance with prior HSD17B13_junction calls
j_chk <- merge(junc[, .(sample_id, new_j = HSD17B13_junction)], old_sg[, .(sample_id, old_j = HSD17B13_junction)], by = "sample_id")
stopifnot(identical(is.na(j_chk$new_j), is.na(j_chk$old_j)),
          all(j_chk$new_j[!is.na(j_chk$new_j)] == j_chk$old_j[!is.na(j_chk$old_j)]))

junc_qc <- junc[, .(
  n = .N,
  n_called = sum(!is.na(HSD17B13_junction)),
  call_rate = mean(!is.na(HSD17B13_junction)),
  n0 = sum(HSD17B13_junction == 0, na.rm = TRUE),
  n1 = sum(HSD17B13_junction == 1, na.rm = TRUE),
  n2 = sum(HSD17B13_junction == 2, na.rm = TRUE),
  alt_freq = (sum(HSD17B13_junction == 1, na.rm = TRUE) + 2 * sum(HSD17B13_junction == 2, na.rm = TRUE)) / (2 * sum(!is.na(HSD17B13_junction)))
)]
fwrite(junc_qc, file.path(OUT, "hsd17b13_junction_qc.tsv"), sep = "\t")

# Join state scores + g_* + metadata + ancestry PCs
g_cols <- old_sg[, .(sample_id, g_FCAMR, g_TRIM31, g_CXCL10, g_IL32, g_CCL20)]
anc_cols <- anc[, .(sample_id = run, nearest_superpop, ancestry_call, passes_min_sites, cohort_contaminated,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]
df <- Reduce(function(a, b) merge(a, b, by = "sample_id", all.x = TRUE), list(
  scores[, .(sample_id, dataset, group_binary, fibrosis_stage, nas_score, steatosis_grade,
             INFLAMED_HEP, IMMUNE, FIBRO, DUCTULAR, HEPMASS)],
  g_cols,
  meta_dt[, .(sample_id, diagnosis_harmonized, sex_final)],
  v10[, c("sample_id", vmap$variant), with = FALSE],
  v20[, c("sample_id", paste0(vmap$variant, "_dp20")), with = FALSE],
  junc[, .(sample_id, HSD17B13_junction)],
  anc_cols
))
stopifnot(nrow(df) == 774L, sum(is.na(df$PC1)) == 0L)

# Clean metadata fibrosis_stage: set GSE213621 Control (68 non-biopsy/un-staged controls coded 0) to NA
# Note: in meta_matched.rds for GSE213621, diagnosis_harmonized is NA on Disease rows while group_binary is Control/Disease
is_ctrl_213621 <- (df$dataset == "GSE213621") & (
  (!is.na(df$group_binary) & df$group_binary == "Control") |
  (!is.na(df$diagnosis_harmonized) & df$diagnosis_harmonized == "Control")
)
stopifnot(sum(is_ctrl_213621) == 68L)
df[, fibrosis_stage_raw := as.numeric(fibrosis_stage)]
df[, fibrosis_stage_clean := ifelse(is_ctrl_213621, NA_real_, as.numeric(fibrosis_stage))]
fwrite(df, file.path(OUT, "sample_scores_genotypes_pc.tsv"), sep = "\t")

# Save clean metadata audit counts
meta_audit <- df[, .(
  n_total = .N,
  n_fib_raw_nonNA = sum(!is.na(fibrosis_stage_raw)),
  n_fib_clean_nonNA = sum(!is.na(fibrosis_stage_clean)),
  n_gse213621_ctrl_set_NA = sum(dataset == "GSE213621" & !is.na(group_binary) & group_binary == "Control" & !is.na(fibrosis_stage_raw)),
  n_eur = sum(nearest_superpop == "EUR")
), by = dataset]
fwrite(meta_audit, file.path(OUT, "metadata_and_ancestry_audit.tsv"), sep = "\t")

fit_term <- function(d_sub, formula_str, term_name, analysis_label, outcome_name, variant_name, cohort_val = NA_character_) {
  f_in <- as.formula(formula_str)
  vars <- all.vars(f_in)
  cc <- d_sub[complete.cases(d_sub[, ..vars])]
  if (nrow(cc) < 5L) return(NULL)
  rhs <- attr(terms(f_in), "term.labels")
  keep_rhs <- rhs[vapply(rhs, function(tr) {
    v <- gsub("^factor\\((.*)\\)$", "\\1", tr)
    if (v %in% names(cc)) length(unique(cc[[v]])) >= 2L else TRUE
  }, logical(1L))]
  if (!(term_name %in% keep_rhs)) return(NULL)
  f_use <- reformulate(keep_rhs, response = vars[1L])
  m <- tryCatch(lm(f_use, data = cc), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- coef(m)
  if (!(term_name %in% names(b)) || is.na(b[term_name])) return(NULL)
  X <- model.matrix(m)
  if (ncol(X) >= nrow(X)) return(NULL)
  e <- residuals(m)
  h <- hatvalues(m)
  if (any(!is.finite(h)) || any(1 - h <= 1e-8)) return(NULL)
  XtXi <- tryCatch(solve(crossprod(X)), error = function(err) NULL)
  if (is.null(XtXi)) return(NULL)
  V <- XtXi %*% crossprod(X * (e / (1.0 - h))) %*% XtXi
  se_hc3 <- sqrt(diag(V))
  if (!(term_name %in% names(se_hc3)) || !is.finite(se_hc3[term_name]) || se_hc3[term_name] <= 0) return(NULL)
  est <- as.numeric(b[term_name])
  se_v <- as.numeric(se_hc3[term_name])
  dfr <- as.integer(df.residual(m))
  t_v <- est / se_v
  p_v <- 2.0 * pt(-abs(t_v), df = dfr)
  tcrit <- qt(0.975, df = dfr)
  data.table(
    term = variant_name,
    estimate = est,
    se = se_v,
    t_stat = t_v,
    p = p_v,
    ci_low = est - tcrit * se_v,
    ci_high = est + tcrit * se_v,
    n = nobs(m),
    df_resid = dfr,
    analysis = analysis_label,
    outcome = outcome_name,
    variant = variant_name,
    cohort = cohort_val
  )
}

# Variants to test: PNPLA3_148M (primary in DIG1, audit strict HWE status), all legacy-usable variants, and SERPINA1_Z (strict HWE>=0.05)
TEST_VARIANTS <- c("PNPLA3_148M", "PNPLA3_434E", "TM6SF2_167K", "GCKR_446L", "MTARC1_165T", "SERPINA1_Z", "HSD17B13_junction")
OUTCOMES_STATE <- c("INFLAMED_HEP", "IMMUNE", "FIBRO", "DUCTULAR", "HEPMASS", "g_FCAMR", "g_TRIM31", "g_CXCL10", "g_IL32", "g_CCL20")
HIST_OUTCOMES <- c("nas_score", "steatosis_grade", "fibrosis_stage_raw", "fibrosis_stage_clean")

pc4_str <- "PC1 + PC2 + PC3 + PC4"
pc10_str <- paste0("PC", 1:10, collapse = " + ")

qc_lookup <- qc_audit[, .(variant, min_cohort_hwe_p, pass_qc, min_cohort_hwe_p_dp20, pass_qc_dp20, min_cohort_hwe_p_eur, pass_qc_eur_hwe05, pass_qc_hwe_1e4_legacy)]

model_rows <- list()

for (vn in TEST_VARIANTS) {
  v_dp20 <- if (vn == "HSD17B13_junction") "HSD17B13_junction" else paste0(vn, "_dp20")

  # 1) State outcomes across raw, PC4, PC10, dp>=20, EUR-only, histology-raw vs clean, and disease-only
  for (outc in OUTCOMES_STATE) {
    # G_raw
    d0 <- df[!is.na(get(outc)) & !is.na(get(vn))]
    model_rows[[length(model_rows) + 1L]] <- fit_term(d0, sprintf("%s ~ %s + factor(dataset) + factor(sex_final)", outc, vn), vn, "G_raw", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(d0, sprintf("%s ~ %s + factor(dataset) + factor(sex_final) + %s", outc, vn, pc4_str), vn, "G_pc4", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(d0, sprintf("%s ~ %s + factor(dataset) + factor(sex_final) + %s", outc, vn, pc10_str), vn, "G_pc10", outc, vn)

    # dp>=20
    d20 <- df[!is.na(get(outc)) & !is.na(get(v_dp20))]
    model_rows[[length(model_rows) + 1L]] <- fit_term(d20, sprintf("%s ~ %s + factor(dataset) + factor(sex_final)", outc, v_dp20), v_dp20, "G_dp20_raw", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(d20, sprintf("%s ~ %s + factor(dataset) + factor(sex_final) + %s", outc, v_dp20, pc4_str), v_dp20, "G_dp20_pc4", outc, vn)

    # EUR-only
    deur <- d0[nearest_superpop == "EUR"]
    model_rows[[length(model_rows) + 1L]] <- fit_term(deur, sprintf("%s ~ %s + factor(dataset) + factor(sex_final)", outc, vn), vn, "G_eur_raw", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(deur, sprintf("%s ~ %s + factor(dataset) + factor(sex_final) + %s", outc, vn, pc4_str), vn, "G_eur_pc4", outc, vn)

    # Histology-adjusted (raw vs clean GSE213621 Control=NA)
    dh_raw <- d0[!is.na(fibrosis_stage_raw)]
    model_rows[[length(model_rows) + 1L]] <- fit_term(dh_raw, sprintf("%s ~ %s + fibrosis_stage_raw + factor(dataset) + factor(sex_final)", outc, vn), vn, "G_plus_histology_raw_meta", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(dh_raw, sprintf("%s ~ %s + fibrosis_stage_raw + factor(dataset) + factor(sex_final) + %s", outc, vn, pc4_str), vn, "G_plus_histology_raw_meta_pc4", outc, vn)

    dh_cln <- d0[!is.na(fibrosis_stage_clean)]
    model_rows[[length(model_rows) + 1L]] <- fit_term(dh_cln, sprintf("%s ~ %s + fibrosis_stage_clean + factor(dataset) + factor(sex_final)", outc, vn), vn, "G_plus_histology_clean_meta", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(dh_cln, sprintf("%s ~ %s + fibrosis_stage_clean + factor(dataset) + factor(sex_final) + %s", outc, vn, pc4_str), vn, "G_plus_histology_clean_meta_pc4", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(d0[!is.na(fibrosis_stage_clean) & nearest_superpop == "EUR"],
                                                      sprintf("%s ~ %s + fibrosis_stage_clean + factor(dataset) + factor(sex_final) + %s", outc, vn, pc4_str),
                                                      vn, "G_plus_histology_clean_meta_eur_pc4", outc, vn)

    # Disease-only
    ddis <- d0[group_binary == "Disease"]
    model_rows[[length(model_rows) + 1L]] <- fit_term(ddis, sprintf("%s ~ %s + factor(dataset) + factor(sex_final)", outc, vn), vn, "G_disease_only", outc, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(ddis, sprintf("%s ~ %s + factor(dataset) + factor(sex_final) + %s", outc, vn, pc4_str), vn, "G_disease_only_pc4", outc, vn)
  }

  # 2) Positive-control histology models
  for (hout in HIST_OUTCOMES) {
    dh <- df[!is.na(get(hout)) & !is.na(get(vn))]
    ds_cov <- if (uniqueN(dh$dataset) > 1L) "+ factor(dataset)" else ""
    model_rows[[length(model_rows) + 1L]] <- fit_term(dh, sprintf("%s ~ %s %s + factor(sex_final)", hout, vn, ds_cov), vn, "control_histology_raw", hout, vn)
    model_rows[[length(model_rows) + 1L]] <- fit_term(dh, sprintf("%s ~ %s %s + factor(sex_final) + %s", hout, vn, ds_cov, pc4_str), vn, "control_histology_pc4", hout, vn)
  }

  # 3) Per-cohort INFLAMED_HEP (raw, PC4, EUR_raw, EUR_pc4)
  for (ds in sort(unique(df$dataset))) {
    dsub <- df[dataset == ds & !is.na(INFLAMED_HEP) & !is.na(get(vn))]
    if (nrow(dsub) >= 15L && uniqueN(dsub[[vn]]) >= 2L) {
      sx <- if (uniqueN(dsub$sex_final) > 1L) "+ factor(sex_final)" else ""
      model_rows[[length(model_rows) + 1L]] <- fit_term(dsub, sprintf("INFLAMED_HEP ~ %s %s", vn, sx), vn, "per_cohort", "INFLAMED_HEP", vn, ds)
      model_rows[[length(model_rows) + 1L]] <- fit_term(dsub, sprintf("INFLAMED_HEP ~ %s %s + %s", vn, sx, pc4_str), vn, "per_cohort_pc4", "INFLAMED_HEP", vn, ds)
      de_c <- dsub[nearest_superpop == "EUR"]
      if (nrow(de_c) >= 15L && uniqueN(de_c[[vn]]) >= 2L) {
        sxe <- if (uniqueN(de_c$sex_final) > 1L) "+ factor(sex_final)" else ""
        model_rows[[length(model_rows) + 1L]] <- fit_term(de_c, sprintf("INFLAMED_HEP ~ %s %s + %s", vn, sxe, pc4_str), vn, "per_cohort_eur_pc4", "INFLAMED_HEP", vn, ds)
      }
    }
  }
}

mods <- rbindlist(Filter(Negate(is.null), model_rows), fill = TRUE)
mods <- merge(mods, qc_lookup, by = "variant", all.x = TRUE, sort = FALSE)
mods[, strict_hwe05_gate := ifelse(variant == "HSD17B13_junction", "pass_splice_junction",
                                   ifelse(pass_qc, "pass_hwe05", "fail_hwe05_sensitivity"))]
fwrite(mods, file.path(OUT, "genotype_models.tsv"), sep = "\t")

# Primary PNPLA3 summary across all sensitivity axes
prim_pnpla3 <- mods[variant == "PNPLA3_148M" & outcome == "INFLAMED_HEP" &
                    analysis %in% c("G_raw", "G_pc4", "G_pc10", "G_dp20_raw", "G_dp20_pc4",
                                    "G_eur_raw", "G_eur_pc4",
                                    "G_plus_histology_raw_meta", "G_plus_histology_raw_meta_pc4",
                                    "G_plus_histology_clean_meta", "G_plus_histology_clean_meta_pc4",
                                    "G_plus_histology_clean_meta_eur_pc4",
                                    "G_disease_only", "G_disease_only_pc4")]
fwrite(prim_pnpla3, file.path(OUT, "primary_pnpla3.tsv"), sep = "\t")

# Secondary Holm comparisons:
# 1) Legacy secondary family (PNPLA3_434E, TM6SF2_167K, GCKR_446L, MTARC1_165T, HSD17B13_junction)
# 2) Strict HWE>=0.05 family (SERPINA1_Z + HSD17B13_junction)
sec_vars_legacy <- c("PNPLA3_434E", "TM6SF2_167K", "GCKR_446L", "MTARC1_165T", "HSD17B13_junction")
sec_vars_strict <- c("SERPINA1_Z", "HSD17B13_junction")
sec_list <- list()
for (an in c("G_raw", "G_pc4", "G_pc10", "G_dp20_pc4", "G_eur_pc4")) {
  sub_leg <- copy(mods[outcome == "INFLAMED_HEP" & analysis == an & variant %in% sec_vars_legacy])
  sub_leg[, `:=`(holm_family = "legacy_5var", p_holm = p.adjust(p, method = "holm"))]
  sub_str <- copy(mods[outcome == "INFLAMED_HEP" & analysis == an & variant %in% sec_vars_strict])
  sub_str[, `:=`(holm_family = "strict_hwe05_only", p_holm = p.adjust(p, method = "holm"))]
  sec_list[[paste0(an, "_leg")]] <- sub_leg
  sec_list[[paste0(an, "_str")]] <- sub_str
}
fwrite(rbindlist(sec_list, fill = TRUE), file.path(OUT, "secondary_holm.tsv"), sep = "\t")

# Per-cohort heterogeneity across raw, pc4, and eur_pc4
het_rows <- list()
for (an in c("per_cohort", "per_cohort_pc4", "per_cohort_eur_pc4")) {
  for (vn in TEST_VARIANTS) {
    sub <- mods[analysis == an & variant == vn & outcome == "INFLAMED_HEP"]
    k <- nrow(sub)
    if (k < 2L) next
    w <- 1.0 / (sub$se^2)
    fe <- sum(w * sub$estimate) / sum(w)
    fe_se <- sqrt(1.0 / sum(w))
    Q <- sum(w * (sub$estimate - fe)^2)
    I2 <- if (Q > k - 1L) (Q - (k - 1L)) / Q else 0.0
    het_rows[[paste0(an, "_", vn)]] <- data.table(analysis = an, variant = vn, k = k, fe = fe, fe_se = fe_se,
                                                  fe_p = 2 * pnorm(-abs(fe / fe_se)), Q = Q,
                                                  Q_p = pchisq(Q, df = k - 1L, lower.tail = FALSE), I2 = I2)
  }
}
fwrite(rbindlist(het_rows, fill = TRUE), file.path(OUT, "per_cohort_heterogeneity.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("DONE 63_genotype_state_20261118\n")
