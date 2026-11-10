# Independent re-derivation of Figure S4D (13-pair stage screen: 6 testable, q = 1, 7 untestable; 69 donors;
# 55,402 groups, 349 fitted, 55,053 untestable; 27 BH families, Bonferroni over 9) and Figure S4E
# (26 planned program x pair models, 22 estimable, 2 supported in 58 donors, betas and q values).
#
# Written from METHODS.md lines 387-391 and the S4D/S4E legends, starting from the raw inputs:
#   LIANA donor-collapsed scores, the 2026-10-14 source-diagnosis donor roster, donor-level program scores.
# Methods do not state the stage-count gate or the stage-model form; those two choices follow the values
# the S4D source table reports (>= 30 donors and >= 10 donors in each of the four stages; random dataset
# intercept with GSE174748 and GSE189600 pooled, log10 sender and receiver cell counts, REML,
# Satterthwaite df). HC3 is computed here from the QR decomposition, not from the producer's code.
#
# Usage: Rscript --vanilla rederive_s4de.R OUTDIR
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest)})
set.seed(20260923)  # no random numbers are drawn; recorded for completeness
args <- commandArgs(trailingOnly = TRUE); out <- args[1]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
f_ccc <- file.path(root, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/all_donor_lr_scores_v2_dc.tsv.gz")
f_roster <- file.path(root, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv")
f_prog <- file.path(root, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv")
f_s4d <- file.path(root, "figures/main/fig4_singlecell_programs/source_tables/current_candidate/source_qualified_ccc_20261014/fixed_13_current_stage_context.tsv")
f_s4e <- file.path(root, "figures/main/fig4_singlecell_programs/source_tables/current_candidate/source_qualified_ccc_20261014/fig4d_program_ccc_associations.tsv")
f_full <- file.path(root, "figures/main/fig4_singlecell_programs/source_tables/current_candidate/source_qualified_ccc_20261014/full_donor_lr_coarse_family.tsv")
inputs <- c(f_ccc, f_roster, f_prog, f_s4d, f_s4e, f_full)
fwrite(data.table(path = inputs, sha256 = vapply(inputs, function(p) strsplit(system2("sha256sum", p, stdout = TRUE), " +")[[1]][1], "")),
       file.path(out, "input_sha256.tsv"), sep = "\t")
log <- function(...) { cat(sprintf(...), "\n"); flush.console() }

STAGES4 <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
roster <- fread(f_roster)
stopifnot(!anyDuplicated(roster$donor))
roster[, excl := tolower(as.character(exclude)) == "true"]
lr <- fread(f_ccc)
log("LR rows %d; units %d; datasets: %s", nrow(lr), uniqueN(lr$sample), paste(names(table(lr$dataset)), collapse = ","))
log("sample == biological_donor in every row: %s", all(lr$sample == lr$biological_donor))
lr[, `:=`(ct_pair = paste(source, target, sep = "->"), lr_pair = paste(ligand_complex, receptor_complex, sep = "__"),
          score = -log10(pmax(magnitude_rank, 1e-4)))]
n_groups_all <- uniqueN(lr[, .(ct_pair, lr_pair)])
dup <- lr[, .N, by = .(sample, ct_pair, lr_pair)][N > 1]
log("all (ct_pair, lr_pair) groups: %d; duplicated donor x group rows: %d", n_groups_all, nrow(dup))

# ---- join to the source-diagnosis roster (rows without a roster entry are not qualified donors)
units <- unique(lr[, .(sample, dataset)])
units <- merge(units, roster[, .(sample = donor, stage = disease_stage_coarse, excl)], by = "sample", all.x = TRUE)
units[, why := fifelse(is.na(excl), "no_roster_entry", fifelse(excl, "roster_exclude", fifelse(stage %in% STAGES4, "kept", "stage_not_in_4")))]
fwrite(units, file.path(out, "s4d_unit_dispositions.tsv"), sep = "\t")
print(units[, .N, by = .(dataset, why)][order(dataset, why)])
keep <- units[why == "kept", .(sample, stage)]
d <- merge(lr, keep, by = "sample")
d[, stage := factor(stage, levels = STAGES4)]
d[, pool := ifelse(dataset %in% c("GSE174748", "GSE189600"), "Other_small", dataset)]
log("qualified donors %d; datasets %d (%s); pools %d", uniqueN(d$sample), uniqueN(d$dataset),
    paste(sort(unique(d$dataset)), collapse = ","), uniqueN(d$pool))
print(unique(d[, .(sample, stage)])[, .N, by = stage])

# ---- stage-count gate over every group in the native input
cnt <- d[, .(n = .N, nH = sum(stage == "Healthy"), nS = sum(stage == "Steatosis"), nSH = sum(stage == "Steatohepatitis"),
             nC = sum(stage == "Cirrhosis")), by = .(ct_pair, lr_pair)]
grp <- merge(unique(lr[, .(ct_pair, lr_pair)]), cnt, by = c("ct_pair", "lr_pair"), all.x = TRUE)
for (v in c("n", "nH", "nS", "nSH", "nC")) set(grp, which(is.na(grp[[v]])), v, 0L)
grp[, eligible := n >= 30 & nH >= 10 & nS >= 10 & nSH >= 10 & nC >= 10]
log("groups %d; eligible %d; not eligible %d", nrow(grp), sum(grp$eligible), sum(!grp$eligible))

# ---- stage models for every eligible group
el <- grp[eligible == TRUE, .(ct_pair, lr_pair)]
t0 <- Sys.time()
del <- split(d[el, on = .(ct_pair, lr_pair)], by = c("ct_pair", "lr_pair"), keep.by = TRUE)
fits <- rbindlist(lapply(seq_len(nrow(el)), function(i) {
  g <- copy(del[[paste(el$ct_pair[i], el$lr_pair[i], sep = ".")]])
  g[, `:=`(ls = log10(pmax(n_source_cells, 1)), lt = log10(pmax(n_target_cells, 1)))]
  m <- try(lmerTest::lmer(score ~ stage + ls + lt + (1 | pool), data = g, REML = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) return(data.table(ct_pair = el$ct_pair[i], lr_pair = el$lr_pair[i], fit_ok = FALSE))
  co <- coef(summary(m))
  rr <- grep("^stage", rownames(co))
  data.table(ct_pair = el$ct_pair[i], lr_pair = el$lr_pair[i], fit_ok = TRUE, term = rownames(co)[rr],
             est = co[rr, "Estimate"], se = co[rr, "Std. Error"], df = co[rr, "df"], p = co[rr, "Pr(>|t|)"],
             n = nrow(g), n_pools = uniqueN(g$pool), singular = isSingular(m))
}), fill = TRUE)
log("fits took %.1f s; groups fitted ok %d; stage coefficients %d; singular fits %d",
    as.numeric(Sys.time() - t0, units = "secs"), uniqueN(fits[fit_ok == TRUE, .(ct_pair, lr_pair)]),
    nrow(fits[fit_ok == TRUE]), uniqueN(fits[fit_ok == TRUE & singular == TRUE, .(ct_pair, lr_pair)]))
ff <- fits[fit_ok == TRUE]
ff[, q_bh := p.adjust(p, "BH"), by = .(ct_pair, term)]
n_fam_bh <- uniqueN(ff[, .(ct_pair, term)]); n_ct <- uniqueN(ff$ct_pair)
ff[, q := pmin(1, q_bh * n_ct)]
log("BH families (ct_pair x term): %d; Bonferroni multiplier (ct_pairs): %d; min q any term %.4g", n_fam_bh, n_ct, min(ff$q))
fwrite(ff, file.path(out, "s4d_full_family_independent.tsv"), sep = "\t")
fwrite(grp, file.path(out, "s4d_group_gate_independent.tsv"), sep = "\t")

# compare with the producer's full family
pf <- fread(f_full)
pf[, term := sub("^disease_stage_coarse", "stage", term)]
cmp <- merge(ff, pf[, .(ct_pair, lr_pair, term, Estimate, pval, q_bonferroni_family)], by = c("ct_pair", "lr_pair", "term"), all = TRUE)
log("full family rows: independent %d, producer %d, joined %d; max |est diff| %.2e; max |p diff| %.2e; max |q diff| %.2e",
    nrow(ff), nrow(pf), nrow(cmp), max(abs(cmp$est - cmp$Estimate), na.rm = TRUE), max(abs(cmp$p - cmp$pval), na.rm = TRUE),
    max(abs(cmp$q - cmp$q_bonferroni_family), na.rm = TRUE))
fwrite(cmp, file.path(out, "s4d_full_family_vs_producer.tsv"), sep = "\t")

# ---- the fixed 13-pair roster
r13 <- data.table(ct_pair = c("Cholangiocytes->Endothelial cells", "Cholangiocytes->Fibroblasts", "Cholangiocytes->Hepatocytes",
                              "Endothelial cells->Hepatocytes", "Endothelial cells->Macrophages", "Fibroblasts->Hepatocytes",
                              "Fibroblasts->Hepatocytes", "Hepatocytes->Endothelial cells", "Hepatocytes->Endothelial cells",
                              "Hepatocytes->Endothelial cells", "Hepatocytes->Fibroblasts", "Hepatocytes->Fibroblasts",
                              "Hepatocytes->Macrophages"),
                  lr_pair = c("C3__NRP1", "COL4A1__ITGA1_ITGB1", "ADAM10__IL6R", "ST6GAL1__EGFR", "COL4A2__CD44", "LAMC3__DAG1",
                              "THBS1__LRP5", "C4BPA__BMPR2", "CDH1__PTPRM", "COL18A1__KDR", "C4BPA__BMPR2", "COL18A1__ITGA1_ITGB1",
                              "C3__CR1"))
s13 <- merge(r13, grp, by = c("ct_pair", "lr_pair"), all.x = TRUE)
s13 <- merge(s13, dcast(ff, ct_pair + lr_pair ~ term, value.var = "q"), by = c("ct_pair", "lr_pair"), all.x = TRUE)
fwrite(s13, file.path(out, "s4d_fixed13_independent.tsv"), sep = "\t")
print(s13)
log("S4D fixed 13: testable %d; untestable %d; max donors in a testable pair %d", sum(s13$eligible), sum(!s13$eligible), max(s13$n[s13$eligible]))
p13 <- fread(f_s4d)
c13 <- merge(s13[, .(ct_pair, lr_pair, eligible, n, q_SH = stageSteatohepatitis)], p13[, .(ct_pair, lr_pair, eligible_stage_counts, n_donors, dc_q)],
             by = c("ct_pair", "lr_pair"))
log("S4D fixed 13 vs producer: eligibility agrees %d/13; donor n agrees %d/13; SH q agrees %d/%d",
    sum(c13$eligible == c13$eligible_stage_counts), sum(c13$n == c13$n_donors), sum(abs(c13$q_SH - c13$dc_q) < 1e-8, na.rm = TRUE), sum(!is.na(c13$dc_q)))

# ---- S4E: 2 programs x 13 pairs from raw donor program scores
focus <- c("hotspot_hepatocytes_f05c535ae5bbc0b9", "hotspot_hepatocytes_48f39dd4d817a10e")
ps <- fread(f_prog)[program_uid %in% focus, .(donor, program_uid, pscore = score, pdataset = dataset, old_stage = disease_stage_coarse,
                                             old_excl = exclude_stage_analysis)]
ps <- merge(ps, roster[, .(donor, stage = disease_stage_coarse, excl)], by = "donor", all.x = TRUE)
log("program score donors per program: %s", paste(ps[, .N, by = program_uid]$N, collapse = ","))
miss <- roster[excl == FALSE & disease_stage_coarse %in% STAGES4[1:3] & !donor %in% ps$donor]
log("roster-qualified H/S/SH donors absent from the program score file: %d (%s)", nrow(miss), paste(miss$donor, collapse = ","))
print(ps[program_uid == focus[2], .N, by = .(old_stage, old_excl, stage, excl)][order(old_stage, stage)])
ps <- ps[!is.na(excl) & excl == FALSE & stage %in% STAGES4[1:3]]
lr13 <- merge(lr, r13, by = c("ct_pair", "lr_pair"))
m <- merge(lr13[, .(donor = sample, ldataset = dataset, ct_pair, lr_pair, score, n_source_cells, n_target_cells)], ps,
           by = "donor", allow.cartesian = TRUE)
stopifnot(all(m$ldataset == m$pdataset))
m[, `:=`(stage = factor(stage, levels = STAGES4[1:3]), dataset = factor(ldataset),
         ls = log10(pmax(n_source_cells, 1)), lt = log10(pmax(n_target_cells, 1)))]

hc3 <- function(fit) {
  X <- model.matrix(fit); e <- residuals(fit)
  Q <- qr(X); h <- rowSums(qr.Q(Q)^2)
  XtXi <- chol2inv(qr.R(Q)); piv <- Q$pivot; XtXi[piv, piv] <- XtXi
  meat <- crossprod(X * (e / (1 - h)))
  XtXi %*% meat %*% XtXi
}
one <- function(g) {
  g <- g[is.finite(score) & is.finite(pscore)]
  res <- data.table(n = uniqueN(g$donor), n_ds = uniqueN(g$dataset), beta = NA_real_, p = NA_real_, p_hc3 = NA_real_, estimable = FALSE)
  if (res$n < 30) return(res)
  g[, pz := (pscore - mean(pscore)) / sd(pscore)]
  fit <- lm(score ~ pz + stage + dataset + ls + lt, data = g)
  if (fit$rank < ncol(model.matrix(fit))) return(res)
  b <- coef(fit)[["pz"]]; V <- hc3(fit); j <- which(colnames(model.matrix(fit)) == "pz")
  res[, `:=`(beta = b, p = coef(summary(fit))["pz", 4], p_hc3 = 2 * pt(abs(b / sqrt(V[j, j])), fit$df.residual, lower.tail = FALSE),
             estimable = TRUE)]
  res
}
e26 <- m[, one(.SD), by = .(program_uid, ct_pair, lr_pair)]
e26 <- merge(CJ(program_uid = focus, k = seq_len(13))[, `:=`(ct_pair = r13$ct_pair[k], lr_pair = r13$lr_pair[k])][, k := NULL],
             e26, by = c("program_uid", "ct_pair", "lr_pair"), all.x = TRUE)
e26[is.na(estimable), estimable := FALSE]
e26[, `:=`(q = p.adjust(p, "BH", n = 26), q_hc3 = p.adjust(p_hc3, "BH", n = 26),
           q22 = p.adjust(p, "BH"), q22_hc3 = p.adjust(p_hc3, "BH"))]
e26[, state := fifelse(!estimable, "untestable", fifelse(q < .05 & q_hc3 < .05, "supported", "indeterminate"))]
e26[, state22 := fifelse(!estimable, "untestable", fifelse(q22 < .05 & q22_hc3 < .05, "supported", "indeterminate"))]
log("S4E: rows %d; estimable %d; supported (BH n=26) %d; supported (BH over estimable only) %d",
    nrow(e26), sum(e26$estimable), sum(e26$state == "supported"), sum(e26$state22 == "supported"))
print(e26[state == "supported" | state22 == "supported"])

# LODO direction for supported rows
lodo <- rbindlist(lapply(which(e26$state == "supported"), function(i) {
  g0 <- m[program_uid == e26$program_uid[i] & ct_pair == e26$ct_pair[i] & lr_pair == e26$lr_pair[i]]
  rbindlist(lapply(sort(unique(as.character(g0$dataset))), function(ds) {
    g <- g0[as.character(dataset) != ds]; g[, pz := (pscore - mean(pscore)) / sd(pscore)]
    g[, `:=`(stage = droplevels(stage), dataset = droplevels(dataset))]
    f <- if (uniqueN(g$dataset) > 1) lm(score ~ pz + stage + dataset + ls + lt, g) else lm(score ~ pz + stage + ls + lt, g)
    data.table(program_uid = e26$program_uid[i], ct_pair = e26$ct_pair[i], lr_pair = e26$lr_pair[i], left_out = ds,
               n = uniqueN(g$donor), beta = coef(f)[["pz"]], beta_full = e26$beta[i])
  }))
}))
print(lodo)
fwrite(lodo, file.path(out, "s4e_lodo_independent.tsv"), sep = "\t")
fwrite(e26, file.path(out, "s4e_26_independent.tsv"), sep = "\t")
pe <- fread(f_s4e)
ce <- merge(e26, pe[, .(program_uid, ct_pair, lr_pair, n_donors, beta_p = beta, qvalue, hc3_qvalue, testable, evidence_state)],
            by = c("program_uid", "ct_pair", "lr_pair"))
log("S4E vs producer: estimable agrees %d/26; state agrees %d/26; donor n agrees %d/26; max |beta diff| %.2e; max |q diff| %.2e; max |hc3 q diff| %.2e",
    sum(ce$estimable == ce$testable), sum(ce$state == ce$evidence_state), sum(ce$n == ce$n_donors, na.rm = TRUE),
    max(abs(ce$beta - ce$beta_p), na.rm = TRUE), max(abs(ce$q - ce$qvalue), na.rm = TRUE), max(abs(ce$q_hc3 - ce$hc3_qvalue), na.rm = TRUE))
fwrite(ce, file.path(out, "s4e_vs_producer.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
