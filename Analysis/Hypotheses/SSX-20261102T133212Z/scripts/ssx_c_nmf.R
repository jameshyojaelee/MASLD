# SSX-C (PRESPEC_SSX_C.md; DEVIATIONS_SSX.md D-C-1). Fixed k6 NMF axis activity by stage in SteatoSITE.
# Release transport (fixed_nmf_count_transport.R): pairwise TMM to the stored reference, log2 CPM, cohort shift,
# row minimum, clip at 0, NNLS on the 4,517 retained rows of basis k6. Gate: the same code path reproduces the release
# discovery scores (P1-P5, max |difference| < 1e-6) before SteatoSITE is scored.
# SteatoSITE shift per batch b and gene g: mean log2 CPM over the batch's included libraries minus the mean over the
# 774 discovery units of their shifted log2 CPM. Sensitivity: zero shift.
# Testable axes fixed before stage is read (D-C-1). Model per axis: z (within 493 staged units) ~ batch + inferred_sex
# + stage4 (F0-F1 reference); family: bin model F3-F4 vs F0-F1 coefficient, one-sided in the discovery direction,
# 1,000 permutations of stage within batch x sex (seed 20261104), Holm over testable axes.
# Usage: Rscript ssx_c_nmf.R --out <new dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(nnls)})
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k) { i <- match(k, a); if (is.na(i)) NA_character_ else a[i + 1] }
out <- opt("--out"); stopifnot(!is.na(out), !dir.exists(out)); dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SS <- file.path(ROOT, "Analysis/Hypotheses/SS-20261025T135237Z")
SSX <- file.path(ROOT, "Analysis/Hypotheses/SSX-20261102T133212Z")
CAND <- file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1")
stopifnot(system(paste("cd", SSX, "&& sha256sum -c --quiet PRESPEC_SSX.sha256")) == 0)
source(file.path(ROOT, "scripts/analysis/histology_anchored_continuum/molecular_layers/fixed_nmf_count_transport.R"))
as_count_matrix <- function(x) { if (inherits(x, "DGEList")) x <- x$counts; if (inherits(x, "Matrix")) x <- as.matrix(x); x }
PR <- paste0("P", 1:5)

protocol <- readRDS(file.path(CAND, "bulk_execution_sources/revised_fixed_nmf_projection_r7_review/revised_fixed_nmf_projection_protocol_with_transport_checks_r2.rds"))
state <- protocol$historical_preprocessing
w <- protocol$basis[["6"]]; stopifnot(nrow(w) == 4517L, qr(w)$rank == 6L)
ref_ids <- protocol$TMM_reference_gene_ids
transport <- function(counts, cohorts, st) {
  sizes <- colSums(counts)
  f <- fixed_reference_tmm(counts[ref_ids, , drop = FALSE], sizes, protocol$TMM_reference_counts,
                           protocol$TMM_reference_library_size, protocol$TMM_reference_norm_factor)
  list(sizes = sizes, factors = f,
       x = fixed_nmf_transport_expression(counts[rownames(w), , drop = FALSE], sizes, f, st, cohorts)) }
logcpm <- function(counts, sizes, f) log2(sweep(counts[rownames(w), , drop = FALSE] + .5, 2L, sizes * f + 1, "/") * 1e6)
recon_rho <- function(x, h) vapply(seq_len(ncol(x)), function(i) cor(x[, i], drop(w %*% h[, i]), method = "spearman"), numeric(1))

# ---- gate: discovery 774 units through this code path
dc <- as_count_matrix(readRDS(file.path(CAND, "verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")))
roster <- fread(file.path(CAND, "identity_resource_qualified/bulk_execution_inputs/primary_bulk_metadata.tsv"))
roster <- roster[match(colnames(dc), sample_id)]; stopifnot(!anyNA(roster$dataset))
TD <- transport(dc, roster$dataset, state)
hd <- fixed_nmf_nnls_scores(w, TD$x)
REL <- fread(file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/fixed-nmf-corrected-scores-20261009T141312Z/fixed_nmf_scores_full_stage_k6.tsv"))
REL <- REL[match(colnames(hd), sample_id)]
gate_diff <- max(abs(as.matrix(REL[, ..PR]) - t(hd[PR, ])))
fwrite(data.table(n_units = ncol(hd), max_abs_diff_P1_P5 = gate_diff, pass = gate_diff < 1e-6), file.path(out, "gate_reproduction.tsv"), sep = "\t")
cat("GATE max |diff|:", gate_diff, "\n"); stopifnot(gate_diff < 1e-6)
disc_shifted_mean <- rowMeans(logcpm(dc, TD$sizes, TD$factors) - state$cohort_batch_subtraction[rownames(w), roster$dataset, drop = FALSE])
disc_recon <- recon_rho(TD$x, hd); rm(dc); invisible(gc())

# ---- SteatoSITE 515 included libraries
U <- fread(file.path(SS, "out_s2_qc/units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
sc <- as_count_matrix(readRDS(file.path(SS, "out_s2_qc/counts_all_genes.rds")))[, U$sample_alias]
stopifnot(all(ref_ids %in% rownames(sc)), all(rownames(w) %in% rownames(sc)))
sizes_s <- colSums(sc)
f_s <- fixed_reference_tmm(sc[ref_ids, , drop = FALSE], sizes_s, protocol$TMM_reference_counts,
                           protocol$TMM_reference_library_size, protocol$TMM_reference_norm_factor)
L <- logcpm(sc, sizes_s, f_s)
bcol <- paste0("SS_", U$batch)
shift <- sapply(split(seq_len(ncol(L)), bcol), function(i) rowMeans(L[, i, drop = FALSE]) - disc_shifted_mean)
st1 <- state; st1$cohort_batch_subtraction <- cbind(state$cohort_batch_subtraction[rownames(w), , drop = FALSE], shift[rownames(w), , drop = FALSE])
st0 <- state; st0$cohort_batch_subtraction <- cbind(state$cohort_batch_subtraction[rownames(w), , drop = FALSE], SS0 = 0)
x1 <- fixed_nmf_transport_expression(sc[rownames(w), ], sizes_s, f_s, st1, bcol)
x0 <- fixed_nmf_transport_expression(sc[rownames(w), ], sizes_s, f_s, st0, rep("SS0", ncol(sc)))
h1 <- fixed_nmf_nnls_scores(w, x1); h0 <- fixed_nmf_nnls_scores(w, x0)
ss_recon <- recon_rho(x1, h1)

# ---- D-C-1 testable axes (no stage read yet)
big <- names(which(table(bcol) >= 20))
AX <- rbindlist(lapply(PR, function(p) {
  r <- vapply(big, function(b) { i <- bcol == b; suppressWarnings(cor(h1[p, i], h0[p, i], method = "spearman")) }, numeric(1))
  data.table(program = p, n_batches_ge20 = length(big), min_within_batch_rho_shift_vs_zero = min(r),
             fraction_zero_scores = mean(h1[p, ] == 0)) }))
AX[, testable := is.finite(min_within_batch_rho_shift_vs_zero) & min_within_batch_rho_shift_vs_zero >= 0.9]
fwrite(AX, file.path(out, "axis_testable.tsv"), sep = "\t")
REC <- data.table(set = c("discovery_774", "steatosite_515"), median_recon_spearman = c(median(disc_recon), median(ss_recon)),
                  p05 = c(quantile(disc_recon, .05), quantile(ss_recon, .05)), min = c(min(disc_recon), min(ss_recon)))
fwrite(REC, file.path(out, "reconstruction_fit.tsv"), sep = "\t")
fwrite(data.table(sample_alias = colnames(h1), batch = U$batch, t(h1[PR, ]), recon_spearman = ss_recon),
       file.path(out, "steatosite_scores_shift.tsv"), sep = "\t")
fwrite(data.table(sample_alias = colnames(h0), t(h0[PR, ])), file.path(out, "steatosite_scores_zero_shift.tsv"), sep = "\t")
print(AX); print(REC)

# ---- stage (read only after the testable set is written)
st <- fread(cmd = paste("tr -d '\\r' <", shQuote(file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"))),
            colClasses = "character")
U[, stage_num := suppressWarnings(as.numeric(st$nash_crn_kleiner_fibrosis_stage[match(sample_alias, st$sample_alias)]))]
S <- U[!is.na(stage_num), .(sample_alias, batch = factor(batch), inferred_sex = factor(inferred_sex), stage_num)]
S[, stage4 := factor(c("F0F1", "F2", "F3", "F4")[stage_num], levels = c("F0F1", "F2", "F3", "F4"))]
S[, bin := factor(c("F0F1", "F2", "F3F4", "F3F4")[stage_num], levels = c("F0F1", "F2", "F3F4"))]
S[, stratum := paste(batch, inferred_sex)]
stopifnot(nrow(S) == 493L)

# discovery: release scores (z over all 774, as the Fig 3E producer), 387-unit roster, F0+F1 reference
DM <- fread(file.path(ROOT, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_sample_manifest.tsv"))
RZ <- copy(REL); for (p in PR) set(RZ, j = p, value = as.numeric(scale(RZ[[p]])))
DM <- merge(DM, RZ[, c("sample_id", PR), with = FALSE], by = "sample_id")
stopifnot(nrow(DM) == 387L)
DM[, stage := factor(fifelse(fib_group %in% c("F0", "F1"), "F0F1", fib_group), levels = c("F0F1", "F2", "F3", "F4"))]
DM[, bin := factor(fifelse(fib_group %in% c("F0", "F1"), "F0F1", fifelse(fib_group == "F2", "F2", "F3F4")), levels = c("F0F1", "F2", "F3F4"))]
DM[, fib0 := factor(fib_group, levels = paste0("F", 0:4))]
coefs <- function(fit, pre, lev) vapply(lev, function(l) unname(coef(fit)[paste0(pre, l)]), numeric(1))
DISC <- rbindlist(lapply(PR, function(p) {
  f4 <- lm(reformulate(c("dataset", "inferred_sex", "stage"), p), data = DM)
  fb <- lm(reformulate(c("dataset", "inferred_sex", "bin"), p), data = DM)
  f0 <- lm(reformulate(c("dataset", "inferred_sex", "fib0"), p), data = DM)
  data.table(program = p, F2 = coefs(f4, "stage", "F2"), F3 = coefs(f4, "stage", "F3"), F4 = coefs(f4, "stage", "F4"),
             F3F4_bin = coefs(fb, "bin", "F3F4"), check_F4_vs_F0_release = coefs(f0, "fib0", "F4")) }))

ss_fit <- function(H, perm = FALSE) {
  d <- cbind(S, as.data.table(t(H[PR, S$sample_alias])))
  for (p in PR) set(d, j = p, value = as.numeric(scale(d[[p]])))
  rbindlist(lapply(PR, function(p) {
    f4 <- lm(reformulate(c("batch", "inferred_sex", "stage4"), p), data = d)
    fb <- lm(reformulate(c("batch", "inferred_sex", "bin"), p), data = d)
    se <- sqrt(diag(vcov(fb)))[["binF3F4"]]
    data.table(program = p, F2 = coefs(f4, "stage4", "F2"), F3 = coefs(f4, "stage4", "F3"), F4 = coefs(f4, "stage4", "F4"),
               F3F4_bin = coefs(fb, "bin", "F3F4"), F3F4_bin_se = se) })) }
R1 <- ss_fit(h1); R0 <- ss_fit(h0)

# permutation null for the bin coefficient (stage permuted within batch x sex)
set.seed(20261104L)
d0 <- cbind(S, as.data.table(t(h1[PR, S$sample_alias]))); for (p in PR) set(d0, j = p, value = as.numeric(scale(d0[[p]])))
X0 <- model.matrix(~ batch + inferred_sex, d0)
perm_b <- replicate(1000, {
  sp <- d0[, .(i = .I, s = sample(stage_num)), by = stratum]
  sn <- integer(nrow(d0)); sn[sp$i] <- sp$s
  b <- factor(c("F0F1", "F2", "F3F4", "F3F4")[sn], levels = c("F0F1", "F2", "F3F4"))
  X <- cbind(X0, model.matrix(~ b)[, -1, drop = FALSE])
  qr.coef(qr(X), as.matrix(d0[, ..PR]))["bF3F4", ] })
FAM <- merge(R1[, .(program, ss_F3F4 = F3F4_bin, ss_se = F3F4_bin_se)], DISC[, .(program, disc_F3F4 = F3F4_bin)], by = "program")
FAM <- merge(FAM, AX[, .(program, testable)], by = "program")
FAM[, direction := sign(disc_F3F4)]
FAM[, p_perm_one_sided := vapply(seq_len(.N), function(i)
  (1 + sum(perm_b[program[i], ] * direction[i] >= ss_F3F4[i] * direction[i])) / 1001, numeric(1))]
FAM[testable == TRUE, holm := p.adjust(p_perm_one_sided, "holm")]
FAM[, zero_shift_F3F4 := R0$F3F4_bin[match(program, R0$program)]]

prof <- function(a, b) cor(c(0, a), c(0, b), method = "spearman")
CMP <- merge(DISC, R1, by = "program", suffixes = c("_disc", "_ss"))
CMP[, sign_agree_F2_F3_F4 := (sign(F2_disc) == sign(F2_ss)) + (sign(F3_disc) == sign(F3_ss)) + (sign(F4_disc) == sign(F4_ss))]
CMP[, profile_spearman := mapply(function(a2, a3, a4, b2, b3, b4) prof(c(a2, a3, a4), c(b2, b3, b4)), F2_disc, F3_disc, F4_disc, F2_ss, F3_ss, F4_ss)]
CMP[, label := c("Inflammatory-EMT", "Hepatocyte-Metabolic", "Fibrotic-ECM", "Unresolved", "Noncoding (protocol-sensitive)")[match(program, PR)]]

fwrite(DISC, file.path(out, "discovery_refit_F0F1_reference.tsv"), sep = "\t")
fwrite(R1, file.path(out, "steatosite_stage_activity.tsv"), sep = "\t")
fwrite(R0, file.path(out, "steatosite_stage_activity_zero_shift.tsv"), sep = "\t")
fwrite(CMP, file.path(out, "discovery_vs_steatosite.tsv"), sep = "\t")
fwrite(FAM, file.path(out, "family_F3F4_holm.tsv"), sep = "\t")
fwrite(data.table(program = PR, perm_b), file.path(out, "perm_null_F3F4.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200); print(DISC, digits = 3); print(CMP, digits = 3); print(FAM, digits = 3)
