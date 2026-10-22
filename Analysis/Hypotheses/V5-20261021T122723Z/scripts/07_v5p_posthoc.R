# V5P (post hoc; V3 known; not a test). Rescore V3 (GSE256398 snRNA) and V5 (GSE298719 Flex)
# hepatocyte pseudobulks on the 18 Ductular injury genes shared by both, fit the same model in each
# cohort, report a fixed-effect inverse-variance summary and a pooled donor model with a cohort
# term, with leave-S35-out. Confidence intervals only. PRESPEC_V5.md, section V5P.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z")
OUT <- file.path(RUN, "outputs/v5p_posthoc")
if (dir.exists(OUT)) stop("refusing to overwrite")
dir.create(OUT)
source(file.path(ROOT, "scripts/manuscript/program_observability_map/analysis_lib.R"))
HOT <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
V3 <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk")
PRIMARY <- "hotspot_hepatocytes_48f39dd4d817a10e"
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")

membership <- fread(file.path(HOT, "program_membership_v2.tsv"))
registry <- fread(file.path(HOT, "program_registry_v2.tsv"))[program_uid == PRIMARY]
annotation <- fread(ANNOT, select = c("gene_id", "gene_name"))
shared <- strsplit(fread(file.path(RUN, "outputs/job3/realised_gene_set_check.tsv"))$realised_genes, ";")[[1]]
mem <- membership[program_uid == PRIMARY & mapped_symbol %in% shared]
stopifnot(nrow(mem) == 18)

score <- function(counts, donors, cohort, keep = NULL) {
  lc <- cpm(calcNormFactors(DGEList(counts[, donors])), log = TRUE, prior.count = 1)
  if (!is.null(keep)) lc <- lc[rownames(lc) %in% keep, , drop = FALSE]
  S <- collapse_symbols(lc, rownames(lc), annotation)
  sc <- score_programs(S, data.table(sample_id = donors, dataset = cohort), mem, registry, 0.80)
  as.numeric(sc$primary[PRIMARY, donors])
}
# V5: detection-restricted rows, as in the primary
pb5 <- fread(file.path(RUN, "outputs/job3/pseudobulk_hep.tsv.gz"))
det5 <- fread(file.path(RUN, "outputs/job3/detection_by_set.tsv.gz"))
dq <- fread(file.path(RUN, "outputs/donor_qc.tsv"), colClasses = list(character = "id"))[included == TRUE][order(gsm)]
c5 <- as.matrix(pb5[, -1]); rownames(c5) <- pb5$gene_id
d5 <- data.table(donor = dq$gsm, cohort = "GSE298719", disease = dq$disease, sex_f = dq$sex_f,
                 group = dq$group, score = score(c5, dq$gsm, "GSE298719", det5$gene_id[det5$hep >= 0.01]))
# V3: V3 route (all genes on the axis)
pb3 <- fread(file.path(V3, "pseudobulk_hepatocyte.tsv.gz"))
m3 <- fread(file.path(V3, "donor_metadata.tsv"))
c3 <- as.matrix(pb3[, m3$donor, with = FALSE]); rownames(c3) <- pb3$gene_id
d3 <- data.table(donor = m3$donor, cohort = "GSE256398", disease = m3$disease, sex_f = as.integer(m3$sex == "female"),
                 group = m3$stage, score = score(c3, m3$donor, "GSE256398"))
d <- rbind(d5, d3)
wr(d, "donor_scores_shared18.tsv")

hc3 <- function(y, X, j = 2) {
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y); e <- as.vector(y - X %*% b)
  h <- rowSums((X %*% XtXi) * X); C <- XtXi %*% t(X)
  se <- sqrt(sum(C[j, ]^2 * e^2 / (1 - h)^2)); df <- nrow(X) - ncol(X)
  data.table(beta = b[j], hc3_se = se, df = df, ci_low = b[j] - qt(0.975, df) * se, ci_high = b[j] + qt(0.975, df) * se)
}
fit_cohort <- function(x, label) cbind(data.table(analysis = label, n = nrow(x), n_disease = sum(x$disease)),
                                       hc3(x$score, cbind(1, x$disease, x$sex_f)))
fe <- function(rows, label) {
  w <- 1 / rows$hc3_se^2; b <- sum(w * rows$beta) / sum(w); se <- sqrt(1 / sum(w))
  data.table(analysis = label, n = sum(rows$n), n_disease = sum(rows$n_disease), beta = b, hc3_se = se, df = NA,
             ci_low = b - 1.96 * se, ci_high = b + 1.96 * se)
}
pooled <- function(x, label) cbind(data.table(analysis = label, n = nrow(x), n_disease = sum(x$disease)),
                                   hc3(x$score, cbind(1, x$disease, x$sex_f, as.integer(x$cohort == "GSE256398"))))
v3_f2 <- d3[group %in% c("Healthy", "MASH_F2_F3", "MASH_cirrhosis")]
r <- list()
r$v5 <- fit_cohort(d5, "GSE298719, score ~ disease + sex")
r$v3 <- fit_cohort(d3, "GSE256398 (V3 C2, all disease stages), score ~ disease + sex")
r$v3_noS35 <- fit_cohort(d3[donor != "S35"], "GSE256398 without S35")
r$v3_f2 <- fit_cohort(v3_f2, "GSE256398 F2-F4 only (matches V5 disease range)")
r$fe <- fe(rbind(r$v5, r$v3), "fixed-effect IV summary (V5 + V3 C2)")
r$fe_noS35 <- fe(rbind(r$v5, r$v3_noS35), "fixed-effect IV summary, V3 without S35")
r$fe_f2 <- fe(rbind(r$v5, r$v3_f2), "fixed-effect IV summary, F2-F4 in both")
r$pool <- pooled(d, "pooled donors, score ~ disease + sex + cohort")
r$pool_noS35 <- pooled(d[donor != "S35"], "pooled donors without S35")
res <- rbindlist(r, fill = TRUE)[, label := "post hoc; V3 known; not a test"]
wr(res, "v5p_summary.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(res)
