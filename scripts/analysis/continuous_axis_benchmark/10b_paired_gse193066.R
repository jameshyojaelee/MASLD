#!/usr/bin/env Rscript
# 10b: the within-person test. 54 donors, two biopsies each.
#
# Cross-sectional data cannot distinguish "people at different stages differ along
# one axis" from "people move along that axis". Paired biopsies can. The powered
# question is P2: is the WITHIN-PERSON delta vector parallel to the cross-sectional
# F0->F4 direction?
#
# If yes, the amplification structure is a property of change within people. If no,
# it is a property of differences between people and must be described that way.
#
# This script refuses to run unless the sealed predictions still hash to the value
# 10a recorded, so the predictions cannot be edited after seeing these data.

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
})
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")

# ------------------------------------------------------------- seal check ----
pred_path <- file.path(OUT, "paired", "predictions.json")
sha_path <- file.path(OUT, "paired", "predictions.sha256")
assert_true(file.exists(pred_path) && file.exists(sha_path),
            "Sealed predictions are absent; run 10a before any paired data is loaded")
recorded <- readLines(sha_path)[1]
actual <- sha256_file(pred_path)
assert_true(identical(recorded, actual),
            paste0("PREDICTION SEAL BROKEN. recorded=", recorded, " actual=", actual,
                   ". The paired analysis is not interpretable and will not run."))
log_step("prediction seal verified: ", substr(actual, 1, 16), "...")

pre <- read_prespec()
set.seed(pre$seeds$master)

# ------------------------------------------------------------- build data ----
cw <- fread(file.path(project_root(),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates",
  "program-context-v2-candidate-2026-08-07/fibrosis-adjacent-true-kleiner-lvqw-v1",
  "audits/gse193066_donor_biopsy_crosswalk.tsv"))
cw <- cw[pass_technical == TRUE]
paired_donors <- cw[, .N, by = donor_id][N == 2L, donor_id]
cw <- cw[donor_id %in% paired_donors]
log_step("paired donors with both biopsies passing QC: ", length(paired_donors))

counts_path <- file.path(project_root(), "results/remediation/bg001",
  "bg001-fragment-v211-gencode49-20260807T194845Z/frozen_sets/read_counts/GSE193066/gene_counts.txt")
fc <- fread(counts_path, skip = 1)
gid <- fc[[1]]
mat <- as.matrix(fc[, 7:ncol(fc)])
# featureCounts names its columns with the full BAM path, e.g.
#   .../alignments/star/SRR17442769/SRR17442769.Aligned.sortedByCoord.out.bam
# so the run accession has to be pulled out explicitly. Stripping only the .bam
# suffix leaves "SRR17442769.Aligned.sortedByCoord.out", which matches nothing.
raw_cols <- colnames(mat)
acc <- sub("^.*?(SRR[0-9]+).*$", "\\1", basename(raw_cols))
assert_true(all(grepl("^SRR[0-9]+$", acc)),
            paste0("Could not parse run accessions from featureCounts columns; first: ",
                   raw_cols[1]))
assert_true(!anyDuplicated(acc), "Duplicate run accessions among count columns")
colnames(mat) <- acc
rownames(mat) <- gid
log_step("count matrix: ", nrow(mat), " genes x ", ncol(mat), " runs")

ref_genes <- rownames(readRDS(substrate_path("dge")))
shared <- intersect(ref_genes, rownames(mat))
log_step("genes shared with the 23,370 reference universe: ", length(shared))
assert_true(length(shared) > 15000L, "Too few shared genes; the annotation build does not match")

ids <- intersect(cw$sample_id, colnames(mat))
assert_true(length(ids) >= 80L, paste0(
  "Only ", length(ids), " of ", nrow(cw), " crosswalk samples are present in the count matrix. ",
  "Crosswalk examples: ", paste(head(cw$sample_id, 3), collapse = ", "),
  " | count-matrix examples: ", paste(head(colnames(mat), 3), collapse = ", ")))
cw <- cw[sample_id %in% ids]
paired_donors <- cw[, .N, by = donor_id][N == 2L, donor_id]
cw <- cw[donor_id %in% paired_donors][order(donor_id, biopsy)]
log_step("final paired donors: ", length(paired_donors), " (", nrow(cw), " samples)")

dge <- DGEList(mat[shared, cw$sample_id])
dge <- calcNormFactors(dge, method = "RLE")
cw[, biopsy_idx := factor(fifelse(grepl("1st", biopsy), "first", "second"),
                          levels = c("first", "second"))]
cw[, donor := factor(donor_id)]
X <- model.matrix(~ donor + biopsy_idx, data = cw)
v <- voomWithQualityWeights(dge, X, plot = FALSE)
E <- v$E

# Per-donor delta on the log scale: second minus first biopsy.
first <- cw[biopsy_idx == "first"]; second <- cw[biopsy_idx == "second"]
setkey(first, donor_id); setkey(second, donor_id)
D <- E[, second[.(paired_donors), sample_id]] - E[, first[.(paired_donors), sample_id]]
colnames(D) <- paired_donors
log_step("delta matrix: ", nrow(D), " genes x ", ncol(D), " donors")

# ------------------------------------------------------------------- P1 -----
dstage <- second[.(paired_donors), harmonized_fibrosis_stage] -
          first[.(paired_donors), harmonized_fibrosis_stage]
# Within-person axis change, projected on the cross-sectional direction (below).
em <- fread(file.path(OUT, "q2", "A4_histology_M_gene.tsv.gz"))
f4 <- em[window == "F4"][match(rownames(D), feature)]
ok <- is.finite(f4$estimate)
w <- f4$estimate[ok]
proj <- as.numeric(t(D[ok, , drop = FALSE]) %*% w) / sqrt(sum(w^2))
p1 <- suppressWarnings(cor.test(proj, dstage, method = "spearman",
                                alternative = "greater", exact = FALSE))

# ------------------------------------------------------------------- P2 -----
# Cosine between each donor's delta vector and the cross-sectional F0->F4 direction.
cosines <- apply(D[ok, , drop = FALSE], 2, function(d)
  sum(d * w) / (sqrt(sum(d^2)) * sqrt(sum(w^2))))
mean_delta <- rowMeans(D[ok, , drop = FALSE])
cos_mean <- sum(mean_delta * w) / (sqrt(sum(mean_delta^2)) * sqrt(sum(w^2)))
sv <- svd(scale(D[ok, , drop = FALSE], center = TRUE, scale = FALSE))
pve1_delta <- sv$d[1]^2 / sum(sv$d^2)

boot_cos <- replicate(2000, {
  i <- sample(ncol(D), replace = TRUE)
  md <- rowMeans(D[ok, i, drop = FALSE])
  sum(md * w) / (sqrt(sum(md^2)) * sqrt(sum(w^2)))
})

thr <- pre$paired_gse193066$predictions$P2$threshold_disattenuated_cosine
res <- data.table(
  n_donors = ncol(D), n_genes = sum(ok),
  P1_spearman_rho = unname(p1$estimate), P1_p_one_sided = p1$p.value,
  P1_predicted_positive = TRUE, P1_supported = unname(p1$estimate) > 0 && p1$p.value < 0.05,
  P2_cosine_mean_delta = cos_mean,
  P2_cosine_boot_lo = quantile(boot_cos, 0.025), P2_cosine_boot_hi = quantile(boot_cos, 0.975),
  P2_median_donor_cosine = median(cosines),
  P2_frac_donors_positive_cosine = mean(cosines > 0),
  P2_threshold = thr, P2_supported = cos_mean > thr,
  delta_pve1 = pve1_delta,
  n_stage_increased = sum(dstage > 0, na.rm = TRUE),
  n_stage_decreased = sum(dstage < 0, na.rm = TRUE),
  n_stage_unchanged = sum(dstage == 0, na.rm = TRUE))

write_tsv_once(res, file.path(OUT, "paired", "10b_paired_verdict.tsv"))
write_tsv_once(data.table(donor = colnames(D), cosine = cosines,
                          projection = proj, delta_stage = dstage),
               file.path(OUT, "paired", "10b_per_donor.tsv"))
for (n in names(res)) cat(sprintf("%-34s %s\n", n, format(res[[n]][1])))
log_step("PAIRED_COMPLETE P1=", res$P1_supported, " P2=", res$P2_supported)
