# V6 job 2: Stromal ECM within-ROI contrast (Vinculin+ minus Albumin+), exact sign-flip over
# 24 ROIs, matched random-set null (NC3), slide and model sensitivities, descriptive arms.
# Follows PRESPEC_V6.md. Unit: ROI (donor replication not established).
# Usage: Rscript 02_outcome.R
suppressPackageStartupMessages({library(data.table); library(lme4)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/V6-20261021T234432Z")
OUT <- file.path(RUN, "outputs")
J2 <- file.path(OUT, "job2")
MEMBERSHIP <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv")
PROGRAMS <- c(stromal_ecm = "hotspot_hepatocytes_f05c535ae5bbc0b9", ductular = "hotspot_hepatocytes_48f39dd4d817a10e")
SEED_SETS <- 6042; N_SETS <- 1000

chk <- readLines(file.path(RUN, "CHECKED_BEFORE_OUTCOME.txt"))
if (!any(grepl("PROCEED: yes", chk, fixed = TRUE))) stop("STOP: CHECKED_BEFORE_OUTCOME.txt not PROCEED: yes")
sha <- function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1]
if (sha(file.path(RUN, "PRESPEC_V6.md")) != strsplit(readLines(file.path(RUN, "PRESPEC_V6.md.sha256")), " ")[[1]][1])
  stop("STOP: prespec hash differs")
if (dir.exists(J2)) stop("outputs/job2 exists; refusing to overwrite")
dir.create(J2)
wr <- function(x, f) fwrite(as.data.table(x), file.path(J2, f), sep = "\t", na = "NA")

d <- readRDS(file.path(RUN, "derived/geomx_rebuilt.rds"))
L <- d$L; aoi <- as.data.table(d$aoi)
det <- fread(file.path(OUT, "K4_gene_detection_alb_vin.tsv"))
cov <- fread(file.path(OUT, "K4_coverage.tsv"))
mem <- fread(MEMBERSHIP)[mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
rois <- sort(unique(aoi$roi))

# score per AOI: genes z-scored across the AOIs of the contrast, weighted sum, standardized
score_aois <- function(genes, w, aois) {
  z <- t(scale(t(L[genes, aois, drop = FALSE])))
  z[!is.finite(z)] <- 0
  s <- as.numeric(crossprod(w / sum(w), z))
  setNames(as.numeric(scale(s)), aois)
}
roi_diff <- function(s, target, ref, roi_set = rois) {
  a <- aoi[roi %in% roi_set]
  vapply(roi_set, function(r) s[a[roi == r & stain == target, dcc]] - s[a[roi == r & stain == ref, dcc]], 0)
}
# exact two-sided sign-flip p for the mean of x (meet in the middle over two halves)
signflip_p <- function(x) {
  n <- length(x); h <- n %/% 2
  sums <- function(v) {
    m <- length(v); s <- numeric(2^m)
    for (i in 0:(2^m - 1)) s[i + 1] <- sum(ifelse(bitwAnd(i, bitwShiftL(1L, 0:(m - 1))) > 0, -v, v))
    s
  }
  a <- sums(x[seq_len(h)]); b <- sort(sums(x[(h + 1):n]))
  obs <- abs(sum(x)) - 1e-9
  hits <- sum(vapply(a, function(ai) sum(b >= obs - ai) + sum(b <= -obs - ai), 0))
  hits / 2^n
}
weights_for <- function(uid, genes_ok) {
  w <- mem[program_uid == uid & mapped_symbol %in% genes_ok, .(gene = mapped_symbol, weight = original_l1_weight)]
  w[order(gene)]
}

# matched random sets: non-members of the tested program, passing K4, 10 x 10 deciles of mean
# log2 Q3 and detection fraction over the 48 Albumin+/Vinculin+ AOIs
d_ok <- det[detected == TRUE]
d_ok[, cell := paste(cut(mean_log2q3, unique(quantile(mean_log2q3, 0:10 / 10)), include.lowest = TRUE, labels = FALSE),
                     cut(frac_above_loq, unique(quantile(frac_above_loq, 0:10 / 10)), include.lowest = TRUE, labels = FALSE))]
matched_null <- function(w, uid, target, ref, seed) {
  aois <- aoi[stain %in% c(target, ref), dcc]
  members <- mem[program_uid == uid, mapped_symbol]
  pool <- d_ok[!gene %in% members]
  cells <- d_ok[match(w$gene, gene), cell]
  set.seed(seed)
  vapply(seq_len(N_SETS), function(k) {
    g <- vapply(cells, function(cl) { cand <- pool[cell == cl, gene]; if (length(cand)) sample(cand, 1) else NA_character_ }, "")
    keep <- !is.na(g) & !duplicated(g)
    mean(roi_diff(score_aois(g[keep], sample(w$weight)[keep], aois), target, ref))
  }, 0)
}
contrast <- function(p, target, ref, seed) {
  uid <- PROGRAMS[[p]]
  w <- weights_for(uid, d_ok$gene)
  aois <- aoi[stain %in% c(target, ref), dcc]
  s <- score_aois(w$gene, w$weight, aois)
  D <- roi_diff(s, target, ref)
  null <- matched_null(w, uid, target, ref, seed)
  list(w = w, s = s, D = D, null = null,
       row = data.table(program = p, contrast = paste(target, "minus", ref), n_roi = length(D),
                        n_genes = nrow(w), weight = sum(w$weight), testable_k4 = cov[program == p, testable],
                        mean_D = mean(D), sd_D = sd(D), n_roi_positive = sum(D > 0),
                        signflip_p = signflip_p(D), null_median = median(null),
                        null_q95 = quantile(null, 0.95), percentile = mean(null <= mean(D)),
                        n_sets = sum(is.finite(null))))
}

# ---- primary --------------------------------------------------------------------------
prim <- contrast("stromal_ecm", "Vinculin", "Albumin", SEED_SETS)
r <- prim$row
r[, verdict := if (!testable_k4) "untestable: program not measured" else
  if (mean_D > 0 && signflip_p < 0.05 && mean_D >= null_q95) "pass" else "not supported"]
wr(r, "results_primary.tsv")
wr(data.table(roi = rois, D = prim$D, slide = substr(rois, 1, 1)), "primary_roi_differences.tsv")
wr(data.table(set = seq_along(prim$null), mean_D = prim$null), "NC3_matched_null_primary.tsv")
wr(data.table(dcc = names(prim$s), score = prim$s)[aoi, on = "dcc", nomatch = 0], "primary_aoi_scores.tsv")

# ---- sensitivities ---------------------------------------------------------------------
sens <- list()
for (sl in c("A", "B", "C")) {
  Ds <- prim$D[substr(rois, 1, 1) == sl]
  sens[[length(sens) + 1]] <- data.table(analysis = paste("within slide", sl), n_roi = length(Ds), mean_D = mean(Ds), signflip_p = NA_real_)
  Dl <- prim$D[substr(rois, 1, 1) != sl]
  sens[[length(sens) + 1]] <- data.table(analysis = paste("leave out slide", sl), n_roi = length(Dl), mean_D = mean(Dl), signflip_p = signflip_p(Dl))
}
Dn <- prim$D[rois != "A9"]
sens[[length(sens) + 1]] <- data.table(analysis = "without ROI A9", n_roi = length(Dn), mean_D = mean(Dn), signflip_p = signflip_p(Dn))
md <- merge(data.table(dcc = names(prim$s), score = prim$s), aoi, by = "dcc")
md[, compartment := factor(stain, levels = c("Albumin", "Vinculin"))][, log10_aligned := log10(aligned)]
for (f in list(score ~ compartment + (1 | slide / roi), score ~ compartment + log10_aligned + (1 | slide / roi))) {
  fit <- lmer(f, data = md, REML = TRUE)
  cf <- summary(fit)$coefficients
  vc <- as.data.table(VarCorr(fit))
  sens[[length(sens) + 1]] <- data.table(analysis = paste("lmer", deparse(f)), n_roi = length(rois),
    mean_D = cf["compartmentVinculin", "Estimate"], signflip_p = NA_real_,
    se = cf["compartmentVinculin", "Std. Error"], t = cf["compartmentVinculin", "t value"],
    var_slide = vc[grp == "slide", vcov], var_roi = vc[grp == "roi:slide", vcov], singular = isSingular(fit))
}
wr(rbindlist(sens, fill = TRUE), "sensitivities_primary.tsv")
aois_av <- aoi[stain %in% c("Vinculin", "Albumin"), dcc]
zg <- t(scale(t(L[prim$w$gene, aois_av, drop = FALSE])))
wr(data.table(gene = prim$w$gene, weight = prim$w$weight,
              mean_D = vapply(prim$w$gene, function(g) mean(roi_diff(setNames(zg[g, ], aois_av), "Vinculin", "Albumin")), 0)),
   "primary_per_gene_D.tsv")

# ---- descriptive arms ------------------------------------------------------------------
e1 <- contrast("ductular", "Vinculin", "Albumin", SEED_SETS + 1)
e2a <- contrast("stromal_ecm", "CD45", "Albumin", SEED_SETS + 2)
e2b <- contrast("ductular", "CD45", "Albumin", SEED_SETS + 3)
desc <- rbindlist(list(e1$row[, arm := "E1"], e2a$row[, arm := "E2"], e2b$row[, arm := "E2"]))
wr(desc, "E_descriptive_arms.tsv")
e3 <- rbindlist(lapply(names(PROGRAMS), function(p) {
  w <- weights_for(PROGRAMS[[p]], d_ok$gene)
  s <- score_aois(w$gene, w$weight, aoi$dcc)
  merge(data.table(dcc = names(s), score = s), aoi, by = "dcc")[, .(program = p, mean_score = mean(score), n_aoi = .N), by = .(slide, disease, stain)]
}))
wr(e3, "E3_scores_by_slide_stain.tsv")
writeLines(capture.output(sessionInfo()), file.path(J2, "sessionInfo.txt"))
print(r); print(rbindlist(sens, fill = TRUE)); print(desc)
