#!/usr/bin/env Rscript
# U10 re-derivation of the within-stage permutation tests for the 43 molecular systems
# (435 stage-complete units; Figure S4U strict maxT calls and the joint two-cohort maxT).
#
# Test as described in the legends and run parameters: residualize the standardized axis
# and each system score on fibrosis stage (factor) and inferred sex within cohort; the
# statistic is the residual Pearson correlation; permute the axis residual within
# stage x sex strata; two-sided; one maximum |r| over the 43 systems per draw for maxT;
# p = (1 + exceedances) / (1 + B). Strict call: maxT p < 0.05 in both cohorts on both
# co-primary axes and the same sign in all four.
# Joint test: Fisher-z weighted sum over the two cohorts (weight n - (k - 1) - 3, k = nuisance
# columns), one permutation stream per axis, both cohorts permuted in each draw.
#
# Seeds: producer contract seed 20260817 with offsets axis*100000 + cohort*10000 (per cohort)
# and axis*700000 (joint), then alternate base seeds 20261101 and 20261202 with the same offsets.
# Row order within cohort is sample_id ascending (C collation), strata split in factor-level
# order; both are needed only to consume the random stream in the same order.
#
# Usage: u10_maxT.R <out_dir>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
hac <- file.path(root, "RNA-seq/results/histology_anchored_continuum")
core <- file.path(hac, "candidates/hac-continuum-corrected774-20261009T171017Z")
sys_dir <- file.path(hac, "molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables")
stress_dir <- file.path(hac, "molecular_systems_stress/hac-system-stress-corrected774-20261009T171017Z/tables")
B <- 10000L
cohorts <- c("GSE162694", "GSE213621")
axes <- c("signature_pc1", "fixed_projection")
base_seeds <- c(20260817L, 20261101L, 20261202L)

read_tsv <- function(path) read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
inputs <- c(
  projection = file.path(core, "projection/participant_scores.tsv"),
  unsupervised = file.path(core, "unsupervised/participant_scores.tsv"),
  system_scores = file.path(sys_dir, "participant_system_scores.tsv.gz"),
  stored_permutations = file.path(sys_dir, "system_within_stage_permutations.tsv"),
  stored_joint = file.path(stress_dir, "system_joint_permutations.tsv")
)
proj <- read_tsv(inputs[["projection"]])
unsup <- read_tsv(inputs[["unsupervised"]])
cov <- proj[proj$axis_id == "fixed_projection", c("sample_id", "dataset", "inferred_sex", "fibrosis_stage")]
axis_tab <- rbind(
  data.frame(proj[proj$axis_id == "fixed_projection", c("sample_id", "dataset")], axis_id = "fixed_projection",
             axis_raw = proj$axis_raw[proj$axis_id == "fixed_projection"]),
  data.frame(unsup[unsup$axis_id == "signature_pc1", c("sample_id", "dataset")], axis_id = "signature_pc1",
             axis_raw = unsup$axis_raw[unsup$axis_id == "signature_pc1"])
)
ss <- read_tsv(gzfile(inputs[["system_scores"]]))
ss <- ss[ss$aggregation_method == "collection_balanced" & ss$dataset %in% cohorts, ]
systems <- sort(unique(ss$community_id))
stopifnot(length(systems) == 43L)

zstd <- function(x) (x - mean(x)) / sd(x)

# One cohort x axis frame: rows sample_id ascending, residualized x and Y.
build <- function(cohort, axis) {
  cv <- cov[cov$dataset == cohort & !is.na(cov$fibrosis_stage) & !is.na(cov$inferred_sex), ]
  cv <- cv[order(cv$sample_id, method = "radix"), ]
  ax <- axis_tab[axis_tab$dataset == cohort & axis_tab$axis_id == axis, ]
  x_raw <- ax$axis_raw[match(cv$sample_id, ax$sample_id)]
  sub <- ss[ss$dataset == cohort, ]
  Y <- sapply(systems, function(s) {
    t <- sub[sub$community_id == s, ]
    t$system_score_z[match(cv$sample_id, t$sample_id)]
  })
  stopifnot(all(is.finite(x_raw)), all(is.finite(Y)))
  nuis <- model.matrix(~ factor(fibrosis_stage) + factor(inferred_sex), data = cv)
  qrn <- qr(nuis)
  xr <- as.numeric(qr.resid(qrn, zstd(x_raw))); xr <- xr - mean(xr)
  Yr <- qr.resid(qrn, Y); Yr <- sweep(Yr, 2L, colMeans(Yr), "-")
  strata <- interaction(factor(cv$fibrosis_stage), factor(cv$inferred_sex), drop = TRUE)
  list(n = nrow(cv), xr = xr, Yr = Yr, syy = colSums(Yr^2), groups = split(seq_len(nrow(cv)), strata),
       k = ncol(nuis), strata_sizes = table(strata))
}

permute_x <- function(x, groups) {
  for (g in groups) x[g] <- x[g][sample.int(length(g), length(g), replace = FALSE)]
  x
}
corr_vec <- function(x, f) as.numeric(crossprod(x, f$Yr)) / sqrt(sum(x^2) * f$syy)

frames <- list()
for (a in axes) for (cc in cohorts) frames[[paste(a, cc, sep = "|")]] <- build(cc, a)
for (cc in cohorts) {
  cat(cc, "n =", frames[[paste("fixed_projection", cc, sep = "|")]]$n, "strata sizes:\n")
  print(frames[[paste("fixed_projection", cc, sep = "|")]]$strata_sizes)
}

per_cohort <- list()
joint <- list()
for (base in base_seeds) {
  for (ai in seq_along(axes)) {
    for (ci in seq_along(cohorts)) {
      f <- frames[[paste(axes[ai], cohorts[ci], sep = "|")]]
      obs <- corr_vec(f$xr, f)
      seed <- base + ai * 100000L + ci * 10000L
      set.seed(seed)
      exceed <- integer(43); max_exceed <- integer(43)
      for (b in seq_len(B)) {
        nul <- corr_vec(permute_x(f$xr, f$groups), f)
        exceed <- exceed + (abs(nul) >= abs(obs))
        max_exceed <- max_exceed + (max(abs(nul)) >= abs(obs))
      }
      per_cohort[[length(per_cohort) + 1L]] <- data.frame(
        base_seed = base, seed = seed, axis_id = axes[ai], dataset = cohorts[ci], community_id = systems,
        n = f$n, residual_pearson = obs, perm_p = (1 + exceed) / (1 + B), maxT_p = (1 + max_exceed) / (1 + B))
    }
    # joint two-cohort Fisher-z statistic
    fl <- lapply(cohorts, function(cc) frames[[paste(axes[ai], cc, sep = "|")]])
    w <- sapply(fl, function(f) f$n - (f$k - 1L) - 3L)
    comb <- function(xs) {
      acc <- 0
      for (i in seq_along(fl)) {
        r <- pmin(pmax(corr_vec(xs[[i]], fl[[i]]), -0.999999999), 0.999999999)
        acc <- acc + w[i] * atanh(r)
      }
      acc / sqrt(sum(w))
    }
    obs_j <- comb(lapply(fl, function(f) f$xr))
    seed_j <- base + ai * 700000L
    set.seed(seed_j)
    ex <- integer(43); mx <- integer(43)
    for (b in seq_len(B)) {
      nul <- comb(lapply(fl, function(f) permute_x(f$xr, f$groups)))
      ex <- ex + (abs(nul) >= abs(obs_j))
      mx <- mx + (max(abs(nul)) >= abs(obs_j))
    }
    joint[[length(joint) + 1L]] <- data.frame(
      base_seed = base, seed = seed_j, axis_id = axes[ai], community_id = systems, fisher_weights = paste(w, collapse = ";"),
      joint_z = obs_j, joint_p = (1 + ex) / (1 + B), joint_maxT_p = (1 + mx) / (1 + B))
  }
}
per_cohort <- do.call(rbind, per_cohort)
joint <- do.call(rbind, joint)
per_cohort$perm_q <- ave(per_cohort$perm_p, per_cohort$base_seed, per_cohort$axis_id, per_cohort$dataset,
                         FUN = function(p) p.adjust(p, "BH", n = 43))
joint$joint_q <- ave(joint$joint_p, joint$base_seed, joint$axis_id, FUN = function(p) p.adjust(p, "BH", n = 43))
write.table(per_cohort, file.path(out_dir, "maxT_per_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(joint, file.path(out_dir, "maxT_joint.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# calls per base seed
calls <- do.call(rbind, lapply(base_seeds, function(base) {
  pc <- per_cohort[per_cohort$base_seed == base, ]
  jt <- joint[joint$base_seed == base, ]
  do.call(rbind, lapply(systems, function(s) {
    r <- pc[pc$community_id == s, ]
    j <- jt[jt$community_id == s, ]
    same4 <- length(unique(sign(r$residual_pearson))) == 1L
    data.frame(base_seed = base, community_id = s, same_sign_four = same4,
               max_maxT_p = max(r$maxT_p), strict_maxT = same4 && all(r$maxT_p < 0.05),
               perm_bh_all_four = same4 && all(r$perm_q < 0.05),
               max_joint_maxT_p = max(j$joint_maxT_p),
               joint_maxT = length(unique(sign(j$joint_z))) == 1L && all(j$joint_maxT_p < 0.05),
               joint_bh = length(unique(sign(j$joint_z))) == 1L && all(j$joint_q < 0.05))
  }))
}))
write.table(calls, file.path(out_dir, "maxT_calls_by_seed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nCall counts by base seed\n")
print(aggregate(cbind(strict_maxT, perm_bh_all_four, joint_maxT, joint_bh) ~ base_seed, data = calls, FUN = sum))
wide <- reshape(calls[, c("base_seed", "community_id", "strict_maxT")], idvar = "community_id",
                timevar = "base_seed", direction = "wide")
cat("\nSystems whose strict maxT call changes across seeds\n")
print(wide[apply(wide[, -1], 1, function(v) length(unique(v)) > 1), ])
wj <- reshape(calls[, c("base_seed", "community_id", "joint_maxT")], idvar = "community_id",
              timevar = "base_seed", direction = "wide")
cat("\nSystems whose joint maxT call changes across seeds\n")
print(wj[apply(wj[, -1], 1, function(v) length(unique(v)) > 1), ])
cat("\nStrict-call systems with max maxT p in [0.03, 0.07] at the producer seed\n")
print(calls[calls$base_seed == 20260817L & calls$max_maxT_p >= 0.03 & calls$max_maxT_p <= 0.07, ])

# same-seed comparison with stored producer tables
st <- read_tsv(inputs[["stored_permutations"]])
m <- merge(per_cohort[per_cohort$base_seed == 20260817L, ], st,
           by = c("community_id", "axis_id", "dataset"), suffixes = c("", ".stored"))
cat("\nSame-seed per-cohort comparison: rows", nrow(m),
    "| max|d r|", max(abs(m$residual_pearson - m$residual_pearson.stored)),
    "| p identical", sum(abs(m$perm_p - m$permutation_p_value) < 1e-12),
    "| maxT identical", sum(abs(m$maxT_p - m$maxT_fwer_p_value) < 1e-12), "\n")
if (file.exists(inputs[["stored_joint"]])) {
  sj <- read_tsv(inputs[["stored_joint"]])
  mj <- merge(joint[joint$base_seed == 20260817L, ], sj, by = c("community_id", "axis_id"))
  cat("Same-seed joint comparison: rows", nrow(mj), "| max|d z|", max(abs(mj$joint_z - mj$joint_fisher_z)),
      "| p identical", sum(abs(mj$joint_p - mj$joint_permutation_p_value) < 1e-12),
      "| maxT identical", sum(abs(mj$joint_maxT_p - mj$joint_maxT_fwer_p_value) < 1e-12), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_R.txt"))
sha <- sapply(inputs, function(p) strsplit(system2("sha256sum", p, stdout = TRUE), " ")[[1]][1])
write.table(data.frame(role = names(inputs), path = inputs, sha256 = sha),
            file.path(out_dir, "inputs_sha256_R.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("B", B), paste("base_seeds", paste(base_seeds, collapse = ",")),
             "per-cohort seed = base + axis_index*100000 + cohort_index*10000 (axis order signature_pc1, fixed_projection; cohort order GSE162694, GSE213621)",
             "joint seed = base + axis_index*700000", "family = 43 systems per axis x cohort; BH n = 43"),
           file.path(out_dir, "run_parameters_R.txt"))
