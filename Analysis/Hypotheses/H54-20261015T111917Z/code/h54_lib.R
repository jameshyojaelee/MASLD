# H54 shared helpers. Sourced by stage1.R and stage2.R (rnaseq env).
# The T4 library and the program-map helpers are sourced read-only from the repo.
# Nothing outside the H54 run directory is written.

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(Matrix)
})
options(digits = 15, scipen = 999)

H54_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H54_RUN <- file.path(H54_ROOT, "Analysis/Hypotheses/H54-20261015T111917Z")
H54_PRESPEC_SHA <- "047480f979ad5ed29d6c22e6ae70063001a10635a5788e8a492e686aa94f661f"
H54_SEED <- 20260923L
SCRIPT_DIR <- file.path(H54_ROOT, "scripts/manuscript/program_observability_map")

# MASLD_PROJECT_ROOT pins config.R to the same tree the hashes were taken on.
Sys.setenv(MASLD_PROJECT_ROOT = H54_ROOT)
source(file.path(SCRIPT_DIR, "config.R"))
source(file.path(SCRIPT_DIR, "analysis_lib.R"))
source(file.path(SCRIPT_DIR, "t4_interaction_lib.R"))

T4_DIR <- file.path(H54_ROOT, "RNA-seq/results/manuscript_release/candidates",
                    "resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9",
                    "systems/T4_interaction")

# Prespecified constants (PRESPEC_H54.md).
PC1_S1 <- 0.897361912613355
PC1_BASE <- 0.709472008490316
N_PERM <- 1000L          # stored within-lineage label permutations per cohort
N_PERM_E <- 200L         # permutations 1..200 define the null mean inside E
N_BOOT <- 1000L          # stored participant bootstrap draws per set
N_NC1 <- 200L            # gene-permutation draws per cohort
N_G1 <- 400L             # holdout power draws per scenario
MDE_MULT <- 2.80
EXPECTED_N <- c(A91 = 91L, B87 = 87L, A66 = 66L, B53 = 53L, A94 = 94L, A37 = 37L, H105 = 105L)
SET_ROLE <- c(A91 = "primary", B87 = "primary", A66 = "secondary", B53 = "secondary",
              A94 = "descriptive", A37 = "descriptive", H105 = "positive_control")
SET_COHORT <- c(A91 = "A", B87 = "B", A66 = "A", B53 = "B", A94 = "A", A37 = "A", H105 = "H")
EXPECTED_TESTABLE <- c(A = 113L, B = 110L, H = 113L)

msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), sprintf(...))

sha256_of <- function(p) sub("[[:space:]].*$", "", system2("sha256sum", shQuote(p), stdout = TRUE)[[1]])

# Prespec hash and all 28 input hashes. Every stage calls this first.
verify_prespec_and_inputs <- function() {
  got <- sha256_of(file.path(H54_RUN, "PRESPEC_H54.md"))
  if (!identical(got, H54_PRESPEC_SHA)) fail("PRESPEC_H54.md sha256 changed: ", got)
  rec <- sub("[[:space:]].*$", "", readLines(file.path(H54_RUN, "PRESPEC_H54.md.sha256"))[[1]])
  if (!identical(rec, H54_PRESPEC_SHA)) fail("recorded prespec sha256 differs")
  man <- fread(file.path(H54_RUN, "code/inputs_sha256.tsv"), colClasses = "character")
  if (nrow(man) != 28L) fail("input manifest must list 28 inputs")
  prespec_text <- paste(readLines(file.path(H54_RUN, "PRESPEC_H54.md")), collapse = "\n")
  out <- man[, {
    p <- file.path(H54_ROOT, path)
    if (!file.exists(p)) fail("missing input ", path)
    if (!grepl(sha256, prespec_text, fixed = TRUE)) fail("manifest hash not in prespec: ", path)
    s <- sha256_of(p)
    sz <- file.size(p)
    if (!identical(s, sha256)) fail("sha256 mismatch ", path, ": ", s)
    if (!identical(as.character(sz), bytes)) fail("byte size mismatch ", path)
    .(ok = TRUE)
  }, by = .(path, bytes, sha256)]
  msg("prespec sha256 %s and %d of 28 input hashes verified", substr(got, 1, 16), sum(out$ok))
  invisible(out)
}

refuse_if_exists <- function(dir) {
  if (file.exists(dir)) fail("Refusing to run: output directory exists: ", dir)
  dir.create(dir, recursive = FALSE)
  invisible(dir)
}

# ------------------------------------------------------------- frozen objects
load_frozen <- function() {
  communities <- fread(file.path(T4_DIR, "program_communities.tsv"))
  coord <- fread(file.path(T4_DIR, "coordinated_edges.tsv"))
  allz <- fread(file.path(T4_DIR, "program_interaction_edges.tsv.gz"),
                select = c("program_a", "program_b", "z_partial_standardized", "r_partial", "support"))
  membership <- fread(program_membership)
  registry <- fread(program_registry)
  features <- communities$feature_id
  weights <- t4_membership_weights(membership, features)
  pairs <- t4_pair_table(features, weights)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
  pairs[, ekey := key(program_a, program_b)]
  coord[, ekey := key(program_a, program_b)]
  allz[, ekey := key(program_a, program_b)]
  if (anyDuplicated(coord$ekey) || anyDuplicated(allz$ekey)) fail("duplicate edge keys")
  pairs[, is_coord := ekey %in% coord$ekey]
  pairs[, z_disc_coord := coord$z_partial_standardized[match(ekey, coord$ekey)]]
  pairs[, z_disc_all := allz$z_partial_standardized[match(ekey, allz$ekey)]]
  pairs[, r_disc_all := allz$r_partial[match(ekey, allz$ekey)]]
  pairs[, support := allz$support[match(ekey, allz$ekey)]]
  pairs[, is_c3 := is_coord & coord$r_partial_k3[match(ekey, coord$ekey)] > 0 &
          coord$q_partial_k3[match(ekey, coord$ekey)] < 0.05]
  pairs[is.na(is_c3), is_c3 := FALSE]
  pairs[, is_noncoord := !is_coord & is.finite(z_disc_all)]
  pairs[, is_infnull := support %in% "informative_null"]
  lineage <- registry$cell_type[match(features, registry$program_uid)]
  if (anyNA(lineage)) fail("lineage missing for a frozen program")
  names(lineage) <- features
  pairs[, lineage_a := lineage[program_a]]
  pairs[, lineage_b := lineage[program_b]]
  pairs[, same_lineage := lineage_a == lineage_b]
  pairs[, community_a := communities$community[match(program_a, communities$feature_id)]]
  pairs[, community_b := communities$community[match(program_b, communities$feature_id)]]
  # Frozen-count assertions from the prespec.
  stopifnot(nrow(coord) == 2426L, sum(pairs$is_coord) == 2426L, sum(pairs$is_c3) == 1448L,
            sum(pairs$is_infnull) == 869L, sum(pairs$is_noncoord) == 3769L,
            all(pairs$z_disc_coord[pairs$is_coord] > 0),
            sum(pairs$is_coord & pairs$same_lineage) == 650L)
  ngenes <- membership[program_uid %in% features & as.numeric(original_l1_weight) != 0 &
                         !is.na(mapped_symbol) & trimws(mapped_symbol) != "",
                       .(n_genes = uniqueN(toupper(trimws(mapped_symbol)))), by = .(feature_id = program_uid)]
  list(communities = communities, coord = coord, registry = registry, features = features,
       weights = weights, pairs = pairs, lineage = lineage,
       n_genes = setNames(ngenes$n_genes, ngenes$feature_id)[features],
       annotation = fread(gene_annotation, select = c("gene_id", "gene_name")),
       universe = sub("\\..*$", "", rownames(readRDS(holdout_dge))))
}

# ------------------------------------------------------------------ cohorts
read_char <- function(path) fread(path, colClasses = "character", na.strings = "")

num_or_na <- function(x) {
  x <- trimws(x)
  bad <- is.na(x) | x == "" | toupper(x) %in% c("NA", "NAN", "NONE", "NULL")
  out <- rep(NA_real_, length(x))
  out[!bad] <- as.numeric(x[!bad])
  if (anyNA(out[!bad])) fail("non-numeric value in a numeric field")
  out
}

# Set definitions. Columns resolved by name; literal "NA"/"None" handled explicitly.
build_set_meta <- function() {
  a <- read_char(file.path(H54_ROOT, "data/GSE281797/metadata/sample_phenotypes_resolved.tsv"))
  a[, `:=`(fib = num_or_na(fibrosis), nasv = num_or_na(nas))]
  if (nrow(a) != 94L || uniqueN(a$sample_id) != 94L) fail("GSE281797 must have 94 unique rows")
  if (anyNA(a$fib) || anyNA(a$nasv)) fail("GSE281797 fibrosis/NAS missing")
  if (!all(a$sex_concordant %in% c("True", "False"))) fail("unexpected sex_concordant values")
  disc <- sort(a[sex_concordant == "False", sample_id])
  if (!identical(disc, c("S371", "S373", "S385"))) fail("discordant set differs: ", paste(disc, collapse = ","))
  mkA <- function(d, sexcol) d[, .(sample_id, matrix_column, dataset = "GSE281797",
                                    fibrosis_stage = as.integer(fib), nas_score = nasv,
                                    sex_final = get(sexcol), diagnosis)]
  A94 <- mkA(a, "sex_deposited")
  A91 <- mkA(a[sex_concordant == "True"], "sex_expressed")
  stopifnot(all(a[sex_concordant == "True", sex_expressed == sex_deposited]))
  A66 <- A91[fibrosis_stage <= 1L]
  A37 <- A66[diagnosis != "No Pathology"]

  b <- read_char(file.path(H54_ROOT, "data/GSE225740/metadata/sample_phenotypes_resolved.tsv"))
  b[, `:=`(fib = num_or_na(fibrosis), nasv = num_or_na(nas))]
  if (nrow(b) != 93L || uniqueN(b$pt_id) != 93L) fail("GSE225740 must have 93 unique rows")
  if (sum(is.na(b$fib)) != 6L || anyNA(b$nasv)) fail("GSE225740 must have exactly 6 NA fibrosis, 0 NA NAS")
  if (any(!is.na(b$sex_deposited) & b$sex_deposited != "")) fail("GSE225740 has deposited sex; prespec says none")
  # Cross-check against the deposited phenotype CSV (literal NA strings there).
  dep <- fread(file.path(H54_ROOT, "data/GSE225740/GSE225740_NAFLD_phenotypes.csv.gz"),
               colClasses = "character", na.strings = "")
  setnames(dep, "pt. ID", "pt_id")
  dep <- dep[match(b$pt_id, pt_id)]
  if (!identical(num_or_na(dep$fibrosis), b$fib) || !identical(num_or_na(dep$nas), b$nasv))
    fail("GSE225740 resolved phenotypes disagree with the deposited CSV")
  b[, sample_id := paste0("pt", pt_id)]
  mkB <- function(d) d[, .(sample_id, matrix_column, dataset = "GSE225740",
                            fibrosis_stage = as.integer(fib), nas_score = nasv,
                            sex_final = sex_expressed, diagnosis = NA_character_)]
  B87 <- mkB(b[!is.na(fib)])
  B53 <- B87[fibrosis_stage <= 1L]
  B93 <- mkB(b)

  sets <- list(A91 = A91, B87 = B87, A66 = A66, B53 = B53, A94 = A94, A37 = A37)
  for (s in names(sets)) {
    if (nrow(sets[[s]]) != EXPECTED_N[[s]]) fail(s, " n = ", nrow(sets[[s]]), " != ", EXPECTED_N[[s]])
    if (anyNA(sets[[s]]$sex_final) || any(sets[[s]]$sex_final == "")) fail(s, " missing sex")
  }
  list(sets = sets, A_all = A94, B_all = B93)
}

read_counts <- function(path) {
  cts <- fread(path, header = TRUE)
  ids <- sub("\\..*$", "", cts[[1]])
  if (anyDuplicated(ids)) fail("duplicate gene ids in ", path)
  M <- as.matrix(cts[, -1, with = FALSE]); rownames(M) <- ids
  storage.mode(M) <- "double"
  M
}

# Expression path (prespec "Modality and method", steps 1-3), per analysis set.
build_set_expression <- function(counts, meta, frozen) {
  M <- counts[, meta$matrix_column, drop = FALSE]
  colnames(M) <- meta$sample_id
  keep <- rownames(M) %in% frozen$universe
  Mu <- M[keep, , drop = FALSE]
  dge <- calcNormFactors(DGEList(Mu), method = "TMM")
  lib <- cpm(dge, log = TRUE, prior.count = 1)
  sm <- collapse_symbols(lib, rownames(lib), frozen$annotation)
  panel <- intersect(sort(unique(frozen$weights$gene_symbol)), rownames(sm))
  list(PM = sm[panel, meta$sample_id, drop = FALSE], panel = panel,
       n_universe_present = sum(keep), qc = library_qc(Mu),
       norm_factors = dge$samples$norm.factors)
}

library_qc <- function(Mu) {
  lib <- colSums(Mu)
  top100 <- apply(Mu, 2L, function(v) sum(sort(v, decreasing = TRUE)[seq_len(100L)])) / lib
  data.table(sample_id = colnames(Mu), library_size_universe = lib, top100_share = top100)
}

# H105: the holdout path. pooled = TRUE uses the stored pooled TMM factors (PC1);
# pooled = FALSE recomputes TMM on the 105 first-biopsy columns (PC1b).
build_holdout <- function(frozen, pooled = TRUE) {
  dge <- readRDS(holdout_dge)
  meta <- as.data.table(readRDS(holdout_meta))
  cw <- fread(gse193066_crosswalk)
  meta <- merge(meta, cw[, .(sample_id = run_id, participant_token, is_first_biopsy)],
                by = "sample_id", all.x = TRUE)
  meta <- meta[is_first_biopsy == TRUE & !is.na(fibrosis_stage) & !is.na(nas_score)]
  if (anyDuplicated(meta$participant_token)) fail("holdout participant duplicated")
  if (nrow(meta) != 105L) fail("H105 n = ", nrow(meta))
  if (pooled) {
    lib <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    lib <- lib[, meta$sample_id, drop = FALSE]
  } else {
    sub <- DGEList(counts = dge$counts[, meta$sample_id, drop = FALSE])
    if (max(abs(sub$samples$lib.size - dge$samples[meta$sample_id, "lib.size"])) > 1e-6)
      fail("holdout lib.size is not the column sum")
    sub <- calcNormFactors(sub, method = "TMM")
    lib <- edgeR::cpm(sub, log = TRUE, prior.count = 1)
  }
  sm <- collapse_symbols(lib, rownames(lib), frozen$annotation)
  panel <- intersect(sort(unique(frozen$weights$gene_symbol)), rownames(sm))
  m <- meta[, .(sample_id, dataset, fibrosis_stage = as.integer(fibrosis_stage),
                nas_score = as.numeric(nas_score), sex_final, diagnosis = diagnosis_harmonized)]
  list(meta = m, PM = sm[panel, m$sample_id, drop = FALSE], panel = panel,
       n_universe_present = nrow(dge), qc = library_qc(dge$counts[, m$sample_id, drop = FALSE]))
}

# --------------------------------------------------------------- statistics
# Block + operators for a set of rows (idx into meta). Duplicated participants get
# unique column names, as in gate job 3131.
make_block <- function(PM, meta, idx, frozen) {
  m <- copy(meta[idx])
  m[, orig := sample_id]
  m[, sample_id := paste0(orig, "__", seq_len(.N))]
  X <- PM[, m$orig, drop = FALSE]; colnames(X) <- m$sample_id
  tryCatch(t4_build_block(X, m, frozen$weights, frozen$features, program_weight_coverage),
           error = function(e) NULL)
}

ops_for <- function(block, frozen, panel, cache = NULL) {
  if (is.null(cache)) return(t4_disjoint_operators(block, frozen$pairs, frozen$features, 0.8))
  gi <- as.numeric(match(block$genes, panel))
  if (anyNA(gi)) fail("block gene outside panel")
  gkey <- paste(length(gi), sum(gi), sum(gi^2), sum(gi * seq_along(gi)),
                paste(as.integer(block$coverage$cohort_testable), collapse = ""))
  ops <- cache[[gkey]]
  if (is.null(ops)) {
    ops <- t4_disjoint_operators(block, frozen$pairs, frozen$features, 0.8)
    assign(gkey, ops, envir = cache)
  }
  ops
}

# Which pairs CAN carry a finite statistic given testability and leave-shared-out
# eligibility. Uses no correlation, so it is safe before the outcome stage.
eligible_pairs <- function(block, ops, frozen) {
  tst <- block$coverage$cohort_testable[match(frozen$features, block$coverage$feature_id)]
  ok <- tst[frozen$pairs$i] & tst[frozen$pairs$j]
  inel <- ops$pair_index[!ops$eligible]
  ok[inel] <- FALSE
  ok
}

# Pair-index lookup matrix: pidx[i, j] = row of the pair table.
pair_index_matrix <- function(frozen) {
  nf <- length(frozen$features)
  P <- matrix(NA_integer_, nf, nf)
  P[cbind(frozen$pairs$i, frozen$pairs$j)] <- frozen$pairs$pair_index
  P[cbind(frozen$pairs$j, frozen$pairs$i)] <- frozen$pairs$pair_index
  P
}

# S1 for one edge mask and z vector.
s1_frac <- function(z, mask) {
  ok <- mask & is.finite(z)
  c(s1 = mean(z[ok] > 0), n = sum(ok))
}

# Label-permutation S1 for every stored permutation (rows of perms).
# Edge set = mask & finite observed z. Mapped pairs that are not finite are dropped.
perm_s1 <- function(z, mask, perms, pidx, frozen) {
  e <- which(mask & is.finite(z))
  ia <- frozen$pairs$i[e]; ib <- frozen$pairs$j[e]
  vapply(seq_len(nrow(perms)), function(r) {
    p <- perms[r, ]
    zz <- z[pidx[cbind(p[ia], p[ib])]]
    ok <- is.finite(zz)
    mean(zz[ok] > 0)
  }, numeric(1))
}

perm_nonfinite_fraction <- function(okpair, mask, perms, pidx, frozen) {
  e <- which(mask & okpair)
  ia <- frozen$pairs$i[e]; ib <- frozen$pairs$j[e]
  mean(vapply(seq_len(nrow(perms)), function(r) {
    p <- perms[r, ]
    mean(!okpair[pidx[cbind(p[ia], p[ib])]])
  }, numeric(1)))
}

# Everything about one realisation: S1_k1, S1_k3, label-null means over perms 1..200,
# E_k1, E_k3, and the holdout-style non-coordinated concordance.
evaluate_z <- function(z1, z3, perms, pidx, frozen, n_perm_e = N_PERM_E, all_perms = FALSE) {
  pr <- frozen$pairs
  a1 <- s1_frac(z1, pr$is_coord); a3 <- s1_frac(z3, pr$is_c3)
  P <- if (all_perms) perms else perms[seq_len(n_perm_e), , drop = FALSE]
  n1 <- perm_s1(z1, pr$is_coord, P, pidx, frozen)
  n3 <- perm_s1(z3, pr$is_c3, P, pidx, frozen)
  m1 <- mean(n1[seq_len(n_perm_e)]); m3 <- mean(n3[seq_len(n_perm_e)])
  okn <- is.finite(z1) & pr$is_noncoord
  base <- mean(sign(z1[okn]) == sign(pr$z_disc_all[okn]))
  okc <- is.finite(z1) & pr$is_coord
  s1_sign <- mean(sign(z1[okc]) == sign(pr$z_disc_coord[okc]))
  out <- list(s1_k1 = a1[["s1"]], s1_k3 = a3[["s1"]], n_c1 = a1[["n"]], n_c3 = a3[["n"]],
              null_mean_k1 = m1, null_mean_k3 = m3, e_k1 = a1[["s1"]] - m1, e_k3 = a3[["s1"]] - m3,
              base_noncoord = base, n_noncoord = sum(okn), s1_k1_signmatch = s1_sign)
  if (all_perms) { out$null_k1 <- n1; out$null_k3 <- n3 }
  out
}

p_perm_upper <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))

# Within-lineage bijections of the testable programs. Untestable programs map to
# themselves; their pairs are non-finite and drop out.
draw_lineage_perms <- function(n, testable, strata, features) {
  out <- matrix(seq_along(features), n, length(features), byrow = TRUE)
  idx_t <- which(features %in% testable)
  groups <- split(idx_t, strata[idx_t])
  for (r in seq_len(n)) {
    for (g in groups) {
      out[r, g] <- if (length(g) > 1L) g[sample.int(length(g))] else g
    }
  }
  colnames(out) <- features
  out
}

# Bootstrap: participants resampled with replacement within fibrosis stage.
draw_stratified_boot <- function(n_draw, stage) {
  lv <- sort(unique(stage))
  t(vapply(seq_len(n_draw), function(b) {
    out <- integer(length(stage))
    for (l in lv) {
      w <- which(stage == l)
      out[w] <- w[sample.int(length(w), length(w), replace = TRUE)]
    }
    out
  }, integer(length(stage))))
}

design_ok <- function(meta_rows, formula = T4_RESIDUAL_FORMULA) {
  d <- copy(meta_rows)
  d[, `:=`(fibrosis_factor = droplevels(factor(fibrosis_stage)),
           nas_factor = droplevels(t4_nas_bin(nas_score)),
           sex_final = droplevels(factor(sex_final)))]
  X <- tryCatch(stats::model.matrix(formula, data = d), error = function(e) NULL)
  if (is.null(X)) return(c(ok = 0, cols = NA, rank = NA, df = NA))
  r <- qr(X)$rank
  c(ok = as.numeric(r == ncol(X) && nrow(X) > ncol(X) + 3L), cols = ncol(X), rank = r,
    df = nrow(X) - ncol(X))
}

# Leading residual components without computing any edge statistic.
residual_components <- function(block, k = 3L) {
  U <- t4_residual_gene_matrix(block)
  P <- U %*% block$W
  Ps <- scale(P); Ps[!is.finite(Ps)] <- 0
  sv <- svd(Ps, nu = k, nv = 0L)
  list(share = sv$d[[1L]]^2 / sum(sv$d^2), u = sv$u[, seq_len(k), drop = FALSE])
}
