#!/usr/bin/env Rscript
# validity-L10-01: discovery-cohort-only checks of the proposed S1 statistic.
# GSE213621 is NOT read for any edge statistic (only discovery disease units).
# Reuses job 3811 helper copies (read only).
#  1. reproduce job 3811 result_D (coarse LOCO S1)
#  2. split S1 and its edge-draw null by edge category
#  3. add library-size and library-complexity covariates to the residual model
#  4. donor bootstrap of S1 and of the excess in one held-out cohort
suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix)})
options(warn = 1)
SEED <- 20261015L
set.seed(SEED)
NBOOT <- as.integer(Sys.getenv("NBOOT", "0"))
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
j3811 <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3811")
out <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-01/out")
dir.create(out, showWarnings = FALSE)
source(file.path(j3811, "analysis_lib.R"))
source(file.path(j3811, "t4_interaction_lib.R"))

cand <- file.path(root, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
T4 <- file.path(cand, "workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction")
rel <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
prog_root <- file.path(root, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
annotation <- fread(file.path(cand, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
                    select = c("gene_id", "gene_name"))
communities <- fread(file.path(T4, "program_communities.tsv"))
features <- communities$feature_id
membership <- fread(file.path(prog_root, "program_membership_v2.tsv"))
weights <- t4_membership_weights(membership, features)
pairs0 <- t4_pair_table(features, weights)
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
pairs0[, ekey := key(program_a, program_b)]
fe <- fread(file.path(T4, "program_interaction_edges.tsv.gz"),
            select = c("program_a", "program_b", "cell_type_a", "cell_type_b", "module_name_a", "module_name_b", "same_cell_type"))
fe[, ekey := key(program_a, program_b)]
# [Inference] module-name rule for technical / housekeeping axes in bulk
hk_rx <- "Translation|ribosome|Oxidative phosphorylation|Mitochondrial|Immediate-early|AP-1|housekeeping"
fe[, hk := grepl(hk_rx, module_name_a) | grepl(hk_rx, module_name_b)]
fe[, cat2 := ifelse(same_cell_type, "same_lineage", "cross_lineage")]
cat("programs matching housekeeping rule:",
    uniqueN(c(fe[grepl(hk_rx, module_name_a), program_a], fe[grepl(hk_rx, module_name_b), program_b])), "\n")

dge <- readRDS(file.path(rel, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
disc <- c("GSE130970", "GSE135251", "GSE162694")
md <- meta[dataset %in% disc & group_binary == "Disease"]
cnt <- dge$counts[, md$sample_id]
libsize <- colSums(cnt)
pr <- sweep(cnt, 2, libsize, "/")
ent <- apply(pr, 2, function(p) { p <- p[p > 0]; -sum(p * log(p)) })
top100 <- apply(pr, 2, function(p) sum(sort(p, decreasing = TRUE)[1:100]))
md[, `:=`(log_lib = log10(libsize[sample_id]), entropy = ent[sample_id], top100 = top100[sample_id])]
rm(pr); gc()
lc <- edgeR::cpm(dge[, md$sample_id], log = TRUE, prior.count = 1)
sm <- collapse_symbols(lc, rownames(lc), annotation); rm(lc, dge, cnt); gc()
panel <- intersect(sort(unique(weights$gene_symbol)), rownames(sm))
sm <- sm[panel, , drop = FALSE]
cat("discovery disease units", nrow(md), "; panel", length(panel), "\n")
print(md[, .(n = .N, med_log_lib = median(log_lib), med_entropy = median(entropy), med_top100 = median(top100)), by = dataset])

bin3 <- function(f) ifelse(f <= 1, 1L, ifelse(f == 2, 2L, 3L))
md[, fibrosis_stage := bin3(fibrosis_stage)]
COARSE <- ~ sex_final + fibrosis_factor
COARSE_LIB <- ~ sex_final + fibrosis_factor + log_lib + entropy + top100
make_block <- function(mc, formula) { T4_RESIDUAL_FORMULA <<- formula; t4_build_block(sm, mc, weights, features, 0.80) }

approx_discovery <- function(blocks) {
  pairs <- copy(pairs0)
  ops <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features, min_retained_weight = 0.80)
  elig <- Reduce(`&`, lapply(ops, function(op) { e <- rep(TRUE, nrow(pairs)); e[op$pair_index] <- op$eligible; e }))
  obs <- t4_one_pass(blocks, pairs, ops, want_matrix = TRUE)
  w_c <- vapply(blocks, function(b) b$df - 3L, numeric(1))
  wsum <- function(zc) rowSums(matrix(w_c, nrow(zc), ncol(zc), byrow = TRUE) * is.finite(zc))
  pairs[, `:=`(entangled = !elig, z_disj = obs$z_disjoint, z_part = obs$z_partial_k1,
               zs_disj = obs$z_disjoint * sqrt(wsum(obs$zc_disjoint)),
               zs_part = obs$z_partial_k1 * sqrt(wsum(obs$zc_partial[[1L]])))]
  fam <- pairs[, !entangled & is.finite(zs_disj) & is.finite(zs_part)]
  pairs[, `:=`(q_disj = NA_real_, q_part = NA_real_)]
  pairs[fam, q_disj := p.adjust(2 * pnorm(-abs(zs_disj)), "BH")]
  pairs[fam, q_part := p.adjust(2 * pnorm(-abs(zs_part)), "BH")]
  pairs[, coordinated := fam & q_disj < 0.05 & q_part < 0.05 & z_part > 0]
  pairs[, in_family := fam]
  pairs
}
replicate_z <- function(block) {
  pairs <- copy(pairs0)
  ops <- list(t4_disjoint_operators(block, pairs, features, min_retained_weight = 0.80))
  obs <- t4_one_pass(list(block), pairs, ops)
  pairs[, z_rep := obs$z_partial_k1]
  pairs[, .(ekey, z_rep)]
}
s1_by <- function(d, r, n_draws = 1000L) {
  m <- merge(d[in_family == TRUE, .(ekey, coordinated, z_part)], r, by = "ekey")
  m <- merge(m[is.finite(z_rep) & is.finite(z_part)], fe[, .(ekey, cat2, hk)], by = "ekey")
  m[, agree := sign(z_rep) == sign(z_part)]
  one <- function(mm, label) {
    co <- mm[coordinated == TRUE]; nc <- mm[coordinated == FALSE]; ncpos <- nc[z_part > 0]
    if (nrow(co) < 20 || nrow(nc) < nrow(co)) {
      nul <- rep(NA_real_, 2)
    } else {
      nd <- replicate(n_draws, mean(nc$agree[sample.int(nrow(nc), nrow(co))]))
      nul <- c(mean(nd), sd(nd))
    }
    data.table(subset = label, n_coord = nrow(co), n_noncoord = nrow(nc), s1 = mean(co$agree),
               null_mean = nul[1], null_sd = nul[2], excess = mean(co$agree) - nul[1],
               noncoord_pos_concord = mean(ncpos$agree), n_noncoord_pos = nrow(ncpos),
               coord_mean_abs_r_rep = mean(abs(tanh(co$z_rep))))
  }
  rbind(one(m, "all"), one(m[cat2 == "same_lineage"], "same_lineage"), one(m[cat2 == "cross_lineage"], "cross_lineage"),
        one(m[hk == TRUE], "housekeeping_involved"), one(m[hk == FALSE], "no_housekeeping"),
        one(m[cat2 == "cross_lineage" & hk == FALSE], "cross_lineage_no_housekeeping"))
}

# ---------------- composition proxies from lineage markers OUTSIDE program membership [Inference: marker choice]
mk <- list(chol = c("KRT7", "EPCAM", "SOX9", "KRT23", "MMP7", "HNF1B"),
           fib = c("COL1A1", "COL1A2", "COL3A1", "LUM"),
           endo = c("VWF", "CDH5"), bcell = c("MS4A1", "CD79A"))
stopifnot(!any(unlist(mk) %in% weights$gene_symbol))
dge2 <- readRDS(file.path(rel, "merged_dge.rds"))
lc2 <- edgeR::cpm(dge2[, md$sample_id], log = TRUE, prior.count = 1)
sym2 <- annotation$gene_name[match(sub("\\..*$", "", rownames(lc2)), sub("\\..*$", "", annotation$gene_id))]
rm(dge2); gc()
for (n in names(mk)) {
  rows <- which(sym2 %in% mk[[n]])
  cat(n, "marker rows found:", length(rows), "\n")
  md[, (paste0("prx_", n)) := 0]
  for (cc in disc) {
    ids <- md[dataset == cc, sample_id]
    z <- t(scale(t(lc2[rows, ids, drop = FALSE])))
    md[dataset == cc, (paste0("prx_", n)) := colMeans(z, na.rm = TRUE)]
  }
}
rm(lc2); gc()
COARSE_COMP <- ~ sex_final + fibrosis_factor + prx_chol + prx_fib + prx_endo + prx_bcell

part <- setNames(communities$community, communities$feature_id)
lin <- unique(rbind(fe[, .(f = program_a, l = cell_type_a)], fe[, .(f = program_b, l = cell_type_b)]))
lin <- setNames(lin$l, lin$f)[names(part)]
c3 <- names(part)[part == 3]; c3_lin <- table(lin[c3])
pk <- function(a, b) { k <- key(a, b); k }
s23 <- function(r) {
  g <- merge(pairs0[, .(ekey, program_a, program_b)], r, by = "ekey")[is.finite(z_rep) & z_rep > 0]
  mod_of <- function(assign) {
    ca <- assign[g$program_a]; cb <- assign[g$program_b]; w <- g$z_rep; tot <- sum(w)
    ki <- tapply(c(w, w), c(ca, cb), sum)[as.character(sort(unique(assign)))]; ki[is.na(ki)] <- 0
    sum(w[ca == cb]) / tot - sum((ki / (2 * tot))^2)
  }
  s2 <- mod_of(part)
  nul_g <- replicate(200L, mod_of(setNames(sample(part), names(part))))
  nul_l <- replicate(200L, { p2 <- part; for (L in unique(lin)) { ix <- which(lin == L); p2[ix] <- p2[ix][sample.int(length(ix))] }; mod_of(p2) })
  dens <- function(set) nrow(g[program_a %in% set & program_b %in% set]) / (length(set) * (length(set) - 1) / 2)
  s3 <- dens(c3)
  nul3 <- replicate(1000L, dens(unlist(lapply(names(c3_lin), function(L) { pool <- names(lin)[lin == L]; sample(pool, c3_lin[[L]]) }))))
  data.table(s2 = s2, s2_null_label = mean(nul_g), s2_p_label = (1 + sum(nul_g >= s2)) / 201,
             s2_null_within_lineage = mean(nul_l), s2_p_within_lineage = (1 + sum(nul_l >= s2)) / 201,
             s3 = s3, s3_null = mean(nul3), s3_null_sd = sd(nul3), s3_p = (1 + sum(nul3 >= s3)) / 1001)
}
res <- list()
for (des in c("COARSE", "COARSE_COMP")) {
  bl <- lapply(split(md, by = "dataset"), function(mc) make_block(mc, get(des)))
  cat(des, ": df", vapply(bl, `[[`, numeric(1), "df"), "\n")
  for (h in disc) {
    set.seed(SEED)
    d <- approx_discovery(bl[setdiff(disc, h)])
    r <- replicate_z(bl[[h]])
    s1 <- s1_by(d, r)[subset %in% c("all", "cross_lineage", "cross_lineage_no_housekeeping")]
    s1w <- dcast(s1[, .(k = 1, subset, s1, excess)], k ~ subset, value.var = c("s1", "excess"))[, k := NULL]
    set.seed(SEED)
    res[[paste(des, h)]] <- cbind(data.table(design = des, held_out = h, n = bl[[h]]$n,
                                             n_coord = s1[subset == "all", n_coord]), s1w, s23(r))
  }
}
R <- rbindlist(res, fill = TRUE)
options(width = 250)
print(R[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
fwrite(R, file.path(out, "loco_composition_proxy.tsv"), sep = "\t")
# how much do proxies relate to fibrosis bin within cohort (residual severity carriers)
print(md[, .(r_chol_bin = cor(prx_chol, fibrosis_stage, method = "spearman"),
             r_fib_bin = cor(prx_fib, fibrosis_stage, method = "spearman")), by = dataset])
cat("DONE\n")
