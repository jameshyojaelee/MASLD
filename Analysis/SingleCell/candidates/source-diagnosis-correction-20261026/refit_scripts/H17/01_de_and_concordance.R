#!/usr/bin/env Rscript
# H17 step 1: limma-voom in both systems, positive control, non-responders,
# S1/S2 point/S3, core, matched-set and sign-permutation nulls, Hotspot and
# bulk placement. Usage: Rscript 01_de_and_concordance.R <HYP_DIR> <OUT_DIR>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
hyp_dir <- args[1]; out_dir <- args[2]
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "h17_lib.R"))
stopifnot(file.exists(file.path(hyp_dir, "PRESPEC_H17.md")), dir.exists(out_dir))
wt <- function(x, f) fwrite(x, file.path(out_dir, f), sep = "\t")
N_PERM <- 10000L; N_SETS <- 10000L

# ---------------------------------------------------------------- iPSC
A <- load_ipsc(hyp_dir)
wt(A$drop_log, "ipsc_symbol_drop_log.tsv")
wt(A$sheet, "ipsc_sample_sheet.tsv")
arm_A <- A$sheet[, .(n_participants = uniqueN(participant)), by = treatment]
wt(arm_A, "ipsc_arm_counts.tsv")
if (any(arm_A$n_participants < 3)) {
  writeLines("iPSC system has < 3 participants in an arm", file.path(out_dir, "STOP_untestable.txt"))
  quit(status = 0)
}
genes_A <- filter_genes(A$counts, A$sheet$treatment)
des_A <- ipsc_design(A$sheet$participant, A$sheet$treatment)
de_A <- fit_full(A$counts, genes_A, des_A, "treatmentHepFA")
wt(de_A, "ipsc_de.tsv")
message("iPSC tested genes: ", nrow(de_A), "; q<0.05: ", sum(de_A$q < Q_THR))

# positive control
pc <- data.table(gene = POS_CONTROLS)
pc[, expressed := gene %in% de_A$gene]
pc <- merge(pc, de_A[, .(gene, logFC, q)], by = "gene", all.x = TRUE, sort = FALSE)
pc[, up_sig := expressed & !is.na(q) & q < Q_THR & logFC > 0]
pc_pass <- isTRUE(pc[gene == "PLIN2", up_sig]) &&
  sum(pc$up_sig) >= ceiling(sum(pc$expressed) / 2)
pc[, overall_pass := pc_pass]
wt(pc, "positive_control_ipsc.tsv")
if (!pc_pass) {
  writeLines("Positive control failed in iPSC system: not informative",
             file.path(out_dir, "STOP_positive_control_failed.txt"))
  quit(status = 0)
}

# non-responders (PLIN2 within-pair change in TMM log2-CPM)
dgeA <- calcNormFactors(DGEList(counts = A$counts[genes_A, ]))
lcpm <- cpm(dgeA, log = TRUE, prior.count = 2)
sh <- A$sheet
plin <- data.table(participant = sh$participant, treatment = sh$treatment,
                   v = lcpm["PLIN2", sh$library])
nr <- dcast(plin, participant ~ treatment, value.var = "v")
nr[, delta_PLIN2_log2cpm := HepFA - HepC]
nr[, non_responder := delta_PLIN2_log2cpm <= 0]
wt(nr, "ipsc_nonresponders.tsv")
rm(dgeA, lcpm)

# ---------------------------------------------------------------- in vivo
B <- load_invivo("steatotic")
wt(B$roster_log, "invivo_roster_log.tsv")
wt(B$meta[, .(sample, dataset, disease_stage_coarse, group, n_hep, n_runs, sex)], "invivo_roster.tsv")
arm_B <- B$meta[, .N, by = .(dataset, group)]
wt(arm_B, "invivo_arm_counts.tsv")
if (any(B$meta[, .N, by = group]$N < 3) || uniqueN(B$meta$group) < 2) {
  writeLines("In vivo system has < 3 donors in an arm", file.path(out_dir, "STOP_untestable.txt"))
  quit(status = 0)
}
genes_B <- filter_genes(B$counts, B$meta$group)
des_B <- invivo_design(B$meta)
de_B <- fit_full(B$counts, genes_B, des_B, "groupsteatotic")
wt(de_B, "invivo_de.tsv")
message("in vivo tested genes: ", nrow(de_B), "; q<0.05: ", sum(de_B$q < Q_THR))
pcB <- de_B[gene %in% POS_CONTROLS, .(gene, logFC, q)]
wt(pcB, "positive_control_genes_invivo_descriptive.tsv")

# ---------------------------------------------------------------- joint
U <- intersect(de_A$gene, de_B$gene)
j <- merge(de_A[gene %in% U, .(gene, logFC_A = logFC, AveExpr_A = AveExpr, t_A = t, q_A = q)],
           de_B[gene %in% U, .(gene, logFC_B = logFC, AveExpr_B = AveExpr, t_B = t, q_B = q)],
           by = "gene")
setorder(j, gene)
j[, sig_A := q_A < Q_THR][, sig_B := q_B < Q_THR]
j[, sA := sign(logFC_A)][, sB := sign(logFC_B)]
j[, core := sig_A & sig_B & sA == sB]

# S1 (A = iPSC) and reverse
idx_A <- which(j$sig_A); idx_B <- which(j$sig_B)
s1_obs <- sign_agreement(j$sA[idx_A], j$sB[idx_A])
s1_null <- sign_perm_null(j$sA[idx_A], j$sB, idx_A, N_PERM)
s1r_obs <- if (length(idx_B)) sign_agreement(j$sB[idx_B], j$sA[idx_B]) else NA_real_
s1r_null <- if (length(idx_B)) sign_perm_null(j$sB[idx_B], j$sA, idx_B, N_PERM) else NA_real_
# S2, S3
rho <- cor(j$logFC_A, j$logFC_B, method = "spearman")
nU <- nrow(j); nA <- length(idx_A); nB <- length(idx_B); k <- sum(j$sig_A & j$sig_B)
s3_p <- phyper(k - 1, nA, nU - nA, nB, lower.tail = FALSE)
core <- j[core == TRUE]

# matched-set negative control for S1
qbin <- function(x) cut(x, breaks = unique(quantile(x, 0:5 / 5)), include.lowest = TRUE, labels = FALSE)
j[, bin := paste(qbin(AveExpr_A), qbin(AveExpr_B), sep = "_")]
same <- j$sA == j$sB
sets_GA <- matched_sets(j$bin, idx_A, N_SETS)
mset_agree <- vapply(sets_GA, function(ix) mean(same[ix]), numeric(1))

conc <- data.table(
  n_U = nU, n_sig_ipsc_in_U = nA, n_sig_invivo_in_U = nB, n_sig_both = k,
  n_core = nrow(core), n_core_up = sum(core$sA > 0), n_core_down = sum(core$sA < 0),
  invivo_up_fraction_U = mean(j$sB > 0),
  S1_agreement = s1_obs, S1_signperm_null_mean = mean(s1_null),
  S1_p_upper = p_upper(s1_obs, s1_null), S1_p_lower = p_lower(s1_obs, s1_null),
  S1_matchedset_null_mean = mean(mset_agree), S1_matchedset_q95 = unname(quantile(mset_agree, 0.95)),
  S1_matchedset_p_upper = p_upper(s1_obs, mset_agree),
  S1rev_agreement = s1r_obs,
  S1rev_signperm_null_mean = if (length(idx_B)) mean(s1r_null) else NA_real_,
  S1rev_p_upper = if (length(idx_B)) p_upper(s1r_obs, s1r_null) else NA_real_,
  S2_spearman = rho, S3_hypergeom_p = s3_p,
  S3_expected_overlap = nA * nB / nU)
wt(conc, "concordance_primary.tsv")
wt(data.table(perm = seq_len(N_PERM), S1_null = s1_null), "S1_signperm_null.tsv.gz")
wt(data.table(set = seq_len(N_SETS), agreement = mset_agree), "S1_matchedset_null.tsv.gz")

# ---------------------------------------------------------------- Hotspot programs
mem <- fread(PATHS$membership)[!is.na(mapped_symbol) & nzchar(mapped_symbol)]
reg <- fread(PATHS$registry, select = c("program_uid", "cell_type", "module", "module_name"))
progs <- reg$program_uid
stopifnot(length(progs) == 117L, setequal(progs, unique(fread(PATHS$membership)$program_uid)))
W <- matrix(0, nrow = nU, ncol = length(progs), dimnames = list(j$gene, progs))
m_in <- mem[mapped_symbol %in% j$gene]
m_in <- m_in[, .(w = sum(original_l1_weight)), by = .(program_uid, mapped_symbol)]
W[cbind(match(m_in$mapped_symbol, j$gene), match(m_in$program_uid, progs))] <- m_in$w
n_mem_U <- colSums(W > 0)

program_test <- function(target_idx, label) {
  if (length(target_idx) == 0L) return(NULL)
  obs <- colMeans(W[target_idx, , drop = FALSE])
  sets <- matched_sets(j$bin, target_idx, N_SETS)
  nullm <- t(vapply(sets, function(ix) colMeans(W[ix, , drop = FALSE]), numeric(ncol(W))))
  p <- vapply(seq_along(obs), function(i) p_upper(obs[i], nullm[, i]), numeric(1))
  data.table(set = label, program_uid = progs, n_members_in_U = n_mem_U,
             n_target_members = colSums(W[target_idx, , drop = FALSE] > 0),
             mean_weight_target = obs, null_mean = colMeans(nullm),
             p_emp = p, q_BH_117 = p.adjust(p, "BH"))
}
idx_core <- which(j$core)
prog <- rbindlist(list(
  program_test(idx_core, "core"),
  program_test(which(j$core & j$sA > 0), "core_up"),
  program_test(which(j$core & j$sA < 0), "core_down")))
if (nrow(prog)) {
  prog <- merge(reg, prog, by = "program_uid")
  setorder(prog, set, p_emp)
}
wt(prog, "hotspot_program_tests.tsv")

# ---------------------------------------------------------------- bulk signature
bulk <- fread(PATHS$bulk)
bulk_n_rows <- nrow(bulk)
bulk <- bulk[!is.na(symbol) & nzchar(symbol)]
bulk <- bulk[!(symbol %in% symbol[duplicated(symbol)])]
bulk[, sig := !is.na(padj) & padj < 0.05 & abs(logFC) > 0.5]
j[, in_bulk := gene %in% bulk$symbol]
j <- merge(j, bulk[, .(gene = symbol, logFC_bulk = logFC, padj_bulk = padj, sig_bulk = sig)],
           by = "gene", all.x = TRUE, sort = TRUE)
j <- j[match(rownames(W), gene)]
stopifnot(identical(j$gene, rownames(W)), identical(which(j$core), idx_core))
j[is.na(sig_bulk), sig_bulk := FALSE]
Ub <- j[in_bulk == TRUE]
k_b <- Ub[core == TRUE & sig_bulk == TRUE, .N]
p_b <- phyper(k_b - 1, Ub[, sum(sig_bulk)], nrow(Ub) - Ub[, sum(sig_bulk)],
              Ub[, sum(core)], lower.tail = FALSE)
cs <- Ub[core == TRUE & sig_bulk == TRUE]
agree_b <- if (nrow(cs)) mean(sign(cs$logFC_bulk) == cs$sA) else NA_real_
bulk_null <- if (length(idx_core)) {
  sets <- matched_sets(j$bin, idx_core, N_SETS)
  rbindlist(lapply(sets, function(ix) {
    s <- j[ix][in_bulk == TRUE & sig_bulk == TRUE]
    data.table(n_overlap = nrow(s),
               agree = if (nrow(s)) mean(sign(s$logFC_bulk) == s$sA) else NA_real_)
  }))
} else NULL
bulk_res <- data.table(
  bulk_rows = bulk_n_rows, n_Ub = nrow(Ub), n_bulk_sig_in_Ub = Ub[, sum(sig_bulk)],
  n_core_in_Ub = Ub[, sum(core)], n_core_and_bulk_sig = k_b,
  expected_overlap = Ub[, sum(core)] * Ub[, sum(sig_bulk)] / nrow(Ub),
  hypergeom_p = p_b, core_bulk_sign_agreement = agree_b,
  matchedset_overlap_mean = if (!is.null(bulk_null)) mean(bulk_null$n_overlap) else NA_real_,
  matchedset_overlap_p_upper = if (!is.null(bulk_null)) p_upper(k_b, bulk_null$n_overlap) else NA_real_,
  matchedset_agree_mean = if (!is.null(bulk_null)) mean(bulk_null$agree, na.rm = TRUE) else NA_real_,
  matchedset_agree_p_upper = if (!is.null(bulk_null) && !is.na(agree_b))
    p_upper(agree_b, bulk_null$agree[!is.na(bulk_null$agree)]) else NA_real_)
wt(bulk_res, "bulk_signature_comparison.tsv")

# ---------------------------------------------------------------- outputs
prog_core <- if (nrow(prog)) prog[set == "core"] else NULL
core_tab <- j[core == TRUE, .(gene, direction = ifelse(sA > 0, "up", "down"),
                              logFC_ipsc = logFC_A, q_ipsc = q_A,
                              logFC_invivo = logFC_B, q_invivo = q_B,
                              logFC_bulk, padj_bulk, sig_bulk)]
if (nrow(core_tab)) {
  memb <- mem[mapped_symbol %in% core_tab$gene, .(programs = paste(sort(unique(program_uid)), collapse = ";")),
              by = .(gene = mapped_symbol)]
  core_tab <- merge(core_tab, memb, by = "gene", all.x = TRUE)
}
wt(core_tab, "core_genes.tsv")
wt(j[, !c("bin"), with = FALSE], "joint_table.tsv.gz")
saveRDS(list(genes_A = genes_A, genes_B = genes_B, U = j$gene,
             idx_A = idx_A, sA = j$sA, sB = j$sB, q_A = j$q_A, logFC_A = j$logFC_A,
             s1_obs = s1_obs, rho = rho, nonresp = nr[non_responder == TRUE, as.character(participant)]),
        file.path(out_dir, "state_step1.rds"))
writeLines("ok", file.path(out_dir, "STEP1_DONE"))
message("step 1 done")
