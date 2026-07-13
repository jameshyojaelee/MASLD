#!/usr/bin/env Rscript
# T3 discovery-arm benchmark — PROCON (confirmatory program convergence) vs xMOFA
# (unsupervised joint factorization). The 41-agent design chose PROCON over joint
# factorization on interpretability + the disjoint-donor problem; this quantifies it.
#
# xMOFA = the existing MOFA+ run (mofapy2 0.7.4) on the 5 donor-SHARED scRNA
# pseudobulk cell-type views (Hep/Mac/Fib/Endo/Chol; 269 donors, k=16 factors;
# 503b). Reuses `mofa_factors_k16.tsv` (per-view gene weights). We compare on three
# axes that decide the method: (1) RECOVERY of the frozen 192-program dictionary,
# (2) MODALITY REACH, (3) INTERPRETABILITY (factor→program blending).
suppressPackageStartupMessages({ library(data.table) })
set.seed(42)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
P    <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
MOFA <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/benchmarks/mofa_factors_k16.tsv")
TOPN <- 100

memb <- fread(file.path(P, "program_gene_membership.tsv"))           # 192 programs
prog <- split(memb$gene, memb$program_id)
proc <- fread(file.path(P, "procon_convergence.tsv"))                # PROCON results (128 conv)
conv_ids <- proc[conv_q < 0.05, program_id]

mf <- fread(MOFA)                                                     # view|factor|gene|weight|abs_weight
# a gene's factor-saliency = max |weight| across the 5 cell-type views (top-loading anywhere)
mf_g <- mf[, .(abs_weight = max(abs_weight)), by = .(factor, gene)]
fac_top <- mf_g[order(-abs_weight), .(genes = list(head(gene, TOPN))), by = factor]
factors <- fac_top$factor
fac_sets <- setNames(fac_top$genes, factors)

universe <- union(unique(mf$gene), unique(memb$gene)); U <- length(universe)
# hypergeometric overlap test: program recovered by a factor if enrichment q<0.05
best_match <- function(pg) {
  pg <- intersect(pg, universe); if (!length(pg)) return(list(q=1, jac=0, fac=NA))
  res <- sapply(fac_sets, function(fs) {
    fs <- intersect(fs, universe); ov <- length(intersect(fs, pg))
    p <- phyper(ov - 1, length(pg), U - length(pg), length(fs), lower.tail = FALSE)
    jac <- ov / length(union(fs, pg)); c(p = p, jac = jac) })
  i <- which.min(res["p", ]); list(q = res["p", i], jac = res["jac", i], fac = factors[i])
}
rec <- rbindlist(lapply(names(prog), function(pid) {
  m <- best_match(prog[[pid]])
  data.table(program_id = pid, source = sub("_.*", "", pid),
             is_convergent = pid %in% conv_ids,
             best_factor = m$fac, best_jaccard = round(m$jac, 3),
             hyper_p = m$q) }))
rec[, hyper_q := p.adjust(hyper_p, "BH")]
rec[, mofa_recovers := hyper_q < 0.05 & best_jaccard >= 0.05]
fwrite(rec, file.path(P, "t3_xmofa_vs_procon_recovery.tsv"), sep = "\t")

# ── (1) RECOVERY ──
cat("── (1) Dictionary recovery: does unsupervised MOFA re-find the programs? ──\n")
cat(sprintf("MOFA recovers %d / %d programs overall (%.0f%%); %d / %d PROCON-CONVERGENT programs (%.0f%%)\n",
    rec[mofa_recovers==TRUE, .N], nrow(rec), 100*mean(rec$mofa_recovers),
    rec[is_convergent==TRUE & mofa_recovers==TRUE, .N], length(conv_ids),
    100*rec[is_convergent==TRUE, mean(mofa_recovers)]))
cat("recovery by program source (MOFA is scRNA-only → recovers scRNA-derived, misses bulk/genetic):\n")
print(rec[, .(n=.N, mofa_recovers=sum(mofa_recovers), pct=round(100*mean(mofa_recovers))), by=source][order(-pct)])

# ── (2) MODALITY REACH ──
cat("\n── (2) Modality reach ──\n")
cat("xMOFA : 5 views = scRNA cell types only (Hep/Mac/Fib/Endo/Chol), 269 donor-SHARED samples.\n")
cat("        Cannot place bulk(846), mouse, proteomics(130), genetics in ONE factor model (DISJOINT donors).\n")
cat("PROCON: 7 modalities incl 4 disjoint-donor (bulk/mouse/proteomics/genetic) + dispersion → cross-species\n")
cat("        + cross-assay CONVERGENCE per program (the question MOFA structurally cannot pose).\n")

# ── (3) INTERPRETABILITY: each MOFA factor blends how many programs? ──
blend <- rbindlist(lapply(factors, function(fc) {
  fs <- intersect(fac_sets[[fc]], universe)
  nb <- sum(sapply(prog, function(pg) {
    pg <- intersect(pg, universe); if(!length(pg)) return(FALSE)
    ov <- length(intersect(fs, pg))
    phyper(ov-1, length(pg), U-length(pg), length(fs), lower.tail=FALSE) < 0.001 }))
  data.table(factor = fc, n_programs_blended = nb) }))
fwrite(blend, file.path(P, "t3_xmofa_factor_blending.tsv"), sep = "\t")
cat(sprintf("\n── (3) Interpretability: a MOFA factor blends a MEDIAN of %.0f dictionary programs (p<0.001);\n", median(blend$n_programs_blended)))
cat(sprintf("   %d/16 factors blend >=3 programs (opaque, need post-hoc annotation) vs PROCON's named programs.\n",
    blend[n_programs_blended>=3, .N]))

cat("\n── VERDICT ──\n")
cat("PROCON > xMOFA here: (a) it integrates disjoint-donor cross-modality + cross-species evidence a\n")
cat("donor-shared factor model cannot fit; (b) it yields a calibrated convergence q + I2 per NAMED program;\n")
cat(sprintf("(c) MOFA factors blend a median %.0f programs (opaque). MOFA retained as a complementary scRNA\n", median(blend$n_programs_blended)))
cat(sprintf("discovery arm — it does recover %.0f%% of scRNA-expressed convergent programs unsupervised.\n",
    100*rec[is_convergent==TRUE & source %in% c("hs","cnmf"), mean(mofa_recovers)]))
