# DIG1 job 61 (exploratory, post hoc). Do the CCL20+ inflamed hepatocytes carry a p53 / cell-cycle-arrest
# program that is separate from their NF-kB response?
# Disclosure: before this script, gene-level excess was read for 34 cell-cycle, senescence and death genes
# (CDKN1A, MDM2, TNFRSF10B, SERPINE1 up in all four cohorts; BBC3, CDKN2A up in three; PCNA, CASP3 flat).
# Input: genome-wide gene excess of jobs 53 (Flex GSE298719, snRNA GSE256398; MASLD donors) and 56 (atlas
# GSE202379, GSE244832): anchor+ hepatocytes vs depth(+injury)-matched anchor- hepatocytes of the same donor,
# minus 20 detection-matched random anchors, mean over eligible donors. Unit of the statistic is the gene;
# donors were the unit when each gene's excess was formed.
# Test: limma cameraPR on the gene excess (competitive; inter-gene correlation 0.01 primary, 0.05 and a
# rank version as sensitivities). Sets from MSigDB v2025.1 symbols (sha256 recorded in the sbatch log).
#   Family H: the 50 Hallmarks, BH within cohort x anchor x mode.
#   Family S (named here): FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS, HALLMARK_P53_PATHWAY, FRIDMAN_SENESCENCE_UP,
#     FRIDMAN_SENESCENCE_DN, SAUL_SEN_MAYO, REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP,
#     KAMMINGA_SENESCENCE, HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT; Holm within cohort x anchor x mode.
#   Each family S set is also tested "clean": genes removed if in HALLMARK_TNFA_SIGNALING_VIA_NFKB,
#   HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_INTERFERON_ALPHA_RESPONSE,
#   the job 56 NFKB / IFN / CORE / BILIARY sets, or the anchors.
# Main readout (named here): CCL20 anchor, depth mode, clean FISCHER p53 set, direction and Holm p in each of
# four cohorts. Contrast: the same for SOX9 and CXCL10 anchors. Depth+injury mode is a sensitivity.
# Usage: Rscript 61_state_p53_senescence_sets.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1)
out <- args[1]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
GMT <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt")
g53 <- fread(file.path(RUN, "out_chemokine_hepatocyte_state/gene_lfc_matched.tsv.gz"))[donors == "masld"][, donors := NULL]
g56 <- fread(file.path(RUN, "out_atlas_state_replication/gene_lfc_matched.tsv.gz"))
G <- rbind(g53, g56, use.names = TRUE)
G <- G[is.finite(mean_excess)]
gl <- strsplit(readLines(GMT), "\t")
all_sets <- setNames(lapply(gl, function(x) unique(x[-(1:2)])), sapply(gl, `[`, 1))
HALL <- grep("^HALLMARK_", names(all_sets), value = TRUE); stopifnot(length(HALL) == 50)
SEN <- c("FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS", "HALLMARK_P53_PATHWAY", "FRIDMAN_SENESCENCE_UP", "FRIDMAN_SENESCENCE_DN",
         "SAUL_SEN_MAYO", "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP", "KAMMINGA_SENESCENCE",
         "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
stopifnot(all(SEN %in% names(all_sets)))
st56 <- fread(file.path(RUN, "out_atlas_state_replication/gene_sets.tsv"))
EXCL <- unique(c(unlist(all_sets[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                                    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")]),
                 unlist(strsplit(st56$genes, ",")), "CCL20", "CXCL10", "SOX9"))
sets_S <- all_sets[SEN]
sets_C <- setNames(lapply(sets_S, setdiff, EXCL), paste0(SEN, "__clean"))
writeLines(sprintf("%s\t%d\t%d", SEN, lengths(sets_S), lengths(sets_C)), file.path(out, "family_S_sizes_before_matching.tsv"))
run <- function(stat, sets, fam) {
  idx <- ids2indices(sets, names(stat)); idx <- idx[lengths(idx) >= 10]
  a <- cameraPR(stat, idx, inter.gene.cor = 0.01, sort = FALSE)
  b <- cameraPR(stat, idx, inter.gene.cor = 0.05, sort = FALSE)
  r <- cameraPR(stat, idx, use.ranks = TRUE, sort = FALSE)
  data.table(set = names(idx), family = fam, n_genes = lengths(idx),
             mean_excess_set = sapply(idx, function(i) mean(stat[i])), mean_excess_rest = sapply(idx, function(i) mean(stat[-i])),
             frac_pos = sapply(idx, function(i) mean(stat[i] > 0)),
             direction = a$Direction, p = a$PValue, p_cor05 = b$PValue, dir_rank = r$Direction, p_rank = r$PValue)
}
res <- list(); lead <- list()
for (k in unique(G[, .(cohort, anchor, mode)])[, paste(cohort, anchor, mode, sep = "|")]) {
  kk <- strsplit(k, "\\|")[[1]]
  d <- G[cohort == kk[1] & anchor == kk[2] & mode == kk[3]]
  d <- d[!duplicated(gene)]
  stat <- setNames(d$mean_excess, d$gene)
  h <- run(stat, all_sets[HALL], "H"); h[, q_bh := p.adjust(p, "BH")]
  s <- run(stat, sets_S, "S"); s[, p_holm := p.adjust(p, "holm")]
  cl <- run(stat, sets_C, "S_clean"); cl[, p_holm := p.adjust(p, "holm")]
  x <- rbind(h, s, cl, fill = TRUE); x[, `:=`(cohort = kk[1], anchor = kk[2], mode = kk[3], n_genes_tested = length(stat))]
  res[[k]] <- x
  for (nm in c("FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS", "HALLMARK_P53_PATHWAY")) {
    gg <- intersect(setdiff(all_sets[[nm]], EXCL), names(stat))
    lead[[paste(k, nm)]] <- data.table(cohort = kk[1], anchor = kk[2], mode = kk[3], set = paste0(nm, "__clean"),
                                       gene = gg, excess = stat[gg], n_pos = d$n_pos_excess[match(gg, d$gene)],
                                       n_donors = d$n_donors[match(gg, d$gene)])
  }
}
R <- rbindlist(res, fill = TRUE); L <- rbindlist(lead)
setcolorder(R, c("cohort", "anchor", "mode", "family", "set"))
fwrite(R, file.path(out, "camera_state_sets.tsv"), sep = "\t")
fwrite(L, file.path(out, "p53_clean_gene_excess.tsv"), sep = "\t")
A <- R[family == "H" & mode == "depth", .(n_up_q05 = sum(direction == "Up" & q_bh < 0.05), n_down_q05 = sum(direction == "Down" & q_bh < 0.05),
                                          n_cohorts = .N), by = .(anchor, set)]
fwrite(A[order(anchor, -n_up_q05, -n_down_q05)], file.path(out, "hallmark_cross_cohort_agreement_depth.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 250)
print(R[family %in% c("S", "S_clean"), .(cohort, anchor, mode, set = sub("_META_ANALYSIS|REACTOME_", "", set), n_genes,
                                         m = round(mean_excess_set, 2), rest = round(mean_excess_rest, 2), direction,
                                         p = signif(p, 2), p_holm = signif(p_holm, 2), p_cor05 = signif(p_cor05, 2), p_rank = signif(p_rank, 2))][order(anchor, mode, set, cohort)])
print(A[order(anchor, -n_up_q05, -n_down_q05)][n_up_q05 + n_down_q05 >= 3])
