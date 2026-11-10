# DIG1 job 48b (exploratory, descriptive; follows jobs 46 and 47). Which Hallmark pathways rise or fall with
# the KRT8/SQSTM1 score inside hepatocytes? Per-gene statistic = mean within-donor excess (job 47 A2: depth,
# zonation and long-gene share removed, bin-centred on mean x detection x span tertile; job 46: depth only).
# limma cameraPR (inter-gene correlation 0.01, default), Hallmark GMT used in job 05. BH over 50 sets per
# cohort and adjustment. Genes are correlated within cells, so these p values are a ranking aid.
# Usage: Rscript 48b_hallmark_injury_axis.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
gl <- strsplit(readLines(file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt")), "\t")
HM <- setNames(lapply(gl, function(x) x[-(1:2)]), vapply(gl, `[`, "", 1))
a2 <- fread(file.path(RUN, "out_length_confound/genome_A2_cross_cohort.tsv.gz"))
f46 <- fread(file.path(RUN, "out_injury_coupling_genome/flex_gene_coupling.tsv.gz"))
s46 <- fread(file.path(RUN, "out_injury_coupling_genome/snrna_gene_coupling.tsv.gz"))
stats_list <- list(flex_A2 = setNames(a2$mean_excess_A2_flex, a2$gene), snrna_A2 = setNames(a2$mean_excess_A2_sn, a2$gene),
                   flex_depth = setNames(f46$mean_excess, f46$gene), snrna_depth = setNames(s46$mean_excess, s46$gene))
res <- rbindlist(lapply(names(stats_list), function(k) {
  st <- stats_list[[k]]; st <- st[is.finite(st)]
  idx <- lapply(HM, function(g) which(names(st) %in% g)); idx <- idx[lengths(idx) >= 10]
  r <- cameraPR(st, idx, sort = FALSE)
  data.table(analysis = k, set = rownames(r), n = r$NGenes, direction = r$Direction, p = r$PValue, bh = p.adjust(r$PValue, "BH"),
             mean_member = vapply(idx, function(i) mean(st[i]), 0), mean_all = mean(st))
}))
res[, signed := ifelse(direction == "Up", 1, -1) * -log10(p)]
fwrite(res, file.path(out, "hallmark_cameraPR.tsv"), sep = "\t")
w <- dcast(res, set ~ analysis, value.var = "signed"); fwrite(w, file.path(out, "hallmark_signed_log10p_wide.tsv"), sep = "\t")
agree <- data.table(pair = c("flex_A2 vs snrna_A2", "flex_depth vs snrna_depth"),
                    spearman = c(cor(w$flex_A2, w$snrna_A2, method = "spearman", use = "complete"), cor(w$flex_depth, w$snrna_depth, method = "spearman", use = "complete")))
fwrite(agree, file.path(out, "hallmark_cross_cohort_agreement.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(agree)
w2 <- merge(w, dcast(res, set ~ analysis, value.var = "bh"), by = "set", suffixes = c("", "_bh"))
print(w2[order(-(flex_A2 + snrna_A2))], digits = 3)
