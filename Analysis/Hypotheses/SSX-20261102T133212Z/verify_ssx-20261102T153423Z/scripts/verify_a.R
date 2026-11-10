# Independent recompute of SSX-A primary statistics from raw ASE counts.
suppressPackageStartupMessages(library(data.table))
wd <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z"
out <- file.path(wd, "verify_ssx-20261102T153423Z")
ss <- file.path(wd, "../SS-20261025T135237Z")
pg <- fread(file.path(wd, "out_a-20261102T133643Z/primary_genes.tsv"))
units <- fread(file.path(ss, "out_s2_qc/units.tsv"))
inc <- units[include == TRUE, sample_alias]
cat("included units:", length(inc), "\n")
ase <- fread(file.path(ss, "out_s3fgh-20261026T233706Z/counts/ase_counts.tsv.gz"))
ase <- ase[sample_alias %in% inc]
cat("dup sample x site rows:", sum(duplicated(ase[, .(sample_alias, chrom, pos)])), "\n")
x <- merge(ase, unique(pg[, .(chrom, pos, gene = id)]), by = c("chrom", "pos"))
x[, dp := nref + nalt]
het <- x[dp >= 10 & nalt / dp > 0.1 & nalt / dp < 0.9]
res <- het[, .(n_het_v = .N, med_dp_v = as.numeric(median(dp)), pooled_f_v = sum(nalt) / sum(dp)), by = .(chrom, pos, gene)]
res <- merge(pg[, .(chrom, pos, gene = id, n_het, med_dp, pooled_f, null_med_logit, slope, shift, agree)], res,
             by = c("chrom", "pos", "gene"), all.x = TRUE)
res[, shift_v := qlogis(pooled_f_v) - null_med_logit]
res[, agree_v := sign(shift_v) == sign(slope)]
res[, rel_f := abs(pooled_f_v - pooled_f) / abs(pooled_f)]
res[, rel_shift := abs(shift_v - shift) / abs(shift)]
fwrite(res, file.path(out, "a_per_gene.tsv"), sep = "\t")
cat("genes:", nrow(res), "; n_het mismatches:", sum(res$n_het != res$n_het_v, na.rm = TRUE),
    "; med_dp mismatches:", sum(res$med_dp != res$med_dp_v, na.rm = TRUE), "\n")
cat("max rel diff pooled_f:", max(res$rel_f), "; max rel diff shift:", max(res$rel_shift), "\n")
k <- sum(res$agree_v); n <- nrow(res)
p <- binom.test(k, n, 0.5, alternative = "greater")$p.value
cat("share:", k, "of", n, "=", k / n, "; exact one-sided binomial p:", format(p, digits = 6), "\n")
cat("reported agree column count:", sum(res$agree), "\n")
print(res[1:6, .(gene, n_het, n_het_v, pooled_f, pooled_f_v, shift, shift_v)])
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_a.txt"))
