suppressPackageStartupMessages(library(data.table))
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/Methods"
F <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/Methods"
set.seed(20260923)
st <- fread(file.path(B, "7130/null_full_stellate.tsv")); st[, observed := 15L]
hp <- rbindlist(lapply(list.files(file.path(B, "7131"), "^null_chunk_hepatocyte", full.names = TRUE), fread))
setorder(hp, relabeling)
stopifnot(nrow(st) == 2002, nrow(hp) == 2002, all(hp$relabeling == 1:2002), sum(hp$is_observed) == 1, sum(st$is_observed) == 1)
nul <- rbind(st, hp)
s <- nul[, .(n_disc = c(stellate = 15L, hepatocyte = 1016L)[lineage[1]], observed = observed[1],
             identity_value = same_sign[is_observed], n_relabelings = .N,
             null_median = median(same_sign), null_q95 = quantile(same_sign, 0.95), null_max = max(same_sign),
             n_ge_observed = sum(same_sign >= observed[1])), by = lineage]
s[, p_exact_donor_relabel := n_ge_observed / n_relabelings]
s[, binomial_p_peaks_independent := mapply(function(k, n) binom.test(k, n, 0.5, alternative = "greater")$p.value, observed, n_disc)]
s[, BH_q_donor_relabel_2_lineages := p.adjust(p_exact_donor_relabel, "BH")]
s[, saved_binomial_BH_q := c(stellate = 4.58e-5, hepatocyte = 2.64e-110)[lineage]]
print(s)
fwrite(s, file.path(F, "R11-04_atac_directional_donor_relabel_summary.tsv"), sep = "\t")
fwrite(nul, file.path(F, "R11-04_atac_directional_donor_relabel_null.tsv.gz"), sep = "\t")
library(ggplot2)
nul[, lineage_lab := factor(paste0(tools::toTitleCase(lineage), " (", c(stellate = 15, hepatocyte = 1016)[lineage], " GSE281367 discovery peaks)"))]
obs <- unique(nul[, .(lineage_lab, observed)])
p <- ggplot(nul, aes(same_sign)) +
  geom_histogram(bins = 40, fill = "#9E9E9E", colour = NA) +
  geom_vline(data = obs, aes(xintercept = observed), linewidth = 0.4, colour = "black") +
  facet_wrap(~lineage_lab, scales = "free", nrow = 1) +
  labs(x = "Discovery peaks with the same sign in GSE244832", y = "Donor relabelings (of 2,002)") +
  theme_classic(base_size = 6) +
  theme(text = element_text(size = 6, colour = "black"), axis.text = element_text(colour = "black"),
        strip.background = element_blank(), strip.text = element_text(size = 6))
ggsave(file.path(F, "R11-04_atac_directional_donor_relabel_null.pdf"), p, width = 4.2, height = 1.7, useDingbats = FALSE)
writeLines(capture.output(sessionInfo()), file.path(F, "R11-04_sessionInfo.txt"))
