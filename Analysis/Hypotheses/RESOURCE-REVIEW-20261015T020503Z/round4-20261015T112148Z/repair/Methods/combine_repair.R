# KEY MESSAGE: under GSE244832 donor relabeling the directional-replication counts vary widely, so the
# observed counts (hepatocyte 849 of 1,016; stellate 15 of 15; 4 stellate peaks significant in both cohorts)
# sit in the upper tail but are not extreme; the binomial q values that treated peaks as independent are withdrawn.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set.seed(20260923)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/Methods"
F <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/Methods/repair"
nul <- rbindlist(lapply(list.files(B, "^null_(stellate|hepatocyte)_[0-9]+_[0-9]+\\.tsv$", full.names = TRUE), fread))
setorder(nul, lineage, relabeling)
idn <- rbindlist(lapply(list.files(B, "^identity_.*\\.tsv$", full.names = TRUE), fread))
for (L in c("stellate", "hepatocyte")) {
  x <- nul[lineage == L]
  stopifnot(nrow(x) == 2002, identical(x$relabeling, 1:2002), sum(x$is_observed) == 1)
  i <- idn[lineage == L]
  stopifnot(i$max_abs_logFC_diff_vs_release < 1e-10,
            x[is_observed == TRUE, same_sign_disc] == i$saved_same_sign_disc,
            x[is_observed == TRUE, n_supported] == i$saved_n_supported)
}
ndisc <- setNames(as.list(idn$n_disc_gse281367), idn$lineage)
obs <- nul[is_observed == TRUE, .(lineage, obs_same = same_sign_disc, obs_supported = n_supported, obs_sig244 = n_sig_gse244832)]
s <- nul[obs, on = "lineage"][, .(
  n_disc_gse281367 = ndisc[[.BY$lineage]],
  observed_same_sign = obs_same[1], n_relabelings = .N,
  null_mean_same_sign = mean(same_sign_disc), null_sd_same_sign = sd(same_sign_disc),
  n_ge_observed_same_sign = sum(same_sign_disc >= obs_same[1]),
  observed_supported = obs_supported[1], null_mean_supported = mean(n_supported),
  n_ge_observed_supported = sum(n_supported >= obs_supported[1]),
  null_median_sig_gse244832 = median(n_sig_gse244832), observed_sig_gse244832 = obs_sig244[1]
), by = lineage]
s[, p_same_sign := n_ge_observed_same_sign / n_relabelings]
s[, binomial_sd_independent_peaks := sqrt(n_disc_gse281367 / 4)]
s[, p_supported := n_ge_observed_supported / n_relabelings]
# Release family (04_run_da.R:346): BH over three lineage rows, macrophage untestable at p = 1.
p3 <- c(s$p_same_sign, 1); names(p3) <- c(s$lineage, "macrophage")
q3 <- p.adjust(p3, "BH")
s[, q_BH_three_rows_release_family := q3[lineage]]
s[, q_BH_two_testable_rows_post_hoc := p.adjust(p_same_sign, "BH")]
s[, saved_binomial_q := c(stellate = 4.57763671875e-05, hepatocyte = 2.63941524908429e-110)[lineage]]
print(s)
fwrite(s, file.path(F, "R11-04_R11-08_donor_relabel_summary_repair.tsv"), sep = "\t")
fwrite(nul, file.path(F, "R11-04_R11-08_donor_relabel_null_repair.tsv.gz"), sep = "\t")

fmt <- function(v) formatC(v, format = "d", big.mark = ",")
int_breaks <- function(lim) { b <- pretty(lim); b[b == round(b)] }
panel <- function(v, obs_v, n_ge, title, discrete, sqrt_y = FALSE) {
  sub <- paste0("observed ", fmt(obs_v), "; ", fmt(n_ge), " of 2,002 relabelings\nat or above observed (p = ",
                formatC(n_ge / 2002, format = "f", digits = 3), ")")
  g <- ggplot()
  if (discrete) {
    g <- g + geom_col(data = data.table(value = v)[, .N, by = value], aes(value, N), width = 0.8, fill = "#9E9E9E")
  } else {
    g <- g + geom_histogram(data = data.table(value = v), aes(value), binwidth = 25, boundary = 0, fill = "#9E9E9E", colour = "#9E9E9E", linewidth = 0.05)
  }
  ylab <- if (sqrt_y) "Relabelings (of 2,002), square root scale" else "Relabelings (of 2,002)"
  yb <- if (sqrt_y) c(0, 10, 50, 200, 500, 1000, 2000) else waiver()
  g + geom_vline(xintercept = obs_v, linewidth = 0.3, colour = "black") +
    scale_x_continuous(breaks = int_breaks, labels = fmt) +
    scale_y_continuous(trans = if (sqrt_y) "sqrt" else "identity", breaks = yb, labels = fmt,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, subtitle = sub, x = "Count under GSE244832 donor relabeling", y = ylab) +
    theme_classic(base_size = 6) +
    theme(text = element_text(size = 6, colour = "black"), axis.text = element_text(size = 6, colour = "black"),
          plot.title = element_text(size = 6, face = "plain"), plot.subtitle = element_text(size = 6),
          axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
}
sh <- s[lineage == "hepatocyte"]; ss <- s[lineage == "stellate"]
p1 <- panel(nul[lineage == "hepatocyte", same_sign_disc], sh$observed_same_sign, sh$n_ge_observed_same_sign,
            "Hepatocyte: same sign among 1,016\nGSE281367 discovery peaks", FALSE)
p2 <- panel(nul[lineage == "stellate", same_sign_disc], ss$observed_same_sign, ss$n_ge_observed_same_sign,
            "Stellate: same sign among 15\nGSE281367 discovery peaks", TRUE)
p3 <- panel(nul[lineage == "stellate", n_supported], ss$observed_supported, ss$n_ge_observed_supported,
            "Stellate: peaks with q < 0.05 in both\ncohorts and the same sign", TRUE, sqrt_y = TRUE)
ggsave(file.path(F, "R11-04_hepatocyte_same_sign_relabel_null_repair.pdf"), p1, width = 2.1, height = 2, useDingbats = FALSE)
ggsave(file.path(F, "R11-04_stellate_same_sign_relabel_null_repair.pdf"), p2, width = 2.1, height = 2, useDingbats = FALSE)
ggsave(file.path(F, "R11-08_stellate_supported_relabel_null_repair.pdf"), p3, width = 2.1, height = 2, useDingbats = FALSE)
message("Caption: Exact GSE244832 donor-relabeling nulls for ATAC v3 directional replication. Each of the 2,002 ",
        "assignments of 9 MASH labels among 14 donors refits the release model (condition-blind filterByExpr, TMM, ",
        "voomWithQualityWeights, robust eBayes); GSE281367 fits are fixed. Gray, null counts; line, observed. ",
        "Exact p = relabelings at or above observed / 2,002, identity included.")
writeLines(capture.output(sessionInfo()), file.path(F, "R11-04_R11-08_sessionInfo_repair.txt"))
