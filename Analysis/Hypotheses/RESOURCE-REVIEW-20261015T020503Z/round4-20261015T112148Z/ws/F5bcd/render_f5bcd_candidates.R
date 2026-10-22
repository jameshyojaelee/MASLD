#!/usr/bin/env Rscript
# Candidate panels for review round 3, workstream F5bcd. Writes only to the candidate dir.
# 5D-gene:  KEY MESSAGE: counted per gene, most colocalized genes place under 5% of their
#           colocalization posterior in any lineage's open chromatin; both-source overlap is a minority.
# 5D-bg:    KEY MESSAGE: averaged over liver-enzyme genes, posterior mass in peaks is about 2.5-3x
#           the locus-matched expectation in every lineage; the 22 direct-trait genes are too few to say.
# 5B-bar:   KEY MESSAGE: protein-RNA sign agreement is near complete where RNA is significant and
#           near chance where it is not.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
grDevices::pdf.options(useDingbats = FALSE)
set.seed(20260923)
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd"
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
th <- theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(text = element_text(size = 6, face = "plain", colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        strip.text = element_text(size = 6, hjust = 0), strip.background = element_blank(),
        legend.text = element_text(size = 6), legend.title = element_blank(),
        legend.key.size = grid::unit(0.18, "cm"), legend.position = "bottom",
        panel.spacing.x = grid::unit(9, "pt"), plot.margin = margin(1, 5, 1, 1))
lin_lab <- c(hepatocyte = "Hepatocyte", stellate = "Stellate", macrophage = "Macrophage",
             cholangiocyte = "Cholangiocyte", t_nk = "T/NK")
trait_lab <- c(direct_MASLD = "MASLD/liver fat", liver_enzyme = "Liver enzyme")

# ---------- 5D gene-level, 4 display states ----------
r <- fread(file.path(OUT, "fig5d_pair_rows_with_background.tsv"))
ORDER <- c("replicated_accessible", "source_dependent", "partial_ge_floor", "trace_lt_floor", "indeterminate")
g <- r[, .(best = min(rank_f05)), by = .(trait_class, lineage, ensembl, gene)]
g[, state := ORDER[best + 1L]]
g[, disp := fifelse(state %in% c("trace_lt_floor", "indeterminate"), "under_5pct", state)]
disp_lab <- c(replicated_accessible = "Both sources ≥50%", source_dependent = "One source ≥50%",
              partial_ge_floor = "5% to <50%", under_5pct = "Under 5% (incl. none)")
cnt <- g[, .N, by = .(trait_class, lineage, disp)]
grid <- CJ(trait_class = names(trait_lab), lineage = names(lin_lab), disp = names(disp_lab))
cnt <- merge(grid, cnt, all.x = TRUE); cnt[is.na(N), N := 0L]
den <- unique(g[, .(trait_class, ensembl)])[, .(denominator = .N), by = trait_class]
cnt <- merge(cnt, den, by = "trait_class")
cnt[, `:=`(proportion = N / denominator, unit = "gene (best state over its gene-study pairs)",
           partial_floor = 0.05)]
fwrite(cnt, file.path(OUT, "fig5d_gene_level_display_counts.tsv"), sep = "\t")
cnt[, lineage_f := factor(lin_lab[lineage], levels = rev(lin_lab))]
cnt[, disp_f := factor(disp_lab[disp], levels = disp_lab)]
cnt[, facet := factor(sprintf("%s (n = %d genes)", trait_lab[trait_class], denominator),
                      levels = sprintf("%s (n = %d genes)", trait_lab, den[match(names(trait_lab), trait_class), denominator]))]
fills <- c("Both sources ≥50%" = "#00695C", "One source ≥50%" = "#7B1FA2",
           "5% to <50%" = "#8EC7E2", "Under 5% (incl. none)" = "#D9D9D9")
p <- ggplot(cnt, aes(proportion, lineage_f, fill = disp_f)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.12) +
  geom_text(aes(label = ifelse(proportion >= 0.065 & N > 0, N, "")),
            position = position_stack(vjust = 0.5), size = 6 / .pt, family = "Helvetica") +
  facet_wrap(~facet, nrow = 1) +
  scale_fill_manual(values = fills, drop = FALSE) +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "50", "100"), expand = expansion(mult = c(0, .04))) +
  labs(x = "Genes (%)", y = NULL) + th +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(file.path(OUT, "fig5d_candidate_gene_level_mass_floor.pdf"), p, width = 2.96, height = 1.88,
       device = grDevices::cairo_pdf)

# ---------- 5D background: observed vs locus-matched expected ----------
b <- fread(file.path(OUT, "fig5d_locus_background_summary.tsv"))
bl <- rbind(
  b[, .(trait_class, lineage, n_genes, what = "Observed", est = mean_observed, lo = obs_lo, hi = obs_hi)],
  b[, .(trait_class, lineage, n_genes, what = "Locus-matched expectation", est = mean_expected, lo = exp_lo, hi = exp_hi)])
bl[, lineage_f := factor(lin_lab[lineage], levels = rev(lin_lab))]
bl[, facet := factor(sprintf("%s (n = %d genes)", trait_lab[trait_class], n_genes),
                     levels = unique(sprintf("%s (n = %d genes)", trait_lab[b$trait_class], b$n_genes)))]
bl[, what := factor(what, levels = c("Observed", "Locus-matched expectation"))]
pb <- ggplot(bl, aes(100 * est, lineage_f, colour = what)) +
  geom_errorbar(aes(xmin = 100 * lo, xmax = 100 * hi), width = 0, orientation = "y", linewidth = 0.3,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 0.9, position = position_dodge(width = 0.5)) +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = c("Observed" = "#00695C", "Locus-matched expectation" = "#9E9E9E")) +
  labs(x = "Colocalization posterior in peaks, gene mean (%)", y = NULL) + th
ggsave(file.path(OUT, "fig5d_candidate_locus_matched_background.pdf"), pb, width = 2.96, height = 1.88,
       device = grDevices::cairo_pdf)

# ---------- 5B bar split by RNA call status ----------
d <- fread(file.path(BASE, "figures/main/fig5_molecular_context/panels/data/fig5b_protein_triage_genes.tsv"))
s <- d[protein_padj < 0.05]
s[, rna := fifelse(bulk_padj < 0.05, "RNA BH padj <0.05", "RNA BH padj ≥0.05")]
s[, agree := fifelse(sign(bulk_logFC) == sign(protein_logFC), "Same sign", "Opposite sign")]
bc <- s[, .N, by = .(rna, agree)]
bc[, n_rna := sum(N), by = rna][, proportion := N / n_rna]
bc[, binom_p_vs_half := sapply(rna, function(x) binom.test(s[rna == x & agree == "Same sign", .N], s[rna == x, .N])$p.value)]
fwrite(bc, file.path(OUT, "fig5b_candidate_sign_agreement_by_rna_status.tsv"), sep = "\t")
bc[, rna_f := factor(sprintf("%s (n = %d)", rna, n_rna))]
bc[, agree := factor(agree, levels = c("Opposite sign", "Same sign"))]
pbb <- ggplot(bc, aes(proportion, rna_f, fill = agree)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.12) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 6 / .pt, family = "Helvetica") +
  scale_fill_manual(values = c("Same sign" = "#00695C", "Opposite sign" = "#9E9E9E")) +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "50", "100"), expand = expansion(mult = c(0, .04))) +
  labs(x = "Protein-significant genes (%)", y = NULL) + th +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
ggsave(file.path(OUT, "fig5b_candidate_sign_agreement_by_rna_status.pdf"), pbb, width = 2.6, height = 1.0,
       device = grDevices::cairo_pdf)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("CAPTION (5D candidate): genes, best state over pairs; partial floor 0.05 of normalized SNP.PP.H4 posterior.")
