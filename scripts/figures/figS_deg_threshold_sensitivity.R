#!/usr/bin/env Rscript
# Supplementary panels — DEG-threshold sensitivity of the orthogonality result.
#
# Panel A: % of jointly-testable SuSiE genes that are NOT differentially expressed,
#          across |logFC| floors paired with padj<0.05, with the canonical TREAT
#          gate marked.
# Panel B: Fisher odds ratio for overlap against the jointly-testable background,
#          same sweep. Shows the enrichment/depletion direction across the ladder.
#
# Point: the headline is stable for any defensible effect-size floor. Publishing
# the ladder pre-empts the "you chose a conservative threshold" objection.
#
# 2026-08-13: Panel B changed materially when the Fisher background was corrected
# from all 27,638 bulk-tested genes to the 14,931 jointly testable ones (see
# figS_deg_threshold_sensitivity_data.py). The two significant ENRICHMENTS at the
# permissive end (no floor, 1.612, p=7.8e-07; >0.10, 1.296, p=0.0075) were
# artifacts of the inflated denominator and are gone (1.104 and 1.082, both n.s.).
# The corrected ladder decreases monotonically across all eleven floors and is
# never significantly above 1. Panel A is unchanged -- it is computed on the 447
# genetic genes alone, which the background correction does not touch.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

DIR <- file.path(BASE, "figures/supplementary/figS_deg_threshold_sensitivity")
s <- fread(file.path(DIR, "figS_deg_threshold_sensitivity_source.tsv"))

conv  <- s[gate == "conventional" & lfc_floor > 0]
nofl  <- s[gate == "conventional" & lfc_floor == 0]
treat <- s[gate == "treat"]

stopifnot(nrow(treat) == 1,
          treat$n_deg_all == 1918, treat$n_deg_joint == 1261,
          abs(treat$pct_non_de - 92.3937) < 0.01,
          abs(treat$fisher_or - 0.8895) < 1e-3,   # matches orthogonality_audit.tsv
          all(s$n_background == 14931), all(s$n_genetic == 447))

GREY <- "#9E9E9E"; ACC <- "#C9265E"

# --- Panel A: % non-DE -------------------------------------------------------
pA <- ggplot(conv, aes(lfc_floor, pct_non_de)) +
  annotate("rect", xmin = -Inf, xmax = 0.15, ymin = -Inf, ymax = Inf,
           fill = GREY, alpha = 0.18) +
  geom_line(linewidth = 0.4, colour = "black") +
  geom_point(size = 0.9, colour = "black") +
  geom_point(data = treat, aes(lfc_floor, pct_non_de), size = 2.0,
             shape = 18, colour = ACC) +
  geom_hline(yintercept = treat$pct_non_de, linewidth = 0.3,
             linetype = "dashed", colour = ACC) +
  scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "|log2FC| floor (paired with padj < 0.05)",
       y = "SuSiE genes not differentially expressed (%)") +
  theme_masld(base_size = 6) +
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"))

ggsave(file.path(DIR, "figS_deg_threshold_pct_nonoverlap.pdf"), pA,
       width = 3.4, height = 2.3, useDingbats = FALSE)

# --- Panel B: Fisher OR ------------------------------------------------------
pB <- ggplot(conv, aes(lfc_floor, fisher_or)) +
  annotate("rect", xmin = -Inf, xmax = 0.15, ymin = -Inf, ymax = Inf,
           fill = GREY, alpha = 0.18) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = GREY) +
  geom_line(linewidth = 0.4, colour = "black") +
  geom_point(size = 0.9, colour = "black") +
  geom_point(data = treat, aes(lfc_floor, fisher_or), size = 2.0,
             shape = 18, colour = ACC) +
  scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "|log2FC| floor (paired with padj < 0.05)",
       y = "Overlap odds ratio vs jointly-testable background") +
  theme_masld(base_size = 6) +
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"))

ggsave(file.path(DIR, "figS_deg_threshold_enrichment_or.pdf"), pB,
       width = 3.4, height = 2.3, useDingbats = FALSE)

m    <- conv[which.min(abs(n_deg_joint - treat$n_deg_joint))]
# Max/min over EVERY conventional row including the unplotted no-floor point,
# so the caption's range and its extreme value cannot disagree.
allc <- s[gate == "conventional"]
hi   <- allc[which.max(fisher_or)]
lo   <- allc[which.min(fisher_or)]
sigd <- allc[fisher_p < 0.05][order(lfc_floor)]

message(sprintf(paste0(
"DEG-threshold sensitivity of the orthogonality result. All tests use the %s-gene ",
"jointly testable background -- genes assayed in both the bulk contrast and ",
"colocalization -- and the %d colocalized genes within it. Across |log2FC| floors ",
"paired with padj<0.05 (black), the percentage of colocalized genes that are not ",
"differentially expressed runs %.1f-%.1f%% for any floor >=0.30 (panel A). The ",
"interval-null comparator gate (red diamond; %s genes genome-wide, %s within the ",
"background, %.1f%%) sits within this range, and the closest conventional definition ",
"(%s; %s genes) gives %.1f%%. The overlap odds ratio decreases monotonically across ",
"the whole ladder, from %.2f with no floor at all (not plotted) to %.2f at ",
"|log2FC|>1.00, and is not significantly above 1 anywhere -- the largest ",
"value, %.2f, carries p = %.2f. Depletion reaches nominal significance only at floors ",
"%s (p = %s), and no floor survives correction across the eleven tested. The direction ",
"is therefore stable across every threshold, while the magnitude tracks effect-size ",
"amplitude: colocalized genes are progressively less represented among genes with ",
"larger disease-state responses. The shaded region (<0.15) admits %s-%s of the %s ",
"background genes and is shown for completeness rather than as a defensible DEG set. ",
"The decline in panel B is an expression artifact and not a property of inherited ",
"risk: colocalized genes are higher-expressed than the background (AveExpr 3.70 vs ",
"2.59, SMD 0.36), and when each is matched to an expression-matched control the ",
"ladder flattens to unity at every floor (matched odds ratios 0.67-1.17, all McNemar ",
"p >= 0.27; 23 colocalized genes strongly differentially expressed against 24 matched ",
"controls at the 0.50 floor). Colocalized genes respond as much as comparable genes ",
"and must not be described as buffered."),
format(unique(s$n_background), big.mark=","), unique(s$n_genetic),
min(conv[lfc_floor>=0.30]$pct_non_de), max(conv[lfc_floor>=0.30]$pct_non_de),
format(treat$n_deg_all, big.mark=","), format(treat$n_deg_joint, big.mark=","),
treat$pct_non_de,
m$label, format(m$n_deg_joint, big.mark=","), m$pct_non_de,
nofl$fisher_or, lo$fisher_or, hi$fisher_or, hi$fisher_p,
paste(sprintf("%.2f", sigd$lfc_floor), collapse=" and "),
paste(sprintf("%.3f", sigd$fisher_p), collapse=" and "),
format(nofl$n_deg_joint, big.mark=","),
format(min(conv[lfc_floor<0.15]$n_deg_joint), big.mark=","),
format(unique(s$n_background), big.mark=",")))
