#!/usr/bin/env Rscript
# Supplementary panels — DEG-threshold sensitivity of the orthogonality result.
#
# Panel A: % of jointly-testable SuSiE genes that are NOT differentially expressed,
#          across |logFC| floors paired with padj<0.05, with the canonical TREAT
#          gate marked.
# Panel B: Fisher odds ratio for overlap against the jointly-testable background,
#          same sweep. Shows the enrichment/depletion direction and where it flips.
#
# Point: the headline is stable for any defensible effect-size floor and only
# breaks where "DEG" would mean a third of the transcriptome. Publishing the ladder
# pre-empts the "you chose a conservative threshold" objection.

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

stopifnot(nrow(treat) == 1, treat$n_deg == 1918, abs(treat$pct_non_de - 92.3937) < 0.01)

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

m <- conv[which.min(abs(n_deg - treat$n_deg))]
message(sprintf(paste0(
"DEG-threshold sensitivity of the orthogonality result. Across |log2FC| floors paired ",
"with padj<0.05 (black), the percentage of jointly-testable SuSiE genes that are not ",
"differentially expressed is stable at %.0f-%.0f%% for any floor >=0.30. The canonical ",
"TREAT interval-null gate (red diamond; %d genes, %.1f%%) sits within this range, and ",
"the closest conventional definition (%s; %d genes) gives %.1f%%. Overlap is not ",
"enriched at any defensible floor; TREAT's odds ratio (%.2f) is in fact HIGHER than ",
"every conventional variant, so the canonical gate is not the choice most favourable ",
"to the orthogonality claim. The shaded region (<0.15) is where the definition admits ",
"%s-%s genes, a third of the tested transcriptome, and is shown for completeness ",
"rather than as a defensible DEG set."),
min(conv[lfc_floor>=0.30]$pct_non_de), max(conv[lfc_floor>=0.30]$pct_non_de),
treat$n_deg, treat$pct_non_de, m$label, m$n_deg, m$pct_non_de, treat$fisher_or,
format(max(conv$n_deg), big.mark=","), format(nofl$n_deg, big.mark=",")))
