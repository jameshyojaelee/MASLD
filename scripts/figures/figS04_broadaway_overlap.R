#!/usr/bin/env Rscript
##############################################################################
# figS04_broadaway_overlap.R
#
# Supplementary coloc panels - our MASLD SuSiE/ABF-COLOC effector genes vs the
# 747 colocalized eGenes of Broadaway et al. 2024 (AJHG, PMID 39173627), the
# source liver-eQTL meta-analysis (N=1,183) whose panel we reuse.
#
# KEY MESSAGE: Despite near-identical counts (751 vs 747), the two sets share
# only ~36% of genes. Their 747 is a broad *cardiometabolic* effector list
# (lipid-dominated); ours is *MASLD-anchored*. Two-thirds of the genes novel to
# us are reachable only through data Broadaway's design lacked - EUR AST,
# non-European ancestry (EAS/AFR/SAS), and dedicated NAFLD/PDFF/NASH GWAS.
#
# Input : GWAS/finemapping/results/susie_coloc/broadaway_coloc_comparison.csv
#         (per-gene: in_broadaway_747, novel_source, our_best_pp4, ...)
# Output (PDF only, flat in FIGS04_DIR):
#   broadaway_overlap_venn.pdf     area-proportional 2-set Venn
#   broadaway_novel_source.pdf     decomposition of our 751 genes
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

COMP <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/broadaway_coloc_comparison.csv")
stopifnot(file.exists(COMP))
d <- fread(COMP)

# ---------------------------------------------------------------------------
# Counts (derived, not hardcoded)
# ---------------------------------------------------------------------------
n_ours       <- nrow(d)
n_shared     <- d[, sum(in_broadaway_747)]
n_ours_only  <- n_ours - n_shared
n_their      <- 747L                       # Broadaway 747 effector eGenes (validated extraction)
n_their_only <- n_their - n_shared

src <- d[, .N, by = novel_source]
get <- function(k) { v <- src[novel_source == k, N]; if (length(v)) v else 0L }
n_newdata <- get("new_data(AST/nonEUR/MASLD-GWAS)")
n_method  <- get("method_diff(EUR_ALT_GGT)")

cat(sprintf("ours=%d their=%d shared=%d ours_only=%d their_only=%d | new_data=%d method=%d\n",
            n_ours, n_their, n_shared, n_ours_only, n_their_only, n_newdata, n_method))
stopifnot(n_shared + n_newdata + n_method == n_ours)

# colours: Broadaway = neutral gray (external reference); ours = signature magenta
COL_OURS  <- masld_colors$masld    # #C9265E deep magenta
COL_THEIR <- "#9E9E9E"             # Liang neutral gray (external reference)
COL_NEW   <- masld_colors$masld    # magenta - genuinely new data
COL_METH  <- masld_colors$masl     # #F4A674 warm peach - method/threshold
COL_SHARE <- "#9E9E9E"             # gray - shared with prior work

# ===========================================================================
# Panel A - area-proportional 2-set Venn (equal radii; 751 ~ 747)
# ===========================================================================
r <- 1
frac <- mean(c(n_shared / n_ours, n_shared / n_their))   # lens area / circle area
lens_area <- function(d) 2 * r^2 * acos(d / (2 * r)) - (d / 2) * sqrt(4 * r^2 - d^2)
dd <- uniroot(function(x) lens_area(x) - frac * pi * r^2,
              lower = 1e-4, upper = 2 * r - 1e-4)$root

circ <- data.table(
  x0   = c(-dd/2, dd/2),
  set  = c("ours", "their"),
  fill = c(COL_OURS, COL_THEIR)
)
lab <- data.table(
  x   = c(-dd/2 - r*0.55, 0, dd/2 + r*0.55),
  y   = c(0, 0, 0),
  n   = c(n_ours_only, n_shared, n_their_only),
  txt = c(sprintf("%d\nnovel to us", n_ours_only),
          sprintf("%d\nshared", n_shared),
          sprintf("%d\nBroadaway\nonly", n_their_only))
)
set_lab <- data.table(
  x = c(-dd/2, dd/2), y = c(r*1.18, r*1.18),
  txt = c(sprintf("Ours - %d\nMASLD 23-GWAS", n_ours),
          sprintf("Broadaway - %d\ncardiometabolic", n_their))
)

pA <- ggplot() +
  geom_circle(data = circ, aes(x0 = x0, y0 = 0, r = r, fill = set),
              colour = NA, alpha = 0.55, show.legend = FALSE) +
  scale_fill_manual(values = c(ours = COL_OURS, their = COL_THEIR)) +
  geom_text(data = lab, aes(x, y, label = txt),
            size = PUB_GEOM_TEXT + 0.4, lineheight = 0.9, fontface = "plain",
            colour = "black") +
  geom_text(data = set_lab, aes(x, y, label = txt),
            size = PUB_GEOM_TEXT, lineheight = 0.9, colour = "black") +
  coord_fixed(xlim = c(-dd/2 - r*1.2, dd/2 + r*1.2),
              ylim = c(-r*1.25, r*1.45)) +
  theme_void(base_family = "Helvetica") +
  theme(plot.margin   = margin(4, 4, 4, 4))

message(sprintf("[caption] MASLD COLOC effector genes vs Broadaway 747: overlap %d genes (%.0f%% of ours, %.0f%% of theirs); PP.H4 > 0.5",
                n_shared, 100*n_shared/n_ours, 100*n_shared/n_their))

ggsave(file.path(FIGS04_DIR, "broadaway_overlap_venn.pdf"),
       pA, width = 75, height = 70, units = "mm", useDingbats = FALSE)

# ===========================================================================
# Panel B - decomposition of our 751 genes (mutually exclusive, sums to 751)
# ===========================================================================
bar <- data.table(
  cat = factor(c("Shared with\nBroadaway 747",
                 "Novel - data they lacked\n(AST / non-EUR / MASLD GWAS)",
                 "Novel - method / threshold\n(shared EUR ALT·GGT)"),
               levels = c("Novel - method / threshold\n(shared EUR ALT·GGT)",
                          "Novel - data they lacked\n(AST / non-EUR / MASLD GWAS)",
                          "Shared with\nBroadaway 747")),
  n   = c(n_shared, n_newdata, n_method),
  col = c(COL_SHARE, COL_NEW, COL_METH)
)
bar[, pct := 100 * n / sum(n)]

pB <- ggplot(bar, aes(cat, n, fill = cat)) +
  geom_col(width = 0.72, colour = NA) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct)),
            hjust = -0.08, size = PUB_GEOM_TEXT + 0.3, colour = "black") +
  scale_fill_manual(values = setNames(bar$col, bar$cat), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "eGenes") +
  theme_masld() + theme_pub() +
  theme(axis.text.y    = element_text(size = PUB_AXIS_TEXT, lineheight = 0.85))

message("[caption] What distinguishes our 751 COLOC genes: MASLD SuSiE/ABF-COLOC (PP.H4 > 0.5) vs Broadaway 747")

ggsave(file.path(FIGS04_DIR, "broadaway_novel_source.pdf"),
       pB, width = 108, height = 55, units = "mm", useDingbats = FALSE)

cat("Wrote:\n  ", file.path(FIGS04_DIR, "broadaway_overlap_venn.pdf"),
    "\n  ", file.path(FIGS04_DIR, "broadaway_novel_source.pdf"), "\n")
