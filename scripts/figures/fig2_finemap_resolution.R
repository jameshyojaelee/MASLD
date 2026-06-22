#!/usr/bin/env Rscript
# fig2_finemap_resolution.R  (2026-06-17)
# Redesign of the cluttered single-bar "fine-mapping cascade" (fig2_finemap_cascade.R,
# left untouched). ONE unified figure (single title, single shared x-axis, no panel
# tags) with two stacked bar blocks on a common "Loci" axis:
#
#   Consolidation block :  Total independent / Shared (>=2 GWAS) / Unique (1 GWAS)
#   Fine-mapper block   :  SuSiE, CARMA (within-ancestry) ; SuSiEx, meSuSiE (cross-
#                          ancestry). Each tool's locus set is segmented by the BEST
#                          variant PIP in the locus:
#                            variant PIP >= 0.9   (resolved / high-confidence)
#                            best PIP 0.5-0.9     (moderate)
#                            no variant > 0.5     (unresolved)
#   NB: each tool ran on its OWN locus set (different denominators) — SuSiE on all 184
#   independent loci, CARMA on the subset it covered, SuSiEx on trait-specific loci,
#   meSuSiE on EUR-EAS converged loci — so the bar-end totals differ by design.
#
# All counts read live from frozen on-disk sources so they always match the text:
#   - merged_loci_map.csv          -> consolidation + SuSiE max_pip per independent locus
#   - carma_all_results.csv        -> CARMA per-variant PIP (mapped to independent loci)
#   - susiex_4way/susiex_cs_4way.csv-> SuSiEx OVRL_PIP per trait-specific locus
#   - mesusie_polyfun/mesusie_locus_summary.csv -> meSuSiE max_pip per converged locus
#
# Out: figures/main/fig2_genetics/panels/fig2B_finemap_cascade.pdf
#      (replaces the old single-bar cascade; supersedes fig2_finemap_cascade.R)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggnewscale)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")
FM        <- file.path(BASE, "GWAS/finemapping/results")

# Single shared y-axis: row order top -> bottom (one blank spacer separates blocks)
Y_T2B <- c("Total independent", "Unique (1 GWAS)", "Shared (≥2 GWAS)", " ",
           "SuSiE", "CARMA", "SuSiEx", "meSuSiE")
yfac  <- function(v) factor(v, levels = rev(Y_T2B))

# ---------------------------------------------------------------------------
# Consolidation block (loci)
# ---------------------------------------------------------------------------
lm <- fread(file.path(FM, "merged_loci_map.csv"))
raw_screened <- nrow(lm)
ind  <- lm[kept == TRUE]
n_total  <- nrow(ind); n_shared <- sum(ind$group_size >= 2); n_unique <- sum(ind$group_size == 1)

dA <- data.table(cat = c("Total independent", "Shared (≥2 GWAS)", "Unique (1 GWAS)"),
                 n   = c(n_total, n_shared, n_unique))
dA[, `:=`(catf = factor(cat, levels = c("Total independent", "Shared (≥2 GWAS)", "Unique (1 GWAS)")),
          y    = yfac(cat))]
colA <- c("Total independent" = "#9E9E9E", "Shared (≥2 GWAS)" = "#00695C", "Unique (1 GWAS)" = "#1565C0")

# ---------------------------------------------------------------------------
# Fine-mapper block (loci classified by best-variant PIP)
# ---------------------------------------------------------------------------
cls <- function(mx) c(hi = sum(mx >= 0.9), mod = sum(mx >= 0.5 & mx < 0.9), lo = sum(mx < 0.5))
susie_b <- cls(ind$max_pip)

ca  <- fread(file.path(FM, "carma_all_results.csv"))
key <- lm[, .(study, original_locus, merged_locus, kept)]
ca  <- merge(ca, key, by.x = c("study", "locus"), by.y = c("study", "original_locus"), all.x = TRUE)
carma_b <- cls(ca[kept == TRUE, .(mx = max(PIP, na.rm = TRUE)), by = merged_locus]$mx)

sx <- fread(file.path(FM, "susiex_4way/susiex_cs_4way.csv"))
susiex_b <- cls(sx[, .(mx = max(OVRL_PIP, na.rm = TRUE)), by = locus_id]$mx)

me <- fread(file.path(FM, "mesusie_polyfun/mesusie_locus_summary.csv"))
mesusie_b <- cls(me[converged == TRUE]$max_pip)

tools <- rbindlist(lapply(list(
  list("SuSiE", susie_b), list("CARMA", carma_b),
  list("SuSiEx", susiex_b), list("meSuSiE", mesusie_b)),
  function(x) data.table(tool = x[[1]], hi = x[[2]]["hi"], mod = x[[2]]["mod"], lo = x[[2]]["lo"])))

dB <- melt(tools, id.vars = "tool", measure.vars = c("hi", "mod", "lo"),
           variable.name = "pip", value.name = "n")
dB[, `:=`(pip = factor(fcase(pip == "hi", "variant PIP ≥ 0.9",
                             pip == "mod", "best PIP 0.5–0.9",
                             default = "no variant > 0.5"),
                       levels = c("variant PIP ≥ 0.9", "best PIP 0.5–0.9", "no variant > 0.5")),
          y   = yfac(tool))]
tot <- tools[, .(tot = hi + mod + lo, y = yfac(tool)), tool]
colB <- c("variant PIP ≥ 0.9" = "#0D3B66", "best PIP 0.5–0.9" = "#7FB3D5", "no variant > 0.5" = "#D9DCE0")

# ---------------------------------------------------------------------------
# One figure, one axis
# ---------------------------------------------------------------------------
XMAX <- n_total * 1.14

p <- ggplot() +
  # consolidation bars (own fill scale, no legend — labels are self-explanatory)
  geom_col(data = dA, aes(n, y, fill = catf), width = 0.64) +
  geom_text(data = dA, aes(n, y, label = n), hjust = -0.20, size = 2.8,
            fontface = "bold", color = "grey15") +
  scale_fill_manual(values = colA, guide = "none") +
  ggnewscale::new_scale_fill() +
  # fine-mapper bars (PIP fill scale -> the only legend)
  geom_col(data = dB, aes(n, y, fill = pip), width = 0.64) +
  geom_text(data = tot, aes(tot, y, label = tot), hjust = -0.28, size = 2.8,
            fontface = "bold", color = "grey15") +
  scale_fill_manual(values = colB, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, XMAX),
                     breaks = seq(0, 200, 50)) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Loci", y = NULL,
       title = "GWAS fine-mapping: consolidation and resolution") +
  theme_masld(base_size = 9) +
  theme(plot.title    = element_text(size = 9.5, face = "bold"),
        axis.text.y   = element_text(size = 8),
        legend.position = c(0.99, 0.04), legend.justification = c(1, 0),
        legend.key.size = unit(0.32, "cm"), legend.text = element_text(size = 7),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA))

save_fig(p, file.path(PANEL_DIR, "fig2B_finemap_cascade.pdf"),
         width = fig_col_width * 1.15, height = 3.9)
cat(sprintf("[finemap_resolution] loci %d/%d/%d (tot/shared/uniq) | per-tool loci by best PIP (hi/mod/lo):\n",
            n_total, n_shared, n_unique))
print(tools[, .(tool, hi, mod, lo, total = hi + mod + lo)])
