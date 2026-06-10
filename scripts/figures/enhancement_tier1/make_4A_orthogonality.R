#!/usr/bin/env Rscript
# Tier 4A — convergence REFRAME to evidence orthogonality (Fig 5 lead).
# Uses existing (2026-06-10) multi-evidence outputs; no recompute.
#   (i)  modality correlation heatmap  -> sources are near-orthogonal (|rho|<0.2)
#   (ii) complementarity bars          -> each source adds independent validation
# Individual PDFs, house style. Output: figures/supplementary/figS_convergence_orthogonality/
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
ME  <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUT <- file.path(BASE, "figures/supplementary/figS_convergence_orthogonality")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
style <- function(p) p + theme_masld() + theme_pub()

# (i) modality correlation heatmap
mc <- fread(file.path(ME, "convergence_evidence_modality_correlations.csv"))
lvl <- unique(c(mc$Mod_A, mc$Mod_B))
mc[, `:=`(Mod_A = factor(Mod_A, levels = lvl), Mod_B = factor(Mod_B, levels = lvl))]
p_i <- style(ggplot(mc[Mod_A != Mod_B], aes(Mod_A, Mod_B, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Spearman_rho)), size = 2) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, limits = c(-0.25, 0.25), name = "Spearman rho") +
  labs(x = NULL, y = NULL,
       title = "Evidence sources are near-orthogonal (all |rho| < 0.21)"))
save_fig(p_i, file.path(OUT, "panel_4A_modality_orthogonality.pdf"), width = 4.6, height = 3.8)

# (ii) complementarity: independent validation gain per source (V1 DGIdb, V3 Conserved)
cp <- fread(file.path(ME, "complementarity_pvalues.csv"))
cp <- cp[validation %in% c("V1_DGIdb", "V3_Conserved")]
cp[, validation := factor(validation, labels = c("DGIdb drug targets", "Conserved core"),
                          levels = c("V1_DGIdb", "V3_Conserved"))]
cp[, sig := ifelse(perm_pvalue < 0.05, "perm p < 0.05", "n.s.")]
cp[, source := factor(source, levels = sort(unique(source)))]
p_ii <- style(ggplot(cp, aes(improvement, source, fill = sig)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  facet_wrap(~ validation, nrow = 1) +
  scale_fill_manual(values = c("perm p < 0.05" = "#2E7D32", "n.s." = "#9E9E9E"), name = NULL) +
  labs(x = "Validation enrichment gain (all - without source)", y = "Evidence source",
       title = "Each source adds independent validation signal (complementarity)"))
save_fig(p_ii, file.path(OUT, "panel_4A_complementarity.pdf"), width = 5.4, height = 3)

cat("Wrote 4A orthogonality panels to", OUT, "\n"); print(list.files(OUT))
