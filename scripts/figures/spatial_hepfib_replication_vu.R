#!/usr/bin/env Rscript
# KEY MESSAGE: the GSE192741 spatial immune-exclusion finding -- hepatocytes and
# fibroblasts occupy mutually exclusive spatial territory -- REPLICATES in the
# independent Vu et al. 2025 Visium cohort. The hepatocyte<->fibroblast per-spot
# co-occurrence Spearman rho is negative in BOTH cohorts (GSE192741 = -0.562,
# Vu = -0.301), and every one of the 5 Vu donors is negative. Computed WITHIN
# each cohort (never pooled, batch rule).
#
# Panel A: hep-fib co-occurrence rho, GSE192741 vs Vu (Vu shows 5 per-donor points).
# Panel B: cell-type composition of immune-excluded-fibrotic spots in Vu
#          (fibroblast-rich, hepatocyte-poor -- same as GSE192741).
#
# Sources:
#   GSE : Analysis/Spatial/results/immune_exclusion/celltype_cooccurrence_spearman.csv (Hep-Fib = -0.562)
#   Vu  : Analysis/Spatial/results/immune_exclusion/vu/{hepfib_cooccurrence_per_section,celltype_mean_by_spot_class}.csv
# Output: figures/main/fig4_validation/spatial_hepfib_replication_vu.pdf
# Env: rnaseq
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
IE <- file.path(BASE, "Analysis/Spatial/results/immune_exclusion")

# ── hep-fib co-occurrence rho per cohort ──────────────────────────────────────
gse <- read.csv(file.path(IE, "celltype_cooccurrence_spearman.csv"), stringsAsFactors = FALSE)
names(gse) <- c("ct1", "ct2", "rho")
gse_hf <- gse$rho[(gse$ct1 == "Hepatocytes" & gse$ct2 == "Fibroblasts")][1]

vu_sec <- read.csv(file.path(IE, "vu/hepfib_cooccurrence_per_section.csv"), stringsAsFactors = FALSE)
vu_donor <- vu_sec %>% group_by(donor) %>% summarise(rho = mean(hepfib_rho), .groups = "drop")
vu_hf <- mean(vu_donor$rho)

pts <- data.frame(cohort = "Vu (n=5)", rho = vu_donor$rho)
bars <- data.frame(
  cohort = factor(c("GSE192741 (n=5)", "Vu (n=5)"), levels = c("GSE192741 (n=5)", "Vu (n=5)")),
  rho = c(gse_hf, vu_hf))

pA <- ggplot(bars, aes(rho, cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_segment(aes(x = 0, xend = rho, yend = cohort), linewidth = 1.1, color = masld_colors$conserved) +
  geom_point(size = 3, color = masld_colors$conserved) +
  geom_point(data = pts, aes(rho, cohort), color = "gray35", size = 1.4,
             position = position_nudge(y = 0.18), alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rho)), nudge_y = -0.22,
            size = 3.1, fontface = "bold", color = masld_colors$conserved) +
  scale_x_continuous(limits = c(-0.74, 0.04), breaks = c(-0.6, -0.4, -0.2, 0)) +
  labs(x = "Hepatocyte-Fibroblast spatial\nco-occurrence (Spearman rho)", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10, color = "black"),
        panel.grid.major.y = element_blank())

# ── Vu composition: fibroblast enrichment in excluded vs non-excluded spots ───
comp <- read.csv(file.path(IE, "vu/celltype_mean_by_spot_class.csv"), stringsAsFactors = FALSE)
keep <- c("Hepatocytes", "Fibroblasts", "Macrophages")
cl_levels <- c("non_excluded", "immune_excluded_fibrotic")
cdf <- comp %>% filter(spot_class %in% cl_levels) %>%
  select(spot_class, all_of(keep)) %>%
  tidyr::pivot_longer(-spot_class, names_to = "ct", values_to = "abund") %>%
  mutate(spot_class = factor(ifelse(spot_class == "immune_excluded_fibrotic",
                                     "Immune-excluded", "Non-excluded"),
                             levels = c("Non-excluded", "Immune-excluded")),
         ct = factor(ct, levels = keep))
pB <- ggplot(cdf, aes(ct, abund, fill = spot_class)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = c("Non-excluded" = masld_colors$ns,
                               "Immune-excluded" = masld_colors$up), name = NULL) +
  labs(x = NULL, y = "Mean cell abundance, Vu spots") +
  theme_masld() + theme_pub() +
  theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = c(0.72, 0.86), legend.text = element_text(size = 9),
        legend.key.size = unit(0.32, "cm"))

p_out <- (pA | pB) + plot_layout(widths = c(1, 1))
out <- file.path(FIG4_DIR, "_supp", "spatial_hepfib_replication_vu.pdf")
cairo_pdf(out, width = fig_full_width * 0.86, height = 2.8)
print(p_out); dev.off()
message(sprintf("hep-fib rho: GSE192741 = %.3f | Vu = %.3f (per-donor: %s)",
                gse_hf, vu_hf, paste(sprintf("%.2f", vu_donor$rho), collapse = ", ")))
message("Saved: ", out)
