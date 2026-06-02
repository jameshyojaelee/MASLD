#!/usr/bin/env Rscript
# KEY MESSAGE: Cross-species conservation enriches for genes validated at the
# protein level (Conserved_Core rho vs overall DEG rho), showing the
# multi-evidence filter is not arbitrary — it preferentially marks robust biology.
# NOTE (G9-012): the specific rho values are NOT hardcoded here — they are
# computed at runtime from the concordance data read below. (Prior stale
# header numbers were removed after the rebuild.)
#
# Output: figures/main/fig4_validation/panels/fig4a_proteomics_overview.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggrastr)  # rasterize background points; keeps PDF small
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

df <- atlas %>%
  select(human_symbol, dream_logFC, dream_padj,
         best_protein_logFC, best_protein_padj, is_conserved) %>%
  filter(!is.na(best_protein_logFC), !is.na(best_protein_padj),
         !is.na(dream_logFC), !is.na(dream_padj),
         best_protein_padj < 0.05) %>%
  mutate(
    is_deg       = dream_padj < 0.05 & abs(dream_logFC) > 0.5,
    is_conserved = as.logical(is_conserved),
    point_class  = case_when(
      is_conserved & is_deg  ~ "conserved",
      is_deg                 ~ "deg",
      TRUE                   ~ "ns"
    )
  )

# Spearman ρ
r_all  <- round(cor(df$dream_logFC[df$is_deg],
                    df$best_protein_logFC[df$is_deg],
                    method = "spearman", use = "complete.obs"), 3)
r_cons <- round(cor(df$dream_logFC[df$is_deg & df$is_conserved],
                    df$best_protein_logFC[df$is_deg & df$is_conserved],
                    method = "spearman", use = "complete.obs"), 3)
n_deg  <- sum(df$is_deg, na.rm = TRUE)
n_cons <- sum(df$is_deg & df$is_conserved, na.rm = TRUE)

# Genes to label — split HKDC1 out for manual nudging
label_genes <- c("THRB", "SERPINE1", "CHI3L1", "GSN", "CAPN2", "HKDC1", "SOD2", "FADS2")
df_label       <- df %>% filter(human_symbol %in% label_genes, human_symbol != "HKDC1", is_deg)
df_label_hkdc1 <- df %>% filter(human_symbol == "HKDC1", is_deg)

# ── Point colors ─────────────────────────────────────────────────────────────
pt_colors <- c(
  ns        = "#E8E8E8",
  deg       = "#E91E63",   # bright magenta (Tier2 disease signal)
  conserved = "#7B1FA2"    # violet/purple (cross-species confirmed)
)
pt_alpha  <- c(ns = 0.3, deg = 0.55, conserved = 0.90)
pt_size   <- c(ns = 0.25, deg = 0.55, conserved = 0.9)

# ── Annotation positions ──────────────────────────────────────────────────────
xlim_r <- max(abs(df$dream_logFC), na.rm = TRUE) * 1.05
ylim_r <- max(abs(df$best_protein_logFC), na.rm = TRUE) * 1.05

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df %>% arrange(point_class),  # ns → deg → conserved painted last
            aes(x = dream_logFC, y = best_protein_logFC,
                color = point_class, alpha = point_class, size = point_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, color = "gray50") +
  rasterise(geom_point(data = filter(df, point_class == "ns"),   shape = 16), dpi = 300) +
  rasterise(geom_point(data = filter(df, point_class == "deg"),  shape = 16), dpi = 300) +
  geom_point(data = filter(df, point_class == "conserved"), shape = 16) +
  scale_color_manual(values = pt_colors,
                     labels  = c(ns = "Non-DEG", deg = "DEG", conserved = "Cross-species DEG"),
                     name    = NULL) +
  scale_alpha_manual(values = pt_alpha, guide = "none") +
  scale_size_manual(values  = pt_size,  guide = "none") +
  geom_text_repel(data = df_label,
                  aes(label = human_symbol),
                  size = PUB_GEOM_TEXT + 0.5, fontface = "bold.italic",
                  box.padding = 0.5, point.padding = 0.4,
                  force = 3, force_pull = 0.5,
                  min.segment.length = 0,
                  segment.size = 0.4, segment.color = "gray30",
                  max.overlaps = Inf,
                  bg.color = "white", bg.r = 0.12,
                  color = "black") +
  # HKDC1 nudged away from CHI3L1 cluster (upper-right)
  geom_text_repel(data = df_label_hkdc1,
                  aes(label = human_symbol),
                  size = PUB_GEOM_TEXT + 0.5, fontface = "bold.italic",
                  nudge_x = -0.55, nudge_y = 0.45,
                  point.padding = 0.5,
                  min.segment.length = 0,
                  segment.size = 0.4, segment.color = "gray30",
                  bg.color = "white", bg.r = 0.12,
                  color = "black") +
  # Correlation annotations — two lines in bottom-right
  annotate("text", x = xlim_r, y = -ylim_r,
           label = sprintf("All DEGs:              ρ = %.3f (n = %d)\nCross-species DEGs: ρ = %.3f (n = %d)",
                           r_all, n_deg, r_cons, n_cons),
           hjust = 1, vjust = 0, size = PUB_GEOM_TEXT + 0.3,
           color = "gray20", fontface = "italic", lineheight = 1.3) +
  coord_cartesian(xlim = c(-xlim_r, xlim_r), ylim = c(-ylim_r, ylim_r)) +
  labs(x = "mRNA log₂FC (dream mega-analysis)",
       y = "Protein log₂FC",
       title = "mRNA–protein concordance") +
  theme_masld() + theme_pub() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_blank())

# ── Save ──────────────────────────────────────────────────────────────────────
out <- file.path(FIG4_DIR, "panels", "fig4a_proteomics_overview.pdf")
pdf(out, width = fig_half_width, height = fig_half_width, useDingbats = FALSE)
print(p)
dev.off()
message("Saved: ", out)
