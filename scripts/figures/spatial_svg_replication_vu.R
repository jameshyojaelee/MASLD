#!/usr/bin/env Rscript
# KEY MESSAGE: the GSE192741 spatial-structure findings REPLICATE in an
# independent Visium cohort (Vu et al. 2025, 5 donors / 10 slides), killing the
# single-cohort (n=5) vulnerability. Of the GSE192741 steatotic SVGs that are
# detectable in Vu, 96% are also spatially variable in Vu, and per-gene Moran's I
# correlates at Spearman rho = 0.55 (healthy 97% / 0.60). The disease-EMERGENT
# SVGs (the novel gain-of-structure set) replicate at the same rate.
#
# Numbers are recomputed live from the CURRENT canonical SVG sets (the stale
# concordance_gse192741.csv used a retired 451-SVG steatotic set; current = 365).
# Within-cohort metric per cohort; NEVER pooled (Vu VLP* vs GSE192741 JBO* batch).
# Vu carries no per-donor F-stage, so this is SVG-PRESENCE replication, not a
# stage-stratified or disease-direction claim.
#
# Sources:
#   Vu  : Analysis/Spatial/results/validation_vu/consensus_svgs_vu.csv (mean_morans_i, consensus_svg)
#   GSE : Analysis/Spatial/results/svg/svgs_{Steatotic,Healthy}.csv (morans_i_mean, svg)
#         Analysis/Spatial/results/svg/differential_svgs.csv (category == disease_emergent_SVG)
#
# Output: figures/main/fig5_molecular_context/spatial_svg_replication_vu.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(ggrepel)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
SP <- file.path(BASE, "Analysis/Spatial/results")

# ── Load ──────────────────────────────────────────────────────────────────────
vu <- read.csv(file.path(SP, "validation_vu/consensus_svgs_vu.csv"), stringsAsFactors = FALSE) %>%
  transmute(gene, vu_moran = mean_morans_i, vu_svg = tolower(as.character(consensus_svg)) == "true")
load_gse <- function(f) {
  d <- read.csv(file.path(SP, "svg", f), stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(gene = d[[1]], gse_moran = as.numeric(d$morans_i_mean),
             gse_svg = tolower(as.character(d$svg)) == "true", stringsAsFactors = FALSE)
}
steat <- load_gse("svgs_Steatotic.csv"); healthy <- load_gse("svgs_Healthy.csv")
diff  <- read.csv(file.path(SP, "svg/differential_svgs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
emergent <- diff[[1]][diff$category == "disease_emergent_SVG"]

# ── Concordance (recomputed; within shared tested universe) ──────────────────
repl_stats <- function(gse) {
  m <- inner_join(vu, gse, by = "gene")                 # shared tested universe
  svg_shared <- m %>% filter(gse_svg)
  n_repl <- sum(svg_shared$vu_svg)
  rho <- suppressWarnings(cor(m$gse_moran, m$vu_moran, method = "spearman"))
  ft  <- fisher.test(table(factor(m$gse_svg, c(FALSE, TRUE)), factor(m$vu_svg, c(FALSE, TRUE))))
  list(m = m, n_shared = nrow(m), n_svg = nrow(svg_shared), n_repl = n_repl,
       pct = 100 * n_repl / nrow(svg_shared), rho = rho, OR = unname(ft$estimate), p = ft$p.value)
}
S <- repl_stats(steat); H <- repl_stats(healthy)
# disease-emergent replication
em_shared <- intersect(emergent, vu$gene)
em_repl   <- sum(vu$vu_svg[match(em_shared, vu$gene)])
message(sprintf("[svg_repl] STEATOTIC: %d/%d (%.0f%%) GSE SVGs (testable in Vu) replicate; Moran rho=%.2f; Fisher OR=%.1f p=%.1e",
                S$n_repl, S$n_svg, S$pct, S$rho, S$OR, S$p))
message(sprintf("[svg_repl] HEALTHY:   %d/%d (%.0f%%); Moran rho=%.2f; OR=%.1f",
                H$n_repl, H$n_svg, H$pct, H$rho, H$OR))
message(sprintf("[svg_repl] DISEASE-EMERGENT: %d/%d (%.0f%%) of the 94 gain-of-structure SVGs detectable in Vu also Vu-SVG",
                em_repl, length(em_shared), 100*em_repl/length(em_shared)))
message(sprintf("[svg_repl] CAVEAT: within-cohort (not pooled); denom = GSE SVGs detectable in Vu (%d/%d steatotic GSE SVGs are in Vu's tested universe; the rest undetected, FFPE/depth). Vu has no per-donor F-stage -> SVG-presence replication.",
                S$n_svg, sum(steat$gse_svg)))

# ── Panel A: replication-rate bars ───────────────────────────────────────────
rate_df <- data.frame(
  lab = factor(c("Healthy SVGs", "Steatotic SVGs", "Disease-emergent"),
               levels = c("Disease-emergent", "Steatotic SVGs", "Healthy SVGs")),
  pct = c(H$pct, S$pct, 100*em_repl/length(em_shared)),
  num = c(H$n_repl, S$n_repl, em_repl),
  den = c(H$n_svg, S$n_svg, length(em_shared)))
pA <- ggplot(rate_df, aes(pct, lab)) +
  geom_segment(aes(x = 0, xend = pct, yend = lab), linewidth = 1.1, color = masld_colors$conserved) +
  geom_point(size = 2.8, color = masld_colors$conserved) +
  geom_text(aes(label = sprintf("%.0f%% (%d/%d)", pct, num, den)), hjust = 0, nudge_x = 2.5,
            size = 3.2, fontface = "plain", color = "black") +
  scale_x_continuous(limits = c(0, 118), breaks = c(0, 50, 100),
                     labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0))) +
  labs(x = "GSE192741 SVGs replicating in Vu", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10, color = "black"),
        panel.grid.major.y = element_blank(), plot.margin = margin(6, 40, 6, 6))

# ── Panel B: per-gene Moran's I scatter (steatotic) ──────────────────────────
m <- S$m %>% mutate(cls = ifelse(gse_svg & vu_svg, "replicating SVG",
                          ifelse(gse_svg & !vu_svg, "GSE-only", "other")))
lab_genes <- m %>% filter(gse_svg & vu_svg) %>% arrange(desc(gse_moran + vu_moran)) %>% slice_head(n = 8)
pB <- ggplot(m, aes(gse_moran, vu_moran)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_point(data = filter(m, cls == "other"), color = masld_colors$ns, size = 0.5, alpha = 0.3) +
  geom_point(data = filter(m, cls == "GSE-only"), color = "#E8A33D", size = 1.1, alpha = 0.8) +
  geom_point(data = filter(m, cls == "replicating SVG"), color = masld_colors$conserved, size = 1.1, alpha = 0.8) +
  geom_text_repel(data = lab_genes, aes(label = gene), size = 3.0, fontface = "italic",
                  color = "gray15", max.overlaps = Inf, seed = 1, box.padding = 0.4) +
  annotate("text", x = max(m$gse_moran, na.rm = TRUE) * 0.98, y = 0.02,
           label = sprintf("Spearman rho = %.2f", S$rho), hjust = 1, vjust = 0,
           size = 3.4, fontface = "plain", color = "gray20") +
  labs(x = "GSE192741 Moran's I (steatotic)", y = "Vu Moran's I") +
  theme_masld() + theme_pub() +
  theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10, color = "black"))

p_out <- (pA | pB) + plot_layout(widths = c(1, 1.25))
out <- file.path(FIG4_DIR, "_supp", "spatial_svg_replication_vu.pdf")
cairo_pdf(out, width = fig_full_width * 0.92, height = 2.7)
print(p_out); dev.off()
message("Saved: ", out)
