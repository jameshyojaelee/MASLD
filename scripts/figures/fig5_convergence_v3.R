##############################################################################
# fig5_convergence_v3.R — Fig 5 REDESIGN (2026-04-15, updated 2026-04-17)
#
# Binding spec: docs/manuscript/FIGURE_PLAN_REVISED.md §"Figure 5"
# Retractions enforced: D5 numbers banned (OR=1.67, 3.6%, Jaccard=0.012,
# "0/5,605 at 6/6 layers"). No "true NASH core" or "discrete F2 switch"
# language. Continuous 1.25x / Wilcoxon p=4.3e-43 only permitted in text.
#
# Panels:
#   5a  therapeutic-axis multi-evidence matrix heatmap
#       (written directly by fig5_convergence.R → panels/fig5a_therapeutic_axes.pdf)
#   5b  atlas convergence modality-count histogram
#
# Output: $FIG5_DIR/panels/fig5b.pdf + $FIG5_DIR/fig5b_convergence_histogram.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUTDIR  <- FIG5_DIR
PANDIR  <- file.path(OUTDIR, "panels")
dir.create(PANDIR, recursive = TRUE, showWarnings = FALSE)

# panels/fig5a_therapeutic_axes.pdf is written directly by fig5_convergence.R — no copy needed here.

# ═══════════════════════════════════════════════════════════════════════════
# Load atlas
# ═══════════════════════════════════════════════════════════════════════════
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat(sprintf("Atlas: %d genes x %d cols\n", nrow(atlas), ncol(atlas)))

# ═══════════════════════════════════════════════════════════════════════════
# Fig 5b uses the canonical 46d convergence table for modality counts.
# Do not maintain an independent local count here: this previously drifted by
# keeping mouse DEG as a local channel after mouse/cross-species had been dropped
# from the canonical Fig 5 convergence score.
# ═══════════════════════════════════════════════════════════════════════════

conv_file <- file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence.csv")
if (!file.exists(conv_file)) {
  stop("Missing canonical convergence table: ", conv_file,
       "\nRun RNA-seq/46d_convergence_evidence.R first.")
}
conv <- fread(conv_file, select = c("human_symbol", "n_modalities_active",
                                    "excluded_from_ranking"))
atlas <- merge(atlas, conv, by = "human_symbol", all.x = TRUE)
atlas[is.na(n_modalities_active), n_modalities_active := 0L]
atlas[is.na(excluded_from_ranking), excluded_from_ranking := TRUE]

sa_dist <- atlas[excluded_from_ranking == FALSE,
                 .N, by = .(sa = as.integer(n_modalities_active))][order(sa)]
sa_dist[, sa := as.integer(sa)]
cat("\nCanonical 46d convergence distribution (mouse/cross-species excluded):\n")
print(sa_dist)

max_sa         <- max(sa_dist$sa, na.rm = TRUE)
top_sa_n       <- sa_dist[sa == max_sa, N]
convergent_ge3 <- sum(sa_dist[sa >= 3, N])
convergent_ge4 <- sum(sa_dist[sa >= 4, N])
convergent_ge5 <- sum(sa_dist[sa >= 5, N])
cat(sprintf("\nHeadline: %d genes at max %d/6; %d at >=3/6; %d at >=4/6; %d at >=5/6\n",
            top_sa_n, max_sa, convergent_ge3, convergent_ge4, convergent_ge5))

# ═══════════════════════════════════════════════════════════════════════════
# PANEL 5b — atlas convergence winnowing funnel (cumulative >=k modalities)
# ═══════════════════════════════════════════════════════════════════════════
# Compact horizontal funnel replacing the former 0-6 histogram (2026-07-08):
# the linear 0/1-dominated histogram wasted a full panel and hid the >=3/>=4
# tail that Fig 5a is about. Cumulative counts (genes reaching AT LEAST k active
# modalities) reframe the same canonical 46d numbers as the winnowing that
# isolates Fig 5a's high-convergence genes (>=3). Bars for >=3 are highlighted
# in the conserved-green used elsewhere in Fig 5.

n_total <- sum(sa_dist$N)
cum <- data.table(k = seq_len(max_sa))
cum[, n_ge := vapply(k, function(kk) sum(sa_dist[sa >= kk, N]), numeric(1))]
cum[, lab := factor(sprintf("≥%d", k),
                    levels = rev(sprintf("≥%d", k)))]  # >=1 at top
cum[, hi := k >= 3]                                          # >=3 = Fig 5a genes
cat("\nCumulative convergence funnel (genes with >=k active modalities):\n")
print(cum[, .(k, n_ge)])

# Horizontal cumulative bars. Linear x keeps the winnowing honest (area = count);
# tip labels carry the exact counts so the tiny >=3/>=4 tail stays legible.
p5b <- ggplot(cum, aes(x = n_ge, y = lab)) +
  geom_col(aes(fill = hi), width = 0.72, alpha = 0.95) +
  geom_text(aes(label = comma(n_ge)),
            hjust = -0.18, size = GEOM_TEXT_6PT, family = "Helvetica",
            fontface = "plain") +
  scale_fill_manual(values = c(`FALSE` = "grey65",
                               `TRUE`  = masld_colors$conserved),
                    guide = "none") +
  scale_x_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.30))) +
  # No in-panel title (PI directive 2026-06-11). Total tested (n_total) and the
  # >=3 = Fig 5a linkage belong in the caption (see "[fig5b caption]" below).
  labs(x = "Genes (≥k active modalities)", y = NULL) +
  theme_masld_compact() +
  theme(text = element_text(family = "Helvetica", size = 6, face = "plain"),
        axis.title = element_text(size = 6, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(3, 6, 3, 3))

ggsave(file.path(PANDIR, "fig5b.pdf"), p5b,
       width = 3.24, height = 2.22, device = cairo_pdf)
cat("Saved: panels/fig5b.pdf\n")
message(sprintf(paste0("[fig5b caption] Cumulative convergence funnel: genes ",
  "reaching >=k active evidence modalities (canonical 46d; mouse/cross-species ",
  "excluded), out of %s non-excluded genes tested. Green = >=3 modalities ",
  "(the Fig 5a therapeutic-axis genes)."), comma(n_total)))

# ═══════════════════════════════════════════════════════════════════════════
# Sidecar CSVs: per-bucket histogram + cumulative funnel (as plotted)
# ═══════════════════════════════════════════════════════════════════════════
fwrite(sa_dist, file.path(OUTDIR, "fig5b_convergence_histogram.csv"))
fwrite(cum[, .(k, genes_ge_k = n_ge, is_fig5a_tier = hi)],
       file.path(OUTDIR, "fig5b_convergence_funnel.csv"))

cat("\nDONE.\n")
