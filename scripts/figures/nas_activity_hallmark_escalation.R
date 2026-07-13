#!/usr/bin/env Rscript
# ============================================================================
# nas_activity_hallmark_escalation.R
# Supp Fig S02 — Hallmark program escalation/suppression across NAS activity.
#
# Dotplot: NAS-activity bin (x) x Hallmark pathway (y). Restricted to pathways
# whose NES trends monotonically across NAS 1-2 / 3-4 / 5-8 vs NAS 0.
# NAS-activity analog of fibrosis_stage_hallmark_escalation.R (the two read as a
# pair: fibrosis-axis vs activity-axis dose-response).
#
# Source: nas_grouped_gsea.csv (fgsea on Hallmark, per grouped-NAS contrast; 14g),
#   itself built on the TREAT-canonical grouped-NAS LVQW contrast (14e).
#
# Output: figures/supplementary/figS02_progression/nas_activity_hallmark_escalation.pdf
#   sized 100 x 130 mm
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS02_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DATA_DIR <- file.path(OUT_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF  <- file.path(OUT_DIR, "nas_activity_hallmark_escalation.pdf")

# ---------------------------------------------------------------------------
# Load GSEA results
# ---------------------------------------------------------------------------
g <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_grouped_gsea.csv"))

bin_levels <- c("NAS1-2_vs_NAS0", "NAS3-4_vs_NAS0", "NAS5-8_vs_NAS0")
g <- g[contrast %in% bin_levels]
g[, nas := factor(contrast, levels = bin_levels, labels = c("1–2", "3–4", "5–8"))]
g[, nas_idx := as.integer(nas)]

# ---------------------------------------------------------------------------
# Monotonic-trend filter: keep pathways present in all 3 bins whose NES is
# monotonic (all consecutive diffs same sign, ties allowed) across bins.
# ---------------------------------------------------------------------------
g3 <- g[, if (.N == 3) .SD, by = pathway]
setorder(g3, pathway, nas_idx)

mono_check <- g3[, {
  d <- diff(NES)
  inc <- all(d >= 0)
  dec <- all(d <= 0)
  rho <- suppressWarnings(cor(NES, nas_idx, method = "spearman"))
  .(monotonic = inc | dec, direction = ifelse(inc, "up", "down"),
    span = abs(NES[3] - NES[1]), rho = rho)
}, by = pathway]

keep <- mono_check[monotonic == TRUE]
# Rank by trend magnitude; cap to top ~24 for legibility
setorder(keep, -span)
if (nrow(keep) > 24) keep <- keep[1:24]

dat <- g3[pathway %in% keep$pathway]

# Order pathways by NES in the highest bin (suppressed at bottom, induced at top)
ord <- merge(keep, dat[nas == "5–8", .(pathway, NES_hi = NES)], by = "pathway")
setorder(ord, NES_hi)
dat[, pathway := factor(pathway, levels = ord$pathway)]

# Pretty labels: strip HALLMARK_ + underscores -> Title Case
pretty <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}
lab_map <- setNames(pretty(levels(dat$pathway)), levels(dat$pathway))

dat[, neglog10padj := -log10(pmax(padj, 1e-300))]

fwrite(dat[, .(pathway, nas, NES, padj, neglog10padj)],
       file.path(DATA_DIR, "nas_activity_hallmark_escalation.csv"))

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
cat(sprintf("[hero] n pathways shown = %d (monotonic across NAS 1-2..5-8)\n",
            length(levels(dat$pathway))))
for (k in c("EPITHELIAL_MESENCHYMAL", "TGF_BETA", "INFLAMMATORY",
            "OXIDATIVE_PHOSPHORYLATION", "FATTY_ACID_METABOLISM")) {
  pk <- grep(k, levels(dat$pathway), value = TRUE)
  if (length(pk) == 1) {
    pp <- dat[pathway == pk]
    cat(sprintf("[hero] %s NES NAS1-2=%.2f -> NAS5-8=%.2f\n", k,
                pp[nas == "1–2", NES], pp[nas == "5–8", NES]))
  } else {
    cat(sprintf("[hero] %s NOT in monotonic set\n", k))
  }
}

# ---------------------------------------------------------------------------
# Dotplot
# ---------------------------------------------------------------------------
nmax <- max(abs(dat$NES))
p <- ggplot(dat, aes(x = nas, y = pathway)) +
  geom_point(aes(size = neglog10padj, fill = NES), shape = 21, color = "grey30", stroke = 0.2) +
  scale_fill_gradient2(low = "#1565C0", mid = "#F5F5F5", high = "#C9265E",
                       midpoint = 0, limits = c(-nmax, nmax), name = "NES") +
  scale_size_continuous(range = c(0.6, 4), name = expression(-log[10]~padj)) +
  scale_y_discrete(labels = lab_map) +
  labs(x = "NAS activity (vs NAS 0)", y = NULL,
       title = "Hallmark program escalation across NAS activity") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "plain", margin = margin(b = 6)),
    axis.text.x     = element_text(size = 6.5, face = "plain"),
    axis.text.y     = element_text(size = 5.6),
    legend.position = "right",
    legend.text     = element_text(size = 5.2),
    legend.title    = element_text(size = 5.8),
    legend.key.size = unit(0.22, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 100 / 25.4,
       height = 130 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
