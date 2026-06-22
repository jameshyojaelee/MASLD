#!/usr/bin/env Rscript
# ============================================================================
# fibrosis_stage_hallmark_escalation.R
# Supp Fig S02 — Hallmark program escalation/suppression across fibrosis stage.
#
# Dotplot: stage (x) x Hallmark pathway (y). Restricted to pathways whose NES
# trends monotonically across F1..F4 vs F0. Inflammatory/EMT/TGFbeta programs
# intensify; OXPHOS/fatty-acid metabolism deepen suppression.
#
# Source: fibrosis_stage_gsea.csv (fgsea on Hallmark, per stage-vs-F0 contrast)
#
# Output: figures/supplementary/figS02_progression/fibrosis_stage_hallmark_escalation.pdf
#   sized 110 x 130 mm
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
OUT_PDF  <- file.path(OUT_DIR, "fibrosis_stage_hallmark_escalation.pdf")

# ---------------------------------------------------------------------------
# Load GSEA results
# ---------------------------------------------------------------------------
g <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_gsea.csv"))

stage_levels <- c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")
g <- g[contrast %in% stage_levels]
g[, stage := factor(contrast, levels = stage_levels, labels = c("F1", "F2", "F3", "F4"))]
g[, stage_idx := as.integer(stage)]

# ---------------------------------------------------------------------------
# Monotonic-trend filter: keep pathways present in all 4 stages whose NES is
# monotonic (all consecutive diffs same sign, ties allowed) across stages.
# ---------------------------------------------------------------------------
g4 <- g[, if (.N == 4) .SD, by = pathway]
setorder(g4, pathway, stage_idx)

mono_check <- g4[, {
  d <- diff(NES)
  inc <- all(d >= 0)
  dec <- all(d <= 0)
  rho <- suppressWarnings(cor(NES, stage_idx, method = "spearman"))
  .(monotonic = inc | dec, direction = ifelse(inc, "up", "down"),
    span = abs(NES[4] - NES[1]), rho = rho)
}, by = pathway]

keep <- mono_check[monotonic == TRUE]
# Rank by trend magnitude; cap to top ~24 for legibility
setorder(keep, -span)
if (nrow(keep) > 24) keep <- keep[1:24]

dat <- g4[pathway %in% keep$pathway]

# Order pathways by direction then span (suppressed at bottom, induced at top)
ord <- merge(keep, dat[stage == "F4", .(pathway, NES_f4 = NES)], by = "pathway")
setorder(ord, NES_f4)
dat[, pathway := factor(pathway, levels = ord$pathway)]

# Pretty labels: strip HALLMARK_ + underscores -> Title Case
pretty <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}
lab_map <- setNames(pretty(levels(dat$pathway)), levels(dat$pathway))

dat[, neglog10padj := -log10(pmax(padj, 1e-300))]

fwrite(dat[, .(pathway, stage, NES, padj, neglog10padj)],
       file.path(DATA_DIR, "fibrosis_stage_hallmark_escalation.csv"))

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
emt_key <- grep("EPITHELIAL_MESENCHYMAL", levels(dat$pathway), value = TRUE)
cat(sprintf("[hero] n pathways shown = %d (monotonic across F1..F4)\n", length(levels(dat$pathway))))
if (length(emt_key) == 1) {
  emt <- dat[pathway == emt_key]
  cat(sprintf("[hero] EMT NES F1=%.2f -> F4=%.2f\n",
              emt[stage == "F1", NES], emt[stage == "F4", NES]))
}
for (k in c("TGF_BETA", "OXIDATIVE_PHOSPHORYLATION", "FATTY_ACID_METABOLISM")) {
  pk <- grep(k, levels(dat$pathway), value = TRUE)
  if (length(pk) == 1) {
    pp <- dat[pathway == pk]
    cat(sprintf("[hero] %s NES F1=%.2f -> F4=%.2f\n", k,
                pp[stage == "F1", NES], pp[stage == "F4", NES]))
  } else {
    cat(sprintf("[hero] %s NOT in monotonic set\n", k))
  }
}

# ---------------------------------------------------------------------------
# Dotplot
# ---------------------------------------------------------------------------
nmax <- max(abs(dat$NES))
p <- ggplot(dat, aes(x = stage, y = pathway)) +
  geom_point(aes(size = neglog10padj, fill = NES), shape = 21, color = "grey30", stroke = 0.2) +
  scale_fill_gradient2(low = "#1565C0", mid = "#F5F5F5", high = "#C9265E",
                       midpoint = 0, limits = c(-nmax, nmax), name = "NES") +
  scale_size_continuous(range = c(0.6, 4), name = expression(-log[10]~padj)) +
  scale_y_discrete(labels = lab_map) +
  labs(x = "Fibrosis stage (vs F0)", y = NULL,
       title = "Hallmark program escalation across fibrosis stage") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    axis.text.x     = element_text(size = 6.5, face = "bold"),
    axis.text.y     = element_text(size = 5.6),
    legend.position = "right",
    legend.text     = element_text(size = 5.2),
    legend.title    = element_text(size = 5.8),
    legend.key.size = unit(0.22, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 110 / 25.4,
       height = 130 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
