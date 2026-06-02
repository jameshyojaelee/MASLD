#!/usr/bin/env Rscript
# ============================================================================
# fig2_panel_hotspot_cascade.R
#
# KEY MESSAGE: Two opposing hepatocyte module trajectories across MASLD severity:
#   GAIN — oxidative stress (hep-24: SOD2/HKDC1), AP-1 injury (hep-26: JUN/ATF3),
#           metabolic reprogramming (hep-20: BICC1/GLS/ITGAV)
#   LOSS — hepatocyte identity (hep-19: HNF4A-AS1/EHHADH), hepatic synthesis
#           (hep-27: CFH/NR1H4)
# Shows only the 5 modules actively discussed in the manuscript.
#
# Cascade restricted to 3 stages (Healthy/ST/SH) after protocol contamination
# remediation 2026-05-22 — all Cirrhosis donors come from GSE136103 NPC-enriched
# protocol and are excluded. Cirrhosis biology is presented as a supplementary
# panel with explicit n=19 Andrews/Gribben single-dataset caveat (separate task).
#
# Output: FIG2_DIR/panels/fig2_hotspot_cascade.pdf
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

HS_RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
META_FILE <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")

# Liang disease-stage gradient (publication_theme.R: masld_colors)
# 3-stage cascade (Cirrhosis removed 2026-05-22 — see header note)
STAGE_LEVELS <- c("Healthy", "Steatosis", "Steatohepatitis")
STAGE_COLORS <- c(Healthy        = masld_colors$control,   # #9E9E9E
                  Steatosis      = masld_colors$masl,      # #F4A674
                  Steatohepatitis= masld_colors$mash)      # #C9265E

# Five hepatocyte modules ordered to trace the directional cascade:
# gain-of-stress, then loss-of-mature-function (effect direction informs order).
# Audit 2026-05-17: hep__24 actually matches the bulk-NMF "Progression-
# Inflammatory" program (3-gene overlap HKDC1/AKR1B10/DEFB1, hyper p=0.007 —
# NOT the labelled "P1 vascular-smooth-muscle"; the original P1 tag was the W
# column index, not the canonical program label). The (bulkNMF) marker is
# retained but documented honestly in the README.
mods <- fread(file.path(HS_RES, "all_modules.tsv"))
mods <- mods[cell_type == "hepatocytes"]

# Five hepatocyte modules actively discussed in the manuscript results section.
# Ordered: gaining modules first (increasing with severity), then losing modules.
FOCAL_MODS <- c(20L, 24L, 26L, 19L, 27L)
FOCAL_LABELS <- c(
  "20"  = "Hep-20\nBICC1 · GLS · ITGAV\n(gain with severity)",
  "24"  = "Hep-24\nSOD2 · HKDC1 · SQSTM1\n(oxidative stress, gain)",
  "26"  = "Hep-26\nJUN · ATF3 · SERPINE1\n(AP-1 injury, gain)",
  "19"  = "Hep-19\nHNF4A-AS1 · EHHADH · BDH1\n(hepatocyte identity, loss)",
  "27"  = "Hep-27\nCFH · NR1H4 · C8A\n(hepatic synthesis, loss)"
)

SELECT <- mods[module %in% FOCAL_MODS,
               .(module_int = module, disease_stage_beta, disease_stage_q,
                 stability_score, bulk_replicated)]
SELECT[, tag         := FOCAL_LABELS[as.character(module_int)]]
SELECT[, facet_label := factor(tag, levels = FOCAL_LABELS)]
setorder(SELECT, facet_label)

# Donor module scores
ds <- fread(file.path(HS_RES, "donor_scores_all.tsv"))
ds <- ds[cell_type == "hepatocytes"]
ds[, module_int := as.integer(sub(".*__", "", module))]
ds <- ds[module_int %in% SELECT$module_int]

# Pull exclude_stage_analysis flag so GSE136103 + Liver_Atlas NPC-enriched
# donors are dropped from the cascade (set by protocol-contamination
# remediation 2026-05-22). Read header first to decide whether the column
# exists in this version of donor_metadata_extended.tsv.
meta_header <- names(fread(META_FILE, nrows = 0))
meta_cols   <- intersect(
  c("sample", "disease_stage_coarse", "dataset", "exclude_stage_analysis"),
  meta_header)
meta <- fread(META_FILE, select = meta_cols)
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
ds <- ds[meta, on = "sample", nomatch = 0]
ds <- ds[disease_stage_coarse %in% STAGE_LEVELS]
ds[, disease_stage_coarse := factor(disease_stage_coarse, levels = STAGE_LEVELS)]
ds <- ds[SELECT[, c("module_int", "facet_label")], on = "module_int", nomatch = 0]
ds[, facet_label := factor(facet_label, levels = SELECT$facet_label)]

# Dataset-center per (module, dataset) — matches the (1|dataset) random intercept
# that 505 uses for the β estimate; raw pooled scores are confounded by cohort
# composition (Wang heavy in Steatosis, Liver_Atlas all Healthy, etc.).
ds[, score_centered := score - mean(score, na.rm = TRUE),
   by = .(module_int, dataset)]

# Per-donor n for caption
n_by_stage <- ds[, .(n = uniqueN(sample)), by = disease_stage_coarse][
  order(disease_stage_coarse)]

# Median trend per facet (for line overlay)
trend <- ds[, .(med = median(score_centered, na.rm = TRUE)),
            by = .(facet_label, disease_stage_coarse)]

p <- ggplot(ds, aes(disease_stage_coarse, score_centered,
                    fill = disease_stage_coarse)) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.2,
              color = "black", alpha = 0.85) +
  geom_line(data = trend, aes(x = disease_stage_coarse, y = med, group = 1),
            inherit.aes = FALSE, color = "gray25", linewidth = 0.35) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.1,
               color = "black", fill = "white", stroke = 0.3) +
  facet_wrap(~ facet_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = STAGE_COLORS, name = NULL, guide = "none") +
  labs(x = NULL,
       y = "Donor module score (dataset-centered)") +
  theme_masld(base_size = 11) +
  theme(axis.text.x  = element_text(size = 10, angle = 30, hjust = 1),
        axis.text.y  = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text   = element_text(size = 11, face = "bold", lineheight = 1.1))

save_fig(p, file.path(PANEL_DIR, "fig2_hotspot_cascade.pdf"),
         width = fig_full_width * 2.5, height = 3.8)

# Side data for caption transparency
fwrite(SELECT, file.path(DATA_DIR, "fig2_hotspot_cascade_modules.csv"))
fwrite(ds[, .(sample, module_int, score, disease_stage_coarse, dataset)],
       file.path(DATA_DIR, "fig2_hotspot_cascade_donor_scores.csv"))
cat("Wrote fig2_hotspot_cascade.pdf\n")
