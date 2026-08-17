#!/usr/bin/env Rscript
# ============================================================================
# hotspot_cascade.R
#
# KEY MESSAGE: Two opposing hepatocyte module trajectories across MASLD severity:
#   GAIN — NRF2 antioxidant (hep-24: TXNRD1), AP-1 stress (hep-26: ATF3/GDF15),
#           Ductular injury (hep-20: BICC1)
#   LOSS — Fatty-acid / peroxisomal metab (hep-19), Complement (hep-27: CFH)
# Shows only the 5 modules actively discussed in the manuscript.
#
# Cascade restricted to 3 stages (Healthy/ST/SH) after protocol contamination
# remediation 2026-05-22 — all Cirrhosis donors come from GSE136103 NPC-enriched
# protocol and are excluded. Cirrhosis biology is presented as a supplementary
# panel with explicit n=19 Andrews/Gribben single-dataset caveat (separate task).
#
# Output: FIG4_SC_DIR/panels/figs3_hotspot_cascade.pdf
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

PANEL_DIR <- file.path(FIG4_SC_DIR, "panels", "supplementary")
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
# Directional arm per focal module: gain-of-stress first, then loss-of-function.
# Direction is shown IN-panel by the signed β + significance stars, so it is not
# repeated in the strip title — that keeps the longer canonical names from
# overflowing (and clipping in) the narrow facets.
FOCAL_DIR <- c("20" = "gain", "24" = "gain", "26" = "gain",
               "19" = "loss", "27" = "loss")
# Canonical module names (single source of truth: module_names.tsv), joined on
# cell_type + module so the facet labels stay in sync with Fig 3I + the composite.
mod_names <- fread(file.path(HS_RES, "module_names.tsv"))
hep_names <- mod_names[cell_type == "hepatocytes",
                       setNames(module_name, as.character(module))]
# Two-line facet labels: "Hep-N · <canonical name>" word-wrapped at ~22 chars so
# the longer canonical names fit the narrow facets with no left/right clipping.
wrap_label <- function(s, width = 22)
  vapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
FOCAL_LABELS <- setNames(
  wrap_label(sprintf("Hep-%s · %s",
                     names(FOCAL_DIR), hep_names[names(FOCAL_DIR)])),
  names(FOCAL_DIR))

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

# Faint violins show the central 90% of donors (5-95%): the module-score
# distributions have fat tails (e.g. Hep-27 q95 = 1.95 but q99 = 6.86) that
# otherwise stretch the free y-axes and flatten the stage-to-stage shift. The
# bold MEDIAN trajectory and the annotated mixed-model effect (disease_stage_
# beta/q) use the FULL, untrimmed data.
ds[, q_lo := quantile(score_centered, 0.05, na.rm = TRUE), by = module_int]
ds[, q_hi := quantile(score_centered, 0.95, na.rm = TRUE), by = module_int]
ds_violin <- ds[score_centered >= q_lo & score_centered <= q_hi]

# Per-facet effect-size annotation: mixed-model disease-stage beta + significance
# stars. This carries the true magnitude (incl. Hep-26's tail-driven gain, whose
# median barely moves) so severity is communicated honestly, not by axis tricks.
star <- function(q) fifelse(q < 1e-3, "***", fifelse(q < 1e-2, "**",
                       fifelse(q < 0.05, "*", "ns")))
beta_lab <- SELECT[, .(facet_label,
  txt = sprintf("β=%+.2f%s", disease_stage_beta, star(disease_stage_q)))]
beta_lab[, facet_label := factor(facet_label, levels = SELECT$facet_label)]

p <- ggplot(ds_violin, aes(disease_stage_coarse, score_centered,
                    fill = disease_stage_coarse)) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.15,
              color = "gray55", alpha = 0.32) +
  geom_line(data = trend, aes(x = disease_stage_coarse, y = med, group = 1),
            inherit.aes = FALSE, color = "gray15", linewidth = 0.9) +
  stat_summary(data = ds, fun = median, geom = "point", shape = 21, size = 1.7,
               color = "black", fill = "white", stroke = 0.45) +
  geom_text(data = beta_lab, aes(x = 0.55, y = Inf, label = txt),
            inherit.aes = FALSE, hjust = 0, vjust = 1.35, size = 2.5,
            fontface = "plain", color = "gray10") +
  facet_wrap(~ facet_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = STAGE_COLORS, name = NULL, guide = "none") +
  labs(x = NULL,
       y = "Module score (centered)") +
  theme_masld(base_size = 9) +
  theme(axis.text.x   = element_text(size = 7.5, angle = 30, hjust = 1),
        axis.text.y   = element_text(size = 7),
        axis.title.y  = element_text(size = 8.5),
        strip.text    = element_text(size = 7.5, face = "plain", lineheight = 1.05),
        panel.spacing = unit(0.5, "lines"))

# ── PANEL CUT 2026-07-08 (user request): figs3_hotspot_cascade is permanently removed — do NOT re-enable this save. ──
message("[CUT 2026-07-08] figs3_hotspot_cascade panel removed per user request; no PDF written.")
# save_fig(p, file.path(PANEL_DIR, "figs3_hotspot_cascade.pdf"),
#          width = fig_full_width * 1.3, height = 2.3)

# Side data for caption transparency
fwrite(SELECT, file.path(DATA_DIR, "hotspot_cascade_modules.csv"))
fwrite(ds[, .(sample, module_int, score, disease_stage_coarse, dataset)],
       file.path(DATA_DIR, "hotspot_cascade_donor_scores.csv"))
cat("figs3_hotspot_cascade panel is CUT (2026-07-08) — no PDF written.\n")

# ============================================================================
# Non-parenchymal cascade — SAME method (dataset-center + median-per-stage) on
# the three disease-GAINING NPC Hotspot modules. Supports the "carrier hands off
# from hepatocytes to non-parenchymal cells" thesis: as the hepatocyte identity
# arm falls, the fibroblast/macrophage/cholangiocyte injury arms rise.
#   Fibroblasts   module 14 — activated stellate (PDGFRA), β=+0.63
#   Macrophages   module  3 — glycolysis/hypoxia (HK2),   β=+1.23
#   Cholangiocytes module 14 — EMT/ductal injury,          β=+0.38
# T cells carry no significant disease module and are omitted by design.
# Each per-ct donor_scores.tsv uses plain integer `module` (not the ct__n form
# of the pooled donor_scores_all.tsv) and SRR-style sample ids that join the
# same donor_metadata_extended.tsv used above.
# ============================================================================
NPC_FOCAL <- data.table(
  cell_type    = c("fibroblasts", "macrophages", "cholangiocytes"),
  module       = c(14L, 3L, 14L),
  module_label = c("Fib-14 · activated stellate (PDGFRA)",
                   "Mac-3 · glycolysis/hypoxia (HK2)",
                   "Chol-14 · EMT/ductal injury"))

# Reuse the SAME stage/exclusion filter already applied to `meta` above
# (exclude_stage_analysis != TRUE). Restrict to the shared 3-stage cascade.
npc_list <- lapply(seq_len(nrow(NPC_FOCAL)), function(i) {
  ct  <- NPC_FOCAL$cell_type[i]
  mod <- NPC_FOCAL$module[i]
  dsi <- fread(file.path(HS_RES, ct, "donor_scores.tsv"))   # sample, module, score
  dsi <- dsi[module == mod]
  dsi <- dsi[meta, on = "sample", nomatch = 0]              # join stage/dataset
  dsi <- dsi[disease_stage_coarse %in% STAGE_LEVELS]
  # dataset-center per (module, dataset) — identical to the hepatocyte track.
  dsi[, score_centered := score - mean(score, na.rm = TRUE),
      by = .(module, dataset)]
  dsi[, .(cell_type = ct,
          module = mod,
          module_label = NPC_FOCAL$module_label[i],
          stage = disease_stage_coarse,
          score_centered)]
})
npc <- rbindlist(npc_list)

# Median score_centered per (module, stage), long form matching the task schema.
npc_med <- npc[, .(med = median(score_centered, na.rm = TRUE)),
               by = .(cell_type, module, module_label, stage)]
setnames(npc_med, "stage", "stage")   # keep explicit
npc_med[, stage := factor(stage, levels = STAGE_LEVELS)]
setorder(npc_med, cell_type, stage)

fwrite(npc_med[, .(cell_type, module, module_label, stage, med)],
       file.path(DATA_DIR, "singlecell_disease_cascade_npc.csv"))
cat("Wrote singlecell_disease_cascade_npc.csv\n")
print(npc_med[, .(cell_type, module, stage, med = round(med, 3))])
