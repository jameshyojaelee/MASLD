#!/usr/bin/env Rscript
# figS_integration_value.R
# Supplementary figure: why the dream mega-analysis matters
# Four panels argue that per-study DEG calls are fragmented and dream recovers
# the consistent cross-cohort signal that no individual study can see.
#
#   dream_venn.pdf       — 5-petal flower (dream-validated + total counts)
#   cohort_venn.pdf      — 5-petal flower (totals only, no dream)
#   twin_chords.pdf     — twin chord: validated + cohort attribution
#   overlap_chord.pdf — corrected chord: cohort × replication tier
#   overlap_alluvial.pdf — alluvial: study → replication tier
#   concordance_volcano.pdf — % cohort-sign concordance vs dream sig
#   pairwise_lfc_grid.pdf  — 5x5 LFC scatter + Spearman + sign agreement
#
# Thresholds match fig1d: padj<0.05 & |logFC|>0.5 (primary DEG framework).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(patchwork)
  library(scales)
  library(circlize)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---- Output dir (new supplementary figure) ----
FIGS_INT_DIR <- file.path(FIG_SUPP, "figS_integration_value")
PANEL_DIR <- file.path(FIGS_INT_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Constants ----
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
COHORT_SHORT <- c(GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
                  GSE162694 = "Bril",   GSE213621 = "Chen")
PADJ_THR <- 0.05
LFC_THR  <- 0.5

# ---- Data ----
message("Loading data...")
dream <- load_dream_results()
per_study <- load_per_study_de()

# Clean Ensembl gene IDs (strip version)
dream[, gene_clean := sub("\\..*", "", gene)]
per_study[, gene_clean := sub("\\..*", "", gene)]

is_deg <- function(p, l) !is.na(p) & p < PADJ_THR & !is.na(l) & abs(l) > LFC_THR

# Dream DEG set + direction
dream[, dream_deg := is_deg(dream_padj, dream_logFC)]
dream_degs <- dream[dream_deg == TRUE, gene_clean]
n_dream <- length(dream_degs)

# Per-cohort DEG sets + direction (sign(logFC))
cohort_de <- per_study[dataset %in% COHORTS,
                       .(gene_clean, dataset, logFC, padj)]
cohort_de[, cohort_deg := is_deg(padj, logFC)]

deg_lists <- lapply(COHORTS, function(ds) {
  unique(cohort_de[dataset == ds & cohort_deg == TRUE, gene_clean])
})
names(deg_lists) <- COHORT_SHORT[COHORTS]
study_sizes <- sapply(deg_lists, length)
message(sprintf("Dream DEGs: %d", n_dream))
for (nm in names(deg_lists))
  message(sprintf("  %s DEGs: %d", nm, length(deg_lists[[nm]])))

# Per-cohort overlap with dream
overlap_with_dream <- sapply(deg_lists, function(s) length(intersect(s, dream_degs)))
message("Per-cohort overlap with dream:")
for (nm in names(overlap_with_dream))
  message(sprintf("  %s ∩ Dream = %d (%.1f%% of cohort DEGs)",
                  nm, overlap_with_dream[nm],
                  100 * overlap_with_dream[nm] / study_sizes[nm]))

# ============================================================================
# PANEL A — Classical 5-set Venn diagram (Grünbaum 5-ellipse layout)
#   + dream as a sixth set overlay (central magenta circle)
#   Every one of the 31 non-empty intersection regions is labelled with
#   total DEG count and the dream-validated subset in parentheses.
# ============================================================================
message("\n── Panel A: 5-set Venn (Grünbaum) with dream centerpiece ──")

cohort_names <- names(deg_lists)             # 5 cohorts in COHORT order
stopifnot(length(cohort_names) == 5)

# --- 1. Compute 31 intersection region counts ----------------------------
all_genes <- unique(unlist(deg_lists))
mask_mat <- sapply(deg_lists, function(s) all_genes %in% s)   # rows = genes
if (is.null(dim(mask_mat)))
  mask_mat <- matrix(mask_mat, nrow = length(all_genes))
mask_key <- apply(mask_mat, 1, function(x) paste0(as.integer(x), collapse = ""))
dream_hit <- all_genes %in% dream_degs

region_dt <- data.table(gene = all_genes,
                        region_key = mask_key,
                        in_dream = dream_hit)
region_counts <- region_dt[, .(N = .N, N_dream = sum(in_dream)), by = region_key]
region_counts <- region_counts[region_key != "00000"]

message("All ", nrow(region_counts), " non-empty intersection regions (sorted desc):")
print(region_counts[order(-N)])

# --- 2. Grünbaum 5-ellipse layout ---------------------------------------
# Classical 5-fold rotationally symmetric Venn (Grünbaum 1975; same family
# used by the VennDiagram and eulerr packages for k=5).
# Parameters chosen by exhaustive sweep so that all 31 non-empty intersection
# regions exist geometrically:
#   semi-major a = 0.40, semi-minor b = 0.15, center offset r = 0.14,
#   major axis rotated by (center_angle + 18°). The +18° tilt off the radial
#   is what creates the staircase pattern that yields 31 distinct regions.
ell_a <- 0.40    # semi-major
ell_b <- 0.15    # semi-minor
ell_r <- 0.14    # center offset from origin
ROT_OFFSET_DEG <- 18

# Center angles — i-th center at angle (90° - (i-1)*72°), so first ellipse at top
center_angles_deg <- 90 - (seq_len(5) - 1) * 72
center_angles_rad <- center_angles_deg * pi / 180
ell_rot_deg <- center_angles_deg + ROT_OFFSET_DEG
ell_rot_rad <- ell_rot_deg * pi / 180

ellipses <- data.table(
  cohort = cohort_names,
  i      = seq_len(5),
  x0     = ell_r * cos(center_angles_rad),
  y0     = ell_r * sin(center_angles_rad),
  a      = ell_a,
  b      = ell_b,
  angle  = ell_rot_rad
)

# --- 3. Point-in-ellipse helper for region centroid computation ----------
point_in_ellipse <- function(px, py, x0, y0, a, b, angle) {
  ca <- cos(-angle); sa <- sin(-angle)
  dx <- px - x0; dy <- py - y0
  xr <- dx * ca - dy * sa
  yr <- dx * sa + dy * ca
  (xr / a)^2 + (yr / b)^2 <= 1
}

# Fine grid covering the diagram, classify each grid point by which ellipses
# contain it, then take the geometric centroid per region.
n_grid <- 800
g_lim  <- 0.58
gx <- seq(-g_lim, g_lim, length.out = n_grid)
gy <- seq(-g_lim, g_lim, length.out = n_grid)
G <- as.data.table(expand.grid(x = gx, y = gy))

memb_mat <- matrix(FALSE, nrow = nrow(G), ncol = 5)
for (i in seq_len(5)) {
  memb_mat[, i] <- point_in_ellipse(G$x, G$y,
                                    ellipses$x0[i], ellipses$y0[i],
                                    ellipses$a[i], ellipses$b[i],
                                    ellipses$angle[i])
}
G[, region_key := apply(memb_mat, 1, function(x) paste0(as.integer(x), collapse = ""))]
G <- G[region_key != "00000"]
centroids <- G[, .(cx = mean(x), cy = mean(y), area = .N), by = region_key]
message(sprintf("Geometric regions found: %d (expected 31)", nrow(centroids)))

missing_geom <- setdiff(region_counts$region_key, centroids$region_key)
if (length(missing_geom) > 0)
  warning("Data regions with no matching geometry: ",
          paste(missing_geom, collapse = ", "))

# Leader positions for tiny regions (area < threshold): push label outward
# along the centroid radius and draw a thin guide line.
plot_regions <- merge(region_counts, centroids, by = "region_key", all.x = TRUE)
plot_regions[is.na(cx), `:=`(cx = 0, cy = 0, area = 0)]
plot_regions[, n_active := nchar(gsub("0", "", region_key))]

# Split labels into two layers so the dream-validated subcount can be drawn
# in teal (vs total in gray). Big regions get two-line stacked text; tiny
# regions get a compact one-line "N (N_dream)" with a small offset for the
# dream subcount.
plot_regions[, lbl_total := comma(N)]
plot_regions[, lbl_dream := sprintf("(%s)", comma(N_dream))]
plot_regions[, is_compact := N < 100]
plot_regions[N < 100, lbl_total := as.character(N)]
plot_regions[N < 100, lbl_dream := sprintf("(%d)", N_dream)]

# Detect tiny geometric regions that need a leader line. Threshold tuned so
# only the smallest regions get leaders.
area_q <- quantile(plot_regions$area[plot_regions$area > 0],
                   probs = 0.35, na.rm = TRUE)
plot_regions[, use_leader := area > 0 & area < area_q]
plot_regions[, c_r := sqrt(cx^2 + cy^2)]
plot_regions[, c_ang := atan2(cy, cx)]
# Default leader push distance
plot_regions[, leader_push := 0.12]
# 5-way intersection stays at centroid
plot_regions[region_key == "11111", `:=`(use_leader = FALSE,
                                          leader_push = 0)]
# Compute initial leader endpoints
plot_regions[, leader_x := cx + ifelse(use_leader,
                                       leader_push * cos(c_ang), 0)]
plot_regions[, leader_y := cy + ifelse(use_leader,
                                       leader_push * sin(c_ang), 0)]

# --- Anti-collision: cluster regions by angular bucket and stagger them ---
# Sort leader-using regions by angle, then within each cluster (angle delta
# < 25°) push successive labels further out on alternating radii. This
# prevents bunching where multiple tiny regions cluster along a similar
# angular direction (e.g. the "401/217" / "440/971" stacks).
ldr_idx <- which(plot_regions$use_leader)
if (length(ldr_idx) > 0) {
  ord <- ldr_idx[order(plot_regions$c_ang[ldr_idx])]
  prev_ang <- NA_real_
  push_pad <- 0  # additional outward push within a cluster
  for (k in seq_along(ord)) {
    i <- ord[k]
    th <- plot_regions$c_ang[i]
    if (!is.na(prev_ang)) {
      dtheta <- abs(prev_ang - th)
      dtheta <- min(dtheta, 2 * pi - dtheta)
      if (dtheta < (25 * pi / 180)) {
        push_pad <- push_pad + 0.06   # stagger further out
      } else {
        push_pad <- 0
      }
    }
    new_push <- plot_regions$leader_push[i] + push_pad
    plot_regions$leader_x[i] <- plot_regions$cx[i] + new_push * cos(th)
    plot_regions$leader_y[i] <- plot_regions$cy[i] + new_push * sin(th)
    prev_ang <- th
  }
}

# --- 4. Cohort labels — placed along radial direction outside each ellipse
# (radial placement reads more naturally than along the tilted major axis).
ellipses[, c_ang_rad := center_angles_rad]
ellipses[, `:=`(
  lab_x = (ell_r + ell_a + 0.06) * cos(c_ang_rad),
  lab_y = (ell_r + ell_a + 0.06) * sin(c_ang_rad)
)]

# --- 5. Color scheme — cool teal→violet from palette3 --------------------
cohort_palette <- c(
  Suppli  = "#40b499",   # palette3[1] — teal
  Hoang   = "#4aa2c2",   # palette3[3] — cyan-blue
  Govaere = "#518dc9",   # palette3[4] — medium blue
  Bril    = "#6470d0",   # palette3[5] — periwinkle
  Chen    = "#9b75d6"    # palette3[6] — soft violet
)

# --- 6. Dream sixth set ---------------------------------------------------
# Strategy: the count in parentheses in EVERY region IS the dream-validated
# subset, so dream is omnipresent. We add a small magenta circle at the
# diagram center as a visual anchor that "dream sits inside the consensus
# region" and a corner legend ("N (X) = total (dream-validated)" + dream
# headline count). This keeps the inner 11111 region's "80" label visible.
r_dream_overlay <- 0.05

extent <- 0.82
p_A <- ggplot() +
  # Five cohort ellipses (drawn first, semi-transparent)
  geom_ellipse(data = ellipses,
               aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle,
                   fill = cohort),
               color = "white", linewidth = 0.4, alpha = 0.32) +
  scale_fill_manual(values = cohort_palette, guide = "none") +
  # Leader lines for tiny regions
  geom_segment(data = plot_regions[use_leader == TRUE],
               aes(x = cx, y = cy, xend = leader_x, yend = leader_y),
               color = "gray45", linewidth = 0.22) +
  # Region count labels — total in gray, dream-validated subcount in teal
  # so dream coverage reads at a glance across every region (no misleading
  # ring at the center: dream DEGs are distributed across regions, not
  # concentrated in any single one).
  # Two-line stacked labels for big regions
  geom_text(data = plot_regions[is_compact == FALSE],
            aes(x = leader_x, y = leader_y + 0.012, label = lbl_total),
            size = 1.45, lineheight = 0.85, color = "gray15") +
  geom_text(data = plot_regions[is_compact == FALSE],
            aes(x = leader_x, y = leader_y - 0.014, label = lbl_dream),
            size = 1.30, lineheight = 0.85,
            color = masld_colors$conserved, fontface = "bold") +
  # Compact single-line labels for tiny regions: "N" then "(N_dream)"
  # placed just to the right
  geom_text(data = plot_regions[is_compact == TRUE],
            aes(x = leader_x - 0.008, y = leader_y, label = lbl_total),
            size = 1.40, color = "gray15", hjust = 1) +
  geom_text(data = plot_regions[is_compact == TRUE],
            aes(x = leader_x + 0.008, y = leader_y, label = lbl_dream),
            size = 1.25, color = masld_colors$conserved,
            fontface = "bold", hjust = 0) +
  # Cohort labels at the outer tips
  geom_text(data = ellipses,
            aes(x = lab_x, y = lab_y,
                label = sprintf("%s\n(%s)", cohort,
                                comma(study_sizes[cohort]))),
            size = 2.1, fontface = "bold", lineheight = 0.9,
            color = "gray15") +
  # Corner legend — explains the (X) notation and shows the dream total
  annotate("text",
           x = -extent + 0.02, y = -extent + 0.10,
           label = sprintf("Dream DEGs total: %s (%s in ≥1 cohort + %d dream-only)",
                           comma(n_dream),
                           comma(sum(dream_degs %in% unique(unlist(deg_lists)))),
                           n_dream - sum(dream_degs %in% unique(unlist(deg_lists)))),
           color = masld_colors$conserved, fontface = "bold",
           size = 1.95, hjust = 0) +
  coord_fixed(xlim = c(-extent, extent), ylim = c(-extent, extent),
              clip = "off") +
  theme_void() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5,
                                  margin = margin(b = 4)),
        plot.margin = margin(6, 6, 6, 6)) +
  labs(title = "5-cohort DEG overlap")

save_fig(p_A, file.path(PANEL_DIR, "dream_venn.pdf"),
         width = fig_half_width, height = fig_half_width)
message("  Saved dream_venn.pdf")

# ============================================================================
# PANEL B — Chord: replication tier × dream direction
# ============================================================================
message("\n── Panel B: chord replication tier × direction ──")

# Tier definition (sign-agnostic, matches fig1d / Panel A Venn):
#   cohort-DEG = padj<0.05 AND |logFC|>0.5  (no direction filter)
# Single unified definition so dream-only counts cross-reference with Panel A.
# Note: this gives Dream-only = 42 (under fig1d convention) instead of 35
# (fig1f's direction-aware convention) — see figure caption / docstring.

n_cohort_called <- cohort_de[cohort_deg == TRUE,
                              .(n_cohorts_called = uniqueN(dataset)),
                              by = gene_clean]

all_tested <- unique(c(dream$gene_clean, cohort_de$gene_clean))
gene_anno <- data.table(gene_clean = all_tested)
gene_anno <- merge(gene_anno, n_cohort_called, by = "gene_clean", all.x = TRUE)
gene_anno[is.na(n_cohorts_called),  n_cohorts_called  := 0L]

# Dream direction
dream_dir <- dream[, .(gene_clean,
                       dream_class = fcase(
                         dream_deg == TRUE & dream_logFC > 0, "Dream-Up",
                         dream_deg == TRUE & dream_logFC < 0, "Dream-Down",
                         default = "Dream-NS"))]
gene_anno <- merge(gene_anno, dream_dir, by = "gene_clean", all.x = TRUE)
gene_anno[is.na(dream_class), dream_class := "Dream-NS"]

# Unified tier (same definition for dream and non-dream genes)
gene_anno[, n_cohorts_sig := n_cohorts_called]
gene_anno <- gene_anno[n_cohorts_sig > 0 | dream_class != "Dream-NS"]

tier_labels <- c("0" = "Dream-only", "1" = "1 cohort", "2" = "2 cohorts",
                 "3" = "3 cohorts", "4" = "4 cohorts", "5" = "All 5")
gene_anno[, tier := tier_labels[as.character(n_cohorts_sig)]]

chord_dt <- gene_anno[, .N, by = .(tier, dream_class)]
chord_dt <- chord_dt[N > 0]
message("Chord matrix (all):")
print(dcast(chord_dt, tier ~ dream_class, value.var = "N", fill = 0))

# Save the Dream-NS counts as a sidebar (study-private noise per tier)
tier_order <- c("Dream-only", "1 cohort", "2 cohorts", "3 cohorts", "4 cohorts", "All 5")
ns_per_tier <- chord_dt[dream_class == "Dream-NS",
                        .(tier = factor(tier, levels = tier_order), N)]
ns_per_tier <- ns_per_tier[order(tier)]

# Chord shows only dream-validated genes (Up/Down) — noise filtered out
chord_validated <- chord_dt[dream_class %in% c("Dream-Up", "Dream-Down")]
mat <- dcast(chord_validated, tier ~ dream_class, value.var = "N", fill = 0)
mat_m <- as.matrix(mat[, -1, with = FALSE])
rownames(mat_m) <- mat$tier

dream_order <- c("Dream-Up", "Dream-Down")
mat_m <- mat_m[intersect(tier_order, rownames(mat_m)),
               intersect(dream_order, colnames(mat_m)), drop = FALSE]
message("Chord matrix (dream-validated only):")
print(mat_m)

# Colors — palette ramps from light (1 cohort) → dark (All 5) so replication
# depth reads at a glance
tier_cols <- c("Dream-only"  = "#C2185B",
               "1 cohort"    = "#F8BBD0",
               "2 cohorts"   = "#CE93D8",
               "3 cohorts"   = "#9575CD",
               "4 cohorts"   = "#5C6BC0",
               "All 5"       = "#1A237E")
dir_cols  <- c("Dream-Up"   = fig1_colors$up,
               "Dream-Down" = fig1_colors$down)
sector_cols <- c(tier_cols[rownames(mat_m)], dir_cols[colnames(mat_m)])



# ----------------------------------------------------------------------------
# B variant 3: twin chords — top = tier → Dream-Up/Down, bottom = tier → 5 cohorts
# Bottom chord shows which cohort each tier's rejected hits come from.
# ----------------------------------------------------------------------------
message("\n── Panel B v3: twin chords (validated above, cohort attribution below) ──")

# Top matrix = same as the canonical (Up/Down only)
mat_top <- mat_m   # validated chord matrix from earlier

# Bottom matrix: for each tier, # rejected hits attributable to each cohort
# (a hit is "in" a cohort if that cohort individually called it at padj<0.05).
# Genes can be in multiple cohorts → we count per (tier × cohort) cell.
rejected_genes <- gene_anno[dream_class == "Dream-NS",
                            .(gene_clean, tier)]
ps_called_long <- cohort_de[cohort_deg == TRUE,
                              .(gene_clean, dataset)]
bot <- merge(rejected_genes, ps_called_long, by = "gene_clean",
             allow.cartesian = TRUE)
bot[, cohort := COHORT_SHORT[dataset]]
bot_dt <- bot[, .N, by = .(tier, cohort)]
mat_bot <- dcast(bot_dt, tier ~ cohort, value.var = "N", fill = 0)
mat_bot_m <- as.matrix(mat_bot[, -1, with = FALSE])
rownames(mat_bot_m) <- mat_bot$tier
mat_bot_m <- mat_bot_m[intersect(tier_order, rownames(mat_bot_m)),
                        intersect(unname(COHORT_SHORT), colnames(mat_bot_m)),
                        drop = FALSE]

cohort_palette <- c(Suppli = "#518dc9", Hoang = "#9b75d6",
                    Govaere = "#d980dc", Bril = "#ea868d", Chen = "#f7bf87")
sector_cols_bot <- c(tier_cols[rownames(mat_bot_m)],
                     cohort_palette[colnames(mat_bot_m)])
# Wider gaps between cohorts so small cohort (Hoang ~3K) labels don't crowd
gap_dirs_bot <- setNames(rep(8, ncol(mat_bot_m) - 1), head(colnames(mat_bot_m), -1))
gap_dirs_bot <- c(gap_dirs_bot, setNames(26, tail(colnames(mat_bot_m), 1)))

pdf_v3 <- file.path(PANEL_DIR, "twin_chords.pdf")
cairo_pdf(pdf_v3, width = fig_half_width + 0.4, height = fig_half_width * 2 + 0.4)
layout(matrix(c(1, 2), nrow = 2), heights = c(1, 1))

# --- Top chord: validated (same as canonical B but wider bottom gap) ---
par(mar = c(1, 3.5, 3.5, 3.5))
circos.clear()
gap_tiers_top <- c(`Dream-only` = 8, `1 cohort` = 3, `2 cohorts` = 3,
                   `3 cohorts` = 8, `4 cohorts` = 12, `All 5` = 30)
gap_dirs_top  <- c(`Dream-Up` = 6, `Dream-Down` = 24)
circos.par(start.degree = 90,
           gap.after = c(gap_tiers_top[rownames(mat_top)],
                         gap_dirs_top[colnames(mat_top)]))
chordDiagram(
  mat_top, grid.col = sector_cols,
  transparency = 0.35, annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08),
  directional = -1, direction.type = "diffHeight",
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrackHeight = c(0.04, 0.04)
)
sector_totals_top <- c(rowSums(mat_top), colSums(mat_top))
circos.trackPlotRegion(
  track.index = 1, panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
    n_sec <- sector_totals_top[sector]
    label <- sprintf("%s (%s)", sector, format(n_sec, big.mark = ","))
    circos.text(mean(xlim), ylim[1] + 2.0, label,
                facing = "bending.outside", niceFacing = TRUE,
                cex = 0.45)
  }, bg.border = NA)
title(main = "Dream DEGs validation",
      cex.main = 0.85, line = 1.5)
circos.clear()

# --- Bottom chord: rejected by tier → cohort of origin ---
par(mar = c(1, 3.5, 3.5, 3.5))
gap_tiers_bot <- c(`Dream-only` = 8, `1 cohort` = 3, `2 cohorts` = 3,
                   `3 cohorts` = 8, `4 cohorts` = 20, `All 5` = 30)
circos.par(start.degree = 90,
           gap.after = c(gap_tiers_bot[rownames(mat_bot_m)],
                         gap_dirs_bot[colnames(mat_bot_m)]))
chordDiagram(
  mat_bot_m, grid.col = sector_cols_bot,
  transparency = 0.4, annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08),
  directional = -1, direction.type = "diffHeight",
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrackHeight = c(0.04, 0.04)
)
sector_totals_bot <- c(rowSums(mat_bot_m), colSums(mat_bot_m))
circos.trackPlotRegion(
  track.index = 1, panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
    n_sec <- sector_totals_bot[sector]
    label <- sprintf("%s (%s)", sector, format(n_sec, big.mark = ","))
    circos.text(mean(xlim), ylim[1] + 2.0, label,
                facing = "bending.outside", niceFacing = TRUE,
                cex = 0.45)
  }, bg.border = NA)
title(main = "Dream-NS validation",
      cex.main = 0.85, line = 1.5)
circos.clear()
dev.off()
message("  Saved twin_chords.pdf")

# ----------------------------------------------------------------------------
# B variant 4: cohort × replication-tier chord (corrected — no double-counting)
# Each DEG is assigned to a tier bucket (0–4 other cohorts) exactly once per
# cohort it belongs to. Arc length = true per-cohort DEG count.
# Pairwise intersection matrices are NOT used here because a gene in k cohorts
# would contribute to C(k,2) pairwise entries, inflating arc lengths.
# ----------------------------------------------------------------------------
message("\n── Panel B v4: per-cohort replication-tier chord ──")

# For each cohort, for each DEG, count how many OTHER cohorts also have it.
all_genes_v4 <- unique(unlist(deg_lists))
n_coh_per_gene_v4 <- setNames(
  rowSums(sapply(deg_lists, function(s) all_genes_v4 %in% s)),
  all_genes_v4)

# Tier labels = total cohort count (no repeated "others")
tier_names_v4 <- c("Unique to study", "2 cohorts", "3 cohorts", "4 cohorts", "5 cohorts")
mat_tier_v4 <- matrix(0L, length(deg_lists), 5,
                      dimnames = list(names(deg_lists), tier_names_v4))
for (nm in names(deg_lists)) {
  n_others <- n_coh_per_gene_v4[deg_lists[[nm]]] - 1
  for (k in 0:4) mat_tier_v4[nm, k + 1] <- sum(n_others == k)
}
stopifnot(all(rowSums(mat_tier_v4) == sapply(deg_lists, length)))

cohort_palette_v4 <- c(Suppli  = "#518dc9", Hoang   = "#9b75d6",
                       Govaere = "#d980dc", Bril    = "#ea868d", Chen    = "#f7bf87")
tier_palette_v4   <- c(`Unique to study` = "#BDBDBD", `2 cohorts` = "#CE93D8",
                       `3 cohorts` = "#9575CD", `4 cohorts` = "#5C6BC0",
                       `5 cohorts` = "#1A237E")
n_no_overlap_v4 <- as.integer(table(n_coh_per_gene_v4)["1"])

tier_order_right_v4  <- c("5 cohorts", "4 cohorts", "3 cohorts", "2 cohorts", "Unique to study")
cohort_order_left_v4 <- names(deg_lists)
mat_tier_v4_T <- t(mat_tier_v4)[tier_order_right_v4, cohort_order_left_v4]
sector_cols_v4 <- c(tier_palette_v4[tier_order_right_v4],
                    cohort_palette_v4[cohort_order_left_v4])

gap_v4 <- c(setNames(c(12, 12, 10, 14), tier_order_right_v4[-5]),
            setNames(30, tier_order_right_v4[5]),
            setNames(c(4, 4, 4, 4),    cohort_order_left_v4[-5]),
            setNames(30, cohort_order_left_v4[5]))

draw_label_v4 <- function(n_uniq, r_near = 1.08, r_far = 1.12, clash_deg = 22, cex = 0.70,
                           x_nudge = c(`5 cohorts` = -0.22, `4 cohorts` = -0.06)) {
  all_sec <- c(tier_order_right_v4, cohort_order_left_v4)
  si <- lapply(all_sec, function(s) {
    set.current.cell(sector.index = s, track.index = 1)
    th <- circlize(mean(get.cell.meta.data("xlim")), 1,
                   sector.index = s, track.index = 1)[1, "theta"]
    list(s = s, th = (th + 360) %% 360)
  })
  si <- si[order(sapply(si, `[[`, "th"), decreasing = TRUE)]
  rails <- rep("near", length(si))
  for (i in seq_along(si)[-1]) {
    d <- min(abs(si[[i]]$th - si[[i-1]]$th), 360 - abs(si[[i]]$th - si[[i-1]]$th))
    if (d < clash_deg) rails[i] <- if (rails[i-1] == "near") "far" else "near"
  }
  for (i in seq_along(si)) {
    s <- si[[i]]$s; th <- si[[i]]$th
    r  <- if (rails[i] == "near") r_near else r_far
    tr <- th * pi / 180
    lbl <- if (s == "Unique to study")
      sprintf("Unique to study\n(%s)", format(n_uniq, big.mark = ","))
    else s
    ah <- if (th > 90 && th < 270) 1 else 0
    nx <- if (s %in% names(x_nudge)) x_nudge[s] else 0
    text(r * cos(tr) + (if (ah == 0) 0.02 else -0.02) + nx, r * sin(tr),
         lbl, cex = cex, adj = c(ah, 0.5), col = "black", xpd = TRUE)
  }
}

pdf_v4 <- file.path(PANEL_DIR, "overlap_chord.pdf")
cairo_pdf(pdf_v4, width = fig_half_width + 1.4, height = fig_half_width + 1.4)
par(mar = c(1, 1, 3, 1))
circos.clear()
circos.par(start.degree = 90, canvas.xlim = c(-1.35, 1.35),
           canvas.ylim = c(-1.35, 1.35), gap.after = gap_v4)
chordDiagram(
  mat_tier_v4_T, grid.col = sector_cols_v4, transparency = 0.35,
  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05),
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrackHeight = c(0.03, 0.03)
)
draw_label_v4(n_no_overlap_v4)
title(main = "Per-cohort DEG replication tier",
      cex.main = 0.85, line = 1.5)
circos.clear()
dev.off()
message("  Saved overlap_chord.pdf")

# ----------------------------------------------------------------------------
# B alluvial: same data as B_v4, alluvial layout (Liang et al lncRNA style)
# Studies on the left, replication tier on the right (4 others → Unique to study,
# top to bottom). Labels outside the strata, soft pastel palette, no numbers
# inside the alluvial to keep the visual clean.
# ----------------------------------------------------------------------------
message("\n── Panel B alluvial: study × replication tier ──")

suppressPackageStartupMessages(library(ggalluvial))

df_alluv <- as.data.table(as.table(mat_tier_v4))
setnames(df_alluv, c("Study", "Tier", "Freq"))
df_alluv[, Freq := as.numeric(Freq)]

cohort_order_top_down <- names(sort(sapply(deg_lists, length), decreasing = TRUE))
tier_order_top_down   <- c("5 cohorts", "4 cohorts", "3 cohorts",
                           "2 cohorts", "Unique to study")
df_alluv[, Study := factor(Study, levels = cohort_order_top_down)]
df_alluv[, Tier  := factor(Tier,  levels = tier_order_top_down)]

# Soft pastel palette (Liang et al lncRNA-style)
cohort_pal_pastel <- c(Suppli = "#a8c7e5", Hoang   = "#c9b8e8",
                       Govaere = "#e8bde9", Bril    = "#f3bec1",
                       Chen    = "#fbdebd")
tier_pal_pastel   <- c(`Unique to study` = "#D5D5D5",
                       `2 cohorts`         = "#e3bee8",
                       `3 cohorts`        = "#beace0",
                       `4 cohorts`        = "#a3acd3",
                       `5 cohorts`        = "#7f87b6")

# Stratum mid-y positions for outside labels
left_dt  <- df_alluv[, .(Freq = sum(Freq)), by = Study]
setorder(left_dt, Study)
total_y  <- sum(left_dt$Freq)
left_dt[, y_top := total_y - cumsum(c(0, head(Freq, -1)))]
left_dt[, y_bot := y_top - Freq]
left_dt[, y_mid := (y_top + y_bot) / 2]

right_dt <- df_alluv[, .(Freq = sum(Freq)), by = Tier]
setorder(right_dt, Tier)
right_dt[, y_top := total_y - cumsum(c(0, head(Freq, -1)))]
right_dt[, y_bot := y_top - Freq]
right_dt[, y_mid := (y_top + y_bot) / 2]

p_alluv <- ggplot(df_alluv, aes(y = Freq, axis1 = Study, axis2 = Tier)) +
  geom_alluvium(aes(fill = Study), width = 1/10, alpha = 0.55,
                knot.pos = 0.42, curve_type = "sigmoid", linewidth = 0) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/10,
               color = "white", linewidth = 0.45) +
  annotate("text", x = 1 - 0.085, y = left_dt$y_mid,
           label = as.character(left_dt$Study), hjust = 1,
           size = 2.2, fontface = "bold", color = "gray15") +
  annotate("text", x = 2 + 0.085, y = right_dt$y_mid,
           label = ifelse(as.character(right_dt$Tier) == "Unique to study",
                          sprintf("Unique to study\n(%s)", format(n_no_overlap_v4, big.mark = ",")),
                          as.character(right_dt$Tier)), hjust = 0,
           size = 2.2, fontface = "bold", color = "gray15") +
  scale_x_discrete(limits = c("Study", "Overlap"),
                   expand = c(0.30, 0.30), position = "top") +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  scale_fill_manual(values = c(cohort_pal_pastel, tier_pal_pastel),
                    guide = "none") +
  labs(title = "Per-cohort DEG overlap",
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title      = element_text(size = 8.5, face = "bold",
                                       hjust = 0.5, margin = margin(b = 6)),
        axis.text.x.top = element_text(size = 7.5, face = "bold",
                                       color = "gray15",
                                       margin = margin(b = 4)),
        plot.margin     = margin(8, 14, 6, 14))

save_fig(p_alluv,
         file.path(PANEL_DIR, "overlap_alluvial.pdf"),
         width = fig_half_width + 1.6, height = fig_half_width + 0.5)
message("  Saved overlap_alluvial.pdf")

# ============================================================================
# PANEL C — Direction-concordance volcano (Idea 3)
# ============================================================================
message("\n── Panel C: direction-concordance volcano ──")

# Per gene: sign of per-cohort logFC vs dream sign; count concordant cohorts
# (only over cohorts where the gene was tested).
ps_wide <- dcast(cohort_de, gene_clean ~ dataset,
                 value.var = c("logFC", "padj"),
                 fill = NA_real_)
ps_wide <- merge(ps_wide,
                 dream[, .(gene_clean, dream_logFC, dream_padj, dream_deg)],
                 by = "gene_clean", all = TRUE)

lfc_cols <- paste0("logFC_", COHORTS)
padj_cols <- paste0("padj_", COHORTS)

# Direction concordance: fraction of cohorts where sign(per-cohort LFC) == sign(dream LFC)
ps_wide[, sign_dream := sign(dream_logFC)]
sign_mat <- as.matrix(ps_wide[, ..lfc_cols])
sign_mat <- sign(sign_mat)
n_tested <- rowSums(!is.na(sign_mat))
n_concord <- rowSums(sign_mat == ps_wide$sign_dream, na.rm = TRUE)
ps_wide[, n_tested_cohorts := n_tested]
ps_wide[, n_concordant     := n_concord]
ps_wide[, pct_concordant   := ifelse(n_tested_cohorts > 0,
                                     100 * n_concordant / n_tested_cohorts, NA_real_)]

# Individual-cohort DEG count
padj_mat <- as.matrix(ps_wide[, ..padj_cols])
lfc_mat  <- as.matrix(ps_wide[, ..lfc_cols])
ind_sig <- (padj_mat < PADJ_THR) & (abs(lfc_mat) > LFC_THR)
ind_sig[is.na(ind_sig)] <- FALSE
ps_wide[, n_cohorts_indiv_sig := rowSums(ind_sig)]

# Universe: genes called DEG in dream OR in ≥1 cohort
volcano_dt <- ps_wide[!is.na(dream_padj) & n_tested_cohorts >= 3 &
                      (dream_deg == TRUE | n_cohorts_indiv_sig > 0)]
volcano_dt[, neglog10_p := -log10(pmax(dream_padj, 1e-30))]
volcano_dt[, indiv_cat := factor(pmin(n_cohorts_indiv_sig, 5),
                                  levels = 0:5,
                                  labels = c("0 (dream-rescued)", "1", "2", "3", "4", "5"))]

# Rescue quadrant: dream-significant, ≥80% concordant, individually-significant in 0 or 1 cohorts
rescue_dt <- volcano_dt[dream_deg == TRUE & pct_concordant >= 80 &
                        n_cohorts_indiv_sig <= 1]
n_rescue <- nrow(rescue_dt)
message(sprintf("  Rescue quadrant (dream-sig, ≥80%% concordant, ≤1 individual): %d genes",
                n_rescue))

# Hex/raster scatter
indiv_palette <- c("0 (dream-rescued)" = "#C2185B",
                   "1" = "#E91E63", "2" = "#7B1FA2",
                   "3" = "#3F51B5", "4" = "#0288D1", "5" = "#00695C")

set.seed(42)
volcano_dt[, x_jit := pct_concordant + runif(.N, -1.5, 1.5)]

p_C <- ggplot(volcano_dt,
              aes(x = x_jit, y = neglog10_p, color = indiv_cat)) +
  rasterize_layer(geom_point(size = 0.25, alpha = 0.55, shape = 16)) +
  geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed",
             color = "gray30", linewidth = 0.3) +
  geom_vline(xintercept = 80, linetype = "dashed",
             color = "gray30", linewidth = 0.3) +
  annotate("rect", xmin = 80, xmax = 102, ymin = -log10(PADJ_THR), ymax = Inf,
           fill = NA, color = "#880E4F", linewidth = 0.5, linetype = "solid") +
  annotate("text", x = 81, y = max(volcano_dt$neglog10_p) * 0.97,
           label = sprintf("Dream rescued: %s", comma(n_rescue)),
           hjust = 0, vjust = 1, size = 2.2, color = "#880E4F", fontface = "bold") +
  scale_color_manual(values = indiv_palette,
                     name = "Individually\nDEG (cohorts)") +
  scale_x_continuous(limits = c(-2, 102), breaks = c(0, 20, 40, 60, 80, 100),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "% of cohorts with same LFC sign as dream",
       y = expression(-log[10]("dream padj")),
       title = "Dream rescues consistent-direction genes individual cohorts miss") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        legend.position = "right")

save_fig(p_C, file.path(PANEL_DIR, "concordance_volcano.pdf"),
         width = fig_col_width, height = 3.2)
message("  Saved concordance_volcano.pdf")

# ============================================================================
# PANEL D — 5×5 pairwise LFC grid
# ============================================================================
message("\n── Panel D: 5×5 pairwise LFC grid ──")

# Use per-cohort LFC matrix from earlier (lfc_mat is over all genes; merge dream sig)
pair_dt <- ps_wide[, c("gene_clean", "dream_deg", lfc_cols), with = FALSE]
pair_dt[is.na(dream_deg), dream_deg := FALSE]
setnames(pair_dt, lfc_cols, COHORT_SHORT[COHORTS])

# Add per-cohort DEG flags (padj < PADJ_THR & |logFC| > LFC_THR) for filtering in make_upper
deg_flags <- dcast(cohort_de[dataset %in% COHORTS], gene_clean ~ dataset,
                   value.var = "cohort_deg", fun.aggregate = any, fill = FALSE)
# rename dataset IDs → short names with "deg_" prefix
setnames(deg_flags, COHORTS, paste0("deg_", COHORT_SHORT[COHORTS]))
pair_dt <- merge(pair_dt, deg_flags, by = "gene_clean", all.x = TRUE)
deg_cols <- paste0("deg_", COHORT_SHORT[COHORTS])
for (col in deg_cols) pair_dt[is.na(get(col)), (col) := FALSE]

cohort_names <- unname(COHORT_SHORT[COHORTS])

# Diagonal counts
diag_dt <- data.table(cohort = cohort_names, n_degs = study_sizes)

# Helper: pairwise panel
make_lower <- function(coh_x, coh_y) {
  dat <- pair_dt[!is.na(get(coh_x)) & !is.na(get(coh_y)),
                 .(x = get(coh_x), y = get(coh_y), dream_deg)]
  if (nrow(dat) > 30000) dat <- dat[sample(.N, 30000)]
  lim <- range(c(dat$x, dat$y), na.rm = TRUE)
  ggplot(dat, aes(x = x, y = y)) +
    rasterize_layer(geom_point(aes(color = dream_deg, alpha = dream_deg),
                               size = 0.15, shape = 16)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray60", linewidth = 0.25) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "gray80") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "gray80") +
    scale_color_manual(values = c(`TRUE` = masld_colors$conserved,
                                  `FALSE` = "gray75"), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 0.90, `FALSE` = 0.30), guide = "none") +
    coord_fixed(xlim = lim, ylim = lim) +
    theme_void() +
    theme(panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.25),
          plot.margin = margin(1, 1, 1, 1))
}

make_upper <- function(coh_x, coh_y) {
  # Filter to genes called DEG in at least one of the two cohorts (union)
  # so ρ and % same dir reflect replication among biologically relevant genes
  deg_x <- paste0("deg_", coh_x)
  deg_y <- paste0("deg_", coh_y)
  dat <- pair_dt[!is.na(get(coh_x)) & !is.na(get(coh_y)) &
                   (get(deg_x) | get(deg_y))]
  rho <- suppressWarnings(cor(dat[[coh_x]], dat[[coh_y]],
                              method = "spearman", use = "pairwise.complete.obs"))
  pct_same <- 100 * mean(sign(dat[[coh_x]]) == sign(dat[[coh_y]]),
                         na.rm = TRUE)
  n_genes <- nrow(dat)
  # Tile background: rho mapped to white (0) → teal (#00695C at 1)
  # Clamp rho to [0, 1] for fill (all pairs have rho >= 0 in practice)
  rho_clamped <- max(0, min(1, rho))
  bg_col  <- colorRamp(c("white", "#B2DFDB", "#00695C"))(rho_clamped)
  bg_fill <- rgb(bg_col[1], bg_col[2], bg_col[3], maxColorValue = 255)
  txt_col <- if (rho_clamped > 0.55) "white" else "gray15"
  sub_col <- if (rho_clamped > 0.55) "#E0F2F1" else "gray40"
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = bg_fill, color = NA) +
    annotate("text", x = 0.5, y = 0.68,
             label = sprintf("ρ = %.2f", rho),
             size = 3.2, color = txt_col, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.42,
             label = sprintf("%.0f%% same dir", pct_same),
             size = 3.2, color = sub_col) +
    annotate("text", x = 0.5, y = 0.17,
             label = sprintf("n = %s", comma(n_genes)),
             size = 3.2, color = sub_col) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.25),
          plot.margin = margin(1, 1, 1, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
}

make_diag <- function(coh) {
  n_d <- diag_dt[cohort == coh, n_degs]
  ggplot() +
    annotate("text", x = 0.5, y = 0.62, label = coh,
             size = 3.2, color = "gray10", fontface = "bold") +
    annotate("text", x = 0.5, y = 0.35,
             label = sprintf("%s DEGs", comma(n_d)),
             size = 3.2, color = masld_colors$conserved) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.4),
          plot.margin = margin(1, 1, 1, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
}

# Assemble 5x5 patchwork
panels <- list()
for (i in seq_along(cohort_names)) {
  for (j in seq_along(cohort_names)) {
    if (i == j) {
      panels[[length(panels) + 1]] <- make_diag(cohort_names[i])
    } else if (i > j) {
      # Lower triangle (row > col): scatter, y = cohort i, x = cohort j
      panels[[length(panels) + 1]] <- make_lower(cohort_names[j], cohort_names[i])
    } else {
      # Upper triangle (row < col): stats summary
      panels[[length(panels) + 1]] <- make_upper(cohort_names[j], cohort_names[i])
    }
  }
}

# Row/column labels via patchwork::wrap_plots
grid_p <- wrap_plots(panels, nrow = 5, ncol = 5)

# Add a header and footer caption via patchwork title
grid_p_titled <- grid_p +
  plot_annotation(
    title = "Pairwise per-cohort LFC concordance · ρ and % same dir computed among per-study DEGs (union)",
    theme = theme(plot.title = element_text(size = 8, face = "bold"))
  )

save_fig(grid_p_titled, file.path(PANEL_DIR, "pairwise_lfc_grid.pdf"),
         width = fig_full_width, height = fig_full_width)
message("  Saved pairwise_lfc_grid.pdf")

message("\nAll 4 panels saved under: ", PANEL_DIR)
