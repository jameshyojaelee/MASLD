#!/usr/bin/env Rscript
#
# *** lfsr VARIANT (auto-derived from the raw-LFC script) ***
# DEG set = ashr lfsr < 0.05 (local false-sign-rate; sign-confidence); effect = shrunk_logFC.
# Output -> figures/supplementary/figS_methods_validation/lfc_sensitivity/patient_concordance_lfsr/
# KEY MESSAGE: Dream DEGs are concordantly dysregulated in the majority of
# individual MASLD patients; patient-level consistency improves with higher
# LFC cutoffs, confirming robustness of the population-level signature.
#
# Output (4 individual panels — flat in figS_lfc_sensitivity/):
#   A_patient_concordance_violin.pdf   — per-patient % DEGs in correct direction, violin per cutoff
#   B_pct_patients_threshold_lines.pdf — % patients exceeding 50/60/75/90% concordance
#   C_deg_counts_lollipop.pdf          — N DEGs (up + down) at each LFC cutoff
#   D_spearman_rho_violin.pdf          — Spearman rho (dream LFC vs patient LFC) per cutoff
#
# Data:
#   dream_results_ashr.csv  — STAR-based dream DEGs
#   patient_lfc_matrix.csv.gz    — per-patient log2FC relative to within-dataset control mean

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
  library(patchwork)
  library(ggdist)
  library(edgeR)
  library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# Output directory — all panels from this script live in the patient_concordance/ subfolder
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/lfc_sensitivity/patient_concordance_lfsr")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

INT_RESULTS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

# ── Publication theme (matches existing panels in this project) ──────────────
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(colour = "black", linewidth = 0.3),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9),
    plot.margin        = margin(8, 10, 8, 8)
  )

theme_set(theme_pub)

# ── Colors ───────────────────────────────────────────────────────────────────
col_up       <- masld_colors$up       # "#C9265E" magenta
col_down     <- masld_colors$down     # "#1565C0" deep blue
col_combined <- "#00695C"             # teal for combined concordance

# Threshold line palette (colorblind-safe, increasing saturation)
thresh_colors <- c("50%" = "#9E9E9E", "60%" = "#FB8C00",
                   "75%" = "#C9265E", "90%" = "#4527A0")

# Dataset colors (Set2 — colorblind-safe)
ds_palette <- c(
  GSE126848 = "#66C2A5", GSE130970 = "#FC8D62",
  GSE135251 = "#8DA0CB", GSE162694 = "#E78AC3",
  GSE213621 = "#A6D854"
)

# ── Load data ────────────────────────────────────────────────────────────────
message("Loading dream STAR results...")
dream <- fread(file.path(INT_RESULTS, "dream_results_ashr.csv"))
setnames(dream, "logFC", "dream_logFC", skip_absent = FALSE)
setnames(dream, "padj",  "dream_padj",  skip_absent = FALSE)

# === lfsr VARIANT: alias ashr columns so every downstream "dream_padj < 0.05" gate
# becomes lfsr < 0.05 and every "dream_logFC" effect becomes the ashr-shrunk posterior,
# WITHOUT editing the analytical body. (sign is preserved by ashr; magnitude shrinks.) ===
stopifnot(all(c("shrunk_logFC","lfsr") %in% names(dream)))
dream[, dream_logFC := shrunk_logFC]
dream[, dream_padj  := lfsr]

message("Loading patient LFC matrix (34k genes × 691 patients)...")
lfc_mat <- fread(file.path(INT_RESULTS, "patient_lfc_matrix.csv.gz"))
setkey(lfc_mat, gene)

# Sample IDs (all columns except 'gene')
sample_ids <- setdiff(colnames(lfc_mat), "gene")
n_patients <- length(sample_ids)
message(sprintf("  %d genes × %d patients", nrow(lfc_mat), n_patients))

# Load metadata for dataset coloring
meta <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))
meta <- meta[sample_id %in% sample_ids]
sid_to_ds <- setNames(meta$dataset, meta$sample_id)

# Align dream and LFC matrix on gene IDs
common_genes <- intersect(dream$gene, lfc_mat$gene)
message(sprintf("  Genes in both: %d", length(common_genes)))
dream_sub <- dream[gene %in% common_genes]
lfc_sub   <- lfc_mat[gene %in% common_genes]
setkey(lfc_sub, gene)

# Convert LFC matrix to a standard matrix (rows = genes, cols = samples)
# Order rows to match dream_sub
lfc_sub <- lfc_sub[dream_sub$gene]
lfc_m <- as.matrix(lfc_sub[, sample_ids, with = FALSE])  # genes × patients
rownames(lfc_m) <- dream_sub$gene

# ── Sweep over LFC cutoffs ───────────────────────────────────────────────────
lfc_cutoffs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0)

message("Computing per-patient concordance across LFC cutoffs...")

patient_results <- rbindlist(lapply(lfc_cutoffs, function(cut) {
  up_idx <- which(dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC > cut)
  dn_idx <- which(dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -cut)
  n_up <- length(up_idx)
  n_dn <- length(dn_idx)
  n_total <- n_up + n_dn

  if (n_total == 0) return(NULL)

  # Per-patient concordance fractions
  # Up: fraction of up-DEGs where patient LFC > 0
  pct_up <- if (n_up > 0)
    colSums(lfc_m[up_idx, , drop = FALSE] > 0) / n_up * 100
  else
    rep(NA_real_, n_patients)

  # Down: fraction of down-DEGs where patient LFC < 0
  pct_dn <- if (n_dn > 0)
    colSums(lfc_m[dn_idx, , drop = FALSE] < 0) / n_dn * 100
  else
    rep(NA_real_, n_patients)

  # Combined: all DEGs in the correct direction
  n_concordant <- numeric(n_patients)
  if (n_up > 0) n_concordant <- n_concordant + colSums(lfc_m[up_idx, , drop = FALSE] > 0)
  if (n_dn > 0) n_concordant <- n_concordant + colSums(lfc_m[dn_idx, , drop = FALSE] < 0)
  pct_combined <- n_concordant / n_total * 100

  # Spearman rho per patient between dream_logFC and patient LFC (for this gene set)
  all_idx <- c(up_idx, dn_idx)
  dream_lfc_sub <- dream_sub$dream_logFC[all_idx]
  pat_lfc_sub   <- lfc_m[all_idx, , drop = FALSE]
  # rank dream LFC once (same for all patients)
  rk_dream <- rank(dream_lfc_sub)
  # Vectorised Spearman via Pearson on ranks
  rk_pat <- apply(pat_lfc_sub, 2, rank)  # genes × patients (rank per col)
  rho <- cor(rk_dream, rk_pat, method = "pearson")  # 1 × n_patients

  dt <- data.table(
    sample_id  = sample_ids,
    lfc_cutoff = cut,
    n_up_degs  = n_up,
    n_dn_degs  = n_dn,
    n_degs     = n_total,
    pct_up     = pct_up,
    pct_dn     = pct_dn,
    pct_combined = pct_combined,
    spearman_rho = as.numeric(rho)
  )
  dt[, dataset := sid_to_ds[sample_id]]
  dt
}))

# Cutoff label (x-axis): show as "|LFC| > X"
patient_results[, cutoff_label := factor(
  sprintf("%.1f", lfc_cutoff),
  levels = sprintf("%.1f", lfc_cutoffs)
)]

# DEG count table (one row per cutoff)
deg_counts <- patient_results[, .(
  n_up   = n_up_degs[1],
  n_dn   = n_dn_degs[1],
  n_total = n_degs[1]
), by = .(lfc_cutoff, cutoff_label)]

message(sprintf("  Patients: %d  |  Cutoffs: %d", n_patients,
                length(unique(patient_results$lfc_cutoff))))
message("DEG counts at each cutoff:")
print(deg_counts[, .(lfc_cutoff, n_up, n_dn, n_total)])

# ── Panel A: Violin — per-patient combined concordance at each LFC cutoff ─────
# One half-eye per cutoff, combined (up + down) concordance only.

pA <- ggplot(patient_results,
             aes(x = cutoff_label, y = pct_combined)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_hline(yintercept = 75, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  stat_halfeye(
    fill = col_combined, color = col_combined, alpha = 0.65,
    adjust = 0.8, width = 0.95,
    point_interval = median_qi, .width = c(0.5, 0.9),
    point_size = 1.0, interval_size = 0.6,
    slab_linewidth = 0.2
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = paste0(c(0, 25, 50, 75, 100), "%")) +
  labs(
    x     = expression("|shrunk log"[2]*"FC| cutoff  (dream DEGs, lfsr < 0.05)"),
    y     = "% of DEGs concordant per patient",
    title = "A  Patient-level concordance with population DEG direction",
    subtitle = sprintf("n = %d patients (5 cohorts with healthy controls)", n_patients)
  ) +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(OUTDIR, "patient_concordance_violins.pdf"),
       pA, width = 9, height = 4, device = cairo_pdf)
message("Saved: A_patient_concordance_violin.pdf")

# ── Panels B+C: Dual-panel — N DEGs (top) + % patients ≥75% concordant (bottom)
# Shares x-axis to show the tradeoff directly: stricter LFC → fewer genes but
# more patients satisfy the signature at the individual level.

ref_line <- geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey55",
                       linewidth = 0.35)
x_scale  <- scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0))

# ── Helper: build BC dual-panel (shared across standard + strict05 versions) ──
MAIN_THRESH  <- c(50, 75, 90)
INTER_THRESH <- c(55, 60, 65, 70, 80, 85)

make_bc_panels <- function(pat_results, pct_col, title_tag, file_suffix) {
  # Compute % patients exceeding each threshold
  all_thresholds <- sort(c(MAIN_THRESH, INTER_THRESH))
  pct_vals <- pat_results[!is.na(get(pct_col)), {
    row_vals <- lapply(all_thresholds, function(thr)
      mean(get(pct_col) >= thr, na.rm = TRUE) * 100)
    setNames(as.list(row_vals), paste0("p", all_thresholds))
  }, by = .(lfc_cutoff)]

  bc_wide_loc <- merge(
    pat_results[, .(n_total = n_degs[1]), by = lfc_cutoff],
    pct_vals, by = "lfc_cutoff")

  # Main 3 lines — colored, with points and end-labels
  bc_main <- melt(bc_wide_loc,
    id.vars   = c("lfc_cutoff", "n_total"),
    measure.vars = paste0("p", MAIN_THRESH),
    variable.name = "threshold", value.name = "pct_patients")
  bc_main[, label := paste0(sub("p", "", threshold), "%")]
  bc_main[, label := factor(label, levels = paste0(MAIN_THRESH, "%"))]

  # Faint intermediate lines — no points, drawn below main lines
  bc_inter <- melt(bc_wide_loc,
    id.vars   = c("lfc_cutoff", "n_total"),
    measure.vars = paste0("p", INTER_THRESH),
    variable.name = "threshold", value.name = "pct_patients")
  bc_inter[, inter_label := paste0(sub("p", "", threshold), "%")]

  thresh_colors3   <- c("50%" = "#9E9E9E", "75%" = col_combined, "90%" = "#4527A0")
  bc_end_labels    <- bc_main[lfc_cutoff == max(lfc_cutoff)]
  bc_inter_labels  <- bc_inter[lfc_cutoff == max(lfc_cutoff)]

  x_expand <- expansion(mult = c(0.02, 0.12))

  pTop_loc <- ggplot(bc_wide_loc, aes(x = lfc_cutoff, y = n_total)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_segment(aes(xend = lfc_cutoff, y = 0, yend = n_total),
                 color = col_combined, linewidth = 0.6, alpha = 0.8) +
    geom_point(color = col_combined, size = 2.2) +
    scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0), expand = x_expand) +
    scale_y_continuous(labels = label_comma()) +
    annotate("text", x = 0.52, y = Inf, label = "|LFC| > 0.5",
             hjust = 0, vjust = 1.4, size = 2.8, color = "grey40") +
    labs(x = NULL, y = "Number of DEGs (lfsr < 0.05)",
         title = paste("B ", title_tag)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(6, 10, 2, 8))

  pBot_loc <- ggplot(bc_main, aes(x = lfc_cutoff, y = pct_patients,
                                   color = label, group = label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    # Faint intermediate lines drawn first (behind main lines)
    geom_line(data = bc_inter,
              aes(x = lfc_cutoff, y = pct_patients, group = threshold),
              color = "grey82", linewidth = 0.32, inherit.aes = FALSE) +
    geom_text(data = bc_inter_labels,
              aes(x = lfc_cutoff, y = pct_patients, label = inter_label),
              hjust = -0.2, size = 2.3, color = "grey72", inherit.aes = FALSE) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.6) +
    geom_text(data = bc_end_labels, aes(label = label),
              hjust = -0.2, size = 2.8, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = thresh_colors3, guide = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0), expand = x_expand) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       labels = paste0(seq(0, 100, 25), "%")) +
    labs(x = expression("|shrunk log"[2]*"FC| cutoff  (lfsr < 0.05)"),
         y = "% patients concordant") +
    theme(plot.margin = margin(2, 10, 6, 8))

  # Align panel widths via gtable so the vline at 0.5 stays in sync
  library(grid); library(gridExtra)
  g1 <- ggplotGrob(pTop_loc)
  g2 <- ggplotGrob(pBot_loc)
  g1$widths <- g2$widths <- unit.pmax(g1$widths, g2$widths)

  outfile <- file.path(OUTDIR, paste0("deg_count_patient_concordance", file_suffix, ".pdf"))
  cairo_pdf(outfile, width = 5, height = 5.5)
  grid.arrange(g1, g2, ncol = 1, heights = c(1.1, 1))
  dev.off()
  message("Saved: ", basename(outfile))
}

# ── Standard BC (sign-only directional concordance) ──────────────────────────
make_bc_panels(patient_results[lfc_cutoff != 0.4],
               pct_col    = "pct_combined",
               title_tag  = "DEG count and patient-level replication",
               file_suffix = "")

# ── Strict BC (patient |LFC| ≥ 0.5 required per DEG) ─────────────────────────
PAT_THRESH <- 0.5

message("Computing strict concordance (patient |LFC| >= 0.5 required)...")
strict_results <- rbindlist(lapply(lfc_cutoffs, function(cut) {
  up_idx <- which(dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  cut)
  dn_idx <- which(dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -cut)
  n_up <- length(up_idx); n_dn <- length(dn_idx)
  n_total <- n_up + n_dn
  if (n_total == 0) return(NULL)

  mat_up <- if (n_up > 0) lfc_m[up_idx, , drop = FALSE] else matrix(0, 0, n_patients)
  mat_dn <- if (n_dn > 0) lfc_m[dn_idx, , drop = FALSE] else matrix(0, 0, n_patients)

  n_concordant <- numeric(n_patients)
  if (n_up > 0) n_concordant <- n_concordant + colSums(mat_up >  PAT_THRESH)
  if (n_dn > 0) n_concordant <- n_concordant + colSums(mat_dn < -PAT_THRESH)

  n_expressed <- numeric(n_patients)
  if (n_up > 0) n_expressed <- n_expressed + colSums(abs(mat_up) >= PAT_THRESH)
  if (n_dn > 0) n_expressed <- n_expressed + colSums(abs(mat_dn) >= PAT_THRESH)

  pct_strict <- ifelse(n_expressed > 0, n_concordant / n_expressed * 100, NA_real_)

  data.table(sample_id = sample_ids, lfc_cutoff = cut,
             n_up_degs = n_up, n_dn_degs = n_dn, n_degs = n_total,
             pct_strict = pct_strict)
}))

make_bc_panels(strict_results[lfc_cutoff != 0.4],
               pct_col    = "pct_strict",
               title_tag  = "DEG count and patient-level replication (strict)",
               file_suffix = "_strict")

# ── Panel A strict: violin using same |LFC| ≥ 0.5 per-patient filter ─────────
strict_results[, cutoff_label := factor(
  sprintf("%.1f", lfc_cutoff),
  levels = sprintf("%.1f", lfc_cutoffs)
)]

pA_strict <- ggplot(strict_results[!is.na(pct_strict)],
                    aes(x = cutoff_label, y = pct_strict)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_hline(yintercept = 75, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  stat_halfeye(
    fill = col_combined, color = col_combined, alpha = 0.65,
    adjust = 0.8, width = 0.95,
    point_interval = median_qi, .width = c(0.5, 0.9),
    point_size = 1.0, interval_size = 0.6,
    slab_linewidth = 0.2
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = paste0(c(0, 25, 50, 75, 100), "%")) +
  labs(
    x     = expression("|shrunk log"[2]*"FC| cutoff  (dream DEGs, lfsr < 0.05)"),
    y     = "% of DEGs concordant per patient",
    title = "A  Patient-level concordance (patient |LFC| \u2265 0.5 required)"
  ) +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(OUTDIR, "patient_concordance_violins_strict.pdf"),
       pA_strict, width = 9, height = 4, device = cairo_pdf)
message("Saved: A_patient_concordance_violin_strict05.pdf")

# ── Panel D: Violin — Spearman ρ per patient at each LFC cutoff ──────────────
# Rho = rank correlation between dream LFC and patient LFC across the DEG set.
# Higher = patient's personal effect sizes align better with population direction.

pd_dat <- patient_results[!is.na(spearman_rho)]

pD <- ggplot(pd_dat, aes(x = cutoff_label, y = spearman_rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  stat_halfeye(
    fill = col_combined, color = col_combined, alpha = 0.65,
    adjust = 0.8, width = 0.95,
    point_interval = median_qi, .width = c(0.5, 0.9),
    point_size = 1.0, interval_size = 0.6,
    slab_linewidth = 0.2
  ) +
  scale_y_continuous(breaks = seq(-0.2, 1.0, 0.2),
                     labels = function(x) sprintf("%.1f", x)) +
  labs(
    x     = expression("|shrunk log"[2]*"FC| cutoff  (lfsr < 0.05)"),
    y     = expression("Spearman "*rho),
    title = expression("D  Per-patient rank concordance  "*rho)
  ) +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(OUTDIR, "patient_spearman_rho_violins.pdf"),
       pD, width = 9, height = 4, device = cairo_pdf)
message("Saved: D_spearman_rho_violin.pdf")

# ── Summary table ─────────────────────────────────────────────────────────────
summary_tbl <- patient_results[, .(
  n_degs          = n_degs[1],
  median_pct_combined = round(median(pct_combined), 1),
  pct75_combined      = round(mean(pct_combined >= 75) * 100, 1),
  median_rho          = round(median(spearman_rho, na.rm = TRUE), 3)
), by = .(lfc_cutoff)]

message("\n===== SUMMARY =====")
print(summary_tbl)
fwrite(summary_tbl, file.path(OUTDIR, "patient_concordance_summary.csv"))
message("\nAll panels saved to: ", OUTDIR)

# ═══════════════════════════════════════════════════════════════════════════════
# SANITY CHECK: Repeat concordance for HEALTHY CONTROLS
# Controls are used as the reference mean → their LFC ~ 0 → concordance ~ 50%
# ═══════════════════════════════════════════════════════════════════════════════

message("\n── Loading STAR DGE to compute control LFC matrix ──")
INT_RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
dge_s <- readRDS(file.path(INT_RDIR, "merged_dge.rds"))

ycfg_d     <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_c     <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg_d))
dge_s      <- dge_s[, dge_s$samples$dataset %in% mega_c]
lcpm_s     <- edgeR::cpm(dge_s, log = TRUE, prior.count = 1)

sinfo <- data.table(sample_id = colnames(dge_s),
                    dataset   = dge_s$samples$dataset,
                    group     = as.character(dge_s$samples$group_binary))
valid_ds <- sinfo[, .(nc = sum(group == "Control"), nd = sum(group == "Disease")),
                  by = dataset][nc > 0 & nd > 0, dataset]

# Within-dataset control means (same reference as disease patients)
cmeans <- lapply(setNames(valid_ds, valid_ds), function(ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  if (length(ids) == 1) lcpm_s[, ids] else rowMeans(lcpm_s[, ids, drop = FALSE])
})

ctrl_sids <- sinfo[group == "Control" & dataset %in% valid_ds, sample_id]
ctrl_lfc  <- matrix(NA_real_, nrow = nrow(lcpm_s), ncol = length(ctrl_sids),
                    dimnames = list(rownames(lcpm_s), ctrl_sids))
for (ds in valid_ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  for (sid in ids) ctrl_lfc[, sid] <- lcpm_s[, sid] - cmeans[[ds]]
}

# Align to dream_sub gene order
common_ctrl <- intersect(dream_sub$gene, rownames(ctrl_lfc))
ctrl_lfc    <- ctrl_lfc[dream_sub$gene[dream_sub$gene %in% common_ctrl], , drop = FALSE]
n_ctrl      <- ncol(ctrl_lfc)
message(sprintf("  Control LFC matrix: %d genes × %d controls", nrow(ctrl_lfc), n_ctrl))

ctrl_sid_ds <- setNames(sinfo[sample_id %in% ctrl_sids, dataset],
                        sinfo[sample_id %in% ctrl_sids, sample_id])

# ── Compute concordance for controls ─────────────────────────────────────────
ctrl_concordance <- rbindlist(lapply(lfc_cutoffs, function(cut) {
  up_idx <- which(dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  cut &
                    dream_sub$gene %in% rownames(ctrl_lfc))
  dn_idx <- which(dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -cut &
                    dream_sub$gene %in% rownames(ctrl_lfc))
  n_up <- length(up_idx); n_dn <- length(dn_idx)
  n_total <- n_up + n_dn
  if (n_total == 0) return(NULL)

  # Standard concordance (sign only)
  n_conc <- numeric(n_ctrl)
  if (n_up > 0) n_conc <- n_conc + colSums(ctrl_lfc[up_idx, , drop = FALSE] > 0)
  if (n_dn > 0) n_conc <- n_conc + colSums(ctrl_lfc[dn_idx, , drop = FALSE] < 0)
  pct_combined <- n_conc / n_total * 100

  # Strict concordance (|ctrl_lfc| >= 0.5 required)
  n_conc_s <- numeric(n_ctrl); n_expr_s <- numeric(n_ctrl)
  if (n_up > 0) {
    m <- ctrl_lfc[up_idx, , drop = FALSE]
    n_conc_s <- n_conc_s + colSums(m >  PAT_THRESH)
    n_expr_s <- n_expr_s + colSums(abs(m) >= PAT_THRESH)
  }
  if (n_dn > 0) {
    m <- ctrl_lfc[dn_idx, , drop = FALSE]
    n_conc_s <- n_conc_s + colSums(m < -PAT_THRESH)
    n_expr_s <- n_expr_s + colSums(abs(m) >= PAT_THRESH)
  }
  pct_strict <- ifelse(n_expr_s > 0, n_conc_s / n_expr_s * 100, NA_real_)

  data.table(sample_id = ctrl_sids, lfc_cutoff = cut,
             n_up_degs = n_up, n_dn_degs = n_dn, n_degs = n_total,
             pct_combined = pct_combined, pct_strict = pct_strict)
}))

ctrl_concordance[, cutoff_label := factor(sprintf("%.1f", lfc_cutoff),
                                          levels = sprintf("%.1f", lfc_cutoffs))]
ctrl_concordance[, group := "Control"]
patient_results[,    group := "Disease"]
strict_results[,     group := "Disease"]

# ── Helper: violin with disease + control overlaid ────────────────────────────
group_colors <- c("Disease" = col_combined, "Control" = "#757575")

make_a_with_ctrl <- function(dis_dat, ctrl_dat, pct_col, file_name) {
  combined <- rbind(
    dis_dat[,  .(sample_id, cutoff_label, pct = get(pct_col), group)],
    ctrl_dat[, .(sample_id, cutoff_label, pct = get(pct_col), group)]
  )
  combined <- combined[!is.na(pct)]

  ggplot(combined, aes(x = cutoff_label, y = pct, fill = group, color = group)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey60", linewidth = 0.35) +
    geom_hline(yintercept = 75, linetype = "dashed", color = "grey60", linewidth = 0.35) +
    stat_halfeye(
      alpha = 0.55, adjust = 0.8, width = 0.95,
      point_interval = median_qi, .width = c(0.5, 0.9),
      point_size = 0.7, interval_size = 0.45, slab_linewidth = 0.18,
      position = position_dodge(width = 0.85)
    ) +
    scale_fill_manual(values = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                       labels = paste0(c(0, 25, 50, 75, 100), "%")) +
    labs(x = expression("|shrunk log"[2]*"FC| cutoff  (dream DEGs, lfsr < 0.05)"),
         y = "% of DEGs concordant per patient") +
    theme(axis.text.x = element_text(size = 8),
          legend.position = "top",
          legend.key.size = unit(0.35, "cm"))
}

pA_ctrl <- make_a_with_ctrl(patient_results, ctrl_concordance,
                             "pct_combined", NULL)
ggsave(file.path(OUTDIR, "patient_concordance_violins_with_controls.pdf"),
       pA_ctrl, width = 13, height = 4, device = cairo_pdf)
message("Saved: A_patient_concordance_violin_with_controls.pdf")

pA_ctrl_s <- make_a_with_ctrl(strict_results, ctrl_concordance,
                               "pct_strict", NULL)
ggsave(file.path(OUTDIR, "patient_concordance_violins_strict_with_controls.pdf"),
       pA_ctrl_s, width = 13, height = 4, device = cairo_pdf)
message("Saved: A_patient_concordance_violin_strict05_with_controls.pdf")

# ── Helper: BC dual-panel with disease (solid) + control (dashed) ─────────────
make_bc_with_ctrl <- function(dis_results, ctrl_results, pct_col,
                               title_tag, file_suffix) {
  all_thresholds <- sort(c(MAIN_THRESH, INTER_THRESH))

  compute_thresh <- function(dat, col) {
    dat[!is.na(get(col)), {
      row_vals <- lapply(all_thresholds, function(thr)
        mean(get(col) >= thr, na.rm = TRUE) * 100)
      setNames(as.list(row_vals), paste0("p", all_thresholds))
    }, by = .(lfc_cutoff)]
  }

  dis_pct  <- compute_thresh(dis_results,  pct_col)
  ctrl_pct <- compute_thresh(ctrl_results, pct_col)

  bc_wide_loc <- merge(dis_results[, .(n_total = n_degs[1]), by = lfc_cutoff],
                       dis_pct, by = "lfc_cutoff")

  # Main disease lines
  bc_main <- melt(bc_wide_loc,
    id.vars = c("lfc_cutoff", "n_total"), measure.vars = paste0("p", MAIN_THRESH),
    variable.name = "threshold", value.name = "pct_patients")
  bc_main[, label := paste0(sub("p", "", threshold), "%")]
  bc_main[, label := factor(label, levels = paste0(MAIN_THRESH, "%"))]

  # Intermediate disease lines
  bc_inter <- melt(bc_wide_loc,
    id.vars = c("lfc_cutoff", "n_total"), measure.vars = paste0("p", INTER_THRESH),
    variable.name = "threshold", value.name = "pct_patients")
  bc_inter[, inter_label := paste0(sub("p", "", threshold), "%")]

  # Control main lines (dashed overlay)
  ctrl_wide <- merge(ctrl_results[, .(n_total = n_degs[1]), by = lfc_cutoff],
                     ctrl_pct, by = "lfc_cutoff")
  bc_ctrl <- melt(ctrl_wide,
    id.vars = c("lfc_cutoff", "n_total"), measure.vars = paste0("p", MAIN_THRESH),
    variable.name = "threshold", value.name = "pct_patients")
  bc_ctrl[, label := paste0(sub("p", "", threshold), "%")]
  bc_ctrl[, label := factor(label, levels = paste0(MAIN_THRESH, "%"))]

  thresh_colors3 <- c("50%" = "#9E9E9E", "75%" = col_combined, "90%" = "#4527A0")
  x_expand <- expansion(mult = c(0.02, 0.14))
  lfc_max  <- max(dis_results$lfc_cutoff)
  bc_end_labels  <- bc_main[ lfc_cutoff == lfc_max]
  bc_inter_labels <- bc_inter[lfc_cutoff == lfc_max]
  bc_ctrl_labels  <- bc_ctrl[ lfc_cutoff == lfc_max]

  pTop_loc <- ggplot(bc_wide_loc, aes(x = lfc_cutoff, y = n_total)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_segment(aes(xend = lfc_cutoff, y = 0, yend = n_total),
                 color = col_combined, linewidth = 0.6, alpha = 0.8) +
    geom_point(color = col_combined, size = 2.2) +
    scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0), expand = x_expand) +
    scale_y_continuous(labels = label_comma()) +
    annotate("text", x = 0.52, y = Inf, label = "|LFC| > 0.5",
             hjust = 0, vjust = 1.4, size = 2.8, color = "grey40") +
    labs(x = NULL, y = "Number of DEGs (lfsr < 0.05)",
         title = paste("B ", title_tag)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(6, 10, 2, 8))

  pBot_loc <- ggplot(bc_main, aes(x = lfc_cutoff, y = pct_patients,
                                   color = label, group = label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    # Faint intermediate disease lines
    geom_line(data = bc_inter,
              aes(x = lfc_cutoff, y = pct_patients, group = threshold),
              color = "grey82", linewidth = 0.32, inherit.aes = FALSE) +
    geom_text(data = bc_inter_labels,
              aes(x = lfc_cutoff, y = pct_patients, label = inter_label),
              hjust = -0.2, size = 2.3, color = "grey72", inherit.aes = FALSE) +
    # Control lines (dashed, same palette, no points)
    geom_line(data = bc_ctrl,
              aes(x = lfc_cutoff, y = pct_patients, color = label, group = label),
              linetype = "dashed", linewidth = 0.5, alpha = 0.6, inherit.aes = FALSE) +
    geom_text(data = bc_ctrl_labels,
              aes(x = lfc_cutoff, y = pct_patients, label = paste0(label, "\u00B7ctrl"),
                  color = label),
              hjust = -0.15, size = 2.2, inherit.aes = FALSE) +
    # Disease lines (solid, with points)
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.6) +
    geom_text(data = bc_end_labels, aes(label = label),
              hjust = -0.15, size = 2.8, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = thresh_colors3, guide = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0), expand = x_expand) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       labels = paste0(seq(0, 100, 25), "%")) +
    labs(x = expression("|shrunk log"[2]*"FC| cutoff  (lfsr < 0.05)"),
         y = "% patients concordant\n(solid=disease, dashed=control)") +
    theme(plot.margin = margin(2, 10, 6, 8))

  library(grid); library(gridExtra)
  g1 <- ggplotGrob(pTop_loc); g2 <- ggplotGrob(pBot_loc)
  g1$widths <- g2$widths <- unit.pmax(g1$widths, g2$widths)
  outfile <- file.path(OUTDIR, paste0("deg_count_patient_concordance",
                                      file_suffix, "_with_controls.pdf"))
  cairo_pdf(outfile, width = 5.5, height = 5.5)
  grid.arrange(g1, g2, ncol = 1, heights = c(1.1, 1))
  dev.off()
  message("Saved: ", basename(outfile))
}

make_bc_with_ctrl(patient_results[lfc_cutoff != 0.4], ctrl_concordance[lfc_cutoff != 0.4],
                  pct_col = "pct_combined",
                  title_tag = "DEG count and patient-level replication",
                  file_suffix = "")

make_bc_with_ctrl(strict_results[lfc_cutoff != 0.4], ctrl_concordance[lfc_cutoff != 0.4],
                  pct_col = "pct_strict",
                  title_tag = "DEG count and patient-level replication (strict)",
                  file_suffix = "_strict")

# ════════════════════════════════════════════════════════════════════════════
# NEW PANELS (2026-05-29): drive CONTROL baseline to ~0 + show per-patient robustness
#
# WHY the existing strict panels keep controls at ~50%:
#   pct_strict = n_concordant / n_expressed is a CONDITIONAL (precision-like) ratio.
#   n_expressed = #DEGs where the individual personally clears |LFC| >= 0.5. For a
#   control (LFC = deviation from its own cohort's control mean, ~0), only a few genes
#   cross 0.5 by noise, and among those the sign is ~50/50 -> control median pins ~50%.
#   ANY sign-only fraction floors at ~50% for controls. To send controls to 0 we either
#   (A1) require magnitude AND divide by ALL DEGs ("strict recall"), or use a
#   (A2) signed score / (A3) control-standardized z that is centered at 0 by design.
#
# All control computations below index ctrl_lfc BY GENE NAME (robust to the disease/
# control gene-set difference), not by dream_sub position.
# ════════════════════════════════════════════════════════════════════════════

message("\n========== NEW controls-to-zero + robustness panels ==========")
# All panels (original + new) share OUTDIR = .../figS_lfc_sensitivity/patient_concordance.
NEWDIR <- OUTDIR
grp_levels <- c("Disease", "Control")
mat_for    <- list(Disease = lfc_m, Control = ctrl_lfc)

# ── Shared per-sample metric sweep (strict recall, conditional ratio, sig score) ──
metric_dt <- rbindlist(lapply(lfc_cutoffs, function(cut) {
  rbindlist(lapply(grp_levels, function(grp) {
    mat <- mat_for[[grp]]
    gp  <- rownames(mat)
    up_g <- intersect(dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  cut], gp)
    dn_g <- intersect(dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -cut], gp)
    n_up <- length(up_g); n_dn <- length(dn_g); n_total <- n_up + n_dn
    if (n_total == 0) return(NULL)
    ns <- ncol(mat)

    mat_up <- if (n_up) mat[up_g, , drop = FALSE] else NULL
    mat_dn <- if (n_dn) mat[dn_g, , drop = FALSE] else NULL

    # concordant at magnitude: correct sign AND |LFC| >= PAT_THRESH
    n_conc <- numeric(ns)
    if (n_up) n_conc <- n_conc + colSums(mat_up >  PAT_THRESH)
    if (n_dn) n_conc <- n_conc + colSums(mat_dn < -PAT_THRESH)
    # discordant at magnitude: WRONG sign AND |LFC| >= PAT_THRESH (controls cross
    # |LFC|>=0.5 by noise ~symmetrically, so n_disc ~ n_conc -> net cancels to ~0)
    n_disc <- numeric(ns)
    if (n_up) n_disc <- n_disc + colSums(mat_up < -PAT_THRESH)
    if (n_dn) n_disc <- n_disc + colSums(mat_dn >  PAT_THRESH)
    # expressed-at-magnitude (the conditional-ratio denominator)
    n_expr <- numeric(ns)
    if (n_up) n_expr <- n_expr + colSums(abs(mat_up) >= PAT_THRESH)
    if (n_dn) n_expr <- n_expr + colSums(abs(mat_dn) >= PAT_THRESH)
    # threshold-free signature score = mean over DEGs of sign(dream)*patientLFC
    su <- if (n_up) colSums(mat_up) else rep(0, ns)
    sd_ <- if (n_dn) colSums(mat_dn) else rep(0, ns)

    data.table(sample_id = colnames(mat), group = grp, lfc_cutoff = cut,
               n_total = n_total, n_expressed = n_expr,
               strict_recall = n_conc / n_total * 100,
               disc_recall   = n_disc / n_total * 100,
               net_recall    = (n_conc - n_disc) / n_total * 100,
               cond_ratio    = ifelse(n_expr > 0, n_conc / n_expr * 100, NA_real_),
               sig_score     = (su - sd_) / n_total)
  }))
}))
metric_dt[, group := factor(group, levels = grp_levels)]
metric_dt[, cutoff_label := factor(sprintf("%.1f", lfc_cutoff),
                                   levels = sprintf("%.1f", lfc_cutoffs))]

# ── Control-derived per-gene null for the z-score (A3), with leave-one-out for ctrls ──
n_c       <- ncol(ctrl_lfc)
g_ctrl    <- rownames(ctrl_lfc)
sx        <- rowSums(ctrl_lfc)
sx2       <- rowSums(ctrl_lfc^2)
mean_full <- sx / n_c
sd_full   <- pmax(sqrt(pmax((sx2 - n_c * mean_full^2) / (n_c - 1), 0)), 0.1)  # sd floor 0.1
loo_mean  <- (sx - ctrl_lfc) / (n_c - 1)                                       # genes x n_c
loo_sd    <- pmax(sqrt(pmax((sx2 - ctrl_lfc^2 - (n_c - 1) * loo_mean^2) / (n_c - 2), 0)), 0.1)
dsign_v   <- setNames(sign(dream_sub$dream_logFC), dream_sub$gene)[g_ctrl]     # +1/-1 per ctrl gene

# Oriented z: disease vs full control null; controls vs leave-one-out null
dis_z  <- ((lfc_m[g_ctrl, , drop = FALSE] - mean_full) / sd_full) * dsign_v
ctrl_z <- ((ctrl_lfc - loo_mean) / loo_sd) * dsign_v

z_dt <- rbindlist(lapply(lfc_cutoffs, function(cut) {
  deg_g <- intersect(dream_sub$gene[dream_sub$dream_padj < 0.05 & abs(dream_sub$dream_logFC) > cut], g_ctrl)
  if (length(deg_g) < 5) return(NULL)
  rbind(
    data.table(sample_id = colnames(lfc_m),    group = "Disease", lfc_cutoff = cut,
               value = colMeans(dis_z[deg_g, , drop = FALSE])),
    data.table(sample_id = colnames(ctrl_lfc), group = "Control", lfc_cutoff = cut,
               value = colMeans(ctrl_z[deg_g, , drop = FALSE]))
  )
}))
z_dt[, group := factor(group, levels = grp_levels)]
z_dt[, cutoff_label := factor(sprintf("%.1f", lfc_cutoff), levels = sprintf("%.1f", lfc_cutoffs))]

# Console diagnostics so verification is readable from the log
.diag <- metric_dt[lfc_cutoff == 0.5, .(med_cond_ratio    = round(median(cond_ratio, na.rm = TRUE), 1),
                                        med_strict_recall = round(median(strict_recall), 1),
                                        med_net_recall    = round(median(net_recall), 1),
                                        med_sig_score     = round(median(sig_score), 3)), by = group]
message("At dream |LFC|>0.5 (PAT_THRESH=0.5) — controls should drop toward 0 for net/sig:")
print(.diag)
.zd <- z_dt[lfc_cutoff == 0.5, .(med_z = round(median(value), 3)), by = group]
message("Control-standardized z at |LFC|>0.5 (control should be ~0):"); print(.zd)

# ── Generic disease-vs-control violin across the dream-cutoff sweep ───────────
make_metric_violin <- function(dat, ylab, title, outfile,
                               pct = FALSE, yref = NULL, width = 13) {
  dat <- dat[!is.na(value)]
  p <- ggplot(dat, aes(x = cutoff_label, y = value, fill = group, color = group))
  if (!is.null(yref))
    p <- p + geom_hline(yintercept = yref, linetype = "dashed", color = "grey55", linewidth = 0.4)
  p <- p +
    stat_halfeye(alpha = 0.55, adjust = 0.8, width = 0.95,
                 point_interval = median_qi, .width = c(0.5, 0.9),
                 point_size = 0.7, interval_size = 0.45, slab_linewidth = 0.18,
                 position = position_dodge(width = 0.85)) +
    scale_fill_manual(values = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(x = expression("|shrunk log"[2]*"FC| cutoff  (dream DEGs, lfsr < 0.05)"),
         y = ylab, title = title) +
    theme(axis.text.x = element_text(size = 8), legend.position = "top",
          legend.key.size = unit(0.35, "cm"))
  if (pct) {
    p <- p + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                                labels = paste0(seq(0, 100, 25), "%"))
  } else {
    yl <- as.numeric(quantile(dat$value, c(0.005, 0.995), na.rm = TRUE))
    p <- p + coord_cartesian(ylim = yl)
  }
  ggsave(outfile, p, width = width, height = 4, device = cairo_pdf)
  message("Saved: ", basename(outfile))
}

# ── Generic disease-vs-control median line + IQR ribbon (continuous x) ────────
make_metric_lines <- function(dat, ylab, title, outfile,
                              pct = FALSE, yref = NULL, width = 6) {
  dat <- dat[!is.na(value)]
  p <- ggplot(dat, aes(x = lfc_cutoff, y = value, color = group, fill = group))
  if (!is.null(yref))
    p <- p + geom_hline(yintercept = yref, linetype = "dashed", color = "grey55", linewidth = 0.4)
  p <- p +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey75", linewidth = 0.3) +
    stat_summary(fun = median, fun.min = function(x) quantile(x, 0.25),
                 fun.max = function(x) quantile(x, 0.75),
                 geom = "ribbon", alpha = 0.16, color = NA) +
    stat_summary(fun = median, geom = "line", linewidth = 0.9) +
    stat_summary(fun = median, geom = "point", size = 1.6) +
    scale_color_manual(values = group_colors, name = NULL) +
    scale_fill_manual(values = group_colors, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0)) +
    labs(x = expression("|shrunk log"[2]*"FC| cutoff  (lfsr < 0.05)"), y = ylab, title = title) +
    theme(legend.position = "top", legend.key.size = unit(0.35, "cm"))
  if (pct)
    p <- p + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                                labels = paste0(seq(0, 100, 25), "%"))
  ggsave(outfile, p, width = width, height = 4.2, device = cairo_pdf)
  message("Saved: ", basename(outfile))
}

# ── A1: Strict recall (global denominator) ───────────────────────────────────
# NOTE: requiring magnitude halves the control baseline (~50% -> ~31%) but does
# NOT reach 0: controls still cross |LFC|>=0.5 in the concordant direction by
# noise. Net recall (A1b) and the signed score (A2) subtract/cancel that noise.
a1_dat <- metric_dt[, .(sample_id, group, lfc_cutoff, cutoff_label, value = strict_recall)]
make_metric_violin(
  a1_dat,
  ylab    = "% of DEGs recapitulated per individual",
  title   = "A1  Strict recall — fraction of signature reproduced (correct sign & |LFC| ≥ 0.5 ÷ all DEGs)",
  outfile = file.path(NEWDIR, "A1_strict_recall.pdf"),
  pct = TRUE)
make_metric_lines(
  a1_dat,
  ylab    = "Median % of DEGs recapitulated",
  title   = "A1  Strict recall (median ± IQR): disease high/rising; controls retain a noise floor (~30%)",
  outfile = file.path(NEWDIR, "A1_strict_recall_lines.pdf"),
  pct = TRUE)

# ── A1b: Net recall (concordant − discordant) — count-based metric that ───────
#         reaches ~0 for controls (their symmetric |LFC|>=0.5 noise cancels) ───
a1b_dat <- metric_dt[, .(sample_id, group, lfc_cutoff, cutoff_label, value = net_recall)]
make_metric_violin(
  a1b_dat,
  ylab    = "Net recall (%)",
  title   = "A1b  Net recall = (concordant − discordant at |LFC| ≥ 0.5) ÷ all DEGs — controls fall to ~0",
  outfile = file.path(NEWDIR, "A1b_net_recall.pdf"),
  yref = 0)
make_metric_lines(
  a1b_dat,
  ylab    = "Median net recall (%)",
  title   = "A1b  Net recall (median ± IQR): controls ~0, disease high and rising",
  outfile = file.path(NEWDIR, "A1b_net_recall_lines.pdf"),
  yref = 0)

# ── A2: Threshold-free signature score ───────────────────────────────────────
a2_dat <- metric_dt[, .(sample_id, group, lfc_cutoff, cutoff_label, value = sig_score)]
make_metric_violin(
  a2_dat,
  ylab    = "Signature score",
  title   = "How strongly each patient reproduces the disease signature",
  outfile = file.path(NEWDIR, "A2_signature_score.pdf"),
  yref = 0)

# ── A3: Control-standardized z (leave-one-out null) ──────────────────────────
make_metric_violin(
  z_dt[, .(sample_id, group, lfc_cutoff, cutoff_label, value)],
  ylab    = "Per-individual mean z vs control null",
  title   = "A3  Control-standardized z (oriented; LOO for controls) — healthy controls are the zero line",
  outfile = file.path(NEWDIR, "A3_control_zscore.pdf"),
  yref = 0)

# ── B: Per-patient threshold (tau) sweep, incl. tau = 1.0 (user request) ──────
tau_set <- c(0, 0.5, 1.0, 1.5, 2.0, 3.0)
tau_dt <- rbindlist(lapply(tau_set, function(tau) {
  rbindlist(lapply(lfc_cutoffs, function(cut) {
    rbindlist(lapply(grp_levels, function(grp) {
      mat <- mat_for[[grp]]; gp <- rownames(mat)
      up_g <- intersect(dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  cut], gp)
      dn_g <- intersect(dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -cut], gp)
      nt <- length(up_g) + length(dn_g); if (nt == 0) return(NULL)
      nc <- numeric(ncol(mat))
      if (length(up_g)) nc <- nc + colSums(mat[up_g, , drop = FALSE] >  tau)
      if (length(dn_g)) nc <- nc + colSums(mat[dn_g, , drop = FALSE] < -tau)
      data.table(sample_id = colnames(mat), group = grp, lfc_cutoff = cut,
                 tau = tau, value = nc / nt * 100)
    }))
  }))
}))
tau_dt[, group := factor(group, levels = grp_levels)]
tau_dt[, tau_label := factor(sprintf("patient |LFC| ≥ %.1f", tau),
                             levels = sprintf("patient |LFC| ≥ %.1f", tau_set))]

pB <- ggplot(tau_dt, aes(x = lfc_cutoff, y = value, color = group, fill = group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey75", linewidth = 0.3) +
  stat_summary(fun = median, fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75),
               geom = "ribbon", alpha = 0.16, color = NA) +
  stat_summary(fun = median, geom = "line", linewidth = 0.8) +
  stat_summary(fun = median, geom = "point", size = 1.3) +
  facet_wrap(~ tau_label, ncol = 3) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), "%")) +
  labs(x = expression("shrunk dream |log"[2]*"FC| cutoff  (lfsr < 0.05)"),
       y = "% of DEGs concordant per patient",
       title = "B  Patient concordance by per-patient LFC requirement") +
  theme(legend.position = "top", legend.key.size = unit(0.35, "cm"),
        axis.text.x = element_text(size = 8))
ggsave(file.path(NEWDIR, "B_patient_threshold_sweep.pdf"),
       pB, width = 12, height = 8, device = cairo_pdf)
message("Saved: B_patient_threshold_sweep.pdf")

# ── C1: Didactic "why controls sat at ~50%" ──────────────────────────────────
# Top: in controls, DEGs that cross |LFC|>=0.5 split ~evenly into concordant vs
# discordant (random direction), so the conditional ratio = conc/(conc+disc) ~ 50%.
cd  <- metric_dt[, .(concordant = median(strict_recall), discordant = median(disc_recall)),
                 by = .(group, lfc_cutoff)]
cdm <- melt(cd, id.vars = c("group", "lfc_cutoff"),
            variable.name = "direction", value.name = "pct")
pC1top <- ggplot(cdm, aes(x = lfc_cutoff, y = pct, color = group, linetype = direction)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey75", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_linetype_manual(values = c("concordant" = "solid", "discordant" = "dotted"), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "% of all DEGs at |LFC| ≥ 0.5",
       title = "C1  Why the old strict ratio kept controls at ~50%",
       subtitle = "Top: in controls, concordant ≈ discordant (|LFC|≥0.5 crossed in random directions)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7), plot.margin = margin(6, 10, 2, 8))

# Bottom: same data, two metrics. Conditional ratio (~50% for ctrl) vs strict recall (~0 for ctrl).
cmp <- metric_dt[, .(`1 · conditional ratio (sign only ÷ n_expressed)` = median(cond_ratio, na.rm = TRUE),
                     `2 · strict recall (sign & |LFC|≥0.5 ÷ all DEGs)`  = median(strict_recall),
                     `3 · net recall (concordant − discordant ÷ all)`   = median(net_recall)),
                 by = .(group, lfc_cutoff)]
cmpm <- melt(cmp, id.vars = c("group", "lfc_cutoff"),
             variable.name = "metric", value.name = "pct")
pC1bot <- ggplot(cmpm, aes(x = lfc_cutoff, y = pct, color = group, linetype = metric)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey75", linewidth = 0.3) +
  geom_hline(yintercept = c(0, 50), linetype = "dotted", color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%")) +
  coord_cartesian(ylim = c(-5, 100)) +
  labs(x = expression("shrunk dream |log"[2]*"FC| cutoff  (lfsr < 0.05)"),
       y = "median % concordant",
       subtitle = "Bottom: each refinement lowers the control floor (50% → ~31% → ~0); disease stays high") +
  theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7), plot.margin = margin(2, 10, 6, 8))

g1 <- ggplotGrob(pC1top); g2 <- ggplotGrob(pC1bot)
g1$widths <- g2$widths <- unit.pmax(g1$widths, g2$widths)
cairo_pdf(file.path(NEWDIR, "C1_why_controls_50pct.pdf"), width = 6.5, height = 6)
grid.arrange(g1, g2, ncol = 1, heights = c(1, 1.1))
dev.off()
message("Saved: C1_why_controls_50pct.pdf")

# ── C2: Per-patient robustness scatter grid (dream LFC vs each individual's LFC) ─
# Representative individuals spanning the strict-recall distribution at |LFC|>0.5.
deg05  <- dream_sub[dream_padj < 0.05 & abs(dream_logFC) > 0.5 & gene %in% g_ctrl]
gsel   <- deg05$gene
dx     <- setNames(deg05$dream_logFC, gsel)
rec_d  <- metric_dt[group == "Disease" & lfc_cutoff == 0.5]
rec_c  <- metric_dt[group == "Control" & lfc_cutoff == 0.5]
pick_by_q <- function(tab, qs) {
  vapply(qs, function(q) tab$sample_id[which.min(abs(tab$strict_recall - quantile(tab$strict_recall, q, na.rm = TRUE)))],
         character(1))
}
qgrid  <- seq(0.10, 0.90, length.out = 5)   # 5 individuals spanning recall (low -> high)
pick_d <- unique(pick_by_q(rec_d, qgrid))
pick_c <- unique(pick_by_q(rec_c, qgrid))

scat <- rbindlist(c(
  lapply(seq_along(pick_d), function(i) {
    s <- pick_d[i]; rec <- rec_d[sample_id == s, strict_recall]
    data.table(facet = sprintf("Disease %d (recall %2.0f%%)", i, rec), ord = i,
               gene = gsel, dream = dx, patient = lfc_m[gsel, s])
  }),
  lapply(seq_along(pick_c), function(i) {
    s <- pick_c[i]; rec <- rec_c[sample_id == s, strict_recall]
    data.table(facet = sprintf("Control %d (recall %2.0f%%)", i, rec), ord = 100 + i,
               gene = gsel, dream = dx, patient = ctrl_lfc[gsel, s])
  })
))
# Magnitude-aware per-gene status (matches the strict/net recall logic):
scat[, status := fifelse(abs(patient) < 0.5, "sub-threshold (|LFC| < 0.5)",
                 fifelse(sign(dream) == sign(patient), "concordant (|LFC| ≥ 0.5)",
                                                        "discordant (|LFC| ≥ 0.5)"))]
scat[, status := factor(status, levels = c("concordant (|LFC| ≥ 0.5)",
                                            "discordant (|LFC| ≥ 0.5)",
                                            "sub-threshold (|LFC| < 0.5)"))]
scat[, facet := factor(facet, levels = unique(facet[order(ord)]))]
scat <- scat[order(-as.integer(status))]  # draw concordant on top

pC2 <- ggplot(scat, aes(x = dream, y = patient, color = status)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~ facet, ncol = 5) +
  scale_color_manual(values = c("concordant (|LFC| ≥ 0.5)"   = col_combined,
                                "discordant (|LFC| ≥ 0.5)"    = "#C0392B",
                                "sub-threshold (|LFC| < 0.5)" = "#D9D9D9"), name = NULL) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = expression("population (dream) log"[2]*"FC of DEGs  (|LFC| > 0.5)"),
       y = expression("individual's log"[2]*"FC vs own-cohort control mean"),
       title = "C2  Per-individual robustness: each disease patient reproduces the population signature") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme(legend.position = "top", legend.key.size = unit(0.35, "cm"),
        strip.text = element_text(size = 8))
ggsave(file.path(NEWDIR, "C2_robustness_scatter.pdf"),
       pC2, width = 11, height = 5.2, device = cairo_pdf)
message("Saved: C2_robustness_scatter.pdf")

# ── D: Density of per-patient LFC for the signature (disease vs control) ──────
# x = patient log2FC = (sample − own-cohort control mean), for the canonical
# |LFC| > 0.5 DEG set. Genes indexed by name (robust to disease/control gene-set diff).
DTHR <- 0.5
deg_u_all <- dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  DTHR]
deg_d_all <- dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -DTHR]
du_d <- intersect(deg_u_all, rownames(lfc_m));    dd_d <- intersect(deg_d_all, rownames(lfc_m))
du_c <- intersect(deg_u_all, rownames(ctrl_lfc)); dd_c <- intersect(deg_d_all, rownames(ctrl_lfc))

# D1 — oriented density: sign(dream) × patient LFC (concordant = positive)
dens_or <- rbindlist(list(
  data.table(group = "Disease", value = c( as.vector(lfc_m[du_d, ]),    -as.vector(lfc_m[dd_d, ]))),
  data.table(group = "Control", value = c( as.vector(ctrl_lfc[du_c, ]), -as.vector(ctrl_lfc[dd_c, ])))
))
dens_or[, group := factor(group, levels = grp_levels)]
pD1 <- ggplot(dens_or, aes(x = value, fill = group, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_density(alpha = 0.45, linewidth = 0.5, adjust = 1.1) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = group_colors, name = NULL) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = "patient log2FC oriented to dream direction  =  sign(dream) × (sample − cohort control mean)",
       y = "density",
       title = "Per-patient LFC distribution of dream DEGs, oriented to dream direction") +
  theme(legend.position = "top", legend.key.size = unit(0.35, "cm"))
ggsave(file.path(NEWDIR, "D1_density_oriented.pdf"), pD1, width = 8, height = 4.2, device = cairo_pdf)
message("Saved: D1_density_oriented.pdf")

# D2 — raw patient LFC (sample − cohort control mean), faceted by DEG direction
dens_raw <- rbindlist(list(
  data.table(group = "Disease", direction = "up-regulated DEGs",   value = as.vector(lfc_m[du_d, ])),
  data.table(group = "Disease", direction = "down-regulated DEGs", value = as.vector(lfc_m[dd_d, ])),
  data.table(group = "Control", direction = "up-regulated DEGs",   value = as.vector(ctrl_lfc[du_c, ])),
  data.table(group = "Control", direction = "down-regulated DEGs", value = as.vector(ctrl_lfc[dd_c, ]))
))
dens_raw[, group := factor(group, levels = grp_levels)]
dens_raw[, direction := factor(direction, levels = c("up-regulated DEGs", "down-regulated DEGs"))]
pD2 <- ggplot(dens_raw, aes(x = value, fill = group, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_density(alpha = 0.45, linewidth = 0.5, adjust = 1.1) +
  facet_wrap(~ direction) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = group_colors, name = NULL) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = "patient log2FC  (sample − own-cohort control mean)", y = "density",
       title = "Per-patient LFC distribution of dream DEGs") +
  theme(legend.position = "top", legend.key.size = unit(0.35, "cm"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(NEWDIR, "D2_density_raw_updown.pdf"), pD2, width = 9, height = 4.2, device = cairo_pdf)
message("Saved: D2_density_raw_updown.pdf")

# D3 — per-individual RAW LFC HISTOGRAM over the dream DEG signature (lfsr < 0.05,
# |LFC| > 0.5). 20 random disease + 20 random control. No recall metric or label.
allu   <- dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  0.5]
alld   <- dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -0.5]
allu_d <- intersect(allu, rownames(lfc_m));    alld_d <- intersect(alld, rownames(lfc_m))
allu_c <- intersect(allu, rownames(ctrl_lfc)); alld_c <- intersect(alld, rownames(ctrl_lfc))
all_d  <- c(allu_d, alld_d); all_c <- c(allu_c, alld_c)
ds_vec <- c(sid_to_ds, ctrl_sid_ds)
flab   <- function(s) sprintf("%s\n%s", s, ds_vec[s])
set.seed(42)
dis_samp  <- sample(colnames(lfc_m),    20)
ctrl_samp <- sample(colnames(ctrl_lfc), 20)

d3 <- rbindlist(c(
  lapply(seq_along(dis_samp),  function(i) data.table(facet = flab(dis_samp[i]),  group = "Disease",
                                                      ord = i,       value = lfc_m[all_d, dis_samp[i]])),
  lapply(seq_along(ctrl_samp), function(i) data.table(facet = flab(ctrl_samp[i]), group = "Control",
                                                      ord = 100 + i, value = ctrl_lfc[all_c, ctrl_samp[i]]))
))
d3[, group := factor(group, levels = grp_levels)]
d3[, facet := factor(facet, levels = unique(d3[order(ord), facet]))]
pD3 <- ggplot(d3, aes(x = value, fill = group, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 bins = 50, alpha = 0.55, linewidth = 0.1) +
  facet_wrap(~ facet, ncol = 10) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = group_colors, name = NULL) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Patient-level log2FC  (patient − cohort control mean)", y = "Density",
       title = "Per-individual LFC distribution over the disease signature") +
  theme(legend.position = "top", legend.key.size = unit(0.35, "cm"),
        strip.text = element_text(size = 6))
ggsave(file.path(NEWDIR, "D3_per_patient_density_grid.pdf"), pD3,
       width = 16, height = 9, device = cairo_pdf)
message("Saved: D3_per_patient_density_grid.pdf")

# D4 — same individuals as D3, ALL significant DEGs (lfsr < 0.05, no |LFC| cutoff) split
# by Dream direction (Up = magenta, Down = blue), overlaid per box. No recall.
d4 <- rbindlist(c(
  lapply(dis_samp, function(s) rbind(
    data.table(facet = flab(s), direction = "Up",   value = lfc_m[allu_d, s]),
    data.table(facet = flab(s), direction = "Down", value = lfc_m[alld_d, s]))),
  lapply(ctrl_samp, function(s) rbind(
    data.table(facet = flab(s), direction = "Up",   value = ctrl_lfc[allu_c, s]),
    data.table(facet = flab(s), direction = "Down", value = ctrl_lfc[alld_c, s])))
))
d4[, direction := factor(direction, levels = c("Up", "Down"))]
d4[, facet := factor(facet, levels = levels(d3$facet))]   # same individual order as D3
pD4 <- ggplot(d4, aes(x = value, fill = direction, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 bins = 50, alpha = 0.5, linewidth = 0.1) +
  facet_wrap(~ facet, ncol = 10) +
  scale_fill_manual(values = c("Up" = col_up, "Down" = col_down), name = "Dream direction") +
  scale_color_manual(values = c("Up" = col_up, "Down" = col_down), name = "Dream direction") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Patient-level log2FC  (patient − cohort control mean)", y = "Density",
       title = "Per-individual LFC distribution by Dream direction") +
  theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
        strip.text = element_text(size = 6))
ggsave(file.path(NEWDIR, "D4_per_patient_density_updown.pdf"), pD4,
       width = 16, height = 9, device = cairo_pdf)
message("Saved: D4_per_patient_density_updown.pdf")

# D5 — fibrosis-severity-stratified per-individual histogram (Up/Down Dream direction),
# dream DEGs |LFC| > 0.5: 10 severe (F3/F4) + 10 mild (F0/F1) disease + 10 control.
# No recall. (NAS too sparse to stratify on — 800 NA; using fibrosis_stage.)
stage_vec <- meta[, setNames(fibrosis_stage, sample_id)]
sev_pool  <- intersect(names(stage_vec)[which(stage_vec %in% c(3, 4))], colnames(lfc_m))
mild_pool <- intersect(names(stage_vec)[which(stage_vec %in% c(0, 1))], colnames(lfc_m))
set.seed(42)
sev_s  <- sample(sev_pool,  min(10, length(sev_pool)))
mild_s <- sample(mild_pool, min(10, length(mild_pool)))
ctrl_s <- sample(colnames(ctrl_lfc), 10)
message(sprintf("D5 pools: severe F3/F4=%d, mild F0/F1=%d, controls=%d",
                length(sev_pool), length(mild_pool), ncol(ctrl_lfc)))

mk5 <- function(samp, mat, gu, gd, sevlab) rbindlist(lapply(samp, function(s) rbind(
  data.table(srr = s, gse = ds_vec[s], sevgroup = sevlab, direction = "Up",   value = mat[gu, s]),
  data.table(srr = s, gse = ds_vec[s], sevgroup = sevlab, direction = "Down", value = mat[gd, s]))))
d5 <- rbindlist(list(
  mk5(sev_s,  lfc_m,    allu_d, alld_d, "F3–F4 (severe)"),
  mk5(mild_s, lfc_m,    allu_d, alld_d, "F0–F1 (mild)"),
  mk5(ctrl_s, ctrl_lfc, allu_c, alld_c, "Control")
))
d5[, sevgroup := factor(sevgroup, levels = c("F3–F4 (severe)", "F0–F1 (mild)", "Control"))]
d5[, direction := factor(direction, levels = c("Up", "Down"))]
# severity = left row strip; individuals = columns within each row; SRR/GSE labelled in-panel
idx5 <- unique(d5[, .(sevgroup, srr)])[, col_idx := seq_len(.N), by = sevgroup]
d5   <- merge(d5, idx5, by = c("sevgroup", "srr"))
d5[, col_idx := factor(col_idx, levels = sort(unique(col_idx)))]
lab5 <- unique(d5[, .(sevgroup, col_idx, srr, gse)])
lab5[, lab := paste0(srr, "\n", gse)]
pD5 <- ggplot(d5, aes(x = value, fill = direction, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 bins = 50, alpha = 0.5, linewidth = 0.1) +
  geom_text(data = lab5, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.15, size = 1.8, lineheight = 0.85, color = "grey25") +
  facet_grid(sevgroup ~ col_idx, switch = "y") +
  scale_fill_manual(values = c("Up" = col_up, "Down" = col_down), name = "Dream direction") +
  scale_color_manual(values = c("Up" = col_up, "Down" = col_down), name = "Dream direction") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Patient-level log2FC  (patient − cohort control mean)", y = "Density",
       title = "Per-individual LFC by fibrosis severity",
       subtitle = "dream DEGs (lfsr < 0.05, |log2FC| > 0.5); Up/Down = dream-up / dream-down; 10 individuals per group") +
  theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
        strip.text.x = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),
        axis.title.y = element_text(margin = margin(r = 8)))
ggsave(file.path(NEWDIR, "D5_per_patient_density_by_fibrosis.pdf"), pD5,
       width = 16, height = 7.5, device = cairo_pdf)
message("Saved: D5_per_patient_density_by_fibrosis.pdf")

message("========== NEW panels complete -> ", NEWDIR, " ==========\n")
