#!/usr/bin/env Rscript
# figS_cas13_library_tpm_by_diet.R
# ---------------------------------------------------------------------------
# Median TPM in DISEASE mouse liver for library gene sets defined by 4 human
# ashr shrunk_logFC thresholds (0.40 / 0.30 / 0.20 / 0.15), across 4 diet
# groups (MCD, CDAHFD, Western, HFD).
#
# Each option = human spine (ashr > threshold) UNION mouse-confirmed tier
# (>=3 diets, human-concordant).  Mirrors figS_cas13_library_options.R.
# TPM computed from featureCounts raw counts + gene Length (gencode.vM38).
# Per-gene value = median TPM across disease samples in that diet.
#
# Output: S_lib_10_tpm_by_diet.pdf  ->  Cas13_Library_Design/figures/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR    <- FIGS_CAS13LIB_DIR
FC_DIR     <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/counts/featurecounts")
PUB_FC     <- file.path(BASE, "RNA-seq/Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
WD_DIR     <- file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets")
MAIN_META  <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv")
ASHR       <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                        "results/integration/dream_results_ashr.csv")
ORTHO      <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
MOUSE_META <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
PERDIET    <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")

DIETS         <- c("MCD", "CDAHFD", "Western", "HFD")
KEEP_BIOTYPES <- c("protein_coding", "lncRNA")
LFSR_THR      <- 0.05
SHRUNK_MOUSE  <- 0.5

# Four thresholds to compare (strict -> loose)
OPT_THRS  <- c(0.40, 0.30, 0.20, 0.15)
OPT_LABELS <- sprintf("ashr > %.2f", OPT_THRS)   # e.g. "ashr > 0.40"

strip_v <- function(x) sub("[.][0-9]+$", "", x)

# =============================================================================
# 1. Build gene sets for each threshold (human spine UNION mouse-confirmed tier)
# =============================================================================
message("Building library gene sets ...")
ash <- fread(ASHR, select = c("gene", "shrunk_logFC", "lfsr", "logFC"))
ash[, hb := strip_v(gene)]
human_logfc <- ash[, .(hb, hlfc = logFC)]

ortho <- fread(cmd = paste0("zcat ", ORTHO),
               select = c("mouse_ensembl", "human_ensembl", "confidence_tier", "is_one2one"))
ortho[, mouse_ensembl := strip_v(mouse_ensembl)]
ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H", "M")]
ortho[, trank    := match(confidence_tier, c("H", "M"))]
ortho[, o2o_rank := fifelse(is_one2one %in% c(TRUE, "True", "TRUE", "true"), 0L, 1L)]

setorder(ortho, human_ensembl, trank, o2o_rank, mouse_ensembl)
ortho_byhuman <- unique(ortho, by = "human_ensembl")   # best mouse per human

setorder(ortho, mouse_ensembl, trank, o2o_rank, human_ensembl)
ortho_m2h <- unique(ortho, by = "mouse_ensembl")        # best human per mouse

mouse_meta_dt <- fread(MOUSE_META)
mouse_meta_dt[grepl("protein_coding", mouse_biotype), bt2 := "protein_coding"]
mouse_meta_dt[grepl("lncRNA|lincRNA",  mouse_biotype), bt2 := "lncRNA"]
pc_lnc_ids <- mouse_meta_dt[bt2 %in% KEEP_BIOTYPES, mouse_ensembl_base]

# Mouse-confirmed tier (constant across thresholds)
up_by_diet <- lapply(DIETS, function(d) {
  dt <- fread(file.path(PERDIET, paste0(d, "_de_results.csv")),
              select = c("gene", "shrunk_logFC", "lfsr"))
  dt[, gene_base := strip_v(gene)]
  dt[lfsr < LFSR_THR & shrunk_logFC > SHRUNK_MOUSE, gene_base]
})
mouse_any <- unique(unlist(up_by_diet))
n_diets_up <- sapply(mouse_any, function(g) sum(sapply(up_by_diet, function(u) g %in% u)))
mouse_3plus <- mouse_any[n_diets_up >= 3]
mc_map <- merge(ortho_m2h[mouse_ensembl %in% mouse_3plus, .(mouse_ensembl, hb = human_ensembl)],
                human_logfc, by = "hb", all.x = TRUE)
mouse_conf <- intersect(mc_map[!is.na(hlfc) & hlfc > 0, mouse_ensembl], pc_lnc_ids)

get_spine <- function(thr) {
  hup <- ash[!is.na(lfsr) & lfsr < LFSR_THR & shrunk_logFC > thr, hb]
  intersect(ortho_byhuman[human_ensembl %in% hup, mouse_ensembl], pc_lnc_ids)
}

gene_sets <- lapply(OPT_THRS, function(thr) sort(union(get_spine(thr), mouse_conf)))
names(gene_sets) <- OPT_LABELS

for (lbl in OPT_LABELS)
  cat(sprintf("  %s: %d genes\n", lbl, length(gene_sets[[lbl]])))
cat(sprintf("  mouse-confirmed (constant): %d genes\n", length(mouse_conf)))

# =============================================================================
# 2. TPM computation helpers
# =============================================================================
load_fc <- function(path) {
  dt        <- fread(path, skip = "Geneid")
  samp_cols <- setdiff(names(dt), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  srr_ids   <- basename(dirname(samp_cols))
  lengths   <- setNames(dt[["Length"]], strip_v(dt[["Geneid"]]))
  mat       <- as.matrix(dt[, ..samp_cols])
  rownames(mat) <- strip_v(dt[["Geneid"]])
  colnames(mat) <- srr_ids
  list(counts = mat, lengths = lengths)
}

compute_tpm <- function(counts, lengths) {
  genes   <- intersect(rownames(counts), names(lengths))
  counts  <- counts[genes, , drop = FALSE]
  rpk     <- sweep(counts, 1, lengths[genes] / 1e3, "/")
  sweep(rpk, 2, colSums(rpk) / 1e6, "/")
}

row_medians <- function(mat, cols) {
  cols <- intersect(colnames(mat), cols)
  if (length(cols) == 0L) return(setNames(rep(NA_real_, nrow(mat)), rownames(mat)))
  if (length(cols) == 1L) return(setNames(mat[, cols], rownames(mat)))
  apply(mat[, cols, drop = FALSE], 1, median)
}

# =============================================================================
# 3. Compute per-diet median TPM across disease samples
# =============================================================================
meta <- fread(MAIN_META)
if ("group_binary" %in% names(meta)) {
  disease_meta <- meta[group_binary != "Control"]
} else {
  disease_meta <- meta[!grepl("Control|LFD|Chow|Ctrl", condition)]
}

message("Computing TPM: MCD ...")
mcd_fcs    <- lapply(file.path(FC_DIR, c("gene_counts_inhouse.txt",
                                          "gene_counts_gse156918.txt",
                                          "gene_counts_gse205974.txt")), load_fc)
gene_lengths <- mcd_fcs[[1]]$lengths   # reference lengths (same GTF throughout)
mcd_tpm    <- compute_tpm(do.call(cbind, lapply(mcd_fcs, `[[`, "counts")), gene_lengths)
mcd_med    <- row_medians(mcd_tpm, disease_meta[diet_model == "MCD", sample_id])
cat(sprintf("  MCD disease samples: %d\n",
            length(intersect(colnames(mcd_tpm), disease_meta[diet_model == "MCD", sample_id]))))

message("Computing TPM: CDAHFD ...")
pub_tpm    <- compute_tpm(load_fc(PUB_FC)$counts, gene_lengths)
cdahfd_med <- row_medians(pub_tpm, disease_meta[diet_model == "CDAHFD", sample_id])
cat(sprintf("  CDAHFD disease samples: %d\n",
            length(intersect(colnames(pub_tpm), disease_meta[diet_model == "CDAHFD", sample_id]))))

message("Computing TPM: Western ...")
gse220_fc  <- load_fc(file.path(WD_DIR, "GSE220575/counts/featurecounts/gene_counts.txt"))
gse220_tpm <- compute_tpm(gse220_fc$counts, gene_lengths)
gse220_m   <- fread(file.path(WD_DIR, "GSE220575/metadata/sample_metadata.csv"))
gse220_g2s <- fread(file.path(WD_DIR, "GSE220575/metadata/gsm_to_srr.tsv"))
gse220_dis <- merge(gse220_m, gse220_g2s, by.x = "gsm", by.y = "gsm")[
  condition %in% c("MASH", "HCC"), srr]

gse246_fc  <- load_fc(file.path(WD_DIR, "GSE246088/counts/featurecounts/gene_counts.txt"))
gse246_tpm <- compute_tpm(gse246_fc$counts, gene_lengths)
gse246_m   <- fread(file.path(WD_DIR, "GSE246088/metadata/sample_metadata.csv"))
gse246_g2s <- fread(file.path(WD_DIR, "GSE246088/metadata/gsm_to_srr.tsv"))
gse246_mrg <- merge(gse246_m, gse246_g2s, by.x = "geo_accession", by.y = "gsm")
gse246_wd_dis  <- gse246_mrg[genotype == "Plvap_Control" & diet == "Western_Diet",  srr]
gse246_hfd_dis <- gse246_mrg[genotype == "Plvap_Control" & diet == "High_Fat_Diet", srr]

gse305_fc  <- load_fc(file.path(WD_DIR, "GSE305484/counts/featurecounts/gene_counts.txt"))
gse305_tpm <- compute_tpm(gse305_fc$counts, gene_lengths)
gse305_m   <- fread(file.path(WD_DIR, "GSE305484/metadata/sample_metadata.csv"))
gse305_g2s <- fread(file.path(WD_DIR, "GSE305484/metadata/gsm_to_srr.tsv"))
gse305_dis <- merge(gse305_m, gse305_g2s, by.x = "gsm_accession", by.y = "gsm")[condition != "Chow", srr]

wd_genes <- Reduce(intersect, list(rownames(gse220_tpm), rownames(gse246_tpm), rownames(gse305_tpm)))
wd_mat   <- cbind(
  gse220_tpm[wd_genes, intersect(colnames(gse220_tpm), gse220_dis),   drop = FALSE],
  gse246_tpm[wd_genes, intersect(colnames(gse246_tpm), gse246_wd_dis), drop = FALSE],
  gse305_tpm[wd_genes, intersect(colnames(gse305_tpm), gse305_dis),   drop = FALSE])
western_med <- apply(wd_mat, 1, median)
cat(sprintf("  Western disease samples: GSE220575=%d, GSE246088-WD=%d, GSE305484=%d\n",
            length(intersect(colnames(gse220_tpm), gse220_dis)),
            length(intersect(colnames(gse246_tpm), gse246_wd_dis)),
            length(intersect(colnames(gse305_tpm), gse305_dis))))

message("Computing TPM: HFD ...")
hfd_genes <- intersect(rownames(pub_tpm), rownames(gse246_tpm))
hfd_mat   <- cbind(
  pub_tpm[hfd_genes,   intersect(colnames(pub_tpm),   disease_meta[diet_model == "HFD", sample_id]), drop = FALSE],
  gse246_tpm[hfd_genes, intersect(colnames(gse246_tpm), gse246_hfd_dis), drop = FALSE])
hfd_med <- apply(hfd_mat, 1, median)
cat(sprintf("  HFD disease samples: main=%d, GSE246088-HFD=%d\n",
            length(intersect(colnames(pub_tpm),   disease_meta[diet_model == "HFD", sample_id])),
            length(intersect(colnames(gse246_tpm), gse246_hfd_dis))))

# Named list: diet -> per-gene median TPM vector
diet_medians <- list(MCD = mcd_med, CDAHFD = cdahfd_med, Western = western_med, HFD = hfd_med)

# =============================================================================
# 4. Assemble long-format table
# =============================================================================
message("Assembling plot data ...")
plot_rows <- lapply(OPT_LABELS, function(lbl) {
  gset <- gene_sets[[lbl]]
  lapply(DIETS, function(d) {
    med_vec <- diet_medians[[d]]
    g       <- intersect(gset, names(med_vec))
    data.table(gene = g, med_tpm = med_vec[g], diet = d, cutoff = lbl)
  }) |> rbindlist()
}) |> rbindlist()

plot_rows[, log2tpm1 := log2(med_tpm + 1)]
plot_rows[, diet   := factor(diet,   levels = DIETS)]
plot_rows[, cutoff := factor(cutoff, levels = OPT_LABELS)]

detect_dt <- plot_rows[, .(
  n       = .N,
  pct_ge1 = round(mean(med_tpm >= 1) * 100, 0)
), by = .(cutoff, diet)]

cat("\nDetection at median TPM >= 1:\n"); print(detect_dt)

# =============================================================================
# 5. Violin panel (one facet per threshold, 4 diets per facet)
# =============================================================================
diet_colors <- c(MCD = "#E64B35", CDAHFD = "#4DBBD5", Western = "#00A087", HFD = "#F39B7F")

# n-label positions: slightly below 0 on log2 scale
n_label_dt <- detect_dt[, .(
  diet, cutoff,
  label = paste0(pct_ge1, "%")
)]

panel_violin <- ggplot(plot_rows, aes(x = diet, y = log2tpm1, fill = diet)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.75, linewidth = 0.25) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               linewidth = 0.30, color = "grey20") +
  geom_hline(yintercept = log2(2), linetype = "dashed",   # TPM = 1
             linewidth = 0.28, color = "grey45") +
  geom_text(data = n_label_dt,
            aes(x = diet, y = -0.15, label = label),
            size = 1.7, color = "grey35", vjust = 1, inherit.aes = FALSE) +
  facet_wrap(~cutoff, nrow = 1) +
  scale_fill_manual(values = diet_colors, guide = "none") +
  scale_y_continuous(
    name   = "log2(median TPM + 1)",
    breaks = c(0, 2, 4, 6, 8, 10),
    expand = expansion(add = c(0.55, 0.2))
  ) +
  scale_x_discrete(name = NULL) +
  theme_masld() + theme_pub() +
  theme(
    strip.text         = element_text(size = 6.5, face = "bold"),
    strip.background   = element_rect(fill = "grey94", color = NA),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.22),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 5.5),
    axis.title.y       = element_text(size = 6)
  )

# =============================================================================
# 6. CDF panel (one line per threshold, averaged across diets)
# =============================================================================
cdf_dt <- plot_rows[, {
  thrs <- seq(0, 10, by = 0.1)
  .(thr = thrs, frac = sapply(thrs, function(t) mean(log2tpm1 >= t)))
}, by = .(cutoff)]
cdf_dt[, cutoff := factor(cutoff, levels = OPT_LABELS)]

# Assign a color per cutoff (use a sequential blue-to-red ramp)
cutoff_colors <- setNames(
  colorRampPalette(c("#1565C0", "#C9265E"))(length(OPT_THRS)),
  OPT_LABELS
)

panel_cdf <- ggplot(cdf_dt, aes(x = thr, y = frac * 100, color = cutoff)) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = log2(2), linetype = "dashed",
             linewidth = 0.28, color = "grey45") +
  scale_color_manual(values = cutoff_colors, name = "Cutoff") +
  scale_x_continuous(name = "log2(median TPM + 1)", breaks = 0:10) +
  scale_y_continuous(name = "% genes above threshold",
                     limits = c(0, 100), breaks = seq(0, 100, 25)) +
  theme_masld() + theme_pub() +
  theme(
    legend.position   = "right",
    legend.key.height = unit(0.35, "cm"),
    legend.text       = element_text(size = 5.5),
    legend.title      = element_text(size = 6),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.22)
  )

# =============================================================================
# 7. Combine and save — landscape, shorter
# =============================================================================
gene_counts_str <- paste(
  sprintf("%s: %d", OPT_LABELS, sapply(gene_sets, length)),
  collapse = " | "
)

fig <- panel_violin / panel_cdf +
  plot_annotation(
    title    = "Mouse disease liver expression of Cas13 library targets by cutoff threshold",
    subtitle = paste0(gene_counts_str,
                      "\nMedian TPM across disease samples; dashed line = TPM 1; diet averaged in CDF"),
    theme    = theme(
      plot.title    = element_text(size = 7, face = "bold"),
      plot.subtitle = element_text(size = 5, color = "grey40", lineheight = 1.3)
    )
  ) +
  plot_layout(heights = c(1.5, 1))

out_file <- file.path(OUT_DIR, "S_lib_10_tpm_by_diet.pdf")
ggsave(out_file, fig, width = 11, height = 5.5, useDingbats = FALSE)
message("Saved: ", out_file)
