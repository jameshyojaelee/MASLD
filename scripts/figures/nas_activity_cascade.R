#!/usr/bin/env Rscript

# NAS analogue of progression_cascade.R.
# Tracks: grouped DEGs vs NAS0, disease-only NMF composition, hepatocyte and
# non-parenchymal cell fractions. NAS0 NMF composition is not estimated because
# the canonical NMF was fitted on disease samples and has n=1 disease NAS0.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

GROUPS <- c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8")
# CANONICAL = treat @ lfc=0.25 (repo-wide TREAT definition; McCarthy & Smyth
# 2009). ashr / raw / sig remain as env-selectable sensitivity arms only.
# Canonical 2026-08-12: conventional padj + effect-size floor ("raw" branch).
# Was "treat" at lfc=0.25 from 2026-06-29 until the gate migration.
DEG_METHOD <- tolower(Sys.getenv("NAS_CASCADE_DEG_METHOD", "raw"))
LFC_CUTOFF <- as.numeric(Sys.getenv("NAS_CASCADE_LFC", "0.5"))
FDR_CUTOFF <- as.numeric(Sys.getenv("NAS_CASCADE_FDR", "0.05"))
stopifnot(DEG_METHOD %in% c("ashr", "raw", "treat", "sig"))
X_EXPAND <- expansion(mult = c(0.04, 0.04))
SHARED_MARGIN <- margin(3, 3, 3, 3)

bin_nas <- function(x) {
  fcase(x == 0, "NAS0", x <= 2, "NAS1-2", x <= 4, "NAS3-4",
        x <= 8, "NAS5-8", default = NA_character_)
}

# Track 1: grouped NAS contrasts vs strict NAS0.
deg <- fread(file.path(BASE, paste0(
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/",
  "disease_signatures/nas_grouped_vs_nas0_lvqw.csv")))
if (DEG_METHOD == "ashr") {
  deg[, effect_for_direction := shrunk_logFC]
  deg_sig <- deg[lfsr < FDR_CUTOFF & abs(shrunk_logFC) > LFC_CUTOFF]
} else if (DEG_METHOD == "raw") {
  deg[, effect_for_direction := logFC]
  deg_sig <- deg[padj < FDR_CUTOFF & abs(logFC) > LFC_CUTOFF]
} else if (DEG_METHOD == "sig") {
  # Significance only: ordinary moderated-t FDR with NO effect-size floor.
  # (Equivalent to treat() at lfc=0, which collapses to the plain eBayes test;
  #  LFC_CUTOFF is intentionally ignored here.) Direction from raw logFC sign.
  deg[, effect_for_direction := logFC]
  deg_sig <- deg[padj < FDR_CUTOFF]
} else {
  # Reconstruct limma::treat exactly from the saved moderated coefficient, SE,
  # t statistic and P value. df.total is common within each fitted contrast.
  infer_df <- function(t_stat, p_value) {
    ok <- which(is.finite(t_stat) & t_stat != 0 & is.finite(p_value) &
                p_value > 0 & p_value < 1)[1]
    uniroot(function(df) 2 * pt(-abs(t_stat[ok]), df) - p_value[ok],
            c(1, 1e7), tol = 1e-10)$root
  }
  deg[, treat_df := infer_df(t, P.Value), by = contrast]
  deg[, treat_p := {
    a <- abs(logFC)
    pt((a - LFC_CUTOFF) / SE, df = treat_df, lower.tail = FALSE) +
      pt((a + LFC_CUTOFF) / SE, df = treat_df, lower.tail = FALSE)
  }]
  deg[, treat_fdr := p.adjust(treat_p, method = "BH"), by = contrast]
  deg[, effect_for_direction := logFC]
  deg_sig <- deg[treat_fdr < FDR_CUTOFF]
}
deg_counts <- deg_sig[, .(up = sum(effect_for_direction > 0),
                          down = sum(effect_for_direction < 0)), by = contrast]
deg_counts[, nas_group := sub("_vs_NAS0$", "", contrast)]
deg_counts <- rbind(data.table(contrast = "reference", nas_group = "NAS0",
                               up = 0L, down = 0L), deg_counts, fill = TRUE)
deg_counts[, nas_group := factor(nas_group, levels = GROUPS)]
deg_long <- melt(deg_counts, id.vars = c("contrast", "nas_group"),
                 measure.vars = c("up", "down"),
                 variable.name = "direction", value.name = "n")
deg_long[, n_signed := fifelse(direction == "down", -n, n)]
deg_long[, direction := factor(direction, levels = c("up", "down"),
                               labels = c("Up vs NAS0", "Down vs NAS0"))]

deg_colors <- c("Up vs NAS0" = masld_colors$mash,
                "Down vs NAS0" = masld_colors$down)
p1 <- ggplot(deg_long, aes(nas_group, n_signed, fill = direction)) +
  geom_col(width = 0.66) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  geom_text(data = deg_long[n > 0],
            aes(label = n, vjust = ifelse(n_signed > 0, -0.35, 1.25)),
            size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = deg_colors, name = NULL) +
  scale_y_continuous(labels = function(x) abs(x),
                     expand = expansion(mult = c(0.12, 0.18))) +
  scale_x_discrete(drop = FALSE, expand = X_EXPAND) +
  labs(x = NULL, y = "Activity DEGs\n(vs NAS0)") +
  theme_masld_compact() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = SHARED_MARGIN)

# Track 2: existing disease-only NMF assignments. NAS0 is intentionally blank.
nmf <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
labs <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
meta <- clean_nas_meta()
meta[, nas_group := factor(bin_nas(nas_num), levels = GROUPS)]
nmf <- merge(nmf[, .(sample_id, dominant_program_code)],
             meta[, .(sample_id, nas_group)], by = "sample_id")
nmf <- nmf[nas_group != "NAS0"]
lab_map <- setNames(labs$biological_label, labs$program_code)
nmf[, prog_label := lab_map[dominant_program_code]]
nmf_comp <- nmf[, .N, by = .(nas_group, dominant_program_code, prog_label)]
nmf_comp[, frac := N / sum(N), by = nas_group]
prog_order <- c("P1", "P3", "P2", "P4", "P5", "P6")
prog_labels_ordered <- unname(lab_map[prog_order])
nmf_comp[, prog_label := factor(prog_label, levels = prog_labels_ordered)]

nmf_colors <- c(
  "Inflammatory-EMT" = cat_palette[3], "Fibrotic-ECM" = cat_palette[6],
  "Hepatocyte-Metabolic" = "#BDBDBD", "Unresolved" = "#D6D6D6",
  "Noncoding" = "#9E9E9E", "Skeletal-muscle" = "#7D7D7D")
p2 <- ggplot(nmf_comp, aes(nas_group, frac, fill = prog_label)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = nmf_colors, name = "NMF program",
                    breaks = prog_labels_ordered) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(drop = FALSE, limits = GROUPS, expand = X_EXPAND) +
  labs(x = NULL, y = "NMF program\n(% disease samples)") +
  # 2026-07-09: MUST match progression_cascade.R's fib_p2 guide exactly (ncol,
  # keyheight) -- guides="collect" only merges legends with identical guide
  # specs; a mismatch here renders two separate "NMF program" legends.
  guides(fill = guide_legend(ncol = 2, keyheight = unit(0.22, "cm"))) +
  theme_masld_compact() +
  theme(legend.position = "right", legend.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = SHARED_MARGIN)

# Tracks 3-4: all cell-proportion samples with an observed NAS score.
ct <- fread(file.path(BASE,
  "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"))
ct <- merge(ct, meta[, .(sample_id, nas_group)], by = "sample_id")
ct <- ct[!is.na(nas_group)]
ct_keep <- c("Hepatocytes", "Fibroblasts", "Macrophages",
             "T cells", "Endothelial cells")
ct_mean <- ct[, lapply(.SD, mean, na.rm = TRUE), by = nas_group, .SDcols = ct_keep]
ct_long <- melt(ct_mean, id.vars = "nas_group",
                variable.name = "celltype", value.name = "frac")
ct_long[, nas_group := factor(nas_group, levels = GROUPS)]
ct_colors <- c(
  "Hepatocytes" = "#9E9E9E", "Fibroblasts" = ct_palette[["Fibroblasts"]],
  "Macrophages" = ct_palette[["Macrophages"]], "T cells" = ct_palette[["T cells"]],
  "Endothelial cells" = ct_palette[["Endothelial cells"]])

ct_hep <- ct_long[celltype == "Hepatocytes"]
ct_npc <- ct_long[celltype != "Hepatocytes"]
p3a <- ggplot(ct_hep, aes(nas_group, frac, group = 1)) +
  geom_line(linewidth = 0.8, color = "#9E9E9E") +
  geom_point(size = 1.4, color = "#9E9E9E") +
  geom_text(data = ct_hep[nas_group %in% c("NAS0", "NAS5-8")],
            aes(label = scales::percent(frac, accuracy = 1)),
            vjust = -0.9, size = GEOM_TEXT_6PT, color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0.08, 0.14))) +
  scale_x_discrete(drop = FALSE, limits = GROUPS, expand = X_EXPAND) +
  labs(x = NULL, y = "Hepatocyte\n(%)") +
  theme_masld_compact() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = SHARED_MARGIN)

p3b <- ggplot(ct_npc, aes(nas_group, frac, color = celltype, group = celltype)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  scale_color_manual(values = ct_colors, name = "Cell type",
                     breaks = setdiff(ct_keep, "Hepatocytes")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0.05, 0.10))) +
  scale_x_discrete(drop = FALSE, limits = GROUPS, expand = X_EXPAND) +
  labs(x = "NAS activity group", y = "Non-parenchymal\n(%)") +
  # 2026-07-09: MUST match progression_cascade.R's fib_p3b guide exactly.
  guides(color = guide_legend(ncol = 2, keyheight = unit(0.22, "cm"))) +
  theme_masld_compact() +
  theme(legend.position = "right",
        legend.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain"),
        plot.margin = SHARED_MARGIN)

cascade <- p1 / p2 / p3a / p3b +
  plot_layout(heights = c(1.05, 1.05, 0.7, 0.9), guides = "collect") +
  plot_annotation(title = "Molecular and cellular remodeling across NAS activity",
                  theme = theme(plot.title = element_text(size = 6, face = "plain")))

out_name <- if (DEG_METHOD == "treat" && LFC_CUTOFF == 0.25 && FDR_CUTOFF == 0.05) {
  "figs3_nas_activity_cascade.pdf"                          # canonical
} else if (DEG_METHOD == "sig") {
  # Significance-only: no LFC floor, so omit the meaningless lfc tag.
  "figs3_nas_activity_cascade_sig.pdf"
} else {
  sprintf("figs3_nas_activity_cascade_%s_lfc%s.pdf", DEG_METHOD,
          gsub("[.]", "p", format(LFC_CUTOFF, trim = TRUE)))
}
out <- file.path(PANEL_DIR, out_name)
write_output <- tolower(Sys.getenv("CASCADE_WRITE_OUTPUT", "true")) %in%
  c("true", "1", "yes")
if (write_output) {
  save_fig(cascade, out, width = 95 / 25.4, height = 120 / 25.4)
  deg_counts[, `:=`(deg_method = DEG_METHOD,
                    lfc_cutoff = if (DEG_METHOD == "sig") NA_real_ else LFC_CUTOFF,
                    fdr_cutoff = FDR_CUTOFF)]
  fwrite(deg_counts, file.path(DATA_DIR,
                              sprintf("nas_activity_cascade_deg_counts_%s.csv", DEG_METHOD)))
  fwrite(nmf_comp, file.path(DATA_DIR, "nas_activity_cascade_nmf_composition.csv"))
  fwrite(ct_mean, file.path(DATA_DIR, "nas_activity_cascade_celltype_means.csv"))
  cat("Wrote", out, "\n")
  print(deg_counts)
}
