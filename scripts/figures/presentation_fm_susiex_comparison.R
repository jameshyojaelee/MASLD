#!/usr/bin/env Rscript
# presentation_fm_susiex_comparison.R
# Presentation-grade "1KG ≈ UKBB" comparison for standalone fine-mapping (SuSiE)
# and cross-ancestry SuSiEX joint fine-mapping.
#
# Outputs (figures/presentation/):
#   fm_concordance_presentation.pdf        — SuSiE PIP scatter + CS-size comparison
#   fm_concordance_scatter.pdf             — SuSiE PIP scatter alone
#   fm_concordance_cs_size.pdf             — CS-size bar
#   susiex_concordance_presentation.pdf    — SuSiEX PIP scatter + CS-count comparison
#   susiex_concordance_scatter.pdf         — SuSiEX PIP scatter alone
#   fm_susiex_summary.txt                  — plain-text stats

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

FM_DIR <- file.path(PROJ, "GWAS/finemapping")
# Consolidated into the 3-way folder (UKBB-vs-1KG PIP scatter + SuSiEX top-PIP
# scatter are already produced by presentation_3way_full.R as facets; this
# script now only emits the unique 2-way deliverables: CS-size scatter,
# reproducibility bars, top-SNP-agreement bar, summary text).
PROJ_OUT <- file.path(PROJ, "figures/supplementary/figS04_coloc/ld_panel_3way_eur")
OUT_DIR <- PROJ_OUT
dir.create(file.path(OUT_DIR, "fm"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "susiex"), recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Standalone SuSiE fine-mapping (Phase 4b)
#    - 1KG EUR: output/{gwas}/1kg_eur_0.5Mb/susie/
#    - UKBB v1: output/{gwas}/EUR_0.5Mb/susie/
#    Both pipelines produce *.tsv per locus with (chromosome, position, PIP, CS, ...)
# ===========================================================================

# Studies to include (17 EUR GWAS)
registry <- fread(file.path(FM_DIR, "config/gwas_registry.tsv"))
eur_studies <- registry[ancestry == "EUR"]$study_name

load_fm <- function(study, panel_dir) {
  susie_dir <- file.path(FM_DIR, "output", study, panel_dir, "susie")
  if (!dir.exists(susie_dir)) return(NULL)
  files <- list.files(susie_dir, pattern = "\\.tsv$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    dt <- try(fread(f), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    dt[, study := study]
    dt[, locus_id := gsub("_cov0\\.95_.*", "", sub(".*/", "", f))]
    # Don't trust PIP from non-converged (locus_id has _notconverged suffix is a filename thing)
    dt[, converged_flag := !grepl("_notconverged", basename(f))]
    dt
  })
  dt <- rbindlist(rows, fill = TRUE)
  dt[]
}

cat("Loading SuSiE fine-mapping results...\n")
fm_1kg <- rbindlist(lapply(eur_studies, load_fm, panel_dir = "1kg_eur_0.5Mb"), fill = TRUE)
fm_v1  <- rbindlist(lapply(eur_studies, load_fm, panel_dir = "EUR_0.5Mb"),  fill = TRUE)
cat("  1KG EUR:  ", nrow(fm_1kg), " variant-locus rows\n", sep = "")
cat("  UKBB v1:  ", nrow(fm_v1),  " variant-locus rows\n", sep = "")

# Merge by (study, locus, variant) — match per-variant PIP
fm_merge <- merge(
  fm_1kg[, .(study, locus_id, SNP, PIP_1kg = PIP, CS_1kg = CS,
             conv_1kg = converged_flag)],
  fm_v1[, .(study, locus_id, SNP, PIP_v1 = PIP, CS_v1 = CS,
            conv_v1 = converged_flag)],
  by = c("study", "locus_id", "SNP"),
  all = FALSE
)
# Keep only converged pairs
fm_merge <- fm_merge[conv_1kg & conv_v1 & !is.na(PIP_1kg) & !is.na(PIP_v1)]
cat("  paired variant rows (both converged): ", nrow(fm_merge), "\n\n", sep = "")

pearson_fm  <- cor(fm_merge$PIP_1kg, fm_merge$PIP_v1)
spearman_fm <- cor(fm_merge$PIP_1kg, fm_merge$PIP_v1, method = "spearman")

# High-PIP agreement
cs_agree   <- nrow(fm_merge[PIP_1kg > 0.5 & PIP_v1 > 0.5])
only_1kg   <- nrow(fm_merge[PIP_1kg > 0.5 & PIP_v1 <= 0.5])
only_v1    <- nrow(fm_merge[PIP_1kg <= 0.5 & PIP_v1 > 0.5])
cat(sprintf("SuSiE PIP: r = %.3f, ρ = %.3f\n", pearson_fm, spearman_fm))
cat(sprintf("  PIP>0.5: both = %d, 1KG-only = %d, UKBB-only = %d\n\n",
            cs_agree, only_1kg, only_v1))

# Per-variant PIP scatter (UKBB v1 vs 1KG) is now produced by
# presentation_3way_full.R as the UKBB-vs-1KG facet of fm/scatter.pdf — not
# duplicated here.

# --- Panel B: PIP>0.5 reproducibility bar ---
bar_fm <- data.frame(
  category = factor(c("Both panels agree", "UKBB only", "1KG only"),
                    levels = c("Both panels agree", "UKBB only", "1KG only")),
  n = c(cs_agree, only_v1, only_1kg)
)
p_fm_b <- ggplot(bar_fm, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), vjust = -0.3, size = 9, fontface = "plain") +
  scale_fill_manual(values = c("#2CA02C", "#FF7F0E", "#9467BD"), guide = "none") +
  labs(x = NULL, y = "# variants at PIP > 0.5",
       title = sprintf("High-PIP variant reproducibility (%.1f%%)",
                       100 * cs_agree / (cs_agree + only_v1))) +
  theme_masld(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "plain"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        panel.grid.major.x = element_blank()) +
  expand_limits(y = max(bar_fm$n) * 1.15)

# --- Panel C: CS-size comparison (per locus) ---
cs_size <- fm_merge[, .(
  cs_size_1kg = sum(CS_1kg > 0, na.rm = TRUE),
  cs_size_v1  = sum(CS_v1  > 0, na.rm = TRUE)
), by = .(study, locus_id)]
cs_size <- cs_size[cs_size_1kg > 0 | cs_size_v1 > 0]

p_fm_c <- ggplot(cs_size, aes(x = cs_size_v1, y = cs_size_1kg)) +
  geom_abline(slope = 1, intercept = 0, color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
  geom_point(alpha = 0.5, size = 3, color = "#8c564b") +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x = 1.2, y = max(cs_size$cs_size_1kg) * 0.8,
           label = sprintf("n = %d loci", nrow(cs_size)),
           size = 7, hjust = 0) +
  labs(x = "CS size — UKBB EUR", y = "CS size — 1KG EUR",
       title = "Credible set sizes per locus") +
  theme_masld(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "plain"),
        axis.title = element_text(size = 22),
        axis.text  = element_text(size = 18))

fm_unique <- (p_fm_b | p_fm_c) +
  plot_annotation(
    title = "FM — UKBB vs 1KG (CS size + PIP>0.5 reproducibility)",
    theme = theme(plot.title = element_text(size = 24, face = "plain"))
  )

ggsave(file.path(OUT_DIR, "fm/cs_size_and_reproducibility.pdf"), fm_unique,
       width = 18, height = 8, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "fm/cs_size.pdf"), p_fm_c,
       width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "fm/pip05_reproducibility_bar.pdf"), p_fm_b,
       width = 8, height = 7, device = cairo_pdf)

cat("Wrote: cs_size_and_reproducibility.pdf + 2 sub-panels\n\n")

# ===========================================================================
# 2. SuSiEX cross-ancestry joint fine-mapping (Phase 4c)
#    140 shared loci (UKBB EUR × BBJ EAS) for ALT/AST/GGT
#    Compare MAX_PIP + CS_LENGTH per locus between v1 and 1kg outputs.
# ===========================================================================
cat("Loading SuSiEX results...\n")
load_susiex <- function(root) {
  if (!dir.exists(root)) return(NULL)
  files <- list.files(root, pattern = "\\.summary$", recursive = TRUE, full.names = TRUE)
  rows <- lapply(files, function(f) {
    ln <- try(readLines(f), silent = TRUE)
    if (inherits(ln, "try-error") || length(ln) < 2) return(NULL)
    # skip comment lines
    data_lines <- ln[!grepl("^#", ln) & nchar(ln) > 0]
    if (length(data_lines) < 2) return(NULL)
    dt <- try(fread(text = paste(data_lines, collapse = "\n")), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    locus_id <- sub("\\.summary$", "", basename(f))
    trait    <- basename(dirname(f))
    dt[, locus_id := locus_id]
    dt[, trait := trait]
    dt
  })
  if (length(rows) == 0) return(NULL)
  rbindlist(rows, fill = TRUE)
}

sx_v1  <- load_susiex(file.path(FM_DIR, "results/susiex"))
sx_1kg <- load_susiex(file.path(FM_DIR, "results/susiex_1kg"))
cat("  UKBB v1:  ", nrow(sx_v1),  " CS rows across ",
    uniqueN(sx_v1$locus_id),  " loci\n", sep = "")
cat("  1KG EUR:  ", nrow(sx_1kg), " CS rows across ",
    uniqueN(sx_1kg$locus_id), " loci\n", sep = "")

# Aggregate per locus → top MAX_PIP + n_CS
agg_sx <- function(dt) {
  dt[, .(top_pip = max(MAX_PIP, na.rm = TRUE),
         n_cs   = .N,
         total_cs_length = sum(CS_LENGTH, na.rm = TRUE),
         top_pip_snp = MAX_PIP_SNP[which.max(MAX_PIP)]),
     by = .(trait, locus_id)]
}
ag_v1  <- agg_sx(sx_v1)
ag_1kg <- agg_sx(sx_1kg)

ag_m <- merge(
  ag_v1[,  .(trait, locus_id, top_pip_v1 = top_pip, n_cs_v1 = n_cs,
             total_cs_len_v1 = total_cs_length, top_snp_v1 = top_pip_snp)],
  ag_1kg[, .(trait, locus_id, top_pip_1kg = top_pip, n_cs_1kg = n_cs,
             total_cs_len_1kg = total_cs_length, top_snp_1kg = top_pip_snp)],
  by = c("trait", "locus_id"))
cat("  paired loci: ", nrow(ag_m), "\n\n", sep = "")

pearson_sx  <- cor(ag_m$top_pip_1kg, ag_m$top_pip_v1)
spearman_sx <- cor(ag_m$top_pip_1kg, ag_m$top_pip_v1, method = "spearman")
same_top_snp <- sum(ag_m$top_snp_v1 == ag_m$top_snp_1kg, na.rm = TRUE)
cat(sprintf("SuSiEX top MAX_PIP: r = %.3f, ρ = %.3f; same top-PIP SNP: %d / %d (%.1f%%)\n\n",
            pearson_sx, spearman_sx,
            same_top_snp, nrow(ag_m),
            100 * same_top_snp / nrow(ag_m)))

# SuSiEX top-PIP scatter (UKBB v1 vs 1KG) is now produced by
# presentation_3way_full.R as the UKBB-vs-1KG facet of susiex/scatter.pdf —
# not duplicated here.

# --- SuSiEX Panel B: CS-count comparison ---
cs_cmp <- data.table(
  category = factor(c("Same top-PIP SNP", "Different"),
                    levels = c("Same top-PIP SNP", "Different")),
  n = c(same_top_snp, nrow(ag_m) - same_top_snp)
)
p_sx_b <- ggplot(cs_cmp, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), vjust = -0.3, size = 9, fontface = "plain") +
  scale_fill_manual(values = c("#2CA02C", "#9467BD"), guide = "none") +
  labs(x = NULL, y = "# loci",
       title = sprintf("Top-PIP SNP reproducibility (%.1f%%)",
                       100 * same_top_snp / nrow(ag_m))) +
  theme_masld(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "plain"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        panel.grid.major.x = element_blank()) +
  expand_limits(y = max(cs_cmp$n) * 1.15)

# --- SuSiEX Panel C: CS-count scatter ---
p_sx_c <- ggplot(ag_m, aes(x = n_cs_v1, y = n_cs_1kg, color = trait)) +
  geom_abline(slope = 1, intercept = 0, color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
  geom_jitter(alpha = 0.7, size = 3, width = 0.1, height = 0.1) +
  scale_color_manual(values = c(ALT = "#1F77B4", AST = "#2CA02C", GGT = "#D62728"),
                     guide = "none") +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  scale_y_continuous(breaks = seq(1, 10, 1)) +
  labs(x = "# CS — UKBB EUR", y = "# CS — 1KG EUR",
       title = "Credible-set counts per locus") +
  theme_masld(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "plain"),
        axis.title = element_text(size = 22),
        axis.text  = element_text(size = 18))

sx_unique <- (p_sx_b | p_sx_c) +
  plot_annotation(
    title = "SuSiEX — UKBB vs 1KG (top-SNP agreement + CS counts)",
    theme = theme(plot.title = element_text(size = 24, face = "plain"))
  )

ggsave(file.path(OUT_DIR, "susiex/top_snp_agreement_and_cs_counts.pdf"), sx_unique,
       width = 18, height = 8, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "susiex/top_snp_agreement_bar.pdf"), p_sx_b,
       width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "susiex/cs_counts.pdf"), p_sx_c,
       width = 9, height = 7, device = cairo_pdf)

cat("Wrote: top_snp_agreement_and_cs_counts.pdf + 2 sub-panels\n\n")

# ===========================================================================
# Summary txt
# ===========================================================================
sink(file.path(OUT_DIR, "fm_susiex_summary.txt"))
cat("Fine-mapping + SuSiEX concordance: UKBB v1 (sghatan, N≈337K) vs 1KG EUR (N=379)\n")
cat("Generated: ", as.character(Sys.time()), "\n\n", sep = "")
cat("=== Standalone SuSiE fine-mapping (Phase 4b) ===\n")
cat("Paired variant rows (both converged): ", format(nrow(fm_merge), big.mark = ","), "\n", sep = "")
cat("Pearson r:  ", sprintf("%.4f", pearson_fm),  "\n", sep = "")
cat("Spearman ρ: ", sprintf("%.4f", spearman_fm), "\n", sep = "")
cat("At PIP > 0.5:\n")
cat("  Both agree: ", cs_agree,  "\n", sep = "")
cat("  UKBB-only:  ", only_v1,   "\n", sep = "")
cat("  1KG-only:   ", only_1kg,  "\n", sep = "")
cat("  Reproducibility: ", sprintf("%.1f%%", 100*cs_agree/(cs_agree+only_v1)), "\n\n", sep = "")

cat("=== SuSiEX cross-ancestry joint fine-mapping (Phase 4c) ===\n")
cat("Paired loci: ", nrow(ag_m), "\n", sep = "")
cat("Pearson r (top PIP):  ", sprintf("%.4f", pearson_sx),  "\n", sep = "")
cat("Spearman ρ (top PIP): ", sprintf("%.4f", spearman_sx), "\n", sep = "")
cat("Same top-PIP SNP: ",  same_top_snp, " / ", nrow(ag_m),
    " (", sprintf("%.1f%%", 100*same_top_snp/nrow(ag_m)), ")\n", sep = "")
sink()

cat("All outputs in: ", OUT_DIR, "\n")
