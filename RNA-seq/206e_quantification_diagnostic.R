#!/usr/bin/env Rscript
# 206e_quantification_diagnostic.R
# Diagnostic figure proving that Hoang concordance gap is driven by
# Salmon vs STAR+featureCounts quantification differences at low expression.
# Contrasts with Govaere (HT-Seq vs featureCounts) where the gap is minimal.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
INTB <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PUB  <- file.path(BASE, "data/published_degs")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_sensitivity")

cat("=== 206e: Quantification Diagnostic ===\n")

# ── Gene mapping ────────────────────────────────────────────────────────────
annot <- fread(file.path(INT, "gene_annotation/human_ensg_to_symbol.tsv"))
annot_pc <- annot[gene_type == "protein_coding"]
sym2ens <- setNames(annot_pc$gene_base, annot_pc$symbol)
extra <- annot[!gene_base %in% annot_pc$gene_base & !symbol %in% names(sym2ens)]
sym2ens <- c(sym2ens, setNames(extra[!duplicated(symbol)]$gene_base,
                                extra[!duplicated(symbol)]$symbol))

# ── Load counts + metadata ──────────────────────────────────────────────────
cat("Loading data...\n")
counts <- readRDS(file.path(INTB, "results/integration/merged_counts_raw.rds"))
meta <- readRDS(file.path(INTB, "results/integration/meta_matched.rds"))
qc <- fread(file.path(INTB, "qc/sample_qc_report.csv"))
meta <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# ── Hoang: published CLM + our CLM (same method) ───────────────────────────
cat("Loading Hoang data...\n")
hoang_fib <- as.data.table(read_excel(file.path(PUB, "GSE130970/MOESM2.xlsx"),
                                       sheet = "fibrosis ordinal regression"))
setnames(hoang_fib, names(hoang_fib), tolower(gsub("\\s+", "_", names(hoang_fib))))
hoang_fib[, gene_base := sym2ens[gene_symbol]]
clm_fib <- fread(file.path(OUTDIR, "clm_fib_results.csv"))

h <- merge(hoang_fib[!is.na(gene_base), .(gene_base, pub = range_log2fc)],
           clm_fib[!is.na(gene_base), .(gene_base, ours = range_log2FC)],
           by = "gene_base")
h <- h[is.finite(pub) & is.finite(ours)]

# Get expression levels
m130 <- meta[dataset == "GSE130970"]
idx <- colnames(counts) %in% m130$sample_id
dge130 <- DGEList(counts = counts[, idx])
dge130 <- calcNormFactors(dge130)
l2cpm130 <- cpm(dge130, log = TRUE, prior.count = 1)
expr130 <- data.table(gene_base = sub("\\.\\d+$", "", rownames(l2cpm130)),
                      mean_expr = rowMeans(l2cpm130))
h <- merge(h, expr130, by = "gene_base")

# ── Govaere: published limma + our limma (same method) ──────────────────────
cat("Loading Govaere data and computing reproduced contrast...\n")
gov <- as.data.table(read_excel(file.path(PUB,
  "GSE135251/aba4448_supplementary_tables.xlsx"), sheet = "Table S5", skip = 1))
setnames(gov, names(gov), tolower(gsub("\\s+", "_", names(gov))))
ens_col <- grep("ensembl|gene_id", names(gov), value = TRUE)[1]
lfc_col <- grep("log2fc|logfc", names(gov), value = TRUE)[1]
gov[, gene_base := sub("\\.\\d+$", "", as.character(gov[[ens_col]]))]
gov[, pub := as.numeric(gov[[lfc_col]])]

# Reproduce NASH F4 vs NAFL
m135 <- meta[dataset == "GSE135251" & !is.na(fibrosis_stage)]
m135[, grp := fcase(condition == "NAFL", "NAFL",
                     condition == "NASH_Fibrosis" & fibrosis_stage == 4, "NASH_F4",
                     default = NA_character_)]
ms <- m135[!is.na(grp)]
ms[, grp := factor(grp, levels = c("NAFL", "NASH_F4"))]
d <- model.matrix(~ grp + inferred_sex, data = ms)
idx2 <- colnames(counts) %in% ms$sample_id
dge135 <- DGEList(counts = counts[, idx2])
dge135 <- calcNormFactors(dge135)
keep <- filterByExpr(dge135, design = d)
dge135 <- dge135[keep, , keep.lib.sizes = FALSE]
v <- voom(dge135, d, plot = FALSE)
fit <- lmFit(v, d); fit <- eBayes(fit)
tt <- as.data.table(topTable(fit, coef = "grpNASH_F4", number = Inf, sort.by = "none"))
tt[, gene_base := sub("\\.\\d+$", "", rownames(topTable(fit, coef = "grpNASH_F4",
                                                         number = Inf, sort.by = "none")))]
g <- merge(gov[!is.na(gene_base), .(gene_base, pub)],
           tt[, .(gene_base, ours = logFC)], by = "gene_base")
g <- g[is.finite(pub) & is.finite(ours)]

# Govaere expression levels
l2cpm135 <- cpm(DGEList(counts = counts[, colnames(counts) %in% m135$sample_id]) |>
                  calcNormFactors(), log = TRUE, prior.count = 1)
expr135 <- data.table(gene_base = sub("\\.\\d+$", "", rownames(l2cpm135)),
                      mean_expr = rowMeans(l2cpm135))
g <- merge(g, expr135, by = "gene_base")

# ══════════════════════════════════════════════════════════════════════════
# ANALYSIS
# ══════════════════════════════════════════════════════════════════════════
cat("\nComputing expression-stratified metrics...\n")

stratify <- function(dt, n_bins = 10) {
  dt[, expr_decile := cut(mean_expr, quantile(mean_expr, 0:n_bins/n_bins),
                          include.lowest = TRUE, labels = paste0("D", 1:n_bins))]
  dt[, residual := ours - pub]
  dt[, concordant := sign(pub) == sign(ours)]

  stats <- dt[, .(
    n = .N,
    r = cor(pub, ours),
    rho = cor(pub, ours, method = "spearman"),
    sd_ratio = sd(ours) / sd(pub),
    mean_abs_res = mean(abs(residual)),
    pct_disc = 100 * mean(!concordant),
    median_expr = median(mean_expr)
  ), by = expr_decile][order(expr_decile)]
  stats
}

h_stats <- stratify(h)
g_stats <- stratify(g)

cat("\nHoang (Salmon vs featureCounts):\n")
print(h_stats[, .(expr_decile, n, r = round(r, 3), sd_ratio = round(sd_ratio, 3),
                   mean_abs_res = round(mean_abs_res, 4), pct_disc = round(pct_disc, 1))])
cat("\nGovaere (HT-Seq vs featureCounts):\n")
print(g_stats[, .(expr_decile, n, r = round(r, 3), sd_ratio = round(sd_ratio, 3),
                   mean_abs_res = round(mean_abs_res, 4), pct_disc = round(pct_disc, 1))])

# ══════════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════════
cat("\nBuilding figure...\n")

# Combine stats for plotting
h_stats[, study := "Hoang (Salmon vs featureCounts)"]
g_stats[, study := "Govaere (HT-Seq vs featureCounts)"]
both <- rbind(h_stats, g_stats)
both[, decile_num := as.integer(gsub("D", "", expr_decile))]

# --- Panel a: Hoang scatter colored by expression decile ---
h[, expr_decile := cut(mean_expr, quantile(mean_expr, 0:10/10),
                       include.lowest = TRUE, labels = paste0("D", 1:10))]
h[, decile_num := as.integer(gsub("D", "", expr_decile))]
h[, expr_group := fcase(decile_num <= 3, "Low (D1-D3)",
                         decile_num <= 7, "Mid (D4-D7)",
                         default = "High (D8-D10)")]
h[, expr_group := factor(expr_group, levels = c("Low (D1-D3)", "Mid (D4-D7)", "High (D8-D10)"))]

pa <- ggplot(h, aes(x = pub, y = ours, color = expr_group)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
  scale_color_manual(values = c("Low (D1-D3)" = masld_colors$up,
                                "Mid (D4-D7)" = masld_colors$ns,
                                "High (D8-D10)" = masld_colors$down),
                     name = "Expression") +
  labs(x = "Hoang (Salmon) range_log2FC", y = "Our (featureCounts) range_log2FC",
       title = "Hoang: Salmon vs featureCounts") +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"))

# --- Panel b: Govaere scatter colored by expression ---
g[, expr_decile := cut(mean_expr, quantile(mean_expr, 0:10/10),
                       include.lowest = TRUE, labels = paste0("D", 1:10))]
g[, decile_num := as.integer(gsub("D", "", expr_decile))]
g[, expr_group := fcase(decile_num <= 3, "Low (D1-D3)",
                         decile_num <= 7, "Mid (D4-D7)",
                         default = "High (D8-D10)")]
g[, expr_group := factor(expr_group, levels = c("Low (D1-D3)", "Mid (D4-D7)", "High (D8-D10)"))]

pb <- ggplot(g, aes(x = pub, y = ours, color = expr_group)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
  scale_color_manual(values = c("Low (D1-D3)" = masld_colors$up,
                                "Mid (D4-D7)" = masld_colors$ns,
                                "High (D8-D10)" = masld_colors$down),
                     name = "Expression") +
  labs(x = "Govaere (HT-Seq) logFC", y = "Our (featureCounts) logFC",
       title = "Govaere: HT-Seq vs featureCounts") +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"))

# --- Panel c: Pearson r by expression decile ---
pc <- ggplot(both, aes(x = decile_num, y = r, color = study)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Hoang (Salmon vs featureCounts)" = masld_colors$up,
                                "Govaere (HT-Seq vs featureCounts)" = masld_colors$down),
                     name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  labs(x = "Expression decile (D1=lowest)", y = "Pearson r",
       title = "Correlation by expression level") +
  ylim(0, 1) +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"),
                                      legend.position = "bottom")

# --- Panel d: SD ratio by decile ---
pd <- ggplot(both, aes(x = decile_num, y = sd_ratio, color = study)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 1.0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Hoang (Salmon vs featureCounts)" = masld_colors$up,
                                "Govaere (HT-Seq vs featureCounts)" = masld_colors$down),
                     name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  labs(x = "Expression decile (D1=lowest)", y = "SD ratio (ours / published)",
       title = "Variance inflation by expression") +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"),
                                      legend.position = "bottom")

# --- Panel e: Direction discordance by decile ---
pe <- ggplot(both, aes(x = decile_num, y = pct_disc, color = study)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Hoang (Salmon vs featureCounts)" = masld_colors$up,
                                "Govaere (HT-Seq vs featureCounts)" = masld_colors$down),
                     name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  labs(x = "Expression decile (D1=lowest)", y = "Direction discordance (%)",
       title = "Sign disagreement by expression") +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"),
                                      legend.position = "bottom")

# --- Panel f: Mean |residual| by decile ---
pf <- ggplot(both, aes(x = decile_num, y = mean_abs_res, color = study)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Hoang (Salmon vs featureCounts)" = masld_colors$up,
                                "Govaere (HT-Seq vs featureCounts)" = masld_colors$down),
                     name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  labs(x = "Expression decile (D1=lowest)", y = "Mean |residual|",
       title = "Absolute error by expression") +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"),
                                      legend.position = "bottom")

# --- Panel g: Hoang residual distribution by expression group ---
pg <- ggplot(h, aes(x = residual, fill = expr_group)) +
  geom_density(alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
  scale_fill_manual(values = c("Low (D1-D3)" = masld_colors$up,
                               "Mid (D4-D7)" = masld_colors$ns,
                               "High (D8-D10)" = masld_colors$down),
                    name = "Expression") +
  labs(x = "Residual (ours - published)", y = "Density",
       title = "Hoang: residual distribution") +
  xlim(-1.5, 1.5) +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"))

# --- Panel h: Cumulative r — what r would we get if we excluded low-expression genes ---
cum_r_h <- sapply(1:10, function(d) {
  sub <- h[decile_num >= d]
  if (nrow(sub) < 10) return(NA)
  cor(sub$pub, sub$ours)
})
cum_r_g <- sapply(1:10, function(d) {
  sub <- g[decile_num >= d]
  if (nrow(sub) < 10) return(NA)
  cor(sub$pub, sub$ours)
})
cum_dt <- data.table(
  decile_cutoff = rep(1:10, 2),
  r = c(cum_r_h, cum_r_g),
  study = rep(c("Hoang (Salmon vs featureCounts)",
                "Govaere (HT-Seq vs featureCounts)"), each = 10)
)

ph <- ggplot(cum_dt, aes(x = decile_cutoff, y = r, color = study)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Hoang (Salmon vs featureCounts)" = masld_colors$up,
                                "Govaere (HT-Seq vs featureCounts)" = masld_colors$down),
                     name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = paste0("≥D", 1:10)) +
  labs(x = "Include genes from decile ≥ X",
       y = "Pearson r (cumulative)",
       title = "Correlation after excluding low-expression genes") +
  ylim(0.5, 1) +
  theme_masld(base_size = 7) + theme(legend.key.size = unit(0.35, "cm"),
                                      legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ══════════════════════════════════════════════════════════════════════════
cat("Assembling figure...\n")

fig <- (pa | pb) /
       (pc | pd) /
       (pe | pf) /
       (pg | ph) +
  plot_annotation(tag_levels = "a",
                  title = "Quantification method drives concordance gap at low expression") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, file.path(OUTDIR, "quantification_diagnostic.pdf"),
         width = fig_full_width, height = 10)
cat(sprintf("  Saved %s\n", file.path(OUTDIR, "quantification_diagnostic.pdf")))

# Save stats table
fwrite(both, file.path(OUTDIR, "quantification_decile_stats.csv"))
cat("  Saved decile stats\n")

cat("\n=== Done ===\n")
