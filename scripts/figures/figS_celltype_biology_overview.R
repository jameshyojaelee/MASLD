#!/usr/bin/env Rscript
# figS_celltype_biology_overview.R
#
# Integrative overview figure for the 17-analysis cell-type-resolved MASLD
# biology pipeline (docs/superpowers/specs/2026-04-17-cell-type-resolved-
# masld-biology-design.md).
#
# Panels:
#   A — B2 bulk reverse validation summary (43.3% vs 25% null; fgsea NES)
#   B — A2 composition shift forest (Mac/Endo/T up, Hep down)
#   C — A1 primary cell-type attribution counts
#   D — G2 enrichment: GWAS-ATAC motif-disrupted ligands (OR=4.11)
#   E — F1 secretome funnel (448 LIANA -> 87 triple-concordant biomarkers)
#   F — L1 sex enrichment (Female receptors OR=2.30)
#   G — D1 Kupffer<->LAM: # LR pairs by Mac-pseudotime association
#   H — Per-cell D1_LAM score by cell type × condition (top deltas)
#   I — L3 cross-species conservation axes (Fully_Conserved MASLD-enriched)

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

# ---- load all ---------------------------------------------------------------
b2    <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/liana_bulk_concordance_perLR.csv"))
b2_fg <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/liana_fgsea_in_bulk.csv"))
a2    <- fread(file.path(BASE, "RNA-seq/results/celltype_attribution/composition_shifts.csv"))
a1    <- fread(file.path(BASE, "RNA-seq/results/celltype_attribution/celltype_primary_attribution.csv"))
g2    <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/chromatin_ccc/gwas_atac_ccc_enrichment.csv"))
f1    <- fread(file.path(BASE, "RNA-seq/results/secretome_chain/secretome_chain.csv"))
f1_trip <- fread(file.path(BASE, "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"))
l1    <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/sex_ccc/sex_ccc_enrichment.csv"))
d1    <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/macrophage_trajectory/hep_to_mac_ligands_x_pseudotime_receptors.csv"))
percell <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/percell_signatures/disease_minus_healthy_per_celltype.csv"))
l3    <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/conserved_masld_ccc_priority.csv"))

# ---- Panel A: B2 bulk reverse ----------------------------------------------
strong <- b2[abs(score_diff) >= 0.10 & !is.na(lig_lfc) & !is.na(rec_lfc)]
rates <- data.table(
  level = factor(c("Ligand concordant","Receptor concordant","Both concordant","Null (25%)"),
                 levels = c("Ligand concordant","Receptor concordant","Both concordant","Null (25%)")),
  rate  = c(mean(strong$lig_concordant), mean(strong$rec_concordant),
            mean(strong$both_concordant), 0.25)
)
b2_masld_nes <- b2_fg[pathway == "LIANA_MASLD_up_LR", NES]
pA <- ggplot(rates, aes(level, rate, fill = level)) +
  geom_col(width = 0.55, color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", rate*100)), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = c("#27AE60","#2980B9","#C0392B","grey60"), guide = "none") +
  scale_y_continuous(labels = percent_format(1), limits = c(0, .55), expand = c(0,0)) +
  labs(x = NULL, y = "Concordance rate",
       title = sprintf("B2: LIANA scRNA CCC reverse-validated in bulk\n(fgsea MASLD-up NES=%.2f, p<1e-24)", b2_masld_nes)) +
  theme_masld() + theme(axis.text.x = element_text(size = 6.5))

# ---- Panel B: A2 composition shifts ----------------------------------------
a2_mv <- a2[contrast == "Control_vs_Disease"][order(-abs(t))][1:12]
a2_mv[, celltype := factor(celltype, levels = rev(celltype))]
a2_mv[, se := abs(logit_diff / t)]
a2_mv[, ci_lo := logit_diff - 1.96 * se]
a2_mv[, ci_hi := logit_diff + 1.96 * se]
a2_mv[, sig := fifelse(padj < 0.05, "sig", "ns")]
pB <- ggplot(a2_mv, aes(logit_diff, celltype, color = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.4) +
  geom_point(size = 1.6) +
  scale_color_manual(values = c("sig" = "#C0392B","ns" = "grey60"),
                     labels = c("sig" = "padj<0.05", "ns" = "ns"), name = NULL) +
  labs(x = "Logit shift (Disease - Control)", y = NULL,
       title = "A2: Composition shifts") +
  theme_masld() + theme(axis.text.y = element_text(size = 6))

# ---- Panel C: A1 attribution counts ----------------------------------------
deg <- a1[attribution_class != "NS_bulk"]
deg[, primary_celltype_label := sub("_intrinsic_(strong|likely)$", "", attribution_class)]
deg[attribution_class %in% c("bulk_only","multi_celltype"),
    primary_celltype_label := attribution_class]
c_counts <- deg[, .N, by = primary_celltype_label][order(-N)][1:10]
c_counts[, primary_celltype_label := factor(primary_celltype_label, levels = rev(primary_celltype_label))]
pC <- ggplot(c_counts, aes(N, primary_celltype_label, fill = primary_celltype_label)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = N), hjust = -0.2, size = 2.3) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "DEGs", y = NULL, title = "A1: Bulk DEGs by primary cell type") +
  theme_masld() + theme(axis.text.y = element_text(size = 6))

# ---- Panel D: G2 enrichment ------------------------------------------------
g2_sub <- g2[grepl("ligands|receptors|DEGs", test)]
g2_sub[, test_short := gsub("enriched for ", "", test)]
g2_sub[, sig := fifelse(pvalue < 0.05, "sig", "ns")]
g2_sub[, test_short := factor(test_short, levels = rev(test_short[order(odds_ratio)]))]
pD <- ggplot(g2_sub, aes(odds_ratio, test_short, fill = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("OR=%.2f\np=%.2g", odds_ratio, pvalue)),
            hjust = -0.05, size = 2) +
  scale_fill_manual(values = c("sig" = "#C0392B","ns" = "grey65"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(x = "Odds ratio", y = NULL,
       title = "G2: GWAS-ATAC × LIANA ligand/receptor") +
  theme_masld() + theme(axis.text.y = element_text(size = 5.5))

# ---- Panel E: F1 secretome funnel ------------------------------------------
funnel <- data.table(
  stage = factor(c("LIANA\nligands","MASLD-up\n(LIANA)","In Olink\npanel","Triple\nconcordant"),
                 levels = c("LIANA\nligands","MASLD-up\n(LIANA)","In Olink\npanel","Triple\nconcordant")),
  N = c(nrow(f1),
        sum(f1$liana_direction == "MASLD_up", na.rm = TRUE),
        sum(f1$in_olink, na.rm = TRUE),
        nrow(f1_trip)))
pE <- ggplot(funnel, aes(stage, N, fill = stage)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = N), vjust = -0.3, size = 2.5) +
  scale_fill_brewer(palette = "YlOrRd", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Ligand count",
       title = "F1: Tissue → plasma biomarker funnel") +
  theme_masld() + theme(axis.text.x = element_text(size = 6.5))

# ---- Panel F: L1 sex enrichment --------------------------------------------
l1_sub <- l1[grepl("ligands|receptors", test)]
l1_sub[, sig := fifelse(pvalue < 0.05, "sig", "ns")]
l1_sub[, test_s := factor(gsub(" enriched.*", "", test),
                           levels = rev(unique(gsub(" enriched.*", "", test))))]
l1_sub[, what := ifelse(grepl("ligand", test), "ligand", "receptor")]
pF <- ggplot(l1_sub, aes(odds_ratio, test_s, fill = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("OR=%.2f\np=%.2g (%s)",
                                odds_ratio, pvalue, what)),
            hjust = -0.05, size = 2) +
  scale_fill_manual(values = c("sig" = "#E91E63","ns" = "grey65"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(x = "Odds ratio", y = NULL, title = "L1: Sex-biased DEGs × CCC") +
  theme_masld() + theme(axis.text.y = element_text(size = 6))

# ---- Panel G: D1 Kupffer/LAM LR count --------------------------------------
d1_valid <- d1[!is.na(receptor_rho) & receptor_padj < 0.05]
d1_cnt <- data.table(
  cls  = c("LAM-associated\n(receptor_rho>0)","Kupffer-associated\n(receptor_rho<0)"),
  N    = c(nrow(d1_valid[receptor_rho > 0]), nrow(d1_valid[receptor_rho < 0])))
d1_cnt[, cls := factor(cls, levels = cls)]
pG <- ggplot(d1_cnt, aes(cls, N, fill = cls)) +
  geom_col(width = 0.55, color = "white") +
  geom_text(aes(label = N), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("#C0392B","#2980B9"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Hep→Mac LR pairs",
       title = "D1: Kupffer↔LAM Hep→Mac axes") +
  theme_masld() + theme(axis.text.x = element_text(size = 7))

# ---- Panel H: per-cell top deltas ------------------------------------------
key_sigs <- c("H2_senescence_core","H2_SASP_profibrotic","H2_SenMayo",
              "D1_Kupffer","D1_LAM","D1_M1","Metab_FAO","Metab_Lipogenic",
              "Metab_Ferroptotic")
pc_top <- percell[signature %in% key_sigs][order(-abs(delta))][1:15]
pc_top[, label := paste0(cell_type, ": ", signature)]
pc_top[, label := factor(label, levels = rev(label))]
pc_top[, dir := fifelse(delta > 0, "disease-enriched","healthy-enriched")]
pH <- ggplot(pc_top, aes(delta, label, fill = dir)) +
  geom_col(width = 0.65, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  scale_fill_manual(values = c("disease-enriched" = "#C0392B",
                                "healthy-enriched" = "#2980B9"), name = NULL) +
  labs(x = "Per-cell mean delta (1.2M cells, scanpy)", y = NULL,
       title = "H2+H3+D1-deep: per-cell signature disease-delta") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5.8),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ---- Panel I: L3 conserved MASLD axes --------------------------------------
ct_pri <- l3[, .N, by = .(source, target)][order(-N)][1:10]
ct_pri[, pair := paste0(source, " -> ", target)]
ct_pri[, pair := factor(pair, levels = rev(pair))]
pI <- ggplot(ct_pri, aes(N, pair)) +
  geom_col(fill = "#27AE60", width = 0.65, color = "white") +
  geom_text(aes(label = N), hjust = -0.2, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Fully-conserved MASLD-enriched LR pairs", y = NULL,
       title = "L3: Cross-species conserved CCC (top pairs)") +
  theme_masld() + theme(axis.text.y = element_text(size = 6))

# ---- Assemble 3×3 -----------------------------------------------------------
fig <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)
fig <- fig + plot_annotation(
  tag_levels = "A",
  title = "Cell-type-resolved MASLD biology — integrative overview",
  subtitle = "17 supplementary analyses across 12 themes (A-L); see figS_celltype_biology/ for individual panels",
  theme = theme(plot.title = element_text(size = 11, face = "bold"),
                plot.subtitle = element_text(size = 8, color = "grey35"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_celltype_biology_OVERVIEW.pdf")
ggsave(out_path, fig, width = 18, height = 16)
message("Saved overview figure: ", out_path)
