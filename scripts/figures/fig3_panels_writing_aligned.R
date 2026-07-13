#!/usr/bin/env Rscript
##############################################################################
# fig3_panels_writing_aligned.R
# Journal-quality individual panels for Fig 3, in writing order.
# No titles / subtitles — data speaks through design.
#
# Panels:
#   scatter_coloc_vs_deg.pdf    — PP.H4 vs logFC
#   atac_rora_thrb.pdf          — GWAS-ATAC motif disruption
#   cross_ancestry_labeled.pdf  — EUR × EAS concordance
#   cyp26a1_locus.pdf           — CYP26A1 EAS locus zoom
#   regulon_tf_lollipop.pdf     — disease-regulon TF lollipop
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(ggrastr)
  library(patchwork)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
ATAC_DIR  <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")

save_panel <- function(p, name, width = fig_half_width, height = 3.2) {
  out <- file.path(PANEL_DIR, name)
  save_fig(p, out, width = width, height = height)
  cat("[fig3] Saved:", basename(out), "\n")
}

# Shared palette across all panels
COL_THRB   <- "#880E4F"   # dark magenta  — THRB / drug target
COL_RORA   <- "#00695C"   # deep teal     — RORA
COL_HNF4A  <- "#1565C0"   # deep blue     — HNF4A
COL_SHARED <- "#E6A817"   # amber         — shared variants
COL_EUR    <- "#3B7DA5"   # steel blue    — European
COL_EAS    <- "#E07B39"   # burnt orange  — East Asian
COL_BOTH   <- "#7B3294"   # purple        — cross-ancestry
COL_MUTED  <- "#90A4AE"   # blue-gray     — secondary


# ═══════════════════════════════════════════════════════════════════════════
# 1. COLOC PP.H4 vs RNA-seq logFC
#    Visual story: THRB + RORA sit in the top half (genetically supported +
#    transcriptionally suppressed). PPARG / GLP1R / PPARA are dysregulated
#    but have no GWAS-eQTL support — two distinct drug-target phenotypes.
# ═══════════════════════════════════════════════════════════════════════════
cat("[1] Scatter: COLOC PP.H4 vs logFC ...\n")

atlas <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol","bulk_logFC","bulk_padj",
             "coloc_susie_best_pp4","coloc_abf_best_pp4"))
setnames(atlas, "human_symbol", "gene")
atlas[, pp4 := fcoalesce(coloc_susie_best_pp4, coloc_abf_best_pp4)]
atlas[is.na(pp4), pp4 := 0]

# Gene sets
key_genes <- c("THRB","RORA","PNPLA3","NR1H4","PPARA","PPARG","GLP1R",
               "SCD","HSD17B13","DGAT2","CYP26A1","GOT2","ADH4","ALDH2","MARC1")

gene_role <- function(g) fcase(
  g %in% c("THRB","RORA"),                          "causal_supported",
  g %in% c("PPARG","GLP1R","PPARA","NR1H4","SCD",
            "HSD17B13","DGAT2","PNPLA3"),            "drug_unsupported",
  g %in% c("CYP26A1","GOT2","ADH4","ALDH2","MARC1"), "novel_coloc",
  default                                            = "other"
)

fg <- atlas[gene %in% key_genes]
fg[, role := gene_role(gene)]
bg <- atlas[!gene %in% key_genes & !is.na(bulk_padj) & bulk_padj < 0.05]
bg[, pp4 := fcoalesce(coloc_susie_best_pp4, coloc_abf_best_pp4)]
bg[is.na(pp4), pp4 := 0]

role_cols <- c(
  causal_supported  = COL_THRB,
  drug_unsupported  = COL_MUTED,
  novel_coloc       = COL_RORA,
  other             = "gray82"
)
role_sizes <- c(causal_supported = 2.8, drug_unsupported = 2.0,
                novel_coloc = 2.2, other = 1.0)

p1 <- ggplot() +
  # Subtle shading: PP4 > 0.5 region
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.05,
           fill = "#FFF8E1", alpha = 0.6) +
  # Reference lines
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, color = "gray55") +
  geom_hline(yintercept = 0.9, linetype = "dotted",
             linewidth = 0.25, color = "gray65") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  # Background DEG cloud
  rasterize(
    geom_point(data = bg, aes(x = bulk_logFC, y = pp4),
               color = "gray85", size = 0.3, alpha = 0.5, shape = 16),
    dpi = 600
  ) +
  # Annotated genes
  geom_point(data = fg[role == "other"],
             aes(x = bulk_logFC, y = pp4),
             color = role_cols["other"], size = role_sizes["other"],
             shape = 16, alpha = 0.7) +
  geom_point(data = fg[role == "drug_unsupported"],
             aes(x = bulk_logFC, y = pp4),
             color = role_cols["drug_unsupported"],
             size = role_sizes["drug_unsupported"], shape = 16) +
  geom_point(data = fg[role == "novel_coloc"],
             aes(x = bulk_logFC, y = pp4),
             color = role_cols["novel_coloc"],
             size = role_sizes["novel_coloc"], shape = 16) +
  geom_point(data = fg[role == "causal_supported"],
             aes(x = bulk_logFC, y = pp4),
             color = role_cols["causal_supported"],
             size = role_sizes["causal_supported"], shape = 16) +
  # Labels — key genes only
  geom_text_repel(
    data = fg[role %in% c("causal_supported","drug_unsupported","novel_coloc")],
    aes(x = bulk_logFC, y = pp4, label = gene,
        color = role),
    size = GEOM_TEXT_6PT, fontface = "italic",
    box.padding = 0.4, point.padding = 0.2,
    segment.size = 0.22, segment.color = "gray55",
    min.segment.length = 0, max.overlaps = Inf,
    seed = 42, force = 5, bg.color = "white", bg.r = 0.1,
    show.legend = FALSE
  ) +
  # Axis annotations (replace subtitle)
  annotate("text", x = Inf, y = 0.515, label = "PP.H4 = 0.5",
           hjust = 1.05, size = GEOM_TEXT_6PT, color = "gray50", fontface = "plain") +
  annotate("text", x = Inf, y = 0.915, label = "PP.H4 = 0.9",
           hjust = 1.05, size = GEOM_TEXT_6PT, color = "gray55", fontface = "plain") +
  scale_color_manual(values = role_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_y_continuous(limits = c(-0.01, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = expression(log[2]~"fold change (MASLD vs. control)"),
       y = "SuSiE colocalization PP.H4") +
  theme_masld() +
  theme(panel.border = element_rect(color = "gray70", fill = NA,
                                    linewidth = 0.3))

# RETIRED 2026-07-07 (not a Fig 2 / FigS2 panel — stale leftover in fig2_genetics/panels/):
# scatter_coloc_vs_deg.pdf. Panel p1 kept computed above for provenance; save disabled so it
# is never regenerated. The other panels in this script (GWAS-ATAC, cross-ancestry, etc.) are
# unaffected.
# save_panel(p1, "scatter_coloc_vs_deg.pdf",
#            width = fig_half_width + 0.4, height = 3.6)


# ═══════════════════════════════════════════════════════════════════════════
# 2. GWAS-ATAC: PIP vs motif disruption strength
#    Visual story: RORA and THRB loci carry high-PIP variants that
#    strongly perturb their own (and each other's) binding motifs.
#    One shared chr5 variant disrupts both simultaneously.
# ═══════════════════════════════════════════════════════════════════════════
cat("[2] GWAS-ATAC motif disruption scatter ...\n")

motif_raw <- fread(file.path(ATAC_DIR, "motif_disruption_scores.csv"))
ms <- unique(motif_raw[, .(SNP_id, tf_name, alleleDiff, max_pip,
                            motif_in_disease_regulon)])

va <- fread(file.path(ATAC_DIR, "gwas_atac_variant_annotation.csv"),
            select = c("variant_id","nearest_gene"))
va <- va[, .(nearest_gene = nearest_gene[1L]), by = variant_id]

merged <- merge(ms, va, by.x = "SNP_id", by.y = "variant_id", all.x = TRUE)
credible <- merged[max_pip >= 0.2]
credible[, abs_diff := abs(alleleDiff)]
credible[, gene_label := fcoalesce(nearest_gene, "")]

rora_snps   <- unique(credible[tf_name == "RORA"]$SNP_id)
thrb_snps   <- unique(credible[tf_name == "THRB"]$SNP_id)
shared_snps <- intersect(rora_snps, thrb_snps)

# Per-variant summary (max disruption across TF hits)
per_var <- credible[, .(
  abs_diff  = max(abs_diff),
  max_pip   = first(max_pip),
  gene_label = first(gene_label),
  in_regulon = any(motif_in_disease_regulon)
), by = SNP_id]

per_var[, type := fcase(
  SNP_id %in% shared_snps, "shared",
  SNP_id %in% rora_snps,   "RORA",
  SNP_id %in% thrb_snps,   "THRB",
  in_regulon == TRUE,       "regulon",
  default                   = "other"
)]
per_var[, type := factor(type, levels = c("shared","RORA","THRB","regulon","other"))]

type_cols  <- c(shared = COL_SHARED, RORA = COL_RORA, THRB = COL_THRB,
                regulon = COL_MUTED, other = "gray82")
type_sizes <- c(shared = 3.2, RORA = 2.4, THRB = 2.4, regulon = 1.6, other = 0.9)
type_alpha <- c(shared = 1.0, RORA = 0.9, THRB = 0.9, regulon = 0.7, other = 0.4)

# Labels for named loci
label_dt <- per_var[type %in% c("RORA","THRB","shared")]
label_dt[, label := gene_label]
# deduplicate same gene at similar positions
label_dt <- label_dt[!duplicated(label)]

p2 <- ggplot(per_var, aes(x = max_pip, y = abs_diff)) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.25, color = "gray65") +
  geom_vline(xintercept = 0.8, linetype = "dotted",
             linewidth = 0.22, color = "gray70") +
  # Plot layers from bottom to top
  geom_point(data = per_var[type == "other"],
             color = type_cols["other"], size = type_sizes["other"],
             alpha = type_alpha["other"], shape = 16) +
  geom_point(data = per_var[type == "regulon"],
             color = type_cols["regulon"], size = type_sizes["regulon"],
             alpha = type_alpha["regulon"], shape = 16) +
  geom_point(data = per_var[type == "RORA"],
             color = type_cols["RORA"], size = type_sizes["RORA"],
             alpha = type_alpha["RORA"], shape = 16) +
  geom_point(data = per_var[type == "THRB"],
             color = type_cols["THRB"], size = type_sizes["THRB"],
             alpha = type_alpha["THRB"], shape = 16) +
  geom_point(data = per_var[type == "shared"],
             color = type_cols["shared"], size = type_sizes["shared"],
             alpha = type_alpha["shared"], shape = 18) +
  geom_text_repel(
    data = label_dt,
    aes(label = label, color = type),
    size = GEOM_TEXT_6PT, fontface = "italic",
    segment.size = 0.2, min.segment.length = 0.05,
    box.padding = 0.3, point.padding = 0.15,
    max.overlaps = Inf, seed = 42, bg.color = "white", bg.r = 0.1,
    show.legend = FALSE
  ) +
  # Inline legend annotations
  annotate("point", x = 0.98, y = max(per_var$abs_diff) * 0.98,
           color = COL_RORA, size = 2.4, shape = 16) +
  annotate("text",  x = 0.95, y = max(per_var$abs_diff) * 0.98,
           label = "RORA motif", hjust = 1, size = GEOM_TEXT_6PT,
           color = COL_RORA, fontface = "plain") +
  annotate("point", x = 0.98, y = max(per_var$abs_diff) * 0.90,
           color = COL_THRB, size = 2.4, shape = 16) +
  annotate("text",  x = 0.95, y = max(per_var$abs_diff) * 0.90,
           label = "THRB motif", hjust = 1, size = GEOM_TEXT_6PT,
           color = COL_THRB, fontface = "plain") +
  annotate("point", x = 0.98, y = max(per_var$abs_diff) * 0.82,
           color = COL_SHARED, size = 3.0, shape = 18) +
  annotate("text",  x = 0.95, y = max(per_var$abs_diff) * 0.82,
           label = "RORA + THRB", hjust = 1, size = GEOM_TEXT_6PT,
           color = COL_SHARED, fontface = "italic") +
  scale_color_manual(values = type_cols, guide = "none") +
  scale_x_continuous(limits = c(0.18, 1.05),
                     breaks = c(0.2, 0.5, 0.8, 1.0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "Max credible-set PIP",
       y = "|Motif alleleDiff|") +
  theme_masld()

# RETIRED 2026-06-12 (not a Fig 2 panel): atac_rora_thrb.pdf
# save_panel(p2, "atac_rora_thrb.pdf",
#            width = fig_half_width + 0.6, height = 3.2)


# ═══════════════════════════════════════════════════════════════════════════
# 3. Cross-ancestry colocalization: EUR × EAS
#    Visual story: corner of gold = shared regulatory architecture across
#    populations. RORA also replicates in South Asian GWAS.
# ═══════════════════════════════════════════════════════════════════════════
cat("[3] Cross-ancestry EUR × EAS scatter ...\n")

EUR_17 <- c("2019_31311600_NAFLD_EUR","2020_32298765_NAFLD_EUR",
            "2021_34128465_PDFF_EUR","2021_34841290_NAFLD_EUR",
            "2021_34957434_PDFF_EUR","2022_36402844_PDFF_EUR",
            "2023_36280732_NAFLD_deCode_EUR","2023_36280732_NAFLD_Intermountain_EUR",
            "2023_36280732_NAFLD_UKBB_EUR","FinnGen_NAFLD","FinnGen_NASH",
            "UKBB_ALT","UKBB_AST","UKBB_GGT")
BBJ_5  <- c("BBJ_ALT","BBJ_AST","BBJ_GGT")
SAS_3  <- c("PanUKBB_CSA_ALT","PanUKBB_CSA_AST","PanUKBB_CSA_GGT")

sc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  select = c("gene","gwas_name","PP.H4.abf","PP.H4.susie"))
mhc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
  select = c("gene","is_mhc"))
sc <- merge(sc, mhc, by = "gene", all.x = TRUE)
sc <- sc[is.na(is_mhc) | is_mhc == FALSE]; sc[, is_mhc := NULL]

anc <- function(g) fcase(g %in% EUR_17,"EUR", g %in% BBJ_5,"EAS",
                          g %in% SAS_3,"SAS", default=NA_character_)
sc[, ancestry := anc(gwas_name)]
sc[, pp4 := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]

per_anc <- sc[!is.na(ancestry) & !is.na(pp4),
              .(max_pp4 = max(pp4)), by = .(gene, ancestry)]
gene_xa <- dcast(per_anc[ancestry %in% c("EUR","EAS")],
                 gene ~ ancestry, value.var = "max_pp4", fill = 0)
n_phen  <- sc[!is.na(pp4) & pp4 >= 0.5, .(n_phen = uniqueN(gwas_name)), by = gene]
gene_xa <- merge(gene_xa, n_phen, by = "gene", all.x = TRUE)
gene_xa[is.na(n_phen), n_phen := 0L]

gene_xa[, tier := fcase(EUR >= 0.5 & EAS >= 0.5, "both",
                         EUR >= 0.5, "EUR-only",
                         EAS >= 0.5, "EAS-only",
                         default   = "bg")]
n_both <- sum(gene_xa$tier == "both")

# RORA SAS value for annotation
rora_sas <- max(per_anc[gene == "RORA" & ancestry == "SAS"]$max_pp4, na.rm = TRUE)

# Jitter at ceiling
gene_xa[EUR > 0.995, EUR := 1 - runif(.N, 0, 0.016)]
gene_xa[EAS > 0.995, EAS := 1 - runif(.N, 0, 0.016)]

label_genes <- c("RORA","THRB","GOT2","CYP26A1","ADH4","ALDH2","HNF1A","MARC1","ALDOB")
label_dt <- gene_xa[gene %in% label_genes & tier != "bg"]

tier_cols <- c(both = COL_BOTH, "EUR-only" = COL_EUR,
               "EAS-only" = COL_EAS, bg = "#C8C2BC")

p3 <- ggplot() +
  # Gold quadrant
  annotate("rect", xmin = 0.5, xmax = 1.02, ymin = 0.5, ymax = 1.02,
           fill = "#FFF8E1", alpha = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.22, color = "gray60") +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.22, color = "gray60") +
  # Count annotation in the corner
  annotate("text", x = 0.76, y = 0.54,
           label = sprintf("n = %d genes", n_both),
           size = GEOM_TEXT_6PT, color = COL_BOTH, fontface = "plain") +
  # RORA SAS callout
  annotate("text", x = 1.01, y = 0.535,
           label = sprintf("RORA: SAS PP.H4 = %.3f", rora_sas),
           size = GEOM_TEXT_6PT, color = COL_RORA, hjust = 1, fontface = "plain") +
  # Points
  rasterize(
    geom_point(data = gene_xa[tier == "bg"],
               aes(x = EUR, y = EAS),
               color = "#C8C2BC", size = 0.35, alpha = 0.2, stroke = 0),
    dpi = 600
  ) +
  geom_point(data = gene_xa[tier %in% c("EUR-only","EAS-only")],
             aes(x = EUR, y = EAS, color = tier),
             size = 1.2, alpha = 0.75, stroke = 0) +
  geom_point(data = gene_xa[tier == "both"],
             aes(x = EUR, y = EAS, size = n_phen),
             color = COL_BOTH, alpha = 0.85, stroke = 0) +
  # White ring on "both" dots for crispness
  geom_point(data = gene_xa[tier == "both"],
             aes(x = EUR, y = EAS, size = n_phen),
             color = "white", shape = 1, stroke = 0.5, alpha = 0.4) +
  geom_label_repel(
    data = label_dt,
    aes(x = EUR, y = EAS, label = gene),
    color = "gray15", fill = alpha("white", 0.90),
    size = GEOM_TEXT_6PT, fontface = "italic",
    label.size = 0.15, label.padding = unit(0.09, "lines"),
    label.r = unit(0.06, "lines"),
    box.padding = 0.5, point.padding = 0.25,
    segment.size = 0.18, segment.color = "gray45",
    min.segment.length = 0, force = 5,
    max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = tier_cols, guide = "none") +
  scale_size_continuous(range = c(1.1, 3.2), name = "GWAS\n(PP.H4≥0.5)",
                        breaks = c(1,3,6,10)) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1),
                     expand = c(0, 0)) +
  labs(x = "Max PP.H4  (14 EUR GWAS)",
       y = "Max PP.H4  (3 EAS GWAS)") +
  theme_masld() +
  theme(panel.border = element_rect(color = "gray55", fill = NA,
                                    linewidth = 0.35),
        legend.position  = c(0.08, 0.78),
        legend.background = element_blank(),
        legend.key.size  = unit(0.22, "cm"),
        legend.title = element_text(size = 6, face = "plain"),
        legend.text  = element_text(size = 6)) +
  guides(size = guide_legend(override.aes = list(color = COL_BOTH)))

# RETIRED 2026-06-12 (no longer a Fig 2 panel): cross_ancestry_labeled.pdf
# save_panel(p3, "cross_ancestry_labeled.pdf",
#            width = fig_half_width + 0.8, height = 3.6)

# Sidecar CSV
fwrite(gene_xa[tier == "both"][order(-(EUR + EAS)),
               .(gene, EUR_max_pp4 = round(EUR,4),
                 EAS_max_pp4 = round(EAS,4), n_phen)],
       file.path(FIG3_DIR, "cross_ancestry_pp4_table_aligned.csv"))


# ═══════════════════════════════════════════════════════════════════════════
# 4. CYP26A1 EAS-specific locus (BBJ GGT, chr10)
#    Visual story: this East Asian-specific signal is invisible in EUR data —
#    ancestry-matched LD references unlock it.
# ═══════════════════════════════════════════════════════════════════════════
cat("[4] CYP26A1 locus zoom ...\n")

CYP_CHR   <- 10L
CYP_LEAD  <- 94839724L
CYP_WIN   <- 600000L

bbj <- fread(
  file.path(BASE, "GWAS/finemapping/data/sumstats/BBJ_GGT_reformatted_hg19.tsv"),
  select = c("chromosome","position","pval"))
lz  <- bbj[chromosome == CYP_CHR &
            position  >= CYP_LEAD - CYP_WIN &
            position  <= CYP_LEAD + CYP_WIN]
lz[, mlog10p := -log10(pmax(as.numeric(pval), 1e-300))]

pip_dt <- fread(file.path(ATAC_DIR, "gwas_atac_variant_annotation.csv"),
                select = c("chromosome","position","max_pip",
                           "nearest_gene","in_susie_cs"))
cyp_pip <- pip_dt[chromosome == CYP_CHR &
                  position  >= CYP_LEAD - CYP_WIN &
                  position  <= CYP_LEAD + CYP_WIN]

lead <- lz[which.max(mlog10p)]
GWS  <- -log10(5e-8)

# Top GWAS panel
p4a <- ggplot(lz, aes(x = position, y = mlog10p)) +
  rasterize(
    geom_point(color = "gray78", size = 0.5, alpha = 0.55, shape = 16),
    dpi = 600
  ) +
  geom_hline(yintercept = GWS, linetype = "longdash",
             linewidth = 0.3, color = "gray30") +
  geom_point(data = lead, fill = COL_EAS, color = "gray10",
             size = 2.4, shape = 21, stroke = 0.4) +
  annotate("text", x = CYP_LEAD, y = GWS + 0.3,
           label = expression(italic(CYP26A1)), size = GEOM_TEXT_6PT,
           color = COL_EAS, fontface = "plain") +
  scale_x_continuous(limits = c(CYP_LEAD - CYP_WIN, CYP_LEAD + CYP_WIN),
                     labels = function(x) paste0(round(x/1e6, 1), " Mb"),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(y = expression(-log[10](italic(p)))) +
  theme_masld() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank())

# PIP panel
p4b <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             linewidth = 0.25, color = "gray60") +
  {if (nrow(cyp_pip) > 0)
    geom_segment(data = cyp_pip,
                 aes(x = position, xend = position, y = 0, yend = max_pip,
                     color = in_susie_cs),
                 linewidth = 0.5, alpha = 0.65)
  } +
  {if (nrow(cyp_pip) > 0)
    geom_point(data = cyp_pip,
               aes(x = position, y = max_pip, color = in_susie_cs),
               size = 1.3, alpha = 0.85)
  } +
  scale_color_manual(values = c("TRUE" = COL_EAS, "FALSE" = "gray75"),
                     na.value = "gray75", guide = "none") +
  scale_x_continuous(limits = c(CYP_LEAD - CYP_WIN, CYP_LEAD + CYP_WIN),
                     labels = function(x) paste0(round(x/1e6, 1), " Mb"),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1),
                     expand = c(0, 0)) +
  labs(x = sprintf("chr%d (hg19)", CYP_CHR), y = "PIP") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6))

p4 <- p4a / p4b + plot_layout(heights = c(2.4, 1))

# RETIRED 2026-06-12 (no longer a Fig 2 panel): cyp26a1_locus.pdf
# save_panel(p4, "cyp26a1_locus.pdf",
#            width = fig_half_width, height = 3.4)


# ═══════════════════════════════════════════════════════════════════════════
# 5. Disease-regulon TF lollipop
#    Visual story: 15 SCENIC+ disease-active TFs carry disrupting GWAS
#    credible-set variants. HNF4A, RORA, THRB and KLF15 top the list.
# ═══════════════════════════════════════════════════════════════════════════
cat("[5] Disease-regulon TF lollipop ...\n")

mot <- fread(file.path(ATAC_DIR, "motif_disruption_scores.csv"))
enrich <- fread(file.path(ATAC_DIR, "enrichment_statistics.csv"))
hep_e  <- enrich[grepl("Hepatocyte", cell_type)]

reg_tfs <- mot[motif_in_disease_regulon == TRUE,
               .(n_vars    = uniqueN(SNP_id),
                 n_vars_hi = uniqueN(SNP_id[max_pip >= 0.5])),
               by = tf_name]
setorder(reg_tfs, n_vars)
reg_tfs[, tf_name := factor(tf_name, levels = tf_name)]

reg_tfs[, color := fcase(
  tf_name == "HNF4A", COL_HNF4A,
  tf_name == "RORA",  COL_RORA,
  tf_name == "THRB",  COL_THRB,
  default             = COL_MUTED
)]

p5 <- ggplot(reg_tfs, aes(y = tf_name, x = n_vars)) +
  # Stem
  geom_segment(aes(x = 0, xend = n_vars, yend = tf_name,
                   color = color),
               linewidth = 0.7, alpha = 0.8) +
  # High-PIP inner ring
  geom_point(aes(x = n_vars_hi, color = color),
             shape = 21, fill = "white", size = 2.8, stroke = 0.5,
             alpha = 0.5) +
  # Main dot
  geom_point(aes(color = color), size = 3.8, shape = 16) +
  # Count label
  geom_text(aes(label = n_vars, x = n_vars + 0.08),
            hjust = 0, size = GEOM_TEXT_6PT, fontface = "plain",
            color = "gray25") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, max(reg_tfs$n_vars) + 1.2),
                     breaks = 0:max(reg_tfs$n_vars),
                     expand = c(0, 0)) +
  labs(x = "Credible-set variants disrupting TF motif  (n)",
       y = NULL) +
  theme_masld() +
  theme(
    axis.text.y = element_text(
      size = 6, face = "italic",
      color = ifelse(levels(reg_tfs$tf_name) %in% c("HNF4A","RORA","THRB"),
                     "gray10", "gray40")
    ),
    panel.grid.major.x = element_line(color = "gray93", linewidth = 0.3),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(4, 8, 4, 4)
  )

# RETIRED 2026-06-12 (no longer a Fig 2 panel): regulon_tf_lollipop.pdf
# save_panel(p5, "regulon_tf_lollipop.pdf",
#            width = fig_half_width + 0.4, height = 3.2)


cat("\n[fig3] Done. Panels in:", PANEL_DIR, "\n")
