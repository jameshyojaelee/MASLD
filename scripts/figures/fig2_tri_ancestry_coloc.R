#!/usr/bin/env Rscript
# fig2_tri_ancestry_coloc.R  (2026-06-17; redesigned 2026-07-01; MVP 5-ancestry 2026-07-05)  — Fig 2 (defends para 6)
# Cross-ancestry colocalization of the pan-ancestry-portable effector genes.
# Scoped 2026-07-06 to the Tier-1/2 (liver-specific) MAIN strata only.
# Eligibility: PP.H4 > 0.5 in >=3 of the 5 tested ancestry panels (EUR/AFR/AMR/EAS/SAS,
# 35 Tier-1/2 GWAS incl. MVP NAFLD/ALT/AST), ranked by overall convergence evidence
# (RNA-seq/results/multi_evidence/convergence_evidence.csv) so the headline genes are the
# most manuscript-relevant, not just the 3 liver-enzyme genes (RORA/EPHA2/GGT1) that
# happened to be hand-picked originally. Renamed 2026-07-07 (Fig2F_crossancestry_coloc)
# from the stale *_rora_tri_* name: the panel shows 9 genes across 5 ancestries, so
# neither "rora" nor "tri" described it (repo grep confirmed no external ref used the
# old name). 2026-07-07 also: outline now encodes the estimator (SuSiE vs ABF) so the
# estimator/ancestry confound is visible; caption discloses EUR-only eQTL + per-ancestry
# trait heterogeneity; rows reordered by portability; "portable" overclaim softened.
# All values computed from disk so they always match the manuscript text.
#
# 2026-07-01 redesign: the original 3x4 tile heatmap wasted the panel on color alone
# (weakest channel for magnitude) and repeated the same 4-ancestry axis across 3
# facets. Now: PP.H4 on one shared continuous x-axis, gene on y (ranked by
# convergence importance, one axis instance total), ancestry as a dodged shape
# within each gene's row, color flags whether the colocalizing signal sits at the
# EUR discovery locus's lead variant (<10kb) or a distinct nearby signal.
#
# Out: figures/main/fig2_genetics/panels/Fig2F_crossancestry_coloc.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

sc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
# MAIN (Tier-1/2, liver-specific) restriction (2026-07-06): placement=="main" strata only
# (NAFLD/NASH/PDFF + ALT/AST/GGT); Tier-3/4 supp strata move to a supplementary figure.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
sc <- sc[gwas_name %in% MAIN_STUDIES]
# ancestry + trait from the GWAS registry (Tier-1/2 portfolio incl. MVP NAFLD/ALT/AST) —
# NOT the retired grepl() heuristics. The old anc() dropped every MVP stratum into EUR (and
# had no AMR bin); the old trait() collapsed ChronLiver/Cirrhosis/Albumin/Platelet
# all into "NAFLD". gwas_ancestry()/gwas_trait() are registry-driven (load_figure_data.R).
sc[, ancestry := as.character(gwas_ancestry(gwas_name))]
sc[, trait := gwas_trait(gwas_name)]
sc[, best := fifelse(!is.na(PP.H4.susie), PP.H4.susie, PP.H4.abf)]
sc[, method_used := fifelse(!is.na(PP.H4.susie), "SuSiE", "ABF")]

ANC <- GWAS_ANCESTRY_LEVELS   # c("EUR","AFR","AMR","EAS","SAS")
N_GENES <- 9

# per gene x ancestry, the TRUE max PP.H4 (any value, not pre-filtered) -- needed both
# to decide eligibility and to show real (low) values for ancestries that fail to
# colocalize, rather than blanking them out
d <- sc[ancestry %in% ANC & !is.na(best)]
best_any <- d[, .SD[which.max(best)], by = .(gene, ancestry)]

# Eligibility: PP.H4 > 0.5 in >= 3 of the 5 ancestry panels (EUR/AFR/AMR/EAS/SAS)
elig <- best_any[best > 0.5, .(n_anc = uniqueN(ancestry)), by = gene]
eligible <- elig[n_anc >= 3, gene]

# Rank eligible genes by overall multi-evidence convergence (most manuscript-relevant
# first), dropping the bottom "4_Weak" evidence tier
ce <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence.csv"),
            select = c("human_symbol", "convergence_rank", "tier"))
rank <- ce[human_symbol %in% eligible & tier != "4_Weak"]
setorder(rank, convergence_rank)
GENES <- head(rank$human_symbol, N_GENES)

# full gene x ancestry grid (keep non-colocalizing cells at their real PP.H4 so
# failures still show, rather than dropping to NA)
grid <- CJ(gene = GENES, ancestry = ANC)
m <- merge(grid, best_any[gene %in% GENES, .(gene, ancestry, best, trait, method_used, top_snp)],
           by = c("gene", "ancestry"), all.x = TRUE)
m[, colocalizes := !is.na(best) & best > 0.5]

# EUR = each gene's discovery locus; classify other ancestries by distance from the
# EUR lead variant. <10kb ~ same LD block / credible set (still the same causal
# signal); >=10kb = a distinct secondary signal at the same gene.
eur_lead <- m[ancestry == "EUR", .(gene, eur_snp = top_snp)]
m <- merge(m, eur_lead, by = "gene")
m[, dist_kb := fifelse(!is.na(top_snp) & !is.na(eur_snp),
                        abs(as.numeric(sub(".*:", "", top_snp)) -
                            as.numeric(sub(".*:", "", eur_snp))) / 1000, NA_real_)]
m[, status := fifelse(!colocalizes, "no_coloc",
             fifelse(dist_kb < 10, "same_signal", "secondary_signal"))]
m[, status := factor(status, levels = c("same_signal", "secondary_signal", "no_coloc"))]

# Selection is convergence-ranked (above); DISPLAY order = portability so the visual
# message (how many ancestries replicate) reads top-to-bottom. Order by n ancestries
# colocalizing, then mean PP.H4; keep convergence_rank as the final tie-break.
crank <- rank[, .(human_symbol, convergence_rank)]
port <- m[, .(n_coloc = sum(colocalizes), mean_best = mean(best, na.rm = TRUE)), by = gene]
port <- merge(port, crank, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
setorder(port, -n_coloc, -mean_best, convergence_rank)
GENE_ORDER <- as.character(port$gene)
m[, gene := factor(gene, levels = rev(GENE_ORDER))]
m[, ancestry := factor(ancestry, levels = ANC)]
m[, method_used := factor(method_used, levels = c("SuSiE", "ABF"))]  # outline aes
dodge <- c(EUR = 0.30, AFR = 0.15, AMR = 0.0, EAS = -0.15, SAS = -0.30)
m[, y_num := as.numeric(gene) + dodge[as.character(ancestry)]]
# No per-point value labels: x-position already encodes PP.H4 for every point, so a
# floating number just restates the coordinate, and labeling only a couple of points
# would falsely flag them as special. The 0.5 gridline shows near-misses directly.

status_cols <- c(same_signal      = "#00695C",   # dark teal: colocalizes, EUR's own lead signal (<10kb)
                  secondary_signal = "#F5A623",  # amber (cat_palette): colocalizes, distinct signal (>=10kb)
                  no_coloc         = "#9E9E9E")  # house ns/control gray: does not colocalize
# Fill = genomic proximity of the ancestry's lead variant to the EUR discovery lead.
# Wording is deliberately "lead variant" not "signal": the colocalizing TRAIT can differ
# across ancestries at the same variant (e.g. RORA ALT in EUR vs GGT in EAS/SAS), so
# "same variant" must NOT be read as "same trait" (trait is in the source CSV + caption).
# Compact labels so the 3-item status legend fits the 3.26in width without clipping;
# the full "same lead variant as the EUR discovery lead" wording lives in the caption.
status_labs <- c(same_signal      = "same variant (<10 kb)",
                  secondary_signal = "distinct (≥10 kb)",
                  no_coloc          = "no coloc")
anc_shapes <- c(EUR = 21, AFR = 22, AMR = 25, EAS = 23, SAS = 24)  # all fillable (21-25)
# Outline color encodes the estimator. The estimator is unevenly distributed across
# ancestries (AFR/SAS almost entirely ABF; EUR/EAS mostly SuSiE), and SuSiE PP.H4
# (multi-causal) vs ABF PP.H4 (single-causal) are different estimands -- making the
# estimator visible stops "AFR colocalizes less" from being read as pure ancestry.
method_cols <- c(SuSiE = "grey15", ABF = "#C0C0C0")

# zebra background band per gene row (helps track 9 rows without axis clutter)
band <- data.table(y = seq_along(GENES))[y %% 2 == 0]

p <- ggplot(m, aes(x = best, y = y_num)) +
  { if (nrow(band)) geom_rect(data = band, inherit.aes = FALSE,
        aes(ymin = y - 0.5, ymax = y + 0.5), xmin = -Inf, xmax = Inf,
        fill = "grey96") } +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.35, color = "grey60") +
  geom_point(aes(fill = status, shape = ancestry, colour = method_used),
             size = 1.8, stroke = 0.5) +
  scale_fill_manual(values = status_cols, labels = status_labs[levels(m$status)],
                    name = NULL, drop = FALSE,
                    guide = guide_legend(override.aes = list(shape = 21, colour = "grey15"), order = 1)) +
  scale_shape_manual(values = anc_shapes, name = "ancestry",
                     guide = guide_legend(override.aes = list(fill = "grey70", colour = "grey15"), order = 2)) +
  scale_colour_manual(values = method_cols, name = "estimator",
                      guide = guide_legend(override.aes = list(shape = 21, fill = "grey85"), order = 3)) +
  scale_x_continuous(limits = c(-0.02, 1.12), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                      expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(GENES), labels = rev(GENE_ORDER),
                      limits = c(0.5, length(GENES) + 0.5), expand = c(0, 0)) +
  labs(x = "Colocalization posterior (PP.H4)", y = NULL) +
  theme_masld() + theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 6, face = "plain"),
        axis.text.x  = element_text(size = 6, face = "plain"),
        axis.text.y  = element_text(size = 6, face = "italic"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.spacing = unit(0.12, "cm"),
        legend.spacing.y = unit(0.12, "cm"),
        legend.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        legend.key.size = unit(0.24, "cm"),
        legend.margin = margin(t = 2, b = 0))

source(file.path(BASE, "figures/layout_specs/regenerate_panels.R"))   # save_panel(): exact contract size + cairo_pdf
save_panel(p, "main/fig2_genetics/panels/Fig2E_crossancestry_coloc.pdf",
           read_sizes(file.path(BASE, "figures/layout_specs/figure2_panel_sizes.tsv")),
           file.path(BASE, "figures"))

# Caption (house style: no in-plot title/subtitle) -> stdout
n_secondary <- m[status == "secondary_signal", .N]
n_no_coloc  <- m[status == "no_coloc", .N]
n_max_anc   <- max(port$n_coloc)   # honest ceiling: no gene is 5/5
message(sprintf(
  "CAPTION (Fig2E): Cross-ancestry colocalization of the %d top convergence-ranked effector genes that show cross-ancestry portability (eQTL-GWAS PP.H4 > 0.5 in >=3 of 5 ancestry panels). Portability is PARTIAL: %d-%d of 5 panels replicate; no gene colocalizes in all five, so this is a replicated-subset panel, not a pan-ancestry claim. Point shape = ancestry (EUR/AFR/AMR/EAS/SAS across the 35 Tier-1/2 liver-specific GWAS); x = max PP.H4 across the tested liver traits (ALT/AST/GGT/NAFLD). OUTLINE = estimator (dark = SuSiE multi-causal; light gray = ABF single-causal): the estimator is unevenly distributed across ancestries (AFR/SAS almost all ABF), so apparent ancestry differences partly reflect the estimator. eQTL panel = EUR liver (Broadaway, N=1,183) for ALL ancestries -- no ancestry-matched liver eQTL exists -- so non-EUR cells test whether the EUR-defined effector's signal ports to the non-EUR GWAS, NOT ancestry-matched colocalization. FILL = genomic proximity of the ancestry's lead variant to the EUR discovery lead (dark teal <10 kb = same variant; amber >=10 kb = distinct signal, n=%d; gray = does not colocalize, n=%d). NOTE the colocalizing trait can differ across ancestries at the SAME variant (e.g. RORA: ALT in EUR, GGT in EAS/SAS) -- 'same variant' does not imply 'same trait' (per-cell trait in the source CSV column 'colocalizing_trait'). Rows ordered by number of replicating ancestries.",
  length(GENES), min(port$n_coloc), n_max_anc, n_secondary, n_no_coloc))

fwrite(m[order(gene, ancestry), .(gene, ancestry, PP_H4 = round(best, 4),
        colocalizing_trait = trait, method = method_used, lead_variant = top_snp,
        dist_from_eur_kb = round(dist_kb, 2), status, colocalizes)],
       file.path(PANEL_DIR, "Fig2E_crossancestry_coloc_source.csv"))
cat("[fig2 cross-ancestry] wrote Fig2E_crossancestry_coloc.pdf  (genes, portability order:",
    paste(GENE_ORDER, collapse = ", "), ")\n")
print(m[order(gene, ancestry), .(gene, ancestry, best = round(best, 3), status, dist_kb = round(dist_kb,1))])
