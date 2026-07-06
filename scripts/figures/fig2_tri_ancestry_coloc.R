#!/usr/bin/env Rscript
# fig2_tri_ancestry_coloc.R  (2026-06-17; redesigned 2026-07-01; MVP 5-ancestry 2026-07-05)  — Fig 2 (defends para 6)
# Cross-ancestry colocalization of the pan-ancestry-portable effector genes.
# Eligibility: PP.H4 > 0.5 in >=3 of the 5 tested ancestry panels (EUR/AFR/AMR/EAS/SAS,
# 50-GWAS portfolio incl. MVP), ranked by overall convergence evidence
# (RNA-seq/results/multi_evidence/convergence_evidence.csv) so the headline genes are the
# most manuscript-relevant, not just the 3 liver-enzyme genes (RORA/EPHA2/GGT1) that
# happened to be hand-picked originally. Filename kept as *_tri_* for assembly-ref
# stability though the panel is now up to 5-ancestry.
# All values computed from disk so they always match the manuscript text.
#
# 2026-07-01 redesign: the original 3x4 tile heatmap wasted the panel on color alone
# (weakest channel for magnitude) and repeated the same 4-ancestry axis across 3
# facets. Now: PP.H4 on one shared continuous x-axis, gene on y (ranked by
# convergence importance, one axis instance total), ancestry as a dodged shape
# within each gene's row, color flags whether the colocalizing signal sits at the
# EUR discovery locus's lead variant (<10kb) or a distinct nearby signal.
#
# Out: figures/main/fig2_genetics/panels/Fig2F_rora_tri_ancestry_coloc.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

sc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
# ancestry + trait from the GWAS registry (50-GWAS portfolio incl. MVP) — NOT the
# retired grepl() heuristics. The old anc() dropped every MVP stratum into EUR (and
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

# gene rows ordered by convergence importance (top gene at the top of the plot)
m[, gene := factor(gene, levels = rev(GENES))]
m[, ancestry := factor(ancestry, levels = ANC)]
dodge <- c(EUR = 0.30, AFR = 0.15, AMR = 0.0, EAS = -0.15, SAS = -0.30)
m[, y_num := as.numeric(gene) + dodge[as.character(ancestry)]]
m[, lab := fifelse(status == "no_coloc", sprintf("%.3f", best), "")]

status_cols <- c(same_signal      = "#00695C",   # dark teal: colocalizes, EUR's own lead signal (<10kb)
                  secondary_signal = "#F5A623",  # amber (cat_palette): colocalizes, distinct signal (>=10kb)
                  no_coloc         = "#9E9E9E")  # house ns/control gray: does not colocalize
status_labs <- c(same_signal       = "same signal as EUR (<10kb)",
                  secondary_signal = "distinct signal (≥10kb)",
                  no_coloc          = "does not colocalize")
anc_shapes <- c(EUR = 21, AFR = 22, AMR = 25, EAS = 23, SAS = 24)  # all fillable (21-25)

# zebra background band per gene row (helps track 9 rows without axis clutter)
band <- data.table(y = seq_along(GENES))[y %% 2 == 0]

p <- ggplot(m, aes(x = best, y = y_num)) +
  { if (nrow(band)) geom_rect(data = band, inherit.aes = FALSE,
        aes(ymin = y - 0.5, ymax = y + 0.5), xmin = -Inf, xmax = Inf,
        fill = "grey96") } +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.35, color = "grey60") +
  geom_point(aes(fill = status, shape = ancestry), size = 1.8, color = "grey20", stroke = 0.25) +
  geom_text(aes(label = lab), hjust = -0.3, size = 2.1, color = "grey30") +
  scale_fill_manual(values = status_cols, labels = status_labs[levels(m$status)],
                    name = NULL, drop = FALSE,
                    guide = guide_legend(override.aes = list(shape = 21), order = 1)) +
  scale_shape_manual(values = anc_shapes, name = "ancestry", guide = guide_legend(order = 2)) +
  scale_x_continuous(limits = c(-0.02, 1.12), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                      expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(GENES), labels = rev(GENES),
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

save_fig(p, file.path(PANEL_DIR, "Fig2F_rora_tri_ancestry_coloc.pdf"),
         width = fig_col_width * 0.95, height = 3.15)

# Caption (house style: no in-plot title/subtitle) -> stdout
n_secondary <- m[status == "secondary_signal", .N]
n_no_coloc  <- m[status == "no_coloc", .N]
message(sprintf(
  "CAPTION (Fig2F): Cross-ancestry colocalization of the %d most convergence-relevant pan-ancestry-portable effector genes (PP.H4 > 0.5 in >=3/5 ancestry panels; ranked by overall multi-evidence convergence score). Point shape = ancestry panel (EUR/AFR/AMR/EAS/SAS across the 50-GWAS portfolio); position = max eQTL-GWAS colocalization posterior across tested liver/enzyme traits. Dark teal = colocalizes within 10kb of the EUR discovery lead variant (same causal signal); amber = colocalizes but >=10kb away (distinct secondary signal, n=%d); gray = does not colocalize (n=%d).",
  length(GENES), n_secondary, n_no_coloc))

fwrite(m[order(gene, ancestry), .(gene, ancestry, PP_H4 = round(best, 4),
        colocalizing_trait = trait, method = method_used, lead_variant = top_snp,
        dist_from_eur_kb = round(dist_kb, 2), status, colocalizes)],
       file.path(PANEL_DIR, "Fig2F_rora_tri_ancestry_coloc_source.csv"))
cat("[fig2 cross-ancestry] wrote Fig2F_rora_tri_ancestry_coloc.pdf  (genes:",
    paste(GENES, collapse = ", "), ")\n")
print(m[order(gene, ancestry), .(gene, ancestry, best = round(best, 3), status, dist_kb = round(dist_kb,1))])
