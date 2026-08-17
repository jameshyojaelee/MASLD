#!/usr/bin/env Rscript
# KEY MESSAGE: Cross-ancestry portability is partial and depends on which signal
# model is evaluable; absent multi-signal output is not a negative COLOC result.
# fig2_tri_ancestry_coloc.R  (2026-06-17; redesigned 2026-07-01; MVP 5-ancestry 2026-07-05)  — Fig 2 (defends para 6)
# Cross-ancestry colocalization of selected, partially portable candidate genes.
# Scoped 2026-07-06 to the Tier-1/2 (liver-specific) MAIN strata only.
# Eligibility: PP.H4 > 0.5 in >=3 of the 5 tested ancestry panels (EUR/AFR/AMR/EAS/SAS,
# 35 Tier-1/2 GWAS incl. MVP NAFLD/ALT/AST). Rows are ranked using genetics alone:
# number of ancestry panels with multi-signal support, total supported panels, mean
# multi-signal PP.H4, then gene symbol. The artwork encodes the four evidence/evaluability
# states; exact methods, traits, lead variants, LD, and ancestry provenance stay in the
# source table and legend.
# All values computed from disk so they always match the manuscript text.
#
# 2026-08-12 revision: genetics-only gene-by-ancestry support matrix. Cell fill
# reports multi-signal support, single-signal-only support, evaluated/no support,
# or not evaluable. Exact posteriors, traits, lead variants, and ancestry/LD
# provenance remain in the source table rather than competing in the artwork.
#
# Out: figures/main/fig2_genetics/panels/Fig2F_crossancestry_coloc.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

COLOC_INPUT <- Sys.getenv(
  "FIG2_COLOC_INPUT",
  file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
FROZEN_SOURCE <- Sys.getenv("FIG2_CROSSANCESTRY_SOURCE")
OUT_DIR <- Sys.getenv("FIG2_CANDIDATE_DIR", file.path(FIG3_DIR, "panels"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ANC <- GWAS_ANCESTRY_LEVELS   # c("EUR","AFR","AMR","EAS","SAS")
N_GENES <- 9

if (nzchar(FROZEN_SOURCE)) {
  m <- fread(FROZEN_SOURCE)
  required <- c("gene", "ancestry", "evidence_state", "multi_pp4", "single_pp4")
  stopifnot(all(required %in% names(m)))
  GENE_ORDER <- unique(m$gene)
  stopifnot(length(GENE_ORDER) == N_GENES, nrow(m) == N_GENES * length(ANC))
} else {
sc <- fread(COLOC_INPUT)
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

# Collapse each gene-by-ancestry cell separately for the two signal models. The
# displayed state follows a fixed hierarchy, but both maxima and their provenance
# are retained in the source table.
d <- sc[ancestry %in% ANC]
best_susie <- d[is.finite(PP.H4.susie), .SD[which.max(PP.H4.susie)],
  by = .(gene, ancestry)][, .(
    gene, ancestry, multi_pp4 = PP.H4.susie, multi_study = gwas_name,
    multi_trait = trait, multi_lead_variant = top_snp,
    multi_ld_panel = ld_panel, multi_ld_panel_n = ld_panel_n,
    multi_ld_reliability = ld_reliability)]
best_abf <- d[is.finite(PP.H4.abf), .SD[which.max(PP.H4.abf)],
  by = .(gene, ancestry)][, .(
    gene, ancestry, single_pp4 = PP.H4.abf, single_study = gwas_name,
    single_trait = trait, single_lead_variant = top_snp,
    single_ld_panel = ld_panel, single_ld_panel_n = ld_panel_n,
    single_ld_reliability = ld_reliability)]
cells <- merge(best_susie, best_abf, by = c("gene", "ancestry"), all = TRUE)
cells[, evidence_state := fcase(
  is.finite(multi_pp4) & multi_pp4 > 0.5, "multi_signal",
  is.finite(single_pp4) & single_pp4 > 0.5, "single_signal_only",
  is.finite(multi_pp4) | is.finite(single_pp4), "evaluated_no_support",
  default = "not_evaluable")]

# Genetics-only row selection. Require support in at least three ancestry panels,
# then rank by multi-signal breadth, total support, mean multi-signal posterior,
# and gene symbol. No disease-state or convergence output enters this decision.
rank <- cells[, .(
  n_multi = sum(evidence_state == "multi_signal"),
  n_supported = sum(evidence_state %in% c("multi_signal", "single_signal_only")),
  mean_multi_pp4 = if (any(evidence_state == "multi_signal"))
    mean(multi_pp4[evidence_state == "multi_signal"]) else NA_real_), by = gene]
rank <- rank[n_supported >= 3]
setorder(rank, -n_multi, -n_supported, -mean_multi_pp4, gene, na.last = TRUE)
GENES <- head(rank$gene, N_GENES)
if (!length(GENES)) stop("No genes meet the prespecified >=3-ancestry support rule")

grid <- CJ(gene = GENES, ancestry = ANC)
m <- merge(grid, cells[gene %in% GENES], by = c("gene", "ancestry"), all.x = TRUE)
m[is.na(evidence_state), evidence_state := "not_evaluable"]
port <- rank[gene %in% GENES]
setorder(port, -n_multi, -n_supported, -mean_multi_pp4, gene, na.last = TRUE)
GENE_ORDER <- as.character(port$gene)
}
GENES <- GENE_ORDER
m[, gene := factor(gene, levels = rev(GENE_ORDER))]
m[, ancestry := factor(ancestry, levels = ANC)]
# The main artwork shows only the evidence state. Exact posterior, trait,
# method availability, lead variant, LD source, and ancestry provenance remain
# in the source table below.
state_levels <- c("multi_signal", "single_signal_only",
                  "evaluated_no_support", "not_evaluable")
m[, evidence_state := factor(evidence_state, levels = state_levels)]
m[, gene := factor(as.character(gene), levels = rev(GENE_ORDER))]
m[, ancestry := factor(as.character(ancestry), levels = ANC)]
state_cols <- c(
  multi_signal = "#1565C0",
  single_signal_only = "#64B5F6",
  evaluated_no_support = "#B0BEC5",
  not_evaluable = "white")
state_labs <- c(
  multi_signal = "Multi-signal COLOC",
  single_signal_only = "Single-signal COLOC only",
  evaluated_no_support = "Evaluated, PP.H4 ≤ 0.5",
  not_evaluable = "□ Not evaluable")

# Circle area reports the posterior used for the displayed evidence state.
# For evaluated/no-support cells, use the larger available posterior so the
# gray mark shows the strongest evaluated model without implying support.
m[, display_pp4 := fcase(
  evidence_state == "multi_signal", multi_pp4,
  evidence_state == "single_signal_only", single_pp4,
  evidence_state == "evaluated_no_support",
    pmax(multi_pp4, single_pp4, na.rm = TRUE),
  default = NA_real_)]
m[!is.finite(display_pp4), display_pp4 := NA_real_]
m[, display_method := fcase(
  evidence_state == "multi_signal", "multi_signal",
  evidence_state == "single_signal_only", "single_signal",
  evidence_state == "evaluated_no_support", "maximum_available",
  default = "not_evaluable")]
stopifnot(all(is.na(m$display_pp4) | between(m$display_pp4, 0, 1)))

# Diagonal marks distinguish not-evaluable cells from supported or evaluated
# states without making absence look like a negative result.
hatch <- m[evidence_state == "not_evaluable", .(
  gene, ancestry, x = as.numeric(ancestry), y = as.numeric(gene))]

p <- ggplot(m, aes(x = ancestry, y = gene)) +
  geom_tile(width = 0.92, height = 0.88, fill = "white",
            colour = "#E0E0E0", linewidth = 0.25) +
  geom_point(data = m[is.finite(display_pp4)],
             aes(size = display_pp4, fill = evidence_state),
             shape = 21, colour = "#4D4D4D", stroke = 0.25) +
  geom_segment(data = hatch, inherit.aes = FALSE,
               aes(x = x - 0.34, xend = x + 0.34,
                   y = y - 0.32, yend = y + 0.32),
               colour = "#9E9E9E", linewidth = 0.25) +
  scale_fill_manual(values = state_cols, labels = state_labs, name = NULL,
                    drop = TRUE,
                    guide = guide_legend(
                      order = 1, ncol = 1,
                      override.aes = list(size = 2.5))) +
  scale_size_area(name = "PP.H4", max_size = 3.6, limits = c(0, 1),
                  breaks = c(0.25, 0.50, 0.75, 1.00),
                  labels = c("0.25", "0.50", "0.75", "1.00"),
                  guide = guide_legend(order = 2, nrow = 1)) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.08)) +
  scale_y_discrete(expand = expansion(add = 0.08)) +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(panel.grid = element_blank(),
        axis.text.x  = element_text(size = 6, face = "plain", colour = "black"),
        axis.text.y  = element_text(size = 6, face = "italic"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        legend.key.width = unit(0.22, "cm"),
        legend.key.height = unit(0.16, "cm"),
        legend.spacing.y = unit(0.01, "cm"),
        legend.margin = margin(t = 1, b = 0),
        plot.margin = margin(2, 3, 2, 3))

source(file.path(BASE, "figures/layout_specs/regenerate_panels.R"))   # save_panel(): exact contract size + cairo_pdf
out_pdf <- file.path(OUT_DIR, "Fig2F_crossancestry_coloc.pdf")
ggsave(out_pdf, p, width = 2.76, height = 2.44, units = "in",
       device = cairo_pdf, family = "Helvetica")

# Caption (house style: no in-plot title/subtitle) -> stdout
n_supported_by_gene <- m[, .(n = sum(evidence_state %in%
  c("multi_signal", "single_signal_only"))), by = gene]
message(sprintf(paste0(
  "CAPTION (Fig2F): Signal-model-aware portability of EUR liver-eQTL colocalization across five GWAS ancestry panels. Rows are the %d genes selected using genetics alone: number of ancestry panels with multi-signal support, total supported panels, mean multi-signal PP.H4, and gene symbol. Circle area denotes the displayed PP.H4. Dark blue uses multi-signal SuSiE PP.H4; light blue uses single-signal ABF PP.H4 when it alone exceeds 0.5; gray uses the larger available posterior when neither model exceeds 0.5. Portability is partial: %d-%d of %d panels support each displayed gene. Exact PP.H4, trait, method availability, lead variant, LD source, GWAS ancestry, and eQTL ancestry are retained in the source table. The eQTL resource is European in every column, so non-European GWAS columns test portability of a European-defined regulatory relationship rather than ancestry-matched molecular regulation. COLOC nominates shared-signal candidates; it does not establish mediation or biological differences between ancestry groups."),
  length(GENES), min(n_supported_by_gene$n), max(n_supported_by_gene$n), length(ANC)))

source_out <- copy(m)
source_out[, gene := as.character(gene)]
source_out[, ancestry := as.character(ancestry)]
source_out[, evidence_state := as.character(evidence_state)]
source_out[, `:=`(
  gwas_ancestry = ancestry,
  eqtl_ancestry = "EUR",
  eqtl_source = "Broadaway liver eQTL (1,183 donors)",
  multi_signal_available = is.finite(multi_pp4),
  single_signal_available = is.finite(single_pp4)
)]
setcolorder(source_out, c(
  "gene", "ancestry", "gwas_ancestry", "eqtl_ancestry", "evidence_state",
  "display_pp4", "display_method",
  "multi_pp4", "multi_study", "multi_trait", "multi_lead_variant",
  "multi_ld_panel", "multi_ld_panel_n", "multi_ld_reliability",
  "single_pp4", "single_study", "single_trait", "single_lead_variant",
  "single_ld_panel", "single_ld_panel_n", "single_ld_reliability"))
source_out <- source_out[order(match(gene, GENE_ORDER), match(ancestry, ANC))]
fwrite(source_out, file.path(OUT_DIR, "Fig2F_crossancestry_coloc_source.tsv"), sep = "\t")
fwrite(source_out, file.path(OUT_DIR, "Fig2F_crossancestry_coloc_source.csv"))
cat("[fig2 cross-ancestry] wrote Fig2F_crossancestry_coloc.pdf  (genes, portability order:",
    paste(GENE_ORDER, collapse = ", "), ")\n")
print(m[order(gene, ancestry), .(gene, ancestry, evidence_state,
                                 multi_pp4 = round(multi_pp4, 3),
                                 single_pp4 = round(single_pp4, 3))])
