#!/usr/bin/env Rscript
# figS_integrated_vs_perstudy_venn.R
# Cas13-library view of the integrated-vs-per-study DEG overlap, restricted to
# CAS13-ELIGIBLE targets, i.e. the library's targeting filters applied to each
# human DEG set:
#   (1) UPREGULATED in disease  : logFC > 0.5  (positive sign + Tier-1 magnitude)
#   (2) biotype                 : protein_coding or lncRNA (human GENCODE v49)
#   (3) mouse-ortholog hepatic expression (the screen target is the MOUSE
#       ortholog): map human -> mouse via the library's master_ortholog_table
#       (confidence tier H/M, one2one-preferred), then require mouse in-house
#       MCD mean TPM >= 1 (protein_coding) / >= 0.5 (lncRNA).
#
# Same 3-panel area-proportional Euler layout as
# figS_methods_validation/integration_value/panels/integrated_vs_perstudy_venn.pdf,
# but every set (integrated + 5 cohorts) is the Cas13-eligible subset.
#
# Output: Cas13_Library_Design/figures/integrated_vs_perstudy_venn.pdf  (FIGS_CAS13LIB_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(FIGS_CAS13LIB_DIR, recursive = TRUE, showWarnings = FALSE)

COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
COHORT_SHORT <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
                  GSE162694 = "GSE162694", GSE213621 = "GSE213621")
PADJ_THR  <- 0.05
LFC_THR   <- 0.5                       # UP: logFC > +0.5
TPM_PC    <- 1.0                       # protein_coding mouse-ortholog MCD TPM gate
TPM_LNC   <- 0.5                       # lncRNA           mouse-ortholog MCD TPM gate
strip_v   <- function(x) sub("\\..*", "", as.character(x))

# ---------------------------------------------------------------------------
# 1. Human biotype (GENCODE v49) — classify protein_coding vs lncRNA
# ---------------------------------------------------------------------------
meta <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"),
              select = c("ensembl_base", "gene_biotype"))
meta <- unique(meta, by = "ensembl_base")
meta[, biotype := fcase(gene_biotype == "protein_coding", "protein_coding",
                        gene_biotype == "lncRNA",         "lncRNA",
                        default = NA_character_)]
biotype_of <- setNames(meta$biotype, meta$ensembl_base)

# ---------------------------------------------------------------------------
# 2. Human -> mouse ortholog (library master_ortholog_table; H/M, one2one-pref)
# ---------------------------------------------------------------------------
ortho <- fread(cmd = paste0("zcat ", file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")),
               select = c("mouse_ensembl", "human_ensembl", "confidence_tier", "is_one2one"))
ortho <- ortho[confidence_tier %in% c("H", "M")]
ortho[, `:=`(hb = strip_v(human_ensembl), mb = strip_v(mouse_ensembl))]
ortho[, trank := match(confidence_tier, c("H", "M"))]
ortho[, o2o_rank := ifelse(is_one2one %in% c(TRUE, "True", "TRUE", "true"), 0L, 1L)]
# best mouse per human (deterministic: tier, one2one, then ensembl id)
setorder(ortho, hb, trank, o2o_rank, mb)
h2m <- unique(ortho, by = "hb")[, .(hb, mb)]
mouse_of <- setNames(h2m$mb, h2m$hb)

# ---------------------------------------------------------------------------
# 3. Mouse in-house MCD hepatic mean TPM (the library expression reference)
# ---------------------------------------------------------------------------
mtpm <- fread(file.path(BASE, "RNA-seq/Mouse/InHouse_MCD/results/mean_tpm_mcd.csv"),
              select = c("gene_id", "mean_tpm"))
mtpm[, mb := strip_v(gene_id)]
mtpm <- mtpm[, .(mean_tpm = max(mean_tpm)), by = mb]   # collapse any dup base ids
tpm_of <- setNames(mtpm$mean_tpm, mtpm$mb)

# ---------------------------------------------------------------------------
# Cas13-eligibility for a vector of human ENSG (base): UP is handled per-set;
# here = biotype PC/lncRNA AND mouse ortholog AND mouse MCD TPM >= biotype gate.
# ---------------------------------------------------------------------------
cas13_eligible <- function(genes_base) {
  bt  <- biotype_of[genes_base]
  mb  <- mouse_of[genes_base]
  tpm <- tpm_of[mb]
  keep <- !is.na(bt) & !is.na(mb) & !is.na(tpm) &
          ((bt == "protein_coding" & tpm >= TPM_PC) |
           (bt == "lncRNA"          & tpm >= TPM_LNC))
  genes_base[keep]
}

# ---------------------------------------------------------------------------
# 4. DEG sets — UP only (logFC>0.5), then Cas13-eligible
# ---------------------------------------------------------------------------
integ <- load_dream_results()      # canonical_deg_results.csv (emits bulk_* cols)
integ[, hb := strip_v(gene)]
integ_up <- unique(integ[!is.na(bulk_padj) & !is.na(bulk_logFC) &
                         bulk_padj < PADJ_THR & bulk_logFC > LFC_THR, hb])
integ_set <- cas13_eligible(integ_up)

ps <- load_per_study_de()
ps[, hb := strip_v(gene)]
deg_sets <- lapply(COHORTS, function(ds) {
  up <- unique(ps[dataset == ds & !is.na(padj) & !is.na(logFC) &
                  padj < PADJ_THR & logFC > LFC_THR, hb])
  cas13_eligible(up)
})
names(deg_sets) <- COHORT_SHORT[COHORTS]

# diagnostics (attrition)
cat(sprintf("Integrated: UP=%d -> Cas13-eligible=%d\n", length(integ_up), length(integ_set)))
for (nm in names(deg_sets)) {
  ds <- COHORTS[match(nm, COHORT_SHORT[COHORTS])]
  up <- unique(ps[dataset == ds & padj < PADJ_THR & logFC > LFC_THR & !is.na(padj) & !is.na(logFC), hb])
  cat(sprintf("  %-8s UP=%d -> Cas13-eligible=%d\n", nm, length(up), length(deg_sets[[nm]])))
}

per_study_union <- unique(unlist(deg_sets))
all5_core       <- Reduce(intersect, deg_sets)
gene_cohort_n   <- table(unlist(deg_sets))
two_plus        <- names(gene_cohort_n)[gene_cohort_n >= 2]
n_int <- length(integ_set)
ov_union <- length(intersect(integ_set, per_study_union))
ov_2plus <- length(intersect(integ_set, two_plus))
ov_all5  <- length(intersect(integ_set, all5_core))
cat(sprintf("\nCas13-eligible: integrated=%d | union=%d (∩=%d, int-only=%d) | 2+=%d (∩=%d) | all5=%d (∩=%d, all5-not-int=%d)\n",
            n_int, length(per_study_union), ov_union, n_int - ov_union,
            length(two_plus), ov_2plus, length(all5_core), ov_all5, length(all5_core) - ov_all5))

# ---------------------------------------------------------------------------
# Geometry helpers (same as the methods-validation venn)
# ---------------------------------------------------------------------------
integrated_color <- "#D4834A"; union_color <- "#5480C2"
twoplus_color    <- "#2E9A86"; all5_color  <- "#00695C"
darken <- function(col, f = 0.68) { m <- col2rgb(col)/255; rgb(m[1]*f, m[2]*f, m[3]*f) }
circle_intersect_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- acos(pmin(1, pmax(-1, (d^2 + r1^2 - r2^2)/(2*d*r1))))
  a2 <- acos(pmin(1, pmax(-1, (d^2 + r2^2 - r1^2)/(2*d*r2))))
  r1^2*(a1 - sin(2*a1)/2) + r2^2*(a2 - sin(2*a2)/2)
}
find_center_dist <- function(r1, r2, target) {
  max_ov <- pi * min(r1, r2)^2
  if (target <= 1e-10) return(r1 + r2 + 0.01)
  if (target >= max_ov - 1e-10) return(max(0, abs(r1 - r2)))
  uniroot(function(d) circle_intersect_area(r1, r2, d) - target,
          lower = abs(r1 - r2) + 1e-8, upper = r1 + r2 - 1e-8, tol = 1e-9)$root
}
two_set_venn <- function(nA, nB, nOverlap, colA, colB, titleA, titleB) {
  nA_only <- nA - nOverlap; nB_only <- nB - nOverlap
  r_B <- 1.0; r_A <- sqrt(nA / nB); target_A <- (nOverlap / nB) * pi
  d <- find_center_dist(r_A, r_B, target_A)
  cx_A <- -d/2; cx_B <- d/2
  tx_A <- -(r_A + r_B)/2; tx_ovlp <- (r_A - r_B)/2; tx_B <- (r_A + r_B)/2
  r_max <- max(r_A, r_B)
  circ <- data.frame(x0 = c(cx_A, cx_B), y0 = c(0, 0), r = c(r_A, r_B), grp = c("A", "B"))
  ggplot(circ) +
    geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = grp), color = NA, alpha = 0.32) +
    scale_fill_manual(values = c(A = colA, B = colB), guide = "none") +
    annotate("text", x = tx_A,    y = 0, label = comma(nA_only),  size = 2.8, color = "black", fontface = "bold") +
    annotate("text", x = tx_ovlp, y = 0, label = comma(nOverlap), size = 2.8, color = "black", fontface = "bold") +
    annotate("text", x = tx_B,    y = 0, label = comma(nB_only),  size = 2.8, color = "black", fontface = "bold") +
    annotate("text", x = cx_A, y = r_A + 0.20, label = titleA, size = 2.3, color = "black", fontface = "bold") +
    annotate("text", x = cx_B, y = r_B + 0.20, label = titleB, size = 2.3, color = "black", fontface = "bold") +
    coord_fixed(xlim = c(min(tx_A, cx_A - r_A) - 0.55, max(tx_B, cx_B + r_B) + 0.55),
                ylim = c(-r_max - 0.12, r_max + 0.55), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(-2, 6, -2, 6))
}

p_union <- two_set_venn(n_int, length(per_study_union), ov_union,
  integrated_color, union_color, "Integrated", "Per-study union")
p_2plus <- two_set_venn(n_int, length(two_plus), ov_2plus,
  integrated_color, twoplus_color, "Integrated", "DE in 2+ cohorts")
p_all5  <- two_set_venn(length(all5_core), n_int, ov_all5,
  all5_color, integrated_color, "DE in all 5 cohorts", "Integrated")

combined <- (p_union / p_2plus / p_all5) +
  plot_annotation(
    title = "Integrated vs per-study DEGs",
    theme = theme(plot.title  = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 1)),
                  plot.margin = margin(6, 6, 4, 6)))

# caption detail (was the on-figure subtitle) -> stdout for the manuscript caption
cat("Eligibility: UP (logFC>0.5) + MCD mouse TPM >= 1 (PC) / 0.5 (lncRNA)\n")

out <- file.path(FIGS_CAS13LIB_DIR, "integrated_vs_perstudy_venn.pdf")
save_fig(combined, out, width = fig_half_width + 0.8, height = (fig_half_width + 0.4) * 3 * 0.62)
cat("Saved: ", out, "\n")
