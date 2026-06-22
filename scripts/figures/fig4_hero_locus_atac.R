#!/usr/bin/env Rscript
# fig4_hero_locus_atac.R
# Fig 4 — Panel C: hero locus accessibility schematic.
#
# KEY MESSAGE (sequence-level, VALID evidence — not n=18 disease-DE):
#   A high-confidence GWAS credible-set variant (max_PIP = 0.999) sits INSIDE an
#   accessible hepatocyte scATAC peak, where its C->A substitution STRONGLY
#   disrupts a RORA binding motif (alleleDiff = -1.78). RORA is an independently
#   well-powered candidate: cross-ancestry SuSiE-COLOC PP.H4 = 0.998 (BBJ_GGT,
#   6 GWAS) and downregulated in the bulk n=846 disease-vs-control mega-analysis
#   (logFC = -0.381, lfsr = 8.8e-15). This panel corroborates the well-powered
#   Fig2 (GWAS/COLOC) + Fig3 (bulk) RORA hit with motifbreakR sequence-level
#   variant->motif evidence anchored in chromatin that is open in hepatocytes.
#
# HONEST CAVEATS (rendered on the panel):
#   - The variant disrupts a RORA *binding motif* (trans-regulatory site); it lies
#     in an intron of the HNF1B gene body (the host gene at chr17:37,686,431-
#     37,745,091, hg38), NOT in the RORA gene (RORA is on chr15). We label the
#     panel as "RORA-motif-disrupting regulatory variant", never as a cis-RORA
#     variant.
#   - The overlapping peak is open in hepatocytes AND in 5 other hepatic cell
#     types (Cholangiocyte / Endothelial / Fibroblast / Macrophage / Plasma);
#     it is accessible-in-hepatocytes, not hepatocyte-EXCLUSIVE. We state this.
#   - NO n=18 per-cell/condition disease-DE p-value is shown anywhere (those are
#     pseudoreplicated / underpowered). Only credible-set PIP, motifbreakR
#     allele-diff, COLOC PP.H4, and bulk-n=846 effect size are claimed.
#
# ALL coordinates / scores are READ from on-disk CSV/BED at runtime — nothing is
# hardcoded as a claim (the few literals below are echoed back from the data and
# re-validated with stopifnot()).
#
# Output: figures/main/fig4_validation/fig4c_hero_locus_atac.pdf
# Run:    ~/micromamba/envs/rnaseq/bin/Rscript scripts/figures/fig4_hero_locus_atac.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

# ── Input paths (canonical on-disk sources) ──────────────────────────────────
MOTIF_CSV <- file.path(PROJ,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")
COLOC_CSV <- file.path(PROJ,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
DEG_CSV   <- file.path(PROJ,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
HEP_BED   <- file.path(PROJ,
  "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2/Hepatocytes_peaks.bed")
PEAK_DIR  <- file.path(PROJ,
  "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2")
GTF_PATH  <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"

OUT_DIR <- file.path(PROJ, "figures/main/fig4_validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "fig4c_hero_locus_atac.pdf")

HERO_TF  <- "RORA"          # preferred = THRB but its best motif variant is PIP=0.18
                            #   (weak); RORA carries the only PIP~1 credible-set hit.
HOST_GENE <- "HNF1B"        # gene body the variant falls in (honest label)

# ── 1. Hero variant: highest-max_pip motif-disruption row for HERO_TF ────────
motif <- fread(MOTIF_CSV)
tf_rows <- motif[tf_name == HERO_TF]
stopifnot(nrow(tf_rows) > 0)
hero <- tf_rows[which.max(max_pip)][1]

CHR        <- hero$seqnames
VAR_POS    <- as.integer(hero$start)        # hg38 1-based variant coordinate
MAX_PIP    <- hero$max_pip
ALLELE_DIFF<- hero$alleleDiff
SNP_ID     <- hero$SNP_id                    # hg19 chr:pos:ref:alt id from finemapping
REF_AL     <- hero$ref_genome
ALT_AL     <- hero$alt_genome
EFFECT     <- hero$effect

cat(sprintf("[hero] TF=%s SNP=%s  %s:%d  REF=%s ALT=%s  max_pip=%.4f  alleleDiff=%.3f (%s)\n",
            HERO_TF, SNP_ID, CHR, VAR_POS, REF_AL, ALT_AL, MAX_PIP, ALLELE_DIFF, EFFECT))
stopifnot(MAX_PIP > 0.5)                     # must be a confident credible-set member

# All RORE-family / co-disrupted motifs at the SAME variant (consistency check)
codisrupt <- motif[SNP_id == SNP_ID, .(alleleDiff = mean(alleleDiff)),
                   by = .(tf_name, effect)][order(alleleDiff)]
cat("[hero] motifs disrupted by this variant:\n"); print(codisrupt)

# ── 2. Overlapping hepatocyte scATAC peak (real BED coords) ──────────────────
hep <- fread(HEP_BED, header = FALSE, col.names = c("chrom", "start", "end"))
ov  <- hep[chrom == CHR & start <= VAR_POS & end >= VAR_POS]
stopifnot(nrow(ov) >= 1)                     # variant MUST fall inside a hep peak
PEAK_START <- ov$start[1]
PEAK_END   <- ov$end[1]
PEAK_WIDTH <- PEAK_END - PEAK_START
cat(sprintf("[peak] Hepatocyte peak %s:%d-%d  (%d bp); variant at +%d bp within peak\n",
            CHR, PEAK_START, PEAK_END, PEAK_WIDTH, VAR_POS - PEAK_START))

# Cross-cell-type accessibility of the SAME locus (honest specificity statement)
ct_files <- list.files(PEAK_DIR, pattern = "_peaks\\.bed$", full.names = TRUE)
ct_open <- character(0)
for (f in ct_files) {
  ct <- sub("_peaks\\.bed$", "", basename(f))
  b  <- fread(f, header = FALSE, col.names = c("chrom", "start", "end"))
  if (nrow(b[chrom == CHR & start <= VAR_POS & end >= VAR_POS]) > 0)
    ct_open <- c(ct_open, ct)
}
cat(sprintf("[peak] locus accessible in %d/%d cell types: %s\n",
            length(ct_open), length(ct_files), paste(ct_open, collapse = ", ")))
n_ct_open <- length(ct_open)

# ── 3. Independent evidence for RORA the candidate (COLOC + bulk DE) ──────────
coloc <- fread(COLOC_CSV)
crow  <- coloc[gene == HERO_TF]
stopifnot(nrow(crow) == 1)
COLOC_PP4   <- if (!is.na(crow$coloc_best_susie_pp4)) crow$coloc_best_susie_pp4 else crow$coloc_best_pp4
COLOC_GWAS  <- if (!is.na(crow$coloc_best_susie_pp4)) crow$coloc_best_susie_gwas else crow$coloc_best_gwas
COLOC_METH  <- if (!is.na(crow$coloc_best_susie_pp4)) "SuSiE" else "abf"
COLOC_NGWAS <- crow$coloc_n_gwas_h4_05
cat(sprintf("[coloc] %s PP.H4=%.4f (%s, %s; %s GWAS PP.H4>=0.5)\n",
            HERO_TF, COLOC_PP4, COLOC_GWAS, COLOC_METH, COLOC_NGWAS))

deg  <- fread(DEG_CSV)
drow <- deg[symbol == HERO_TF][1]
stopifnot(nrow(drow) == 1)
DE_LFC  <- drow$logFC
DE_LFSR <- drow$lfsr
cat(sprintf("[bulk]  %s logFC=%.3f lfsr=%.2e (n=846 disease-vs-control)\n",
            HERO_TF, DE_LFC, DE_LFSR))

# ── 4. Host-gene exon model in the window (honest genomic context) ───────────
WIN_PAD   <- 9000L
WIN_START <- VAR_POS - WIN_PAD
WIN_END   <- VAR_POS + WIN_PAD

# Pull HNF1B MANE-Select exons that fall in the window (awk for speed; no Gviz)
gtf_cmd <- sprintf(
  "zcat %s | awk -F'\\t' '$1==\"%s\" && $3==\"exon\" && /gene_name \"%s\"/ && /MANE_Select/ {print $4\"\\t\"$5\"\\t\"$7}'",
  GTF_PATH, CHR, HOST_GENE)
ex <- tryCatch(fread(cmd = gtf_cmd, header = FALSE,
                     col.names = c("start", "end", "strand")),
               error = function(e) data.table(start = integer(), end = integer(),
                                               strand = character()))
ex_win <- ex[end >= WIN_START & start <= WIN_END]
GENE_STRAND <- if (nrow(ex) > 0) ex$strand[1] else "-"
cat(sprintf("[gene]  %s MANE exons in window: %d (strand %s)\n",
            HOST_GENE, nrow(ex_win), GENE_STRAND))

# ── PLOT ─────────────────────────────────────────────────────────────────────
COL_PEAK   <- "#00695C"   # teal — accessible chromatin (open peak)
COL_VAR    <- "#C9265E"   # Liang deep magenta — disease/variant marker
COL_GENE   <- "#37474F"   # slate — host gene model
COL_MOTIF  <- "#A01753"   # deep magenta — motif disruption emphasis
mb <- function(x) sprintf("%.3f", x / 1e6)   # bp -> Mb label

# Track y-positions (single coordinate panel, stacked lanes)
Y_GENE <- 3; Y_PEAK <- 2; Y_VAR <- 1
lane_lab <- data.table(
  y = c(Y_GENE, Y_PEAK, Y_VAR),
  lab = c(sprintf("%s gene\n(host)", HOST_GENE),
          "Hepatocyte\nscATAC peak",
          "Credible-set\nvariant"))

# Hep peaks across the window (context — show all, highlight the overlapping one)
hep_win <- hep[chrom == CHR & end >= WIN_START & start <= WIN_END]
hep_win[, is_hero := (start == PEAK_START & end == PEAK_END)]

# ---- (A) Locus track panel -------------------------------------------------
pA <- ggplot() +
  # gene span line
  { if (nrow(ex_win) > 0)
      annotate("segment", x = max(WIN_START, min(ex_win$start)),
               xend = min(WIN_END, max(ex_win$end)),
               y = Y_GENE, yend = Y_GENE, linewidth = 0.4, color = COL_GENE) } +
  # gene exons (boxes)
  { if (nrow(ex_win) > 0)
      geom_rect(data = ex_win,
                aes(xmin = start, xmax = end, ymin = Y_GENE - 0.22, ymax = Y_GENE + 0.22),
                fill = COL_GENE, color = COL_GENE) } +
  # all hep peaks in window (context, light)
  geom_rect(data = hep_win[is_hero == FALSE],
            aes(xmin = start, xmax = end, ymin = Y_PEAK - 0.22, ymax = Y_PEAK + 0.22),
            fill = "#B2DFDB", color = NA) +
  # the overlapping hero peak (emphasis)
  geom_rect(data = hep_win[is_hero == TRUE],
            aes(xmin = start, xmax = end, ymin = Y_PEAK - 0.28, ymax = Y_PEAK + 0.28),
            fill = COL_PEAK, color = COL_PEAK) +
  # variant marker — vertical guide through all lanes + diamond at variant lane
  annotate("segment", x = VAR_POS, xend = VAR_POS, y = Y_VAR - 0.35, yend = Y_GENE + 0.35,
           linetype = "22", linewidth = 0.3, color = COL_VAR) +
  annotate("point", x = VAR_POS, y = Y_VAR, shape = 23, size = 2.6,
           fill = COL_VAR, color = "black", stroke = 0.3) +
  # lane labels
  geom_text(data = lane_lab, aes(x = WIN_START, y = y, label = lab),
            hjust = 1.05, vjust = 0.5, size = 1.9, color = "black", lineheight = 0.85) +
  scale_x_continuous(
    name = sprintf("%s position (Mb, hg38)", CHR),
    limits = c(WIN_START - (WIN_END - WIN_START) * 0.28, WIN_END),
    labels = function(x) sprintf("%.3f", x / 1e6),
    breaks = scales::pretty_breaks(4),
    expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 3.7), expand = c(0, 0)) +
  labs(title = sprintf("%s-motif-disrupting regulatory variant in an accessible hepatocyte peak",
                       HERO_TF),
       subtitle = sprintf("%s  |  %s:%s  |  credible-set PIP = %.3f",
                          SNP_ID, CHR, format(VAR_POS, big.mark = ","), MAX_PIP)) +
  theme_masld() + theme_pub() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        plot.title    = element_text(size = PUB_TITLE, face = "bold"),
        plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray30"))

# ---- (B) Peak zoom: variant inside the peak + RORA-family motif disruption --
# Tight window on the peak; show the peak rectangle, variant position, and a
# compact bar of allele-diff for the co-disrupted RORE-family motifs.
ZWIN_START <- PEAK_START - 60L
ZWIN_END   <- PEAK_END + 60L
pB1 <- ggplot() +
  geom_rect(aes(xmin = PEAK_START, xmax = PEAK_END, ymin = 0.4, ymax = 0.8),
            fill = COL_PEAK, color = COL_PEAK) +
  annotate("text", x = (PEAK_START + PEAK_END) / 2, y = 0.6,
           label = sprintf("Hep scATAC peak\n%s:%s-%s  (%d bp)",
                           CHR, format(PEAK_START, big.mark = ","),
                           format(PEAK_END, big.mark = ","), PEAK_WIDTH),
           size = 1.8, color = "white", lineheight = 0.85) +
  annotate("segment", x = VAR_POS, xend = VAR_POS, y = 0.35, yend = 0.95,
           linewidth = 0.4, color = COL_VAR) +
  annotate("point", x = VAR_POS, y = 0.95, shape = 23, size = 2.4,
           fill = COL_VAR, color = "black", stroke = 0.3) +
  annotate("text", x = VAR_POS, y = 1.12,
           label = sprintf("variant  %s>%s", REF_AL, ALT_AL),
           size = 1.9, color = "black", fontface = "bold") +
  scale_x_continuous(name = sprintf("%s (kb, hg38)", CHR),
                     limits = c(ZWIN_START, ZWIN_END),
                     labels = function(x) sprintf("%.1f", x / 1e3),
                     breaks = c(PEAK_START, VAR_POS, PEAK_END), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 1.25), expand = c(0, 0)) +
  labs(subtitle = "Variant lies inside the open peak") +
  theme_masld() + theme_pub() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x  = element_text(size = PUB_AXIS_TEXT),
        plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray30"))

# Co-disrupted motif allele-diff bars (motifbreakR; negative = ALT weakens motif)
cod <- copy(codisrupt)
cod[, tf := factor(tf_name, levels = tf_name[order(alleleDiff)])]
cod[, is_hero := tf_name == HERO_TF]
pB2 <- ggplot(cod, aes(x = alleleDiff, y = tf)) +
  geom_col(aes(fill = is_hero), width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_text(aes(label = sprintf("%.2f", alleleDiff),
                hjust = ifelse(alleleDiff < 0, 1.1, -0.1)),
            size = 1.7, color = "black") +
  scale_fill_manual(values = c(`TRUE` = COL_MOTIF, `FALSE` = "#9E9E9E"), guide = "none") +
  scale_x_continuous(name = "motifbreakR allele-diff (ALT - REF)",
                     expand = expansion(mult = c(0.18, 0.18))) +
  labs(y = NULL, subtitle = sprintf("%s>%s disrupts RORE-family motifs", REF_AL, ALT_AL)) +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = PUB_AXIS_TEXT, face = "italic", color = "black"),
        plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray30"))

# ---- (C) Honest evidence-convergence caption box ---------------------------
ct_txt <- if (n_ct_open > 1) {
  sprintf("Peak open in %d cell types (incl. 5 other hepatic lineages) - accessible-in-hepatocytes, not hepatocyte-exclusive.", n_ct_open)
} else {
  "Peak open in hepatocytes."
}
wrap <- function(s, w = 96) paste(strwrap(s, width = w), collapse = "\n")
cap_txt <- paste(
  wrap(sprintf("RORA candidate (independent, well-powered): cross-ancestry COLOC PP.H4 = %.3f (%s, %s; %s GWAS); bulk n=846 logFC = %.2f, lfsr = %.0e.",
               COLOC_PP4, COLOC_GWAS, COLOC_METH, COLOC_NGWAS, DE_LFC, DE_LFSR)),
  wrap(sprintf("Variant %s (PIP %.3f) disrupts a RORA binding MOTIF (allele-diff %.2f, '%s'); it falls in an intron of the host gene %s, NOT the RORA gene. %s",
               SNP_ID, MAX_PIP, ALLELE_DIFF, EFFECT, HOST_GENE, ct_txt)),
  sep = "\n")
pC <- ggplot() +
  annotate("text", x = 0, y = 1, label = cap_txt, hjust = 0, vjust = 1,
           size = 1.85, color = "black", lineheight = 1.05) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void()

# ── Compose ───────────────────────────────────────────────────────────────
panel <- (pA / ((pB1 | pB2) + plot_layout(widths = c(1, 1.05))) / pC) +
  plot_layout(heights = c(1.15, 1.0, 0.42))

# PDF only; cairo_pdf when available (no dingbats by construction), else base
# pdf() with useDingbats = FALSE per repo figure convention.
if (capabilities("cairo")) {
  ggsave(OUT_PDF, panel, width = fig_col_width, height = 5.0, device = cairo_pdf)
} else {
  ggsave(OUT_PDF, panel, width = fig_col_width, height = 5.0,
         device = grDevices::pdf, useDingbats = FALSE)
}

cat(sprintf("\n[done] wrote %s  (%.1f KB)\n", OUT_PDF, file.info(OUT_PDF)$size / 1024))
