##############################################################################
# Fig 3 Panel K: GWAS-ATAC regulatory chain exemplar (Liang refinement)
#
# Layout (2 x 2 grid, ~6.5 x 7.5 in):
#
#   a  Locus zoom (GWAS PIP) + gene-model strip below      (top, full-width)
#   b  Motif disruption — HNF4A + THRB PWM logos           (middle-left)
#   c  Forest of effect sizes (motif + caQTL + eQTL)       (middle-right)
#   d  Compact schematic flow (variant -> drug)            (bottom, full-width)
#
# Style: Liang manuscript palette (#C9265E magenta, #1565C0 blue,
# #9E9E9E gray); bold lowercase panel letters top-left; no in-plot titles;
# hollow circles for points; dashed gray zero/reference lines; minimal
# whitespace. PWM logos via ggseqlogo (variant column highlighted in magenta).
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(grid)
  library(ggrepel)
  library(MotifDb)
  library(ggseqlogo)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
theme_set(theme_masld())

OUT_DIR <- file.path(FIG_SUPP, "figS_cyfip2_locus")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# Individual per-panel PDFs (no composite).
LOCUS_ZOOM_DIR <- file.path(FIG3_DIR, "locus_zoom")  # FIG3_DIR == figures/main/fig2_genetics
dir.create(LOCUS_ZOOM_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_A <- file.path(LOCUS_ZOOM_DIR, "GGT_CYFIP2.pdf")  # 2026-07-01: moved out of OUT_DIR into the consolidated Fig2 locus_zoom/ dir
OUT_B <- file.path(OUT_DIR, "cyfip2_b_motif_logos.pdf")
OUT_C <- file.path(OUT_DIR, "cyfip2_c_effect_forest.pdf")
OUT_D <- file.path(OUT_DIR, "cyfip2_d_schematic.pdf")

# ---------------------------------------------------------------------------
# Exemplar locus constants
# ---------------------------------------------------------------------------
EX <- list(
  variant_id    = "5:156753817:G:T",
  ref           = "G",
  alt           = "T",
  hg19_pos      = 156753817L,
  hg38_chr      = "chr5",
  hg38_pos      = 157326809L,
  window        = 100000L,
  locus_name    = "chr5q33 / CYFIP2",
  gene          = "CYFIP2",
  primary_tf    = "THRB",
  parallel_tf   = "RORA",
  gwas          = "UKBB_GGT",
  rec_pip       = 0.0003,     # low individual PIP; variant is in CARMA CS
  coloc_pp4     = 0.958,
  primary_eff   = -1.234,     # ALT LOSES THRB motif
  parallel_eff  = -1.698,     # ALT LOSES RORA motif
  eqtl_beta     = -0.329,
  eqtl_se       = 0.041,      # approx from p=6.6e-15
  eqtl_pval     = 6.6e-15,
  caqtl_beta    = -0.338,
  caqtl_se      = 0.097,
  caqtl_pval    = 7.3e-4,
  drug          = "resmetirom"
)

PURPLE <- "#C9265E"
COOL   <- "#1565C0"
GRAY   <- "#9E9E9E"
INK    <- "#222222"

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
GWAS_ATAC_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
FM_FILE       <- file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv")
PEAK_DIR      <- file.path(BASE,
                          "Analysis/ATAC/Human_Multiome/results/label_transfer",
                          "cell_type_peak_sets_v2")
GTF_FILE      <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"

cat("Loading data...\n")

fm_all <- fread(FM_FILE)
fm <- fm_all[
  chromosome == 10 &
  position >= EX$hg19_pos - EX$window &
  position <= EX$hg19_pos + EX$window
]
cat("Fine-mapped variants in +/-500kb:", nrow(fm), "\n")

# Cell-type peak overlap (hg38)
ct_files <- list.files(PEAK_DIR, pattern = "_peaks\\.bed$", full.names = TRUE)
ct_overlap <- rbindlist(lapply(ct_files, function(f) {
  bed <- tryCatch(fread(f, header = FALSE, sep = "\t", select = 1:3,
                        col.names = c("chr", "start", "end")),
                  error = function(e) NULL)
  if (is.null(bed) || nrow(bed) == 0) {
    return(data.table(ct = sub("_peaks\\.bed$", "", basename(f)),
                      overlap = FALSE))
  }
  has <- any(bed$chr == EX$hg38_chr &
             bed$start <= EX$hg38_pos &
             bed$end   >= EX$hg38_pos)
  data.table(ct = sub("_peaks\\.bed$", "", basename(f)), overlap = has)
}))
cat("CTs with peak overlap:", sum(ct_overlap$overlap), "/", nrow(ct_overlap), "\n")

# ---------------------------------------------------------------------------
# Panel a: GWAS PIP locus zoom + gene-model strip
# ---------------------------------------------------------------------------
fm_plot <- copy(fm)
fm_plot[, is_lead := position == EX$hg19_pos]
# Use the most informative PIP across methods (rec/max/susie/carma) so that
# weak-rec-PIP / strong-CARMA-CS variants still show structure on the y-axis.
fm_plot[, plot_pip := pmax(recommended_pip, max_pip, susie_pip, carma_pip,
                           na.rm = TRUE)]
fm_plot[plot_pip < 0, plot_pip := 0]
# Adaptive y-axis: cap at observed max + 20% headroom (no fixed 0-1).
y_max_locus <- max(0.05, max(fm_plot$plot_pip, na.rm = TRUE) * 1.25)
# Mark CARMA credible-set members for context
fm_plot[, in_cs := (!is.na(carma_cs) & carma_cs > 0) |
                   (!is.na(susie_cs) & susie_cs > 0)]

xlims_mb <- c((EX$hg19_pos - EX$window) / 1e6,
              (EX$hg19_pos + EX$window) / 1e6)

p_locus <- ggplot(fm_plot, aes(x = position / 1e6, y = plot_pip)) +
  geom_point(data = fm_plot[is_lead == FALSE & in_cs == FALSE],
             shape = 1, color = "grey55", stroke = 0.3, size = 1.0) +
  geom_point(data = fm_plot[is_lead == FALSE & in_cs == TRUE],
             shape = 1, color = INK, stroke = 0.45, size = 1.4) +
  geom_segment(data = fm_plot[is_lead == TRUE],
               aes(x = position / 1e6, xend = position / 1e6,
                   y = 0, yend = plot_pip),
               color = PURPLE, linewidth = 0.5) +
  geom_point(data = fm_plot[is_lead == TRUE],
             shape = 21, fill = PURPLE, color = "black",
             stroke = 0.5, size = 3.4) +
  geom_text_repel(
    data = fm_plot[is_lead == TRUE],
    aes(label = sprintf("%s:%s  %s>%s\nCARMA CS  PP4 %.2f",
                        EX$hg38_chr,
                        format(position, big.mark = ","),
                        EX$ref, EX$alt, EX$coloc_pp4)),
    nudge_y = 0.20, nudge_x = -0.05, size = 2.1,
    segment.size = 0.3, segment.color = "grey45",
    family = "Helvetica", lineheight = 0.95
  ) +
  scale_y_continuous(limits = c(0, y_max_locus),
                     breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(limits = xlims_mb,
                     labels = function(x) sprintf("%.2f", x),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "PIP") +
  theme_masld() + theme_pub() +
  theme(
    axis.title.y  = element_text(size = 6, face = "plain", color = "black",
                                  angle = 90),
    axis.text     = element_text(size = 6, color = "black"),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid    = element_blank(),
    plot.margin   = margin(2, 6, 0, 6)
  )

# Gene-model strip below
# Liftover-aware: combined_finemapping is hg19; variant hg19=94839724 vs
# hg38=93079967, delta ~ -1,759,757. Apply same delta to gene coords from
# hg38 GTF so genes align with the hg19 x-axis above.
HG38_TO_HG19 <- EX$hg19_pos - EX$hg38_pos  # +1,759,757

gtf <- tryCatch({
  cmd <- sprintf(
    "zcat %s | awk '$1==\"%s\" && $3==\"gene\" && $4>%d && $4<%d'",
    GTF_FILE, EX$hg38_chr,
    EX$hg38_pos - EX$window - 50000, EX$hg38_pos + EX$window + 50000)
  fread(cmd = cmd, header = FALSE, sep = "\t")
}, error = function(e) NULL)

gene_strip <- data.table()
if (!is.null(gtf) && nrow(gtf) > 0) {
  attrs <- gtf$V9
  gene_strip <- data.table(
    start_hg19 = gtf$V4 + HG38_TO_HG19,
    end_hg19   = gtf$V5 + HG38_TO_HG19,
    strand     = gtf$V7,
    gene_type  = sub('.*gene_type "([^"]+)".*', '\\1', attrs),
    gene_name  = sub('.*gene_name "([^"]+)".*', '\\1', attrs)
  )
  gene_strip <- gene_strip[gene_type %in% c("protein_coding")]
  gene_strip <- gene_strip[
    end_hg19 >= EX$hg19_pos - EX$window &
    start_hg19 <= EX$hg19_pos + EX$window
  ]
  gene_strip[, mid_hg19 := (start_hg19 + end_hg19) / 2]
  # Label the locus gene + 5 nearest flanking protein-coding genes; drop
  # readthroughs / numeric Ensembl-only IDs
  gene_strip <- gene_strip[!grepl("^ENSG", gene_name)]
  gene_strip[, dist_to_var := abs(mid_hg19 - EX$hg19_pos)]
  gene_strip[, label := ifelse(rank(dist_to_var) <= 6 |
                                gene_name == EX$gene, gene_name, NA_character_)]
  gene_strip[, is_locus := gene_name == EX$gene]
}

# Build gene strip plot (always returns a small ggplot to keep layout stable)
p_genes <- ggplot() +
  scale_x_continuous(limits = xlims_mb, expand = c(0, 0),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = NULL, expand = c(0, 0)) +
  labs(x = sprintf("%s position (Mb, hg19)", EX$hg38_chr), y = NULL) +
  theme_masld() + theme_pub() +
  theme(
    axis.title.x = element_text(size = 6, face = "plain", color = "black"),
    axis.text.x  = element_text(size = 6, color = "black"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(0, 6, 2, 6)
  )

if (nrow(gene_strip) > 0) {
  p_genes <- p_genes +
    # variant guide line spanning gene strip
    geom_vline(xintercept = EX$hg19_pos / 1e6, color = PURPLE,
               linewidth = 0.4, linetype = "dotted") +
    geom_rect(
      data = gene_strip,
      aes(xmin = start_hg19 / 1e6, xmax = end_hg19 / 1e6,
          ymin = -0.18, ymax = 0.18,
          fill = is_locus, color = is_locus),
      linewidth = 0.25
    ) +
    geom_text_repel(
      data = gene_strip[!is.na(label)],
      aes(x = mid_hg19 / 1e6, y = 0.25, label = label,
          fontface = "italic"),
      size = 2.1, family = "Helvetica",
      color = ifelse(gene_strip[!is.na(label)]$is_locus, PURPLE, "grey20"),
      box.padding = 0.18, point.padding = 0.05,
      segment.size = 0.2, segment.color = "grey55",
      max.overlaps = 30,
      direction = "x", nudge_y = 0.20, ylim = c(0.22, 0.58)
    ) +
    scale_fill_manual(values = c("TRUE" = PURPLE, "FALSE" = "#9E9E9E"),
                      guide = "none") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey45"),
                       guide = "none")
}

# Stack locus + gene strip (height ratio 4 : 1.5). Wrap so patchwork tagging
# treats this as a single panel-a element rather than two children.
panel_a <- wrap_elements(
  full = (p_locus / p_genes) + plot_layout(heights = c(4, 1.5))
)

# ---------------------------------------------------------------------------
# Panel b: PWM logos for HNF4A and THRB via ggseqlogo
# ---------------------------------------------------------------------------
get_motif <- function(tf) {
  q <- query(MotifDb, c(tf, "Hsapiens", "jaspar"))
  if (length(q) == 0) return(NULL)
  prefer <- c("jaspar2022", "jaspar2018", "jaspar2016", "JASPAR_2014")
  for (p in prefer) {
    hit <- names(q)[grepl(p, names(q), ignore.case = TRUE)]
    if (length(hit) > 0) return(q[[hit[1]]])
  }
  q[[1]]
}

# motifbreakR reports the SNP position (no motif-match offset stored in our
# summary CSVs), so we visually mark the column where the variant sits by
# centering the motif on the variant: column ceiling(L/2). This is consistent
# with motifbreakR's "best window" reporting heuristic for both REF and ALT
# and is unambiguous in the figure caption.
build_ppm <- function(tf) {
  ppm <- get_motif(tf)
  if (is.null(ppm)) stop("No motif for ", tf)
  ppm <- pmax(ppm, 1e-6)
  ppm <- sweep(ppm, 2, colSums(ppm), "/")
  rownames(ppm) <- c("A", "C", "G", "T")
  ppm
}

ppm_h <- build_ppm(EX$primary_tf)
ppm_t <- build_ppm(EX$parallel_tf)

# Variant column = center of motif (motifbreakR centers strong matches).
# Use the highest-IC column nearest the center as the marked SNP column.
pick_var_col <- function(ppm) {
  ic <- 2 - (-colSums(ppm * log2(ppm)))
  ic[ic < 0] <- 0
  L <- ncol(ppm)
  center <- (L + 1) / 2
  # IC-weighted distance to center; pick the IC-rich column closest to mid
  scores <- ic - 0.4 * abs(seq_len(L) - center)
  which.max(scores)
}

var_col_h <- pick_var_col(ppm_h)
var_col_t <- pick_var_col(ppm_t)

logo_data <- setNames(list(ppm_h, ppm_t),
                      c(EX$primary_tf, EX$parallel_tf))
ann_data <- data.frame(
  seq_group  = factor(c(EX$primary_tf, EX$parallel_tf),
                      levels = c(EX$primary_tf, EX$parallel_tf)),
  variant_x  = c(var_col_h, var_col_t),
  effect_lab = c(sprintf("alleleDiff %+.2f", EX$primary_eff),
                 sprintf("alleleDiff %+.2f", EX$parallel_eff)),
  motif_len  = c(ncol(ppm_h), ncol(ppm_t))
)

# Compute facet-specific ymax (IC ceiling = 2) for label placement
p_logo <- ggplot() +
  geom_rect(data = ann_data,
            aes(xmin = variant_x - 0.5, xmax = variant_x + 0.5,
                ymin = 0, ymax = 2.05),
            fill = PURPLE, alpha = 0.18, inherit.aes = FALSE) +
  geom_logo(logo_data, method = "bits", seq_type = "dna") +
  geom_text(data = ann_data,
            aes(x = 0.6, y = 2.30, label = effect_lab),
            hjust = 0, vjust = 1, inherit.aes = FALSE,
            size = GEOM_TEXT_6PT, fontface = "plain", color = PURPLE,
            family = "Helvetica") +
  geom_text(data = ann_data,
            aes(x = variant_x, y = -0.35,
                label = sprintf("%s>%s", EX$ref, EX$alt)),
            inherit.aes = FALSE, size = GEOM_TEXT_6PT, fontface = "plain",
            color = PURPLE, family = "Helvetica") +
  facet_wrap(~ seq_group, ncol = 1, scales = "free_x",
             strip.position = "left") +
  scale_y_continuous(limits = c(-0.5, 2.45), breaks = c(0, 1, 2),
                     name = "bits", expand = c(0, 0)) +
  scale_x_continuous(breaks = NULL, name = NULL) +
  theme_masld() + theme_pub() +
  theme(
    axis.title.y      = element_text(size = 6, face = "plain", color = "black",
                                      angle = 90),
    axis.text.y       = element_text(size = 6, color = "black"),
    axis.line.x       = element_blank(),
    strip.text.y.left = element_text(size = 6, face = "italic",
                                      color = "black", angle = 0,
                                      hjust = 1),
    strip.background  = element_blank(),
    strip.placement   = "outside",
    panel.spacing.y   = unit(2.5, "mm"),
    panel.grid        = element_blank(),
    plot.margin       = margin(2, 4, 2, 2)
  )

# ---------------------------------------------------------------------------
# Panel c: Forest plot of 4 measurements
# ---------------------------------------------------------------------------
# Standardize motif alleleDiffs into a comparable "ALT effect" scale (motif
# alleleDiff is in PWM-score units; we report it as-is on a secondary x range
# using a single panel with a vertical separator).
#
# For visual concordance the message is "all four > 0". We split the forest
# into two strips (motif gain | QTL effect) so units are not conflated.

forest_motif <- data.table(
  measurement = c(sprintf("%s motif", EX$primary_tf),
                  sprintf("%s motif", EX$parallel_tf)),
  beta = c(EX$primary_eff, EX$parallel_eff),
  se   = c(NA_real_, NA_real_),
  units = "motif alleleDiff",
  pval  = c(NA_real_, NA_real_)
)
forest_qtl <- data.table(
  measurement = c(sprintf("caQTL  (hepatocyte, %s peak)", EX$gene),
                  sprintf("eQTL  (%s, liver Broadaway)", EX$gene)),
  beta = c(EX$caqtl_beta, EX$eqtl_beta),
  se   = c(EX$caqtl_se, EX$eqtl_se),
  units = "QTL beta (ALT)",
  pval  = c(EX$caqtl_pval, EX$eqtl_pval)
)
forest <- rbind(forest_motif, forest_qtl)
forest[, lo := ifelse(is.na(se), beta, beta - 1.96 * se)]
forest[, hi := ifelse(is.na(se), beta, beta + 1.96 * se)]
forest[, label := ifelse(is.na(pval),
                          sprintf("%+.2f", beta),
                          sprintf("%+.2f  p=%.0e", beta, pval))]
forest[, units := factor(units, levels = c("motif alleleDiff",
                                            "QTL beta (ALT)"))]
forest[, measurement := factor(measurement, levels = rev(measurement))]

p_forest <- ggplot(forest, aes(x = beta, y = measurement)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35,
             linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 color = INK, linewidth = 0.6, na.rm = TRUE) +
  geom_point(shape = 21, fill = "white", color = PURPLE,
             size = 3.6, stroke = 0.9) +
  geom_text(aes(label = label),
            hjust = -0.10, vjust = -0.8, size = GEOM_TEXT_6PT,
            color = INK, family = "Helvetica", fontface = "plain") +
  facet_wrap(~ units, ncol = 1, scales = "free",
             strip.position = "top") +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.06))) +
  labs(x = "ALT-allele effect", y = NULL) +
  theme_masld() + theme_pub() +
  theme(
    axis.title.x      = element_text(size = 6, face = "plain", color = "black"),
    axis.text         = element_text(size = 6, color = "black"),
    strip.text        = element_text(size = 6, face = "plain", color = "black"),
    strip.background  = element_blank(),
    panel.spacing.y   = unit(3, "mm"),
    panel.grid        = element_blank(),
    plot.margin       = margin(2, 6, 2, 2)
  )

# ---------------------------------------------------------------------------
# Panel d: Schematic flow — boxes + arrows
# ---------------------------------------------------------------------------
nodes <- data.table(
  i = 1:6,
  x = c(0.06, 0.22, 0.40, 0.58, 0.76, 0.94),
  y = 0.55,
  primary = c("chr5:156753817\nG>T  (CARMA CS)",
              "THRB + RORA\nmotif loss",
              "chromatin closing\ncaQTL  beta -0.34",
              "CYFIP2 down\neQTL  beta -0.33",
              "GGT GWAS\n(UKBB)  PP4 0.96",
              "resmetirom\nTHRB agonism\n(compensates)"),
  layer   = c("variant", "motif loss", "chromatin", "expression",
              "trait", "drug rescue"),
  hue     = c(GRAY, COOL, COOL, PURPLE, PURPLE, "#1B7F3A")
)

p_flow <- ggplot() +
  # arrows between consecutive nodes
  geom_segment(
    data = data.table(x = nodes$x[-6] + 0.055,
                      xend = nodes$x[-1] - 0.055,
                      y = 0.55, yend = 0.55),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
    color = "grey35", linewidth = 0.55
  ) +
  geom_point(data = nodes,
             aes(x = x, y = y, fill = hue),
             shape = 21, color = "black", stroke = 0.6, size = 7) +
  geom_text(data = nodes,
            aes(x = x, y = y - 0.28, label = primary),
            size = GEOM_TEXT_6PT, family = "Helvetica",
            fontface = "plain", color = INK,
            lineheight = 0.95) +
  geom_text(data = nodes,
            aes(x = x, y = y + 0.20, label = layer),
            size = 2.1, family = "Helvetica",
            fontface = "plain", color = "grey35") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(2, 6, 2, 6))

# ---------------------------------------------------------------------------
# Save each panel as its own PDF at a natural aspect ratio (no composite).
# Sizes chosen so each panel reads cleanly on its own without composite-grid
# squeezing.
# ---------------------------------------------------------------------------
message(sprintf("[caption] Panel a: %s locus zoom (GGT GWAS)", EX$gene))
cat("Saving panel a:", OUT_A, "\n")
ggsave(OUT_A, panel_a,        width = 4.6, height = 2.8, device = cairo_pdf)

cat("Saving panel b:", OUT_B, "\n")
ggsave(OUT_B, p_logo,         width = 5.0, height = 3.0, device = cairo_pdf)

cat("Saving panel c:", OUT_C, "\n")
ggsave(OUT_C, p_forest,       width = 5.0, height = 3.0, device = cairo_pdf)

cat("Saving panel d:", OUT_D, "\n")
ggsave(OUT_D, p_flow,         width = 8.0, height = 1.8, device = cairo_pdf)

# ---------------------------------------------------------------------------
# Source CSV for caption transparency
# ---------------------------------------------------------------------------
chain_csv <- file.path(OUT_DIR, "cyfip2_chain_values.csv")
chain_dt <- data.table(
  layer = c("GWAS variant (hg19)", "GWAS variant (hg38)", "GWAS",
            "recommended PIP",
            sprintf("COLOC PP.H4 (%s)", EX$gene), "CT peak overlap",
            sprintf("%s motif alleleDiff", EX$primary_tf),
            sprintf("%s motif alleleDiff", EX$parallel_tf),
            sprintf("Broadaway eQTL beta (%s)", EX$gene), "Broadaway eQTL p",
            "Currin caQTL beta", "Currin caQTL p",
            "Linked gene", "Primary TF (drug)", "Parallel TF", "Drug"),
  value = c(EX$variant_id,
            sprintf("%s:%d %s>%s", EX$hg38_chr, EX$hg38_pos, EX$ref, EX$alt),
            EX$gwas,
            sprintf("%.4f", EX$rec_pip),
            sprintf("%.3f", EX$coloc_pp4),
            paste(ct_overlap[overlap == TRUE]$ct, collapse = ", "),
            sprintf("%+.3f", EX$primary_eff),
            sprintf("%+.3f", EX$parallel_eff),
            sprintf("%+.4f", EX$eqtl_beta),
            sprintf("%.2e", EX$eqtl_pval),
            sprintf("%+.4f", EX$caqtl_beta),
            sprintf("%.2e", EX$caqtl_pval),
            EX$gene, EX$primary_tf, EX$parallel_tf, EX$drug)
)
fwrite(chain_dt, chain_csv)
cat("Wrote chain CSV:", chain_csv, "\n")
cat("Done.\n")
