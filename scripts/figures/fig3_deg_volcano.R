#!/usr/bin/env Rscript
# fig3_deg_volcano.R — Bulk DEG volcanoes for Fig 3 (RNA-seq).
#
# Left: canonical limma-voom quality-weighted C2 disease-vs-control call.
# Right: sensitivity call on the exact definitive samples plotted in Fig S3U.
#
# Tier 1 DEG definition (canonical, 2026-06-29): TREAT FDR < 0.05 at lfc = 0.25
# (= 1,918 pooled DEGs; 1,419 up / 499 down). treat() tests H0:|log2FC| <= 0.25,
# so the effect-size floor is folded INTO the test — there is NO separate
# |logFC|/|shrunk| filter. The volcano X = raw log2FC (treat tests the raw effect),
# Y = -log10(treat FDR); up/down by sign(logFC); dashed guides mark the ±0.25 lfc.
# A raw-scale companion (padj<0.05 & |logFC|>0.5) is emitted as a sensitivity
# supplement: *_raw_sensitivity.pdf.
#
# Output: figures/main/fig3_RNAseq/panels/figs3a_deg_volcano.pdf
# (moved from fig1_atlas_overview 2026-06-12 — this is an RNA-seq DEG panel.)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")   # FIG2_DIR = figures/main/fig3_RNAseq
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# Canonical TREAT gate
TREAT_FDR_CUT <- 0.05
TREAT_LFC     <- 0.25   # treat() effect-size offset (H0:|log2FC|<=lfc), folded into the test
# Raw sensitivity gate
PADJ_CUT  <- 0.05
LFC_CUT   <- 0.5
N_TOP_DIR <- 6

# --- Analytical TREAT (identical to limma::treat; reconstructed from the
# moderated-t stat). Mirrors Cas13 rebuild_cas13_library.R::add_treat_fdr.
# Used for tables (e.g. the definitive-sample arm) that ship logFC+SE+t+P.Value
# but not a precomputed treat_fdr column. ----------------------------------
.infer_df <- function(dt) {
  pr <- dt[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
  idx <- unique(round(seq(1, nrow(pr), length.out = min(nrow(pr), 12))))
  median(vapply(idx, function(i)
    uniroot(function(df) 2 * pt(-abs(pr$t[i]), df = df) - pr$P.Value[i], c(0.1, 1e6))$root,
    numeric(1)))
}
treat_fdr_vec <- function(d, lfc = TREAT_LFC, se_col = "SE", df_col = NULL) {
  se <- if (!is.null(se_col) && se_col %in% names(d)) d[[se_col]] else abs(d$logFC / d$t)
  se[!is.finite(se) | se <= 0] <- NA_real_
  dfu <- if (!is.null(df_col) && df_col %in% names(d)) d[[df_col]] else .infer_df(d)
  p <- pt((abs(d$logFC) - lfc) / se, df = dfu, lower.tail = FALSE) +
       pt((abs(d$logFC) + lfc) / se, df = dfu, lower.tail = FALSE)
  p.adjust(p, method = "BH")
}

# DEG_VOLCANO_VECTOR=TRUE renders point clouds as native vector geometry and
# writes to deg_volcano_vector.pdf (heavier file, but fully editable).
VECTOR_MODE <- identical(toupper(Sys.getenv("DEG_VOLCANO_VECTOR", "FALSE")), "TRUE")

CANONICAL_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
DEFINITIVE_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sensitivity/figs3u_plotted_points_de.csv")

if (!file.exists(DEFINITIVE_PATH)) {
  stop("Missing definitive-sample DE table: ", DEFINITIVE_PATH,
       "\nRun RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05j_figs3u_plotted_points_de.R first.")
}

# Curated MASLD-relevant anchor genes — labelled if Tier 1.
CURATED <- c(
  "THRB",     # resmetirom target
  "HNF4A",    # master hepatocyte TF
  "NR1H4",    # FXR / obeticholic acid
  "PPARA",    # elafibranor / fenofibrate
  "RORA",     # cross-ancestry COLOC
  "AKR1B10",  # largest bulk effect, MASLD biomarker
  "TREM2",    # LAM macrophage marker
  "LPL",      # lipoprotein lipase, lipid metabolism
  "GSTM1",    # antioxidant
  "POSTN",    # fibrosis ECM
  "GSN",      # advanced fibrosis
  "FABP4",    # macrophage / lipid
  "SPP1",     # scar-associated macrophage
  "COL1A1",   # fibrosis
  "HKDC1"     # F1->F2 progression driver
)

volc_colors <- c(
  "Up"   = masld_colors$up,
  "Down" = masld_colors$down,
  "n.s." = "#D8D8D8"
)

# mode = "treat" (canonical): x = raw logFC, y = -log10(treat FDR); significance
#   = treat_fdr < 0.05 (the effect floor is folded into the test); up/down by
#   sign(logFC); dashed guides mark the ±0.25 lfc. If the table lacks a
#   precomputed treat_fdr it is reconstructed analytically (definitive arm).
# mode = "raw" (sensitivity):  x = logFC, y = -log10(padj), gate |x|>0.5 & padj<0.05
prep_volcano <- function(dt, mode = "treat") {
  out <- copy(dt)

  if (mode == "treat") {
    if (!"treat_fdr" %in% names(out)) {
      req <- c("symbol", "logFC", "t", "P.Value")
      missing <- setdiff(req, names(out))
      if (length(missing) > 0) stop("Missing cols for TREAT: ", paste(missing, collapse = ", "))
      out[, treat_fdr := treat_fdr_vec(out, lfc = TREAT_LFC,
            se_col = if ("SE" %in% names(out)) "SE" else NULL,
            df_col = if ("df.total" %in% names(out)) "df.total" else NULL)]
    }
    out <- out[!is.na(treat_fdr) & !is.na(logFC),
               .(symbol, xx = logFC, pp = treat_fdr)]
    xcut <- TREAT_LFC; pcut <- TREAT_FDR_CUT
  } else {
    req <- c("symbol", "logFC", "padj")
    missing <- setdiff(req, names(out))
    if (length(missing) > 0) stop("Missing expected columns: ", paste(missing, collapse = ", "))
    out <- out[!is.na(padj) & !is.na(logFC),
               .(symbol, xx = logFC, pp = padj)]
    xcut <- LFC_CUT; pcut <- PADJ_CUT
  }
  min_pos <- min(out[pp > 0, pp], na.rm = TRUE)
  out[pp <= 0, pp := min_pos]
  out[, neglog10p := -log10(pp)]
  if (mode == "treat") {
    # significance by treat FDR alone; direction by sign(logFC) (no |x| filter)
    out[, status := fcase(
      pp < pcut & xx > 0, "Up",
      pp < pcut & xx < 0, "Down",
      default = "n.s."
    )]
  } else {
    out[, status := fcase(
      pp < pcut &  xx >  xcut, "Up",
      pp < pcut &  xx < -xcut, "Down",
      default = "n.s."
    )]
  }
  out[, status := factor(status, levels = c("n.s.", "Down", "Up"))]
  setorder(out, status)
  attr(out, "xcut") <- xcut; attr(out, "pcut") <- pcut; attr(out, "mode") <- mode
  out
}

label_genes <- function(volc) {
  sig_up   <- volc[status == "Up"   & symbol != "" & !grepl("^ENSG", symbol)]
  sig_down <- volc[status == "Down" & symbol != "" & !grepl("^ENSG", symbol)]
  setorder(sig_up, pp)
  setorder(sig_down, pp)

  top_up   <- head(sig_up,   N_TOP_DIR)
  top_down <- head(sig_down, N_TOP_DIR)
  curated_lbl <- volc[symbol %in% CURATED & status != "n.s." & !grepl("^ENSG", symbol)]

  unique(rbind(top_up, top_down, curated_lbl), by = "symbol")
}

make_volcano <- function(volc, title) {
  label_df <- label_genes(volc)
  xcut <- attr(volc, "xcut"); pcut <- attr(volc, "pcut"); mode <- attr(volc, "mode")
  n_up   <- sum(volc$status == "Up")
  n_down <- sum(volc$status == "Down")
  n_ns   <- sum(volc$status == "n.s.")
  gate_lab <- if (mode == "treat") sprintf("adjusted FDR<%.2g (min fold-change lfc=%.2f)", pcut, xcut)
              else sprintf("padj<%.2g, |log2FC|>%.1f", pcut, xcut)
  message(sprintf("%s: %s up, %s down, %s n.s. at %s",
                  title, comma(n_up), comma(n_down), comma(n_ns), gate_lab))
  message(sprintf("  Labelling %d genes (%d top per direction + curated anchors)",
                  nrow(label_df), N_TOP_DIR))

  x_axis_lab <- expression("log"[2]*" fold change")
  y_axis_lab <- if (mode == "treat") expression(-log[10]~italic("FDR"))
                else expression(-log[10]~"adjusted "*italic(P))

  x_lim <- max(abs(volc$xx), na.rm = TRUE) * 1.04
  y_lim <- max(volc$neglog10p, na.rm = TRUE) * 1.05

  ggplot(volc, aes(x = xx, y = neglog10p, color = status)) +
    geom_hline(yintercept = -log10(pcut),
               linetype = "dashed", color = "gray70", linewidth = 0.25) +
    geom_vline(xintercept = c(-xcut, xcut),
               linetype = "dashed", color = "gray70", linewidth = 0.25) +
    (if (VECTOR_MODE) {
      geom_point(data = volc[status == "n.s."],
                 size = 0.32, alpha = 0.45, shape = 16)
    } else {
      rasterize_layer(geom_point(data = volc[status == "n.s."],
                                 size = 0.32, alpha = 0.45, shape = 16))
    }) +
    (if (VECTOR_MODE) {
      geom_point(data = volc[status != "n.s."],
                 size = 0.50, alpha = 0.85, shape = 16)
    } else {
      rasterize_layer(geom_point(data = volc[status != "n.s."],
                                 size = 0.50, alpha = 0.85, shape = 16))
    }) +
    geom_point(data = label_df,
               aes(x = xx, y = neglog10p, fill = status),
               color = "black", shape = 21, size = 1.15,
               stroke = 0.25, inherit.aes = FALSE) +
    geom_text_repel(data = label_df,
                    aes(x = xx, y = neglog10p, label = symbol),
                    inherit.aes = FALSE,
                    size = GEOM_TEXT_6PT, color = "black", fontface = "italic",
                    segment.size = 0.2, segment.color = "gray45",
                    box.padding = 0.32, point.padding = 0.18,
                    min.segment.length = 0,
                    max.overlaps = Inf, force = 4, seed = 42,
                    show.legend = FALSE) +
    annotate("text", x = -x_lim * 0.98, y = y_lim * 0.97,
             label = sprintf("%s down", comma(n_down)),
             hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, fontface = "plain",
             color = masld_colors$down) +
    annotate("text", x =  x_lim * 0.98, y = y_lim * 0.97,
             label = sprintf("%s up", comma(n_up)),
             hjust = 1, vjust = 1, size = GEOM_TEXT_6PT, fontface = "plain",
             color = masld_colors$up) +
    scale_color_manual(values = volc_colors, guide = "none") +
    scale_fill_manual(values = volc_colors, guide = "none") +
    scale_x_continuous(limits = c(-x_lim, x_lim),
                       expand = expansion(mult = 0),
                       breaks = pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(0, y_lim),
                       expand = expansion(mult = c(0, 0)),
                       breaks = pretty_breaks(n = 5)) +
    labs(
      x = x_axis_lab,
      y = y_axis_lab
    ) +
    theme_masld(base_size = 6) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, color = "gray92"),
      legend.position  = "none",
      axis.title       = element_text(size = 6),
      plot.margin      = margin(4, 6, 2, 4)
    )
}

message("Loading canonical and definitive-sample DEG results...")
canonical_raw  <- fread(CANONICAL_PATH)
definitive_raw <- fread(DEFINITIVE_PATH)

# ── CANONICAL TREAT (raw logFC + treat FDR) — the headline panels ────────────
canonical_volc  <- prep_volcano(canonical_raw,  mode = "treat")
definitive_volc <- prep_volcano(definitive_raw, mode = "treat")
n_can <- sum(canonical_volc$status != "n.s.")
n_def <- sum(definitive_volc$status != "n.s.")
message("CAPTION: Disease-vs-control bulk RNA-seq (limma-voom quality-weighted C2, ",
        "~ dataset + sex + group). DEGs gated at adjusted FDR < 0.05 with a minimum ",
        "fold-change threshold (lfc = 0.25; effect floor folded into the significance test — ",
        "no separate |logFC| cut). ",
        "Left: pooled 5-cohort canonical call (n=846; 156 Ctrl / 690 Dis). ",
        "Right: definitive-sample sensitivity arm (same gate; Ctrl=NAS0/F0 or ",
        "labeled Control; Dis=NAS>=5 or F>=3). Gene labels italic; control series gray.")

p <- make_volcano(
  canonical_volc,
  sprintf("Canonical C2 (%s DEGs; adjusted FDR < 0.05, lfc = 0.25)", comma(n_can))
) |
  make_volcano(
    definitive_volc,
    sprintf("Definitive samples (%s DEGs; adjusted FDR < 0.05, lfc = 0.25)", comma(n_def))
  )

out_name <- if (VECTOR_MODE) "deg_volcano_vector.pdf" else "figs3a_deg_volcano.pdf"
out_pdf  <- file.path(PANEL_DIR, out_name)
save_fig(p, out_pdf, width = fig_full_width, height = 3.65)

# ── RAW sensitivity volcano (padj + raw logFC) — clearly-labelled supplement ──
canonical_raw_volc  <- prep_volcano(canonical_raw,  mode = "raw")
definitive_raw_volc <- prep_volcano(definitive_raw, mode = "raw")
n_can_r <- sum(canonical_raw_volc$status != "n.s.")
n_def_r <- sum(definitive_raw_volc$status != "n.s.")
message("CAPTION (sensitivity): same data on the RAW scale (unshrunk logFC, padj). ",
        "DEGs gated at raw |log2FC| > 0.5 & padj < 0.05.")
p_raw <- make_volcano(
  canonical_raw_volc,
  sprintf("Canonical C2 (%s DEGs; raw |log2FC| > 0.5, padj < 0.05; sensitivity)", comma(n_can_r))
) |
  make_volcano(
    definitive_raw_volc,
    sprintf("Definitive samples (%s DEGs; raw |log2FC| > 0.5, padj < 0.05; sensitivity)", comma(n_def_r))
  )
raw_name <- if (VECTOR_MODE) "deg_volcano_vector_raw_sensitivity.pdf" else "figs3a_deg_volcano_raw_sensitivity.pdf"
raw_pdf  <- file.path(PANEL_DIR, raw_name)
save_fig(p_raw, raw_pdf, width = fig_full_width, height = 3.65)

# Companion CSV: every labelled gene (canonical TREAT panels).
if (!VECTOR_MODE) {
  labels_out <- rbindlist(list(
    cbind(panel = "Canonical C2",       label_genes(canonical_volc)),
    cbind(panel = "Definitive samples", label_genes(definitive_volc))
  ), use.names = TRUE, fill = TRUE)
  setnames(labels_out, c("xx", "pp"), c("logFC", "treat_fdr"))
  fwrite(labels_out[, .(panel, symbol, logFC, treat_fdr, status)],
         file.path(PANEL_DIR, "deg_volcano_labels.csv"))
}

for (f in c(out_pdf, raw_pdf)) {
  if (file.exists(f)) {
    message(sprintf("\nOutput: %s (%s)", f,
                    utils:::format.object_size(file.size(f), "auto")))
  }
}
