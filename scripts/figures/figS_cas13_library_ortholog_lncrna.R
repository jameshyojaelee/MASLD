#!/usr/bin/env Rscript
# =============================================================================
# Figure S_lib_4 (ortholog mapping) and S_lib_7 (lncRNA targets)
# Cas13 library presentation supplementary figures
#
# KEY MESSAGE (S_lib_4): The 12-layer ortholog bridge recovers mappings across
#   biotypes and confidence tiers, with complementary coverage across methods.
# KEY MESSAGE (S_lib_7): lncRNA library targets are supported by multi-diet
#   mouse evidence, with canonical lncRNAs recovered by sequence-based methods.
#
# Outputs:
#   figures/supplementary/figS_cas13_library/04_ortholog_mapping.pdf
#   figures/supplementary/figS_cas13_library/S_lib_7_lncrna_targets.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(UpSetR)
  library(patchwork)
  library(grid)
  library(gridExtra)
})

# -- Load theme and paths -----------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

# -- Load data -----------------------------------------------------------------
ortho <- fread(file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz"))
lib   <- fread(file.path(BASE, "Cas13_Library_Design/data/cas13_library.csv"))
canon <- fread(file.path(BASE, "Cas13_Library_Design/results/ortholog_validation/canonical_lncrna_resolution.csv"))

# miRNA is OUT of scope for the Cas13 library: RfxCas13d cannot knock down mature
# ~22 nt miRNAs, and the miRNA tier was retired at v7 (library v8 = PC + lncRNA
# only). Drop miRNA orthologs so this bridge figure reflects the in-scope biotypes;
# the dedicated miRNA mapping layers (miRBase, MirGeneDB) are removed below.
ortho <- ortho[mouse_biotype != "miRNA" & human_biotype != "miRNA"]

# =============================================================================
# FIGURE S_lib_4: Ortholog mapping
# =============================================================================

# ---- Define layers -----------------------------------------------------------
# Two layers dropped from figure display 2026-05-28 (kept in master table for
# recall — display only):
#   - Plain (non-reciprocal) BLAST: tens of thousands of spurious single-method
#     hits at ~90 Myr divergence.
#   - Phase J synteny: 40k pairs, ~90% Gm-predicted, up to 976 mouse partners
#     per human anchor; 92% land at Tier L. Positional-only, not a reliable
#     orthology caller — ortho2align supersedes it for syntenic lncRNAs.
#   - Seekr k-mer: dropped 2026-05-28 (only 3 Tier-M pairs; near-zero recall on
#     canonical lncRNAs — compositional similarity too weak to be useful here).
layer_cols <- c(
  "tier_H_toga", "tier_H_biomart", "tier_H_orthofinder",
  "tier_M_blast_rbh",
  "tier_M_mmseqs2_rbh", "tier_M_ortho2align",
  "tier_M_noncode"
)

layer_labels <- c(
  tier_H_toga           = "TOGA",
  tier_H_biomart        = "biomaRt",
  tier_H_orthofinder    = "OrthoFinder",
  tier_M_blast_rbh      = "BLAST RBH",
  tier_M_mmseqs2_rbh    = "MMseqs2 RBH",
  tier_M_ortho2align    = "ortho2align",
  tier_M_noncode        = "NONCODE"
)

# Deduplicate to unique ortholog pairs
ortho_dedup <- unique(ortho, by = c("mouse_ensembl", "human_ensembl"))

# Ensure layer columns are integer (0/1) for UpSetR
for (col in layer_cols) {
  ortho_dedup[, (col) := as.integer(get(col))]
}

# ---- Panel A: UpSet plot (data prep) -----------------------------------------
upset_df <- as.data.frame(ortho_dedup[, ..layer_cols])
colnames(upset_df) <- layer_labels[layer_cols]

# ---- Panel B: Tier x biotype stacked bar -------------------------------------
ortho_dedup[, biotype_group := fcase(
  mouse_biotype == "protein_coding", "Protein-coding",
  mouse_biotype == "lncRNA",         "lncRNA",
  default = "Other"
)]
ortho_dedup[, biotype_group := factor(biotype_group,
                                       levels = c("Protein-coding", "lncRNA",
                                                  "Other"))]

tier_bio <- ortho_dedup[, .N, by = .(confidence_tier, biotype_group)]
tier_labels <- c("H" = "High", "M" = "Medium", "L" = "Low")
tier_bio[, confidence_tier := factor(tier_labels[confidence_tier],
                                     levels = tier_labels)]
tier_totals <- tier_bio[, .(total = sum(N)), by = confidence_tier]

biotype_pal <- c(
  "Protein-coding" = "#1565C0",
  "lncRNA"         = "#C9265E",
  "Other"          = "#9E9E9E"
)

# Only label segments large enough to be legible (>3% of tier total)
tier_bio <- merge(tier_bio, tier_totals, by = "confidence_tier")
tier_bio[, show_label := (N / total) > 0.03]

pB <- ggplot(tier_bio, aes(x = confidence_tier, y = N, fill = biotype_group)) +
  geom_col(width = 0.65) +
  geom_text(data = tier_bio[show_label == TRUE],
            aes(label = scales::comma(N)),
            position = position_stack(vjust = 0.5),
            size = PUB_GEOM_TEXT, color = "white") +
  geom_text(data = tier_totals,
            aes(x = confidence_tier, y = total, label = scales::comma(total)),
            vjust = -0.4, size = PUB_GEOM_TEXT + 0.3, fontface = "plain",
            inherit.aes = FALSE) +
  scale_fill_manual(values = biotype_pal, name = "Biotype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = scales::comma) +
  labs(x = NULL, y = "Ortholog pairs", tag = "B") +
  theme_masld() + theme_pub() +
  theme(legend.position = "right",
        plot.tag = element_text(size = 9, face = "plain"))

# ---- Panel C: Canonical lncRNA detection tile --------------------------------
lnc_method_cols   <- c("tier_M_blast_rbh", "tier_M_mmseqs2_rbh",
                       "tier_M_ortho2align", "tier_M_noncode")
lnc_method_labels <- c("BLAST RBH", "MMseqs2 RBH",
                       "ortho2align", "NONCODE")

canonical_genes <- c("MALAT1", "NEAT1", "MEG3", "XIST", "H19",
                     "HOTAIR", "KCNQ1OT1", "TERC", "PVT1", "NORAD", "HOTTIP")

canon_detect <- ortho_dedup[human_symbol %in% canonical_genes,
                            lapply(.SD, function(x) as.integer(any(x == 1, na.rm = TRUE))),
                            by = human_symbol,
                            .SDcols = lnc_method_cols]
setnames(canon_detect, lnc_method_cols, lnc_method_labels)

canon_melt <- melt(canon_detect, id.vars = "human_symbol",
                   variable.name = "method", value.name = "detected")
canon_melt[, human_symbol := factor(human_symbol, levels = rev(canonical_genes))]
canon_melt[, method := factor(method, levels = lnc_method_labels)]

pC <- ggplot(canon_melt, aes(x = method, y = human_symbol,
                              fill = factor(detected))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "#E0E0E0", "1" = "#C9265E"),
                    labels = c("Not detected", "Detected"),
                    name = "") +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(face = "italic", size = 6),
        legend.position = "bottom",
        panel.grid = element_blank(),
        plot.tag = element_text(size = 9, face = "plain"))

# ---- Compose S_lib_4 --------------------------------------------------------
# UpSetR draws via base graphics and always emits a blank first page through
# grid.arrange. We capture just the UpSet page and the B+C ggplot page, then
# strip the blanks with a targeted approach: open device with onefile=TRUE,
# render UpSet, then render B+C on a new page.
p_BC <- pB + pC + plot_layout(widths = c(1.2, 1))

out_path_slib4 <- file.path(OUT_DIR, "ortholog_mapping.pdf")
pdf_device(out_path_slib4, width = fig_full_width, height = 4.2, onefile = TRUE)

# Page 1: UpSet (UpSetR internally calls grid.newpage + grid.arrange)
upset(upset_df,
      sets = rev(unname(layer_labels)),
      nintersects = 20,
      nsets = length(layer_labels),
      order.by = "freq",
      decreasing = TRUE,
      show.numbers = "yes",
      number.angles = 0,
      point.size = 1.8,
      line.size = 0.5,
      mb.ratio = c(0.6, 0.4),
      text.scale = c(1.2, 1, 0.9, 0.9, 1.1, 0.8),
      mainbar.y.label = "Intersection size",
      sets.x.label = "Pairs per layer")

# Page 2: Panels B + C
grid.newpage()
print(p_BC)

invisible(dev.off())

# Strip blank pages: UpSetR always emits one blank page before the actual plot.
# Use gs (ghostscript) to extract only the content pages.
gs_bin <- Sys.which("gs")
if (nzchar(gs_bin)) {
  # Detect page count via gs
  n_pages <- tryCatch({
    info <- system2(gs_bin,
      c("-q", "-dNODISPLAY", "-dNOSAFER",
        paste0("-c \"(", out_path_slib4, ") (r) file runpdfbegin pdfpagecount = quit\"")),
      stdout = TRUE, stderr = FALSE)
    as.integer(trimws(info[length(info)]))
  }, error = function(e) NA_integer_)

  # Determine which pages to keep (UpSet content = page 2; B+C = last page)
  if (!is.na(n_pages) && n_pages >= 3L) {
    keep_pages <- c(2L, n_pages)  # page 2 = UpSet, last page = B+C
    tmp_out <- tempfile(fileext = ".pdf")
    # gs page selection: -dFirstPage / -dLastPage for each page, combined via
    # a two-pass extraction (gs only supports contiguous ranges natively).
    # Simplest: extract page 2, extract last page, merge.
    tmp_p1 <- tempfile(fileext = ".pdf")
    tmp_p2 <- tempfile(fileext = ".pdf")
    system2(gs_bin, c("-sDEVICE=pdfwrite", "-dNOPAUSE", "-dBATCH", "-dQUIET",
                      paste0("-dFirstPage=", keep_pages[1]),
                      paste0("-dLastPage=", keep_pages[1]),
                      paste0("-sOutputFile=", tmp_p1), out_path_slib4),
            stdout = FALSE, stderr = FALSE)
    system2(gs_bin, c("-sDEVICE=pdfwrite", "-dNOPAUSE", "-dBATCH", "-dQUIET",
                      paste0("-dFirstPage=", keep_pages[2]),
                      paste0("-dLastPage=", keep_pages[2]),
                      paste0("-sOutputFile=", tmp_p2), out_path_slib4),
            stdout = FALSE, stderr = FALSE)
    # Merge
    system2(gs_bin, c("-sDEVICE=pdfwrite", "-dNOPAUSE", "-dBATCH", "-dQUIET",
                      paste0("-sOutputFile=", tmp_out), tmp_p1, tmp_p2),
            stdout = FALSE, stderr = FALSE)
    if (file.exists(tmp_out) && file.size(tmp_out) > 0) {
      file.copy(tmp_out, out_path_slib4, overwrite = TRUE)
      cat("Stripped blank pages (", n_pages, "-> 2) via gs\n")
    }
    unlink(c(tmp_p1, tmp_p2, tmp_out))
  } else {
    cat("Page count is", n_pages, "-- skipping blank-page strip\n")
  }
} else {
  cat("gs not found -- blank UpSetR pages retained\n")
}

cat("Saved 04_ortholog_mapping.pdf\n")
cat("Done. (S_lib_7 lncRNA-targets figure removed from the library design per PI 2026-06-01.)\n")
quit(save = "no")   # stop here -- S_lib_7 section below is retained for provenance but not executed

# =============================================================================
# FIGURE S_lib_7: lncRNA targets  [REMOVED 2026-06-01 -- not generated]
# =============================================================================

# ---- Panel A: lncRNA evidence depth bar chart --------------------------------
lib_lnc <- lib[biotype == "lncRNA"]

# Evidence source for each lncRNA, given the v4 library definition
# (human ashr DEGs lfsr<0.05 & shrunk_logFC>0.2 UNION mouse cross-diet >=3 diets):
#   Human + Mouse DE : human DE AND mouse-replicated (>=3 diets)
#   Human DE only    : in via human ashr backbone; mouse DE absent / sub-threshold
#   Mouse DE only    : in via mouse cross-diet (>=3); no human DE
lib_lnc[, evidence_cat := fcase(
  has_human_de == TRUE & n_diets_up >= 3, "Human + Mouse DE",
  has_human_de == TRUE,                   "Human DE only",
  default =                               "Mouse DE only"
)]
lib_lnc[, evidence_cat := factor(evidence_cat,
                                  levels = c("Mouse DE only",
                                             "Human DE only",
                                             "Human + Mouse DE"))]

evid_counts <- lib_lnc[, .N, by = evidence_cat]
# Drop empty levels for clean display
evid_counts <- evid_counts[!is.na(evidence_cat) & N > 0]
total_lnc <- nrow(lib_lnc)

evid_pal <- c(
  "Mouse DE only"    = "#9E9E9E",
  "Human DE only"    = "#1565C0",
  "Human + Mouse DE" = "#C9265E"
)

p7A <- ggplot(evid_counts, aes(x = N, y = evidence_cat, fill = evidence_cat)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = N), hjust = -0.2, size = PUB_GEOM_TEXT + 0.5,
            fontface = "plain") +
  scale_fill_manual(values = evid_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Gene count", y = NULL) +
  theme_masld() + theme_pub()
message("[caption] lncRNA library targets (n = ", total_lnc, ")")

# ---- Panel B: Canonical lncRNA expanded annotation tile ----------------------
# Get tier for each canonical lncRNA
canon_tier <- ortho_dedup[human_symbol %in% canonical_genes,
                          .(confidence_tier = confidence_tier[1]),
                          by = human_symbol]

# Library inclusion (by human symbol or mouse symbol fallback)
canon_lib <- data.table(human_symbol = canonical_genes)
canon_lib[, in_library := human_symbol %in% lib[biotype == "lncRNA", gene_symbol_human]]

# n_diets_up from library (direct match)
lib_diets_h <- lib[gene_symbol_human %in% canonical_genes,
                   .(n_diets_up = max(n_diets_up)),
                   by = gene_symbol_human]
canon_lib <- merge(canon_lib, lib_diets_h,
                   by.x = "human_symbol", by.y = "gene_symbol_human", all.x = TRUE)
canon_lib[is.na(n_diets_up), n_diets_up := 0L]

# Fallback via mouse ortholog symbol
for (i in seq_len(nrow(canon_lib))) {
  if (canon_lib$n_diets_up[i] == 0L) {
    hsym <- canon_lib$human_symbol[i]
    msyms <- ortho_dedup[human_symbol == hsym, unique(mouse_symbol)]
    match_row <- lib[gene_symbol_mouse %in% msyms]
    if (nrow(match_row) > 0) {
      canon_lib[i, in_library := TRUE]
      canon_lib[i, n_diets_up := max(match_row$n_diets_up)]
    }
  }
}

# Merge all annotation sources
canon_annot <- merge(canon_tier, canon_detect, by = "human_symbol", all = TRUE)
canon_annot <- merge(canon_annot, canon_lib, by = "human_symbol", all = TRUE)
canon_annot[is.na(in_library), in_library := FALSE]
canon_annot[is.na(n_diets_up), n_diets_up := 0L]

# Build tile data: methods (binary) + annotations
method_melt <- melt(canon_annot, id.vars = "human_symbol",
                    measure.vars = lnc_method_labels,
                    variable.name = "feature", value.name = "value")
method_melt[, type := "method"]

tier_dt <- canon_annot[, .(human_symbol,
                           feature = "Tier",
                           value = ifelse(confidence_tier == "H", 1,
                                   ifelse(confidence_tier == "M", 0.5, 0)),
                           type = "annotation")]

lib_dt <- canon_annot[, .(human_symbol,
                          feature = "In library",
                          value = as.numeric(in_library),
                          type = "annotation")]

max_diets <- 4L   # 4 diet groups (NASH+Western merged 2026-05-29)
diets_dt <- canon_annot[, .(human_symbol,
                            feature = "n_diets_up",
                            value = n_diets_up / max_diets,
                            type = "annotation")]

all_tile <- rbindlist(list(method_melt, tier_dt, lib_dt, diets_dt))
all_tile[, human_symbol := factor(human_symbol, levels = rev(canonical_genes))]

feat_order  <- c(lnc_method_labels, "Tier", "In library", "n_diets_up")
feat_labels <- c(lnc_method_labels, "Tier", "In library", "Diets (n)")
all_tile[, feature := factor(feature, levels = feat_order, labels = feat_labels)]

# Discrete fill categories
all_tile[, fill_cat := fcase(
  type == "method" & value == 1,                              "detected",
  type == "method" & value == 0,                              "not_detected",
  type == "annotation" & feature == "Tier" & value == 1,      "tier_H",
  type == "annotation" & feature == "Tier" & value == 0.5,    "tier_M",
  type == "annotation" & feature == "Tier" & value == 0,      "tier_L",
  type == "annotation" & feature == "In library" & value == 1, "detected",
  type == "annotation" & feature == "In library" & value == 0, "not_detected",
  type == "annotation" & feature == "Diets (n)" & value > 0,  "diets_pos",
  type == "annotation" & feature == "Diets (n)" & value == 0, "not_detected",
  default = "not_detected"
)]

tile_pal <- c(
  "detected"     = "#C9265E",
  "not_detected" = "#E0E0E0",
  "tier_H"       = "#0D47A1",
  "tier_M"       = "#42A5F5",
  "tier_L"       = "#BDBDBD",
  "diets_pos"    = "#00695C"
)

# Text labels for tier and diets columns
all_tile[, text_label := fcase(
  type == "annotation" & feature == "Tier" & value == 1,   "H",
  type == "annotation" & feature == "Tier" & value == 0.5, "M",
  type == "annotation" & feature == "Tier" & value == 0,   "L",
  type == "annotation" & feature == "Diets (n)",
    as.character(as.integer(round(value * max_diets))),
  default = ""
)]

p7B <- ggplot(all_tile, aes(x = feature, y = human_symbol, fill = fill_cat)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = text_label), size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_manual(values = tile_pal,
                    labels = c(detected = "Detected / Yes",
                               not_detected = "Not detected / No",
                               tier_H = "Tier H",
                               tier_M = "Tier M",
                               tier_L = "Tier L",
                               diets_pos = "Diets > 0"),
                    name = "") +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(face = "italic", size = 6),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        panel.grid = element_blank())
message("[caption] Canonical lncRNA validation (expanded)")

# ---- Compose S_lib_7 --------------------------------------------------------
p_slib7 <- p7A / p7B +
  plot_layout(heights = c(0.8, 1.2)) +
  plot_annotation(tag_levels = "A")

save_fig(p_slib7,
         file.path(OUT_DIR, "S_lib_7_lncrna_targets.pdf"),
         width = fig_full_width, height = 6)

cat("Saved S_lib_7_lncrna_targets.pdf\n")
cat("Done.\n")
