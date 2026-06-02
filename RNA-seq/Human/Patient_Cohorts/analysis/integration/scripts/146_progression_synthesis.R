#!/usr/bin/env Rscript
# 146_progression_synthesis.R
# ---------------------------------------------------------------------------
# Progression Synthesis: multi-panel publication figure + summary tables.
#
# Generates 8 figure panels from all progression analysis results
# (Scripts 130-145) and combines into a composite figure.
#
# Panels:
#   A — UpSet plot: DEG overlap across C1-C9 contrasts
#   B — Pathway heatmap: top 30 Hallmark pathways x contrasts (NES fill)
#   C — Progression volcano: C2 with progression-specific genes highlighted
#   D — Drug Venn: onset-reversing vs progression-reversing CGP signatures
#   E — Deconv comparison: attribution class proportions C1 vs C2
#   F — Mouse concordance heatmap: human contrasts x 5 mouse diets
#   G — Subtype validation: S1/S2 proportion across NASH vs NAFL
#   H — Top genes table: top 50 progression-specific genes
#
# All inputs are loaded with graceful skip-on-missing.  If a required file
# is absent, the corresponding panel renders as a placeholder.
#
# Output:
#   figures/supplementary/progression/panel_a_upset.pdf  ..  panel_h_top_genes.pdf
#   figures/supplementary/progression/progression_synthesis.pdf  (20 x 15 in, landscape)
#   results/progression/progression_synthesis_summary.csv
#
# Usage: Rscript 146_progression_synthesis.R
# SLURM: cpu, 4 CPU, 16GB RAM, ~10min
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(gridExtra)
})

set.seed(42)

cat("=== Script 146: Progression Synthesis ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
#  Paths
# ============================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R"))

# Also source publication_theme.R from figures for masld_colors / theme_masld
pub_theme_path <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(pub_theme_path)) {
  source(pub_theme_path)
}

INT      <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PROG_DIR <- file.path(INT, "results/progression")
INT_RES  <- file.path(INT, "results/integration")
SIG_DIR  <- file.path(INT, "results/disease_signatures")
ATLAS_DIR <- file.path(BASE, "RNA-seq/results/multi_evidence")

FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(PROG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
#  Colors (from publication_theme.R / theme_publication.R)
# ============================================================
# Use sanjana_colors if available from theme_publication.R, otherwise define
if (!exists("sanjana_colors")) {
  sanjana_colors <- c(
    Magenta  = "#e14b9d", Pink  = "#e35070", Purple = "#d358c7",
    Blue     = "#4baeef", Orange = "#e1b172", Green  = "#30d796",
    Grey     = "#b0b0b0", DarkGrey = "#606060"
  )
}

# Contrast display labels (canonical ordering)
CONTRAST_LABELS <- c(
  C1  = "MASLD vs Ctrl",
  C2  = "NASH vs NAFL",
  C3  = "F3-F4 vs F0-F2",
  C4  = "NAFL vs Ctrl",
  C5  = "NAS>=5 vs <5",
  C6  = "Extreme",
  C7a = "Steatosis",
  C7b = "Inflammation",
  C7c = "Ballooning",
  C8  = "Cirrhosis",
  C9  = "F2 inflection",
  C11 = "NASH vs Ctrl",
  C12 = "Early vs Late NASH",
  C13 = "NASH vs NAFL (fib-adj)",
  C17 = "Fibrosis ordinal"
)

# Onset contrasts vs progression contrasts (for classification)
ONSET_IDS <- c("C1", "C4", "C11")
PROGRESSION_IDS <- c("C2", "C3", "C5", "C12", "C13")

# Attribution class colors
attr_colors <- c(
  Hepatocyte_intrinsic = "#0D47A1",
  Composition_driven   = "#C2185B",
  Unmasked             = "#7B1FA2",
  Not_significant      = "#BDBDBD"
)

# ============================================================
#  Helper: safe file loader
# ============================================================
safe_load <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    cat(sprintf("  [SKIP] %s — file not found\n", label))
    return(NULL)
  }
  dt <- fread(path)
  cat(sprintf("  [OK]   %s: %d rows x %d cols\n", label, nrow(dt), ncol(dt)))
  dt
}

# ============================================================
#  Helper: placeholder panel
# ============================================================
make_placeholder <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             size = 2.5, color = "gray50") +
    theme_void()
}

# ============================================================
# 1. Load all results
# ============================================================
cat("=== Loading results ===\n")

gene_class  <- safe_load(file.path(PROG_DIR, "gene_progression_classification.csv"))
consensus   <- safe_load(file.path(PROG_DIR, "progression_consensus_matrix.csv"))
nes_matrix  <- safe_load(file.path(PROG_DIR, "progression_gsea_nes_matrix.csv"))
gsea_all    <- safe_load(file.path(PROG_DIR, "progression_gsea_all.csv"))
drug_all    <- safe_load(file.path(PROG_DIR, "progression_drug_reversal_all.csv"))
drug_comp   <- safe_load(file.path(PROG_DIR, "progression_drug_comparison.csv"))
deconv_c2   <- safe_load(file.path(PROG_DIR, "progression_deconv_attribution_c2.csv"))
concordance <- safe_load(file.path(PROG_DIR, "progression_concordance_summary.csv"))
subtype_int <- safe_load(file.path(PROG_DIR, "subtype_progression_interaction.csv"))
subtype_stg <- safe_load(file.path(PROG_DIR, "subtype_stage_enrichment.csv"))
tf_summary  <- safe_load(file.path(PROG_DIR, "progression_tf_activity_summary.csv"))
deconv_c1   <- safe_load(file.path(BASE, "RNA-seq/results/causal_inference/deconv_attribution_scores.csv"))
atlas       <- safe_load(file.path(ATLAS_DIR, "multi_evidence_atlas.csv"))

# Also try to load raw dream results for C1 and C2 (for volcano / upset)
dream_c1 <- safe_load(file.path(INT_RES, "dream_results.csv"))
dream_c2 <- safe_load(file.path(SIG_DIR, "nafl_vs_nash_dream.csv"))
if (is.null(dream_c2)) {
  dream_c2 <- safe_load(file.path(PROG_DIR, "c2_nafl_vs_nash_dream.csv"))
}

# Load additional contrast dream files for UpSet
dream_contrasts <- list()
contrast_files <- list(
  C1  = file.path(INT_RES, "dream_results.csv"),
  C2  = c(file.path(SIG_DIR, "nafl_vs_nash_dream.csv"),
          file.path(PROG_DIR, "c2_nafl_vs_nash_dream.csv")),
  C3  = file.path(PROG_DIR, "c3_adv_vs_early_fib_dream.csv"),
  C4  = file.path(PROG_DIR, "c4_nafl_vs_ctrl_dream.csv"),
  C5  = file.path(PROG_DIR, "c5_nas_ge5_vs_lt5_dream.csv"),
  C6  = file.path(PROG_DIR, "c6_extreme_endpoints_dream.csv"),
  C7a = file.path(PROG_DIR, "c7a_steatosis_ordinal_dream.csv"),
  C7b = file.path(PROG_DIR, "c7b_inflammation_ordinal_dream.csv"),
  C7c = file.path(PROG_DIR, "c7c_ballooning_ordinal_dream.csv"),
  C8  = file.path(PROG_DIR, "c8_cirrhosis_dream.csv"),
  C9  = file.path(PROG_DIR, "c9_f2_inflection_dream.csv"),
  C11 = file.path(PROG_DIR, "c11_nash_vs_ctrl_dream.csv"),
  C12 = file.path(PROG_DIR, "c12_early_vs_late_nash_dream.csv"),
  C13 = file.path(PROG_DIR, "c13_nash_vs_nafl_fib_adj_dream.csv"),
  C17 = file.path(PROG_DIR, "c17_fibrosis_ordinal_dream.csv")
)

for (cid in names(contrast_files)) {
  paths <- contrast_files[[cid]]
  loaded_ok <- FALSE
  for (p in paths) {
    if (file.exists(p)) {
      dt <- fread(p)
      # Harmonise padj column
      if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt)) {
        setnames(dt, "adj.P.Val", "padj")
      }
      dream_contrasts[[cid]] <- dt
      loaded_ok <- TRUE
      break
    }
  }
  if (!loaded_ok) {
    cat(sprintf("  [SKIP] %s dream file not found\n", cid))
  }
}
cat(sprintf("\n  Loaded %d / %d contrast dream files\n\n",
  length(dream_contrasts), length(contrast_files)))

# ============================================================
# Gene symbol lookup
# ============================================================
gene_map <- NULL
if (!is.null(atlas)) {
  gene_map <- atlas[!is.na(human_symbol) & human_symbol != "",
    .(ensembl_clean = sub("\\..*", "", ensembl_id), symbol = human_symbol)]
  gene_map <- unique(gene_map, by = "ensembl_clean")
}

add_symbols <- function(dt, gene_col = "gene") {
  if (is.null(dt) || !gene_col %in% names(dt)) return(dt)
  # If dt already has a symbol column, use it directly

  if ("symbol" %in% names(dt)) {
    dt[is.na(symbol) | symbol == "", symbol := sub("\\..*", "", get(gene_col))]
    return(dt)
  }
  if (is.null(gene_map)) return(dt)
  dt[, ensembl_clean := sub("\\..*", "", get(gene_col))]
  dt <- merge(dt, gene_map, by = "ensembl_clean", all.x = TRUE)
  dt[is.na(symbol), symbol := ensembl_clean]
  dt
}

# ============================================================
# 2. Generate panels
# ============================================================
cat("=== Generating panels ===\n\n")

# Pre-initialise all panels as placeholders
p_a <- make_placeholder("(a) UpSet: DEG overlap — data not available")
p_b <- make_placeholder("(b) Pathway heatmap — data not available")
p_c <- make_placeholder("(c) Progression volcano — data not available")
p_d <- make_placeholder("(d) Drug Venn — data not available")
p_e <- make_placeholder("(e) Deconv comparison — data not available")
p_f <- make_placeholder("(f) Mouse concordance — data not available")
p_g <- make_placeholder("(g) Subtype validation — data not available")
p_h <- make_placeholder("(h) Top progression genes — data not available")

# ------------------------------------------------------------------
# Panel A: UpSet plot — DEG overlap across contrasts
# ------------------------------------------------------------------
tryCatch({
  cat("Panel A: UpSet plot...\n")

  if (length(dream_contrasts) >= 3) {
    # Build binary membership lists
    deg_lists <- list()
    for (cid in names(dream_contrasts)) {
      dt <- dream_contrasts[[cid]]
      if ("padj" %in% names(dt) && "gene" %in% names(dt)) {
        sig_genes <- dt[padj < 0.1, gene]
        if (length(sig_genes) > 0) {
          deg_lists[[cid]] <- sig_genes
        }
      }
    }

    if (length(deg_lists) >= 3) {
      # Use UpSetR if available, else manual bar plot
      if (requireNamespace("UpSetR", quietly = TRUE)) {
        library(UpSetR)

        # Save upset directly to PDF (UpSetR uses base graphics)
        upset_pdf <- file.path(FIG_DIR, "panels", "panel_a_upset.pdf")
        cairo_pdf(upset_pdf, width = 7, height = 4.5)
        upset(fromList(deg_lists),
              sets = rev(names(deg_lists)),
              keep.order = TRUE,
              nsets = length(deg_lists),
              nintersects = 30,
              order.by = "freq",
              mainbar.y.label = "Intersection size",
              sets.x.label = "DEGs per contrast",
              text.scale = c(1.1, 1, 0.9, 0.9, 1.1, 0.8),
              point.size = 2,
              line.size = 0.5,
              mb.ratio = c(0.6, 0.4),
              main.bar.color = sanjana_colors["Magenta"],
              sets.bar.color = sanjana_colors["Blue"])
        dev.off()
        cat(sprintf("  Saved UpSet to %s\n", upset_pdf))

        # Create a ggplot placeholder for the composite
        p_a <- make_placeholder("(a) UpSet plot — see panel_a_upset.pdf")
      } else {
        # Fallback: bar plot of DEG counts per contrast
        deg_counts <- data.table(
          contrast = names(deg_lists),
          n_degs = sapply(deg_lists, length)
        )
        deg_counts[, contrast := factor(contrast, levels = names(CONTRAST_LABELS))]
        p_a <- ggplot(deg_counts, aes(x = contrast, y = n_degs, fill = contrast)) +
          geom_col(show.legend = FALSE) +
          scale_fill_manual(values = rep(sanjana_colors["Blue"], nrow(deg_counts))) +
          labs(title = "(a) DEGs per contrast", x = NULL, y = "N DEGs (padj < 0.1)") +
          theme_publication() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
  } else if (!is.null(consensus)) {
    # Fallback from consensus matrix
    sig_cols <- grep("^sig_", names(consensus), value = TRUE)
    if (length(sig_cols) >= 3) {
      deg_counts <- data.table(
        contrast = sub("^sig_", "", sig_cols),
        n_degs = sapply(sig_cols, function(col) sum(consensus[[col]] == 1, na.rm = TRUE))
      )
      deg_counts[, contrast := toupper(contrast)]
      p_a <- ggplot(deg_counts, aes(x = contrast, y = n_degs)) +
        geom_col(fill = sanjana_colors["Blue"]) +
        labs(title = "(a) DEGs per contrast", x = NULL, y = "N DEGs (padj < 0.1)") +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }
}, error = function(e) {
  cat(sprintf("  Panel A error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel B: Pathway heatmap — top 30 Hallmark by significance
# ------------------------------------------------------------------
tryCatch({
  cat("Panel B: Pathway heatmap...\n")

  if (!is.null(nes_matrix)) {
    # nes_matrix: rows = pathways, first col = pathway name, rest = contrasts (NES values)
    nm <- copy(nes_matrix)
    pathway_col <- names(nm)[1]
    setnames(nm, pathway_col, "pathway")

    # Clean pathway names
    nm[, pathway := gsub("^HALLMARK_", "", pathway)]
    nm[, pathway := gsub("_", " ", pathway)]

    contrast_cols <- setdiff(names(nm), "pathway")

    # Melt to long
    nm_long <- melt(nm, id.vars = "pathway", variable.name = "contrast",
                    value.name = "NES")
    nm_long[, NES := as.numeric(NES)]

    # Rank pathways by max absolute NES across contrasts
    pw_rank <- nm_long[, .(max_abs_nes = max(abs(NES), na.rm = TRUE)), by = pathway]
    pw_rank <- pw_rank[order(-max_abs_nes)]
    top30 <- head(pw_rank$pathway, 30)

    nm_plot <- nm_long[pathway %in% top30]
    nm_plot[, pathway := factor(pathway, levels = rev(top30))]

    # Map contrast names to labels
    nm_plot[, contrast_label := CONTRAST_LABELS[as.character(contrast)]]
    nm_plot[is.na(contrast_label), contrast_label := as.character(contrast)]

    p_b <- ggplot(nm_plot, aes(x = contrast_label, y = pathway, fill = NES)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = sanjana_colors["Blue"],
                           mid = "white",
                           high = sanjana_colors["Magenta"],
                           midpoint = 0, na.value = "grey90",
                           name = "NES") +
      labs(title = "(b) Hallmark pathway NES across contrasts",
           x = NULL, y = NULL) +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5),
            legend.key.height = unit(0.4, "cm"),
            legend.key.width  = unit(0.2, "cm"))
  } else if (!is.null(gsea_all)) {
    # Build NES matrix from gsea_all (has columns: pathway, contrast, NES, padj, collection)
    hallmark_gsea <- gsea_all[grepl("Hallmark", collection, ignore.case = TRUE) |
                              grepl("^HALLMARK_", pathway)]
    if (nrow(hallmark_gsea) > 0) {
      # Rank by minimum padj across contrasts
      pw_rank <- hallmark_gsea[, .(min_padj = min(padj, na.rm = TRUE)), by = pathway]
      pw_rank <- pw_rank[order(min_padj)]
      top30 <- head(pw_rank$pathway, 30)

      hm_plot <- hallmark_gsea[pathway %in% top30]
      hm_plot[, pathway := gsub("^HALLMARK_", "", pathway)]
      hm_plot[, pathway := gsub("_", " ", pathway)]
      hm_plot[, pathway := factor(pathway, levels = rev(unique(pathway[order(-NES)])))]

      hm_plot[, contrast_label := CONTRAST_LABELS[as.character(contrast)]]
      hm_plot[is.na(contrast_label), contrast_label := as.character(contrast)]

      p_b <- ggplot(hm_plot, aes(x = contrast_label, y = pathway, fill = NES)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(low = sanjana_colors["Blue"],
                             mid = "white",
                             high = sanjana_colors["Magenta"],
                             midpoint = 0, na.value = "grey90",
                             name = "NES") +
        labs(title = "(b) Hallmark pathway NES across contrasts",
             x = NULL, y = NULL) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
              axis.text.y = element_text(size = 5),
              legend.key.height = unit(0.4, "cm"),
              legend.key.width  = unit(0.2, "cm"))
    }
  }
}, error = function(e) {
  cat(sprintf("  Panel B error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel C: Progression volcano — C2 (NASH vs NAFL)
# ------------------------------------------------------------------
tryCatch({
  cat("Panel C: Progression volcano...\n")

  if (!is.null(dream_c2)) {
    dt <- copy(dream_c2)
    # Harmonise padj
    if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt)) {
      setnames(dt, "adj.P.Val", "padj")
    }

    # Identify progression-specific genes: sig in C2 but NOT in C1
    c1_sig <- character(0)
    if (!is.null(dream_c1)) {
      c1_dt <- copy(dream_c1)
      if ("adj.P.Val" %in% names(c1_dt) && !"padj" %in% names(c1_dt)) {
        setnames(c1_dt, "adj.P.Val", "padj")
      }
      c1_sig <- c1_dt[padj < 0.1, gene]
    }

    dt[, nlog10p := pmin(-log10(padj), 50)]  # cap at 50

    # Classification: progression-specific (sig in C2, not in C1)
    dt[, category := "Not significant"]
    dt[padj < 0.1 & gene %in% c1_sig, category := "Shared (C1 & C2)"]
    dt[padj < 0.1 & !gene %in% c1_sig, category := "Progression-specific"]

    vol_colors <- c(
      "Progression-specific" = unname(sanjana_colors["Magenta"]),
      "Shared (C1 & C2)"     = unname(sanjana_colors["Blue"]),
      "Not significant"      = unname(sanjana_colors["Grey"])
    )

    # Add symbols for labeling top genes
    dt <- add_symbols(dt)

    # Label top 15 progression-specific genes by significance
    prog_top <- dt[category == "Progression-specific"][order(padj)]
    if (nrow(prog_top) == 0) {
      # Fallback: label top C2-significant genes
      prog_top <- dt[padj < 0.1][order(padj)][1:min(15, sum(dt$padj < 0.1, na.rm = TRUE))]
    } else {
      prog_top <- head(prog_top, 15)
    }

    # Only include categories present in the data
    present_cats <- unique(dt$category)
    vol_colors <- vol_colors[names(vol_colors) %in% present_cats]

    p_c <- ggplot(dt, aes(x = logFC, y = nlog10p, color = category)) +
      geom_point(size = 0.3, alpha = 0.6) +
      scale_color_manual(values = vol_colors, name = NULL) +
      labs(title = "(c) Progression volcano: NASH vs NAFL",
           x = expression(log[2]~FC),
           y = expression(-log[10]~padj)) +
      theme_publication() +
      theme(legend.position = c(0.02, 0.98),
            legend.justification = c(0, 1),
            legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
            legend.key.size = unit(0.3, "cm"))

    if (nrow(prog_top) > 0 && "symbol" %in% names(prog_top)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p_c <- p_c +
          ggrepel::geom_text_repel(
            data = prog_top,
            aes(label = symbol),
            size = 1.8, color = "black",
            max.overlaps = 12,
            segment.size = 0.2,
            segment.color = "grey60",
            min.segment.length = 0)
      }
    }

    # Annotate counts
    n_prog <- sum(dt$category == "Progression-specific", na.rm = TRUE)
    n_shared <- sum(dt$category == "Shared (C1 & C2)", na.rm = TRUE)
    p_c <- p_c +
      annotate("text", x = max(dt$logFC, na.rm = TRUE) * 0.7,
               y = max(dt$nlog10p, na.rm = TRUE) * 0.15,
               label = sprintf("Prog.-specific: %s\nShared: %s",
                               formatC(n_prog, format = "d", big.mark = ","),
                               formatC(n_shared, format = "d", big.mark = ",")),
               size = 2, hjust = 0.5, color = "grey30")
  }
}, error = function(e) {
  cat(sprintf("  Panel C error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel D: Drug Venn — onset vs progression CGP signatures
# ------------------------------------------------------------------
tryCatch({
  cat("Panel D: Drug Venn...\n")

  if (!is.null(drug_all)) {
    # Identify reversal hits (NES < 0, padj < 0.05) per contrast
    drug_all_copy <- copy(drug_all)

    # Standardise column names — look for contrast/contrast_id column
    cid_col <- intersect(c("contrast", "contrast_id", "short_label"), names(drug_all_copy))[1]
    if (!is.na(cid_col)) {
      setnames(drug_all_copy, cid_col, "cid", skip_absent = TRUE)
    }

    nes_col <- intersect(c("NES", "nes"), names(drug_all_copy))[1]
    padj_col <- intersect(c("padj", "adj.P.Val", "pval"), names(drug_all_copy))[1]
    pw_col <- intersect(c("pathway", "gene_set", "gs_name"), names(drug_all_copy))[1]

    if (!any(is.na(c(nes_col, padj_col, pw_col))) && "cid" %in% names(drug_all_copy)) {
      reversal <- drug_all_copy[get(nes_col) < 0 & get(padj_col) < 0.05]

      # Onset: C1; Progression: C2, C3, C5
      onset_drugs <- unique(reversal[grepl("C1|onset", cid, ignore.case = TRUE), get(pw_col)])
      prog_drugs  <- unique(reversal[grepl("C2|C3|C5|nash|fib|progression", cid, ignore.case = TRUE), get(pw_col)])

      n_onset_only <- length(setdiff(onset_drugs, prog_drugs))
      n_prog_only  <- length(setdiff(prog_drugs, onset_drugs))
      n_both       <- length(intersect(onset_drugs, prog_drugs))

      # Build a simple Venn as a ggplot (avoids VennDiagram/ggVennDiagram deps)
      venn_data <- data.table(
        category = c("Onset only", "Both", "Progression only"),
        count = c(n_onset_only, n_both, n_prog_only),
        x = c(-1.2, 0, 1.2),
        y = c(0, 0, 0)
      )

      if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
        venn_list <- list(
          "Onset (C1)" = onset_drugs,
          "Progression (C2/C3/C5)" = prog_drugs
        )
        p_d <- ggVennDiagram::ggVennDiagram(venn_list,
                 label = "count", label_alpha = 0,
                 edge_size = 0.5) +
          scale_fill_gradient(low = "white", high = sanjana_colors["Magenta"]) +
          labs(title = "(d) Onset vs progression reversal drugs") +
          theme_void() +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0),
                legend.position = "none")
      } else {
        # Manual Euler-style plot
        p_d <- ggplot() +
          # Left circle (onset)
          annotate("path",
            x = -0.4 + 1.2 * cos(seq(0, 2*pi, length.out = 100)),
            y = 1.2 * sin(seq(0, 2*pi, length.out = 100)),
            color = sanjana_colors["Blue"], linewidth = 0.8) +
          # Right circle (progression)
          annotate("path",
            x = 0.4 + 1.2 * cos(seq(0, 2*pi, length.out = 100)),
            y = 1.2 * sin(seq(0, 2*pi, length.out = 100)),
            color = sanjana_colors["Magenta"], linewidth = 0.8) +
          # Labels
          annotate("text", x = -1.1, y = 0, label = n_onset_only,
                   size = 4, fontface = "bold") +
          annotate("text", x = 0, y = 0, label = n_both,
                   size = 4, fontface = "bold") +
          annotate("text", x = 1.1, y = 0, label = n_prog_only,
                   size = 4, fontface = "bold") +
          annotate("text", x = -1.1, y = 1.5, label = "Onset (C1)",
                   size = 2.5, color = sanjana_colors["Blue"]) +
          annotate("text", x = 1.1, y = 1.5, label = "Progression",
                   size = 2.5, color = sanjana_colors["Magenta"]) +
          coord_fixed(xlim = c(-2, 2), ylim = c(-1.8, 2)) +
          labs(title = "(d) Onset vs progression reversal drugs") +
          theme_void() +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0))
      }
    }
  } else if (!is.null(drug_comp)) {
    # Use pre-computed comparison table
    p_d <- make_placeholder("(d) Drug comparison — see progression_drug_comparison.csv")
  }
}, error = function(e) {
  cat(sprintf("  Panel D error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel E: Deconv comparison — C1 vs C2 attribution proportions
# ------------------------------------------------------------------
tryCatch({
  cat("Panel E: Deconv comparison...\n")

  if (!is.null(deconv_c2) || !is.null(deconv_c1)) {
    # Try to classify genes for each contrast
    build_attr_props <- function(dt, label) {
      # Look for attribution_class column
      class_col <- intersect(c("attribution_class", "class", "deconv_class", "category"), names(dt))
      if (length(class_col) == 0) return(NULL)
      class_col <- class_col[1]

      # Harmonise class names
      dt[, attr_class := get(class_col)]
      dt[grepl("intrinsic|hepatocyte", attr_class, ignore.case = TRUE),
         attr_class := "Hepatocyte_intrinsic"]
      dt[grepl("composition|driven", attr_class, ignore.case = TRUE),
         attr_class := "Composition_driven"]
      dt[grepl("unmasked", attr_class, ignore.case = TRUE),
         attr_class := "Unmasked"]
      dt[grepl("not_sig|NS|not_significant", attr_class, ignore.case = TRUE),
         attr_class := "Not_significant"]

      counts <- dt[, .N, by = attr_class]
      counts[, contrast := label]
      counts[, prop := N / sum(N)]
      counts
    }

    prop_list <- list()
    if (!is.null(deconv_c1)) {
      c1_props <- build_attr_props(copy(deconv_c1), "C1: MASLD vs Ctrl")
      if (!is.null(c1_props)) prop_list[["C1"]] <- c1_props
    }
    if (!is.null(deconv_c2)) {
      c2_props <- build_attr_props(copy(deconv_c2), "C2: NASH vs NAFL")
      if (!is.null(c2_props)) prop_list[["C2"]] <- c2_props
    }

    if (length(prop_list) > 0) {
      props <- rbindlist(prop_list, fill = TRUE)
      props[, attr_class := factor(attr_class,
        levels = c("Hepatocyte_intrinsic", "Composition_driven",
                   "Unmasked", "Not_significant"))]

      p_e <- ggplot(props, aes(x = contrast, y = prop, fill = attr_class)) +
        geom_col(width = 0.6, color = "white", linewidth = 0.3) +
        scale_fill_manual(values = attr_colors, name = "Attribution",
                          drop = FALSE) +
        scale_y_continuous(labels = scales::percent_format(),
                           expand = expansion(mult = c(0, 0.02))) +
        labs(title = "(e) Deconvolution attribution: C1 vs C2",
             x = NULL, y = "Proportion of DEGs") +
        theme_publication() +
        theme(legend.position = "right",
              legend.key.size = unit(0.3, "cm"))
    }
  }
}, error = function(e) {
  cat(sprintf("  Panel E error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel F: Mouse concordance heatmap
# ------------------------------------------------------------------
tryCatch({
  cat("Panel F: Mouse concordance heatmap...\n")

  if (!is.null(concordance)) {
    conc <- copy(concordance)

    # Expect columns like: human_contrast, mouse_diet, rho_all, rho_sig,
    #   direction_concordance, jaccard, etc.
    # Identify key columns
    hc_col <- intersect(c("human_contrast", "contrast", "contrast_id"), names(conc))[1]
    md_col <- intersect(c("mouse_diet", "diet", "diet_model"), names(conc))[1]
    rho_col <- intersect(c("rho_all", "spearman_rho", "rho"), names(conc))[1]

    if (!any(is.na(c(hc_col, md_col, rho_col)))) {
      setnames(conc, c(hc_col, md_col, rho_col),
               c("human_contrast", "mouse_diet", "rho"), skip_absent = TRUE)

      conc[, rho := as.numeric(rho)]

      # Map contrast labels
      conc[, contrast_label := CONTRAST_LABELS[as.character(human_contrast)]]
      conc[is.na(contrast_label), contrast_label := as.character(human_contrast)]

      p_f <- ggplot(conc, aes(x = mouse_diet, y = contrast_label, fill = rho)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.2f", rho)), size = 2, color = "black") +
        scale_fill_gradient2(low = sanjana_colors["Blue"],
                             mid = "white",
                             high = sanjana_colors["Magenta"],
                             midpoint = 0,
                             limits = c(-1, 1),
                             name = expression(rho)) +
        labs(title = "(f) Cross-species concordance (Spearman rho)",
             x = "Mouse diet model", y = NULL) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.key.height = unit(0.4, "cm"),
              legend.key.width  = unit(0.2, "cm"))
    }
  }
}, error = function(e) {
  cat(sprintf("  Panel F error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel G: Subtype validation — S1/S2 in NASH vs NAFL
# ------------------------------------------------------------------
tryCatch({
  cat("Panel G: Subtype validation...\n")

  if (!is.null(subtype_stg)) {
    stg <- copy(subtype_stg)

    # Expect: stage/diagnosis, subtype/nmf_subtype, n/count
    stg_col <- intersect(c("diagnosis", "stage", "group"), names(stg))[1]
    sub_col <- intersect(c("nmf_subtype", "subtype", "Subtype"), names(stg))[1]
    n_col   <- intersect(c("n", "count", "N"), names(stg))[1]

    if (!any(is.na(c(stg_col, sub_col, n_col)))) {
      setnames(stg, c(stg_col, sub_col, n_col),
               c("diagnosis", "subtype", "n"), skip_absent = TRUE)
      stg[, n := as.numeric(n)]

      # Keep only NAFL and NASH rows (or equivalent)
      stg_filt <- stg[grepl("NAFL|NASH|nafl|nash", diagnosis, ignore.case = TRUE)]
      if (nrow(stg_filt) == 0) stg_filt <- stg  # fallback: use all

      stg_filt[, total := sum(n), by = diagnosis]
      stg_filt[, prop := n / total]

      subtype_cols <- c("S1" = sanjana_colors["Blue"], "S2" = sanjana_colors["Magenta"])

      p_g <- ggplot(stg_filt, aes(x = diagnosis, y = prop, fill = subtype)) +
        geom_col(width = 0.6, color = "white", linewidth = 0.3) +
        scale_fill_manual(values = subtype_cols, name = "Subtype") +
        scale_y_continuous(labels = scales::percent_format(),
                           expand = expansion(mult = c(0, 0.02))) +
        labs(title = "(g) NMF subtype by diagnosis",
             x = NULL, y = "Proportion") +
        theme_publication() +
        theme(legend.position = "right",
              legend.key.size = unit(0.3, "cm"))

      # Add chi-squared p-value if subtype_int has it
      if (!is.null(subtype_int)) {
        chi_row <- subtype_int[grepl("chi.*nafl|nafl.*chi|diagnosis.*chi",
                                     paste(names(subtype_int), collapse = "|"),
                                     ignore.case = TRUE)]
        # Try extracting p-value from any column named pvalue/p_value/p.value
        pval_col <- intersect(c("p_value", "pvalue", "p.value", "p"), names(subtype_int))
        test_col <- intersect(c("test", "analysis", "comparison"), names(subtype_int))

        if (length(pval_col) > 0 && length(test_col) > 0) {
          chi_rows <- subtype_int[grepl("chi|diagnosis|nafl", get(test_col[1]),
                                        ignore.case = TRUE)]
          if (nrow(chi_rows) > 0) {
            chi_p <- chi_rows[[pval_col[1]]][1]
            p_g <- p_g +
              annotate("text",
                x = 1.5, y = max(stg_filt$prop) * 1.08,
                label = sprintf("Chi-sq p = %s",
                  formatC(chi_p, format = "e", digits = 2)),
                size = 2.2, color = "grey30")
          }
        }
      }
    }
  }
}, error = function(e) {
  cat(sprintf("  Panel G error: %s\n", conditionMessage(e)))
})

# ------------------------------------------------------------------
# Panel H: Top 50 progression-specific genes (table)
# ------------------------------------------------------------------
tryCatch({
  cat("Panel H: Top genes table...\n")

  if (!is.null(gene_class)) {
    gc <- copy(gene_class)

    # Filter to progression-specific or progression_only genes
    class_col <- intersect(c("class", "gene_class", "classification", "progression_class"),
                           names(gc))[1]
    if (!is.na(class_col)) {
      setnames(gc, class_col, "gene_class", skip_absent = TRUE)

      prog_genes <- gc[grepl("progression", gene_class, ignore.case = TRUE)]
      if (nrow(prog_genes) == 0) {
        # Fallback: use all genes, sort by tau or n_contrasts
        prog_genes <- gc
      }
    } else {
      prog_genes <- gc
    }

    # Add C2 stats if available
    if (!is.null(dream_c2)) {
      c2_cols <- intersect(c("gene", "logFC", "padj"), names(dream_c2))
      if (length(c2_cols) == 3) {
        c2_sub <- dream_c2[, ..c2_cols]
        if ("adj.P.Val" %in% names(dream_c2)) {
          c2_sub <- dream_c2[, .(gene, logFC, padj = adj.P.Val)]
        }
        setnames(c2_sub, c("logFC", "padj"), c("logFC_C2", "padj_C2"))
        prog_genes <- merge(prog_genes, c2_sub, by = "gene", all.x = TRUE)
      }
    }

    # Add TWAS hit from atlas
    if (!is.null(atlas) && "twas_pval" %in% names(atlas)) {
      twas_lookup <- atlas[, .(ensembl_id, twas_pval)]
      twas_lookup[, ensembl_clean := sub("\\..*", "", ensembl_id)]
      twas_lookup[, TWAS_hit := fifelse(!is.na(twas_pval) & twas_pval < 0.05, "Yes", "No")]
      prog_genes[, ensembl_clean := sub("\\..*", "", gene)]
      prog_genes <- merge(prog_genes, twas_lookup[, .(ensembl_clean, TWAS_hit)],
                          by = "ensembl_clean", all.x = TRUE)
      prog_genes[is.na(TWAS_hit), TWAS_hit := ""]
    }

    # Add symbols
    prog_genes <- add_symbols(prog_genes)

    # Sort: prefer padj_C2, then tau, then n_contrasts
    sort_col <- intersect(c("padj_C2", "tau", "n_contrasts_sig", "n_sig"), names(prog_genes))
    if ("padj_C2" %in% sort_col) {
      prog_genes <- prog_genes[order(padj_C2)]
    } else if ("tau" %in% sort_col) {
      prog_genes <- prog_genes[order(-tau)]
    } else if (any(c("n_contrasts_sig", "n_sig") %in% sort_col)) {
      nc <- intersect(c("n_contrasts_sig", "n_sig"), names(prog_genes))[1]
      prog_genes <- prog_genes[order(-get(nc))]
    }

    top50 <- head(prog_genes, 50)

    # Select display columns
    display_cols <- intersect(
      c("symbol", "logFC_C2", "padj_C2", "n_contrasts_sig", "n_sig",
        "gene_class", "deconv_class", "attribution_class", "TWAS_hit"),
      names(top50)
    )
    if (length(display_cols) == 0) {
      display_cols <- names(top50)[1:min(6, ncol(top50))]
    }
    # Always include symbol if available
    if ("symbol" %in% names(top50) && !"symbol" %in% display_cols) {
      display_cols <- c("symbol", display_cols)
    }

    tbl <- top50[, ..display_cols]

    # Format numeric columns
    for (col in names(tbl)) {
      if (is.numeric(tbl[[col]])) {
        tbl[, (col) := signif(get(col), 3)]
      }
    }

    # Rename for display
    pretty_names <- c(
      symbol = "Gene",
      logFC_C2 = "logFC (C2)",
      padj_C2 = "padj (C2)",
      n_contrasts_sig = "N contrasts",
      n_sig = "N contrasts",
      gene_class = "Class",
      deconv_class = "Deconv",
      attribution_class = "Deconv",
      TWAS_hit = "TWAS"
    )
    for (i in seq_along(display_cols)) {
      if (display_cols[i] %in% names(pretty_names)) {
        setnames(tbl, display_cols[i], pretty_names[display_cols[i]],
                 skip_absent = TRUE)
      }
    }

    # Create table grob
    tt <- ttheme_minimal(
      base_size = 5,
      core    = list(fg_params = list(fontsize = 5)),
      colhead = list(fg_params = list(fontsize = 5.5, fontface = "bold"))
    )
    tbl_grob <- tableGrob(tbl, rows = NULL, theme = tt)
    p_h <- wrap_elements(tbl_grob)
  }
}, error = function(e) {
  cat(sprintf("  Panel H error: %s\n", conditionMessage(e)))
})

# ============================================================
# 3. Save individual panels as PDFs
# ============================================================
cat("\n=== Saving individual panel PDFs ===\n")

save_panel <- function(p, name, width = 5, height = 4) {
  out_path <- file.path(FIG_DIR, "panels", paste0(name, ".pdf"))
  tryCatch({
    pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
    ggsave(out_path, plot = p, device = pdf_device,
           width = width, height = height)
    cat(sprintf("  Saved %s\n", out_path))
  }, error = function(e) {
    cat(sprintf("  Failed to save %s: %s\n", name, conditionMessage(e)))
  })
}

save_panel(p_b, "panel_b_pathway_heatmap", width = 7, height = 5)
save_panel(p_c, "panel_c_volcano", width = 5, height = 4)
save_panel(p_d, "panel_d_drug_venn", width = 4.5, height = 3.5)
save_panel(p_e, "panel_e_deconv_comparison", width = 4.5, height = 3.5)
save_panel(p_f, "panel_f_concordance_heatmap", width = 5, height = 4)
save_panel(p_g, "panel_g_subtype_validation", width = 4, height = 3.5)
save_panel(p_h, "panel_h_top_genes", width = 7, height = 6)

# Panel A was saved by UpSetR directly; save placeholder/fallback
if (!inherits(p_a, "ggplot") || identical(p_a$data, data.frame())) {
  # Already saved by UpSetR
} else {
  save_panel(p_a, "panel_a_upset", width = 7, height = 4.5)
}

# ============================================================
# 4. Composite figure (patchwork)
# ============================================================
cat("\n=== Assembling composite figure ===\n")

tryCatch({
  # Layout: 4 rows x 2 columns
  # Row 1: A (upset) | B (heatmap)
  # Row 2: C (volcano) | D (venn)
  # Row 3: E (deconv) | F (concordance)
  # Row 4: G (subtype) | H (table)

  composite <- (p_a | p_b) /
               (p_c | p_d) /
               (p_e | p_f) /
               (p_g | p_h) +
    plot_annotation(
      title = "Progression Synthesis",
      tag_levels = "a",
      theme = theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        plot.tag   = element_text(size = 8, face = "bold")
      )
    )

  composite_path <- file.path(FIG_DIR, "progression_synthesis.pdf")
  pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
  ggsave(composite_path, plot = composite, device = pdf_device,
         width = 20, height = 15)
  cat(sprintf("  Saved composite: %s\n", composite_path))
}, error = function(e) {
  cat(sprintf("  Composite figure error: %s\n", conditionMessage(e)))
  cat("  Individual panels were saved separately.\n")
})

# ============================================================
# 5. Summary CSV
# ============================================================
cat("\n=== Generating summary CSV ===\n")

summary_rows <- list()

# Count progression-specific / onset-specific from gene_class
if (!is.null(gene_class)) {
  class_col <- intersect(c("class", "gene_class", "classification", "progression_class"),
                         names(gene_class))[1]
  if (!is.na(class_col)) {
    class_counts <- gene_class[, .N, by = class_col]
    setnames(class_counts, class_col, "category")
    class_counts[, metric := "gene_class_count"]
    summary_rows <- c(summary_rows, list(class_counts))
  }
}

# DEG counts per contrast
if (length(dream_contrasts) > 0) {
  deg_summary <- data.table(
    category = paste0(names(dream_contrasts), "_DEGs"),
    N = sapply(dream_contrasts, function(dt) {
      if ("padj" %in% names(dt)) sum(dt$padj < 0.1, na.rm = TRUE) else NA_integer_
    }),
    metric = "contrast_deg_count"
  )
  summary_rows <- c(summary_rows, list(deg_summary))
}

# Drug reversal counts
if (!is.null(drug_all)) {
  nes_col <- intersect(c("NES", "nes"), names(drug_all))[1]
  padj_col <- intersect(c("padj", "adj.P.Val", "pval"), names(drug_all))[1]
  cid_col <- intersect(c("contrast", "contrast_id", "short_label", "cid"), names(drug_all))[1]

  if (!any(is.na(c(nes_col, padj_col, cid_col)))) {
    rev_counts <- drug_all[get(nes_col) < 0 & get(padj_col) < 0.05,
                           .N, by = get(cid_col)]
    if (nrow(rev_counts) > 0) {
      setnames(rev_counts, "get", "category")
      rev_counts[, category := paste0(category, "_reversal_hits")]
      rev_counts[, metric := "drug_reversal_count"]
      summary_rows <- c(summary_rows, list(rev_counts))
    }
  }
}

# Concordance stats
if (!is.null(concordance)) {
  rho_col <- intersect(c("rho_all", "spearman_rho", "rho"), names(concordance))[1]
  if (!is.na(rho_col)) {
    concordance[, rho_val := as.numeric(get(rho_col))]
    conc_summary <- data.table(
      category = c("concordance_mean_rho", "concordance_max_rho"),
      N = c(mean(concordance$rho_val, na.rm = TRUE),
            max(concordance$rho_val, na.rm = TRUE)),
      metric = "concordance_stat"
    )
    summary_rows <- c(summary_rows, list(conc_summary))
  }
}

# Deconv attribution breakdown
if (!is.null(deconv_c2)) {
  class_col <- intersect(c("attribution_class", "class", "deconv_class"), names(deconv_c2))
  if (length(class_col) > 0) {
    attr_counts <- deconv_c2[, .N, by = class_col[1]]
    setnames(attr_counts, class_col[1], "category")
    attr_counts[, category := paste0("C2_deconv_", category)]
    attr_counts[, metric := "deconv_attribution"]
    summary_rows <- c(summary_rows, list(attr_counts))
  }
}

if (length(summary_rows) > 0) {
  summary_dt <- rbindlist(summary_rows, fill = TRUE)
  summary_path <- file.path(PROG_DIR, "progression_synthesis_summary.csv")
  fwrite(summary_dt, summary_path)
  cat(sprintf("  Saved summary: %s (%d rows)\n", summary_path, nrow(summary_dt)))
} else {
  cat("  No summary data available (upstream results not yet generated)\n")
}

# ============================================================
# Done
# ============================================================
cat(sprintf("\n=== Script 146 complete: %s ===\n", as.character(Sys.time())))
cat(sprintf("  Figure directory: %s\n", FIG_DIR))
cat(sprintf("  Results directory: %s\n", PROG_DIR))
cat("  NOTE: Panels with missing upstream data render as placeholders.\n")
cat("  Re-run after Scripts 130-145 complete for full figure.\n")
