#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity.R  (new 2026-04-27; dual-scale 2026-06-27)
# Supp figure: downstream conclusions are stable across |log2FC| ∈ {0.2, 0.3, 0.5}.
#
# Produced on BOTH effect-size scales (looped):
#   raw     : bulk_padj < 0.05 & |bulk_logFC|        > cut   (suffix _raw)
#   shrunk  : bulk_lfsr < 0.05 & |bulk_shrunk_logFC| > cut   (suffix _shrunk, PRIMARY)
#
# Five panels (this script builds A..E only; F_loo / I_loo are owned by the
# figS_lfc_sensitivity_loo*.R scripts and are not touched here):
#   A  COLOC × DEG overlap (counts + Fisher OR + p-value at each cutoff)
#   B  Conserved (human DEG ∩ mouse meta-DEG via ortholog) count
#   C  Triple convergence: DEG ∩ COLOC ∩ Conserved gene count
#   D  Pathway enrichment top-20 Jaccard (3x3) -- Hallmark ORA at each cutoff
#   E  Atlas "high-confidence multi-evidence" gene count at each cutoff
#
# Outputs (per scale, scale suffix before extension):
#   figS_lfc_sensitivity_{raw,shrunk}.pdf
#   panels/{coloc_overlap,conserved,triple_convergence,pathway_jaccard,
#           multi_evidence}_{raw,shrunk}.pdf
#   sensitivity_data_{raw,shrunk}.csv
#   pathway_jaccard_{raw,shrunk}.csv
#   top20_pathways_per_cutoff_{raw,shrunk}.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR  <- FIGS_LFCSENS_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

LFC_TIERS <- c(0.2, 0.3, 0.5)
PADJ <- 0.05
COLOC_PP4 <- 0.5

# -----------------------------------------------------------------------------
# Load atlas
# -----------------------------------------------------------------------------
cat("Loading multi-evidence atlas...\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat(sprintf("  atlas: %d genes x %d cols\n", nrow(atlas), ncol(atlas)))

# COLOC support flag (scale-independent; SuSiE PP4, fall back to ABF)
coloc_support <- {
  pp4 <- atlas$coloc_susie_best_pp4
  if (is.null(pp4) || all(is.na(pp4))) pp4 <- atlas$coloc_best_susie_pp4
  if (is.null(pp4) || all(is.na(pp4))) pp4 <- atlas$coloc_best_pp4
  !is.na(pp4) & pp4 > COLOC_PP4
}
cat(sprintf("  COLOC PP4 > %.2f: %d genes\n", COLOC_PP4, sum(coloc_support)))

# Mouse DEG flag (mouse meta has no ashr-shrunk arm; mirror human grid on the
# mouse meta effect size in both scales)
is_mouse_deg <- function(dt, lfc) {
  !is.na(dt$mouse_meta_padj) & dt$mouse_meta_padj < PADJ &
    !is.na(dt$mouse_meta_logFC) & abs(dt$mouse_meta_logFC) > lfc
}

# Hallmark via msigdbr (load once)
have_msig <- requireNamespace("msigdbr", quietly = TRUE)
if (have_msig) {
  msig <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"))
  if (!"gene_symbol" %in% names(msig)) setnames(msig, "human_gene_symbol", "gene_symbol")
} else {
  warning("msigdbr unavailable; pathway panel will be skipped")
}

pathway_top20 <- function(deg_genes, bg_genes, msig) {
  if (length(deg_genes) == 0) return(character(0))
  bg_size <- length(bg_genes)
  res <- msig[, {
    set <- intersect(.SD$gene_symbol, bg_genes)
    if (length(set) < 5) {
      .(p = NA_real_, n_overlap = 0L, n_set = length(set))
    } else {
      ov <- intersect(deg_genes, set)
      m <- length(set); k <- length(deg_genes); q <- length(ov); n <- bg_size - m
      p_val <- phyper(q - 1, m, n, k, lower.tail = FALSE)
      .(p = p_val, n_overlap = q, n_set = m)
    }
  }, by = gs_name]
  res <- res[!is.na(p) & n_overlap > 0]
  res <- res[order(p)][1:min(.N, 20), gs_name]
  res
}

# Scale configurations (raw + ashr-shrunk)
SCALES <- list(
  raw = list(suffix = "_raw", pcol = "bulk_padj", lcol = "bulk_logFC",
             title_lfc = "raw |log2FC|"),
  shrunk = list(suffix = "_shrunk", pcol = "bulk_lfsr", lcol = "bulk_shrunk_logFC",
                title_lfc = "ashr-shrunk |log2FC|")
)

# -----------------------------------------------------------------------------
# Per-scale builder
# -----------------------------------------------------------------------------
build_scale <- function(cfg) {
  sfx <- cfg$suffix; pcol <- cfg$pcol; lcol <- cfg$lcol; tlfc <- cfg$title_lfc
  cat("\n========== scale:", pcol, "/", lcol, "(suffix", sfx, ") ==========\n")

  is_human_deg <- function(dt, lfc) {
    !is.na(dt[[pcol]]) & dt[[pcol]] < PADJ &
      !is.na(dt[[lcol]]) & abs(dt[[lcol]]) > lfc
  }

  # ----- Panel A: COLOC x DEG overlap --------------------------------------
  panelA <- rbindlist(lapply(LFC_TIERS, function(lfc) {
    is_deg <- is_human_deg(atlas, lfc)
    tab <- table(deg = is_deg, coloc = coloc_support)
    fish <- fisher.test(tab)
    data.table(lfc = lfc, n_deg = sum(is_deg), n_coloc = sum(coloc_support),
               n_overlap = sum(is_deg & coloc_support),
               OR = unname(fish$estimate), p = fish$p.value,
               ci_lo = fish$conf.int[1], ci_hi = fish$conf.int[2])
  }))
  print(panelA)

  p_a <- ggplot(panelA, aes(x = factor(lfc), y = n_overlap)) +
    geom_col(fill = "#1565C0", width = 0.55) +
    geom_col(data = panelA[lfc == 0.5], aes(x = factor(lfc), y = n_overlap),
             fill = "#0D47A1", width = 0.55) +
    geom_text(aes(label = sprintf("%d\n(OR=%.2f, p=%.1e)", n_overlap, OR, p)),
              vjust = -0.2, size = GEOM_TEXT_6PT, lineheight = 0.95) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = expression("|log"[2]*"FC| cutoff"),
         y = "Genes (DEG ∩ COLOC PP4 > 0.5)") +
    theme_masld()

  # ----- Panel B: Conserved ------------------------------------------------
  panelB <- rbindlist(lapply(LFC_TIERS, function(lfc) {
    hd <- is_human_deg(atlas, lfc); md <- is_mouse_deg(atlas, lfc)
    data.table(lfc = lfc, n_human = sum(hd, na.rm = TRUE),
               n_mouse = sum(md, na.rm = TRUE),
               n_conserved = sum(hd & md, na.rm = TRUE))
  }))
  print(panelB)

  p_b <- ggplot(panelB, aes(x = factor(lfc), y = n_conserved)) +
    geom_col(aes(fill = factor(lfc)), width = 0.55) +
    geom_text(aes(label = format(n_conserved, big.mark = ",")),
              vjust = -0.3, size = GEOM_TEXT_6PT, fontface = "plain") +
    scale_fill_manual(values = c(`0.2` = "#90CAF9", `0.3` = "#1976D2",
                                 `0.5` = "#0D47A1"), guide = "none") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    labs(x = expression("|log"[2]*"FC| cutoff (both species)"),
         y = "Conserved genes") +
    theme_masld()

  # ----- Panel C: Triple convergence ---------------------------------------
  panelC <- rbindlist(lapply(LFC_TIERS, function(lfc) {
    hd <- is_human_deg(atlas, lfc); md <- is_mouse_deg(atlas, lfc)
    data.table(lfc = lfc, n_triple = sum(hd & md & coloc_support, na.rm = TRUE))
  }))
  print(panelC)

  p_c <- ggplot(panelC, aes(x = factor(lfc), y = n_triple)) +
    geom_col(aes(fill = factor(lfc)), width = 0.55) +
    geom_text(aes(label = n_triple), vjust = -0.3, size = GEOM_TEXT_6PT, fontface = "plain") +
    scale_fill_manual(values = c(`0.2` = "#90CAF9", `0.3` = "#1976D2",
                                 `0.5` = "#0D47A1"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = expression("|log"[2]*"FC| cutoff"), y = "Triple-convergent genes") +
    theme_masld()

  # ----- Panel D: Top-20 pathway Jaccard -----------------------------------
  if (have_msig) {
    bg <- atlas[!is.na(human_symbol) & human_symbol != "", unique(human_symbol)]
    top20_by_lfc <- lapply(LFC_TIERS, function(lfc) {
      deg_idx <- is_human_deg(atlas, lfc)
      deg_genes <- atlas[deg_idx & !is.na(human_symbol) & human_symbol != "",
                         unique(human_symbol)]
      pathway_top20(deg_genes, bg, msig)
    })
    names(top20_by_lfc) <- as.character(LFC_TIERS)

    jacc <- CJ(a = LFC_TIERS, b = LFC_TIERS)
    jacc[, jaccard := mapply(function(x, y) {
      sa <- top20_by_lfc[[as.character(x)]]; sb <- top20_by_lfc[[as.character(y)]]
      if (length(sa) == 0 || length(sb) == 0) return(NA_real_)
      length(intersect(sa, sb)) / length(union(sa, sb))
    }, a, b)]
    print(dcast(jacc, a ~ b, value.var = "jaccard"))

    p_d <- ggplot(jacc, aes(x = factor(a), y = factor(b), fill = jaccard)) +
      geom_tile(color = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.2f", jaccard)), size = GEOM_TEXT_6PT, color = "white") +
      scale_fill_viridis_c(option = "mako", limits = c(0, 1), name = "Top-20\nJaccard") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(x = expression("|log"[2]*"FC| cutoff (rows)"),
           y = expression("|log"[2]*"FC| cutoff (cols)")) +
      theme_masld() +
      theme(panel.grid = element_blank(), axis.ticks = element_blank())
  } else {
    p_d <- ggplot() + theme_masld()
    message("[caption] Pathway panel skipped (msigdbr unavailable)")
    jacc <- data.table(); top20_by_lfc <- list()
  }

  # ----- Panel E: Multi-evidence high-confidence count ---------------------
  panelE <- rbindlist(lapply(LFC_TIERS, function(lfc) {
    hd <- is_human_deg(atlas, lfc); md <- is_mouse_deg(atlas, lfc)
    is_cc <- !is.na(atlas$is_conserved) & atlas$is_conserved == TRUE
    mev <- hd & (coloc_support | md | is_cc)
    data.table(lfc = lfc, n_multi_evidence = sum(mev, na.rm = TRUE),
               pct_of_deg = if (sum(hd) > 0)
                 round(100 * sum(mev, na.rm = TRUE) / sum(hd), 1) else NA_real_)
  }))
  print(panelE)

  p_e <- ggplot(panelE, aes(x = factor(lfc), y = n_multi_evidence)) +
    geom_col(aes(fill = factor(lfc)), width = 0.55) +
    geom_text(aes(label = sprintf("%s\n(%.1f%% of DEGs)",
                                  format(n_multi_evidence, big.mark = ","), pct_of_deg)),
              vjust = -0.2, size = GEOM_TEXT_6PT, lineheight = 0.95) +
    scale_fill_manual(values = c(`0.2` = "#90CAF9", `0.3` = "#1976D2",
                                 `0.5` = "#0D47A1"), guide = "none") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.20))) +
    labs(x = expression("|log"[2]*"FC| cutoff"), y = "Multi-evidence DEGs") +
    theme_masld()

  # ----- Save individual panels --------------------------------------------
  ggsave(file.path(PANEL_DIR, paste0("coloc_overlap", sfx, ".pdf")), p_a,
         width = 5, height = 4, device = cairo_pdf)
  ggsave(file.path(PANEL_DIR, paste0("conserved", sfx, ".pdf")), p_b,
         width = 5, height = 4, device = cairo_pdf)
  ggsave(file.path(PANEL_DIR, paste0("triple_convergence", sfx, ".pdf")), p_c,
         width = 5, height = 4, device = cairo_pdf)
  ggsave(file.path(PANEL_DIR, paste0("pathway_jaccard", sfx, ".pdf")), p_d,
         width = 5, height = 4, device = cairo_pdf)
  ggsave(file.path(PANEL_DIR, paste0("multi_evidence", sfx, ".pdf")), p_e,
         width = 5, height = 4, device = cairo_pdf)

  # ----- Combined composite ------------------------------------------------
  combined <- (p_a | p_b) / (p_c | p_d) / (p_e + plot_spacer()) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 6, face = "plain"))
  ggsave(file.path(OUT_DIR, paste0("figS_lfc_sensitivity", sfx, ".pdf")), combined,
         width = 12, height = 12, device = cairo_pdf)

  # ----- Companion data ----------------------------------------------------
  companion <- merge(panelA[, .(lfc, n_deg, n_overlap, fisher_OR = OR, fisher_p = p)],
                     panelB[, .(lfc, n_human_deg = n_human, n_mouse_deg = n_mouse,
                                n_conserved)], by = "lfc")
  companion <- merge(companion, panelC, by = "lfc")
  companion <- merge(companion, panelE, by = "lfc")
  fwrite(companion, file.path(OUT_DIR, paste0("sensitivity_data", sfx, ".csv")))

  if (length(top20_by_lfc) > 0) {
    fwrite(jacc, file.path(OUT_DIR, paste0("pathway_jaccard", sfx, ".csv")))
    top_long <- rbindlist(lapply(names(top20_by_lfc), function(n) {
      data.table(lfc = as.numeric(n), rank = seq_along(top20_by_lfc[[n]]),
                 pathway = top20_by_lfc[[n]])
    }))
    fwrite(top_long, file.path(OUT_DIR, paste0("top20_pathways_per_cutoff", sfx, ".csv")))
  }

  message("Caption: LFC-cutoff sensitivity on the ", tlfc, " scale (sig column ",
          pcol, " < 0.05, effect column ", lcol, "). Downstream convergence ",
          "(COLOC overlap, cross-species conservation, triple convergence, ",
          "Hallmark pathway Jaccard, multi-evidence count) across ",
          "|log2FC| in {0.2, 0.3, 0.5}; canonical cutoff 0.5 highlighted.")
  cat("Saved scale", sfx, "to:", OUT_DIR, "\n")
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Run both scales
# -----------------------------------------------------------------------------
for (nm in names(SCALES)) build_scale(SCALES[[nm]])

cat("\nDone (raw + shrunk).\n")
