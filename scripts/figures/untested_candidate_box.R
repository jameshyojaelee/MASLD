#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# KEY MESSAGE: The unsupervised convergence score nominates a set of genetically
# convergent targets that are NOT YET DRUGGED — the prospective end of the
# drug-development gradient. After approved (THRB/GLP1R/SLC5A2), clinical
# (RORA, Phase 1 TB-840) and preclinically validated (HKDC1, AKR1B10) anchors
# validate the prioritizer, these 9 candidates are the forward-looking
# nominations: high convergence rank, drug_dev_status == "discovery"
# (no drug / probe / MASLD-liver perturbation in adversarial PubMed mining).
#
# This is the PROSPECTIVE Cas13 target list — drug-development-gradient
# validation framing, NOT a novel-discovery claim.
#
# DESIGN: a compact horizontal lollipop / table, ranked by convergence_rank
# (best at top). Each gene shows its convergence rank (lollipop position +
# label), concordance state (point colour), and a verified "not-yet-drugged"
# status badge. Minimal/clean, no 3D.
#
# DATA (verified, traced):
#   * convergence_rank + concordance_state <- RNA-seq/results/multi_evidence/
#       convergence_evidence.csv (column human_symbol)
#   * drug_dev_status (== "discovery" for all 9) <- data/external/drug_targets/
#       drug_target_classification.tsv (column symbol)
# RE-PULLED 2026-06-19 after the INTACT->COLOC genetic-gate swap re-scored the
# convergence ranking. Verified-clean candidates now in the convergence top:
#   GCAT(2), TMEM184B(5), RAPH1(8), EPB41L4B(9), COL25A1(14), OCEL1(17),
#   SDC2(18; DGIdb-druggable), RHOBTB3(19; Tier-1), IKZF5(22).
#   OCEL1/RHOBTB3/IKZF5 are NEW (0 liver-disease PubMed hits). PCOLCE2 (now #325)
#   and EFHD1 (#134) remain verified-untested but fell out of the top under the
#   COLOC re-score, so they no longer headline.
# SCREENED OUT on the same PubMed deep-mine (substantial prior literature):
#   FOXN3 (now #1; hepatocyte-KO alleviates NAFLD), CIITA(10), ESRP2(15),
#   SH3YL1(188) — hidden MASH perturbation; plus ARHGEF39(13; cancer oncogene),
#   ADAMTS6(23; NASH expression panel), TUSC3(24; HCC-TS + XMEN therapy).
#
# Output: figures/main/fig5_convergence/panels/untested_candidate_box.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── 1. The 9 verified genuinely-untested candidates (re-pulled 2026-06-19 ──────
#       after the COLOC gate swap; all top convergence + PubMed-clean) ─────────
cand_genes <- c("GCAT", "TMEM184B", "RAPH1", "EPB41L4B",
                "COL25A1", "OCEL1", "SDC2", "RHOBTB3", "IKZF5")

# ── 2. Pull convergence rank + concordance state ──────────────────────────────
conv <- read.csv(file.path(BASE,
  "RNA-seq/results/multi_evidence/convergence_evidence.csv"),
  stringsAsFactors = FALSE)
conv_sub <- conv[match(cand_genes, conv$human_symbol),
                 c("human_symbol", "convergence_rank", "concordance_state")]
stopifnot(!any(is.na(conv_sub$convergence_rank)))

# ── 3. Pull drug_dev_status (must be "discovery" = not-yet-drugged) ───────────
drug <- read.delim(file.path(BASE,
  "data/external/drug_targets/drug_target_classification.tsv"),
  stringsAsFactors = FALSE)
drug_sub <- drug[match(cand_genes, drug$symbol),
                 c("symbol", "drug_dev_status", "pharos_tdl", "dgidb_n_drugs")]
# Sanity: the spec defines these as genuinely untested (discovery)
if (!all(drug_sub$drug_dev_status == "discovery")) {
  warning("Some candidate genes are no longer drug_dev_status == 'discovery': ",
          paste(drug_sub$symbol[drug_sub$drug_dev_status != "discovery"],
                collapse = ", "))
}

df <- merge(conv_sub, drug_sub,
            by.x = "human_symbol", by.y = "symbol", sort = FALSE)
df <- df[order(df$convergence_rank), ]

# Display direction from concordance_state (Concordant-up vs everything-down)
df$direction <- ifelse(grepl("up", df$concordance_state, ignore.case = TRUE),
                       "Up in MASLD", "Down in MASLD")
# DGIdb-druggable footnote flag (SDC2 has a known small molecule)
df$dgidb_n_drugs[is.na(df$dgidb_n_drugs)] <- 0
df$druggable_flag <- ifelse(df$dgidb_n_drugs > 0, " *", "")

# Order factor so best rank sits at the TOP of the lollipop
df$gene_lab <- paste0(df$human_symbol, df$druggable_flag)
df$gene_lab <- factor(df$gene_lab, levels = rev(df$gene_lab))

# ── 4. Colours: concordance direction (reuse semantic disease/control hues) ───
dir_cols <- c("Up in MASLD"   = unname(masld_colors$up),    # Liang deep magenta
              "Down in MASLD" = unname(masld_colors$down))  # deep blue

# ── 5. Plot — horizontal lollipop ─────────────────────────────────────────────
xmax <- max(df$convergence_rank) * 1.18

p <- ggplot(df, aes(x = convergence_rank, y = gene_lab)) +
  geom_segment(aes(x = 0, xend = convergence_rank,
                   y = gene_lab, yend = gene_lab,
                   colour = direction),
               linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(colour = direction), size = 2.4) +
  geom_text(aes(label = convergence_rank),
            colour = "white", size = PUB_GEOM_TEXT - 0.1, fontface = "bold") +
  scale_colour_manual(values = dir_cols, name = "Concordance") +
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0)),
                     name = "Convergence rank (lower = stronger)") +
  labs(
    title = "Genetically convergent, not-yet-drugged candidates",
    subtitle = "Prospective Cas13 targets — drug_dev_status = discovery\n(no drug / probe / MASLD-liver perturbation)",
    y = NULL,
    caption = paste0(
      "Rank + concordance: convergence_evidence.csv (COLOC-gated re-score).\n",
      "Drug-dev status: drug_target_classification.tsv (all 9 = \"discovery\").\n",
      "* = DGIdb small molecule exists (SDC2) but no MASLD program.\n",
      "Screened out by adversarial PubMed deep-mine (prior perturbation/\n",
      "characterisation): FOXN3 (top-ranked #1), CIITA, ESRP2, SH3YL1,\n",
      "ARHGEF39, ADAMTS6, TUSC3 — the verification check works.")
  ) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(0.25, "cm"),
    axis.text.y = element_text(size = PUB_AXIS_TITLE, face = "italic"),
    plot.caption = element_text(size = PUB_SUBTITLE - 1.5, colour = "gray35",
                                hjust = 0, lineheight = 1.05),
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )

# ── 6. Save ───────────────────────────────────────────────────────────────────
out <- file.path(FIG5_DIR, "panels", "untested_candidate_box.pdf")
save_fig(p, out, width = fig_half_width, height = 2.9)
cat("Wrote:", out, "\n")
cat("Genes (rank):",
    paste(df$human_symbol, df$convergence_rank, sep = "=", collapse = ", "), "\n")
