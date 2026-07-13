#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# KEY MESSAGE: After the convergence prioritiser is calibrated by drugged anchors
# (approved THRB/GLP1R/SLC5A2; clinical RORA; preclinical HKDC1/AKR1B10), it
# nominates a set of genetically-convergent MASLD targets that carry NO existing
# therapeutic program (drug_dev_status == "discovery": no approved/clinical drug,
# no chemical probe, no documented MASLD-liver perturbation in adversarial PubMed
# mining). These are the prioritised, MODALITY-AGNOSTIC therapeutic opportunities.
#
# Each candidate's Pharos Target Development Level (TDL) gives generalisable
# therapeutic guidance independent of any one platform:
#   Tbio  = biologically characterised but no chemical matter yet
#           (ligand-discovery / genetic-modality entry point)
#   Tdark = understudied "dark" target (highest novelty, highest risk)
# and the disease direction states the intervention logic (inhibit up-genes /
# restore down-genes).
#
# DESIGN (house style: no title/subtitle/in-plot text, no lollipop, black text,
# control = #9E9E9E, compact/dense): a horizontal CONVERGENCE-RANK BAR per gene
# (fill = disease direction) paired with a Pharos-TDL druggability strip. All
# narrative lives in the caption emitted to stdout, never on the panel.
#
# DATA (verified, traced):
#   * convergence_rank + concordance_state  <- RNA-seq/results/multi_evidence/
#       convergence_evidence.csv (column human_symbol)
#   * drug_dev_status (== "discovery") + pharos_tdl + dgidb_n_drugs
#       <- data/external/drug_targets/drug_target_classification.tsv (column symbol)
# RE-PULLED 2026-06-19 after the INTACT->COLOC genetic-gate swap re-scored the
# convergence ranking. Verified-clean undrugged candidates now in the top:
#   GCAT(2), TMEM184B(5), RAPH1(8), EPB41L4B(9), COL25A1(14), OCEL1(17),
#   SDC2(18; one DGIdb small molecule, no MASLD program), RHOBTB3(19), IKZF5(22).
#   OCEL1/RHOBTB3/IKZF5 are NEW (0 liver-disease PubMed hits). PCOLCE2 (#325) and
#   EFHD1 (#134) remain verified-untested but fell out of the top under the re-score.
# SCREENED OUT on the same PubMed deep-mine (substantial prior literature):
#   FOXN3 (now #1; hepatocyte-KO alleviates NAFLD), CIITA(10), ESRP2(15),
#   SH3YL1(188), ARHGEF39(13), ADAMTS6(23), TUSC3(24) — the verification works.
#
# Output: figures/main/fig5_convergence/panels/untested_candidate_box.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── 1. Verified genetically-convergent, therapeutically-unexploited candidates ─
cand_genes <- c("GCAT", "TMEM184B", "RAPH1", "EPB41L4B",
                "COL25A1", "OCEL1", "SDC2", "RHOBTB3", "IKZF5")

# ── 2. Convergence rank + concordance state ───────────────────────────────────
conv <- read.csv(file.path(BASE,
  "RNA-seq/results/multi_evidence/convergence_evidence.csv"),
  stringsAsFactors = FALSE)
conv_sub <- conv[match(cand_genes, conv$human_symbol),
                 c("human_symbol", "convergence_rank", "concordance_state")]
stopifnot(!any(is.na(conv_sub$convergence_rank)))

# ── 3. Therapeutic status: drug_dev_status + Pharos druggability tier ──────────
drug <- read.delim(file.path(BASE,
  "data/external/drug_targets/drug_target_classification.tsv"),
  stringsAsFactors = FALSE)
drug_sub <- drug[match(cand_genes, drug$symbol),
                 c("symbol", "drug_dev_status", "pharos_tdl", "dgidb_n_drugs")]
if (!all(drug_sub$drug_dev_status == "discovery")) {
  warning("Candidate no longer drug_dev_status == 'discovery': ",
          paste(drug_sub$symbol[drug_sub$drug_dev_status != "discovery"], collapse = ", "))
}

df <- merge(conv_sub, drug_sub, by.x = "human_symbol", by.y = "symbol", sort = FALSE)
df <- df[order(df$convergence_rank), ]
df$dgidb_n_drugs[is.na(df$dgidb_n_drugs)] <- 0

# Disease direction = intervention logic (restore down-genes / inhibit up-genes)
df$direction <- factor(
  ifelse(grepl("up", df$concordance_state, ignore.case = TRUE), "Up in MASLD", "Down in MASLD"),
  levels = c("Down in MASLD", "Up in MASLD"))

# Pharos druggability tier (modality-agnostic tractability), descriptive labels
tdl_levels <- c("Tdark", "Tbio", "Tchem", "Tclin")
tdl_labels <- c(Tdark = "Tdark (understudied)", Tbio = "Tbio (characterised)",
                Tchem = "Tchem (chemical matter)", Tclin = "Tclin (drugged)")
df$tdl <- factor(df$pharos_tdl, levels = tdl_levels)

# Best rank at the TOP
df$gene <- factor(df$human_symbol, levels = rev(df$human_symbol))

# ── 4. Colours ────────────────────────────────────────────────────────────────
dir_cols <- c("Down in MASLD" = unname(masld_colors$down),   # deep blue
              "Up in MASLD"   = unname(masld_colors$up))     # deep magenta
tdl_cols <- c("Tdark" = gray_gradient[10], "Tbio" = gray_gradient[5],
              "Tchem" = gray_gradient[3],  "Tclin" = gray_gradient[1])

# ── 5a. Druggability-tier strip + gene names (left) ───────────────────────────
p_tdl <- ggplot(df, aes(x = 1, y = gene)) +
  geom_tile(aes(fill = tdl), colour = "white", width = 0.96, height = 0.82) +
  scale_fill_manual(values = tdl_cols, breaks = tdl_levels, labels = tdl_labels,
                    name = "Druggability (Pharos TDL)", drop = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme(axis.text.y  = element_text(face = "italic", colour = "black", size = PUB_AXIS_TITLE),
        axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(), panel.border = element_blank(),
        plot.margin  = margin(2, 1, 2, 2))

# ── 5b. Convergence-rank bars, filled by disease direction (right) ────────────
xmax <- max(df$convergence_rank) * 1.10
p_bar <- ggplot(df, aes(x = convergence_rank, y = gene, fill = direction)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = convergence_rank), hjust = -0.35, colour = "black",
            size = PUB_GEOM_TEXT + 0.4) +
  scale_fill_manual(values = dir_cols, name = "Disease direction") +
  scale_x_continuous(name = "Convergence rank (lower = stronger)",
                     limits = c(0, xmax), expand = expansion(mult = c(0, 0))) +
  labs(y = NULL) +
  theme_masld(base_size = 7) +
  theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin  = margin(2, 2, 2, 1))

# ── 6. Compose + save ─────────────────────────────────────────────────────────
p <- p_tdl + p_bar +
  plot_layout(widths = c(0.9, 6), guides = "collect") &
  theme(legend.position = "top", legend.justification = "left",
        legend.direction = "vertical", legend.box = "vertical", legend.box.just = "left",
        legend.title = element_text(size = PUB_AXIS_TITLE, face = "plain"),
        legend.text = element_text(size = PUB_AXIS_TITLE - 0.5),
        legend.key.size = unit(0.22, "cm"), legend.spacing.y = unit(0.02, "cm"),
        legend.margin = margin(0, 0, 0, 0))

out <- file.path(FIG5_DIR, "panels", "untested_candidate_box.pdf")
save_fig(p, out, width = fig_half_width, height = 3.2)
cat("Wrote:", out, "\n")

# ── Caption (house rule: narrative to stdout, never on the panel) ─────────────
cat("\nCAPTION: Genetically-convergent MASLD targets with no existing therapeutic\n",
    "program (drug_dev_status = 'discovery': no approved/clinical drug, no chemical\n",
    "probe, no documented MASLD-liver perturbation). The calibrated convergence score\n",
    "(ranks from convergence_evidence.csv) nominates these as the prioritised,\n",
    "modality-agnostic therapeutic opportunities. Bar = convergence rank (lower =\n",
    "stronger); fill = disease direction (restore down-regulated / inhibit up-regulated);\n",
    "left strip = Pharos Target Development Level (Tbio = biologically characterised but\n",
    "chemically unexploited; Tdark = understudied). SDC2 carries one DGIdb small molecule\n",
    "but no MASLD program. Genes with prior MASLD literature (e.g. FOXN3, ESRP2) were\n",
    "removed by adversarial PubMed mining.\n", sep = "")
cat("\nGenes (rank | direction | TDL):\n")
for (i in seq_len(nrow(df)))
  cat(sprintf("  %-9s %2d  %-14s %s\n", df$human_symbol[i], df$convergence_rank[i],
              as.character(df$direction[i]), df$pharos_tdl[i]))
