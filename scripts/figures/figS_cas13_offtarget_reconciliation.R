#!/usr/bin/env Rscript
# =============================================================================
# Fig S_cas13_offtarget_reconciliation — Why the simple bowtie QC "off-targets"
#   are NOT real off-targets (vs the columba+oahmed production pipeline).
# KEY MESSAGE: every apparent QC flag resolves to a non-off-target category;
#   0 real off-targets remain — reproducing columba's clean call.
# Data: A8 REAL bowtie1 -v run on the SAME gencode.vM38.transcripts.filt.fa
#   columba used (job 18493837), reconciled against columba offtarget_calls.tsv.
#   Source CSVs written by scratchpad/prep_offtarget_fig_data.py.
# Panels (two individual PDFs; assemble in Illustrator):
#   A  Filter cascade: 3,695 raw flags -> 0 real off-targets (by off-target gene-pair)
#   B  Named survivors at <=2 mm: all are gene-family paralogs or columba
#      indel-filtered near-cognates (none real)
# Output: Cas13_Library_Design/figures/guide_qc/offtarget_reconciliation_{A,B}_*.pdf
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
QC_DIR <- file.path(FIGS_CAS13LIB_DIR, "guide_qc")

casc <- fread(file.path(QC_DIR, "offtarget_reconciliation_cascade_source.csv"))
pairs <- fread(file.path(QC_DIR, "offtarget_reconciliation_named_source.csv"))

# ---------------------------------------------------------------------------
# Panel A — filter cascade (horizontal bars, sqrt length axis, counts on bars)
# ---------------------------------------------------------------------------
cat_cols <- c(raw = "#455A64", threshold = "#4baeef", control = "#9E9E9E",
              indel = "#00b3a4", paralog = "#C9265E")
casc[, step_label := factor(step_label, levels = rev(step_label))]  # raw at top
casc[, end_lab := ifelse(category == "paralog", "0 real off-targets",
                         formatC(remaining, format = "d", big.mark = ","))]

pA <- ggplot(casc, aes(remaining, step_label, fill = category)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = end_lab), hjust = -0.08, size = GEOM_TEXT_6PT, colour = "black") +
  scale_fill_manual(values = cat_cols, guide = "none") +
  scale_x_sqrt(breaks = c(0, 100, 1000, 3695),
               labels = scales::label_comma(),
               expand = expansion(mult = c(0, 0.42))) +
  labs(x = "Off-target gene-pairs remaining (√ scale)", y = NULL) +
  theme_masld_compact() +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
        axis.text.y = element_text(size = 6, colour = "black", hjust = 0))
save_fig(pA, file.path(QC_DIR, "offtarget_reconciliation_A_cascade.pdf"),
         width = 6.4, height = 2.3)

# ---------------------------------------------------------------------------
# Panel B — named survivors at <=2 mm, faceted by class, coloured by biotype
# ---------------------------------------------------------------------------
pairs[, class := factor(class, levels = c("Gene-family paralog",
                                          "Near-cognate (columba indel-filtered)"))]
pairs[, biotype := factor(biotype, levels = c("protein-coding", "lncRNA", "pseudogene"))]
# order rows within each facet by mismatch then label
setorder(pairs, class, min_mismatch, pair_label)
pairs[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
bt_cols <- c("protein-coding" = "#1565C0", "lncRNA" = "#00b3a4", "pseudogene" = "#9E9E9E")

pB <- ggplot(pairs, aes(min_mismatch, pair_label, colour = biotype)) +
  geom_point(size = 1.5) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = bt_cols, name = "Off-target biotype") +
  scale_x_continuous(breaks = 0:2, limits = c(-0.3, 2.3)) +
  labs(x = "Mismatches to off-target (bowtie −v, ungapped)", y = NULL) +
  theme_masld_compact() +
  theme(axis.text.y = element_text(size = 6, colour = "black", face = "italic"),
        panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey92"),
        legend.position = "top", strip.text.y = element_text(angle = 90))
save_fig(pB, file.path(QC_DIR, "offtarget_reconciliation_B_named.pdf"),
         width = 3.9, height = 6.4)

message("── CAPTIONS (add in layout; not baked into panels) ─────────────────")
message("Panel A. Off-target QC on the v9 library (bowtie1 -v3, target_seq, mouse ",
        "mature transcriptome gencode.vM38.filt) flags 3,695 guide→gene pairs. Each ",
        "is resolved against the columba+oahmed pipeline: 3,629 are ≥3-mismatch ",
        "(20/23 nt) hits below the Cas13 activity floor and columba's edit≤2 cap; 18 ",
        "are non-targeting/essential control guides absent from the screen roster; 33 ",
        "are 2-edit sites columba represents as a mid-seed indel and drops via its ",
        "CIGAR filter (bowtie shows the ungapped 2-substitution form of the same site ",
        "at the same position); 15 are gene-family paralogs (≥98% identity). 0 real ",
        "off-targets remain, reproducing columba's independent clean call.")
message("Panel B. The 39 distinct ≤2-mm off-target pairs (columba's 8,654-guide input, ",
        "apples-to-apples). All 0-mm perfect matches are same-family paralogs of the ",
        "guide's own target (e.g. Hba-a2→Hba-a1, Dynlt1c→Dynlt1f), which no 23-mer can ",
        "spare; the remainder are near-cognates columba filtered as mid-seed indels. ",
        "None is a novel off-target to an unrelated gene.")
message("Wrote: ", file.path(QC_DIR, "offtarget_reconciliation_A_cascade.pdf"))
message("Wrote: ", file.path(QC_DIR, "offtarget_reconciliation_B_named.pdf"))
