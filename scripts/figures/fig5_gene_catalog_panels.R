#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 5 — MASLD Gene Catalog panels, STEP 2 of 2 (rendering).
#
# KEY MESSAGE: what an absence MEANS depends on whether the gene was testable.
# Two panels, no score, no rank, no leaderboard, no default ordering by count.
#
#   fig6_gene_catalog_state.pdf
#       Joint assay-call matrix (genetics call x transcriptomics call) over the
#       whole accepted gene universe. Fill encodes only whether an inferential
#       result EXISTS (both layers testable) or not. The `tested_negative` row
#       and column are drawn EMPTY, at n = 0, so the reader sees that this
#       release contains no informative negative anywhere.
#
#   fig6_gene_catalog_next_experiment.pdf
#       The one deterministic discriminating experiment each class routes to,
#       with the outcome that would falsify it. Bars are gene counts, ordered by
#       evidence-class semantics, NOT by magnitude.
#
# Inputs are written by scripts/figures/fig5_gene_catalog_panels.py from the sealed
# PASS06 candidate. No scientific call is recomputed here.
#
# Env: micromamba activate rnaseq
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

DATA_DIR  <- file.path(BASE, "figures/main/fig6_gene_catalog/data")
PANEL_DIR <- file.path(BASE, "figures/main/fig6_gene_catalog/panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

PT <- 6 / ggplot2::.pt   # geom_text size for 6 pt

# Neutral, non-ranked fills. Nothing here reads as "better".
FILL_RESULT   <- "#E8EDF3"   # both layers testable -> an inferential result exists
FILL_NOMEAS   <- "#EDEDED"   # >=1 layer untestable -> no inferential result
FILL_ABSENT   <- "#FFFFFF"   # informative-negative band, empty in this release
BAR_FILL      <- "#C7D3DF"

CALL_LEVELS <- c("supported", "indeterminate", "tested_negative", "untestable")
CALL_LABS   <- c(supported = "supported", indeterminate = "indeterminate",
                 tested_negative = "tested negative", untestable = "untestable")

# ── Panel 1 ─────────────────────────────────────────────────────────────────
m <- read.csv(file.path(DATA_DIR, "fig5_gene_catalog_call_matrix.csv"),
              stringsAsFactors = FALSE)
prov <- read.csv(file.path(DATA_DIR, "fig5_gene_catalog_provenance.csv"),
                 stringsAsFactors = FALSE)
marg <- read.csv(file.path(DATA_DIR, "fig5_gene_catalog_domain_margin.csv"),
                 stringsAsFactors = FALSE)

n_total <- sum(m$n_genes)
stopifnot(n_total > 0)

m$genetics_call        <- factor(m$genetics_call, levels = rev(CALL_LEVELS))
m$transcriptomics_call <- factor(m$transcriptomics_call, levels = CALL_LEVELS)
m$cell_fill <- ifelse(m$cell_kind == "informative_negative_absent", FILL_ABSENT,
               ifelse(m$cell_kind == "not_measurable", FILL_NOMEAS, FILL_RESULT))
m$cell_lab  <- format(m$n_genes, big.mark = ",", trim = TRUE)

# Class name printed under the count only where the frozen index gives one
# unambiguous class for the cell.
wrap_lab <- function(x, w) vapply(x, function(s)
  if (!nzchar(s)) "" else paste(strwrap(s, width = w), collapse = "\n"),
  character(1), USE.NAMES = FALSE)

m$class_lab <- ifelse(m$n_genes > 0 & !grepl(";", m$primary_evidence_class) &
                        m$primary_evidence_class != "",
                      gsub("_", " ", m$primary_evidence_class), "")
m$class_lab <- wrap_lab(m$class_lab, 12)

p1 <- ggplot(m, aes(transcriptomics_call, genetics_call)) +
  geom_tile(aes(fill = cell_fill), colour = "black", linewidth = 0.3) +
  geom_tile(data = subset(m, cell_kind == "informative_negative_absent"),
            fill = NA, colour = "black", linewidth = 0.3, linetype = "22") +
  scale_fill_identity() +
  geom_text(aes(label = cell_lab), size = PT, colour = "black",
            vjust = ifelse(m$class_lab == "", 0.5, -0.6)) +
  geom_text(aes(label = class_lab), size = PT, colour = "black",
            vjust = 0.9, lineheight = 1.05) +
  scale_x_discrete(labels = CALL_LABS, position = "top", expand = c(0, 0)) +
  scale_y_discrete(labels = CALL_LABS, expand = c(0, 0)) +
  labs(x = "Established disease-state transcriptome (assay call)",
       y = "Liver-trait regulatory genetics (assay call)") +
  coord_equal() +
  theme_masld() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 0, vjust = 0),
        plot.margin = margin(3, 3, 3, 3))

save_fig(p1, file.path(PANEL_DIR, "fig6_gene_catalog_state.pdf"),
         width = 4.2, height = 3.6)

# CAPTION via message() — never a plot subtitle.
gen_unt <- subset(marg, evidence_domain == "genetics" & call_state == "untestable")
txn_unt <- subset(marg, evidence_domain == "transcriptomics" & call_state == "untestable")
one_sided <- sum(subset(m, cell_kind == "not_measurable" &
                          (genetics_call == "supported" |
                             transcriptomics_call == "supported"))$n_genes)
message(sprintf(
  paste0("[fig5_gene_catalog_state] n = %s genes, %s gene-evidence rows. ",
         "Cells are gene counts; fill marks only whether both layers were testable. ",
         "Genetics untestable n = %s (reason: %s); transcriptome untestable n = %s ",
         "(reason: %s). %s genes carry a supported call in one layer while the other ",
         "layer could not be evaluated. The informative-negative band is empty: ",
         "%s rows in the release meet the frozen tested-negative rule, so no cell ",
         "in this matrix is evidence of absence. Provenance: %s of %s gene-evidence ",
         "rows are '%s'; none is independent."),
  format(n_total, big.mark = ","), format(sum(prov$n_rows), big.mark = ","),
  format(sum(gen_unt$n_genes), big.mark = ","), paste(unique(gen_unt$testability_reason), collapse = "/"),
  format(sum(txn_unt$n_genes), big.mark = ","), paste(unique(txn_unt$testability_reason), collapse = "/"),
  format(one_sided, big.mark = ","),
  sum(m$n_genes[m$cell_kind == "informative_negative_absent"]),
  format(sum(prov$n_rows), big.mark = ","), format(sum(prov$n_rows), big.mark = ","),
  paste(unique(prov$provenance_state), collapse = "/")))

# ── Panel 2 ─────────────────────────────────────────────────────────────────
nx <- read.csv(file.path(DATA_DIR, "fig5_gene_catalog_next_experiment.csv"),
               stringsAsFactors = FALSE)
# Descriptive frequency order only. No merit hierarchy is implied or encoded.
nx <- nx[order(nx$display_order_by_count), ]

# Row label = the evidence class that routes here + the perturbation it asks for.
nx$row_lab <- sprintf("%s\n%s", gsub("_", " ", nx$primary_evidence_class),
                      wrap_lab(nx$perturbation, 34))
nx$row_lab <- factor(nx$row_lab, levels = rev(nx$row_lab))
nx$count_lab <- format(nx$n_genes, big.mark = ",", trim = TRUE)

p2 <- ggplot(nx, aes(n_genes, row_lab)) +
  geom_col(fill = BAR_FILL, colour = "black", linewidth = 0.3, width = 0.62) +
  geom_text(aes(label = count_lab), hjust = -0.15, size = PT, colour = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.34)),
                     labels = scales::comma) +
  labs(x = "Genes routed (count, not priority)", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(hjust = 0, lineheight = 1.15),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_fig(p2, file.path(PANEL_DIR, "fig6_gene_catalog_next_experiment.pdf"),
         width = 5.0, height = 3.0)

message(sprintf(
  paste0("[fig5_gene_catalog_next_experiment] One deterministic experiment per gene ",
         "from the frozen rulebook; %s genes routed across %d rules, ordered by ",
         "evidence-class semantics, not by count. Falsifiers: %s. Each proposal ",
         "is a discriminating research design, NOT a calibrated probability of ",
         "therapeutic success and not a clinical recommendation."),
  format(sum(nx$n_genes), big.mark = ","), nrow(nx),
  paste(sprintf("%s -> %s", gsub("_", " ", nx$primary_evidence_class),
                nx$falsifying_outcome), collapse = "; ")))

message("[done] wrote 2 panels to ", PANEL_DIR)
