#!/usr/bin/env Rscript
# crossspecies_conservation.R — Fig 4g (cross-species conserved core)
#
# Message: a cross-species conserved core (1,355 genes concordant in >=3 of 4
# mouse models) is a translatability filter enriched for druggable targets
# (61 drug targets; OR 4.92 vs a biotype-matched null).
#
# DESIGN (2026-06-19 redesign): LEAN composition only — a single compact
# conserved-core composition element (3-of-4 vs all-4 models) with hero genes
# named and the druggable count annotated. The drug-target ENRICHMENT (OR 4.92)
# and the per-gene drug annotations are NOT plotted as a one-dot panel; they are
# written to a SUPPLEMENTARY TABLE and cited. (No lollipop, no single-data-point
# panel — see memory/feedback-no-lollipop.)
#
# Output: figures/main/fig4_validation/fig4g_crossspecies_conservation.pdf
#         + supp table figures/main/fig4_validation/_supp/data/conserved_core_drug_targets.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Conserved core composition ───────────────────────────────────────────────
conc <- fread(file.path(CONCORDANCE, "concordance_atlas_unified.csv"))
core <- conc[n_concordant >= 3]
n_core  <- nrow(core); n_three <- sum(core$n_concordant == 3); n_four <- sum(core$n_concordant >= 4)
stopifnot(n_core == 1355L)

# ── Drug-target annotation (atlas) for the supp table + hero-gene picking ─────
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "dgidb_druggable", "opentargets_drug", "drug_dev_status",
             "clintrial_n_drugs", "clintrial_drugs"))
core_a <- merge(core[, .(gene = human_symbol, n_models = n_concordant,
                         human_logFC = mean_h_lfc, translatability = translatability_score)],
                atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
# A real drug target = actual drug-development evidence (MASLD-specific stage,
# drugged for another indication, or an active clinical trial) — NOT the broad
# "druggable genome" (dgidb) or "discovery/undetermined" status.
core_a[, ctn := suppressWarnings(as.integer(clintrial_n_drugs))]
core_a[, is_drug_target := drug_dev_status %in% c("masld_approved", "masld_clinical",
        "masld_preclinical", "masld_discontinued", "drugged_other_indi") |
        (!is.na(ctn) & ctn > 0)]
n_drug <- core_a[is_drug_target == TRUE, .N]

# Canonical enrichment number (H6 positive-control): 61 drug targets, OR 4.92.
h6 <- fread(file.path(AUDIT, "H6_inv_cc_positive_control.csv"))
hr <- h6[comparator == "drug_target" & test == "B_biotypematched" & universe == "full_atlas"]
obs_or <- hr$observed_OR[1]; n_dt_canon <- 61L

# Hero genes: recognizable druggable genes, prefer all-4-model depth, then |logFC|.
priority <- c("THRB","FGF21","RORA","PPARA","DGAT2","SCD","ACACA","LPL","SPP1","HMGCR","NR1H4","CIDEC")
heroes4 <- core_a[n_models >= 4 & (is_drug_target | gene %in% priority)][
  order(match(gene, priority), -abs(human_logFC))][, head(gene, 6)]
heroes3 <- core_a[n_models == 3 & (is_drug_target | gene %in% priority) & !(gene %in% heroes4)][
  order(match(gene, priority), -abs(human_logFC))][, head(gene, 6)]

# ── Lean composition: ONE horizontal segmented bar (geom_rect) + hero callouts ─
# Compact (Liang style): black header, counts inside bar, hero genes right under it.
YB <- 0.80; YT <- 1.14   # bar vertical band
seg <- data.table(
  lab  = c("3 of 4 models", "all 4 models"),
  xmin = c(0, n_three), xmax = c(n_three, n_core),
  n    = c(n_three, n_four),
  fill = c("#9CC9C2", masld_colors$conserved),
  tcol = c("grey15", "white"),
  heroes = c(paste(heroes3, collapse = ", "), paste(heroes4, collapse = ", ")))
seg[, xmid := (xmin + xmax) / 2]

p <- ggplot(seg) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = YB, ymax = YT, fill = fill),
            color = "white", linewidth = 0.5) +
  geom_text(aes(x = xmid, y = (YB + YT) / 2,
                label = sprintf("%s\n%s genes", lab, formatC(n, big.mark = ",", format = "d")),
                color = tcol), size = PUB_GEOM_TEXT, lineheight = 0.9) +
  # hero-gene callouts just beneath each segment
  geom_text(aes(x = xmid, y = YB - 0.04, label = heroes), vjust = 1, size = PUB_GEOM_TEXT - 0.3,
            fontface = "italic", color = "grey30") +
  scale_fill_identity() + scale_color_identity() +
  scale_x_continuous(limits = c(0, n_core * 1.01),
                     breaks = c(0, 500, 1000, 1355),
                     labels = scales::comma, expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(0.55, 1.22), expand = c(0, 0)) +
  labs(x = "number of conserved-core genes", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text = element_text(color = "black"),
        plot.margin = margin(2, 6, 2, 3))

out_pdf <- file.path(FIG4_DIR, "fig4g_crossspecies_conservation.pdf")
save_fig(p, out_pdf, width = fig_col_width + 0.6, height = 1.35)
message("Saved: ", out_pdf)

# ── Supplementary table (the enrichment + per-gene drug annotations) ─────────
supp_dir <- file.path(FIG4_DIR, "_supp", "data")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
supp <- core_a[order(-n_models, -is_drug_target, -abs(human_logFC)),
               .(gene, n_models_concordant = n_models,
                 human_logFC = round(human_logFC, 3),
                 translatability = round(translatability, 3),
                 is_drug_target,
                 opentargets_drug, drug_dev_status, clintrial_drugs)]
supp_path <- file.path(supp_dir, "conserved_core_drug_targets.csv")
# header line documents the canonical enrichment
writeLines(sprintf("# Conserved core (>=3/4 mouse models): %d genes; %d drug targets (canonical, H6); biotype-matched-null OR %.2f (emp_p<0.001). Per-gene atlas drug annotations below.",
                   n_core, n_dt_canon, obs_or), supp_path)
fwrite(supp, supp_path, append = TRUE, col.names = TRUE)
message("Supp table: ", supp_path)
message(sprintf("[crossspecies] core=%d (3of4=%d, all4=%d); atlas-druggable=%d (canonical 61); OR=%.2f",
                n_core, n_three, n_four, n_drug, obs_or))
message("[crossspecies] heroes all-4: ", paste(heroes4, collapse=", "),
        " | 3-of-4: ", paste(heroes3, collapse=", "))
