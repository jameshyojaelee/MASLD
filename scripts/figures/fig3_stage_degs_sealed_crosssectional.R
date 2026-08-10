#!/usr/bin/env Rscript
# ============================================================================
# fig3_stage_degs_sealed_crosssectional.R
#
# Fig 3 — CROSS-SECTIONAL adjacent-stage differential expression.
#
# NEW script (2026-08-08). Does NOT edit or replace any existing figure script.
# Supersedes, for manuscript use, the retired-source panels `fib_stage_degs.pdf`
# and `fib_stage_upset.pdf`.
#
# SINGLE ADMISSIBLE SOURCE — the sealed candidate bundle
#   .../candidates/program-context-v2-candidate-2026-08-07/
#       fibrosis-adjacent-true-kleiner-lvqw-v1/
# whose READY manifest SHA256 is verified in-script before anything is read.
#
# Model (from the bundle's own design audit, re-asserted below):
#   ~ dataset + inferred_sex + fib_group      coefficient `fib_grouphigh`
#   limma-voom with quality weights; BH padj < 0.05; NO logFC filter.
#   Six exact-Kleiner cohorts; one pass-QC first biopsy per GSE193066 donor;
#   GSE213621 excluded (its F0F1/F3F4 labels are BINNED, not exact stages).
#
# SEMANTICS — these are CROSS-SECTIONAL contrasts between donors at adjacent
# Kleiner stages. They are NOT within-donor change and NOT a time course.
# Permitted wording: "F1 versus F0", "stage-associated", "stage-ordered
# remodeling". The words transition / progression / cascade / therapeutic
# window are forbidden and appear nowhere in this script's output.
#
# Every headline number is RE-DERIVED from the raw per-gene results here and
# asserted with stopifnot(). Any drift fails the job loudly rather than
# rendering a wrong panel.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# FIG2_DIR is the back-compat constant that points at figures/main/fig3_RNAseq
PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 0. Sealed bundle — locate and verify BEFORE reading any result
# ---------------------------------------------------------------------------
BUNDLE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates",
  "program-context-v2-candidate-2026-08-07",
  "fibrosis-adjacent-true-kleiner-lvqw-v1")

READY_SHA256 <- "f6a6d274124749b5c8bbdbf204629c847a33b45bc9d2a570ba98352d895fbcc9"

f_ready    <- file.path(BUNDLE, "READY")
f_results  <- file.path(BUNDLE, "results/fibrosis_consecutive_lvqw_true_kleiner.csv")
f_design   <- file.path(BUNDLE, "audits/design_audit.tsv")
f_manifest <- file.path(BUNDLE, "audits/sample_manifest.tsv")

for (f in c(f_ready, f_results, f_design, f_manifest)) {
  if (!file.exists(f)) stop("SEALED BUNDLE FILE MISSING: ", f)
}

observed_sha <- sub("\\s.*$", "", system2("sha256sum", shQuote(f_ready), stdout = TRUE))
if (!identical(observed_sha, READY_SHA256)) {
  stop("SEALED BUNDLE SHA256 MISMATCH.\n  expected: ", READY_SHA256,
       "\n  observed: ", observed_sha,
       "\nRefusing to render from an unverified source.")
}
message("[verify] READY SHA256 OK: ", observed_sha)

# ---------------------------------------------------------------------------
# 1. Re-derive every count from the raw per-gene results
# ---------------------------------------------------------------------------
res <- fread(f_results)

FAMILIES <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
# Display labels state the comparison, never a direction of travel.
LABELS   <- c(F0_to_F1 = "F1 vs F0", F1_to_F2 = "F2 vs F1",
              F2_to_F3 = "F3 vs F2", F3_to_F4 = "F4 vs F3")

PADJ_THR <- 0.05        # BH; NO logFC filter, by the bundle's own contract

counts <- res[, .(
  n_tested  = .N,
  n_deg     = sum(padj < PADJ_THR, na.rm = TRUE),
  n_up      = sum(padj < PADJ_THR & logFC > 0, na.rm = TRUE),
  n_down    = sum(padj < PADJ_THR & logFC < 0, na.rm = TRUE),
  n_donors  = unique(n_samples),
  n_cohorts = unique(n_cohorts),
  n_low     = unique(n_low),
  n_high    = unique(n_high)
), by = transition]

setkey(counts, transition)
counts <- counts[FAMILIES]
counts[, label := LABELS[transition]]

# ---------------------------------------------------------------------------
# 2. ASSERTIONS — fail loudly on ANY mismatch
# ---------------------------------------------------------------------------
EXPECT <- data.table(
  transition = FAMILIES,
  n_tested   = c(27638L, 27638L, 27638L, 27638L),
  n_deg      = c(  798L,  2190L,  3490L,   783L),
  n_up       = c(  442L,  1510L,  2131L,   504L),
  n_down     = c(  356L,   680L,  1359L,   279L),
  n_donors   = c(  313L,   361L,   307L,   152L),
  n_cohorts  = c(    6L,     6L,     6L,     5L)
)

chk <- merge(counts, EXPECT, by = "transition", suffixes = c("", "_exp"))
for (v in c("n_tested", "n_deg", "n_up", "n_down", "n_donors", "n_cohorts")) {
  bad <- chk[get(v) != get(paste0(v, "_exp"))]
  if (nrow(bad)) {
    stop("ASSERTION FAILED on ", v, ":\n",
         paste(capture.output(print(
           bad[, c("transition", v, paste0(v, "_exp")), with = FALSE])),
           collapse = "\n"))
  }
}
# up + down must exhaust the DEG set (no logFC == 0 significant gene)
stopifnot(all(counts$n_up + counts$n_down == counts$n_deg))
message("[assert] DEG counts OK: 798 / 2,190 / 3,490 / 783")
message("[assert] up/down OK: 442-356 / 1510-680 / 2131-1359 / 504-279")
message("[assert] donor denominators OK: 313 / 361 / 307 / 152")

# --- design contract ---------------------------------------------------------
design <- fread(f_design)
stopifnot(all(design$formula == "~ dataset + inferred_sex + fib_group"))
stopifnot(all(design$coefficient == "fib_grouphigh"))
stopifnot(all(design$design_rank == design$n_design_columns))   # no aliasing
stopifnot(all(design$cohort_set == design$expected_cohort_set))
stopifnot(all(design$n_unique_donors == design$n_samples))      # 1 sample = 1 donor
stopifnot(!any(grepl("GSE213621", design$cohort_set)))
message("[assert] design contract OK (formula, coefficient, full rank, ",
        "1 sample = 1 donor, GSE213621 absent)")

# --- donor census ------------------------------------------------------------
man <- fread(f_manifest)
stopifnot(all(man$pass_technical))                     # every row passed QC
stopifnot(all(man$biopsy %in% c("single deposited biopsy", "1st biopsy")))
stopifnot(!any(man$dataset == "GSE213621"))

census <- man[, .(n_donors = uniqueN(donor_id)), by = fibrosis_stage][
  order(fibrosis_stage)]
stopifnot(identical(census$n_donors, c(126L, 187L, 174L, 133L, 44L)))
stopifnot(uniqueN(man$donor_id) == 664L)               # no donor at two stages
stopifnot(sum(census$n_donors) == 664L)
message("[assert] donor census OK: F0-F4 = 126 / 187 / 174 / 133 / 44, ",
        "664 unique donors, no donor contributes two stages")

# NOTE the one asymmetry worth stating in the legend: F3 has 133 donors overall,
# but only 108 of them sit in the F4-vs-F3 comparison because GSE193066 has no
# F4 arm and is dropped from that comparison.
f3_in_f4vsf3 <- man[fibrosis_stage == 3L & transition == "F3_to_F4",
                    uniqueN(donor_id)]
stopifnot(f3_in_f4vsf3 == 108L)

# ---------------------------------------------------------------------------
# 3. Panel 1 — cross-sectional adjacent-stage DEG counts (diverging bars)
# ---------------------------------------------------------------------------
BASE_SIZE <- 6
LBL_SIZE  <- 6 / ggplot2::.pt
col_up    <- masld_colors$up      # "#C9265E"
col_down  <- masld_colors$down    # "#1565C0"

plot_dt <- rbindlist(list(
  counts[, .(label, direction = "Higher at the later stage",
             n = n_up,   y =  n_up)],
  counts[, .(label, direction = "Lower at the later stage",
             n = n_down, y = -n_down)]
))
plot_dt[, label := factor(label, levels = LABELS[FAMILIES])]
plot_dt[, direction := factor(direction,
  levels = c("Higher at the later stage", "Lower at the later stage"))]

# donor denominator annotation, placed under the axis
den <- counts[, .(label = factor(label, levels = LABELS[FAMILIES]),
                  txt = sprintf("n = %d donors", n_donors))]

ymax <- max(counts$n_up); ymin <- -max(counts$n_down)
pad  <- 0.16 * (ymax - ymin)

p1 <- ggplot(plot_dt, aes(x = label, y = y, fill = direction)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "black") +
  geom_text(aes(label = scales::comma(n),
                vjust = ifelse(y >= 0, -0.45, 1.35)),
            size = LBL_SIZE, colour = "black") +
  geom_text(data = den, inherit.aes = FALSE,
            aes(x = label, y = ymin - pad * 0.62, label = txt),
            size = LBL_SIZE, colour = "black") +
  scale_fill_manual(values = c("Higher at the later stage" = col_up,
                               "Lower at the later stage"  = col_down),
                    name = NULL) +
  scale_y_continuous(labels = function(v) scales::comma(abs(v)),
                     expand = expansion(mult = c(0.14, 0.10))) +
  labs(x = "Cross-sectional comparison of donors at adjacent Kleiner stages",
       y = "Stage-associated genes (BH FDR < 0.05)") +
  theme_masld_compact(base_size = BASE_SIZE) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank())

pdf(file.path(PANEL_DIR, "fig3_stage_degs_crosssectional.pdf"),
    width = 3.4, height = 2.6, useDingbats = FALSE)
print(p1)
invisible(dev.off())

# ---------------------------------------------------------------------------
# 4. Panel 2 — donor denominators behind each comparison
# ---------------------------------------------------------------------------
arm_dt <- rbindlist(list(
  counts[, .(label, arm = "Earlier stage", n = n_low)],
  counts[, .(label, arm = "Later stage",   n = n_high)]
))
arm_dt[, label := factor(label, levels = LABELS[FAMILIES])]
arm_dt[, arm := factor(arm, levels = c("Earlier stage", "Later stage"))]

p2 <- ggplot(arm_dt, aes(x = label, y = n, fill = arm)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.62) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.68),
            vjust = -0.45, size = LBL_SIZE, colour = "black") +
  scale_fill_manual(values = c("Earlier stage" = "#9E9E9E",
                               "Later stage"   = col_up),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = "Cross-sectional comparison of donors at adjacent Kleiner stages",
       y = "Donors") +
  theme_masld_compact(base_size = BASE_SIZE) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank())

pdf(file.path(PANEL_DIR, "fig3_stage_donor_census_crosssectional.pdf"),
    width = 3.4, height = 2.4, useDingbats = FALSE)
print(p2)
invisible(dev.off())

# ---------------------------------------------------------------------------
# 5. Sidecar data + captions
# ---------------------------------------------------------------------------
out <- counts[, .(comparison = label, n_genes_tested = n_tested,
                  n_stage_associated = n_deg, n_higher = n_up, n_lower = n_down,
                  n_donors, n_donors_earlier = n_low, n_donors_later = n_high,
                  n_cohorts)]
fwrite(out, file.path(DATA_DIR, "fig3_stage_degs_crosssectional_data.csv"))
fwrite(census[, .(kleiner_stage = fibrosis_stage, n_donors)],
       file.path(DATA_DIR, "fig3_stage_donor_census_crosssectional_data.csv"))

message("")
message("CAPTION (stage DEG panel): Cross-sectional differential expression ",
        "between donors at adjacent Kleiner fibrosis stages. Genes higher ",
        "(magenta) or lower (blue) at the later of the two stages, ",
        "limma-voom with quality weights, ~ dataset + inferred_sex + ",
        "fib_group, BH FDR < 0.05, no fold-change filter, 27,638 genes ",
        "tested per comparison. F1 vs F0 798 (442/356), F2 vs F1 2,190 ",
        "(1,510/680), F3 vs F2 3,490 (2,131/1,359), F4 vs F3 783 (504/279). ",
        "Six cohorts with exact Kleiner staging (GSE130970, GSE135251, ",
        "GSE162694, GSE174478, GSE193066, GSE240729); one pass-QC first ",
        "biopsy per donor. F4 vs F3 uses five cohorts because GSE193066 ",
        "contributes no F4 donors. These are between-donor contrasts, not ",
        "within-donor change.")
message("")
message("CAPTION (donor census panel): Donors contributing to each ",
        "cross-sectional comparison. Unique donors per stage F0-F4 = ",
        "126 / 187 / 174 / 133 / 44 (664 total, each donor at one stage ",
        "only). Of the 133 F3 donors, 108 enter the F4 vs F3 comparison; ",
        "GSE193066 is dropped there for lack of an F4 arm.")
message("")
message("Wrote:")
message("  ", file.path(PANEL_DIR, "fig3_stage_degs_crosssectional.pdf"))
message("  ", file.path(PANEL_DIR, "fig3_stage_donor_census_crosssectional.pdf"))
message("  ", file.path(DATA_DIR,  "fig3_stage_degs_crosssectional_data.csv"))
message("  ", file.path(DATA_DIR,  "fig3_stage_donor_census_crosssectional_data.csv"))
message("DONE — all assertions passed.")
