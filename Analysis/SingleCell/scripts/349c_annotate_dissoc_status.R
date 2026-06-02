#!/usr/bin/env Rscript
# ============================================================================
# 349c_annotate_dissoc_status.R
#
# Companion to Script 349b. Cross-references the paralog-dedup top-30
# independent CCC signals with the GSE136103-holdout dissociation
# sensitivity test (Script 347b output) and writes an annotated TSV +
# heatmap PDF flagging dissociation-suspect entries.
#
# A dedup row collapses 1+ ligand symbols into one family-level signal
# (e.g. rank-22 collapses COL4A1+COL4A2 under COL_typeIV). The dissoc
# test, by contrast, is per-LR-pair. We therefore explode each dedup row
# back into its constituent (ct_pair, ligand, receptor) members, join
# against 347b, and pick the WORST individual-member status as the
# row-level call (precautionary: if ANY collapsed paralog weakens or
# reverses, the whole family-level signal carries a caveat).
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/dedup_top30_independent_signals.tsv
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/cirrhosis_dissoc_sensitivity.tsv
#
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/dedup_top30_with_dissoc.tsv
#   figures/supplementary/stage_ccc/figS_stage_ccc_trajectory_dedup_annotated.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
SUPP_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

DEDUP_TSV  <- file.path(OUT_DIR, "dedup_top30_independent_signals.tsv")
# Prefer the combined SH+Cirrhosis dissoc table (Script 347c) when it exists;
# fall back to the Cirrhosis-only 347b output for backward compatibility.
DISSOC_COMBINED <- file.path(OUT_DIR, "all_dissoc_sensitivity.tsv")
DISSOC_CIRR     <- file.path(OUT_DIR, "cirrhosis_dissoc_sensitivity.tsv")
DISSOC_TSV <- if (file.exists(DISSOC_COMBINED)) DISSOC_COMBINED else DISSOC_CIRR
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")
META_EXT   <- file.path(OUT_DIR, "donor_metadata_extended.tsv")

stopifnot(file.exists(DEDUP_TSV), file.exists(DISSOC_TSV),
          file.exists(MERGED_TSV), file.exists(META_EXT))

dedup  <- fread(DEDUP_TSV)
dissoc <- fread(DISSOC_TSV)
cat(sprintf("[input] dissoc table: %s (%d rows)\n", DISSOC_TSV, nrow(dissoc)))
if ("contrast" %in% names(dissoc)) {
  cat("[input] contrast breakdown:\n")
  print(table(dissoc$contrast, useNA = "ifany"))
}

# The dedup top-30 figure is built on the SH-vs-Healthy contrast (Script
# 349b). When the combined table carries both contrasts, prefer the
# SH-emergent rows for any (ct_pair, lr_pair) tested under both. If only
# Cirrhosis-tested rows exist for a pair, keep them as a secondary witness.
if ("contrast" %in% names(dissoc)) {
  dissoc[, contrast_pref :=
           ifelse(contrast == "SH_vs_Healthy", 1L,
           ifelse(contrast == "Cirrhosis_vs_Healthy", 2L, 3L))]
  setorder(dissoc, ct_pair, lr_pair, contrast_pref)
  dissoc <- dissoc[, .SD[1L], by = .(ct_pair, lr_pair)]
  dissoc[, contrast_pref := NULL]
}

# ---------------------------------------------------------------------------
# Precedence order for the family-level "worst" call. Suspect = anything
# that does not survive (weakened / reversed / lost). untestable is below
# survives_dissoc but above not_in_top100.
# ---------------------------------------------------------------------------
status_levels <- c("reversed", "weakened", "lost",
                   "untestable", "survives_dissoc", "not_tested")
worst_status <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return("not_tested")
  for (lvl in status_levels) if (lvl %in% v) return(lvl)
  v[1]
}

# ---------------------------------------------------------------------------
# Explode dedup -> (ct_pair, ligand_member) rows, join against 347b
# ---------------------------------------------------------------------------
exploded <- dedup[, .(
    rank = rank,
    ct_pair = ct_pair,
    ligand_family = ligand_family,
    representative_ligand = representative_ligand,
    receptor = receptor,
    member = unlist(strsplit(paralogs_collapsed, ","))
  ), by = seq_len(nrow(dedup))]
exploded[, member_lr_pair := paste(member, receptor, sep = "__")]

# 347b/c is per (ct_pair, lr_pair). Direct join — if a paralog member
# isn't in either top-100 (SH or Cirrhosis), it's marked not_tested.
diss_cols <- c("ct_pair", "lr_pair", "full_pval", "holdout_pval",
               "full_estimate", "holdout_estimate", "classification")
if ("contrast" %in% names(dissoc)) diss_cols <- c(diss_cols, "contrast")
diss_key <- dissoc[, ..diss_cols]
setnames(diss_key,
         old = c("full_pval", "holdout_pval",
                 "full_estimate", "holdout_estimate", "classification"),
         new = c("member_full_pval", "member_holdout_pval",
                 "member_full_estimate", "member_holdout_estimate",
                 "member_classification"))
if ("contrast" %in% names(diss_key))
  setnames(diss_key, "contrast", "member_contrast")
setnames(diss_key, "lr_pair", "member_lr_pair")

annot <- merge(exploded, diss_key,
               by = c("ct_pair", "member_lr_pair"),
               all.x = TRUE)
annot[is.na(member_classification), member_classification := "not_tested"]

# Map dissoc-test class -> row-level dissoc_status. Family-level pick =
# worst across collapsed members.
row_status <- annot[, .(
    dissoc_status = worst_status(member_classification),
    n_members_tested = sum(member_classification != "not_in_top100"),
    n_members_total = .N,
    members_weakened = paste(member[member_classification == "weakened"],
                             collapse = ","),
    members_reversed = paste(member[member_classification == "reversed"],
                             collapse = ","),
    members_lost     = paste(member[member_classification == "lost"],
                             collapse = ",")
  ),
  by = .(rank, ct_pair, ligand_family)]

out_dt <- merge(dedup, row_status,
                by = c("rank", "ct_pair", "ligand_family"),
                sort = FALSE)
setorder(out_dt, rank)

out_tsv <- file.path(OUT_DIR, "dedup_top30_with_dissoc.tsv")
fwrite(out_dt, out_tsv, sep = "\t")
cat(sprintf("[output] %s\n", out_tsv))

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
cat("\n========== DISSOC ANNOTATION SUMMARY ==========\n")
cat("dissoc_status table:\n")
print(table(out_dt$dissoc_status, useNA = "ifany"))
cat("\nSuspect rows (weakened / reversed):\n")
print(out_dt[dissoc_status %in% c("weakened", "reversed"),
             .(rank, ct_pair, ligand_family, representative_ligand,
               receptor, dissoc_status, members_weakened, members_reversed)])
cat("===============================================\n")

# ---------------------------------------------------------------------------
# Heatmap re-draw with right-margin "star dissoc" annotation
# ---------------------------------------------------------------------------
lr_long <- fread(MERGED_TSV)
meta    <- fread(META_EXT)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long <- merge(lr_long,
                 meta[, .(sample, disease_stage_coarse)],
                 by = "sample", all.x = TRUE)
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]

top_keys <- unique(out_dt[, .(ct_pair, lr_pair)])
hm <- merge(lr_long, top_keys, by = c("ct_pair", "lr_pair"))
hm_mean <- hm[, .(score_mean = mean(score, na.rm = TRUE),
                  n_donors   = .N),
              by = .(ct_pair, lr_pair, disease_stage_coarse)]
hm_mean[, z := scale(score_mean)[, 1], by = .(ct_pair, lr_pair)]
hm_mean[, disease_stage_coarse := factor(disease_stage_coarse,
        levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]

out_dt[, family_tag := ifelse(n_members >= 2,
                              sprintf(" [+%d paralogs:%s]",
                                      n_members - 1,
                                      ligand_family),
                              "")]
out_dt[, suspect_mark := ifelse(dissoc_status %in% c("weakened", "reversed"),
                                " * dissoc", "")]
out_dt[, pair_id := paste0(ct_pair, " | ", lr_pair, family_tag, suspect_mark)]
hm_mean <- merge(hm_mean,
                 out_dt[, .(ct_pair, lr_pair, pair_id)],
                 by = c("ct_pair", "lr_pair"))

ord <- out_dt[order(-Estimate)]$pair_id
hm_mean[, pair_id := factor(pair_id, levels = ord)]

# MASLD project palette (publication_theme.R `masld_colors`):
#   diverging: down (blue) -> white -> up (magenta)
DIV_LOW  <- masld_colors$down   # "#1565C0"
DIV_HIGH <- masld_colors$up     # "#C2185B"

# Build a small dissoc strip column on the right.
# Dissoc-status palette: green = survives, magenta/dark-magenta = caveats,
# matching MASLD project semantics (gray = neutral/missing, magenta = disease/risk).
strip <- out_dt[, .(pair_id, dissoc_status)]
strip[, pair_id := factor(pair_id, levels = ord)]
status_color <- c(survives_dissoc = masld_colors$conserved,  # teal
                  weakened        = masld_colors$mash,        # magenta
                  reversed        = masld_colors$fibrosis,    # dark magenta
                  lost            = "#7B1FA2",                # violet
                  untestable      = masld_colors$ns,          # gray
                  not_tested      = "#F5F5F5")                # off-white

p_hm <- ggplot(hm_mean,
               aes(x = disease_stage_coarse, y = pair_id, fill = z)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                       midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish,
                       name = "z(-log10 mag.rank)") +
  scale_x_discrete(labels = c("Healthy", "Steat.", "SH", "Cirr.")) +
  labs(x = NULL, y = NULL,
       title = sprintf("Top %d independent LR signals with dissociation annotation",
                       nrow(out_dt))) +
  theme_masld(base_size = 7) +
  theme(axis.text.y      = element_text(size = 5.5),
        axis.text.x      = element_text(angle = 30, hjust = 1, size = 7),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "top",
        legend.key.width = unit(0.6, "cm"),
        legend.key.height= unit(0.22, "cm"),
        legend.text      = element_text(size = 5.5),
        legend.title     = element_text(size = 6, face = "bold"),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 9, hjust = 0,
                                        margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(3, 3, 3, 3))

p_strip <- ggplot(strip,
                  aes(x = 1, y = pair_id, fill = dissoc_status)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_text(data = strip[dissoc_status %in% c("weakened", "reversed")],
            aes(label = "*"), color = "white", size = 2.6) +
  scale_fill_manual(values = status_color, name = "dissoc.",
                    breaks = c("survives_dissoc", "weakened", "reversed",
                               "lost", "untestable", "not_tested")) +
  labs(x = NULL, y = NULL, title = "GSE136103 holdout") +
  theme_masld(base_size = 6) +
  theme(axis.text.y      = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "right",
        legend.key.size  = unit(0.25, "cm"),
        legend.text      = element_text(size = 5.5),
        legend.title     = element_text(size = 6, face = "bold"),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 7, hjust = 0,
                                        margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(3, 3, 3, 3))

panel <- (p_hm | p_strip) + plot_layout(widths = c(4, 1.4))

out_pdf <- file.path(SUPP_DIR, "figS_stage_ccc_trajectory_dedup_annotated.pdf")
save_fig(panel, out_pdf, width = 8.5, height = 7)
cat(sprintf("[output] %s\n", out_pdf))

cat("\n[done] dissoc-annotated dedup figure + TSV written\n")
