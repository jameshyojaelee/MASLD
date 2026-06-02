#!/usr/bin/env Rscript
# ============================================================================
# 349b_dedup_trajectory_figure.R
#
# Companion to Script 349. Builds a paralog/dedup-aware stage-trajectory CCC
# heatmap by collapsing the top-60 LR pairs into independent biological signals:
#
#   Step 1: per-ligand dedup — keep only the top-ranked LR pair (by
#           |Estimate| x -log10(p) on SH-vs-Healthy) for each ligand symbol.
#   Step 2: per-(ct_pair, ligand_family) collapse — collapse known paralog
#           families (HLA Class I, HLA Class II, collagens, FGFs, IGFs,
#           CCLs, CXCLs, APOA, complement, ADAM, etc.). Anything outside
#           the curated map falls back to the first-4-char prefix.
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/stage_lr_lmm_coarse.tsv
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/all_donor_lr_scores.tsv.gz
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
#
# Outputs:
#   figures/supplementary/stage_ccc/figS_stage_ccc_trajectory_dedup.pdf
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/dedup_top30_independent_signals.tsv
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

TOP_N_RAW   <- 60   # match Script 349 selection pool
TOP_N_FINAL <- 30   # dedup target

# ---------------------------------------------------------------------------
# Paralog family map (curated for the obvious cases). Anything not present
# falls back to a 4-char prefix.
# ---------------------------------------------------------------------------
paralog_map <- function(sym) {
  # Vectorised; sym is a character vector of ligand symbols
  fam <- character(length(sym))
  for (i in seq_along(sym)) {
    s <- sym[i]
    if (is.na(s) || s == "") { fam[i] <- NA_character_; next }
    # HLA Class I
    if (grepl("^HLA-[ABCEFG]$", s)) { fam[i] <- "HLA_classI"; next }
    # HLA Class II
    if (grepl("^HLA-D", s))         { fam[i] <- "HLA_classII"; next }
    # Collagens — split by TYPE (different chains, different biology).
    # Type I (fibrillar): COL1A1/A2; Type III commonly co-deposits and is
    # grouped with Type I for fibrotic-deposition signal.
    if (s %in% c("COL1A1", "COL1A2", "COL3A1")) {
      fam[i] <- "COL_typeI"; next
    }
    # Type IV (basement membrane)
    if (grepl("^COL4A[1-6]$", s))   { fam[i] <- "COL_typeIV"; next }
    # Type V (regulator of fibrillogenesis)
    if (grepl("^COL5A[1-3]$", s))   { fam[i] <- "COL_typeV"; next }
    # Type VI (microfibrillar)
    if (grepl("^COL6A[1-3]$", s))   { fam[i] <- "COL_typeVI"; next }
    # Other collagens (II, VII, VIII, IX, ...) — keep type-specific, do NOT
    # collapse across types. Fall through to per-gene "COL_typeX" tag.
    if (grepl("^COL([0-9]+)A[0-9]+$", s)) {
      ctype <- sub("^COL([0-9]+)A[0-9]+$", "\\1", s)
      fam[i] <- paste0("COL_type", as.roman(as.integer(ctype)))
      next
    }
    # FGFs
    if (grepl("^FGF[0-9]+$", s))    { fam[i] <- "FGF_family"; next }
    # IGFs (ligands, not BPs)
    if (s %in% c("IGF1", "IGF2"))   { fam[i] <- "IGF_family"; next }
    # CCL / CXCL chemokines
    if (grepl("^CCL[0-9]+$", s))    { fam[i] <- "CCL_chemokines"; next }
    if (grepl("^CXCL[0-9]+$", s))   { fam[i] <- "CXCL_chemokines"; next }
    # Apolipoproteins
    if (grepl("^APOA[0-9]?$", s))   { fam[i] <- "APOA_apolipoproteins"; next }
    if (grepl("^APOC[0-9]?$", s))   { fam[i] <- "APOC_apolipoproteins"; next }
    # Complement
    if (s %in% c("C3", "C4A", "C4B", "C4BPA", "C4BPB", "C5"))
                                    { fam[i] <- "Complement"; next }
    # ADAM metalloproteases
    if (grepl("^ADAM[0-9]+$", s))   { fam[i] <- "ADAM_metalloproteases"; next }
    # TGFB
    if (grepl("^TGFB[0-9]?$", s))   { fam[i] <- "TGFB_family"; next }
    # BMPs
    if (grepl("^BMP[0-9]+$", s))    { fam[i] <- "BMP_family"; next }
    # WNTs
    if (grepl("^WNT[0-9A-Z]+$", s)) { fam[i] <- "WNT_family"; next }
    # Notch ligands
    if (s %in% c("JAG1","JAG2","DLL1","DLL3","DLL4"))
                                    { fam[i] <- "Notch_ligands"; next }
    # Thrombospondins
    if (grepl("^THBS[0-9]+$", s))   { fam[i] <- "THBS_family"; next }
    # Serpins
    if (grepl("^SERPIN", s))        { fam[i] <- "SERPIN_family"; next }
    # Laminins
    if (grepl("^LAM[ABC][0-9]+$", s)) { fam[i] <- "LAM_family"; next }
    # Integrins (mostly receptors but treated for completeness)
    if (grepl("^ITG[AB][0-9A-Z]+$", s)) { fam[i] <- "ITG_family"; next }
    # Annexins
    if (grepl("^ANXA[0-9]+$", s))   { fam[i] <- "ANX_family"; next }
    # Fallback: 4-char prefix marker so reviewers can audit
    fam[i] <- paste0("prefix:", substr(s, 1, 4))
  }
  fam
}

# ---------------------------------------------------------------------------
# Load and rank
# ---------------------------------------------------------------------------
COARSE   <- fread(file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv"))
META_EXT <- fread(file.path(OUT_DIR, "donor_metadata_extended.tsv"))
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")

CONTINUOUS <- if (file.exists(file.path(OUT_DIR, "stage_lr_lmm_continuous.tsv"))) {
  fread(file.path(OUT_DIR, "stage_lr_lmm_continuous.tsv"))
} else NULL
BULK_CONC  <- if (file.exists(file.path(OUT_DIR, "lr_bulk_concordance.tsv"))) {
  fread(file.path(OUT_DIR, "lr_bulk_concordance.tsv"))
} else NULL
LOO_RATE   <- if (file.exists(file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv"))) {
  fread(file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv"))
} else NULL

sh_term <- "disease_stage_coarseSteatohepatitis"
sh <- COARSE[term == sh_term]
sh[, neglog10p := -log10(pmax(pval, 1e-30))]
sh[, rank_score := abs(Estimate) * neglog10p]

top_raw <- sh[order(-rank_score)][1:TOP_N_RAW]
cat(sprintf("[dedup] raw top-%d ligands: %d unique symbols\n",
            TOP_N_RAW, length(unique(top_raw$ligand_complex))))

# How bad is the redundancy?
dup_tbl <- top_raw[, .N, by = ligand_complex][N >= 2][order(-N)]
cat(sprintf("[dedup] %d ligands appear >=2x in raw top-%d (worst: %s = %d rows)\n",
            nrow(dup_tbl), TOP_N_RAW,
            ifelse(nrow(dup_tbl) > 0, dup_tbl$ligand_complex[1], "none"),
            ifelse(nrow(dup_tbl) > 0, dup_tbl$N[1], 0L)))

# ---------------------------------------------------------------------------
# Step 1: per-ligand dedup over the FULL SH table (not just the top-60),
# so we can backfill independent signals if dedup empties the pool.
# Keep best row per ligand_complex by rank_score.
# ---------------------------------------------------------------------------
setorder(sh, -rank_score)
per_ligand_best <- sh[, .SD[1], by = ligand_complex]
setorder(per_ligand_best, -rank_score)

# ---------------------------------------------------------------------------
# Step 2: assign paralog family and dedup per (ct_pair, family).
# Within each (ct_pair, family) keep only the single top-ranked row, and
# track the collapsed members for the side TSV.
# ---------------------------------------------------------------------------
per_ligand_best[, ligand_family := paralog_map(ligand_complex)]

# Collapse memo so we can list paralogs grouped into each kept row.
collapsed <- per_ligand_best[, .(
    members  = paste(sort(unique(ligand_complex)), collapse = ","),
    n_members = uniqueN(ligand_complex)
  ),
  by = .(ct_pair, ligand_family)]

setorder(per_ligand_best, ct_pair, ligand_family, -rank_score)
fam_best <- per_ligand_best[, .SD[1], by = .(ct_pair, ligand_family)]
fam_best <- merge(fam_best, collapsed, by = c("ct_pair", "ligand_family"),
                  all.x = TRUE)

setorder(fam_best, -rank_score)
final <- fam_best[1:min(TOP_N_FINAL, nrow(fam_best))]

cat(sprintf("[dedup] after per-ligand + paralog collapse: %d independent signals\n",
            nrow(final)))
cat(sprintf("[dedup] families with >=2 collapsed members: %d\n",
            sum(final$n_members >= 2)))

# ---------------------------------------------------------------------------
# Side TSV: collapsed mapping
# ---------------------------------------------------------------------------
out_tsv <- file.path(OUT_DIR, "dedup_top30_independent_signals.tsv")
fwrite(final[, .(rank      = .I,
                 ct_pair, ligand_family,
                 representative_ligand = ligand_complex,
                 receptor = receptor_complex,
                 lr_pair, source, target,
                 Estimate, pval, padj_within_ct,
                 rank_score,
                 paralogs_collapsed = members,
                 n_members)],
       out_tsv, sep = "\t")
cat(sprintf("[output] %s\n", out_tsv))

# ---------------------------------------------------------------------------
# Build heatmap (mirrors Script 349 structure)
# ---------------------------------------------------------------------------
lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long <- merge(lr_long,
                 META_EXT[, .(sample, disease_stage_coarse)],
                 by = "sample", all.x = TRUE)
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]

top_keys <- unique(final[, .(ct_pair, lr_pair)])
hm <- merge(lr_long, top_keys, by = c("ct_pair", "lr_pair"))
hm_mean <- hm[, .(score_mean = mean(score, na.rm = TRUE),
                  n_donors   = .N),
              by = .(ct_pair, lr_pair, disease_stage_coarse)]

hm_mean[, z := scale(score_mean)[, 1], by = .(ct_pair, lr_pair)]
hm_mean[, disease_stage_coarse := factor(disease_stage_coarse,
        levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]

# Label rows with an explicit family tag so the heatmap is self-documenting
final[, family_tag := ifelse(n_members >= 2,
                             sprintf(" [+%d paralogs:%s]",
                                     n_members - 1,
                                     ligand_family),
                             "")]
final[, pair_id := paste0(ct_pair, " | ", lr_pair, family_tag)]
hm_mean <- merge(hm_mean,
                 final[, .(ct_pair, lr_pair, pair_id)],
                 by = c("ct_pair", "lr_pair"))

ord <- final[order(-Estimate)]$pair_id
hm_mean[, pair_id := factor(pair_id, levels = ord)]

# Right-margin annotation
annot <- final[, .(pair_id, Estimate)]
setnames(annot, "Estimate", "coarse_SH_beta")

if (!is.null(CONTINUOUS)) {
  cont_sh <- CONTINUOUS[term == "macrophage_pseudotime_mean",
                        .(ct_pair = paste(source, target, sep = "->"),
                          lr_pair = paste(ligand_complex, receptor_complex,
                                          sep = "__"),
                          continuous_beta = Estimate,
                          continuous_p = pval)]
  cont_sh <- merge(cont_sh, final[, .(ct_pair, lr_pair, pair_id)],
                   by = c("ct_pair", "lr_pair"))
  annot <- merge(annot, cont_sh[, .(pair_id, continuous_beta)],
                 by = "pair_id", all.x = TRUE)
}
if (!is.null(BULK_CONC)) {
  bc <- unique(BULK_CONC[, .(ct_pair = paste(source, target, sep = "->"),
                              lr_pair = paste(ligand_complex, receptor_complex,
                                              sep = "__"),
                              bulk_both = both_concordant)])
  bc <- merge(bc, final[, .(ct_pair, lr_pair, pair_id)],
              by = c("ct_pair", "lr_pair"))
  annot <- merge(annot, bc[, .(pair_id, bulk_both)],
                 by = "pair_id", all.x = TRUE)
}
if (!is.null(LOO_RATE)) {
  lr <- LOO_RATE[, .(ct_pair, lr_pair, replication_rate)]
  lr <- merge(lr, final[, .(ct_pair, lr_pair, pair_id)],
              by = c("ct_pair", "lr_pair"))
  annot <- merge(annot, lr[, .(pair_id, replication_rate)],
                 by = "pair_id", all.x = TRUE)
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
# MASLD project palette (publication_theme.R `masld_colors`):
#   diverging: down (blue) -> white -> up (magenta)
DIV_LOW  <- masld_colors$down   # "#1565C0"
DIV_HIGH <- masld_colors$up     # "#C2185B"

p_hm <- ggplot(hm_mean,
               aes(x = disease_stage_coarse, y = pair_id, fill = z)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                       midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish,
                       name = "z(-log10 mag.rank)") +
  scale_x_discrete(labels = c("Healthy", "Steat.", "SH", "Cirr.")) +
  labs(x = NULL, y = NULL,
       title = sprintf("Stage-trajectory CCC: top %d independent LR signals",
                       nrow(final))) +
  theme_masld(base_size = 7) +
  theme(axis.text.y      = element_text(size = 5),
        axis.text.x      = element_text(angle = 30, hjust = 1, size = 6),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "top",
        legend.key.width = unit(0.6, "cm"),
        legend.key.height= unit(0.22, "cm"),
        legend.text      = element_text(size = 5),
        legend.title     = element_text(size = 6, face = "bold"),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 8,
                                        hjust = 0, margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(2, 2, 2, 2))

annot[, pair_id := factor(pair_id, levels = ord)]
annot_long <- melt(annot, id.vars = "pair_id",
                   measure.vars = intersect(c("continuous_beta",
                                               "bulk_both",
                                               "replication_rate"),
                                            names(annot)),
                   variable.name = "track", value.name = "value")
annot_long[, value := as.numeric(value)]
# Shorten track labels for compactness
annot_long[, track := factor(track,
  levels = c("continuous_beta", "bulk_both", "replication_rate"),
  labels = c("PT β", "bulk", "LOO"))]

p_annot <- ggplot(annot_long,
                  aes(x = track, y = pair_id, fill = value)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                       midpoint = 0, na.value = "grey90",
                       name = "value") +
  labs(x = NULL, y = NULL, title = "Tracks") +
  theme_masld(base_size = 6) +
  theme(axis.text.y      = element_blank(),
        axis.text.x      = element_text(angle = 30, hjust = 1, size = 5),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "top",
        legend.key.width = unit(0.35, "cm"),
        legend.key.height= unit(0.18, "cm"),
        legend.text      = element_text(size = 5),
        legend.title     = element_text(size = 5.5, face = "bold"),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 7,
                                        hjust = 0, margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(2, 2, 2, 2))

panel <- (p_hm | p_annot) + plot_layout(widths = c(5, 1.2))

out_pdf <- file.path(SUPP_DIR, "figS_stage_ccc_trajectory_dedup.pdf")
save_fig(panel, out_pdf, width = 7, height = 6.5)
cat(sprintf("[output] %s\n", out_pdf))

# ---------------------------------------------------------------------------
# Console summary for the reviewer-facing report
# ---------------------------------------------------------------------------
cat("\n========== DEDUP SUMMARY ==========\n")
cat(sprintf("Raw top-%d:           %d rows, %d unique ligands\n",
            TOP_N_RAW, nrow(top_raw), uniqueN(top_raw$ligand_complex)))
cat(sprintf("Ligands repeated >=2x: %d\n", nrow(dup_tbl)))
if (nrow(dup_tbl) > 0) {
  cat("  Top offenders:\n")
  for (i in seq_len(min(10, nrow(dup_tbl)))) {
    cat(sprintf("    %s : %d rows\n",
                dup_tbl$ligand_complex[i], dup_tbl$N[i]))
  }
}
cat(sprintf("After per-ligand dedup: %d signals\n",
            nrow(per_ligand_best[ligand_complex %in% top_raw$ligand_complex])))
cat(sprintf("After paralog-family collapse (final): %d independent signals\n",
            nrow(final)))
cat(sprintf("Families with collapsed paralogs (>=2 members): %d\n",
            sum(final$n_members >= 2)))
if (any(final$n_members >= 2)) {
  cat("  Collapsed families:\n")
  for (i in which(final$n_members >= 2)) {
    cat(sprintf("    %s | %s -> rep=%s, paralogs={%s}\n",
                final$ct_pair[i], final$ligand_family[i],
                final$ligand_complex[i], final$members[i]))
  }
}
cat("===================================\n")

cat("\n[done] dedup figure + TSV written\n")
