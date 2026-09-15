#!/usr/bin/env Rscript
# Candidate supplementary panels for the AlphaGenome Atlas extension (Universe A, 29 direct-trait signals).
# KEY MESSAGE (panel A): at expression-colocalized direct-trait signals, predicted molecular effects are
#   channel-specific and incompletely covered; no single channel dominates and coverage is shown beside magnitude.
# KEY MESSAGE (panel B): posterior mass in measured accessible chromatin differs by lineage, and the few
#   disease-changed peaks under genetic variants are recorded with their native effects, not enriched.
# KEY MESSAGE (panel D): predicted-vs-measured eQTL direction is mostly unresolved or split, matching the
#   prior in-house benchmark.
# Individual PDFs only; 6 pt; control gray #9E9E9E; captions written to captions.txt, not drawn on panels.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT <- Sys.getenv("AGA_OUT_ROOT", ""); if (!nzchar(OUT)) stop("AGA_OUT_ROOT is unset")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)
T <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures"); if (dir.exists(FIG)) stop("Refusing to overwrite: ", FIG); dir.create(FIG)
GRAY <- "#9E9E9E"
save_panel <- function(p, name, w, h) ggsave(file.path(FIG, name), p, width = w, height = h, units = "in", device = cairo_pdf)
caption <- file(file.path(FIG, "captions.txt"), "w"); cap <- function(...) { message(...); writeLines(paste0(...), caption) }

sig <- fread(file.path(T, "eligible_signals.tsv"))[universe == "A_direct"]
sig[, label := paste0(gene, " | ", gsub("_EUR", "", gwas_name))]
sig[, label := ifelse(duplicated(label) | duplicated(label, fromLast = TRUE), paste0(label, " #", signal_pair_index), label)]
chan <- fread(file.path(T, "signal_channel_profiles.tsv"))[signal_uid %in% sig$signal_uid]
prof <- fread(file.path(T, "signal_profiles.tsv"))[signal_uid %in% sig$signal_uid]
keep <- c("expression", "splice_site_usage", "polyadenylation", "accessibility_atac", "accessibility_dnase", "histone", "tf_binding", "cage")
chan <- chan[channel %in% keep]
chan[, channel := factor(channel, levels = keep, labels = c("expression", "splice-site usage", "polyadenylation", "ATAC", "DNase", "histone", "TF ChIP", "CAGE"))]
avi <- prof[, .(signal_uid, channel = factor("AVI (top variant)"), value = avi_top_quantile, C_liver = avi_C, liver_tracks_available = TRUE)]
hm <- rbind(chan[, .(signal_uid, channel, value = A_quantile_over_C, C_liver, liver_tracks_available)], avi)
hm <- merge(hm, sig[, .(signal_uid, label)], by = "signal_uid")
hm[, label := factor(label, levels = sig[order(anchor_hg38_chrom, anchor_hg38_position)]$label)]
hm[, low_cov := !is.na(C_liver) & C_liver < 0.5]

pA <- ggplot(hm, aes(channel, label)) +
  geom_tile(aes(fill = value), colour = "white", linewidth = 0.2) +
  geom_point(data = hm[low_cov == TRUE], shape = 4, size = 0.8, colour = "black") +
  scale_fill_gradientn(colours = purple_gradient, limits = c(0, 1), na.value = GRAY, name = "posterior-weighted |quantile| / coverage") +
  theme_masld() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.key.height = unit(2, "mm")) +
  labs(x = NULL, y = NULL)
save_panel(pA, "panelA_mechanism_landscape.pdf", 3.4, 4.6)
cap("Panel A. Predicted molecular effects at the 29 expression-colocalized direct-trait signals (Universe A), ordered by genomic position. ",
    "Cell colour: posterior-weighted mean absolute Atlas calibrated score (0-1) over primary-liver and hepatocyte tracks, divided by covered posterior mass; ",
    "gray: no liver track for that channel; x: covered mass below 0.5. Weights are COLOC SNP.PP.H4 conditional on colocalization. AVI is the top-weight variant's calibrated AVI score (not liver-specific). ",
    "Values are within-scorer calibrated scores, not probabilities; no channel is selected as the mechanism.")

cov <- hm[, .(signal_uid, label, channel, C_liver)]
pAc <- ggplot(cov, aes(channel, label)) + geom_tile(aes(fill = C_liver), colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(colours = gray_gradient, limits = c(0, 1), na.value = "white", name = "covered posterior mass") +
  theme_masld() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.key.height = unit(2, "mm")) + labs(x = NULL, y = NULL)
save_panel(pAc, "panelA_coverage.pdf", 3.4, 4.6)
cap("Panel A coverage strip. Posterior mass of each signal carried by variants with a liver-track Atlas score for the channel (mapping, query floor, Atlas availability and liver-track availability all reduce it).")

# Panel B: posterior mass in accessible peaks by lineage
cen <- fread(file.path(T, "chromatin_census.tsv"))[signal_uid %in% sig$signal_uid]
cen <- merge(cen, sig[, .(signal_uid, label)], by = "signal_uid")
cols <- grep("_any_mass$", names(cen), value = TRUE)
long <- melt(cen[, c("label", cols), with = FALSE], id.vars = "label", variable.name = "lineage", value.name = "mass")
long[, lineage := gsub("_any_mass", "", lineage)]
long[, lineage := factor(lineage, levels = c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"))]
long[, label := factor(label, levels = levels(hm$label))]
pB <- ggplot(long, aes(mass, label, colour = lineage)) + geom_point(size = 1.2, alpha = 0.9) +
  scale_colour_manual(values = cat_palette[c(2, 3, 1, 4, 6)], name = NULL) + xlim(0, 1) +
  theme_masld() + theme(legend.position = "bottom") + labs(x = "posterior mass in accessible peaks (either snATAC cohort)", y = NULL)
save_panel(pB, "panelB_accessible_mass_by_lineage.pdf", 3.4, 4.6)
cap("Panel B. Posterior mass of each Universe A signal that falls inside consensus snATAC peaks accessible in GSE244832 or GSE281367, by lineage (mass counted once per variant per lineage). ",
    "Static accessibility only; disease change is available for hepatocyte and stellate cells alone (panel B2).")

vp <- fread(file.path(T, "variant_peak_overlaps.tsv.gz"))[signal_uid %in% sig$signal_uid & weight >= 0.01]
if (nrow(vp)) {
  vp[, evidence_state := factor(evidence_state, levels = c("supported", "discordant", "source_dependent", "indeterminate"))]
  pal <- c(supported = "#C9265E", discordant = "#1565C0", source_dependent = "#F4A674", indeterminate = GRAY)
  pB2 <- ggplot(vp, aes(logFC_gse244832, logFC_gse281367, colour = evidence_state, shape = lineage)) +
    geom_hline(yintercept = 0, colour = GRAY, linewidth = 0.2) + geom_vline(xintercept = 0, colour = GRAY, linewidth = 0.2) +
    geom_errorbar(aes(ymin = logFC_gse281367 - se_gse281367, ymax = logFC_gse281367 + se_gse281367), width = 0, linewidth = 0.2) +
    geom_errorbarh(aes(xmin = logFC_gse244832 - se_gse244832, xmax = logFC_gse244832 + se_gse244832), height = 0, linewidth = 0.2) +
    geom_point(size = 1.2) + scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    theme_masld() + theme(legend.position = "bottom") + labs(x = "log2FC MASH vs normal, GSE244832", y = "log2FC MASH vs normal, GSE281367")
  save_panel(pB2, "panelB2_variant_peak_disease_effects.pdf", 3.4, 3.4)
  cap("Panel B2. Native differential-accessibility effects (log2 fold change MASH vs normal, +/- SE, donor-level models) of hepatocyte and stellate peaks that contain a Universe A variant with posterior weight >= 0.01, ",
      "coloured by the existing peak evidence state (BH within jointly testable peaks per cohort). The allele-effect sign is not expected to equal the case-control sign.")
}

# Panel D: predicted vs measured eQTL direction mass
dc <- fread(file.path(T, "direction_comparisons.tsv"))[comparison == "predicted_vs_measured_eqtl" & signal_uid %in% sig$signal_uid]
dc <- merge(dc, sig[, .(signal_uid, label)], by = "signal_uid")
dl <- melt(dc[, .(label, concordant = concordant_mass, discordant = discordant_mass, unresolved = unresolved_mass)], id.vars = "label")
dl[, label := factor(label, levels = levels(hm$label))]
pD <- ggplot(dl, aes(value, label, fill = variable)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(concordant = "#00695C", discordant = "#C9265E", unresolved = GRAY), name = NULL) +
  theme_masld() + theme(legend.position = "bottom") + labs(x = "posterior mass", y = NULL)
save_panel(pD, "panelD_direction_mass.pdf", 3.4, 4.6)
cap("Panel D. For each Universe A signal, posterior mass of variants whose Atlas liver RNA prediction for the colocalized gene agrees (concordant) or disagrees (discordant) in sign with the measured liver eQTL effect of the same allele, ",
    "or is unresolved (palindromic, unmapped, unscored, or zero). A direction is called only at >= 0.95 of the original mass; nothing is renormalised.")
close(caption)
message("panels written to ", FIG)
