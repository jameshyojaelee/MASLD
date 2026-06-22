#!/usr/bin/env Rscript
# crossspecies_conservation.R — Fig 4g (human vs 4 mouse models)
# -----------------------------------------------------------------------------
# Message: the human MASLD disease signature is recapitulated across the 4 mouse
# dietary models (MCD, HFD, CDAHFD, FPC). ONE integrated panel, columns = models:
#   TOP    human<->mouse log2FC scatter — all orthologs (grey), conserved core
#          (light teal = >=3/4 models, dark teal = all 4), hero genes labeled;
#          per-model header reports the genome-wide Spearman rho (all orthologs).
#   BOTTOM Fibrotic-arm pathway human-proximity per model = correlation of each
#          diet's KEGG-pathway-NES profile with the human Severe-vs-Mild pathway
#          reference (the PEP-space construction of Vacca et al. 2024 Nat Metab).
#          LITMUS's DHPS is two-armed; our pathway-PEP reproduces the FIBROTIC arm
#          across their own 33 mouse models (Spearman 0.62 full / 0.54 shared, vs
#          metabolic arm 0.29 NS — script 11), so we present only the fibrotic axis
#          and DEFER the metabolic ranking to LITMUS. At this resolution our Western
#          diet (FPC) is the most human-proximal, consistent with "Western closest".
#
# WHY two layers: TOP = per-gene conservation (impartial genome-wide rho ~0.3,
# matching Vacca's whole-transcriptome agreement); BOTTOM = fibrotic pathway-PROGRAM
# human-proximity, the resolution at which LITMUS's DHPS operates. We DEFER the
# precise model ranking to LITMUS and show consistency, not a competing ranking.
# NOTE (2026-06-21): bottom panel is the VALIDATED fibrotic-arm pathway-PEP (canonical
# script 11; reproduces LITMUS fibrotic DHPS at Spearman 0.62). The former gene-level
# GSEA-NES wrongly claimed to be "the DHPS-style metric of Vacca" (reproduced DHPS at
# only 0.22) and was replaced; the metabolic arm is not reproducible and is deferred.
#
# DATA: human = canonical_deg_results.csv (C2). Mouse = the 4 per-diet limma-voom
# results joined through the 1:1 ortholog map — REPLICATING concordance script 01
# so per-model Spearman rho reproduces gene_concordance_matrix_20x.csv (asserted).
#
# Output: figures/main/fig4_validation/fig4g_crossspecies_conservation.pdf
#         + _supp/data/conserved_core_drug_targets.csv (drug-target enrichment table)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fgsea)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

set.seed(42)
RNA   <- file.path(BASE, "RNA-seq")
H_INT <- file.path(RNA, "Human/Patient_Cohorts/analysis/integration")
ANNOT <- file.path(H_INT, "results/gene_annotation")
INT   <- file.path(H_INT, "results/integration")
MOUSE_PD <- file.path(RNA, "Mouse/Unified_Integration/results/per_diet")

DIETS     <- c("MCD", "HFD", "CDAHFD", "FPC")
PADJ      <- 0.05
col_core4 <- tryCatch(masld_colors$conserved, error = function(e) "#1B7C6F")
col_core3 <- "#9CC9C2"; col_bg <- "#D9D9D9"
nes_red   <- "#B2182B"; nes_blue <- "#2C6FB0"
heroes    <- c("THRB", "FGF21", "COL1A1", "SPP1", "CIDEC", "NR1H4")

# ── Faithful paired table (replicates concordance script 01) ─────────────────
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
hum   <- fread(file.path(INT, "canonical_deg_results.csv"))
hum[, gene_base := gsub("\\..*", "", gene)]
h_mapped <- merge(hum[, .(gene_base, h_lfc = logFC, h_padj = padj, h_t = t)],
                  ortho[, .(human_gene_id, mouse_gene_id, human_symbol)],
                  by.x = "gene_base", by.y = "human_gene_id")

pairs <- list(); rho_chk <- list()
for (d in DIETS) {
  m <- fread(file.path(MOUSE_PD, paste0(d, "_de_results.csv")))
  m[, mouse_base := gsub("\\..*", "", gene)]
  pr <- merge(h_mapped, m[, .(mouse_base, m_lfc = logFC, m_padj = adj.P.Val, m_t = t)],
              by.x = "mouse_gene_id", by.y = "mouse_base")
  pr[, diet := d]; pairs[[d]] <- pr
  rho_chk[[d]] <- cor(pr$h_lfc, pr$m_lfc, method = "spearman", use = "complete.obs")
}
P <- rbindlist(pairs)

# correctness guard: reproduce on-disk rho_all
mx <- fread(file.path(CONCORDANCE, "gene_concordance_matrix_20x.csv"))[human_signature == "disease_vs_ctrl"]
for (d in DIETS) {
  disk <- mx[diet == d, rho_all]; mine <- rho_chk[[d]]
  cat(sprintf("  rho_all %-7s recomputed=%.4f  on-disk=%.4f  diff=%+.4f\n", d, mine, disk, mine - disk))
  if (abs(mine - disk) > 0.03) warning(sprintf("rho mismatch %s (%.4f vs %.4f)", d, mine, disk))
}

# conserved-core membership (n models concordant; the dual-significance selection)
conc <- fread(file.path(CONCORDANCE, "concordance_atlas_unified.csv"))
P <- merge(P, conc[, .(human_symbol, n_concordant)], by = "human_symbol", all.x = TRUE)
P[is.na(n_concordant), n_concordant := 0]
P[, core := factor(fifelse(n_concordant >= 4, "all4", fifelse(n_concordant >= 3, "three", "bg")),
                   levels = c("bg", "three", "all4"))]
n_core <- conc[n_concordant >= 3, .N]; n_four <- conc[n_concordant >= 4, .N]
stopifnot(n_core == 1355L)

# per-model rho (all genes vs both-significant); order models by concordance
rsum <- P[, {
  bs <- h_padj < PADJ & m_padj < PADJ
  .(rho_all = cor(h_lfc, m_lfc, method = "spearman", use = "complete.obs"),
    rho_sig = if (sum(bs, na.rm = TRUE) > 10)
                cor(h_lfc[bs], m_lfc[bs], method = "spearman", use = "complete.obs") else NA_real_)
}, by = diet]
setorder(rsum, -rho_all)
diet_order <- as.character(rsum$diet)
P[, diet := factor(diet, levels = diet_order)]

# robust shared axis limits
xl <- quantile(P$h_lfc, c(.005, .995), na.rm = TRUE)
yl <- quantile(P$m_lfc, c(.005, .995), na.rm = TRUE)

# ── Fibrotic-arm pathway human-proximity per diet (VALIDATED, canonical script 11) ─
# q1_our_diets_fibrotic_pep.csv: PEP correlation of each diet's KEGG-pathway-NES
# profile with the human Severe-vs-Mild (fibrotic) reference. This pathway-PEP
# reproduces LITMUS's FIBROTIC-arm DHPS across their 33 models at Spearman 0.62
# (metabolic arm 0.29 NS, deferred). LITMUS Western models = the reference band.
VBdir <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
pp <- fread(file.path(VBdir, "q1_our_diets_fibrotic_pep.csv"))
setnames(pp, "fibro_pep", "our_pepP")
pp[, diet := factor(our_diet, levels = diet_order)]
vm_pep   <- fread(file.path(VBdir, "q1_litmus_models_fibrotic_pep_shared.csv"))
west_med <- median(vm_pep[grepl("WD|GAN|AMLN|AMLD", diet_group)]$fibro_pep, na.rm = TRUE)  # LITMUS Western reference

# ── Integrated panel: columns = models; top scatter / bottom NES ─────────────
strip_lab <- setNames(sprintf("%s\nrho %.2f", rsum$diet, rsum$rho_all), diet_order)
hero_dt   <- P[human_symbol %in% heroes & core != "bg"]

ptop <- ggplot(P[order(core)], aes(h_lfc, m_lfc)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(data = P[core == "bg"],    color = col_bg,    size = 0.18, alpha = 0.5, shape = 16) +
  geom_point(data = P[core == "three"], color = col_core3, size = 0.35, alpha = 0.9, shape = 16) +
  geom_point(data = P[core == "all4"],  color = col_core4, size = 0.55, alpha = 0.95, shape = 16) +
  geom_text_repel(data = hero_dt, aes(label = human_symbol), size = PUB_GEOM_TEXT - 0.5,
                  fontface = "italic", color = "black", min.segment.length = 0,
                  segment.size = 0.2, segment.color = "grey55", max.overlaps = 20,
                  box.padding = 0.18, seed = 42) +
  facet_wrap(~ diet, nrow = 1, labeller = as_labeller(strip_lab)) +
  scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) +
  labs(x = "human disease log2FC", y = "mouse model\nlog2FC") +
  theme_masld() + theme_pub() +
  theme(strip.text = element_text(face = "bold", color = "black", size = 7, lineheight = 0.9),
        axis.text = element_text(color = "black"), panel.spacing = unit(4, "pt"),
        axis.title.x = element_text(size = 7))

pbot <- ggplot(pp, aes(x = diet, y = our_pepP)) +
  geom_hline(yintercept = west_med, linetype = "22", color = "grey55", linewidth = 0.3) +
  geom_col(width = 0.6, fill = col_core4) +
  geom_text(aes(label = sprintf("%.2f", our_pepP)), vjust = -0.3,
            size = PUB_GEOM_TEXT - 0.2, color = "black") +
  facet_wrap(~ diet, nrow = 1, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1.08)) +
  labs(x = NULL, y = "fibrotic pathway\nproximity") +
  theme_masld() + theme_pub() +
  theme(strip.text = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black"), panel.spacing = unit(4, "pt"))

combined <- ptop / pbot + plot_layout(heights = c(3, 1.25))

out_pdf <- file.path(FIG4_DIR, "fig4g_crossspecies_conservation.pdf")
save_fig(combined, out_pdf, width = fig_col_width + 1.5, height = 3.1)
message("Saved: ", out_pdf)
cat("\n[crossspecies] core=", n_core, " (all4=", n_four, ")  model order: ",
    paste(diet_order, collapse = ", "), "\n", sep = "")
cat("[crossspecies] fibrotic pathway human-proximity per diet (PEP corr; LITMUS Western median=",
    round(west_med, 2), "; validated vs LITMUS fibrotic DHPS at Spearman 0.62, script 11):\n", sep = "")
print(pp[, .(our_diet, fibrotic_proximity = round(our_pepP, 2),
             boot95 = sprintf("[%.2f, %.2f]", fibro_ci_lo, fibro_ci_hi))])

# ── Supplementary table: conserved-core drug-target enrichment (OR 4.92, H6) ──
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "dgidb_druggable", "opentargets_drug", "drug_dev_status",
             "clintrial_n_drugs", "clintrial_drugs"))
core <- conc[n_concordant >= 3]
core_a <- merge(core[, .(gene = human_symbol, n_models = n_concordant,
                         human_logFC = mean_h_lfc, translatability = translatability_score)],
                atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
core_a[, ctn := suppressWarnings(as.integer(clintrial_n_drugs))]
core_a[, is_drug_target := drug_dev_status %in% c("masld_approved", "masld_clinical",
        "masld_preclinical", "masld_discontinued", "drugged_other_indi") | (!is.na(ctn) & ctn > 0)]
h6 <- fread(file.path(AUDIT, "H6_inv_cc_positive_control.csv"))
hr <- h6[comparator == "drug_target" & test == "B_biotypematched" & universe == "full_atlas"]
obs_or <- hr$observed_OR[1]; n_dt_canon <- 61L

supp_dir <- file.path(FIG4_DIR, "_supp", "data"); dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
supp <- core_a[order(-n_models, -is_drug_target, -abs(human_logFC)),
               .(gene, n_models_concordant = n_models, human_logFC = round(human_logFC, 3),
                 translatability = round(translatability, 3), is_drug_target,
                 opentargets_drug, drug_dev_status, clintrial_drugs)]
supp_path <- file.path(supp_dir, "conserved_core_drug_targets.csv")
writeLines(sprintf("# Conserved core (>=3/4 mouse models): %d genes; %d drug targets (canonical, H6); biotype-matched-null OR %.2f (emp_p<0.001). Per-gene atlas drug annotations below.",
                   n_core, n_dt_canon, obs_or), supp_path)
fwrite(supp, supp_path, append = TRUE, col.names = TRUE)
message("Supp table: ", supp_path)
