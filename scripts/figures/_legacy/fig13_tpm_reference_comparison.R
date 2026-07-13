#!/usr/bin/env Rscript
# =============================================================================
# Compare final library size vs expression cutoff using two mouse references:
#   1) in-house MCD bulk liver TPM (genuine length-normalized TPM)
#   2) mouse single-cell hepatocyte pseudobulk CPM
#
# Reuses the exact library tier logic from figS_cas13_library_options.R so only
# the expression reference changes.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(GenomicFeatures)
  library(rtracklayer)
  library(GenomicRanges)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/figS_cas13_library_options.R"), local = FALSE)

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REF_TPM_PC  <- 1.0
REF_TPM_LNC <- 0.1
GTF <- "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz"
PBDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mouse_sc/pseudobulk")

cat("Building single-cell hepatocyte CPM reference...\n")
gtf <- rtracklayer::import(GTF)
exons <- gtf[gtf$type == "exon"]
gene_len <- sum(width(reduce(split(exons, exons$gene_name))))
gene_len <- as.numeric(gene_len)
names(gene_len) <- names(sum(width(reduce(split(exons, exons$gene_name)))))

hep <- fread(file.path(PBDIR, "Hepatocytes_pseudobulk.csv"))
setnames(hep, 1, "gene_symbol")
hep <- hep[gene_symbol %in% names(gene_len)]
mat <- as.matrix(hep[, -1])
rownames(mat) <- hep$gene_symbol
len_kb <- gene_len[rownames(mat)] / 1e3
rpk <- sweep(mat, 1, len_kb, "/")
sf  <- colSums(rpk) / 1e6
sf[sf == 0] <- 1
tpm <- sweep(rpk, 2, sf, "/")
hep_tpm_symbol <- rowMeans(tpm)

meta <- fread(file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv"))
setnames(meta, c("mouse_ensembl_base", "mouse_symbol_gtf"), c("gene_id_mouse", "gene_symbol_mouse"), skip_absent = TRUE)
sc_tpm_dt <- unique(meta[, .(gene_id_mouse, gene_symbol_mouse)])
sc_tpm_dt[, sc_hep_tpm := unname(hep_tpm_symbol[gene_symbol_mouse])]
sc_tpm_dt[is.na(sc_hep_tpm), sc_hep_tpm := 0]
sc_tpm_of <- setNames(sc_tpm_dt$sc_hep_tpm, sc_tpm_dt$gene_id_mouse)

gate_tpm_custom <- function(g, pc_thr, lnc_thr, ref_vec) {
  pc <- intersect(g, pc_ids)
  lnc <- intersect(g, lnc_ids)
  pc_tpm <- unname(ref_vec[pc]); pc_tpm[is.na(pc_tpm)] <- 0
  lnc_tpm <- unname(ref_vec[lnc]); lnc_tpm[is.na(lnc_tpm)] <- 0

  exempt_g <- function(ids) {
    (ids %in% posctrl_mouse) | (!is.na(coloc_pp4_of[ids]) & coloc_pp4_of[ids] >= COLOC_EXEMPT)
  }
  pc_ok <- pc[pc_tpm >= pc_thr | exempt_g(pc)]
  lnc_ok <- lnc[lnc_tpm >= lnc_thr | exempt_g(lnc)]
  union(pc_ok, lnc_ok)
}

final_count <- function(ref_vec, pc_thr = REF_TPM_PC, lnc_thr = REF_TPM_LNC) {
  core_set <- core_mouse_for(LFC_LNC)
  coh <- cohort_mouse_for(LFC_LNC, CHOSEN_COHORT)
  mo  <- mouse_conf_for(LFC_LNC)
  co  <- coloc_mouse
  ct  <- posctrl_mouse
  ct_g   <- ct
  core_g <- setdiff(core_set, ct_g)
  coh_g  <- setdiff(coh, union(ct_g, core_g))
  mo_g   <- setdiff(mo,  Reduce(union, list(ct_g, core_g, coh_g)))
  co_g   <- setdiff(co,  Reduce(union, list(ct_g, core_g, coh_g, mo_g)))
  ungated <- Reduce(union, list(ct_g, core_g, coh_g, mo_g, co_g))
  final <- setdiff(gate_tpm_custom(ungated, pc_thr, lnc_thr, ref_vec), ung_ids)
  list(total = length(final), pc = length(intersect(final, pc_ids)), lnc = length(intersect(final, lnc_ids)))
}

sweep_ref <- function(ref_name, ref_vec, cuts, mode = c("lnc", "pc")) {
  mode <- match.arg(mode)
  rbindlist(lapply(cuts, function(x) {
    r <- if (mode == "lnc") final_count(ref_vec, pc_thr = REF_TPM_PC, lnc_thr = x)
         else final_count(ref_vec, pc_thr = x, lnc_thr = REF_TPM_LNC)
    data.table(reference = ref_name, cutoff = x, n_lnc = r$lnc, n_pc = r$pc, n_total = r$total)
  }))
}

cuts_tpm_lnc <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
cuts_tpm_pc  <- c(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0)

s_lnc <- rbindlist(list(
  sweep_ref("MCD bulk TPM", tpm_of, cuts_tpm_lnc, "lnc"),
  sweep_ref("sc hepatocyte CPM", sc_tpm_of, cuts_tpm_lnc, "lnc")
))
s_pcg <- rbindlist(list(
  sweep_ref("MCD bulk TPM", tpm_of, cuts_tpm_pc, "pc"),
  sweep_ref("sc hepatocyte CPM", sc_tpm_of, cuts_tpm_pc, "pc")
))

fwrite(rbindlist(list(
  s_lnc[, parameter := "lncRNA_TPM_compare"],
  s_pcg[, parameter := "PCG_TPM_compare"]
), fill = TRUE), file.path(BASE, "Cas13_Library_Design/data/13_tpm_reference_comparison.csv"))

cols <- c("MCD bulk TPM" = "#C0143C", "sc hepatocyte CPM" = "#2271B2")
pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

p_lnc <- ggplot(s_lnc, aes(cutoff, n_lnc, color = reference)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = REF_TPM_LNC, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = cuts_tpm_lnc) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Expression cutoff", y = "lncRNA genes in final library",
       title = "lncRNA: MCD bulk TPM vs sc hepatocyte CPM") +
  theme_masld(base_size = 11) + theme_pub() +
  theme(legend.position = "top")

p_pcg <- ggplot(s_pcg, aes(cutoff, n_pc, color = reference)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = REF_TPM_PC, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = cuts_tpm_pc) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Expression cutoff", y = "PCG genes in final library",
       title = "PCG: MCD bulk TPM vs sc hepatocyte CPM") +
  theme_masld(base_size = 11) + theme_pub() +
  theme(legend.position = "top")

ggsave(file.path(OUT_DIR, "13a_lncrna_count_vs_tpm_cutoff_compare.pdf"), p_lnc,
       width = 6.5, height = 4.7, device = pdf_dev)
ggsave(file.path(OUT_DIR, "13b_pcg_count_vs_tpm_cutoff_compare.pdf"), p_pcg,
       width = 6.5, height = 4.7, device = pdf_dev)

cat("fig13_tpm_reference_comparison.R complete.\n")
