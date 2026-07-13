#!/usr/bin/env Rscript
# =============================================================================
# Figure 13 TPM sensitivity -- final library size after all current library
# rules, sweeping MCD TPM thresholds while reusing the exact same tier-building
# logic as figS_cas13_library_options.R.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Reuse the canonical Cas13 library option builder so the TPM sensitivity curves
# and the rebuilt library stay on the same exact rule set.
source(file.path(BASE, "scripts/figures/figS_cas13_library_options.R"), local = FALSE)

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REF_TPM_PC  <- 1.0
REF_TPM_LNC <- 0.1

gate_tpm_custom <- function(g, pc_thr, lnc_thr) {
  pc <- intersect(g, pc_ids)
  lnc <- intersect(g, lnc_ids)

  pc_tpm <- unname(tpm_of[pc]);  pc_tpm[is.na(pc_tpm)] <- 0
  lnc_tpm <- unname(tpm_of[lnc]); lnc_tpm[is.na(lnc_tpm)] <- 0

  exempt_g <- function(ids) {
    (ids %in% posctrl_mouse) | (!is.na(coloc_pp4_of[ids]) & coloc_pp4_of[ids] >= COLOC_EXEMPT)
  }

  pc_ok <- pc[pc_tpm >= pc_thr | exempt_g(pc)]
  lnc_ok <- lnc[lnc_tpm >= lnc_thr | exempt_g(lnc)]
  union(pc_ok, lnc_ok)
}

final_count <- function(pc_thr = REF_TPM_PC, lnc_thr = REF_TPM_LNC) {
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

  final <- setdiff(gate_tpm_custom(ungated, pc_thr, lnc_thr), ung_ids)
  list(
    total = length(final),
    pc = length(intersect(final, pc_ids)),
    lnc = length(intersect(final, lnc_ids))
  )
}

sweep_lnc <- function(cuts) {
  rbindlist(lapply(cuts, function(x) {
    r <- final_count(pc_thr = REF_TPM_PC, lnc_thr = x)
    data.table(parameter = "lncRNA_TPM", cutoff = x, n_lnc = r$lnc, n_pc = r$pc, n_total = r$total)
  }))
}

sweep_pc <- function(cuts) {
  rbindlist(lapply(cuts, function(x) {
    r <- final_count(pc_thr = x, lnc_thr = REF_TPM_LNC)
    data.table(parameter = "PCG_TPM", cutoff = x, n_lnc = r$lnc, n_pc = r$pc, n_total = r$total)
  }))
}

cuts_tpm_lnc <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
cuts_tpm_pc  <- c(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0)

cat("Sweeping lncRNA TPM cutoff...\n")
s_lnc <- sweep_lnc(cuts_tpm_lnc)
print(s_lnc)

cat("Sweeping PCG TPM cutoff...\n")
s_pcg <- sweep_pc(cuts_tpm_pc)
print(s_pcg)

fwrite(rbindlist(list(s_lnc, s_pcg), fill = TRUE),
       file.path(BASE, "Cas13_Library_Design/data/13_tpm_cutoff_sensitivity.csv"))

pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
col_lnc <- "#C0143C"
col_pc  <- "#2271B2"

p13a <- ggplot(s_lnc, aes(cutoff, n_lnc)) +
  geom_vline(xintercept = REF_TPM_LNC, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = REF_TPM_LNC + 0.02, y = max(s_lnc$n_lnc) * 0.98,
           label = paste0("reference\n(", REF_TPM_LNC, ")"),
           hjust = 0, vjust = 1, size = 2.8, color = "grey35") +
  geom_line(linewidth = 1, color = col_lnc) +
  geom_point(size = 2.5, color = col_lnc) +
  geom_text(data = s_lnc[cutoff == REF_TPM_LNC], aes(label = n_lnc),
            vjust = -1, hjust = 0.5, size = 3.5, fontface = "bold", color = col_lnc) +
  scale_x_continuous(breaks = cuts_tpm_lnc) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Mouse MCD TPM cutoff", y = "lncRNA genes in final library",
       title = "lncRNA: library size vs TPM cutoff") +
  theme_masld(base_size = 11) + theme_pub()

ggsave(file.path(OUT_DIR, "13a_lncrna_count_vs_tpm_cutoff.pdf"), p13a,
       width = 6, height = 4.5, device = pdf_dev)

p13b <- ggplot(s_pcg, aes(cutoff, n_pc)) +
  geom_vline(xintercept = REF_TPM_PC, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = REF_TPM_PC + 0.05, y = max(s_pcg$n_pc) * 0.98,
           label = paste0("reference\n(", REF_TPM_PC, ")"),
           hjust = 0, vjust = 1, size = 2.8, color = "grey35") +
  geom_line(linewidth = 1, color = col_pc) +
  geom_point(size = 2.5, color = col_pc) +
  geom_text(data = s_pcg[cutoff == REF_TPM_PC], aes(label = n_pc),
            vjust = -1, hjust = 0.5, size = 3.5, fontface = "bold", color = col_pc) +
  scale_x_continuous(breaks = cuts_tpm_pc) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Mouse MCD TPM cutoff", y = "PCG genes in final library",
       title = "PCG: library size vs TPM cutoff") +
  theme_masld(base_size = 11) + theme_pub()

ggsave(file.path(OUT_DIR, "13b_pcg_count_vs_tpm_cutoff.pdf"), p13b,
       width = 6, height = 4.5, device = pdf_dev)

cat("fig13_tpm_cutoff_sensitivity.R complete.\n")
