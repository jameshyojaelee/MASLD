#!/usr/bin/env Rscript
# 324_hep_stromal_circuits.R
#
# Analysis D2 + D3 (v1) — Hepatocyte-stromal + stromal-stromal CCC circuits.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   Uses LIANA differential LR pairs (already computed; 43,566 pairs).
#   1. D2 Hep-HSC bidirectional feedback (Hepatocytes <-> Fibroblasts as HSC proxy):
#        - Hep -> Fib ("activation signals"): LR enriched MASLD
#        - Fib -> Hep ("lipotoxic/metabolic return signals"): LR enriched MASLD
#        - Bidirectional gene pairs (genes acting as ligand in one dir + receptor
#          in the other)
#   2. D2 Hep-Macrophage bidirectional (parallel).
#   3. D3 Stromal-stromal (Fibroblasts, Endothelial, Macrophages among themselves).
#   4. Cross-reference with bulk reverse validation (B2 output) — prioritize
#      LR pairs that are ALSO bulk-concordant.
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/ccc/

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CCC_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
LIANA_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")

message("[1] Loading LIANA + B2 bulk-validated concordance table...")
liana <- fread(file.path(LIANA_DIR, "liana_differential_interactions.csv"))
b2    <- fread(file.path(CCC_DIR, "liana_bulk_concordance_perLR.csv"))
# Inner join on the LR keys (source, target, ligand, receptor)
key_cols <- c("source","target","ligand_complex","receptor_complex")
lr_ann <- merge(liana, b2[, c(key_cols,
                              "lig_lfc","lig_padj","rec_lfc","rec_padj",
                              "lig_concordant","rec_concordant","both_concordant"),
                          with = FALSE],
               by = key_cols, all.x = TRUE)
message(sprintf("  Annotated LR pairs: %d", nrow(lr_ann)))

SCORE_DIFF_MIN <- 0.10

# ---------------------------------------------------------------------------
# Helper: extract directional LR set between two cell-type labels
# ---------------------------------------------------------------------------
extract_dir <- function(dt, sender, receiver, direction_label) {
  dt[source == sender & target == receiver & abs(score_diff) >= SCORE_DIFF_MIN,
     .(direction = direction_label, source, target,
       ligand = ligand_complex, receptor = receptor_complex,
       score_masld, score_control, score_diff,
       lig_lfc, lig_padj, rec_lfc, rec_padj,
       lig_concordant, rec_concordant, both_concordant)][
     order(-score_diff)]
}

# ---------------------------------------------------------------------------
# D2 Hep <-> Fib (stellate proxy); Hep <-> Macrophage
# ---------------------------------------------------------------------------
message("[2] D2 Hep <-> Fib bidirectional...")
hep_fib_fwd <- extract_dir(lr_ann, "Hepatocytes", "Fibroblasts",      "Hep_to_Fib")
fib_hep_ret <- extract_dir(lr_ann, "Fibroblasts",  "Hepatocytes",      "Fib_to_Hep")
hep_fib <- rbind(hep_fib_fwd, fib_hep_ret, fill = TRUE)
fwrite(hep_fib, file.path(CCC_DIR, "D2_hep_HSC_bidirectional.csv"))
message(sprintf("  Hep->Fib: %d LR, Fib->Hep: %d LR", nrow(hep_fib_fwd), nrow(fib_hep_ret)))

message("[3] D2 Hep <-> Macrophage bidirectional...")
hep_mac_fwd <- extract_dir(lr_ann, "Hepatocytes", "Macrophages", "Hep_to_Mac")
mac_hep_ret <- extract_dir(lr_ann, "Macrophages", "Hepatocytes", "Mac_to_Hep")
hep_mac <- rbind(hep_mac_fwd, mac_hep_ret, fill = TRUE)
fwrite(hep_mac, file.path(CCC_DIR, "D2_hep_Mac_bidirectional.csv"))
message(sprintf("  Hep->Mac: %d LR, Mac->Hep: %d LR", nrow(hep_mac_fwd), nrow(mac_hep_ret)))

# Hep <-> Endothelial (LSEC proxy)
hep_end_fwd <- extract_dir(lr_ann, "Hepatocytes", "Endothelial cells", "Hep_to_LSEC")
end_hep_ret <- extract_dir(lr_ann, "Endothelial cells", "Hepatocytes", "LSEC_to_Hep")
hep_end <- rbind(hep_end_fwd, end_hep_ret, fill = TRUE)
fwrite(hep_end, file.path(CCC_DIR, "D2_hep_LSEC_bidirectional.csv"))
message(sprintf("  Hep->LSEC: %d LR, LSEC->Hep: %d LR", nrow(hep_end_fwd), nrow(end_hep_ret)))

# ---------------------------------------------------------------------------
# D2 Reciprocal gene pairs: genes that are ligand in one direction + receptor
# in the opposite direction. These represent TRUE feedback loops.
# ---------------------------------------------------------------------------
reciprocal_loops <- function(dt_fwd, dt_ret, fwd_name, ret_name) {
  # Find genes appearing as ligand in fwd AND as receptor in ret (and vice versa)
  fwd_lig <- unique(dt_fwd$ligand)
  fwd_rec <- unique(dt_fwd$receptor)
  ret_lig <- unique(dt_ret$ligand)
  ret_rec <- unique(dt_ret$receptor)

  # Gene X is ligand fwd and receptor ret
  genes_dual <- intersect(fwd_lig, ret_rec)
  # Gene Y is ligand ret and receptor fwd
  genes_dual2 <- intersect(ret_lig, fwd_rec)

  data.table(
    direction = c(rep(paste0(fwd_name, "->", ret_name, " (gene-is-ligand-fwd, receptor-ret)"),
                      length(genes_dual)),
                  rep(paste0(ret_name, "->", fwd_name, " (gene-is-ligand-ret, receptor-fwd)"),
                      length(genes_dual2))),
    gene = c(genes_dual, genes_dual2)
  )
}
message("[4] Reciprocal-gene feedback loops...")
recip_hepfib <- reciprocal_loops(hep_fib_fwd, fib_hep_ret, "Hep", "Fib")
recip_hepmac <- reciprocal_loops(hep_mac_fwd, mac_hep_ret, "Hep", "Mac")
recip_heplsec <- reciprocal_loops(hep_end_fwd, end_hep_ret, "Hep", "LSEC")
recip_all <- rbind(
  recip_hepfib[, axis := "Hep<->Fib"],
  recip_hepmac[, axis := "Hep<->Mac"],
  recip_heplsec[, axis := "Hep<->LSEC"])
fwrite(recip_all, file.path(CCC_DIR, "D2_reciprocal_feedback_genes.csv"))
message(sprintf("  Reciprocal genes: Hep<->Fib=%d, Hep<->Mac=%d, Hep<->LSEC=%d",
                nrow(recip_hepfib), nrow(recip_hepmac), nrow(recip_heplsec)))

# ---------------------------------------------------------------------------
# D3 Stromal-stromal (Fibroblasts, Endothelial cells, Macrophages among themselves)
# ---------------------------------------------------------------------------
stromal_cts <- c("Fibroblasts", "Endothelial cells", "Macrophages")
stromal_stromal <- lr_ann[source %in% stromal_cts & target %in% stromal_cts &
                          source != target & abs(score_diff) >= SCORE_DIFF_MIN,
                          .(source, target,
                            ligand = ligand_complex, receptor = receptor_complex,
                            score_masld, score_control, score_diff,
                            lig_lfc, lig_padj, rec_lfc, rec_padj,
                            lig_concordant, rec_concordant, both_concordant)][
                          order(-score_diff)]
fwrite(stromal_stromal, file.path(CCC_DIR, "D3_stromal_stromal.csv"))
message(sprintf("  Stromal-stromal LR pairs: %d", nrow(stromal_stromal)))

# ---------------------------------------------------------------------------
# Summary with top pairs for each axis
# ---------------------------------------------------------------------------
top_masld <- function(dt, n = 12) {
  if (nrow(dt) == 0) return(NULL)
  dt[order(-score_diff)][1:min(n, nrow(dt)),
     .(direction, source, target, ligand, receptor, score_diff,
       both_concordant)]
}

summary_lines <- c(
  sprintf("Hep <-> Fib (HSC proxy): %d fwd, %d ret, %d bulk-concordant total",
          nrow(hep_fib_fwd), nrow(fib_hep_ret), sum(hep_fib$both_concordant, na.rm = TRUE)),
  "",
  "Top MASLD-enriched Hep->Fib + Fib->Hep (bulk-concordant highlighted):",
  capture.output(print(top_masld(hep_fib), nrows = 20)),
  "",
  sprintf("Hep <-> Mac: %d fwd, %d ret", nrow(hep_mac_fwd), nrow(mac_hep_ret)),
  "Top Hep->Mac + Mac->Hep:",
  capture.output(print(top_masld(hep_mac), nrows = 20)),
  "",
  sprintf("Hep <-> LSEC: %d fwd, %d ret", nrow(hep_end_fwd), nrow(end_hep_ret)),
  "Top Hep->LSEC + LSEC->Hep:",
  capture.output(print(top_masld(hep_end), nrows = 15)),
  "",
  sprintf("Reciprocal feedback genes (acting as ligand one way + receptor other): %d",
          nrow(recip_all)),
  capture.output(print(recip_all, nrows = 30)),
  "",
  sprintf("D3 Stromal-stromal LR pairs: %d", nrow(stromal_stromal)),
  "Top stromal-stromal (MASLD-enriched):",
  capture.output(print(stromal_stromal[, .(source,target,ligand,receptor,score_diff,both_concordant)][1:20], nrows = 20))
)
writeLines(summary_lines, file.path(CCC_DIR, "D2_D3_circuits_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", CCC_DIR)
