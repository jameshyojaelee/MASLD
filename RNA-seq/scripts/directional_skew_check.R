#!/usr/bin/env Rscript
# directional_skew_check.R
# Tests whether the up-skew in the STAR -s 2 mega DEGs is driven by cell-type
# composition vs hepatocyte-intrinsic biology.
#   Part A: directional ORA (up vs down DEGs) — Hallmark + curated liver cell-type markers.
#   Part B: composition decomposition from the -s 2 deconvolution attribution.
# Output: RNA-seq/results/audit_sensitivity/directional_skew/

suppressPackageStartupMessages({
  library(data.table); library(msigdbr); library(clusterProfiler)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/directional_skew")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PADJ <- 0.05; LFC <- 0.5

# ---- gene-ID -> symbol map (from atlas) ----
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas[, ens := sub("\\..*", "", ensembl_id)]
ens2sym <- setNames(atlas$human_symbol, atlas$ens)
sym2ens <- setNames(atlas$ens, atlas$human_symbol)

# ---- 1. canonical -s 2 dream: up / down / universe ----
dream <- fread(file.path(INT, "dream_results.csv"))
dream[, ens := sub("\\..*", "", gene)]
up    <- dream[padj < PADJ & logFC >  LFC, ens]
down  <- dream[padj < PADJ & logFC < -LFC, ens]
univ  <- dream$ens
cat(sprintf("Up: %d  Down: %d  (ratio %.2f:1)  Universe: %d\n",
            length(up), length(down), length(up)/length(down), length(univ)))

# ---- 2. Hallmark ORA (msigdbr ensembl) ----
h <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
t2g <- unique(as.data.table(h)[, .(gs_name, ensembl_gene)])
run_ora <- function(genes) {
  r <- tryCatch(as.data.table(enricher(genes, TERM2GENE = t2g, universe = univ,
                  pvalueCutoff = 1.1, qvalueCutoff = 1.1)), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) return(data.table())
  r[order(p.adjust)][, .(ID, GeneRatio, pvalue, padj = p.adjust, Count)]
}
hm_up   <- run_ora(up);   hm_up[,   direction := "up"]
hm_down <- run_ora(down); hm_down[, direction := "down"]
hm <- rbind(hm_up, hm_down, fill = TRUE)
fwrite(hm, file.path(OUT, "hallmark_ora_up_vs_down.csv"))
cat("\n== Top Hallmark (UP) ==\n");   print(head(hm_up[,   .(ID, padj, Count)], 8))
cat("\n== Top Hallmark (DOWN) ==\n"); print(head(hm_down[, .(ID, padj, Count)], 8))

# ---- 3. curated liver cell-type marker enrichment (Fisher) ----
markers <- list(
  Hepatocyte = c("ALB","APOA1","APOB","CYP2E1","CYP3A4","CYP2A6","ASGR1","TTR","SERPINA1",
                 "HP","FGA","FGB","FGG","CPS1","ASS1","GLUL","HPX","APOC3","TF","AHSG"),
  Macrophage = c("CD68","CD163","MARCO","LYZ","CSF1R","ITGAM","CD14","C1QA","C1QB","C1QC","AIF1","VSIG4"),
  Stellate_Fibroblast = c("COL1A1","COL1A2","COL3A1","ACTA2","PDGFRB","DCN","LUM","TIMP1","TAGLN","PDGFRA","COL6A3","THY1"),
  Cholangiocyte = c("KRT19","KRT7","SOX9","EPCAM","CFTR","ANXA4","SPP1","HNF1B","CLDN4"),
  Endothelial = c("PECAM1","VWF","CLEC4G","CLEC4M","STAB2","LYVE1","CD34","ENG","KDR"),
  T_NK = c("CD3D","CD3E","CD2","TRAC","NKG7","GZMB","CCL5","IL7R","CD8A","GZMK")
)
fisher_dir <- function(genes, set_ens) {
  inset <- univ %in% set_ens
  indeg <- univ %in% genes
  tab <- table(factor(indeg, c(TRUE,FALSE)), factor(inset, c(TRUE,FALSE)))
  ft <- fisher.test(tab, alternative = "greater")
  data.table(n_overlap = sum(indeg & inset), set_size = sum(inset),
             OR = unname(ft$estimate), p = ft$p.value)
}
ct <- rbindlist(lapply(names(markers), function(ct_name) {
  set_ens <- unique(sym2ens[markers[[ct_name]]]); set_ens <- set_ens[!is.na(set_ens)]
  rbind(
    cbind(celltype = ct_name, direction = "up",   fisher_dir(up,   set_ens)),
    cbind(celltype = ct_name, direction = "down", fisher_dir(down, set_ens))
  )
}))
ct[, padj := p.adjust(p, "BH")]
fwrite(ct, file.path(OUT, "celltype_marker_enrichment_up_vs_down.csv"))
cat("\n== Cell-type marker enrichment (up vs down) ==\n")
print(ct[order(direction, -OR), .(celltype, direction, n_overlap, set_size, OR = round(OR,2), padj = signif(padj,2))])

# ---- 4. composition decomposition (deconv -s 2) ----
decf <- file.path(BASE, "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")
if (file.exists(decf)) {
  dec <- fread(decf)
  dec[, dir := ifelse(logFC_unadj > 0, "up", "down")]
  sig <- dec[category %in% c("Composition_driven","Hepatocyte_intrinsic","Unmasked")]
  cat("\n== Composition decomposition (deconv -s 2) ==\n")
  tab <- dcast(sig[, .N, by = .(category, dir)], category ~ dir, value.var = "N", fill = 0)
  tab[, ratio_up_down := round(up / pmax(down,1), 2)]
  print(tab)
  # bulk vs hepatocyte-intrinsic up:down among bulk-significant genes
  bs <- dec[sig_unadj == TRUE]
  bulk_ratio <- sum(bs$logFC_unadj > 0) / max(sum(bs$logFC_unadj < 0), 1)
  hi <- dec[category == "Hepatocyte_intrinsic"]
  hi_ratio <- sum(hi$logFC_adj > 0) / max(sum(hi$logFC_adj < 0), 1)
  cat(sprintf("\nUp:down ratio  BULK (sig_unadj): %.2f:1   HEPATOCYTE-INTRINSIC (logFC_adj): %.2f:1\n",
              bulk_ratio, hi_ratio))
  cat(sprintf("Composition_driven genes are %.0f%% up; Hepatocyte_intrinsic are %.0f%% up\n",
              100*mean(sig[category=="Composition_driven"]$dir=="up"),
              100*mean(sig[category=="Hepatocyte_intrinsic"]$dir=="up")))
  fwrite(tab, file.path(OUT, "composition_decomposition.csv"))
} else {
  cat("\nWARNING: deconv_attribution_scores.csv not found — Part B skipped.\n")
}
cat("\n[done] outputs in", OUT, "\n")
