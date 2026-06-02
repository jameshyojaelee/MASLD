#!/usr/bin/env Rscript
# 86_cellstate_signatures_bulk.R
#
# Analyses H2 + J1 + J2 + K1 (v1) — Cell-state signature scoring in bulk.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Scores curated cell-state signatures in bulk dream t-stat ranking via fgsea:
#   H2 Senescence: SenMayo (Saul et al.); core senescence (CDKN2A, CDKN1A, TP53, LMNB1 loss)
#   H2 SASP subtypes: pro-fibrotic, pro-inflammatory, growth-arrest
#   J1 T-cell exhaustion: PDCD1/HAVCR2/LAG3/TOX/TIGIT/CTLA4/BTLA/CD244
#   J1 NK dysfunction: GZMB/PRF1/NKG7/IFNG loss markers
#   J2 Ductular reaction: HNF1B/KRT19/KRT7/SOX9/EpCAM; HNF4A-downregulated Hep
#   K1 LSEC capillarization: VWF/CD34/PLVAP gain; LYVE1/STAB2/FCGR2B loss
#   K1 Angiocrine signals: HGF, WNT2, RSPO3
#
# Env: rnaseq
# Outputs: RNA-seq/results/celltype_attribution/

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR  <- file.path(BASE, "RNA-seq/results/celltype_attribution")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading bulk dream DE...")
bulk <- fread(file.path(INT_RES, "dream_results_ashr.csv"),
              select = c("gene","symbol","logFC","t","padj"))
bulk <- bulk[!is.na(symbol) & symbol != ""]
bulk <- bulk[order(-abs(t))][!duplicated(symbol)]
ranks <- setNames(bulk$t, bulk$symbol)
ranks <- ranks[!is.na(ranks)]

# ---------------------------------------------------------------------------
# Signatures
# ---------------------------------------------------------------------------
signatures <- list(
  # H2 Senescence
  "H2_senescence_core_up"   = c("CDKN2A","CDKN1A","TP53","GDF15","SERPINE1","IGFBP3",
                                 "SERPINB2","MMP3","MMP10","MMP12","GLB1","TNFRSF10D"),
  "H2_senescence_core_down" = c("LMNB1","MKI67","HMGB1","BIRC5"),
  "H2_SASP_profibrotic"     = c("TIMP1","SERPINE1","TGFB1","TGFB2","CCN2","CTGF","PAI1",
                                 "IL11","IGFBP3","IGFBP6","SPP1"),
  "H2_SASP_proinflammatory" = c("IL6","IL8","CXCL1","CXCL2","CXCL3","CXCL8",
                                 "IL1A","IL1B","TNF","CCL2","CCL20"),
  "H2_SASP_growth_arrest"   = c("SPP1","IGFBP3","IGFBP6","IGFBP7","CDKN1A","CDKN2A",
                                 "TNFRSF10D","TNFRSF10B"),
  "H2_SenMayo_subset"       = c("ACVR1B","ANG","AREG","AXL","BEX3","BMP2","BMP6","C3",
                                 "CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26",
                                 "CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8","CD55","CD9",
                                 "CSF1","CSF2","CSF2RB","CST10","CTNNB1","CTSB","CXCL1",
                                 "CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8","CXCR2",
                                 "DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2","FAS",
                                 "FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1",
                                 "ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4",
                                 "IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18",
                                 "IL1A","IL1B","IL2","IL6","IL6ST","IL7","INHA","IQGAP2",
                                 "ITGA2","ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10",
                                 "MMP12","MMP13","MMP14","MMP2","MMP3","MMP9","NAP1L4",
                                 "NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU","PLAUR",
                                 "PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F",
                                 "SERPINB3","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX",
                                 "TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B",
                                 "TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2"),

  # J1 T / NK dysfunction
  "J1_Tcell_exhaustion"     = c("PDCD1","HAVCR2","LAG3","TOX","TIGIT","CTLA4","BTLA",
                                 "CD244","CXCL13","LAYN","ENTPD1","NR4A1","NR4A2","NR4A3",
                                 "BATF","EOMES","IRF4","TNFRSF9"),
  "J1_NK_cytotoxicity"      = c("GZMA","GZMB","GZMH","PRF1","NKG7","KLRD1","KLRF1","NCR1",
                                 "NCR3","FCGR3A","GNLY","IFNG","TNF","FASLG"),
  "J1_NK_inhibitory_down"   = c("KLRC1","KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2",
                                 "NCR2","CD96"),

  # J2 Ductular reaction
  "J2_ductular_reaction"    = c("HNF1B","KRT19","KRT7","SOX9","EPCAM","PROM1","ANXA4",
                                 "CD24","SPP1","OPN","LGR5","AFP","DLK1","EpCAM","NOTCH2"),
  "J2_cholangiocyte"        = c("KRT7","KRT8","KRT18","KRT19","SOX9","HNF1B","EPCAM",
                                 "JAG1","CFTR","SLC4A4","TFF1","MUC5B","SCTR"),

  # K1 LSEC capillarization
  "K1_LSEC_capillarization_up" = c("VWF","CD34","PLVAP","PECAM1","CDH5","EFNB2"),
  "K1_LSEC_healthy_markers_down" = c("LYVE1","STAB2","STAB1","FCGR2B","CLEC4G","CLEC4M",
                                      "CLEC1B","CD32B"),
  "K1_angiocrine_signals"   = c("HGF","WNT2","RSPO3","ANGPT1","ANGPT2","VEGFA","DLL4",
                                 "JAG1","JAG2","BMP2","BMP6")
)

message("[2] fgsea scoring of cell-state signatures in bulk dream t-stat...")
# Limit to signatures with enough genes present in ranks
min_size <- 3
sigs_filtered <- lapply(signatures, function(s) intersect(s, names(ranks)))
sigs_filtered <- sigs_filtered[sapply(sigs_filtered, length) >= min_size]
message(sprintf("  Signatures passing min_size=%d: %d / %d",
                min_size, length(sigs_filtered), length(signatures)))

set.seed(42)
fr <- fgsea(pathways = sigs_filtered, stats = ranks, eps = 0,
            nPermSimple = 10000)
fr[, direction := fifelse(NES > 0, "MASLD_up", "MASLD_down")]
fr[, theme := sub("_.*", "", pathway)]
fr[, abs_NES := abs(NES)]
setorder(fr, -abs_NES)
fr[, leading_edge_str := sapply(leadingEdge, paste, collapse = ",")]

message("[3] Per-signature gene-level view (bulk DE for each gene)...")
gene_lvl <- rbindlist(lapply(names(sigs_filtered), function(s) {
  g <- intersect(signatures[[s]], bulk$symbol)
  b <- bulk[symbol %in% g, .(symbol, bulk_lfc = logFC, bulk_t = t, bulk_padj = padj)]
  b[, signature := s]
  b
}), fill = TRUE)

fwrite(fr[, !"leadingEdge"], file.path(OUTDIR, "cellstate_signatures_fgsea.csv"))
fwrite(gene_lvl,             file.path(OUTDIR, "cellstate_signatures_gene_level.csv"))

summary_lines <- c(
  "=== fgsea of cell-state signatures in bulk dream t-stat ranking ===",
  capture.output(print(fr[, .(pathway, theme, size, NES, pval, padj, direction)],
                       nrows = 30)),
  "",
  "=== Gene-level bulk DE for key signature genes (top 5 up + top 5 down per sig) ===",
  "",
  capture.output({
    for (s in unique(gene_lvl$signature)) {
      sub <- gene_lvl[signature == s][order(-bulk_t)]
      cat(sprintf("\n--- %s ---\n", s))
      print(rbind(head(sub, 5), tail(sub, 5)))
    }
  })
)
writeLines(summary_lines, file.path(OUTDIR, "cellstate_signatures_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
