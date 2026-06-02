# 243_multi_anchor_signature_concordance.R
# Phase 1.3 — Multi-anchor signature concordance for F3 sub-states.
# Li et al. F3a/F3b is treated as ONE hypothesis among several anchors
# (Reactome ECM, SenMayo senescence, NMF P5/P6, MASH bulk signatures).
#
# Inputs:
#   results/integration/merged_dge.rds, meta_matched.rds
#   results/granular_staging/f3_substate_per_cohort_labels.csv (from 241)
#   results/granular_staging/f3_substate_pooled_labels.csv (from 241)
#   results/subtypes/nmf_results_cache_clean.rds (k=6 NMF programs P5/P6 as anchors)
#
# Outputs (to results/granular_staging/):
#   f3_substate_signature_scores.csv  — per-sample ssGSEA scores per signature
#   f3_substate_signature_concordance.csv — cluster-mean + Wilcoxon + effect-size per signature
#   f3_substate_signature_summary.md  — human-readable summary

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(GSVA)
  library(msigdbr)
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
NMF_CACHE <- file.path(PROJECT_ROOT, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")

# Sanity: 241 must have produced labels
labels_path <- file.path(OUT_DIR, "f3_substate_pooled_labels.csv")
stopifnot("241_f3_substate_discovery.R must run first" = file.exists(labels_path))

cat("== Phase 1.3 Multi-anchor signature concordance ==\n")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]

pooled_labels <- read.csv(labels_path, stringsAsFactors = FALSE)
per_cohort_labels <- read.csv(file.path(OUT_DIR, "f3_substate_per_cohort_labels.csv"),
                              stringsAsFactors = FALSE)

f3_idx <- which(!is.na(meta$fibrosis_stage) & meta$fibrosis_stage == 3)
stopifnot(length(f3_idx) == nrow(pooled_labels))

# ---------------------------------------------------------------------------
# Build signature catalogue
# ---------------------------------------------------------------------------
cat("Building signature catalogue...\n")

# Convert ENSG to symbols using multi-evidence atlas as the source of truth.
atlas_sym <- read.csv(file.path(PROJECT_ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
                      stringsAsFactors = FALSE)[, c("ensembl_id", "human_symbol")]
atlas_sym$ensembl_clean <- sub("\\..*", "", atlas_sym$ensembl_id)
ensg_to_sym <- function(ensg) {
  ensg_clean <- sub("\\..*", "", ensg)
  m <- atlas_sym$human_symbol[match(ensg_clean, atlas_sym$ensembl_clean)]
  out <- ifelse(is.na(m) | m == "", ensg_clean, m)
  setNames(out, ensg)
}

all_genes_sym <- ensg_to_sym(rownames(dge))

# --- Anchor 1: Li et al. F3a (UPR/SREBP/cholesterol/lipid) ---
li_f3a_genes <- c(
  # UPR / ER stress
  "HSPA5", "DDIT3", "ATF4", "ATF6", "EIF2AK3", "ERN1", "XBP1", "DNAJB9",
  # SREBP / cholesterol biosynthesis
  "SREBF1", "SREBF2", "HMGCR", "HMGCS1", "SQLE", "DHCR7", "DHCR24",
  "MVD", "MVK", "FASN", "ACACA", "SCD", "INSIG1",
  # Lipid metabolism
  "FABP1", "ACOX1", "CPT1A", "PLIN2", "DGAT1"
)

# --- Anchor 2: Li et al. F3b (ECM/senescence/Ig) ---
li_f3b_genes <- c(
  # Top reported in Li et al. letter
  "IGFBP7", "BGN", "COL1A2", "COL3A1", "TIMP1",
  # Additional ECM
  "COL1A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
  "FN1", "VCAN", "DCN", "LUM", "LOX", "LOXL1", "LOXL2",
  "MMP2", "MMP9", "TIMP2", "ACTA2", "TAGLN",
  # Senescence-adjacent
  "CDKN1A", "CDKN2A", "GLB1", "SERPINE1"
)

# --- Anchor 3: Reactome ECM / Integrin / Collagen ---
msig_reactome <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  error = function(e) {
    # fallback to old API
    suppressWarnings(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
  }
)
react_terms <- c(
  "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_COLLAGEN_FORMATION",
  "REACTOME_COLLAGEN_DEGRADATION",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF",
  "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
  "REACTOME_FATTY_ACID_METABOLISM"
)
gene_col <- intersect(c("gene_symbol", "human_gene_symbol", "gs_symbol"), colnames(msig_reactome))[1]
set_col <- intersect(c("gs_name", "gs_id"), colnames(msig_reactome))[1]

reactome_sigs <- lapply(react_terms, function(term) {
  unique(msig_reactome[[gene_col]][msig_reactome[[set_col]] == term])
})
names(reactome_sigs) <- react_terms

# --- Anchor 4: Hallmark UPR / EMT ---
msig_hallmark <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) suppressWarnings(msigdbr(species = "Homo sapiens", category = "H"))
)
hallmark_terms <- c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "HALLMARK_FATTY_ACID_METABOLISM",
                    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                    "HALLMARK_INFLAMMATORY_RESPONSE",
                    "HALLMARK_CHOLESTEROL_HOMEOSTASIS")
hallmark_sigs <- lapply(hallmark_terms, function(term) {
  unique(msig_hallmark[[gene_col]][msig_hallmark[[set_col]] == term])
})
names(hallmark_sigs) <- hallmark_terms

# --- Anchor 5: SenMayo senescence ---
# SenMayo from Saul et al. 2022 Nat Commun. 125 genes.
senmayo_genes <- c(
  "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6",
  "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3",
  "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2",
  "CSF2RB", "CST10", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16",
  "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG",
  "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG",
  "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3",
  "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18",
  "IL1A", "IL1B", "IL2", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2",
  "ITPKA", "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13",
  "MMP14", "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1",
  "PGF", "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES",
  "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB2", "SERPINB3A",
  "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C",
  "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC",
  "VGF", "WNT16", "WNT2"
)

# --- Anchor 6: Govaere/Hardy MASH bulk signatures (use existing dream DEGs) ---
dream_path <- file.path(INT_DIR, "dream_results_ashr.csv")
dream_sigs <- list()
if (file.exists(dream_path)) {
  dr <- read.csv(dream_path, stringsAsFactors = FALSE)
  if ("gene_symbol" %in% colnames(dr) || "symbol" %in% colnames(dr)) {
    sym_col <- intersect(c("gene_symbol", "symbol"), colnames(dr))[1]
    lfc_col <- intersect(c("logFC", "logFC_disease", "shrunk_logFC"), colnames(dr))[1]
    p_col <- intersect(c("padj", "lfsr", "adj.P.Val"), colnames(dr))[1]
    if (!is.null(sym_col) && !is.null(lfc_col) && !is.null(p_col)) {
      dr$lfc <- dr[[lfc_col]]; dr$padj <- dr[[p_col]]
      up <- dr[[sym_col]][!is.na(dr$padj) & dr$padj < 0.05 & dr$lfc > 0.5]
      dn <- dr[[sym_col]][!is.na(dr$padj) & dr$padj < 0.05 & dr$lfc < -0.5]
      dream_sigs <- list(
        DREAM_DISEASE_UP = up[!is.na(up) & up != ""],
        DREAM_DISEASE_DN = dn[!is.na(dn) & dn != ""]
      )
    }
  }
}

# --- Anchor 7: NMF k=6 P5/P6 top loadings ---
nmf_sigs <- list()
if (file.exists(NMF_CACHE)) {
  nmf_obj <- tryCatch(readRDS(NMF_CACHE), error = function(e) NULL)
  if (!is.null(nmf_obj)) {
    # nmf_results_cache_clean.rds should be a list keyed by k
    if (is.list(nmf_obj) && !is.null(nmf_obj[["6"]])) {
      W <- tryCatch(NMF::basis(nmf_obj[["6"]]$res), error = function(e) NULL)
      if (!is.null(W) && nrow(W) > 0) {
        # P5 and P6 (Hepatic-metabolic, Stellate-myofibroblast)
        for (p in c("P5", "P6")) {
          col_idx <- as.integer(sub("P", "", p))
          if (col_idx <= ncol(W)) {
            top_genes <- rownames(W)[order(W[, col_idx], decreasing = TRUE)[1:200]]
            # If genes in W are ENSG, convert
            if (any(grepl("^ENSG", top_genes))) {
              top_sym <- ensg_to_sym(top_genes)
              nmf_sigs[[paste0("NMF_K6_", p, "_top200")]] <- as.character(top_sym)
            } else {
              nmf_sigs[[paste0("NMF_K6_", p, "_top200")]] <- top_genes
            }
          }
        }
      }
    }
  }
}

# Aggregate all signatures
all_sigs <- c(
  list(LI_F3A_HYPOTHESIS = li_f3a_genes,
       LI_F3B_HYPOTHESIS = li_f3b_genes,
       SENMAYO_SENESCENCE = senmayo_genes),
  reactome_sigs,
  hallmark_sigs,
  dream_sigs,
  nmf_sigs
)
all_sigs <- lapply(all_sigs, function(g) g[!is.na(g) & g != ""])
all_sigs <- all_sigs[lengths(all_sigs) >= 5]

cat(sprintf("Anchor signatures: %d (sizes %d-%d)\n",
            length(all_sigs), min(lengths(all_sigs)), max(lengths(all_sigs))))
print(data.frame(signature = names(all_sigs), n_genes = lengths(all_sigs)))

# ---------------------------------------------------------------------------
# Compute ssGSEA on F3 samples (batch-corrected log-CPM)
# ---------------------------------------------------------------------------
cat("\nComputing ssGSEA scores for F3 samples...\n")
dge_f3 <- dge[, f3_idx]
keep <- filterByExpr(dge_f3, group = factor(rep("F3", ncol(dge_f3))), min.count = 5)
dge_f3 <- dge_f3[keep, , keep.lib.sizes = FALSE]
dge_f3 <- calcNormFactors(dge_f3)
v_f3 <- voom(dge_f3, design = NULL)
expr_f3 <- v_f3$E
expr_f3_bc <- removeBatchEffect(expr_f3, batch = factor(meta$dataset[f3_idx]))
rownames(expr_f3_bc) <- ensg_to_sym(rownames(expr_f3_bc))

# Drop duplicated symbols (keep first)
expr_f3_bc <- expr_f3_bc[!duplicated(rownames(expr_f3_bc)), ]

gsva_param <- ssgseaParam(expr_f3_bc, all_sigs, normalize = TRUE)
ssgsea_scores <- gsva(gsva_param)
cat(sprintf("ssGSEA scores: %d signatures × %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))

write.csv(t(ssgsea_scores), file.path(OUT_DIR, "f3_substate_signature_scores.csv"))

# ---------------------------------------------------------------------------
# Cluster-level concordance (POOLED labels)
# ---------------------------------------------------------------------------
cat("\nComputing cluster-level concordance...\n")
samples <- colnames(ssgsea_scores)
labels_match <- pooled_labels$cluster_pooled[match(samples, pooled_labels$sample_id)]
stopifnot(all(!is.na(labels_match)))

unique_clusters <- sort(unique(labels_match))
cat(sprintf("Unique pooled clusters: %s\n", paste(unique_clusters, collapse = " / ")))
cat(sprintf("ssGSEA matrix dim: %d sigs × %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))
if (length(unique_clusters) < 2 || nrow(ssgsea_scores) < 1) {
  cat("Insufficient clusters or signatures — writing empty concordance and exiting cleanly.\n")
  empty <- data.frame(signature = character(), cluster_a = character(), cluster_b = character(),
                       n_a = integer(), n_b = integer(), mean_a = numeric(), mean_b = numeric(),
                       wilcox_p = numeric(), cohens_d = numeric(), mw_r_effect = numeric(),
                       wilcox_padj = numeric(), stringsAsFactors = FALSE)
  write.csv(empty, file.path(OUT_DIR, "f3_substate_signature_concordance.csv"), row.names = FALSE)
  sink(file.path(OUT_DIR, "f3_substate_signature_summary.md"))
  cat("# Phase 1.3 — Multi-anchor signature concordance\n\n")
  cat("Insufficient signatures or clusters; concordance not computed.\n")
  sink()
  cat("== Phase 1.3 complete (empty). ==\n")
  quit(status = 0)
}
concord_rows <- list()
for (sig_name in rownames(ssgsea_scores)) {
  scores_v <- ssgsea_scores[sig_name, ]
  for (i in seq_along(unique_clusters)[-length(unique_clusters)]) {
    for (j in (i+1):length(unique_clusters)) {
      c1 <- unique_clusters[i]; c2 <- unique_clusters[j]
      a <- scores_v[labels_match == c1]
      b <- scores_v[labels_match == c2]
      if (length(a) < 3 || length(b) < 3) next
      wt <- suppressWarnings(wilcox.test(a, b))
      cohens_d <- (mean(a) - mean(b)) / sqrt((var(a) * (length(a)-1) + var(b) * (length(b)-1)) /
                                              (length(a) + length(b) - 2))
      mw_r <- abs(qnorm(wt$p.value / 2)) / sqrt(length(a) + length(b))
      concord_rows[[paste(sig_name, c1, c2, sep = "_")]] <- data.frame(
        signature = sig_name,
        cluster_a = c1, cluster_b = c2,
        n_a = length(a), n_b = length(b),
        mean_a = mean(a), mean_b = mean(b),
        wilcox_p = wt$p.value,
        cohens_d = cohens_d,
        mw_r_effect = mw_r,
        stringsAsFactors = FALSE
      )
    }
  }
}
concord <- if (length(concord_rows) > 0) do.call(rbind, concord_rows) else NULL
if (is.null(concord) || nrow(concord) == 0) {
  cat("No concordance rows produced (likely due to imbalanced cluster sizes — skipping summary).\n")
  empty <- data.frame(signature = character(), cluster_a = character(), cluster_b = character(),
                       n_a = integer(), n_b = integer(), mean_a = numeric(), mean_b = numeric(),
                       wilcox_p = numeric(), cohens_d = numeric(), mw_r_effect = numeric(),
                       wilcox_padj = numeric(), stringsAsFactors = FALSE)
  write.csv(empty, file.path(OUT_DIR, "f3_substate_signature_concordance.csv"), row.names = FALSE)
  sink(file.path(OUT_DIR, "f3_substate_signature_summary.md"))
  cat("# Phase 1.3 — Multi-anchor signature concordance\n\n")
  cat(sprintf("Concordance NOT computed: pooled clusters too imbalanced (sizes %s).\n",
              paste(sapply(unique_clusters, function(c) sum(labels_match == c)), collapse = ", ")))
  sink()
  cat("== Phase 1.3 complete (empty due to imbalanced clusters). ==\n")
  quit(status = 0)
}
concord$wilcox_padj <- p.adjust(concord$wilcox_p, method = "BH")
write.csv(concord, file.path(OUT_DIR, "f3_substate_signature_concordance.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Summary markdown
# ---------------------------------------------------------------------------
top_sig <- concord[order(-abs(concord$cohens_d)), ][1:15, ]
sink(file.path(OUT_DIR, "f3_substate_signature_summary.md"))
cat("# Phase 1.3 — Multi-anchor signature concordance\n\n")
cat(sprintf("- F3 samples scored: %d\n", ncol(ssgsea_scores)))
cat(sprintf("- Pooled clusters: %s\n", paste(unique_clusters, collapse = " / ")))
cat(sprintf("- Anchor signatures: %d\n\n", nrow(ssgsea_scores)))
cat("## Top 15 cluster-discriminating signatures (by |Cohen's d|)\n\n")
cat("| Signature | Clusters | Cohen's d | Wilcox p | BH-padj |\n")
cat("|---|---|---:|---:|---:|\n")
for (i in seq_len(nrow(top_sig))) {
  cat(sprintf("| %s | %s vs %s | %.2f | %.2e | %.2e |\n",
              top_sig$signature[i], top_sig$cluster_a[i], top_sig$cluster_b[i],
              top_sig$cohens_d[i], top_sig$wilcox_p[i], top_sig$wilcox_padj[i]))
}
cat("\n## Li et al. anchor concordance\n\n")
li <- concord[concord$signature %in% c("LI_F3A_HYPOTHESIS", "LI_F3B_HYPOTHESIS"), ]
print(li[, c("signature", "cluster_a", "cluster_b", "cohens_d", "wilcox_padj")], row.names = FALSE)
cat("\n## NMF P5/P6 anchor concordance\n\n")
nmf_c <- concord[grepl("^NMF_K6", concord$signature), ]
if (nrow(nmf_c) > 0) print(nmf_c[, c("signature", "cluster_a", "cluster_b", "cohens_d", "wilcox_padj")], row.names = FALSE)
sink()

cat(sprintf("\nWrote %d concordance rows to %s\n", nrow(concord),
            file.path(OUT_DIR, "f3_substate_signature_concordance.csv")))
cat("== Phase 1.3 complete. ==\n")
