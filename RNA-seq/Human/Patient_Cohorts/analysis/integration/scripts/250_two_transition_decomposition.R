# 250_two_transition_decomposition.R
# Two-transition: Test whether F1→F2 and F2→F3 activate distinct molecular programs.
#
# Per Phase 1 finding: F1→F2 and F2→F3 both have substantial mean shifts
# (median |LFC| top1000: 0.245 and 0.400 respectively). Tests whether consecutive
# transitions activate distinct molecular programs across a multi-step cascade.
# NOTE: The "metabolic-to-inflammatory switch at F2" framing was retired 2026-05-09;
# current thesis is multi-step cellular cascade across 4 CRN transitions (see paper_outline.md).
#
# Tests:
#   (1) Pathway-level decomposition — Hallmark + Reactome enrichment per transition;
#       which pathways are unique to F1→F2 vs unique to F2→F3 vs shared.
#   (2) NMF k=6 program activity — mean P1..P6 score at F0/F1/F2/F3/F4 stages.
#   (3) Cell-type proportion shifts per stage — using existing deconvolution.
#   (4) Sex-stratified effect per transition.
#
# Inputs:
#   results/integration/merged_dge.rds, meta_matched.rds
#   results/celltype_attribution/persample_celltype_proportions.csv
#   results/subtypes/nmf_results_cache_clean.rds
#   results/multi_evidence/multi_evidence_atlas.csv (for ENSG → symbol)
#
# Outputs (results/granular_staging/):
#   two_transition_pathway_decomposition.csv   — pathway × transition NES table
#   two_transition_pathway_unique.md           — unique up/down pathways per transition
#   two_transition_nmf_program_by_stage.csv    — k=6 NMF program scores per stage
#   two_transition_celltype_by_stage.csv       — cell-type proportion shifts per transition
#   two_transition_sex_stratified.csv          — sex × transition effect-size table
#   two_transition_summary.md                   — narrative summary

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
NMF_CACHE <- file.path(PROJECT_ROOT, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")
DECONV_PATH <- file.path(PROJECT_ROOT, "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv")
ATLAS_PATH <- file.path(PROJECT_ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

cat("== Two-transition program decomposition ==\n")

# Symbol mapping
atlas_sym <- read.csv(ATLAS_PATH, stringsAsFactors = FALSE)[, c("ensembl_id", "human_symbol")]
atlas_sym$ensembl_clean <- sub("\\..*", "", atlas_sym$ensembl_id)
ensg_to_sym <- function(ensg) {
  ensg_clean <- sub("\\..*", "", ensg)
  m <- atlas_sym$human_symbol[match(ensg_clean, atlas_sym$ensembl_clean)]
  ifelse(is.na(m) | m == "", ensg_clean, m)
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]

cat(sprintf("Samples: %d\n", ncol(dge)))
cat("Stage distribution:\n"); print(table(meta$fibrosis_stage, useNA = "ifany"))

# ---------------------------------------------------------------------------
# (1) Pathway decomposition per transition
# ---------------------------------------------------------------------------
cat("\n--- (1) Pathway decomposition per transition ---\n")

run_contrast <- function(idx_a, idx_b) {
  group <- factor(c(rep("A", length(idx_a)), rep("B", length(idx_b))), levels = c("A", "B"))
  cohort <- factor(meta$dataset[c(idx_a, idx_b)])
  dge_sub <- dge[, c(idx_a, idx_b)]
  keep <- filterByExpr(dge_sub, group = group, min.count = 5)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub)
  design <- if (length(unique(cohort)) > 1) model.matrix(~ group + cohort) else model.matrix(~ group)
  v <- voom(dge_sub, design)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = "groupB", number = Inf, sort.by = "none")
  res$gene_id <- rownames(res)
  res$gene_symbol <- ensg_to_sym(res$gene_id)
  res
}

f0 <- which(meta$fibrosis_stage == 0)
f1 <- which(meta$fibrosis_stage == 1)
f2 <- which(meta$fibrosis_stage == 2)
f3 <- which(meta$fibrosis_stage == 3)
f4 <- which(meta$fibrosis_stage == 4)

cat(sprintf("F0: %d  F1: %d  F2: %d  F3: %d  F4: %d\n",
            length(f0), length(f1), length(f2), length(f3), length(f4)))

cat("\nComputing F0→F1...\n"); F0F1 <- run_contrast(f0, f1)
cat("Computing F1→F2...\n");   F1F2 <- run_contrast(f1, f2)
cat("Computing F2→F3...\n");   F2F3 <- run_contrast(f2, f3)
cat("Computing F3→F4...\n");   F3F4 <- run_contrast(f3, f4)

transitions <- list(F0F1 = F0F1, F1F2 = F1F2, F2F3 = F2F3, F3F4 = F3F4)

# Hallmark + Reactome
msig_h <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) suppressWarnings(msigdbr(species = "Homo sapiens", category = "H"))
)
msig_r <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  error = function(e) suppressWarnings(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
)
gene_col <- intersect(c("gene_symbol", "human_gene_symbol", "gs_symbol"), colnames(msig_h))[1]
set_col <- intersect(c("gs_name", "gs_id"), colnames(msig_h))[1]

hallmark_sets <- split(msig_h[[gene_col]], msig_h[[set_col]])
hallmark_sets <- lapply(hallmark_sets, unique)
# Top Reactome by relevance — use a curated list focused on transition biology
react_terms <- c(
  "REACTOME_ECM_PROTEOGLYCANS", "REACTOME_COLLAGEN_FORMATION", "REACTOME_COLLAGEN_DEGRADATION",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX1",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF",
  "REACTOME_FATTY_ACID_METABOLISM", "REACTOME_BETA_OXIDATION_OF_OCTANOYL_COA_TO_HEXANOYL_COA",
  "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
  "REACTOME_NEUTROPHIL_DEGRANULATION",
  "REACTOME_COMPLEMENT_CASCADE",
  "REACTOME_BILE_ACID_AND_BILE_SALT_METABOLISM"
)
react_sets <- lapply(react_terms, function(term) {
  unique(msig_r[[gene_col]][msig_r[[set_col]] == term])
})
names(react_sets) <- react_terms
react_sets <- react_sets[lengths(react_sets) >= 5]

all_pathways <- c(hallmark_sets, react_sets)
cat(sprintf("Pathway sets: %d Hallmark + %d Reactome = %d total\n",
            length(hallmark_sets), length(react_sets), length(all_pathways)))

run_fgsea <- function(res, label) {
  ranks <- res$logFC
  names(ranks) <- res$gene_symbol
  ranks <- ranks[!duplicated(names(ranks)) & !is.na(names(ranks)) & names(ranks) != ""]
  ranks <- sort(ranks, decreasing = TRUE)
  fgr <- suppressWarnings(fgsea(pathways = all_pathways, stats = ranks,
                                minSize = 10, maxSize = 500, nproc = 4))
  fgr$transition <- label
  fgr
}

pwy_rows <- lapply(names(transitions), function(n) run_fgsea(transitions[[n]], n))
pwy_df <- do.call(rbind, pwy_rows)
# Compact representation
pwy_compact <- data.frame(
  pathway = pwy_df$pathway,
  transition = pwy_df$transition,
  NES = pwy_df$NES,
  padj = pwy_df$padj,
  stringsAsFactors = FALSE
)
write.csv(pwy_compact, file.path(OUT_DIR, "two_transition_pathway_decomposition.csv"), row.names = FALSE)
cat(sprintf("Wrote pathway × transition: %d rows\n", nrow(pwy_compact)))

# Identify pathways UP / DOWN at each transition
sig_at <- function(transition, direction = "up") {
  pw <- pwy_compact[pwy_compact$transition == transition & pwy_compact$padj < 0.05, ]
  if (direction == "up") pw <- pw[pw$NES > 0, ] else pw <- pw[pw$NES < 0, ]
  pw[order(-abs(pw$NES)), ]
}

f1f2_up <- sig_at("F1F2", "up")
f2f3_up <- sig_at("F2F3", "up")
f1f2_dn <- sig_at("F1F2", "down")
f2f3_dn <- sig_at("F2F3", "down")

unique_f1f2_up <- setdiff(f1f2_up$pathway, f2f3_up$pathway)
unique_f2f3_up <- setdiff(f2f3_up$pathway, f1f2_up$pathway)
shared_up <- intersect(f1f2_up$pathway, f2f3_up$pathway)

cat(sprintf("Pathways UP at F1→F2 only: %d\n", length(unique_f1f2_up)))
cat(sprintf("Pathways UP at F2→F3 only: %d\n", length(unique_f2f3_up)))
cat(sprintf("Pathways UP at BOTH:        %d\n", length(shared_up)))

# Markdown summary
sink(file.path(OUT_DIR, "two_transition_pathway_unique.md"))
cat("# Two-transition: Pathways unique to each transition\n\n")
cat(sprintf("- F0→F1: %d significant (padj<0.05)\n", sum(pwy_compact$transition=="F0F1" & pwy_compact$padj<0.05)))
cat(sprintf("- F1→F2: %d significant\n", sum(pwy_compact$transition=="F1F2" & pwy_compact$padj<0.05)))
cat(sprintf("- F2→F3: %d significant\n", sum(pwy_compact$transition=="F2F3" & pwy_compact$padj<0.05)))
cat(sprintf("- F3→F4: %d significant\n\n", sum(pwy_compact$transition=="F3F4" & pwy_compact$padj<0.05)))

cat("## Pathways UP at F1→F2 ONLY (not F2→F3)\n\n")
for (i in seq_len(min(20, length(unique_f1f2_up)))) {
  pw <- unique_f1f2_up[i]
  nes <- f1f2_up$NES[f1f2_up$pathway == pw]
  cat(sprintf("- **%s** (NES = %.2f)\n", pw, nes))
}
cat("\n## Pathways UP at F2→F3 ONLY (not F1→F2)\n\n")
for (i in seq_len(min(20, length(unique_f2f3_up)))) {
  pw <- unique_f2f3_up[i]
  nes <- f2f3_up$NES[f2f3_up$pathway == pw]
  cat(sprintf("- **%s** (NES = %.2f)\n", pw, nes))
}
cat("\n## Pathways UP at BOTH transitions (shared)\n\n")
for (i in seq_len(min(20, length(shared_up)))) {
  pw <- shared_up[i]
  nes_a <- f1f2_up$NES[f1f2_up$pathway == pw]
  nes_b <- f2f3_up$NES[f2f3_up$pathway == pw]
  cat(sprintf("- **%s** (F1→F2 NES = %.2f, F2→F3 NES = %.2f)\n", pw, nes_a, nes_b))
}
sink()

# ---------------------------------------------------------------------------
# (2) NMF k=6 program activity per stage
# ---------------------------------------------------------------------------
cat("\n--- (2) NMF k=6 program activity per stage ---\n")
nmf_obj <- tryCatch(readRDS(NMF_CACHE), error = function(e) NULL)

nmf_per_stage <- NULL
if (!is.null(nmf_obj)) {
  if (is.list(nmf_obj) && !is.null(nmf_obj[["6"]])) {
    obj6 <- nmf_obj[["6"]]
    H <- if (!is.null(obj6$res)) {
      tryCatch(NMF::coef(obj6$res), error = function(e) NULL)
    } else NULL
    if (!is.null(H)) {
      cat(sprintf("NMF H matrix: %d programs × %d samples\n", nrow(H), ncol(H)))
      cat("First 5 sample names in H:", paste(head(colnames(H), 5), collapse = ", "), "\n")
      # Match samples
      common_samples <- intersect(colnames(H), meta$sample_id)
      cat(sprintf("Samples in NMF cache matched to current metadata: %d / %d\n",
                  length(common_samples), ncol(H)))
      meta_match <- meta[match(common_samples, meta$sample_id), ]
      H_match <- H[, common_samples, drop = FALSE]

      # Average program score per stage
      stages <- 0:4
      prog_names <- paste0("P", seq_len(nrow(H_match)))
      nmf_rows <- list()
      for (s in stages) {
        s_idx <- which(meta_match$fibrosis_stage == s)
        if (length(s_idx) >= 5) {
          for (p in seq_len(nrow(H_match))) {
            score <- H_match[p, s_idx]
            nmf_rows[[paste(s, p)]] <- data.frame(
              stage = s, program = prog_names[p],
              n_samples = length(s_idx),
              mean_score = mean(score),
              median_score = median(score),
              sd_score = sd(score),
              stringsAsFactors = FALSE
            )
          }
        }
      }
      nmf_per_stage <- do.call(rbind, nmf_rows)
      write.csv(nmf_per_stage, file.path(OUT_DIR, "two_transition_nmf_program_by_stage.csv"),
                row.names = FALSE)
      cat("\nNMF program × stage:\n")
      print(reshape(nmf_per_stage[, c("stage", "program", "mean_score")],
                    idvar = "program", timevar = "stage", direction = "wide"))
    }
  }
}

# ---------------------------------------------------------------------------
# (3) Cell-type proportion shifts per stage
# ---------------------------------------------------------------------------
cat("\n--- (3) Cell-type proportion shifts per stage ---\n")
deconv <- read.csv(DECONV_PATH, stringsAsFactors = FALSE)
ct_cols <- setdiff(colnames(deconv),
                   c("sample_id", "dataset", "group_binary", "condition",
                     "fibrosis_stage", "inferred_sex", "diagnosis_harmonized",
                     "f2_group"))
cat(sprintf("Cell-type columns: %d\n", length(ct_cols)))
cat(sprintf("Samples with deconv: %d\n", nrow(deconv)))

ct_shifts <- list()
for (ct in ct_cols) {
  for (transition in c("F0vF1", "F1vF2", "F2vF3", "F3vF4")) {
    s_a <- as.integer(sub("F", "", strsplit(transition, "v")[[1]][1]))
    s_b <- as.integer(sub("F", "", strsplit(transition, "v")[[1]][2]))
    a <- deconv[[ct]][deconv$fibrosis_stage == s_a]
    b <- deconv[[ct]][deconv$fibrosis_stage == s_b]
    a <- a[!is.na(a) & is.finite(a)]; b <- b[!is.na(b) & is.finite(b)]
    if (length(a) < 5 || length(b) < 5) next
    wt <- suppressWarnings(wilcox.test(a, b))
    ct_shifts[[paste(ct, transition)]] <- data.frame(
      celltype = ct, transition = transition,
      n_a = length(a), n_b = length(b),
      mean_a = mean(a), mean_b = mean(b),
      log2FC_proportion = log2((mean(b) + 1e-6) / (mean(a) + 1e-6)),
      wilcox_p = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}
ct_df <- do.call(rbind, ct_shifts)
ct_df$padj <- p.adjust(ct_df$wilcox_p, method = "BH")
write.csv(ct_df, file.path(OUT_DIR, "two_transition_celltype_by_stage.csv"), row.names = FALSE)

# Top cell-type shifts per transition
cat("\nTop cell-type shifts per transition (BH-significant only):\n")
for (transition in c("F0vF1", "F1vF2", "F2vF3", "F3vF4")) {
  sig <- ct_df[ct_df$transition == transition & ct_df$padj < 0.05, ]
  sig <- sig[order(-abs(sig$log2FC_proportion)), ]
  if (nrow(sig) > 0) {
    cat(sprintf("\n%s (%d significant):\n", transition, nrow(sig)))
    for (i in seq_len(min(8, nrow(sig)))) {
      cat(sprintf("  %-30s log2FC=%6.2f  p=%.2e\n",
                  sig$celltype[i], sig$log2FC_proportion[i], sig$padj[i]))
    }
  }
}

# ---------------------------------------------------------------------------
# (4) Sex-stratified effect per transition
# ---------------------------------------------------------------------------
cat("\n--- (4) Sex-stratified effect per transition ---\n")

run_sex_contrast <- function(idx_a, idx_b, sex) {
  m_a <- intersect(idx_a, which(meta$inferred_sex == sex))
  m_b <- intersect(idx_b, which(meta$inferred_sex == sex))
  if (length(m_a) < 5 || length(m_b) < 5) {
    return(list(effect = NA_real_, n_a = length(m_a), n_b = length(m_b),
                n_padj_05 = NA_integer_))
  }
  res <- run_contrast(m_a, m_b)
  res_o <- res[order(res$P.Value), ][1:1000, ]
  list(effect = median(abs(res_o$logFC), na.rm = TRUE),
       n_a = length(m_a), n_b = length(m_b),
       n_padj_05 = sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05))
}

sex_rows <- list()
for (transition in c("F1F2", "F2F3", "F3F4")) {
  parts <- regmatches(transition, regexpr("F\\d", transition, perl = TRUE))
  pair <- regmatches(transition, gregexpr("F\\d", transition, perl = TRUE))[[1]]
  sa <- as.integer(sub("F", "", pair[1]))
  sb <- as.integer(sub("F", "", pair[2]))
  idx_a <- which(meta$fibrosis_stage == sa)
  idx_b <- which(meta$fibrosis_stage == sb)
  for (sex in c("F", "M")) {
    res <- run_sex_contrast(idx_a, idx_b, sex)
    sex_rows[[paste(transition, sex)]] <- data.frame(
      transition = transition, sex = sex,
      n_a = res$n_a, n_b = res$n_b,
      median_abs_logFC_top1000 = res$effect,
      n_padj_05 = res$n_padj_05,
      stringsAsFactors = FALSE
    )
  }
}
sex_df <- do.call(rbind, sex_rows)
write.csv(sex_df, file.path(OUT_DIR, "two_transition_sex_stratified.csv"), row.names = FALSE)
cat("\nSex × transition effect sizes:\n"); print(sex_df)

# ---------------------------------------------------------------------------
# (5) Summary
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "two_transition_summary.md"))
cat("# Two-transition program decomposition\n\n")
cat("## Effect-size context\n\n")
cat("Two transitions both show substantial mean shifts (median |log2FC| top 1000 DEGs):\n")
for (transition in c("F0F1", "F1F2", "F2F3", "F3F4")) {
  res <- transitions[[transition]]
  res_o <- res[order(res$P.Value), ][1:1000, ]
  cat(sprintf("- %s: %.3f (%d padj<0.05 DEGs)\n",
              transition, median(abs(res_o$logFC), na.rm = TRUE),
              sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05)))
}

cat("\n## Pathway-level distinctness\n\n")
cat(sprintf("- Pathways UP F1→F2 only: %d\n", length(unique_f1f2_up)))
cat(sprintf("- Pathways UP F2→F3 only: %d\n", length(unique_f2f3_up)))
cat(sprintf("- Pathways UP both:       %d\n\n", length(shared_up)))

if (!is.null(nmf_per_stage)) {
  cat("## NMF k=6 program activity per stage\n\n")
  cat("| Program | F0 | F1 | F2 | F3 | F4 |\n|---|---:|---:|---:|---:|---:|\n")
  prog_table <- reshape(nmf_per_stage[, c("stage", "program", "mean_score")],
                        idvar = "program", timevar = "stage", direction = "wide")
  for (i in seq_len(nrow(prog_table))) {
    cat(sprintf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f |\n",
                prog_table$program[i],
                prog_table[["mean_score.0"]][i] %||% NA,
                prog_table[["mean_score.1"]][i] %||% NA,
                prog_table[["mean_score.2"]][i] %||% NA,
                prog_table[["mean_score.3"]][i] %||% NA,
                prog_table[["mean_score.4"]][i] %||% NA))
  }
}

cat("\n## Cell-type proportion shifts (top 5 per transition, BH<0.05)\n\n")
for (transition in c("F0vF1", "F1vF2", "F2vF3", "F3vF4")) {
  sig <- ct_df[ct_df$transition == transition & ct_df$padj < 0.05, ]
  sig <- sig[order(-abs(sig$log2FC_proportion)), ]
  cat(sprintf("\n### %s\n\n", transition))
  if (nrow(sig) > 0) {
    cat("| Cell type | log2FC proportion | BH-padj |\n|---|---:|---:|\n")
    for (i in seq_len(min(5, nrow(sig)))) {
      cat(sprintf("| %s | %.2f | %.2e |\n",
                  sig$celltype[i], sig$log2FC_proportion[i], sig$padj[i]))
    }
  } else {
    cat("(none significant)\n")
  }
}

cat("\n## Sex-stratified effect sizes per transition\n\n")
cat("| Transition | Sex | n_a | n_b | median |LFC| top1000 | n padj<0.05 |\n|---|---|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(sex_df))) {
  cat(sprintf("| %s | %s | %d | %d | %s | %s |\n",
              sex_df$transition[i], sex_df$sex[i],
              sex_df$n_a[i], sex_df$n_b[i],
              ifelse(is.na(sex_df$median_abs_logFC_top1000[i]), "NA",
                     sprintf("%.3f", sex_df$median_abs_logFC_top1000[i])),
              ifelse(is.na(sex_df$n_padj_05[i]), "NA", as.character(sex_df$n_padj_05[i]))))
}
sink()

cat("\n== Two-transition complete. ==\n")

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
