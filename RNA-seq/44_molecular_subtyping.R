#!/usr/bin/env Rscript
# =============================================================================
# 44_molecular_subtyping.R — MASLD k-Program Molecular Subtyping
#
# REFACTORED 2026-04-20: replaced hardcoded k=2 (S1/S2) with k-program
# decomposition. The chosen k comes from Script 92's biology-rubric diagnostic
# (nmf_ksweep_rubric_scores.csv) unless overridden via env var NMF_CHOSEN_K.
#
# Emits:
#   nmf_assignments.csv       — sample_id, dominant_program, dominant_program_code,
#                               P1..Pk continuous scores, legacy nmf_subtype alias
#   program_labels.csv        — program_code, biological_label, top_pathway, top_genes
#   subtype_markers.csv       — long: gene x program x (logFC, t, padj, direction)
#   subtype_pathways.csv      — per-program Hallmark fgsea
#   subtype_clinical_associations.csv — per-program demographics + clinical
#   nmf_metrics.csv           — k-sweep stability metrics (unchanged)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(edgeR)
  library(limma)
  library(NMF)
  library(ConsensusClusterPlus)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

select <- dplyr::select
filter <- dplyr::filter

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

dge_path   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
qc_path    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

out_dir <- file.path(BASE, "RNA-seq/results/subtypes")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module A1: k-Program Molecular Subtyping ===\n")
cat("Start:", format(Sys.time()), "\n\n")

K_MIN <- 2
K_MAX <- 10
NMF_RUNS <- 50
CCP_REPS <- 1000
N_TOP_GENES <- 5000

# ============================================================
# Canonical MASLD biology gene sets for tagging
# ============================================================
FIBROTIC_CANON  <- c(
  "COL1A1","COL1A2","COL3A1","COL5A1","COL6A3","COL15A1",
  "LUM","DCN","FBN1","LOX","SPARC","ELN","FN1",
  "ACTA2","DES","MYOCD","CNN1","TAGLN","MYL9","MYH11","FLNC",
  "PDGFRB","THY1","VIM","CDH11","LTBP2","CTSK",
  "TIMP1","PDPN","SMOC2","EDIL3","POSTN","RCN3","BGN",
  "LRRC15","SFRP4","PTPRQ")
HEPATOCYTE_CANON <- c("ALB","APOB","APOA1","APOA2","TF","SERPINA1","TTR",
                      "HP","HRG","HPX","AHSG","CP","C3","FGA","FGB","FGG",
                      "PCK1","G6PC","HNF4A","HNF1A","CPS1","FABP1")
NEURAL_CANON     <- c("NCAM1","NCAM2","ASCL1","NLGN4X","NLGN4Y","CNTN2",
                      "STMN2","NRXN1","NEFL","GAP43","SYN1","GABRG2",
                      "RBFOX1","RBFOX3","MAP2","TUBB3")
CHOL_DUCTULAR    <- c("KRT19","KRT7","SOX9","EPCAM","CFTR","HNF1B","ANXA4",
                      "SPP1","MUC1","MUC5B","MUC6","MUC17","CLDN4","CLDN7")
METABOLIC_CANON <- c("CYP3A4","CYP2E1","ALB","APOB","APOA1","G6PC","PCK1",
                     "HMGCS2","CPT1A","ACOX1","SLC10A1","ABCB11","FABP1",
                     "ALDOB","GLUD1","TAT","SERPINA1","CPS1","HNF4A")
INFLAM_CANON    <- c("CCL2","CXCL10","IL6","TNF","NFKB1","STAT1","IRF1",
                     "CD14","CD68","IL1B","NLRP3","S100A8","S100A9","TLR4",
                     "CCL20","CXCL8")
YCHR_GENES      <- c("KDM5D","USP9Y","UTY","RPS4Y1","DDX3Y","EIF1AY","ZFY",
                     "NLGN4Y","TMSB4Y","PRKY","TTTY14","TTTY15","TTTY10",
                     "TXLNGY","XKRY","AMELY","SRY","TBL1Y")
XIST_SET        <- c("XIST","TSIX")
HCC_CANON       <- c("AKR1B10","GPC3","AFP","HKDC1","IGF2BP3","SPINK1")
CHOLESTATIC     <- CHOL_DUCTULAR  # alias
CELLCYCLE_CANON <- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","AURKB","BIRC5")

# Labels are DIRECTIONAL along the MASLD clinical spectrum (Option A),
# anchored to stage-level dominance-rate Spearman ρ (not per-sample H ρ).
# Only "Fibrogenic" uses gene-level biology (HSC activation markers).
#
# stage_rho[j] = cor(p_fstage[j, 0:4], 0:4, method="spearman")
#   where p_fstage[j, s] = Pr(dominant_program == j | F == s)
tag_program <- function(top50_syms, fgsea_tbl,
                        stage_rho = NA_real_,
                        switch_ratio = NA_real_, pct_dom = NA_real_) {
  sig <- as.data.table(fgsea_tbl)[padj < 0.05 & NES > 0][order(-NES)]
  emt_NES <- sig[grepl("EPITHELIAL_MESENCHYMAL|TGF_BETA|MATRISOME|ANGIOGENESIS|COAGULATION",
                       pathway), max(NES, na.rm = TRUE)]
  inf_NES <- sig[grepl("TNFA|INFLAMMATORY|INTERFERON|IL6|IL2|COMPLEMENT",
                       pathway), max(NES, na.rm = TRUE)]
  if (!is.finite(emt_NES)) emt_NES <- 0
  if (!is.finite(inf_NES)) inf_NES <- 0
  n_fib <- sum(top50_syms %in% FIBROTIC_CANON)

  # 1. Fibrogenic (biology override): HSC/myofibroblast + EMT + rising.
  if (n_fib >= 3 && emt_NES >= 2 && !is.na(stage_rho) && stage_rho > 0.5) {
    return("Fibrogenic")
  }

  # 2. Directional labels based on stage-level dominance trend
  if (is.na(stage_rho)) return("Unlabeled")

  if (stage_rho > 0.6) {
    if (inf_NES >= 2) return("Progression-Inflammatory")
    return("Progression-Other")
  }
  if (stage_rho < -0.6) return("Quiescent-Parenchyma")
  if (!is.na(switch_ratio) && switch_ratio > 1.3) return("Transition-Late")
  if (!is.na(switch_ratio) && switch_ratio < 0.77) return("Transition-Early")
  "Stable"
}

# ============================================================
# 1. Determine chosen k (env override, else from rubric)
# ============================================================
cat("--- Choosing k ---\n")
chosen_k_env <- Sys.getenv("NMF_CHOSEN_K", unset = "")
rubric_path  <- file.path(out_dir, "nmf_ksweep_rubric_scores.csv")

if (chosen_k_env != "") {
  chosen_k <- as.integer(chosen_k_env)
  cat(sprintf("  chosen_k = %d (from NMF_CHOSEN_K env var)\n", chosen_k))
} else if (file.exists(rubric_path)) {
  rubric <- fread(rubric_path)
  full_pass <- rubric[score == 6][order(k)]
  if (nrow(full_pass) > 0) {
    chosen_k <- full_pass$k[1]
    cat(sprintf("  chosen_k = %d (smallest k with score=6/6)\n", chosen_k))
  } else {
    best <- rubric[order(-score, k)]
    chosen_k <- best$k[1]
    cat(sprintf("  chosen_k = %d (best score=%d/6, no full pass)\n",
                chosen_k, best$score[1]))
  }
} else {
  stop("No NMF_CHOSEN_K env var and no ", rubric_path,
       "\nRun Script 92 first (sbatch RNA-seq/run_92_ksweep.sbatch).")
}

# ============================================================
# 2. Load data (disease samples only, same as before)
# ============================================================
cat("\n--- Loading data ---\n")
dge <- readRDS(dge_path)
meta <- fread(meta_path)
qc <- fread(qc_path)
meta <- meta %>% left_join(qc %>% select(sample_id, pass_technical), by = "sample_id")

disease_samples <- meta %>%
  filter(group_binary == "Disease" & pass_technical == TRUE) %>%
  pull(sample_id)
available_samples <- intersect(disease_samples, colnames(dge))
cat(sprintf("  Disease + QC-passing + in DGEList: %d samples\n", length(available_samples)))

dge_sub <- dge[, available_samples]
meta_sub <- meta %>% filter(sample_id %in% available_samples)

# ============================================================
# 3. Normalize + select top variable genes (reuse cached if possible)
# ============================================================
# Cache precedence: NMF_CACHE_PATH env var > clean (Script 95) > nosex (Script 94) > full (Script 92).
env_cache   <- Sys.getenv("NMF_CACHE_PATH", "")
clean_cache <- file.path(out_dir, "nmf_results_cache_clean.rds")
nosex_cache <- file.path(out_dir, "nmf_results_cache_nosex.rds")
full_cache  <- file.path(out_dir, "nmf_results_cache.rds")
cache_path  <- if (nzchar(env_cache)) {
  env_cache
} else if (file.exists(clean_cache)) {
  clean_cache
} else if (file.exists(nosex_cache)) {
  nosex_cache
} else {
  full_cache
}
if (!file.exists(cache_path)) {
  stop("No NMF cache at ", clean_cache, ", ", nosex_cache, ", or ", full_cache,
       " — run Script 95 (confounder-stripped) / 94 / 92 first.")
}
cat("\n--- Loading cached NMF fits from ", basename(cache_path), " ---\n")
cached <- readRDS(cache_path)
nmf_results <- cached$nmf_results
mat_nn <- cached$mat_nn
mat <- cached$mat
metrics <- cached$metrics
cat(sprintf("  Cache keys: %s\n", paste(names(nmf_results), collapse=",")))

if (!as.character(chosen_k) %in% names(nmf_results)) {
  cat(sprintf("  k=%d not cached — refitting (nrun=%d)\n", chosen_k, NMF_RUNS))
  set.seed(42 + chosen_k)
  fit <- nmf(mat_nn, rank = chosen_k, nrun = NMF_RUNS,
             method = "brunet", seed = "random",
             .options = list(verbose = FALSE))
  nmf_results[[as.character(chosen_k)]] <- fit
  saveRDS(list(nmf_results = nmf_results, metrics = metrics,
               mat_nn = mat_nn, mat = mat), cache_path)
}

res_k <- nmf_results[[as.character(chosen_k)]]
W <- basis(res_k)   # genes x k
H <- coef(res_k)    # k x samples
rownames(W) <- rownames(mat_nn)
colnames(H) <- colnames(mat_nn)

# ============================================================
# 4. Ensembl -> symbol map for labeling
# ============================================================
atlas <- fread(atlas_path, select = c("ensembl_id","human_symbol"))
id2sym <- setNames(atlas$human_symbol, atlas$ensembl_id)
strip_ver <- function(x) sub("\\.\\d+$", "", x)
sym_W <- id2sym[strip_ver(rownames(W))]; names(sym_W) <- rownames(W)

# ============================================================
# 5. Tag programs biologically + top genes + fgsea
# ============================================================
cat("\n--- Tagging programs + fgsea ---\n")
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hallmark$gene_symbol, hallmark$gs_name)

W_centered <- sapply(seq_len(chosen_k), function(j) W[, j] - rowMeans(W[, -j, drop = FALSE]))
rownames(W_centered) <- rownames(W)

# Pre-compute per-program clinical behavior so labels can be directional.
fib_int <- as.integer(meta_sub$fibrosis_stage[match(colnames(H), meta_sub$sample_id)])
dom_idx_pre <- apply(H, 2, which.max)

# Stage-level dominance rate per program (cleaner signal than per-sample H ρ)
stage_rho_vec <- vapply(seq_len(chosen_k), function(j) {
  prop <- vapply(0:4, function(s) {
    n_s <- sum(fib_int == s, na.rm = TRUE)
    if (n_s == 0) NA_real_ else mean(dom_idx_pre[fib_int == s] == j, na.rm = TRUE)
  }, numeric(1))
  if (sum(!is.na(prop)) < 3) return(NA_real_)
  suppressWarnings(cor(prop, 0:4, method = "spearman", use = "pairwise"))
}, numeric(1))

program_labels <- character(chosen_k)
top_genes_per_prog <- list()
fgsea_per_prog     <- list()
clin_pre <- data.frame(program = paste0("P", seq_len(chosen_k)),
                       stage_rho = stage_rho_vec,
                       rho_fib = NA_real_, switch_ratio = NA_real_,
                       pct_dom = NA_real_)
for (j in seq_len(chosen_k)) {
  w_j <- W_centered[, j]
  ord <- order(w_j, decreasing = TRUE)
  top50_ids <- rownames(W)[ord][1:50]
  top50_syms <- unname(sym_W[top50_ids])
  top50_syms <- top50_syms[!is.na(top50_syms) & top50_syms != ""]
  top_genes_per_prog[[j]] <- top50_syms

  has_sym <- !is.na(sym_W) & sym_W != ""
  ranked <- w_j[has_sym]; names(ranked) <- sym_W[has_sym]
  ranked <- ranked[!duplicated(names(ranked))]
  ranked <- sort(ranked, decreasing = TRUE)
  set.seed(42)
  fg <- suppressWarnings(fgsea(pathways = hall_list, stats = ranked,
                               minSize = 10, maxSize = 500,
                               nPermSimple = 10000, scoreType = "std"))
  fg$program <- paste0("P", j)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, FUN.VALUE=character(1), collapse=";")
  fgsea_per_prog[[j]] <- as.data.table(fg)

  # Clinical behavior
  h_j <- H[j, ]
  rho_fib <- suppressWarnings(cor(h_j, fib_int, method = "spearman",
                                  use = "pairwise.complete"))
  ok_f <- !is.na(fib_int)
  m_lo <- mean(h_j[ok_f & fib_int <= 2], na.rm = TRUE)
  m_hi <- mean(h_j[ok_f & fib_int >= 3], na.rm = TRUE)
  sw_ratio <- m_hi / (m_lo + 1e-10)
  pct_dom <- mean(dom_idx_pre == j) * 100
  clin_pre[j, c("rho_fib","switch_ratio","pct_dom")] <-
    list(rho_fib, sw_ratio, pct_dom)

  program_labels[j] <- tag_program(top50_syms, fg,
                                    stage_rho = stage_rho_vec[j],
                                    switch_ratio = sw_ratio,
                                    pct_dom = pct_dom)
  cat(sprintf("  P%d: %-26s stage_ρ=%+.2f (H ρ_fib=%+.2f) sw=%.2f cov=%4.1f%%  top: %s\n",
              j, program_labels[j], stage_rho_vec[j], rho_fib, sw_ratio, pct_dom,
              paste(head(top50_syms, 6), collapse=", ")))
}

# Make unique biological names: if two programs share the same label, append Roman numerals
lab_counts <- table(program_labels)
program_codes <- paste0("P", seq_len(chosen_k))
program_bio_names <- program_labels
for (lbl in names(lab_counts)) {
  if (lab_counts[lbl] > 1) {
    idx <- which(program_labels == lbl)
    for (ii in seq_along(idx)) {
      program_bio_names[idx[ii]] <- paste0(lbl, "_", ii)
    }
  }
}

program_labels_dt <- data.table(
  program_code  = program_codes,
  biological_label = program_bio_names,
  label_category = program_labels,
  top_pathway = vapply(seq_len(chosen_k), function(j) {
    s <- fgsea_per_prog[[j]][padj < 0.05 & NES > 0][order(-NES)]
    if (nrow(s)) sprintf("%s(NES=%.2f)", s$pathway[1], s$NES[1]) else "ns"
  }, character(1)),
  top_genes = vapply(top_genes_per_prog,
                     function(x) paste(head(x, 15), collapse=","),
                     character(1))
)
fwrite(program_labels_dt, file.path(out_dir, "program_labels.csv"))
cat("\nProgram labels:\n"); print(program_labels_dt[, .(program_code, biological_label, top_pathway)], row.names = FALSE)

# ============================================================
# 6. Per-sample assignments (dominant + continuous)
# ============================================================
cat("\n--- Per-sample assignments ---\n")
dom_idx <- apply(H, 2, which.max)
dom_program_code <- program_codes[dom_idx]
dom_program_bio  <- program_bio_names[dom_idx]

# Legacy S1/S2 alias: Fibrotic program -> S2, everything else -> S1
# (preserves the S2 = "fibrogenic progressor" semantics of prior work)
fibrotic_codes <- program_codes[program_labels %in% c("Fibrogenic")]
legacy_subtype <- ifelse(dom_program_code %in% fibrotic_codes, "S2", "S1")

H_scores <- as.data.table(t(H))
setnames(H_scores, old = names(H_scores), new = program_codes)
# Also add biological-label-aliased score columns
bio_score_cols <- paste0("P_", gsub("[-/ ]", "_", program_bio_names), "_score")
for (j in seq_len(chosen_k)) H_scores[[bio_score_cols[j]]] <- H_scores[[program_codes[j]]]

final_classes <- data.table(
  sample_id = colnames(H),
  dominant_program_code = dom_program_code,
  dominant_program      = dom_program_bio,
  nmf_subtype           = legacy_subtype  # S1 / S2 legacy for transitional downstream
)
final_classes <- cbind(final_classes, H_scores)
fwrite(final_classes, file.path(out_dir, "nmf_assignments.csv"))
cat(sprintf("  nmf_assignments.csv: %d samples\n", nrow(final_classes)))
cat("  Dominant program distribution:\n"); print(table(final_classes$dominant_program))
cat("  Legacy S1/S2 mapping:\n"); print(table(final_classes$nmf_subtype))

# ============================================================
# 7. Per-program marker DE (limma one-vs-rest)
# ============================================================
cat("\n--- Per-program one-vs-rest limma DE ---\n")
meta_sub <- meta_sub %>% mutate(
  dominant_program_code = dom_program_code[match(sample_id, colnames(H))],
  dominant_program      = dom_program_bio[match(sample_id, colnames(H))],
  nmf_subtype           = legacy_subtype[match(sample_id, colnames(H))]
)

design <- model.matrix(~ 0 + factor(dominant_program_code, levels = program_codes) + factor(dataset),
                       data = meta_sub)
colnames(design) <- gsub("factor\\(dominant_program_code, levels = program_codes\\)", "",
                         colnames(design))
colnames(design) <- gsub("factor\\(dataset\\)", "dataset_", colnames(design))

v2 <- voom(dge_sub, design = design, plot = FALSE)
fit <- lmFit(v2, design)

all_markers <- list()
for (j in seq_len(chosen_k)) {
  pc <- program_codes[j]
  others <- setdiff(program_codes, pc)
  contr_vec <- rep(0, ncol(design)); names(contr_vec) <- colnames(design)
  contr_vec[pc] <- 1
  contr_vec[others] <- -1 / length(others)
  contr_mat <- matrix(contr_vec, ncol = 1); colnames(contr_mat) <- pc
  fit2 <- eBayes(contrasts.fit(fit, contr_mat))
  tt <- topTable(fit2, number = Inf, sort.by = "t")
  tt$gene <- rownames(tt)
  tt$program_code <- pc
  tt$program <- program_bio_names[j]
  tt$direction <- ifelse(tt$logFC > 0, "up", "down")
  # Keep all rows (not truncate at 50) so downstream can filter by padj/logFC
  all_markers[[j]] <- as.data.table(tt[, c("gene","program_code","program","logFC","AveExpr","t","P.Value","adj.P.Val","direction")])
  setnames(all_markers[[j]], c("P.Value","adj.P.Val"), c("pvalue","padj"))
  # Legacy S1/S2 alias: S2 = Fibrotic program, S1 = any other program.
  # Kept for transitional compat with Scripts 209/215/figS08.
  all_markers[[j]][, subtype := ifelse(program_code %in% fibrotic_codes, "S2", "S1")]
  # Legacy: adj.P.Val column name expected by some consumers
  all_markers[[j]][, adj.P.Val := padj]
}
all_markers_dt <- rbindlist(all_markers)
fwrite(all_markers_dt, file.path(out_dir, "subtype_markers.csv"))
cat(sprintf("  subtype_markers.csv: %d rows (%d genes x %d programs)\n",
            nrow(all_markers_dt), nrow(all_markers_dt)/chosen_k, chosen_k))

# ============================================================
# 8. Per-program Hallmark pathways (rewrite output)
# ============================================================
cat("\n--- Per-program pathway enrichment ---\n")
pathway_results <- rbindlist(fgsea_per_prog, fill = TRUE)
pathway_results[, program_bio := program_bio_names[match(program, program_codes)]]
fwrite(pathway_results, file.path(out_dir, "subtype_pathways.csv"))
cat(sprintf("  subtype_pathways.csv: %d rows\n", nrow(pathway_results)))

# ============================================================
# 9. Per-program clinical associations
# ============================================================
cat("\n--- Per-program clinical associations ---\n")
clinical_assoc <- meta_sub %>%
  filter(!is.na(dominant_program_code)) %>%
  group_by(dominant_program_code, dominant_program) %>%
  summarise(
    n_samples = n(),
    n_datasets = n_distinct(dataset),
    datasets = paste(unique(dataset), collapse = ";"),
    n_male = sum(sex == "M", na.rm = TRUE),
    n_female = sum(sex == "F", na.rm = TRUE),
    pct_male = mean(sex == "M", na.rm = TRUE) * 100,
    mean_age = mean(as.numeric(age), na.rm = TRUE),
    n_nafl = sum(diagnosis_harmonized == "NAFL", na.rm = TRUE),
    n_nash = sum(diagnosis_harmonized %in% c("NASH","Borderline"), na.rm = TRUE),
    pct_nash = mean(diagnosis_harmonized %in% c("NASH","Borderline"), na.rm = TRUE) * 100,
    mean_fibrosis = mean(as.numeric(fibrosis_stage), na.rm = TRUE),
    mean_nas = mean(as.numeric(nas_score), na.rm = TRUE),
    .groups = "drop"
  )
fwrite(clinical_assoc, file.path(out_dir, "subtype_clinical_associations.csv"))
cat("\nClinical associations by program:\n")
print(clinical_assoc, row.names = FALSE)

# Multi-way contingency: NASH x program
if (nrow(meta_sub) > 0) {
  nash_tab <- meta_sub %>%
    filter(!is.na(dominant_program_code) & !is.na(diagnosis_harmonized)) %>%
    mutate(is_nash = diagnosis_harmonized %in% c("NASH","Borderline")) %>%
    with(table(dominant_program_code, is_nash))
  if (all(dim(nash_tab) >= 2)) {
    cs <- chisq.test(nash_tab)
    cat(sprintf("\nNASH x program chi-sq: X2=%.2f df=%d p=%.2e\n",
                cs$statistic, cs$parameter, cs$p.value))
  }
}

# ============================================================
# 10. Figures: metrics + marker heatmap
# ============================================================
cat("\n--- Figures ---\n")

pdf(file.path(fig_dir, "nmf_metrics.pdf"), width = 10, height = 4)
par(mfrow = c(1, 3))
plot(metrics$k, metrics$cophenetic, type = "b", pch = 19,
     xlab = "k", ylab = "Cophenetic", main = "NMF Cophenetic")
abline(h = 0.9, lty = 2, col = "red")
plot(metrics$k, metrics$dispersion, type = "b", pch = 19,
     xlab = "k", ylab = "Dispersion", main = "NMF Dispersion")
plot(metrics$k, metrics$silhouette, type = "b", pch = 19,
     xlab = "k", ylab = "Silhouette", main = "NMF Silhouette")
dev.off()

# Marker heatmap
top_markers <- all_markers_dt[direction == "up" & padj < 0.05][order(program_code, -t)]
top_markers <- top_markers[, head(.SD, 20), by = program_code]
tm_genes <- unique(top_markers$gene)
tm_genes <- intersect(tm_genes, rownames(mat))
if (length(tm_genes) > 5) {
  m_scaled <- t(scale(t(mat[tm_genes, , drop = FALSE])))
  ann_col <- data.frame(
    Program = meta_sub$dominant_program[match(colnames(m_scaled), meta_sub$sample_id)],
    Dataset = meta_sub$dataset[match(colnames(m_scaled), meta_sub$sample_id)],
    row.names = colnames(m_scaled)
  )
  k_palette <- setNames(brewer.pal(max(3, chosen_k), "Set2")[1:chosen_k], program_bio_names)
  pdf(file.path(fig_dir, "program_marker_heatmap.pdf"), width = 14, height = 10)
  pheatmap(m_scaled, annotation_col = ann_col, show_colnames = FALSE,
           clustering_method = "ward.D2",
           annotation_colors = list(Program = k_palette),
           main = sprintf("Top 20/program up-regulated markers (k=%d)", chosen_k))
  dev.off()
}

# ============================================================
# 11. Summary
# ============================================================
cat("\n--- Final Summary ---\n")
cat(sprintf("Chosen k: %d\n", chosen_k))
cat(sprintf("Cophenetic at k=%d: %.3f\n",
            chosen_k, metrics$cophenetic[metrics$k == chosen_k]))
cat("Program → biological label → dominant-sample count:\n")
print(data.table(
  code = program_codes,
  label = program_bio_names,
  n = as.integer(table(factor(dom_program_code, levels = program_codes)))
), row.names = FALSE)
cat(sprintf("\n=== Module A1 k-program refactor complete (%s) ===\n",
            format(Sys.time())))
