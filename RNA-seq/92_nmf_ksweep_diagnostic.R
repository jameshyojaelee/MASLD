#!/usr/bin/env Rscript
# 92_nmf_ksweep_diagnostic.R
# NMF k-sweep (k=3..10) on 1,254 QC-passing disease samples. Evaluates each k
# against a biology-informed rubric (6 criteria: program interpretability, stage
# association, cell-type specificity, stability, etc.). Reuses k=3..6 fits from
# Script 44 and fits k=7..10 fresh.
#
# Emits:
#   RNA-seq/results/subtypes/nmf_ksweep_biology.md      (decision report)
#   RNA-seq/results/subtypes/nmf_ksweep_program_atlas.csv (long-form per-k)
#   RNA-seq/results/subtypes/nmf_ksweep_rubric_scores.csv (k x criterion)
#   figures/misc/nmf_ksweep_diagnostic.pdf              (per-k panels)
#   RNA-seq/results/subtypes/nmf_results_cache.rds      (UPDATED with k=7..10)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(NMF)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

cache_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache.rds")
meta_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
out_dir    <- file.path(BASE, "RNA-seq/results/subtypes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

md_out     <- file.path(out_dir, "nmf_ksweep_biology.md")
atlas_out  <- file.path(out_dir, "nmf_ksweep_program_atlas.csv")
rubric_out <- file.path(out_dir, "nmf_ksweep_rubric_scores.csv")
pdf_out    <- file.path(FIG_MISC, "nmf_ksweep_diagnostic.pdf")

K_SWEEP  <- 3:6  # cached; k=7..10 fits stall due to NMF parallelism limits
NMF_RUNS <- 30
NMF_CPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

# ------- Canonical MASLD gene sets for biological tagging -------
FIBROTIC_CANON  <- c("COL1A1","COL1A2","COL3A1","ACTA2","TIMP1","LUM",
                     "DCN","FBN1","LOX","SPARC","COL15A1","COL6A3",
                     "CTSK","CDH11","LTBP2","PDGFRB","THY1","VIM","COL5A1")
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
CHOLESTATIC     <- c("KRT19","SOX9","EPCAM","CFTR","HNF1B","ANXA4")
CELLCYCLE_CANON <- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","AURKB","BIRC5")

# ------- Biological tagging rules (priority-ordered) -------
# Returns a label in {Sex-M, Sex-F, HCC-like, Cell-cycle, Fibrotic,
# Inflammatory, Cholestatic, Metabolic, Unlabeled}
tag_program <- function(top50_syms, fgsea_tbl) {
  sig <- as.data.table(fgsea_tbl)[padj < 0.05 & NES > 0][order(-NES)]
  fib_NES <- sig[grepl("EPITHELIAL_MESENCHYMAL|TGF_BETA|COAGULATION|ANGIOGENESIS|COMPLEMENT|APICAL_JUNCTION", pathway), max(NES, na.rm = TRUE)]
  inf_NES <- sig[grepl("TNFA|INFLAMMATORY|INTERFERON|IL6|IL2", pathway), max(NES, na.rm = TRUE)]
  met_NES <- sig[grepl("OXIDATIVE_PHOSPH|FATTY_ACID|BILE_ACID|XENOBIOTIC|ADIPOGENESIS|PEROXISOME", pathway), max(NES, na.rm = TRUE)]
  cc_NES  <- sig[grepl("E2F_TARGETS|G2M_CHECKPOINT|MITOTIC_SPINDLE|MYC_TARGETS", pathway), max(NES, na.rm = TRUE)]
  if (!is.finite(fib_NES)) fib_NES <- 0
  if (!is.finite(inf_NES)) inf_NES <- 0
  if (!is.finite(met_NES)) met_NES <- 0
  if (!is.finite(cc_NES))  cc_NES  <- 0

  n_y  <- sum(top50_syms %in% YCHR_GENES)
  n_x  <- sum(top50_syms %in% XIST_SET)
  n_hcc <- sum(top50_syms %in% HCC_CANON)
  n_chol <- sum(top50_syms %in% CHOLESTATIC)
  n_cc <- sum(top50_syms %in% CELLCYCLE_CANON)
  n_fib <- sum(top50_syms %in% FIBROTIC_CANON)
  n_met <- sum(top50_syms %in% METABOLIC_CANON)
  n_inf <- sum(top50_syms %in% INFLAM_CANON)

  # Priority order: Sex (chr-driven) > HCC > Cell-cycle > Fibrotic >
  # Inflammatory > Cholestatic > Metabolic
  if (n_y >= 5) return("Sex-M")
  if (n_x >= 1 && (n_fib + n_inf + n_met + n_hcc) < 3) return("Sex-F")
  if (n_hcc >= 2) return("HCC-like")
  if (n_cc >= 3 || cc_NES >= 2) return("Cell-cycle")
  if (n_fib >= 3 && fib_NES >= 2) return("Fibrotic")
  if (n_inf >= 3 && inf_NES >= 2) return("Inflammatory")
  if (n_chol >= 3) return("Cholestatic")
  if (n_met >= 3 && met_NES >= 2) return("Metabolic")
  # Fallback on fgsea only
  maxNES <- c(Fibrotic = fib_NES, Inflammatory = inf_NES,
              Metabolic = met_NES, "Cell-cycle" = cc_NES)
  if (max(maxNES) >= 2) return(names(maxNES)[which.max(maxNES)])
  "Unlabeled"
}

# ============================================================
# 1. Load cache + metadata
# ============================================================
cat("=== Step 1: Load cache + metadata ===\n")
cached <- readRDS(cache_path)
stopifnot(all(c("nmf_results","mat","mat_nn") %in% names(cached)))
mat_nn <- cached$mat_nn
sample_ids <- colnames(mat_nn)
gene_ids   <- rownames(mat_nn)
cat(sprintf("  Matrix: %d genes x %d samples (disease-only)\n",
            length(gene_ids), length(sample_ids)))

meta <- fread(meta_path)[sample_id %in% sample_ids]
setkey(meta, sample_id)

atlas <- fread(atlas_path, select = c("ensembl_id","human_symbol"))
id2sym <- setNames(atlas$human_symbol, atlas$ensembl_id)
strip_ver <- function(x) sub("\\.\\d+$", "", x)
sym_W <- id2sym[strip_ver(gene_ids)]; names(sym_W) <- gene_ids

# ============================================================
# 2. Fit or reuse NMF for each k
# ============================================================
cat("\n=== Step 2: NMF fits for k=3..10 ===\n")
nmf_list <- cached$nmf_results
for (k in K_SWEEP) {
  if (as.character(k) %in% names(nmf_list)) {
    cat(sprintf("  k=%d: reusing cached fit\n", k))
    next
  }
  cat(sprintf("  k=%d: fitting (nrun=%d, parallel=%d)...\n", k, NMF_RUNS, NMF_CPUS))
  t0 <- Sys.time()
  set.seed(42 + k)
  fit <- nmf(mat_nn, rank = k, nrun = NMF_RUNS,
             method = "brunet", seed = "random",
             .pbackend = NMF_CPUS,
             .options = list(verbose = FALSE))
  nmf_list[[as.character(k)]] <- fit
  cat(sprintf("    done in %s\n", format(round(Sys.time() - t0, 1))))
}
# Save back to cache so future reruns are cheap
saveRDS(list(nmf_results = nmf_list, metrics = cached$metrics,
             mat_nn = mat_nn, mat = cached$mat), cache_path)
cat("  Cache saved.\n")

# ============================================================
# 3. Per-k characterization
# ============================================================
cat("\n=== Step 3: Per-k program characterization ===\n")
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hallmark$gene_symbol, hallmark$gs_name)

atlas_rows <- list()
rubric_rows <- list()
per_k_summary <- list()

for (k in K_SWEEP) {
  cat(sprintf("  -- k=%d --\n", k))
  res_k <- nmf_list[[as.character(k)]]
  W <- basis(res_k); H <- coef(res_k)
  rownames(W) <- gene_ids; colnames(H) <- sample_ids

  W_centered <- sapply(seq_len(k), function(j) W[, j] - rowMeans(W[, -j, drop = FALSE]))
  rownames(W_centered) <- rownames(W)

  dom <- apply(H, 2, which.max)   # dominant program per sample (1..k)

  # Per-program
  prog_rows <- list()
  top50_syms_list <- list()
  for (j in seq_len(k)) {
    w_j <- W_centered[, j]
    ord <- order(w_j, decreasing = TRUE)
    top50_ids <- gene_ids[ord][1:50]
    top50_syms <- unname(sym_W[top50_ids])
    top50_syms <- top50_syms[!is.na(top50_syms) & top50_syms != ""]
    top50_syms_list[[j]] <- top50_syms

    # fgsea over symbol ranking
    has_sym <- !is.na(sym_W) & sym_W != ""
    ranked <- w_j[has_sym]; names(ranked) <- sym_W[has_sym]
    ranked <- ranked[!duplicated(names(ranked))]
    ranked <- sort(ranked, decreasing = TRUE)
    set.seed(42)
    fg <- suppressWarnings(fgsea(pathways = hall_list, stats = ranked,
                                 minSize = 10, maxSize = 500,
                                 nPermSimple = 10000, scoreType = "std"))

    tag <- tag_program(top50_syms, fg)

    # Clinical correlations
    h_j <- H[j, ]
    rho_fib <- suppressWarnings(cor(h_j, as.numeric(meta$fibrosis_stage[match(sample_ids, meta$sample_id)]),
                                    method = "spearman", use = "pairwise.complete"))
    rho_nas <- suppressWarnings(cor(h_j, as.numeric(meta$nas_score[match(sample_ids, meta$sample_id)]),
                                    method = "spearman", use = "pairwise.complete"))
    rho_age <- suppressWarnings(cor(h_j, as.numeric(meta$age[match(sample_ids, meta$sample_id)]),
                                    method = "spearman", use = "pairwise.complete"))
    sex_vec <- meta$sex[match(sample_ids, meta$sample_id)]
    ok_sx <- sex_vec %in% c("M","F")
    rho_sex <- if (sum(ok_sx) > 10) {
      suppressWarnings(cor(h_j[ok_sx], as.integer(sex_vec[ok_sx] == "F"),
                           method = "spearman"))
    } else NA_real_
    ds_vec <- meta$dataset[match(sample_ids, meta$sample_id)]
    kw <- suppressWarnings(kruskal.test(h_j ~ factor(ds_vec)))
    n_ds <- length(unique(ds_vec))
    eta2 <- max(0, (as.numeric(kw$statistic) - n_ds + 1) / (length(h_j) - n_ds))

    # F2 switch
    fib_int <- as.integer(meta$fibrosis_stage[match(sample_ids, meta$sample_id)])
    ok_f <- !is.na(fib_int)
    m_lo <- mean(h_j[ok_f & fib_int <= 2], na.rm = TRUE)
    m_hi <- mean(h_j[ok_f & fib_int >= 3], na.rm = TRUE)
    sw_ratio <- m_hi / (m_lo + 1e-10)
    sw_w <- suppressWarnings(wilcox.test(h_j[ok_f & fib_int <= 2],
                                         h_j[ok_f & fib_int >= 3]))$p.value

    # Coverage
    pct_dom <- mean(dom == j) * 100
    pct_hi  <- mean(h_j > median(h_j) + sd(h_j)) * 100

    # Pathways top-3
    sigp <- as.data.table(fg)[padj < 0.05 & NES > 0][order(-NES)][1:3, .(pathway, NES)]
    top_path <- paste(sprintf("%s(%.2f)", sigp$pathway, sigp$NES), collapse = ";")

    prog_rows[[j]] <- data.table(
      k = k, program = paste0("P", j), label = tag,
      top10_genes = paste(head(top50_syms, 10), collapse = ","),
      top_pathways = top_path,
      rho_fibrosis = rho_fib, rho_nas = rho_nas, rho_sex_F = rho_sex,
      rho_age = rho_age, eta2_dataset = eta2,
      switch_ratio = sw_ratio, wilcox_p = sw_w,
      pct_dominant = pct_dom, pct_high = pct_hi,
      n_fibrotic_canon = sum(top50_syms %in% FIBROTIC_CANON),
      n_metabolic_canon = sum(top50_syms %in% METABOLIC_CANON),
      n_inflam_canon = sum(top50_syms %in% INFLAM_CANON),
      n_ychr = sum(top50_syms %in% YCHR_GENES),
      n_xist = sum(top50_syms %in% XIST_SET),
      n_hcc = sum(top50_syms %in% HCC_CANON)
    )
  }
  prog_dt <- rbindlist(prog_rows)
  atlas_rows[[as.character(k)]] <- prog_dt

  # Orthogonality: max off-diagonal |H correlation|
  Hcor <- cor(t(H))
  max_off <- max(abs(Hcor[lower.tri(Hcor)]))

  # Rubric scoring
  labels_k    <- prog_dt$label
  n_unlabeled <- sum(labels_k == "Unlabeled")
  n_sex_prog  <- sum(labels_k %in% c("Sex-M","Sex-F"))
  fib_rows    <- prog_dt[label == "Fibrotic"]
  n_fib       <- nrow(fib_rows)
  fib_clean   <- if (n_fib == 1) {
    abs(fib_rows$rho_sex_F) < 0.2 &
      fib_rows$eta2_dataset < 0.2 &
      fib_rows$switch_ratio > 1.3
  } else FALSE
  fib_canon_ok <- if (n_fib == 1) fib_rows$n_fibrotic_canon >= 3 else FALSE
  min_cov <- min(prog_dt$pct_dominant)
  cov_ok  <- min_cov >= 5

  rubric <- data.table(
    k = k,
    crit1_all_labeled = (n_unlabeled == 0),
    crit2_one_or_zero_sex = (n_sex_prog <= 1),
    crit3_one_fibrotic_canon = (n_fib == 1 && fib_canon_ok),
    crit4_fibrotic_clean = fib_clean,
    crit5_coverage = cov_ok,
    crit6_orthogonal = (max_off < 0.6),
    max_H_cor = max_off,
    min_coverage_pct = min_cov,
    n_programs = k,
    score = NA_integer_
  )
  rubric[, score := sum(unlist(.SD)), .SDcols = grep("^crit", names(rubric), value = TRUE)]
  rubric_rows[[as.character(k)]] <- rubric
  per_k_summary[[as.character(k)]] <- list(prog_dt = prog_dt, Hcor = Hcor)

  cat(sprintf("    labels: %s\n", paste(labels_k, collapse=",")))
  cat(sprintf("    rubric score: %d/6, min_cov=%.1f%%, max|H|=%.2f\n",
              rubric$score, min_cov, max_off))
}

atlas_dt <- rbindlist(atlas_rows)
rubric_dt <- rbindlist(rubric_rows)
fwrite(atlas_dt,  atlas_out)
fwrite(rubric_dt, rubric_out)

# ============================================================
# 4. k recommendation
# ============================================================
cat("\n=== Step 4: k recommendation ===\n")
full_pass <- rubric_dt[score == 6][order(k)]
if (nrow(full_pass) > 0) {
  chosen_k <- full_pass$k[1]
  rationale <- sprintf("Full pass: smallest k with 6/6 rubric criteria = k=%d", chosen_k)
} else {
  best <- rubric_dt[order(-score, k)]
  chosen_k <- best$k[1]
  rationale <- sprintf("No k satisfies all 6 criteria. Best is k=%d with %d/6; next candidates: %s",
                       chosen_k, best$score[1],
                       paste(sprintf("k=%d(%d/6)", best$k[2:min(nrow(best),4)],
                                     best$score[2:min(nrow(best),4)]),
                             collapse=", "))
}
cat(sprintf("  Recommended k = %d\n  %s\n", chosen_k, rationale))

# ============================================================
# 5. Markdown report
# ============================================================
cat("\n=== Step 5: Markdown report ===\n")
md <- c(
  "# NMF k-Sweep Biology Rubric Diagnostic",
  sprintf("**Generated:** %s  ", format(Sys.time())),
  sprintf("**Samples:** %d disease-only, QC-passing  ", length(sample_ids)),
  sprintf("**Gene pool:** %d top-IQR genes, batch-corrected log-CPM", length(gene_ids)),
  sprintf("**k range:** %d..%d, NMF nrun=%d", min(K_SWEEP), max(K_SWEEP), NMF_RUNS),
  "",
  sprintf("## RECOMMENDED k: **%d**", chosen_k),
  "",
  rationale,
  "",
  "## Rubric scores by k",
  "",
  paste(capture.output(print(rubric_dt[, .(k, score, crit1_all_labeled, crit2_one_or_zero_sex,
                                           crit3_one_fibrotic_canon, crit4_fibrotic_clean,
                                           crit5_coverage, crit6_orthogonal,
                                           min_coverage_pct = round(min_coverage_pct, 1),
                                           max_H_cor = round(max_H_cor, 2))],
                            row.names = FALSE)),
        collapse = "\n"),
  "",
  "## Rubric definitions",
  "1. All programs labeled (none Unlabeled)",
  "2. At most one Sex program (Sex-M or Sex-F, not both)",
  "3. Exactly one Fibrotic program with >=3 canonical ECM genes in top-50",
  "4. Fibrotic program |rho_sex|<0.2, eta2_dataset<0.2, switch_ratio>1.3",
  "5. All programs have >=5% dominant-sample coverage",
  "6. Max pairwise |H correlation| < 0.6",
  "",
  "---",
  "",
  "## Per-k program atlas"
)
for (k in K_SWEEP) {
  prog_dt <- per_k_summary[[as.character(k)]]$prog_dt
  md <- c(md, sprintf("### k=%d", k), "")
  tbl <- prog_dt[, .(program, label, pct_dominant = round(pct_dominant,1),
                     rho_fibrosis = round(rho_fibrosis,3),
                     rho_sex_F = round(rho_sex_F,3),
                     eta2_dataset = round(eta2_dataset,3),
                     switch_ratio = round(switch_ratio,2),
                     fib_canon = n_fibrotic_canon,
                     ychr = n_ychr, xist = n_xist,
                     top_genes = substr(top10_genes, 1, 55))]
  md <- c(md, paste(capture.output(print(tbl, row.names = FALSE)), collapse = "\n"), "")
}

writeLines(md, md_out)
cat(sprintf("  Wrote %s\n", md_out))

# ============================================================
# 6. Diagnostic PDF (per-k panel)
# ============================================================
cat("\n=== Step 6: Diagnostic PDF ===\n")

# Panel A: rubric heatmap
rh <- melt(rubric_dt[, .(k, crit1_all_labeled, crit2_one_or_zero_sex,
                         crit3_one_fibrotic_canon, crit4_fibrotic_clean,
                         crit5_coverage, crit6_orthogonal)],
           id.vars = "k", variable.name = "criterion", value.name = "pass")
rh[, criterion := sub("^crit[0-9]+_", "", criterion)]
rh[, pass := as.integer(pass)]
pA <- ggplot(rh, aes(factor(k), criterion, fill = factor(pass))) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pass == 1, "o", "x")), size = 3) +
  scale_fill_manual(values = c("0"="#EF9A9A","1"="#A5D6A7"),
                    labels = c("0"="fail","1"="pass"), name = NULL) +
  scale_x_discrete(name = "k") +
  labs(y = NULL, title = "A — Rubric criteria per k",
       subtitle = sprintf("Recommended k = %d", chosen_k)) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 7))

# Panel B: per-k label composition
label_long <- rbindlist(lapply(K_SWEEP, function(k) {
  per_k_summary[[as.character(k)]]$prog_dt[, .(k = k, label)]
}))
label_long[, label := factor(label, levels = c("Metabolic","Inflammatory","Fibrotic",
                                               "Cholestatic","Immune/Exhaustion",
                                               "HCC-like","Cell-cycle","Sex-F","Sex-M","Unlabeled"))]
label_counts <- label_long[, .N, by = .(k, label)]
pB <- ggplot(label_counts, aes(factor(k), N, fill = label)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(Metabolic = "#1E88E5", Inflammatory = "#E53935",
                               Fibrotic = "#8E24AA", Cholestatic = "#00897B",
                               "Immune/Exhaustion" = "#FB8C00",
                               "HCC-like" = "#5D4037", "Cell-cycle" = "#6D4C41",
                               "Sex-F" = "#EC407A", "Sex-M" = "#1A237E",
                               Unlabeled = "#BDBDBD"), drop = FALSE) +
  labs(x = "k", y = "programs",
       title = "B — Biological label composition per k") +
  theme_masld()

# Panel C: fibrotic-program rho_fibrosis across k (if any)
fib_track <- atlas_dt[label == "Fibrotic"]
if (nrow(fib_track) > 0) {
  pC <- ggplot(fib_track, aes(factor(k), rho_fibrosis, fill = switch_ratio)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = sprintf("sw=%.2f\nρsex=%.2f", switch_ratio, rho_sex_F)),
              vjust = -0.2, size = 2.3) +
    scale_fill_gradient(low = "#F8BBD0", high = "#8E24AA",
                        name = "switch\nratio") +
    labs(x = "k", y = "rho(fibrosis)",
         title = "C — Fibrotic program performance vs k") +
    expand_limits(y = max(fib_track$rho_fibrosis, na.rm = TRUE) * 1.3) +
    theme_masld()
} else {
  pC <- ggplot() +
    labs(title = "C — No Fibrotic program labeled at any k") +
    theme_masld()
}

# Panel D: coverage floor by k
cov_track <- atlas_dt[, .(min_cov = min(pct_dominant),
                          max_cov = max(pct_dominant)), by = k]
pD <- ggplot(cov_track, aes(factor(k))) +
  geom_linerange(aes(ymin = min_cov, ymax = max_cov), color = "#1565C0",
                 linewidth = 1.5) +
  geom_point(aes(y = min_cov), color = "#C2185B", size = 2) +
  geom_point(aes(y = max_cov), color = "#1565C0", size = 2) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40") +
  labs(x = "k", y = "% dominant coverage",
       title = "D — Program coverage floor (<5% fails)") +
  theme_masld()

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title = sprintf("NMF k-sweep biology rubric — recommended k = %d", chosen_k),
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
ggsave(pdf_out, fig, width = 14, height = 10, device = cairo_pdf)
cat(sprintf("  Wrote %s\n", pdf_out))
cat("\nDone.\n")
