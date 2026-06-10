#!/usr/bin/env Rscript
# figS_disease_signature_sweep.R
# Extends panelK_supervised_disease (limma-voom, N=200) across:
#   (1) 7 DE method families from the multimethod validation battery
#   (2) Gene number sweep: N = 25, 50, 100, 200, 500, 1000, 2000
#   (3) Random null baseline: 100 replicates per (fold x N)
#   (4) Train-on-1-cohort, project-to-4 (within-test z-scores, limma-voom)
#
# Outputs → figures/supplementary/figS_methods_validation/disease_signature_sweep/panels/

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(pROC)
  library(BiocParallel); library(RColorBrewer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS_METHVAL_DIR, "disease_signature_sweep", "panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_CPUS          <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
BPPARAM_PAR     <- MulticoreParam(N_CPUS)
BPPARAM_SER     <- SerialParam()

CTRL            <- "#9E9E9E"
set.seed(42)

MEGA         <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
                  GSE162694 = "Bril",   GSE213621 = "Chen")

N_VALS      <- c(25L, 50L, 100L, 200L, 500L, 1000L, 2000L)
N_RAND_REPS <- 100L

# ── Data loading ────────────────────────────────────────────────────────────────
cat("Loading merged DGE ...\n")
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep_s <- samp$dataset %in% MEGA
dge  <- dge[, keep_s]
samp <- samp[keep_s]

# dge$samples$sex is 573/846 NA for the mega cohorts (GSE135251 + GSE213621 are
# sex-inferred only). The canonical COMPLETE inferred_sex (zero NA) lives in
# meta_matched.rds keyed by sample_id — the same source the dream mega-analysis
# uses. Join it here, or the model-based methods (DESeq2/dream/ComBat-seq) would
# drop the NA rows and crash on a count/design dimension mismatch.
meta_m <- as.data.frame(readRDS(file.path(INT_RESULTS, "meta_matched.rds")))
samp[, inferred_sex := meta_m$inferred_sex[match(sample_id, meta_m$sample_id)]]
n_na_sex <- sum(is.na(samp$inferred_sex))
if (n_na_sex > 0)
  stop(sprintf("inferred_sex join left %d NA — sample_id mismatch with meta_matched.rds", n_na_sex))
cat(sprintf("inferred_sex joined (complete): %s\n",
            paste(names(table(samp$inferred_sex)), table(samp$inferred_sex),
                  sep = "=", collapse = ", ")))

grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
cat(sprintf("Dataset: %d samples, %d genes, %d cohorts\n",
            ncol(logcpm), nrow(logcpm), length(unique(samp$dataset))))

# ── Shared helpers ───────────────────────────────────────────────────────────────
project_score <- function(logcpm, te_idx, genes, weights) {
  Xte <- logcpm[genes, te_idx, drop = FALSE]
  Zte <- t(scale(t(Xte)))
  Zte[!is.finite(Zte)] <- 0
  as.numeric(crossprod(weights, Zte)) / length(weights)
}

auc_from_score <- function(labels, scores) {
  if (length(unique(labels)) < 2 || all(is.na(scores))) return(NA_real_)
  suppressMessages(
    as.numeric(pROC::auc(pROC::roc(
      response = labels, predictor = scores,
      levels = c("Control", "Disease"), direction = "<", quiet = TRUE)))
  )
}

# ── Method fitting functions ─────────────────────────────────────────────────────
# Each returns data.table(gene, stat) — stat is SIGNED (Disease vs Control direction).
# Genes are ranked by |stat|; sign(stat) is used as projection weights.

fit_limma_voom <- function(dge, samp, tr) {
  grp_tr <- droplevels(factor(samp$group_binary[tr], levels = c("Control", "Disease")))
  ds_tr  <- droplevels(factor(samp$dataset[tr]))
  sex_tr <- droplevels(factor(samp$inferred_sex[tr]))
  terms  <- "ds_tr"
  if (nlevels(sex_tr) >= 2) terms <- c(terms, "sex_tr")
  terms  <- c(terms, "grp_tr")
  des_tr <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
  dge_tr <- edgeR::calcNormFactors(edgeR::DGEList(dge$counts[, tr]))
  v_tr   <- limma::voom(dge_tr, des_tr)
  fit_tr <- limma::eBayes(limma::lmFit(v_tr, des_tr))
  tt     <- limma::topTable(fit_tr, coef = "grp_trDisease", n = Inf, sort.by = "none")
  data.table(gene = rownames(tt), stat = tt$t)
}

fit_dream <- function(dge, samp, tr) {
  suppressPackageStartupMessages(library(variancePartition))
  meta_tr <- data.frame(
    group_binary = factor(samp$group_binary[tr], levels = c("Control", "Disease")),
    inferred_sex = factor(samp$inferred_sex[tr]),
    dataset      = factor(samp$dataset[tr]),
    row.names    = colnames(dge)[tr]
  )
  has_sex <- length(unique(meta_tr$inferred_sex)) >= 2
  form    <- if (has_sex)
    ~ group_binary + inferred_sex + (1 | dataset) else
    ~ group_binary + (1 | dataset)
  dge_tr  <- edgeR::calcNormFactors(edgeR::DGEList(dge$counts[, tr]))
  v_obj   <- variancePartition::voomWithDreamWeights(dge_tr, form, meta_tr,
                                                      BPPARAM = BPPARAM_SER)
  fit_obj <- variancePartition::dream(v_obj, form, meta_tr, BPPARAM = BPPARAM_SER)
  fit_obj <- variancePartition::eBayes(fit_obj)
  tt      <- limma::topTable(fit_obj, coef = "group_binaryDisease", n = Inf, sort.by = "none")
  data.table(gene = rownames(tt), stat = tt$t)
}

fit_deseq2 <- function(dge, samp, tr) {
  suppressPackageStartupMessages(library(DESeq2))
  counts_tr <- dge$counts[, tr]
  meta_tr <- data.frame(
    group_binary = factor(samp$group_binary[tr], levels = c("Control", "Disease")),
    inferred_sex = factor(samp$inferred_sex[tr]),
    dataset      = factor(samp$dataset[tr]),
    row.names    = colnames(counts_tr)
  )
  has_sex     <- length(unique(meta_tr$inferred_sex)) >= 2
  design_form <- if (has_sex) ~ dataset + inferred_sex + group_binary else ~ dataset + group_binary
  dds <- DESeq2::DESeqDataSetFromMatrix(counts_tr, meta_tr, design = design_form)
  suppressWarnings(
    dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_PAR, quiet = TRUE)
  )
  res <- DESeq2::results(dds, contrast = c("group_binary", "Disease", "Control"),
                         independentFiltering = FALSE)
  data.table(gene = rownames(res), stat = ifelse(is.na(res$stat), 0, res$stat))
}

fit_edger_qlf <- function(dge, samp, tr) {
  grp_tr <- droplevels(factor(samp$group_binary[tr], levels = c("Control", "Disease")))
  ds_tr  <- droplevels(factor(samp$dataset[tr]))
  sex_tr <- droplevels(factor(samp$inferred_sex[tr]))
  has_sex <- nlevels(sex_tr) >= 2
  terms   <- "ds_tr"
  if (has_sex) terms <- c(terms, "sex_tr")
  terms   <- c(terms, "grp_tr")
  des_tr  <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
  dge_tr  <- edgeR::calcNormFactors(edgeR::DGEList(dge$counts[, tr], group = grp_tr))
  dge_tr  <- edgeR::estimateDisp(dge_tr, des_tr)
  fit_tr  <- edgeR::glmQLFit(dge_tr, des_tr)
  qlf     <- edgeR::glmQLFTest(fit_tr, coef = "grp_trDisease")
  tt      <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  # sign(logFC) * sqrt(|F|) gives a signed t-like statistic for ranking
  data.table(gene = rownames(tt), stat = sign(tt$logFC) * sqrt(pmax(tt$F, 0)))
}

fit_metafor_re <- function(dge, samp, tr) {
  # Vectorized DL random-effects meta-analysis across training cohorts.
  # Equivalent to per-gene metafor::rma(method="DL") for gene ranking purposes.
  train_cohorts <- unique(samp$dataset[tr])

  per_cohort_dt <- rbindlist(lapply(train_cohorts, function(coh) {
    idx   <- tr[samp$dataset[tr] == coh]
    if (length(idx) < 4) return(NULL)
    grp_c <- droplevels(factor(samp$group_binary[idx], levels = c("Control", "Disease")))
    if (length(unique(grp_c)) < 2) return(NULL)
    sex_c <- droplevels(factor(samp$inferred_sex[idx]))
    has_sex <- nlevels(sex_c) >= 2
    terms   <- if (has_sex) c("sex_c", "grp_c") else "grp_c"
    des_c   <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
    dge_c   <- edgeR::calcNormFactors(edgeR::DGEList(dge$counts[, idx]))
    v_c     <- limma::voom(dge_c, des_c)
    fit_c   <- limma::eBayes(limma::lmFit(v_c, des_c))
    tt_c    <- limma::topTable(fit_c, coef = "grp_cDisease", n = Inf, sort.by = "none")
    # SE = |logFC / t|; guard against t ≈ 0
    se_est  <- abs(tt_c$logFC / pmax(abs(tt_c$t), 1e-10))
    data.table(gene = rownames(tt_c), cohort = coh, logFC = tt_c$logFC, se = se_est)
  }))

  if (is.null(per_cohort_dt) || nrow(per_cohort_dt) == 0) {
    return(data.table(gene = rownames(dge$counts), stat = 0))
  }

  # Pivot to genes × cohorts matrices
  dt_lfc  <- dcast(per_cohort_dt, gene ~ cohort, value.var = "logFC")
  dt_se   <- dcast(per_cohort_dt, gene ~ cohort, value.var = "se")
  coh_cols <- setdiff(names(dt_lfc), "gene")

  LFC_mat <- as.matrix(dt_lfc[, ..coh_cols])
  SE_mat  <- as.matrix(dt_se[, ..coh_cols])

  bad            <- !is.finite(SE_mat) | SE_mat <= 0 | !is.finite(LFC_mat)
  LFC_mat[bad]   <- NA_real_; SE_mat[bad] <- NA_real_
  K_valid        <- rowSums(!is.na(SE_mat))

  W_mat          <- 1 / SE_mat^2; W_mat[is.na(W_mat)] <- 0
  W_sum          <- rowSums(W_mat)
  wLFC           <- W_mat * LFC_mat; wLFC[is.na(wLFC)] <- 0
  mu_fe          <- rowSums(wLFC) / pmax(W_sum, 1e-10)
  Q              <- rowSums(W_mat * (LFC_mat - mu_fe)^2, na.rm = TRUE)
  C_             <- W_sum - rowSums(W_mat^2) / pmax(W_sum, 1e-10)
  tau2           <- pmax(0, (Q - pmax(K_valid - 1, 0)) / pmax(C_, 1e-10))

  W_re           <- 1 / (SE_mat^2 + tau2); W_re[is.na(W_re)] <- 0
  W_re_s         <- rowSums(W_re)
  wLFC_r         <- W_re * LFC_mat; wLFC_r[is.na(wLFC_r)] <- 0
  mu_re          <- rowSums(wLFC_r) / pmax(W_re_s, 1e-10)
  se_re          <- sqrt(1 / pmax(W_re_s, 1e-10))
  z_re           <- mu_re / pmax(se_re, 1e-10)
  z_re[K_valid < 2] <- 0

  data.table(gene = dt_lfc$gene, stat = z_re)
}

fit_combatseq_deseq2 <- function(dge, samp, tr) {
  suppressPackageStartupMessages({ library(sva); library(DESeq2) })
  counts_tr  <- dge$counts[, tr]
  batch_tr   <- samp$dataset[tr]
  group_tr   <- samp$group_binary[tr]
  cat("    ComBat-seq ...\n")
  counts_corr <- sva::ComBat_seq(counts_tr, batch = batch_tr, group = group_tr, full_mod = TRUE)
  counts_corr <- pmax(counts_corr, 0L)  # guard against rare negative edge cases
  meta_tr <- data.frame(
    group_binary = factor(samp$group_binary[tr], levels = c("Control", "Disease")),
    inferred_sex = factor(samp$inferred_sex[tr]),
    row.names    = colnames(counts_tr)
  )
  has_sex     <- length(unique(meta_tr$inferred_sex)) >= 2
  design_form <- if (has_sex) ~ inferred_sex + group_binary else ~ group_binary
  dds <- DESeq2::DESeqDataSetFromMatrix(counts_corr, meta_tr, design = design_form)
  suppressWarnings(
    dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_PAR, quiet = TRUE)
  )
  res <- DESeq2::results(dds, contrast = c("group_binary", "Disease", "Control"),
                         independentFiltering = FALSE)
  data.table(gene = rownames(res), stat = ifelse(is.na(res$stat), 0, res$stat))
}

# SVA+limma is intentionally excluded: SVA removes latent components that may
# include disease signal, producing a more aggressively deconfounded view than
# any method that explicitly includes batch as a fixed covariate. This makes
# AUROC comparisons across methods misleading. All included methods use dataset
# as an explicit term in their design/formula.

# Single-cohort fit for Part 2 (no dataset covariate)
fit_limma_voom_single <- function(dge, samp, tr) {
  grp_tr  <- droplevels(factor(samp$group_binary[tr], levels = c("Control", "Disease")))
  sex_tr  <- droplevels(factor(samp$inferred_sex[tr]))
  has_sex <- nlevels(sex_tr) >= 2
  terms   <- if (has_sex) c("sex_tr", "grp_tr") else "grp_tr"
  des_tr  <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
  dge_tr  <- edgeR::calcNormFactors(edgeR::DGEList(dge$counts[, tr]))
  v_tr    <- limma::voom(dge_tr, des_tr)
  fit_tr  <- limma::eBayes(limma::lmFit(v_tr, des_tr))
  tt      <- limma::topTable(fit_tr, coef = "grp_trDisease", n = Inf, sort.by = "none")
  data.table(gene = rownames(tt), stat = tt$t)
}

# ── Method registry ──────────────────────────────────────────────────────────────
FIT_METHODS <- list(
  limma_voom       = fit_limma_voom,
  dream            = fit_dream,
  deseq2           = fit_deseq2,
  edger_qlf        = fit_edger_qlf,
  metafor_re       = fit_metafor_re,
  combatseq_deseq2 = fit_combatseq_deseq2
)
METHOD_NAMES <- names(FIT_METHODS)

method_labels <- c(
  limma_voom       = "limma-voom",
  dream            = "dream",
  deseq2           = "DESeq2",
  edger_qlf        = "edgeR QLF",
  metafor_re       = "metafor RE",
  combatseq_deseq2 = "ComBat-seq+DESeq2"
)
method_colors <- setNames(RColorBrewer::brewer.pal(6, "Set2"), METHOD_NAMES)

# ── Part 1: LOCO sweep ───────────────────────────────────────────────────────────
cat("\n=== Part 1: LOCO sweep (5 folds x 7 methods x 7 N + 100 null reps) ===\n")

# Fixed gene universe for the random null (method-agnostic: random genes, no sign).
GENE_UNIV <- rownames(logcpm)

results_loco <- rbindlist(lapply(MEGA, function(test_cohort) {
  cat(sprintf("\n-- Held-out: %s (%s) --\n", test_cohort, cohort_short[test_cohort]))
  tr <- which(samp$dataset != test_cohort)
  te <- which(samp$dataset == test_cohort)
  if (length(unique(grp[te])) < 2) { cat("  Skipped (single class in held-out)\n"); return(NULL) }

  # --- Random null: computed ONCE per fold (does not depend on method) ---------
  cat("  [null] sampling random gene sets ...\n")
  null_rows <- rbindlist(lapply(N_VALS, function(N) {
    n_draw <- min(N, length(GENE_UNIV))
    rbindlist(lapply(seq_len(N_RAND_REPS), function(r) {
      set.seed(r * 10000L + N)               # collision-free: max(N)=2000 < 10000
      rand_g <- sample(GENE_UNIV, n_draw)
      sc     <- project_score(logcpm, te, rand_g, rep(1L, n_draw))
      data.table(method = NA_character_, held_out = test_cohort, N = N, type = "random",
                 auroc = auc_from_score(grp[te], sc), rep = r)
    }))
  }))

  # --- Real signatures: top-N by |stat| per method -----------------------------
  real_rows <- rbindlist(lapply(METHOD_NAMES, function(m) {
    cat(sprintf("  [%s] fitting ...\n", m))
    t0 <- proc.time()["elapsed"]
    stats_dt <- tryCatch(
      FIT_METHODS[[m]](dge, samp, tr),
      error = function(e) {
        message("  ERROR in ", m, " (", test_cohort, "): ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(stats_dt)) return(NULL)
    cat(sprintf("  [%s] done in %.1f s, %d genes\n",
                m, proc.time()["elapsed"] - t0, nrow(stats_dt)))

    n_genes   <- nrow(stats_dt)
    stats_ord <- stats_dt[order(-abs(stat))]  # pre-sorted by |stat|

    rbindlist(lapply(N_VALS, function(N) {
      top_genes <- head(stats_ord$gene, min(N, n_genes))
      w   <- sign(stats_ord$stat[seq_along(top_genes)])
      sc  <- project_score(logcpm, te, top_genes, w)
      data.table(method = m, held_out = test_cohort, N = N, type = "real",
                 auroc = auc_from_score(grp[te], sc), rep = NA_integer_)
    }))
  }))

  rbind(real_rows, null_rows)
}))

fwrite(results_loco, file.path(OUT, "auroc_loco_sweep_data.csv"))
cat("\nSaved auroc_loco_sweep_data.csv\n")

# ── Part 2: Train-on-1, project-to-4 ─────────────────────────────────────────────
cat("\n=== Part 2: Train-on-1 cohort, project-to-4 (limma-voom) ===\n")

results_train1 <- rbindlist(lapply(MEGA, function(train_cohort) {
  cat(sprintf("\n-- Training: %s (%s) --\n", train_cohort, cohort_short[train_cohort]))
  tr <- which(samp$dataset == train_cohort)
  if (length(unique(grp[tr])) < 2) { cat("  Skipped (single class)\n"); return(NULL) }

  stats_dt <- tryCatch(
    fit_limma_voom_single(dge, samp, tr),
    error = function(e) { message("ERROR: ", conditionMessage(e)); return(NULL) }
  )
  if (is.null(stats_dt)) return(NULL)
  stats_ord <- stats_dt[order(-abs(stat))]

  rbindlist(lapply(N_VALS, function(N) {
    top_genes <- head(stats_ord$gene, min(N, nrow(stats_ord)))
    w         <- sign(stats_ord$stat[seq_along(top_genes)])
    rbindlist(lapply(setdiff(MEGA, train_cohort), function(test_cohort) {
      te <- which(samp$dataset == test_cohort)
      if (length(unique(grp[te])) < 2) return(NULL)
      sc <- project_score(logcpm, te, top_genes, w)
      data.table(train_cohort = train_cohort, test_cohort = test_cohort, N = N,
                 auroc = auc_from_score(grp[te], sc))
    }))
  }))
}))

fwrite(results_train1, file.path(OUT, "auroc_train1_project4_data.csv"))
cat("\nSaved auroc_train1_project4_data.csv\n")

# ── Figures ───────────────────────────────────────────────────────────────────────
cat("\n=== Generating figures ===\n")

# Ordered factor levels for methods and cohorts
m_lvls  <- unname(method_labels[METHOD_NAMES])
c_lvls  <- unname(cohort_short[MEGA])
col_map <- setNames(method_colors[METHOD_NAMES], m_lvls)

loco_real <- results_loco[type == "real"]
loco_null <- results_loco[type == "random"]

loco_real[, method_lbl := factor(method_labels[method], levels = m_lvls)]
loco_real[, cohort_lbl := factor(cohort_short[held_out], levels = c_lvls)]

# Summaries
loco_sum <- loco_real[, .(auroc_mean = mean(auroc, na.rm = TRUE),
                           auroc_sd  = sd(auroc,   na.rm = TRUE)),
                      by = .(method, method_lbl, N)]
null_sum <- loco_null[, .(null_mean = mean(auroc, na.rm = TRUE),
                           null_sd  = sd(auroc,   na.rm = TRUE)), by = N]

# ── Panel 1: AUROC vs N line plot ─────────────────────────────────────────────────
p1 <- ggplot() +
  geom_ribbon(data = null_sum,
              aes(x = N, ymin = null_mean - 2 * null_sd, ymax = null_mean + 2 * null_sd),
              fill = "grey80", alpha = 0.5) +
  geom_line(data = null_sum, aes(x = N, y = null_mean),
            color = "grey55", linewidth = 0.4, linetype = "dashed") +
  geom_line(data = loco_sum,
            aes(x = N, y = auroc_mean, color = method_lbl), linewidth = 0.6) +
  geom_point(data = loco_sum,
             aes(x = N, y = auroc_mean, color = method_lbl), size = 1.4) +
  geom_errorbar(data = loco_sum,
                aes(x = N, ymin = auroc_mean - auroc_sd, ymax = auroc_mean + auroc_sd,
                    color = method_lbl),
                width = 0.06, linewidth = 0.4) +
  scale_x_log10(breaks = N_VALS, labels = N_VALS) +
  scale_color_manual(values = col_map, name = NULL) +
  coord_cartesian(ylim = c(0.45, 1.0)) +
  labs(x = "Gene signature size (N, log scale)",
       y = "AUROC — 5-fold LOCO mean ± SD",
       title = "Supervised disease separation: method × signature size") +
  theme_masld(base_size = 7) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "auroc_vs_N.pdf"), p1,
       width = 8.4, height = 4.2, device = cairo_pdf)
cat("  auroc_vs_N.pdf\n")

# ── Panel 2: Method x N heatmap ──────────────────────────────────────────────────
loco_sum[, N_lbl := factor(N, levels = N_VALS)]

p2 <- ggplot(loco_sum, aes(x = N_lbl, y = method_lbl, fill = auroc_mean)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", auroc_mean)), size = 2.2, color = "grey15") +
  scale_fill_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#08306b",
                       midpoint = 0.75, limits = c(0.5, 1.0), na.value = "grey90",
                       name = "Mean\nAUROC") +
  labs(x = "Gene signature size (N)", y = NULL,
       title = "Mean AUROC (5-fold LOCO) by method and N") +
  theme_masld(base_size = 7) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "method_N_heatmap.pdf"), p2,
       width = 7.0, height = 3.2, device = cairo_pdf)
cat("  method_N_heatmap.pdf\n")

# ── Panel 3: Per-cohort breakdown at N=200 ────────────────────────────────────────
p3_dt <- loco_real[N == 200]
p3 <- ggplot(p3_dt, aes(x = method_lbl, y = auroc, color = method_lbl)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_point(size = 2.2) +
  facet_wrap(~ cohort_lbl, nrow = 1) +
  scale_color_manual(values = col_map) +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(x = NULL, y = "AUROC at N = 200", color = NULL,
       title = "Per held-out cohort AUROC at N = 200") +
  theme_masld(base_size = 7) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 6.5))
ggsave(file.path(OUT, "per_cohort_N200.pdf"), p3,
       width = 9.0, height = 3.5, device = cairo_pdf)
cat("  per_cohort_N200.pdf\n")

# ── Panel 4: Null distribution vs real methods ───────────────────────────────────
# Null is method-agnostic (computed once per fold over the full gene universe).
null_p4 <- loco_null   # all 5 folds x 7 N x 100 reps

p4 <- ggplot() +
  geom_boxplot(data = null_p4,
               aes(x = factor(N), y = auroc),
               fill = "grey88", color = "grey50",
               outlier.size = 0.3, linewidth = 0.4, width = 0.55) +
  geom_jitter(data = loco_real,
              aes(x = factor(N), y = auroc, color = method_lbl),
              position = position_jitter(width = 0.18, seed = 42),
              size = 0.9, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = col_map, name = NULL) +
  coord_cartesian(ylim = c(0.35, 1.0)) +
  labs(x = "Gene signature size (N)",
       y = "AUROC",
       title = "Random null (grey box) vs. real method AUROCs (coloured points)") +
  theme_masld(base_size = 7) +
  theme(legend.position = "right")
ggsave(file.path(OUT, "null_vs_real.pdf"), p4,
       width = 9.0, height = 4.0, device = cairo_pdf)
cat("  null_vs_real.pdf\n")

# ── Panel 5: Train-on-1, project-to-4 heatmap (N=200) ────────────────────────────
p5_dt <- results_train1[N == 200]
p5_dt[, train_lbl := factor(cohort_short[train_cohort], levels = c_lvls)]
p5_dt[, test_lbl  := factor(cohort_short[test_cohort],  levels = c_lvls)]

p5 <- ggplot(p5_dt, aes(x = test_lbl, y = train_lbl, fill = auroc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", auroc)), size = 2.4, color = "grey15") +
  scale_fill_gradient(low = "#fff5eb", high = "#7f2704",
                      limits = c(0.5, 1.0), na.value = "grey92",
                      name = "AUROC") +
  labs(x = "Test cohort (held-out)", y = "Training cohort",
       title = "Train-on-1-cohort, project-to-4 (limma-voom, N = 200)") +
  theme_masld(base_size = 7) +
  theme(panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "train1_project4.pdf"), p5,
       width = 4.5, height = 4.0, device = cairo_pdf)
cat("  train1_project4.pdf\n")

cat(sprintf("\nDone. All outputs written to:\n  %s\n", OUT))
