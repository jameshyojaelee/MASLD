#!/usr/bin/env Rscript
# figS_multimethod_disease_separation.R
# Supervised disease separation — the companion to the batch-model PCA
# (batch_model_pca.pdf), which does NOT separate Control/Disease because disease
# is <1% of variance.
# This panel shows the disease contrast IS recoverable on a SUPERVISED axis, with
# NO double-dipping: a leave-one-cohort-out (LOCO) disease signature is trained on
# 4 cohorts and projected onto the held-out 5th; the held-out samples never inform
# their own signature or scaling.
#   supervised_disease.pdf  (was panelK_supervised_disease.pdf)
#   left  : held-out LOCO disease score per cohort, coloured by true label
#   right : pooled held-out ROC (+ per-cohort AUROC)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats); library(pROC)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"; set.seed(42)

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
                  GSE162694 = "GSE162694", GSE213621 = "GSE213621")
TOPK <- 200L            # signature size (genes by |t| on the training cohorts)

dge <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep_s <- samp$dataset %in% MEGA
dge <- dge[, keep_s]; samp <- samp[keep_s]
grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds  <- factor(samp$dataset)
sex <- factor(samp$inferred_sex)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
cat(sprintf("LOCO disease separation: %d samples, %d genes, %d cohorts\n",
            ncol(logcpm), nrow(logcpm), nlevels(ds)))

# --- leave-one-cohort-out: train signature on 4, project onto held-out 5th ----
scores <- rbindlist(lapply(MEGA, function(test) {
  tr <- which(samp$dataset != test); te <- which(samp$dataset == test)
  # need both classes in the held-out cohort to evaluate separation
  if (length(unique(grp[te])) < 2) return(NULL)
  # train DE (cohort-adjusted), supervised on TRAIN ONLY
  ds_tr <- droplevels(ds[tr]); sex_tr <- droplevels(sex[tr]); grp_tr <- droplevels(grp[tr])
  terms <- c("ds_tr", if (nlevels(sex_tr) >= 2) "sex_tr", "grp_tr")
  des_tr <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
  v_tr   <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(dge$counts[, tr])), des_tr)
  fit_tr <- limma::eBayes(limma::lmFit(v_tr, des_tr))
  tt     <- limma::topTable(fit_tr, coef = "grp_trDisease", n = Inf, sort.by = "none")
  tt     <- tt[order(-abs(tt$t)), ]
  sig    <- head(rownames(tt), TOPK); w <- sign(tt$t[match(sig, rownames(tt))])
  # project onto held-out cohort, scaled WITHIN the held-out cohort (no leakage)
  Xte <- logcpm[sig, te, drop = FALSE]
  Zte <- t(scale(t(Xte)))                       # z within held-out cohort
  Zte[!is.finite(Zte)] <- 0
  score <- as.numeric(crossprod(w, Zte)) / length(w)
  data.table(sample_id = samp$sample_id[te], cohort = cohort_short[test],
             disease = grp[te], score = score)
}))
scores[, cohort := factor(cohort, levels = unname(cohort_short))]

# per-cohort + pooled AUROC (held-out)
auc_by <- scores[, .(auc = as.numeric(pROC::auc(pROC::roc(
  response = disease, predictor = score, levels = c("Control", "Disease"),
  direction = "<", quiet = TRUE)))), by = cohort]
roc_all <- pROC::roc(response = scores$disease, predictor = scores$score,
                     levels = c("Control", "Disease"), direction = "<", quiet = TRUE)
auc_all <- as.numeric(pROC::auc(roc_all))
cat("Pooled held-out AUROC =", round(auc_all, 3), "\n"); print(auc_by)

# --- Panel K-left: held-out score per cohort ---------------------------------
lab_auc <- auc_by[, setNames(sprintf("%s\nAUROC %.2f", cohort, auc), cohort)]
pL <- ggplot(scores, aes(cohort, score, fill = disease, colour = disease)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_violin(position = position_dodge(width = 0.8), width = 0.78, alpha = 0.35,
              linewidth = 0.3, scale = "width", draw_quantiles = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             size = 0.3, alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c(Control = CTRL, Disease = masld_colors$nash), name = NULL) +
  scale_colour_manual(values = c(Control = CTRL, Disease = masld_colors$nash), name = NULL) +
  scale_x_discrete(labels = lab_auc) +
  labs(x = NULL, y = "held-out LOCO disease score") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(size = 6, lineheight = 0.85),
        legend.position = "top")

# --- Panel K-right: pooled held-out ROC --------------------------------------
roc_dt <- data.table(fpr = 1 - roc_all$specificities, tpr = roc_all$sensitivities)
setorder(roc_dt, fpr, tpr)
pR <- ggplot(roc_dt, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_path(colour = masld_colors$nash, linewidth = 0.7) +
  annotate("text", x = 0.6, y = 0.18, size = 2.8, colour = "grey15",
           label = sprintf("pooled held-out\nAUROC = %.2f", auc_all)) +
  coord_equal() +
  labs(x = "false positive rate", y = "true positive rate") +
  theme_masld(base_size = 7)

fig <- (pL | pR) + plot_layout(widths = c(1.5, 1)) +
  plot_annotation(
    title = "Leave-one-cohort-out supervised disease score",
    theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(OUT, "supervised_disease.pdf"), fig,
       width = 8.4, height = 4.2, device = cairo_pdf)
fwrite(scores, file.path(OUT, "supervised_disease_data.csv"))
fwrite(auc_by, file.path(OUT, "supervised_disease_auc.csv"))
cat("Wrote supervised_disease.pdf\n")
