#!/usr/bin/env Rscript
# figS_raw_vs_shrunk_patient_concordance.R
# ---------------------------------------------------------------------------
# Raw-log2FC vs ashr-shrunk-log2FC comparison at the PATIENT level, recomputed
# from source over an effect-size cutoff grid extending to 3.0, WITH a healthy-
# control baseline (the negative control: controls are the within-cohort reference
# so their LFC ~ 0 -> rho ~ 0, directional concordance ~ 50%).
#
# 2026-06-24: the effect-size source is the LIBRARY's canonical limma-voom QW C2
# DEGs (canonical_deg_results.csv), swapped off the retired metafor meta_results_ashr.csv
# so human AND mouse DE both use the same limma-voom engine. This validates the
# raw-vs-ashr shrinkage on the SAME instrument the v8 library is built from.
#
# Inputs:
#   canonical_deg_results.csv  (gene, padj, logFC, lfsr, shrunk_logFC; limma-voom QW C2 + ashr)
#   patient_lfc_matrix.csv.gz  (690 disease samples; log2FC vs within-cohort control mean)
#   merged_dge.rds             (rebuild the 156 healthy-control LFC the same way)
# Conditions (significance fixed at 0.05; effect-size cutoff x swept 0 -> 3.0):
#   RAW   : padj < 0.05 & |log2FC|        > x
#   ASHR  : lfsr < 0.05 & |shrunk_log2FC| > x
#   TREAT : limma TREAT FDR < 0.05 with lfc = x
# Metrics (median over samples), two-sided DEG sets:
#   - median per-sample directional concordance %  (Disease vs Healthy)
#   - median per-sample Spearman rho (sample LFC vs dream effect)  (Disease vs Healthy)
#   - n_DEGs (two-sided) and raw/ashr DEG-set Jaccard (DEG-set level, no group split)
#
# Output (library dir): Cas13_Library_Design/figures/
#   08a_patient_concordance_rho.pdf, _pct.pdf, 08c_deg_counts_by_cutoff.pdf, 08d_raw_ashr_jaccard.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(edgeR); library(limma); library(yaml) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
FIGLIB <- file.path(BASE, "Cas13_Library_Design/figures")
dir.create(FIGLIB, showWarnings = FALSE, recursive = TRUE)
strip_v <- function(x) sub("[.][0-9]+$", "", x)

dream <- fread(file.path(INT, "canonical_deg_results.csv"))  # canonical limma-voom QW C2 (2026-06-24: was metafor meta_results_ashr.csv); gene,padj,logFC,lfsr,shrunk_logFC
lfc_mat <- fread(file.path(INT, "patient_lfc_matrix.csv.gz")) # disease samples

# ---- healthy-control LFC matrix (same within-cohort reference as disease) ----
message("Building healthy-control LFC matrix from merged_dge.rds ...")
dge <- readRDS(file.path(INT, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge  <- dge[, dge$samples$dataset %in% mega]
lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
sinfo <- data.table(sample_id = colnames(dge), dataset = dge$samples$dataset,
                    group = as.character(dge$samples$group_binary))
valid_ds <- sinfo[, .(nc = sum(group == "Control"), nd = sum(group == "Disease")),
                  by = dataset][nc > 0 & nd > 0, dataset]
cmeans <- lapply(setNames(valid_ds, valid_ds), function(ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  if (length(ids) == 1) lcpm[, ids] else rowMeans(lcpm[, ids, drop = FALSE])
})
ctrl_sids <- sinfo[group == "Control" & dataset %in% valid_ds, sample_id]
ctrl_lfc <- matrix(NA_real_, nrow = nrow(lcpm), ncol = length(ctrl_sids),
                   dimnames = list(rownames(lcpm), ctrl_sids))
for (ds in valid_ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  for (sid in ids) ctrl_lfc[, sid] <- lcpm[, sid] - cmeans[[ds]]
}

# ---- common gene set across dream + disease + control matrices ---------------
sample_ids <- setdiff(colnames(lfc_mat), "gene")
common <- Reduce(intersect, list(dream$gene, lfc_mat$gene, rownames(ctrl_lfc)))
if (length(common) < 10000) {
  dream[, gene := strip_v(gene)]; lfc_mat[, gene := strip_v(gene)]
  rownames(ctrl_lfc) <- strip_v(rownames(ctrl_lfc))
  common <- Reduce(intersect, list(dream$gene, lfc_mat$gene, rownames(ctrl_lfc)))
}
dream_c <- dream[match(common, gene)]
setkey(lfc_mat, gene)
lfc_dis <- as.matrix(lfc_mat[common, sample_ids, with = FALSE])  # disease genes x samples
ctrl    <- ctrl_lfc[common, , drop = FALSE]                       # control genes x samples
cat(sprintf("common genes: %d | disease: %d | healthy: %d\n", length(common), ncol(lfc_dis), ncol(ctrl)))

# per-sample metrics on matrix M for one (significance, effect) pair at cutoff x
eval_grp <- function(cut, sigvec, effvec, M) {
  up <- which(sigvec < 0.05 & effvec >  cut)
  dn <- which(sigvec < 0.05 & effvec < -cut)
  idx <- c(up, dn); n <- length(idx)
  if (n == 0) return(list(n_degs = 0L, median_rho = NA_real_, median_pct = NA_real_))
  nc <- numeric(ncol(M))
  if (length(up)) nc <- nc + colSums(M[up, , drop = FALSE] > 0)
  if (length(dn)) nc <- nc + colSums(M[dn, , drop = FALSE] < 0)
  pct <- nc / n * 100
  if (n >= 2) rho <- as.numeric(cor(rank(effvec[idx]), apply(M[idx, , drop = FALSE], 2, rank)))
  else        rho <- rep(NA_real_, ncol(M))
  list(n_degs = n, median_rho = median(rho, na.rm = TRUE), median_pct = median(pct))
}

CUTS <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0)
RAW   <- "raw:   padj<0.05 & |log2FC| > x"
ASHR  <- "ashr:  lfsr<0.05 & |shrunk| > x"
TREAT <- "treat: FDR<0.05 with lfc = x"
DZ <- sprintf("Disease (n=%d)", ncol(lfc_dis))
HZ <- sprintf("Healthy (n=%d)", ncol(ctrl))
build <- function(cond, sigvec, effvec) rbindlist(lapply(CUTS, function(x) {
  dz <- eval_grp(x, sigvec, effvec, lfc_dis); hz <- eval_grp(x, sigvec, effvec, ctrl)
  rbind(data.table(condition = cond, cutoff = x, group = DZ, n_degs = dz$n_degs,
                   median_rho = dz$median_rho, median_pct = dz$median_pct),
        data.table(condition = cond, cutoff = x, group = HZ, n_degs = hz$n_degs,
                   median_rho = hz$median_rho, median_pct = hz$median_pct))
}))
df_probe <- dream_c[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
df_idx <- unique(round(seq(1, nrow(df_probe), length.out = min(nrow(df_probe), 12))))
treat_df <- median(vapply(
  df_idx,
  function(i) {
    uniroot(function(df) 2 * pt(-abs(df_probe$t[i]), df = df) - df_probe$P.Value[i],
            c(0.1, 1e6))$root
  },
  numeric(1)
))
treat_counts <- rbindlist(lapply(CUTS, function(x) {
  tp <- pt((abs(dream_c$logFC) - x) / dream_c$SE, df = treat_df, lower.tail = FALSE) +
    pt((abs(dream_c$logFC) + x) / dream_c$SE, df = treat_df, lower.tail = FALSE)
  sig <- p.adjust(tp, method = "BH") < 0.05
  data.table(
    condition = TREAT,
    cutoff = x,
    group = DZ,
    n_degs = sum(sig),
    median_rho = NA_real_,
    median_pct = NA_real_
  )
}))
res <- rbind(
  build(RAW, dream_c$padj, dream_c$logFC),
  build(ASHR, dream_c$lfsr, dream_c$shrunk_logFC),
  treat_counts
)
res[, condition := factor(condition, levels = c(RAW, ASHR, TREAT))]
res[, group := factor(group, levels = c(DZ, HZ))]
pal <- setNames(c("#B2182B", "#3B4CC0", "#1B9E77"), c(RAW, ASHR, TREAT))
fwrite(res, file.path(BASE, "Cas13_Library_Design/data/patient_concordance_by_cutoff_to3.csv"))
fwrite(
  treat_counts[, .(treat_lfc = cutoff, fdr = 0.05, total = n_degs)],
  file.path(BASE, "Cas13_Library_Design/data/treat_deg_counts.csv")
)
cat("\n=== sanity (cutoff 0.5; limma-voom C2 source -- expect Disease rho>Healthy rho, ashr>=raw concordance) ===\n")
print(res[cutoff == 0.5])

bt <- theme_bw(base_size = 6) + theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "top", legend.box = "vertical")
ltys <- setNames(c("solid", "dotted"), c(DZ, HZ))
mk <- function(y, ytitle, title, ylab_pct = FALSE) {
  g <- ggplot(res[!is.na(get(y))], aes(cutoff, get(y), color = condition, linetype = group)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.7) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = ltys) +
    scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
    labs(x = "Effect-size cutoff", y = ytitle) + bt
  if (ylab_pct) g <- g + scale_y_continuous(labels = function(z) paste0(z, "%"))
  g
}
pRho <- mk("median_rho", "Median per-sample Spearman rho",
           "Per-sample concordance: disease vs healthy")
pPct <- mk("median_pct", "Median per-sample concordance (%)",
           "Per-sample directional concordance: disease vs healthy", ylab_pct = TRUE)

# ---- DEG count + Jaccard (DEG-set level; no disease/healthy split) -----------
ndeg <- unique(res[group == DZ, .(condition, cutoff, n_degs)])
pNdeg <- ggplot(ndeg, aes(cutoff, n_degs, color = condition)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.0) +
  scale_color_manual(values = pal) + scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Effect-size cutoff", y = "Total DEGs (up + down)") + bt
jac <- rbindlist(lapply(CUTS, function(x) {
  r <- dream_c[padj < 0.05 & abs(logFC) > x, gene]
  a <- dream_c[lfsr < 0.05 & abs(shrunk_logFC) > x, gene]
  u <- length(union(r, a)); i <- length(intersect(r, a))
  data.table(cutoff = x, jaccard = if (u > 0) i / u else NA_real_)
}))
fwrite(jac, file.path(BASE, "Cas13_Library_Design/data/raw_ashr_deg_jaccard_by_cutoff.csv"))
pJac <- ggplot(jac[!is.na(jaccard)], aes(cutoff, jaccard)) +
  geom_line(linewidth = 0.8, color = "#333333") + geom_point(size = 2.0, color = "#333333") +
  scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
  scale_y_continuous(labels = function(z) paste0(round(z * 100), "%"), limits = c(0, 1)) +
  labs(x = "Effect-size cutoff", y = "Jaccard overlap of DEG sets") + bt

message("[caption] Per-sample concordance: disease vs healthy (rho)")
message("[caption] Per-sample directional concordance: disease vs healthy (%)")
message("[caption] Total DEG count by cutoff")
message("[caption] DEG-set overlap between raw and ashr")
ggsave(file.path(FIGLIB, "patient_concordance_rho.pdf"), pRho, width = 7.0, height = 4.6, useDingbats = FALSE)
ggsave(file.path(FIGLIB, "patient_concordance_pct.pdf"), pPct, width = 7.0, height = 4.6, useDingbats = FALSE)
ggsave(file.path(FIGLIB, "deg_counts_by_cutoff.pdf"),              pNdeg, width = 6.8, height = 4.2, useDingbats = FALSE)
ggsave(file.path(FIGLIB, "raw_ashr_jaccard.pdf"),        pJac, width = 6.8, height = 4.2, useDingbats = FALSE)
cat("\nWrote S_lib_9 panels (rho/pct now with disease vs healthy) to", FIGLIB, "\n")
