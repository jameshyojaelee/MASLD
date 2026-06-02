#!/usr/bin/env Rscript
# 05_concordance.R
# ---------------------------------------------------------------------------
# Aggregates dream (canonical), edgeR-QL, voomLmFit, metafor-HKSJ results
# and computes:
#   - Spearman rho / Pearson r matrices on disease t-stats and logFC
#   - Pairwise Jaccard on padj<0.05 DEG sets
#   - Direction concordance %
#   - Per-gene robustness flag (5 levels)
# Emits atlas_columns.tsv (gene x 13 cols) for 27a_assemble_evidence_atlas.R
# and panel_data_A_to_F.rds for the figure renderer.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(yaml)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR_IN  <- file.path(INT, "results/mega_validation")
RDIR_OUT <- file.path(INT, "results/mega_validation/concordance")
dir.create(RDIR_OUT, recursive = TRUE, showWarnings = FALSE)

# --- Load canonical dream + all sensitivity arms ---
dream <- fread(file.path(INT, "results/integration/dream_results.csv"))
setnames(dream, c("logFC","t","P.Value","padj"),
                c("dream_logFC","dream_t","dream_P","dream_padj"))

eql <- fread(file.path(RDIR_IN, "edgeRql/ql_results.csv"))
setnames(eql, c("logFC","F","P.Value","padj"),
              c("eqlogFC","eqF","eqP","eqpadj"))
# edgeR returns F; convert to signed-sqrt(F) as a t-stat proxy (1-df contrast).
eql[, eqt := sign(eqlogFC) * sqrt(eqF)]

vlm <- fread(file.path(RDIR_IN, "voomLmFit/vlm_results.csv"))
setnames(vlm, c("logFC","t","P.Value","padj"),
              c("vlmlogFC","vlmt","vlmP","vlmpadj"))

mfr <- fread(file.path(RDIR_IN, "metafor/per_gene_meta.csv"))
setnames(mfr, c("meta_logFC","meta_z","HKSJ_p","HKSJ_padj"),
              c("mfrlogFC","mfrt","mfrP","mfrpadj"))

# DESeq2 (Arm 5; added 2026-05-19). Tolerate absence for older runs.
deseq_path <- file.path(RDIR_IN, "DESeq2/deseq2_results.csv")
deseq_present <- file.exists(deseq_path)
if (deseq_present) {
  deseq <- fread(deseq_path)
  setnames(deseq, c("logFC","t","P.Value","padj"),
                  c("dslogFC","dst","dsP","dspadj"))
  cat("DESeq2 arm present:", nrow(deseq), "genes\n")
} else {
  cat("DESeq2 arm absent (skipping)\n")
}

# --- Extended arms (added 2026-05-19): trend, DESeq2-LRT, edgeR-LRT, voom+SVA
# These enter rho/Jaccard/UpSet panels but NOT the robustness-flag voting,
# since they are robustness CHECKS on the 3 canonical voters (eql/vlm/deseq).
load_extended <- function(path, prefix) {
  if (!file.exists(path)) {
    cat("  ", basename(path), " absent\n"); return(NULL)
  }
  dt <- fread(path)
  cn <- names(dt)
  # Cover the four output schemas: limma (t), DESeq2 (stat), edgeR-LRT (LR
  # statistic → convert to signed-sqrt(LR) as a t-stat proxy)
  if ("t" %in% cn && "P.Value" %in% cn) {
    setnames(dt, c("logFC","t","P.Value","padj"),
                  paste0(prefix, c("logFC","t","P","padj")))
  } else if ("stat" %in% cn && "pvalue" %in% cn) {
    setnames(dt, c("log2FoldChange","stat","pvalue","padj"),
                  paste0(prefix, c("logFC","t","P","padj")))
  } else if ("LR" %in% cn && "P.Value" %in% cn) {
    # edgeR LRT: 1-df LRT → signed-sqrt(LR) is a t-stat proxy with sign of logFC
    dt[, signed_sqrt_LR := sign(logFC) * sqrt(LR)]
    setnames(dt, c("logFC","signed_sqrt_LR","P.Value","padj"),
                  paste0(prefix, c("logFC","t","P","padj")))
  }
  cat("  ", prefix, " arm present:", nrow(dt), "genes\n")
  dt
}
trend  <- load_extended(file.path(RDIR_IN, "limma_trend/trend_results.csv"),  "tr")
dlrt   <- load_extended(file.path(RDIR_IN, "DESeq2_LRT/deseq2_lrt_results.csv"), "dl")
elrt   <- load_extended(file.path(RDIR_IN, "edgeR_LRT/edger_lrt_results.csv"),  "el")
sva_a  <- load_extended(file.path(RDIR_IN, "voom_sva/voom_sva_results.csv"),    "sv")
trend_present <- !is.null(trend)
dlrt_present  <- !is.null(dlrt)
elrt_present  <- !is.null(elrt)
sva_present   <- !is.null(sva_a)

# --- Left-join: dream is the anchor; arm columns become NA where genes filter out ---
merged <- dream[, .(gene, dream_logFC, dream_t, dream_P, dream_padj)]
merged <- merge(merged, eql[, .(gene, eqlogFC, eqt, eqP, eqpadj)],
                by = "gene", all.x = TRUE)
merged <- merge(merged, vlm[, .(gene, vlmlogFC, vlmt, vlmP, vlmpadj)],
                by = "gene", all.x = TRUE)
merged <- merge(merged, mfr[, .(gene, mfrlogFC, mfrt, mfrP, mfrpadj,
                                 meta_SE, tau2, I2, n_studies)],
                by = "gene", all.x = TRUE)
if (deseq_present) {
  merged <- merge(merged, deseq[, .(gene, dslogFC, dst, dsP, dspadj)],
                  by = "gene", all.x = TRUE)
}
merge_ext <- function(merged, ext, prefix) {
  if (is.null(ext)) return(merged)
  cols <- c("gene", paste0(prefix, c("logFC","t","P","padj")))
  merge(merged, ext[, ..cols], by = "gene", all.x = TRUE)
}
merged <- merge_ext(merged, trend, "tr")
merged <- merge_ext(merged, dlrt,  "dl")
merged <- merge_ext(merged, elrt,  "el")
merged <- merge_ext(merged, sva_a, "sv")
cat("Genes from dream anchor:", nrow(merged), "\n")
cat("  with edgeR-QL:",  sum(!is.na(merged$eqpadj)),  "\n")
cat("  with voomLmFit:", sum(!is.na(merged$vlmpadj)), "\n")
cat("  with metafor:",   sum(!is.na(merged$mfrpadj)), "\n")
if (deseq_present) cat("  with DESeq2:",  sum(!is.na(merged$dspadj)), "\n")
if (trend_present) cat("  with limma-trend:", sum(!is.na(merged$trpadj)), "\n")
if (dlrt_present)  cat("  with DESeq2-LRT:",  sum(!is.na(merged$dlpadj)), "\n")
if (elrt_present)  cat("  with edgeR-LRT:",   sum(!is.na(merged$elpadj)), "\n")
if (sva_present)   cat("  with voom+SVA:",    sum(!is.na(merged$svpadj)), "\n")

# --- Concordance matrices (grows with available arms; 4-8+) ---
tstat_cols <- c("dream_t","eqt","vlmt","mfrt")
lfc_cols   <- c("dream_logFC","eqlogFC","vlmlogFC","mfrlogFC")
if (deseq_present) { tstat_cols <- c(tstat_cols, "dst"); lfc_cols <- c(lfc_cols, "dslogFC") }
if (trend_present) { tstat_cols <- c(tstat_cols, "trt"); lfc_cols <- c(lfc_cols, "trlogFC") }
if (dlrt_present)  { tstat_cols <- c(tstat_cols, "dlt"); lfc_cols <- c(lfc_cols, "dllogFC") }
if (elrt_present)  { tstat_cols <- c(tstat_cols, "elt"); lfc_cols <- c(lfc_cols, "ellogFC") }
if (sva_present)   { tstat_cols <- c(tstat_cols, "svt"); lfc_cols <- c(lfc_cols, "svlogFC") }
tstat_mat <- as.matrix(merged[, ..tstat_cols])
lfc_mat   <- as.matrix(merged[, ..lfc_cols])
rho_t <- cor(tstat_mat, method = "spearman", use = "pairwise.complete.obs")
r_lfc <- cor(lfc_mat,   method = "pearson",  use = "pairwise.complete.obs")
fwrite(as.data.table(rho_t, keep.rownames = "method"),
       file.path(RDIR_OUT, "tstat_spearman_rho.csv"))
fwrite(as.data.table(r_lfc, keep.rownames = "method"),
       file.path(RDIR_OUT, "logFC_pearson_r.csv"))

# --- Jaccard on padj < 0.05 DEG sets ---
deg_sets <- list(
  dream = merged[!is.na(dream_padj) & dream_padj < 0.05, gene],
  eq    = merged[!is.na(eqpadj)     & eqpadj   < 0.05, gene],
  vlm   = merged[!is.na(vlmpadj)    & vlmpadj  < 0.05, gene],
  mfr   = merged[!is.na(mfrpadj)    & mfrpadj  < 0.05, gene]
)
if (deseq_present) deg_sets$deseq <- merged[!is.na(dspadj) & dspadj < 0.05, gene]
if (trend_present) deg_sets$trend <- merged[!is.na(trpadj) & trpadj < 0.05, gene]
if (dlrt_present)  deg_sets$dlrt  <- merged[!is.na(dlpadj) & dlpadj < 0.05, gene]
if (elrt_present)  deg_sets$elrt  <- merged[!is.na(elpadj) & elpadj < 0.05, gene]
if (sva_present)   deg_sets$sva   <- merged[!is.na(svpadj) & svpadj < 0.05, gene]
K <- length(deg_sets)
jacc <- matrix(NA_real_, K, K, dimnames = list(names(deg_sets), names(deg_sets)))
for (i in seq_len(K)) for (j in seq_len(K)) {
  a <- deg_sets[[i]]; b <- deg_sets[[j]]
  jacc[i,j] <- length(intersect(a,b)) / max(length(union(a,b)), 1)
}
fwrite(as.data.table(jacc, keep.rownames = "method"),
       file.path(RDIR_OUT, "jaccard_padj005.csv"))

# --- Direction concordance % among union of padj<0.05 calls ---
union_sig <- Reduce(union, deg_sets)
ds <- merged[gene %in% union_sig]
dir_conc <- c(
  dream_vs_eq  = mean(sign(ds$dream_logFC) == sign(ds$eqlogFC),  na.rm = TRUE),
  dream_vs_vlm = mean(sign(ds$dream_logFC) == sign(ds$vlmlogFC), na.rm = TRUE),
  dream_vs_mfr = mean(sign(ds$dream_logFC) == sign(ds$mfrlogFC), na.rm = TRUE)
)
if (deseq_present) {
  dir_conc <- c(dir_conc,
    dream_vs_deseq = mean(sign(ds$dream_logFC) == sign(ds$dslogFC), na.rm = TRUE))
}
fwrite(as.data.table(dir_conc, keep.rownames = "pair"),
       file.path(RDIR_OUT, "direction_concordance.csv"))

# --- Per-gene robustness flag ---
# Voting arms: edgeR-QL + voomLmFit + DESeq2 (NB/voom-family). metafor-HKSJ is
# a heterogeneity diagnostic (τ²/I²), NOT a DEG voter — its small-K HKSJ p
# is systematically too conservative to participate in voting (0 DEGs at
# K=5). Total voters n_voters = 3 if DESeq2 present, else 2.
n_voters <- 2L + as.integer(deseq_present)
sens_sig <- with(merged,
                 as.integer(!is.na(eqpadj)  & eqpadj  < 0.05) +
                 as.integer(!is.na(vlmpadj) & vlmpadj < 0.05))
if (deseq_present) {
  sens_sig <- sens_sig +
              with(merged, as.integer(!is.na(dspadj) & dspadj < 0.05))
}
dream_sig <- !is.na(merged$dream_padj) & merged$dream_padj < 0.05

flag <- ifelse(dream_sig & sens_sig == n_voters,           "confirmed_all",
        ifelse(dream_sig & sens_sig >= n_voters - 1L,      "confirmed_majority",
        ifelse(dream_sig & sens_sig == 1L,                 "confirmed_minority",
        ifelse(dream_sig & sens_sig == 0L,                 "primary_only",
        ifelse(!dream_sig & sens_sig >= ceiling(n_voters/2),
                                                            "discordant",
                                                            "ns")))))
merged[, dream_robustness_flag := factor(flag,
       levels = c("confirmed_all","confirmed_majority","confirmed_minority",
                  "primary_only","discordant","ns"))]
cat("Robustness voters: edgeR-QL + voomLmFit",
    if (deseq_present) "+ DESeq2" else "",
    " (n =", n_voters, "); metafor excluded (heterogeneity-only)\n")

# --- Bring in mashr sharing classes ---
mash_path <- file.path(RDIR_IN, "mashr/sharing_classes.csv")
if (file.exists(mash_path)) {
  mash_sh <- fread(mash_path)
  merged <- merge(merged,
                  mash_sh[, .(gene, n_sig_cohorts, sharing_class, pan_cohort)],
                  by = "gene", all.x = TRUE)
  setnames(merged, c("n_sig_cohorts","sharing_class","pan_cohort"),
                   c("mashr_n_sig_cohorts","mashr_sharing_class","mashr_pan_cohort"))
} else {
  merged[, `:=`(mashr_n_sig_cohorts = NA_integer_,
                mashr_sharing_class = NA_character_,
                mashr_pan_cohort    = NA)]
  mash_sh <- data.table()
}

# --- Build the atlas_columns.tsv (13 cols, or 15 with DESeq2) ---
atlas_cols <- merged[, .(
  gene,
  sensitivity_edgeRql_logFC   = eqlogFC,
  sensitivity_edgeRql_padj    = eqpadj,
  sensitivity_voomLmFit_logFC = vlmlogFC,
  sensitivity_voomLmFit_padj  = vlmpadj,
  metafor_logFC               = mfrlogFC,
  metafor_SE                  = meta_SE,
  metafor_tau2                = tau2,
  metafor_I2                  = I2,
  metafor_HKSJ_padj           = mfrpadj,
  mashr_n_sig_cohorts,
  mashr_sharing_class,
  mashr_pan_cohort,
  dream_robustness_flag
)]
if (deseq_present) {
  atlas_cols[, sensitivity_DESeq2_logFC := merged$dslogFC]
  atlas_cols[, sensitivity_DESeq2_padj  := merged$dspadj]
}
fwrite(atlas_cols, file.path(RDIR_OUT, "atlas_columns.tsv"),
       sep = "\t", na = "NA")

# --- Robustness summary for manuscript ---
tier1 <- merged[dream_sig & abs(dream_logFC) > 0.5]
robust_summary <- data.table(
  flag = c("confirmed_all","confirmed_majority","confirmed_minority",
           "primary_only","discordant"),
  n = sapply(c("confirmed_all","confirmed_majority","confirmed_minority",
               "primary_only","discordant"),
             function(f) tier1[dream_robustness_flag == f, .N]),
  pct_tier1 = NA_real_
)
robust_summary[, pct_tier1 := 100 * n / max(nrow(tier1), 1)]
fwrite(robust_summary, file.path(RDIR_OUT, "robustness_summary.csv"))

# --- Manuscript placeholder values ---
vlm_rho_dt <- if (file.exists(file.path(RDIR_IN, "voomLmFit/intra_rho.csv"))) {
  fread(file.path(RDIR_IN, "voomLmFit/intra_rho.csv"))
} else { data.table(intra_rho_round2 = NA_real_) }
# arm 2 emits intra_rho_round2 (final consensus) since 2026-05-19 rewrite to
# manual dupCor idiom; older bundles emitted intra_rho only.
intra_rho_value <- if ("intra_rho_round2" %in% names(vlm_rho_dt)) {
  vlm_rho_dt$intra_rho_round2
} else {
  vlm_rho_dt$intra_rho
}

summary_rows <- list(
  c("n_genes_dream_anchor",                 nrow(merged)),
  c("rho_dream_eq",                         round(rho_t["dream_t","eqt"],   3)),
  c("rho_dream_vlm",                        round(rho_t["dream_t","vlmt"],  3)),
  c("rho_dream_mfr",                        round(rho_t["dream_t","mfrt"],  3)),
  c("jaccard_dream_eq",                     round(jacc["dream","eq"],   3)),
  c("jaccard_dream_vlm",                    round(jacc["dream","vlm"],  3)),
  c("jaccard_dream_mfr",                    round(jacc["dream","mfr"],  3)),
  c("direction_conc_dream_eq",              round(dir_conc["dream_vs_eq"],  3)),
  c("direction_conc_dream_vlm",             round(dir_conc["dream_vs_vlm"], 3)),
  c("direction_conc_dream_mfr",             round(dir_conc["dream_vs_mfr"], 3))
)
if (deseq_present) {
  summary_rows <- c(summary_rows, list(
    c("rho_dream_deseq",            round(rho_t["dream_t","dst"], 3)),
    c("jaccard_dream_deseq",        round(jacc["dream","deseq"],   3)),
    c("direction_conc_dream_deseq", round(dir_conc["dream_vs_deseq"], 3)),
    c("n_deseq_padj005",            length(deg_sets$deseq))
  ))
}
summary_rows <- c(summary_rows, list(
  c("n_tier1_dream",                      nrow(tier1)),
  c("pct_confirmed_majority_or_better",
    round(sum(tier1$dream_robustness_flag %in%
              c("confirmed_all","confirmed_majority")) /
          max(nrow(tier1), 1) * 100, 1)),
  c("n_voters",                           n_voters),
  c("intra_cohort_rho_vlm",               round(intra_rho_value, 3))
))
summary_tsv <- data.table(
  metric = sapply(summary_rows, `[`, 1),
  value  = sapply(summary_rows, `[`, 2)
)
fwrite(summary_tsv, file.path(RDIR_OUT, "summary.tsv"), sep = "\t")

# --- Panel data for figure renderer ---
saveRDS(list(
  merged           = merged,
  rho_t            = rho_t,
  r_lfc            = r_lfc,
  jaccard          = jacc,
  direction_conc   = dir_conc,
  deg_sets         = deg_sets,
  sharing_classes  = mash_sh,
  robust_summary   = robust_summary,
  tier1            = tier1
), file.path(RDIR_OUT, "panel_data_A_to_F.rds"))

cat("\nSaved aggregated outputs to ", RDIR_OUT, "\n")
print(summary_tsv)
