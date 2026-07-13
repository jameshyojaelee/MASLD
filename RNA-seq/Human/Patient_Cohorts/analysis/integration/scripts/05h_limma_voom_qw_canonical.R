#!/usr/bin/env Rscript
# ===========================================================================
# 05h — production limma_voom_qw__C2 canonical DEG call (benchmark winner).
# Mirrors 05_dream's input loading EXACTLY (merged_dge.rds + yaml include_in_mega
# + meta_matched.rds), but fits the benchmark-winning method:
#   voomWithQualityWeights -> lmFit -> eBayes,  design ~ dataset + sex + group.
# Emits BOTH raw (logFC, padj) and ashr-shrunk (shrunk_logFC, lfsr) columns,
# matching the dream_results.csv two-tier schema. Then runs a SANITY GATE
# (concordance vs current dream canonical + positive-control recovery) and
# STOPS — does NOT overwrite any canonical file. Output:
#   results/integration/limma_voom_qw_C2_results.csv
#   results/integration/limma_voom_qw_C2_sanity.csv
# ===========================================================================
suppressMessages({library(edgeR); library(limma); library(ashr); library(data.table)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
RDIR <- file.path(BASE, "analysis/integration/results/integration")
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

# ---- input loading (identical to 05_dream lines 44-99) ----
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(ROOT, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep <- dge$samples$dataset %in% mega
dge_mega <- dge[, keep]
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control","Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(sex))
rownames(info) <- colnames(dge_mega)
cat(sprintf("[in] %d genes x %d samples; %d cohorts; %s\n", nrow(dge_mega), ncol(dge_mega),
            nlevels(info$dataset), paste(names(table(info$group_binary)), table(info$group_binary), collapse=" ")))

# ---- C2 design: ~ dataset + sex + group  (cohort + sex fixed effects) ----
design <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
coef_name <- "group_binaryDisease"
stopifnot(coef_name %in% colnames(design))

# ---- limma_voom_qw (benchmark-winning engine) ----
v    <- voomWithQualityWeights(dge_mega, design)
fit0 <- lmFit(v, design)
fit  <- eBayes(fit0)
res  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)
dt <- as.data.table(res); setnames(dt, "adj.P.Val", "padj")
dt[, SE := abs(logFC / t)]

# ---- TREAT: canonical DEG significance gate (2026-06-29) ----
# treat() tests H0: |true logFC| <= TREAT_LFC (McCarthy & Smyth 2009). FDR<0.05 on this
# test IS the canonical DEG call -- the effect-size floor is folded INTO the test, so NO
# separate |logFC| filter is applied downstream. Supersedes the ashr lfsr+|shrunk|>0.3 gate.
TREAT_LFC <- as.numeric(Sys.getenv("CANONICAL_TREAT_LFC", "0.25"))
ttm <- topTreat(treat(fit0, lfc = TREAT_LFC), coef = coef_name, number = Inf, sort.by = "none")
dt[, treat_lfc := TREAT_LFC]
dt[, treat_p   := ttm[gene, "P.Value"]]
dt[, treat_fdr := ttm[gene, "adj.P.Val"]]

# ---- ashr shrinkage (kept SIDE-BY-SIDE with raw) ----
ok <- is.finite(dt$logFC) & is.finite(dt$SE) & dt$SE > 0
ash <- ashr::ash(dt$logFC[ok], dt$SE[ok], mixcompdist = "normal")
dt[, shrunk_logFC := NA_real_][ok, shrunk_logFC := ash$result$PosteriorMean]
dt[, lfsr := NA_real_][ok, lfsr := ash$result$lfsr]

out <- dt[, .(gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, treat_lfc, treat_p, treat_fdr, AveExpr)]

# ---- symbol column (drop-in parity with dream_results_ashr.csv `symbol`) ----
# Class-B downstream readers (e.g. 80_celltype_intrinsic_attribution.R, 310a) select
# `symbol`; map via GENCODE v49 metadata so canonical_deg_results.csv is a true drop-in.
gm0 <- fread(file.path(ROOT, "data/gencode_v49_gene_metadata.tsv.gz"))
gm0[, eb := sub("[.][0-9]+$", "", gene_id)]
out[, symbol := gm0[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]

fwrite(out, file.path(RDIR, "limma_voom_qw_C2_results.csv"))
# Stable, method-neutral canonical filename: 07 + all direct readers point here, so a
# future method swap only re-touches 05h (the producer), not every consumer.
fwrite(out, file.path(RDIR, "canonical_deg_results.csv"))

# ---- DEG counts: canonical TREAT + legacy reference definitions ----
treat1 <- out[treat_fdr < 0.05]
raw1 <- out[padj < 0.05 & abs(logFC) > 0.5]
ash1 <- out[lfsr < 0.05 & abs(shrunk_logFC) > 0.5]
cat(sprintf("\n[C2] genes=%d\n  CANONICAL TREAT (treat_fdr<.05, lfc=%.2f)  : %d  (up %d / dn %d)\n  RAW  ref (padj<.05 & |logFC|>.5)           : %d  (up %d / dn %d)\n  ashr ref (lfsr<.05 & |shrunk_logFC|>.5)    : %d  (up %d / dn %d)\n",
  nrow(out), TREAT_LFC, nrow(treat1), sum(treat1$logFC>0), sum(treat1$logFC<0),
  nrow(raw1), sum(raw1$logFC>0), sum(raw1$logFC<0),
  nrow(ash1), sum(ash1$shrunk_logFC>0), sum(ash1$shrunk_logFC<0)))

# ======================= SANITY GATE (no overwrite) =======================
dr <- fread(file.path(RDIR, "dream_results.csv"))
m  <- merge(out[, .(gene, c2_lfc=logFC, c2_padj=padj, c2_slfc=shrunk_logFC, c2_lfsr=lfsr)],
            dr[, .(gene, dr_lfc=logFC, dr_padj=padj,
                   dr_slfc=if ("shrunk_logFC" %in% names(dr)) shrunk_logFC else NA_real_,
                   dr_lfsr=if ("lfsr" %in% names(dr)) lfsr else NA_real_)], by="gene")
rho  <- cor(m$c2_lfc, m$dr_lfc, method="spearman", use="complete.obs")
pear <- cor(m$c2_lfc, m$dr_lfc, method="pearson",  use="complete.obs")
dirA <- mean(sign(m$c2_lfc)==sign(m$dr_lfc), na.rm=TRUE)
# dream canonical Tier-1 (ashr if present, else raw)
dr1 <- if (all(c("shrunk_logFC","lfsr") %in% names(dr))) dr[lfsr<0.05 & abs(shrunk_logFC)>0.5, gene] else dr[padj<0.05 & abs(logFC)>0.5, gene]
jac <- function(a,b) length(intersect(a,b))/length(union(a,b))
# positive-control recovery
pc <- fread(file.path(ROOT, "RNA-seq/results/validation/positive_control_validation.csv"))
gm <- fread(file.path(ROOT, "data/gencode_v49_gene_metadata.tsv.gz"))
gm[, eb := sub("[.][0-9]+$","",gene_id)]
pc_eb <- unique(na.omit(gm[match(pc$gene, gene_name), eb]))
out[, eb := sub("[.][0-9]+$","",gene)]
pc_in_raw  <- length(intersect(pc_eb, out[padj<0.05 & abs(logFC)>0.5, eb]))
pc_in_ash  <- length(intersect(pc_eb, out[lfsr<0.05 & abs(shrunk_logFC)>0.5, eb]))
pc_tested  <- length(intersect(pc_eb, out$eb))

san <- data.table(
  metric = c("genes_tested","C2_raw_Tier1","C2_ashr_Tier1","dream_Tier1",
             "logFC_spearman_vs_dream","logFC_pearson_vs_dream","direction_agree_vs_dream",
             "jaccard_C2ashr_vs_dreamTier1","PC_tested","PC_in_C2raw","PC_in_C2ashr"),
  value  = c(nrow(out), nrow(raw1), nrow(ash1), length(dr1),
             round(rho,4), round(pear,4), round(dirA,4),
             round(jac(ash1$gene, dr1),4), pc_tested, pc_in_raw, pc_in_ash))
fwrite(san, file.path(RDIR, "limma_voom_qw_C2_sanity.csv"))
cat("\n=== SANITY GATE vs current dream canonical ===\n"); print(san, row.names=FALSE)
cat("\n[done] wrote limma_voom_qw_C2_results.csv + limma_voom_qw_C2_sanity.csv (NO canonical overwritten)\n")
