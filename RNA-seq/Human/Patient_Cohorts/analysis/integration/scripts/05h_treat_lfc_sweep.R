#!/usr/bin/env Rscript
# ===========================================================================
# 05h_treat_lfc_sweep.R  (audit P2#18 / R-T1, 2026-06-29)
#
# Sensitivity artifact for the canonical TREAT effect-size floor. Re-fits the
# EXACT C2 limma_voom_qw engine from 05h_limma_voom_qw_canonical.R (same input
# loading: merged_dge.rds subset to include_in_mega cohorts + meta_matched.rds;
# same design ~ dataset + inferred_sex + group_binary), then evaluates the
# TREAT interval-null gate (treat_fdr < 0.05) at FOUR lfc thresholds
#   CANONICAL_TREAT_LFC in {0.15, 0.20, 0.25, 0.30}
# and writes a single combined table of DEG counts per lfc. This shows that the
# canonical 1,918 at lfc=0.25 is not a knife-edge.
#
# IMPORTANT: this script does NOT overwrite canonical_deg_results.csv (unlike
# 05h, which rewrites it on every run). It writes ONLY:
#   RNA-seq/results/audit_sensitivity/treat_lfc_sweep.csv
#
# Env: rnaseq.  SLURM: cpu, 16 CPU, ~64G, 48h.
# ===========================================================================
suppressMessages({library(edgeR); library(limma); library(data.table); library(yaml)})
set.seed(42)

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BASE <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts")
RDIR <- file.path(BASE, "analysis/integration/results/integration")
OUTDIR <- file.path(ROOT, "RNA-seq/results/audit_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- input loading (identical to 05h lines 20-33) ----
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

# ---- C2 design + fit (identical to 05h lines 36-43) ----
design <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
coef_name <- "group_binaryDisease"
stopifnot(coef_name %in% colnames(design))
v    <- voomWithQualityWeights(dge_mega, design)
fit0 <- lmFit(v, design)
fit  <- eBayes(fit0)
res  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)
dt <- as.data.table(res); setnames(dt, "adj.P.Val", "padj")

# ---- TREAT gate at each lfc in the sweep ----
LFC_GRID <- as.numeric(strsplit(Sys.getenv("TREAT_LFC_GRID", "0.15,0.20,0.25,0.30"), ",")[[1]])
rows <- lapply(LFC_GRID, function(L) {
  tt <- topTreat(treat(fit0, lfc = L), coef = coef_name, number = Inf, sort.by = "none")
  fdr <- tt[dt$gene, "adj.P.Val"]
  lfc <- dt$logFC
  sig <- !is.na(fdr) & fdr < 0.05
  data.table(
    treat_lfc      = L,
    genes_tested   = nrow(dt),
    n_deg          = sum(sig),
    n_up           = sum(sig & lfc > 0),
    n_down         = sum(sig & lfc < 0),
    min_abs_logFC  = if (any(sig)) round(min(abs(lfc[sig])), 4) else NA_real_)
})
sweep <- rbindlist(rows)
cat("\n=== TREAT lfc sweep ===\n"); print(sweep, row.names = FALSE)

fwrite(sweep, file.path(OUTDIR, "treat_lfc_sweep.csv"))
cat(sprintf("\n[done] wrote %s (canonical_deg_results.csv NOT touched)\n",
            file.path(OUTDIR, "treat_lfc_sweep.csv")))
