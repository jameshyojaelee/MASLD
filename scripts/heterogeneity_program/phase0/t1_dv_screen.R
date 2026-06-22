#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# T1 — DV-MASLD-CN genome-wide differential-variability screen (disease vs control)
#
# Primary axis = disease-vs-control differential variability on the control-bearing
# mega cohorts (best-powered). limma Brown-Forsythe DV statistic (validated:
# test_dv_spikein.R), edgeR NB confirmer, LOCO PASS/FAIL keep rule, genome-wide
# cohort-blocked permutation null (perm_null engine; refit voom+DV per replicate),
# and the headline 2×2 variance-only vs mean-and-variance partition vs Tier-1 DEGs.
#
# Modes (env N_PERM): 0 = observed+LOCO+confirmer only (smoke, fast);
#                     >0 = + permutation null (heavy → sbatch, N_PERM=1000).
# Substrate: merged_dge.rds (RAW counts; never corrected_logcpm — GATE-N8).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/heterogeneity_program/lib/dv_stat.R"))
source(file.path(BASE, "scripts/heterogeneity_program/lib/perm_null.R"))
N_PERM  <- as.integer(Sys.getenv("N_PERM", "0"))
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
OUT <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

# ── Load substrate ────────────────────────────────────────────────────────────
dge  <- readRDS(file.path(INT, "merged_dge.rds"))
cat(sprintf("merged_dge: %d genes × %d samples\n", nrow(dge), ncol(dge)))
cat("samples columns:", paste(colnames(dge$samples), collapse=", "), "\n")
gb  <- as.character(dge$samples$group_binary)         # "Control" / "Disease"
ds  <- as.character(dge$samples$dataset)
# disease-vs-control DV needs cohorts with BOTH groups (≥3 each) so the
# block×group median is estimable in every cell.
ct <- table(ds, gb)
cb_cohorts <- rownames(ct)[ct[, "Control"] >= 3 & ct[, "Disease"] >= 3]
keep_s <- ds %in% cb_cohorts & gb %in% c("Control", "Disease")
cat(sprintf("control-bearing cohorts (both groups, ≥3 each): %s\n", paste(cb_cohorts, collapse=", ")))
cat(sprintf("mega (disease-vs-control) subset: %d samples across %d cohorts\n",
            sum(keep_s), length(cb_cohorts)))
print(table(dataset = ds[keep_s], group = gb[keep_s]))

dge <- dge[, keep_s, keep.lib.sizes = FALSE]
group <- factor(as.character(dge$samples$group_binary), levels = c("Control", "Disease"))
block <- droplevels(factor(as.character(dge$samples$dataset)))

# ── Observed DV (fixed gene universe) ─────────────────────────────────────────
keep_g <- filterByExpr(dge, group = group)
gene_set <- rownames(dge)[keep_g]
cat(sprintf("filterByExpr universe: %d genes\n", length(gene_set)))

# ── Strip confound genes from the universe (repo convention; gencode_v49) ──────
# chrY/chrM/XIST + ribosomal (RP[SL]/MRP[SL]) + IG/TR + hemoglobin carry
# composition/library-prep VARIANCE (sex, technical) — not disease heterogeneity.
# Same strip the NMF pipeline applies (CLAUDE.md: chrY/chrM/RPS-RPL/IG/HB/XIST).
gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
gms <- gm[match(sub("[.].*$", "", gene_set), ensembl_base)]
hb  <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
strip <- (gms$chromosome %in% c("chrY","chrM")) | (gms$gene_name == "XIST") |
         (grepl("^(RP[SL]|MRP[SL])", gms$gene_name) & !grepl("K[ABCL]?[0-9]", gms$gene_name)) |
         grepl("^IG_|^TR_", gms$gene_biotype) | (gms$gene_name %in% hb)
strip[is.na(strip)] <- FALSE
cat(sprintf("stripping %d confound genes (chrY/chrM/XIST/ribosomal/IG/HB) → universe %d → %d\n",
            sum(strip), length(gene_set), sum(!strip)))
gene_set <- gene_set[!strip]
obs <- as.data.table(dv_screen(dge, group, block, gene_set = gene_set))
cat(sprintf("observed DV computed for %d genes (gene id eg: %s)\n",
            nrow(obs), paste(head(obs$gene, 2), collapse=", ")))

# ── NB confirmer (edgeR per-group dispersion log-ratio) ───────────────────────
nb <- as.data.table(dv_nb_confirmer(dge, group, block, gene_set))
obs <- merge(obs, nb, by = "gene", all.x = TRUE)
obs[, nb_concordant := sign(nb_logratio) == sign(dv_t)]

# ── Join canonical bulk DEG (mean-shift axis) — strip ENSG version ────────────
deg <- fread(file.path(INT, "canonical_deg_results.csv"))
deg[, gene_base := sub("\\..*$", "", gene)]
obs[, gene_base := sub("\\..*$", "", gene)]
obs <- merge(obs, deg[, .(gene_base, symbol, bulk_logFC = logFC, bulk_tstat = t,
                          bulk_shrunk_logFC = shrunk_logFC, bulk_lfsr = lfsr)],
             by = "gene_base", all.x = TRUE)
obs[, is_mean_deg := !is.na(bulk_lfsr) & bulk_lfsr < 0.05 & abs(bulk_shrunk_logFC) > 0.5]
low_lfc <- obs[!is.na(bulk_shrunk_logFC) & abs(bulk_shrunk_logFC) < 0.5]
cat(sprintf("\nDecoupling: corr(dv_t, bulk_tstat) ALL = %+.3f | among low-|logFC| genes = %+.3f (target ≈0)\n",
            obs[!is.na(bulk_tstat), cor(dv_t, bulk_tstat)],
            low_lfc[!is.na(bulk_tstat), cor(dv_t, bulk_tstat)]))

# ── LOCO PASS/FAIL keep rule (≥4/5 cohorts, sign-stable) ──────────────────────
cohorts <- levels(block)
loco <- matrix(NA_real_, nrow(obs), length(cohorts), dimnames = list(obs$gene, cohorts))
for (co in cohorts) {
  ti <- which(block != co)
  r <- as.data.table(dv_screen(dge, group, block, train_idx = ti, gene_set = gene_set))
  loco[r$gene, co] <- r$dv_t
}
obs[, loco_signfrac := rowMeans(sign(loco[gene, , drop=FALSE]) == sign(dv_t), na.rm = TRUE)]
obs[, loco_pass := loco_signfrac >= 0.8]

# ── Permutation null (heavy; N_PERM>0) ────────────────────────────────────────
if (N_PERM > 0) {
  cat(sprintf("\nPermutation null: N_PERM=%d, cores=%d (refit voom+DV per replicate)\n", N_PERM, N_CORES))
  stat_fn <- function(g) {
    g <- factor(as.character(g), levels = c("Control","Disease"))
    r <- dv_screen(dge, g, block, gene_set = gene_set)
    setNames(r$dv_t, r$gene)[obs$gene]
  }
  t0 <- Sys.time()
  pn <- perm_null(stat_fn, group = group, strata = block, n_perm = N_PERM,
                  alternative = "two.sided", n_cores = N_CORES, return_null = TRUE)
  cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
  obs[, dv_p_perm := pn$p_perm]; obs[, dv_q_perm := pn$q_perm]   # per-gene (coarse, 1/1001 floor)

  # ── Decile-matched POOLED null (fine resolution; 46d idiom) ─────────────────
  # 1000 perms floor per-gene p at 1e-3 → BH can't reach q<0.05 even for genuine
  # hits (0/1000, LOCO-stable). Pool permuted |stat| within AveExpr deciles for
  # ~1/(2.6M) resolution under within-decile exchangeability of the moderated stat.
  dec <- cut(obs$AveExpr, quantile(obs$AveExpr, seq(0, 1, .1), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)
  # Z-STANDARDISED pooling (review C1 fix): each gene's permutation null has its
  # own spread (eBayes moderation + per-gene df → heteroskedastic within a decile),
  # so pooling raw |null| is anti-conservative for tight-null genes. Standardise
  # each gene's obs+null by its OWN null SD → unit-variance nulls → valid to pool.
  null_sd <- apply(pn$null, 1L, sd, na.rm = TRUE)
  null_sd[!is.finite(null_sd) | null_sd == 0] <- NA_real_
  z_obs  <- abs(pn$observed) / null_sd
  z_null <- abs(pn$null) / null_sd
  p_pool <- rep(NA_real_, nrow(obs))
  for (dd in sort(unique(dec))) {
    idx  <- which(dec == dd & is.finite(z_obs))
    if (length(idx) == 0L) next
    pool <- sort(as.numeric(z_null[idx, , drop = FALSE])); pool <- pool[is.finite(pool)]
    np   <- length(pool)
    cnt_ge <- np - findInterval(z_obs[idx] - 1e-9, pool)         # # pooled standardized |null| ≥ |z_obs|
    p_pool[idx] <- (1 + cnt_ge) / (1 + np)
  }
  obs[, dv_p_pool := p_pool]; obs[, dv_q_pool := p.adjust(p_pool, "BH")]
  cat(sprintf("  z-pooled-null calibration: median p_pool=%.3f (≈0.5 if null; per-gene perm median=%.3f same ⇒ real enrichment); cor(p_perm,p_pool) mid=%.3f\n",
              median(p_pool, na.rm = TRUE), median(obs$dv_p_perm, na.rm = TRUE),
              obs[dv_p_perm > 0.1 & dv_p_perm < 0.9, cor(dv_p_perm, dv_p_pool, use = "complete.obs")]))
  # DV-sig = z-pooled FDR + LOCO-stable (rigorous gates). NB-concordance is a
  # confidence FLAG (review M3/M4: NB tagwise dispersion is noisy at small control
  # cells — must not hard-veto genuine hits; doc says "concordance, not significance").
  obs[, dv_sig := !is.na(dv_q_pool) & dv_q_pool < 0.05 & loco_pass]
  obs[, dv_confidence := fcase(dv_sig & nb_concordant, "high", dv_sig & !nb_concordant, "moderate", default = "ns")]
  cat(sprintf("  DV-sig (z-pool q<0.05 & LOCO) = %d  [of which NB-concordant 'high' = %d]  [per-gene q_perm<0.05 = %d]\n",
              sum(obs$dv_sig), sum(obs$dv_confidence == "high"), sum(obs$dv_q_perm < 0.05, na.rm = TRUE)))
} else {
  obs[, dv_q_eb := p.adjust(dv_p, "BH")]
  obs[, dv_sig := dv_q_eb < 0.05 & loco_pass]   # provisional (eBayes) until perms
}

# ── Headline 2×2: variance-only vs mean-and-variance ──────────────────────────
obs[, partition := fcase(
  dv_sig & !is_mean_deg & (is.na(bulk_shrunk_logFC) | abs(bulk_shrunk_logFC) < 0.5), "variance_only",
  dv_sig &  is_mean_deg, "mean_and_variance",
  !dv_sig & is_mean_deg, "mean_only",
  default = "neither")]
cat("\n── 2×2 partition (DV-sig × mean-DEG) ──\n"); print(obs[, .N, by = partition][order(-N)])
cat(sprintf("Tier-1 DEGs (mean): %d | DV-sig: %d | variance-ONLY (novel): %d\n",
            obs[is_mean_deg == TRUE, .N], obs[dv_sig == TRUE, .N], obs[partition=="variance_only", .N]))

fwrite(obs, file.path(OUT, sprintf("dv_results%s.csv", if (N_PERM>0) "" else "_smoke")))
cat(sprintf("\nWrote dv_results%s.csv\n", if (N_PERM>0) "" else "_smoke"))
