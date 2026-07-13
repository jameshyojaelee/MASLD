#!/usr/bin/env Rscript
# 07_consensus_degs.R
# ---------------------------------------------------------------------------
# Identify primary disease-associated DEGs from the CANONICAL pooled
# (cohort-adjusted) analysis: limma-voom quality-weighted, C2 design
# (~ dataset + inferred_sex + group_binary; benchmark winner 2026-06-07,
#  replaces dream). Reads the stable canonical table produced by
# 05h_limma_voom_qw_canonical.R.
#
# Canonical threshold (ashr): lfsr < 0.05 AND |shrunk_logFC| > 0.5.
# Emits method-neutral `bulk_*` columns consumed by 27a and downstream
# (hard cutover from the legacy `dream_*` names, 2026-06-07).
#
# Output: results/integration/consensus_degs.csv
#         results/integration/consensus_significant_degs.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load the canonical DEG table (limma_voom_qw C2) ---
canon_file <- file.path(RDIR, "canonical_deg_results.csv")
if (!file.exists(canon_file)) {
  stop("canonical_deg_results.csv not found — run 05h_limma_voom_qw_canonical.R first")
}
deg <- fread(canon_file)
stopifnot(all(c("shrunk_logFC", "lfsr", "treat_fdr") %in% names(deg)))
cat("Loaded canonical (limma_voom_qw C2) results:", nrow(deg), "genes\n")

# --- Significance: TREAT FDR (canonical 2026-06-29) ---
# treat_fdr = real limma::treat() FDR at lfc=0.25 (computed in 05h). FDR<0.05 IS the DEG
# call: treat() tests H0:|true logFC|<=lfc, so the effect-size floor is folded INTO the
# test and NO separate |logFC| filter is applied. Supersedes the ashr lfsr+|shrunk|>0.3 gate
# (which is retained in the table as a reference column).
TREAT_FDR_THRESH <- 0.05
deg[, bulk_sig := treat_fdr < TREAT_FDR_THRESH]
deg[, bulk_dir := sign(logFC)]   # TREAT tests raw logFC; direction = logFC sign

cat("\n===== CANONICAL DE RESULTS (TREAT FDR, lfc =", deg$treat_lfc[1], ") =====\n")
cat("Threshold: treat_fdr <", TREAT_FDR_THRESH, "\n")
cat("Total genes:", nrow(deg), "\n\n")

sig_genes <- deg[bulk_sig == TRUE]
cat("Significant DEGs:", nrow(sig_genes), "\n")
cat("  Up:", sum(sig_genes$bulk_dir == 1, na.rm = TRUE), "\n")
cat("  Down:", sum(sig_genes$bulk_dir == -1, na.rm = TRUE), "\n")

# --- Rename to method-neutral bulk_* schema (hard cutover from dream_*) ---
# bulk_logFC        = raw logFC (reference; thresholding uses shrunk effect)
# bulk_shrunk_logFC = ashr shrunk effect size (canonical)
# bulk_padj         = BH-adjusted p (reference; thresholding uses lfsr)
# bulk_lfsr         = ashr local false sign rate (canonical significance)
# se                = standard error (from SE); t / P.Value / AveExpr / symbol kept.
# Dropped vs legacy dream schema: z.std, shrunk_se, svalue (unused past consensus).
setnames(deg, "logFC", "bulk_logFC")
setnames(deg, "shrunk_logFC", "bulk_shrunk_logFC")
setnames(deg, "padj", "bulk_padj")
setnames(deg, "lfsr", "bulk_lfsr")
if ("treat_fdr" %in% names(deg)) setnames(deg, "treat_fdr", "bulk_treat_fdr")
if ("treat_p" %in% names(deg))   setnames(deg, "treat_p", "bulk_treat_p")
if ("SE" %in% names(deg)) setnames(deg, "SE", "se")

# --- Save ---
deg <- deg[order(bulk_padj)]
fwrite(deg, file.path(RDIR, "consensus_degs.csv"))
cat("\nSaved: consensus_degs.csv\n")

fwrite(deg[bulk_sig == TRUE], file.path(RDIR, "consensus_significant_degs.csv"))
cat("Saved: consensus_significant_degs.csv (", nrow(sig_genes), "genes)\n")

cat("\n===== COMPLETE =====\n")
