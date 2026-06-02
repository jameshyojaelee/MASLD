#!/usr/bin/env python3
"""
160b_fix_hcc_score.py
Fix HCC molecular score gene mapping bug.

Problem: Script 150 computed HCC score from 3K most-variable genes in
prepared_data.h5. Only 15/59 COLOC genes mapped. High-PP.H4 genes like
ATF4 (0.952), HKDC1 (0.896), STOX2 (0.787), KIF7 (0.795) were missed.

Fix: Load full-genome expression (rank_expression_full.rds, ~34K genes),
map all 59 HCC COLOC genes via Ensembl annotation, z-score per-gene
across ALL samples, and recompute the weighted HCC molecular score.

Input:
  - results/staging_classifier/rank_expression_full.rds   (34K genes x 1,444 samples)
  - results/prognosis/hcc_coloc_genes.csv                 (59 genes with PP.H4 weights)
  - results/gene_annotation/human_ensg_to_symbol.tsv      (Ensembl -> symbol)
  - results/staging_classifier/modeling_metadata.csv       (sample metadata)
  - results/prognosis/prognosis_pseudo_labels.csv          (P(F4), CPS, s2_binary)
  - RNA-seq/results/subtypes/nmf_assignments.csv           (NMF S1/S2 labels)

Output (results/prognosis_v2/):
  - hcc_molecular_score_v2.csv      (per-sample scores with full gene set)
  - hcc_validation_v2.csv           (6 validation tests)
  - hcc_gene_mapping_report.csv     (mapping details for all 59 genes)

Env:   micromamba activate spatial
SLURM: cpu partition, 4 CPUs, 16G, 48h
"""

import os
import sys
import time
import subprocess
import tempfile
import logging
import warnings

import numpy as np
import pandas as pd
from scipy import stats
import h5py

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Paths ───────────────────────────────────────────────────────────────────
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INTEGRATION = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")

RANK_FULL_RDS = os.path.join(INTEGRATION, "results/staging_classifier/rank_expression_full.rds")
HCC_COLOC_CSV = os.path.join(INTEGRATION, "results/prognosis/hcc_coloc_genes.csv")
GENE_ANNOT    = os.path.join(INTEGRATION, "results/gene_annotation/human_ensg_to_symbol.tsv")
MODEL_META    = os.path.join(INTEGRATION, "results/staging_classifier/modeling_metadata.csv")
PSEUDO_LABELS = os.path.join(INTEGRATION, "results/prognosis/prognosis_pseudo_labels.csv")
NMF_ASSIGN    = os.path.join(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")

OUT_DIR = os.path.join(INTEGRATION, "results/prognosis_v2")
os.makedirs(OUT_DIR, exist_ok=True)


# ── Step 0: Extract full-genome rank expression via R subprocess ────────────
def extract_rank_expression_to_h5(rds_path, out_h5):
    """Use R to load rank_expression_full.rds and write to HDF5 for Python."""
    log.info("Extracting full-genome rank expression from RDS via R subprocess...")

    r_script = f"""
suppressPackageStartupMessages({{
  library(rhdf5)
}})

mat <- readRDS("{rds_path}")
cat(sprintf("Loaded rank expression: %d genes x %d samples\\n", nrow(mat), ncol(mat)))

gene_ids <- rownames(mat)
sample_ids <- colnames(mat)

out_file <- "{out_h5}"
if (file.exists(out_file)) file.remove(out_file)
h5createFile(out_file)

h5createDataset(out_file, "expression",
                dims = dim(mat),
                chunk = c(min(1000, nrow(mat)), ncol(mat)),
                level = 4)
h5write(mat, out_file, "expression")
h5write(gene_ids, out_file, "gene_ids")
h5write(sample_ids, out_file, "sample_ids")

cat(sprintf("Saved %d genes x %d samples to %s\\n", nrow(mat), ncol(mat), out_file))
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
        f.write(r_script)
        r_tmp = f.name

    try:
        result = subprocess.run(
            ["micromamba", "run", "-n", "rnaseq", "Rscript", r_tmp],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            log.error("R extraction failed:\n%s", result.stderr)
            sys.exit(1)
        log.info("R output:\n%s", result.stdout.strip())
    finally:
        os.unlink(r_tmp)


def load_h5_expression(h5_path):
    """Load expression matrix from HDF5 written by R (rhdf5 = column-major)."""
    with h5py.File(h5_path, "r") as f:
        raw = f["expression"][:]  # R col-major -> appears transposed in h5py
        gene_ids = [g.decode() if isinstance(g, bytes) else g for g in f["gene_ids"][:]]
        sample_ids = [s.decode() if isinstance(s, bytes) else s for s in f["sample_ids"][:]]

    # rhdf5 stores (genes, samples) in column-major; h5py reads as (samples, genes)
    # Transpose back to (genes, samples)
    expr = raw.T
    log.info("Loaded H5: %d genes x %d samples", expr.shape[0], expr.shape[1])

    assert expr.shape[0] == len(gene_ids), (
        f"Gene mismatch: {expr.shape[0]} rows vs {len(gene_ids)} gene_ids"
    )
    assert expr.shape[1] == len(sample_ids), (
        f"Sample mismatch: {expr.shape[1]} cols vs {len(sample_ids)} sample_ids"
    )
    return expr, gene_ids, sample_ids


# ── Main ────────────────────────────────────────────────────────────────────
def main():
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Load HCC COLOC genes
    # ------------------------------------------------------------------
    log.info("=== Step 1: Load HCC COLOC genes ===")
    hcc_genes = pd.read_csv(HCC_COLOC_CSV)
    n_total = len(hcc_genes)
    n_old_mapped = hcc_genes["in_expression_matrix"].sum()
    log.info("HCC COLOC genes: %d total, %d previously mapped (3K subset)", n_total, n_old_mapped)

    # ------------------------------------------------------------------
    # 2. Load gene annotation (symbol <-> Ensembl)
    # ------------------------------------------------------------------
    log.info("=== Step 2: Load gene annotation ===")
    annot = pd.read_csv(GENE_ANNOT, sep="\t")
    log.info("Gene annotation: %d entries", len(annot))

    # Build symbol -> list of Ensembl base IDs (some symbols map to multiple)
    # Prefer protein_coding biotype
    symbol_to_ensg = {}
    for _, row in annot.iterrows():
        sym = row["symbol"]
        ensg = row["gene_base"]  # e.g., ENSG00000128272
        gtype = row.get("gene_type", "")
        if sym not in symbol_to_ensg:
            symbol_to_ensg[sym] = []
        symbol_to_ensg[sym].append((ensg, gtype))

    log.info("Unique symbols in annotation: %d", len(symbol_to_ensg))

    # ------------------------------------------------------------------
    # 3. Load full-genome rank expression
    # ------------------------------------------------------------------
    log.info("=== Step 3: Load full-genome expression ===")
    h5_tmp = os.path.join(OUT_DIR, "_tmp_rank_full.h5")

    if not os.path.exists(h5_tmp):
        extract_rank_expression_to_h5(RANK_FULL_RDS, h5_tmp)
    else:
        log.info("Using cached H5: %s", h5_tmp)

    expr, gene_ids, sample_ids = load_h5_expression(h5_tmp)
    # gene_ids are versioned Ensembl IDs (e.g., ENSG00000128272.17)
    log.info("Expression matrix: %d genes x %d samples", len(gene_ids), len(sample_ids))

    # Build gene_id (versioned) -> index mapping
    gene_id_to_idx = {gid: i for i, gid in enumerate(gene_ids)}
    # Build gene_base -> index mapping (strip version)
    gene_base_to_idx = {}
    for gid, idx in gene_id_to_idx.items():
        base = gid.split(".")[0]
        if base not in gene_base_to_idx:
            gene_base_to_idx[base] = idx

    log.info("Unique base gene IDs in expression: %d", len(gene_base_to_idx))

    # Build gene_base -> symbol mapping from annotation (vectorized, not iterating)
    base_to_symbol = dict(zip(annot["gene_base"], annot["symbol"]))

    # Build symbol -> expression index mapping
    symbol_to_expr_idx = {}
    for gid, idx in gene_id_to_idx.items():
        base = gid.split(".")[0]
        sym = base_to_symbol.get(base)
        if sym and sym not in symbol_to_expr_idx:
            symbol_to_expr_idx[sym] = idx

    log.info("Symbols mapped to expression indices: %d", len(symbol_to_expr_idx))

    # ------------------------------------------------------------------
    # 4. Map HCC COLOC genes to expression matrix
    # ------------------------------------------------------------------
    log.info("=== Step 4: Map HCC genes to full expression ===")
    mapping_rows = []
    mapped_indices = []
    mapped_weights = []
    mapped_symbols = []

    for _, row in hcc_genes.iterrows():
        sym = row["gene"]
        pp_h4 = row["PP.H4"]
        was_mapped = row["in_expression_matrix"]

        # Try direct symbol lookup
        idx = symbol_to_expr_idx.get(sym)

        # If direct failed, try via annotation: symbol -> Ensembl base -> expression
        if idx is None and sym in symbol_to_ensg:
            candidates = symbol_to_ensg[sym]
            # Prefer protein_coding
            pc = [c for c in candidates if c[1] == "protein_coding"]
            search = pc if pc else candidates
            for ensg_base, _ in search:
                if ensg_base in gene_base_to_idx:
                    idx = gene_base_to_idx[ensg_base]
                    break

        now_mapped = idx is not None
        mapping_rows.append({
            "gene": sym,
            "PP.H4": pp_h4,
            "was_mapped_3k": was_mapped,
            "now_mapped_full": now_mapped,
            "newly_recovered": now_mapped and not was_mapped,
            "expression_index": idx if now_mapped else np.nan,
        })

        if now_mapped:
            mapped_indices.append(idx)
            mapped_weights.append(pp_h4)
            mapped_symbols.append(sym)

    mapping_df = pd.DataFrame(mapping_rows)
    n_now_mapped = mapping_df["now_mapped_full"].sum()
    n_recovered = mapping_df["newly_recovered"].sum()
    n_still_missing = n_total - n_now_mapped

    log.info("Gene mapping results:")
    log.info("  Total HCC COLOC genes:     %d", n_total)
    log.info("  Previously mapped (3K):    %d", int(n_old_mapped))
    log.info("  Now mapped (full 34K):     %d", n_now_mapped)
    log.info("  Newly recovered:           %d", n_recovered)
    log.info("  Still unmapped:            %d", n_still_missing)

    # Log the unmapped genes
    unmapped = mapping_df[~mapping_df["now_mapped_full"]]
    if len(unmapped) > 0:
        log.info("Unmapped genes (%d):", len(unmapped))
        for _, row in unmapped.iterrows():
            log.info("  %s (PP.H4=%.3f)", row["gene"], row["PP.H4"])

    # Log the newly recovered genes
    recovered = mapping_df[mapping_df["newly_recovered"]]
    if len(recovered) > 0:
        log.info("Newly recovered genes (%d):", len(recovered))
        for _, row in recovered.sort_values("PP.H4", ascending=False).iterrows():
            log.info("  %s (PP.H4=%.3f)", row["gene"], row["PP.H4"])

    if len(mapped_indices) == 0:
        log.error("No HCC genes mapped. Cannot compute score.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 5. Compute HCC molecular score
    # ------------------------------------------------------------------
    log.info("=== Step 5: Compute HCC score ===")
    weights = np.array(mapped_weights)
    gene_exprs = expr[mapped_indices, :]  # (n_mapped_genes, n_samples)

    # Z-score per gene across ALL samples
    gene_z = np.zeros_like(gene_exprs, dtype=np.float64)
    for j in range(gene_exprs.shape[0]):
        row = gene_exprs[j, :]
        mu = np.nanmean(row)
        sd = np.nanstd(row)
        if sd > 0:
            gene_z[j, :] = (row - mu) / sd
        else:
            gene_z[j, :] = 0.0

    # HCC score = weighted sum across genes for each sample
    hcc_score = gene_z.T @ weights  # (n_samples,)
    log.info("HCC score: mean=%.4f, std=%.4f, n=%d",
             np.mean(hcc_score), np.std(hcc_score), len(hcc_score))

    # ------------------------------------------------------------------
    # 6. Load metadata for validation
    # ------------------------------------------------------------------
    log.info("=== Step 6: Load metadata ===")
    meta = pd.read_csv(MODEL_META)
    log.info("Modeling metadata: %d samples", len(meta))

    pseudo = pd.read_csv(PSEUDO_LABELS)
    log.info("Pseudo labels: %d samples", len(pseudo))

    nmf = pd.read_csv(NMF_ASSIGN)
    log.info("NMF assignments: %d samples", len(nmf))

    # Build per-sample output DataFrame
    # sample_ids from expression matrix
    score_df = pd.DataFrame({
        "sample_id": sample_ids,
        "hcc_molecular_score": hcc_score,
    })

    # Merge metadata
    meta_cols = meta[["sample_id", "dataset", "fibrosis_stage", "nas_score"]].copy()
    score_df = score_df.merge(meta_cols, on="sample_id", how="left")

    # Merge pseudo labels (P(F4), CPS, s2_binary)
    pseudo_cols = pseudo[["sample_id", "s2_binary", "p_f4", "cps"]].copy()
    score_df = score_df.merge(pseudo_cols, on="sample_id", how="left")

    # Merge NMF subtype
    nmf_cols = nmf[["sample_id", "nmf_subtype"]].copy()
    score_df = score_df.merge(nmf_cols, on="sample_id", how="left")

    # ------------------------------------------------------------------
    # 7. Validation tests (reproduce Script 152 tests)
    # ------------------------------------------------------------------
    log.info("=== Step 7: Validation tests ===")
    validation_rows = []

    hcc = score_df["hcc_molecular_score"].values.astype(float)
    fib = pd.to_numeric(score_df["fibrosis_stage"], errors="coerce").values
    pf4 = pd.to_numeric(score_df["p_f4"], errors="coerce").values
    cps = pd.to_numeric(score_df["cps"], errors="coerce").values
    nas = pd.to_numeric(score_df["nas_score"], errors="coerce").values
    subtype = score_df["nmf_subtype"].values

    # Test 1: HCC score vs fibrosis stage (Spearman)
    valid = ~np.isnan(hcc) & ~np.isnan(fib) & (fib >= 0)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc[valid], fib[valid])
        validation_rows.append({
            "test": "HCC_score_vs_fibrosis_stage",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "positive = HCC risk increases with fibrosis",
        })
        log.info("1. vs fibrosis:  rho=%.4f, p=%.2e, n=%d", rho, pval, valid.sum())

    # Test 2: HCC score vs P(F4) (Spearman)
    valid = ~np.isnan(hcc) & ~np.isnan(pf4)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc[valid], pf4[valid])
        validation_rows.append({
            "test": "HCC_score_vs_P_F4",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "positive = HCC risk co-tracks with cirrhosis fate probability",
        })
        log.info("2. vs P(F4):     rho=%.4f, p=%.2e, n=%d", rho, pval, valid.sum())

    # Test 3: HCC score vs CPS (Spearman)
    valid = ~np.isnan(hcc) & ~np.isnan(cps)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc[valid], cps[valid])
        validation_rows.append({
            "test": "HCC_score_vs_CPS",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "positive = HCC risk tracks with composite progression",
        })
        log.info("3. vs CPS:       rho=%.4f, p=%.2e, n=%d", rho, pval, valid.sum())

    # Test 4: Kruskal-Wallis across fibrosis stages
    groups = []
    group_labels = []
    for stage in sorted(np.unique(fib[~np.isnan(fib) & (fib >= 0)])):
        mask = (fib == stage) & ~np.isnan(hcc)
        if mask.sum() > 0:
            groups.append(hcc[mask])
            group_labels.append(f"F{int(stage)}")
    if len(groups) >= 2:
        h_stat, h_pval = stats.kruskal(*groups)
        validation_rows.append({
            "test": "HCC_score_Kruskal_Wallis_by_fibrosis",
            "method": "Kruskal-Wallis",
            "statistic": h_stat,
            "p_value": h_pval,
            "n_samples": sum(len(g) for g in groups),
            "interpretation": f"Groups: {group_labels}; significant = HCC score differs by fibrosis",
        })
        log.info("4. KW fibrosis:  H=%.2f, p=%.2e, n=%d", h_stat, h_pval,
                 sum(len(g) for g in groups))
        for label, group in zip(group_labels, groups):
            log.info("   %s: mean=%.3f, n=%d", label, np.mean(group), len(group))

    # Test 5: S2 vs S1 (Mann-Whitney U)
    s1_mask = (subtype == "S1") & ~np.isnan(hcc)
    s2_mask = (subtype == "S2") & ~np.isnan(hcc)
    if s1_mask.sum() > 5 and s2_mask.sum() > 5:
        s1_scores = hcc[s1_mask]
        s2_scores = hcc[s2_mask]
        u_stat, u_pval = stats.mannwhitneyu(s2_scores, s1_scores, alternative="two-sided")
        validation_rows.append({
            "test": "HCC_score_S2_vs_S1",
            "method": "Mann-Whitney U",
            "statistic": u_stat,
            "p_value": u_pval,
            "n_samples": int(s1_mask.sum() + s2_mask.sum()),
            "interpretation": (
                f"S2 mean={np.mean(s2_scores):.3f} vs S1 mean={np.mean(s1_scores):.3f}; "
                f"S2 is female-enriched progressor subtype"
            ),
        })
        log.info("5. S2 vs S1:     U=%.1f, p=%.2e (S2=%.3f, S1=%.3f)",
                 u_stat, u_pval, np.mean(s2_scores), np.mean(s1_scores))

    # Test 6: HCC score vs NAS (Spearman)
    valid = ~np.isnan(hcc) & ~np.isnan(nas)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc[valid], nas[valid])
        validation_rows.append({
            "test": "HCC_score_vs_NAS",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "HCC risk vs NAFLD Activity Score",
        })
        log.info("6. vs NAS:       rho=%.4f, p=%.2e, n=%d", rho, pval, valid.sum())

    # ------------------------------------------------------------------
    # 8. Compare with old score
    # ------------------------------------------------------------------
    log.info("=== Step 8: Compare old vs new HCC score ===")
    old_score_path = os.path.join(INTEGRATION, "results/prognosis/hcc_molecular_score.csv")
    if os.path.exists(old_score_path):
        old_df = pd.read_csv(old_score_path)
        merged = score_df.merge(
            old_df[["sample_id", "hcc_molecular_score"]].rename(
                columns={"hcc_molecular_score": "hcc_old"}
            ),
            on="sample_id",
            how="inner",
        )
        if len(merged) > 10:
            rho_old_new, _ = stats.spearmanr(merged["hcc_molecular_score"], merged["hcc_old"])
            log.info("Old vs new score correlation: rho=%.4f (n=%d)", rho_old_new, len(merged))
            log.info("Old score: %d genes (from 3K subset)", int(n_old_mapped))
            log.info("New score: %d genes (from full 34K)", n_now_mapped)

    # ------------------------------------------------------------------
    # 9. Save outputs
    # ------------------------------------------------------------------
    log.info("=== Step 9: Saving outputs ===")

    # HCC molecular score v2
    score_path = os.path.join(OUT_DIR, "hcc_molecular_score_v2.csv")
    score_df.to_csv(score_path, index=False)
    log.info("Saved: %s (%d samples)", os.path.basename(score_path), len(score_df))

    # Validation results v2
    val_df = pd.DataFrame(validation_rows)
    val_path = os.path.join(OUT_DIR, "hcc_validation_v2.csv")
    val_df.to_csv(val_path, index=False)
    log.info("Saved: %s (%d tests)", os.path.basename(val_path), len(val_df))

    # Gene mapping report
    mapping_path = os.path.join(OUT_DIR, "hcc_gene_mapping_report.csv")
    mapping_df.to_csv(mapping_path, index=False)
    log.info("Saved: %s (%d genes)", os.path.basename(mapping_path), len(mapping_df))

    # Cleanup temp H5
    if os.path.exists(h5_tmp):
        os.remove(h5_tmp)
        log.info("Cleaned up temp H5")

    elapsed = time.time() - t0
    log.info("=" * 70)
    log.info("DONE in %.1f seconds", elapsed)
    log.info("  Genes mapped: %d / %d (was %d)", n_now_mapped, n_total, int(n_old_mapped))
    log.info("  Genes recovered: %d", n_recovered)
    log.info("  Genes still unmapped: %d", n_still_missing)
    if n_still_missing > 0:
        for _, r in unmapped.iterrows():
            log.info("    %s (PP.H4=%.3f)", r["gene"], r["PP.H4"])


if __name__ == "__main__":
    main()
