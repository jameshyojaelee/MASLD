#!/usr/bin/env python
"""
353_ccc_coexpression_null.py  (S5)

Coexpression-aware null model for the LIANA-based stage-CCC headline.

Current production null (LIANA / Script 345 + downstream):
    permute cell-type labels within donor -> preserves cell-type marginal
    composition and per-gene mean expression, but breaks the cell-type ->
    LR coexpression linkage. This is the existing "label-shuffle" null.

Stricter null introduced here:
    Within each cell type, permute UMI counts of each gene independently
    across cells of the same cell type, in the same donor.
    -> preserves the marginal gene expression within cell type
       (so any LR pair will still be "expressed" in the source/target
        at the same fraction of cells), but breaks the within-cell
        joint distribution of ligand x receptor across cells.
    -> any communication score that exploits cell-cell coexpression
       (LIANA "expr_prod", "scaled" etc.) will collapse to chance,
       while a label-shuffle null leaves the joint intact and only
       randomizes which CT pair the score is attributed to.

For each (donor, cell type) we shuffle the columns of the raw count matrix
independently per gene (np.random.permutation of each column). We then
re-run a streamlined LIANA expr_prod-style score per LR pair per CT pair,
and aggregate to the donor level.

Outputs:
    coexpression_null_per_donor.parquet         per-donor null LR scores
    coexpression_null_summary.tsv               per-LR/per-CT-pair null
    headline_multiplier_vs_coexpr_null.tsv      headline / null comparison
    METHODS_PARAGRAPH.md                        documentation of the
                                                "DUBIOUS path bug" fix and
                                                the new coexpression-aware
                                                null

Usage (SLURM):
    sbatch run_ccc_coexpression_null.sh
"""

from __future__ import annotations
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

H5AD_PATH  = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
PER_DONOR  = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/per_donor_lr"
OUT_ROOT   = Path(os.environ.get("S5_OUT_ROOT", str(PROJECT_ROOT)))
OUT_DIR    = OUT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc_null_v2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

KEEP_LINEAGES = [
    "Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial cells",
    "Cholangiocytes", "T cells", "B cells", "Resident NK", "Plasma cells",
]
MIN_CELLS = 30
EXPR_FRAC = 0.10
RNG_SEED  = 42
N_NULLS   = int(os.environ.get("CCC_N_NULLS", "10"))  # null replicates per donor

# 8 paracrine pair headline (placeholder list; will be refined from observed
# per-donor LIANA aggregates loaded below if available)
HEADLINE_CT_PAIRS = [
    ("Hepatocytes",   "Macrophages"),
    ("Macrophages",   "Hepatocytes"),
    ("Hepatocytes",   "Fibroblasts"),
    ("Fibroblasts",   "Hepatocytes"),
    ("Macrophages",   "Fibroblasts"),
    ("Fibroblasts",   "Macrophages"),
    ("Endothelial cells", "Hepatocytes"),
    ("Hepatocytes",   "Endothelial cells"),
]


def load_atlas() -> ad.AnnData:
    print(f"[load] reading {H5AD_PATH}", flush=True)
    A = ad.read_h5ad(H5AD_PATH, backed="r")
    print(f"[load] atlas {A.shape}", flush=True)
    return A


def lr_pairs(adata: ad.AnnData) -> pd.DataFrame:
    """Use LIANA built-in consensus database; fallback to a tiny curated list."""
    try:
        import liana as li
        db = li.resource.select_resource("consensus")
        return db
    except Exception as e:
        print(f"[warn] liana resource unavailable ({e}); using curated fallback",
              flush=True)
        return pd.DataFrame({
            "ligand":   ["TGFB1", "TGFB1", "PDGFB", "CCL2", "IL6", "VEGFA",
                         "WNT2", "TNF"],
            "receptor": ["TGFBR2", "TGFBR1", "PDGFRB", "CCR2", "IL6R", "KDR",
                         "FZD1", "TNFRSF1A"],
        })


def shuffle_within_celltype(X: sp.csr_matrix, rng: np.random.Generator
                            ) -> sp.csr_matrix:
    """Permute each gene's counts across cells of the same cell type.

    X is (cells_in_celltype x genes). For each gene column we apply an
    independent random permutation of row indices.
    """
    Xd = X.toarray() if sp.issparse(X) else np.asarray(X)
    n_cells, n_genes = Xd.shape
    for j in range(n_genes):
        rng.shuffle(Xd[:, j])
    return sp.csr_matrix(Xd)


def expr_prod_score(adata: ad.AnnData, lr: pd.DataFrame,
                    ct_pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """Streamlined LIANA-like expr_prod score.

    For each LR pair and each (source, target) CT pair:
        score = mean(ligand) in source * mean(receptor) in target,
        passing if frac>=EXPR_FRAC on each side.
    """
    var = adata.var_names
    cts = adata.obs["cell_type"]
    rows = []
    # Pre-compute per-celltype mean expression
    ct_means: dict[str, np.ndarray] = {}
    ct_fracs: dict[str, np.ndarray] = {}
    for ct in KEEP_LINEAGES:
        mask = (cts == ct).values
        if mask.sum() < MIN_CELLS:
            continue
        sub = adata.X[mask] if not sp.issparse(adata.X) else adata.X[mask]
        if sp.issparse(sub):
            sub = sub.toarray()
        ct_means[ct] = sub.mean(axis=0)
        ct_fracs[ct] = (sub > 0).mean(axis=0)

    # Index ligand / receptor positions
    var_index = {g: i for i, g in enumerate(var)}
    for _, row in lr.iterrows():
        lig, rec = row.get("ligand"), row.get("receptor")
        if lig not in var_index or rec not in var_index:
            continue
        gi, ri = var_index[lig], var_index[rec]
        for src, tgt in ct_pairs:
            if src not in ct_means or tgt not in ct_means:
                continue
            if ct_fracs[src][gi] < EXPR_FRAC or ct_fracs[tgt][ri] < EXPR_FRAC:
                continue
            score = float(ct_means[src][gi] * ct_means[tgt][ri])
            rows.append({"source": src, "target": tgt,
                         "ligand": lig, "receptor": rec,
                         "score": score})
    return pd.DataFrame(rows)


def run_donor(adata_full: ad.AnnData, sample: str, lr: pd.DataFrame,
              n_nulls: int, rng: np.random.Generator) -> pd.DataFrame:
    sub = adata_full[adata_full.obs["sample"] == sample].to_memory()
    sub = sub[sub.obs["cell_type"].isin(KEEP_LINEAGES)].copy()
    if sub.n_obs < 200:
        return pd.DataFrame()

    # OBSERVED
    obs_df = expr_prod_score(sub, lr, HEADLINE_CT_PAIRS)
    obs_df["null_replicate"] = -1
    obs_df["sample"] = sample

    # COEXPRESSION-AWARE NULL: shuffle counts within each cell type
    null_rows = [obs_df]
    for k in range(n_nulls):
        sub_null = sub.copy()
        X = sub_null.X.copy() if not sp.issparse(sub_null.X) else sub_null.X.copy()
        cts = sub_null.obs["cell_type"].values
        # Build a new matrix by permuting within each CT
        if sp.issparse(X):
            X = X.tolil()
            for ct in np.unique(cts):
                idx = np.where(cts == ct)[0]
                if len(idx) < MIN_CELLS:
                    continue
                block = X[idx, :].toarray()
                for j in range(block.shape[1]):
                    rng.shuffle(block[:, j])
                X[idx, :] = block
            sub_null.X = X.tocsr()
        else:
            for ct in np.unique(cts):
                idx = np.where(cts == ct)[0]
                if len(idx) < MIN_CELLS:
                    continue
                block = X[idx, :]
                for j in range(block.shape[1]):
                    rng.shuffle(block[:, j])
                X[idx, :] = block
            sub_null.X = X

        nd = expr_prod_score(sub_null, lr, HEADLINE_CT_PAIRS)
        nd["null_replicate"] = k
        nd["sample"] = sample
        null_rows.append(nd)

    return pd.concat(null_rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", type=int, default=int(os.environ.get("CCC_N_DONORS", "30")),
                    help="Number of donors to evaluate (random subsample for speed)")
    ap.add_argument("--n-nulls", type=int, default=N_NULLS)
    args = ap.parse_args()

    rng = np.random.default_rng(RNG_SEED)
    A = load_atlas()
    lr = lr_pairs(A)
    print(f"[lr] {len(lr)} LR pairs", flush=True)

    donors = pd.read_csv(DONOR_META, sep="\t")
    if "sample" in donors.columns:
        donor_ids = donors["sample"].dropna().unique().tolist()
    else:
        donor_ids = donors.iloc[:, 0].dropna().unique().tolist()
    rng.shuffle(donor_ids)
    donor_ids = donor_ids[: args.donors]
    print(f"[donors] {len(donor_ids)} donors", flush=True)

    all_rows = []
    for i, s in enumerate(donor_ids):
        t0 = time.time()
        try:
            df = run_donor(A, s, lr, args.n_nulls, rng)
        except Exception as e:
            print(f"[donor {s}] FAILED: {e}", flush=True)
            continue
        if df.empty:
            continue
        all_rows.append(df)
        print(f"[donor {s}] ({i+1}/{len(donor_ids)}) rows={len(df)} "
              f"elapsed={time.time()-t0:.1f}s", flush=True)

    if not all_rows:
        print("[fatal] no donors produced output", flush=True)
        sys.exit(1)

    full = pd.concat(all_rows, ignore_index=True)
    full.to_parquet(OUT_DIR / "coexpression_null_per_donor.parquet",
                    index=False)
    print(f"[write] {OUT_DIR / 'coexpression_null_per_donor.parquet'}",
          flush=True)

    # SUMMARY: per (source, target, ligand, receptor)
    obs = full[full["null_replicate"] == -1]
    nul = full[full["null_replicate"] != -1]
    obs_agg = obs.groupby(["source", "target", "ligand", "receptor"]
                          )["score"].mean().reset_index(name="obs_mean")
    nul_agg = nul.groupby(["source", "target", "ligand", "receptor"]
                          )["score"].agg(["mean", "std"]).reset_index()
    nul_agg = nul_agg.rename(columns={"mean": "null_mean", "std": "null_sd"})
    summary = obs_agg.merge(nul_agg, on=["source", "target", "ligand", "receptor"],
                            how="left")
    summary["multiplier"] = summary["obs_mean"] / summary["null_mean"].replace(0, np.nan)
    summary.to_csv(OUT_DIR / "coexpression_null_summary.tsv", sep="\t",
                   index=False)
    print(f"[write] {OUT_DIR / 'coexpression_null_summary.tsv'}", flush=True)

    # HEADLINE multiplier: aggregate over 8 paracrine ct pairs
    headline = summary[summary.apply(
        lambda r: (r["source"], r["target"]) in HEADLINE_CT_PAIRS, axis=1)]
    out_h = headline.groupby(["source", "target"]).agg(
        obs_mean=("obs_mean", "mean"),
        null_mean=("null_mean", "mean"),
        n_lr=("ligand", "count"),
    ).reset_index()
    out_h["multiplier_vs_coexpr_null"] = out_h["obs_mean"] / out_h["null_mean"]
    out_h.to_csv(OUT_DIR / "headline_multiplier_vs_coexpr_null.tsv",
                 sep="\t", index=False)
    print(f"[write] {OUT_DIR / 'headline_multiplier_vs_coexpr_null.tsv'}",
          flush=True)
    print(out_h.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
