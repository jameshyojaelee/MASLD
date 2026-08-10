#!/usr/bin/env python3
"""Evidence-class spatial autocorrelation across the 35 HMSMA Visium samples.

QUESTION
--------
Do genetic-only and disease-state-only genes occupy different tissue contexts?
The paper's central result — that liver-trait genetics and the disease-state
transcriptome nominate largely non-overlapping genes — is currently a counting
statistic. If the classes are also *physically* organised differently in tissue,
the orthogonality becomes a mechanism rather than an overlap count.

WHY THIS NEEDS NO SAMPLE LABELS
-------------------------------
Evidence classes are defined by OUR data (SuSiE colocalization + bulk DE), not by
HMSMA phenotype. Spatial autocorrelation is a within-tissue property. So this runs
on unlabelled samples. What labels would later add is stage-dependence, not the
base comparison.

ESTIMAND
--------
Deliberately the same as the existing spatial arm (GSE192741, Vu):
`residual_spatial_autocorrelation_against_matched_gene_null`. That makes this a
direct replication at 35 samples rather than a new measure — the existing result
rests on 4 donors plus Vu with an unresolved donor count.

TWO DELIBERATE DIVERGENCES from `05c_spatially_variable_genes.py`, both required:
  1. It tests only `highly_variable` genes. That is an asymmetric denominator here:
     if one class is systematically less variable it would be under-tested and the
     comparison rigged. We test every adequately-detected gene.
  2. We match the null on expression AND detection rate. Moran's I rises with
     expression, so an unmatched comparison would recover expression differences
     between classes, not spatial biology.
`n_neighs=6` is kept (Visium hex lattice) for comparability.

UNIT OF REPLICATION
-------------------
The sample. Per-sample class-minus-matched-null deltas are computed first, then
tested across samples. Spots are not replicates.

CAVEATS carried into the output
-------------------------------
- No H&E: on-tissue status is a UMI-threshold proxy from the upstream build.
- No phenotype labels: no disease contrast, no stage stratification.
- Coordinates are array positions from the Visium v1 barcode map, not
  image-registered.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
from scipy import stats

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
IN_DIR = PROJ / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo"
CLASSES = (PROJ / "Analysis/Multimodal_Program_Projection/candidates"
           / "program-context-v2-candidate-2026-08-07/genetics_context"
           / "frozen_evidence_classes.tsv")
OUT_DIR = PROJ / "Analysis/Spatial/results/hmsma_class_spatial"

TEST_CLASSES = ["genetic_only", "disease_state_only", "convergent"]
NULL_CLASS = "neither"


def load_classes() -> pd.DataFrame:
    df = pd.read_csv(CLASSES, sep="\t")
    df = df[df["ensembl_bulk"].notna() & (df["ensembl_bulk"] != "")]
    df = df[["ensembl_bulk", "static_class"]].drop_duplicates("ensembl_bulk")
    return df.set_index("ensembl_bulk")["static_class"]


def per_sample_moran(a: ad.AnnData, min_frac: float) -> pd.Series:
    """Moran's I per gene for one sample. Returns Series indexed by unversioned ENSG."""
    a = a.copy()
    a.var_names = [g.split(".")[0] for g in a.var_names]
    a.var_names_make_unique()

    # detection filter — power for Moran's I, applied identically to every class
    det = np.asarray((a.X > 0).sum(axis=0)).ravel()
    keep = det >= max(3, int(min_frac * a.n_obs))
    a = a[:, keep].copy()
    if a.n_vars < 100:
        return pd.Series(dtype=float)

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    sq.gr.spatial_neighbors(a, coord_type="generic", n_neighs=6)
    sq.gr.spatial_autocorr(a, mode="moran", genes=a.var_names.tolist(),
                           n_perms=None, n_jobs=1)
    mi = a.uns["moranI"]["I"]

    # covariates for matching, on the same gene set
    expr = np.asarray(a.X.mean(axis=0)).ravel()
    detf = np.asarray((a.X > 0).mean(axis=0)).ravel()
    cov = pd.DataFrame({"expr": expr, "detf": detf}, index=a.var_names)
    out = pd.concat([mi.rename("I"), cov], axis=1).dropna()
    return out


def matched_delta(df: pd.DataFrame, cls: pd.Series, n_bins: int, seed: int) -> dict:
    """Class mean Moran's I minus expression/detection-matched null mean."""
    rng = np.random.default_rng(seed)
    d = df.join(cls.rename("cls"), how="inner")
    d["ebin"] = pd.qcut(d["expr"].rank(method="first"), n_bins, labels=False)
    d["dbin"] = pd.qcut(d["detf"].rank(method="first"), n_bins, labels=False)

    null_pool = d[d["cls"] == NULL_CLASS]
    res = {}
    for c in TEST_CLASSES:
        tgt = d[d["cls"] == c]
        if len(tgt) < 5:
            res[c] = dict(n=len(tgt), delta=np.nan, obs=np.nan, null=np.nan)
            continue
        drawn = []
        for (eb, db), grp in tgt.groupby(["ebin", "dbin"], observed=True):
            pool = null_pool[(null_pool["ebin"] == eb) & (null_pool["dbin"] == db)]
            if len(pool) == 0:
                continue
            k = min(len(pool), 20 * len(grp))
            drawn.append(pool["I"].to_numpy()[rng.choice(len(pool), k, replace=False)])
        null_vals = np.concatenate(drawn) if drawn else np.array([np.nan])
        res[c] = dict(n=int(len(tgt)), obs=float(tgt["I"].mean()),
                      null=float(np.nanmean(null_vals)),
                      delta=float(tgt["I"].mean() - np.nanmean(null_vals)))
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-detect-frac", type=float, default=0.05)
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cls = load_classes()
    print(f"class table: {len(cls)} genes; "
          f"{cls.value_counts().to_dict()}", flush=True)

    samples = sorted(p.stem for p in IN_DIR.glob("HRA_*.h5ad"))
    print(f"samples: {len(samples)}", flush=True)

    rows, per_gene = [], []
    for s in samples:
        a = ad.read_h5ad(IN_DIR / f"{s}.h5ad")
        df = per_sample_moran(a, args.min_detect_frac)
        if df.empty:
            print(f"  {s}: SKIP (too few genes)", flush=True)
            continue
        r = matched_delta(df, cls, args.n_bins, args.seed)
        for c, v in r.items():
            rows.append(dict(sample=s, cls=c, **v))
        g = df.join(cls.rename("cls"), how="inner")
        g["sample"] = s
        per_gene.append(g.reset_index().rename(columns={"index": "gene"}))
        print(f"  {s}: {df.shape[0]} genes tested | "
              + " ".join(f"{c}Δ={r[c]['delta']:+.4f}(n={r[c]['n']})" for c in TEST_CLASSES),
              flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(args.out_dir / "per_sample_class_delta.tsv", sep="\t", index=False)
    pd.concat(per_gene).to_csv(args.out_dir / "per_gene_morans_i.tsv.gz",
                               sep="\t", index=False, compression="gzip")

    # sample-level inference: one-sample t-test on per-sample deltas
    summary = []
    for c in TEST_CLASSES:
        d = res.loc[res["cls"] == c, "delta"].dropna()
        if len(d) < 3:
            continue
        t, p = stats.ttest_1samp(d, 0.0)
        w = stats.wilcoxon(d)[1] if len(d) >= 6 else np.nan
        summary.append(dict(cls=c, n_samples=int(len(d)),
                            mean_delta=float(d.mean()),
                            sd=float(d.std(ddof=1)),
                            ci_lo=float(d.mean() - 1.96 * d.std(ddof=1) / np.sqrt(len(d))),
                            ci_hi=float(d.mean() + 1.96 * d.std(ddof=1) / np.sqrt(len(d))),
                            t=float(t), p_ttest=float(p), p_wilcoxon=float(w),
                            n_samples_positive=int((d > 0).sum())))
    sm = pd.DataFrame(summary)
    sm.to_csv(args.out_dir / "class_summary.tsv", sep="\t", index=False)

    meta = dict(
        n_samples=len(samples), estimand="residual_spatial_autocorrelation_vs_matched_null",
        unit_of_replication="sample", n_neighs=6,
        min_detect_frac=args.min_detect_frac, match_bins=args.n_bins,
        hvg_restriction=False,
        phenotype_labels_available=False,
        tissue_mask="UMI>=500 proxy (no H&E)",
        coordinates="Visium v1 array positions (not image-registered)",
    )
    (args.out_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    print("\n=== CLASS SUMMARY ===", flush=True)
    print(sm.to_string(index=False), flush=True)
    print(f"\nWROTE {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
