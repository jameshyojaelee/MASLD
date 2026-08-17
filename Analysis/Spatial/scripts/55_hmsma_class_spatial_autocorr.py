#!/usr/bin/env python3
"""Evidence-class spatial autocorrelation across the 35 HMSMA Visium arrays.

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
on unlabelled arrays. What labels would later add is stage-dependence, not the
base comparison.

ESTIMAND AND UNIT OF INFERENCE
------------------------------
`residual_spatial_autocorrelation_against_matched_gene_null`, the same estimand as
the existing spatial arm (GSE192741, Vu).

The unit of inference is the MATCHED GENE SET, with the 35 arrays held fixed. It
is NOT the array. HMSMA has no array-to-donor key (see
`53_build_hmsma_anndata.py`), so we cannot tell whether 35 arrays represent 35
donors or fewer. Any test that treats arrays as independent replicates — a t-test
or Wilcoxon over the 35 per-array values — assumes exactly the donor independence
this dataset cannot establish, so no such test is computed here. Instead the
observed cohort statistic is compared against a null built by redrawing
expression- and detection-matched gene sets from the `neither` class and repeating
the identical per-array collapse. That null is over genes, not over donors.

Consequence, and it is a real limit: the resulting p-values are conditional on
these 35 arrays and do not license generalization to new donors. Per-array values
are reported descriptively only.

TWO DELIBERATE DIVERGENCES from `05c_spatially_variable_genes.py`, both required:
  1. It tests only `highly_variable` genes. That is an asymmetric denominator here:
     if one class is systematically less variable it would be under-tested and the
     comparison rigged. We test every adequately-detected gene.
  2. We match the null on expression AND detection rate. Moran's I rises with
     expression, so an unmatched comparison would recover expression differences
     between classes, not spatial biology.
`n_neighs=6` is kept (Visium hex lattice) for comparability.

CAVEATS carried into the output
-------------------------------
- No H&E: on-tissue status is a UMI-threshold proxy from the upstream build.
- No phenotype labels: no disease contrast, no stage stratification.
- No donor key: arrays are not donors and are never counted as such.
- Coordinates are array positions from the Visium v1 barcode map, not
  image-registered.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
IN_DIR = PROJ / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo"
CLASSES = (PROJ / "Analysis/Multimodal_Program_Projection/candidates"
           / "program-context-v2-candidate-2026-08-07/genetics_context"
           / "frozen_evidence_classes.tsv")
OUT_DIR = PROJ / "Analysis/Spatial/results/hmsma_class_spatial_v2"
PER_GENE_CACHE = (PROJ / "Analysis/Spatial/results/hmsma_class_spatial"
                  / "per_gene_morans_i.tsv.gz")

TEST_CLASSES = ["genetic_only", "disease_state_only", "convergent"]
NULL_CLASS = "neither"
CONTRAST = ("disease_state_only", "genetic_only")
MIN_CLASS_GENES = 5


def load_classes() -> pd.Series:
    df = pd.read_csv(CLASSES, sep="\t")
    df = df[df["ensembl_bulk"].notna() & (df["ensembl_bulk"] != "")]
    df = df[["ensembl_bulk", "static_class"]].drop_duplicates("ensembl_bulk")
    return df.set_index("ensembl_bulk")["static_class"]


def per_sample_moran(a, min_frac: float) -> pd.DataFrame:
    """Moran's I per gene for one array. Returns frame indexed by unversioned ENSG."""
    import scanpy as sc
    import squidpy as sq

    a = a.copy()
    a.var_names = [g.split(".")[0] for g in a.var_names]
    a.var_names_make_unique()

    # detection filter — power for Moran's I, applied identically to every class
    det = np.asarray((a.X > 0).sum(axis=0)).ravel()
    keep = det >= max(3, int(min_frac * a.n_obs))
    a = a[:, keep].copy()
    if a.n_vars < 100:
        return pd.DataFrame()

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
    return pd.concat([mi.rename("I"), cov], axis=1).dropna()


def matched_cells(df: pd.DataFrame, cls: pd.Series, n_bins: int) -> dict:
    """Per-array matched-cell structure for one array.

    Genes are stratified into expression x detection quantile cells. For each test
    class we keep only the genes whose cell contains at least one `neither` gene,
    so the observed mean and the matched expectation are computed over exactly the
    same genes. Retaining a class gene whose cell has no control would compare a
    filtered class against an unfiltered background — the asymmetric-denominator
    trap this analysis exists to avoid.
    """
    d = df.join(cls.rename("cls"), how="inner")
    d = d[d["cls"].isin([*TEST_CLASSES, NULL_CLASS])].copy()
    if len(d) < 10 * n_bins:
        return {}
    d["ebin"] = pd.qcut(d["expr"].rank(method="first"), n_bins,
                        labels=False, duplicates="drop")
    d["dbin"] = pd.qcut(d["detf"].rank(method="first"), n_bins,
                        labels=False, duplicates="drop")
    d = d.dropna(subset=["ebin", "dbin"])
    d["ebin"] = d["ebin"].astype(int)
    d["dbin"] = d["dbin"].astype(int)

    pools = {
        key: grp["I"].to_numpy(float)
        for key, grp in d[d["cls"] == NULL_CLASS].groupby(["ebin", "dbin"])
        if len(grp)
    }

    out = {}
    for c in TEST_CLASSES:
        tgt = d[d["cls"] == c]
        cells, observed, dropped = [], [], 0
        for key, grp in tgt.groupby(["ebin", "dbin"]):
            pool = pools.get(key)
            if pool is None or not len(pool):
                dropped += len(grp)
                continue
            cells.append((len(grp), pool))
            observed.append(grp["I"].to_numpy(float))
        n_total = int(sum(n for n, _ in cells))
        if n_total < MIN_CLASS_GENES:
            continue
        obs_values = np.concatenate(observed)
        # E[matched draw] equals this expectation exactly, so the null delta is
        # centered at zero by construction and needs no further recentering.
        expectation = float(
            sum(n * pool.mean() for n, pool in cells) / n_total
        )
        out[c] = {
            "cells": cells,
            "n_total": n_total,
            "n_dropped_no_control": int(dropped),
            "n_control_pool": int(sum(len(pool) for _, pool in cells)),
            "observed_mean": float(obs_values.mean()),
            "expectation": expectation,
            "delta": float(obs_values.mean() - expectation),
        }
    return out


def null_deltas(cells, n_total: int, expectation: float,
                n_null: int, rng: np.random.Generator) -> np.ndarray:
    """Matched-gene-set null for one array and one class.

    Each class gene is replaced by a control drawn uniformly from the `neither`
    genes in its own expression/detection cell. Draws are independent across
    genes, so a control may repeat within a set; that induces positive
    within-set correlation and mildly inflates the null spread, which is
    conservative for an upper-tail test.
    """
    total = np.zeros(n_null, dtype=float)
    for n_k, pool in cells:
        idx = rng.integers(0, len(pool), size=(n_null, n_k))
        total += pool[idx].sum(axis=1)
    return total / n_total - expectation


def upper_p(observed: float, null: np.ndarray) -> float:
    return float((1.0 + np.sum(null >= observed)) / (len(null) + 1.0))


def two_sided_p(observed: float, null: np.ndarray) -> float:
    return float((1.0 + np.sum(np.abs(null) >= abs(observed))) / (len(null) + 1.0))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-detect-frac", type=float, default=0.05)
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--n-null", type=int, default=9999)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--in-dir", type=Path, default=IN_DIR,
                    help="directory of per-array h5ad objects to recompute from")
    ap.add_argument("--from-per-gene", type=Path, default=None,
                    help="reuse a stored per_gene_morans_i.tsv.gz instead of "
                         "recomputing Moran's I from the h5ad objects")
    args = ap.parse_args()
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite populated results directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cls = load_classes()
    print(f"class table: {len(cls)} genes; {cls.value_counts().to_dict()}", flush=True)

    if args.from_per_gene is not None:
        cached = pd.read_csv(args.from_per_gene, sep="\t")
        per_gene = {
            s: part.set_index("gene")[["I", "expr", "detf"]]
            for s, part in cached.groupby("sample")
        }
        source = f"cached per-gene Moran values: {args.from_per_gene}"
    else:
        import anndata as ad

        per_gene = {}
        for path in sorted(args.in_dir.glob("HRA_*.h5ad")):
            df = per_sample_moran(ad.read_h5ad(path), args.min_detect_frac)
            if df.empty:
                print(f"  {path.stem}: SKIP (too few genes)", flush=True)
                continue
            per_gene[path.stem] = df
        source = f"recomputed from {args.in_dir}"
    samples = sorted(per_gene)
    print(f"arrays: {len(samples)} ({source})", flush=True)

    rng = np.random.default_rng(args.seed)
    per_array_rows = []
    obs_delta: dict[str, dict[str, float]] = {c: {} for c in TEST_CLASSES}
    null_delta: dict[str, dict[str, np.ndarray]] = {c: {} for c in TEST_CLASSES}

    for s in samples:
        cells_by_class = matched_cells(per_gene[s], cls, args.n_bins)
        if not cells_by_class:
            print(f"  {s}: SKIP (too few matchable genes)", flush=True)
            continue
        for c, info in cells_by_class.items():
            obs_delta[c][s] = info["delta"]
            null_delta[c][s] = null_deltas(
                info["cells"], info["n_total"], info["expectation"], args.n_null, rng
            )
            per_array_rows.append({
                "sample": s,
                "cls": c,
                "n": info["n_total"],
                "n_dropped_no_control": info["n_dropped_no_control"],
                "n_control_pool": info["n_control_pool"],
                "obs": info["observed_mean"],
                "null": info["expectation"],
                "delta": info["delta"],
            })
        print(f"  {s}: " + " ".join(
            f"{c}Δ={cells_by_class[c]['delta']:+.4f}(n={cells_by_class[c]['n_total']})"
            for c in TEST_CLASSES if c in cells_by_class), flush=True)

    per_array = pd.DataFrame(per_array_rows)
    per_array.to_csv(args.out_dir / "per_array_class_delta.tsv", sep="\t", index=False)

    # Cohort collapse: arrays are held fixed and averaged. The null repeats the
    # identical collapse on redrawn matched gene sets, so the sampling object is
    # the gene set and no donor independence is assumed anywhere.
    summary, cohort_null = [], {}
    for c in TEST_CLASSES:
        used = sorted(obs_delta[c])
        if len(used) < 2:
            continue
        observed = float(np.mean([obs_delta[c][s] for s in used]))
        null = np.mean(np.vstack([null_delta[c][s] for s in used]), axis=0)
        cohort_null[c] = null
        deltas = np.array([obs_delta[c][s] for s in used])
        summary.append(dict(
            cls=c,
            n_arrays=len(used),
            observed_delta=observed,
            null_mean=float(null.mean()),
            null_sd=float(null.std(ddof=1)),
            null_q050=float(np.quantile(null, 0.05)),
            null_q500=float(np.quantile(null, 0.50)),
            null_q950=float(np.quantile(null, 0.95)),
            p_matched_gene_null=upper_p(observed, null),
            n_null=args.n_null,
            sidedness="one_sided_upper",
            n_arrays_delta_positive=int((deltas > 0).sum()),
            median_array_delta=float(np.median(deltas)),
            iqr_array_delta=float(np.subtract(*np.percentile(deltas, [75, 25]))),
        ))
    sm = pd.DataFrame(summary)
    sm.to_csv(args.out_dir / "class_summary.tsv", sep="\t", index=False)

    # Within-array paired contrast, same gene-set null, arrays still fixed.
    hi, lo = CONTRAST
    contrast = pd.DataFrame()
    shared = sorted(set(obs_delta[hi]) & set(obs_delta[lo]))
    if len(shared) >= 2:
        paired = np.array([obs_delta[hi][s] - obs_delta[lo][s] for s in shared])
        observed = float(paired.mean())
        null = np.mean(np.vstack(
            [null_delta[hi][s] - null_delta[lo][s] for s in shared]), axis=0)
        contrast = pd.DataFrame([dict(
            contrast=f"{hi}_minus_{lo}",
            n_arrays=len(shared),
            observed_difference=observed,
            ratio_of_class_deltas=float(
                np.mean([obs_delta[hi][s] for s in shared])
                / np.mean([obs_delta[lo][s] for s in shared])),
            null_mean=float(null.mean()),
            null_q050=float(np.quantile(null, 0.05)),
            null_q950=float(np.quantile(null, 0.95)),
            p_matched_gene_null=two_sided_p(observed, null),
            n_null=args.n_null,
            sidedness="two_sided",
            n_arrays_difference_positive=int((paired > 0).sum()),
            median_array_difference=float(np.median(paired)),
        )])
        contrast.to_csv(args.out_dir / "contrast_summary.tsv", sep="\t", index=False)

    meta = dict(
        n_arrays=len(samples),
        estimand="residual_spatial_autocorrelation_vs_matched_gene_null",
        inference_unit="matched_gene_set_with_arrays_held_fixed",
        uncertainty_semantics="matched_gene_null_quantiles_not_sampling_confidence_interval",
        arrays_are_biological_donors=False,
        donor_identity_resolved=False,
        across_array_parametric_test_performed=False,
        across_array_test_note=(
            "No t-test or Wilcoxon over per-array values is computed: HMSMA has no "
            "array-to-donor key, so array independence cannot be established. "
            "Per-array values are descriptive only."
        ),
        generalization_note=(
            "p-values are conditional on these arrays and do not support "
            "generalization to new donors."
        ),
        n_null=args.n_null,
        seed=args.seed,
        n_neighs=6,
        min_detect_frac=args.min_detect_frac,
        match_bins=args.n_bins,
        match_covariates=["expression", "detection_rate"],
        hvg_restriction=False,
        phenotype_labels_available=False,
        tissue_mask="UMI>=500 proxy (no H&E)",
        coordinates="Visium v1 array positions (not image-registered)",
        moran_source=source,
    )
    (args.out_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    print("\n=== CLASS SUMMARY (matched-gene-set null) ===", flush=True)
    print(sm.to_string(index=False), flush=True)
    if not contrast.empty:
        print("\n=== CONTRAST ===", flush=True)
        print(contrast.to_string(index=False), flush=True)
    print(f"\nWROTE {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
