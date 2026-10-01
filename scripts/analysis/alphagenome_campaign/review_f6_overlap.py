#!/usr/bin/env python3
"""Independent reconstruction of the F6 genetic-overlap contrasts.

Re-derives GC, promoter status and donor signal from their sources by a separate
path, rebuilds the group assignment from the overlap table, and recomputes the
matched-draw contrast without importing the producer's measured_contrast. It also
reports the UNMATCHED difference beside the matched one, because a matched
contrast that equals its unmatched version did no matching, and it reports how
far the group and its matched controls sit apart on each matching covariate
before and after matching.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

import data_p4 as p4
import data_p4_families as fam

OVERLAPS = p4.ROOT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/variant_peak_overlaps.tsv.gz"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
N_DRAWS = 10000
SEED = 20260917


def gc_of(frame: pd.DataFrame) -> np.ndarray:
    import pysam
    genome = pysam.FastaFile(FASTA)
    try:
        return np.array([(s.count("G") + s.count("C")) / max(len(s), 1) for s in
                         (genome.fetch(str(r.chrom), int(r.start0), int(r.end)).upper()
                          for r in frame.itertuples())])
    finally:
        genome.close()


def signal_of(cohort: str, lineage: str, keys: set) -> pd.DataFrame:
    """logCPM mean/sd per peak, library size from every source peak, re-derived here."""
    base = fam.CTX / "counts" / cohort / lineage
    peaks = pd.read_csv(f"{base}.peaks.tsv", sep="\t")
    donors = pd.read_csv(f"{base}.donors.tsv", sep="\t")
    m = mmread(f"{base}.mtx.gz").tocsr()
    lib = np.asarray(m.sum(axis=1)).ravel()
    keep = donors.contrast_eligible.astype(str).str.upper().eq("TRUE").to_numpy() & (lib > 0)
    take = peaks.peak_coordinate.isin(keys).to_numpy()
    x = m[keep][:, take].toarray().astype(float)
    lc = np.log2(1 + x / lib[keep, None] * 1e6)
    return pd.DataFrame({"region_key": peaks.loc[take, "peak_coordinate"].to_numpy(),
                         "signal_mean": lc.mean(axis=0), "signal_sd": lc.std(axis=0, ddof=1)})


def deciles(x: np.ndarray, k: int) -> np.ndarray:
    e = np.quantile(x[np.isfinite(x)], np.linspace(0, 1, k + 1))[1:-1]
    return np.digitize(x, e, right=False)


def contrast(values: np.ndarray, strata: np.ndarray, group: np.ndarray, covars: dict) -> dict:
    pool = {lab: np.where((~group) & (strata == lab))[0] for lab in np.unique(strata[~group])}
    pool = {k: v for k, v in pool.items() if v.size}
    usable = group & np.isfinite(values) & np.array([lab in pool for lab in strata])
    obs = float(np.nanmedian(values[usable]))
    rng = np.random.default_rng(SEED)
    labs = strata[usable]
    order = sorted(pool)
    flat = np.concatenate([pool[lab] for lab in order])
    start, size = {}, {}
    at = 0
    for lab in order:
        start[lab], size[lab] = at, pool[lab].size
        at += pool[lab].size
    st = np.array([start[lab] for lab in labs])
    sz = np.array([size[lab] for lab in labs])
    draws = np.empty(N_DRAWS)
    covar_sum = {name: 0.0 for name in covars}
    for d in range(N_DRAWS):
        picks = flat[st + (rng.random(st.size) * sz).astype(int)]
        draws[d] = float(np.nanmedian(values[picks]))
        for name, v in covars.items():
            covar_sum[name] += float(np.nanmean(v[picks]))
    centre = float(np.median(draws))
    extreme = int(np.sum(np.abs(draws - centre) >= abs(obs - centre)))
    unmatched = float(np.nanmedian(values[(~group) & np.isfinite(values)]))
    out = {"n_group_usable": int(usable.sum()), "n_dropped_no_match": int(group.sum() - usable.sum()),
           "observed_group": obs, "matched_null_median": centre,
           "matched_difference": obs - centre,
           "unmatched_control_median": unmatched,
           "unmatched_difference": obs - unmatched,
           "p_nominal": (extreme + 1) / (N_DRAWS + 1),
           "p_floor": 1.0 / (N_DRAWS + 1)}
    for name, v in covars.items():
        out[f"{name}_group"] = float(np.nanmean(v[usable]))
        out[f"{name}_unmatched_control"] = float(np.nanmean(v[~group]))
        out[f"{name}_matched_control"] = covar_sum[name] / N_DRAWS
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reduced", type=Path, nargs="+", required=True)
    ap.add_argument("--producer", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    feats = pd.concat([pd.read_csv(r / "saturation_region_features.tsv.gz", sep="\t",
                                   keep_default_na=False, low_memory=False)
                       for r in args.reduced], ignore_index=True).drop_duplicates("region_key")
    overlap_peaks = set(pd.read_csv(OVERLAPS, sep="\t", usecols=["peak"], dtype=str).peak)
    da = pd.read_csv(fam.CTX / "da/da_peak_results.tsv.gz", sep="\t", keep_default_na=False)
    manifest = pd.read_csv(fam.CTX / "consensus_peak_manifest.tsv", sep="\t", dtype=str)
    manifest = manifest.loc[~manifest.blacklist_overlap.str.upper().eq("TRUE")]
    keys = set(manifest.lineage + "|" + manifest.chrom + ":" + manifest.start0 + "-" + manifest.end)
    da = da.loc[(da.lineage + "|" + da.peak_coordinate).isin(keys)]

    rows = []
    for lineage in ("hepatocyte", "stellate"):
        sub = da.loc[da.lineage.eq(lineage)]
        d0 = feats.loc[feats.region_key.isin(set(sub.peak_coordinate))].copy().reset_index(drop=True)
        d0["gc"] = gc_of(d0)
        d0["width"] = d0["end"].astype(int) - d0["start0"].astype(int)
        promo = sub[["peak_coordinate", "promoter_genes"]].drop_duplicates("peak_coordinate")
        d0 = d0.merge(promo, left_on="region_key", right_on="peak_coordinate", how="left")
        d0["promoter"] = d0.promoter_genes.fillna("").apply(
            lambda v: int(str(v) not in ("", "NA", "none", "None", "nan")))
        for cohort in ("GSE244832", "GSE281367"):
            d = d0.merge(signal_of(cohort, lineage, set(d0.region_key)),
                         on="region_key", validate="one_to_one").reset_index(drop=True)
            strata = np.array(["|".join(map(str, t)) for t in zip(
                deciles(d.gc.to_numpy(), 10), deciles(d.width.to_numpy(), 10),
                deciles(d.signal_mean.to_numpy(), 5), deciles(d.signal_sd.to_numpy(), 5),
                d.promoter.to_numpy())])
            group = d.region_key.isin(overlap_peaks).to_numpy()
            covars = {"gc": d.gc.to_numpy(), "width": d.width.to_numpy().astype(float),
                      "signal_mean": d.signal_mean.to_numpy(), "signal_sd": d.signal_sd.to_numpy(),
                      "promoter": d.promoter.to_numpy().astype(float)}
            for tg in ("liver_atac", "liver_dnase_adult_verified", "liver_h3k27ac_adult_verified"):
                col = f"{tg}_rel_mean"
                if col not in d:
                    continue
                r = contrast(pd.to_numeric(d[col], errors="coerce").to_numpy(), strata, group, covars)
                rows.append({"cohort": cohort, "lineage": lineage, "track_group": tg,
                             "scale": "rel", **r})
    review = pd.DataFrame(rows)
    review.to_csv(args.out / "F6_independent_contrasts.tsv", sep="\t", index=False)

    prod = pd.read_csv(args.producer / "F6_genetic_overlap.tsv", sep="\t")
    prod = prod.loc[prod.contrast.eq("overlap_vs_background") & prod.scale.eq("rel")]
    merged = review.merge(prod[["cohort", "lineage", "track_group", "observed_group",
                                "null_median", "observed_minus_null", "p_nominal"]],
                          on=["cohort", "lineage", "track_group"], suffixes=("_review", "_producer"))
    merged["abs_gap_observed"] = (merged.observed_group_review - merged.observed_group_producer).abs()
    merged["abs_gap_difference"] = (merged.matched_difference - merged.observed_minus_null).abs()
    merged.to_csv(args.out / "F6_review_vs_producer.tsv", sep="\t", index=False)

    matching_bites = bool((
        (review.matched_difference - review.unmatched_difference).abs() > 1e-9).all())
    sign_agrees = bool((np.sign(merged.matched_difference) ==
                        np.sign(merged.observed_minus_null)).all())
    # This reconstruction derives its own region pool, so its decile edges are its
    # own and stratum membership differs at the bin edges; the group sizes differ by
    # a few regions and exact equality is not the right test. The test is that the
    # signs agree and that the remaining gap is small beside the effects themselves.
    gap = float(merged.abs_gap_observed.max())
    effect = float(review.matched_difference.abs().median())
    reversed_by_matching = int((np.sign(review.unmatched_difference)
                                != np.sign(review.matched_difference)).sum())
    checks = {
        "status": "pass" if (sign_agrees and matching_bites and gap < 1e-3) else "fail",
        "agreement_rule": ("sign agrees on every contrast, matching changes every contrast, and the "
                           "largest reconstruction gap stays below 1e-3; exact equality is not required "
                           "because the reconstruction builds its own pool and its own decile edges"),
        "largest_gap_over_median_effect": gap / effect if effect else None,
        "contrasts_whose_sign_matching_reverses": reversed_by_matching,
        "sign_reversal_meaning": ("the unmatched control pool is promoter-defined background, so before "
                                  "matching the overlap group looks MORE sensitive; holding promoter "
                                  "status, GC, width and measured signal fixed reverses that"),
        "contrasts_reconstructed": int(len(review)),
        "max_abs_gap_observed_group": float(merged.abs_gap_observed.max()),
        "max_abs_gap_matched_difference": float(merged.abs_gap_difference.max()),
        "sign_agrees_on_every_contrast": sign_agrees,
        "matching_bites": matching_bites,
        "matching_bites_meaning": ("the matched difference differs from the unmatched difference on every "
                                   "contrast; equality would mean the strata did nothing"),
        "median_abs_matched_difference": float(review.matched_difference.abs().median()),
        "median_abs_unmatched_difference": float(review.unmatched_difference.abs().median()),
        "n_draws": N_DRAWS, "seed": SEED, "p_floor": 1.0 / (N_DRAWS + 1),
        "covariates_rederived_here": ["gc from the reference FASTA", "width from the keyed interval",
                                      "promoter from the source DA annotation",
                                      "signal_mean/signal_sd from the cohort count matrix"],
        "adopted": False,
    }
    (args.out / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2), flush=True)


if __name__ == "__main__":
    main()
