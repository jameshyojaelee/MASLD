#!/usr/bin/env python3
"""C1 endpoint 2 v2, step 6: two checks the main recompute did not close.

1. The deposited-calibration comparison. `c1e2_05_recompute_v2.py` wrote an EMPTY
   `deposited_calibration_keys` because `currin_external_calibration.json` nests its Spearman and
   sign concordance under `point_estimates.direction`, and the filter only looked at top-level
   scalars. The comparison itself is made here, value against value.

2. A second code path for the two Holm-family C1 statistics. The headline contrasts are
   re-derived from the DEPOSITED row-level tables (`endpoint1_matched_tierA4.tsv.gz`,
   `endpoint2_matched.tsv`) rather than from the in-memory frames, with the same seed and the same
   block definition. A statistic that does not survive a reload of its own deposit is not a
   statistic anyone else can reproduce.

Writes `tables/qc_addendum.json`. Adds nothing to and overwrites nothing in the main deposit.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CALIB_JSON = os.path.join(
    ROOT, "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/"
          "currin_calibration/currin_external_calibration.json")
BOOT_SEED = 20260914
N_BOOT = 10000


def _pearson(a, b):
    a = a - a.mean()
    b = b - b.mean()
    da, db = np.sqrt((a * a).sum()), np.sqrt((b * b).sum())
    return np.nan if da <= 0 or db <= 0 else float((a * b).sum() / (da * db))


def signed_spearman(s, y):
    return _pearson(stats.rankdata(s), stats.rankdata(y))


def weighted_spearman(s, y, w):
    rs, rl, sw = stats.rankdata(s), stats.rankdata(y), w.sum()
    ms, ml = (w * rs).sum() / sw, (w * rl).sum() / sw
    cov = (w * (rs - ms) * (rl - ml)).sum() / sw
    vs = (w * (rs - ms) ** 2).sum() / sw
    vl = (w * (rl - ml) ** 2).sum() / sw
    return float(cov / np.sqrt(vs * vl))


def block_draws(blocks, n_boot, seed):
    uniq, starts, sizes = np.unique(blocks, return_index=True, return_counts=True)
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        s, st = sizes[pick], starts[pick]
        ends = np.cumsum(s)
        offs = np.arange(int(s.sum())) - np.repeat(ends - s, s)
        yield np.repeat(st, s) + offs


def contrast(df, a, b, label_col, weight_col, seed=BOOT_SEED):
    df = df.sort_values("block_1mb", kind="mergesort").reset_index(drop=True)
    blocks = df.block_1mb.values
    y = df[label_col].values.astype(float)
    w = df[weight_col].values.astype(float) if weight_col else None
    sa, sb = df[a].values.astype(float), df[b].values.astype(float)
    f = (lambda s, yy, ww: weighted_spearman(s, yy, ww)) if w is not None else \
        (lambda s, yy, ww: signed_spearman(s, yy))
    point = f(sa, y, w) - f(sb, y, w)
    d = np.empty(N_BOOT)
    for i, idx in enumerate(block_draws(blocks, N_BOOT, seed)):
        ww = w[idx] if w is not None else None
        d[i] = f(sa[idx], y[idx], ww) - f(sb[idx], y[idx], ww)
    p = min(1.0, 2.0 * min((d <= 0).mean(), (d >= 0).mean()))
    floor = 1.0 / N_BOOT
    return dict(point=float(point), ci_lo=float(np.percentile(d, 2.5)),
                ci_hi=float(np.percentile(d, 97.5)),
                boot_p_two_sided=float(max(p, floor)),
                p_floored=bool(p < floor), n=int(len(df)),
                n_blocks=int(pd.unique(blocks).size))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out, tab = args.outdir, os.path.join(args.outdir, "tables")
    res = {}

    # ---- 1. deposited calibration, value against value ----------------------
    with open(CALIB_JSON) as fh:
        calib = json.load(fh)
    dep = calib["point_estimates"]["direction"]
    t = pd.read_csv(os.path.join(tab, "endpoint1_matched_tierA4.tsv.gz"), sep="\t")
    # the deposit's number is over ALL Currin leads, so it is recomputed from the full
    # endpoint-2 score table's Currin join instead of the matched subset
    e2 = pd.read_csv(os.path.join(tab, "endpoint2_ase_union_scores.tsv"), sep="\t")
    res["deposited_direction"] = {k: dep[k] for k in ("spearman", "sign_concordance",
                                                      "n_oriented_positive")}
    res["note_scope"] = ("the deposited direction values are over all 32,336 Currin leads; the "
                         "loader reproduction of them is in qc_checks.json as "
                         "loader_reproduction_full_currin_*")
    qc = json.load(open(os.path.join(out, "qc_checks.json")))
    res["loader_vs_deposit_spearman_abs_diff"] = abs(
        qc["loader_reproduction_full_currin_spearman"] - dep["spearman"])
    res["loader_vs_deposit_sign_concordance_abs_diff"] = abs(
        qc["loader_reproduction_full_currin_sign_concordance"] - dep["sign_concordance"])
    res["loader_reproduction_full_currin_spearman"] = qc[
        "loader_reproduction_full_currin_spearman"]
    res["loader_reproduction_full_currin_sign_concordance"] = qc[
        "loader_reproduction_full_currin_sign_concordance"]
    res["endpoint2_union_rows_loaded"] = int(len(e2))

    # ---- 2. second code path for both Holm C1 members -----------------------
    res["endpoint1_tierA4_reload"] = contrast(
        t, "chrombpnet_adult_hep", "alphagenome_atac_liver", "beta_alt", None)
    m2 = pd.read_csv(os.path.join(tab, "endpoint2_matched.tsv"), sep="\t")
    res["endpoint2_weighted_reload"] = contrast(
        m2, "chrombpnet_adult_hep", "alphagenome_atac_liver",
        "mean_log2_alt_over_ref", "precision_weight")
    res["endpoint2_borzoi_atac_vs_atlas_weighted_reload"] = contrast(
        m2, "borzoi_atac", "alphagenome_atac_liver",
        "mean_log2_alt_over_ref", "precision_weight")

    # compare against what the main run deposited
    con = pd.read_csv(os.path.join(tab, "c1e2_contrasts.tsv"), sep="\t")

    def dep_row(ep, c, kind):
        s = con[(con.endpoint == ep) & (con.contrast == c) & (con.weighting == kind)]
        return s.iloc[0]

    prim = "chrombpnet_adult_hep - alphagenome_atac_liver"
    pairs = [("endpoint1_tierA4_reload", "endpoint1_currin_tierA4", prim, "unweighted"),
             ("endpoint2_weighted_reload", "endpoint2_allelic_imbalance", prim, "weighted"),
             ("endpoint2_borzoi_atac_vs_atlas_weighted_reload", "endpoint2_allelic_imbalance",
              "borzoi_atac - alphagenome_atac_liver", "weighted")]
    agree = {}
    for key, ep, c, kind in pairs:
        r = dep_row(ep, c, kind)
        agree[key] = dict(
            point_abs_diff=abs(res[key]["point"] - float(r.point)),
            ci_lo_abs_diff=abs(res[key]["ci_lo"] - float(r.ci_lo)),
            ci_hi_abs_diff=abs(res[key]["ci_hi"] - float(r.ci_hi)),
            p_abs_diff=abs(res[key]["boot_p_two_sided"] - float(r.boot_p_two_sided)))
    res["reload_vs_deposit_abs_diffs"] = agree
    res["reload_all_agree_below_1e-12"] = bool(all(
        v <= 1e-12 for a in agree.values() for v in a.values()))

    with open(os.path.join(tab, "qc_addendum.json"), "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print(json.dumps(res, indent=2, default=float), flush=True)


if __name__ == "__main__":
    main()
