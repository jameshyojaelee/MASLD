#!/usr/bin/env python3
"""Reproduce the reference application's per-cohort call from its STORED matches.

NOT part of the tool. This is an adapter for one specific application (a BulkFormer /
PreBULK audit) kept as provenance, and `exposure_scan.py` neither imports nor knows about
it. It skips cleanly when the producing lane's artifacts are not on disk.

WHAT IT VERIFIES, AND WHAT IT DOES NOT. The 581,503 x 20,010 PreBULK corpus is no longer
on disk, so the SCAN cannot be re-run here. What the producing lane did keep is every
query sample's matched corpus row and max correlation. Feeding those into this release's
own statistics re-derives the entire DECISION layer -- window, injectivity, aligned run,
constant offset, permutation null, the same-sample band, and the eight per-cohort states --
from code written independently of the producing lane's scripts.

  verified here            everything downstream of the scan, against a real 46 GB corpus
  verified by the fixture  the scan itself, including chunked streaming and the argmax

Run:  python3 reproduce_reference_call.py [--bench <MASLD_Model_Benchmark dir>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import exposure_scan as X  # noqa: E402

LANE = "executions/bulkformer-feasibility-20260905T000000Z/results"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default=str(HERE.parents[2]),
                    help="the MASLD_Model_Benchmark directory holding executions/")
    ap.add_argument("--out", default=str(HERE / "reference_call_reproduction.json"))
    a = ap.parse_args()
    lane = Path(a.bench) / LANE
    raw, inp = lane / "overlap_raw.npz", lane / "bf_inputs.npz"
    published = lane / "bijection_test.json"
    for p in (raw, inp, published):
        if not p.exists():
            print(f"SKIP: producing-lane artifact not present at {p}")
            return 0

    d = np.load(raw)
    bi = np.load(inp, allow_pickle=True)
    ref = json.loads(published.read_text())
    coh = np.asarray(bi["pool_cohort"]).astype(str)

    groups = {}
    for c in sorted(set(coh.tolist())):
        m = coh == c
        groups[c] = (d["POOL_V2_topk_idx"][m, 0].astype(np.int64), d["POOL_V2_max"][m])
    for t in ("GSE268273", "GSE213621", "GSE276114"):
        groups[t] = (d[f"{t}_topk_idx"][:, 0].astype(np.int64), d[f"{t}_max"])

    crit = X.Criterion()
    rng = np.random.default_rng(X.DEFAULT_SEED)
    per = {}
    for name in sorted(groups):
        idx, mx = groups[name]
        n = len(idx)
        run, off = X.largest_constant_offset_run(idx)
        null = X.permutation_null(idx, crit.n_permutations, rng)
        rec = {"n": n, "largest_aligned_run": run, "aligned_fraction": round(run / n, 4),
               "constant_offset": off, "n_permutations": crit.n_permutations,
               "n_exceedances": int((null >= run).sum()),
               "perm_null_mean_run": float(null.mean()),
               "analytic_expected_max_run": round(X.analytic_expected_max_run(idx), 4),
               "max_correlation_median": float(np.median(mx))}
        rec.update(X.block_stats(idx))
        per[name] = rec

    band = X.calibrate_band(per, groups, crit)
    lo = band["lower_bound_applied"]
    state = {}
    for c, v in per.items():
        leg1 = (v["distinct_fraction"] >= crit.injectivity_min
                and v["window_width_over_n"] <= crit.window_over_n_max)
        leg2 = v["max_correlation_median"] >= lo
        v["leg1_structure_passes"], v["leg2_magnitude_passes"] = bool(leg1), bool(leg2)
        state[c] = "encoder_seen" if (leg1 and leg2) else "unknown"

    # ---- compare against what the producing lane published --------------------------
    diffs, exact, tol = [], 0, 0
    for c, v in per.items():
        r = ref["cohorts"][c]
        for k, rk in (("largest_aligned_run", "largest_aligned_run"),
                      ("constant_offset", "constant_offset"),
                      ("n_distinct_corpus_rows", "n_distinct_corpus_rows"),
                      ("narrowest_window_holding_90pct", "narrowest_window_holding_90pct"),
                      ("span", "span")):
            exact += 1
            if v[k] != r[rk]:
                diffs.append(f"{c}.{k}: {v[k]} vs published {r[rk]}")
        for k, rk, t in (("max_correlation_median", "max_correlation_median", 1e-12),
                         ("window_width_over_n", "window_width_over_n", 5e-4),
                         ("distinct_fraction", "distinct_fraction", 5e-5),
                         # a different RNG stream order: sampling noise on 2000 draws only
                         ("perm_null_mean_run", "perm_null_mean_run", 0.30)):
            tol += 1
            if abs(v[k] - r[rk]) > t:
                diffs.append(f"{c}.{k}: {v[k]} vs published {r[rk]} (tol {t})")
    for c, s in state.items():
        exact += 1
        if s != ref["exposure_state_measured"][c]:
            diffs.append(f"{c}.state: {s} vs published {ref['exposure_state_measured'][c]}")
    rb = ref["same_sample_cross_pipeline_band"]
    for k, rk in (("min", "min"), ("median", "median"), ("max", "max")):
        exact += 1
        if abs(band[k] - rb[rk]) > 1e-12:
            diffs.append(f"band.{k}: {band[k]} vs published {rb[rk]}")

    out = {
        "what": "the reference application's decision layer, re-derived by this release "
                "from the producing lane's stored per-sample matches",
        "not_verified_here": "the corpus scan itself; PreBULK.h5ad is no longer on disk. "
                             "The scan is verified on the synthetic fixture instead.",
        "n_exact_comparisons": exact, "n_tolerance_comparisons": tol,
        "disagreements": diffs,
        "reproduces_published_call": not diffs,
        "same_sample_band": band,
        "cohorts": per,
        "exposure_state_measured": state,
        "n_encoder_seen": sum(s == "encoder_seen" for s in state.values()),
        "n_unknown": sum(s == "unknown" for s in state.values()),
        "n_clean": 0,
    }
    Path(a.out).write_text(json.dumps(out, indent=2) + "\n")
    for c in sorted(state):
        v = per[c]
        print(f"  {c:<12} {state[c]:<13} n={v['n']:>3} run {v['largest_aligned_run']:>3}"
              f" off {v['constant_offset']:>7} win {v['window_width_over_n']:.3f}x"
              f" inj {v['distinct_fraction']:.4f} med r {v['max_correlation_median']:.4f}"
              f" exc {v['n_exceedances']}/{v['n_permutations']}")
    print(f"  band  {band['lower_bound_applied']:.6f} lower bound from "
          f"{band['from_cohorts']}, n={band['n_samples']}")
    print(f"  {exact} exact + {tol} tolerance comparisons against the published call")
    if diffs:
        print("  DISAGREEMENTS:")
        for x in diffs:
            print(f"    {x}")
        return 1
    print("  reproduces the published call with 0 disagreements")
    return 0


if __name__ == "__main__":
    sys.exit(main())
