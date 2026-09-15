#!/usr/bin/env python3
"""Step 75 (P3): is the allelic concordance robust to the het-calling filters?

The het call is the headline's main free choice: a site enters only if enough reads carry each allele in
enough donors. This sweeps the three thresholds over a fixed grid and reports the concordance, the
marginal-skew expectation it has to beat, and the LD-respecting block sign test at every combination. No
threshold is selected here; the grid is reported whole so the reader can see whether the result depends on
the filter or survives it.

Outputs (tables/): allelic_filter_sensitivity.tsv
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import math

import numpy as np
import pandas as pd

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
GRID_TOTAL = (10, 15, 20)
GRID_ALLELE = (3, 4)
GRID_DONORS = (3, 4, 5)
COHORTS = (("gse281367", "allelic_donor_counts.tsv.gz"), ("gse244832", "gse244832_allelic_donor_counts.tsv.gz"))
MIN_SITES = 20


def _step73():
    spec = importlib.util.spec_from_file_location("s73", la.SCRIPT_DIR / "73_atac_allelic_analysis.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    m = _step73()
    preds = {}
    with gzip.open(m.LIVER_SUMMARIES, "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["scorer"] == "ATAC" and r["liver_median_quantile"] not in ("", "nan"):
                preds[r["variant_uid"]] = float(r["liver_median_quantile"])
    rows = []
    for cohort, cfile in COHORTS:
        path = TABLES / cfile
        if not path.exists():
            continue
        raw = pd.read_csv(path, sep="\t")
        for min_total in GRID_TOTAL:
            for min_allele in GRID_ALLELE:
                d = raw.copy()
                d["het"] = [m.call_het(a, b, c, min_total=min_total, min_allele=min_allele)
                            for a, b, c in zip(d.n_ref, d.n_alt, d.n_other)]
                d = d[d.het]
                for min_don in GRID_DONORS:
                    pr, me, bl = [], [], []
                    for uid, x in d.groupby("uid"):
                        if len(x) < min_don:
                            continue
                        q = preds.get(uid)
                        if q is None or q == 0:
                            continue
                        st = m.site_statistic(list(zip(x.n_ref, x.n_alt)))
                        if np.sign(st["mean_log2"]) == 0:
                            continue
                        c, p, _, _ = uid.split(":")
                        pr.append(np.sign(q)); me.append(np.sign(st["mean_log2"])); bl.append(f"{c}~{int(p) // 1_000_000}")
                    row = {"cohort": cohort, "min_reads_total": min_total, "min_reads_per_allele": min_allele,
                           "min_het_donors": min_don, "n_sites": len(pr)}
                    if len(pr) < MIN_SITES:
                        row["verdict"] = f"too few sites (< {MIN_SITES})"
                        rows.append(row); continue
                    pr, me, bl = np.array(pr), np.array(me), np.array(bl)
                    b = m.block_level_concordance(pr, me, bl)
                    row.update({"concordance": float(np.mean(pr == me)),
                                "marginal_expected": la.marginal_expected_concordance(pr, me),
                                "n_blocks_above_half": b["n_blocks_above_half"], "n_blocks_decided": b["n_blocks_decided"],
                                "block_sign_test_p": b["sign_test_p"], "per_block_mean": b["mean_per_block"],
                                "verdict": "reported"})
                    rows.append(row)
    la.write_tsv_once(TABLES / "allelic_filter_sensitivity.tsv", rows, sorted({k for r in rows for k in r}))
    ok = [r for r in rows if r.get("verdict") == "reported"]
    for cohort, _ in COHORTS:
        sub = [r for r in ok if r["cohort"] == cohort]
        if not sub:
            continue
        la.log(f"{cohort}: {len(sub)} filter combinations; concordance "
               f"{min(r['concordance'] for r in sub):.3f}-{max(r['concordance'] for r in sub):.3f}; "
               f"block sign-test p {min(r['block_sign_test_p'] for r in sub):.1e}-{max(r['block_sign_test_p'] for r in sub):.1e}")


if __name__ == "__main__":
    main()
