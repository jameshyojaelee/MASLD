#!/usr/bin/env python3
"""Step 64b: re-derive step 64's paired indel-vs-SNV statistics by a separate code path.

Nothing is imported from steps 63 or 64. The effects table is read with pandas, 1-Mb blocks are recomputed
from coordinates, the sign test is scipy's exact binomial test and BH is scipy's false_discovery_control.
Every count, median, p and q must reproduce exactly. The bootstrap interval is re-drawn with its own
generator, so it is compared within Monte Carlo error rather than for equality.

It also compares the rerun with the run it replaces: a target or control scored on the same alleles in both
runs must carry identical predictions if the model API is deterministic, so any difference there is either
nondeterminism or a changed input, and is reported either way.

Outputs (tables/ of the rescue run): indel_vs_snv_paired_rederived.json
"""

from __future__ import annotations

import json
import os
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

RUN = pathlib.Path(os.environ["AGA_OUT_ROOT"])
RESULTS = RUN.parent
PREVIOUS = RESULTS / "p6f-indel-rescue-20260914T135052Z" / "tables" / "indel_rescue_effects.tsv"
CHANNELS = ["rna", "atac", "dnase", "h3k27ac", "splice"]
STRATA = ["all", "dbsnp_oriented", "dbsnp_unique_record", "reference_resolved"]
N_BOOT = 2000


def load(path: pathlib.Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for ch in CHANNELS:
        for part in ("log2", "ref", "alt"):
            col = f"{ch}_{part}"
            if col in df:
                df[col] = pd.to_numeric(df[col].replace("", np.nan), errors="coerce")
    return df


def reference_resolved(results: pathlib.Path) -> set:
    tables = sorted(p for p in results.glob("p0-indel-uid-repair-*/tables/deferred_indel_uid_repair.tsv")
                    if not (p.parents[1] / "SUPERSEDED.txt").exists())
    rep = pd.read_csv(tables[-1], sep="\t", dtype=str, keep_default_na=False)
    keep = rep[(rep["repair_state"] == "recovered") & (rep["orientation_ambiguous"] == "False")]
    return set(keep["source_variant_id"])


def paired(ind: pd.DataFrame, ctl: pd.DataFrame, ch: str) -> pd.DataFrame:
    c = ctl.assign(v=ctl[f"{ch}_log2"].abs()).dropna(subset=["v"])
    med = c.groupby("control_for")["v"].median().rename("ctl_median")
    i = ind.assign(v=ind[f"{ch}_log2"].abs()).dropna(subset=["v"])
    out = i.merge(med, left_on="source_variant_id", right_index=True, how="inner")
    pos = out["pos_hg38"].astype(int)
    return out.assign(diff=out["v"] - out["ctl_median"],
                      block=out["chrom"] + ":" + (pos // 1_000_000).astype(str))


def block_bootstrap(d: pd.DataFrame, seed: int) -> tuple:
    """Blocks resampled with replacement, one value drawn within each chosen block, median of the draw."""
    groups = [g.to_numpy() for _, g in d.groupby("block")["diff"]]
    rng = np.random.default_rng(seed)
    meds = np.empty(N_BOOT)
    for b in range(N_BOOT):
        pick = rng.integers(0, len(groups), size=len(groups))
        meds[b] = np.median([groups[j][rng.integers(0, groups[j].size)] for j in pick])
    return float(np.quantile(meds, 0.025)), float(np.quantile(meds, 0.975))


def main() -> None:
    eff = load(RUN / "tables" / "indel_rescue_effects.tsv")
    reported = json.loads((RUN / "tables" / "indel_vs_snv_paired.json").read_text())
    ind_all = eff[(eff["arm"] == "indel") & (eff["state"] == "scored")]
    ctl = eff[eff["arm"] == "snv_control"]
    ref_ok = reference_resolved(RESULTS)
    subsets = {
        "all": ind_all,
        "dbsnp_oriented": ind_all[ind_all["orientation_source"] == "dbsnp"],
        "dbsnp_unique_record": ind_all[(ind_all["orientation_source"] == "dbsnp")
                                       & (ind_all["dbsnp_rule"] == "unique_record")],
        "reference_resolved": ind_all[ind_all["source_variant_id"].isin(ref_ok)],
    }
    site = ind_all["chrom"] + ":" + ind_all["pos_hg38"] + ":" + [
        "/".join(sorted((a, b))) for a, b in zip(ind_all["ref"], ind_all["alt"])]
    out = {"n_scored_indels": int(len(ind_all)), "n_distinct_variants": int(site.nunique()),
           "n_controls": int(len(ctl)), "mismatches": [], "strata": {}}
    for name in STRATA:
        ind = subsets[name]
        rows, ps = {}, []
        for ch in CHANNELS:
            d = paired(ind, ctl, ch)
            bm = d.groupby("block")["diff"].median()
            n_nonzero, k = int((bm != 0).sum()), int((bm > 0).sum())
            p = stats.binomtest(k, n_nonzero, 0.5).pvalue if n_nonzero else float("nan")
            lo, hi = block_bootstrap(d, seed=7) if len(d) else (None, None)
            rows[ch] = {"n_paired_indels": int(len(d)), "n_blocks": int(bm.size),
                        "median_paired_diff_log2": float(d["diff"].median()) if len(d) else None,
                        "blocks_positive": f"{k}/{n_nonzero}", "sign_test_p": float(p), "lo": lo, "hi": hi}
            ps.append(p)
        finite = [x for x in ps if x == x]
        qs = iter(stats.false_discovery_control(finite, method="bh")) if finite else iter([])
        for ch, p in zip(CHANNELS, ps):
            rows[ch]["bh_q"] = float(next(qs)) if p == p else float("nan")
        out["strata"][name] = {"n_indels_scored": int(len(ind)), **rows}

        rep = reported[name]
        if rep["n_indels_scored"] != len(ind):
            out["mismatches"].append(f"{name}: n {rep['n_indels_scored']} vs {len(ind)}")
        for ch in CHANNELS:
            a, b = rep[ch], rows[ch]
            for key in ("n_paired_indels", "n_blocks", "blocks_positive"):
                if a[key] != b[key]:
                    out["mismatches"].append(f"{name}/{ch}/{key}: {a[key]} vs {b[key]}")
            for key in ("median_paired_diff_log2", "sign_test_p", "bh_q"):
                x, y = a[key], b[key]
                if not ((x is None and y is None) or (x != x and y != y) or abs(x - y) <= 1e-12 * max(1.0, abs(x))):
                    out["mismatches"].append(f"{name}/{ch}/{key}: {x} vs {y}")
            if a["lo"] is not None and b["lo"] is not None:
                width = max(a["hi"] - a["lo"], 1e-12)
                shift = max(abs(a["lo"] - b["lo"]), abs(a["hi"] - b["hi"])) / width
                b["ci_shift_over_width"] = shift
                if shift > 0.25:
                    out["mismatches"].append(f"{name}/{ch}/bootstrap interval moved {shift:.2f} of its width")

    if PREVIOUS.exists():
        prev = load(PREVIOUS)
        cols = [f"{ch}_{part}" for ch in CHANNELS for part in ("ref", "alt", "log2")]
        det = {}
        for arm, key in (("snv_control", ["control_for", "source_variant_id"]),
                         ("indel", ["source_variant_id", "ref", "alt"])):
            a = eff[(eff["arm"] == arm) & (eff["state"] == "scored")][key + cols]
            b = prev[(prev["arm"] == arm) & (prev["state"] == "scored")][key + cols]
            m = a.merge(b, on=key, suffixes=("_new", "_old"))
            # both missing counts as identical; one missing is a difference, never silently equal
            same = np.column_stack([((m[f"{c}_new"] == m[f"{c}_old"])
                                     | (m[f"{c}_new"].isna() & m[f"{c}_old"].isna())).to_numpy() for c in cols]) \
                if len(m) else np.zeros((0, 1), bool)
            diffs = np.column_stack([(m[f"{c}_new"] - m[f"{c}_old"]).abs().to_numpy() for c in cols]) \
                if len(m) else np.zeros((0, 1))
            det[arm] = {"n_matched": int(len(m)), "n_identical_on_all_channels": int(same.all(axis=1).sum()),
                        "max_abs_difference": float(np.nanmax(diffs)) if diffs.size and np.isfinite(diffs).any() else None}
        out["same_alleles_in_previous_run"] = det

    (RUN / "tables" / "indel_vs_snv_paired_rederived.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({"mismatches": out["mismatches"], "determinism": out.get("same_alleles_in_previous_run"),
                      "n_scored": out["n_scored_indels"], "n_distinct": out["n_distinct_variants"]}, indent=1))
    sys.exit(1 if out["mismatches"] else 0)


if __name__ == "__main__":
    main()
