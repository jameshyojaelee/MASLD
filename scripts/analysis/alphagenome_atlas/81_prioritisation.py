#!/usr/bin/env python3
"""Step 81 (P6e): does the validated accessibility channel prioritise the functional variant?

The Atlas accessibility channel is validated on three assays, so the translational question is whether it
IDENTIFIES the right variant. Ground truth is measured allelic imbalance in GSE281367 snATAC: within a
credible set, is the allelically ACTIVE variant ranked above the inactive ones?

Every variant compared is peak-resident and read-covered, because only those can be allelically tested. That
is what makes the comparison fair: ranking an allelic site against the whole credible set would be won by
peak membership alone, since the Atlas scores peak-resident variants higher and only peak-resident variants
can be tested. Both classes here were measured in the same donors under the same filters.

Posterior weight is the competing ranker, and the honest baseline: fine-mapping says which variant is
statistically credible, not which is chromatin-active.

Prespecification: 81_prioritisation_prespec.json.
Outputs (tables/): prioritisation_sets.tsv, prioritisation.json
"""

from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import json
import shutil

import numpy as np
import pandas as pd
from scipy.stats import rankdata

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
PRESPEC = la.SCRIPT_DIR / "81_prioritisation_prespec.json"
ALLELIC = la.PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/allelic_sites.tsv"
TRACK0 = la.track0_root() / "tables"
ACTIVE_THRESHOLD = 0.5
MIN_TESTED = 5
DRAWS = 2000
SEED = 20260911
BLOCK_BP = 1_000_000


def mean_active_rank_percentile(score: np.ndarray, active: np.ndarray) -> float:
    """Mean within-set rank percentile of the ACTIVE variants. 0.5 is chance.

    Ties get the average rank, so a set of identical scores puts the active variant at 0.5 rather than the
    top. NaN scores are dropped rather than ranked last, which would otherwise reward the active variant
    whenever a comparator happened to be unscored.
    """
    s = np.asarray(score, float); a = np.asarray(active, bool)
    ok = ~np.isnan(s)
    s, a = s[ok], a[ok]
    if a.sum() == 0 or (~a).sum() == 0:
        return float("nan")
    r = rankdata(s)                       # average ranks for ties
    pct = (r - 1) / (len(r) - 1) if len(r) > 1 else np.array([0.5])
    return float(np.mean(pct[a]))


def permuted_active(sets: np.ndarray, active: np.ndarray, draws: int, seed: int) -> list:
    """Active/inactive labels reshuffled WITHIN each set, preserving how many are active per set."""
    rng = np.random.default_rng(seed)
    sets = np.asarray(sets); active = np.asarray(active, bool)
    idx = {s: np.flatnonzero(sets == s) for s in np.unique(sets)}
    out = []
    for _ in range(draws):
        lab = active.copy()
        for members in idx.values():
            lab[members] = rng.permutation(active[members])
        out.append(lab)
    return out


def combined_rank(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Product of the two rankers' within-set rank percentiles."""
    def pct(x):
        x = np.asarray(x, float)
        r = rankdata(np.where(np.isnan(x), -np.inf, x))
        return (r - 1) / (len(r) - 1) if len(r) > 1 else np.array([0.5] * len(x))
    return pct(a) * pct(b)


def block_bootstrap(values: dict, blocks: dict, draws: int, seed: int) -> dict:
    """One set per 1-Mb block per draw: sets at one locus are not independent."""
    rng = np.random.default_rng(seed)
    by_block = collections.defaultdict(list)
    for s, v in values.items():
        if v == v:
            by_block[blocks[s]].append(v)
    uniq = sorted(by_block)
    if len(uniq) < 2:
        return {"mean": float("nan"), "ci95": [float("nan"), float("nan")], "n_blocks": len(uniq)}
    vals = []
    for _ in range(draws):
        picked = rng.choice(len(uniq), size=len(uniq), replace=True)
        vals.append(float(np.mean([rng.choice(by_block[uniq[i]]) for i in picked])))
    v = np.asarray(vals)
    return {"mean": float(v.mean()), "ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
            "n_blocks": len(uniq)}


def main() -> None:
    spec = json.loads(PRESPEC.read_text())
    digest = hashlib.sha256(PRESPEC.read_bytes()).hexdigest()
    shutil.copy(PRESPEC, ROOT / "81_prioritisation_prespec.json")

    al = pd.read_csv(ALLELIC, sep="\t")
    meas = dict(zip(al.uid, al.mean_log2_alt_over_ref.abs()))
    atac = dict(zip(al.uid, al.atac_quantile.abs()))

    xw = {}
    with gzip.open(TRACK0 / "variant_crosswalk.tsv.gz", "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r.get("variant_uid"):
                xw[r["source_variant_id"]] = r["variant_uid"]
    mem = collections.defaultdict(list)
    with gzip.open(TRACK0 / "signal_variant_weights.tsv.gz", "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            u = xw.get(r["source_variant_id"])
            if u in meas:
                mem[r["signal_uid"]].append((u, float(r["weight"])))

    rows, per_set = [], {}
    for sig, vs in mem.items():
        vs = list({u: w for u, w in vs}.items())
        if len(vs) < MIN_TESTED:
            continue
        uids = [u for u, _ in vs]
        act = np.array([meas[u] > ACTIVE_THRESHOLD for u in uids])
        if act.sum() == 0 or (~act).sum() == 0:
            continue
        a_score = np.array([atac.get(u, np.nan) for u in uids])
        p_score = np.array([w for _, w in vs], float)
        c_score = combined_rank(a_score, p_score)
        chrom, pos = uids[0].split(":")[0], int(uids[0].split(":")[1])
        per_set[sig] = {
            "atlas": mean_active_rank_percentile(a_score, act),
            "posterior": mean_active_rank_percentile(p_score, act),
            "combined": mean_active_rank_percentile(c_score, act),
            "block": f"{chrom}~{pos // BLOCK_BP}", "n_tested": len(uids), "n_active": int(act.sum()),
            "uids": uids, "active": act, "a": a_score, "p": p_score, "c": c_score}
        rows.append({"signal_uid": sig, "n_tested": len(uids), "n_active": int(act.sum()),
                     "block": per_set[sig]["block"], "rank_atlas": per_set[sig]["atlas"],
                     "rank_posterior": per_set[sig]["posterior"], "rank_combined": per_set[sig]["combined"]})

    summary = {"prespec_sha256": digest, "n_eligible_sets": len(per_set),
               "active_threshold_abs_log2": ACTIVE_THRESHOLD, "min_tested_per_set": MIN_TESTED}
    if len(per_set) < 30:
        summary["verdict"] = f"cannot conclude: only {len(per_set)} eligible sets (minimum 30)"
        json.dump(summary, (TABLES / "prioritisation.json").open("w"), indent=1, default=float)
        la.log(f"P6e: {len(per_set)} eligible sets, cannot conclude")
        return

    la.write_tsv_once(TABLES / "prioritisation_sets.tsv", rows, list(rows[0].keys()))
    blocks = {s: v["block"] for s, v in per_set.items()}
    sets_arr = np.concatenate([[s] * len(per_set[s]["uids"]) for s in per_set])
    act_arr = np.concatenate([per_set[s]["active"] for s in per_set])
    score_arrs = {k: np.concatenate([per_set[s][k[0]] for s in per_set]) for k in ("atlas", "posterior", "combined")}
    key = {"atlas": "a", "posterior": "p", "combined": "c"}
    score_arrs = {k: np.concatenate([per_set[s][key[k]] for s in per_set]) for k in key}

    perms = permuted_active(sets_arr, act_arr, DRAWS, SEED)
    for ranker in ("atlas", "posterior", "combined"):
        obs = {s: per_set[s][ranker] for s in per_set}
        pooled = float(np.nanmean(list(obs.values())))
        boot = block_bootstrap(obs, blocks, DRAWS, SEED)
        null = []
        sc = score_arrs[ranker]
        for lab in perms:
            vals = []
            for s in per_set:
                m = sets_arr == s
                v = mean_active_rank_percentile(sc[m], lab[m])
                if v == v:
                    vals.append(v)
            null.append(float(np.mean(vals)) if vals else np.nan)
        null = np.asarray(null, float)
        summary[ranker] = {"mean_rank_percentile": pooled, "block_bootstrap": boot,
                           "null_mean": float(np.nanmean(null)),
                           "p_empirical": float((np.sum(null >= pooled) + 1) / (np.sum(~np.isnan(null)) + 1)),
                           "beats_chance": bool(boot["ci95"][0] > 0.5)}
    summary["atlas_minus_posterior"] = summary["atlas"]["mean_rank_percentile"] - summary["posterior"]["mean_rank_percentile"]
    summary["combined_minus_atlas"] = summary["combined"]["mean_rank_percentile"] - summary["atlas"]["mean_rank_percentile"]
    diag = []
    if summary["atlas"]["block_bootstrap"]["n_blocks"] < 15:
        diag.append(f"only {summary['atlas']['block_bootstrap']['n_blocks']} independent blocks")
    if abs(summary["atlas"]["null_mean"] - 0.5) > 0.01:
        diag.append(f"permutation null centres at {summary['atlas']['null_mean']:.3f}, not 0.5")
    summary["diagnostics_triggered"] = diag
    summary["predictions"] = spec["predictions_before_results"]
    summary["claim_boundary"] = spec["claim_boundary"]
    summary["verdict"] = ("cannot conclude" if diag else
                          "the accessibility channel prioritises above chance" if summary["atlas"]["beats_chance"]
                          else "no prioritisation above chance")
    json.dump(summary, (TABLES / "prioritisation.json").open("w"), indent=1, default=float)
    la.log(f"P6e: {len(per_set)} sets; atlas {summary['atlas']['mean_rank_percentile']:.3f} "
           f"CI{[round(x,3) for x in summary['atlas']['block_bootstrap']['ci95']]}; "
           f"posterior {summary['posterior']['mean_rank_percentile']:.3f}; delta {summary['atlas_minus_posterior']:+.3f}")


if __name__ == "__main__":
    main()
