#!/usr/bin/env python3
"""Step 80 (P1-E1b): do inherited liver-trait risk variants converge on disruption of specific TFs?

The accessibility channel is validated on three independent assays, so the Atlas is used here as an
INSTRUMENT rather than benchmarked. With 21 adult-liver and 539 HepG2 TF ChIP tracks covering ~501 factors,
the question is which factors' predicted binding is preferentially disrupted by the credible variants of
liver-trait signals.

Two design choices carry the analysis:

  * **Within-signal pairing.** Credible variants (posterior weight >= 0.05) are compared with non-credible
    variants of the SAME signal (weight < 0.005). Same locus, same LD block, same GC and peak context. No
    matched draw has to be constructed because the control is the rest of the credible set.
  * **Factor-specificity deviation.** A variant in a strong regulatory element perturbs many factors at once.
    Subtracting each variant's own mean across factors separates WHICH factor is preferentially hit from HOW
    MUCH the variant does overall, which is the question actually being asked.

Prespecification: 80_tf_convergence_prespec.json (filters, seeds, family, written predictions).
Outputs (tables/): tf_convergence_factors.tsv, tf_convergence_signals.tsv.gz, tf_convergence.json
"""

from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import json
import math
import shutil
import sys

import numpy as np
import pandas as pd

import atlas_archive as aa
import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
PRESPEC = la.SCRIPT_DIR / "80_tf_convergence_prespec.json"
W_TEST = 0.05
W_CTRL = 0.005
MIN_TEST = 1
MIN_CTRL = 5
MIN_SIGNALS_PER_FACTOR = 20
DRAWS = 2000
SEED = 20260911
BH_Q = 0.05
ADULT_LIVER = {"UBERON:0002107", "UBERON:0001114", "UBERON:0001115"}
HEPG2 = {"EFO:0001187"}


def factor_deviation(v: np.ndarray) -> np.ndarray:
    """Per-variant factor-specificity deviation: each factor's |quantile| minus the variant's own mean.

    Removes the variant's overall sequence impact, so a uniformly high-impact variant contributes zero to
    every factor and only PREFERENTIAL disruption survives. NaN factors stay NaN rather than counting as 0.
    """
    a = np.asarray(v, float)
    m = np.nanmean(a) if np.any(~np.isnan(a)) else np.nan
    return a - m


def signal_contrast(dev: np.ndarray, is_test: np.ndarray, weight: np.ndarray) -> float:
    """Posterior-weighted mean deviation over test variants minus the unweighted mean over controls."""
    dev = np.asarray(dev, float); is_test = np.asarray(is_test, bool); w = np.asarray(weight, float)
    t, c = is_test & ~np.isnan(dev), (~is_test) & ~np.isnan(dev)
    if not t.any() or not c.any():
        return float("nan")
    wt = w[t]
    wt = wt / wt.sum() if wt.sum() > 0 else np.full(t.sum(), 1.0 / t.sum())
    return float(np.sum(dev[t] * wt) - np.mean(dev[c]))


def permuted_labels(signal: np.ndarray, is_test: np.ndarray, draws: int, seed: int) -> list:
    """Test/control labels reshuffled WITHIN each signal; across-signal shuffling would break the pairing."""
    rng = np.random.default_rng(seed)
    signal = np.asarray(signal); is_test = np.asarray(is_test, bool)
    idx = {s: np.flatnonzero(signal == s) for s in np.unique(signal)}
    out = []
    for _ in range(draws):
        lab = is_test.copy()
        for members in idx.values():
            lab[members] = rng.permutation(is_test[members])
        out.append(lab)
    return out


def stable_seed(name: str) -> int:
    """Deterministic per-signal seed from a sha256 digest.

    Python randomises str.__hash__ per process, so seeding a subsample with hash(signal_uid) draws a
    different control set on every run - two runs of this job reported 38,523 and 38,456 variants before this
    was fixed. A digest is stable across processes and machines.
    """
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def signal_groups(signal: np.ndarray) -> list:
    """Row indices per signal, computed ONCE. Recomputing a mask per factor per draw is what made the first
    implementation unusable (about 1.4e10 element comparisons per permutation draw)."""
    signal = np.asarray(signal)
    order = np.argsort(signal, kind="mergesort")
    s_sorted = signal[order]
    bounds = np.flatnonzero(np.r_[True, s_sorted[1:] != s_sorted[:-1], True])
    return [order[bounds[i]:bounds[i + 1]] for i in range(len(bounds) - 1)]


def per_factor_fast(dev: np.ndarray, groups: list, labels: np.ndarray, weight: np.ndarray,
                    min_signals: int) -> np.ndarray:
    """Mean signal-level contrast per factor, signal-major and vectorised across factors.

    Numerically identical to looping signal_contrast over every factor x signal pair, which the fixture
    tests assert; it is simply ordered so each signal's slice is touched once per draw instead of once per
    factor per draw.
    """
    n_fac = dev.shape[1]
    total = np.zeros(n_fac)
    count = np.zeros(n_fac, dtype=int)
    for idx in groups:
        lab = labels[idx]
        if not lab.any() or lab.all():
            continue
        d = dev[idx]
        t_rows, c_rows = d[lab], d[~lab]
        w = weight[idx][lab].astype(float)
        ok_t = ~np.isnan(t_rows)
        ok_c = ~np.isnan(c_rows)
        wt = np.where(ok_t, w[:, None], 0.0)
        wsum = wt.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            # a factor with no posterior mass on its scored test rows falls back to an unweighted mean
            flat = np.where(ok_t, 1.0, 0.0)
            wt = np.where(wsum[None, :] > 0, wt, flat)
            wsum = np.where(wsum > 0, wsum, flat.sum(axis=0))
            t_mean = np.nansum(np.where(ok_t, t_rows, 0.0) * wt, axis=0) / wsum
            c_mean = np.nansum(np.where(ok_c, c_rows, 0.0), axis=0) / ok_c.sum(axis=0)
        good = (ok_t.any(axis=0)) & (ok_c.any(axis=0)) & ~np.isnan(t_mean) & ~np.isnan(c_mean)
        total[good] += (t_mean - c_mean)[good]
        count[good] += 1
    out = np.full(n_fac, np.nan)
    enough = count >= min_signals
    out[enough] = total[enough] / count[enough]
    return out


def permuted_label_weight(signal: np.ndarray, is_test: np.ndarray, weight: np.ndarray,
                          draws: int, seed: int) -> list:
    """Permute the (label, weight) PAIR within each signal.

    The statistic weights the test side by posterior and leaves the control side unweighted, so shuffling
    labels alone is not exchangeable: a permuted test set would be built from low-weight rows and the null
    drifts off zero. Moving each row's weight with its label keeps the permuted data structurally like the
    observed data.
    """
    rng = np.random.default_rng(seed)
    signal = np.asarray(signal); is_test = np.asarray(is_test, bool); weight = np.asarray(weight, float)
    idx = {s: np.flatnonzero(signal == s) for s in np.unique(signal)}
    out = []
    for _ in range(draws):
        lab = is_test.copy(); w = weight.copy()
        for members in idx.values():
            perm = rng.permutation(len(members))
            lab[members] = is_test[members][perm]
            w[members] = weight[members][perm]
        out.append((lab, w))
    return out


def rank_enrichment(effects: np.ndarray, names: np.ndarray, wanted: set, draws: int, seed: int) -> dict:
    """Do the named factors rank above chance by effect size?

    One test instead of 502. The per-factor family cannot carry BH because the deviations are compositional -
    they sum to zero within each variant - so a real shift in profile shape necessarily lights up many
    factors at once. Rank 0 is the strongest effect.
    """
    e = np.asarray(effects, float); n = np.asarray(names)
    ok = ~np.isnan(e)
    e, n = e[ok], n[ok]
    order = np.argsort(-e)
    rank = {str(nm): i for i, nm in enumerate(n[order])}
    hit = sorted({str(x) for x in n} & set(wanted))
    if not hit:
        return {"n_in_set": 0, "mean_rank": float("nan"), "p_empirical": float("nan")}
    obs = float(np.mean([rank[h] for h in hit]))
    rng = np.random.default_rng(seed)
    null = np.array([float(np.mean(rng.choice(len(e), size=len(hit), replace=False))) for _ in range(draws)])
    return {"n_in_set": len(hit), "members": hit, "mean_rank": obs,
            "expected_mean_rank": (len(e) - 1) / 2.0, "n_tested": int(len(e)),
            "null_mean": float(null.mean()),
            "p_empirical": float((np.sum(null <= obs) + 1) / (draws + 1)),
            "ranks": {h: rank[h] for h in hit}}


def bh_reject(p: np.ndarray, q: float) -> np.ndarray:
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    out = np.zeros(p.shape, bool)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx], kind="mergesort")]
    m = order.size
    passing = np.flatnonzero(p[order] <= q * (np.arange(1, m + 1) / m))
    if passing.size:
        out[order[: passing[-1] + 1]] = True
    return out


def load_factor_panel() -> tuple[dict, dict]:
    """track index -> factor name, restricted to adult-liver and HepG2 biosamples."""
    md = pd.read_csv(RAW / "scorer_metadata" / "CHIP_TF.track_metadata.tsv", sep="\t")
    keep = md["ontology_curie"].isin(ADULT_LIVER | HEPG2)
    fac = md.loc[keep, "transcription_factor"].astype(str).str.strip()
    src = np.where(md.loc[keep, "ontology_curie"].isin(ADULT_LIVER), "adult_liver", "HepG2")
    return {i: f for i, f in zip(np.flatnonzero(keep.to_numpy()), fac)}, \
           {i: s for i, s in zip(np.flatnonzero(keep.to_numpy()), src)}


def main() -> None:
    spec = json.loads(PRESPEC.read_text())
    digest = hashlib.sha256(PRESPEC.read_bytes()).hexdigest()
    shutil.copy(PRESPEC, ROOT / "80_tf_convergence_prespec.json")

    weights = collections.defaultdict(dict)
    with gzip.open(TABLES / "signal_variant_weights.tsv.gz", "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            weights[r["signal_uid"]][r["source_variant_id"]] = float(r["weight"])
    xw = {}
    with gzip.open(TABLES / "variant_crosswalk.tsv.gz", "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r.get("variant_uid"):
                xw[r["source_variant_id"]] = r["variant_uid"]

    # membership: uid -> list of (signal, weight, is_test)
    member = collections.defaultdict(list)
    contributing = set()
    for sig, vw in weights.items():
        test = [(v, w) for v, w in vw.items() if w >= W_TEST and v in xw]
        ctrl = [(v, w) for v, w in vw.items() if w < W_CTRL and v in xw]
        if len(test) < MIN_TEST or len(ctrl) < MIN_CTRL:
            continue
        contributing.add(sig)
        rng = np.random.default_rng(stable_seed(sig))
        if len(ctrl) > 50:                      # cap controls per signal so one huge set cannot dominate
            ctrl = [ctrl[i] for i in rng.choice(len(ctrl), 50, replace=False)]
        for v, w in test:
            member[xw[v]].append((sig, w, True))
        for v, w in ctrl:
            member[xw[v]].append((sig, w, False))
    la.log(f"TF convergence: {len(contributing)} contributing signals, {len(member)} distinct variants")

    fac_of, src_of = load_factor_panel()
    cols = sorted(set(fac_of.values()))
    col_idx = {f: i for i, f in enumerate(cols)}
    tracks_by_factor = collections.defaultdict(list)
    for ti, f in fac_of.items():
        tracks_by_factor[f].append(ti)

    rows = []
    seen = set()
    for chunk in sorted((RAW / "atlas_direct").glob("chunk_*")) + sorted((RAW / "atlas_enzyme").glob("chunk_*")):
        res = aa.load_archive(chunk)
        a = res.get("CHIP_TF")
        if a is None or a.shape[0] == 0 or "variant" not in a.obs:
            continue
        obs = a.obs.reset_index(drop=True)
        uids = [aa.variant_uid_from_str(str(x)) for x in obs["variant"]]
        want = [i for i, u in enumerate(uids) if u in member and u not in seen]
        if not want:
            continue
        q = np.abs(np.asarray(a.layers["quantiles"], np.float32))
        for i in want:
            u = uids[i]
            seen.add(u)
            v = np.full(len(cols), np.nan)
            for f, tis in tracks_by_factor.items():
                vals = q[i, tis]
                if np.any(~np.isnan(vals)):
                    v[col_idx[f]] = float(np.nanmean(vals))
            dev = factor_deviation(v)
            for sig, w, is_t in member[u]:
                rows.append((sig, u, w, is_t, dev))
        if len(seen) % 5000 < 50:
            la.log(f"  {len(seen)} variants read")
    la.log(f"TF convergence: {len(seen)} variants had CHIP_TF rows; {len(rows)} signal-variant pairs")

    if len(rows) < 100:
        json.dump({"verdict": "cannot conclude: too few variants with CHIP_TF rows", "n_rows": len(rows),
                   "prespec_sha256": digest}, (TABLES / "tf_convergence.json").open("w"), indent=1)
        return

    sig_arr = np.array([r[0] for r in rows])
    w_arr = np.array([r[2] for r in rows], float)
    t_arr = np.array([r[3] for r in rows], bool)
    dev_mat = np.vstack([r[4] for r in rows])

    groups = signal_groups(sig_arr)
    np.savez_compressed(TABLES / "tf_convergence_matrix.npz", dev=dev_mat.astype(np.float32),
                        sig=sig_arr, w=w_arr, t=t_arr, cols=np.array(cols))
    la.log("deviation matrix cached; a rerun of the scoring need not re-read the archive")

    def per_factor(labels: np.ndarray) -> np.ndarray:
        return per_factor_fast(dev_mat, groups, labels, w_arr, MIN_SIGNALS_PER_FACTOR)

    obs_stat = per_factor(t_arr)
    null = np.vstack([per_factor_fast(dev_mat, groups, lab, w, MIN_SIGNALS_PER_FACTOR)
                      for lab, w in permuted_label_weight(sig_arr, t_arr, w_arr, DRAWS, SEED)])
    pvals = np.array([((np.sum(np.abs(null[:, j]) >= abs(obs_stat[j])) + 1) / (null.shape[0] + 1))
                      if obs_stat[j] == obs_stat[j] else np.nan for j in range(len(cols))])
    # Per-factor BH is NOT a valid family: the deviations sum to zero within every variant, so any real
    # shift in profile shape pushes some factors up and an equal mass down and many cross together (247 of
    # 502 did, on a correctly centred null). The per-factor p-values are kept as descriptive output and the
    # inferential question is the prespecified rank enrichment below - one test, with its own negative control.
    rej = bh_reject(pvals, BH_Q)

    contrib = np.zeros(len(cols), dtype=int)
    for idx in groups:
        lab = t_arr[idx]
        if not lab.any() or lab.all():
            continue
        d = dev_mat[idx]
        contrib += ((~np.isnan(d[lab])).any(axis=0) & (~np.isnan(d[~lab])).any(axis=0)).astype(int)
    n_sig = {f: int(contrib[i]) for i, f in enumerate(cols)}
    out_rows = [{"factor": f, "n_signals": n_sig[f], "effect": obs_stat[i],
                 "null_mean": float(np.nanmean(null[:, i])), "p_empirical": pvals[i], "bh_reject": bool(rej[i]),
                 "track_source": "adult_liver" if all(src_of[t] == "adult_liver" for t in tracks_by_factor[f])
                 else ("HepG2" if all(src_of[t] == "HepG2" for t in tracks_by_factor[f]) else "both"),
                 "n_tracks": len(tracks_by_factor[f])}
                for i, f in enumerate(cols)]
    out_rows.sort(key=lambda r: (-(r["effect"] if r["effect"] == r["effect"] else -9), r["factor"]))
    la.write_tsv_once(TABLES / "tf_convergence_factors.tsv", out_rows, list(out_rows[0].keys()))

    tested = [r for r in out_rows if r["effect"] == r["effect"]]
    hits = [r for r in tested if r["bh_reject"]]
    HEP = {"HNF4A", "HNF1A", "HNF4G", "CEBPA", "CEBPB", "FOXA1", "FOXA2", "NR2F2", "ONECUT1", "HNF1B"}
    GENERAL = {"POLR2A", "TAF1", "CTCF", "RAD21", "YY1", "EP300"}
    top10 = [r["factor"] for r in tested[:10]]
    name_arr = np.array(cols)
    hep_enr = rank_enrichment(obs_stat, name_arr, HEP, draws=20000, seed=SEED)
    gen_enr = rank_enrichment(obs_stat, name_arr, GENERAL, draws=20000, seed=SEED)
    summary = {
        "prespec_sha256": digest, "n_contributing_signals": len(contributing),
        "n_variants_with_chip_tf": len(seen), "n_factors_tested": len(tested), "n_factors_bh": len(hits),
        "bh_q": BH_Q, "draws": DRAWS, "seed": SEED,
        "null_centre_max_abs": float(np.nanmax(np.abs(np.nanmean(null, axis=0)))),
        "top10_by_effect": top10,
        "hepatocyte_regulators_among_hits": sorted({r["factor"] for r in hits} & HEP),
        "hepatocyte_share_of_hits": (len({r["factor"] for r in hits} & HEP) / len(hits)) if hits else None,
        "hepatocyte_share_of_tested": len({r["factor"] for r in tested} & HEP) / max(len(tested), 1),
        "general_machinery_in_top10": sorted(set(top10) & GENERAL),
        "per_factor_bh_is_not_a_valid_family": ("the factor deviations sum to zero within every variant, so the "
                                                "502 factors are compositional; per-factor p-values are descriptive only"),
        "effects_sum": float(np.nansum(obs_stat)),
        "n_factors_crossing_bh_descriptive_only": len(hits),
        "rank_enrichment_hepatocyte_regulators": hep_enr,
        "rank_enrichment_general_machinery_CONTROL": gen_enr,
        "control_beats_hypothesis": bool(gen_enr["mean_rank"] == gen_enr["mean_rank"] and
                                         hep_enr["mean_rank"] == hep_enr["mean_rank"] and
                                         gen_enr["mean_rank"] < hep_enr["mean_rank"]),
        "predictions": spec["predictions_before_results"],
        "claim_boundary": spec["claim_boundary"],
    }
    diag = []
    if len(contributing) < 50:
        diag.append("fewer than 50 contributing signals")
    if summary["null_centre_max_abs"] > 0.002:
        diag.append(f"permutation null does not centre at zero (max |mean| {summary['null_centre_max_abs']:.4f})")
    summary["diagnostics_triggered"] = diag
    if diag:
        summary["verdict"] = "cannot conclude"
    elif summary["control_beats_hypothesis"]:
        summary["verdict"] = ("NO MASLD-specific TF convergence: the prespecified negative control (general "
                              "machinery) ranks higher than the hepatocyte regulators it was meant to be "
                              "distinguished from, so the enrichment reflects credible variants sitting in more "
                              "active regions where ubiquitously bound factors move most")
    elif hep_enr["p_empirical"] < 0.05:
        summary["verdict"] = "hepatocyte regulators rank above chance and above the general-machinery control"
    else:
        summary["verdict"] = "no convergence detected"
    json.dump(summary, (TABLES / "tf_convergence.json").open("w"), indent=1, default=float)
    la.log(f"TF convergence: {len(tested)} factors tested, {len(hits)} survive BH q<{BH_Q}; top {top10[:5]}")


if __name__ == "__main__":
    main()
