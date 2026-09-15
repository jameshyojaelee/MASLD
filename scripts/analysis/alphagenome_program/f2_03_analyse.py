#!/usr/bin/env python3
"""F2/F3 step 3: verdicts from the scored rows. No network use.

Reads raw/scored_rows.jsonl (written by f2_02_score.py) and emits:
  tables/f2_gnmt_effects.tsv        the four-arm effects, additive expectation and residual per channel
  tables/f2_gnmt_verdict.json       F2.1 against the 100-draw chr6 matched null
  tables/f3_observed_effects.tsv    the 200 rescored pairs with empirical p and BH flags per channel
  tables/f3_null_effects.tsv        the 1,000 corrected null draws
  tables/f3_summary.json            per-channel verdicts beside the existing P5 deposit's
"""

from __future__ import annotations

import csv
import json
import math
import os
import pathlib

import numpy as np

import f2_recipe as R

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
P5 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables"
CHANNELS = R.CHANNELS


def load_rows() -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {"gnmt": [], "gnmtnull": [], "obs": [], "null": []}
    seen = set()
    with open(OUT / "raw" / "scored_rows.jsonl") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r["tag"] in seen:
                continue
            seen.add(r["tag"])
            t = r["tag"]
            key = "gnmtnull" if t.startswith("gnmtnull") else ("gnmt" if t.startswith("gnmt") else ("null" if t.startswith("null") else "obs"))
            groups[key].append(r)
    return groups


def write_tsv(path: pathlib.Path, rows: list[dict]) -> None:
    if not rows:
        return
    cols = sorted({k for r in rows for k in r})
    lead = [c for c in ("tag", "chrom", "pos1", "pos2", "separation_bp", "gene", "signal_uid", "rna_readout") if c in cols]
    cols = lead + [c for c in cols if c not in lead]
    with open(path, "w") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def channel_verdict(ch: str, obs: list[dict], nul: list[dict], family_label: str) -> dict:
    o = np.array([r[f"{ch}_residual"] for r in obs], float)
    n = np.array([r[f"{ch}_residual"] for r in nul], float)
    v1 = np.array([r.get(f"{ch}_v1_log2", math.nan) for r in obs], float)
    v2 = np.array([r.get(f"{ch}_v2_log2", math.nan) for r in obs], float)
    nv1 = np.array([r.get(f"{ch}_v1_log2", math.nan) for r in nul], float)
    nv2 = np.array([r.get(f"{ch}_v2_log2", math.nan) for r in nul], float)
    thr = R.effect_threshold(nv1, nv2) if nv1.size else math.nan
    keep = R.effect_stratum(v1, v2, thr) if thr == thr else np.zeros(o.shape, bool)
    null_keep = R.effect_stratum(nv1, nv2, thr) if (thr == thr and nv1.size) else np.zeros(n.shape, bool)
    p = np.array([R.empirical_two_sided_p(x, n) for x in o], float)
    q95 = float(np.nanquantile(np.abs(n), 0.95)) if n.size else math.nan
    out = {"family": family_label, "n_pairs": int(o.size), "n_null_draws": int(n.size),
           "p_floor": float(1.0 / (int(np.sum(~np.isnan(n))) + 1)) if n.size else math.nan,
           "median_abs_residual": float(np.nanmedian(np.abs(o))) if o.size else math.nan,
           "null_median_abs_residual": float(np.nanmedian(np.abs(n))) if n.size else math.nan,
           "null_q95_abs_residual": q95,
           "n_pairs_exceeding_null_q95_UNCORRECTED": int(np.nansum(np.abs(o) > q95)) if q95 == q95 else None,
           "median_abs_single_variant_effect": float(np.nanmedian(np.abs(np.concatenate([v1, v2])))) if o.size else math.nan,
           "null_median_abs_single_variant_effect": float(np.nanmedian(np.abs(np.concatenate([nv1, nv2])))) if n.size else math.nan,
           "single_variant_effect_threshold_null_q95": thr,
           "null_pairs_meeting_active_stratum": int(null_keep.sum()),
           "bh_q": R.BH_Q, "min_p_empirical": float(np.nanmin(p)) if np.any(~np.isnan(p)) else math.nan}
    for label, mask in (("all_pairs", np.ones(o.shape, bool)), ("both_variants_active", keep)):
        gate = R.stratum_is_testable(int(mask.sum()), int(null_keep.sum()) if label == "both_variants_active" else int(n.size))
        if not gate["testable"]:
            out[label] = {"family_size_pairs": int(mask.sum()), **gate}
            continue
        pm = np.where(mask, p, np.nan)
        rej = R.bh_reject(pm, R.BH_Q)
        for i, r in enumerate(obs):
            r[f"{ch}_p_empirical"] = float(p[i])
            r[f"{ch}_both_active"] = bool(keep[i])
            r[f"{ch}_bh_reject_{'active' if label == 'both_variants_active' else 'all'}"] = bool(rej[i])
        out[label] = {"family_size_pairs": int(np.sum(~np.isnan(pm))), **gate,
                      "n_pairs_bh": int(rej.sum()),
                      "n_signals_bh": len(sorted({obs[i].get("signal_uid", obs[i]["tag"]) for i in np.flatnonzero(rej)})),
                      "signals_bh": sorted({obs[i].get("signal_uid", obs[i]["tag"]) for i in np.flatnonzero(rej)}),
                      "min_p_empirical": float(np.nanmin(pm)) if np.any(~np.isnan(pm)) else math.nan}
    return out


def p5_null_reference() -> dict:
    """The existing deposit's 300-draw null, recomputed here from its own table (a relayed number is a claim)."""
    rows = []
    with open(P5 / "haplotype_null_effects.tsv") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            rows.append(r)
    out = {"n_null": len(rows), "p_floor": 1.0 / (len(rows) + 1)}
    for ch in CHANNELS:
        a = np.array([float(r[f"{ch}_residual"]) if r[f"{ch}_residual"] not in ("", "nan") else math.nan for r in rows], float)
        s1 = np.array([float(r[f"{ch}_v1_log2"]) if r[f"{ch}_v1_log2"] not in ("", "nan") else math.nan for r in rows], float)
        s2 = np.array([float(r[f"{ch}_v2_log2"]) if r[f"{ch}_v2_log2"] not in ("", "nan") else math.nan for r in rows], float)
        out[ch] = {"median_abs_residual": float(np.nanmedian(np.abs(a))),
                   "q95_abs_residual": float(np.nanquantile(np.abs(a), 0.95)),
                   "median_abs_single_variant_effect": float(np.nanmedian(np.abs(np.concatenate([s1, s2]))))}
    out["rna_readout_counts"] = {}
    return out


def main() -> None:
    g = load_rows()
    # stage 1 (1,000 draws) and the amendment-03 extension (2,500) each keep their own tables
    tag = f"n{len(g['null'])}"
    obs, nul = g["obs"], g["null"]
    gnmt, gnull = g["gnmt"], g["gnmtnull"]
    summary: dict = {"n_obs": len(obs), "n_null": len(nul), "n_gnmt_pairs": len(gnmt), "n_gnmt_null": len(gnull),
                     "seed": 20260914}

    # ---------------------------------------------------------------- determinism
    prim = [r for r in gnmt if r.get("is_primary")]
    det = {}
    if prim and f"rna_ref_sum_repeat" in prim[0]:
        r = prim[0]
        det = {ch: {"ref_sum": r[f"{ch}_ref_sum"], "ref_sum_repeat": r[f"{ch}_ref_sum_repeat"],
                    "relative_difference": (abs(r[f"{ch}_ref_sum_repeat"] - r[f"{ch}_ref_sum"]) / abs(r[f"{ch}_ref_sum"])
                                            if r[f"{ch}_ref_sum"] else math.nan)} for ch in CHANNELS}
        det["verdict"] = ("deterministic at this window (every channel's repeated REF sum is identical)"
                          if all(det[ch]["relative_difference"] == 0 for ch in CHANNELS)
                          else "NOT bit-identical on a repeated REF call; residuals are differences of non-reproducible sums")
    summary["determinism_check_gnmt_primary_ref_arm"] = det

    # ---------------------------------------------------------------- F2
    f2 = {"pairs": [], "null_draws": len(gnull)}
    for r in gnmt:
        entry = {"tag": r["tag"], "rsid1": r.get("rsid1"), "rsid2": r.get("rsid2"),
                 "pos1": r["pos1"], "ref1": r["ref1"], "alt1": r["alt1"],
                 "pos2": r["pos2"], "ref2": r["ref2"], "alt2": r["alt2"],
                 "palindromic1": r.get("palindromic1"), "palindromic2": r.get("palindromic2"),
                 "separation_bp": r["separation_bp"], "rna_readout": r["rna_readout"], "is_primary": bool(r.get("is_primary"))}
        for ch in CHANNELS:
            nn = np.array([x[f"{ch}_residual"] for x in gnull], float)
            lo, hi = (float(np.nanquantile(nn, 0.025)), float(np.nanquantile(nn, 0.975))) if nn.size else (math.nan, math.nan)
            res = r[f"{ch}_residual"]
            # native scale: the summed predicted signal in the readout window, before any log
            for arm in ("ref", "v1", "v2", "joint"):
                entry[f"{ch}_{arm}_sum"] = r.get(f"{ch}_{arm}_sum")
            entry[f"{ch}_additive_expectation_sum"] = (r.get(f"{ch}_ref_sum", math.nan)
                                                       * 2.0 ** (r[f"{ch}_v1_log2"] + r[f"{ch}_v2_log2"]))
            entry[f"{ch}_v1_log2"] = r[f"{ch}_v1_log2"]
            entry[f"{ch}_v2_log2"] = r[f"{ch}_v2_log2"]
            entry[f"{ch}_joint_log2"] = r[f"{ch}_joint_log2"]
            entry[f"{ch}_additive_expectation_log2"] = r[f"{ch}_v1_log2"] + r[f"{ch}_v2_log2"]
            entry[f"{ch}_residual"] = res
            entry[f"{ch}_null_central95_lo"] = lo
            entry[f"{ch}_null_central95_hi"] = hi
            entry[f"{ch}_inside_central95"] = bool(lo <= res <= hi) if (lo == lo and res == res) else None
            entry[f"{ch}_p_empirical_vs_gnmt_null"] = R.empirical_two_sided_p(res, nn)
        f2["pairs"].append(entry)
    prim_entry = next((e for e in f2["pairs"] if e["is_primary"]), None)
    if prim_entry:
        f2["F2.1_verdict"] = {ch: prim_entry[f"{ch}_inside_central95"] for ch in CHANNELS}
        f2["F2.1_met"] = all(prim_entry[f"{ch}_inside_central95"] for ch in CHANNELS)
    # the six GNMT-region pairs x four channels are 24 tests; the central-95% read is uncorrected, so BH the family
    keys = [(e["tag"], ch) for e in f2["pairs"] for ch in CHANNELS]
    pv = np.array([next(e for e in f2["pairs"] if e["tag"] == t)[f"{c}_p_empirical_vs_gnmt_null"] for t, c in keys], float)
    rj = R.bh_reject(pv, R.BH_Q)
    nn = len(gnull)
    f2["gnmt_family_bh"] = {
        "family": "6 GNMT-region pairs x 4 channels = 24 tests against the 100-draw chr6 matched null",
        "q": R.BH_Q, "n_tests": int(pv.size), "p_floor": 1.0 / (nn + 1),
        "pairs_that_must_sit_AT_the_floor_together_for_any_rejection": int(math.ceil((1.0 / (nn + 1)) * pv.size / R.BH_Q)),
        "n_at_the_floor": int(np.sum(np.abs(pv - 1.0 / (nn + 1)) < 1e-12)),
        "n_outside_central95_UNCORRECTED": int(sum(1 for e in f2["pairs"] for ch in CHANNELS if e[f"{ch}_inside_central95"] is False)),
        "n_bh_rejected": int(rj.sum()),
        "rejected": [f"{keys[i][0]}:{keys[i][1]}" for i in np.flatnonzero(rj)],
        "min_p": float(np.nanmin(pv))}
    for i, (t, c) in enumerate(keys):
        for e in f2["pairs"]:
            if e["tag"] == t:
                e[f"{c}_bh_reject_gnmt_family"] = bool(rj[i])
    json.dump(f2, (OUT / "tables" / f"f2_gnmt_verdict_{tag}.json").open("w"), indent=1, default=float)
    write_tsv(OUT / "tables" / f"f2_gnmt_effects_{tag}.tsv", f2["pairs"])
    write_tsv(OUT / "tables" / f"f2_gnmt_null_effects_{tag}.tsv", gnull)

    # ---------------------------------------------------------------- F3
    per_channel = {}
    for ch in CHANNELS:
        per_channel[ch] = channel_verdict(ch, obs, nul, "200 highest weight-product P5 pairs, corrected 1,000-draw matched null")
    summary["f3_corrected"] = per_channel

    # readout-matched RNA diagnostic (amendment A1)
    obs_c20 = [r for r in obs if r.get("rna_readout") == "central_20kb"]
    summary["rna_readout_counts"] = {"observed_gene_span": sum(1 for r in obs if r.get("rna_readout") == "gene_span"),
                                     "observed_central_20kb": len(obs_c20),
                                     "null_gene_span": sum(1 for r in nul if r.get("rna_readout") == "gene_span"),
                                     "null_central_20kb": sum(1 for r in nul if r.get("rna_readout") == "central_20kb")}
    summary["rna_readout_matched_diagnostic"] = (channel_verdict("rna", [dict(r) for r in obs_c20], nul,
                                                                 "observed pairs whose RNA readout is also central_20kb")
                                                 if len(obs_c20) >= 5 else {"n_pairs": len(obs_c20), "testable": False})

    # is pooling the null across decile bins safe? (amendment 02)
    from scipy import stats as sstats
    pooling = {}
    for ch in CHANNELS:
        nb = {}
        for b in sorted({int(r.get("decile_bin", -1)) for r in nul}):
            a = np.array([abs(r[f"{ch}_residual"]) for r in nul if int(r.get("decile_bin", -1)) == b], float)
            nb[b] = {"n": int(a.size), "median_abs_residual": float(np.nanmedian(a)) if a.size else math.nan}
        sn = np.array([r["separation_bp"] for r in nul], float)
        an = np.array([abs(r[f"{ch}_residual"]) for r in nul], float)
        so = np.array([r["separation_bp"] for r in obs], float)
        ao = np.array([abs(r[f"{ch}_residual"]) for r in obs], float)
        ok_n, ok_o = ~np.isnan(an), ~np.isnan(ao)
        rn = sstats.spearmanr(sn[ok_n], an[ok_n]) if ok_n.sum() > 2 else None
        ro = sstats.spearmanr(so[ok_o], ao[ok_o]) if ok_o.sum() > 2 else None
        pooling[ch] = {"null_by_decile_bin": nb,
                       "null_spearman_sep_vs_abs_residual": [float(rn.statistic), float(rn.pvalue)] if rn else None,
                       "observed_spearman_sep_vs_abs_residual": [float(ro.statistic), float(ro.pvalue)] if ro else None}
    summary["null_pooling_diagnostic"] = {
        "why": ("each draw is matched to one sampled observed pair's chromosome and separation decile, so the pooled null is "
                "matched to the family's joint distribution but not conditionally to each pair. If |residual| does not vary "
                "with separation decile, pooling costs nothing; if it does, the per-pair p-values are mis-matched and that is "
                "reported. This also tests the P5 prespecification's H4 (residual magnitude does not grow with separation)."),
        "per_channel": pooling}

    # blocks: the resampling unit
    blocks = sorted({r.get("analysis_block", "") for r in obs})
    summary["blocks"] = {"n_distinct_analysis_blocks": len([b for b in blocks if b]),
                         "n_distinct_signals": len({r.get("signal_uid", "") for r in obs}),
                         "note": "pairs inside one 1-Mb analysis block share sequence and are not independent tests; the BH family of 200 pairs therefore overstates the number of independent tests"}

    # POST HOC (amendment 07): every stage-1 rejection sits below 2 kb of separation, where the chromatin
    # readout's two +/-1 kb windows OVERLAP. Split the family at 2 kb and keep the null in the same stratum.
    SEP_SPLIT = 2 * R.FLANK
    strat = {}
    for label, keep_o, keep_n in (("separation_under_2kb_readout_windows_overlap",
                                   [r for r in obs if r["separation_bp"] < SEP_SPLIT],
                                   [r for r in nul if r["separation_bp"] < SEP_SPLIT]),
                                  ("separation_2kb_or_more_readout_windows_disjoint",
                                   [r for r in obs if r["separation_bp"] >= SEP_SPLIT],
                                   [r for r in nul if r["separation_bp"] >= SEP_SPLIT])):
        blk = {}
        for ch in CHANNELS:
            blk[ch] = channel_verdict(ch, [dict(r) for r in keep_o], keep_n, label)
        strat[label] = {"n_observed": len(keep_o), "n_null": len(keep_n), "per_channel": blk}
    bands = [(1, 3), (4, 27), (28, 127), (128, 305), (306, 1085), (1086, 1999), (2000, 15134), (15135, 10 ** 9)]
    feas = []
    for lo, hi in bands:
        feas.append({"band_bp": f"{lo}-{hi if hi < 10 ** 9 else 'max'}",
                     "observed_pairs": int(sum(1 for r in obs if lo <= r["separation_bp"] <= hi)),
                     "null_draws": int(sum(1 for r in nul if lo <= r["separation_bp"] <= hi)),
                     "observed_share": round(sum(1 for r in obs if lo <= r["separation_bp"] <= hi) / max(1, len(obs)), 4),
                     "null_share": round(sum(1 for r in nul if lo <= r["separation_bp"] <= hi) / max(1, len(nul)), 4)})
    summary["separation_stratified_POST_HOC"] = {
        "status": "POST HOC. Added after every stage-1 rejection turned out to lie below 2 kb of separation.",
        "mechanism": ("the ATAC, DNASE and H3K27ac readout is the union of +/-1 kb around each variant, so below "
                      "2 kb of separation the two windows overlap and below a few bp they coincide. Two edits read "
                      "out over the same bases have no reason to add on a log2 of a summed signal, so a residual "
                      "appears from the readout definition rather than from any interaction."),
        "split_bp": SEP_SPLIT, "strata": strat,
        "why_the_decile_match_cannot_fix_it": ("the separation decile bin at the bottom is [1, 27.7] bp, and inside it the "
                                               "observed pairs sit at 1-3 bp while common EUR SNV pairs that close are rare. "
                                               "The band table below shows the residual mismatch the decile match leaves."),
        "separation_band_feasibility": feas}

    # POST HOC (amendment 06): the residual is a function of the two arms, and the observed family's arms are
    # not the same size as the null's. Match the null on arm magnitude and see whether the verdict survives.
    from scipy import stats as sstats2
    LO, HI, MIN_MATCHED = 2.0 / 3.0, 3.0 / 2.0, 50
    mag = {}
    for ch in CHANNELS:
        om = np.array([abs(r[f"{ch}_v1_log2"]) + abs(r[f"{ch}_v2_log2"]) for r in obs], float)
        nm = np.array([abs(r[f"{ch}_v1_log2"]) + abs(r[f"{ch}_v2_log2"]) for r in nul], float)
        ores = np.array([abs(r[f"{ch}_residual"]) for r in obs], float)
        nres = np.array([abs(r[f"{ch}_residual"]) for r in nul], float)
        ko, kn = ~np.isnan(om) & ~np.isnan(ores), ~np.isnan(nm) & ~np.isnan(nres)
        so = sstats2.spearmanr(om[ko], ores[ko]) if ko.sum() > 2 else None
        sn = sstats2.spearmanr(nm[kn], nres[kn]) if kn.sum() > 2 else None
        order = np.argsort(nm)
        nm_s, nres_s = nm[order], nres[order]
        p = np.full(om.shape, np.nan)
        nmatch = np.zeros(om.shape, int)
        for i in range(om.size):
            if not ko[i] or om[i] <= 0:
                continue
            a = int(np.searchsorted(nm_s, om[i] * LO, side="left"))
            b = int(np.searchsorted(nm_s, om[i] * HI, side="right"))
            sub = nres_s[a:b]
            sub = sub[~np.isnan(sub)]
            nmatch[i] = sub.size
            if sub.size >= MIN_MATCHED:
                p[i] = (np.sum(sub >= ores[i]) + 1) / (sub.size + 1)
        rej = R.bh_reject(p, R.BH_Q)
        for i, r in enumerate(obs):
            r[f"{ch}_p_magnitude_matched"] = float(p[i])
            r[f"{ch}_n_magnitude_matched_null"] = int(nmatch[i])
            r[f"{ch}_bh_reject_magnitude_matched"] = bool(rej[i])
        mag[ch] = {"observed_median_arm_magnitude": float(np.nanmedian(om)),
                   "null_median_arm_magnitude": float(np.nanmedian(nm)),
                   "observed_over_null_arm_magnitude": float(np.nanmedian(om) / np.nanmedian(nm)) if np.nanmedian(nm) else math.nan,
                   "spearman_arm_magnitude_vs_abs_residual_observed": [float(so.statistic), float(so.pvalue)] if so else None,
                   "spearman_arm_magnitude_vs_abs_residual_null": [float(sn.statistic), float(sn.pvalue)] if sn else None,
                   "n_pairs_testable": int(np.sum(~np.isnan(p))), "n_pairs_untestable_too_few_matched_draws": int(np.sum(np.isnan(p) & ko)),
                   "median_matched_null_draws": float(np.median(nmatch[nmatch > 0])) if np.any(nmatch > 0) else 0.0,
                   "n_pairs_bh": int(rej.sum()), "min_p": float(np.nanmin(p)) if np.any(~np.isnan(p)) else math.nan,
                   "p_floor_given_matched_draws": float(1.0 / (np.median(nmatch[nmatch > 0]) + 1)) if np.any(nmatch > 0) else math.nan}
    summary["magnitude_matched_diagnostic_POST_HOC"] = {
        "status": "POST HOC. Added after the stage-1 result showed the observed family's median single-variant effect is "
                  "larger than the null's in three of four channels. It bounds the primary verdict; it does not replace it.",
        "why": "the residual is joint - v1 - v2, a function of the two single-variant arms. If the observed pairs carry "
               "larger arms than the null pairs, a larger residual follows from arm size and not from departure from "
               "additivity. The matched comparison asks whether a pair's residual is extreme among null draws whose arms "
               "are the same size.",
        "rule": f"per observed pair, the null draws whose |v1|+|v2| lies in [{LO:.3f}, {HI:.3f}] times the pair's own; "
                f"at least {MIN_MATCHED} matched draws or the pair is untestable; then BH at q=0.10 over the testable pairs.",
        "per_channel": mag}

    # family-level contrast with an interval (amendment 04): are Resource pairs MORE additive than matched
    # null pairs? The existing deposit asserts this from two medians with no interval.
    rngb = np.random.default_rng(20260914)
    blocks_of = {}
    for i, r in enumerate(obs):
        blocks_of.setdefault(r.get("analysis_block", f"_nb{i}"), []).append(i)
    bkeys = sorted(blocks_of)
    contrast = {}
    for ch in CHANNELS:
        ao = np.array([abs(r[f"{ch}_residual"]) for r in obs], float)
        an = np.array([abs(r[f"{ch}_residual"]) for r in nul], float)
        point = float(np.nanmedian(ao) - np.nanmedian(an))
        draws = np.empty(10000)
        for b in range(10000):
            bi = rngb.integers(0, len(bkeys), len(bkeys))
            idx = np.concatenate([blocks_of[bkeys[j]] for j in bi])
            ni = rngb.integers(0, an.size, an.size)
            draws[b] = np.nanmedian(ao[idx]) - np.nanmedian(an[ni])
        contrast[ch] = {"observed_median_abs_residual": float(np.nanmedian(ao)),
                        "null_median_abs_residual": float(np.nanmedian(an)),
                        "difference": point,
                        "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
                        "excludes_zero": bool(np.percentile(draws, 2.5) > 0 or np.percentile(draws, 97.5) < 0)}
    summary["family_level_contrast"] = {
        "statistic": "median |residual| of the 200 observed pairs minus median |residual| of the matched null draws",
        "resampling": ("10,000 draws, seed 20260914; observed side resampled over its 1-Mb analysis blocks (pairs inside "
                       "one block are not independent), null side resampled over its independent draws"),
        "n_blocks": len(bkeys), "per_channel": contrast,
        "why": ("the existing deposit states 'Resource variant pairs are if anything more additive than random pairs at the "
                "same separation' from two medians with no interval and against a null that matched separation only "
                "marginally. This is that comparison with an interval and a jointly matched null.")}

    # p-floor arithmetic (amendment 03): can either BH test reject at all?
    def floor_block(n_draws: int, m: int, q: float = R.BH_Q) -> dict:
        # BH is a step-up: it rejects k pairs when p_(k) <= q*k/m. A p-value floored at 1/(draws+1) can still
        # be rejected, but only if enough pairs sit AT the floor together. k_min is that number.
        floor = 1.0 / (n_draws + 1)
        k_min = int(math.ceil(floor * m / q))
        return {"null_draws": int(n_draws), "family_size": int(m), "q": q,
                "p_floor": floor, "bh_threshold_at_rank_1": q / m,
                "pairs_that_must_sit_AT_the_floor_together_for_any_rejection": k_min,
                "any_rejection_possible_in_principle": bool(k_min <= m),
                "draws_needed_so_a_SINGLE_pair_could_reject": int(math.ceil(m / q) - 1)}
    summary["p_floor_arithmetic"] = {
        "existing_p5_deposit_300_draws_1572_pairs": floor_block(300, 1572),
        "this_package": floor_block(len(nul), len(obs)),
        "meaning": ("an empirical p from N draws cannot fall below 1/(N+1). Benjamini-Hochberg is a step-up test, so a "
                    "floored p is still rejectable, but only when at least k_min pairs reach the floor together, where "
                    "k_min = ceil(p_floor * m / q). A single floored pair can reject only if the floor is below q/m."),
        "correction_note": ("amendment 03 stated this as 'no pair could have been rejected whatever the model returned', "
                            "comparing the floor only with the rank-1 threshold q/m. That is wrong for a step-up test. "
                            "The correct bound is k_min: the deposit's 300-draw null over 1,572 pairs could have rejected "
                            "only if at least 53 pairs simultaneously exceeded all 300 null draws; this package's "
                            "1,000-draw null over the scored family needs 2. The correction is recorded in "
                            "prespec/f2_00_amendment_05.json and does not change any decision rule."),
        "n_pairs_at_the_p_floor_per_channel": {ch: int(sum(1 for r in obs if r.get(f"{ch}_p_empirical") is not None
                                                           and abs(r.get(f"{ch}_p_empirical", 1.0) - 1.0 / (len(nul) + 1)) < 1e-12))
                                               for ch in CHANNELS}}

    summary["existing_p5_deposit"] = {"family_size_pairs": 1572, "n_null": 300,
                                      "p_floor": 1 / 301,
                                      "reported": "0 pairs and 0 signals survive BH q<0.10 in every channel; the both_variants_active stratum was declared untestable (1 to 3 matched null pairs)",
                                      "recomputed_null_here": p5_null_reference()}
    summary["claim_boundary"] = ("model prediction only. Every number is the model's response to substitutions in one reference "
                                 "sequence. Levels 2 to 4 of the plan's evidence ladder (measured reporter interaction, measured "
                                 "endogenous interaction, measured state-dependent interaction) have no substrate for these pairs, "
                                 "so nothing here is an epistasis claim.")
    json.dump(summary, (OUT / "tables" / f"f3_summary_{tag}.json").open("w"), indent=1, default=float)
    write_tsv(OUT / "tables" / f"f3_observed_effects_{tag}.tsv", obs)
    write_tsv(OUT / "tables" / f"f3_null_effects_{tag}.tsv", nul)
    print(json.dumps({"f2_met": f2.get("F2.1_met"),
                      "f3_bh_rejections": {ch: per_channel[ch]["all_pairs"].get("n_pairs_bh") for ch in CHANNELS}}, indent=1))


if __name__ == "__main__":
    main()
