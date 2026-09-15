#!/usr/bin/env python3
"""f2-haplotype-v2 step 3: the repaired verdicts and every accounting correction. No network use.

Reads
  this package : raw/scored_genespan_null.jsonl          the 100 gene-span-matched chr6 null draws
  v1 read-only : raw/scored_rows.jsonl                   the six GNMT pairs, the 199 observed pairs,
                                                         the 100 central-20kb chr6 draws, the 2,500 draws
Writes
  tables/f2v2_gnmt_verdict.json          F2.1 per channel against the repaired null, beside v1's
  tables/f2v2_genespan_null_effects.tsv  the repaired null's rows
  tables/f2v2_both_active_stratum.json   defect 2: p-values with the null INSIDE the stratum
  tables/f2v2_accounting.json            call counts, checkpoint deposits, prespec and amendment timestamps
  tables/f2v2_readout_geometry.json      the RNA readout is ONE window at every separation
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import math
import os
import pathlib
import re

import numpy as np
from scipy import stats as sstats

import f2_recipe as R

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
V1 = PROJECT / "GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z"
CH = R.CHANNELS
GNMT_SPAN = (42960690, 42963883)
GNMT_REF_RNA_SUM = None          # read from v1's own table, never hard-coded


def load_jsonl(path: pathlib.Path) -> list[dict]:
    rows, seen = [], set()
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r["tag"] in seen:
                continue
            seen.add(r["tag"])
            rows.append(r)
    return rows


def group_v1(rows: list[dict]) -> dict[str, list[dict]]:
    g = {"gnmt": [], "gnmtnull": [], "obs": [], "null": []}
    for r in rows:
        t = r["tag"]
        key = ("gnmtnull" if t.startswith("gnmtnull") else
               "gnmt" if t.startswith("gnmt") else
               "null" if t.startswith("null") else "obs")
        g[key].append(r)
    return g


def write_tsv(path: pathlib.Path, rows: list[dict]) -> None:
    if not rows:
        return
    cols = sorted({k for r in rows for k in r})
    lead = [c for c in ("tag", "chrom", "pos1", "pos2", "separation_bp", "ensembl", "gene_name",
                        "span_width_bp", "rna_readout", "rna_readout_bp") if c in cols]
    cols = lead + [c for c in cols if c not in lead]
    with open(path, "w") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def arr(rows: list[dict], key: str) -> np.ndarray:
    return np.array([float(r.get(key, math.nan)) if r.get(key, "") not in ("", None) else math.nan for r in rows], float)


# ---------------------------------------------------------------- F2.1 against a given null
def f2_block(gnmt: list[dict], nul: list[dict], label: str) -> dict:
    out = {"null_label": label, "n_null_draws": len(nul), "p_floor": 1.0 / (len(nul) + 1), "pairs": []}
    for r in gnmt:
        e = {"tag": r["tag"], "rsid1": r.get("rsid1"), "rsid2": r.get("rsid2"),
             "separation_bp": r["separation_bp"], "is_primary": bool(r.get("is_primary")),
             "rna_readout": r["rna_readout"]}
        for ch in CH:
            nn = arr(nul, f"{ch}_residual")
            lo, hi = float(np.nanquantile(nn, 0.025)), float(np.nanquantile(nn, 0.975))
            res = float(r[f"{ch}_residual"])
            e[f"{ch}_residual"] = res
            e[f"{ch}_null_central95_lo"] = lo
            e[f"{ch}_null_central95_hi"] = hi
            e[f"{ch}_inside_central95"] = bool(lo <= res <= hi)
            e[f"{ch}_p_empirical"] = R.empirical_two_sided_p(res, nn)
            e[f"{ch}_null_q95_abs_residual"] = float(np.nanquantile(np.abs(nn), 0.95))
        out["pairs"].append(e)
    prim = next(e for e in out["pairs"] if e["is_primary"])
    out["F2.1_verdict_per_channel"] = {ch: prim[f"{ch}_inside_central95"] for ch in CH}
    out["F2.1_met"] = all(prim[f"{ch}_inside_central95"] for ch in CH)
    out["F2.1_p_per_channel"] = {ch: prim[f"{ch}_p_empirical"] for ch in CH}
    out["null_q95_abs_residual"] = {ch: float(np.nanquantile(np.abs(arr(nul, f"{ch}_residual")), 0.95)) for ch in CH}
    out["null_median_abs_residual"] = {ch: float(np.nanmedian(np.abs(arr(nul, f"{ch}_residual")))) for ch in CH}
    keys = [(e["tag"], ch) for e in out["pairs"] for ch in CH]
    pv = np.array([next(e for e in out["pairs"] if e["tag"] == t)[f"{c}_p_empirical"] for t, c in keys], float)
    rj = R.bh_reject(pv, R.BH_Q)
    floor = 1.0 / (len(nul) + 1)
    out["gnmt_family_bh"] = {
        "family": f"6 GNMT-region pairs x 4 channels = 24 tests against the {len(nul)}-draw {label}",
        "q": R.BH_Q, "n_tests": int(pv.size), "p_floor": floor,
        "pairs_that_must_sit_AT_the_floor_together_for_any_rejection": int(math.ceil(floor * pv.size / R.BH_Q)),
        "n_at_the_floor": int(np.sum(np.abs(pv - floor) < 1e-12)),
        "n_outside_central95_UNCORRECTED": int(sum(1 for e in out["pairs"] for ch in CH
                                                   if e[f"{ch}_inside_central95"] is False)),
        "n_bh_rejected": int(rj.sum()),
        "rejected": [f"{keys[i][0]}:{keys[i][1]}" for i in np.flatnonzero(rj)],
        "min_p": float(np.nanmin(pv))}
    return out


# ---------------------------------------------------------------- defect 2
def both_active_block(obs: list[dict], nul: list[dict], stage: str) -> dict:
    res = {"stage": stage, "n_obs": len(obs), "n_null": len(nul), "per_channel": {}}
    for ch in CH:
        o = arr(obs, f"{ch}_residual")
        n = arr(nul, f"{ch}_residual")
        v1, v2 = arr(obs, f"{ch}_v1_log2"), arr(obs, f"{ch}_v2_log2")
        nv1, nv2 = arr(nul, f"{ch}_v1_log2"), arr(nul, f"{ch}_v2_log2")
        thr = R.effect_threshold(nv1, nv2)
        keep_o = R.effect_stratum(v1, v2, thr)
        keep_n = R.effect_stratum(nv1, nv2, thr)
        n_in = n[keep_n]
        floor_pooled = 1.0 / (int(np.sum(~np.isnan(n))) + 1)
        floor_in = 1.0 / (int(np.sum(~np.isnan(n_in))) + 1) if n_in.size else math.nan

        p_pooled = np.full(o.shape, np.nan)
        p_in = np.full(o.shape, np.nan)
        for i in np.flatnonzero(keep_o):
            p_pooled[i] = R.empirical_two_sided_p(o[i], n)
            p_in[i] = R.empirical_two_sided_p(o[i], n_in) if n_in.size else math.nan
        rej_pooled = R.bh_reject(p_pooled, R.BH_Q)
        rej_in = R.bh_reject(p_in, R.BH_Q)
        gate = R.stratum_is_testable(int(keep_o.sum()), int(keep_n.sum()))
        # marginal rates behind the "by construction" claim
        a = np.abs(np.concatenate([nv1, nv2]))
        a = a[~np.isnan(a)]
        frac1 = float(np.mean(np.abs(nv1[~np.isnan(nv1)]) > thr))
        frac2 = float(np.mean(np.abs(nv2[~np.isnan(nv2)]) > thr))
        ok = ~np.isnan(nv1) & ~np.isnan(nv2)
        rho = sstats.spearmanr(np.abs(nv1[ok]), np.abs(nv2[ok]))
        res["per_channel"][ch] = {
            "single_variant_effect_threshold_null_q95": float(thr),
            "n_observed_in_stratum": int(keep_o.sum()),
            "n_null_draws_in_stratum": int(keep_n.sum()),
            "null_draws_in_stratum_share_percent": round(100.0 * float(keep_n.sum()) / max(1, len(nul)), 4),
            "share_expected_if_the_two_arms_were_independent_percent": round(100.0 * frac1 * frac2, 4),
            "share_if_both_arms_were_independent_at_exactly_5_percent": 0.25,
            "fraction_of_null_v1_arms_above_threshold": round(frac1, 4),
            "fraction_of_null_v2_arms_above_threshold": round(frac2, 4),
            "spearman_abs_v1_vs_abs_v2_in_null": [float(rho.statistic), float(rho.pvalue)],
            "p_floor_pooled_null": floor_pooled,
            "p_floor_null_inside_the_stratum": floor_in,
            "second_smallest_attainable_p_inside_the_stratum": 2.0 * floor_in if floor_in == floor_in else math.nan,
            "bh_threshold_at_rank_1_in_this_stratum": R.BH_Q / max(1, int(keep_o.sum())),
            "bh_threshold_at_the_last_rank": R.BH_Q,
            "pairs_that_must_sit_AT_the_in_stratum_floor_together_for_any_rejection":
                int(math.ceil(floor_in * int(keep_o.sum()) / R.BH_Q)) if floor_in == floor_in and keep_o.sum() else None,
            "recipe_min_matched_null": R.MIN_MATCHED_NULL, "testable_under_the_recipe_gate": gate["testable"],
            "V1_AS_PUBLISHED_p_against_the_pooled_null": {
                "n_bh_rejected": int(rej_pooled.sum()),
                "min_p": float(np.nanmin(p_pooled)) if np.any(~np.isnan(p_pooled)) else math.nan,
                "n_pairs_at_the_pooled_floor": int(np.sum(np.abs(p_pooled - floor_pooled) < 1e-12))},
            "REPAIRED_p_against_the_null_inside_the_stratum": {
                "n_bh_rejected_gated": int(rej_in.sum()) if gate["testable"] else None,
                "n_bh_rejected_if_the_20_draw_gate_is_ignored": int(rej_in.sum()),
                "min_p": float(np.nanmin(p_in)) if np.any(~np.isnan(p_in)) else math.nan,
                "n_pairs_at_the_in_stratum_floor": int(np.sum(np.abs(p_in - floor_in) < 1e-12)) if floor_in == floor_in else None,
                "p_values": [None if np.isnan(x) else round(float(x), 6) for x in p_in[keep_o]],
                "bh_family_size": int(keep_o.sum())},
        }
    return res


# ---------------------------------------------------------------- readout geometry
def readout_geometry(obs: list[dict], spans: dict) -> dict:
    rows = []
    for r in obs:
        p1, p2 = int(r["pos1"]), int(r["pos2"])
        mid = (p1 + p2) // 2
        start0 = max(0, mid - R.WINDOW // 2)
        end = start0 + R.WINDOW
        span = spans.get(r.get("ensembl", ""))
        use_span = bool(span and span[0] == r["chrom"] and span[2] > start0 and span[1] < end)
        rna_regions = [(span[1], span[2])] if use_span else [(mid - 10_000, mid + 10_000)]
        w_rna = int(R.window_mask(R.WINDOW, start0, rna_regions).sum())
        w_loc = int(R.window_mask(R.WINDOW, start0,
                                  [(p1 - R.FLANK, p1 + R.FLANK), (p2 - R.FLANK, p2 + R.FLANK)]).sum())
        a, b = rna_regions[0]
        in1, in2 = bool(a <= p1 < b), bool(a <= p2 < b)
        rows.append({"tag": r["tag"], "chrom": r["chrom"], "pos1": p1, "pos2": p2,
                     "ref1": r.get("ref1", ""), "alt1": r.get("alt1", ""),
                     "ref2": r.get("ref2", ""), "alt2": r.get("alt2", ""),
                     "signal_uid": r.get("signal_uid", ""), "analysis_block": r.get("analysis_block", ""),
                     "separation_bp": int(r["separation_bp"]),
                     "rna_readout": "gene_span" if use_span else "central_20kb",
                     "rna_readout_regions": len(rna_regions), "rna_readout_bp": w_rna,
                     "rna_readout_contains_variant_1": in1, "rna_readout_contains_variant_2": in2,
                     "rna_readout_contains_neither_variant": bool(not in1 and not in2),
                     "bp_from_variant_1_to_the_rna_readout": 0 if in1 else int(min(abs(p1 - a), abs(p1 - b))),
                     "local_readout_regions_before_union": 2, "local_readout_bp": w_loc,
                     "local_windows_overlap": bool(w_loc < 2 * 2 * R.FLANK),
                     "rna_bh_reject_all": str(r.get("rna_bh_reject_all", "")).lower() in ("true", "1"),
                     "rna_residual": float(r["rna_residual"]) if r.get("rna_residual", "") not in ("", None) else math.nan})
    sep = np.array([x["separation_bp"] for x in rows], float)
    wl = np.array([x["local_readout_bp"] for x in rows], float)
    wr = np.array([x["rna_readout_bp"] for x in rows], float)
    rej = [x for x in rows if x["rna_bh_reject_all"]]
    span_rows = [x for x in rows if x["rna_readout"] == "gene_span"]
    ws = np.array([x["rna_readout_bp"] for x in span_rows], float)
    return {
        "what_this_measures": ("the RNA readout is ONE region at every separation (the gene span, or the central "
                              "20 kb); the ATAC, DNASE and H3K27ac readout is the UNION of two +/-1 kb windows, so "
                              "its width shrinks as separation falls below 2 kb and the two arms come to be read out "
                              "over shared bases. The overlap mechanism therefore cannot explain an RNA rejection."),
        "n_pairs": len(rows),
        "rna_readout_regions_per_pair": sorted({x["rna_readout_regions"] for x in rows}),
        "local_readout_bp_vs_separation_spearman": list(map(float, sstats.spearmanr(sep, wl)[:2])),
        "rna_readout_bp_vs_separation_spearman": list(map(float, sstats.spearmanr(sep, wr)[:2])),
        "n_pairs_with_overlapping_local_windows": int(sum(1 for x in rows if x["local_windows_overlap"])),
        "n_pairs_with_overlapping_rna_windows": 0,
        "rna_readout_bp": {"min": float(wr.min()), "median": float(np.median(wr)), "max": float(wr.max())},
        "n_gene_span": len(span_rows), "n_central_20kb": len(rows) - len(span_rows),
        "gene_span_readout_bp": ({"min": float(ws.min()), "q25": float(np.percentile(ws, 25)),
                                  "median": float(np.median(ws)), "q75": float(np.percentile(ws, 75)),
                                  "max": float(ws.max())} if ws.size else None),
        "v1_null_rna_readout_bp_for_2470_of_2500_draws": 20000,
        "why_that_matters": ("in the RNA channel the observed family's readout width spans three orders of magnitude, "
                            "because it is the target gene's own span, while 2,470 of v1's 2,500 null draws read out "
                            "over a fixed 20,000 bp. The RNA arm of the 199-pair test is therefore unmatched on "
                            "readout width pair by pair, not only for the GNMT pair."),
        "n_pairs_whose_rna_readout_contains_neither_variant": int(sum(1 for x in rows if x["rna_readout_contains_neither_variant"])),
        "n_pairs_whose_rna_readout_contains_both_variants": int(sum(1 for x in rows if x["rna_readout_contains_variant_1"]
                                                                   and x["rna_readout_contains_variant_2"])),
        "rna_bh_rejections_v1": [{k: x[k] for k in ("tag", "signal_uid", "analysis_block", "separation_bp",
                                                    "rna_readout", "rna_readout_bp", "local_readout_bp",
                                                    "rna_readout_contains_neither_variant",
                                                    "bp_from_variant_1_to_the_rna_readout", "rna_residual")}
                                 for x in rej],
        "all_rna_rejections_below_2kb": bool(all(x["separation_bp"] < 2 * R.FLANK for x in rej)) if rej else None,
        "n_rna_rejections": len(rej),
        "duplicate_rows_in_the_199_pair_family": duplicates(rows),
        "rejections_as_rows_and_as_distinct_variant_pairs": rejection_multiplicity(obs),
    }


def _pair_key(r: dict) -> tuple:
    p1, p2 = int(r["pos1"]), int(r["pos2"])
    lo, hi = min(p1, p2), max(p1, p2)
    return (r["chrom"], lo, hi)


def duplicates(rows: list[dict]) -> dict:
    from collections import Counter
    c = Counter(_pair_key(r) for r in rows)
    return {"rows": len(rows), "distinct_variant_pairs": len(c),
            "variant_pairs_appearing_more_than_once": int(sum(1 for v in c.values() if v > 1)),
            "rows_they_account_for": int(sum(v for v in c.values() if v > 1)),
            "multiplicity_histogram": {str(k): int(v) for k, v in sorted(Counter(c.values()).items())},
            "why_that_matters": ("the same two positions enter the 199-test BH family once per COLOC signal that "
                                "names them. A count of rejected ROWS is therefore not a count of variant pairs, and "
                                "the BH family size of 199 counts the same sequence more than once.")}


def rejection_multiplicity(obs_tsv: list[dict]) -> dict:
    out = {}
    for ch in CH:
        rej = [r for r in obs_tsv if str(r.get(f"{ch}_bh_reject_all", "")).lower() in ("true", "1")]
        out[ch] = {"rejected_rows": len(rej),
                   "distinct_variant_pairs": len({_pair_key(r) for r in rej}),
                   "distinct_analysis_blocks": len({r.get("analysis_block", "") for r in rej})}
    return out


# ---------------------------------------------------------------- accounting
def accounting() -> dict:
    ck = V1 / "raw" / "scored_rows.jsonl"
    tags = [json.loads(l)["tag"] for l in open(ck) if l.strip()]
    logs = sorted((V1 / "logs").glob("agp_f2-haplotype-2*.out"))
    per_job = []
    for p in logs:
        txt = p.read_text()
        jid = re.search(r"-(\d+)\.out$", p.name).group(1)
        done = re.search(r"DONE: (\d+) calls in ([\d.]+) min", txt)
        already = re.search(r"checkpoint: (\d+) pairs already scored", txt)
        worklist = re.search(r"work list: (\d+) pairs, (\d+) to score", txt)
        last_calls = re.findall(r"\((\d+) calls, ([\d.]+)/min\)", txt)
        cancelled = "CANCELLED" in txt
        per_job.append({"jobid": jid, "cancelled": cancelled,
                        "checkpoint_pairs_at_start": int(already.group(1)) if already else None,
                        "work_list_pairs": int(worklist.group(1)) if worklist else None,
                        "pairs_to_score": int(worklist.group(2)) if worklist else None,
                        "calls_reported_at_DONE": int(done.group(1)) if done else None,
                        "minutes_reported_at_DONE": float(done.group(2)) if done else None,
                        "last_progress_line_calls": int(last_calls[-1][0]) if last_calls else None,
                        "max_rate_per_min_printed": max(float(x[1]) for x in last_calls) if last_calls else None})
    cancel = next(j for j in per_job if j["cancelled"])
    successor = min((j for j in per_job if j["checkpoint_pairs_at_start"]), key=lambda j: j["jobid"])
    deposited = successor["checkpoint_pairs_at_start"]
    # the primary GNMT pair costs one extra call (the determinism probe), every other pair costs four
    primary_in_deposit = any(t == "gnmt|rs2296805|rs2296804" for t in tags[:deposited])
    calls_cancelled = 4 * deposited + (1 if primary_in_deposit else 0)
    completed = [j for j in per_job if j["calls_reported_at_DONE"]]
    total = calls_cancelled + sum(j["calls_reported_at_DONE"] for j in completed)
    return {
        "per_job": per_job,
        "cancelled_job": {
            "jobid": cancel["jobid"],
            "pairs_deposited_in_the_checkpoint": deposited,
            "evidence": (f"job {successor['jobid']}'s own log line 'checkpoint: {deposited} pairs already scored', "
                         f"and the first {deposited} tags of raw/scored_rows.jsonl are "
                         f"{tags[:deposited]}"),
            "predict_sequence_calls": calls_cancelled,
            "how_that_is_derived": (f"4 arms per pair x {deposited} deposited pairs, plus 1 for the repeated REF arm "
                                    f"on the primary GNMT pair, which is in the deposit: "
                                    f"4 x {deposited} + 1 = {calls_cancelled}. The job's last progress line printed "
                                    f"{cancel['last_progress_line_calls']} calls at pair 10 of 1,306, and pair 11 "
                                    f"cost the remaining 4."),
            "previously_reported": {"calls": 86, "pairs": 21},
            "correction": (f"86 calls and 21 pairs were reported for this job; the log and the checkpoint give "
                           f"{calls_cancelled} calls and {deposited} pairs."),
        },
        "total_predict_sequence_calls": total,
        "total_previously_reported_about": 11222,
        "total_how_derived": (" + ".join([str(calls_cancelled)] + [str(j["calls_reported_at_DONE"]) for j in completed])
                             + f" = {total}"),
        "note_on_the_counter": ("the scorer increments its counter once per successful predict_sequence return, so a "
                               "request retried under the quota wrapper is counted once. The total is successful "
                               "returns, not HTTP requests."),
        "checkpoint_rows": len(tags),
        "checkpoint_expected_rows": "6 GNMT pairs + 100 chr6 draws + 199 scored observed pairs + 2,500 draws = 2,805",
        "rate_ceiling": {
            "v1_setting_seconds_between_call_starts": 2.1,
            "v1_implied_ceiling_per_min": round(60 / 2.1, 2),
            "v1_max_rate_printed": max(j["max_rate_per_min_printed"] or 0 for j in per_job),
            "v2_setting_seconds_between_call_starts": 2.15,
            "v2_implied_ceiling_per_min": round(60 / 2.15, 2),
            "note": "v1's 2.1 s spacing implies 28.6 calls/min and its logs printed up to 29.2/min, above the 28/min "
                    "ceiling this repair was asked to hold. v2 uses 2.15 s = 27.9/min."},
    }


def job_start(jobid: str) -> str:
    """sacct Start for one job, so a log's time-of-day gets its date from the scheduler."""
    import subprocess
    try:
        r = subprocess.run(["sacct", "-j", jobid, "-n", "-P", "-o", "Start", "-X"],
                           capture_output=True, text=True, timeout=60)
        return r.stdout.strip().splitlines()[0] if r.stdout.strip() else ""
    except Exception:  # noqa: BLE001 - sacct absent on a node is not a reason to fail the analysis
        return ""


def timestamps() -> dict:
    def mt(p: pathlib.Path) -> str:
        return dt.datetime.fromtimestamp(p.stat().st_mtime, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    rows = []
    for p in sorted((V1 / "prespec").glob("*.json")):
        d = json.loads(p.read_text())
        rows.append({"file": p.name, "mtime_utc": mt(p), "written_utc_field": d.get("written_utc"),
                     "written_utc_is_after_its_own_mtime": bool(d.get("written_utc") and d["written_utc"] > mt(p))})
    sha = V1 / "prespec" / "f2_00_prespec.sha256"
    # the first API call is the first "API key from" line of the EARLIEST job, so order by jobid, not by the
    # time-of-day string: the later jobs' 02:xx times sort before the first job's 23:xx.
    per_log = []
    for p in sorted((V1 / "logs").glob("agp_f2-haplotype-2*.out"),
                    key=lambda q: int(re.search(r"-(\d+)\.out$", q.name).group(1))):
        m = re.search(r"\[(\d\d:\d\d:\d\d)\] API key from", p.read_text())
        if m:
            per_log.append({"jobid": re.search(r"-(\d+)\.out$", p.name).group(1),
                            "first_api_call_time_of_day_utc": m.group(1),
                            "job_start_local": job_start(re.search(r"-(\d+)\.out$", p.name).group(1))})
    first_call = per_log[0]["first_api_call_time_of_day_utc"] if per_log else None
    order_ok = [r["written_utc_field"] for r in rows if r["file"].startswith("f2_00_amendment")]
    return {
        "v1_prespec_files": rows,
        "v1_sha256_file": {"file": sha.name, "mtime_utc": mt(sha)},
        "v1_first_api_call_per_job": per_log,
        "v1_first_api_call_utc_time_of_day": first_call,
        "v1_claim_in_RESULTS_md": "prespecification, sha256 stamped before the first API call",
        "correction": (f"the prespecification JSON's own mtime {mt(V1 / 'prespec' / 'f2_00_prespec.json')} does "
                       f"precede the first API call at {first_call} UTC on 2026-09-14, so the prespecification "
                       f"itself is pre-registered. The digest file that is supposed to evidence it was written at "
                       f"{mt(sha)}, after stage 1 finished, so the digest is not evidence of anything: it was "
                       f"computed over files that already existed at the end of the run. Only the JSON's mtime "
                       f"carries the pre-registration."),
        "written_utc_defect": ("every one of v1's seven amendment files carries a written_utc later than its own file "
                               "mtime, by 4 to 23 minutes, so the field is not a record of when the bytes were "
                               "written. In mtime order the amendments are 01, 02, 03, 04, 05, 06, 07, but amendment "
                               f"03's written_utc {order_ok[2]} is later than amendment 04's {order_ok[3]}, so the "
                               "recorded order contradicts the file order for 03 and 04."),
        "v2_practice": ("v2 writes written_utc from the machine clock in the same call that writes the bytes, and "
                        "rewrites the digest file in the same second, so neither field can post-date its file."),
    }


def main() -> None:
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    v1_rows = group_v1(load_jsonl(V1 / "raw" / "scored_rows.jsonl"))
    gnmt, gnull_c20 = v1_rows["gnmt"], v1_rows["gnmtnull"]
    obs, nul2500 = v1_rows["obs"], v1_rows["null"]
    nul1000 = [r for r in nul2500 if int(re.sub(r"\D", "", r["tag"]) or 0) < 1000]
    spans = R.gene_spans()
    global GNMT_REF_RNA_SUM
    prim = next(r for r in gnmt if r.get("is_primary"))
    GNMT_REF_RNA_SUM = float(prim["rna_ref_sum"])

    out: dict = {"package": "f2-haplotype-v2", "supersedes": V1.name,
                 "v1_row_counts": {k: len(v) for k, v in v1_rows.items()},
                 "gnmt_primary_ref_rna_sum": GNMT_REF_RNA_SUM}

    # ---- defect 1
    gs_path = OUT / "raw" / "scored_genespan_null.jsonl"
    if gs_path.exists():
        gs = load_jsonl(gs_path)
        # the scored rows carry the design fields the scorer copied; join the draw table for the rest
        drawn = {d["tag"]: d for d in csv.DictReader(open(OUT / "tables" / "gnmt_genespan_null_pairs.tsv"),
                                                     delimiter="\t")}
        for r in gs:
            for k, v in drawn.get(r["tag"], {}).items():
                r.setdefault(k, v)
        design = json.loads((OUT / "tables" / "gnmt_genespan_null_design.json").read_text())
        rep = f2_block(gnmt, gs, "gene-span-matched chr6 null (REPAIRED, primary)")
        old = f2_block(gnmt, gnull_c20, "central-20kb chr6 null (v1, SUPERSEDED for the RNA channel)")
        first_only = [r for r in gs if str(r.get("is_first_pair_from_this_gene", "")).lower() in ("true", "1")]
        sub_first = f2_block(gnmt, first_only, "one pair per host gene sub-null (sensitivity)")
        refs = np.array([float(r["rna_ref_sum"]) for r in gs], float)
        sig = [r for r in gs if 0.5 * GNMT_REF_RNA_SUM <= float(r["rna_ref_sum"]) <= 2.0 * GNMT_REF_RNA_SUM]
        sub_sig = (f2_block(gnmt, sig, "signal-matched sub-null, REF RNA sum within 2x of GNMT's (POST HOC)")
                   if len(sig) >= R.MIN_MATCHED_NULL else {"n_null_draws": len(sig), "testable": False,
                                                           "reason": f"only {len(sig)} draws within 2x of GNMT's REF RNA sum"})
        widths = np.array([int(r["rna_readout_bp"]) for r in gs], float)
        out["defect_1_repair"] = {
            "design": design,
            "repaired_null_readout": {
                "all_draws_use_the_gene_span": bool(all(r["rna_readout"] == "gene_span" for r in gs)),
                "rna_readout_bp": {"min": float(widths.min()), "median": float(np.median(widths)),
                                   "max": float(widths.max())},
                # the recipe's mask covers [span_start, span_end) of the 1-based GTF span, so end - start bases
                "gnmt_rna_readout_bp": GNMT_SPAN[1] - GNMT_SPAN[0],
                "gnmt_span_width_bp_1based_inclusive": GNMT_SPAN[1] - GNMT_SPAN[0] + 1,
                "v1_null_rna_readout_bp": 20000,
                "width_ratio_v1_null_over_observed": round(20000 / (GNMT_SPAN[1] - GNMT_SPAN[0]), 3),
                "width_ratio_repaired_null_over_observed": round(float(np.median(widths)) / (GNMT_SPAN[1] - GNMT_SPAN[0]), 3),
                "ref_rna_sum_vs_gnmt": {"gnmt": GNMT_REF_RNA_SUM, "null_min": float(refs.min()),
                                        "null_median": float(np.median(refs)), "null_max": float(refs.max()),
                                        "n_within_2x_of_gnmt": len(sig)}},
            "REPAIRED_primary": rep, "V1_SUPERSEDED": old,
            "sensitivity_one_pair_per_host_gene": sub_first,
            "POST_HOC_signal_matched": sub_sig,
            "moved": {ch: {"v1_p": old["F2.1_p_per_channel"][ch], "v2_p": rep["F2.1_p_per_channel"][ch],
                           "v1_inside": old["F2.1_verdict_per_channel"][ch],
                           "v2_inside": rep["F2.1_verdict_per_channel"][ch],
                           "v1_null_q95_abs_residual": old["null_q95_abs_residual"][ch],
                           "v2_null_q95_abs_residual": rep["null_q95_abs_residual"][ch],
                           "relative_change_in_null_q95_percent":
                               round(100.0 * (rep["null_q95_abs_residual"][ch] / old["null_q95_abs_residual"][ch] - 1), 2)}
                      for ch in CH}}
        write_tsv(OUT / "tables" / "f2v2_genespan_null_effects.tsv", gs)
        json.dump(out["defect_1_repair"], (OUT / "tables" / "f2v2_gnmt_verdict.json").open("w"), indent=1, default=float)
    else:
        out["defect_1_repair"] = {"status": "NOT RUN: raw/scored_genespan_null.jsonl absent"}

    # ---- loader equivalence (step 2b), if it has been run
    le = OUT / "tables" / "f2v2_loader_equivalence.json"
    out["loader_equivalence"] = json.loads(le.read_text()) if le.exists() else {"status": "NOT RUN"}

    # ---- defect 2
    d2 = {"what_v1_did": ("channel_verdict computed the empirical p once against the full pooled null and then masked "
                          "it for the both-active stratum, so the stratum's p-values were against 2,500 pooled draws "
                          "rather than the draws that meet the stratum. The stratum selects pairs whose two arms are "
                          "both above the null's 95th percentile of |arm|, and the residual is a function of those "
                          "arms, so this is the arm-magnitude confound v1's own amendment 06 exists to avoid."),
          "n2500": both_active_block(obs, nul2500, "2,500 draws (v1's published stage)"),
          "n1000": both_active_block(obs, nul1000, "1,000 draws (v1's stage 1)")}
    out["defect_2_repair"] = d2
    json.dump(d2, (OUT / "tables" / "f2v2_both_active_stratum.json").open("w"), indent=1, default=float)

    # ---- numbers carried forward unchanged, so v2 is self-contained once v1 is marked superseded
    v1s = json.loads((V1 / "tables" / "f3_summary_n2500.json").read_text())
    v1gv = json.loads((V1 / "tables" / "f2_gnmt_verdict_n2500.json").read_text())
    v1prim = next(e for e in v1gv["pairs"] if e.get("is_primary"))
    out["carried_forward_unchanged_from_v1"] = {
        "source": str(V1 / "tables/f3_summary_n2500.json"),
        "n_obs": v1s["n_obs"], "n_null": v1s["n_null"],
        "per_channel": {ch: {k: v1s["f3_corrected"][ch][k] for k in
                            ("median_abs_residual", "null_median_abs_residual", "null_q95_abs_residual",
                             "n_pairs_exceeding_null_q95_UNCORRECTED", "min_p_empirical",
                             "median_abs_single_variant_effect", "null_median_abs_single_variant_effect")}
                       | {"n_pairs_bh": v1s["f3_corrected"][ch]["all_pairs"]["n_pairs_bh"],
                          "n_signals_bh": v1s["f3_corrected"][ch]["all_pairs"]["n_signals_bh"]}
                       for ch in CH},
        "family_level_contrast": v1s["family_level_contrast"],
        "magnitude_matched": v1s["magnitude_matched_diagnostic_POST_HOC"]["per_channel"],
        "separation_stratified": {k: {"n_observed": v["n_observed"], "n_null": v["n_null"],
                                      "per_channel": {ch: {"n_pairs_bh": v["per_channel"][ch]["all_pairs"].get("n_pairs_bh"),
                                                           "min_p": v["per_channel"][ch]["all_pairs"].get("min_p_empirical")}
                                                      for ch in CH}}
                                  for k, v in v1s["separation_stratified_POST_HOC"]["strata"].items()},
        "p_floor_arithmetic": v1s["p_floor_arithmetic"],
        "determinism": v1s["determinism_check_gnmt_primary_ref_arm"],
        "rna_readout_counts": v1s["rna_readout_counts"],
        "gnmt_primary_arms": {ch: {"v1_log2": v1prim[f"{ch}_v1_log2"], "v2_log2": v1prim[f"{ch}_v2_log2"],
                                   "joint_log2": v1prim[f"{ch}_joint_log2"],
                                   "additive_expectation_log2": v1prim[f"{ch}_additive_expectation_log2"],
                                   "ref_sum": v1prim[f"{ch}_ref_sum"]} for ch in CH},
        "blocks": v1s["blocks"],
    }

    # ---- readout geometry and accounting
    # the BH flags v1 published live in its analysed TSV, not in the checkpoint
    obs_tsv = list(csv.DictReader(open(V1 / "tables" / "f3_observed_effects_n2500.tsv"), delimiter="\t"))
    geom = readout_geometry(obs_tsv, spans)
    json.dump(geom, (OUT / "tables" / "f2v2_readout_geometry.json").open("w"), indent=1, default=float)
    out["readout_geometry"] = geom
    acc = {"calls_and_checkpoints": accounting(), "prespec_timestamps": timestamps()}
    json.dump(acc, (OUT / "tables" / "f2v2_accounting.json").open("w"), indent=1, default=float)
    out["accounting"] = acc

    json.dump(out, (OUT / "tables" / "f2v2_summary.json").open("w"), indent=1, default=float)
    print(json.dumps({"f2.1_v2": out["defect_1_repair"].get("REPAIRED_primary", {}).get("F2.1_verdict_per_channel"),
                      "f2.1_v1": out["defect_1_repair"].get("V1_SUPERSEDED", {}).get("F2.1_verdict_per_channel"),
                      "both_active_repaired": {ch: d2["n2500"]["per_channel"][ch]["REPAIRED_p_against_the_null_inside_the_stratum"]["n_bh_rejected_if_the_20_draw_gate_is_ignored"] for ch in CH},
                      "both_active_v1": {ch: d2["n2500"]["per_channel"][ch]["V1_AS_PUBLISHED_p_against_the_pooled_null"]["n_bh_rejected"] for ch in CH},
                      "total_calls": acc["calls_and_checkpoints"]["total_predict_sequence_calls"]}, indent=1))


if __name__ == "__main__":
    main()
