#!/usr/bin/env python3
"""Step 66 (P6d follow-up): read the dossier against the step-65 prespecification's written predictions.

Separate from step 65 so that the dossier is written before any prediction is scored, and so that the
block-aware statistics here never end up inside the exported table. Two of the five predictions need a unit
of inference the dossier itself does not impose: 204 rows carrying both a measured and a predicted liver RNA
direction sit in only 19 independent 1-Mb blocks, and a single locus contributes 45 of them.

The peak-overlap comparison here is POST HOC. The prespecified R2 was a calibration check with an interval
that turned out to be wrong by an order of magnitude; the within-locus high- vs low-weight contrast is what
that failure pointed at, and it is labelled as post hoc wherever it is reported.

Outputs (tables/): dossier_prediction_verdicts.json
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
PRESPEC = SCRIPT_DIR / "65_clinical_dossier_prespec.json"
ROOT = la.out_root()
TABLES = ROOT / "tables"
TRACK0 = la.track0_root() / "tables"
WEIGHT_FLOOR = 0.01


def _load_63():
    import importlib.util
    spec = importlib.util.spec_from_file_location("indel_rescue", SCRIPT_DIR / "63_indel_rescue.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


block_sign_test = _load_63().block_sign_test


def concordance_rows(rows) -> list:
    """One entry per dossier row that carries both directions, as a signed deviation from chance.

    The sign test then takes the median within each 1-Mb block, so a locus with 45 rows counts once.
    """
    out = []
    for r in rows:
        c = str(r.get("direction_concordant", ""))
        if c not in ("True", "False"):
            continue
        out.append({"analysis_block": r["analysis_block"], "value": 0.5 if c == "True" else -0.5})
    return out


def overlap_contrast_rows(per_block: dict) -> list:
    """Within-locus difference in measured-peak overlap between high- and low-weight variants.

    A block missing either arm is dropped rather than scored zero: it carries no contrast, and scoring it
    zero would shrink the estimate toward no difference with blocks that never tested anything.
    """
    out = []
    for blk, arms in sorted(per_block.items()):
        hk, hn = arms["high"]
        lk, ln = arms["low"]
        if hn <= 0 or ln <= 0:
            continue
        out.append({"analysis_block": blk, "value": hk / hn - lk / ln})
    return out


def verdict(value, lo: float, hi: float) -> str:
    if value is None:
        return "not_evaluable"
    v = float(value)
    if v < lo:
        return "not_met_below"
    if v > hi:
        return "not_met_above"
    return "met"


def main() -> None:
    prespec = json.load(PRESPEC.open())
    rows = la.read_tsv(TABLES / "clinical_variant_dossier.tsv")
    by_signal = defaultdict(list)
    for r in rows:
        by_signal[r["signal_uid"]].append(r)

    # R1
    signals = {r["signal_uid"] for r in la.read_tsv(TRACK0 / "eligible_signals.tsv")
               if r["universe"] in ("A_direct", "B_direct")}
    r1 = {"n_signals_expected": len(signals), "n_signals_in_dossier": len(by_signal),
          "verdict": "met" if len(by_signal) == len(signals) else "not_met"}

    # R3: does the export itself show the hazard, or does only the coverage column carry it?
    non_covered = [s for s, rs in by_signal.items() if rs[0]["signal_atlas_coverage_state"] != "covered"]
    hidden = [s for s in non_covered if not any(r["atlas_served"] == "False" for r in by_signal[s])]
    r3 = {"n_non_covered_signals": len(non_covered), "n_with_no_unserved_row_in_the_export": len(hidden),
          "hidden_signals": sorted(hidden),
          "hidden_indel_deferred_shares": sorted(float(by_signal[s][0]["signal_excluded_share_indel_deferred"]) for s in hidden),
          "verdict": "met" if not hidden else "not_met",
          "what_it_means": ("the 0.01 inclusion floor does not by itself expose every signal that loses mass to "
                            "deferred indels; the per-row coverage-state column, not the row set, is what gates a reader")}

    # R4
    top = [r for r in rows if int(r["weight_rank"]) == 1]
    unserved_top = sorted({r["signal_uid"] for r in top if r["atlas_served"] == "False"})
    r4 = {"n_signals_top_variant_unserved": len(unserved_top), "threshold": 10,
          "verdict": "met" if len(unserved_top) >= 10 else "not_met", "signals": unserved_top}

    # R5, with the block as the unit and the marginal sign skew as the baseline
    scored = [r for r in rows if r["direction_concordant"] in ("True", "False")]
    obs = sum(r["direction_concordant"] == "True" for r in scored) / len(scored) if scored else None
    base = la.marginal_expected_concordance([float(r["predicted_direction"]) for r in scored],
                                            [float(r["measured_direction"]) for r in scored]) if scored else None
    bst = block_sign_test(concordance_rows(rows)) if scored else {}
    biggest = max(((b, sum(1 for r in scored if r["analysis_block"] == b))
                   for b in {r["analysis_block"] for r in scored}), key=lambda kv: kv[1], default=("", 0))
    r5 = {"observed_concordance": obs, "marginal_expected_baseline": base, "n_rows": len(scored),
          "block_sign_test": bst, "largest_block": {"analysis_block": biggest[0], "n_rows": biggest[1]},
          "verdict": verdict(obs, 0.45, 0.65),
          "what_it_means": ("inside the prespecified interval and consistent with the prior in-house AlphaGenome "
                            "directional benchmark (auROC 0.5607, n 211); the block sign test rests on very few "
                            "independent loci and is not a directional claim")}

    # R2 as written, then the post-hoc within-locus contrast it pointed at
    dossier_uids = {r["variant_uid"] for r in rows if r["variant_uid"]}
    in_peak, queried_w, blk_of = set(), {}, {}
    with la.open_text(TRACK0 / "variant_peak_overlaps.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["universe"] in ("A_direct", "B_direct") and r["peak"]:
                in_peak.add(r["variant_uid"])
    sig_block = {r["signal_uid"]: r["analysis_block"] for r in la.read_tsv(TRACK0 / "eligible_signals.tsv")}
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["signal_uid"] in signals and r["in_query_set"] == "True":
                u = r["variant_uid"]
                queried_w[u] = max(queried_w.get(u, 0.0), float(r["weight"]))
                blk_of[u] = sig_block[r["signal_uid"]]
    per_block = defaultdict(lambda: {"high": [0, 0], "low": [0, 0]})
    for u, w in queried_w.items():
        arm = "high" if w >= WEIGHT_FLOOR else "low"
        per_block[blk_of[u]][arm][1] += 1
        if u in in_peak:
            per_block[blk_of[u]][arm][0] += 1
    hi_k = sum(1 for u, w in queried_w.items() if w >= WEIGHT_FLOOR and u in in_peak)
    hi_n = sum(1 for w in queried_w.values() if w >= WEIGHT_FLOOR)
    lo_k = sum(1 for u, w in queried_w.items() if w < WEIGHT_FLOOR and u in in_peak)
    lo_n = sum(1 for w in queried_w.values() if w < WEIGHT_FLOOR)
    share_rows = sum(1 for r in rows if r["measured_peak"]) / len(rows)
    contrast = overlap_contrast_rows({k: {"high": tuple(v["high"]), "low": tuple(v["low"])}
                                      for k, v in per_block.items()})
    r2 = {"prespecified_interval": [0.20, 0.60],
          "observed_share_of_rows_in_a_measured_peak": share_rows,
          "verdict": verdict(share_rows, 0.20, 0.60),
          "why_the_prediction_was_wrong": ("the interval was anchored on nothing: the measured snATAC consensus "
                                           "peaks cover on the order of 1% of the genome, so a fine-mapped variant "
                                           "landing in one is the exception, not the rule"),
          "POST_HOC_within_locus_contrast": {
              "high_weight_variants_in_a_peak": [hi_k, hi_n, hi_k / hi_n if hi_n else None],
              "low_weight_variants_in_a_peak": [lo_k, lo_n, lo_k / lo_n if lo_n else None],
              "block_sign_test": block_sign_test(contrast) if contrast else {},
              "label": "POST_HOC: not prespecified, reported as an observation and not as a test result",
              "dossier_variants_only": len(dossier_uids)}}

    out = {"prespec_sha256": la.sha256_file(PRESPEC), "dossier_rows": len(rows),
           "R1_every_direct_signal_present": r1, "R2_measured_peak_share": r2,
           "R3_hazard_visible_in_the_export": r3, "R4_top_variant_unserved": r4,
           "R5_direction_concordance": r5,
           "claim_boundary": prespec["claim_boundary"]}
    json.dump(out, (TABLES / "dossier_prediction_verdicts.json").open("w"), indent=1, default=float)
    la.log(f"P6d checks: R1 {r1['verdict']}, R2 {r2['verdict']}, R3 {r3['verdict']}, "
           f"R4 {r4['verdict']}, R5 {r5['verdict']}")


if __name__ == "__main__":
    main()
