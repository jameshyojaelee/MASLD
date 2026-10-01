#!/usr/bin/env python3
"""AD.1: does local zero-shot scoring reproduce the archived hosted Atlas accessibility ranking?

Prespecified bar (ADAPTATION_ARM_PRESPEC.md, prediction AD.1): Spearman >= 0.9 between the local and
the archived score sets on the 3,845 matched Currin leads. If it is below 0.9 the arm stops and
reports what differs.

The diagnostics below are fixed here, before any local score was read, together with a written
prediction of each one's result. Predicting only the decision rule and not the diagnostics is how a
failure gets explained after the fact instead of before it.

  D1 length          which of the five SDK lengths agrees best with the archive.
                     PREDICTED: 1,048,576. The runtime probe's HSD17B13 anchor reproduced the hosted
                     model API's ATAC and DNase log2 to four decimal places at 1 Mb, and 1 Mb is the
                     longest length the released weights were trained at.
  D2 gate            ATAC Spearman on the full matched set, with a 1-Mb-block bootstrap interval.
                     PREDICTED: >= 0.9, i.e. AD.1 holds.
  D3 channel parity  DNase Spearman is within 0.05 of ATAC's.
                     PREDICTED: held. Both are the same scorer at the same width on adjacent assays.
  D4 effect size     sign agreement by quartile of |archived score|. PREDICTED: the bottom quartile
                     is below 0.90 and the top quartile is above 0.95. The runtime probe measured
                     only 5.6x concentration of the ATAC |alt-ref| mass within +/-2 kb of a
                     substitution, so a 501-bp mask on a near-null variant is dominated by the
                     diffuse component and two implementations will disagree in sign there.
  D5 scale           median ratio local/archived and median |local - archived| against the median
                     |archived|. PREDICTED: ratio within 0.9 to 1.1; a high correlation with a ratio
                     far from 1 would mean a scale convention differs, not the model.
  D6 per track       the three liver ATAC tracks correlate individually at least as well as the
                     3-track mean. PREDICTED: held, within 0.03 of the mean's value.
  D7 rejected        windows rejected for a non-ACGT byte or a chromosome edge. PREDICTED: between
                     1 and 60 of 3,845. The label build only ever required a clean 4,096-bp window,
                     so a 1-Mb window can still reach an assembly gap or a chromosome end.
  D8 orientation     the archived column's own Spearman against beta_alt must come out positive and
                     near the 0.7361 that `c1-endpoint2-v2` printed for tier A4. This is the guard
                     that no sign was applied twice; it is not a new result.

Every prediction is reported beside its observed value, including the ones that fail.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import i1_common as C  # noqa: E402

GATE_BAR = 0.9
TIER_A4_ARCHIVED_ATAC_VS_BETA = 0.7361  # c1-endpoint2-v2 RESULTS.md, tier A4 (3,845)


def channel_block(m: pd.DataFrame, chan: str, resamples: int) -> dict:
    loc = m[f"local_{chan}_liver"].to_numpy(float)
    ar = m[f"archived_{chan}_liver"].to_numpy(float)
    blocks = m["block_1mb"].to_numpy(str)
    ok = np.isfinite(loc) & np.isfinite(ar)
    loc, ar, blocks = loc[ok], ar[ok], blocks[ok]

    rho = C.fast_spearman(ar, loc)
    draws = C.block_bootstrap_plain(ar, {"local": loc}, blocks, resamples=resamples)["local"]
    low, high, se = C.interval_of(draws)
    rec = {
        "channel": chan,
        "n": int(ok.sum()),
        "blocks": int(len(set(blocks.tolist()))),
        "spearman": rho,
        "ci_low": low,
        "ci_high": high,
        "bootstrap_se": se,
        "bootstrap_resamples": int(resamples),
        "bootstrap_seed": C.BOOTSTRAP_SEED,
        "bootstrap_unit": "1Mb_block",
        "pearson": float(np.corrcoef(ar, loc)[0, 1]),
        "sign_agreement": float((np.sign(ar) == np.sign(loc)).mean()),
        "median_abs_archived": float(np.median(np.abs(ar))),
        "median_abs_local": float(np.median(np.abs(loc))),
        "median_abs_difference": float(np.median(np.abs(loc - ar))),
        "max_abs_difference": float(np.max(np.abs(loc - ar))),
        "median_ratio_local_over_archived": float(np.median(loc / ar)),
        "passes_bar": bool(rho >= GATE_BAR),
        "ci_low_above_bar": bool(low >= GATE_BAR),
    }
    # D4: agreement by effect-size quartile of |archived|.
    q = np.quantile(np.abs(ar), [0.25, 0.5, 0.75])
    edges = [-np.inf, *q, np.inf]
    strata = []
    for i in range(4):
        sel = (np.abs(ar) > edges[i]) & (np.abs(ar) <= edges[i + 1])
        if sel.sum() < 3:
            continue
        strata.append(
            {
                "quartile_of_abs_archived": i + 1,
                "n": int(sel.sum()),
                "abs_archived_range": [float(np.abs(ar)[sel].min()), float(np.abs(ar)[sel].max())],
                "spearman": C.fast_spearman(ar[sel], loc[sel]),
                "sign_agreement": float((np.sign(ar[sel]) == np.sign(loc[sel])).mean()),
            }
        )
    rec["by_effect_size_quartile"] = strata
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--scores", required=True, help="the full-gate local score TSV")
    ap.add_argument("--resamples", type=int, default=C.RESAMPLES)
    args = ap.parse_args()

    run = pathlib.Path(args.run)
    gate = pd.read_csv(run / "inputs/gate_variants.tsv", sep="\t")
    local = pd.read_csv(args.scores, sep="\t")

    cards = sorted(set(local.card_tag.astype(str)))
    if len(cards) != 1:
        raise C.ContractError(
            f"local scores span {cards}; scores are not reproducible across GPU architectures and a "
            "mixed-card table must be rescored, not merged"
        )
    lengths = sorted(set(local.length_bp.astype(int)))
    if len(lengths) != 1:
        raise C.ContractError(f"local scores span input lengths {lengths}")

    arch_cols = [c for c in gate.columns if c.startswith("archived_")]
    m = local.merge(
        gate[["key", "beta_alt", "heldout_fold", *arch_cols]].rename(
            columns={"beta_alt": "beta_alt_gate", "heldout_fold": "fold_gate"}
        ),
        on="key",
        how="inner",
        validate="one_to_one",
    )
    C.log(f"{len(m)} leads with both a local and an archived score, card {cards[0]}, {lengths[0]} bp")

    out: dict = {
        "card": str(local.card.iloc[0]),
        "card_tag": cards[0],
        "length_bp": lengths[0],
        "mask_width_bp": 501,
        "aggregation": "DIFF_LOG2_SUM",
        "gate_bar_spearman": GATE_BAR,
        "gate_set_rows": int(len(gate)),
        "gate_set_blocks": int(gate.block_1mb.nunique()),
        "local_rows_scored": int(len(local)),
        "paired_rows": int(len(m)),
        "paired_blocks": int(m.block_1mb.nunique()),
    }

    rej = pathlib.Path(args.scores)
    rej = rej.with_name(rej.stem + "_rejected.tsv")
    if rej.exists():
        rd = pd.read_csv(rej, sep="\t")
        out["rejected_windows"] = int(len(rd))
        out["rejected_reasons"] = {k: int(v) for k, v in rd.reason.value_counts().items()}
        out["rejected_keys"] = sorted(rd.key.astype(str).tolist())[:400]
        # Is the exclusion at random with respect to effect size? Dropping near-null leads would
        # make the gate easier, because near-null leads are where two implementations most disagree
        # in sign, so the direction of any skew has to be on the page.
        gi = gate.set_index("key")
        rej_arch = gi.loc[rd.key.astype(str), "archived_atac_liver"].to_numpy(float)
        kept_arch = gi.loc[
            [k for k in gi.index if k not in set(rd.key.astype(str))], "archived_atac_liver"
        ].to_numpy(float)
        out["rejected_vs_kept_effect_size"] = {
            "rejected_n": int(len(rej_arch)),
            "kept_n": int(len(kept_arch)),
            "rejected_median_abs_archived_atac": float(np.median(np.abs(rej_arch))),
            "kept_median_abs_archived_atac": float(np.median(np.abs(kept_arch))),
            "note": (
                "a rejected window is not scored and not replaced, so these leads leave the gate "
                "set; if their archived effects are systematically smaller the remaining set is "
                "easier than the full one"
            ),
        }
        # Would a shorter window have been clean at the same variant? This separates a genome
        # defect from a consequence of choosing a 1-Mb window.
        import pysam

        fa = pysam.FastaFile(C.FASTA_PATH)
        shorter = {}
        for L in (524288, 131072, 16384):
            clean = 0
            for r in rd.itertuples(index=False):
                s0, e0 = C.window_bounds(int(r.pos_hg38), L)
                chk = C.validate_window(fa, str(r.chr), s0, e0, int(r.pos_hg38), str(r.ref))
                clean += int(chk["acgt_ok"] and chk["reference_match"])
            shorter[L] = clean
        out["rejected_would_be_clean_at_shorter_length"] = shorter
    else:
        out["rejected_windows"] = 0
        out["rejected_reasons"] = {}
        out["rejected_keys"] = []
        out["rejected_vs_kept_effect_size"] = {}
        out["rejected_would_be_clean_at_shorter_length"] = {}

    channels = {}
    for chan in ("atac", "dnase"):
        channels[chan] = channel_block(m, chan, args.resamples)
    out["channels"] = channels

    # D6: per liver ATAC / DNase track.
    per_track = []
    for chan, tracks in (("atac", C.AG_ATAC_TRACKS), ("dnase", C.AG_DNASE_TRACKS)):
        for t in tracks:
            tag = t.split()[0]
            lc, ac = f"local_{chan}__{tag}", f"archived_{chan}__{tag}"
            if lc not in m or ac not in m:
                continue
            a, b = m[ac].to_numpy(float), m[lc].to_numpy(float)
            ok = np.isfinite(a) & np.isfinite(b)
            per_track.append(
                {
                    "channel": chan,
                    "track": t,
                    "n": int(ok.sum()),
                    "spearman": C.fast_spearman(a[ok], b[ok]),
                    "sign_agreement": float((np.sign(a[ok]) == np.sign(b[ok])).mean()),
                    "median_abs_difference": float(np.median(np.abs(a[ok] - b[ok]))),
                }
            )
    out["per_track"] = per_track

    # D8: the orientation guard. The archived column's own signed Spearman against beta_alt on this
    # same set must reproduce what c1-endpoint2-v2 printed; if it does not, a sign moved.
    y = m["beta_alt"].to_numpy(float)
    ar = m["archived_atac_liver"].to_numpy(float)
    lo = m["local_atac_liver"].to_numpy(float)
    ok = np.isfinite(y) & np.isfinite(ar) & np.isfinite(lo)
    out["orientation_guard"] = {
        "archived_atac_vs_beta_alt_spearman_here": C.fast_spearman(y[ok], ar[ok]),
        "archived_atac_vs_beta_alt_spearman_c1e2_tierA4": TIER_A4_ARCHIVED_ATAC_VS_BETA,
        "local_atac_vs_beta_alt_spearman_here": C.fast_spearman(y[ok], lo[ok]),
        "n": int(ok.sum()),
        "note": (
            "pooled over all folds, so it is comparable with the c1-endpoint2-v2 tier A4 value and "
            "is NOT the fold-macro statistic the zero-shot endpoint reports"
        ),
        "beta_alt_from_gate_table_matches_score_table": bool(
            np.allclose(
                m["beta_alt"].to_numpy(float),
                m["beta_alt_gate"].to_numpy(float),
                equal_nan=True,
            )
        ),
    }

    # Cross-process reproducibility, which the runtime probe left open: it established that repeat
    # calls WITHIN one process on one card are bit-identical and said plainly that reproducibility
    # across separate processes was not tested. The gate rescored, in its own process, the same 60
    # leads the pilot sweep scored at the same length on the same card, so the check is free.
    sp = run / f"raw/sweep_{lengths[0]}.tsv"
    if sp.exists():
        sd = pd.read_csv(sp, sep="\t")
        both = local.merge(sd, on="key", how="inner", suffixes=("_gate", "_sweep"))
        if len(both):
            rec = {"n": int(len(both)), "length_bp": int(lengths[0])}
            for chan in ("atac", "dnase"):
                a = both[f"local_{chan}_liver_gate"].to_numpy(float)
                b = both[f"local_{chan}_liver_sweep"].to_numpy(float)
                rec[f"{chan}_max_abs_difference"] = float(np.max(np.abs(a - b)))
                rec[f"{chan}_n_bit_identical"] = int(np.sum(a == b))
            rec["same_card"] = bool(
                set(sd.card_tag.astype(str)) == set(local.card_tag.astype(str))
            )
            rec["note"] = (
                "two separate processes, same card, same weights, same sequences; the runtime probe "
                "tested only within-process determinism"
            )
            out["cross_process_reproducibility"] = rec

    sweep_p = run / "tables/pilot_length_sweep.tsv"
    if sweep_p.exists():
        out["pilot_length_sweep"] = pd.read_csv(sweep_p, sep="\t").to_dict(orient="records")
    pv = run / "tables/pilot_verdict.json"
    if pv.exists():
        out["pilot_verdict"] = json.loads(pv.read_text())

    atac = channels["atac"]
    dnase = channels["dnase"]
    q = {s["quartile_of_abs_archived"]: s for s in atac["by_effect_size_quartile"]}
    preds = [
        {
            # Judged on the matched table, which compares lengths on the leads every length scored.
            # The unmatched table compares different n per length and is not a length comparison.
            "id": "D1",
            "prediction": "1,048,576 bp agrees best with the archive",
            "observed": (out.get("pilot_verdict") or {}).get("best_length_on_matched_leads")
            or (out.get("pilot_verdict") or {}).get("chosen_length_bp"),
            "met": (
                (out.get("pilot_verdict") or {}).get("best_length_on_matched_leads")
                or (out.get("pilot_verdict") or {}).get("chosen_length_bp")
            )
            == 1_048_576,
        },
        {
            "id": "D2 (AD.1)",
            "prediction": f"ATAC Spearman >= {GATE_BAR}",
            "observed": atac["spearman"],
            "met": atac["passes_bar"],
        },
        {
            "id": "D3",
            "prediction": "DNase Spearman within 0.05 of ATAC's",
            "observed": abs(dnase["spearman"] - atac["spearman"]),
            "met": bool(abs(dnase["spearman"] - atac["spearman"]) <= 0.05),
        },
        {
            "id": "D4",
            "prediction": "sign agreement below 0.90 in the bottom |archived| quartile and above 0.95 in the top",
            "observed": {
                "q1_sign_agreement": q.get(1, {}).get("sign_agreement"),
                "q4_sign_agreement": q.get(4, {}).get("sign_agreement"),
            },
            "met": bool(
                q.get(1, {}).get("sign_agreement", 1.0) < 0.90
                and q.get(4, {}).get("sign_agreement", 0.0) > 0.95
            ),
        },
        {
            "id": "D5",
            "prediction": "median local/archived ratio between 0.9 and 1.1",
            "observed": atac["median_ratio_local_over_archived"],
            "met": bool(0.9 <= atac["median_ratio_local_over_archived"] <= 1.1),
        },
        {
            "id": "D6",
            "prediction": "every liver ATAC track's own Spearman within 0.03 of the 3-track mean's",
            "observed": [
                {"track": p["track"], "spearman": p["spearman"]}
                for p in per_track
                if p["channel"] == "atac"
            ],
            "met": bool(
                all(
                    abs(p["spearman"] - atac["spearman"]) <= 0.03
                    for p in per_track
                    if p["channel"] == "atac"
                )
            ),
        },
        {
            "id": "D7",
            "prediction": "between 1 and 60 of 3,845 windows rejected",
            "observed": out["rejected_windows"],
            "met": bool(1 <= out["rejected_windows"] <= 60),
        },
        {
            "id": "D8",
            "prediction": (
                "the archived column's own Spearman against beta_alt here reproduces the "
                f"c1-endpoint2-v2 tier A4 value {TIER_A4_ARCHIVED_ATAC_VS_BETA} to within 0.01"
            ),
            "observed": out["orientation_guard"]["archived_atac_vs_beta_alt_spearman_here"],
            "met": bool(
                abs(
                    out["orientation_guard"]["archived_atac_vs_beta_alt_spearman_here"]
                    - TIER_A4_ARCHIVED_ATAC_VS_BETA
                )
                <= 0.01
            ),
        },
    ]
    out["predictions"] = preds
    out["gate_verdict"] = "PASS" if atac["passes_bar"] else "FAIL"
    out["gate_verdict_statement"] = (
        f"AD.1 {'holds' if atac['passes_bar'] else 'fails'}: local-versus-archived ATAC Spearman "
        f"{atac['spearman']:.4f} [{atac['ci_low']:.4f}, {atac['ci_high']:.4f}] on {atac['n']} leads "
        f"in {atac['blocks']} 1-Mb blocks, bar {GATE_BAR}"
    )

    C.write_json(run / "tables/ad1_gate.json", out)
    pd.DataFrame([channels["atac"], channels["dnase"]]).drop(
        columns=["by_effect_size_quartile"]
    ).to_csv(run / "tables/ad1_gate_channels.tsv", sep="\t", index=False)
    pd.DataFrame(per_track).to_csv(run / "tables/ad1_gate_per_track.tsv", sep="\t", index=False)
    pd.DataFrame(preds).to_csv(run / "tables/ad1_predictions.tsv", sep="\t", index=False)
    m.to_csv(run / "tables/ad1_paired_scores.tsv.gz", sep="\t", index=False)

    print(out["gate_verdict_statement"])
    print(json.dumps({k: v for k, v in out.items() if k != "pilot_length_sweep"}, indent=1)[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
