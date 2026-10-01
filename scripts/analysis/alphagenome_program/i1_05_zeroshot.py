#!/usr/bin/env python3
"""Zero-shot local AlphaGenome on the endogenous accessibility allele-effect endpoint.

Zero-shot means no training: the released weights score REF and ALT and the resulting signed number
is ranked against the measured caQTL `beta_alt` directly. It is the reference point arms 2 and 3 have
to beat, so the number that matters is not this arm's absolute Spearman but its paired difference
against the fitted comparators on the variants all arms cover.

Statistic, identical in construction to `c2-endogenous-head-20260914T194000Z`:
  - signed Spearman inside each of the 5 ChromBPNet chromosome-group held-out folds,
  - Fisher-z macro over the 5 folds (`fisher_macro`),
  - 10,000 bootstrap resamples of 1-Mb blocks WITHIN each fold, seed 20260914,
  - every arm recomputed inside the same draw, so arm differences are paired.

Two things this script insists on rather than assumes.

1. **The comparator is regenerated, not quoted.** The published HyenaDNA delta-ridge macro is
   0.1736479364554445 on all 32,322 leads. This loader recomputes it from that package's own
   `oof_predictions.tsv.gz` and fails if it does not reproduce it. A relayed number is a claim until
   it is found in its producing file.

2. **The scored set is a doubly selected subset, so the headline contrast is the PAIRED one
   recomputed there.** Two selections stack. The first is the Atlas: only 3,845 of the 32,322 Currin
   leads carry an archived hosted score, and `c1-endpoint2-v2` already recorded that the same
   ChromBPNet reaches 0.4495 on all 32,322 against 0.6450 on the 3,845, so an Atlas-served subset is
   an easier set for every model. The second is this arm's own window guard: a 1-Mb window is
   rejected outright when it carries an assembly-gap N or runs off a chromosome end, which the
   4,096-bp windows the deposited comparators were built on almost never did, and the leads that
   survive do not have the same archived-effect distribution as the leads that do not.

   So the published 0.1737 and 0.0435 are values on a different population and are never the
   contrast. Every comparator is recomputed from its deposited out-of-fold predictions on exactly
   the surviving variant set, with the same blocks and the same folds, and every difference is taken
   inside one bootstrap draw. The published values are printed beside the recomputed ones so the
   size of the shift is visible, and if the shift is material this script says so in `comparator_
   shift` rather than leaving a reader to infer it. A shift is a fact about the subset, not an
   error in either run.

No sign is applied anywhere. `beta_alt` is the ALT-dosage slope, positive meaning ALT increases
accessibility, and the local score is log2(1 + sum ALT) - log2(1 + sum REF) over the 501-bp mask,
positive meaning ALT increases predicted accessibility. Those already point the same way; applying a
sign here would cancel for half the family.
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

HYENA = "hyenadna__delta_ridge__ensemble"
CONTROL = "available_simple_controls__allele_identity_ridge__ensemble"
CADUCEUS = "caduceus__delta_ridge__ensemble"

PUBLISHED = {
    HYENA: {
        "macro": 0.1736479364554445,
        "ci": [0.16326892437614288, 0.18435287486797755],
        "se": 0.005362261642359096,
    },
    CONTROL: {
        "macro": 0.04349628551774106,
        "ci": [0.03248318878325781, 0.0543378172005105],
        "se": 0.00560903441366127,
    },
    CADUCEUS: {
        "macro": 0.12760125075118817,
        "ci": [0.11700213192842923, 0.13859091786256023],
        "se": 0.005505725817392921,
    },
}
ARM_NAME = {
    HYENA: "hyenadna_delta_ridge",
    CADUCEUS: "caduceus_delta_ridge",
    CONTROL: "control_allele_identity_ridge",
}

# The reproduction tolerance is not uniform, and the reason is measured rather than assumed.
# HyenaDNA and Caduceus predictions are continuous: all 32,322 values are distinct, and this loader
# reproduces their published macro at 0.0 and 2.8e-17. The allele-identity control's prediction is a
# ridge on a 16-level one-hot of (REF, ALT), so it takes only 12 to 13 distinct values in each fold
# of roughly 6,000 variants. With that many ties a last-bit difference in the deposited decimal moves
# which rows are tied and `rankdata`'s average ranks shift with it, which is why the control
# reproduces at 1.3e-05 rather than at machine precision. That is 0.24 percent of the control's own
# published bootstrap standard error of 0.0056, so it is a deposit-precision difference in a
# tie-degenerate arm, not a disagreement about the statistic. It is reported, not absorbed.
TOLERANCE = {HYENA: 1e-12, CADUCEUS: 1e-12, CONTROL: 1e-3}


def summarise(name, macro_draws, pooled_draws, point_macro, point_pooled, n, blocks) -> dict:
    lo, hi, se = C.interval_of(macro_draws)
    plo, phi, pse = C.interval_of(pooled_draws)
    return {
        "arm": name,
        "variants": int(n),
        "blocks": int(blocks),
        "macro_spearman": point_macro,
        "macro_ci_low": lo,
        "macro_ci_high": hi,
        "macro_bootstrap_se": se,
        "pooled_spearman": point_pooled,
        "pooled_ci_low": plo,
        "pooled_ci_high": phi,
        "pooled_bootstrap_se": pse,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--resamples", type=int, default=C.RESAMPLES)
    args = ap.parse_args()

    run = pathlib.Path(args.run)
    local = pd.read_csv(args.scores, sep="\t")
    cards = sorted(set(local.card_tag.astype(str)))
    if len(cards) != 1:
        raise C.ContractError(f"local scores span cards {cards}")
    lengths = sorted(set(local.length_bp.astype(int)))
    if len(lengths) != 1:
        raise C.ContractError(f"local scores span lengths {lengths}")

    labels = pd.read_csv(C.C2_LABELS, sep="\t").sort_values("row_index").reset_index(drop=True)
    oof = pd.read_csv(C.C2_OOF, sep="\t").sort_values("row_index").reset_index(drop=True)
    if not labels.row_index.equals(oof.row_index):
        raise C.ContractError("label and out-of-fold tables are not row-aligned")
    if not np.allclose(labels.beta_alt.to_numpy(float), oof.beta_alt.to_numpy(float)):
        raise C.ContractError("beta_alt differs between the label table and the out-of-fold table")

    out: dict = {
        "card": str(local.card.iloc[0]),
        "card_tag": cards[0],
        "length_bp": lengths[0],
        "readout": "CenterMask 501 bp, DIFF_LOG2_SUM, mean over the liver tracks; no training",
        "bootstrap_resamples": int(args.resamples),
        "bootstrap_seed": C.BOOTSTRAP_SEED,
        "bootstrap_unit": "1Mb_block_within_fold",
        "comparison_seeds_note": (
            "the comparison seeds 1103/2909/4721/6673/8111 index the fitted heads of the c2 package, "
            "which this arm consumes as a seed-ensemble prediction; a zero-shot arm has nothing to "
            "seed, so no seed enters the local score"
        ),
        "no_sign_applied": True,
    }

    # ---------------------------------------------------------------- loader equivalence, 32,322
    y_all = labels.beta_alt.to_numpy(np.float64)
    folds_all = labels.heldout_fold.to_numpy(np.int64)
    blocks_all = labels.block_1mb.to_numpy(str)
    repro = {}
    failures = []
    for name, pub in PUBLISHED.items():
        v_all = oof[name].to_numpy(np.float64)
        got = C.macro_over_folds(y_all, v_all, folds_all)
        tol = TOLERANCE[name]
        rec = {
            "published_macro": pub["macro"],
            "recomputed_macro": got,
            "abs_difference": abs(got - pub["macro"]),
            "tolerance": tol,
            "within_tolerance": bool(abs(got - pub["macro"]) <= tol),
            "distinct_predicted_values_all": int(np.unique(v_all).size),
            "distinct_predicted_values_per_fold": {
                int(f): int(np.unique(v_all[folds_all == f]).size) for f in range(5)
            },
            "variants_per_fold": {
                int(f): int((folds_all == f).sum()) for f in range(5)
            },
        }
        repro[name] = rec
        if not rec["within_tolerance"]:
            failures.append((name, rec["abs_difference"], tol))
    out["loader_equivalence_32322"] = repro
    out["loader_equivalence_max_abs_difference"] = max(
        v["abs_difference"] for v in repro.values()
    )
    out["loader_equivalence_note"] = (
        "tolerance is per arm and the reason is the arm's tie structure: the two continuous arms "
        "have 32,322 distinct predicted values and reproduce at machine precision, while the "
        "allele-identity control has 12 to 13 distinct values per fold, so a last-bit difference in "
        "the deposited decimal moves which rows are tied and the average ranks shift with it"
    )
    if failures:
        raise C.ContractError(
            "this loader does not reproduce the published c2 macro Spearman within tolerance: "
            + "; ".join(f"{n} differs by {d:.3e} against tolerance {t:.1e}" for n, d, t in failures)
        )
    for name, rec in repro.items():
        C.log(
            f"loader equivalence {name}: abs diff {rec['abs_difference']:.3e} "
            f"(tolerance {rec['tolerance']:.1e}, {rec['distinct_predicted_values_all']} distinct values)"
        )

    # ---------------------------------------------------------------- the shared set
    # `oof.variant_id` is 'chr1:822354:C:T'; the gate key is '1:822354:C:T'.
    oof = oof.copy()
    oof["key"] = [
        C.key_of(*v.split(":")) for v in oof.variant_id.astype(str)
    ]
    keep = ["key", "row_index", "block_1mb", "heldout_fold", "beta_alt", HYENA, CADUCEUS, CONTROL]
    merged = local.merge(
        oof[keep].rename(
            columns={
                "block_1mb": "block_oof",
                "heldout_fold": "fold_oof",
                "beta_alt": "beta_alt_oof",
            }
        ),
        on="key",
        how="inner",
        validate="one_to_one",
    )
    if not (merged.block_1mb.astype(str) == merged.block_oof.astype(str)).all():
        raise C.ContractError("block assignment disagrees between the local table and the c2 deposit")
    if not (merged.heldout_fold.astype(int) == merged.fold_oof.astype(int)).all():
        raise C.ContractError("fold assignment disagrees between the local table and the c2 deposit")
    if not np.allclose(merged.beta_alt.to_numpy(float), merged.beta_alt_oof.to_numpy(float)):
        raise C.ContractError("beta_alt disagrees between the local table and the c2 deposit")

    y = merged.beta_alt.to_numpy(np.float64)
    folds = merged.heldout_fold.to_numpy(np.int64)
    blocks = merged.block_1mb.to_numpy(str)
    arms = {
        "alphagenome_local_zeroshot_atac": merged.local_atac_liver.to_numpy(np.float64),
        "alphagenome_local_zeroshot_dnase": merged.local_dnase_liver.to_numpy(np.float64),
        "hyenadna_delta_ridge": merged[HYENA].to_numpy(np.float64),
        "caduceus_delta_ridge": merged[CADUCEUS].to_numpy(np.float64),
        "control_allele_identity_ridge": merged[CONTROL].to_numpy(np.float64),
    }
    out["shared_set"] = {
        "variants": int(len(merged)),
        "blocks": int(merged.block_1mb.nunique()),
        "per_fold_variants": {
            int(f): int(n) for f, n in merged.heldout_fold.value_counts().sort_index().items()
        },
        "per_fold_blocks": {
            int(f): int(g.block_1mb.nunique()) for f, g in merged.groupby("heldout_fold")
        },
        "note": (
            "the Atlas-served leads, so a selected subset of the 32,322; every arm here is "
            "recomputed on it and only within-set contrasts are interpretable"
        ),
        "distinct_predicted_values_per_arm": {
            k: int(np.unique(v).size) for k, v in arms.items()
        },
        "distinct_predicted_values_per_arm_per_fold": {
            k: {int(f): int(np.unique(v[folds == f]).size) for f in range(5)}
            for k, v in arms.items()
        },
        "tie_structure_note": (
            "the allele-identity control is a ridge on a 16-level one-hot of (REF, ALT), so its "
            "prediction takes only a handful of distinct values per fold; its Spearman is capped by "
            "that tie structure and is not comparable in precision with a continuous arm's"
        ),
    }

    C.log(f"bootstrapping {len(arms)} arms on {len(merged)} leads, {merged.block_1mb.nunique()} blocks")
    macro_draws, pooled_draws = C.block_bootstrap_within_fold(
        y, arms, folds, blocks, resamples=args.resamples, seed=C.BOOTSTRAP_SEED
    )
    points_macro = {k: C.macro_over_folds(y, v, folds) for k, v in arms.items()}
    points_pooled = {k: C.fast_spearman(y, v) for k, v in arms.items()}

    summary = [
        summarise(
            k,
            macro_draws[k],
            pooled_draws[k],
            points_macro[k],
            points_pooled[k],
            len(merged),
            merged.block_1mb.nunique(),
        )
        for k in arms
    ]
    out["shared_set_arms"] = summary

    # per-fold components, so a macro is never reported without the folds it averages
    per_fold = []
    for f in range(5):
        sel = folds == f
        row = {"fold": f, "variants": int(sel.sum()), "blocks": int(len(set(blocks[sel].tolist())))}
        for k, v in arms.items():
            row[k] = C.fast_spearman(y[sel], v[sel])
        per_fold.append(row)
    out["per_fold_spearman"] = per_fold

    # ---------------------------------------------------------------- paired contrasts
    contrasts = []
    pairs = [
        ("alphagenome_local_zeroshot_atac", "hyenadna_delta_ridge"),
        ("alphagenome_local_zeroshot_atac", "caduceus_delta_ridge"),
        ("alphagenome_local_zeroshot_atac", "control_allele_identity_ridge"),
        ("alphagenome_local_zeroshot_dnase", "hyenadna_delta_ridge"),
        ("alphagenome_local_zeroshot_atac", "alphagenome_local_zeroshot_dnase"),
        ("hyenadna_delta_ridge", "control_allele_identity_ridge"),
    ]
    for a, b in pairs:
        for tag, draws, pts in (
            ("macro", macro_draws, points_macro),
            ("pooled", pooled_draws, points_pooled),
        ):
            d = draws[a] - draws[b]
            lo, hi, se = C.interval_of(d)
            p, note = C.bootstrap_p(d)
            contrasts.append(
                {
                    "contrast": f"{a} - {b}",
                    "statistic": tag,
                    "difference": pts[a] - pts[b],
                    "ci_low": lo,
                    "ci_high": hi,
                    "bootstrap_se": se,
                    "bootstrap_p": p,
                    "p_note": note,
                    "excludes_zero": bool(np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0)),
                    "variants": int(len(merged)),
                    "blocks": int(merged.block_1mb.nunique()),
                }
            )
    out["paired_contrasts"] = contrasts

    # ---------------------------------------------------------------- prespecified reference points
    out["prespecified_comparators_on_32322"] = {
        "hyenadna_delta_ridge_macro": PUBLISHED[HYENA]["macro"],
        "hyenadna_delta_ridge_ci": PUBLISHED[HYENA]["ci"],
        "control_allele_identity_macro": PUBLISHED[CONTROL]["macro"],
        "note": (
            "measured on all 32,322 leads, not on the 3,845 this arm covers; reported as the fixed "
            "reference point of ADAPTATION_ARM_PRESPEC.md and NOT as the paired contrast"
        ),
    }
    same_set = {s["arm"]: s["macro_spearman"] for s in summary}
    same_set_ci = {s["arm"]: (s["macro_ci_low"], s["macro_ci_high"]) for s in summary}
    out["subset_selection_effect"] = {
        "hyenadna_delta_ridge_macro_on_32322": PUBLISHED[HYENA]["macro"],
        "hyenadna_delta_ridge_macro_on_shared_set": same_set["hyenadna_delta_ridge"],
        "control_macro_on_32322": PUBLISHED[CONTROL]["macro"],
        "control_macro_on_shared_set": same_set["control_allele_identity_ridge"],
        "note": (
            "how much of any apparent local advantage is the subset rather than the model: the "
            "fitted comparators move on the same leads with no change to their own recipe"
        ),
    }

    # Is the shift material? Judged against the published arm's own bootstrap standard error and
    # against whether the recomputed value falls outside the published interval. A shift is a fact
    # about the subset, not an error in either run, and it is stated rather than left to inference.
    shift = []
    for col, pub in PUBLISHED.items():
        arm = ARM_NAME[col]
        got = same_set[arm]
        lo, hi = same_set_ci[arm]
        d = got - pub["macro"]
        shift.append(
            {
                "arm": arm,
                "published_macro_on_32322": pub["macro"],
                "published_ci_on_32322": pub["ci"],
                "published_bootstrap_se_on_32322": pub["se"],
                "recomputed_macro_on_shared_set": got,
                "recomputed_ci_on_shared_set": [lo, hi],
                "difference": d,
                "difference_in_published_se_units": d / pub["se"],
                "outside_published_ci": bool(got < pub["ci"][0] or got > pub["ci"][1]),
                "material": bool(abs(d) > 2 * pub["se"]),
            }
        )
    out["comparator_shift"] = shift
    out["comparator_shift_note"] = (
        "the published values were measured on all 32,322 Currin leads from 4,096-bp windows, where "
        "essentially nothing was rejected; this arm scores 1-Mb windows on the Atlas-served subset "
        "and rejects any window carrying an assembly-gap N or running off a chromosome end, so the "
        "two are different populations. A comparator marked material moved by more than two of its "
        "own published bootstrap standard errors. That is a statement about the subset, not about "
        "either run's correctness. Only the recomputed values enter a contrast."
    )
    out["comparators_that_moved_materially"] = [s["arm"] for s in shift if s["material"]]
    for s in shift:
        C.log(
            f"comparator shift {s['arm']}: published {s['published_macro_on_32322']:.4f} -> "
            f"shared set {s['recomputed_macro_on_shared_set']:.4f} "
            f"({s['difference_in_published_se_units']:+.1f} published SE, material={s['material']})"
        )

    # ------------------------------------------------- what the window guard removed, beside the result
    rejected_path = pathlib.Path(args.scores)
    rejected_path = rejected_path.with_name(rejected_path.stem + "_rejected.tsv")
    gate_tbl = pd.read_csv(run / "inputs/gate_variants.tsv", sep="\t").set_index("key")
    cov = {
        "atlas_served_leads": int(len(gate_tbl)),
        "currin_leads_total": 32322,
        "leads_scored_locally": int(len(local)),
        "leads_in_shared_set": int(len(merged)),
    }
    if rejected_path.exists():
        rd = pd.read_csv(rejected_path, sep="\t")
        rej_keys = [k for k in rd.key.astype(str) if k in gate_tbl.index]
        rej = gate_tbl.loc[rej_keys, "archived_atac_liver"].to_numpy(float)
        kept = gate_tbl.loc[
            [k for k in gate_tbl.index if k not in set(rej_keys)], "archived_atac_liver"
        ].to_numpy(float)
        rej_b = gate_tbl.loc[rej_keys, "beta_alt"].to_numpy(float)
        kept_b = gate_tbl.loc[
            [k for k in gate_tbl.index if k not in set(rej_keys)], "beta_alt"
        ].to_numpy(float)
        qs = [0.25, 0.5, 0.75]
        cov.update(
            {
                "leads_rejected": int(len(rd)),
                "rejected_fraction": float(len(rd) / len(gate_tbl)),
                "rejected_reasons": {k: int(v) for k, v in rd.reason.value_counts().items()},
                "rejected_abs_archived_atac_quartiles": [
                    float(x) for x in np.quantile(np.abs(rej), qs)
                ]
                if len(rej)
                else [],
                "kept_abs_archived_atac_quartiles": [
                    float(x) for x in np.quantile(np.abs(kept), qs)
                ],
                "rejected_median_abs_archived_atac": float(np.median(np.abs(rej)))
                if len(rej)
                else None,
                "kept_median_abs_archived_atac": float(np.median(np.abs(kept))),
                "rejected_median_abs_beta_alt": float(np.median(np.abs(rej_b)))
                if len(rej_b)
                else None,
                "kept_median_abs_beta_alt": float(np.median(np.abs(kept_b))),
                "survivors_skew": (
                    "towards larger archived effects"
                    if len(rej) and np.median(np.abs(kept)) > np.median(np.abs(rej))
                    else "towards smaller archived effects"
                    if len(rej)
                    else "no leads rejected"
                ),
                "note": (
                    "the leads a 1-Mb window guard removes are not a random sample of the gate set, "
                    "so the survivor set is the population every number in this section describes; "
                    "the published comparators describe a different one"
                ),
            }
        )
    out["coverage_and_survivor_skew"] = cov

    C.write_json(run / "tables/zeroshot_endpoint.json", out)
    pd.DataFrame(shift).to_csv(run / "tables/zeroshot_comparator_shift.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(run / "tables/zeroshot_arms.tsv", sep="\t", index=False)
    pd.DataFrame(contrasts).to_csv(run / "tables/zeroshot_contrasts.tsv", sep="\t", index=False)
    pd.DataFrame(per_fold).to_csv(run / "tables/zeroshot_per_fold.tsv", sep="\t", index=False)
    np.savez_compressed(
        run / "tables/zeroshot_bootstrap_samples.npz",
        **{f"macro__{k}": v for k, v in macro_draws.items()},
        **{f"pooled__{k}": v for k, v in pooled_draws.items()},
    )
    merged.to_csv(run / "tables/zeroshot_rows.tsv.gz", sep="\t", index=False)

    print(pd.DataFrame(summary).to_string(index=False))
    print(pd.DataFrame(contrasts).to_string(index=False))
    print(pd.DataFrame(shift).to_string(index=False))
    print(json.dumps(out["coverage_and_survivor_skew"], indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
