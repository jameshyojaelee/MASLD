#!/usr/bin/env python3
"""Derived-after, NON-DECISIVE addendum to the four-arm aspect decomposition.

The sealed prespecification decides the axis count with a conditional criterion.
This adds a purely descriptive second reading of the same substrate: how many
directions the ASSOCIATION VECTORS span, per arm, measured as the participation
ratio of the aspect-by-aspect Gram matrix of the marginal association vectors.

It is written and run AFTER the decisive result and it changes no verdict. It
exists because the label-side participation ratio was already reported at
Stage 0 and a molecular-side counterpart is the natural check on it: a low-rank
molecular reading is only interesting if the labels are not low rank for
arithmetic reasons.

Null: the participant order is permuted ONCE per draw and applied to every
aspect at the same time, so the aspect-to-aspect label correlation is preserved
exactly and only the feature-to-label link is destroyed. Under that null the
association vectors are near-orthogonal and the participation ratio approaches
the number of aspects, so the test is one sided and it can fail.
"""
from __future__ import annotations
import argparse, json, pathlib, sys
import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from aspect_axis_arms import load_arm_a, load_arm_b, load_arm_c, load_arm_d
from aspect_axis_common import rank_columns, residualise, standardise
from masld_bench.artifacts import freeze_tree, verify_frozen_tree


def participation_ratio(B: np.ndarray) -> float:
    G = B.T @ B / B.shape[0]
    e = np.linalg.eigvalsh(G)
    e = np.clip(e, 0, None)
    if e.sum() <= 0:
        return float("nan")
    return float(e.sum() ** 2 / (e ** 2).sum())


def association_matrix(Zs: np.ndarray, aspects: np.ndarray, C: np.ndarray,
                       order: np.ndarray) -> np.ndarray:
    # The FEATURE rows are permuted, not the labels. Permuting the labels would
    # also break their alignment with the nuisance design; permuting the feature
    # residuals destroys exactly the feature-to-label link and leaves the
    # aspect-to-aspect and aspect-to-covariate structure untouched.
    Zp = Zs[order, :]
    cols = []
    for j in range(aspects.shape[1]):
        a = residualise(rankdata(aspects[:, j]).reshape(-1, 1), C).ravel()
        sd = a.std()
        a = (a - a.mean()) / sd if sd > 0 else a * 0
        cols.append((Zp * a[:, None]).sum(axis=0) / len(a))
    return np.column_stack(cols)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--result", required=True,
                    help="the decisive Stage 1 result this addendum follows")
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260901)
    args = ap.parse_args()
    out = pathlib.Path(args.output)
    if out.exists():
        sys.exit(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    decisive = json.loads((pathlib.Path(args.result) /
                           "aspect_axis_result.json").read_text())

    rows = []
    for key, loader in [("A", load_arm_a), ("B", load_arm_b),
                        ("C", load_arm_c), ("D", load_arm_d)]:
        arm = loader()
        n = len(arm.participants)
        A = arm.aspects.to_numpy(float)
        names = list(arm.aspects.columns)
        C = np.column_stack([np.ones(n)] + [arm.nuisance[c].to_numpy(float)
                                            for c in arm.nuisance.columns])
        Zs, ok = standardise(residualise(rank_columns(arm.features), C))
        Zs = Zs[:, ok]
        obs = participation_ratio(association_matrix(Zs, A, C, np.arange(n)))
        rng = np.random.default_rng(args.seed + ord(key))
        null = np.array([participation_ratio(
            association_matrix(Zs, A, C, rng.permutation(n)))
            for _ in range(args.n_perm)])
        p_low = float((1 + int((null <= obs).sum())) / (1 + args.n_perm))
        rows.append({
            "arm": key, "cohort": arm.cohort, "n": n, "aspects": names,
            "n_aspects": len(names),
            "observed_participation_ratio": round(obs, 4),
            "null_mean": round(float(null.mean()), 4),
            "null_p05": round(float(np.percentile(null, 5)), 4),
            "one_sided_p_observed_is_lower": p_low,
            "reading": ("the association vectors span fewer directions than the "
                        "aspect count" if p_low < 0.05 else
                        "the association vectors are not measurably lower rank "
                        "than the null"),
        })
        print(json.dumps(rows[-1]), flush=True)

    payload = {
        "status": "DERIVED_AFTER_NON_DECISIVE",
        "follows_result": str(args.result),
        "decisive_outcomes_unchanged": {
            a: s["three_aspect_outcome"] for a, s in decisive["per_arm"].items()},
        "what_this_is_not": ("this changes no verdict, authorises no claim and "
                             "was written after the decisive result was read"),
        "n_permutations": args.n_perm, "seed": args.seed,
        "per_arm": rows,
    }
    (out / "dimensionality_addendum.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    freeze_tree(out, {"artifact_class": "aspect_axis_dimensionality_addendum",
                      "stage": "derived_after", "decisive": False})
    verify_frozen_tree(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
