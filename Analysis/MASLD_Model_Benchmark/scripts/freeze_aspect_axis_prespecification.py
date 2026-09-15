#!/usr/bin/env python3
"""Stage 0 of the four-arm aspect decomposition: substrate, ceilings, label
arithmetic, and the sealed prespecification.

NOTHING here correlates a feature with a label. Stage 0 reads the four label
tables, measures the tie structure, measures the exact functional dependence
among the labels, computes the detectable-effect floor from n and the family
size, and writes the prespecification that Stage 1 must verify by digest before
it opens a matrix. Family sizes are a property of the substrate and are
recorded here so the family-adjusted floor is prespecified rather than chosen
after the counts are seen.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from aspect_axis_arms import load_arm_a, load_arm_b, load_arm_c, load_arm_d  # noqa: E402
from aspect_axis_common import (BENCH, canonical_sha256, mde_spearman,  # noqa: E402
                                tie_ceiling, tie_ceiling_from_counts)
from masld_bench.artifacts import freeze_tree, verify_frozen_tree  # noqa: E402

# Two constants published by earlier stages of this campaign. The ceiling
# formula must reproduce BOTH from their marginals alone, and a one-participant
# tamper must move the value. A check that cannot fail is not a check.
PINNED_CEILINGS = {
    "GSE267145 fibrosis 71/15/9/4": ([71, 15, 9, 4], 0.791771),
    "PXD051911 ballooning 41/16/1": ([41, 16, 1], 0.791170),
}

NAMED_ASPECTS = ["metabolism_steatosis", "inflammation_lobular", "fibrosis"]
CONDITIONAL_INFORMATIVENESS_FLOOR = 0.10


def ceiling_guard() -> dict:
    entries, failures = [], []
    for label, (counts, expected) in PINNED_CEILINGS.items():
        got = tie_ceiling_from_counts(counts)
        tampered = list(counts)
        tampered[0] += 1
        moved = tie_ceiling_from_counts(tampered)
        ok = abs(got - expected) < 5e-6
        fires = abs(moved - got) > 1e-4
        entries.append({
            "pinned": label, "expected": expected, "recomputed": round(got, 6),
            "reproduces_the_published_constant": bool(ok),
            "one_participant_tamper": tampered,
            "value_under_tamper": round(moved, 6),
            "tamper_moves_the_value": bool(fires),
        })
        if not ok or not fires:
            failures.append(label)
    return {
        "guard": "G_CEILING_PIN",
        "what_it_prints_if_the_thing_did_NOT_happen":
            "reproduces_the_published_constant false, or tamper_moves_the_value "
            "false, and the stage exits non-zero",
        "entries": entries, "failures": failures,
        "passed": not failures,
    }


def functional_dependence(x: np.ndarray, y: np.ndarray) -> bool:
    """Is x a deterministic function of y? True when x is constant inside every
    level of y. This is what catches an aspect that is a re-encoding of another
    rather than a second reading of the patient."""
    return all(len(set(x[y == lv])) == 1 for lv in np.unique(y))


def label_arithmetic(arm) -> dict:
    A = arm.aspects
    n = len(arm.participants)
    out = {"arm": arm.arm_id, "cohort": arm.cohort, "n_participants": n}

    if arm.components:
        s = A[arm.components].sum(axis=1).to_numpy()
        hold = int((s == A[arm.composite].to_numpy()).sum())
        out["composite_is_the_component_sum"] = {
            "identity": f"{arm.composite} == " + " + ".join(arm.components),
            "holds_on": f"{hold}/{n}",
            "exact": hold == n,
            "consequence": ("conditioning a component on the composite is label "
                            "arithmetic, not adjustment; every direction in this "
                            "arm conditions on the OTHER COMPONENTS instead"
                            if hold == n else "not an exact identity here"),
        }
    else:
        out["composite_is_the_component_sum"] = {
            "identity": "not evaluable: this arm deposits no components",
            "holds_on": None, "exact": None,
        }

    pairs, deps = {}, []
    cols = list(A.columns)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            rho = spearmanr(A[a], A[b]).statistic
            pairs[f"{a}|{b}"] = round(float(rho), 4)
            if functional_dependence(A[a].to_numpy(), A[b].to_numpy()):
                deps.append({"determined": a, "by": b})
            if functional_dependence(A[b].to_numpy(), A[a].to_numpy()):
                deps.append({"determined": b, "by": a})
    out["pairwise_label_spearman"] = pairs
    out["exact_functional_dependence"] = deps
    out["any_aspect_pair_is_a_deterministic_re_encoding"] = bool(deps)

    R = np.apply_along_axis(rankdata, 0, A.to_numpy(float))
    R = (R - R.mean(axis=0)) / R.std(axis=0)
    sv = np.linalg.svd(R / np.sqrt(len(R)), compute_uv=False)
    ev = sv ** 2
    out["label_matrix"] = {
        "aspects": cols,
        "singular_values": [round(float(v), 4) for v in sv],
        "eigenvalue_share": [round(float(v / ev.sum()), 4) for v in ev],
        "participation_ratio": round(float(ev.sum() ** 2 / (ev ** 2).sum()), 4),
        "reading": ("the participation ratio counts how many directions the "
                    "LABELS themselves span. It bounds nothing about the "
                    "transcriptome and is reported so a low-rank molecular "
                    "result can be checked against low-rank labels."),
    }
    return out


def ceiling_table(arm) -> list:
    rows = []
    n = len(arm.participants)
    family = arm.features.shape[1]
    n_nuis = arm.nuisance.shape[1]
    for aspect in arm.aspects.columns:
        v = arm.aspects[aspect].to_numpy(float)
        counts = pd.Series(v).value_counts().sort_index()
        others = [a for a in arm.aspects.columns
                  if a != aspect and a != arm.composite]
        if aspect == arm.composite:
            others = [a for a in arm.aspects.columns
                      if a != aspect and a not in arm.components]
        k_marg = n_nuis
        k_cond = n_nuis + len(others)
        ceil = tie_ceiling(v)
        mde_lib = mde_spearman(n, k_cond, 0.05)
        mde_fam = mde_spearman(n, k_cond, 0.05 / max(family, 1))
        if ceil < mde_lib:
            call = "untestable_ceiling_below_liberal_mde"
        elif ceil < mde_fam:
            call = "low_power_ceiling_below_family_mde"
        else:
            call = "testable"
        rows.append({
            "arm": arm.arm_id, "cohort": arm.cohort, "aspect": aspect,
            "n": n, "levels": int(counts.size),
            "marginal": "/".join(str(int(c)) for c in counts.to_numpy()),
            "tie_ceiling": round(ceil, 6),
            "family_size": int(family),
            "k_conditional": int(k_cond), "k_marginal": int(k_marg),
            "mde_liberal_alpha0.05": round(mde_lib, 4),
            "mde_family_alpha0.05_over_family": round(mde_fam, 4),
            "testability": call,
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--arm-d-lineage", default="Hepatocytes")
    ap.add_argument("--derived-after-directions", default=None,
                    help="directions.tsv from the calibration smoke whose "
                         "numbers were seen before this version was written")
    args = ap.parse_args()
    out = pathlib.Path(args.output)
    if out.exists():
        sys.exit(f"refusing to overwrite {out}")
    out.mkdir(parents=True)

    guard = ceiling_guard()
    print(json.dumps(guard, indent=2))
    if not guard["passed"]:
        sys.exit("G_CEILING_PIN failed; the tie-structure ceiling is not the "
                 "quantity the campaign already published")

    arms = {
        "A": load_arm_a(),
        "B": load_arm_b(),
        "C": load_arm_c(),
        "D": load_arm_d(lineage=args.arm_d_lineage),
    }
    for k, arm in arms.items():
        print(f"arm {k} {arm.cohort}: n={len(arm.participants)} "
              f"features={arm.features.shape[1]} aspects={list(arm.aspects.columns)}")
        for note in arm.notes:
            print(f"    - {note}")

    substrate = {
        "arms": [{
            "arm": a.arm_id, "cohort": a.cohort, "assay": a.assay,
            "n_participants": len(a.participants),
            "unit_of_inference": "participant" if a.arm_id != "D" else "donor",
            "feature_family_size": int(a.features.shape[1]),
            "aspects_present": list(a.aspects.columns),
            "nuisance_covariates": list(a.nuisance.columns),
            "notes": a.notes, "input_digests": a.inputs,
        } for a in arms.values()],
    }
    present = {asp: [a.arm_id for a in arms.values() if asp in a.aspects.columns]
               for asp in NAMED_ASPECTS + ["ballooning", "composite_nas"]}
    substrate["aspect_availability"] = present
    substrate["pooled_n_by_aspect"] = {
        asp: int(sum(len(arms[k].participants) for k in v))
        for asp, v in present.items()}
    substrate["roadmap_arithmetic_correction"] = {
        "roadmap_states": "pooled n across GSE267145 (99), GSE202379 (40), "
                          "PXD051911 (58) and GSE135251 (180) is 377",
        "realised": int(sum(len(a.participants) for a in arms.values())),
        "why": ("GSE202379 contributes the donors that are BOTH SAF-graded and "
                "present in the pseudobulk substrate; P70 is graded but absent "
                "from the expression substrate entirely."),
    }

    arithmetic = [label_arithmetic(a) for a in arms.values()]
    ceilings = [r for a in arms.values() for r in ceiling_table(a)]

    derived_after = None
    if args.derived_after_directions:
        smoke = pd.read_csv(args.derived_after_directions, sep="\t")
        n2 = smoke[smoke.family == "CONTROL_NEGATIVE_N2"]
        fired = n2[n2.verdict == "UNIQUE"]
        low = smoke[(smoke.verdict == "UNIQUE") & (smoke.count_bh05 < 5)]
        derived_after = {
            "what_was_run": ("a 40-draw calibration smoke of the Stage 1 "
                             "evaluator on the real four-arm substrate, written "
                             "to a scratch directory that is not a result"),
            "why_this_version_exists": (
                "the smoke's own negative control fired. Listing the numbers "
                "that were seen by value is what makes this derived-after "
                "exposure inspectable instead of something a reader has to "
                "reconstruct."),
            "n2_shuffled_exposure_control_fired_in": [
                {"arm": r.arm, "count": int(r.count_bh05),
                 "null_count_p95": float(r.null_count_p95),
                 "perm_p": float(r.perm_p)} for _, r in fired.iterrows()],
            "directions_called_unique_on_fewer_than_five_features": [
                {"arm": r.arm, "family": r.family, "tag": r.tag,
                 "count": int(r.count_bh05), "perm_p": float(r.perm_p)}
                for _, r in low.iterrows()],
            "amendments_this_caused": [
                "MIN_DISCOVERY_COUNT = 5: a direction whose BH count clears its "
                "null but is below five features is reported "
                "UNIQUE_BUT_BELOW_MINIMUM_COUNT and reads indeterminate, never "
                "supported. With a null that is a point mass at zero, a count of "
                "one is extreme by construction.",
                "the Stouffer combination clamps each per-arm permutation p to "
                "[1/(B+1), 1-1/(B+1)]; an unclamped p of exactly 1 sent the "
                "combined z to negative infinity in the smoke",
                "the harmonised meta selects on family as well as tag, because "
                "arm D's H2 direction and its full-conditional direction carry "
                "the same tag and were pooled twice",
            ],
            "smoke_directions_sha256": None,
            "what_was_NOT_changed": (
                "no criterion, no adjustment set, no aspect roster and no "
                "decision cell was altered. The three amendments are a "
                "discreteness floor, a numerical clamp and a de-duplication."),
        }
        derived_after["smoke_directions_sha256"] = __import__("hashlib").sha256(
            pathlib.Path(args.derived_after_directions).read_bytes()).hexdigest()

    prespec = {
        "prespec_id": ("aspect_axis_count_four_arm_v2" if derived_after
                       else "aspect_axis_count_four_arm_v1"),
        "derived_after_exposure": derived_after,
        "minimum_discovery_count": 5,
        "question": ("How many aspect axes does the substrate support, and are "
                     "metabolism, inflammation and fibrosis separable as outcome "
                     "OUTCOMES rather than as one severity reading?"),
        "roadmap_dimension": "docs/ROADMAP.md aspect decomposition, evidence "
                             "state indeterminate, no claim attached",
        "claim_boundary": (
            "This resolves how many aspect axes each arm's molecular substrate "
            "resolves. It is not a claim about MASLD biology, it authorises no "
            "aspect-specific gene attribution, and a negative fills the "
            "dimension exactly as a positive does."),
        "arms_are_four_instruments": (
            "GSE267145 is bulk RNA on a 0-3 fibrosis scale with no F4, "
            "GSE135251 is bulk RNA on 0-4 with no NAS components, PXD051911 is "
            "liver protein on F0-F3, GSE202379 is single-cell pseudobulk on the "
            "SAF system whose activity grade conflates ballooning and lobular "
            "inflammation. They enter as SEPARATE ARMS with per-arm "
            "decomposition. No feature space and no participant axis is pooled. "
            "The only pooling is of decision statistics."),
        "aspect_definitions": {
            "metabolism_steatosis": "steatosis grade, the histological reading "
                                    "of the metabolic aspect",
            "inflammation_lobular": "lobular inflammation grade",
            "fibrosis": "Kleiner or SAF fibrosis stage, arm-native scale",
            "ballooning": "hepatocellular ballooning; reported because it is a "
                          "NASH-CRN component, not one of the three named aspects",
            "composite_nas": "the deposited activity sum (NAS, or SAF activity)",
        },
        "statistic": {
            "per_feature": "partial Spearman of the feature against the exposure, "
                           "adjusting for the arm's other testable aspects and "
                           "the arm's nuisance covariates, all rank-transformed",
            "family": "BH 0.05 within each direction's own family",
            "null": "Freedman-Lane: the Z-residualised exposure is permuted, "
                    "which preserves the exposure's dependence on the adjustment "
                    "set. Permuting the RAW exposure does not, and is reported "
                    "as a secondary because the frozen 2026-08-27 runs used it.",
            "n_permutations": 5000,
        },
        "conditioning_rule_that_defeats_label_arithmetic": (
            "A component is NEVER conditioned on the composite that contains it. "
            "Where the composite is exactly the component sum, conditioning a "
            "component on it is arithmetic, not adjustment. Every three-aspect "
            "direction conditions on the OTHER COMPONENTS."),
        "testability_rules": {
            "untestable_no_label": "the aspect has no column in the arm",
            "untestable_ceiling_below_liberal_mde":
                "the tie-structure ceiling is below the alpha 0.05, power 0.80 "
                "MDE, so no molecular feature could reach significance however "
                "perfectly it read the aspect",
            "low_power_ceiling_below_family_mde":
                "the ceiling clears the liberal MDE but not the family-adjusted "
                "one; a null here is indeterminate, never negative",
            "conditionally_uninformative":
                f"the residual rank-variance fraction of the exposure after "
                f"adjustment is below {CONDITIONAL_INFORMATIVENESS_FLOOR}",
        },
        "channel_power_rule": (
            "Before any direction may carry a negative, the same aspect's "
            "MARGINAL direction (nuisance-adjusted only) must itself be UNIQUE "
            "against its own null. A channel whose known-positive yield is at "
            "chance cannot carry a negative, and its conditional null is then "
            "reported indeterminate."),
        "positive_control": (
            "sex adjusted for the composite and fibrosis, run on the identical "
            "instrument in every arm that deposits sex (A, C, D). Sex is a known "
            "strong and known separable axis. An arm whose sex control is not "
            "UNIQUE cannot carry a negative for any aspect."),
        "negative_controls": {
            "N1_duplicate_aspect": "an aspect conditioned on an exact copy of "
                                   "itself must return NOT_APPLICABLE_COLLINEAR, "
                                   "never a count",
            "N2_shuffled_exposure": "a fully permuted composite must not be "
                                    "UNIQUE; the realised perm p is reported",
        },
        "criteria": {
            "d1_an_aspect_is_a_separate_axis": {
                "decisive": True,
                "threshold": "BH 0.05 count above the 95th percentile of its own "
                             "permutation count null AND permutation p < 0.05 "
                             "AND count >= 5",
                "threshold_is_a_judgment_call": True,
                "rationale": "an aspect is a separate outcome axis only if it "
                             "carries features the other aspects do not explain",
            },
            "d2_an_empty_direction_is_classified_against_its_measured_floor": {
                "decisive": False,
                "threshold": "an empty direction whose observed max |partial "
                             "rho| is below the 95th percentile of the null's "
                             "max |rho| is EMPTY_UNDERPOWERED and reads "
                             "indeterminate; at or above it is "
                             "EMPTY_NOT_INDEPENDENT and reads as evidence",
                "applies_only_to_empty_directions": True,
            },
        },
        "decision_rule": {
            "per_aspect_per_arm": {
                "supported": "d1 met",
                "not_independent": "d2 says EMPTY_NOT_INDEPENDENT and the "
                                   "channel-power and positive controls both hold",
                "indeterminate": "everything else that was tested",
                "untestable": "no label, or a testability rule fired before any "
                              "result was read",
            },
            "axes_per_arm": "the number of aspects labelled supported",
            "declared_arm_outcomes": ["THREE_AXES", "TWO_AXES", "ONE_AXIS",
                                      "NO_AXIS_RESOLVED", "INDETERMINATE"],
            "a_one_axis_answer_is_a_result_not_a_failure": True,
        },
        "meta_analysis": {
            "unit": "the per-arm permutation p of the conditional direction, "
                    "never a feature and never a participant",
            "method": "Stouffer with sqrt(n) weights over the arms in which the "
                      "aspect is testable",
            "concentration_diagnostic": "each arm's share of the combined z, and "
                                        "the leave-one-arm-out recombination. A "
                                        "combined result driven by one arm is "
                                        "reported as a single-arm finding.",
            "floor": "with 5000 draws the per-arm p floors at 1/5001, so the "
                     "per-arm z caps near 3.54 and the combined z cannot be "
                     "driven arbitrarily high by one arm",
        },
        "guard_register": [
            "G_CEILING_PIN: the ceiling reproduces both published constants and "
            "moves under a one-participant tamper",
            "G_PRESPEC_DIGEST: Stage 1 exits before opening a matrix unless this "
            "tree's digest matches",
            "G_ZERO_JOIN: every join asserts its expected size and a zero join "
            "raises",
            "G_RAGGED_HEADER: field count is measured on raw bytes before pandas "
            "aligns anything",
            "G_FROZEN_INPUTS: the input digests are re-checked after the run",
            "N1_duplicate_aspect and N2_shuffled_exposure as declared above",
        ],
        "substrate": substrate,
        "label_arithmetic": arithmetic,
        "ceiling_and_power": ceilings,
        "conditional_informativeness_floor": CONDITIONAL_INFORMATIVENESS_FLOOR,
    }
    prespec["prespecification_sha256_of_content"] = canonical_sha256(prespec)

    (out / "aspect_axis_prespecification.json").write_text(
        json.dumps(prespec, indent=2, sort_keys=True), encoding="utf-8")
    pd.DataFrame(ceilings).to_csv(out / "ceiling_and_power.tsv", sep="\t",
                                  index=False)
    (out / "label_arithmetic.json").write_text(
        json.dumps(arithmetic, indent=2, sort_keys=True), encoding="utf-8")
    (out / "substrate_audit.json").write_text(
        json.dumps(substrate, indent=2, sort_keys=True), encoding="utf-8")
    (out / "ceiling_guard.json").write_text(json.dumps(guard, indent=2),
                                            encoding="utf-8")

    print("\n=== ceiling and power ===")
    print(pd.DataFrame(ceilings)[
        ["arm", "aspect", "n", "marginal", "tie_ceiling",
         "mde_liberal_alpha0.05", "mde_family_alpha0.05_over_family",
         "testability"]].to_string(index=False))
    print("\n=== label arithmetic ===")
    for entry in arithmetic:
        print(f"{entry['arm']} {entry['cohort']}: composite==sum "
              f"{entry['composite_is_the_component_sum']['holds_on']}, "
              f"deterministic re-encodings "
              f"{entry['exact_functional_dependence'] or 'none'}, "
              f"label participation ratio "
              f"{entry['label_matrix']['participation_ratio']}")
    digest = freeze_tree(out, {
        "artifact_class": "aspect_axis_count_prespecification",
        "stage": "stage_0_prespecification",
        "prespec_id": prespec["prespec_id"],
        "no_feature_was_correlated_with_any_label": True,
    })
    verify_frozen_tree(out)
    print(f"\nARTIFACTS.json sha256: {digest}")
    print(f"prespecification file sha256: "
          f"{__import__('hashlib').sha256((out / 'aspect_axis_prespecification.json').read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
