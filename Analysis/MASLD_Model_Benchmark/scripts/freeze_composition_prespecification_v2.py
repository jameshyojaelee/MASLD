"""Hold back again the composition criteria after three defects, two of them in the design.

Stage 3 v2, job 1 of 2. Cites ``composition_dependence_v1`` as parent and
records run ``model-check-1055-21183164`` as **not a result**.

⛔ **21183164's outcome, COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY, is not a
result and may not be cited.** Two independent causes:

1.  The training arm silently dropped out. ``rna_feature_axis.tsv`` was indexed
    by position at column 0, which is ``feature_index``; ``stable_gene_id`` is
    column 1. The two external arms read from a different axis whose column 0
    genuinely *is* ``stable_gene_id``, so one line was correct in two arms and
    wrong in the third. Measured: the join is 61,940 of 61,940 on the right
    column and 0 on the wrong one. Both key spaces were ENSG throughout, so this
    was never a key-space mismatch and never needed a crosswalk.
2.  The verdict it produced was a threshold coin-flip. The two arms' drops were
    0.0491 and 0.0483, and the arm that "passed" cleared its null by 0.00188
    against a set-level minimum detectable shift of 0.005. An outcome named for
    a between-arm difference rested on 0.0008.

**The column lesson runs both ways in this tree.** The ragged-header traps
taught *measure column contents, never trust the name*. This one teaches *read
by name, never trust the position*. Neither rule alone is sufficient: a column
must be resolved explicitly, by whichever of name or position the file actually
guarantees, and its contents validated afterwards. Both failures now exist in
the same tree.

Two criteria change, which is why this is a new prespecification rather than a
bug fix. T1's power guard becomes symmetric, because it was written only for the
mediated direction and a pass inside the detection limit is equally
indistinguishable from the null. T2 moves to the one arm where per-gene BH is
meaningful, because Stage 1 established that per-gene BH does not reproduce at
n=78 or n=106 and the frozen v1 criterion asked for exactly that.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


class PrespecificationError(RuntimeError):
    """Raised when the substrate cannot carry the question being frozen."""


PRESPEC_ID = "composition_dependence_v2"

PARENTS = {
    "composition_dependence_v1_prespecification":
        "4e837504d7a7a4968b99f04276a4c232a3c5b84a2c9a6c253725a3f9fe1f83f5",
    "composition_substrate":
        "97f6c12407ac5d093afc978a777983b14804d971963c3291d6510059f9112623",
    "stage_1_substrate":
        "68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab",
    "stage_1_external_result":
        "daf4f02db0585f96b5605c09f4f71d6182e77d2d9bae25e263b34b1ba2e7eece",
    "not_a_result_run":
        "e6e9bb04064fe1f1d9df1f8346b6633f819e6953ccb00d2e4c04ba6886d610cb",
}

#: The only arm where a per-gene BH criterion is meaningful.
CONFIRMATORY_ARM = "GSE135251"
EXTERNAL_ARMS = ("GSE130970", "GSE193066")

FROZEN_THRESHOLDS = {
    "t1_a_composition_independent_component_exists": (
        "after adjusting for the six lineage proportions in addition to the "
        "other axis, the assigned set's median absolute partial association on "
        "its own axis exceeds the 95th percentile of a matched-permutation null "
        "BY AT LEAST that arm's set-level minimum detectable shift; a margin "
        "inside the MDE is indeterminate on either side of the percentile"
    ),
    "t2_per_gene_composition_dependence": (
        "in GSE135251 only, a gene is composition_independent if it retains BH "
        "0.05 on its own axis after the composition adjustment, "
        "composition_mediated if it loses it while the arm's set-level MDE shows "
        "the loss was detectable, and indeterminate otherwise; in the external "
        "arms every gene is indeterminate_by_construction"
    ),
}


def build(substrate: dict) -> dict:
    family = substrate["eligible_family"]["family"]
    if len(family) != 6:
        raise PrespecificationError(f"the family is {len(family)}, not the approved six")

    return {
        "schema_version": "masld-bench-composition-prespec-v2",
        "prespec_id": PRESPEC_ID,
        "supersedes": "composition_dependence_v1",
        "parents": PARENTS,
        "the_prior_run_is_not_a_result": {
            "run": "model-check-1055-21183164-composition-dependence",
            "its_outcome": "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY",
            "may_it_be_cited": False,
            "cause_1_the_training_arm_dropped_out_silently": {
                "what": (
                    "rna_feature_axis.tsv was indexed by position at column 0, "
                    "which is feature_index; stable_gene_id is column 1"),
                "measured": {
                    "join_on_stable_gene_id": "61940 of 61940",
                    "join_on_feature_index": 0,
                    "both_key_spaces_were_ensg": True,
                },
                "why_it_looked_plausible": (
                    "the two external arms read from the Stage 1 gene axis, "
                    "whose column 0 genuinely is stable_gene_id, so the same "
                    "line was correct in two arms and wrong in the third"),
                "it_was_never_a_key_space_mismatch": (
                    "so it never needed a crosswalk; it needed reading the "
                    "column by name"),
                "why_it_reached_a_verdict": (
                    "the cells came back applicable: false rather than raising, "
                    "and the wrapper checked covariate rank, participant count "
                    "and axis scope but never that every arm contributed an "
                    "applicable condition"),
            },
            "cause_2_the_verdict_was_a_threshold_coin_flip": {
                "the_two_arms_drops": [0.0491, 0.0483],
                "difference": 0.0008,
                "the_passing_margin": 0.00188,
                "that_arms_set_level_mde": 0.005,
                "reading": (
                    "the pass was smaller than the smallest shift the arm can "
                    "resolve, and an outcome named for a between-arm difference "
                    "rested on 0.0008"),
            },
        },
        "correction_to_the_frozen_substrate_artifact": {
            "carried_forward_from_v1_because_it_is_still_live": (
                "v2 binds the same substrate artifact, whose prose still "
                "misstates the count. Dropping the correction because the "
                "criteria changed would lose a record deliberately added, so it "
                "travels with v2 rather than being reachable only through the "
                "v1 parent."),
            "field": "variant_sensitivity_on_the_family.what_it_shows",
            "what_it_says": (
                "Only Hepatocytes and Macrophages are stable across variants"),
            "what_is_true": (
                "three are stable, not two: Hepatocytes, Macrophages and pDCs"),
            "how_it_happened": (
                "the prose was propagated from an upstream summary without being "
                "checked against the data beside it in the same object"),
            "why_it_is_not_retrofitted": (
                "the substrate artifact is frozen and its build script is hashed "
                "into that job's source manifest"),
            "the_substantive_point_is_unchanged": (
                "three of six flip and the backup family is half the size, so "
                "the eligible family is variant-dependent"),
        },
        "closure": {
            "no_between_lineage_claims": (
                "a per-lineage association with an external variable is well "
                "posed; a reciprocal shift is arithmetic. Any between-lineage "
                "statement would need a subcomposition check, and the macrophage "
                "result reversed +0.225 to -0.165 under exactly that test."),
            "one_shared_operator_not_independent_measurements": (
                "cross-cohort agreement is one deconvolution operator "
                "transferring, not independent measurements. The 13_bayesprism "
                "reference path does not exist and REF_HUMAN was overridden "
                "without the value being recorded; the byte-identical 16-name "
                "roster is what confirms one label set, not the path."),
            "never_align_the_mouse_roster_by_name": (
                "different ontologies sharing 11 names"),
        },
        "the_column_lesson_runs_both_ways_in_this_tree": {
            "from_the_ragged_headers": "measure column contents, never trust the name",
            "from_this_defect": "read by name, never trust the position",
            "the_unifying_rule": (
                "a column must be resolved explicitly, by whichever of name or "
                "position the file actually guarantees, and its contents "
                "validated afterwards. Neither rule alone is sufficient, and "
                "this tree now contains failures from both directions."),
        },
        "the_freeze_contract_caught_what_no_test_could": {
            "what_happened": (
                "a global sed renaming the analysis script to _v2 also rewrote "
                "the test filename inside the source-manifest line, into a file "
                "that never existed. The extensionless unittest line survived, "
                "so the tests ran and passed."),
            "why_it_matters": (
                "no test could have caught a manifest pointing at a missing "
                "file; the freeze contract failed the job on it and preserved "
                "the tree as .failed"),
        },
        "the_framing_that_governs_every_reading": {
            "a_composition_dependent_gene_is_not_an_artifact": True,
            "why": (
                "fibrosis genuinely involves cell-population change, so "
                "composition is mediation, not confounding"),
            "the_two_labels_are_both_real_findings": [
                "composition_mediated", "composition_independent"],
            "forbidden": (
                "calling a composition-dependent gene spurious, or naming T1's "
                "failure in any way that implies the two-axis result collapsed"),
        },
        "honest_provenance": {
            "derived_after_a_defective_run_whose_outcome_is_void": True,
            "criteria_fixed_before_any_v2_number_was_computed": True,
            "what_changed_from_v1": [
                "T1's power guard is symmetric rather than one-sided",
                "T2 is confirmatory in GSE135251 only",
                "a new outcome cell forbids naming a within-resolution "
                "difference as a between-arm finding",
            ],
            "what_did_not_change": [
                "the eligible family of six",
                "the mediation framing and SET_IS_COMPOSITION_MEDIATED",
                "the covariate construction and df = n - 9",
            ],
        },
        "substrate": {
            "eligible_family": family,
            "family_size": len(family),
            "derivation_rule": substrate["eligible_family"]["rule"],
            "quantification_is_not_uniform": substrate[
                "quantification_is_not_uniform_across_the_arms"],
        },
        "criteria": {
            "t1_a_composition_independent_component_exists": {
                "decisive": True,
                "what_it_asks": (
                    "is there a composition-independent component, not whether "
                    "the axis holds"),
                "threshold": FROZEN_THRESHOLDS[
                    "t1_a_composition_independent_component_exists"],
                "threshold_is_a_judgment_call": True,
                "failure_is_named": "SET_IS_COMPOSITION_MEDIATED",
                "the_guard_is_symmetric": (
                    "v1 applied the MDE only to the mediated direction. A pass "
                    "whose margin over the null percentile is smaller than the "
                    "arm's minimum detectable shift is equally indistinguishable "
                    "from the null, so it is indeterminate too."),
                "six_lineages_is_a_stricter_test_than_five": (
                    "a larger covariate set makes T1 harder to pass, so a "
                    "mediated verdict becomes more likely for a partly "
                    "technical reason; the MDE separates that from a real "
                    "mediation finding"),
            },
            "t2_per_gene_composition_dependence": {
                "decisive": False,
                "is_the_deliverable": True,
                "threshold": FROZEN_THRESHOLDS["t2_per_gene_composition_dependence"],
                "confirmatory_arm": CONFIRMATORY_ARM,
                "external_arms_are_indeterminate_by_construction": {
                    "arms": list(EXTERNAL_ARMS),
                    "reason": (
                        "Stage 1 established that per-gene BH does not "
                        "reproduce at n=78 or n=106, which is why Stage 1 was "
                        "set-level. v1 then asked for per-gene BH in those same "
                        "arms and 98.6% of GSE130970's activity_only cell "
                        "landed in indeterminate by construction."),
                    "reported_not_dropped": (
                        "the record should show the question was asked and why "
                        "the substrate cannot answer it"),
                },
                "why_gse135251": (
                    "n=180, and it is the arm the assignment was derived in, so "
                    "per-gene BH is meaningful there"),
            },
        },
        "decision_rule": {
            "outcome_map_is_complete": True,
            "cells": [
                {"outcome": "COMPOSITION_INDEPENDENT_COMPONENT_EXISTS",
                 "when": "T1 holds with a margin at or above the MDE in every arm serving the cell"},
                {"outcome": "ARMS_DIFFER_WITHIN_RESOLUTION",
                 "when": ("the arms reach different verdicts but their margins "
                          "differ by less than the smallest MDE among them. The "
                          "cell exists so a threshold coin-flip cannot be named "
                          "a between-arm finding.")},
                {"outcome": "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY",
                 "when": ("the arms genuinely differ, by more than either can "
                          "resolve")},
                {"outcome": "SET_IS_COMPOSITION_MEDIATED",
                 "when": ("T1 fails where the MDE shows the drop was "
                          "detectable; the axis acts through composition, which "
                          "is a mechanism finding")},
                {"outcome": "INDETERMINATE",
                 "when": ("every applicable margin sits inside its MDE, or a "
                          "gate reaches zero applicable conditions")},
            ],
            "vacuity_guard": (
                "an arm matching zero assigned genes, or holding an empty "
                "background pool, raises rather than reporting applicable: "
                "false. A silent arm dropout is how v1 reached a verdict."),
        },
        "controls": {
            "every_arm_must_contribute_an_applicable_condition": True,
            "gene_ids_resolved_by_name_and_validated": True,
            "row_sums_use_tolerance_not_equality": True,
            "the_sealed_instrument_is_imported_not_edited": True,
            "no_model_is_fitted": True,
        },
        "claim_boundary": (
            "This is a statement about whether assigned sets retain axis "
            "association after adjusting for six lineage proportions in these "
            "cohorts. It authorises no between-lineage claim, no causal "
            "mediation claim, and no clinical claim. Composition here is one "
            "deconvolution operator's output, not a measured cell count."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--substrate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise PrespecificationError("refusing to overwrite a prespecification")
    payload = build(json.loads(arguments.substrate.read_text(encoding="utf-8")))
    for key, expected in FROZEN_THRESHOLDS.items():
        if payload["criteria"][key]["threshold"] != expected:
            raise PrespecificationError(f"{key} threshold does not match the literal")
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "prespec_id": payload["prespec_id"],
        "supersedes": payload["supersedes"],
        "prior_run_citable": payload["the_prior_run_is_not_a_result"]["may_it_be_cited"],
        "outcomes": [c["outcome"] for c in payload["decision_rule"]["cells"]],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
