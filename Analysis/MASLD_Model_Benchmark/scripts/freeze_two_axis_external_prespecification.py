"""Freeze what would count as external support for the two-axis assignment.

Stage 1, job 2 of 3. Held-back before any outcome-versus-expression association
exists in either external arm.

Provenance. This descends from Stage 0 (STOP), Stage 0b (ONE_AXIS_ACTIVITY at
n=99) and Stage 0c (TWO_AXIS_ACTIVITY_FIBROSIS at n=180), each cited by digest.
The assignment under test is Stage 0c's, reproduced exactly rather than read: the
substrate build recomputed both partial families through the imported Stage 0b
code path and aborted unless six published invariants landed on the nose.

What is being claimed, and what is not. The deliverable is a **per-gene axis
assignment** whose external support is **set-level**. At n=78 and n=106 a
per-gene BH result will not reproduce, so no part of this stage tests, reports
or implies per-gene external replication. The evidence state on every assigned
gene must record that its external support is set-level, and that distinction is
the honesty constraint of the whole deliverable rather than a footnote.

Two criteria. S1 asks whether an assigned set is elevated on its own axis
against an expression-matched background. S2 asks whether it is *not* elevated
on the other axis, and S2 is the half that matters: without it the Catalog field
can claim association but never specificity, which is the only property that
makes an axis label useful.

S2 accepts a null, so it is guarded. A set that shows no elevation on the other
axis in an arm that could not have detected one is `indeterminate`, never
`tested_negative`. The set-level minimum detectable effect decides which, and it
is computed for the median-shift statistic rather than borrowed from a
single-correlation MDE, which is a different estimand with different power.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Sequence


class PrespecificationError(RuntimeError):
    """Raised when the substrate cannot carry the question being frozen."""


PRESPEC_ID = "two_axis_external_support_v1"

PARENT_DIGESTS = {
    "stage_0_result_aspect_separability_json":
        "4840ad08d9b2020fb6161856f4b3cccc054e38fce95b84653ccf7c1e440e14a2",
    "stage_0b_result_artifacts":
        "78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e",
    "stage_0b_derived_after_addendum_artifacts":
        "91b18fbe8cb061f04825beb3d9c7651ade92e6055782eacd06b6ece9c2cab49f",
    "stage_0b_bootstrap_df_correction_artifacts":
        "97db7a398b3b7dbb9a52c9e7599c3650b0c02b1277c5bcdc8431bde170fd267e",
    "stage_0c_prespecification_artifacts":
        "58f72a4522309f9d4e0dfb1630c37fd994a308e714856b88508aa40cb27732f6",
    "stage_0c_result_artifacts":
        "5f26597ab6b6aeeaff867d54112554a33d4d0dc35be5f333e17e9010a7199fed",
    "stage_1_substrate_artifacts":
        "68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab",
    "activation_record":
        "3c5d3cc4e3f47b331c5f8ece04c213c304d49faa3ce3b974cef878c57e4f6ca4",
}

EXCLUSIVE_CELLS = ("activity_only", "fibrosis_only")
ALL_CELLS = ("activity_only", "fibrosis_only", "both", "neither")
OWN_AXIS = {"activity_only": "activity", "fibrosis_only": "fibrosis"}
OTHER_AXIS = {"activity_only": "fibrosis", "fibrosis_only": "activity"}

#: Which arm is permitted to serve which axis, and why.
ARM_SCOPE = {
    "GSE130970": {
        "axes_served": ["activity", "fibrosis"],
        "unit": "bulk_rna_sample",
        "n": 78,
        "donor_key_exists": False,
    },
    "GSE193066": {
        "axes_served": ["activity"],
        "unit": "participant",
        "n": 106,
        "donor_key_exists": True,
    },
}

FROZEN_THRESHOLDS = {
    "s1_the_assigned_set_is_elevated_on_its_own_axis": (
        "the set's median absolute partial association on its own axis exceeds "
        "the 95th percentile of a null that permutes assignment among "
        "expression-matched genes"
    ),
    "s2_the_assigned_set_is_not_elevated_on_the_other_axis": (
        "the set's median absolute partial association on the other axis does "
        "NOT exceed the 95th percentile of the same matched-permutation null, "
        "AND the arm's set-level minimum detectable effect is at or below the "
        "elevation observed on the set's own axis; if the arm could not have "
        "detected an elevation that size, the result is indeterminate rather "
        "than met"
    ),
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build(substrate: dict, activation: dict) -> dict:
    cells = substrate["assignment"]["cell_sizes"]
    if sum(cells.values()) != substrate["reproduction_of_stage_0c"]["realised_universe"]:
        raise PrespecificationError("the four cells do not partition the universe")
    if cells["activity_only"] + cells["both"] != 4388:
        raise PrespecificationError("the activity marginal does not reconstruct")
    if cells["fibrosis_only"] + cells["both"] != 1305:
        raise PrespecificationError("the fibrosis marginal does not reconstruct")
    arms = substrate["external_arms"]
    if arms["GSE193066"]["samples"] != 164:
        raise PrespecificationError(
            "the deposited GSE193066 matrix is expected to carry 164 columns; "
            "the note that 164 is not the analysis unit depends on it"
        )

    return {
        "schema_version": "masld-bench-two-axis-external-prespec-v1",
        "prespec_id": PRESPEC_ID,
        "question": (
            "Does the Stage 0c two-axis gene assignment find external support, "
            "at set level, in cohorts outside the one it was derived in?"
        ),
        "the_deliverable_and_its_honesty_constraint": {
            "deliverable": "a per-gene axis assignment",
            "external_support_is_set_level_not_per_gene": True,
            "why": (
                "At n=78 and n=106 a per-gene BH result will not reproduce. No "
                "criterion, diagnostic or reported number in this stage tests "
                "or implies per-gene external replication."
            ),
            "what_every_assigned_gene_must_record": (
                "its evidence state must say that its external support is "
                "set-level. A per-gene assignment carrying set-level evidence "
                "is a legitimate object; a per-gene assignment presented as "
                "per-gene replicated is not."
            ),
        },
        "honest_provenance": {
            "derived_after_three_completed_stages": True,
            "parents": PARENT_DIGESTS,
            "stage_0_decision": "STOP",
            "stage_0b_outcome": "ONE_AXIS_ACTIVITY",
            "stage_0c_outcome": "TWO_AXIS_ACTIVITY_FIBROSIS",
            "the_assignment_was_reproduced_not_read": (
                "Stage 0c deposited counts and summaries only, so the substrate "
                "build recomputed both partial families through the imported "
                "Stage 0b code path and aborted unless the realised universe, "
                "both family sizes, both BH counts and both maximum absolute "
                "partial correlations reproduced exactly. A digest read would "
                "have proved only that the right file was opened."
            ),
            "a_lesson_carried_forward": (
                "Stage 0 deposited a per-gene table; Stage 0b and Stage 0c did "
                "not, and neither could hand its memberships to a successor. "
                "The same omission recurs here in a smaller form: the substrate "
                "build deposited assignment without the signed partial "
                "associations, so a direction-concordance check has to "
                "re-derive them. Recorded rather than silently repaired."
            ),
            "no_external_outcome_association_existed_at_the_freeze": True,
            "criteria_fixed_before_any_external_number_was_computed": True,
        },
        "substrate": {
            "training_cohort": "GSE135251, n=180, in-cohort for the assignment",
            "assignment_cell_sizes": cells,
            "jackknife_membership_stability": substrate["assignment"][
                "jackknife_membership_stability"]["per_cell"],
            "the_both_cell_is_the_least_stable": (
                "51.7% of its genes keep the assignment in every one of the 180 "
                "leave-one-out fits, against 74.9% for activity_only and 64.8% "
                "for fibrosis_only. It travels with that number wherever the "
                "cell appears, because it is the cell carrying an "
                "association-without-specificity claim."
            ),
            "join_loss_per_arm": {
                name: arm["join_loss_per_assigned_cell"] for name, arm in arms.items()
            },
            "scope_statement": (
                "every external result applies to the joined subset of each "
                "assigned set, with the surviving fraction named, and never to "
                "the full assigned set"
            ),
            "how_the_mappability_confound_is_removed": (
                "by construction, not by a threshold: the expression-matched "
                "background is drawn from the arm's own joined universe, so an "
                "assigned set and its background always share a mappability "
                "regime. A threshold on unequal loss would invite a judgment "
                "call at the moment the answer is visible."
            ),
        },
        "unit_of_inference": {
            "GSE130970": {
                "unit": "bulk_rna_sample",
                "n": 78,
                "distinct_people_assertable": False,
                "why": (
                    "no donor key is deposited and no repeat-specimen indicator "
                    "exists. Absence of an indicator is not proof of one "
                    "specimen per person, so 78 samples may not be restated as "
                    "78 people."
                ),
            },
            "GSE193066": {
                "unit": "participant",
                "n": 106,
                "rows_in_the_deposited_matrix": 164,
                "REQUIRED_NOTE_ON_THE_SUBSTRATE_ARTIFACT": (
                    "The Stage 1 substrate artifact records samples: 164 for "
                    "GSE193066. That is the width of the deposited matrix and "
                    "is NOT the analysis unit. The analysis unit is 106 "
                    "first-biopsy participants. 164 must never be read as the "
                    "n. The 58 additional rows are repeat specimens of people "
                    "already present, and using them would pseudoreplicate the "
                    "activity arm."
                ),
                "selection_rule": (
                    "first biopsy only, one row per participant, keyed from "
                    "!Sample_title. The suffix forms are heterogeneous: 58 "
                    "first biopsies carry _1, 48 carry no suffix, and 58 second "
                    "biopsies carry _2. Stripping _2 pairs zero of 58 and "
                    "silently reports 164 participants, so the key is taken "
                    "from the deposited title rather than derived by "
                    "truncation."
                ),
                "the_analysis_must_assert_exactly_106_distinct_participants": True,
            },
        },
        "arm_scope": {
            "arms": ARM_SCOPE,
            "gse193066_is_restricted_to_the_activity_axis": {
                "restricted": True,
                "ground_that_carries_the_decision": (
                    "Its cross-sectional fibrosis is [6, 37, 37, 26, 0] and "
                    "omits stage 4 entirely. Under the approved standard a "
                    "scale that omits a level is a different measurement "
                    "instrument, reported separately and never meta-analysed "
                    "against a graded scale. The two arms are therefore not "
                    "measuring the same variable, which is categorical and not "
                    "a power claim better data could overturn."
                ),
                "supporting_coverage_caveats_not_independent_grounds": [
                    "fibrosis_only joins at 75.45% against activity_only's "
                    "92.57%. This bounds generalization, it does not bias S2: "
                    "the specificity test evaluates the same gene set on both "
                    "axes, so whatever the join dropped is dropped from both "
                    "evaluations identically.",
                    "its neither pool joins at 33.69%. This is the weakest of "
                    "the three: 19,015 genes is ample for expression matching, "
                    "and the differential runs toward safety because the "
                    "background is attrited harder than the assigned sets "
                    "precisely by being less well expressed and annotated, "
                    "which narrows the set-versus-background mappability gap "
                    "rather than widening it.",
                ],
                "the_fibrosis_arm_is_still_computed": (
                    "Declining to compute an arm because its answer is "
                    "predicted to be inconvenient is the same failure as "
                    "choosing a substrate because it gives a preferred answer. "
                    "It is computed and deposited under diagnostics, enters no "
                    "outcome cell, and the analysis aborts if it appears in any "
                    "gate."
                ),
                "how_the_fibrosis_arm_result_may_be_cited": (
                    "neither as support nor as refutation for the fibrosis "
                    "half. Only as a flag that the fibrosis result may be "
                    "arm-specific and that a third arm is warranted."
                ),
            },
        },
        "criteria": {
            "s1_the_assigned_set_is_elevated_on_its_own_axis": {
                "decisive": True,
                "statistic": (
                    "the median absolute partial Spearman of the assigned set's "
                    "joined genes with its own axis, adjusting for the other "
                    "axis, in the external arm"
                ),
                "comparator": (
                    "an expression-matched background drawn from the arm's own "
                    "joined neither pool"
                ),
                "threshold": FROZEN_THRESHOLDS[
                    "s1_the_assigned_set_is_elevated_on_its_own_axis"],
                "threshold_is_a_judgment_call": True,
                "evaluated_on": list(EXCLUSIVE_CELLS),
            },
            "s2_the_assigned_set_is_not_elevated_on_the_other_axis": {
                "decisive": True,
                "why_it_is_the_half_that_matters": (
                    "without it the Catalog field can claim association but "
                    "never specificity, and specificity is the only property "
                    "that makes an axis label useful"
                ),
                "statistic": (
                    "the same median absolute partial association, computed "
                    "against the OTHER axis"
                ),
                "threshold": FROZEN_THRESHOLDS[
                    "s2_the_assigned_set_is_not_elevated_on_the_other_axis"],
                "threshold_is_a_judgment_call": True,
                "evaluated_on": list(EXCLUSIVE_CELLS),
                "why_the_both_cell_is_excluded": (
                    "genes in the both cell are BH-significant on both axes by "
                    "construction, so a specificity test that included them "
                    "would fail for an arithmetic reason rather than an "
                    "empirical one"
                ),
                "s2_accepts_a_null_so_it_is_power_guarded": (
                    "a set showing no elevation on the other axis in an arm "
                    "that could not have detected one is indeterminate, never "
                    "tested_negative. The set-level MDE decides which."
                ),
            },
        },
        "power": {
            "set_level_mde_is_the_operative_reference": {
                "definition": (
                    "the smallest median-shift in absolute partial association, "
                    "injected into a random expression-matched set of the same "
                    "size, that the matched-permutation null detects at 80% "
                    "power in that arm"
                ),
                "why_it_is_computed_rather_than_borrowed": (
                    "the primary statistic is a set-level median shift, not a "
                    "single correlation. Borrowing a single-correlation MDE for "
                    "it would be a method-to-estimand mismatch: a set of "
                    "thousands of genes can show a decisive median shift where "
                    "no individual gene clears the single-gene floor, and the "
                    "reverse is also possible."
                ),
                "is_a_gate": False,
            },
            "single_correlation_mde_is_context_only": {
                "values": {
                    "GSE130970_bh_family_12": 0.40367682076575373,
                    "GSE193066_bh_family_12": 0.3502,
                    "GSE130970_nominal_0_05": 0.3126674774693899,
                },
                "what_they_are_the_reference_for": (
                    "the per-gene claim this stage explicitly does not make. "
                    "They are recorded as context and no null is judged against "
                    "them."
                ),
            },
            "ceilings": {
                "GSE130970_fibrosis": 0.9558402968719953,
                "GSE193066_fibrosis": 0.948727110241683,
                "caveat": (
                    "a ceiling is a property of the outcome tie structure. A "
                    "ceiling above an MDE means an effect is detectable, not "
                    "that an effect exists."
                ),
            },
        },
        "method": {
            "partial_association": (
                "the rank-linear partial correlation used throughout this lane, "
                "imported from the Stage 0b analysis rather than "
                "reimplemented, so training and external arms share one "
                "instrument"
            ),
            "expression_matching": (
                "genes are binned by their mean CPM decile within the arm's own "
                "joined universe, and the background is drawn to match the "
                "assigned set's decile histogram exactly. An unmatched "
                "background is the failure that produced the "
                "orthogonality-ladder artifact on this project."
            ),
            "null": (
                "assignment labels are permuted among expression-matched genes, "
                "never shuffled globally"
            ),
            "normalization": (
                "counts per million over the sample axis. A per-gene Spearman "
                "across samples is invariant to a monotone per-gene transform "
                "but not to per-sample scaling, which is what library size is. "
                "GSE193066's library-size ratio is 8.22 and GSE130970's is "
                "1.57, so this is not theoretical."
            ),
            "library_size_orthogonality_is_measured_per_arm": True,
            "p_values_are_reported_at_their_resolution_floor": True,
        },
        "diagnostics_never_gates": {
            "direction_concordance": {
                "is_a_gate": False,
                "what_it_is": (
                    "the fraction of an assigned set's joined genes whose "
                    "external partial association carries the same sign as its "
                    "re-derived training partial association"
                ),
                "why_it_is_reported": (
                    "S1 and S2 read magnitude. An elevated median absolute "
                    "association with no sign agreement would be a much weaker "
                    "result than the same number with agreement, and the "
                    "difference should be visible rather than assumed."
                ),
                "why_it_is_not_a_gate": (
                    "it was not in the approved design and adding a third gate "
                    "after the design was settled would be the mid-flight "
                    "change this lane has already paid for once"
                ),
            },
            "gse193066_fibrosis_arm": {
                "is_a_gate": False,
                "enters_no_outcome_cell": True,
                "citable_as": (
                    "a flag that the fibrosis result may be arm-specific and a "
                    "third arm is warranted"
                ),
                "not_citable_as": ["support for the fibrosis half",
                                   "refutation of the fibrosis half"],
            },
            "both_cell_reported_separately": {
                "is_a_gate": False,
                "catalog_field_states": "association, explicitly not specificity",
                "n_genes": cells["both"],
            },
            "neither_cell_size_is_reported": {
                "is_a_gate": False,
                "n_genes": cells["neither"],
                "why": "it is the denominator that makes the other three interpretable",
            },
        },
        "evidence_states": {
            "rule": (
                "a set whose assignment fails to replicate in an arm powered to "
                "detect it is tested_negative for that arm; one in an arm not "
                "so powered is indeterminate. The two are never merged."
            ),
            "decided_by": "the arm's set-level MDE against the observed elevation",
        },
        "decision_rule": {
            "outcome_map_is_complete": True,
            "cells": [
                {"outcome": "EXTERNALLY_SUPPORTED",
                 "when": ("S1 and S2 hold in GSE130970 for both exclusive cells "
                          "AND in GSE193066 for activity_only")},
                {"outcome": "SUPPORTED_SINGLE_ARM",
                 "when": "S1 and S2 hold in GSE130970 only"},
                {"outcome": "TRANSFERS_WITHOUT_SPECIFICITY",
                 "when": ("S1 holds and S2 fails; the field may state "
                          "association and must not state specificity")},
                {"outcome": "NOT_EXTERNALLY_SUPPORTED",
                 "when": ("S1 fails; the two-axis model stays in-cohort and the "
                          "Gene Catalog says so")},
                {"outcome": "INDETERMINATE",
                 "when": "either gate reaches zero applicable conditions"},
            ],
            "vacuity_guard": (
                "any gate reaching zero applicable conditions returns "
                "NO_APPLICABLE_CONDITIONS and never a silent pass; all([]) is "
                "True and is guarded explicitly"
            ),
        },
        "claim_boundary": {
            "the_fibrosis_axis_has_exactly_one_external_arm": (
                "The fibrosis axis has exactly one external arm, n=78 samples, "
                "with no assertable donor key."
            ),
            "training_is_in_cohort": (
                "the two-axis result is a property of GSE135251 until an "
                "external arm reports; if S1 fails it remains a property of "
                "GSE135251"
            ),
            "no_clinical_or_prognostic_claim": True,
            "external_support_is_set_level": True,
            "activation_user_approved_scope_only": activation["approval"][
                "technical_adjudication"]["user_approved"],
            "no_per_cohort_verdict_may_be_presented_as_user_approved": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--substrate", type=Path, required=True)
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise PrespecificationError("refusing to overwrite a prespecification")
    payload = build(read_json(arguments.substrate), read_json(arguments.activation))
    for key, expected in FROZEN_THRESHOLDS.items():
        if payload["criteria"][key]["threshold"] != expected:
            raise PrespecificationError(f"{key} threshold does not match the literal")
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "prespec_id": payload["prespec_id"],
        "criteria": sorted(payload["criteria"]),
        "outcomes": [c["outcome"] for c in payload["decision_rule"]["cells"]],
        "arms": {k: v["axes_served"] for k, v in ARM_SCOPE.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
