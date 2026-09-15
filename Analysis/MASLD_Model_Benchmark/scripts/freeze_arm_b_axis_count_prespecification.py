"""Freeze whether activity and fibrosis are two axes in GSE135251 at n=180.

Stage 0c of the MASLD showcase model, job 1 of 2.

Provenance, stated first. **This question was derived after two completed
results**, and the output file says so rather than implying a clean start:

*   Stage 0 (`model-check-1044`) asked whether steatosis, ballooning and lobular
    inflammation are three separable axes in GSE267145 and returned STOP.
*   Stage 0b (`model-check-1047`) asked how many axes that substrate supports
    and returned ONE_AXIS_ACTIVITY. Its decisive criterion found 805 genes for
    the NAS sum adjusted for fibrosis and **zero** for fibrosis adjusted for the
    NAS sum, with fibrosis's best gene at |partial r| 0.3767 against a measured
    familywise floor of 0.4502 -- below the level at which any gene could have
    been detected. Its addendum reproduced that sub-floor failure under a second
    covariate definition (0.4026 against a floor of 0.4487).

So Stage 0b could not decide whether fibrosis is a separate axis; it could only
establish that GSE267145 cannot answer. That is the question this stage moves to
a substrate that can. Arm A carries fibrosis 71/15/9/4 on a 0-3 scale with 54.3%
of participant pairs tied and a split-half gene-ordering reliability of 0.357.
Arm B carries 35/41/48/44/12 on a 0-4 scale at n=180.

**Nothing has been computed in GSE135251.** No axis-versus-axis correlation and
no axis-versus-expression association exists for this cohort, deliberately, so
that holding back first is worth something. This file inspects the two axes'
marginals and tie structure and **does not compute the correlation between
them**; that number comes into existence only after this tree is frozen and only
inside the analysis job.

The prior observations that *do* exist are all Arm A, and every one of them is
recorded here by value, so the derived-after exposure is on the record rather
than inferred later by a reader.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import t as student_t


class PrespecificationError(RuntimeError):
    """Raised when the substrate cannot carry the question being frozen."""


PRESPEC_ID = "gse135251_axis_count_v1"

#: Every parent this question descends from, cited by digest.
PARENT_DIGESTS = {
    "stage_0_prespecification_artifacts": (
        "5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76"
    ),
    "stage_0_result_aspect_separability_json": (
        "4840ad08d9b2020fb6161856f4b3cccc054e38fce95b84653ccf7c1e440e14a2"
    ),
    "stage_0_post_hoc_observations_json": (
        "8775349ea96abfa9563318f873f360dd6cbe2391ca34f571f9a3100d2851e024"
    ),
    "stage_0b_prespecification_artifacts": (
        "e7397c3038dfbbdd1619048acfdd6cabc0e018a7af161175925130a7a6ffed3e"
    ),
    "stage_0b_result_artifacts": (
        "78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e"
    ),
    "stage_0b_derived_after_addendum_artifacts": (
        "91b18fbe8cb061f04825beb3d9c7651ade92e6055782eacd06b6ece9c2cab49f"
    ),
}

ACTIVITY = "nas_score"
FIBROSIS = "fibrosis_stage"
FAMILIES = ((ACTIVITY, FIBROSIS), (FIBROSIS, ACTIVITY))

BH_ALPHA = 0.05
COUNT_FLOOR_PERCENTILE = 95.0

_CRITICAL_T = float(student_t.isf(0.025, 177))

FROZEN_THRESHOLDS = {
    "e1_activity_and_fibrosis_are_independent": (
        "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
        "the matched permutation count null"
    ),
    "e2_an_empty_direction_is_classified_against_the_measured_floor": (
        "an empty direction is underpowered if its observed max |partial r| is "
        "below the measured familywise floor p95, and not_independent if it "
        "reaches or clears that floor"
    ),
}


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def tie_profile(values: Sequence[int], name: str, scale: str) -> dict[str, object]:
    counts = Counter(values)
    n = len(list(values))
    tied = sum(count * (count - 1) // 2 for count in counts.values())
    return {
        "axis": name,
        "native_scale": scale,
        "distinct_values_realised": len(counts),
        "value_counts": {str(key): counts[key] for key in sorted(counts)},
        "largest_tied_block": max(counts.values()),
        "fraction_of_pairs_tied": tied / (n * (n - 1) // 2),
    }


#: Every Arm A number this stage descends from, recorded by value.
ARM_A_PRIOR_OBSERVATIONS = {
    "why_they_are_listed": (
        "This stage was designed after these were measured. Listing them by "
        "value is what makes the derived-after exposure inspectable instead of "
        "something a reader has to reconstruct."
    ),
    "substrate": "GSE267145, n=99",
    "label_spearman_sealed_before_expression_was_opened": {
        "steatosis|saf_activity_sum": 0.7409,
        "steatosis|nas_activity_sum": 0.9138,
        "steatosis|fibrosis": 0.3702,
        "saf_activity_sum|nas_activity_sum": 0.9466,
        "saf_activity_sum|fibrosis": 0.6494,
        "nas_activity_sum|fibrosis": 0.5721,
    },
    "partial_families_bh_0_05_counts": {
        "nas_activity_sum|fibrosis": 805,
        "fibrosis|nas_activity_sum": 0,
        "steatosis|saf_activity_sum": 30,
        "saf_activity_sum|steatosis": 0,
        "saf_activity_sum|fibrosis": 141,
        "fibrosis|saf_activity_sum": 0,
    },
    "empty_directions_sat_below_their_measured_floor": {
        "fibrosis|nas_activity_sum": {"observed_max": 0.3767, "floor_p95": 0.4502},
        "fibrosis|saf_activity_sum": {"observed_max": 0.4026, "floor_p95": 0.4487},
        "saf_activity_sum|steatosis": {"observed_max": 0.3959, "floor_p95": 0.4540},
    },
    "marginal_bh_0_05_counts": {
        "nas_activity_sum": 1885, "steatosis": 1819, "lobular_inflammation": 1137,
        "saf_activity_sum": 963, "ballooning": 519, "fibrosis": 8,
    },
    "split_half_gene_ordering_reliability_full_length": {
        "nas_activity_sum": 0.5715, "steatosis": 0.5530,
        "lobular_inflammation": 0.5184, "saf_activity_sum": 0.5164,
        "ballooning": 0.4790, "fibrosis": 0.3574,
    },
}


def build(rows: Sequence[dict[str, str]]) -> dict[str, object]:
    n = len(rows)
    if n != 180:
        raise PrespecificationError(f"expected 180 participants, found {n}")
    if len({row["participant_id"] for row in rows}) != n:
        raise PrespecificationError("participant_id is not unique")

    activity = [int(row[ACTIVITY]) for row in rows]
    fibrosis = [int(row[FIBROSIS]) for row in rows]
    if max(fibrosis) != 4:
        raise PrespecificationError(
            "Arm B fibrosis must realise stage 4; the scale hazard depends on it"
        )
    if max(activity) > 8 or min(activity) < 0:
        raise PrespecificationError("nas_score outside its 0-8 scale")

    return {
        "schema_version": "masld-bench-arm-b-axis-count-prespec-v1",
        "prespec_id": PRESPEC_ID,
        "question": (
            "Does the GSE135251 bulk transcriptome at n=180 support activity "
            "(the deposited NAS sum) and fibrosis stage as two separate outcome "
            "axes, or does one fail to carry genes the other does not explain?"
        ),
        "honest_provenance": {
            "this_question_was_derived_after_two_completed_results": True,
            "parents": PARENT_DIGESTS,
            "stage_0_decision": "STOP",
            "stage_0b_outcome": "ONE_AXIS_ACTIVITY",
            "what_stage_0b_could_not_decide": (
                "whether fibrosis is a separate axis. Its fibrosis direction "
                "returned zero genes with a best |partial r| of 0.3767 against a "
                "measured familywise floor of 0.4502, reproduced at 0.4026 "
                "against 0.4487 under a second covariate. An empty count below "
                "the floor at which anything could be detected is a statement "
                "about the substrate, not about fibrosis."
            ),
            "why_this_substrate": (
                "Arm A carries fibrosis 71/15/9/4 with 54.3% of participant "
                "pairs tied and a split-half gene-ordering reliability of 0.357, "
                "the lowest of any axis measured there. Arm B carries "
                "35/41/48/44/12 at n=180 and deposits both axes natively. The "
                "move is to the substrate that can answer, not to the substrate "
                "that gives a preferred answer, and the Arm A result is reported "
                "rather than discarded."
            ),
            "criteria_fixed_before_any_arm_b_number_was_computed": True,
            "nothing_has_been_computed_in_this_cohort": (
                "No axis-versus-axis correlation and no axis-versus-expression "
                "association exists for GSE135251. The freeze inspects the two "
                "axes' marginals and tie structure only and deliberately does "
                "not compute the correlation between them; that number comes "
                "into existence after this tree is frozen."
            ),
            "in_cohort_and_why_that_is_acceptable": (
                "GSE135251 is the training cohort for the transfer models in "
                "this campaign, so an answer here is in-cohort. That is "
                "acceptable because this gate asks whether the substrate "
                "supports two axes, not whether a model transfers. No model is "
                "fitted and no held-out claim is made. Stated here rather than "
                "left implicit."
            ),
            "derivable_later_is_not_recorded_before": (
                "Nothing here is backdated. The analysis job must verify this "
                "tree's digest before it opens the expression matrix."
            ),
        },
        "arm_a_prior_observations": ARM_A_PRIOR_OBSERVATIONS,
        "cohort": {
            "dataset_id": "gse135251_human_liver",
            "n_participants": n,
            "unit_of_inference": "participant",
            "admission": "model-check-892-21132281-gse135251-admission, admitted",
            "endpoint_table": (
                "executions/model-data-880-21130257-gse135251-source/source/"
                "outcomes/participant_endpoints.tsv"
            ),
            "molecular_fixture": (
                "executions/model-data-880-21130257-gse135251-source/source/molecular/"
            ),
            "n_is_180_not_216": (
                "The source metadata carries full labels for 216 runs and the "
                "36 extra are recorded as a reconciliation table, but the 180 "
                "sealed here underlie already-scored transfer results and are "
                "not re-cut. n=180, deliberately, not 216."
            ),
        },
        "axes": {
            ACTIVITY: "deposited NAS sum, native scale 0-8, all nine levels realised",
            FIBROSIS: "deposited fibrosis stage, native scale 0-4",
            "activity_cannot_be_decomposed_here": (
                "GSE135251 deposits only the NAS sum. Steatosis, ballooning and "
                "lobular inflammation are not separate columns, so the "
                "steatosis-versus-activity question Stage 0b answered in Arm A "
                "cannot be asked here at all. That omission is a property of the "
                "substrate, not a choice, and it is why Stage 0b's D2 is not "
                "repeated."
            ),
        },
        "fibrosis_scale_hazard": {
            "arm_a_native_scale": "0-3, no stage 4 present",
            "arm_b_native_scale": "0-4, all five stages realised",
            "the_hazard": (
                "F3 does not denote the same thing across these two substrates, "
                "so any Arm A to Arm B fibrosis comparison crosses a scale "
                "boundary. This is the GSE213621 trap in a new place: that "
                "cohort's harmonised stage silently pooled F0F1 and F3F4 with no "
                "level 4."
            ),
            "no_recode_is_applied_here": (
                "Stage 0c makes no cross-substrate fibrosis comparison. Every "
                "value carries its native scale and nothing is harmonised. Any "
                "future comparison needs an explicit labelled recode; inventing "
                "one that this stage does not use would be worse than naming the "
                "hazard."
            ),
            "whether_arm_a_lacks_a_level_4_or_sampled_no_f4_participant": (
                "not determinable from the deposit. Recorded as undetermined "
                "rather than guessed."
            ),
        },
        "what_was_inspected_before_freezing": {
            "no_gene_expression_was_opened": True,
            "no_correlation_between_the_two_axes_was_computed": True,
            "marginals_and_tie_structure_only": {
                ACTIVITY: tie_profile(activity, ACTIVITY, "0-8"),
                FIBROSIS: tie_profile(fibrosis, FIBROSIS, "0-4"),
            },
        },
        "criteria": {
            "e1_activity_and_fibrosis_are_independent": {
                "decisive": True,
                "statistic": (
                    "per-gene partial Spearman, counted at BH 0.05 within each "
                    "direction's own family"
                ),
                "directions": [
                    {"exposure": exposure, "adjusted_for": covariate}
                    for exposure, covariate in FAMILIES
                ],
                "threshold": FROZEN_THRESHOLDS["e1_activity_and_fibrosis_are_independent"],
                "threshold_is_a_judgment_call": True,
                "rationale": (
                    "Two axes are independent only if each carries genes the "
                    "other does not explain, and the relation is not symmetric, "
                    "so both directions are required. The count must also clear "
                    "a permutation floor because BH on a rank-based partial "
                    "statistic is not guaranteed calibrated at ordinal tie "
                    "structure."
                ),
                "identical_in_form_to_stage_0b_d1": (
                    "deliberately, so the two substrates are compared on the "
                    "same instrument rather than on two differently tuned ones"
                ),
            },
            "e2_an_empty_direction_is_classified_against_the_measured_floor": {
                "decisive": False,
                "is_a_disambiguation_rule_not_a_pass_fail_gate": True,
                "threshold": FROZEN_THRESHOLDS[
                    "e2_an_empty_direction_is_classified_against_the_measured_floor"
                ],
                "threshold_is_a_judgment_call": True,
                "rationale": (
                    "Stage 0b made this comparison unprompted and it turned out "
                    "to be the number that settled the interpretation, so it is "
                    "stated in advance here rather than discovered again. An "
                    "empty count whose best gene sits below the level at which "
                    "anything could have been detected is underpowered; an empty "
                    "count whose best gene clears that level is genuinely empty."
                ),
                "applies_only_to_empty_directions": (
                    "If no direction is empty this criterion has zero applicable "
                    "conditions and returns NO_APPLICABLE_CONDITIONS. That is "
                    "the correct reading and is not a failure."
                ),
            },
        },
        "diagnostics_never_gates": {
            "split_half_reliability_per_axis": {
                "is_a_gate": False,
                "why_it_is_mandatory_beside_every_count": (
                    "A less reliably estimated axis yields fewer BH-surviving "
                    "genes for a measurement reason rather than a biological "
                    "one, so a count comparison across axes is not interpretable "
                    "without it. Arm A measured fibrosis at 0.357 against 0.5715 "
                    "for the NAS sum; whether Arm B's fibrosis is materially "
                    "better at 35/41/48/44/12 is itself informative, and if it "
                    "is not, that is the answer to why the axis will not resolve."
                ),
                "prediction_recorded_before_the_measurement": (
                    "Arm B fibrosis should exceed Arm A's 0.357, because its "
                    "tie fraction is far lower and n is nearly double. Recorded "
                    "as a prediction about the instrument and reported whichever "
                    "way it falls, as Arm A's prediction was."
                ),
            },
            "marginal_counts_for_both_axes": {
                "is_a_gate": False,
                "why": (
                    "an empty partial direction can mean not independent or "
                    "underpowered, and the marginal count in the same universe "
                    "separates the readings. Arm A's fibrosis marginal was 8 "
                    "genes, which is what made its empty partial unsurprising."
                ),
            },
            "participant_jackknife_not_bootstrap": {
                "is_a_gate": False,
                "why_the_bootstrap_was_dropped": (
                    "Stage 0b's participant bootstrap of a BH count was invalid: "
                    "every observed count sat below its own bootstrap 2.5th "
                    "percentile, because resampling with replacement duplicates "
                    "participants whose identical expression and identical "
                    "labels manufacture exact concordance blocks that inflate "
                    "rank correlation. A leave-one-participant-out jackknife has "
                    "no duplication and is used instead."
                ),
            },
        },
        "method": {
            "partial_spearman_with_ties": (
                "Every variable is average-rank transformed. The partial "
                "association between a gene and an exposure given a covariate is "
                "the Pearson correlation between the residual of the ranked gene "
                "and the residual of the ranked exposure, each after "
                "least-squares projection on [1, ranked covariate]. This is the "
                "rank-linear partial correlation, and it is the identical code "
                "path Stage 0b used, imported rather than reimplemented."
            ),
            "what_it_is_not": (
                "It is not a test of conditional independence and it is not "
                "proof of an independent axis. A surviving partial association "
                "says the exposure carries genes the named covariate does not "
                "explain, not that it is independent of everything. Activity and "
                "fibrosis are correlated in every MASLD cohort, so nothing in "
                "this design can deliver independence in the strong sense, and "
                "a TWO_AXIS outcome names a modelling decision rather than a "
                "claim about orthogonal biology."
            ),
            "per_gene_p_value": (
                "two-sided Student t on the partial correlation with df = n - 3 "
                "= 177, reported beside a permutation p from the same draws"
            ),
            "multiplicity": (
                "Benjamini-Hochberg at 0.05 within each direction's own family, "
                "using masld_bench.evaluators.stats.benjamini_hochberg, with the "
                "family named and its realised size reported beside every count"
            ),
            "count_null": (
                "the exposure residual is permuted and re-residualized on the "
                "covariate and the whole family recomputed, so each draw yields "
                "a complete BH count. ter Braak-type residual permutation: "
                "asymptotically exact rather than exact, which is why the floor "
                "is measured instead of assumed."
            ),
            "primary_reference_is_global_and_why_it_differs_from_arm_a": (
                "Stage 0b used a within-recorded_sex permutation because Arm A "
                "deposits sex and Stage 0 measured it to carry structure. "
                "GSE135251 deposits no sex column and no other biological "
                "stratum. group_in_paper and both outer-fold columns are derived "
                "from the outcome itself, so permuting within them would hold "
                "fixed the very variation the criteria ask about. The global "
                "residual permutation is therefore the primary reference here, "
                "and the difference from Arm A is stated rather than silent."
            ),
            "p_values_are_reported_at_their_resolution_floor": True,
        },
        "power": {
            "n": 180,
            "covariates": 1,
            "residual_df": 177,
            "analytic_partial_rank_null_sd_1_over_sqrt_n_minus_2": 1.0
            / float(np.sqrt(178.0)),
            "two_sided_0_05_critical_abs_partial_r": float(
                _CRITICAL_T / np.sqrt(_CRITICAL_T**2 + 177.0)
            ),
            "note": (
                "Analytic figures are per-gene references. The operative floor "
                "is the measured max-statistic permutation floor over each "
                "realised family, and the criteria read a count rather than a "
                "single correlation, so the measured count null is what E1 is "
                "judged against."
            ),
            "why_this_should_be_better_than_arm_a": (
                "n rises from 99 to 180 and fibrosis moves from 71/15/9/4 with "
                "54.3% of pairs tied to 35/41/48/44/12. If the fibrosis "
                "direction is still empty and still below its measured floor "
                "here, that is a much stronger statement than Arm A's was."
            ),
        },
        "controls": {
            "unit_of_inference_is_the_participant": True,
            "every_count_carries_a_measured_permutation_floor": True,
            "realised_universe_is_re_derived_never_copied": True,
            "vacuity_audit": (
                "every gate condition is classified applicable or not before it "
                "is evaluated and the verdict is reported as met over "
                "applicable, never met over total. A gate with zero applicable "
                "conditions returns NO_APPLICABLE_CONDITIONS and never a silent "
                "pass; all([]) is True and is guarded explicitly."
            ),
            "no_model_is_fitted": True,
            "no_cross_substrate_fibrosis_comparison_is_made": True,
        },
        "decision_rule": {
            "a_one_axis_answer_is_a_result_not_a_failure": True,
            "outcome_map_is_complete": True,
            "cells": [
                {
                    "e1": "both directions survive",
                    "outcome": "TWO_AXIS_ACTIVITY_FIBROSIS",
                    "meaning": (
                        "activity and fibrosis each carry genes the other does "
                        "not explain at n=180. The showcase model keeps two "
                        "outcome axes."
                    ),
                },
                {
                    "e1": "activity survives, fibrosis empty and below its floor",
                    "outcome": "ONE_AXIS_ACTIVITY_FIBROSIS_UNDERPOWERED",
                    "meaning": (
                        "fibrosis is indeterminate, not tested_negative, even at "
                        "n=180. Two substrates would then have failed to power "
                        "it and the honest report is that this design cannot "
                        "resolve fibrosis on bulk transcriptome."
                    ),
                },
                {
                    "e1": "activity survives, fibrosis empty and clears its floor",
                    "outcome": "ONE_AXIS_ACTIVITY_FIBROSIS_NOT_INDEPENDENT",
                    "meaning": (
                        "fibrosis carries no gene beyond activity at a level "
                        "where such genes were detectable. This is the only "
                        "cell in which a fibrosis null is evidence rather than "
                        "an absence of evidence."
                    ),
                },
                {
                    "e1": "fibrosis survives, activity empty and below its floor",
                    "outcome": "ONE_AXIS_FIBROSIS_ACTIVITY_UNDERPOWERED",
                    "meaning": "the mirror case, declared because the gate must be symmetric",
                },
                {
                    "e1": "fibrosis survives, activity empty and clears its floor",
                    "outcome": "ONE_AXIS_FIBROSIS_ACTIVITY_NOT_INDEPENDENT",
                    "meaning": "the mirror case, declared because the gate must be symmetric",
                },
                {
                    "e1": "neither direction survives",
                    "outcome": "NO_AXIS_RESOLVED",
                    "meaning": (
                        "neither axis carries genes the other does not explain. "
                        "The aspect layer is dropped and both axes are recorded "
                        "indeterminate."
                    ),
                },
            ],
            "no_applicable_conditions": (
                "if E1 has zero applicable conditions the outcome is "
                "INDETERMINATE and E1 reports NO_APPLICABLE_CONDITIONS"
            ),
            "e2_qualifies_the_outcome_name_it_never_overrides_e1": True,
        },
        "claim_boundary": (
            "This is a statement about how many outcome axes the GSE135251 bulk "
            "transcriptome resolves in these 180 participants. It is not a claim "
            "about MASLD biology, it authorises no aspect-specific attribution, "
            "and it is in-cohort with respect to the transfer models trained on "
            "this dataset."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise PrespecificationError("refusing to overwrite a prespecification")
    payload = build(read_table(arguments.endpoints))
    for key, expected in FROZEN_THRESHOLDS.items():
        if key not in payload["criteria"]:
            raise PrespecificationError(f"missing criterion {key}")
        if payload["criteria"][key]["threshold"] != expected:
            raise PrespecificationError(f"{key} threshold does not match the literal")
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "prespec_id": payload["prespec_id"],
        "criteria": sorted(payload["criteria"]),
        "decisive": "e1_activity_and_fibrosis_are_independent",
        "outcomes": [cell["outcome"] for cell in payload["decision_rule"]["cells"]],
        "n": payload["cohort"]["n_participants"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
