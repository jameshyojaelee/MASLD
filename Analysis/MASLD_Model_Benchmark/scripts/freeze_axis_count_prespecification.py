"""Freeze what would count as one, two or three histological outcome axes.

Stage 0b of the MASLD showcase model, job 1 of 2.

Provenance, stated first because it is the thing a reader must weigh. **This
question was derived after seeing Stage 0's answer.** Stage 0 asked whether
steatosis, ballooning and lobular inflammation are three separable axes in
GSE267145 and returned STOP: C1 failed (ballooning|inflammation 0.818 above the
0.80 bar), C2 failed (participation ratio 1.42 of 3.0 against a 2.0 floor) and
the decisive C3 failed (gene-ordering correlations 0.825, 0.872, 0.872, none
below 0.80). Its prespecification tree carries ARTIFACTS.json digest
``5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76``.

Stage 0b asks the successor question -- *how many* axes does this substrate
support, and which grouping -- and it is a declared successor rather than a post
hoc rationalisation only because of three things this file records: the
provenance above is stated in the output file in those words, the Stage 0 digest is
cited, and every criterion and literal threshold below was fixed before any
expression value was read. The analysis job verifies this tree's digest with
``sha256sum --check --strict`` before it opens the matrix.

Stage 0 measured *marginal* association-vector correlations. Stage 0b measures
*partial* associations after adjustment, and counts genes. That is a different
estimand, and this file records a measured audit of how much of it is already
recoverable from Stage 0's wrote output files rather than asserting that it is
new. The audit is run on the labels alone, before the freeze, and its numbers
are written into the prespecification verbatim:

*   ``rank(fibrosis)`` sits only 46.5% inside the span of the three Stage 0
    aspect rank vectors, and no fibrosis association vector was deposited, so
    both D1 directions require a measurement Stage 0 never made.
*   ``rank(NAS sum)`` sits 99.6% and ``rank(SAF activity)`` 99.7% inside the
    span of their own constituents' rank vectors. Their association vectors are
    therefore close to reconstructible from Stage 0's wrote per-gene table.
    D2 is kept because no substitute decides two axes against three, but it is
    recorded as a higher-resolution re-reading of the same measurements rather
    than as independent evidence, and D1 rather than D2 is the decisive check.

Nothing here is a re-read of a Stage 0 number: no partial association, no gene
count and no fibrosis or composite association vector appears in the Stage 0
result.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import rankdata
from scipy.stats import t as student_t



class PrespecificationError(RuntimeError):
    """Raised when the substrate cannot carry the question being frozen."""


PRESPEC_ID = "gse267145_axis_count_v1"

STAGE0_PRESPEC_ARTIFACTS_SHA256 = (
    "5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76"
)
STAGE0_RESULT_ARTIFACTS_SHA256 = (
    "682238189f8bd97144eb7a7d39efefea8f17546d60d25ad655498af992f60c48"
)

#: The three NASH-CRN NAS components. lobular_necrosis is not one of them.
COMPONENTS = ("steatosis", "ballooning", "lobular_inflammation")

#: The two composites, kept under distinct names because both are called
#: "activity" in clinical prose and they are not the same vector.
NAS_ACTIVITY = ("steatosis", "ballooning", "lobular_inflammation")
SAF_ACTIVITY = ("ballooning", "lobular_inflammation")

D1_FAMILIES = (
    ("nas_activity_sum", "fibrosis"),
    ("fibrosis", "nas_activity_sum"),
)
D2_FAMILIES = (
    ("steatosis", "saf_activity_sum"),
    ("saf_activity_sum", "steatosis"),
)

BH_ALPHA = 0.05

#: Two-sided 0.05 critical t at the residual df a single covariate leaves.
_critical_t = float(student_t.isf(0.025, 96))
COUNT_FLOOR_PERCENTILE = 95.0

#: Literal threshold prose. The analysis job re-reads these out of the frozen
#: file and aborts if any has drifted.
FROZEN_THRESHOLDS = {
    "d1_activity_and_fibrosis_are_independent": (
        "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
        "the matched within-sex permutation count null"
    ),
    "d2_steatosis_and_saf_activity_are_independent": (
        "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
        "the matched within-sex permutation count null"
    ),
}


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def unit_rank(values: Sequence[float]) -> np.ndarray:
    """Average ranks, centred and scaled so a dot product is a Spearman."""

    ranked = rankdata(np.asarray(values, dtype=float), method="average").astype(float)
    centred = ranked - ranked.mean()
    norm = float(np.sqrt((centred**2).sum()))
    if norm == 0.0:
        raise PrespecificationError("a constant axis cannot carry a rank correlation")
    return centred / norm


def variance_explained(target: np.ndarray, basis: Sequence[np.ndarray]) -> float:
    """Fraction of ``target``'s rank variance inside the span of ``basis``.

    This is the determination audit. A composite whose rank vector sits almost
    entirely inside the span of its constituents' rank vectors carries almost no
    label geometry those constituents do not already carry, so a criterion built
    on it is a re-reading rather than a new measurement. Reported as a number
    rather than argued.
    """

    design = np.column_stack([np.ones(target.size)] + list(basis))
    coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
    residual = target - design @ coefficients
    total = float(((target - target.mean()) ** 2).sum())
    if total == 0.0:
        raise PrespecificationError("a constant target cannot have explained variance")
    return 1.0 - float((residual**2).sum()) / total


def tie_profile(values: Sequence[float], name: str) -> dict[str, object]:
    counts = Counter(int(value) for value in values)
    n = len(list(values))
    tied = sum(count * (count - 1) // 2 for count in counts.values())
    total = n * (n - 1) // 2
    return {
        "axis": name,
        "distinct_values": len(counts),
        "value_counts": {str(key): counts[key] for key in sorted(counts)},
        "largest_tied_block": max(counts.values()),
        "fraction_of_pairs_tied": tied / total,
    }


# --------------------------------------------------------------------------
# the frozen document
# --------------------------------------------------------------------------


def build(rows: Sequence[dict[str, str]]) -> dict[str, object]:
    n = len(rows)
    if n != 99:
        raise PrespecificationError(f"expected 99 participants, found {n}")
    if len({row["participant_id"] for row in rows}) != n:
        raise PrespecificationError("participant_id is not unique")

    def column(name: str) -> np.ndarray:
        return np.asarray([int(row[name]) for row in rows], dtype=float)

    components = {name: column(name) for name in COMPONENTS}
    fibrosis = column("fibrosis")
    deposited_sum = column("nash_crn_component_sum")
    nas_activity = sum(components[name] for name in NAS_ACTIVITY)
    saf_activity = sum(components[name] for name in SAF_ACTIVITY)
    if not np.array_equal(deposited_sum, nas_activity):
        raise PrespecificationError(
            "the deposited component sum is not the three NAS components"
        )

    ranked = {name: unit_rank(values) for name, values in components.items()}
    ranked["fibrosis"] = unit_rank(fibrosis)
    ranked["nas_activity_sum"] = unit_rank(nas_activity)
    ranked["saf_activity_sum"] = unit_rank(saf_activity)

    aspect_basis = [ranked[name] for name in COMPONENTS]
    audit = {
        "what_it_measures": (
            "the fraction of each new axis's rank variance that already sits "
            "inside the span of the three aspect rank vectors Stage 0 "
            "deposited. A value near 1 means the axis adds almost no label "
            "geometry those three do not already carry."
        ),
        "rank_nas_activity_sum_in_the_span_of_its_three_components": (
            variance_explained(ranked["nas_activity_sum"], aspect_basis)
        ),
        "rank_saf_activity_sum_in_the_span_of_its_two_components": (
            variance_explained(
                ranked["saf_activity_sum"],
                [ranked[name] for name in SAF_ACTIVITY],
            )
        ),
        "rank_fibrosis_in_the_span_of_the_three_aspects": (
            variance_explained(ranked["fibrosis"], aspect_basis)
        ),
        "stage0_deposited_per_gene_vectors": list(COMPONENTS),
        "stage0_did_not_deposit_a_fibrosis_or_composite_association_vector": True,
        "reading": (
            "Both D1 directions need a per-gene fibrosis association Stage 0 "
            "never computed, and fibrosis sits less than half inside the aspect "
            "span, so D1 is a new measurement. Both composites sit above 0.99 "
            "inside their own constituents' span, so D2's statistic is close to "
            "reconstructible from Stage 0's deposited per-gene table plus one "
            "label scalar. D2 is kept because nothing else decides two axes "
            "against three, and it is recorded as a higher-resolution re-reading "
            "of the same measurements rather than as independent evidence. D1, "
            "not D2, is the decisive gate."
        ),
        "no_criterion_answer_is_on_record": (
            "No partial association, no gene count, and no fibrosis or "
            "composite association vector appears in the Stage 0 result. "
            "Reconstructible from a deposited table is not the same as already "
            "on record, and the audit above is the measured difference rather "
            "than an assertion."
        ),
    }

    label_axes = {
        "steatosis": components["steatosis"],
        "saf_activity_sum": saf_activity,
        "nas_activity_sum": nas_activity,
        "fibrosis": fibrosis,
    }
    label_spearman = {
        left: {
            right: float(unit_rank(label_axes[left]) @ unit_rank(label_axes[right]))
            for right in label_axes
        }
        for left in label_axes
    }

    sexes = Counter(row["recorded_sex"] for row in rows)

    return {
        "schema_version": "masld-bench-axis-count-prespec-v1",
        "prespec_id": PRESPEC_ID,
        "question": (
            "How many histological outcome axes does the GSE267145 bulk "
            "transcriptome support, and which grouping: one, two "
            "(activity against fibrosis), or three (SAF: steatosis, "
            "activity{ballooning + lobular inflammation}, fibrosis)?"
        ),
        "honest_provenance": {
            "this_question_was_derived_after_seeing_stage_0s_answer": True,
            "stage_0_decision": "STOP",
            "stage_0_prespecification_artifacts_sha256": (
                STAGE0_PRESPEC_ARTIFACTS_SHA256
            ),
            "stage_0_result_artifacts_sha256": STAGE0_RESULT_ARTIFACTS_SHA256,
            "stage_0_what_it_established": (
                "steatosis, ballooning and lobular inflammation are not three "
                "separable axes here: C1 failed at ballooning|inflammation "
                "0.818, C2's participation ratio was 1.42 of 3.0, and the "
                "decisive C3 gene-ordering correlations were 0.825, 0.872 and "
                "0.872 with none below 0.80. Contrast axes through the same "
                "pipeline gave recorded_sex 0.151/0.328/0.258 and fibrosis "
                "0.567/0.749/0.710, so the 0.80 bar was reachable and fibrosis "
                "reaches it."
            ),
            "criteria_fixed_before_any_expression_value_was_read": True,
            "why_the_disclosure_is_not_a_weakness": (
                "A successor question derived from a completed result is "
                "ordinary science. What separates a declared successor from a "
                "post hoc rationalisation is that the derivation is stated, the "
                "parent result is cited by digest, and the successor's criteria "
                "are fixed and frozen before its own data are opened. All three "
                "hold here and the analysis job proves the third by verifying "
                "this tree's digest before it opens the expression matrix."
            ),
            "derivable_later_is_not_recorded_before": (
                "Nothing in this file is backdated. It is frozen now, and the "
                "answer is computed afterwards by a job that must verify this "
                "digest first."
            ),
        },
        "cohort": {
            "dataset_id": "gse267145_znf469_human_liver",
            "n_participants": n,
            "unit_of_inference": "participant",
            "endpoint_table": (
                "executions/model-data-061-21079623/activation/"
                "participant_endpoints.tsv"
            ),
            "molecular_fixture": (
                "executions/model-data-064-21079902/fixture/molecular/"
            ),
        },
        "axes": {
            "steatosis": "deposited NASH-CRN steatosis grade",
            "saf_activity_sum": "ballooning + lobular_inflammation, CRN sum convention",
            "nas_activity_sum": (
                "steatosis + ballooning + lobular_inflammation, the deposited "
                "nash_crn_component_sum, identity verified against the columns"
            ),
            "fibrosis": "deposited fibrosis stage",
            "lobular_necrosis_is_excluded": (
                "deposited but not a NASH-CRN NAS component; it enters no "
                "composite and no criterion"
            ),
            "the_two_composites_are_not_the_same_vector": (
                "clinical prose calls both of them activity. nas_activity_sum "
                "includes steatosis and saf_activity_sum does not, so they are "
                "kept under distinct names throughout."
            ),
        },
        "determination_audit_run_before_the_freeze": audit,
        "what_was_inspected_before_freezing": {
            "no_gene_expression_was_opened": True,
            "stage_0_artifacts_were_read_in_full": True,
            "endpoint_label_marginals": {
                name: tie_profile(values, name)
                for name, values in {
                    **components,
                    "fibrosis": fibrosis,
                    "saf_activity_sum": saf_activity,
                    "nas_activity_sum": nas_activity,
                }.items()
            },
            "label_spearman_among_the_criterion_axes": label_spearman,
            "recorded_sex_counts": dict(sorted(sexes.items())),
            "why_these_were_looked_at": (
                "The determination audit cannot be run without the label rank "
                "geometry, and the audit is what decides whether a criterion is "
                "vacuous. Disclosing the numbers that were seen is stronger "
                "than concealing them: none of the criteria below is a function "
                "of any of them, because every criterion is a count of genes "
                "and no gene was opened."
            ),
        },
        "criteria": {
            "d1_activity_and_fibrosis_are_independent": {
                "decisive": True,
                "statistic": (
                    "per-gene partial Spearman, counted at BH 0.05 within each "
                    "direction's own family"
                ),
                "directions": [
                    {
                        "exposure": exposure,
                        "adjusted_for": covariate,
                        "family": (
                            f"every gene in the realised {exposure} | {covariate} "
                            "partial universe"
                        ),
                    }
                    for exposure, covariate in D1_FAMILIES
                ],
                "threshold": FROZEN_THRESHOLDS["d1_activity_and_fibrosis_are_independent"],
                "threshold_is_a_judgment_call": True,
                "rationale": (
                    "Two axes are independent only if each carries genes the "
                    "other does not explain. A count is required in both "
                    "directions because the relation is not symmetric: fibrosis "
                    "could be entirely explained by NAS without the reverse "
                    "holding. A bare count above zero is not enough, because BH "
                    "on a rank-based partial statistic with this much tie "
                    "structure is not guaranteed calibrated, so the count must "
                    "also clear a permutation floor built from the same "
                    "residuals."
                ),
                "why_the_95th_percentile": (
                    "It is the conventional one-sided 0.05 level and is the "
                    "same bar Stage 0 used for its familywise floors. It is a "
                    "judgment call and is recorded as one."
                ),
            },
            "d2_steatosis_and_saf_activity_are_independent": {
                "decisive": False,
                "statistic": (
                    "per-gene partial Spearman, counted at BH 0.05 within each "
                    "direction's own family"
                ),
                "directions": [
                    {
                        "exposure": exposure,
                        "adjusted_for": covariate,
                        "family": (
                            f"every gene in the realised {exposure} | {covariate} "
                            "partial universe"
                        ),
                    }
                    for exposure, covariate in D2_FAMILIES
                ],
                "threshold": FROZEN_THRESHOLDS[
                    "d2_steatosis_and_saf_activity_are_independent"
                ],
                "threshold_is_a_judgment_call": True,
                "rationale": (
                    "This decides two axes against three. It resolves a real "
                    "tension in the Stage 0 record: at label level steatosis "
                    "sits below the 0.80 bar against both other components "
                    "(0.680 and 0.728), so it is clinically distinguishable, "
                    "while at gene level all three pairs read 0.825 or above, "
                    "so it may not be transcriptionally distinguishable. D2 "
                    "settles which reading governs the model."
                ),
                "recorded_limitation": (
                    "The SAF activity composite sits 99.7% inside the span of "
                    "its two constituents' rank vectors, so D2 is a "
                    "higher-resolution re-reading of measurements Stage 0 "
                    "already deposited rather than independent evidence. Its "
                    "counts and its permutation floor are new; its inputs are "
                    "largely not."
                ),
            },
        },
        "diagnostics_never_gates": {
            "d3_does_a_composite_earn_its_place": {
                "is_a_gate": False,
                "statistic": (
                    "genes reaching BH 0.05 for a composite while reaching BH "
                    "0.05 for none of its constituents alone, in the same "
                    "universe, for nas_activity_sum and for saf_activity_sum"
                ),
                "reported_beside_it": [
                    "the reverse count: genes significant for a constituent "
                    "but not for the composite",
                    "the tie structure of every axis, because a composite has "
                    "more distinct values and fewer tied pairs than any single "
                    "component and therefore more power, so a composite-only "
                    "gene set can be a power difference rather than emergent "
                    "signal",
                    "the matched permutation count null for every family",
                    "the determination audit above, which already says the "
                    "composites add under 1% of new label rank geometry",
                ],
                "why_it_is_never_a_gate": (
                    "The comparison is confounded by power in a direction that "
                    "always favours the composite, and no reweighting removes "
                    "that cleanly at this tie structure. It is informative and "
                    "it is not decisive."
                ),
            },
            "pc1_composite_sensitivity": {
                "is_a_gate": False,
                "definition": (
                    "the first principal component of the correlation matrix of "
                    "the ranked constituents, scored per participant and used "
                    "in place of the CRN sum"
                ),
                "why_it_is_a_sensitivity_and_not_the_gate": (
                    "The CRN sum is the convention a clinician and every other "
                    "cohort in this Resource uses. A data-derived weighting is "
                    "fitted on the same 99 participants the criteria are "
                    "evaluated on, so it cannot be the gate."
                ),
            },
            "marginal_counts_for_every_axis": {
                "is_a_gate": False,
                "why": (
                    "an empty partial direction can mean not independent or "
                    "merely underpowered. The marginal count for the same axis "
                    "in the same universe separates the two readings. Fibrosis "
                    "is the axis at risk: 71 of 99 participants are stage 0 and "
                    "54.3% of participant pairs are tied, and Stage 0 measured "
                    "its split-half gene-ordering reliability at 0.357 against "
                    "0.553, 0.479 and 0.518 for the three aspects."
                ),
            },
            "gain_concentration": {
                "is_a_gate": False,
                "will_be_called_and_expected_to_return_not_applicable": True,
                "why": (
                    "gain_concentration answers where a pooled gain across "
                    "strata comes from. There is one study, one platform and "
                    "one processing run here, so there is no pooled gain, and "
                    "recorded_sex offers two strata against its min_strata of "
                    "four. It is called anyway so the not-applicable verdict is "
                    "a measured return from the package rather than an "
                    "assertion in prose. The question it would answer, whether "
                    "a count is carried by a few participants, is answered "
                    "instead by a participant bootstrap of the count itself."
                ),
            },
        },
        "method": {
            "partial_spearman_with_ties": (
                "Every variable is average-rank transformed (scipy rankdata "
                "method='average', validated against "
                "masld_bench.evaluators.metrics._average_ranks on real "
                "columns). The partial association between a gene and an "
                "exposure given a covariate is the Pearson correlation between "
                "the residual of the ranked gene and the residual of the ranked "
                "exposure, each after least-squares projection on [1, ranked "
                "covariate]. This is the rank-linear partial correlation."
            ),
            "what_it_is_not": (
                "It is not a test of conditional independence and it is not "
                "proof of an independent axis. It removes the linear-on-ranks "
                "component of one covariate and nothing else, so two things "
                "survive it: a nonlinear dependence on that covariate, and a "
                "genuine dependence on some third axis correlated with both. "
                "Fibrosis, steatosis and activity are all correlated here, so a "
                "surviving partial association says the exposure carries genes "
                "the named covariate does not explain, not that the exposure is "
                "an independent cause of them. Reported as a limitation, not "
                "corrected for."
            ),
            "per_gene_p_value": (
                "two-sided Student t on the partial correlation with df = n - 3 "
                "= 96, reported beside a permutation p from the same draws"
            ),
            "multiplicity": (
                "Benjamini-Hochberg at 0.05 within each direction's own family, "
                "using masld_bench.evaluators.stats.benjamini_hochberg. The "
                "family is named and its realised size reported beside every "
                "count. The union across the four criterion families is "
                "reported as a stricter reading and is not the gate."
            ),
            "realised_universe": (
                "re-derived from the matrix, never copied: a gene enters a "
                "family only if it is non-constant across all 99 participants "
                "and its ranked vector is not collinear with the ranked "
                "covariate. Both drop counts are reported per family."
            ),
            "count_null": (
                "the exposure residual is permuted and re-residualized on the "
                "covariate, and the whole family is recomputed, so each draw "
                "yields a complete BH count. This is a ter Braak-type residual "
                "permutation: asymptotically exact rather than exact, which is "
                "precisely why the count floor is measured instead of assumed."
            ),
            "primary_reference_is_stratified": (
                "recorded_sex is an established stratum on this substrate: "
                "Stage 0 measured ballooning +0.353 at p=0.005 and lobular "
                "inflammation +0.281 at p=0.010. The primary count null "
                "permutes within recorded_sex; the global permutation is "
                "reported beside it. outer_fold is not a stratum: Stage 0 ruled "
                "it out because the folds are balanced on histology stage by "
                "construction, so permuting within one would hold fixed part of "
                "the outcome variation being asked about."
            ),
            "p_values_are_reported_at_their_resolution_floor": True,
        },
        "power": {
            "n": 99,
            "covariates": 1,
            "residual_df": 96,
            "analytic_partial_rank_null_sd_1_over_sqrt_n_minus_2": 1.0
            / float(np.sqrt(97.0)),
            "two_sided_0_05_critical_abs_partial_r": float(
                _critical_t / np.sqrt(_critical_t**2 + 96.0)
            ),
            "note": (
                "The analytic figures are per-gene references. The "
                "multiplicity-corrected floor is the measured max-statistic "
                "permutation floor over each realised family, which accounts "
                "for gene-gene correlation as an analytic per-gene threshold "
                "cannot. A count, not a single correlation, is what the "
                "criteria read, so the count null is the operative floor."
            ),
            "the_direction_at_risk": (
                "fibrosis | nas_activity_sum. 71 of 99 participants are stage "
                "0, 54.3% of participant pairs are tied on fibrosis, and only "
                "46.5% of its rank variance sits outside the aspect span. An "
                "empty count there is reported as "
                "not_independent_or_underpowered with the measured floor and "
                "the marginal fibrosis count beside it, never as proof that "
                "fibrosis is not an axis."
            ),
        },
        "controls": {
            "unit_of_inference_is_the_participant": True,
            "every_count_carries_a_measured_permutation_floor": True,
            "participant_bootstrap_is_index_aligned_across_directions": (
                "the two directions of a criterion are compared on the same "
                "resample indices, so the difference has a paired interval "
                "rather than two marginal ones"
            ),
            "vacuity_audit": (
                "every gate condition is classified applicable or not before it "
                "is evaluated, and the verdict is reported as met over "
                "applicable, never met over total. A gate with zero applicable "
                "conditions returns NO_APPLICABLE_CONDITIONS and never a silent "
                "pass; all([]) is True and is guarded against explicitly."
            ),
            "lobular_necrosis_never_enters_a_composite_or_a_criterion": True,
            "no_model_is_fitted": True,
        },
        "decision_rule": {
            "a_one_axis_answer_is_a_result_not_a_failure": True,
            "outcome_map_is_complete_over_the_four_cells": True,
            "cells": [
                {
                    "d1": "both directions survive",
                    "d2": "both directions survive",
                    "outcome": "THREE_AXIS_SAF",
                    "meaning": (
                        "steatosis, activity{ballooning + lobular inflammation} "
                        "and fibrosis"
                    ),
                },
                {
                    "d1": "both directions survive",
                    "d2": "not both",
                    "outcome": "TWO_AXIS_ACTIVITY_FIBROSIS",
                    "meaning": (
                        "the full NAS sum against fibrosis, the standard "
                        "clinical pairing"
                    ),
                },
                {
                    "d1": "not both",
                    "d2": "both directions survive",
                    "outcome": "TWO_AXIS_STEATOSIS_ACTIVITY",
                    "meaning": (
                        "steatosis against activity{ballooning + lobular "
                        "inflammation}, with fibrosis not independent of the "
                        "NAS sum on this substrate. This cell is declared "
                        "because the two criteria can disagree and a gate with "
                        "an undeclared cell is not a gate."
                    ),
                },
                {
                    "d1": "not both",
                    "d2": "not both",
                    "outcome": "ONE_AXIS",
                    "meaning": (
                        "the surviving D1 direction names the axis: "
                        "ONE_AXIS_ACTIVITY if only nas_activity_sum | fibrosis "
                        "survives, ONE_AXIS_FIBROSIS if only fibrosis | "
                        "nas_activity_sum survives, ONE_AXIS_UNDETERMINED if "
                        "neither does. The showcase drops the aspect field and "
                        "becomes lineage-only, and every disease_aspect field "
                        "stays indeterminate rather than tested_negative."
                    ),
                },
            ],
            "no_applicable_conditions": (
                "if either criterion has zero applicable conditions the outcome "
                "is INDETERMINATE and the gate reports "
                "NO_APPLICABLE_CONDITIONS for that criterion"
            ),
        },
        "claim_boundary": (
            "This is a statement about how many outcome axes the GSE267145 "
            "bulk transcriptome resolves in these 99 participants. It is not a "
            "claim about MASLD biology, it authorises no aspect-specific "
            "attribution, and a count of separable axes is not a claim that any "
            "particular gene belongs to any particular axis."
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
    for key in FROZEN_THRESHOLDS:
        if key not in payload["criteria"]:
            raise PrespecificationError(f"missing criterion {key}")
        if payload["criteria"][key]["threshold"] != FROZEN_THRESHOLDS[key]:
            raise PrespecificationError(f"{key} threshold does not match the literal")
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(json.dumps({
        "prespec_id": payload["prespec_id"],
        "criteria": sorted(payload["criteria"]),
        "decisive": "d1_activity_and_fibrosis_are_independent",
        "outcomes": [cell["outcome"] for cell in payload["decision_rule"]["cells"]],
        "determination_audit": {
            key: round(value, 6)
            for key, value in payload["determination_audit_run_before_the_freeze"].items()
            if isinstance(value, float)
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
