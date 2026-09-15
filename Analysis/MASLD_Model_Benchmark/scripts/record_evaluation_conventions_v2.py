"""Freeze the two evaluation conventions adopted from the oFM paper.

Additive, campaign-wide, forward-only for enforcement. Modifies nothing.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


class ConventionError(RuntimeError):
    """Raised when the record would assert more than the evidence supports."""


def build(demonstration: dict[str, object], retrospective: dict[str, object]) -> dict:
    return {
        "schema_version": "masld-bench-evaluation-convention-v1",
        "record_id": "stratified_null_and_gain_concentration_v1",
        "status": "campaign_wide_convention",
        "scope": (
            "CAMPAIGN-WIDE. Enforcement is FORWARD-ONLY: it binds gates and "
            "evaluators registered after this record. It does not revise any "
            "completed lane."
        ),
        "record_is_additive_overlay": True,
        "modifies_nothing": True,
        "external_development_only": True,
        "champion_eligible": False,
        "adopted_from": {
            "citation": (
                "Vorontsov et al., A Multimodal Foundation Model for "
                "Longitudinal Patient Representation and Scalable Insight "
                "Generation in Oncology. Tempus AI. arXiv:2608.24688, "
                "25 Aug 2026."
            ),
            "what_was_taken": (
                "Two evaluation practices only. Their permutation null "
                "reshuffles treatment assignment WITHIN propensity-score bins "
                "rather than globally; and their section 6.2 uses a quantity "
                "that should be symmetric as an internal leakage control - "
                "'that the gain in C-index is similar in both the control and "
                "experimental arms is evidence that the signal is prognostic "
                "rather than a treatment effect leaking into the score'."
            ),
            "what_was_deliberately_not_taken": {
                "architecture": (
                    "The three-stage curriculum, EMA-teacher latent prediction "
                    "and sparse-autoencoder concept extraction are scale-"
                    "dependent at 1.67M patients. Porting them to a "
                    "hundreds-of-donors campaign would be the "
                    "method-across-modalities error our own rules forbid."
                ),
                "longitudinal_framing": (
                    "Their model predicts a future latent conditioned on an "
                    "anchor-day intervention and elapsed time, which requires "
                    "repeated observations per patient. This campaign is "
                    "cross-sectional and OVERALL_PLAN forbids trajectory "
                    "language. They earn longitudinal claims by construction; "
                    "we cannot."
                ),
                "their_average_precision_chance_line": (
                    "DIVERGENCE, recorded so nobody later 'corrects' our "
                    "practice toward the paper's. Their Figure 3 uses positive "
                    "prevalence as the average-precision chance line. This "
                    "project established prevalence is NOT the AP null: at "
                    "n=37 with 15 positives a random scorer averages 0.459 "
                    "against a prevalence of 0.405, and scoring against "
                    "prevalence previously produced six false positives. Their "
                    "strata go down to 20 patients, where that bias is "
                    "material. Keep the measured random-scorer reference."
                ),
            },
        },
        "convention_1_stratified_null": {
            "force": "PROCEDURAL REQUIREMENT",
            "why_procedural_not_a_threshold": (
                "It constrains the method, not the result, so it is "
                "enforceable without inventing a number."
            ),
            "the_rule": (
                "Where the evaluation set has declared strata - study, cohort, "
                "batch, platform, preservation method, genomic or LD block - "
                "the permutation p-value MUST be computed by permuting within "
                "each stratum, and the stratified and global nulls MUST be "
                "reported side by side. A global shuffle destroys the nuisance "
                "structure along with the signal and credits the observed "
                "statistic for structure the null does not have."
            ),
            "implementation": (
                "masld_bench.evaluators.auprc_reference."
                "stratified_permutation_reference for ranking metrics and "
                "stratified_spearman_reference for continuous endpoints."
            ),
            "guards": [
                "A stratum with fewer than two members, or only one distinct "
                "label, destroys no association. Counted and exposed as "
                "non_contributing_strata, never silently folded in.",
                "If every stratum is non-contributing the null collapses onto "
                "the observed value; this raises rather than returning a "
                "meaningless p-value.",
                "The two nulls carry distinct null strings so they cannot be "
                "confused inside a frozen artifact.",
            ],
            "demonstration": demonstration,
        },
        "convention_2_gain_concentration": {
            "force": "MANDATORY REPORTED DIAGNOSTIC, never pass/fail",
            "why_not_a_gate_condition": (
                "There is no defensible threshold for 'too concentrated'. "
                "Inventing one would substitute a judgment call for the "
                "evidence it exists to expose, which is the failure mode this "
                "campaign has already flagged twice."
            ),
            "the_rule": (
                "Any claim of a gain over a baseline must be accompanied by a "
                "leave-one-stratum-out concentration report. Its absence "
                "blocks publishing the gain claim; its value never blocks a "
                "gate."
            ),
            "implementation": "masld_bench.evaluators.stats.gain_concentration",
            "why_leave_one_out_and_not_heterogeneity": (
                "Directly interpretable, assumes no distribution, and is what "
                "actually caught the GSE189600 problem in the 2026-08-25 audit."
            ),
            "min_strata_guard": (
                "Below four strata concentration is not estimable and the "
                "result is returned as not_applicable with a reason and NO "
                "numbers - the same k problem that made between-cohort "
                "heterogeneity unestimable at k=3 in the W2 lane. The vacuity "
                "lesson applied at design time rather than found at review."
            ),
            "retrospective_finding": retrospective,
        },
        "relationship_to_the_zero_applicable_convention": (
            "Both conventions here can produce not_applicable outcomes, so the "
            "zero-applicable rule governs them: a gate must branch on an empty "
            "applicable set before any aggregation, and report met out of "
            "applicable rather than met out of total."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demonstration", type=Path, required=True)
    parser.add_argument("--retrospective", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise ConventionError("refusing to overwrite a convention record")

    demonstration = json.loads(arguments.demonstration.read_text())
    diagnostic = json.loads(arguments.retrospective.read_text())
    if not diagnostic.get("is_post_hoc"):
        raise ConventionError("the retrospective must declare itself post hoc")

    flips = [
        {
            "candidate": row["candidate_model_id"],
            "baseline": row["baseline_model_id"],
            "pooled_gain": row["pooled_gain"],
            "pooled_gain_without_the_carrier": row["pooled_gain_without_it"],
            "carrier_stratum": row["most_influential_stratum"],
        }
        for row in diagnostic["comparisons"]
        if row["applicable"] and row["sign_flips_when_dropped"]
    ]
    retrospective = {
        "is_post_hoc": True,
        "revises_nothing": True,
        "target": "50,000-cell cell-state lane, donor-class-balanced macro-F1",
        "n_studies": diagnostic["n_studies"],
        "n_comparisons": diagnostic["n_comparisons"],
        "comparisons_whose_sign_flips": len(flips),
        "sign_flips": flips,
        "most_influential_study_frequency": diagnostic[
            "most_influential_study_frequency"
        ],
        "reading": (
            "GSE189600 holds 149 cells from 2 donors, roughly 0.3 percent of "
            "the evaluation set, and is the most influential single study in "
            "half of all pairwise model comparisons. The elastic-net advantage "
            "over linear SVM and over logistic regression REVERSES SIGN when "
            "it is dropped. An unweighted macro-average over studies is a "
            "defensible choice for a transfer question and is not the error; "
            "the error would be quoting the resulting margin without saying "
            "what carries it."
        ),
        "what_it_does_not_do": (
            "It does not revise the frozen verdicts, which stand as computed "
            "under the references they declared. A diagnostic computed after a "
            "result is seen can inform how that result is read and can never "
            "be cited as the pre-registered test."
        ),
        "diagnostic_artifact": diagnostic["diagnostic_id"],
    }

    payload = build(demonstration, retrospective)
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"record_id": payload["record_id"], "sign_flips": len(flips)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
