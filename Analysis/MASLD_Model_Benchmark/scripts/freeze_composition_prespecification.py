"""Freeze what would count as composition-independent axis association.

Stage 3, job 2 of 3. Held-back before any composition-adjusted association exists.

**A composition-dependent gene is not an artifact.** Fibrosis genuinely involves
cell-population change, so composition is mediation, not confounding. The two
labels are ``composition_mediated`` and ``composition_independent`` and both are
real findings. The distinction is directly actionable: it tells an experimenter
whether a bulk readout suffices or a sorted assay is required. Nothing in this
stage may call a composition-dependent gene spurious, and T1's failure is named
``SET_IS_COMPOSITION_MEDIATED`` rather than anything implying the axis collapsed.

⛔ **Correction to the frozen substrate output file.** Its
``variant_sensitivity_on_the_family.what_it_shows`` reads "Only Hepatocytes and
Macrophages are stable across variants". **Three are stable, not two** --
Hepatocytes, Macrophages and pDCs -- as the same output file's own
``stable_across_variants`` and ``family_under_the_backup_variant`` lists both
record. The prose was propagated from an upstream summary without being checked
against the data beside it. The substrate output file is frozen and its build
script is hashed into that job's source manifest, so this is a labelled
correction rather than a retrofit. The substantive point is unchanged: three of
six flip and the backup family is half the size, so the family is
variant-dependent.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


class PrespecificationError(RuntimeError):
    """Raised when the substrate cannot carry the question being frozen."""


PRESPEC_ID = "composition_dependence_v1"

PARENT_DIGESTS = {
    "stage_0c_result": "5f26597ab6b6aeeaff867d54112554a33d4d0dc35be5f333e17e9010a7199fed",
    "stage_1_external_result": "daf4f02db0585f96b5605c09f4f71d6182e77d2d9bae25e263b34b1ba2e7eece",
    "stage_1_substrate": "68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab",
    "stage_3_composition_substrate": "97f6c12407ac5d093afc978a777983b14804d971963c3291d6510059f9112623",
    "sealed_instrument_sha256": (
        "5f4e483dc22f1143f1c4447f60e308a5822dd86070540527abbd429e7278e559"),
}

ARM_AXES = {"GSE135251": ("activity", "fibrosis"),
            "GSE130970": ("activity", "fibrosis"),
            "GSE193066": ("activity",)}
EXCLUSIVE_CELLS = ("activity_only", "fibrosis_only")

FROZEN_THRESHOLDS = {
    "t1_a_composition_independent_component_exists": (
        "after adjusting for the six lineage proportions in addition to the "
        "other axis, the assigned set's median absolute partial association on "
        "its own axis still exceeds the 95th percentile of a null that permutes "
        "assignment among expression-matched genes"
    ),
    "t2_per_gene_composition_dependence": (
        "a gene is composition_independent if it retains BH 0.05 significance on "
        "its own axis after the composition adjustment, composition_mediated if "
        "it loses it while the arm's set-level MDE shows the loss was "
        "detectable, and indeterminate otherwise"
    ),
}


def build(substrate: dict) -> dict:
    family = substrate["eligible_family"]["family"]
    if len(family) != 6:
        raise PrespecificationError(f"the family is {len(family)}, not the approved six")
    stable = substrate["variant_sensitivity_on_the_family"]["stable_across_variants"]
    if len(stable) != 3:
        raise PrespecificationError(
            "the correction assumes three variant-stable lineages; the substrate "
            f"records {len(stable)}")
    arms = substrate["arms"]

    return {
        "schema_version": "masld-bench-composition-prespec-v1",
        "prespec_id": PRESPEC_ID,
        "question": (
            "Does each assigned gene set retain an association with its own axis "
            "after adjusting for lineage composition, and which genes act "
            "through composition rather than beside it?"
        ),
        "the_framing_that_governs_every_reading": {
            "a_composition_dependent_gene_is_not_an_artifact": True,
            "why": (
                "fibrosis genuinely involves cell-population change, so "
                "composition is mediation, not confounding. A gene whose axis "
                "association disappears after adjustment is telling you the "
                "axis acts THROUGH composition for that gene."
            ),
            "the_two_labels_are_both_real_findings": [
                "composition_mediated", "composition_independent"],
            "why_it_is_actionable": (
                "it tells an experimenter whether a bulk readout suffices or a "
                "sorted or single-nucleus assay is required. That is the "
                "Catalog value of this stage."
            ),
            "forbidden": (
                "calling a composition-dependent gene spurious, or naming T1's "
                "failure in any way that implies the two-axis result collapsed"
            ),
        },
        "correction_to_the_frozen_substrate_artifact": {
            "field": "variant_sensitivity_on_the_family.what_it_shows",
            "what_it_says": (
                "Only Hepatocytes and Macrophages are stable across variants"),
            "what_is_true": (
                "three are stable, not two: Hepatocytes, Macrophages and pDCs"),
            "the_same_artifact_records_the_correct_list": {
                "stable_across_variants": stable,
                "family_under_the_backup_variant": substrate[
                    "variant_sensitivity_on_the_family"]["family_under_the_backup_variant"],
            },
            "how_it_happened": (
                "the prose was propagated from an upstream summary without being "
                "checked against the data beside it in the same object"),
            "why_it_is_not_retrofitted": (
                "the substrate artifact is frozen and its build script is hashed "
                "into that job's source manifest; editing either would break a "
                "sealed record for a prose error a labelled correction reaches"),
            "the_substantive_point_is_unchanged": (
                "three of six flip and the backup family is half the size, so "
                "the eligible family is variant-dependent"),
        },
        "honest_provenance": {
            "derived_after_stage_0c_and_stage_1": True,
            "parents": PARENT_DIGESTS,
            "no_composition_adjusted_association_existed_at_the_freeze": True,
            "criteria_fixed_before_any_adjusted_number_was_computed": True,
        },
        "substrate": {
            "eligible_family": family,
            "family_size": len(family),
            "derivation_rule": substrate["eligible_family"]["rule"],
            "derived_on_the_two_kallisto_arms_only": True,
            "why_not_pooled": substrate["eligible_family"]["why_not_pooled"],
            "what_the_pooled_rule_would_have_cost": substrate[
                "families_recorded_side_by_side"]["pooled_confounded_rule"],
            "superseded_family": substrate["families_recorded_side_by_side"][
                "superseded_five_from_the_star_backup_table"],
            "quantification_is_not_uniform": substrate[
                "quantification_is_not_uniform_across_the_arms"],
        },
        "covariates": {
            "adjustment_set": list(family) + ["the other axis"],
            "n_covariates": len(family) + 1,
            "residual_df": "n - 2 - 7, i.e. n - 9",
            "rank_transformed": (
                "ranks sidestep the log-ratio question entirely, which matters "
                "because the near-zero values that make log-ratios hazardous "
                "here are exactly what the detection floor excludes"),
            "ranks_do_not_rescue_a_collapsed_part": (
                "a lineage at 1e-18 still receives an ordering set by posterior "
                "noise, which is why the detection floor exists at all"),
            "conditioning_is_a_measured_negative": {
                "what_would_have_broken": (
                    "a rank-deficient covariate block produces partial "
                    "correlations against a smaller subspace than the one "
                    "named, silently"),
                "measured_per_arm": {
                    cohort: {
                        "numerical_rank": arm[
                            "family_conditioning_without_the_other_axis"]["numerical_rank"],
                        "n_covariates": arm[
                            "family_conditioning_without_the_other_axis"]["n_covariates"],
                        "condition_number": arm[
                            "family_conditioning_without_the_other_axis"]["condition_number"],
                    } for cohort, arm in arms.items()},
                "it_did_not_bite": (
                    "every block is full rank at six of six and the worst "
                    "condition number is under five, so the closure-induced "
                    "collinearity from hepatocytes near 0.9 does not degenerate "
                    "the design"),
            },
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
                "evaluated_on": list(EXCLUSIVE_CELLS),
                "failure_is_named": "SET_IS_COMPOSITION_MEDIATED",
                "what_a_failure_means": (
                    "the set acts through composition on this substrate. That "
                    "is a finding about mechanism and is never a failure of the "
                    "two-axis result."),
                "six_lineages_is_a_stricter_test_than_five": (
                    "a larger covariate set makes T1 harder to pass, so a "
                    "SET_IS_COMPOSITION_MEDIATED outcome becomes more likely "
                    "for a partly technical reason. The set-level MDE is what "
                    "separates that from a real mediation finding, and a "
                    "mediated verdict in an arm whose MDE exceeds the observed "
                    "drop is indeterminate rather than mediated."),
            },
            "t2_per_gene_composition_dependence": {
                "decisive": False,
                "is_the_deliverable": True,
                "threshold": FROZEN_THRESHOLDS["t2_per_gene_composition_dependence"],
                "classes": ["composition_independent", "composition_mediated",
                            "indeterminate"],
                "reported_as": "the fraction in each class per assigned set per arm",
                "indeterminate_where_undetectable": (
                    "same MDE guard as Stage 1's S2, because this criterion also "
                    "accepts a null when it calls a gene mediated"),
            },
        },
        "diagnostics_never_gates": {
            "sixteen_lineage_exploratory": {
                "is_a_gate": False,
                "arm": "GSE130970 only",
                "why_only_there": (
                    "it is the one cohort where fibroblasts and endothelium are "
                    "measurable at 1.000 and 0.987 against 0.227 and 0.361 in "
                    "the training cohort"),
                "reported_under_both_variants": (
                    "this is precisely where the variant changes the answer, so "
                    "a single-variant exploratory result would mislead"),
                "enters_no_outcome_cell": True,
            },
            "variant_sensitivity": {
                "is_a_gate": False,
                "arms": ["GSE135251", "GSE130970"],
                "reports_family_change_not_only_conclusion_change": True,
            },
            "gse193066_is_single_variant": (
                "a disclosure, not a defect: it serves activity only, and "
                "activity rests on lineages measured stable across variants"),
        },
        "closure": {
            "no_between_lineage_claims": (
                "a per-lineage association with an external variable is well "
                "posed; a reciprocal shift is arithmetic. Any between-lineage "
                "statement would need a subcomposition check, and the "
                "macrophage result reversed +0.225 to -0.165 under exactly that "
                "test."),
            "one_shared_operator_not_independent_measurements": (
                "cross-cohort agreement is one deconvolution operator "
                "transferring, not independent measurements. The 13_bayesprism "
                "reference path does not exist and REF_HUMAN was overridden "
                "without the value being recorded; the byte-identical 16-name "
                "roster is what confirms one label set, not the path."),
            "never_align_the_mouse_roster_by_name": (
                "different ontologies sharing 11 names"),
        },
        "controls": {
            "vacuity_audit": (
                "met over applicable, never met over total; a gate with zero "
                "applicable conditions returns NO_APPLICABLE_CONDITIONS and "
                "never a silent pass"),
            "row_sums_use_tolerance_not_equality": True,
            "ragged_offsets_are_asserted_not_assumed": True,
            "the_sealed_instrument_is_imported_not_edited": (
                "six sealed jobs hashed it; the multi-covariate extension lives "
                "beside it and is proved to reduce exactly to it at one "
                "covariate"),
            "no_model_is_fitted": True,
        },
        "decision_rule": {
            "outcome_map_is_complete": True,
            "cells": [
                {"outcome": "COMPOSITION_INDEPENDENT_COMPONENT_EXISTS",
                 "when": "T1 holds for both exclusive cells in every arm serving them"},
                {"outcome": "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY",
                 "when": "T1 holds in one arm and not the other"},
                {"outcome": "SET_IS_COMPOSITION_MEDIATED",
                 "when": ("T1 fails where the MDE shows the drop was detectable; "
                          "the axis acts through composition, which is a "
                          "mechanism finding")},
                {"outcome": "INDETERMINATE",
                 "when": ("T1 fails where the MDE does not show the drop was "
                          "detectable, or a gate reaches zero applicable "
                          "conditions")},
            ],
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
        "family": payload["substrate"]["eligible_family"],
        "criteria": sorted(payload["criteria"]),
        "outcomes": [c["outcome"] for c in payload["decision_rule"]["cells"]],
        "variant_stable_corrected_to": len(
            payload["correction_to_the_frozen_substrate_artifact"][
                "the_same_artifact_records_the_correct_list"]["stable_across_variants"]),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
