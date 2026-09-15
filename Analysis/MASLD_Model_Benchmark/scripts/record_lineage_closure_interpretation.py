#!/usr/bin/env python3
"""Bind an interpretive record to the frozen lineage-composition verdict.

The frozen verdict is PASS and it is correct as computed. Its INTERPRETATION is
not what a bare reading of it suggests, because the check's dominance guard was
satisfied by the very output file it was written to exclude.

This record is additive. It rewrites nothing: the verdict, the check and the
per-lineage results stay byte-identical, and their digests are recorded here so
that a later reader can prove this record refers to those exact numbers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(*, evaluation: Path, diagnostic: Path, gate: Path, taskspec: Path,
        output: Path) -> dict[str, Any]:
    from masld_bench.artifacts import freeze_tree

    if output.exists():
        raise RuntimeError(f"refusing to overwrite interpretive record: {output}")

    verdict = json.loads((evaluation / "promotion_gate_verdict.json").read_text(encoding="utf-8"))
    diag = json.loads((diagnostic / "closure_diagnostic_receipt.json").read_text(encoding="utf-8"))
    with (evaluation / "per_cohort_lineage_results.tsv").open(encoding="utf-8", newline="") as h:
        per = [r for r in csv.DictReader(h, delimiter="\t") if r["arm_id"] == "three_cohort_primary"]
    with (evaluation / "macro_average_by_lineage.tsv").open(encoding="utf-8", newline="") as h:
        macro = {r["lineage"]: r for r in csv.DictReader(h, delimiter="\t")
                 if r["arm_id"] == "three_cohort_primary"}

    def cohort_rows(lineage: str) -> list[dict[str, str]]:
        return [r for r in per if r["lineage"] == lineage]

    excluded = []
    for lineage, m in macro.items():
        if m["detection_floor_eligible"] in {"False", "false"}:
            rows = cohort_rows(lineage)
            excluded.append({
                "lineage": lineage,
                "macro_rho_never_tested": float(m["macro_rho"]),
                "above_floor_fraction": {r["cohort"]: float(r["fraction_above_floor"]) for r in rows},
                "minimum_above_floor_fraction": min(float(r["fraction_above_floor"]) for r in rows),
            })
    excluded.sort(key=lambda item: -abs(item["macro_rho_never_tested"]))

    mac = diag["verdicts"]["Macrophages"]
    record = {
        "schema_version": "masld-bench-interpretive-record-v1",
        "record_id": "bulk_lineage_composition_fibrosis_interpretation_v1",
        "status": "interpretive_record_additive_frozen_numbers_unchanged",
        "record_is_additive_overlay": True,
        "frozen_artifacts_this_record_interprets": {
            "evaluation_artifacts_sha256": sha256_file(evaluation / "ARTIFACTS.json"),
            "diagnostic_artifacts_sha256": sha256_file(diagnostic / "ARTIFACTS.json"),
            "gate_sha256": sha256_file(gate),
            "taskspec_sha256": sha256_file(taskspec),
        },
        "frozen_verdict_altered": False,
        "frozen_gate_altered": False,
        "frozen_per_lineage_results_altered": False,

        "the_verdict_as_computed": {
            "verdict": verdict["verdict"],
            "arithmetic": verdict["verdict_arithmetic"],
            "qualifying_lineages": verdict["qualifying_lineages"],
            "statement": (
                "The frozen verdict is PASS, 4 of 4 applicable conditions met, and it "
                "remains so as computed. Nothing in this record changes a number."
            ),
        },

        "the_interpretation": {
            "headline": (
                "Hepatocyte depletion, with the rest following arithmetically. ONE "
                "qualifying lineage, and the dominance requirement genuinely unmet."
            ),
            "qualifying_lineages_after_interpretation": ["Hepatocytes"],
            "macrophages_reclassified_as": "closure artifact, not an independent association",
            "evidence": {
                "hepatocyte_mean_share_by_cohort": diag["dominant_mean_share_by_cohort"],
                "macrophages_full_composition_macro_rho": mac["full_composition_macro_rho"],
                "macrophages_subcomposition_macro_rho": mac["subcomposition_macro_rho"],
                "macrophages_sign_preserved": mac["sign_preserved"],
                "macrophages_per_cohort_subcomposition_rho": mac["per_cohort_subcomposition_rho"],
                "macrophages_per_cohort_full_composition_rho": {
                    r["cohort"]: float(r["rho"]) for r in cohort_rows("Macrophages")
                },
                "reading": (
                    "Macrophage SHARE rises with fibrosis only because hepatocyte share "
                    "falls. Within the non-hepatocyte compartment macrophages do not rise; "
                    "the macro estimate reverses from +0.2246 to -0.1650, and in the one "
                    "cohort whose subcomposition interval excludes zero it is negative."
                ),
            },
        },

        "why_the_dominance_guard_failed": {
            "condition": "qualifying_set_dominance_guard",
            "requirement_as_written": "at least one qualifying lineage must NOT be Hepatocytes",
            "what_happened": (
                "It was satisfied by the artifact it was written to exclude. On a "
                "composition closed to exactly 1, the dominant lineage's arithmetic "
                "complement can satisfy a 'not the dominant lineage' requirement by proxy: "
                "hepatocytes at 82 to 86 percent of the composition mean that any fall in "
                "their share mechanically inflates every other share, so a small lineage "
                "can qualify without any change in its own abundance."
            ),
            "failure_shape": "a condition satisfiable by construction, reading as evidence",
            "the_contract_caught_what_the_gate_condition_did_not": (
                "TaskSpec admission gate 8 forbids reading reciprocal between-lineage shifts "
                "out of the primary, and the reported result was exactly such a shift. The "
                "written contract caught the claim; the scored gate condition did not."
            ),
            "authorship": (
                "The guard was written by this producer and argued as the reason to set the "
                "qualifying threshold at 2 rather than 1. The closure loophole was not seen "
                "at design time."
            ),
        },

        "hepatocyte_depletion_is_a_positive_control": {
            "per_cohort_rho": {r["cohort"]: float(r["rho"]) for r in cohort_rows("Hepatocytes")},
            "per_cohort_ci": {
                r["cohort"]: [float(r["ci_low"]), float(r["ci_high"])]
                for r in cohort_rows("Hepatocytes")
            },
            "all_three_ci_exclude_zero": True,
            "all_three_loco_folds_sign_consistent": True,
            "statement": (
                "Replacement of parenchyma by fibrous tissue and infiltrate is what fibrosis "
                "IS histologically, so a strong negative hepatocyte-fibrosis association is "
                "close to a positive control. Recovering it is evidence that the "
                "reference-plus-BayesPrism projection DETECTS the one compositional change "
                "independently predictable from histology."
            ),
            "this_is_a_statement_about_the_operator_not_a_discovery": True,
            "framing_rule": (
                "Reported as validation of the operator it is genuinely useful. Reported as "
                "a biological finding it would overclaim."
            ),
        },

        "what_is_testable_on_this_substrate": {
            "robustly_testable": ["Hepatocytes"],
            "testable_and_an_artifact": ["Macrophages"],
            "testable_and_null": ["pDCs"],
            "pdcs_per_cohort_rho": {r["cohort"]: float(r["rho"]) for r in cohort_rows("pDCs")},
            "statement": (
                "On this substrate only Hepatocytes is robustly testable. Macrophages is "
                "testable and turns out to be a closure artifact; pDCs is testable and null, "
                "with one cohort negative. The lane's capacity to find anything is bounded by "
                "what the operator can measure in all three cohorts, not by whether the "
                "biology is there."
            ),
        },

        "what_the_detection_floor_cost": {
            "eligible": 3,
            "excluded": len(excluded),
            "of": 16,
            "excluded_lineages": excluded,
            "the_filter_earned_its_place": (
                "Fibroblasts at macro rho 0.488 and Cholangiocytes at 0.426 would have been "
                "the two largest non-hepatocyte effects in the lane, and both rest on cohorts "
                "where the lineage sits below the detection floor in 75.6 and 96.5 percent of "
                "samples respectively. Neutrophils are below floor in 100 percent of "
                "GSE135251."
            ),
            "counterfactual": (
                "Without the filter these would have entered a 16-lineage BH family and "
                "Cholangiocytes at 0.426 would plausibly have been the headline of this lane."
            ),
        },

        "a_field_in_the_frozen_diagnostic_that_misleads_if_read_alone": {
            "field": "verdicts.Macrophages.fraction_of_effect_retained",
            "frozen_value": mac["fraction_of_effect_retained"],
            "why_it_misleads": (
                "It is |subcomposition rho| / |full-composition rho| and is MAGNITUDE ONLY. "
                "It reads as '73 percent of the effect retained' for an effect that CHANGED "
                "DIRECTION. A number that reads as reassuring when the underlying finding is "
                "a reversal is exactly the artifact this campaign keeps getting caught by."
            ),
            "corrected_name": "magnitude_ratio_ignoring_sign",
            "must_be_read_with": "sign_preserved",
            "paired_reading": {
                "magnitude_ratio_ignoring_sign": mac["fraction_of_effect_retained"],
                "sign_preserved": mac["sign_preserved"],
                "combined": (
                    f"magnitude ratio {mac['fraction_of_effect_retained']:.4f} with "
                    f"sign_preserved={mac['sign_preserved']} - the association REVERSED"
                ),
            },
            "note": (
                "The frozen diagnostic's decision logic required BOTH a magnitude ratio at or "
                "above 0.5 AND sign preservation, so it classified Macrophages correctly. The "
                "defect is in the field's name and in its separability from the sign flag, "
                "not in the decision."
            ),
        },

        "shared_estimator_dependency": {
            "required_by": "TaskSpec [shared_estimator_dependency] required_phrasing",
            "statement": (
                "All three cohorts were deconvolved against one REF_HUMAN reference, so any "
                "cross-cohort agreement reported here is evidence that ONE fixed operator "
                "transfers and is NOT evidence that three independent measurements agree. "
                "That applies to the hepatocyte result as much as to any other: the cohorts "
                "are independent in their bulk RNA and are not independent in their lineage "
                "estimates."
            ),
            "forbidden_claims_restated": [
                "the cohorts independently confirm the lineage biology",
                "cross-cohort replication validates the lineage attribution",
                "a lineage's estimated share is its measured abundance",
                "a rise in one lineage's share caused or was caused by a fall in another's",
            ],
        },

        "claim_boundary": {
            "external_development_only": True,
            "champion_eligible": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "clinical_claim_allowed": False,
            "lineage_attribution_claim_allowed": False,
        },
    }
    output.mkdir(parents=True)
    (output / "interpretive_record.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(output, {
        "artifact_class": "bulk_lineage_composition_fibrosis_interpretive_record",
        "frozen_verdict_altered": False,
        "interpreted_verdict": verdict["verdict"],
        "qualifying_lineages_after_interpretation": ["Hepatocytes"],
        "status": "passed",
    })
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", required=True, type=Path)
    parser.add_argument("--diagnostic", required=True, type=Path)
    parser.add_argument("--gate", required=True, type=Path)
    parser.add_argument("--taskspec", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    a = parser.parse_args()
    record = run(evaluation=a.evaluation, diagnostic=a.diagnostic, gate=a.gate,
                 taskspec=a.taskspec, output=a.output)
    print(json.dumps({
        "headline": record["the_interpretation"]["headline"],
        "qualifying_after_interpretation": record["the_interpretation"][
            "qualifying_lineages_after_interpretation"],
        "excluded_by_detection_floor": record["what_the_detection_floor_cost"]["excluded"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
