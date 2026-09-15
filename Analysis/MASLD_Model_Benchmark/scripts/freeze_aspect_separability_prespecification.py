"""Freeze what would count as three separable disease aspects, before looking.

Written and frozen BEFORE any correlation between aspects, or any gene-aspect
association, is computed. A criterion derived after seeing the answer is not a
criterion. Only the marginal label distributions -- which are needed to choose
sensible thresholds at all -- were inspected first, and they are recorded here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


class PrespecError(RuntimeError):
    """Raised when the prespecification would assert more than it should."""


def build() -> dict:
    return {
        "schema_version": "masld-bench-aspect-separability-prespec-v1",
        "prespec_id": "gse267145_aspect_separability_v1",
        "question": (
            "Are steatosis, hepatocyte ballooning and lobular inflammation "
            "separable as three distinct outcome axes in GSE267145, or do they "
            "collapse to one severity axis?"
        ),
        "why_it_gates_everything": (
            "If they collapse, an aspect-resolved model is undefined, every "
            "disease_aspect field is indeterminate, and that is the result. "
            "No later stage can rescue it."
        ),
        "cohort": {
            "dataset_id": "gse267145_znf469_human_liver",
            "n_participants": 99,
            "endpoint_table": (
                "executions/model-data-061-21079623/activation/"
                "participant_endpoints.tsv"
            ),
            "unit_of_inference": "participant",
        },
        "what_was_inspected_before_freezing": {
            "marginal_label_distributions_only": True,
            "steatosis": {"0": 24, "1": 43, "2": 19, "3": 13},
            "ballooning": {"0": 52, "1": 30, "2": 17},
            "lobular_inflammation": {"0": 40, "1": 39, "2": 20},
            "fibrosis": {"0": 71, "1": 15, "2": 9, "3": 4},
            "nash_crn_component_sum_zero_participants": 24,
            "no_cross_tabulation_or_correlation_was_computed": True,
            "no_gene_expression_was_opened": True,
        },
        "criteria": {
            "c1_labels_are_not_redundant": {
                "statistic": "pairwise Spearman among the three aspects",
                "threshold": "every pair |rho| < 0.80",
                "rationale": (
                    "Above 0.80 two graded 3-to-4 level scores carry "
                    "substantially the same ordering and are not separate "
                    "measurements for this purpose."
                ),
                "threshold_is_a_judgment_call": True,
            },
            "c2_effective_dimensionality": {
                "statistic": (
                    "participation ratio of the 3x3 Spearman matrix "
                    "eigenvalues, (sum lambda)^2 / sum lambda^2"
                ),
                "threshold": ">= 2.0 of a possible 3.0",
                "rationale": (
                    "Below 2.0 the three scores span roughly one direction and "
                    "an aspect-resolved model has nothing extra to resolve."
                ),
                "threshold_is_a_judgment_call": True,
            },
            "c3_aspects_rank_genes_differently": {
                "statistic": (
                    "Spearman between the per-gene aspect-association vectors, "
                    "for each of the three aspect pairs"
                ),
                "threshold": "at least one pair with |rho| < 0.80",
                "rationale": (
                    "This is the decisive criterion. Two correlated labels can "
                    "still carry distinguishable gene signatures, and two "
                    "weakly correlated labels can still rank genes "
                    "identically. What the model needs is different gene "
                    "orderings, not merely different scores."
                ),
                "threshold_is_a_judgment_call": True,
                "decisive": True,
            },
        },
        "decision_rule": {
            "go": "c1 AND c2 AND c3 all met",
            "stop_labels_collapse": (
                "c1 or c2 fails: the three aspects are one severity axis. Emit "
                "indeterminate for every disease_aspect field and report the "
                "collapse as the finding."
            ),
            "stop_same_gene_ordering": (
                "c3 fails: the aspects are distinct scores that rank genes the "
                "same way. An aspect-resolved model cannot beat a severity "
                "model; report indeterminate for aspect attribution."
            ),
            "a_stop_is_a_result_not_a_failure": True,
        },
        "power": {
            "n": 99,
            "analytic_rank_null_sd": "1/sqrt(n-1) = 0.10102",
            "two_null_sd": 0.20203,
            "note": (
                "A per-gene association at n=99 is a different estimand from "
                "the donor-level phenotype prediction that required about 150 "
                "donors. The multiplicity-corrected floor over the gene "
                "universe is reported with the result, not assumed here."
            ),
        },
        "controls_that_apply": [
            "measured permutation nulls, never an assumed null",
            "24 participants have an all-zero component sum; the tie structure "
            "is reported and the null is computed from the realised vectors",
            "fibrosis is analysed as a fourth axis and never folded into the "
            "three aspects",
            "lobular_necrosis is deposited but is not a NASH-CRN NAS component "
            "and is excluded from the three aspects",
        ],
        "claim_boundary": (
            "Separability is a statement about this cohort and these deposited "
            "scores. It is not a claim about MASLD biology, and it authorises "
            "no aspect-specific attribution on its own."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise PrespecError("refusing to overwrite a prespecification")
    payload = build()
    for key in ("c1_labels_are_not_redundant", "c2_effective_dimensionality",
                "c3_aspects_rank_genes_differently"):
        if key not in payload["criteria"]:
            raise PrespecError(f"missing criterion {key}")
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"prespec_id": payload["prespec_id"],
                      "criteria": sorted(payload["criteria"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
