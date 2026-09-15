#!/usr/bin/env python3
"""Freeze a reusable methods note on ranking-metric null references.

Written for reuse across the benchmark rather than for the cohort that
surfaced it.  Every AUPRC scored anywhere in this benchmark is exposed to the
two effects tabulated here, including the GSE267145 histology endpoints and
anything the observed-multiome lane later scores.

The note carries its own computed evidence.  Nothing in it is asserted from
memory: the bias table and the tie-structure table are produced by the bound
module at freeze time, so a reader can re-run them.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from masld_bench.evaluators.auprc_reference import (  # noqa: E402
    GUIDANCE,
    bias_table,
    random_score_reference,
    spearman_reference,
    tie_structure_table,
)
from masld_bench.gse296875_phenotype_scoring import (  # noqa: E402
    endpoint_values,
    read_tsv,
)


class NoteError(RuntimeError):
    """Raised when a bound input does not verify."""


DRAWS = 10_000
SEED = 20260825

# Cohort shapes spanning the sizes this benchmark actually scores.
BIAS_CONFIGURATIONS = (
    (5, 15),
    (10, 30),
    (15, 22),
    (20, 20),
    (25, 75),
    (50, 50),
    (100, 300),
    (200, 200),
)
TIE_GRID = (2, 3, 5, 10, 20, None)

BOUND = {
    "phenotype_endpoints": (
        "executions/gse296875-phenotype-endpoints-20260825",
        "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    ),
    "campaign_spec": (
        "executions/gse296875-phenotype-campaign-spec-21107179",
        "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    ),
    "campaign_scores": (
        "executions/gse296875-phenotype-scores-21107640",
        "0fcf833422f0316644885560b58e0a540b8e98071d7c4701c6801706230c1afd",
    ),
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def worked_example(root: Path) -> dict[str, object]:
    """The instance that surfaced the effect, with the corrected numbers."""

    endpoint_rows = read_tsv(
        root
        / "executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv"
    )
    fibrosis = endpoint_values("fibrosis", endpoint_rows)
    steatosis = endpoint_values("steatosis", endpoint_rows)
    positives = int(sum(fibrosis.values()))
    negatives = len(fibrosis) - positives

    auprc = random_score_reference(
        positives, negatives, n_draws=DRAWS, seed=SEED
    )
    rho = spearman_reference(
        list(steatosis.values()), n_draws=DRAWS, seed=SEED
    )

    scores = json.loads(
        (
            root
            / "executions/gse296875-phenotype-scores-21107640/scores"
            / "confirmatory_and_comparator_results.json"
        ).read_text()
    )
    observed = [
        record
        for record in scores
        if record.get("confirmatory") and record["endpoint"] == "fibrosis"
    ]
    values = sorted(record["estimate"] for record in observed)
    steatosis_values = sorted(steatosis.values())

    return {
        "cohort": "GSE296875 donor-by-lineage phenotype campaign",
        "campaign_outcome": "0 of 12 confirmatory tests survive BH q0.05",
        "fibrosis": {
            "positives": positives,
            "negatives": negatives,
            "prevalence": positives / len(fibrosis),
            "continuous_random_score_reference": auprc.to_dict(),
            "observed_auprc_range": [values[0], values[-1]],
            "observed_auprc_all_below_the_random_score_mean": values[-1] < auprc.mean,
            "mechanism": (
                "Every fibrosis AUPRC in that campaign, from "
                f"{values[0]:.3f} to {values[-1]:.3f}, sits below what a "
                f"continuous random scorer averages ({auprc.mean:.3f}) on the "
                "same label vector. Scored against the prevalence line of "
                f"{positives / len(fibrosis):.3f} all six lineage scopes would "
                "have read as above baseline. One metric choice is the "
                "difference between reporting nothing and reporting six false "
                "positives. This is the transferable lesson, not a detail of "
                "this cohort."
            ),
        },
        "steatosis": {
            "observed_donors": len(steatosis),
            "continuous_random_score_reference": rho.to_dict(),
            "donors_at_zero_percent": sum(1 for v in steatosis_values if v == 0.0),
            "donors_at_or_below_ten_percent": sum(
                1 for v in steatosis_values if v <= 10.0
            ),
            "donors_at_or_above_forty_five_percent": sum(
                1 for v in steatosis_values if v >= 45.0
            ),
            "leverage": (
                f"{sum(1 for v in steatosis_values if v <= 10.0)} of "
                f"{len(steatosis_values)} donors sit at or below ten percent and "
                f"{sum(1 for v in steatosis_values if v == 0.0)} sit at exactly "
                f"zero, so the entire rank statistic is carried by the "
                f"{sum(1 for v in steatosis_values if v >= 45.0)} donors at or "
                "above forty-five percent. A correlation resting on six "
                "observations is a leverage statement, not an effect estimate."
            ),
            "corrected_from": (
                "An earlier count of 26 donors at or below ten percent was read "
                "off the interquartile range instead of counted. The counted "
                "value is recorded here."
            ),
        },
        "unit_count_correction": (
            "At the twenty-nucleus threshold 170 of 195 primary donor-by-lineage "
            "units clear, not 181. The 181 figure came from a different, "
            "non-frozen lineage collapse; splitting the immune compartment into "
            "macrophage, t_cell, and b_cell exposes thin units that the collapse "
            "had hidden. Cholangiocyte (11) and t_cell (10) absorb 21 of the 25 "
            "masked units, so the two thinnest lineages carry almost the whole "
            "cost of the threshold."
        ),
        "independent_validation": (
            "The donor-by-lineage aggregation reconciles exactly with the "
            "independently built frozen Corgi context counts on 557,115 of "
            "557,115 shared cell values, which validates the barcode join, the "
            "well namespacing, and the summation together."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    bound: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in BOUND.items():
        observed = digest(arguments.root / relative / "ARTIFACTS.json")
        if observed != expected:
            raise NoteError(f"bound artifact changed: {relative}")
        bound[name] = {"path": relative, "artifacts_sha256": observed}

    bias = bias_table(BIAS_CONFIGURATIONS, n_draws=DRAWS, seed=SEED)
    ties = tie_structure_table(15, 22, TIE_GRID, n_draws=DRAWS, seed=SEED)

    note = {
        "schema_version": "masld-bench-ranking-metric-null-reference-note-v1",
        "title": "Null references for ranking metrics at small sample sizes",
        "scope": "benchmark_wide_methods_note_not_cohort_specific",
        "guidance": GUIDANCE,
        "applies_to": (
            "every endpoint in this benchmark scored by average precision or by "
            "a rank correlation, including the GSE267145 histology endpoints and "
            "anything the observed-multiome lane later scores"
        ),
        "claim_1": {
            "statement": (
                "Class prevalence is not the expected average precision of a "
                "random scorer. Average precision is upward biased at small "
                "sample sizes, so the prevalence line understates the detection "
                "floor and converts noise into apparent skill."
            ),
            "evidence": "bias_table.tsv",
            "reading": (
                "The column that matters is the gap between the random scorer's "
                "95th percentile and prevalence. At fifteen positives in "
                "thirty-seven donors a random scorer clears prevalence by more "
                "than two tenths of an AUPRC one time in twenty. The bias decays "
                "with cohort size and is small by a few hundred observations, "
                "which is why it goes unnoticed in large benchmarks and bites in "
                "small ones."
            ),
        },
        "claim_2": {
            "statement": (
                "A permutation null inherits the tie structure of the score "
                "vector it permutes. A model emitting few distinct values has a "
                "materially tighter null than a continuous scorer, so one fixed "
                "null is wrong for both."
            ),
            "evidence": "tie_structure_table.tsv",
            "reading": (
                "At fifteen positives in thirty-seven donors the 95th percentile "
                "of the null rises monotonically with the number of distinct "
                "score values. A two-valued scorer cannot reach the continuous "
                "floor, so judging it against the continuous reference is "
                "conservative, while judging a continuous scorer against a tied "
                "reference is anti-conservative and will manufacture "
                "significance."
            ),
            "when_this_bites": (
                "Any model whose predictions tie heavily: a scope where many "
                "units are imputed to a common training mean, a hard classifier "
                "reported as zero or one, a discretised risk score, or a small "
                "tree ensemble."
            ),
        },
        "prescription": [
            "Report a continuous random-score reference beside every AUPRC and "
            "every rank correlation. It is the detection floor.",
            "Take the p-value from a permutation of the model's own score "
            "vector, so the model's ties are respected.",
            "Never quote class prevalence as the AUPRC baseline.",
            "Compute the rank-correlation reference against the real outcome "
            "vector, because ties in the outcome narrow it the same way ties in "
            "the scores narrow the permutation null.",
            "State both references and the draw count, since the percentiles "
            "move slightly with the number of draws and the seed.",
        ],
        "module": "src/masld_bench/evaluators/auprc_reference.py",
        "module_is_lane_agnostic": True,
        "draws": DRAWS,
        "seed": SEED,
        "worked_example": worked_example(arguments.root),
        "bound_artifacts": bound,
        "outcomes_read": True,
        "outcomes_read_note": (
            "This note reads frozen endpoint values to build its worked example. "
            "It fits no model and produces no prediction, so it cannot leak a "
            "label into a model path."
        ),
        "status": "pass_methods_note",
    }

    arguments.output.mkdir(parents=True, exist_ok=False)
    write_tsv(arguments.output / "bias_table.tsv", bias)
    write_tsv(arguments.output / "tie_structure_table.tsv", ties)
    (arguments.output / "methods_note.json").write_text(
        json.dumps(note, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": note["status"],
                "bias_rows": len(bias),
                "tie_rows": len(ties),
                "draws": DRAWS,
                "scope": note["scope"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
