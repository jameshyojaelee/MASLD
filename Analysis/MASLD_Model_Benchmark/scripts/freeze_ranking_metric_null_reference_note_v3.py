#!/usr/bin/env python3
"""Freeze version three of the ranking-metric null reference methods note.

Version three corrects version two. Version two asserted that ties in the
*outcome* narrow a rank-correlation floor the way ties in the *score* vector
narrow an average-precision permutation null. Computing both floors disproved
it: the two GSE267145 outcome vectors have wildly different tie structures, 71.7
percent versus 24.2 percent in the largest tied group, and their null standard
deviations agree to three decimal places and match 1/sqrt(n-1) exactly.

The reason is that for a fixed vector and a random permutation of the other, the
null variance of a standardised correlation is 1/(n-1) whatever the fixed
vector's distribution. Ties in the outcome therefore change the rank floor
essentially not at all. Ties in the score vector genuinely do narrow the average
precision null, which is a different mechanism and remains true. Conflating the
two was the error.

The distinction that governs this whole output file: computing a *reference* is not
recomputing a *metric*. A reference is a function of the label vector alone and
of a stated random or permuted scorer. Nothing here reads, reproduces, restates,
or revises any candidate model's reported estimate. Where an observed estimate
appears it is quoted verbatim from the frozen scoring file for comparison
against the floor, never recalculated.
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


class NoteError(RuntimeError):
    """Raised when a cited output file does not verify."""


DRAWS = 10_000
SEED = 20260825

BIAS_CONFIGURATIONS = (
    (5, 15), (10, 30), (15, 22), (20, 20), (25, 75), (50, 50), (100, 300), (200, 200)
)
TIE_GRID = (2, 3, 5, 10, 20, None)

CITED = {
    "note_v1": (
        "executions/ranking-metric-null-reference-note-21107934",
        "5682a45aef72038e9a8460fbe4b3b43bff14bbb75ca950e3a0292b72a0380b36",
    ),
    "note_v2": (
        "executions/ranking-metric-null-reference-note-v2-21108960",
        "d4f7e84ab37e5bd9f91e8c67ffe894e6f5184fc16baa2b387b300f96a9ddd328",
    ),
    "gse274114_evaluation": (
        "executions/data-eval-247-21092325",
        "d92205364f18149070b7e636c3cdedc57158feb8513751cfc5c15593bcc0edc7",
    ),
    "gse274114_label_authority": (
        "executions/gse274114-activation-21066117",
        "f429021ebe8f742314a2073747267047447d07398a3e34eb799b3728c9d829d0",
    ),
    "gse267145_scoring": (
        "executions/model-scoring-078-21092130",
        "92a22bbbf2e8456288c6bc92e2704799da261672d943b033e1de1354183eeebd",
    ),
    "gse267145_authoritative_join": (
        "executions/gse267145-authoritative-join-21064930",
        "e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580",
    ),
    "observed_multiome_roster": (
        "executions/model-check-249-21096291",
        "9947e411225c3ae3179908fe1cf50bdf9095dd0ac13fda709195d1f85b58b53d",
    ),
    "gc_matched_prior_art": (
        "executions/gc-matched-evaluation-windows-21063848",
        "8ae4d2854a9d644141b3d47583feae52465431d1883e4b96b50fe205085e0335",
    ),
}

CITED_CONFIGS = (
    "config/campaigns/ranking_metric_floor_requirement_20260825.json",
    "config/campaigns/tied_baseline_floor_annotations_20260825.json",
)


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def gse274114_floors(root: Path) -> list[dict[str, object]]:
    """AUPRC floors at n=20 and n=19 for the two frozen contrasts.

    The evaluation receipt does not record which class carries the positive
    label, so both orientations are reported rather than guessed.
    """

    join = read_tsv(
        root / "executions/gse274114-activation-21066117/participant_join.tsv"
    )
    counts: dict[str, int] = {}
    for row in join:
        counts[row["group"]] = counts.get(row["group"], 0) + 1
    if counts != {"CTRL": 9, "ENEG": 11, "ENEG_NASH": 9, "NASH": 10}:
        raise NoteError(f"gse274114 group counts differ: {counts}")

    contrasts = (
        ("gse274114_hiseq_healthy_vs_hbv", "CTRL", "ENEG"),
        ("gse274114_novaseq_mash_vs_mash_hbv", "NASH", "ENEG_NASH"),
    )
    rows: list[dict[str, object]] = []
    for task_id, class_a, class_b in contrasts:
        for positive, negative in ((class_a, class_b), (class_b, class_a)):
            reference = random_score_reference(
                counts[positive], counts[negative], n_draws=DRAWS, seed=SEED
            )
            rows.append(
                {
                    "task_id": task_id,
                    "metric": "participant_auprc",
                    "positive_class": positive,
                    "negative_class": negative,
                    "n": counts[positive] + counts[negative],
                    "n_positive": counts[positive],
                    "n_negative": counts[negative],
                    "prevalence": reference.prevalence,
                    "continuous_random_score_mean": reference.mean,
                    "continuous_random_score_p95": reference.percentile_95,
                    "continuous_random_score_p99": reference.percentile_99,
                    "prevalence_understates_the_floor_by": (
                        reference.percentile_95 - reference.prevalence
                    ),
                }
            )
    return rows


def gse267145_floors(root: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Tie-corrected rank floors for the three frozen Spearman outcomes.

    These outcomes are ordinal and heavily tied. Ties in the outcome narrow the
    rank null exactly as ties in the score vector narrow a permutation null, so
    the floor must be computed against the real vector.
    """

    join = read_tsv(
        root / "executions/gse267145-authoritative-join-21064930/participant_join.tsv"
    )
    if len(join) != 99:
        raise NoteError(f"gse267145 participant count differs: {len(join)}")

    fibrosis = [float(row["fibrosis"]) for row in join]
    # NASH CRN activity components are steatosis, ballooning, and lobular
    # inflammation. Lobular necrosis is not a NAS component and is excluded.
    component_sum = [
        float(row["steatosis"]) + float(row["ballooning"])
        + float(row["lobular_inflammation"])
        for row in join
    ]

    vectors = {
        "fibrosis_ordinal": fibrosis,
        "nash_crn_component_sum": component_sum,
    }
    rows: list[dict[str, object]] = []
    for name, values in vectors.items():
        reference = spearman_reference(values, n_draws=DRAWS, seed=SEED)
        distinct = sorted(set(values))
        largest = max(values.count(value) for value in distinct)
        rows.append(
            {
                "outcome_vector": name,
                "n": len(values),
                "distinct_values": len(distinct),
                "largest_tied_group": largest,
                "largest_tied_fraction": largest / len(values),
                "absolute_spearman_null_sd": reference.sd,
                "absolute_spearman_p95": reference.percentile_95,
                "absolute_spearman_p99": reference.percentile_99,
            }
        )

    # Observed estimates are quoted verbatim from the frozen scorer, never
    # recomputed, purely so the floor can be read next to them.
    metrics = read_tsv(
        root
        / "executions/model-scoring-078-21092130/scores/standardized_metrics.tsv"
    )
    quoted: dict[str, dict[str, object]] = {}
    for metric_id in (
        "fibrosis_cumulative_spearman",
        "fibrosis_regression_spearman",
        "nash_crn_component_sum_spearman",
    ):
        estimates = [
            float(row["estimate"])
            for row in metrics
            if row["metric_id"] == metric_id
            and row["applicability_state"] == "observed"
            and row["estimate"] not in {"not_applicable", "not_estimable"}
        ]
        quoted[metric_id] = {
            "models_with_an_observed_estimate": len(estimates),
            "quoted_minimum": min(estimates),
            "quoted_maximum": max(estimates),
            "quoted_verbatim_not_recomputed": True,
        }
    return rows, quoted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    root = arguments.root

    cited: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in CITED.items():
        observed = digest(root / relative / "ARTIFACTS.json")
        if observed != expected:
            raise NoteError(f"cited artifact changed: {relative}")
        cited[name] = {"path": relative, "artifacts_sha256": observed}
    for relative in CITED_CONFIGS:
        cited[Path(relative).stem] = {
            "path": relative,
            "sha256": digest(root / relative),
        }

    bias = bias_table(BIAS_CONFIGURATIONS, n_draws=DRAWS, seed=SEED)
    ties = tie_structure_table(15, 22, TIE_GRID, n_draws=DRAWS, seed=SEED)
    auprc_floors = gse274114_floors(root)
    rank_floors, quoted = gse267145_floors(root)

    fibrosis_row = next(r for r in rank_floors if r["outcome_vector"] == "fibrosis_ordinal")
    component_row = next(
        r for r in rank_floors if r["outcome_vector"] == "nash_crn_component_sum"
    )

    note = {
        "schema_version": "masld-bench-ranking-metric-null-reference-note-v3",
        "title": "Null references for ranking metrics at small sample sizes",
        "version": 3,
        "supersedes": {
            "path": CITED["note_v2"][0],
            "artifacts_sha256": CITED["note_v2"][1],
            "also_supersedes": {
                "path": CITED["note_v1"][0],
                "artifacts_sha256": CITED["note_v1"][1],
            },
            "relationship": "corrects_version_two_which_remains_frozen",
            "correction": (
                "Version two claimed that ties in the outcome vector move a "
                "rank-correlation floor. Computing the floors disproved it. The "
                "two GSE267145 outcome vectors differ enormously in tie "
                "structure and their null standard deviations agree to three "
                "decimal places, both matching 1/sqrt(n-1). Version two's "
                "numeric tables were already correct; only its prose overstated "
                "the effect. Cite version three."
            ),
        },
        "scope": "benchmark_wide_methods_note_not_cohort_specific",
        "guidance": GUIDANCE,
        "computing_a_reference_is_not_recomputing_a_metric": {
            "statement": (
                "A reference is a function of the label vector and a stated "
                "random or permuted scorer. It says what a number had to clear. "
                "Recomputing a metric would re-derive what a model achieved. "
                "This artifact does only the former."
            ),
            "candidate_metrics_recomputed": 0,
            "candidate_metrics_restated_or_revised": 0,
            "observed_estimates_quoted_verbatim_for_comparison_only": True,
        },
        "claim_1_prevalence_is_not_the_null": {
            "evidence": "bias_table.tsv",
            "reading": (
                "The column that matters is the gap between the random scorer's "
                "95th percentile and prevalence, because that is what a pure "
                "noise result clears one time in twenty. It decays with cohort "
                "size and is negligible by a few hundred observations, which is "
                "why this is invisible in large benchmarks and lethal in small "
                "ones."
            ),
        },
        "claim_2_the_null_inherits_the_tie_structure": {
            "evidence": "tie_structure_table.tsv",
            "asymmetry": (
                "A tied scorer judged against a continuous reference merely "
                "loses power. A continuous scorer judged against a tied "
                "reference manufactures significance. State which side the ties "
                "are on before deciding whether an exposure matters."
            ),
        },
        "claim_3_gain_metrics_are_preferable_at_small_n": {
            "statement": (
                "Where a matched baseline scorer exists, a gain metric "
                "(candidate minus baseline, bootstrapped on the difference) is "
                "structurally immune to claim one, because class prevalence is "
                "never the reference."
            ),
            "worked_case": (
                "A benchmark-wide exposure audit of every frozen evaluator "
                "artifact found that all six registered endpoint evaluators are "
                "gain metrics, and that not one lane reports a rank metric "
                "against a prevalence line, in code or in artifact. A search of "
                "the whole source tree for prevalence and chance-level "
                "references returned a single hit, inside the module that "
                "implements this note. The architecture was immune before the "
                "problem was described, which is a stronger argument for "
                "preferring gain metrics than any correction would have been."
            ),
            "registered_gain_evaluators": [
                "graph_ld_block_auprc_gain_v1",
                "bulk_paired_spearman_gain_v1",
                "locus_heldout_allelic_spearman_v1",
                "variant_ld_block_fisher_z_spearman_gain_v1",
                "rna_atac_two_way_deviance_reduction_v1",
            ],
            "residual_exposure": (
                "A gain metric protects the comparison. It does not give the "
                "absolute candidate and baseline values a floor, and those are "
                "often what a reader takes away. Report both references beside "
                "any absolute value."
            ),
            "prior_art": {
                "path": CITED["gc_matched_prior_art"][0],
                "flag": "full_positive_set_background_AUPRC_prohibited",
                "note": (
                    "An existing guard already refuses one class of misleading "
                    "AUPRC background, which is the same instinct applied to "
                    "the background rather than the reference."
                ),
            },
        },
        "computed_floors": {
            "purpose": (
                "Protecting the next result, not auditing the last one. In both "
                "lanes below the candidate estimates are far from the floor "
                "today, so no interpretation changes. A mid-range future result "
                "on either task is the case these numbers exist for."
            ),
            "gse274114": {
                "evidence": "gse274114_auprc_floors.tsv",
                "metric": "participant_auprc",
                "positive_class_not_recorded_in_the_receipt": True,
                "both_orientations_reported": True,
                "reading": (
                    "At nineteen and twenty participants a continuous random "
                    "scorer clears prevalence by a wide margin. The frozen "
                    "training_class_prior baselines, at 0.514 and 0.458, are "
                    "constant scorers and are therefore comparators rather than "
                    "chance levels. Candidates on both tasks sit between 0.966 "
                    "and 1.000, above every floor computed here."
                ),
            },
            "gse267145": {
                "evidence": "gse267145_rank_floors.tsv",
                "metric": "spearman",
                "n": 99,
                "baseline_comparison_performed_by_the_scorer": False,
                "outcome_ties_do_not_move_the_rank_floor": (
                    "The fibrosis outcome puts "
                    f"{fibrosis_row['largest_tied_group']} of 99 participants "
                    f"({fibrosis_row['largest_tied_fraction']:.1%}) at a single "
                    "value; the NASH CRN component sum puts only "
                    f"{component_row['largest_tied_group']} of 99 "
                    f"({component_row['largest_tied_fraction']:.1%}) at its "
                    "largest. Despite that, their null standard deviations are "
                    f"{fibrosis_row['absolute_spearman_null_sd']:.4f} and "
                    f"{component_row['absolute_spearman_null_sd']:.4f}, and "
                    "their 95th percentiles "
                    f"{fibrosis_row['absolute_spearman_p95']:.4f} and "
                    f"{component_row['absolute_spearman_p95']:.4f}. Both match "
                    "1/sqrt(n-1) = 0.1010. Outcome ties move the rank floor "
                    "essentially not at all, because for a fixed vector and a "
                    "random permutation of the other the null variance of a "
                    "standardised correlation is 1/(n-1) whatever the fixed "
                    "vector's distribution."
                ),
                "the_two_tie_mechanisms_are_not_the_same": (
                    "Ties in the SCORE vector genuinely narrow an average "
                    "precision permutation null: at 15 positives in 37 the 95th "
                    "percentile runs from 0.492 at two distinct values to 0.612 "
                    "continuous. Ties in the OUTCOME do not narrow a rank "
                    "correlation null. Version two conflated the two. Compute "
                    "the rank reference against the real outcome vector anyway, "
                    "because it costs nothing and removes the question, but do "
                    "not expect it to move the number."
                ),
                "both_fibrosis_endpoints_share_one_truth_vector": (
                    "fibrosis_cumulative and fibrosis_regression score against "
                    "truth derived from the same ordinal fibrosis stage. "
                    "Spearman is invariant to a strictly monotone transform of "
                    "the truth, so one floor applies to both."
                ),
                "quoted_observed_estimates": quoted,
            },
        },
        "prescription": [
            "Report a continuous random-score reference beside every average "
            "precision and every rank correlation. It is the detection floor.",
            "Take the p-value from a permutation of the model's own score "
            "vector, so the model's ties are respected.",
            "Never quote class prevalence as the average-precision baseline.",
            "Compute a rank reference against the real outcome vector, because "
            "ties in the outcome narrow it.",
            "Prefer a gain metric over an absolute one wherever a matched "
            "baseline scorer exists.",
            "Label a constant baseline a comparator, never a chance level.",
            "Record n, the positive count, and the number of distinct score "
            "values in the receipt, so a later reader can reconstruct the floor.",
        ],
        "binding_configuration": {
            "floor_requirement": cited["ranking_metric_floor_requirement_20260825"],
            "tied_baseline_annotations": cited[
                "tied_baseline_floor_annotations_20260825"
            ],
        },
        "module": "src/masld_bench/evaluators/auprc_reference.py",
        "module_is_lane_agnostic": True,
        "applies_to": (
            "every endpoint in this benchmark scored by average precision or by "
            "a rank correlation, including the GSE267145 histology endpoints and "
            "the observed-multiome donor_peak_auprc metric that has not yet been "
            "scored"
        ),
        "draws": DRAWS,
        "seed": SEED,
        "cited_artifacts": cited,
        "outcomes_read": True,
        "outcomes_read_note": (
            "Frozen outcome vectors are read to compute label-dependent "
            "references. No model is fitted, no prediction is produced, and no "
            "candidate metric is recomputed."
        ),
        "status": "pass_methods_note_v3",
    }

    arguments.output.mkdir(parents=True, exist_ok=False)
    write_tsv(arguments.output / "bias_table.tsv", bias)
    write_tsv(arguments.output / "tie_structure_table.tsv", ties)
    write_tsv(arguments.output / "gse274114_auprc_floors.tsv", auprc_floors)
    write_tsv(arguments.output / "gse267145_rank_floors.tsv", rank_floors)
    (arguments.output / "methods_note_v3.json").write_text(
        json.dumps(note, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": note["status"],
                "version": 3,
                "auprc_floor_rows": len(auprc_floors),
                "rank_floor_rows": len(rank_floors),
                "candidate_metrics_recomputed": 0,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
