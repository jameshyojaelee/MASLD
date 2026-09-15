"""Post-hoc stratified permutation nulls for the completed microarray lanes.

EXPLICITLY POST HOC. The frozen verdicts stand as computed under the global
nulls they declared. Nothing here revises them and no frozen output file is
modified. A stricter null computed after a result is seen can inform how that
result is read and can never be cited as the pre-registered test.

The stratum is SCAN DATE, a genuine technical batch for a microarray assay --
shared reagent lot, instrument state and operator -- not an invented grouping.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from masld_bench.evaluators.auprc_reference import (
    permutation_reference,
    stratified_permutation_reference,
)
from masld_bench.evaluators.metrics import MetricError

DRAWS = 20_000
SEED = 20260827


class LaneDiagnosticError(RuntimeError):
    """Raised when a lane cannot support the diagnostic."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gse49541(root: Path) -> dict[str, object]:
    part = {
        row["participant_id"]: row["sample_accession"]
        for row in read_tsv(root / "executions/model-data-072-21080651-participants/gse49541_participants.tsv")
    }
    qc = {
        row["sample_accession"]: row["scan_date"][:10]
        for row in read_tsv(root / "executions/model-data-086-21110619-summarization/summary/gse49541_qc_metrics.tsv")
    }
    labels = {
        row["row_id"]: 1 if row["fibrosis_stage_group"] == "advanced_f3_f4" else 0
        for row in read_tsv(root / "executions/model-eval-866-21130184/evaluator_only/labels.tsv")
    }
    return {
        "lane": "gse49541_fibrosis_transfer",
        "endpoint": "participant_auprc_advanced_f3_f4",
        "positive_label": "advanced_f3_f4",
        "predictions_root": root / "executions/model-predict-865-21130175/predictions",
        "score_column": "probability_advanced_f3_f4",
        "labels": labels,
        "stratum_of": lambda row_id: qc.get(part.get(row_id, ""), None),
    }


def gse83452(root: Path) -> dict[str, object]:
    rec = {
        row["sample_accession"]: row["scan_date"]
        for row in read_tsv(root / "executions/model-data-072-21080651-participants/gse83452_records.tsv")
    }
    labels = {
        row["row_id"]: 1 if row["nash_status"] == "nash" else 0
        for row in read_tsv(root / "executions/model-data-888-21130317-both_arms-labels/evaluator_only/labels.tsv")
    }
    return {
        "lane": "gse83452_baseline_nash_transfer",
        "endpoint": "participant_auprc_nash_status",
        "positive_label": "nash",
        "predictions_root": root / "executions/model-predict-886-21130315-nash_vs_nafl/predictions",
        "score_column": "probability_nash",
        "labels": labels,
        "stratum_of": lambda row_id: rec.get(row_id.split("::", 1)[1], None),
    }


def score_lane(spec: dict[str, object]) -> dict[str, object]:
    predictions_root = spec["predictions_root"]
    labels = spec["labels"]
    stratum_of = spec["stratum_of"]
    models: list[dict[str, object]] = []
    for directory in sorted(p for p in predictions_root.iterdir() if p.is_dir()):
        table = read_tsv(directory / "predictions.tsv")
        rows = [row for row in table if row["row_id"] in labels]
        dropped = len(table) - len(rows)
        y = [labels[row["row_id"]] for row in rows]
        scores = [float(row[spec["score_column"]]) for row in rows]
        strata = [stratum_of(row["row_id"]) for row in rows]
        if any(value is None for value in strata):
            raise LaneDiagnosticError(f"{directory.name}: a row has no stratum")

        global_ref, observed, global_p = permutation_reference(
            y, scores, n_permutations=DRAWS, seed=SEED
        )
        entry: dict[str, object] = {
            "model_id": directory.name,
            "n_scored": len(rows),
            "rows_without_a_label_excluded": dropped,
            "observed_average_precision": observed,
            "global_null_mean": global_ref.mean,
            "global_null_p95": global_ref.percentile_95,
            "global_p_value": global_p,
        }
        try:
            strat_ref, _, strat_p = stratified_permutation_reference(
                y, scores, strata, n_permutations=DRAWS, seed=SEED
            )
        except MetricError as error:
            entry["stratified"] = "not_applicable"
            entry["stratified_reason"] = str(error)
        else:
            entry.update(
                {
                    "stratified_null_mean": strat_ref.mean,
                    "stratified_null_p95": strat_ref.percentile_95,
                    "stratified_p_value": strat_p,
                    "n_strata": strat_ref.n_strata,
                    "non_contributing_strata": strat_ref.non_contributing_strata,
                    "null_mean_shift": strat_ref.mean - global_ref.mean,
                    "p_value_shift": strat_p - global_p,
                    "conclusion_changes_at_0.05": (global_p < 0.05) != (strat_p < 0.05),
                }
            )
        models.append(entry)
    return {
        "lane": spec["lane"],
        "endpoint": spec["endpoint"],
        "stratum": "scan_date",
        "stratum_justification": (
            "Scan date is a genuine technical batch for a microarray assay - "
            "shared reagent lot, instrument state and operator. It was not "
            "chosen to produce an effect."
        ),
        "models": models,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise LaneDiagnosticError("refusing to overwrite a diagnostic receipt")
    root = arguments.root.resolve(strict=True)

    lanes = [score_lane(gse49541(root)), score_lane(gse83452(root))]
    changed = [
        f"{lane['lane']}::{model['model_id']}"
        for lane in lanes
        for model in lane["models"]
        if model.get("conclusion_changes_at_0.05")
    ]
    payload = {
        "schema_version": "masld-bench-stratified-null-diagnostic-v1",
        "diagnostic_id": "completed_lane_stratified_nulls_v1",
        "is_post_hoc": True,
        "revises_nothing": True,
        "post_hoc_statement": (
            "Computed after these lanes' verdicts were frozen under the global "
            "nulls they declared. It may inform how those p-values are read and "
            "may never be cited as the pre-registered test."
        ),
        "n_permutations": DRAWS,
        "seed": SEED,
        "models_whose_conclusion_changes_at_0.05": changed,
        "lanes": lanes,
    }
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in payload.items() if k != "lanes"}, indent=2))
    for lane in lanes:
        for model in lane["models"]:
            print(
                f"{lane['lane'][:28]:<30}{model['model_id']:<28}"
                f"AP={model['observed_average_precision']:.4f}  "
                f"global p={model['global_p_value']:.5f}  "
                f"strat p={model.get('stratified_p_value', float('nan')):.5f}  "
                f"shift={model.get('p_value_shift', float('nan')):+.5f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
