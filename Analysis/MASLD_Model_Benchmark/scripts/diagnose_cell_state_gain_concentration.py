"""Post-hoc gain-concentration diagnostic for the 50,000-cell cell-state lane.

THIS IS EXPLICITLY POST HOC.  The frozen verdicts stand as computed under the
references they declared; nothing here revises them and no frozen output file is
modified.  A diagnostic computed after a result is seen can inform how that
result is read and can never be cited as the pre-registered test.

It answers one question the lane never asked: when a model's study-balanced
advantage is pooled across studies, does that advantage come from the studies
broadly or from one of them?
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from masld_bench.evaluators.stats import gain_concentration


class DiagnosticError(RuntimeError):
    """Raised when the inputs cannot support the diagnostic."""


def load_study_summaries(path: Path) -> dict[str, dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise DiagnosticError("study summary table is empty")
    table: dict[str, dict[str, float]] = {}
    for row in rows:
        table.setdefault(row["model_id"], {})[row["study"]] = float(
            row["donor_class_balanced_macro_f1"]
        )
    studies = {tuple(sorted(values)) for values in table.values()}
    if len(studies) != 1:
        raise DiagnosticError("models do not share one study set; cannot compare")
    return table


def run(summaries: Path, output: Path) -> dict[str, object]:
    table = load_study_summaries(summaries)
    models = sorted(table)
    studies = sorted(next(iter(table.values())))

    comparisons: list[dict[str, object]] = []
    for candidate in models:
        for baseline in models:
            if candidate == baseline:
                continue
            result = gain_concentration(
                [table[candidate][study] for study in studies],
                [table[baseline][study] for study in studies],
                studies,
            )
            payload = result.to_dict()
            payload["candidate_model_id"] = candidate
            payload["baseline_model_id"] = baseline
            comparisons.append(payload)

    applicable = [row for row in comparisons if row["applicable"]]
    flipping = [row for row in applicable if row["sign_flips_when_dropped"]]
    carriers: dict[str, int] = {}
    for row in applicable:
        stratum = str(row["most_influential_stratum"])
        carriers[stratum] = carriers.get(stratum, 0) + 1

    return {
        "schema_version": "masld-bench-gain-concentration-diagnostic-v1",
        "diagnostic_id": "cell_state_50k_gain_concentration_v1",
        "is_post_hoc": True,
        "post_hoc_statement": (
            "Computed after the lane's verdicts were frozen. It may inform how "
            "those numbers are read and may never be cited as the "
            "pre-registered test. No frozen artifact was modified."
        ),
        "revises_nothing": True,
        "supports_pass_fail_verdict": False,
        "source_table": str(summaries),
        "metric": "donor_class_balanced_macro_f1",
        "unit_of_stratification": "study",
        "n_studies": len(studies),
        "studies": studies,
        "models": models,
        "n_comparisons": len(comparisons),
        "n_applicable": len(applicable),
        "comparisons_whose_sign_flips_when_one_study_is_dropped": len(flipping),
        "most_influential_study_frequency": dict(
            sorted(carriers.items(), key=lambda item: -item[1])
        ),
        "comparisons": comparisons,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-summaries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise DiagnosticError("refusing to overwrite a diagnostic receipt")
    payload = run(arguments.study_summaries, arguments.output)
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in payload.items() if k != "comparisons"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
