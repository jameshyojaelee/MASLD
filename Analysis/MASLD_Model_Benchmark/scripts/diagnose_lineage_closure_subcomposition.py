#!/usr/bin/env python3
"""Test whether the non-hepatocyte associations survive removal of closure.

Hepatocytes are 82 to 86 percent of the composition. Because the composition is
closed to exactly 1, a falling hepatocyte share mechanically inflates every other
share: a few points of hepatocyte loss redistributes across the remaining ~15
percent and lifts a small lineage's share arithmetically, with no change in that
lineage's actual abundance.

This matters for the verdict specifically. The check's dominance guard requires at
least one qualifying lineage that is not Hepatocytes, and it exists to stop a
hepatocyte-only result. If Macrophages qualifies BECAUSE hepatocytes fell, the
guard is satisfied by exactly the output file it was written to exclude -- vacuous
in the same way the seed condition was, but harder to see.

The test: drop Hepatocytes and renormalise the remaining 15 lineages to sum to 1
within each sample, then recompute each lineage's rank association with fibrosis
stage inside that subcomposition. Log-ratio analysis is sub-compositionally
coherent, so a subcomposition is a legitimate object and this is not a second
bite at the same test.

  survives at a comparable effect -> the rise is not merely hepatocyte
                                    displacement and the finding stands
  collapses toward zero           -> the honest headline is hepatocyte depletion
                                    with the rest following arithmetically

This is a DIAGNOSTIC. It never rewrites the frozen verdict, the frozen check or
the frozen per-lineage results; it is reported alongside them.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from scripts.analyze_bulk_lineage_composition_fibrosis import (
    PRIMARY_COHORTS,
    bootstrap_interval,
    load_fixture,
    permutation_null,
    spearman,
)

DOMINANT = "Hepatocytes"


def run(*, fixture: Path, frozen_evaluation: Path, output: Path, seed: int,
        replicates: int) -> dict[str, Any]:
    if output.exists():
        raise RuntimeError(f"refusing to overwrite diagnostic: {output}")
    with (fixture / "lineage_axis.tsv").open(encoding="utf-8", newline="") as handle:
        lineages = [r["lineage"] for r in csv.DictReader(handle, delimiter="\t")]
    dominant_index = lineages.index(DOMINANT)
    kept = [l for i, l in enumerate(lineages) if i != dominant_index]

    with (frozen_evaluation / "macro_average_by_lineage.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        frozen = {
            r["lineage"]: r
            for r in csv.DictReader(handle, delimiter="\t")
            if r["arm_id"] == "three_cohort_primary"
        }

    data = {c: load_fixture(fixture, c) for c in PRIMARY_COHORTS}
    dominant_share = {
        c: float(np.mean(data[c][0][:, dominant_index])) for c in PRIMARY_COHORTS
    }

    rows: list[dict[str, Any]] = []
    macro_sub: dict[str, float] = {}
    for lineage in kept:
        index = lineages.index(lineage)
        per = []
        for cohort in PRIMARY_COHORTS:
            proportions, stage = data[cohort]
            without = np.delete(proportions, dominant_index, axis=1)
            totals = without.sum(axis=1)
            if np.any(totals <= 0):
                raise RuntimeError(f"{cohort} subcomposition has an empty sample")
            renormalised = without / totals[:, None]
            column = renormalised[:, kept.index(lineage)]
            null = permutation_null(column, stage, replicates=replicates, seed=seed + index)
            low, high, _ = bootstrap_interval(
                column, stage, replicates=replicates, seed=seed + 5077 + index
            )
            full_rho = spearman(proportions[:, index], stage)
            per.append(null["rho"])
            rows.append({
                "lineage": lineage,
                "cohort": cohort,
                "full_composition_rho": full_rho,
                "subcomposition_rho": null["rho"],
                "subcomposition_ci_low": low,
                "subcomposition_ci_high": high,
                "subcomposition_permutation_p": null["two_sided_permutation_p"],
                "attenuation": (
                    float(1.0 - abs(null["rho"]) / abs(full_rho))
                    if abs(full_rho) > 0 else float("nan")
                ),
            })
        macro_sub[lineage] = float(np.mean(per))

    output.mkdir(parents=True)
    with (output / "subcomposition_by_cohort.tsv").open("x", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    summary_rows = []
    for lineage in kept:
        full_macro = float(frozen[lineage]["macro_rho"])
        sub_macro = macro_sub[lineage]
        summary_rows.append({
            "lineage": lineage,
            "detection_floor_eligible": frozen[lineage]["detection_floor_eligible"],
            "qualifying_in_frozen_verdict": frozen[lineage]["qualifying"],
            "full_composition_macro_rho": full_macro,
            "subcomposition_macro_rho": sub_macro,
            "attenuation": (
                float(1.0 - abs(sub_macro) / abs(full_macro))
                if abs(full_macro) > 0 else float("nan")
            ),
            "sign_preserved": bool(np.sign(full_macro) == np.sign(sub_macro)),
        })
    with (output / "subcomposition_macro_summary.tsv").open("x", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary_rows)

    qualifying_non_dominant = [
        l for l in kept if frozen[l]["qualifying"] in {"True", "true", True}
    ]
    verdicts = {}
    for lineage in qualifying_non_dominant:
        full_macro = float(frozen[lineage]["macro_rho"])
        sub_macro = macro_sub[lineage]
        retained = abs(sub_macro) / abs(full_macro) if abs(full_macro) > 0 else float("nan")
        per_cohort = [r for r in rows if r["lineage"] == lineage]
        verdicts[lineage] = {
            "full_composition_macro_rho": full_macro,
            "subcomposition_macro_rho": sub_macro,
            "fraction_of_effect_retained": float(retained),
            "sign_preserved": bool(np.sign(full_macro) == np.sign(sub_macro)),
            "per_cohort_subcomposition_rho": {
                r["cohort"]: r["subcomposition_rho"] for r in per_cohort
            },
            "per_cohort_subcomposition_ci_excludes_zero": {
                r["cohort"]: bool(
                    r["subcomposition_ci_low"] > 0 or r["subcomposition_ci_high"] < 0
                )
                for r in per_cohort
            },
            "interpretation": (
                "survives the removal of closure; the association is not merely "
                "hepatocyte displacement"
                if retained >= 0.5 and np.sign(full_macro) == np.sign(sub_macro)
                else "collapses once closure is removed; the association is "
                     "substantially hepatocyte displacement expressed as a share"
            ),
        }

    receipt = {
        "schema_version": "masld-bench-lineage-closure-diagnostic-v1",
        "status": "diagnostic_only_frozen_verdict_unchanged",
        "is_a_diagnostic_not_a_rescore": True,
        "frozen_verdict_altered": False,
        "frozen_gate_altered": False,
        "frozen_evaluation": str(frozen_evaluation),
        "dominant_lineage_removed": DOMINANT,
        "dominant_mean_share_by_cohort": dominant_share,
        "subcomposition_lineages": len(kept),
        "renormalisation": "the 15 non-hepatocyte lineages rescaled to sum to 1 within each sample",
        "legitimacy": (
            "log-ratio analysis is sub-compositionally coherent, so a subcomposition is a "
            "well-defined object and this is not a second test of the same hypothesis"
        ),
        "eligibility_source": (
            "detection-floor eligibility is inherited from the frozen full-composition "
            "decision and is NOT recomputed on renormalised values, whose scale differs"
        ),
        "qualifying_non_dominant_lineages_tested": qualifying_non_dominant,
        "verdicts": verdicts,
        "replicates": replicates,
        "seed": seed,
    }
    (output / "closure_diagnostic_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--frozen-evaluation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--replicates", type=int, default=10000)
    arguments = parser.parse_args()
    receipt = run(
        fixture=arguments.fixture, frozen_evaluation=arguments.frozen_evaluation,
        output=arguments.output, seed=arguments.seed, replicates=arguments.replicates,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
