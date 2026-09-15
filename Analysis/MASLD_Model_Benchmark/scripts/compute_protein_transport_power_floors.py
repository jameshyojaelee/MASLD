"""Compute the prespecified detection floors for the protein_transport task.

This script reads the *design* — the frozen participant join and its label
distributions — and passes plain integers into `masld_bench.protein_transport_power`,
which is itself incapable of I/O. The separation is the point: this script may read a
file, the module that turns counts into thresholds may not, so no threshold can have
been informed by a model output. `tests/contract/test_protein_transport_power_blindness`
enforces that structurally by walking the module's AST.

Nothing here reads, or can read, a fitted model, a prediction, or a metric.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
from hashlib import sha256
from pathlib import Path

from masld_bench import protein_transport_power as power

FIBROSIS_CODE = {"F0": 0, "F1": 1, "F2": 2, "F3": 3, "F4": 4}
ALPHA = 0.05
OUTER_FOLDS = 5
FOLD_SEED = 20260821
DRAW_SEED = 20260825

# The prespecified secondary family. Naming it before any fit is what makes the BH
# bar a design quantity rather than a choice made once the results are visible.
SECONDARY_FAMILY = (
    "fibrosis_F2plus",
    "saf_three_state_pooled",
    "steatosis_ordinal",
    "ballooning_ge1",
    "cross_tissue_S_SCWAT",
    "cross_tissue_S_OWAT",
)


def fold_index(unit_id: str, *, seed: int, outer_folds: int) -> int:
    digest = sha256(f"{seed}\0{unit_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % outer_folds


def compute(join_path: Path, *, n_draws: int) -> dict:
    with join_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    n = len(rows)
    family_size = len(SECONDARY_FAMILY)
    bh_alpha = power.bh_single_discovery_alpha(ALPHA, family_size)

    assignment = {
        row["person_key"]: fold_index(row["person_key"], seed=FOLD_SEED, outer_folds=OUTER_FOLDS)
        for row in rows
    }
    fold_sizes = [sum(1 for k in assignment.values() if k == f) for f in range(OUTER_FOLDS)]

    nas = [int(r["nafld_activity_score"]) for r in rows]
    distinct_per_fold = []
    for f in range(OUTER_FOLDS):
        values = {int(r["nafld_activity_score"]) for r in rows if assignment[r["person_key"]] == f}
        distinct_per_fold.append(len(values))

    binaries = {
        "fibrosis_F2plus": sum(1 for r in rows if FIBROSIS_CODE[r["kleiner_fibrosis_grade"]] >= 2),
        "ballooning_ge1": sum(1 for r in rows if int(r["hepatocellular_ballooning_score"]) >= 1),
        "steatosis_ge2": sum(1 for r in rows if int(r["steatosis_score"]) >= 2),
        "NAS_ge4": sum(1 for r in rows if int(r["nafld_activity_score"]) >= 4),
        "MASLD_any": sum(1 for r in rows if r["saf_diagnosis"] != "No_MASLD"),
    }

    binary_bars = {}
    for name, positives in binaries.items():
        negatives = n - positives
        summary = power.average_precision_null_summary(
            positives, negatives, quantiles=[1 - ALPHA, 1 - bh_alpha],
            n_draws=n_draws, seed=DRAW_SEED,
        )
        binary_bars[name] = {
            "n": n,
            "positive_count": positives,
            "negative_count": negatives,
            "prevalence": round(summary["prevalence"], 4),
            "continuous_random_score_null_mean": round(summary["mean"], 4),
            "prevalence_understates_the_null_mean_by": round(
                summary["mean"] - summary["prevalence"], 4
            ),
            "auprc_operating_bar_alpha_0_05": round(summary["quantiles"][1 - ALPHA], 4),
            "auprc_operating_bar_bh_single_discovery": round(
                summary["quantiles"][1 - bh_alpha], 4
            ),
            "detectable_auroc_alpha_0_05": round(power.detectable_auroc(positives, negatives), 3),
            "detectable_auroc_bh_single_discovery": round(
                power.detectable_auroc(positives, negatives, alpha=bh_alpha), 3
            ),
            "n_draws": n_draws,
            "seed": DRAW_SEED,
        }

    scwat = sum(1 for r in rows if r["scwat_state"] == "observed")
    owat = sum(1 for r in rows if r["owat_state"] == "observed")
    scopes = {}
    for scope, size in (("S-LIVER", n), ("S-SCWAT", scwat), ("S-OWAT", owat)):
        scopes[scope] = {
            "n": size,
            "detectable_spearman_alpha_0_05": round(power.detectable_spearman(size), 3),
            "detectable_spearman_bh_single_discovery": round(
                power.detectable_spearman(size, alpha=bh_alpha), 3
            ),
        }

    kleiner = collections.Counter(r["kleiner_fibrosis_grade"] for r in rows)
    saf_per_fold = {
        level: [
            sum(
                1
                for r in rows
                if r["saf_diagnosis"] == level and assignment[r["person_key"]] == f
            )
            for f in range(OUTER_FOLDS)
        ]
        for level in sorted({r["saf_diagnosis"] for r in rows})
    }

    return {
        "schema_version": "masld-bench-protein-transport-power-floors-v1",
        "status": "PRESPECIFIED_BEFORE_ANY_FITTING",
        "computed_by": "masld_bench.protein_transport_power",
        "blindness": (
            "The threshold module and its whole import closure contain no I/O; "
            "tests/contract/test_protein_transport_power_blindness walks their AST and fails "
            "if any of them can open a file, a socket, or a subprocess."
        ),
        "alpha": ALPHA,
        "multiplicity": {
            "secondary_family": list(SECONDARY_FAMILY),
            "family_size": family_size,
            "bh_single_discovery_alpha": round(bh_alpha, 6),
            "rationale": (
                "Benjamini-Hochberg admits rank k at alpha*k/m. A lane that cannot assume its "
                "neighbours also fire must be powered for rank 1, which is alpha/m."
            ),
        },
        "split": {
            "outer_folds": OUTER_FOLDS,
            "seed": FOLD_SEED,
            "stratified": False,
            "fold_sizes": fold_sizes,
        },
        "primary": {
            "endpoint": "nafld_activity_score_ordinal_0_8",
            "n": n,
            "tested_alone_at_alpha_0_05": True,
            "in_secondary_family": False,
            "detectable_spearman_alpha_0_05": round(power.detectable_spearman(n), 3),
            "distinct_values_per_fold": distinct_per_fold,
            "estimable": power.rank_endpoint_is_estimable(distinct_per_fold),
            "estimability_criterion": "non_constant_within_every_held_out_fold",
            "levels_present": sorted(set(nas)),
        },
        "secondary_binary_bars": binary_bars,
        "cross_tissue_scopes": scopes,
        "kleiner_stratum_estimability": power.smallest_estimable_stratum(
            dict(kleiner), outer_folds=OUTER_FOLDS
        ),
        "saf_per_fold_estimability": power.classification_endpoint_per_fold_estimability(
            saf_per_fold
        ),
        "models_fit": [],
        "models_scored": [],
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--join", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--draws", type=int, default=100_000)
    args = parser.parse_args()
    result = compute(args.join, n_draws=args.draws)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
