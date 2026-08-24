#!/usr/bin/env python3
"""Report structural and QC census for reconstructed GSE281364 oligos."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--constructs", type=Path, required=True)
    arguments = parser.parse_args()
    with arguments.constructs.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_name = {row["construct"]: row for row in rows}
    alternatives = sorted(name for name in by_name if name.endswith("_Mut"))
    paired = [name for name in alternatives if name[:-4] in by_name]
    alternative_only = sorted(set(alternatives) - set(paired))
    paired_names = {name for name in paired} | {name[:-4] for name in paired}
    reference_only = sorted(
        name for name in by_name if not name.endswith("_Mut") and name not in paired_names
    )
    difference_counts: Counter[int] = Counter()
    paired_qc: Counter[str] = Counter()
    for alternative_name in paired:
        reference = by_name[alternative_name[:-4]]
        alternative = by_name[alternative_name]
        ref_sequence = reference["consensus_107bp"]
        alt_sequence = alternative["consensus_107bp"]
        difference_counts[
            sum(left != right for left, right in zip(ref_sequence, alt_sequence, strict=True))
        ] += 1
        minimum_coverage = min(
            int(reference["minimum_position_coverage"]),
            int(alternative["minimum_position_coverage"]),
        )
        minimum_fraction = min(
            float(reference["minimum_consensus_fraction"]),
            float(alternative["minimum_consensus_fraction"]),
        )
        if minimum_coverage < 20:
            paired_qc["below_20x"] += 1
        elif minimum_fraction <= 0.50:
            paired_qc["no_strict_base_majority"] += 1
        else:
            paired_qc["strict_base_majority_at_least_20x"] += 1
    output = {
        "constructs": len(rows),
        "paired_elements": len(paired),
        "reference_only_elements": len(reference_only),
        "alternative_only_elements": len(alternative_only),
        "pair_consensus_difference_count": dict(sorted(difference_counts.items())),
        "paired_consensus_qc": dict(sorted(paired_qc.items())),
        "construct_status": dict(sorted(Counter(row["status"] for row in rows).items())),
        "constructs_at_least_20x_strict_base_majority": sum(
            int(row["minimum_position_coverage"]) >= 20
            and float(row["minimum_consensus_fraction"]) > 0.50
            for row in rows
        ),
        "constructs_at_least_20x_55pct_base_majority": sum(
            int(row["minimum_position_coverage"]) >= 20
            and float(row["minimum_consensus_fraction"]) >= 0.55
            for row in rows
        ),
        "minimum_consensus_fraction_quantiles": {
            str(q): value
            for q, value in zip(
                [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0],
                __import__("numpy").quantile(
                    [float(row["minimum_consensus_fraction"]) for row in rows],
                    [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0],
                ).tolist(),
                strict=True,
            )
        },
    }
    print(json.dumps(output, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
