#!/usr/bin/env python3
"""Independently validate the fragment-native adjacent-stage lncRNA candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


EXPECTED = {
    "F0_to_F1": ("F1_vs_F0", 126, 187, 313, 6),
    "F1_to_F2": ("F2_vs_F1", 187, 174, 361, 6),
    "F2_to_F3": ("F3_vs_F2", 174, 132, 306, 6),
    "F3_to_F4": ("F4_vs_F3", 107, 42, 149, 5),
}
N_GENES = 24_196


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    require(path.is_file() and not path.is_symlink(), f"missing table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    require(all(None not in row for row in rows), f"malformed row width: {path}")
    return rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 1.0
    total = len(values)
    for reversed_rank, index in enumerate(reversed(order), start=1):
        rank = total - reversed_rank + 1
        running = min(running, values[index] * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.candidate.is_dir(), "candidate directory is missing")
    require(not args.candidate.is_symlink(), "candidate directory is symlinked")
    require(not args.output.exists(), "review output already exists")

    manifest = read_tsv(args.candidate / "output_manifest.tsv")
    require(len(manifest) == 9, "artifact manifest cardinality drift")
    for row in manifest:
        path = args.candidate / row["relative_path"]
        require(
            path.stat().st_size == int(row["size_bytes"]), f"size drift: {path.name}"
        )
        require(sha256(path) == row["sha256"], f"checksum drift: {path.name}")

    results = read_tsv(args.candidate / "stage_all_gene_results.tsv")
    require(len(results) == len(EXPECTED) * N_GENES, "complete result family drift")
    by_transition: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in results:
        by_transition[row["transition"]].append(row)
    require(set(by_transition) == set(EXPECTED), "transition family drift")

    recomputed_counts: dict[tuple[str, str], dict[str, int]] = {}
    for transition, rows in by_transition.items():
        contrast, n_low, n_high, n_samples, n_cohorts = EXPECTED[transition]
        require(len(rows) == N_GENES, f"tested family drift: {transition}")
        require(
            len({row["gene_id_versioned"] for row in rows}) == N_GENES,
            f"duplicate gene: {transition}",
        )
        require(
            all(
                row["contrast"] == contrast
                and int(row["n_low"]) == n_low
                and int(row["n_high"]) == n_high
                and int(row["n_samples"]) == n_samples
                and int(row["n_cohorts"]) == n_cohorts
                for row in rows
            ),
            f"contrast census drift: {transition}",
        )
        p_values = [float(row["P.Value"]) for row in rows]
        adjusted = bh(p_values)
        require(
            max(abs(value - float(row["FDR"])) for value, row in zip(adjusted, rows))
            < 1e-10,
            f"BH rederivation failed: {transition}",
        )
        for row in rows:
            require(
                all(
                    math.isfinite(float(row[field]))
                    for field in (
                        "logFC",
                        "SE",
                        "CI_low",
                        "CI_high",
                        "t",
                        "P.Value",
                        "FDR",
                        "AveExpr",
                    )
                ),
                f"non-finite statistic: {transition}",
            )
            key = (transition, row["display_biotype"])
            values = recomputed_counts.setdefault(
                key,
                {"n_tested": 0, "n_fdr_positive": 0, "n_up": 0, "n_down": 0},
            )
            values["n_tested"] += 1
            positive = float(row["FDR"]) < 0.05
            values["n_fdr_positive"] += int(positive)
            values["n_up"] += int(positive and float(row["logFC"]) > 0)
            values["n_down"] += int(positive and float(row["logFC"]) < 0)

    counts = read_tsv(args.candidate / "stage_counts_by_biotype.tsv")
    require(len(counts) == 12, "biotype count table cardinality drift")
    for row in counts:
        observed = {
            field: int(row[field])
            for field in ("n_tested", "n_fdr_positive", "n_up", "n_down")
        }
        require(
            observed == recomputed_counts[(row["transition"], row["display_biotype"])],
            f"biotype count drift: {row['transition']} {row['display_biotype']}",
        )

    lncrna = read_tsv(args.candidate / "stage_lncrna_results.tsv")
    expected_lncrna = [row for row in results if row["gene_type"] == "lncRNA"]
    require(len(lncrna) == len(expected_lncrna), "lncRNA subset cardinality drift")
    require(
        [(row["transition"], row["gene_id_versioned"]) for row in lncrna]
        == [(row["transition"], row["gene_id_versioned"]) for row in expected_lncrna],
        "lncRNA subset identity or order drift",
    )

    samples = read_tsv(args.candidate / "stage_sample_manifest.tsv")
    for transition, (_, n_low, n_high, n_samples, n_cohorts) in EXPECTED.items():
        rows = [row for row in samples if row["transition"] == transition]
        require(len(rows) == n_samples, f"sample count drift: {transition}")
        require(
            len({row["analysis_unit_id"] for row in rows}) == n_samples,
            f"participant duplication: {transition}",
        )
        require(sum(row["arm"] == "low" for row in rows) == n_low, "low-arm drift")
        require(sum(row["arm"] == "high" for row in rows) == n_high, "high-arm drift")
        require(len({row["dataset"] for row in rows}) == n_cohorts, "cohort drift")

    verdict = read_tsv(args.candidate / "stage_lncrna_verdict.tsv")
    require(len(verdict) == 4, "stage verdict cardinality drift")
    require(
        all(
            row["claim"]
            == "cross-sectional_stage-associated_remodeling_not_longitudinal_progression"
            for row in verdict
        ),
        "stage claim wording drift",
    )

    args.output.mkdir(parents=True)
    payload = {
        "status": "pass",
        "candidate": str(args.candidate.resolve()),
        "n_result_rows": len(results),
        "n_genes_per_contrast": N_GENES,
        "n_lncrna_result_rows": len(lncrna),
        "biological_unit": "one_cross-sectional_human_biopsy_per_participant",
        "multiple_testing": "BH_separately_across_complete_24196_gene_family_per_contrast",
    }
    (args.output / "VALIDATED.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
