#!/usr/bin/env python3
"""Diagnose outcome-independent cCRE-free GC matching support."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from scripts.build_gc_matched_evaluation_windows import (
    INPUT_LENGTH,
    PRIMARY_CONTIGS,
    gc_fraction,
    match_gc_targets,
    merge_intervals,
    overlaps,
    percentile,
    stable_integer,
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def summarize(values: Sequence[float]) -> dict[str, float | int]:
    return {
        "pairs": len(values),
        "mean_absolute_gc_difference": sum(values) / len(values),
        "median_absolute_gc_difference": percentile(values, 0.5),
        "p95_absolute_gc_difference": percentile(values, 0.95),
        "maximum_absolute_gc_difference": max(values),
        "fraction_at_most_0_04": sum(value <= 0.04 for value in values) / len(values),
    }


def build_diagnostic(
    *,
    split_contract: Path,
    ccre_bed: Path,
    fasta: Path,
    output: Path,
    grid_steps: Iterable[int],
) -> None:
    import pyfaidx

    if output.exists():
        raise ValueError("refusing to overwrite GC diagnostic")
    windows = read_rows(split_contract / "ccre_evaluation_windows.tsv")
    fold_rows = read_rows(split_contract / "genomic_folds.tsv")
    fold_by_contig = {
        contig: int(row["genomic_fold"])
        for row in fold_rows
        for contig in row["contigs"].split(",")
    }
    forbidden_raw: dict[str, list[tuple[int, int]]] = {
        contig: [] for contig in PRIMARY_CONTIGS
    }
    with gzip.open(ccre_bed, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in forbidden_raw:
                forbidden_raw[fields[0]].append((int(fields[1]), int(fields[2])))
    forbidden = {
        contig: merge_intervals(values) for contig, values in forbidden_raw.items()
    }
    forbidden_starts = {
        contig: tuple(start for start, _end in values)
        for contig, values in forbidden.items()
    }

    reference = pyfaidx.Fasta(
        str(fasta),
        as_raw=True,
        sequence_always_upper=True,
        rebuild=False,
        read_ahead=1_000_000,
    )
    try:
        targets = []
        for row in sorted(
            windows,
            key=lambda value: (
                PRIMARY_CONTIGS.index(value["contig"]),
                int(value["input_start"]),
            ),
        ):
            start, end = int(row["input_start"]), int(row["input_end"])
            fraction = gc_fraction(str(reference[row["contig"]][start:end]))
            targets.append(
                {
                    "id": row["window_id"],
                    "contig": row["contig"],
                    "fold": fold_by_contig[row["contig"]],
                    "ccre_class": row["ccre_class"],
                    "gc": fraction,
                }
            )

        scenarios = {}
        half = INPUT_LENGTH // 2
        for grid_step in grid_steps:
            if grid_step < INPUT_LENGTH:
                raise ValueError("diagnostic grid windows must not overlap")
            candidates_by_contig: dict[str, list[tuple[float, int]]] = {}
            for contig in PRIMARY_CONTIGS:
                offset = stable_integer(
                    "background_grid", contig
                ) % grid_step
                first_center = half + offset
                candidates = []
                for center in range(first_center, len(reference[contig]) - half, grid_step):
                    start = center - half
                    end = start + INPUT_LENGTH
                    if overlaps(
                        forbidden[contig], forbidden_starts[contig], start, end
                    ):
                        continue
                    try:
                        fraction = gc_fraction(str(reference[contig][start:end]))
                    except ValueError:
                        continue
                    candidates.append((fraction, center))
                candidates_by_contig[contig] = candidates

            for matching_unit in ("chromosome", "genomic_fold"):
                targets_by_unit: dict[str, list[tuple[float, str]]] = defaultdict(list)
                candidates_by_unit: dict[str, list[tuple[float, int]]] = defaultdict(list)
                target_lookup = {row["id"]: row for row in targets}
                for row in targets:
                    unit = row["contig"] if matching_unit == "chromosome" else str(row["fold"])
                    targets_by_unit[unit].append((row["gc"], row["id"]))
                for contig, values in candidates_by_contig.items():
                    unit = contig if matching_unit == "chromosome" else str(fold_by_contig[contig])
                    for fraction, center in values:
                        encoded_center = PRIMARY_CONTIGS.index(contig) * 1_000_000_000 + center
                        candidates_by_unit[unit].append((fraction, encoded_center))
                differences = []
                by_class: dict[str, list[float]] = defaultdict(list)
                support = {}
                for unit in sorted(targets_by_unit):
                    unit_targets = targets_by_unit[unit]
                    unit_candidates = candidates_by_unit[unit]
                    support[unit] = {
                        "targets": len(unit_targets),
                        "candidates": len(unit_candidates),
                        "target_gc_min": min(value[0] for value in unit_targets),
                        "target_gc_max": max(value[0] for value in unit_targets),
                        "candidate_gc_min": min(value[0] for value in unit_candidates),
                        "candidate_gc_max": max(value[0] for value in unit_candidates),
                    }
                    matches = match_gc_targets(unit_targets, unit_candidates)
                    for target_id, (background_gc, _center) in matches.items():
                        difference = abs(target_lookup[target_id]["gc"] - background_gc)
                        differences.append(difference)
                        by_class[target_lookup[target_id]["ccre_class"]].append(difference)
                scenarios[f"{matching_unit}_grid_{grid_step}"] = {
                    "overall": summarize(differences),
                    "by_ccre_class": {
                        key: summarize(values) for key, values in sorted(by_class.items())
                    },
                    "support": support,
                }
    finally:
        reference.close()

    output.mkdir(mode=0o750)
    (output / "diagnostic.json").write_text(
        json.dumps(
            {
                "outcome_independent": True,
                "held_atac_used": False,
                "positive_windows": len(windows),
                "background_excludes_all_encode_ccres": True,
                "scenarios": scenarios,
            },
            sort_keys=True,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--ccre-bed", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grid-step", type=int, action="append", required=True)
    arguments = parser.parse_args()
    build_diagnostic(
        split_contract=arguments.split_contract,
        ccre_bed=arguments.ccre_bed,
        fasta=arguments.fasta,
        output=arguments.output,
        grid_steps=arguments.grid_step,
    )
    print(arguments.output.resolve(strict=True))


if __name__ == "__main__":
    main()
