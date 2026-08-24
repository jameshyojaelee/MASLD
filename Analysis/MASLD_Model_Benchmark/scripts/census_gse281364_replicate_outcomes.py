#!/usr/bin/env python3
"""Report the observed schema and topology of exposed GSE281364 count files."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
import gzip
import json
from pathlib import Path


def census(path: Path) -> tuple[dict[str, object], set[str]]:
    rows = 0
    elements: dict[str, set[str]] = {}
    samples: dict[str, set[tuple[str, str]]] = {}
    sample_rows: Counter[str] = Counter()
    element_rows: Counter[str] = Counter()
    fieldnames: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows += 1
            elements.setdefault(row["element_id"], set()).add(row["allele"])
            sample_rows[row["sample_id"]] += 1
            element_rows[row["element_id"]] += 1
            samples.setdefault(row["sample_id"], set()).add(
                (row["condition"], row["replicate"])
            )
    receipt = {
        "path": str(path),
        "fieldnames": fieldnames,
        "rows": rows,
        "elements": len(elements),
        "paired_elements": sum(alleles == {"ref", "alt"} for alleles in elements.values()),
        "reference_only_elements": sum(alleles == {"ref"} for alleles in elements.values()),
        "alternative_only_elements": sum(alleles == {"alt"} for alleles in elements.values()),
        "samples": {key: sorted(value) for key, value in sorted(samples.items())},
        "rows_per_sample": dict(sorted(sample_rows.items())),
        "rows_per_element_distribution": dict(sorted(Counter(element_rows.values()).items())),
    }
    return receipt, set(elements)


def mapping_census(path: Path) -> tuple[dict[str, object], set[str]]:
    constructs: set[str] = set()
    rows = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            construct = line.split("\t", 1)[0]
            if not construct:
                raise ValueError("empty construct in mapping count file")
            constructs.add(construct)
            rows += 1
    alternatives = {value for value in constructs if value.endswith("_Mut")}
    elements = {value[:-4] if value.endswith("_Mut") else value for value in constructs}
    return {
        "path": str(path),
        "barcode_rows": rows,
        "constructs": len(constructs),
        "alternative_constructs": len(alternatives),
        "reference_or_control_constructs": len(constructs - alternatives),
        "elements": len(elements),
    }, elements


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-counts", type=Path)
    parser.add_argument("paths", type=Path, nargs="+")
    arguments = parser.parse_args()
    receipts: list[dict[str, object]] = []
    observed_sets: list[set[str]] = []
    for path in arguments.paths:
        receipt, elements = census(path)
        receipts.append(receipt)
        observed_sets.append(elements)
    output: dict[str, object] = {"observed": receipts}
    if arguments.mapping_counts is not None:
        mapping, expected = mapping_census(arguments.mapping_counts)
        mapping["per_observed_file"] = [
            {
                "observed_path": receipt["path"],
                "expected_elements_missing": len(expected - observed),
                "unexpected_observed_elements": len(observed - expected),
                "missing_examples": sorted(expected - observed)[:30],
                "unexpected_examples": sorted(observed - expected)[:30],
            }
            for receipt, observed in zip(receipts, observed_sets, strict=True)
        ]
        mapping["expected_missing_from_any_observed_file"] = len(
            expected - set.union(*observed_sets)
        )
        mapping["unexpected_in_any_observed_file"] = len(
            set.union(*observed_sets) - expected
        )
        mapping["observed_set_symmetric_difference"] = len(
            set.symmetric_difference(*observed_sets)
        )
        mapping["observed_symmetric_difference_examples"] = sorted(
            set.symmetric_difference(*observed_sets)
        )[:30]
        output["mapping"] = mapping
    print(json.dumps(output, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
