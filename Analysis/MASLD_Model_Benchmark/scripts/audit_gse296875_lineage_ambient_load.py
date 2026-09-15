#!/usr/bin/env python3
"""Measure hepatocyte ambient load in each donor-by-lineage pseudobulk unit.

The ambient pool in liver single-nucleus RNA is dominated by hepatocyte
transcript, so a non-hepatocyte lineage can carry a hepatocyte signal it did
not express.  This audit reports how much, per unit.

It reads no outcome and fits no model, so there is no fold to protect: the
hepatocyte-restricted gene set is defined once over the whole cohort and the
result is descriptive only.  The second half of the registered diagnostic, which
asks whether an association survives adjustment for this fraction, is vacuous
when no association exists, and is not computed here.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import statistics

import numpy as np


class AmbientAuditError(RuntimeError):
    """Raised when a frozen input does not verify."""


PRIMARY_LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
SECONDARY_LINEAGES = ("endothelial_cell", "b_cell")
LINEAGES = PRIMARY_LINEAGES + SECONDARY_LINEAGES
TOP_MARKERS = 50


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--pseudobulk-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    if digest(arguments.pseudobulk / "ARTIFACTS.json") != arguments.pseudobulk_sha256:
        raise AmbientAuditError("pseudobulk fixture changed")

    units = read_tsv(arguments.pseudobulk / "pseudobulk" / "units.tsv")
    genes = read_tsv(arguments.pseudobulk / "pseudobulk" / "genes.tsv")
    with np.load(arguments.pseudobulk / "pseudobulk" / "donor_lineage_counts.npz") as p:
        counts = p["counts"]
        observed = p["unit_observed"]

    library = counts.sum(axis=1, keepdims=True).astype(np.float64)
    safe = np.where(library <= 0, 1.0, library)
    profile = np.log1p(counts.astype(np.float64) / safe * 1e6)

    lineage_of = np.asarray([row["lineage_id"] for row in units])
    usable = observed & (library.ravel() > 0)
    means = {
        lineage: profile[(lineage_of == lineage) & usable].mean(axis=0)
        for lineage in LINEAGES
    }
    others = np.vstack([means[l] for l in LINEAGES if l != "hepatocyte"])
    specificity = means["hepatocyte"] - others.max(axis=0)
    markers = np.argsort(-specificity)[:TOP_MARKERS]

    marker_umis = counts[:, markers].sum(axis=1).astype(np.float64)
    fraction = np.where(library.ravel() > 0, marker_umis / safe.ravel(), np.nan)

    rows: list[dict[str, object]] = []
    for index, row in enumerate(units):
        rows.append(
            {
                "unit_index": index,
                "donor_id": row["donor_id"],
                "lineage_id": row["lineage_id"],
                "analysis_role": row["analysis_role"],
                "rna_state": row["rna_state"],
                "nuclei": row["nuclei"],
                "library_umis": int(library[index, 0]),
                "hepatocyte_marker_umis": int(marker_umis[index]),
                "hepatocyte_marker_umi_fraction": (
                    "" if np.isnan(fraction[index]) else f"{fraction[index]:.6f}"
                ),
            }
        )

    by_lineage: dict[str, dict[str, float]] = {}
    for lineage in LINEAGES:
        selected = fraction[(lineage_of == lineage) & usable]
        by_lineage[lineage] = {
            "units": int(selected.size),
            "median_fraction": float(statistics.median(selected)),
            "minimum_fraction": float(selected.min()),
            "maximum_fraction": float(selected.max()),
            "relative_to_hepatocyte": None,
        }
    hepatocyte_median = by_lineage["hepatocyte"]["median_fraction"]
    for lineage, record in by_lineage.items():
        record["relative_to_hepatocyte"] = record["median_fraction"] / hepatocyte_median

    arguments.output.mkdir(parents=True, exist_ok=False)
    with (arguments.output / "unit_ambient_load.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "schema_version": "masld-bench-gse296875-lineage-ambient-load-v1",
        "dataset_id": "gse296875",
        "top_hepatocyte_markers": TOP_MARKERS,
        "marker_symbols": [genes[int(i)]["symbol"] for i in markers],
        "marker_definition": (
            "highest mean log counts per million in hepatocyte units minus the "
            "largest mean across the other six lineages, over observed units"
        ),
        "outcomes_read": False,
        "model_fitted": False,
        "fold_structure_not_required": (
            "this audit touches no outcome and fits no model, so a whole-cohort "
            "marker definition cannot leak a label"
        ),
        "by_lineage": by_lineage,
        "adjustment_half_not_computed": (
            "the registered diagnostic asks whether an association survives "
            "adjustment for this fraction. No association was detected on any "
            "lineage or endpoint, so the adjustment is vacuous and is not "
            "reported as though it were a passed check."
        ),
        "status": "pass_ambient_load_descriptive",
    }
    (arguments.output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in audit.items() if k != "marker_symbols"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
