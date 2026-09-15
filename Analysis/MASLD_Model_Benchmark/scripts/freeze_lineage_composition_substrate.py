#!/usr/bin/env python3
"""Freeze a hashed copy of the deconvolution lineage-composition substrate.

The BayesPrism proportions and the cohort metadata live outside the benchmark
tree and are read-only there. This copies the exact bytes, records their
SHA-256, and derives the analysable fixture.

Two parsing traps are handled explicitly and fail closed.

1.  Both the metadata files AND the proportions files are R
    ``write.table(row.names = TRUE)`` output. The header carries N fields and
    every data row carries N + 1, with field 0 an unnamed rowname. For the
    proportions files this means the header holds 16 lineage names and NO
    key-column name, so a reader that treats the header as a key plus 15
    lineages is off by one and silently drops a lineage.

2.  Control arms are excluded from the graded-fibrosis fixture and counted.
    A healthy control is not fibrosis stage 0.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

LINEAGES = (
    "Endothelial cells", "Hepatocytes", "Plasma cells", "T cells",
    "Cholangiocytes", "Fibroblasts", "Macrophages", "Circulating NK/NKT",
    "Resident NK", "Mono+mono derived cells", "Basophils", "B cells",
    "cDC1s", "cDC2s", "pDCs", "Neutrophils",
)
PRIMARY_COHORTS = ("GSE135251", "GSE162694", "GSE240729")
REDUCED_COHORT = "GSE213621"

JOIN_KEY = {
    "GSE135251": "Run",
    "GSE162694": "sample_id",
    "GSE240729": "sample_id",
    # Verified: the GSE213621 proportions are keyed by SRR run accession, not by
    # BioSample. BioSample joins 0 of 367; Run joins 367 of 367.
    "GSE213621": "Run",
}
STAGE_FIELD = {
    "GSE135251": "Fibrosis_stage",
    "GSE162694": "condition",
    "GSE240729": "fibrosisscore",
    "GSE213621": "fibrotic_stage",
}
# Source value -> graded stage 0-4, or None for an excluded control arm.
STAGE_MAP: dict[str, dict[str, int | None]] = {
    "GSE135251": {"0": 0, "1": 1, "2": 2, "3": 3, "4": 4},
    "GSE162694": {
        "Control": None, "NASH_F0": 0, "NASH_F1": 1,
        "NASH_F2": 2, "NASH_F3": 3, "NASH_F4": 4,
    },
    "GSE240729": {"F0": 0, "F1": 1, "F2": 2, "F3": 3, "F4": 4},
    "GSE213621": {"Control": None, "F0F1": 0, "F2": 1, "F3F4": 2},
}
CONTROL_ARM_FIELD = {"GSE135251": ("disease", "Control")}

EXPECTED_GRADED = {"GSE135251": 172, "GSE162694": 112, "GSE240729": 66, "GSE213621": 299}
EXPECTED_METADATA_ROWS = {"GSE135251": 180, "GSE162694": 143, "GSE240729": 66, "GSE213621": 367}
EXPECTED_PROPORTION_ROWS = {"GSE135251": 216, "GSE162694": 143, "GSE240729": 66, "GSE213621": 367}


class SubstrateFreezeError(RuntimeError):
    """Raised when the substrate freeze would misparse or mislabel a sample."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t",
            lineterminator="\n", extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_rowname_offset_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read R write.table(row.names = TRUE) output, failing closed on the offset."""

    lines = [l for l in path.read_text(encoding="utf-8", errors="replace").split("\n") if l != ""]
    if len(lines) < 2:
        raise SubstrateFreezeError(f"table has no data rows: {path}")
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header) + 1:
            raise SubstrateFreezeError(
                f"{path} line {index} carries {len(fields)} fields, not {len(header) + 1}; "
                f"the R rowname offset does not hold"
            )
        record = {header[i]: fields[i + 1] for i in range(len(header))}
        record["__rowname"] = fields[0]
        rows.append(record)
    return header, rows


def read_proportions(path: Path) -> tuple[list[str], np.ndarray]:
    """Read a proportions file. Header is 16 lineage names with NO key name."""

    lines = [l for l in path.read_text(encoding="utf-8").split("\n") if l != ""]
    header = lines[0].split("\t")
    if tuple(header) != LINEAGES:
        raise SubstrateFreezeError(
            f"{path} lineage axis differs from the frozen 16-lineage axis; got {len(header)} names"
        )
    keys: list[str] = []
    values: list[list[float]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header) + 1:
            raise SubstrateFreezeError(
                f"{path} line {index} carries {len(fields)} fields, not {len(header) + 1}"
            )
        keys.append(fields[0])
        values.append([float(x) for x in fields[1:]])
    matrix = np.asarray(values, dtype=np.float64)
    if matrix.shape[1] != len(LINEAGES):
        raise SubstrateFreezeError("proportions matrix does not carry 16 lineages")
    if not np.all(np.isfinite(matrix)) or np.any(matrix < 0.0):
        raise SubstrateFreezeError("proportions are not finite and nonnegative")
    sums = matrix.sum(axis=1)
    if not np.allclose(sums, 1.0, rtol=0.0, atol=1e-9):
        raise SubstrateFreezeError("proportions do not sum to 1 within 1e-9")
    if len(set(keys)) != len(keys):
        raise SubstrateFreezeError("a proportions row key repeats")
    return keys, matrix


def build_cohort(
    *, cohort: str, metadata_path: Path, proportions_path: Path
) -> dict[str, Any]:
    header, meta_rows = read_rowname_offset_tsv(metadata_path)
    if len(meta_rows) != EXPECTED_METADATA_ROWS[cohort]:
        raise SubstrateFreezeError(f"{cohort} metadata row census differs")
    keys, matrix = read_proportions(proportions_path)
    if len(keys) != EXPECTED_PROPORTION_ROWS[cohort]:
        raise SubstrateFreezeError(f"{cohort} proportions row census differs")

    join = JOIN_KEY[cohort]
    stage_field = STAGE_FIELD[cohort]
    mapping = STAGE_MAP[cohort]
    position = {key: index for index, key in enumerate(keys)}

    rows: list[dict[str, Any]] = []
    excluded_control = 0
    unmapped: collections.Counter = collections.Counter()
    missing_from_proportions = 0
    for record in meta_rows:
        key = record[join]
        if key not in position:
            missing_from_proportions += 1
            continue
        raw = record[stage_field]
        if raw not in mapping:
            unmapped[raw] += 1
            continue
        stage = mapping[raw]
        control_spec = CONTROL_ARM_FIELD.get(cohort)
        if control_spec is not None:
            field, value = control_spec
            if record.get(field) == value:
                excluded_control += 1
                continue
        if stage is None:
            excluded_control += 1
            continue
        rows.append({"sample_key": key, "row": position[key], "stage": stage, "source_value": raw})
    if unmapped:
        raise SubstrateFreezeError(f"{cohort} carries unmapped stage values: {dict(unmapped)}")
    if len(rows) != EXPECTED_GRADED[cohort]:
        raise SubstrateFreezeError(
            f"{cohort} graded census is {len(rows)}, not {EXPECTED_GRADED[cohort]}"
        )
    rows.sort(key=lambda item: item["sample_key"])
    index = np.asarray([item["row"] for item in rows], dtype=np.int64)
    return {
        "sample_rows": [
            {"sample_key": item["sample_key"], "fibrosis_stage": item["stage"],
             "source_value": item["source_value"]}
            for item in rows
        ],
        "cohort": cohort,
        "samples": len(rows),
        "stage": np.asarray([item["stage"] for item in rows], dtype=np.int64),
        "proportions": matrix[index],
        "sample_keys": [item["sample_key"] for item in rows],
        "excluded_control_arm": excluded_control,
        "metadata_rows_without_proportions": missing_from_proportions,
        "proportion_rows_without_metadata": len(keys) - len(meta_rows) + missing_from_proportions,
        "stage_distribution": {
            str(s): int(np.sum(np.asarray([r["stage"] for r in rows]) == s))
            for s in sorted({r["stage"] for r in rows})
        },
        "metadata_sha256": sha256_file(metadata_path),
        "proportions_sha256": sha256_file(proportions_path),
    }


def run(*, deconvolution_root: Path, output: Path) -> dict[str, Any]:
    from masld_bench.artifacts import freeze_tree

    if output.exists():
        raise SubstrateFreezeError(f"refusing to overwrite substrate freeze: {output}")
    output.mkdir(parents=True)
    raw = output / "raw"
    fixture = output / "fixture"
    raw.mkdir()
    fixture.mkdir()

    audit: dict[str, Any] = {}
    for cohort in PRIMARY_COHORTS + (REDUCED_COHORT,):
        metadata_path = deconvolution_root / "bulk" / cohort / f"{cohort}_metadata.tsv"
        proportions_path = (
            deconvolution_root / "results" / cohort / f"{cohort}_bayesprism_proportions.tsv"
        )
        built = build_cohort(
            cohort=cohort, metadata_path=metadata_path, proportions_path=proportions_path
        )
        shutil.copyfile(metadata_path, raw / metadata_path.name)
        shutil.copyfile(proportions_path, raw / proportions_path.name)
        if sha256_file(raw / metadata_path.name) != built["metadata_sha256"]:
            raise SubstrateFreezeError(f"{cohort} frozen metadata copy differs")
        if sha256_file(raw / proportions_path.name) != built["proportions_sha256"]:
            raise SubstrateFreezeError(f"{cohort} frozen proportions copy differs")

        np.save(fixture / f"{cohort}_proportions.npy", built["proportions"], allow_pickle=False)
        write_tsv(
            fixture / f"{cohort}_samples.tsv",
            ("sample_key", "fibrosis_stage", "source_value"),
            built["sample_rows"],
        )
        audit[cohort] = {
            key: built[key]
            for key in (
                "samples", "excluded_control_arm", "stage_distribution",
                "metadata_rows_without_proportions", "metadata_sha256", "proportions_sha256",
            )
        }

    write_tsv(
        fixture / "lineage_axis.tsv",
        ("lineage_index", "lineage"),
        [{"lineage_index": i, "lineage": name} for i, name in enumerate(LINEAGES)],
    )

    receipt = {
        "schema_version": "masld-bench-lineage-composition-substrate-freeze-v1",
        "status": "pass_substrate_frozen_hashed_copy",
        "deconvolution_root": str(deconvolution_root),
        "source_tree_is_read_only_to_this_step": True,
        "lineages": len(LINEAGES),
        "lineage_axis": list(LINEAGES),
        "rowname_offset_asserted_on_metadata_and_proportions": True,
        "proportions_header_has_no_key_column_name": True,
        "row_sums_verified_to_1_within_1e-9": True,
        "primary_cohorts": list(PRIMARY_COHORTS),
        "reduced_resolution_cohort": REDUCED_COHORT,
        "primary_total_samples": sum(audit[c]["samples"] for c in PRIMARY_COHORTS),
        "cohorts": audit,
        "control_arms_excluded_total": sum(
            audit[c]["excluded_control_arm"] for c in PRIMARY_COHORTS
        ),
        "unit": "bulk_rna_sample",
        "donor_key_exists": False,
        "music_used": False,
    }
    with (output / "substrate_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    for directory, artifact_class in (
        (raw, "lineage_composition_substrate_raw_copy"),
        (fixture, "lineage_composition_substrate_fixture"),
    ):
        freeze_tree(directory, {
            "artifact_class": artifact_class,
            "lineages": len(LINEAGES),
            "primary_total_samples": receipt["primary_total_samples"],
            "status": "passed",
        })
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deconvolution-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run(deconvolution_root=arguments.deconvolution_root, output=arguments.output)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
