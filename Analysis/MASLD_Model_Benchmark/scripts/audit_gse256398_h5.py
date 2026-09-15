#!/usr/bin/env python3
"""Audit GSE256398 human 10x H5 matrices and donor-safe source metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Iterable, Mapping

import h5py
import numpy as np


EXPECTED_GROUPS = {
    "alcohol_associated_cirrhosis": 4,
    "alcohol_associated_hepatitis": 5,
    "healthy_control": 6,
    "masld_f0": 3,
    "mash_cirrhosis": 4,
    "mash_fibrosis": 4,
}
EXPECTED_HUMAN_FILES = 26
EXPECTED_HUMAN_BYTES = 3_392_881_344
FILE_PATTERN = re.compile(
    r"^(GSM809\d{4})_(S\d+)_CB_raw_feature_bc_matrix_filtered\.h5$"
)


class GSE256398AuditError(ValueError):
    """Raised when the source differs from the frozen GSE256398 requirements."""


def digest_values(values: Iterable[str]) -> str:
    value = sha256()
    for item in values:
        value.update(item.encode("utf-8"))
        value.update(b"\n")
    return value.hexdigest()


def decode_many(values: np.ndarray) -> list[str]:
    output = []
    for item in values:
        if isinstance(item, bytes):
            output.append(item.decode("utf-8"))
        else:
            output.append(str(item))
    return output


def disease_group(title: str) -> str:
    lowered = title.lower()
    if "alcohol-associated cirrhosis" in lowered:
        return "alcohol_associated_cirrhosis"
    if "alcohol-associated hepatitis" in lowered:
        return "alcohol_associated_hepatitis"
    if "healthy control" in lowered:
        return "healthy_control"
    if "mash cirrhosis" in lowered:
        return "mash_cirrhosis"
    if "masld f0" in lowered:
        return "masld_f0"
    if "mash fibrosis" in lowered:
        return "mash_fibrosis"
    raise GSE256398AuditError(f"unrecognized human disease group: {title}")


def parse_soft(path: Path) -> list[dict[str, object]]:
    opener = gzip.open if path.suffix == ".gz" else open
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"gsm": line.split(" = ", 1)[1], "characteristics": []}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_organism_ch1 = "):
                current["organism"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                cast = current["characteristics"]
                assert isinstance(cast, list)
                cast.append(value)
    if current is not None:
        samples.append(current)
    humans = [sample for sample in samples if sample.get("organism") == "Homo sapiens"]
    if len(samples) != 30 or len(humans) != EXPECTED_HUMAN_FILES:
        raise GSE256398AuditError("SOFT sample or human donor census differs")
    rows: list[dict[str, object]] = []
    for sample in humans:
        title = str(sample.get("title", ""))
        pieces = [piece.strip() for piece in title.split(",")]
        if len(pieces) < 3 or pieces[1] != "Human":
            raise GSE256398AuditError("human sample title differs")
        metadata: dict[str, str] = {}
        for entry in sample["characteristics"]:
            key, separator, value = str(entry).partition(":")
            if separator:
                metadata[key.strip().lower()] = value.strip()
        group = disease_group(title)
        age_source = metadata.get("age", "")
        age_match = re.fullmatch(r"(\d+) year-old", age_source)
        age_years = int(age_match.group(1)) if age_match else None
        sex = metadata.get("sex", "").lower()
        if sex not in {"female", "male", ""}:
            raise GSE256398AuditError("recorded sex differs")
        rows.append(
            {
                "gsm": sample["gsm"],
                "source_sample_id": pieces[0],
                "person_key": f"gse256398:{pieces[0]}",
                "source_title": title,
                "disease_group": group,
                "development_role": (
                    "etiology_ood" if group.startswith("alcohol_") else "masld_relevant"
                ),
                "age_source": age_source,
                "age_years": "" if age_years is None else age_years,
                "age_state": "observed" if age_years is not None else "structurally_missing",
                "recorded_sex": sex,
                "sex_state": "observed" if sex else "structurally_missing",
            }
        )
    counts = {group: 0 for group in EXPECTED_GROUPS}
    for row in rows:
        counts[str(row["disease_group"])] += 1
    if counts != EXPECTED_GROUPS or len({row["person_key"] for row in rows}) != len(rows):
        raise GSE256398AuditError("disease-group or person-key census differs")
    return rows


def read_filelist(path: Path) -> dict[str, int]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    human: dict[str, int] = {}
    for row in rows:
        name = row["Name"]
        if row["#Archive/File"] == "File" and name.startswith("GSM809"):
            if not FILE_PATTERN.fullmatch(name) or name in human:
                raise GSE256398AuditError("human H5 file-list identity differs")
            human[name] = int(row["Size"])
    if len(human) != EXPECTED_HUMAN_FILES or sum(human.values()) != EXPECTED_HUMAN_BYTES:
        raise GSE256398AuditError("human H5 file-list census differs")
    return human


def audit_h5(path: Path) -> dict[str, object]:
    match = FILE_PATTERN.fullmatch(path.name)
    if match is None:
        raise GSE256398AuditError("H5 filename differs")
    with h5py.File(path, "r") as handle:
        if "matrix" not in handle:
            raise GSE256398AuditError("10x matrix group is absent")
        matrix = handle["matrix"]
        required = {"barcodes", "data", "indices", "indptr", "shape", "features"}
        if not required.issubset(matrix.keys()):
            raise GSE256398AuditError("10x matrix datasets differ")
        shape = tuple(int(value) for value in matrix["shape"][...])
        if len(shape) != 2 or min(shape) <= 0:
            raise GSE256398AuditError("10x matrix shape differs")
        features = matrix["features"]
        if not {"id", "name", "feature_type", "genome"}.issubset(features.keys()):
            raise GSE256398AuditError("10x feature datasets differ")
        ids = decode_many(features["id"][...])
        names = decode_many(features["name"][...])
        feature_types = decode_many(features["feature_type"][...])
        genomes = decode_many(features["genome"][...])
        barcodes = decode_many(matrix["barcodes"][...])
        if (
            len(ids) != shape[0]
            or len(names) != shape[0]
            or len(feature_types) != shape[0]
            or len(genomes) != shape[0]
            or len(barcodes) != shape[1]
            or len(set(barcodes)) != len(barcodes)
            or matrix["indptr"].shape != (shape[1] + 1,)
            or matrix["data"].shape != matrix["indices"].shape
        ):
            raise GSE256398AuditError("10x matrix axis census differs")
        if set(feature_types) != {"Gene Expression"} or set(genomes) != {"GRCh38"}:
            raise GSE256398AuditError("10x feature type or genome differs")
        data = matrix["data"]
        indices = matrix["indices"]
        sample_positions = sorted({0, max(0, data.shape[0] // 2), max(0, data.shape[0] - 1)})
        values = np.asarray([data[position] for position in sample_positions])
        index_values = np.asarray([indices[position] for position in sample_positions])
        if (
            not np.issubdtype(data.dtype, np.integer)
            or np.any(values <= 0)
            or np.any(index_values < 0)
            or np.any(index_values >= shape[0])
        ):
            raise GSE256398AuditError("10x count or sparse-index sample differs")
        return {
            "gsm": match.group(1),
            "source_sample_id": match.group(2),
            "file": path.name,
            "bytes": path.stat().st_size,
            "features": shape[0],
            "nuclei": shape[1],
            "nonzero_entries": int(data.shape[0]),
            "feature_ids_sha256": digest_values(ids),
            "feature_names_sha256": digest_values(names),
            "barcodes_sha256": digest_values(barcodes),
            "genome": "GRCh38",
            "count_dtype": str(data.dtype),
        }


def write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series-soft", type=Path, required=True)
    parser.add_argument("--filelist", type=Path, required=True)
    parser.add_argument("--h5-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise GSE256398AuditError("audit output exists")
    samples = parse_soft(args.series_soft)
    expected = read_filelist(args.filelist)
    observed_paths = {path.name: path for path in args.h5_root.glob("*.h5")}
    if set(observed_paths) != set(expected):
        raise GSE256398AuditError("extracted human H5 roster differs")
    if any(path.stat().st_size != expected[name] for name, path in observed_paths.items()):
        raise GSE256398AuditError("extracted human H5 size differs")
    audits = [audit_h5(observed_paths[name]) for name in sorted(observed_paths)]
    by_gsm = {str(row["gsm"]): row for row in audits}
    if len(by_gsm) != EXPECTED_HUMAN_FILES:
        raise GSE256398AuditError("H5 GSM census differs")
    joined: list[dict[str, object]] = []
    for sample in samples:
        h5 = by_gsm.get(str(sample["gsm"]))
        if h5 is None or h5["source_sample_id"] != sample["source_sample_id"]:
            raise GSE256398AuditError("SOFT-to-H5 donor join differs")
        joined.append({**sample, **h5})
    feature_hashes = {str(row["feature_ids_sha256"]) for row in joined}
    feature_name_hashes = {str(row["feature_names_sha256"]) for row in joined}
    if len(feature_hashes) != 1 or len(feature_name_hashes) != 1:
        raise GSE256398AuditError("feature roster differs across donors")
    args.output.mkdir(parents=True)
    write_tsv(args.output / "donors.tsv", joined)
    receipt = {
        "schema_version": "masld-bench-gse256398-h5-audit-v1",
        "status": "pass_source_topology_and_h5_structure",
        "dataset_id": "gse256398",
        "human_donors": len(joined),
        "masld_relevant_donors": sum(row["development_role"] == "masld_relevant" for row in joined),
        "alcohol_ood_donors": sum(row["development_role"] == "etiology_ood" for row in joined),
        "disease_group_counts": EXPECTED_GROUPS,
        "total_nuclei": sum(int(row["nuclei"]) for row in joined),
        "human_h5_bytes": sum(int(row["bytes"]) for row in joined),
        "feature_ids_sha256": next(iter(feature_hashes)),
        "feature_names_sha256": next(iter(feature_name_hashes)),
        "genome": "GRCh38",
        "exact_cell_ranger_reference": "unresolved",
        "biological_unit": "donor",
        "alcohol_donors_are_masld_negative_controls": False,
        "labels_are_development_visible": True,
        "sealed_or_final_outcomes_read": False,
        "model_training_active": False,
        "champion_eligible": False,
    }
    if receipt["masld_relevant_donors"] != 17 or receipt["alcohol_ood_donors"] != 9:
        raise GSE256398AuditError("development/OOD donor census differs")
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
