#!/usr/bin/env python3
"""Independently recount a seed-fixed donor/peak sample from fragment files."""

from __future__ import annotations

import argparse
import csv
import gzip
import random
from pathlib import Path

import h5py
import pysam
import scipy.io


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
PROJECT = Path(__file__).resolve().parents[4]
EXPECTED = (PROJECT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID).resolve()
LABEL_MAP = {
    "GSE244832": {
        "Hepatocytes": "hepatocyte", "Fibroblasts": "stellate",
        "Macrophages": "macrophage", "Cholangiocytes": "cholangiocyte",
        "T_cells": "t_nk", "Resident_NK": "t_nk", "Circulating_NK_NKT": "t_nk",
    },
    "GSE281367": {
        "Hepatocyte": "hepatocyte", "Stellate_Cell": "stellate",
        "Macrophage": "macrophage", "Kupffer_Cell": "macrophage",
        "Cholangiocyte": "cholangiocyte", "NK_T_Cell": "t_nk",
    },
}
H5AD = {
    "GSE244832": PROJECT / "Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_relabeled.h5ad",
    "GSE281367": PROJECT / "Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=EXPECTED)
    return parser.parse_args()


def decode(values) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def obs_column(group: h5py.Group, name: str) -> list[str]:
    value = group[name]
    if isinstance(value, h5py.Dataset):
        return decode(value[:])
    categories = decode(value["categories"][:])
    return [categories[int(code)] if int(code) >= 0 else "" for code in value["codes"][:]]


def cell_sets(cohort: str) -> dict[tuple[str, str], set[str]]:
    result: dict[tuple[str, str], set[str]] = {}
    with h5py.File(H5AD[cohort], "r") as handle:
        obs = handle["obs"]
        barcodes = decode(obs["_index"][:])
        donors = obs_column(obs, "donor_id")
        labels = obs_column(obs, "cell_type")
    for barcode, donor, label in zip(barcodes, donors, labels):
        lineage = LABEL_MAP[cohort].get(label)
        if lineage:
            result.setdefault((donor, lineage), set()).add(barcode)
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fragment_path(cohort: str, donor: str) -> Path:
    if cohort == "GSE244832":
        return PROJECT / f"Analysis/ATAC/Human_Multiome/results/fragments/{donor}_fragments.tsv.gz"
    return PROJECT / f"Analysis/ATAC/Human_External/cellranger/{donor}/outs/fragments.tsv.gz"


def direct_counts(
    path: Path, chrom: str, start0: int, end: int, barcodes: set[str]
) -> tuple[int, int]:
    records = 0
    multiplicity = 0
    with pysam.TabixFile(str(path)) as tabix:
        for record in tabix.fetch(chrom, start0, end):
            fields = record.split("\t")
            if int(fields[1]) < end and int(fields[2]) > start0 and fields[3] in barcodes:
                records += 1
                multiplicity += int(fields[4]) if len(fields) > 4 else 1
    return records, multiplicity


def main() -> None:
    args = arguments()
    root = args.candidate_root.resolve()
    if root != EXPECTED:
        raise RuntimeError(f"unsafe candidate root: {root}")
    out = root / "recount"
    if out.exists():
        raise RuntimeError(f"refusing to overwrite exact recount: {out}")
    out.mkdir(parents=True)
    rng = random.Random(42)
    audit = []
    for cohort in ("GSE244832", "GSE281367"):
        cells = cell_sets(cohort)
        for lineage in ("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"):
            source = root / "counts" / cohort
            donors = read_tsv(source / f"{lineage}.donors.tsv")
            peaks = read_tsv(source / f"{lineage}.peaks.tsv")
            with gzip.open(source / f"{lineage}.mtx.gz", "rb") as handle:
                matrix = scipy.io.mmread(handle).tocsr()
            eligible = [index for index, row in enumerate(donors) if row["contrast_eligible"] == "TRUE"]
            donor_index = rng.choice(eligible)
            nonzero = matrix.getrow(donor_index).indices.tolist()
            if not nonzero:
                raise RuntimeError(f"no nonzero peaks for recount: {cohort}/{lineage}")
            peak_index = rng.choice(nonzero)
            coordinate = peaks[peak_index]["peak_coordinate"]
            chrom, interval = coordinate.split(":")
            start0, end = map(int, interval.split("-"))
            donor = donors[donor_index]["donor_id"]
            expected = int(matrix[donor_index, peak_index])
            observed_records, observed_multiplicity = direct_counts(
                fragment_path(cohort, donor), chrom, start0, end, cells[(donor, lineage)]
            )
            audit.append({
                "release_id": RELEASE_ID,
                "cohort": cohort,
                "lineage": lineage,
                "donor_id": donor,
                "peak_coordinate": coordinate,
                "counting_unit": "deduplicated_fragment_record",
                "matrix_count": expected,
                "direct_fragment_records": observed_records,
                "direct_fragment_multiplicity": observed_multiplicity,
                "record_count_match": str(expected == observed_records).upper(),
                "multiplicity_count_match": str(expected == observed_multiplicity).upper(),
            })
    path = out / "exact_fragment_recount.tsv"
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(audit[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit)
    if any(row["record_count_match"] != "TRUE" for row in audit):
        raise RuntimeError(
            "at least one direct deduplicated-fragment-record recount disagrees with the matrix; "
            f"diagnostics preserved at {path}"
        )
    print(f"Exact fragment recount passed for {len(audit)} donor-peak combinations")


if __name__ == "__main__":
    main()
