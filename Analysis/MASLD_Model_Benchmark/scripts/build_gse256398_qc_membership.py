#!/usr/bin/env python3
"""Freeze outcome-blind retained-nucleus membership for all GSE256398 donors."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import re
import time

import numpy as np
import scrublet

from scripts.probe_gse256398_scrublet_qc import basic_qc, read_10x_counts


EXPECTED_SOURCE_ARTIFACTS = (
    "f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20"
)
EXPECTED_CROSSWALK_ARTIFACTS = (
    "54ca7af31f6ae223defdee4a7c1d7931f7cba4362f621fe9646fd08530c4c6da"
)
EXPECTED_DONORS = 26
EXPECTED_NUCLEI = 197_942
EXPECTED_FEATURES = 36_601
EXPECTED_VERSIONS = {
    "h5py": "3.15.1",
    "numpy": "2.2.6",
    "scipy": "1.15.3",
    "scrublet": "0.2.3",
}
FILE_PATTERN = re.compile(
    r"^(GSM809\d{4})_(S\d+)_CB_raw_feature_bc_matrix_filtered\.h5$"
)
QC_FIELDS = (
    "row_id",
    "gsm",
    "source_sample_id",
    "barcode",
    "total_counts",
    "detected_genes",
    "percent_mt",
    "scrublet_score",
    "scrublet_predicted_doublet",
    "basic_qc_retained",
    "final_retained",
)
MEMBERSHIP_FIELDS = (
    "row_id",
    "gsm",
    "source_sample_id",
    "barcode",
)


class GSE256398MembershipError(RuntimeError):
    """Raised when the full-cohort QC requirement differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def source_identity(path: Path) -> tuple[str, str]:
    match = FILE_PATTERN.fullmatch(path.name)
    if match is None:
        raise GSE256398MembershipError("GSE256398 H5 filename differs")
    return match.group(1), match.group(2)


def package_versions() -> dict[str, str]:
    observed = {name: importlib.metadata.version(name) for name in EXPECTED_VERSIONS}
    if observed != EXPECTED_VERSIONS:
        raise GSE256398MembershipError(
            f"QC runtime versions differ: expected {EXPECTED_VERSIONS}, observed {observed}"
        )
    return observed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if args.seed != 0:
        raise GSE256398MembershipError("Scrublet seed must equal its frozen default 0")
    if args.output.exists():
        raise GSE256398MembershipError("membership output exists")
    if sha256_file(args.source / "ARTIFACTS.json") != EXPECTED_SOURCE_ARTIFACTS:
        raise GSE256398MembershipError("source ARTIFACTS SHA-256 differs")
    if sha256_file(args.crosswalk / "ARTIFACTS.json") != EXPECTED_CROSSWALK_ARTIFACTS:
        raise GSE256398MembershipError("crosswalk ARTIFACTS SHA-256 differs")
    source_receipt = json.loads(
        (args.source / "audit/receipt.json").read_text(encoding="utf-8")
    )
    crosswalk_receipt = json.loads(
        (args.crosswalk / "summary.json").read_text(encoding="utf-8")
    )
    paths = sorted((args.source / "h5").glob("*.h5"))
    if (
        len(paths) != EXPECTED_DONORS
        or source_receipt["total_nuclei"] != EXPECTED_NUCLEI
        or crosswalk_receipt["source_features"] != EXPECTED_FEATURES
        or crosswalk_receipt["mapping_states"]
        != {
            "absent_or_retired_from_gencode_v49": 1_146,
            "stable_id_exact_gencode_v49": 35_455,
        }
    ):
        raise GSE256398MembershipError("source or crosswalk census differs")
    versions = package_versions()
    args.output.mkdir(parents=True, exist_ok=False)
    qc_path = args.output / "nuclei_qc.tsv"
    membership_path = args.output / "retained_membership.tsv"
    started = time.monotonic()
    total_input = 0
    total_basic = 0
    total_doublets = 0
    total_final = 0
    seen_rows: set[str] = set()
    donor_receipts: list[dict[str, object]] = []
    with qc_path.open("x", encoding="utf-8", newline="") as qc_handle, membership_path.open(
        "x", encoding="utf-8", newline=""
    ) as membership_handle:
        qc_writer = csv.DictWriter(
            qc_handle, fieldnames=QC_FIELDS, delimiter="\t", lineterminator="\n"
        )
        membership_writer = csv.DictWriter(
            membership_handle,
            fieldnames=MEMBERSHIP_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        qc_writer.writeheader()
        membership_writer.writeheader()
        for path in paths:
            donor_started = time.monotonic()
            gsm, sample_id = source_identity(path)
            counts, barcodes, names = read_10x_counts(path)
            metrics = basic_qc(counts, names)
            model = scrublet.Scrublet(counts, random_state=args.seed)
            scores, predicted = model.scrub_doublets()
            scores = np.asarray(scores)
            predicted = np.asarray(predicted)
            if (
                scores.shape != (counts.shape[0],)
                or predicted.shape != (counts.shape[0],)
                or not np.all(np.isfinite(scores))
                or predicted.dtype != np.bool_
                or model.threshold_ is None
            ):
                raise GSE256398MembershipError(f"Scrublet output differs for {sample_id}")
            final = metrics["basic_qc_retained"] & ~predicted
            for index, barcode in enumerate(barcodes):
                row_id = f"gse256398:{sample_id}:{barcode}"
                if row_id in seen_rows:
                    raise GSE256398MembershipError("nucleus row ID is duplicated")
                seen_rows.add(row_id)
                identity = {
                    "row_id": row_id,
                    "gsm": gsm,
                    "source_sample_id": sample_id,
                    "barcode": barcode,
                }
                qc_writer.writerow(
                    {
                        **identity,
                        "total_counts": int(metrics["total_counts"][index]),
                        "detected_genes": int(metrics["detected_genes"][index]),
                        "percent_mt": format(float(metrics["percent_mt"][index]), ".17g"),
                        "scrublet_score": format(float(scores[index]), ".17g"),
                        "scrublet_predicted_doublet": str(bool(predicted[index])).lower(),
                        "basic_qc_retained": str(
                            bool(metrics["basic_qc_retained"][index])
                        ).lower(),
                        "final_retained": str(bool(final[index])).lower(),
                    }
                )
                if final[index]:
                    membership_writer.writerow(identity)
            donor = {
                "gsm": gsm,
                "source_sample_id": sample_id,
                "nuclei_input": counts.shape[0],
                "basic_qc_retained": int(metrics["basic_qc_retained"].sum()),
                "scrublet_predicted_doublets": int(predicted.sum()),
                "final_retained": int(final.sum()),
                "scrublet_threshold": float(model.threshold_),
                "elapsed_seconds": time.monotonic() - donor_started,
            }
            if donor["final_retained"] < 1:
                raise GSE256398MembershipError(f"QC retained no nuclei for {sample_id}")
            donor_receipts.append(donor)
            total_input += int(donor["nuclei_input"])
            total_basic += int(donor["basic_qc_retained"])
            total_doublets += int(donor["scrublet_predicted_doublets"])
            total_final += int(donor["final_retained"])
    if total_input != EXPECTED_NUCLEI or len(seen_rows) != EXPECTED_NUCLEI:
        raise GSE256398MembershipError("full-cohort nucleus census differs")
    receipt = {
        "schema_version": "masld-bench-gse256398-qc-membership-v1",
        "status": "pass",
        "dataset_id": "gse256398",
        "donors": len(donor_receipts),
        "nuclei_input": total_input,
        "basic_qc_retained": total_basic,
        "scrublet_predicted_doublets": total_doublets,
        "final_retained": total_final,
        "features_native": EXPECTED_FEATURES,
        "features_project_allowed": 35_455,
        "features_project_masked": 1_146,
        "seed": args.seed,
        "package_versions": versions,
        "qc_order": "Scrublet_per_donor_then_nFeature_200_to_6499_nCount_below_40000_percent_mt_below_20",
        "source_artifacts_sha256": EXPECTED_SOURCE_ARTIFACTS,
        "crosswalk_artifacts_sha256": EXPECTED_CROSSWALK_ARTIFACTS,
        "nucleus_row_id": "dataset_id:source_sample_id:sample_local_barcode",
        "barcode_level_cell_labels": "structurally_missing",
        "disease_age_sex_or_histology_read": False,
        "normalization_or_variable_feature_selection_fitted": False,
        "model_training_activated": False,
        "donor_receipts": donor_receipts,
        "elapsed_seconds": time.monotonic() - started,
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in receipt.items() if key != "donor_receipts"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
