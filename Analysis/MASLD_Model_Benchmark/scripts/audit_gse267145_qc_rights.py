#!/usr/bin/env python3
"""Freeze assay-native participant QC and conservative reuse terms for GSE267145."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import shlex
from statistics import median


class QCRightsError(RuntimeError):
    """Raised when participant, measurement, QC, or public-use evidence differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_join(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 99 or len({row["participant_id"] for row in rows}) != 99:
        raise QCRightsError("participant join differs")
    return rows


def matrix_qc(
    path: Path,
    participants: set[str],
    *,
    h3: bool,
) -> dict[str, dict[str, float | int]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        if h3:
            header = shlex.split(handle.readline())
            sample_ids = [value.split("_", 1)[0] for value in header]
            row_offset = 1
        else:
            header = handle.readline().rstrip("\n").split("\t")
            sample_ids = header[1:]
            row_offset = 1
        index_by_participant = {
            participant: sample_ids.index(participant) + row_offset
            for participant in participants
            if participant in sample_ids
        }
        if set(index_by_participant) != participants:
            raise QCRightsError("participant matrix axis differs")
        totals = {participant: 0.0 for participant in participants}
        detected = {participant: 0 for participant in participants}
        feature_rows = 0
        for line_number, line in enumerate(handle, start=2):
            values = shlex.split(line) if h3 else line.rstrip("\n").split("\t")
            # The RNA header includes its feature-ID field; the H3 header does not.
            expected = len(header) + (1 if h3 else 0)
            if len(values) != expected:
                raise QCRightsError(f"matrix width differs at {line_number}")
            feature_rows += 1
            for participant, index in index_by_participant.items():
                value = float(values[index])
                if not math.isfinite(value) or value < 0:
                    raise QCRightsError("matrix value is not finite nonnegative")
                totals[participant] += value
                detected[participant] += int(value > 0)
    expected_rows = 96_460 if h3 else 43_285
    if feature_rows != expected_rows:
        raise QCRightsError("matrix feature count differs")
    return {
        participant: {
            "library_sum": totals[participant],
            "detected_features": detected[participant],
            "feature_rows": feature_rows,
        }
        for participant in participants
    }


def robust_z(values: dict[str, float]) -> dict[str, float]:
    center = median(values.values())
    mad = median(abs(value - center) for value in values.values())
    if mad == 0:
        return {key: 0.0 for key in values}
    return {key: 0.67448975 * (value - center) / mad for key, value in values.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--join", type=Path, required=True)
    parser.add_argument("--join-artifacts-sha256", required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--measurement", type=Path, required=True)
    parser.add_argument("--measurement-artifacts-sha256", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reference-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for root, expected in (
        (args.join, args.join_artifacts_sha256),
        (args.matrix, args.matrix_artifacts_sha256),
        (args.measurement, args.measurement_artifacts_sha256),
        (args.reference, args.reference_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise QCRightsError("input ARTIFACTS SHA-256 differs")
    rows = read_join(args.join / "participant_join.tsv")
    participants = {row["participant_id"] for row in rows}
    rna = matrix_qc(
        args.matrix / "raw" / "GSE269412_RNA.txt.gz", participants, h3=False
    )
    h3 = matrix_qc(
        args.matrix / "raw" / "GSE267119_H3K27ac.txt.gz", participants, h3=True
    )
    rna_z = robust_z({key: float(value["library_sum"]) for key, value in rna.items()})
    h3_z = robust_z({key: float(value["library_sum"]) for key, value in h3.items()})
    output_rows: list[dict[str, object]] = []
    for row in rows:
        participant = row["participant_id"]
        hard_failures: list[str] = []
        if float(rna[participant]["library_sum"]) <= 0:
            hard_failures.append("rna_zero_library")
        if float(h3[participant]["library_sum"]) <= 0:
            hard_failures.append("h3k27ac_zero_library")
        output_rows.append(
            {
                "participant_id": participant,
                "stage": row["stage"],
                "sex": row["sex"],
                "outer_fold": row["outer_fold"],
                "rna_library_sum": f"{float(rna[participant]['library_sum']):.12g}",
                "rna_detected_genes": rna[participant]["detected_features"],
                "rna_library_robust_z": f"{rna_z[participant]:.8g}",
                "h3k27ac_library_sum": f"{float(h3[participant]['library_sum']):.12g}",
                "h3k27ac_detected_regions": h3[participant]["detected_features"],
                "h3k27ac_library_robust_z": f"{h3_z[participant]:.8g}",
                "hard_qc_state": "below_qc" if hard_failures else "observed",
                "hard_qc_reasons": ";".join(hard_failures),
                "distributional_outlier_is_automatic_exclusion": "false",
            }
        )
    args.output.mkdir(parents=True, exist_ok=False)
    fields = tuple(output_rows[0])
    with (args.output / "participant_qc.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    hard_failures = sum(row["hard_qc_state"] == "below_qc" for row in output_rows)
    qc = {
        "schema_version": "masld-bench-gse267145-qc-v1",
        "status": "pass" if hard_failures == 0 else "pass_with_below_qc_participants",
        "biological_unit": "participant",
        "participants": len(rows),
        "hard_qc_failures": hard_failures,
        "stage_counts": dict(sorted(Counter(row["stage"] for row in rows).items())),
        "sex_counts": dict(sorted(Counter(row["sex"] for row in rows).items())),
        "hard_qc_rule": "finite_nonnegative_and_positive_assay_library_sum",
        "distributional_outliers_reported_not_excluded": True,
        "cells_or_libraries_as_replicates": False,
        "normalization_fit_outside_outer_training_fold": False,
    }
    rights = {
        "schema_version": "masld-bench-gse267145-rights-v1",
        "access_tier": "public",
        "new_dua_or_controlled_access_required": False,
        "internal_nonclinical_research_use": "allowed_by_project_public_data_definition",
        "source_terms_apply": True,
        "explicit_data_license_detected": False,
        "raw_or_processed_data_redistribution": "prohibited_by_project_policy",
        "released_weight_redistribution": "requires_model_specific_terms_and_legal_review",
        "clinical_use": False,
        "open_champion_eligibility": "blocked_until_derivative_weight_review_if_this_source_materially_trains_the_champion",
        "source_evidence": [
            "official_GEO_SuperSeries_and_subseries_are_openly_downloadable",
            "source_article_data_availability_points_to_GEO",
            "no_new_DUA_or_controlled_access_workflow_encountered",
        ],
    }
    (args.output / "qc_summary.json").write_text(
        json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.output / "rights_contract.json").write_text(
        json.dumps(rights, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"qc": qc, "rights": rights}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
