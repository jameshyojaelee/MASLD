#!/usr/bin/env python3
"""Re-audit GSE106737/GSE83452 aliases without exporting expression values."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


FAMILY_ID = "antwerp_inserm_shared"
EXPECTED_MATRIX_SHA256 = {
    "GSE106737": "d6c54bdb441f53f9fccc6d9592413a5c51a5cdce34900bb0b37ea9099e52cb7d",
    "GSE83452": "f704bf623e4714e0c38317cc9bdf37e064ddae41185d1485b079e6fdc224603d",
}
N_PROBES = 5_000
CORRELATION_THRESHOLD = 0.9999


class CrossAccessionOverlapError(RuntimeError):
    """Raised when the frozen overlap evidence differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_axis(path: Path) -> tuple[list[str], set[str]]:
    samples: list[str] = []
    probes: set[str] = set()
    active = False
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                active = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not active:
                continue
            values = next(csv.reader([line], delimiter="\t"))
            if not samples:
                if values[0] != "ID_REF":
                    raise CrossAccessionOverlapError("series-matrix feature header differs")
                samples = [value.strip('"') for value in values[1:]]
                if len(samples) != len(set(samples)):
                    raise CrossAccessionOverlapError("series-matrix sample IDs are duplicated")
                continue
            probe = values[0].strip('"')
            if not probe or probe in probes:
                raise CrossAccessionOverlapError("series-matrix probe IDs differ")
            probes.add(probe)
    if not samples or not probes:
        raise CrossAccessionOverlapError("series-matrix axis is absent")
    return samples, probes


def _selected_matrix(
    path: Path, expected_samples: list[str], selected: tuple[str, ...]
) -> np.ndarray:
    selected_set = set(selected)
    rows: dict[str, np.ndarray] = {}
    active = False
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                active = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not active:
                continue
            values = next(csv.reader([line], delimiter="\t"))
            if values[0] == "ID_REF":
                observed = [value.strip('"') for value in values[1:]]
                if observed != expected_samples:
                    raise CrossAccessionOverlapError("series-matrix sample order differs")
                continue
            probe = values[0].strip('"')
            if probe not in selected_set:
                continue
            if probe in rows or len(values) != len(expected_samples) + 1:
                raise CrossAccessionOverlapError("selected probe matrix differs")
            try:
                row = np.asarray([float(value) for value in values[1:]], dtype=np.float64)
            except ValueError as error:
                raise CrossAccessionOverlapError("selected expression value is nonnumeric") from error
            if not np.isfinite(row).all():
                raise CrossAccessionOverlapError("selected expression value is nonfinite")
            rows[probe] = row
    if set(rows) != selected_set:
        raise CrossAccessionOverlapError("selected probe rows are incomplete")
    return np.vstack([rows[probe] for probe in selected])


def _correlations(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left_centered = left - left.mean(axis=0, keepdims=True)
    right_centered = right - right.mean(axis=0, keepdims=True)
    left_scale = np.sqrt(np.sum(left_centered * left_centered, axis=0))
    right_scale = np.sqrt(np.sum(right_centered * right_centered, axis=0))
    if np.any(left_scale == 0) or np.any(right_scale == 0):
        raise CrossAccessionOverlapError("constant expression fingerprint found")
    result = (left_centered.T @ right_centered) / np.outer(left_scale, right_scale)
    if not np.isfinite(result).all():
        raise CrossAccessionOverlapError("expression correlation is nonfinite")
    return result


def _read_gse106737_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset_id",
            "sample_id",
            "source_accession",
            "biological_unit_id",
            "timepoint",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise CrossAccessionOverlapError("GSE106737 manifest fields differ")
        rows = [
            row
            for row in reader
            if row["dataset_id"] == "GSE106737"
        ]
    result = {row["sample_id"]: row for row in rows}
    if len(rows) != 111 or len(result) != 111:
        raise CrossAccessionOverlapError("GSE106737 sample census differs")
    if len({row["biological_unit_id"] for row in rows}) != 70:
        raise CrossAccessionOverlapError("GSE106737 participant census differs")
    return result


def _read_gse83452_records(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_accession", "participant_id", "timepoint"}
        if not required.issubset(reader.fieldnames or ()):
            raise CrossAccessionOverlapError("GSE83452 participant fields differ")
        rows = [dict(row) for row in reader]
    result = {row["sample_accession"]: row for row in rows}
    if len(rows) != 231 or len(result) != 231:
        raise CrossAccessionOverlapError("GSE83452 record census differs")
    if len({row["participant_id"] for row in rows}) != 171:
        raise CrossAccessionOverlapError("GSE83452 participant census differs")
    return result


def _evidence_digest(
    left_sha: str,
    right_sha: str,
    probe_sha: str,
    left_sample: str,
    right_sample: str,
    correlation: float,
) -> str:
    payload = "\x1f".join(
        (
            "masld-bench-cross-accession-expression-evidence-v1",
            left_sha,
            right_sha,
            probe_sha,
            left_sample,
            right_sample,
            format(correlation, ".17g"),
        )
    )
    return sha256(payload.encode("utf-8")).hexdigest()


def _read_timepoint_dispositions(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "gse106737_sample_accession",
            "gse83452_sample_accession",
            "gse106737_source_timepoint",
            "gse83452_source_timepoint",
            "timepoint_token_disposition",
            "bad_match_indicated_by_timepoint",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise CrossAccessionOverlapError("timepoint diagnostic fields differ")
        rows = [dict(row) for row in reader]
    result = {
        (row["gse106737_sample_accession"], row["gse83452_sample_accession"]): row
        for row in rows
    }
    allowed = {"exact_same_source_token", "encoding_semantics_punctuation_only"}
    if (
        len(rows) != 78
        or len(result) != 78
        or any(row["timepoint_token_disposition"] not in allowed for row in rows)
        or any(row["bad_match_indicated_by_timepoint"].lower() != "false" for row in rows)
    ):
        raise CrossAccessionOverlapError("timepoint diagnostic disposition differs")
    return result


def audit_overlap(
    gse106737_matrix: Path,
    gse83452_matrix: Path,
    gse106737_manifest: Path,
    gse83452_records: Path,
    timepoint_diagnostic: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    matrix_paths = {
        "GSE106737": gse106737_matrix,
        "GSE83452": gse83452_matrix,
    }
    matrix_hashes = {key: sha256_file(value) for key, value in matrix_paths.items()}
    if matrix_hashes != EXPECTED_MATRIX_SHA256:
        raise CrossAccessionOverlapError("frozen series-matrix SHA-256 differs")

    left_samples, left_probes = _matrix_axis(gse106737_matrix)
    right_samples, right_probes = _matrix_axis(gse83452_matrix)
    if len(left_samples) != 111 or len(right_samples) != 231:
        raise CrossAccessionOverlapError("series-matrix sample count differs")
    common = sorted(left_probes & right_probes)
    if len(common) < N_PROBES:
        raise CrossAccessionOverlapError("insufficient common platform probes")
    selected = tuple(common[:N_PROBES])
    probe_sha = sha256("\n".join(selected).encode("utf-8")).hexdigest()
    left = _selected_matrix(gse106737_matrix, left_samples, selected)
    right = _selected_matrix(gse83452_matrix, right_samples, selected)
    correlations = _correlations(left, right)
    left_best = np.argmax(correlations, axis=1)
    right_best = np.argmax(correlations, axis=0)

    left_manifest = _read_gse106737_manifest(gse106737_manifest)
    right_manifest = _read_gse83452_records(gse83452_records)
    timepoint_dispositions = _read_timepoint_dispositions(
        timepoint_diagnostic / "cross_accession_timepoint_mismatches.tsv"
    )
    if set(left_samples) != set(left_manifest) or set(right_samples) != set(right_manifest):
        raise CrossAccessionOverlapError("matrix and participant axes differ")

    matches: list[dict[str, object]] = []
    participant_pairs: dict[tuple[str, str], list[dict[str, object]]] = {}
    for left_index, right_index in enumerate(left_best):
        correlation = float(correlations[left_index, right_index])
        reciprocal = int(right_best[right_index]) == left_index
        if not reciprocal or correlation < CORRELATION_THRESHOLD:
            continue
        left_sample = left_samples[left_index]
        right_sample = right_samples[int(right_index)]
        left_record = left_manifest[left_sample]
        right_record = right_manifest[right_sample]
        disposition_key = (left_sample, right_sample)
        if disposition_key not in timepoint_dispositions:
            raise CrossAccessionOverlapError("matched alias lacks timepoint disposition")
        timepoint = timepoint_dispositions[disposition_key]
        if (
            timepoint["gse106737_source_timepoint"] != left_record["timepoint"]
            or timepoint["gse83452_source_timepoint"] != right_record["timepoint"]
        ):
            raise CrossAccessionOverlapError("source timepoint differs from frozen diagnostic")
        row: dict[str, object] = {
            "cohort_family_id": FAMILY_ID,
            "gse106737_sample_accession": left_sample,
            "gse83452_sample_accession": right_sample,
            "gse106737_participant_alias": f"gse106737::{left_record['biological_unit_id']}",
            "gse83452_participant_alias": right_record["participant_id"],
            "gse106737_timepoint": left_record["timepoint"],
            "gse83452_timepoint": right_record["timepoint"],
            "timepoint_token_disposition": timepoint["timepoint_token_disposition"],
            "pearson_r": format(correlation, ".15g"),
            "reciprocal_best": True,
            "duplicate_threshold_pass": True,
            "n_fixed_probes": N_PROBES,
            "probe_selection": "first_5000_lexicographic_common_GPL16686_feature_IDs",
            "nonreversible_expression_evidence_sha256": _evidence_digest(
                matrix_hashes["GSE106737"],
                matrix_hashes["GSE83452"],
                probe_sha,
                left_sample,
                right_sample,
                correlation,
            ),
        }
        matches.append(row)
        participant_pairs.setdefault(
            (
                str(row["gse106737_participant_alias"]),
                str(row["gse83452_participant_alias"]),
            ),
            [],
        ).append(row)

    if len(matches) != 78 or len(participant_pairs) != 41:
        raise CrossAccessionOverlapError("cross-accession overlap count differs")
    if len({row["gse83452_sample_accession"] for row in matches}) != 78:
        raise CrossAccessionOverlapError("GSE83452 duplicate match found")
    if any(len({pair[1] for pair in participant_pairs if pair[0] == left}) != 1 for left, _ in participant_pairs):
        raise CrossAccessionOverlapError("GSE106737 participant alias is one-to-many")
    if any(len({pair[0] for pair in participant_pairs if pair[1] == right}) != 1 for _, right in participant_pairs):
        raise CrossAccessionOverlapError("GSE83452 participant alias is one-to-many")
    observed_pair_keys = {
        (str(row["gse106737_sample_accession"]), str(row["gse83452_sample_accession"]))
        for row in matches
    }
    if observed_pair_keys != set(timepoint_dispositions):
        raise CrossAccessionOverlapError("fingerprint and timepoint diagnostic pairs differ")

    aliases: list[dict[str, object]] = []
    for (left_participant, right_participant), rows in sorted(participant_pairs.items()):
        pair_payload = "\x1f".join(
            (
                "masld-bench-cross-accession-participant-alias-v1",
                left_participant,
                right_participant,
                *sorted(str(row["nonreversible_expression_evidence_sha256"]) for row in rows),
            )
        )
        aliases.append(
            {
                "cohort_family_id": FAMILY_ID,
                "family_participant_group_id": (
                    f"{FAMILY_ID}::{sha256(pair_payload.encode('utf-8')).hexdigest()[:24]}"
                ),
                "gse106737_participant_alias": left_participant,
                "gse83452_participant_alias": right_participant,
                "matched_array_aliases": len(rows),
                "matched_timepoints": ";".join(
                    sorted(
                        {
                            f"{row['gse106737_timepoint']}|{row['gse83452_timepoint']}"
                            for row in rows
                        }
                    )
                ),
                "identity_state": "fingerprint_confirmed_cross_accession_alias",
            }
        )

    summary = {
        "schema_version": "masld-bench-cross-accession-overlap-audit-v1",
        "status": "pass_confirmed_cross_accession_family_overlap",
        "cohort_family_id": FAMILY_ID,
        "accession_aliases": ["GSE106737", "GSE83452"],
        "series_matrix_sha256": matrix_hashes,
        "series_matrix_role": "contamination_and_alias_audit_only_prohibited_model_input",
        "common_feature_count": len(common),
        "fixed_probe_count": N_PROBES,
        "fixed_probe_ids_sha256": probe_sha,
        "correlation_threshold": CORRELATION_THRESHOLD,
        "reciprocal_array_aliases": len(matches),
        "cross_accession_participant_aliases": len(aliases),
        "gse106737_arrays": len(left_samples),
        "gse106737_participants": len({row["biological_unit_id"] for row in left_manifest.values()}),
        "gse83452_arrays": len(right_samples),
        "gse83452_participants": len({row["participant_id"] for row in right_manifest.values()}),
        "exact_family_union_participant_count": "join_unresolved_for_unmatched_GSE106737_participants",
        "family_union_participant_range": [171, 200],
        "independent_cohort_claim_allowed_between_accessions": False,
        "same_outer_family_required": True,
        "project_exposure": "project_exposed_external_development",
        "source_timepoints_retained_verbatim": True,
        "timepoint_disposition_registry_required": True,
        "timepoints_coerced_or_overwritten": False,
        "labels_read_for_overlap_audit": False,
        "expression_values_exported": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    return matches, aliases, summary


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gse106737-series-matrix", type=Path, required=True)
    parser.add_argument("--gse83452-series-matrix", type=Path, required=True)
    parser.add_argument("--gse106737-manifest", type=Path, required=True)
    parser.add_argument("--gse83452-records", type=Path, required=True)
    parser.add_argument("--timepoint-diagnostic", type=Path, required=True)
    parser.add_argument("--timepoint-diagnostic-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=False)
    if (
        sha256_file(arguments.timepoint_diagnostic / "ARTIFACTS.json")
        != arguments.timepoint_diagnostic_artifacts_sha256
    ):
        raise CrossAccessionOverlapError("timepoint diagnostic ARTIFACTS SHA-256 differs")
    matches, aliases, summary = audit_overlap(
        arguments.gse106737_series_matrix,
        arguments.gse83452_series_matrix,
        arguments.gse106737_manifest,
        arguments.gse83452_records,
        arguments.timepoint_diagnostic,
    )
    _write_tsv(arguments.output / "cross_accession_sample_aliases.tsv", matches)
    _write_tsv(arguments.output / "cross_accession_participant_aliases.tsv", aliases)
    (arguments.output / "cross_accession_overlap_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
