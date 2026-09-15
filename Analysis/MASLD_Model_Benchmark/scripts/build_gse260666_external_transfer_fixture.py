#!/usr/bin/env python3
"""Build outcome-separated GSE260666 RNA transfer and contamination fixtures."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import hmac
import json
from pathlib import Path
import re
import subprocess
from typing import Any, Mapping, Sequence

import numpy as np


CLASSES = ("NOR", "NAFL", "NASH")
SOURCE_TO_STAGE3 = {
    "healthy_control": "NOR",
    "non-alcoholic_fatty_liver_disease_(nafld)": "NAFL",
    "non-alcoholic_steatohepatitis_(nash)": "NASH",
}
EXPECTED_SOURCE_COUNTS = {
    "healthy_control": 6,
    "non-alcoholic_fatty_liver_disease_(nafld)": 6,
    "non-alcoholic_steatohepatitis_(nash)": 4,
}


class ExternalTransferFixtureError(RuntimeError):
    """Raised when the external transfer fixture would not meet its requirements."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExternalTransferFixtureError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _row_id(sample_accession: str) -> str:
    return "g260_" + hashlib.sha256(
        f"gse260666_external_v1\0{sample_accession}".encode("utf-8")
    ).hexdigest()[:20]


def _salted_signature(key: bytes, label: str, values: np.ndarray) -> str:
    quantized = np.round(np.asarray(values, dtype=np.float64), 6).astype("<f8")
    return hmac.new(key, label.encode("utf-8") + b"\0" + quantized.tobytes(), hashlib.sha256).hexdigest()


def _source_aliases_from_soft(path: Path) -> set[str]:
    aliases: set[str] = set()
    pattern = re.compile(r"\b(?:GSE|GSM|SAMN|SRX|PRJNA)[0-9]+\b")
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            aliases.update(pattern.findall(line))
    if not aliases:
        raise ExternalTransferFixtureError("training SOFT contains no source aliases")
    return aliases


def _prior_project_use(
    project_root: Path,
    benchmark_root: Path,
    aliases: Sequence[str],
    fingerprint_key: bytes,
) -> dict[str, Any]:
    try:
        excluded = benchmark_root.relative_to(project_root).as_posix()
    except ValueError as error:
        raise ExternalTransferFixtureError("benchmark root is outside project root") from error
    command = [
        "rg",
        "-l",
        "-F",
        "--hidden",
        "--glob",
        f"!{excluded}/**",
        "--glob",
        "!.git/**",
    ]
    for alias in sorted(set(aliases)):
        command.extend(("-e", alias))
    command.append(".")
    completed = subprocess.run(
        command,
        cwd=project_root,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode not in {0, 1}:
        raise ExternalTransferFixtureError(
            f"prior project-use search failed with rg exit {completed.returncode}"
        )
    paths = sorted({line for line in completed.stdout.splitlines() if line})
    hashed = [
        hmac.new(fingerprint_key, path.encode("utf-8"), hashlib.sha256).hexdigest()
        for path in paths
    ]
    return {
        "search_scope": "project_root_excluding_benchmark_package_and_git_metadata",
        "aliases_searched": len(set(aliases)),
        "matched_file_count": len(paths),
        "matched_file_hmac_sha256": hashed,
        "raw_paths_retained": False,
    }


def _log_cpm(matrix: np.ndarray) -> np.ndarray:
    value = np.asarray(matrix, dtype=np.float64)
    if value.ndim != 2 or not np.all(np.isfinite(value)) or np.any(value < 0):
        raise ExternalTransferFixtureError("RNA matrix is not finite and nonnegative")
    totals = value.sum(axis=1)
    if np.any(totals <= 0):
        raise ExternalTransferFixtureError("RNA participant library is empty")
    return np.log1p(value * (1_000_000.0 / totals[:, None]))


def _ranks(matrix: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return np.vstack([rankdata(row, method="average") for row in matrix])


def _row_standardize(matrix: np.ndarray) -> np.ndarray:
    centered = matrix - matrix.mean(axis=1, keepdims=True)
    scales = np.sqrt(np.sum(centered * centered, axis=1, keepdims=True))
    if np.any(scales <= np.finfo(np.float64).eps):
        raise ExternalTransferFixtureError("expression fingerprint is constant")
    return centered / scales


def load_external_counts(
    matrix_path: Path,
    crosswalk_path: Path,
    participant_rows: Sequence[Mapping[str, str]],
) -> tuple[np.ndarray, list[str], list[dict[str, Any]], dict[str, int]]:
    crosswalk_fields, crosswalk = read_tsv(crosswalk_path)
    required_crosswalk = {
        "source_feature_id",
        "mapping_state",
        "gencode_v49_stable_id",
    }
    if not required_crosswalk <= set(crosswalk_fields):
        raise ExternalTransferFixtureError("GSE260666 gene crosswalk schema differs")
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_axis = header[1:]
        expected_axis = [row["participant_id"] for row in participant_rows]
        if set(sample_axis) != set(expected_axis) or len(sample_axis) != len(expected_axis):
            raise ExternalTransferFixtureError("raw-count sample axis differs from participant join")
        axis_lookup = {value: index for index, value in enumerate(sample_axis)}
        reorder = [axis_lookup[value] for value in expected_axis]
        stable_to_vectors: dict[str, list[np.ndarray]] = defaultdict(list)
        stable_to_symbols: dict[str, list[str]] = defaultdict(list)
        line_count = 0
        for line_count, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header) or line_count > len(crosswalk):
                raise ExternalTransferFixtureError("raw-count matrix width or feature axis differs")
            cross = crosswalk[line_count - 1]
            if fields[0] != cross["source_feature_id"]:
                raise ExternalTransferFixtureError("raw-count and crosswalk feature order differs")
            try:
                vector = np.asarray([int(value) for value in fields[1:]], dtype=np.int64)
            except ValueError as error:
                raise ExternalTransferFixtureError("raw count is not an integer") from error
            if np.any(vector < 0):
                raise ExternalTransferFixtureError("raw count is negative")
            if cross["mapping_state"] == "unique_symbol":
                stable = cross["gencode_v49_stable_id"]
                if not stable:
                    raise ExternalTransferFixtureError("unique symbol lacks a stable ID")
                stable_to_vectors[stable].append(vector[reorder])
                stable_to_symbols[stable].append(fields[0])
        if line_count != len(crosswalk):
            raise ExternalTransferFixtureError("raw-count and crosswalk row counts differ")
    stable_ids = sorted(stable_to_vectors)
    if not stable_ids:
        raise ExternalTransferFixtureError("no uniquely mapped GSE260666 genes")
    columns: list[np.ndarray] = []
    feature_rows: list[dict[str, Any]] = []
    for index, stable in enumerate(stable_ids):
        combined = np.sum(np.vstack(stable_to_vectors[stable]), axis=0, dtype=np.int64)
        columns.append(combined)
        symbols = sorted(stable_to_symbols[stable])
        feature_rows.append(
            {
                "feature_index": index,
                "stable_gene_id": stable,
                "source_feature_count": len(symbols),
                "source_symbols_sha256": hashlib.sha256("\n".join(symbols).encode()).hexdigest(),
            }
        )
    counts = np.column_stack(columns).astype(np.uint64, copy=False)
    mapping_counts = Counter(row["mapping_state"] for row in crosswalk)
    return counts, stable_ids, feature_rows, dict(sorted(mapping_counts.items()))


def build_fixture(
    *,
    activation: Path,
    matrix: Path,
    training_molecular: Path,
    training_join: Path,
    task_spec: Path,
    project_root: Path,
    benchmark_root: Path,
    fingerprint_salt_file: Path,
    output: Path,
    activation_artifacts_sha256: str,
    matrix_artifacts_sha256: str,
    training_molecular_artifacts_sha256: str,
    training_join_artifacts_sha256: str,
    task_spec_sha256: str,
    expected_external_n: int = 16,
    expected_training_n: int = 99,
) -> dict[str, Any]:
    if output.exists():
        raise ExternalTransferFixtureError(f"refusing to overwrite fixture: {output}")
    expected_inputs = (
        (activation / "ARTIFACTS.json", activation_artifacts_sha256),
        (matrix / "ARTIFACTS.json", matrix_artifacts_sha256),
        (training_molecular / "ARTIFACTS.json", training_molecular_artifacts_sha256),
        (training_join / "ARTIFACTS.json", training_join_artifacts_sha256),
        (task_spec, task_spec_sha256),
    )
    for path, expected in expected_inputs:
        if sha256_file(path) != expected:
            raise ExternalTransferFixtureError(f"input SHA-256 differs: {path}")
    salt_mode = fingerprint_salt_file.stat().st_mode & 0o777
    if salt_mode & 0o077:
        raise ExternalTransferFixtureError("fingerprint salt file is group/world accessible")
    fingerprint_salt = fingerprint_salt_file.read_bytes()
    if len(fingerprint_salt) < 32:
        raise ExternalTransferFixtureError("fingerprint salt is shorter than 32 bytes")
    participant_fields, participant_rows = read_tsv(activation / "participant_join.tsv")
    required_participant = {
        "participant_id",
        "sample_accession",
        "biosample_accession",
        "sra_experiment",
        "group",
        "fibrosis",
        "sex",
        "biological_unit",
    }
    if not required_participant <= set(participant_fields):
        raise ExternalTransferFixtureError("GSE260666 participant join schema differs")
    if len(participant_rows) != expected_external_n:
        raise ExternalTransferFixtureError("GSE260666 participant count differs")
    source_counts = Counter(row["group"] for row in participant_rows)
    if expected_external_n == 16 and source_counts != Counter(EXPECTED_SOURCE_COUNTS):
        raise ExternalTransferFixtureError("GSE260666 source label counts differ")
    if any(
        row["fibrosis"] != "not_reported"
        or row["sex"] != "not_reported"
        or row["biological_unit"] != "participant"
        for row in participant_rows
    ):
        raise ExternalTransferFixtureError("GSE260666 missing metadata or biological unit differs")
    for field in ("participant_id", "sample_accession", "biosample_accession", "sra_experiment"):
        values = [row[field] for row in participant_rows]
        if len(values) != len(set(values)):
            raise ExternalTransferFixtureError(f"GSE260666 identifier is duplicated: {field}")
    counts, stable_ids, feature_rows, mapping_counts = load_external_counts(
        matrix / "raw" / "GSE260nnn" / "GSE260666_raw_counts.txt.gz",
        activation / "gene_crosswalk.tsv",
        participant_rows,
    )
    if counts.shape != (expected_external_n, len(stable_ids)):
        raise ExternalTransferFixtureError("GSE260666 RNA matrix shape differs")
    training_fields, training_axis = read_tsv(training_molecular / "rna_feature_axis.tsv")
    participant_axis_fields, training_participants = read_tsv(
        training_molecular / "participant_axis.tsv"
    )
    if "stable_gene_id" not in training_fields or "participant_id" not in participant_axis_fields:
        raise ExternalTransferFixtureError("GSE267145 molecular axes differ")
    training_stable = [row["stable_gene_id"] for row in training_axis]
    if len(training_stable) != len(set(training_stable)):
        raise ExternalTransferFixtureError("GSE267145 stable-gene axis is duplicated")
    if len(training_participants) != expected_training_n:
        raise ExternalTransferFixtureError("GSE267145 participant count differs")
    training_values = np.load(
        training_molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False
    )
    if training_values.shape != (expected_training_n, len(training_stable)):
        raise ExternalTransferFixtureError("GSE267145 RNA matrix shape differs")
    external_aliases = {
        row[field]
        for row in participant_rows
        for field in (
            "participant_id",
            "sample_accession",
            "biosample_accession",
            "sra_experiment",
        )
    }
    external_aliases.update(("GSE260666", "PRJNA1082656"))
    training_aliases = _source_aliases_from_soft(
        training_join / "raw" / "GSE267145_family.soft.gz"
    )
    training_aliases.update(
        row[field]
        for row in training_participants
        for field in (
            "participant_id",
            "rna_source_sample_accession",
            "h3k27ac_source_sample_accession",
        )
        if field in row
    )
    alias_overlap = sorted(external_aliases & training_aliases)
    common = sorted(set(stable_ids) & set(training_stable))
    if len(common) < 2:
        raise ExternalTransferFixtureError("cross-cohort common gene axis is too small")
    external_lookup = {value: index for index, value in enumerate(stable_ids)}
    training_lookup = {value: index for index, value in enumerate(training_stable)}
    external_common = _log_cpm(counts[:, [external_lookup[value] for value in common]])
    training_common = _log_cpm(
        np.asarray(training_values[:, [training_lookup[value] for value in common]])
    )
    external_z = _row_standardize(external_common)
    training_z = _row_standardize(training_common)
    pearson = external_z @ training_z.T
    external_rank = _row_standardize(_ranks(external_common))
    training_rank = _row_standardize(_ranks(training_common))
    spearman = external_rank @ training_rank.T
    best_flat = int(np.argmax(np.minimum(pearson, spearman)))
    external_index, training_index = np.unravel_index(best_flat, pearson.shape)
    maximum_pearson = float(pearson[external_index, training_index])
    maximum_spearman = float(spearman[external_index, training_index])
    near_duplicate = maximum_pearson >= 0.999 and maximum_spearman >= 0.999
    fingerprint_key = hmac.new(
        fingerprint_salt,
        (
            "gse260666-gse267145-fingerprint-v1\0"
            + activation_artifacts_sha256
            + matrix_artifacts_sha256
            + training_molecular_artifacts_sha256
            + training_join_artifacts_sha256
        ).encode("utf-8"),
        hashlib.sha256,
    ).digest()
    prior_use = _prior_project_use(
        project_root.resolve(strict=True),
        benchmark_root.resolve(strict=True),
        sorted(external_aliases),
        fingerprint_key,
    )
    row_ids = [_row_id(row["sample_accession"]) for row in participant_rows]
    training_row_ids = [
        "g267_"
        + hashlib.sha256(
            f"gse267145_training_v1\0{row['participant_id']}".encode("utf-8")
        ).hexdigest()[:20]
        for row in training_participants
    ]
    model_root = output / "model_input"
    evaluator_root = output / "evaluator_only"
    audit_root = output / "contamination_audit"
    model_root.mkdir(parents=True)
    evaluator_root.mkdir()
    audit_root.mkdir()
    np.save(model_root / "rna_counts.npy", counts, allow_pickle=False)
    np.save(
        model_root / "rna_observed_mask.npy",
        np.ones(expected_external_n, dtype=np.bool_),
        allow_pickle=False,
    )
    write_tsv(
        model_root / "participant_axis.tsv",
        ("row_index", "row_id", "cohort_family_id", "observation_state"),
        [
            {
                "row_index": index,
                "row_id": row_id,
                "cohort_family_id": "gse260666_bulk_rna",
                "observation_state": "observed",
            }
            for index, row_id in enumerate(row_ids)
        ],
    )
    write_tsv(
        model_root / "rna_feature_axis.tsv",
        ("feature_index", "stable_gene_id", "source_feature_count", "source_symbols_sha256"),
        feature_rows,
    )
    model_receipt = {
        "schema_version": "masld-bench-gse260666-model-input-v1",
        "status": "passed_outcome_free_external_query_fixture",
        "participants": expected_external_n,
        "rna_features": len(stable_ids),
        "measurement": "raw_nonnegative_integer_gene_counts",
        "query_preprocessing_fitted": False,
        "query_participants_jointly_normalized": False,
        "labels_included": False,
        "fibrosis_state": "structurally_missing",
        "nas_state": "structurally_missing",
        "recorded_sex_state": "structurally_missing",
        "age_state": "structurally_missing",
        "mapping_counts": mapping_counts,
    }
    (model_root / "receipt.json").write_text(
        json.dumps(model_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    label_rows = []
    for row_id, row in zip(row_ids, participant_rows, strict=True):
        if row["group"] not in SOURCE_TO_STAGE3:
            raise ExternalTransferFixtureError("external source label differs")
        label_rows.append(
            {
                "row_id": row_id,
                "participant_id": row["participant_id"],
                "sample_accession": row["sample_accession"],
                "source_label": row["group"],
                "evaluation_stage3": SOURCE_TO_STAGE3[row["group"]],
            }
        )
    write_tsv(
        evaluator_root / "labels.tsv",
        ("row_id", "participant_id", "sample_accession", "source_label", "evaluation_stage3"),
        label_rows,
    )
    evaluator_receipt = {
        "schema_version": "masld-bench-gse260666-evaluator-outcomes-v1",
        "status": "passed_evaluator_only",
        "participants": expected_external_n,
        "source_class_counts": dict(sorted(source_counts.items())),
        "evaluation_class_counts": dict(
            sorted(Counter(SOURCE_TO_STAGE3[value] for value in source_counts.elements()).items())
        ),
        "mapping": SOURCE_TO_STAGE3,
        "mapping_is_one_to_one": True,
        "nafl_nash_pooled": False,
        "label_definition_mismatch": "source disease-state diagnosis versus GSE267145 histology state",
        "model_environment_access": False,
    }
    (evaluator_root / "receipt.json").write_text(
        json.dumps(evaluator_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    external_signatures = [
        _salted_signature(fingerprint_key, row_ids[index], external_common[index])
        for index in range(expected_external_n)
    ]
    training_signatures = [
        _salted_signature(fingerprint_key, training_row_ids[index], training_common[index])
        for index in range(expected_training_n)
    ]
    contamination = {
        "schema_version": "masld-bench-gse260666-contamination-audit-v1",
        "status": "blocked_duplicate" if alias_overlap or near_duplicate else "passed",
        "cohort_family_overlap": False,
        "accession_families": {
            "development_source": ["GSE267145", "GSE267119", "GSE269412"],
            "external_source": ["GSE260666", "PRJNA1082656"],
        },
        "accession_bioproject_biosample_sra_or_participant_alias_overlap_count": len(alias_overlap),
        "accession_bioproject_biosample_sra_or_participant_alias_overlap_sha256": hashlib.sha256(
            "\n".join(alias_overlap).encode("utf-8")
        ).hexdigest(),
        "external_alias_count": len(external_aliases),
        "training_alias_count": len(training_aliases),
        "prior_project_use": prior_use,
        "prior_project_use_changes_independence_counting": bool(
            prior_use["matched_file_count"]
        ),
        "external_independence_status": (
            "project_exposed_external_development"
            if prior_use["matched_file_count"]
            else "no_prior_project_use_detected"
        ),
        "common_gencode_v49_stable_genes": len(common),
        "fingerprint_method": "per-participant log1p-CPM across common stable genes; Pearson and tied-rank Spearman across genes",
        "near_duplicate_rule": "Pearson >= 0.999 and Spearman >= 0.999 for one cross-cohort pair",
        "maximum_cross_cohort_pearson": maximum_pearson,
        "maximum_cross_cohort_spearman": maximum_spearman,
        "maximum_pair_external_row_id": row_ids[external_index],
        "maximum_pair_training_row_id": training_row_ids[training_index],
        "near_duplicate_detected": near_duplicate,
        "salted_expression_fingerprints_only": True,
        "permissioned_random_salt_outside_release": True,
        "fingerprint_salt_id_sha256": hashlib.sha256(fingerprint_salt).hexdigest(),
        "external_signature_set_sha256": hashlib.sha256(
            "\n".join(sorted(external_signatures)).encode("utf-8")
        ).hexdigest(),
        "training_signature_set_sha256": hashlib.sha256(
            "\n".join(sorted(training_signatures)).encode("utf-8")
        ).hexdigest(),
        "raw_expression_fingerprints_retained": False,
        "checkpoint_exposure": {
            "gse267145_source_trained_task_native_models": "target_label_unexposed_by_construction_pending_selection_lock",
            "pretrained_cell_or_sequence_checkpoints": "not_admitted_to_this_bulk_rna_task",
        },
        "evaluation_role": "external_development_not_project_sealed",
        "champion_claim_eligible": False,
    }
    (audit_root / "audit.json").write_text(
        json.dumps(contamination, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if contamination["status"] != "passed":
        raise ExternalTransferFixtureError("cross-cohort contamination audit failed")
    summary = {
        "schema_version": "masld-bench-gse260666-external-transfer-fixture-v1",
        "status": "passed",
        "participants": expected_external_n,
        "biological_unit": "participant",
        "evaluation_role": "external_development",
        "project_sealed": False,
        "model_input_path": "model_input",
        "evaluator_only_path": "evaluator_only",
        "contamination_audit_path": "contamination_audit",
        "model_input_contains_labels": False,
        "source_fit_or_calibration_performed": False,
        "external_fit_or_calibration_performed": False,
        "external_outcomes_scored": False,
        "full_gse267145_refit_allowed": False,
        "full_refit_blocker": "gse267145_internal_selection_lock_not_supplied",
    }
    (output / "receipt.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", required=True, type=Path)
    parser.add_argument("--activation-artifacts-sha256", required=True)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--training-molecular", required=True, type=Path)
    parser.add_argument("--training-molecular-artifacts-sha256", required=True)
    parser.add_argument("--training-join", required=True, type=Path)
    parser.add_argument("--training-join-artifacts-sha256", required=True)
    parser.add_argument("--task-spec", required=True, type=Path)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--fingerprint-salt-file", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    result = build_fixture(
        activation=arguments.activation,
        matrix=arguments.matrix,
        training_molecular=arguments.training_molecular,
        training_join=arguments.training_join,
        task_spec=arguments.task_spec,
        project_root=arguments.project_root,
        benchmark_root=arguments.benchmark_root,
        fingerprint_salt_file=arguments.fingerprint_salt_file,
        output=arguments.output,
        activation_artifacts_sha256=arguments.activation_artifacts_sha256,
        matrix_artifacts_sha256=arguments.matrix_artifacts_sha256,
        training_molecular_artifacts_sha256=arguments.training_molecular_artifacts_sha256,
        training_join_artifacts_sha256=arguments.training_join_artifacts_sha256,
        task_spec_sha256=arguments.task_spec_sha256,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
