#!/usr/bin/env python3
"""Fail-closed contract for the sealed true-Kleiner fibrosis candidate bundle.

Plan 60 consumes the complete independently validated upstream bundle.  The
primary coefficient table is never accepted as a stand-alone file: READY,
validation report/checks, both manifests, and every producer artifact must be
present, hash-consistent, and semantically self-consistent.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Mapping, Sequence

from release_common import (
    CANDIDATE_ID,
    ReleaseContractError,
    clean_relative_path,
    parse_nonnegative_int,
    project_relative,
    require_sha256,
    resolve_project_path,
    sha256_file,
)


FIBROSIS_ANALYSIS_ID = "fibrosis-adjacent-true-kleiner-lvqw-v1"
FIBROSIS_EXPECTED_N_GENES = 27_638
FIBROSIS_FIXTURE_N_GENES = 8
FIBROSIS_ACCEPTED_READY_SHA256 = (
    "f6a6d274124749b5c8bbdbf204629c847a33b45bc9d2a570ba98352d895fbcc9"
)
FIBROSIS_ACCEPTED_VALIDATED_MANIFEST_SHA256 = (
    "9e8dc350e90bdeac45633007f079bcb5cf5e08f21789078adc50212bcbb64994"
)
FIBROSIS_ACCEPTED_PRIMARY_RESULT_SHA256 = (
    "01a942d33fb321631db0a3cb0579b2b89b96373c399d61f87030a49ed9788f18"
)
FIBROSIS_BUNDLE_PROJECT_REL = (
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates/"
    f"{CANDIDATE_ID}/{FIBROSIS_ANALYSIS_ID}"
)
FIBROSIS_BASE_PREFIX = "inputs/BASE/fibrosis_candidate"
FIBROSIS_PRIMARY_RELATIVE = "results/fibrosis_consecutive_lvqw_true_kleiner.csv"
FIBROSIS_PRIMARY_ARTIFACT_ID = "fibrosis_candidate_primary_result"
FIBROSIS_VALIDATION_ROLE = "fibrosis_validation_artifact"

TRUE_KLEINER_COHORTS = (
    "GSE130970",
    "GSE135251",
    "GSE162694",
    "GSE174478",
    "GSE193066",
    "GSE240729",
)
TRUE_KLEINER_TRANSITIONS = {
    "F1_vs_F0": {
        "transition": "F0_to_F1",
        "low_stage": 0,
        "high_stage": 1,
        "n_low": 126,
        "n_high": 187,
        "n_samples": 313,
        "cohorts": frozenset(TRUE_KLEINER_COHORTS),
    },
    "F2_vs_F1": {
        "transition": "F1_to_F2",
        "low_stage": 1,
        "high_stage": 2,
        "n_low": 187,
        "n_high": 174,
        "n_samples": 361,
        "cohorts": frozenset(TRUE_KLEINER_COHORTS),
    },
    "F3_vs_F2": {
        "transition": "F2_to_F3",
        "low_stage": 2,
        "high_stage": 3,
        "n_low": 174,
        "n_high": 133,
        "n_samples": 307,
        "cohorts": frozenset(TRUE_KLEINER_COHORTS),
    },
    "F4_vs_F3": {
        "transition": "F3_to_F4",
        "low_stage": 3,
        "high_stage": 4,
        "n_low": 108,
        "n_high": 44,
        "n_samples": 152,
        "cohorts": frozenset(set(TRUE_KLEINER_COHORTS) - {"GSE193066"}),
    },
}

FIBROSIS_PRODUCER_ARTIFACTS = frozenset(
    {
        FIBROSIS_PRIMARY_RELATIVE,
        "results/transition_summary.tsv",
        "audits/null_shuffle_results.csv",
        "audits/null_shuffle_summary.tsv",
        "audits/sample_manifest.tsv",
        "audits/cohort_arm_audit.tsv",
        "audits/design_audit.tsv",
        "audits/model_design.tsv",
        "audits/shuffle_audit.tsv",
        "audits/gse193066_donor_biopsy_crosswalk.tsv",
        "manifests/input_manifest.tsv",
        "manifests/code_manifest.tsv",
        "manifests/environment_manifest.tsv",
        "manifests/sessionInfo.txt",
        "manifests/run_contract.tsv",
        "manifests/producer_status.tsv",
    }
)
FIBROSIS_VALIDATED_ARTIFACTS = frozenset(
    {
        *FIBROSIS_PRODUCER_ARTIFACTS,
        "manifests/producer_artifact_manifest.tsv",
        "manifests/validation_checks.tsv",
        "manifests/validation_report.json",
    }
)
FIBROSIS_BUNDLE_FILES = tuple(
    sorted(
        {
            *FIBROSIS_VALIDATED_ARTIFACTS,
            "manifests/validated_artifact_manifest.tsv",
            "READY",
        }
    )
)
FIBROSIS_VALIDATION_CHECK_IDS = tuple(f"V{index:02d}" for index in range(1, 14))


def fibrosis_artifact_id(relative_path: str) -> str:
    """Return the deterministic BASE artifact ID for one bundle file."""

    relative = clean_relative_path(relative_path, "fibrosis bundle relative path")
    if relative == FIBROSIS_PRIMARY_RELATIVE:
        return FIBROSIS_PRIMARY_ARTIFACT_ID
    slug = re.sub(r"[^A-Za-z0-9]+", "_", relative).strip("_")
    return f"fibrosis_candidate_validation__{slug}"


def fibrosis_snapshot_relpath(relative_path: str) -> str:
    relative = clean_relative_path(relative_path, "fibrosis bundle relative path")
    return f"{FIBROSIS_BASE_PREFIX}/{relative}"


def fibrosis_artifact_role(relative_path: str) -> str:
    return (
        "fibrosis_transition_raw"
        if relative_path == FIBROSIS_PRIMARY_RELATIVE
        else FIBROSIS_VALIDATION_ROLE
    )


def _read_rows(
    path: Path, delimiter: str | None = None
) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(
            f"fibrosis bundle table is missing or unsafe: {path}"
        )
    if delimiter is None:
        delimiter = "," if path.suffix == ".csv" else "\t"
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def _read_json(path: Path, context: str) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"{context} is missing or unsafe: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseContractError(f"invalid {context}: {path}") from error
    if not isinstance(value, dict):
        raise ReleaseContractError(f"{context} must be a JSON object")
    return value


def _parse_int(value: object, context: str) -> int:
    try:
        parsed = int(str(value))
    except ValueError as error:
        raise ReleaseContractError(f"{context} is not an integer: {value!r}") from error
    return parsed


def _parse_float(value: object, context: str) -> float:
    try:
        parsed = float(str(value))
    except ValueError as error:
        raise ReleaseContractError(f"{context} is not numeric: {value!r}") from error
    if not math.isfinite(parsed):
        raise ReleaseContractError(f"{context} is not finite")
    return parsed


def _parse_bool(value: object, context: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized not in {"true", "false"}:
        raise ReleaseContractError(f"{context} is not TRUE/FALSE: {value!r}")
    return normalized == "true"


def _bh_adjust(p_values: Sequence[float]) -> list[float]:
    n_values = len(p_values)
    order = sorted(range(n_values), key=p_values.__getitem__)
    adjusted = [1.0] * n_values
    running = 1.0
    for reverse_index in range(n_values - 1, -1, -1):
        original_index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, p_values[original_index] * n_values / rank, 1.0)
        adjusted[original_index] = running
    return adjusted


def _manifest_by_path(
    path: Path,
    expected_paths: frozenset[str],
    fields: tuple[str, ...],
    context: str,
) -> dict[str, dict[str, str]]:
    observed_fields, rows = _read_rows(path, "\t")
    if tuple(observed_fields) != fields:
        raise ReleaseContractError(
            f"{context} schema drift: observed={observed_fields}, expected={list(fields)}"
        )
    by_path = {row["artifact_path"]: row for row in rows}
    if len(by_path) != len(rows) or frozenset(by_path) != expected_paths:
        raise ReleaseContractError(f"{context} artifact universe is not exact")
    return by_path


def _validate_manifests(root: Path) -> None:
    validated_path = root / "manifests/validated_artifact_manifest.tsv"
    validated = _manifest_by_path(
        validated_path,
        FIBROSIS_VALIDATED_ARTIFACTS,
        ("artifact_path", "bytes", "sha256"),
        "fibrosis validated artifact manifest",
    )
    for relative, row in validated.items():
        artifact = root / relative
        expected_hash = require_sha256(
            row["sha256"], f"fibrosis validated artifact {relative} SHA256"
        )
        expected_bytes = parse_nonnegative_int(
            row["bytes"], f"fibrosis validated artifact {relative} bytes"
        )
        if (
            not artifact.is_file()
            or artifact.is_symlink()
            or artifact.stat().st_size != expected_bytes
            or sha256_file(artifact) != expected_hash
        ):
            raise ReleaseContractError(
                f"fibrosis artifact changed after validation: {relative}"
            )

    producer_path = root / "manifests/producer_artifact_manifest.tsv"
    producer = _manifest_by_path(
        producer_path,
        FIBROSIS_PRODUCER_ARTIFACTS,
        ("artifact_path", "artifact_role", "bytes", "sha256"),
        "fibrosis producer artifact manifest",
    )
    for relative, row in producer.items():
        validated_row = validated.get(relative)
        if validated_row is None or any(
            validated_row[field] != row[field] for field in ("bytes", "sha256")
        ):
            raise ReleaseContractError(
                f"fibrosis producer/validated manifest disagreement: {relative}"
            )


def _validate_ready(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    ready = _read_json(root / "READY", "fibrosis READY")
    expected_ready_keys = {
        "status",
        "candidate_id",
        "analysis_id",
        "canonical_promotion_status",
        "validation_report_sha256",
        "validation_checks_sha256",
        "validated_artifact_manifest_sha256",
    }
    if set(ready) != expected_ready_keys:
        raise ReleaseContractError("fibrosis READY schema is not exact")
    expected_identity = {
        "status": "READY",
        "candidate_id": CANDIDATE_ID,
        "analysis_id": FIBROSIS_ANALYSIS_ID,
        "canonical_promotion_status": "not_promoted",
    }
    for key, value in expected_identity.items():
        if ready.get(key) != value:
            raise ReleaseContractError(f"fibrosis READY {key} drift")
    hash_bindings = {
        "validation_report_sha256": "manifests/validation_report.json",
        "validation_checks_sha256": "manifests/validation_checks.tsv",
        "validated_artifact_manifest_sha256": (
            "manifests/validated_artifact_manifest.tsv"
        ),
    }
    for field, relative in hash_bindings.items():
        expected = require_sha256(str(ready[field]), f"fibrosis READY {field}")
        if sha256_file(root / relative) != expected:
            raise ReleaseContractError(f"fibrosis READY hash drift for {relative}")

    report = _read_json(
        root / "manifests/validation_report.json", "fibrosis validation report"
    )
    for key, value in expected_identity.items():
        if report.get(key) != value:
            raise ReleaseContractError(f"fibrosis validation report {key} drift")
    if set(report) != {
        *expected_identity,
        "job_id",
        "n_checks",
        "summary",
        "validated_at_utc",
    }:
        raise ReleaseContractError("fibrosis validation report schema is not exact")
    if not str(report["job_id"]).isdigit():
        raise ReleaseContractError("fibrosis validation report lacks a SLURM job ID")
    fields, checks = _read_rows(root / "manifests/validation_checks.tsv", "\t")
    if fields != ["check_id", "status", "detail"]:
        raise ReleaseContractError("fibrosis validation checks schema drift")
    check_ids = [row["check_id"].split("_", 1)[0] for row in checks]
    if (
        tuple(check_ids) != FIBROSIS_VALIDATION_CHECK_IDS
        or any(row["status"] != "PASS" for row in checks)
        or _parse_int(report["n_checks"], "fibrosis report n_checks") != len(checks)
    ):
        raise ReleaseContractError(
            "fibrosis validation checks are incomplete or failed"
        )
    return ready, report


def _validate_run_contract(
    root: Path,
    expected_n_genes: int,
    project_root: Path | None,
    verify_external_sources: bool,
) -> None:
    fields, rows = _read_rows(root / "manifests/run_contract.tsv", "\t")
    if fields != ["key", "value"]:
        raise ReleaseContractError("fibrosis run contract schema drift")
    contract = {row["key"]: row["value"] for row in rows}
    if len(contract) != len(rows):
        raise ReleaseContractError("fibrosis run contract contains duplicate keys")
    expected = {
        "candidate_id": CANDIDATE_ID,
        "analysis_id": FIBROSIS_ANALYSIS_ID,
        "model": "~ dataset + inferred_sex + fib_group",
        "coefficient": "fib_grouphigh",
        "method": "limma_voom_qw_C2_true_kleiner_first_biopsy",
        "tested_universe": str(expected_n_genes),
        "significance_rule": "BH padj < 0.05 per complete contrast family",
        "lfc_gate": "none",
        "allowlist": ";".join(TRUE_KLEINER_COHORTS),
        "hard_exclusions": "GSE213621;PRJNA512027",
        "gse193066_timepoint": "1st biopsy only",
        "gse193066_donor_key": ("deposited !Sample_title with terminal _1/_2 removed"),
        "shuffle_seed": "20260808",
        "replication_unit": "one first-biopsy biological donor/sample",
        "canonical_write": "false",
    }
    if contract != expected:
        raise ReleaseContractError(
            "fibrosis run contract differs from the locked contract"
        )

    expected_manifests = {
        "input_manifest.tsv": {
            "merged_dge",
            "meta_matched",
            "sample_qc",
            "gse193066_source_metadata",
        },
        "code_manifest.tsv": {
            "producer",
            "validator",
            "slurm_wrapper",
            "lvqw_engine",
        },
    }
    for filename, expected_ids in expected_manifests.items():
        manifest_fields, manifest_rows = _read_rows(root / "manifests" / filename, "\t")
        if manifest_fields != ["artifact_id", "path", "bytes", "sha256"]:
            raise ReleaseContractError(f"fibrosis {filename} schema drift")
        by_id = {row["artifact_id"]: row for row in manifest_rows}
        if len(by_id) != len(manifest_rows) or set(by_id) != expected_ids:
            raise ReleaseContractError(f"fibrosis {filename} artifact IDs drift")
        for artifact_id, row in by_id.items():
            require_sha256(row["sha256"], f"fibrosis {filename} {artifact_id} SHA256")
            parse_nonnegative_int(
                row["bytes"], f"fibrosis {filename} {artifact_id} bytes"
            )
            clean_relative_path(row["path"], f"fibrosis {filename} {artifact_id} path")
            if verify_external_sources:
                if project_root is None:
                    raise ReleaseContractError(
                        "external fibrosis-source verification requires project_root"
                    )
                source = resolve_project_path(
                    project_root, row["path"], f"fibrosis {filename} {artifact_id}"
                )
                if (
                    not source.is_file()
                    or source.is_symlink()
                    or source.stat().st_size != int(row["bytes"])
                    or sha256_file(source) != row["sha256"]
                ):
                    raise ReleaseContractError(
                        f"fibrosis source changed after upstream validation: {artifact_id}"
                    )


def _validate_crosswalk(root: Path) -> dict[str, tuple[str, str]]:
    fields, rows = _read_rows(
        root / "audits/gse193066_donor_biopsy_crosswalk.tsv", "\t"
    )
    required = {
        "sample_id",
        "source_title",
        "donor_id",
        "biopsy",
        "is_first_biopsy",
        "included_primary",
    }
    if not required.issubset(fields) or len(rows) != 164:
        raise ReleaseContractError("fibrosis GSE193066 crosswalk schema/census drift")
    seen_samples: set[str] = set()
    donor_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    included: dict[str, tuple[str, str]] = {}
    first_count = 0
    second_count = 0
    for row in rows:
        sample_id = row["sample_id"]
        title = row["source_title"]
        if not sample_id or sample_id in seen_samples:
            raise ReleaseContractError("fibrosis crosswalk sample IDs are not unique")
        if re.fullmatch(r"HUnafld[0-9]{3}(?:_[12])?", title) is None:
            raise ReleaseContractError("fibrosis crosswalk has an invalid source title")
        expected_donor = "GSE193066::" + re.sub(r"_[12]$", "", title)
        if row["donor_id"] != expected_donor:
            raise ReleaseContractError("fibrosis crosswalk donor key drift")
        is_first = _parse_bool(row["is_first_biopsy"], "crosswalk first-biopsy flag")
        if row["biopsy"] not in {"1st biopsy", "2nd biopsy"} or is_first != (
            row["biopsy"] == "1st biopsy"
        ):
            raise ReleaseContractError("fibrosis crosswalk biopsy flag drift")
        if is_first and title.endswith("_2"):
            raise ReleaseContractError("fibrosis crosswalk first biopsy has _2 suffix")
        if not is_first and not title.endswith("_2"):
            raise ReleaseContractError(
                "fibrosis crosswalk second biopsy lacks _2 suffix"
            )
        if is_first:
            first_count += 1
        else:
            second_count += 1
        if _parse_bool(row["included_primary"], "crosswalk included flag"):
            if not is_first:
                raise ReleaseContractError(
                    "fibrosis crosswalk includes a second biopsy in the primary universe"
                )
            included[sample_id] = (expected_donor, title)
        donor_rows[expected_donor].append(row)
        seen_samples.add(sample_id)
    if first_count != 106 or second_count != 58 or len(donor_rows) != 106:
        raise ReleaseContractError(
            "fibrosis crosswalk first/second-biopsy census drift"
        )
    if sum(len(value) == 2 for value in donor_rows.values()) != 58:
        raise ReleaseContractError("fibrosis crosswalk paired-donor census drift")
    if any(
        sum(row["biopsy"] == "1st biopsy" for row in value) != 1
        or sum(row["biopsy"] == "2nd biopsy" for row in value) > 1
        for value in donor_rows.values()
    ):
        raise ReleaseContractError("fibrosis crosswalk donor/biopsy multiplicity drift")
    if len(included) != 105:
        raise ReleaseContractError(
            "fibrosis crosswalk eligible first-biopsy census drift"
        )
    return included


def _validate_sample_census(
    root: Path, included_gse193066: Mapping[str, tuple[str, str]]
) -> None:
    fields, rows = _read_rows(root / "audits/sample_manifest.tsv", "\t")
    required = {
        "sample_id",
        "donor_id",
        "donor_id_source",
        "dataset",
        "source_title",
        "biopsy",
        "fibrosis_stage",
        "transition",
        "contrast",
        "arm",
        "pass_technical",
    }
    if not required.issubset(fields):
        raise ReleaseContractError("fibrosis sample manifest schema drift")
    expected_by_transition = {
        str(spec["transition"]): (contrast, spec)
        for contrast, spec in TRUE_KLEINER_TRANSITIONS.items()
    }
    by_transition: dict[str, list[dict[str, str]]] = defaultdict(list)
    global_donor_sample: dict[str, str] = {}
    global_sample_state: dict[tuple[str, str], tuple[str, int]] = {}
    observed_gse: set[str] = set()
    for row in rows:
        transition = row["transition"]
        if transition not in expected_by_transition:
            raise ReleaseContractError(
                "fibrosis sample manifest has unknown transition"
            )
        contrast, spec = expected_by_transition[transition]
        if row["contrast"] != contrast:
            raise ReleaseContractError(
                "fibrosis sample transition/contrast mapping drift"
            )
        if row["dataset"] not in spec["cohorts"]:
            raise ReleaseContractError(
                "fibrosis sample cohort is outside true-Kleiner gate"
            )
        if not _parse_bool(row["pass_technical"], "fibrosis pass_technical"):
            raise ReleaseContractError("fibrosis sample failed technical QC")
        stage = _parse_int(row["fibrosis_stage"], "fibrosis stage")
        expected_arm = (
            "low"
            if stage == spec["low_stage"]
            else "high"
            if stage == spec["high_stage"]
            else ""
        )
        if row["arm"] != expected_arm:
            raise ReleaseContractError("fibrosis sample stage/arm drift")
        donor_id = row["donor_id"]
        sample_id = row["sample_id"]
        if not donor_id or not sample_id:
            raise ReleaseContractError("fibrosis sample/donor ID is blank")
        prior_sample = global_donor_sample.setdefault(donor_id, sample_id)
        if prior_sample != sample_id:
            raise ReleaseContractError(
                "fibrosis biological donor maps to multiple samples"
            )
        key = (row["dataset"], sample_id)
        state = (donor_id, stage)
        prior_state = global_sample_state.setdefault(key, state)
        if prior_state != state:
            raise ReleaseContractError(
                "fibrosis sample identity/stage changes across contrasts"
            )
        if row["dataset"] == "GSE193066":
            expected_gse_identity = included_gse193066.get(sample_id)
            if (
                expected_gse_identity is None
                or (row["donor_id"], row["source_title"]) != expected_gse_identity
                or row["biopsy"] != "1st biopsy"
                or row["source_title"].endswith("_2")
                or row["donor_id_source"] != "deposited_title_root"
            ):
                raise ReleaseContractError(
                    "fibrosis GSE193066 sample is not an eligible first biopsy"
                )
            observed_gse.add(sample_id)
        elif (
            row["donor_id"] != f"{row['dataset']}::{sample_id}"
            or row["donor_id_source"] != "single_biopsy_sample_id"
        ):
            raise ReleaseContractError(
                "fibrosis single-biopsy sample/donor identity drift"
            )
        by_transition[transition].append(row)

    expected_total = sum(
        int(spec["n_samples"]) for spec in TRUE_KLEINER_TRANSITIONS.values()
    )
    if len(rows) != expected_total:
        raise ReleaseContractError("fibrosis donor-by-contrast row census drift")
    for transition, (contrast, spec) in expected_by_transition.items():
        subset = by_transition[transition]
        donors = [row["donor_id"] for row in subset]
        samples = [row["sample_id"] for row in subset]
        if len(donors) != len(set(donors)) or len(samples) != len(set(samples)):
            raise ReleaseContractError(
                f"fibrosis duplicate sample/biological donor in {transition}"
            )
        if len(subset) != spec["n_samples"]:
            raise ReleaseContractError(f"fibrosis {transition} sample count drift")
        if {row["dataset"] for row in subset} != spec["cohorts"]:
            raise ReleaseContractError(f"fibrosis {transition} cohort membership drift")
        for cohort in spec["cohorts"]:
            cohort_rows = [row for row in subset if row["dataset"] == cohort]
            if not any(row["arm"] == "low" for row in cohort_rows) or not any(
                row["arm"] == "high" for row in cohort_rows
            ):
                raise ReleaseContractError(
                    f"fibrosis {transition}/{cohort} lacks both adjacent-stage arms"
                )
        if sum(row["arm"] == "low" for row in subset) != spec["n_low"]:
            raise ReleaseContractError(f"fibrosis {transition} low-arm count drift")
        if sum(row["arm"] == "high" for row in subset) != spec["n_high"]:
            raise ReleaseContractError(f"fibrosis {transition} high-arm count drift")
    if len(global_donor_sample) != 664:
        raise ReleaseContractError(
            "fibrosis exact-stage biological-donor universe drift"
        )
    if observed_gse != set(included_gse193066):
        raise ReleaseContractError(
            "fibrosis GSE193066 eligible first-biopsy universe is incomplete"
        )


def _validate_results(
    root: Path,
    expected_n_genes: int,
    report: Mapping[str, object],
) -> dict[str, dict[str, int]]:
    fields, rows = _read_rows(root / FIBROSIS_PRIMARY_RELATIVE, ",")
    required = {
        "gene",
        "logFC",
        "SE",
        "t",
        "P.Value",
        "padj",
        "contrast",
        "transition",
        "method",
        "n_cohorts",
        "n_samples",
        "n_low",
        "n_high",
    }
    if not required.issubset(fields):
        raise ReleaseContractError("fibrosis primary result schema drift")
    if len(rows) != expected_n_genes * len(TRUE_KLEINER_TRANSITIONS):
        raise ReleaseContractError(
            "fibrosis primary result is not the exact complete contrast family"
        )
    genes: dict[str, set[str]] = defaultdict(set)
    p_values: dict[str, list[float]] = defaultdict(list)
    q_values: dict[str, list[float]] = defaultdict(list)
    summary: dict[str, dict[str, int]] = {
        contrast: {"n_genes_tested": 0, "n_deg_05": 0, "n_up": 0, "n_down": 0}
        for contrast in TRUE_KLEINER_TRANSITIONS
    }
    for row in rows:
        contrast = row["contrast"]
        if contrast not in TRUE_KLEINER_TRANSITIONS:
            raise ReleaseContractError("fibrosis primary result has unknown contrast")
        spec = TRUE_KLEINER_TRANSITIONS[contrast]
        expected_pairs = {
            "transition": spec["transition"],
            "method": "limma_voom_qw_C2_true_kleiner_first_biopsy",
            "n_samples": str(spec["n_samples"]),
            "n_cohorts": str(len(spec["cohorts"])),
            "n_low": str(spec["n_low"]),
            "n_high": str(spec["n_high"]),
        }
        if any(row[field] != value for field, value in expected_pairs.items()):
            raise ReleaseContractError(
                f"fibrosis result provenance/census drift: {contrast}"
            )
        gene = row["gene"]
        if not gene or gene in genes[contrast]:
            raise ReleaseContractError(f"fibrosis blank/duplicate gene in {contrast}")
        genes[contrast].add(gene)
        logfc = _parse_float(row["logFC"], f"fibrosis {contrast} logFC")
        se = _parse_float(row["SE"], f"fibrosis {contrast} SE")
        statistic = _parse_float(row["t"], f"fibrosis {contrast} t")
        p_value = _parse_float(row["P.Value"], f"fibrosis {contrast} p")
        q_value = _parse_float(row["padj"], f"fibrosis {contrast} q")
        if se <= 0 or not 0 <= p_value <= 1 or not 0 <= q_value <= 1:
            raise ReleaseContractError("fibrosis coefficient statistics are invalid")
        if abs(logfc - statistic * se) > 5e-10 * (1 + abs(logfc)):
            raise ReleaseContractError("fibrosis logFC/SE/t identity drift")
        p_values[contrast].append(p_value)
        q_values[contrast].append(q_value)
        observed = summary[contrast]
        observed["n_genes_tested"] += 1
        if q_value < 0.05:
            if logfc == 0:
                raise ReleaseContractError("fibrosis significant result has zero logFC")
            observed["n_deg_05"] += 1
            observed["n_up" if logfc > 0 else "n_down"] += 1
    reference_genes: set[str] | None = None
    for contrast in TRUE_KLEINER_TRANSITIONS:
        if len(genes[contrast]) != expected_n_genes:
            raise ReleaseContractError(f"fibrosis {contrast} tested family drift")
        if reference_genes is None:
            reference_genes = genes[contrast]
        elif genes[contrast] != reference_genes:
            raise ReleaseContractError(
                "fibrosis gene universe differs across contrasts"
            )
        rederived = _bh_adjust(p_values[contrast])
        if (
            max(
                abs(observed - expected)
                for observed, expected in zip(
                    q_values[contrast], rederived, strict=True
                )
            )
            > 5e-12
        ):
            raise ReleaseContractError(f"fibrosis {contrast} BH rederivation drift")

    summary_fields, summary_rows = _read_rows(
        root / "results/transition_summary.tsv", "\t"
    )
    if not {
        "contrast",
        "n_genes_tested",
        "n_deg_05",
        "n_up",
        "n_down",
        "n_samples",
        "n_cohorts",
        "significance_rule",
        "lfc_gate",
    }.issubset(summary_fields):
        raise ReleaseContractError("fibrosis transition summary schema drift")
    by_contrast = {row["contrast"]: row for row in summary_rows}
    if len(by_contrast) != len(summary_rows) or set(by_contrast) != set(
        TRUE_KLEINER_TRANSITIONS
    ):
        raise ReleaseContractError(
            "fibrosis transition summary contrast universe drift"
        )
    report_summary = report.get("summary")
    if not isinstance(report_summary, dict) or set(report_summary) != set(
        TRUE_KLEINER_TRANSITIONS
    ):
        raise ReleaseContractError("fibrosis validation report summary universe drift")
    for contrast, observed in summary.items():
        spec = TRUE_KLEINER_TRANSITIONS[contrast]
        row = by_contrast[contrast]
        expected = {
            **observed,
            "n_samples": int(spec["n_samples"]),
            "n_cohorts": len(spec["cohorts"]),
        }
        for field, value in expected.items():
            if _parse_int(row[field], f"fibrosis summary {contrast} {field}") != value:
                raise ReleaseContractError(
                    f"fibrosis transition summary disagrees with raw result: {contrast}"
                )
        if row["significance_rule"] != "BH padj < 0.05" or row["lfc_gate"] != "none":
            raise ReleaseContractError("fibrosis transition summary threshold drift")
        report_row = report_summary[contrast]
        if not isinstance(report_row, dict):
            raise ReleaseContractError(
                "fibrosis report contrast summary is not an object"
            )
        for field, value in expected.items():
            if (
                _parse_int(
                    report_row.get(field, ""), f"fibrosis report {contrast} {field}"
                )
                != value
            ):
                raise ReleaseContractError(
                    f"fibrosis validation report disagrees with raw result: {contrast}"
                )
        if "null_n_deg_05" not in report_row:
            raise ReleaseContractError("fibrosis report omits shuffled-null count")
        _parse_int(report_row["null_n_deg_05"], "fibrosis shuffled-null count")
    return summary


def _validate_status(
    root: Path, expected_n_genes: int, report: Mapping[str, object]
) -> None:
    fields, rows = _read_rows(root / "manifests/producer_status.tsv", "\t")
    if len(rows) != 1 or not {
        "status",
        "candidate_id",
        "analysis_id",
        "n_contrasts",
        "n_result_rows",
        "n_unique_primary_donors",
        "canonical_write",
    }.issubset(fields):
        raise ReleaseContractError("fibrosis producer status schema/census drift")
    row = rows[0]
    expected = {
        "status": "PRODUCER_COMPLETE_PENDING_VALIDATION",
        "candidate_id": CANDIDATE_ID,
        "analysis_id": FIBROSIS_ANALYSIS_ID,
        "n_contrasts": "4",
        "n_result_rows": str(expected_n_genes * 4),
        "n_unique_primary_donors": "664",
        "canonical_write": "FALSE",
    }
    if any(row[field] != value for field, value in expected.items()):
        raise ReleaseContractError("fibrosis producer status identity/census drift")
    fields, environment = _read_rows(root / "manifests/environment_manifest.tsv", "\t")
    if not environment or "slurm_job_id" not in fields:
        raise ReleaseContractError("fibrosis environment manifest schema drift")
    job_ids = {row["slurm_job_id"] for row in environment}
    if len(job_ids) != 1 or next(iter(job_ids)) != str(report["job_id"]):
        raise ReleaseContractError("fibrosis report/environment SLURM identity drift")


def validate_fibrosis_candidate_bundle(
    bundle_root: Path,
    *,
    expected_n_genes: int = FIBROSIS_EXPECTED_N_GENES,
    project_root: Path | None = None,
    verify_external_sources: bool = False,
) -> dict[str, object]:
    """Validate a live or staged sealed candidate bundle without writing."""

    requested_root = Path(bundle_root.absolute())
    if requested_root.is_symlink():
        raise ReleaseContractError("fibrosis candidate bundle root is unsafe")
    root = requested_root.resolve(strict=True)
    if not root.is_dir():
        raise ReleaseContractError("fibrosis candidate bundle root is unsafe")
    if expected_n_genes <= 0:
        raise ReleaseContractError(
            "fibrosis expected gene-family size must be positive"
        )
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ReleaseContractError(
                f"fibrosis candidate bundle contains symlink: {path}"
            )
    observed = {
        path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file()
    }
    if observed != set(FIBROSIS_BUNDLE_FILES):
        raise ReleaseContractError(
            "fibrosis candidate bundle file universe is not exact: "
            f"missing={sorted(set(FIBROSIS_BUNDLE_FILES) - observed)}, "
            f"extra={sorted(observed - set(FIBROSIS_BUNDLE_FILES))}"
        )
    _validate_manifests(root)
    _, report = _validate_ready(root)
    _validate_run_contract(
        root, expected_n_genes, project_root, verify_external_sources
    )
    included = _validate_crosswalk(root)
    _validate_sample_census(root, included)
    summary = _validate_results(root, expected_n_genes, report)
    _validate_status(root, expected_n_genes, report)
    if expected_n_genes == FIBROSIS_EXPECTED_N_GENES:
        accepted_hashes = {
            "READY": FIBROSIS_ACCEPTED_READY_SHA256,
            "manifests/validated_artifact_manifest.tsv": (
                FIBROSIS_ACCEPTED_VALIDATED_MANIFEST_SHA256
            ),
            FIBROSIS_PRIMARY_RELATIVE: FIBROSIS_ACCEPTED_PRIMARY_RESULT_SHA256,
        }
        for relative, accepted_sha256 in accepted_hashes.items():
            if sha256_file(root / relative) != accepted_sha256:
                raise ReleaseContractError(
                    f"fibrosis bundle differs from accepted upstream release: {relative}"
                )
    return {
        "status": "READY",
        "candidate_id": CANDIDATE_ID,
        "analysis_id": FIBROSIS_ANALYSIS_ID,
        "canonical_promotion_status": "not_promoted",
        "bundle_root": str(root),
        "bundle_file_count": len(FIBROSIS_BUNDLE_FILES),
        "validated_manifest_sha256": sha256_file(
            root / "manifests/validated_artifact_manifest.tsv"
        ),
        "ready_sha256": sha256_file(root / "READY"),
        "primary_result_sha256": sha256_file(root / FIBROSIS_PRIMARY_RELATIVE),
        "n_genes_per_contrast": expected_n_genes,
        "summary": summary,
    }


def validate_fibrosis_selection_rows(
    project_root: Path,
    selection_rows: Sequence[Mapping[str, str]],
    *,
    fixture_mode: bool,
    verify_external_sources: bool,
) -> dict[str, object]:
    """Require an exact one-row-per-file selection of the sealed bundle."""

    selected = [
        row
        for row in selection_rows
        if row.get("artifact_role")
        in {"fibrosis_transition_raw", FIBROSIS_VALIDATION_ROLE}
        or str(row.get("snapshot_relpath", "")).startswith(f"{FIBROSIS_BASE_PREFIX}/")
    ]
    if len(selected) != len(FIBROSIS_BUNDLE_FILES):
        raise ReleaseContractError(
            "BASE selection must contain the complete sealed fibrosis bundle"
        )
    expected_ids = {fibrosis_artifact_id(path) for path in FIBROSIS_BUNDLE_FILES}
    observed_ids = {str(row["artifact_id"]) for row in selected}
    if len(observed_ids) != len(selected) or observed_ids != expected_ids:
        raise ReleaseContractError("sealed fibrosis BASE artifact IDs are not exact")
    primary = [
        row for row in selected if row["artifact_role"] == "fibrosis_transition_raw"
    ]
    validation = [
        row for row in selected if row["artifact_role"] == FIBROSIS_VALIDATION_ROLE
    ]
    if len(primary) != 1 or len(validation) != len(FIBROSIS_BUNDLE_FILES) - 1:
        raise ReleaseContractError(
            "sealed fibrosis bundle requires exactly one primary result role"
        )

    source_roots: set[Path] = set()
    by_relative: dict[str, Mapping[str, str]] = {}
    for row in selected:
        snapshot = clean_relative_path(
            str(row["snapshot_relpath"]), "fibrosis BASE snapshot path"
        )
        prefix = f"{FIBROSIS_BASE_PREFIX}/"
        if not snapshot.startswith(prefix):
            raise ReleaseContractError("sealed fibrosis bundle escaped its BASE prefix")
        relative = snapshot[len(prefix) :]
        if relative not in FIBROSIS_BUNDLE_FILES or relative in by_relative:
            raise ReleaseContractError(
                "sealed fibrosis bundle path universe is not exact"
            )
        if row["artifact_id"] != fibrosis_artifact_id(relative) or row[
            "artifact_role"
        ] != fibrosis_artifact_role(relative):
            raise ReleaseContractError("sealed fibrosis artifact identity/role drift")
        source = resolve_project_path(
            project_root,
            str(row["source_path"]),
            f"fibrosis BASE source {relative}",
        )
        relative_parts = Path(relative).parts
        source_root = source
        for _ in relative_parts:
            source_root = source_root.parent
        if source != source_root.joinpath(*relative_parts):
            raise ReleaseContractError("fibrosis source path cannot be reconstructed")
        source_roots.add(source_root)
        by_relative[relative] = row
    if set(by_relative) != set(FIBROSIS_BUNDLE_FILES) or len(source_roots) != 1:
        raise ReleaseContractError("fibrosis BASE sources do not form one exact bundle")
    source_root = next(iter(source_roots))
    if not fixture_mode:
        expected_root = (project_root.resolve() / FIBROSIS_BUNDLE_PROJECT_REL).resolve()
        if source_root != expected_root:
            raise ReleaseContractError(
                "fibrosis BASE source is not the fixed true-Kleiner candidate root"
            )
    expected_n_genes = (
        FIBROSIS_FIXTURE_N_GENES if fixture_mode else FIBROSIS_EXPECTED_N_GENES
    )
    contract = validate_fibrosis_candidate_bundle(
        source_root,
        expected_n_genes=expected_n_genes,
        project_root=project_root,
        verify_external_sources=verify_external_sources,
    )
    for relative, row in by_relative.items():
        source = source_root / relative
        if (
            row["source_sha256"] != sha256_file(source)
            or _parse_int(row["source_bytes"], f"fibrosis {relative} selected bytes")
            != source.stat().st_size
        ):
            raise ReleaseContractError(
                f"fibrosis selected source hash/byte drift: {relative}"
            )
    return contract


def validate_frozen_fibrosis_bundle(
    candidate_root: Path,
    *,
    fixture_mode: bool,
) -> dict[str, object]:
    """Revalidate the complete copied bundle from BASE during REL02/REL05."""

    return validate_fibrosis_candidate_bundle(
        candidate_root / FIBROSIS_BASE_PREFIX,
        expected_n_genes=(
            FIBROSIS_FIXTURE_N_GENES if fixture_mode else FIBROSIS_EXPECTED_N_GENES
        ),
    )


def fixed_live_bundle_root(project_root: Path) -> Path:
    return project_root.resolve() / FIBROSIS_BUNDLE_PROJECT_REL


def fibrosis_source_provenance(project_root: Path, relative_path: str) -> str:
    return project_relative(
        project_root.resolve(), fixed_live_bundle_root(project_root) / relative_path
    )
