#!/usr/bin/env python3
"""Independently validate and seal the no-fit Resource stage preflight."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import stat
import subprocess
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

from resource_contract import (
    CANDIDATE_ID,
    CANDIDATE_ROOT,
    EXPECTED_EXECUTION_LOCALE,
    EXPECTED_EXECUTION_PARTITION,
    EXPECTED_EXECUTION_QOS,
    RNASEQ_PYTHON,
    RNASEQ_PYTHON_SHA256,
    ContractError,
    require,
    require_regular_file,
    require_single_slurm_node,
    sha256,
)
from snapshot_io import (
    SnapshotIOError,
    chmod_nofollow,
    publish_directory_noreplace,
    publish_file_noreplace,
)
from verify_candidate_snapshot import verify as verify_candidate_snapshot


PREFLIGHT_ROOT = CANDIDATE_ROOT / "workstreams/BULK-STAGE/PREFLIGHT-v1"
BG001_ROOT = CANDIDATE_ROOT / "inputs/BG001-DECISION"
ARM_ROOT = BG001_ROOT / "arms/F_legacy"
UNIFIED = (
    BG001_ROOT
    / "source_snapshot/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
)
GSE193066 = (
    CANDIDATE_ROOT
    / "inputs/BULK-NINE-COHORT/source_metadata/GSE193066_metadata.tsv"
)
QC = ARM_ROOT / "qc/sample_qc_report.csv"
RDS_HELPER = CANDIDATE_ROOT / "code/validate_stage_rds_membership.R"
VALIDATION_ROOT = PREFLIGHT_ROOT / "validation"
VALIDATION_MANIFEST = VALIDATION_ROOT / "validated_artifact_manifest.tsv"
READY_PATH = PREFLIGHT_ROOT / "STAGE_PREFLIGHT_READY.json"

COHORTS = (
    "GSE130970",
    "GSE135251",
    "GSE162694",
    "GSE174478",
    "GSE193066",
    "GSE240729",
)
STAGE_COUNTS = {0: 126, 1: 187, 2: 174, 3: 132, 4: 42}
TRANSITIONS = {
    "F0_to_F1": {
        "contrast": "F1_vs_F0",
        "low": 0,
        "high": 1,
        "n_low": 126,
        "n_high": 187,
        "n_samples": 313,
        "cohorts": COHORTS,
        "rank": 8,
    },
    "F1_to_F2": {
        "contrast": "F2_vs_F1",
        "low": 1,
        "high": 2,
        "n_low": 187,
        "n_high": 174,
        "n_samples": 361,
        "cohorts": COHORTS,
        "rank": 8,
    },
    "F2_to_F3": {
        "contrast": "F3_vs_F2",
        "low": 2,
        "high": 3,
        "n_low": 174,
        "n_high": 132,
        "n_samples": 306,
        "cohorts": COHORTS,
        "rank": 8,
    },
    "F3_to_F4": {
        "contrast": "F4_vs_F3",
        "low": 3,
        "high": 4,
        "n_low": 107,
        "n_high": 42,
        "n_samples": 149,
        "cohorts": tuple(item for item in COHORTS if item != "GSE193066"),
        "rank": 7,
    },
}
EXPECTED_HASHES = {
    "merged_counts_raw": "16afc00bf1db706731d07df1dfeabdad0bc423b1caa1b822631685a55ec7f225",
    "merged_dge": "56d1cf97f791a81404ddc6c5ea4c5df9e6a38e7cf011b39d0d7c3279c87afa4d",
    "meta_matched": "90ba6aca68643c4680c08fc982a6c72b76edb740ab726874cb50603711756a83",
    "sample_qc": "1c7fef991c599554369ecc79de2b2f441597cf9b6c1a7ff3ecbe05f1c5969560",
    "unified_metadata": "b9f6afc1e916391951c638a60d8ec95f1b380d20311469c0830f638f16649699",
    "gse193066_metadata": "b66158d20da2127e8c3cc7ebae4f1d85a2d6813b7dd6f2af4bdac60070ff2bc5",
}


def read_table(path: Path, delimiter: str = "\t") -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file() and not path.is_symlink(), f"missing or symlinked table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    require(
        all(None not in row and all(value is not None for value in row.values()) for row in rows),
        f"malformed table row width: {path}",
    )
    return fields, rows


def write_exclusive(path: Path, payload: str) -> None:
    require(not path.exists() and not path.is_symlink(), f"refusing overwrite: {path}")
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def tsv_payload(rows: list[dict[str, object]], fields: list[str]) -> str:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def parse_bool(value: str) -> bool:
    normalized = value.strip().upper()
    require(normalized in {"TRUE", "FALSE"}, f"invalid boolean value: {value}")
    return normalized == "TRUE"


def matrix_rank(matrix: list[list[float]], tolerance: float = 1e-10) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    n_rows = len(work)
    n_cols = len(work[0])
    rank = 0
    for column in range(n_cols):
        pivot = max(range(rank, n_rows), key=lambda row: abs(work[row][column]), default=None)
        if pivot is None or abs(work[pivot][column]) <= tolerance:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(n_rows):
            if row == rank:
                continue
            factor = work[row][column]
            if abs(factor) <= tolerance:
                continue
            work[row] = [
                current - factor * pivot_current
                for current, pivot_current in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == n_rows:
            break
    return rank


def load_json(path: Path) -> dict[str, object]:
    require_regular_file(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def derive_rds_membership() -> dict[str, dict[str, str]]:
    runtime_path = BG001_ROOT / "contract/analysis_runtime_contract.json"
    runtime = load_json(runtime_path)
    require(isinstance(runtime.get("Rscript"), dict), "runtime Rscript schema drift")
    rscript = Path(str(runtime["Rscript"]["path"])).resolve(strict=True)
    require_regular_file(rscript)
    require(
        sha256(rscript) == runtime["Rscript"]["sha256"],
        "stage validation Rscript drift",
    )
    require_regular_file(RDS_HELPER)
    result = subprocess.run(
        [
            str(rscript),
            "--vanilla",
            str(RDS_HELPER),
            str(ARM_ROOT / "results/integration/merged_counts_raw.rds"),
            str(ARM_ROOT / "results/integration/merged_dge.rds"),
            str(ARM_ROOT / "results/integration/meta_matched.rds"),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    require(
        result.returncode == 0,
        f"independent RDS membership audit failed: {result.stderr.strip()}",
    )
    reader = csv.DictReader(result.stdout.splitlines(), delimiter="	")
    require(
        list(reader.fieldnames or [])
        == [
            "sample_id",
            "inferred_sex",
            "raw_present",
            "corrected_dge_present",
        ],
        "independent RDS membership schema drift",
    )
    rows = list(reader)
    require(
        len(rows) == 1281
        and all(None not in row for row in rows),
        "independent RDS membership row drift",
    )
    by_sample = {row["sample_id"]: row for row in rows}
    require(len(by_sample) == 1281, "duplicate independent RDS sample")
    require(
        sum(parse_bool(row["raw_present"]) for row in rows) == 1281
        and sum(parse_bool(row["corrected_dge_present"]) for row in rows) == 1257,
        "independent RDS membership census drift",
    )
    return by_sample


def derive_source_stage_set() -> tuple[
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
    dict[str, dict[str, str]],
]:
    _, unified_rows = read_table(UNIFIED, ",")
    _, qc_rows = read_table(QC, ",")
    _, gse_rows = read_table(GSE193066, "\t")
    require(len(unified_rows) == 1469, "frozen unified metadata row-count drift")
    require(len(qc_rows) == 1281, "corrected QC row-count drift")
    require(len(gse_rows) == 164, "GSE193066 source row-count drift")
    qc_by_sample = {row["sample_id"]: row for row in qc_rows}
    require(len(qc_by_sample) == 1281, "duplicate QC sample")
    require(
        len({row["sample_id"] for row in unified_rows}) == 1469,
        "duplicate unified-metadata sample",
    )
    rds_by_sample = derive_rds_membership()

    gse_by_sample: dict[str, dict[str, object]] = {}
    participants: defaultdict[str, list[str]] = defaultdict(list)
    for row in gse_rows:
        sample = row["sample_id"]
        require(sample == row["Run"], "GSE193066 Run/sample_id mismatch")
        title = row["!Sample_title"]
        require(
            re.fullmatch(r"HUnafld[0-9]{3}(?:_[12])?", title) is not None,
            "GSE193066 title regex drift",
        )
        participant = title[:-2] if title.endswith(("_1", "_2")) else title
        require(row["biopsy"] in {"1st biopsy", "2nd biopsy"}, "GSE193066 biopsy label drift")
        require(row["fibrosis stage"] in {"0", "1", "2", "3", "4"}, "GSE193066 stage drift")
        require(sample not in gse_by_sample, "duplicate GSE193066 sample")
        require(
            (row["biopsy"] == "2nd biopsy" and title.endswith("_2"))
            or (
                row["biopsy"] == "1st biopsy"
                and not title.endswith("_2")
            ),
            "GSE193066 biopsy/title suffix drift",
        )
        gse_by_sample[sample] = {
            "participant": participant,
            "title": title,
            "biopsy": row["biopsy"],
            "stage": row["fibrosis stage"],
        }
        participants[participant].append(row["biopsy"])
    require(
        len({row["title"] for row in gse_by_sample.values()}) == 164,
        "duplicate GSE193066 source title",
    )
    require(len(participants) == 106, "GSE193066 participant count drift")
    require(sum(sorted(value) == ["1st biopsy", "2nd biopsy"] for value in participants.values()) == 58,
            "GSE193066 paired participant count drift")
    require(sum(value["biopsy"] == "1st biopsy" for value in gse_by_sample.values()) == 106,
            "GSE193066 first-biopsy count drift")

    candidates: dict[str, dict[str, object]] = {}
    eligibility_expected: dict[str, dict[str, str]] = {}
    eligible_rows = 0
    exclusion_reasons: Counter[str] = Counter()

    def r_value(value: str) -> str:
        return value if value not in {"", "NA"} else "NA"

    def r_bool(value: bool) -> str:
        return "TRUE" if value else "FALSE"

    for row in unified_rows:
        dataset = row["dataset"]
        if dataset not in COHORTS or row["fibrosis_stage"] not in {"0", "1", "2", "3", "4"}:
            continue
        eligible_rows += 1
        sample = row["sample_id"]
        require(sample not in eligibility_expected, f"duplicate eligible source sample: {sample}")
        if dataset == "GSE193066":
            require(sample in gse_by_sample, f"GSE193066 stage sample absent from source key: {sample}")
            source = gse_by_sample[sample]
            require(source["stage"] == row["fibrosis_stage"], "GSE193066 source/harmonized stage mismatch")
            is_first = source["biopsy"] == "1st biopsy"
            analysis_unit = "GSE193066::{}".format(source["participant"])
            title = str(source["title"])
            biopsy = str(source["biopsy"])
            key_source = "deposited_title_root"
        else:
            source = None
            is_first = True
            analysis_unit = f"{dataset}::{sample}"
            title = sample
            biopsy = "single deposited biopsy"
            key_source = "single_deposited_sample_accession"

        qc_present = sample in qc_by_sample
        rds_present = sample in rds_by_sample
        pass_technical = (
            parse_bool(qc_by_sample[sample]["pass_technical"])
            if qc_present
            else False
        )
        raw_present = (
            parse_bool(rds_by_sample[sample]["raw_present"])
            if rds_present
            else False
        )
        dge_present = (
            parse_bool(rds_by_sample[sample]["corrected_dge_present"])
            if rds_present
            else False
        )
        if dataset == "GSE193066" and not is_first:
            reason = "excluded_gse193066_second_biopsy"
        elif not qc_present:
            reason = "excluded_absent_from_corrected_qc"
        elif not rds_present:
            reason = "excluded_absent_from_matched_metadata"
        elif not raw_present:
            reason = "excluded_absent_from_raw_matrix"
        elif not pass_technical:
            reason = "excluded_failed_corrected_technical_qc"
        elif not dge_present:
            reason = "excluded_absent_from_corrected_dge"
        else:
            reason = "included_primary"
        included = reason == "included_primary"
        exclusion_reasons[reason] += 1

        expected = {
            **{
                # fread preserves empty character fields, but promotes empty
                # numeric fields to NA before write.table serializes them.
                field: (
                    r_value(row[field])
                    if field in {"age", "fibrosis_stage", "nas_score"}
                    else row[field]
                )
                for field in (
                    "sample_id",
                    "dataset",
                    "condition",
                    "group_binary",
                    "sex",
                    "age",
                    "fibrosis_stage",
                    "nas_score",
                    "diagnosis_harmonized",
                )
            },
            "pass_technical": r_bool(pass_technical) if qc_present else "NA",
            "inferred_sex": (
                r_value(rds_by_sample[sample]["inferred_sex"])
                if rds_present
                else "NA"
            ),
            "matched_dataset": dataset if rds_present else "NA",
            "matched_fibrosis_stage": row["fibrosis_stage"] if rds_present else "NA",
            "source_title": title,
            "biopsy": biopsy,
            "analysis_unit_id": analysis_unit,
            "is_first_biopsy": r_bool(is_first),
            "analysis_unit_key_source": key_source,
            "raw_present": r_bool(raw_present),
            "corrected_dge_present": r_bool(dge_present),
            "qc_present": r_bool(qc_present),
            "meta_present": r_bool(rds_present),
            "included_primary": r_bool(included),
            "exclusion_reason": reason,
        }
        eligibility_expected[sample] = expected

        if source is not None:
            require(
                qc_present and rds_present,
                f"GSE193066 sample lacks corrected QC/RDS metadata: {sample}",
            )
            crosswalk_included = (
                is_first and pass_technical and raw_present and dge_present
            )
            if crosswalk_included:
                crosswalk_reason = "included_first_biopsy"
            elif source["biopsy"] == "2nd biopsy":
                crosswalk_reason = "excluded_second_biopsy"
            elif not pass_technical:
                crosswalk_reason = "excluded_failed_corrected_technical_qc"
            elif not raw_present:
                crosswalk_reason = "excluded_absent_from_raw_matrix"
            elif not dge_present:
                crosswalk_reason = "excluded_absent_from_corrected_dge"
            else:
                crosswalk_reason = "excluded_unclassified"
            source.update(
                {
                    "inferred_sex": rds_by_sample[sample]["inferred_sex"],
                    "pass_technical": pass_technical,
                    "raw_present": raw_present,
                    "corrected_dge_present": dge_present,
                    "included_primary": crosswalk_included,
                    "exclusion_reason": crosswalk_reason,
                }
            )

        if included:
            require(sample not in candidates, f"duplicate selected sample: {sample}")
            candidates[sample] = {
                "sample_id": sample,
                "analysis_unit_id": analysis_unit,
                "analysis_unit_key_source": key_source,
                "dataset": dataset,
                "source_title": title,
                "biopsy": biopsy,
                "fibrosis_stage": int(row["fibrosis_stage"]),
                "inferred_sex": rds_by_sample[sample]["inferred_sex"],
            }
    require(eligible_rows == 731, "six-cohort exact-stage source census drift")
    require(
        exclusion_reasons
        == Counter(
            {
                "included_primary": 661,
                "excluded_gse193066_second_biopsy": 58,
                "excluded_absent_from_corrected_qc": 2,
                "excluded_failed_corrected_technical_qc": 10,
            }
        ),
        f"stage eligibility census drift: {dict(exclusion_reasons)}",
    )
    require(len(candidates) == 661, "corrected stage sample count drift")
    require(len({row["analysis_unit_id"] for row in candidates.values()}) == 661,
            "duplicate stage analysis unit")
    require(
        all("included_primary" in row for row in gse_by_sample.values()),
        "GSE193066 source did not join corrected RDS/QC state",
    )
    return candidates, gse_by_sample, eligibility_expected


def verify_producer_manifest() -> None:
    manifest_path = PREFLIGHT_ROOT / "manifests/producer_artifact_manifest.tsv"
    fields, rows = read_table(manifest_path)
    require(fields == ["candidate_path", "size_bytes", "sha256"], "producer manifest schema drift")
    observed: set[str] = set()
    for row in rows:
        relative = PurePosixPath(row["candidate_path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "unsafe producer path")
        text = relative.as_posix()
        require(text not in observed, f"duplicate producer artifact: {text}")
        observed.add(text)
        path = PREFLIGHT_ROOT.joinpath(*relative.parts)
        require(path.is_file() and not path.is_symlink(), f"missing producer artifact: {text}")
        require(path.stat().st_size == int(row["size_bytes"]), f"producer byte drift: {text}")
        require(sha256(path) == row["sha256"], f"producer hash drift: {text}")
    excluded = {
        manifest_path,
        READY_PATH,
    }
    actual = {
        path.relative_to(PREFLIGHT_ROOT).as_posix()
        for path in PREFLIGHT_ROOT.rglob("*")
        if path.is_file()
        and path not in excluded
        and VALIDATION_ROOT not in path.parents
    }
    require(observed == actual, "producer manifest file-set mismatch")


def validate_input_and_code_manifests() -> None:
    expected_inputs = [
        (
            "merged_counts_raw",
            "inputs/BG001-DECISION/arms/F_legacy/results/integration/"
            "merged_counts_raw.rds",
            "primary_stage_native_count_source",
        ),
        (
            "merged_dge",
            "inputs/BG001-DECISION/arms/F_legacy/results/integration/"
            "merged_dge.rds",
            "full_resource_filter_reference_only",
        ),
        (
            "meta_matched",
            "inputs/BG001-DECISION/arms/F_legacy/results/integration/"
            "meta_matched.rds",
            "corrected_sample_metadata",
        ),
        (
            "sample_qc",
            "inputs/BG001-DECISION/arms/F_legacy/qc/sample_qc_report.csv",
            "corrected_fragment_qc",
        ),
        (
            "unified_metadata",
            "inputs/BG001-DECISION/source_snapshot/RNA-seq/Human/"
            "Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv",
            "frozen_harmonized_stage_source",
        ),
        (
            "gse193066_metadata",
            "inputs/BULK-NINE-COHORT/source_metadata/GSE193066_metadata.tsv",
            "authoritative_gse193066_participant_key",
        ),
    ]
    input_fields, input_rows = read_table(
        PREFLIGHT_ROOT / "manifests/input_manifest.tsv"
    )
    require(
        input_fields
        == ["input_id", "candidate_path", "size_bytes", "sha256", "role"],
        "input manifest schema drift",
    )
    require(
        [
            (row["input_id"], row["candidate_path"], row["role"])
            for row in input_rows
        ]
        == expected_inputs,
        "input ID/path/role/order drift",
    )

    _, source_rows = read_table(
        CANDIDATE_ROOT / "manifests/source_selection.tsv"
    )
    source_by_candidate = {
        row["candidate_path"]: row
        for row in source_rows
    }
    require(
        len(source_by_candidate) == len(source_rows),
        "duplicate candidate source-selection path",
    )
    for row in input_rows:
        path = CANDIDATE_ROOT / row["candidate_path"]
        require_regular_file(path)
        source = source_by_candidate.get(row["candidate_path"])
        require(source is not None, f"stage input absent from source selection: {row['input_id']}")
        require(
            path.stat().st_size == int(row["size_bytes"])
            and row["sha256"] == EXPECTED_HASHES[row["input_id"]]
            and source["sha256"] == row["sha256"]
            and sha256(path) == row["sha256"],
            f"stage input sealed identity drift: {row['input_id']}",
        )

    code_fields, code_rows = read_table(
        PREFLIGHT_ROOT / "manifests/producer_code_manifest.tsv"
    )
    require(
        code_fields == ["artifact_id", "path", "role", "sha256"],
        "stage producer-code manifest schema drift",
    )
    expected_code_order = [
        "producer",
        "atomic_publisher",
        "atomic_publisher_io",
        "Rscript",
        "Python",
    ]
    require(
        [row["artifact_id"] for row in code_rows] == expected_code_order,
        "stage producer-code ID/order drift",
    )
    runtime = load_json(
        BG001_ROOT / "contract/analysis_runtime_contract.json"
    )
    require(isinstance(runtime.get("Rscript"), dict), "stage runtime schema drift")
    expected_paths = {
        "producer": CANDIDATE_ROOT / "code/preflight_stage_substrate.R",
        "atomic_publisher": CANDIDATE_ROOT / "code/publish_noreplace.py",
        "atomic_publisher_io": CANDIDATE_ROOT / "code/snapshot_io.py",
    }
    expected_code_roles = {
        "producer": "resource_candidate_code",
        "atomic_publisher": "resource_candidate_code",
        "atomic_publisher_io": "resource_candidate_code",
        "Rscript": "analysis_runtime_rscript",
        "Python": "resource_candidate_python",
    }
    by_id = {row["artifact_id"]: row for row in code_rows}
    for artifact_id, path in expected_paths.items():
        path = path.resolve(strict=True)
        relative = path.relative_to(CANDIDATE_ROOT).as_posix()
        require(
            by_id[artifact_id]["path"] == str(path)
            and by_id[artifact_id]["role"] == expected_code_roles[artifact_id]
            and relative in source_by_candidate
            and source_by_candidate[relative]["role"] == expected_code_roles[artifact_id]
            and by_id[artifact_id]["sha256"] == source_by_candidate[relative]["sha256"]
            and sha256(path) == by_id[artifact_id]["sha256"],
            f"stage producer-code identity drift: {artifact_id}",
        )
    rscript = Path(str(runtime["Rscript"]["path"])).resolve(strict=True)
    require(
        by_id["Rscript"]["path"] == str(rscript)
        and by_id["Rscript"]["role"] == expected_code_roles["Rscript"]
        and by_id["Rscript"]["sha256"] == runtime["Rscript"]["sha256"]
        and sha256(rscript) == runtime["Rscript"]["sha256"],
        "stage Rscript identity drift",
    )
    require(
        by_id["Python"]["path"] == str(RNASEQ_PYTHON)
        and by_id["Python"]["role"] == expected_code_roles["Python"]
        and by_id["Python"]["sha256"] == RNASEQ_PYTHON_SHA256
        and sha256(RNASEQ_PYTHON) == RNASEQ_PYTHON_SHA256,
        "stage Python identity drift",
    )

    environment_fields, environment_rows = read_table(
        PREFLIGHT_ROOT / "manifests/environment_manifest.tsv"
    )
    require(
        environment_fields
        == [
            "component",
            "version",
            "hostname",
            "locale",
            "slurm_job_id",
            "slurm_partition",
            "slurm_qos",
            "slurm_nodelist",
            "slurm_cpus_per_task",
            "omp_num_threads",
            "openblas_num_threads",
            "mkl_num_threads",
            "veclib_maximum_threads",
            "numexpr_num_threads",
        ],
        "stage environment schema drift",
    )
    require(len(environment_rows) == 2, "stage environment row-count drift")
    package_versions = {
        row["name"]: row["version"]
        for row in runtime["R_runtime"]["packages"]
    }
    observed_versions = {
        row["component"]: row["version"]
        for row in environment_rows
    }
    require(
        observed_versions
        == {
            "R": package_versions["base"],
            "data.table": package_versions["data.table"],
        },
        "stage environment version drift",
    )
    hostnames = {row["hostname"] for row in environment_rows}
    locales = {row["locale"] for row in environment_rows}
    job_ids = {row["slurm_job_id"] for row in environment_rows}
    qos_values = {row["slurm_qos"] for row in environment_rows}
    nodelists = {row["slurm_nodelist"] for row in environment_rows}
    require_single_slurm_node(hostnames, nodelists, "stage")
    require(locales == {EXPECTED_EXECUTION_LOCALE}, "stage locale drift")
    require(len(job_ids) == 1 and next(iter(job_ids)).isdigit(), "stage job-ID drift")
    require(qos_values == {EXPECTED_EXECUTION_QOS}, "stage QOS drift")
    require(
        all(row["slurm_partition"] == EXPECTED_EXECUTION_PARTITION for row in environment_rows)
        and all(row["slurm_cpus_per_task"] == "4" for row in environment_rows)
        and all(
            row[field] == "4"
            for row in environment_rows
            for field in (
                "omp_num_threads",
                "openblas_num_threads",
                "mkl_num_threads",
                "veclib_maximum_threads",
                "numexpr_num_threads",
            )
        ),
        "stage SLURM environment drift",
    )

    status_fields, status_rows = read_table(
        PREFLIGHT_ROOT / "manifests/producer_status.tsv"
    )
    require(
        status_fields
        == [
            "status",
            "candidate_id",
            "workstream_id",
            "n_stage_samples",
            "n_transition_rows",
            "gene_filter_status",
            "model_fit_status",
            "manual_review_required",
        ]
        and len(status_rows) == 1,
        "stage producer-status schema drift",
    )
    status = status_rows[0]
    require(
        status
        == {
            "status": "PREFLIGHT_COMPLETE_PENDING_INDEPENDENT_VALIDATION",
            "candidate_id": CANDIDATE_ID,
            "workstream_id": "bulk-stage-preflight-v1",
            "n_stage_samples": "661",
            "n_transition_rows": "1129",
            "gene_filter_status": "not_run",
            "model_fit_status": "not_run",
            "manual_review_required": "TRUE",
        },
        "stage producer-status value drift",
    )


def verify_science_contract() -> dict[str, object]:
    require(PREFLIGHT_ROOT.is_dir() and not PREFLIGHT_ROOT.is_symlink(), "stage preflight root is missing")
    for path in PREFLIGHT_ROOT.rglob("*"):
        require(not path.is_symlink(), f"stage preflight contains symlink: {path}")
    require(not (PREFLIGHT_ROOT / "results").exists(), "no-fit preflight may not contain results")
    verify_producer_manifest()
    validate_input_and_code_manifests()

    input_fields, input_rows = read_table(PREFLIGHT_ROOT / "manifests/input_manifest.tsv")
    require(
        input_fields == ["input_id", "candidate_path", "size_bytes", "sha256", "role"],
        "input manifest schema drift",
    )
    require({row["input_id"] for row in input_rows} == set(EXPECTED_HASHES), "input set drift")
    for row in input_rows:
        relative = PurePosixPath(row["candidate_path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "unsafe input path")
        path = CANDIDATE_ROOT.joinpath(*relative.parts)
        require(path.is_file() and not path.is_symlink(), f"missing frozen input: {relative}")
        require(path.stat().st_size == int(row["size_bytes"]), f"input byte drift: {relative}")
        require(sha256(path) == row["sha256"] == EXPECTED_HASHES[row["input_id"]],
                f"input hash drift: {row['input_id']}")

    candidates, gse_expected, eligibility_expected = derive_source_stage_set()
    stage_fields, stage_rows = read_table(
        PREFLIGHT_ROOT / "audits/stage_sample_manifest.tsv"
    )
    require(
        stage_fields
        == [
            "sample_id",
            "analysis_unit_id",
            "analysis_unit_key_source",
            "dataset",
            "source_title",
            "biopsy",
            "fibrosis_stage",
            "inferred_sex",
            "pass_technical",
            "raw_present",
            "corrected_dge_present",
        ],
        "stage sample-manifest schema drift",
    )
    require(len(stage_rows) == 661, "emitted stage sample count drift")
    emitted = {row["sample_id"]: row for row in stage_rows}
    require(len(emitted) == 661 and set(emitted) == set(candidates), "emitted/source stage sample set mismatch")
    for sample, expected in candidates.items():
        row = emitted[sample]
        for field in (
            "analysis_unit_id",
            "analysis_unit_key_source",
            "dataset",
            "source_title",
            "biopsy",
        ):
            require(row[field] == str(expected[field]), f"stage field mismatch: {sample}/{field}")
        require(int(row["fibrosis_stage"]) == expected["fibrosis_stage"], f"stage mismatch: {sample}")
        require(parse_bool(row["pass_technical"]), f"nonpassing emitted sample: {sample}")
        require(parse_bool(row["raw_present"]), f"raw-absent emitted sample: {sample}")
        require(parse_bool(row["corrected_dge_present"]), f"DGE-absent emitted sample: {sample}")
        require(
            row["inferred_sex"] == expected["inferred_sex"],
            f"inferred-sex mismatch: {sample}",
        )
    require(len({row["analysis_unit_id"] for row in stage_rows}) == 661, "emitted duplicate analysis unit")

    observed_stage = Counter(int(row["fibrosis_stage"]) for row in stage_rows)
    require(dict(sorted(observed_stage.items())) == STAGE_COUNTS, "emitted stage census mismatch")
    stage_census_fields, stage_census = read_table(
        PREFLIGHT_ROOT / "audits/stage_census.tsv"
    )
    require(
        stage_census_fields == ["fibrosis_stage", "n_samples"]
        and len(stage_census) == 5,
        "stage census schema/cardinality drift",
    )
    require(
        {int(row["fibrosis_stage"]): int(row["n_samples"]) for row in stage_census} == STAGE_COUNTS,
        "stage census table mismatch",
    )

    cohort_stage_fields, cohort_stage_rows = read_table(
        PREFLIGHT_ROOT / "audits/cohort_stage_census.tsv"
    )
    require(
        cohort_stage_fields == ["dataset", "fibrosis_stage", "n_samples"]
        and len(cohort_stage_rows) == len(COHORTS) * 5,
        "cohort-stage census schema/cardinality drift",
    )
    cohort_stage_keys = [
        (row["dataset"], int(row["fibrosis_stage"]))
        for row in cohort_stage_rows
    ]
    require(
        cohort_stage_keys
        == [(cohort, stage) for cohort in COHORTS for stage in range(5)],
        "cohort-stage census key/order drift",
    )
    expected_cohort_stage = Counter(
        (str(row["dataset"]), int(row["fibrosis_stage"]))
        for row in candidates.values()
    )
    require(
        all(
            int(row["n_samples"])
            == expected_cohort_stage[(row["dataset"], int(row["fibrosis_stage"]))]
            for row in cohort_stage_rows
        ),
        "cohort-stage census values drift",
    )

    eligibility_fields, eligibility = read_table(
        PREFLIGHT_ROOT / "audits/stage_eligibility_audit.tsv"
    )
    expected_eligibility_fields = [
        "sample_id",
        "dataset",
        "condition",
        "group_binary",
        "sex",
        "age",
        "fibrosis_stage",
        "nas_score",
        "diagnosis_harmonized",
        "pass_technical",
        "inferred_sex",
        "matched_dataset",
        "matched_fibrosis_stage",
        "source_title",
        "biopsy",
        "analysis_unit_id",
        "is_first_biopsy",
        "analysis_unit_key_source",
        "raw_present",
        "corrected_dge_present",
        "qc_present",
        "meta_present",
        "included_primary",
        "exclusion_reason",
    ]
    require(
        eligibility_fields == expected_eligibility_fields
        and len(eligibility) == 731,
        "eligibility audit schema/cardinality drift",
    )
    eligibility_by_sample = {row["sample_id"]: row for row in eligibility}
    require(
        len(eligibility_by_sample) == 731
        and set(eligibility_by_sample) == set(eligibility_expected),
        "eligibility audit sample set drift",
    )
    for sample, expected in eligibility_expected.items():
        require(
            eligibility_by_sample[sample] == expected,
            f"eligibility audit row drift: {sample}",
        )

    crosswalk_fields, crosswalk = read_table(
        PREFLIGHT_ROOT / "audits/gse193066_participant_crosswalk.tsv"
    )
    expected_crosswalk_fields = [
        "run_id",
        "metadata_sample_id",
        "source_title",
        "biopsy",
        "source_fibrosis_stage",
        "participant_token",
        "analysis_unit_id",
        "is_first_biopsy",
        "is_second_biopsy",
        "harmonized_fibrosis_stage",
        "inferred_sex",
        "pass_technical",
        "raw_present",
        "corrected_dge_present",
        "included_primary",
        "exclusion_reason",
    ]
    require(
        crosswalk_fields == expected_crosswalk_fields
        and len(crosswalk) == 164,
        "GSE193066 crosswalk schema/cardinality drift",
    )
    crosswalk_by_sample = {row["run_id"]: row for row in crosswalk}
    require(
        len(crosswalk_by_sample) == 164
        and set(crosswalk_by_sample) == set(gse_expected),
        "GSE193066 crosswalk sample set drift",
    )
    for sample, source in gse_expected.items():
        expected = {
            "run_id": sample,
            "metadata_sample_id": sample,
            "source_title": str(source["title"]),
            "biopsy": str(source["biopsy"]),
            "source_fibrosis_stage": str(source["stage"]),
            "participant_token": str(source["participant"]),
            "analysis_unit_id": "GSE193066::{}".format(source["participant"]),
            "is_first_biopsy": (
                "TRUE" if source["biopsy"] == "1st biopsy" else "FALSE"
            ),
            "is_second_biopsy": (
                "TRUE" if source["biopsy"] == "2nd biopsy" else "FALSE"
            ),
            "harmonized_fibrosis_stage": str(source["stage"]),
            "inferred_sex": str(source["inferred_sex"]),
            "pass_technical": "TRUE" if source["pass_technical"] else "FALSE",
            "raw_present": "TRUE" if source["raw_present"] else "FALSE",
            "corrected_dge_present": (
                "TRUE" if source["corrected_dge_present"] else "FALSE"
            ),
            "included_primary": (
                "TRUE" if source["included_primary"] else "FALSE"
            ),
            "exclusion_reason": str(source["exclusion_reason"]),
        }
        require(
            crosswalk_by_sample[sample] == expected,
            f"GSE193066 crosswalk row drift: {sample}",
        )

    transition_fields, transition_rows = read_table(
        PREFLIGHT_ROOT / "audits/transition_sample_manifest.tsv"
    )
    require(
        transition_fields
        == [
            "sample_id",
            "analysis_unit_id",
            "dataset",
            "fibrosis_stage",
            "inferred_sex",
            "transition",
            "contrast",
            "arm",
        ],
        "transition sample-manifest schema drift",
    )
    require(len(transition_rows) == 1129, "transition sample manifest row-count drift")
    transition_lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in transition_rows:
        key = (row["transition"], row["sample_id"])
        require(key not in transition_lookup, f"duplicate transition sample: {key}")
        transition_lookup[key] = row

    summary_fields, transition_summary = read_table(
        PREFLIGHT_ROOT / "audits/transition_census.tsv"
    )
    require(
        summary_fields
        == [
            "transition",
            "contrast",
            "low_stage",
            "high_stage",
            "n_low",
            "n_high",
            "n_samples",
            "n_cohorts",
        ]
        and len(transition_summary) == 4,
        "transition census schema/cardinality drift",
    )
    summary_by_name = {row["transition"]: row for row in transition_summary}
    require(set(summary_by_name) == set(TRANSITIONS), "transition summary set drift")
    cohort_fields, cohort_rows = read_table(
        PREFLIGHT_ROOT / "audits/transition_cohort_census.tsv"
    )
    require(
        cohort_fields
        == [
            "transition",
            "contrast",
            "dataset",
            "n_low",
            "n_high",
            "n_total",
            "both_arms_present",
        ],
        "transition cohort-census schema drift",
    )
    expected_transition_cohort_rows = sum(
        len(spec["cohorts"]) for spec in TRANSITIONS.values()
    )
    transition_cohort_keys = [
        (row["transition"], row["dataset"]) for row in cohort_rows
    ]
    require(
        len(cohort_rows) == expected_transition_cohort_rows
        and len(set(transition_cohort_keys)) == expected_transition_cohort_rows,
        "transition cohort-census key/cardinality drift",
    )
    require(
        set(row["transition"] for row in cohort_rows) == set(TRANSITIONS)
        and all(
            row["contrast"] == TRANSITIONS[row["transition"]]["contrast"]
            for row in cohort_rows
        ),
        "transition cohort-census transition/contrast drift",
    )
    cohort_by_transition: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cohort_rows:
        cohort_by_transition[row["transition"]].append(row)

    expected_transition_keys: set[tuple[str, str]] = set()
    for transition, spec in TRANSITIONS.items():
        eligible_cohorts = set(spec["cohorts"])
        selected = {
            sample: row
            for sample, row in candidates.items()
            if row["dataset"] in eligible_cohorts
            and row["fibrosis_stage"] in {spec["low"], spec["high"]}
        }
        require(len(selected) == spec["n_samples"], f"source transition n drift: {transition}")
        require(sum(row["fibrosis_stage"] == spec["low"] for row in selected.values()) == spec["n_low"],
                f"source low-arm drift: {transition}")
        require(sum(row["fibrosis_stage"] == spec["high"] for row in selected.values()) == spec["n_high"],
                f"source high-arm drift: {transition}")
        expected_transition_keys.update((transition, sample) for sample in selected)
        for sample, expected_sample in selected.items():
            emitted_transition = transition_lookup[(transition, sample)]
            expected_arm = (
                "high"
                if expected_sample["fibrosis_stage"] == spec["high"]
                else "low"
            )
            require(
                emitted_transition["analysis_unit_id"]
                == expected_sample["analysis_unit_id"]
                and emitted_transition["dataset"] == expected_sample["dataset"]
                and int(emitted_transition["fibrosis_stage"])
                == expected_sample["fibrosis_stage"]
                and emitted_transition["inferred_sex"]
                == expected_sample["inferred_sex"]
                and emitted_transition["contrast"] == spec["contrast"]
                and emitted_transition["arm"] == expected_arm,
                f"transition sample field drift: {transition}/{sample}",
            )
        summary = summary_by_name[transition]
        require(summary["contrast"] == spec["contrast"], f"contrast drift: {transition}")
        require(
            int(summary["low_stage"]) == spec["low"]
            and int(summary["high_stage"]) == spec["high"],
            f"transition stage-label drift: {transition}",
        )
        require(int(summary["n_low"]) == spec["n_low"], f"summary low drift: {transition}")
        require(int(summary["n_high"]) == spec["n_high"], f"summary high drift: {transition}")
        require(int(summary["n_samples"]) == spec["n_samples"], f"summary n drift: {transition}")
        require(int(summary["n_cohorts"]) == len(spec["cohorts"]), f"summary cohort drift: {transition}")
        observed_cohorts = {row["dataset"] for row in cohort_by_transition[transition]}
        require(observed_cohorts == eligible_cohorts, f"cohort set drift: {transition}")
        require(all(parse_bool(row["both_arms_present"]) for row in cohort_by_transition[transition]),
                f"single-arm cohort entered: {transition}")
        for row in cohort_by_transition[transition]:
            expected_subset = [
                sample
                for sample in selected.values()
                if sample["dataset"] == row["dataset"]
            ]
            expected_low = sum(
                sample["fibrosis_stage"] == spec["low"]
                for sample in expected_subset
            )
            expected_high = sum(
                sample["fibrosis_stage"] == spec["high"]
                for sample in expected_subset
            )
            require(
                int(row["n_low"]) == expected_low
                and int(row["n_high"]) == expected_high
                and int(row["n_total"]) == len(expected_subset),
                f"transition cohort count drift: {transition}/{row['dataset']}",
            )
    require(set(transition_lookup) == expected_transition_keys, "transition sample set mismatch")

    design_fields, design_rows = read_table(
        PREFLIGHT_ROOT / "audits/model_design_preflight.tsv"
    )
    require(
        design_fields
        == [
            "sample_id",
            "transition",
            "contrast",
            "(Intercept)",
            "datasetGSE135251",
            "datasetGSE162694",
            "datasetGSE174478",
            "datasetGSE193066",
            "datasetGSE240729",
            "inferred_sexM",
            "fib_grouphigh",
        ],
        "model-design preflight schema drift",
    )
    require(
        len(design_rows) == 1129
        and set(row["transition"] for row in design_rows) == set(TRANSITIONS),
        "model-design row/transition drift",
    )
    audit_fields, design_audit = read_table(
        PREFLIGHT_ROOT / "audits/design_preflight.tsv"
    )
    require(
        audit_fields
        == [
            "transition",
            "contrast",
            "formula",
            "coefficient",
            "design_rank",
            "n_design_columns",
            "design_columns",
            "n_samples",
            "n_cohorts",
            "cohort_set",
            "gene_filter_status",
            "model_fit_status",
        ]
        and len(design_audit) == 4,
        "design-audit schema/cardinality drift",
    )
    audit_by_transition = {row["transition"]: row for row in design_audit}
    require(set(audit_by_transition) == set(TRANSITIONS), "design audit set drift")
    design_by_transition: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in design_rows:
        design_by_transition[row["transition"]].append(row)
    for transition, spec in TRANSITIONS.items():
        rows = design_by_transition[transition]
        require(len(rows) == spec["n_samples"], f"design row count drift: {transition}")
        cohort_levels = sorted(spec["cohorts"])
        expected_columns = ["(Intercept)"] + [f"dataset{item}" for item in cohort_levels[1:]]
        expected_columns += ["inferred_sexM", "fib_grouphigh"]
        require(all(column in design_fields for column in expected_columns), f"design column absent: {transition}")
        matrix: list[list[float]] = []
        seen: set[str] = set()
        for row in rows:
            sample = row["sample_id"]
            require(sample not in seen, f"duplicate design sample: {transition}/{sample}")
            require(
                row["contrast"] == spec["contrast"],
                f"design contrast drift: {transition}/{sample}",
            )
            seen.add(sample)
            source = transition_lookup[(transition, sample)]
            expected_values = [1.0]
            expected_values += [1.0 if source["dataset"] == item else 0.0 for item in cohort_levels[1:]]
            expected_values += [1.0 if source["inferred_sex"] == "M" else 0.0]
            expected_values += [1.0 if source["arm"] == "high" else 0.0]
            observed_values = [float(row[column]) for column in expected_columns]
            require(observed_values == expected_values, f"design value drift: {transition}/{sample}")
            nonapplicable = {
                "datasetGSE135251",
                "datasetGSE162694",
                "datasetGSE174478",
                "datasetGSE193066",
                "datasetGSE240729",
            } - set(expected_columns)
            require(
                all(row[column] == "NA" for column in nonapplicable),
                f"non-applicable design column is populated: {transition}/{sample}",
            )
            matrix.append(observed_values)
        rank = matrix_rank(matrix)
        require(rank == spec["rank"] == len(expected_columns), f"design rank drift: {transition}")
        audit = audit_by_transition[transition]
        require(audit["contrast"] == spec["contrast"], f"audit contrast drift: {transition}")
        require(audit["formula"] == "~ dataset + inferred_sex + fib_group", f"formula drift: {transition}")
        require(audit["coefficient"] == "fib_grouphigh", f"coefficient drift: {transition}")
        require(int(audit["design_rank"]) == rank, f"audit rank drift: {transition}")
        require(int(audit["n_design_columns"]) == len(expected_columns), f"audit column drift: {transition}")
        require(audit["design_columns"].split(";") == expected_columns, f"audit design names drift: {transition}")
        require(int(audit["n_samples"]) == spec["n_samples"], f"audit sample-count drift: {transition}")
        require(int(audit["n_cohorts"]) == len(spec["cohorts"]), f"audit cohort-count drift: {transition}")
        require(
            audit["cohort_set"].split(";") == sorted(spec["cohorts"]),
            f"audit cohort-set drift: {transition}",
        )
        require(audit["gene_filter_status"] == "not_run", f"filter ran in preflight: {transition}")
        require(audit["model_fit_status"] == "not_run", f"model ran in preflight: {transition}")

    contract_fields, contract_rows = read_table(PREFLIGHT_ROOT / "manifests/run_contract.tsv")
    require(contract_fields == ["key", "value"], "run contract schema drift")
    contract = {row["key"]: row["value"] for row in contract_rows}
    require(
        len(contract_rows) == 16
        and len(contract) == 16,
        "run contract key/cardinality drift",
    )
    expected_contract = {
        "candidate_id": CANDIDATE_ID,
        "workstream_id": "bulk-stage-preflight-v1",
        "source_layer": "human_bulk_corrected_ninecohort_candidate",
        "biological_unit": (
            "one deposited biopsy per analysis unit; GSE193066 restricted "
            "to first biopsy per title-root participant"
        ),
        "allowlist": ";".join(COHORTS),
        "gse193066_timepoint": "1st biopsy only",
        "gse193066_participant_key": (
            "deposited !Sample_title with terminal _1/_2 removed"
        ),
        "n_stage_samples": "661",
        "stage_counts": "F0=126;F1=187;F2=174;F3=132;F4=42",
        "gene_filter_status": "not_run",
        "normalization_status": "not_run",
        "model_fit_status": "not_run",
        "hypothesis_test_status": "not_run",
        "canonical_write": "false",
        "figure_write": "false",
        "manual_review_required": "true",
    }
    require(contract == expected_contract, "run contract value drift")

    return {
        "status": "STAGE_PREFLIGHT_VALIDATED",
        "candidate_id": CANDIDATE_ID,
        "workstream_id": "bulk-stage-preflight-v1",
        "n_stage_samples": 661,
        "stage_counts": STAGE_COUNTS,
        "transition_counts": {
            transition: {
                "n_low": spec["n_low"],
                "n_high": spec["n_high"],
                "n_samples": spec["n_samples"],
                "n_cohorts": len(spec["cohorts"]),
                "design_rank": spec["rank"],
            }
            for transition, spec in TRANSITIONS.items()
        },
        "gene_filter_status": "not_run",
        "model_fit_status": "not_run",
        "manual_review_required": True,
        "stage_model_authorized": False,
        "canonical_promotion_authorized": False,
    }


VALIDATION_CHECKS = [
    {"check_id": "source_identity", "status": "PASS"},
    {"check_id": "code_runtime_identity", "status": "PASS"},
    {"check_id": "gse193066_participant_gate", "status": "PASS"},
    {"check_id": "stage_661_census", "status": "PASS"},
    {"check_id": "transition_censuses", "status": "PASS"},
    {"check_id": "no_fit_design_rank", "status": "PASS"},
    {"check_id": "no_model_or_canonical_write", "status": "PASS"},
]


def normalized_verdict(verdict: dict[str, object]) -> dict[str, object]:
    return json.loads(json.dumps(verdict, sort_keys=True))


def build_validation_manifest(
    validation_files: list[Path],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(
        path
        for path in PREFLIGHT_ROOT.rglob("*")
        if path.is_file()
        and path != READY_PATH
        and VALIDATION_ROOT not in path.parents
    ):
        rows.append(
            {
                "candidate_path": path.relative_to(PREFLIGHT_ROOT).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    for path in validation_files:
        rows.append(
            {
                "candidate_path": f"validation/{path.name}",
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    rows.sort(key=lambda row: str(row["candidate_path"]))
    return rows


def publish_validation(verdict: dict[str, object]) -> None:
    require(
        not VALIDATION_ROOT.exists() and not VALIDATION_ROOT.is_symlink(),
        "stage validation bundle already exists",
    )
    require(
        not READY_PATH.exists() and not READY_PATH.is_symlink(),
        "stage READY marker already exists",
    )
    temporary = (
        PREFLIGHT_ROOT.parent
        / f".{PREFLIGHT_ROOT.name}.validation.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    )
    require(
        not temporary.exists() and not temporary.is_symlink(),
        "stage validation temporary collision",
    )
    temporary.mkdir(mode=0o750)
    try:
        checks_path = temporary / "validation_checks.tsv"
        report_path = temporary / "validation_report.json"
        manifest_path = temporary / "validated_artifact_manifest.tsv"
        write_exclusive(
            checks_path,
            tsv_payload(VALIDATION_CHECKS, ["check_id", "status"]),
        )
        report = normalized_verdict(verdict)
        report["validated_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_exclusive(
            report_path,
            json.dumps(report, indent=2, sort_keys=True) + "\n",
        )
        manifest_rows = build_validation_manifest(
            [checks_path, report_path]
        )
        write_exclusive(
            manifest_path,
            tsv_payload(
                manifest_rows,
                ["candidate_path", "size_bytes", "sha256"],
            ),
        )
        verify_candidate_snapshot()
        publication = publish_directory_noreplace(
            temporary,
            VALIDATION_ROOT,
            marker="validated_artifact_manifest.tsv",
        )
        require(
            publication.get("source_retained") is True,
            "stage validation publisher did not report retained source",
        )
        require(
            temporary.is_dir() and not temporary.is_symlink(),
            "stage validation publication source was not retained safely",
        )
        require_regular_file(manifest_path)
        require_regular_file(VALIDATION_MANIFEST)
        require(
            sha256(manifest_path) == sha256(VALIDATION_MANIFEST),
            "published stage-validation commit-marker hash drift",
        )
    except Exception:
        print(
            f"FAILED stage validation temporary retained for audit: {temporary}",
            file=os.sys.stderr,
        )
        raise
    finalize_ready_marker(verdict)


def verify_validation_bundle(
    verdict: dict[str, object],
) -> dict[str, object]:
    require(
        PREFLIGHT_ROOT.is_dir() and not PREFLIGHT_ROOT.is_symlink(),
        "stage preflight root is invalid",
    )
    require(
        VALIDATION_ROOT.is_dir() and not VALIDATION_ROOT.is_symlink(),
        "stage validation bundle is missing or symlinked",
    )
    for path in PREFLIGHT_ROOT.rglob("*"):
        require(not path.is_symlink(), f"symlink in stage validation tree: {path}")
    fields, rows = read_table(VALIDATION_MANIFEST)
    require(
        fields == ["candidate_path", "size_bytes", "sha256"]
        and rows,
        "stage validated-manifest schema drift",
    )
    observed: set[str] = set()
    for row in rows:
        relative = PurePosixPath(row["candidate_path"])
        require(
            not relative.is_absolute() and ".." not in relative.parts,
            "unsafe stage validated path",
        )
        text = relative.as_posix()
        require(text not in observed, f"duplicate stage validated path: {text}")
        observed.add(text)
        path = PREFLIGHT_ROOT.joinpath(*relative.parts)
        require_regular_file(path)
        require(
            path.stat().st_size == int(row["size_bytes"]),
            f"stage validated byte drift: {text}",
        )
        require(
            sha256(path) == row["sha256"],
            f"stage validated hash drift: {text}",
        )
    actual = {
        path.relative_to(PREFLIGHT_ROOT).as_posix()
        for path in PREFLIGHT_ROOT.rglob("*")
        if path.is_file() and path not in {VALIDATION_MANIFEST, READY_PATH}
    }
    require(observed == actual, "stage validated-manifest file-set mismatch")

    check_fields, checks = read_table(
        VALIDATION_ROOT / "validation_checks.tsv"
    )
    require(
        check_fields == ["check_id", "status"]
        and checks == VALIDATION_CHECKS,
        "stage validation check family drift",
    )
    report = load_json(VALIDATION_ROOT / "validation_report.json")
    expected = normalized_verdict(verdict)
    require(
        set(report) == set(expected) | {"validated_at_utc"},
        "stage validation-report schema drift",
    )
    for key, value in expected.items():
        require(report.get(key) == value, f"stage validation-report drift: {key}")
    require(
        isinstance(report.get("validated_at_utc"), str)
        and bool(report["validated_at_utc"]),
        "stage validation timestamp missing",
    )
    return report


def ready_payload(verdict: dict[str, object]) -> dict[str, object]:
    payload = normalized_verdict(verdict)
    _, rows = read_table(VALIDATION_MANIFEST)
    payload.update(
        {
            "validation_report_sha256": sha256(
                VALIDATION_ROOT / "validation_report.json"
            ),
            "validation_checks_sha256": sha256(
                VALIDATION_ROOT / "validation_checks.tsv"
            ),
            "validated_artifact_manifest_sha256": sha256(
                VALIDATION_MANIFEST
            ),
            "validated_file_count": len(rows),
        }
    )
    return payload


def seal_read_only() -> None:
    require(
        PREFLIGHT_ROOT.is_dir() and not PREFLIGHT_ROOT.is_symlink(),
        "stage root invalid before read-only seal",
    )
    descendants = list(PREFLIGHT_ROOT.rglob("*"))
    for path in descendants:
        require(not path.is_symlink(), f"symlink before stage read-only seal: {path}")
    for path in descendants:
        if path.is_file():
            chmod_nofollow(path, 0o440, directory=False)
    for path in sorted(
        (path for path in descendants if path.is_dir()),
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        chmod_nofollow(path, 0o550, directory=True)
    chmod_nofollow(PREFLIGHT_ROOT, 0o550, directory=True)


def finalize_ready_marker(verdict: dict[str, object]) -> None:
    verify_validation_bundle(verdict)
    verify_candidate_snapshot()
    expected = ready_payload(verdict)
    require(not READY_PATH.is_symlink(), "stage READY marker may not be a symlink")
    if READY_PATH.exists():
        require(
            load_json(READY_PATH) == expected,
            "existing stage READY marker is inconsistent",
        )
    else:
        temporary = (
            PREFLIGHT_ROOT.parent
            / f".{PREFLIGHT_ROOT.name}.READY.tmp.{os.getpid()}.{uuid.uuid4().hex}"
        )
        require(
            not temporary.exists() and not temporary.is_symlink(),
            "stage READY temporary collision",
        )
        write_exclusive(
            temporary,
            json.dumps(expected, indent=2, sort_keys=True) + "\n",
        )
        verify_validation_bundle(verdict)
        verify_candidate_snapshot()
        publication = publish_file_noreplace(temporary, READY_PATH)
        require(
            publication.get("source_retained") is True,
            "stage READY publisher did not report retained source",
        )
        require_regular_file(temporary)
        require_regular_file(READY_PATH)
        require(
            sha256(temporary) == sha256(READY_PATH),
            "published stage READY marker hash drift",
        )
    verify_validation_bundle(verdict)
    seal_read_only()


def verify_existing_seal(verdict: dict[str, object]) -> dict[str, object]:
    verify_candidate_snapshot()
    verify_producer_manifest()
    report = verify_validation_bundle(verdict)
    require(
        load_json(READY_PATH) == ready_payload(verdict),
        "stage READY differs from validation bundle",
    )
    for path in PREFLIGHT_ROOT.rglob("*"):
        require(not path.is_symlink(), f"symlink in sealed stage tree: {path}")
        expected_mode = 0o440 if path.is_file() else 0o550
        require(
            stat.S_IMODE(path.stat().st_mode) == expected_mode,
            f"sealed stage mode drift: {path}",
        )
    require(
        stat.S_IMODE(PREFLIGHT_ROOT.stat().st_mode) == 0o550,
        "sealed stage root mode drift",
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    verify_candidate_snapshot()
    verdict = verify_science_contract()
    if args.verify_only:
        print(
            json.dumps(
                verify_existing_seal(verdict),
                indent=2,
                sort_keys=True,
            )
        )
        return
    if READY_PATH.exists() or READY_PATH.is_symlink():
        require(
            READY_PATH.is_file() and not READY_PATH.is_symlink(),
            "invalid stage READY marker",
        )
        require(
            VALIDATION_ROOT.is_dir() and not VALIDATION_ROOT.is_symlink(),
            "stage validation bundle missing",
        )
        finalize_ready_marker(verdict)
    elif VALIDATION_ROOT.exists() or VALIDATION_ROOT.is_symlink():
        require(
            VALIDATION_ROOT.is_dir() and not VALIDATION_ROOT.is_symlink(),
            "invalid stage validation bundle",
        )
        finalize_ready_marker(verdict)
    else:
        publish_validation(verdict)
    print(
        json.dumps(
            verify_existing_seal(verdict),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except (
        ContractError,
        SnapshotIOError,
        KeyError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        raise SystemExit(f"STAGE PREFLIGHT VALIDATION ERROR: {error}") from error
