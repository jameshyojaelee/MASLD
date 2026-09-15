#!/usr/bin/env python3
"""Shared read-only contract for the F_five Resource-paper candidate."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path


PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
).resolve()

CANDIDATE_ID = "resource-f-five-coloc-v6-candidate-2026-08-10"
CANDIDATE_ROOT = (
    PROJECT_ROOT
    / "RNA-seq/results/manuscript_release/candidates"
    / CANDIDATE_ID
)

EXPECTED_EXECUTION_PARTITION = "io"
EXPECTED_EXECUTION_QOS = "nslab"
EXPECTED_EXECUTION_LOCALE = (
    "LC_CTYPE=en_US.UTF-8;LC_NUMERIC=C;LC_TIME=en_US.UTF-8;"
    "LC_COLLATE=en_US.UTF-8;LC_MONETARY=en_US.UTF-8;"
    "LC_MESSAGES=en_US.UTF-8;LC_PAPER=en_US.UTF-8;LC_NAME=C;"
    "LC_ADDRESS=C;LC_TELEPHONE=C;LC_MEASUREMENT=en_US.UTF-8;LC_IDENTIFICATION=C"
)

BG_RUN_ID = "bg001-fragment-v211-gencode49-20260807T195243Z"
BG_ROOT = PROJECT_ROOT / "results/remediation/bg001" / BG_RUN_ID
F_FIVE_ROOT = BG_ROOT / "arms/F_five"
F_LEGACY_ROOT = BG_ROOT / "arms/F_legacy"
CANONICAL_REFERENCE = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
    / "canonical_deg_results.csv"
)
CANONICAL_REFERENCE_SHA256 = "e27fd319033fcd290ae2d6d8eed5b41500e8eb06460eafb4594b7cb20abf7bcd"
GSE193066_METADATA = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata/metadata.tsv"
)
GSE193066_METADATA_SHA256 = "b66158d20da2127e8c3cc7ebae4f1d85a2d6813b7dd6f2af4bdac60070ff2bc5"
UNIFIED_METADATA = (
    BG_ROOT
    / "source_snapshot/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
)
UNIFIED_METADATA_SHA256 = "b9f6afc1e916391951c638a60d8ec95f1b380d20311469c0830f638f16649699"

BASELINE_FROZEN_SHA256 = "4ab863fce3c8129573ed090c0f5c47daa71325a455a3c0123ddc0e5f097bca33"
BASELINE_COMPATIBILITY_SHA256 = "b152a2922d71af7e41f89c2fe172a5513a5269603db92cb18bba89734d359623"
BASELINE_FILES_MANIFEST_SHA256 = "e598ef77c17c2aed8b4f9106d33ec5c83b9ca07b8950e507af293d12a429610a"
BASELINE_MANIFESTS_MANIFEST_SHA256 = "e17cc0e6db3355366105af5fdef6f97f432973b329360618874dc76c39659521"
RECOUNT_ARTIFACTS_MANIFEST_SHA256 = "4932bf96e3d15e2e52ed6e57ffa5c99591f3d1235a60a77a2504436ca6909458"
RECOUNT_COMPLETE_SHA256 = "208e03c6a34785041928a3959babd7ab32644d5cf51d114d19f210518635c58c"
STRUCTURAL_VALIDATION_SHA256 = "9cb45cc751261b00d7c27873734b421680098d33ee8483a0efe3ef99cc41b689"

ACCEPTED_RUN_INPUT = CANDIDATE_ROOT / "inputs/BG001-DECISION"
F_FIVE_INPUT = ACCEPTED_RUN_INPUT / "arms/F_five"
F_LEGACY_INPUT = ACCEPTED_RUN_INPUT / "arms/F_legacy"
POOLED_REPRO_ROOT = CANDIDATE_ROOT / "workstreams/BULK-POOLED-REPRO"
STAGE_ROOT = CANDIDATE_ROOT / "workstreams/BULK-STAGE"

EXPECTED_HASHES = {
    "ANALYSIS_COMPLETE.json": "f8be659440579b510f2945cd43942b6fd617acf7ea34a666561b6e5626dcd389",
    "arms/F_five/ARM_COMPLETE.json": "2e33912cdafc9fa05a692481f2fb02cfe2481573d85eebf8799006ac030a7ff4",
    "arms/F_five/provenance/artifact_manifest.tsv": "b6eeafd82691617c9ea924c0ebf42631d10429592dead6265b47211a81516503",
    "arms/F_five/results/integration/deg_results.csv": "691b09517f42ff24e472fdb9db460b1c339ec734876c40115f8b146293b5b633",
    "arms/F_five/results/integration/merged_counts_raw.rds": "16afc00bf1db706731d07df1dfeabdad0bc423b1caa1b822631685a55ec7f225",
    "arms/F_five/results/integration/merged_dge.rds": "bc3ea19c8b064861339eba862abf58e82ddb5460c622b550dc45aa3b5fcd1f76",
    "arms/F_five/results/integration/meta_matched.rds": "90ba6aca68643c4680c08fc982a6c72b76edb740ab726874cb50603711756a83",
    "arms/F_legacy/ARM_COMPLETE.json": "3b42688402a1db1b2d136286aea705d6c91f98bab1cb06fed5e2fe3e9b46cf50",
    "arms/F_legacy/provenance/artifact_manifest.tsv": "cc9dcf3915945d7f6b6ab4d7f8b616afbd50b3c1a9a2bcbe7a4754d8d4b3f1fe",
    "arms/F_legacy/results/integration/merged_dge.rds": "56d1cf97f791a81404ddc6c5ea4c5df9e6a38e7cf011b39d0d7c3279c87afa4d",
    "arms/F_legacy/qc/sample_qc_report.csv": "1c7fef991c599554369ecc79de2b2f441597cf9b6c1a7ff3ecbe05f1c5969560",
    "comparisons/artifact_manifest.tsv": "f1a43fb7799f27714835e3d5ece8c86b556b3e2fc315e9987a6d89a40bcaadb7",
    "comparisons/verdict.json": "6799d9ddf4d0b3ab3780ecce4248d04bd0ba3e63f3fdf2c6241b3d3e503cf80a",
    "comparisons/dependency_manifest.tsv": "fc9e08e24efadfd1efedfc6f1f2ec3eb992cc65d1c4d8c92216e5164bab7d919",
    "comparisons/escalation_disposition.tsv": "6af7e1ad481012d56d6c2856895b55674ca7955e7acf1e341b83a197e3614a11",
    "promotion/CANDIDATE_ACCEPTANCE.md": "8bbd22663304eeef05fbb79ea666690b0a6a9f5ff4ebb3c7207a7fd8366f1616",
}

EXPECTED_POOLED = {
    "n_samples": 844,
    "n_control": 157,
    "n_disease": 687,
    "n_cohorts": 5,
    "n_genes": 23370,
    "n_treat": 1616,
    "n_up": 1144,
    "n_down": 472,
}

EXPECTED_RESOURCE = {
    "n_matched": 1281,
    "n_pass_technical": 1257,
    "n_raw_genes": 86369,
    "n_f_legacy_genes": 24196,
    "n_cohorts": 9,
}

RNASEQ_PYTHON = Path(
    "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python3.10"
)
RNASEQ_PYTHON_SHA256 = (
    "2e76c2d5e89880deb45cf75abfb81a10029f1d9cbbb155686e036720f72f7502"
)
GENCODE_GTF = Path(
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
    "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)
GENCODE_GTF_SHA256 = (
    "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
)
GENCODE_GTF_BYTES = 97_570_076

REVIEWED_CODE_FILES = (
    "README.md",
    "accepted_run_provenance.py",
    "preflight_bulk_census.R",
    "preflight_stage_substrate.R",
    "prepare_candidate.py",
    "prepare_candidate.sbatch",
    "publish_noreplace.py",
    "reproduce_pooled_fivecohort.R",
    "resource_contract.py",
    "run_pooled_reproduction.sbatch",
    "run_stage_preflight.sbatch",
    "snapshot_io.py",
    "validate_pooled_candidate.py",
    "validate_stage_preflight.py",
    "validate_stage_rds_membership.R",
    "verify_candidate_snapshot.py",
    "verify_mirrored_provenance.py",
)


class ContractError(RuntimeError):
    """Raised when an read-only candidate contract is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _is_dns_label(value: str) -> bool:
    return (
        1 <= len(value) <= 63
        and value[0].isascii()
        and value[0].isalnum()
        and value[-1].isascii()
        and value[-1].isalnum()
        and all(
            character.isascii()
            and (character.isalnum() or character == "-")
            for character in value
        )
    )


def require_single_slurm_node(
    hostnames: set[str],
    nodelists: set[str],
    label: str,
) -> tuple[str, str]:
    require(
        len(hostnames) == 1 and len(nodelists) == 1,
        f"{label} SLURM node cardinality drift",
    )
    hostname = next(iter(hostnames))
    nodelist = next(iter(nodelists))
    require(
        hostname not in {"", "not_slurm"} and nodelist not in {"", "not_slurm"},
        f"{label} SLURM node identity missing",
    )
    normalized_hostname = hostname.removesuffix(".")
    hostname_labels = normalized_hostname.split(".")
    require(
        normalized_hostname
        and len(normalized_hostname) <= 253
        and all(_is_dns_label(part) for part in hostname_labels),
        f"{label} hostname format drift",
    )
    require(
        _is_dns_label(nodelist)
        and nodelist.lower() == hostname_labels[0].lower(),
        f"{label} SLURM node/hostname drift",
    )
    return hostname, nodelist


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def relative_to_project(path: Path) -> str:
    resolved = path.resolve(strict=True)
    try:
        return resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError as error:
        raise ContractError(f"path escapes project root: {resolved}") from error


def require_regular_file(path: Path) -> None:
    require(path.is_file(), f"required file is missing: {path}")
    require(not path.is_symlink(), f"symlinked source is prohibited: {path}")


def verify_manifest(root: Path, manifest_path: Path) -> list[dict[str, str]]:
    require_regular_file(manifest_path)
    rows = read_tsv(manifest_path)
    require(rows, f"empty manifest: {manifest_path}")
    require(
        set(rows[0]) == {"relative_path", "size_bytes", "sha256"},
        f"unexpected artifact-manifest schema: {manifest_path}",
    )
    observed: set[str] = set()
    for row in rows:
        relative = row["relative_path"]
        require(relative not in observed, f"duplicate manifest path: {relative}")
        observed.add(relative)
        source = root / relative
        require_regular_file(source)
        require(source.stat().st_size == int(row["size_bytes"]), f"byte drift: {source}")
        require(sha256(source) == row["sha256"], f"hash drift: {source}")
    return rows
