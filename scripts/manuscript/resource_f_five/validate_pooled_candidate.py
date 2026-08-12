#!/usr/bin/env python3
"""Validate, atomically seal, and reverify the F_five pooled reproduction."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import stat
import uuid
from datetime import datetime, timezone
from pathlib import Path

from resource_contract import (
    CANONICAL_REFERENCE,
    CANONICAL_REFERENCE_SHA256,
    CANDIDATE_ID,
    CANDIDATE_ROOT,
    EXPECTED_POOLED,
    EXPECTED_EXECUTION_LOCALE,
    EXPECTED_EXECUTION_PARTITION,
    EXPECTED_EXECUTION_QOS,
    F_FIVE_INPUT,
    POOLED_REPRO_ROOT,
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


NONPROBABILITY_TOLERANCES = {
    "logFC": 2e-12,
    "SE": 2e-12,
    "t": 2e-10,
    "shrunk_logFC": 2e-10,
    "treat_lfc": 0.0,
    "AveExpr": 2e-12,
}
PROBABILITY_FIELDS = ("P.Value", "padj", "lfsr", "treat_p", "treat_fdr")
PROBABILITY_ABSOLUTE_TOLERANCE = 1e-12
PROBABILITY_RELATIVE_TOLERANCE = 1e-8
VALIDATION_ROOT = POOLED_REPRO_ROOT / "validation"
READY_PATH = POOLED_REPRO_ROOT / "VALIDATED.json"
VALIDATION_MANIFEST = VALIDATION_ROOT / "validated_artifact_manifest.tsv"
DEG_FIELDS = [
    "gene",
    "logFC",
    "SE",
    "t",
    "P.Value",
    "padj",
    "shrunk_logFC",
    "lfsr",
    "treat_lfc",
    "treat_p",
    "treat_fdr",
    "AveExpr",
    "symbol",
]
EXPECTED_CHECK_IDS = {
    "P01_PRODUCER_MANIFEST",
    "P02_FULL_NUMERICAL_REPRODUCTION",
    "P03_DESIGN_IDENTITY",
    "P04_WEIGHT_IDENTITY",
    "P05_RUNTIME_AND_INPUT_IDENTITY",
    "P06_NONCANONICAL_CANARY",
}


def write_exclusive(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def tsv_payload(rows: list[dict[str, object]], fields: list[str]) -> str:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def read_tsv_strict(
    path: Path, expected_fields: list[str] | None = None
) -> tuple[list[str], list[dict[str, str]]]:
    require_regular_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if expected_fields is not None:
        require(fields == expected_fields, f"TSV schema drift: {path}")
    require(
        all(None not in row and all(value is not None for value in row.values()) for row in rows),
        f"malformed TSV row width: {path}",
    )
    return fields, rows


def load_json(path: Path) -> dict[str, object]:
    require_regular_file(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def require_live_canonical_unchanged() -> None:
    require_regular_file(CANONICAL_REFERENCE)
    require(
        sha256(CANONICAL_REFERENCE) == CANONICAL_REFERENCE_SHA256,
        "live canonical DEG table changed during candidate execution",
    )


def validate_producer_manifest() -> None:
    manifest_path = POOLED_REPRO_ROOT / "producer_artifact_manifest.tsv"
    _, rows = read_tsv_strict(manifest_path, ["artifact_path", "bytes", "sha256"])
    require(rows, "empty pooled producer manifest")
    paths: set[str] = set()
    for row in rows:
        relative = row["artifact_path"]
        require(relative not in paths, f"duplicate producer artifact: {relative}")
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts, "unsafe producer path")
        paths.add(relative)
        path = POOLED_REPRO_ROOT / relative
        require_regular_file(path)
        require(path.stat().st_size == int(row["bytes"]), f"producer byte drift: {relative}")
        require(sha256(path) == row["sha256"], f"producer hash drift: {relative}")
    observed = {
        path.relative_to(POOLED_REPRO_ROOT).as_posix()
        for path in POOLED_REPRO_ROOT.rglob("*")
        if path.is_file()
        and path
        not in {
            manifest_path,
            POOLED_REPRO_ROOT / "PRODUCER_COMPLETE.json",
            READY_PATH,
        }
        and VALIDATION_ROOT not in path.parents
    }
    require(paths == observed, "producer manifest file-set mismatch")


def read_gene_table(
    path: Path,
) -> tuple[list[str], list[str], dict[str, dict[str, str]]]:
    require_regular_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        require(fields == DEG_FIELDS, f"exact DEG schema drift: {path}")
        order: list[str] = []
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            require(
                None not in row and all(value is not None for value in row.values()),
                f"malformed DEG row width: {path}",
            )
            gene = row["gene"]
            require(gene and gene not in rows, f"blank or duplicate gene in {path}: {gene!r}")
            order.append(gene)
            rows[gene] = row
    return fields, order, rows


def compare_probability(
    left: float, right: float, field: str, gene: str
) -> tuple[float, float, float | None, bool]:
    require(math.isfinite(left) and math.isfinite(right), f"nonfinite {field}: {gene}")
    require(0 <= left <= 1 and 0 <= right <= 1, f"out-of-range {field}: {gene}")
    absolute_delta = abs(left - right)
    scale = max(abs(left), abs(right))
    relative_delta = absolute_delta / scale if scale else 0.0
    limit = max(
        PROBABILITY_ABSOLUTE_TOLERANCE,
        PROBABILITY_RELATIVE_TOLERANCE * scale,
    )
    require(
        absolute_delta <= limit,
        (
            f"scale-aware {field} reproduction failure for {gene}: "
            f"absolute={absolute_delta}, relative={relative_delta}, limit={limit}"
        ),
    )
    log10_delta = (
        abs(math.log10(left) - math.log10(right))
        if left > 0 and right > 0
        else None
    )
    return absolute_delta, relative_delta, log10_delta, (left == 0) != (right == 0)


def validate_environment() -> dict[str, object]:
    contract = load_json(
        CANDIDATE_ROOT
        / "inputs/BG001-DECISION/contract/analysis_runtime_contract.json"
    )
    package_versions = {
        row["name"]: row["version"] for row in contract["R_runtime"]["packages"]
    }
    _, rows = read_tsv_strict(
        POOLED_REPRO_ROOT / "environment_manifest.tsv",
        [
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
    )
    require(len(rows) == 7, "environment-manifest cardinality drift")
    require(len({row["component"] for row in rows}) == 7, "duplicate environment component")
    observed = {row["component"]: row["version"] for row in rows}
    expected = {
        "R": package_versions["base"],
        **{
            name: package_versions[name]
            for name in ("data.table", "edgeR", "limma", "ashr", "yaml", "jsonlite")
        },
    }
    require(observed == expected, "recorded R environment differs from frozen contract")
    hostnames = {row["hostname"] for row in rows}
    locales = {row["locale"] for row in rows}
    job_ids = {row["slurm_job_id"] for row in rows}
    qos_values = {row["slurm_qos"] for row in rows}
    nodelists = {row["slurm_nodelist"] for row in rows}
    hostname, nodelist = require_single_slurm_node(
        hostnames, nodelists, "pooled"
    )
    require(locales == {EXPECTED_EXECUTION_LOCALE}, "recorded locale drift")
    require(len(job_ids) == 1 and next(iter(job_ids)).isdigit(), "missing SLURM job identity")
    require(qos_values == {EXPECTED_EXECUTION_QOS}, "recorded QOS drift")
    require(all(row["slurm_cpus_per_task"] == "8" for row in rows), "recorded CPU count drift")
    require(
        all(
            row[field] == "8"
            for row in rows
            for field in (
                "omp_num_threads",
                "openblas_num_threads",
                "mkl_num_threads",
                "veclib_maximum_threads",
                "numexpr_num_threads",
            )
        ),
        "recorded thread environment drift",
    )
    require(
        all(row["slurm_partition"] == EXPECTED_EXECUTION_PARTITION for row in rows),
        "recorded partition drift",
    )
    return {
        "packages": observed,
        "hostname": hostname,
        "locale": next(iter(locales)),
        "slurm_job_id": next(iter(job_ids)),
        "slurm_nodelist": nodelist,
        "slurm_partition": EXPECTED_EXECUTION_PARTITION,
        "slurm_qos": EXPECTED_EXECUTION_QOS,
    }


def validate_input_code_manifest() -> None:
    expected_order = [
        "producer",
        "merged_dge",
        "meta_matched",
        "config",
        "gene_metadata",
        "accepted_reference",
        "runtime_contract",
        "Rscript",
        "atomic_publisher",
        "atomic_publisher_io",
        "Python",
    ]
    _, rows = read_tsv_strict(
        POOLED_REPRO_ROOT / "input_code_manifest.tsv",
        ["artifact_id", "path", "role", "sha256"],
    )
    require(
        [row["artifact_id"] for row in rows] == expected_order,
        "input/code-manifest ID/order drift",
    )

    runtime_path = (
        CANDIDATE_ROOT
        / "inputs/BG001-DECISION/contract/analysis_runtime_contract.json"
    )
    runtime = load_json(runtime_path)
    require(isinstance(runtime.get("Rscript"), dict), "runtime Rscript schema drift")
    candidate_paths = {
        "producer": CANDIDATE_ROOT / "code/reproduce_pooled_fivecohort.R",
        "merged_dge": F_FIVE_INPUT / "results/integration/merged_dge.rds",
        "meta_matched": F_FIVE_INPUT / "results/integration/meta_matched.rds",
        "config": (
            CANDIDATE_ROOT
            / "inputs/BULK-F-FIVE/frozen_model_inputs/human_datasets.yaml"
        ),
        "gene_metadata": (
            CANDIDATE_ROOT
            / "inputs/BULK-F-FIVE/frozen_model_inputs/"
            "gencode_v49_gene_metadata.tsv.gz"
        ),
        "accepted_reference": F_FIVE_INPUT / "results/integration/deg_results.csv",
        "runtime_contract": runtime_path,
        "atomic_publisher": CANDIDATE_ROOT / "code/publish_noreplace.py",
        "atomic_publisher_io": CANDIDATE_ROOT / "code/snapshot_io.py",
    }
    external_paths = {
        "Rscript": Path(str(runtime["Rscript"]["path"])),
        "Python": RNASEQ_PYTHON,
    }
    candidate_roles = {
        "producer": "resource_candidate_code",
        "merged_dge": "bg001_f_five_artifact",
        "meta_matched": "bg001_f_five_artifact",
        "config": "pooled_config",
        "gene_metadata": "gene_annotation",
        "accepted_reference": "bg001_f_five_artifact",
        "runtime_contract": (
            "bg001_analysis_bound_input;bg001_baseline_transaction"
        ),
        "atomic_publisher": "resource_candidate_code",
        "atomic_publisher_io": "resource_candidate_code",
    }
    external_roles = {
        "Rscript": "analysis_runtime_rscript",
        "Python": "resource_candidate_python",
    }
    _, source_rows = read_tsv_strict(
        CANDIDATE_ROOT / "manifests/source_selection.tsv",
        [
            "source_scope",
            "source_path",
            "candidate_path",
            "role",
            "size_bytes",
            "sha256",
        ],
    )
    source_by_candidate = {
        row["candidate_path"]: row
        for row in source_rows
    }
    require(
        len(source_by_candidate) == len(source_rows),
        "duplicate candidate path in source selection",
    )
    by_id = {row["artifact_id"]: row for row in rows}
    for artifact_id, expected_path in candidate_paths.items():
        expected_path = expected_path.resolve(strict=True)
        row = by_id[artifact_id]
        require(Path(row["path"]) == expected_path, f"input/code path drift: {artifact_id}")
        candidate_relative = expected_path.relative_to(CANDIDATE_ROOT).as_posix()
        require(
            candidate_relative in source_by_candidate,
            f"input/code absent from source selection: {artifact_id}",
        )
        source_row = source_by_candidate[candidate_relative]
        expected_hash = source_row["sha256"]
        require(
            row["role"] == source_row["role"] == candidate_roles[artifact_id],
            f"input/code sealed role drift: {artifact_id}",
        )
        require(
            row["sha256"] == expected_hash == sha256(expected_path),
            f"input/code sealed hash drift: {artifact_id}",
        )

    rscript = external_paths["Rscript"].resolve(strict=True)
    require(
        by_id["Rscript"]["path"] == str(rscript)
        and by_id["Rscript"]["role"] == external_roles["Rscript"]
        and by_id["Rscript"]["sha256"] == runtime["Rscript"]["sha256"]
        and sha256(rscript) == runtime["Rscript"]["sha256"],
        "Rscript exact identity drift",
    )
    python = external_paths["Python"].resolve(strict=True)
    require(
        by_id["Python"]["path"] == str(python)
        and by_id["Python"]["role"] == external_roles["Python"]
        and by_id["Python"]["sha256"] == RNASEQ_PYTHON_SHA256
        and sha256(python) == RNASEQ_PYTHON_SHA256,
        "Python exact identity drift",
    )


def validate_model_input_manifest() -> None:
    expected_fields = [
        "dge_path",
        "meta_path",
        "config_path",
        "gene_metadata_path",
        "accepted_reference_path",
        "n_genes",
        "n_samples",
        "n_cohorts",
        "n_control",
        "n_disease",
        "treat_lfc",
    ]
    _, rows = read_tsv_strict(
        POOLED_REPRO_ROOT / "model_input_manifest.tsv",
        expected_fields,
    )
    require(len(rows) == 1, "model-input manifest cardinality drift")
    expected = {
        "dge_path": str(
            (F_FIVE_INPUT / "results/integration/merged_dge.rds").resolve(strict=True)
        ),
        "meta_path": str(
            (F_FIVE_INPUT / "results/integration/meta_matched.rds").resolve(strict=True)
        ),
        "config_path": str(
            (
                CANDIDATE_ROOT
                / "inputs/BULK-F-FIVE/frozen_model_inputs/human_datasets.yaml"
            ).resolve(strict=True)
        ),
        "gene_metadata_path": str(
            (
                CANDIDATE_ROOT
                / "inputs/BULK-F-FIVE/frozen_model_inputs/"
                "gencode_v49_gene_metadata.tsv.gz"
            ).resolve(strict=True)
        ),
        "accepted_reference_path": str(
            (F_FIVE_INPUT / "results/integration/deg_results.csv").resolve(strict=True)
        ),
        "n_genes": "23370",
        "n_samples": "844",
        "n_cohorts": "5",
        "n_control": "157",
        "n_disease": "687",
        "treat_lfc": "0.25",
    }
    require(rows[0] == expected, "model-input manifest value drift")


def validate_scientific_reproduction() -> dict[str, object]:
    verify_candidate_snapshot()
    require_live_canonical_unchanged()
    require(POOLED_REPRO_ROOT.is_dir(), "pooled workstream is missing")
    require(not POOLED_REPRO_ROOT.is_symlink(), "pooled workstream is symlinked")
    for path in POOLED_REPRO_ROOT.rglob("*"):
        require(not path.is_symlink(), f"symlink in pooled workstream: {path}")
    require_regular_file(POOLED_REPRO_ROOT / "PRODUCER_COMPLETE.json")

    validate_producer_manifest()
    checks: list[dict[str, object]] = [
        {
            "check_id": "P01_PRODUCER_MANIFEST",
            "status": "PASS",
            "detail": "exact producer file set, bytes, and hashes verified",
        }
    ]

    accepted_path = F_FIVE_INPUT / "results/integration/deg_results.csv"
    reproduced_path = POOLED_REPRO_ROOT / "deg_results.csv"
    accepted_fields, accepted_order, accepted = read_gene_table(accepted_path)
    reproduced_fields, reproduced_order, reproduced = read_gene_table(reproduced_path)
    require(
        accepted_fields == reproduced_fields == DEG_FIELDS,
        "accepted/reproduced DEG schema differs",
    )
    require(accepted_order == reproduced_order, "accepted/reproduced gene order differs")
    require(len(accepted_order) == EXPECTED_POOLED["n_genes"], "tested-gene count drift")

    numeric_maxima = {field: 0.0 for field in NONPROBABILITY_TOLERANCES}
    probability_absolute_maxima = {field: 0.0 for field in PROBABILITY_FIELDS}
    probability_relative_maxima = {field: 0.0 for field in PROBABILITY_FIELDS}
    probability_log10_maxima = {field: 0.0 for field in PROBABILITY_FIELDS}
    probability_zero_nonzero_counts = {field: 0 for field in PROBABILITY_FIELDS}
    accepted_padj: set[str] = set()
    reproduced_padj: set[str] = set()
    accepted_lfsr: set[str] = set()
    reproduced_lfsr: set[str] = set()
    accepted_treat: set[str] = set()
    reproduced_treat: set[str] = set()
    for gene in accepted_order:
        left = accepted[gene]
        right = reproduced[gene]
        require(left["symbol"] == right["symbol"], f"symbol drift for {gene}")
        for field, tolerance in NONPROBABILITY_TOLERANCES.items():
            left_value = float(left[field])
            right_value = float(right[field])
            require(
                math.isfinite(left_value) and math.isfinite(right_value),
                f"nonfinite {field}: {gene}",
            )
            delta = abs(left_value - right_value)
            numeric_maxima[field] = max(numeric_maxima[field], delta)
            require(delta <= tolerance, f"{field} reproduction failed for {gene}: {delta}")
        for field in PROBABILITY_FIELDS:
            absolute_delta, relative_delta, log10_delta, zero_nonzero = compare_probability(
                float(left[field]), float(right[field]), field, gene
            )
            probability_absolute_maxima[field] = max(
                probability_absolute_maxima[field], absolute_delta
            )
            probability_relative_maxima[field] = max(
                probability_relative_maxima[field], relative_delta
            )
            if log10_delta is not None:
                probability_log10_maxima[field] = max(
                    probability_log10_maxima[field], log10_delta
                )
            probability_zero_nonzero_counts[field] += int(zero_nonzero)
        if float(left["padj"]) < 0.05:
            accepted_padj.add(gene)
        if float(right["padj"]) < 0.05:
            reproduced_padj.add(gene)
        if float(left["lfsr"]) < 0.05:
            accepted_lfsr.add(gene)
        if float(right["lfsr"]) < 0.05:
            reproduced_lfsr.add(gene)
        if float(left["treat_fdr"]) < 0.05:
            accepted_treat.add(gene)
        if float(right["treat_fdr"]) < 0.05:
            reproduced_treat.add(gene)
    require(accepted_padj == reproduced_padj, "padj<0.05 membership differs")
    require(accepted_lfsr == reproduced_lfsr, "lfsr<0.05 membership differs")
    require(accepted_treat == reproduced_treat, "TREAT membership differs")
    require(len(reproduced_treat) == EXPECTED_POOLED["n_treat"], "TREAT count drift")
    n_up = sum(float(reproduced[gene]["logFC"]) > 0 for gene in reproduced_treat)
    n_down = sum(float(reproduced[gene]["logFC"]) < 0 for gene in reproduced_treat)
    require(
        n_up == EXPECTED_POOLED["n_up"] and n_down == EXPECTED_POOLED["n_down"],
        "TREAT direction-count drift",
    )
    checks.append(
        {
            "check_id": "P02_FULL_NUMERICAL_REPRODUCTION",
            "status": "PASS",
            "detail": json.dumps(
                {
                    "numeric_absolute_maxima": numeric_maxima,
                    "probability_absolute_maxima": probability_absolute_maxima,
                    "probability_relative_maxima": probability_relative_maxima,
                    "probability_log10_maxima": probability_log10_maxima,
                    "probability_zero_nonzero_counts": probability_zero_nonzero_counts,
                    "probability_tolerances": {
                        "absolute": PROBABILITY_ABSOLUTE_TOLERANCE,
                        "relative": PROBABILITY_RELATIVE_TOLERANCE,
                    },
                    "threshold_membership_counts": {
                        "padj_lt_0.05": len(reproduced_padj),
                        "lfsr_lt_0.05": len(reproduced_lfsr),
                        "treat_fdr_lt_0.05": len(reproduced_treat),
                    },
                },
                sort_keys=True,
            ),
        }
    )

    accepted_fields, accepted_design = read_tsv_strict(
        F_FIVE_INPUT / "results/integration/model_design.tsv"
    )
    reproduced_fields, reproduced_design = read_tsv_strict(
        POOLED_REPRO_ROOT / "model_design.tsv"
    )
    require(accepted_fields == reproduced_fields, "model-design schema differs")
    require(accepted_design == reproduced_design, "model design differs from accepted F_five")
    require(len(reproduced_design) == EXPECTED_POOLED["n_samples"], "model-design row-count drift")
    design_order = [row["sample_id"] for row in reproduced_design]
    require(len(set(design_order)) == len(design_order), "duplicate model-design sample")
    checks.append(
        {
            "check_id": "P03_DESIGN_IDENTITY",
            "status": "PASS",
            "detail": f"{len(reproduced_design)} ordered model rows are identical",
        }
    )

    weight_fields = ["sample_id", "sample_weight"]
    _, accepted_weights = read_tsv_strict(
        F_FIVE_INPUT / "results/integration/lvqw_sample_weights.tsv", weight_fields
    )
    _, reproduced_weights = read_tsv_strict(
        POOLED_REPRO_ROOT / "lvqw_sample_weights.tsv", weight_fields
    )
    require(len(accepted_weights) == len(reproduced_weights) == EXPECTED_POOLED["n_samples"], "weight row-count drift")
    accepted_weight_order = [row["sample_id"] for row in accepted_weights]
    reproduced_weight_order = [row["sample_id"] for row in reproduced_weights]
    require(len(set(accepted_weight_order)) == len(accepted_weight_order), "duplicate accepted weight sample")
    require(len(set(reproduced_weight_order)) == len(reproduced_weight_order), "duplicate reproduced weight sample")
    require(accepted_weight_order == reproduced_weight_order == design_order, "weight sample order drift")
    weight_deltas: list[float] = []
    for left, right in zip(accepted_weights, reproduced_weights, strict=True):
        left_value = float(left["sample_weight"])
        right_value = float(right["sample_weight"])
        require(
            math.isfinite(left_value)
            and math.isfinite(right_value)
            and left_value > 0
            and right_value > 0,
            f"invalid sample weight: {left['sample_id']}",
        )
        weight_deltas.append(abs(left_value - right_value))
    weight_delta = max(weight_deltas)
    require(weight_delta <= 2e-12, f"sample-weight reproduction failed: {weight_delta}")
    checks.append(
        {
            "check_id": "P04_WEIGHT_IDENTITY",
            "status": "PASS",
            "detail": f"max absolute sample-weight delta={weight_delta:.3g}",
        }
    )

    environment = validate_environment()
    validate_input_code_manifest()
    validate_model_input_manifest()
    completion = load_json(POOLED_REPRO_ROOT / "PRODUCER_COMPLETE.json")
    producer_manifest = POOLED_REPRO_ROOT / "producer_artifact_manifest.tsv"
    expected_completion = {
        "status": "PRODUCER_COMPLETE_PENDING_VALIDATION",
        "candidate_id": CANDIDATE_ID,
        "workstream": "BULK-POOLED-REPRO",
        "n_genes": EXPECTED_POOLED["n_genes"],
        "n_samples": EXPECTED_POOLED["n_samples"],
        "n_treat": EXPECTED_POOLED["n_treat"],
        "n_up": EXPECTED_POOLED["n_up"],
        "n_down": EXPECTED_POOLED["n_down"],
        "producer_artifact_manifest_sha256": sha256(producer_manifest),
        "canonical_write": False,
    }
    require(
        set(completion) == {*expected_completion, "created_at_utc"},
        "producer completion schema drift",
    )
    for key, expected in expected_completion.items():
        if key == "canonical_write":
            require(
                completion[key] is False,
                "producer canonical-write flag is not literal false",
            )
        else:
            require(completion[key] == expected, f"producer completion drift: {key}")
    require(
        isinstance(completion["created_at_utc"], str)
        and bool(completion["created_at_utc"]),
        "producer completion timestamp missing",
    )
    checks.extend(
        [
            {
                "check_id": "P05_RUNTIME_AND_INPUT_IDENTITY",
                "status": "PASS",
                "detail": json.dumps(environment, sort_keys=True),
            },
            {
                "check_id": "P06_NONCANONICAL_CANARY",
                "status": "PASS",
                "detail": f"live canonical remained {CANONICAL_REFERENCE_SHA256}; producer canonical_write=false",
            },
        ]
    )
    return {
        "checks": checks,
        "numeric_maxima": numeric_maxima,
        "probability_absolute_maxima": probability_absolute_maxima,
        "probability_relative_maxima": probability_relative_maxima,
        "probability_log10_maxima": probability_log10_maxima,
        "probability_zero_nonzero_counts": probability_zero_nonzero_counts,
        "probability_tolerances": {
            "absolute": PROBABILITY_ABSOLUTE_TOLERANCE,
            "relative": PROBABILITY_RELATIVE_TOLERANCE,
        },
        "threshold_membership_counts": {
            "padj_lt_0.05": len(reproduced_padj),
            "lfsr_lt_0.05": len(reproduced_lfsr),
            "treat_fdr_lt_0.05": len(reproduced_treat),
        },
        "weight_delta": weight_delta,
    }


def build_validation_manifest(
    temporary: Path, validation_files: list[Path]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(
        path
        for path in POOLED_REPRO_ROOT.rglob("*")
        if path.is_file() and path != READY_PATH and VALIDATION_ROOT not in path.parents
    ):
        rows.append(
            {
                "artifact_path": path.relative_to(POOLED_REPRO_ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    for path in validation_files:
        rows.append(
            {
                "artifact_path": f"validation/{path.name}",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    rows.sort(key=lambda row: str(row["artifact_path"]))
    return rows


def publish_validation(result: dict[str, object]) -> None:
    require(
        not READY_PATH.exists() and not READY_PATH.is_symlink(),
        "validation seal already exists",
    )
    require(
        not VALIDATION_ROOT.exists() and not VALIDATION_ROOT.is_symlink(),
        "validation bundle already exists",
    )
    temporary = (
        POOLED_REPRO_ROOT.parent
        / f".{POOLED_REPRO_ROOT.name}.validation.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    )
    require(
        not temporary.exists() and not temporary.is_symlink(),
        "validation temporary path collision",
    )
    temporary.mkdir(mode=0o750)
    try:
        checks_path = temporary / "validation_checks.tsv"
        report_path = temporary / "validation_report.json"
        manifest_path = temporary / "validated_artifact_manifest.tsv"
        write_exclusive(
            checks_path,
            tsv_payload(result["checks"], ["check_id", "status", "detail"]),
        )
        report = {
            "status": "VALIDATED",
            "candidate_id": CANDIDATE_ID,
            "workstream": "BULK-POOLED-REPRO",
            "n_checks": len(result["checks"]),
            "max_numeric_deltas": result["numeric_maxima"],
            "max_probability_absolute_deltas": result["probability_absolute_maxima"],
            "max_probability_relative_deltas": result["probability_relative_maxima"],
            "max_probability_log10_deltas": result["probability_log10_maxima"],
            "probability_zero_nonzero_counts": result["probability_zero_nonzero_counts"],
            "probability_tolerances": result["probability_tolerances"],
            "threshold_membership_counts": result["threshold_membership_counts"],
            "max_sample_weight_delta": result["weight_delta"],
            "live_canonical_sha256": CANONICAL_REFERENCE_SHA256,
            "canonical_promotion_status": "not_promoted",
            "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_exclusive(
            report_path,
            json.dumps(report, indent=2, sort_keys=True) + "\n",
        )
        manifest_rows = build_validation_manifest(
            temporary,
            [checks_path, report_path],
        )
        write_exclusive(
            manifest_path,
            tsv_payload(manifest_rows, ["artifact_path", "bytes", "sha256"]),
        )
        require_live_canonical_unchanged()
        verify_candidate_snapshot()
        publication = publish_directory_noreplace(
            temporary,
            VALIDATION_ROOT,
            marker="validated_artifact_manifest.tsv",
        )
        require(
            publication.get("source_retained") is True,
            "validation publisher did not report retained source",
        )
        require(
            temporary.is_dir() and not temporary.is_symlink(),
            "validation publication source was not retained safely",
        )
        require_regular_file(manifest_path)
        require_regular_file(VALIDATION_MANIFEST)
        require(
            sha256(manifest_path) == sha256(VALIDATION_MANIFEST),
            "published validation commit-marker hash drift",
        )
    except Exception:
        print(
            f"FAILED validation temporary retained for audit: {temporary}",
            file=os.sys.stderr,
        )
        raise
    finalize_ready_marker(result)


def verify_validation_bundle_without_ready(
    expected_result: dict[str, object] | None = None,
) -> dict[str, object]:
    require(
        POOLED_REPRO_ROOT.is_dir() and not POOLED_REPRO_ROOT.is_symlink(),
        "pooled workstream root is invalid",
    )
    require(
        VALIDATION_ROOT.is_dir() and not VALIDATION_ROOT.is_symlink(),
        "validation bundle is missing or symlinked",
    )
    for path in POOLED_REPRO_ROOT.rglob("*"):
        require(not path.is_symlink(), f"symlink in pooled validation tree: {path}")
    _, rows = read_tsv_strict(
        VALIDATION_MANIFEST,
        ["artifact_path", "bytes", "sha256"],
    )
    require(rows, "validated artifact manifest is empty")
    observed_paths: set[str] = set()
    for row in rows:
        relative = row["artifact_path"]
        require(relative not in observed_paths, f"duplicate validated artifact: {relative}")
        require(
            not Path(relative).is_absolute() and ".." not in Path(relative).parts,
            "unsafe validated path",
        )
        observed_paths.add(relative)
        path = POOLED_REPRO_ROOT / relative
        require_regular_file(path)
        require(
            path.stat().st_size == int(row["bytes"]),
            f"validated byte drift: {relative}",
        )
        require(sha256(path) == row["sha256"], f"validated hash drift: {relative}")
    actual_paths = {
        path.relative_to(POOLED_REPRO_ROOT).as_posix()
        for path in POOLED_REPRO_ROOT.rglob("*")
        if path.is_file() and path not in {READY_PATH, VALIDATION_MANIFEST}
    }
    require(observed_paths == actual_paths, "validated manifest file-set mismatch")

    _, checks = read_tsv_strict(
        VALIDATION_ROOT / "validation_checks.tsv",
        ["check_id", "status", "detail"],
    )
    require(
        {row["check_id"] for row in checks} == EXPECTED_CHECK_IDS
        and len(checks) == len(EXPECTED_CHECK_IDS),
        "validation check-ID family drift",
    )
    require(
        all(row["status"] == "PASS" for row in checks),
        "non-PASS validation check",
    )
    if expected_result is not None:
        expected_checks = [
            {
                "check_id": str(row["check_id"]),
                "status": str(row["status"]),
                "detail": str(row["detail"]),
            }
            for row in expected_result["checks"]
        ]
        require(checks == expected_checks, "validation checks differ from recomputation")

    report = load_json(VALIDATION_ROOT / "validation_report.json")
    required_report_keys = {
        "status",
        "candidate_id",
        "workstream",
        "n_checks",
        "max_numeric_deltas",
        "max_probability_absolute_deltas",
        "max_probability_relative_deltas",
        "max_probability_log10_deltas",
        "probability_zero_nonzero_counts",
        "probability_tolerances",
        "threshold_membership_counts",
        "max_sample_weight_delta",
        "live_canonical_sha256",
        "canonical_promotion_status",
        "validated_at_utc",
    }
    require(set(report) == required_report_keys, "validation report schema drift")
    require(report.get("status") == "VALIDATED", "validation report status drift")
    require(report.get("candidate_id") == CANDIDATE_ID, "validation report candidate drift")
    require(
        report.get("workstream") == "BULK-POOLED-REPRO",
        "validation report workstream drift",
    )
    require(report.get("n_checks") == len(checks), "validation report check-count drift")
    require(
        report.get("live_canonical_sha256") == CANONICAL_REFERENCE_SHA256,
        "canonical canary drift",
    )
    require(
        report.get("canonical_promotion_status") == "not_promoted",
        "canonical promotion status drift",
    )
    require(
        isinstance(report.get("validated_at_utc"), str)
        and bool(report["validated_at_utc"]),
        "validation timestamp missing",
    )
    if expected_result is not None:
        require(
            report["max_numeric_deltas"] == expected_result["numeric_maxima"]
            and report["max_probability_absolute_deltas"]
            == expected_result["probability_absolute_maxima"]
            and report["max_probability_relative_deltas"]
            == expected_result["probability_relative_maxima"]
            and report["max_probability_log10_deltas"]
            == expected_result["probability_log10_maxima"]
            and report["probability_zero_nonzero_counts"]
            == expected_result["probability_zero_nonzero_counts"]
            and report["probability_tolerances"]
            == expected_result["probability_tolerances"]
            and report["threshold_membership_counts"]
            == expected_result["threshold_membership_counts"]
            and report["max_sample_weight_delta"] == expected_result["weight_delta"],
            "validation report differs from scientific recomputation",
        )
    return report


def ready_payload() -> dict[str, object]:
    return {
        "status": "VALIDATED",
        "candidate_id": CANDIDATE_ID,
        "workstream": "BULK-POOLED-REPRO",
        "validation_report_sha256": sha256(
            VALIDATION_ROOT / "validation_report.json"
        ),
        "validation_checks_sha256": sha256(
            VALIDATION_ROOT / "validation_checks.tsv"
        ),
        "validated_artifact_manifest_sha256": sha256(VALIDATION_MANIFEST),
        "live_canonical_sha256": CANONICAL_REFERENCE_SHA256,
        "canonical_promotion_status": "not_promoted",
    }


def seal_read_only() -> None:
    require(
        POOLED_REPRO_ROOT.is_dir() and not POOLED_REPRO_ROOT.is_symlink(),
        "pooled workstream root is invalid before sealing",
    )
    descendants = list(POOLED_REPRO_ROOT.rglob("*"))
    for path in descendants:
        require(not path.is_symlink(), f"symlink before read-only seal: {path}")
    for path in descendants:
        if path.is_file():
            chmod_nofollow(path, 0o440, directory=False)
    for path in sorted(
        (path for path in descendants if path.is_dir()),
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        chmod_nofollow(path, 0o550, directory=True)
    chmod_nofollow(POOLED_REPRO_ROOT, 0o550, directory=True)


def finalize_ready_marker(result: dict[str, object]) -> None:
    verify_validation_bundle_without_ready(result)
    verify_candidate_snapshot()
    require_live_canonical_unchanged()
    expected = ready_payload()
    require(not READY_PATH.is_symlink(), "VALIDATED marker may not be a symlink")
    if READY_PATH.exists():
        require(
            load_json(READY_PATH) == expected,
            "existing VALIDATED marker is inconsistent",
        )
    else:
        temporary = (
            POOLED_REPRO_ROOT.parent
            / f".{POOLED_REPRO_ROOT.name}.VALIDATED.tmp.{os.getpid()}.{uuid.uuid4().hex}"
        )
        require(
            not temporary.exists() and not temporary.is_symlink(),
            "VALIDATED temporary path collision",
        )
        write_exclusive(
            temporary,
            json.dumps(expected, indent=2, sort_keys=True) + "\n",
        )
        verify_validation_bundle_without_ready(result)
        verify_candidate_snapshot()
        require_live_canonical_unchanged()
        publication = publish_file_noreplace(temporary, READY_PATH)
        require(
            publication.get("source_retained") is True,
            "VALIDATED publisher did not report retained source",
        )
        require_regular_file(temporary)
        require_regular_file(READY_PATH)
        require(
            sha256(temporary) == sha256(READY_PATH),
            "published VALIDATED marker hash drift",
        )
    verify_validation_bundle_without_ready(result)
    seal_read_only()


def verify_existing_seal(
    expected_result: dict[str, object] | None = None,
) -> dict[str, object]:
    verify_candidate_snapshot()
    require_live_canonical_unchanged()
    validate_producer_manifest()
    report = verify_validation_bundle_without_ready(expected_result)
    ready = load_json(READY_PATH)
    require(
        ready == ready_payload(),
        "VALIDATED marker differs from bound validation bundle",
    )
    for path in POOLED_REPRO_ROOT.rglob("*"):
        require(not path.is_symlink(), f"symlink in validated workstream: {path}")
        expected_mode = 0o440 if path.is_file() else 0o550
        require(
            stat.S_IMODE(path.stat().st_mode) == expected_mode,
            f"validated mode drift: {path}",
        )
    require(
        stat.S_IMODE(POOLED_REPRO_ROOT.stat().st_mode) == 0o550,
        "validated workstream root mode drift",
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    result = validate_scientific_reproduction()
    if args.verify_only:
        print(
            json.dumps(
                verify_existing_seal(result),
                indent=2,
                sort_keys=True,
            )
        )
        return
    if READY_PATH.exists() or READY_PATH.is_symlink():
        require(
            READY_PATH.is_file() and not READY_PATH.is_symlink(),
            "invalid VALIDATED marker",
        )
        require(
            VALIDATION_ROOT.is_dir() and not VALIDATION_ROOT.is_symlink(),
            "missing validation bundle",
        )
        finalize_ready_marker(result)
    elif VALIDATION_ROOT.exists() or VALIDATION_ROOT.is_symlink():
        require(
            VALIDATION_ROOT.is_dir() and not VALIDATION_ROOT.is_symlink(),
            "invalid validation bundle",
        )
        finalize_ready_marker(result)
    else:
        publish_validation(result)
    print(
        json.dumps(
            verify_existing_seal(result),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except (ContractError, SnapshotIOError) as error:
        raise SystemExit(f"VALIDATION ERROR: {error}") from error
