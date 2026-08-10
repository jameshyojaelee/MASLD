#!/usr/bin/env python3
"""Freeze the reviewed phase-C executor before any held-out outcome access."""

from __future__ import annotations

import argparse
import datetime as dt
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    RELEASE_ID,
    atomic_write_text,
    parse_bool,
    read_tsv,
    require_within,
    sha256_file,
    stable_bundle_sha256,
    write_tsv,
)


SPECIFICATION_SHA256 = (
    "1f99b7a7a95d9aebbbb911bc2317860f4a3fada0d1da2a4932d5c1ab190c4cdd"
)
EXPECTED_VERSIONS = {
    "python": "3.10",
    "numpy": "2.2.6",
    "pandas": "2.3.3",
    "scipy": "1.15.2",
    "statsmodels": "0.14.5",
}
SCRIPT_DIR = Path(__file__).resolve().parent
CODE_PATHS = [
    SCRIPT_DIR / "09_freeze_phase_c_code.py",
    SCRIPT_DIR / "10_execute_one_pass.py",
    SCRIPT_DIR / "11_validate_phase_c.py",
    SCRIPT_DIR / "myojin_firewall_common.py",
    SCRIPT_DIR / "phase_c_stats.py",
    SCRIPT_DIR / "phase_c_execution_contract.yaml",
    SCRIPT_DIR / "run_phase_c.sbatch",
    SCRIPT_DIR / "tests/test_phase_c_stats.py",
]
FROZEN_INPUT_NAMES = [
    "SEALED",
    "FIREWALL_READY",
    "DEPMAP_READY",
    "blind_analysis_spec.yaml",
    "specification_sha256.txt",
    "masked_screen_schema.tsv",
    "prediction_manifest.tsv",
    "known_hit_exclusion.tsv",
    "screen_gene_universe.tsv",
    "gene_mapping_audit.tsv",
    "covariate_coverage_audit.tsv",
    "program_testability.tsv",
    "gate_status.tsv",
    "deviations.tsv",
    "firewall_input_manifest.tsv",
    "depmap_source_manifest.tsv",
    "depmap_gene_mapping_audit.tsv",
    "hlf_identity_audit.tsv",
    "source_manifest.tsv",
    "source/raw/Myojin_Supplemental_Table_S1_C218.xlsx",
]
FORBIDDEN_PREEXISTING = [
    "UNSEALED",
    "PHASE_C_COMPLETE",
    "PHASE_C_VALIDATED",
    "screen_results.tsv",
    "class_effects.tsv",
    "class_protective_hits.tsv",
    "program_effects.tsv",
    "fig5_verdict.tsv",
]


def one_row(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one row in {path}; observed {len(rows)}")
    return rows[0]


def validate_environment() -> dict[str, str]:
    observed = {
        "python": ".".join(platform.python_version().split(".")[:2]),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
    }
    if observed != EXPECTED_VERSIONS:
        raise ValueError(
            f"Phase-C environment mismatch: {observed} != {EXPECTED_VERSIONS}"
        )
    return observed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    args = parser.parse_args()
    outdir = require_within(args.outdir, CANDIDATE_ROOT)

    missing_code = [str(path) for path in CODE_PATHS if not path.is_file()]
    missing_inputs = [
        name for name in FROZEN_INPUT_NAMES if not (outdir / name).is_file()
    ]
    if missing_code or missing_inputs:
        raise FileNotFoundError(
            f"Missing phase-C code={missing_code}; inputs={missing_inputs}"
        )
    forbidden = [name for name in FORBIDDEN_PREEXISTING if (outdir / name).exists()]
    if forbidden:
        raise ValueError(f"Cannot freeze after phase C began: {forbidden}")
    if (outdir / "BLOCKED").exists():
        raise ValueError("Core firewall is BLOCKED")
    if read_tsv(outdir / "unseal_record.tsv"):
        raise ValueError("Unseal record is already populated")

    sealed = one_row(outdir / "SEALED")
    ready = one_row(outdir / "FIREWALL_READY")
    if (
        sealed.get("specification_sha256") != SPECIFICATION_SHA256
        or ready.get("specification_sha256") != SPECIFICATION_SHA256
        or sealed.get("external_outcomes_read") != "FALSE"
        or ready.get("external_outcomes_read") != "FALSE"
    ):
        raise ValueError(
            "SEALED/FIREWALL_READY do not authorize the frozen specification"
        )
    spec_line = (
        (outdir / "specification_sha256.txt")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )
    if spec_line != f"specification_bundle_sha256\t{SPECIFICATION_SHA256}":
        raise ValueError("Specification hash file does not match the authorized seal")

    gates = read_tsv(outdir / "gate_status.tsv")
    if len(gates) != 12 or not all(
        parse_bool(row["passed"]) and row["external_outcomes_read"] == "FALSE"
        for row in gates
    ):
        raise ValueError(
            f"Expected all 12 frozen source gates to pass; observed {len(gates)}"
        )
    observed_versions = validate_environment()

    code_bundle_sha256, code_rows = stable_bundle_sha256(CODE_PATHS, PROJECT_ROOT)
    input_bundle_sha256, input_rows = stable_bundle_sha256(
        [outdir / name for name in FROZEN_INPUT_NAMES], PROJECT_ROOT
    )
    manifest_rows = [dict(row, role="phase_c_code") for row in code_rows]
    manifest_rows.extend(dict(row, role="frozen_input") for row in input_rows)
    manifest_rows.extend(
        {
            "relative_path": f"environment::{name}",
            "bytes": "not_applicable",
            "sha256": value,
            "role": "environment_version",
        }
        for name, value in sorted(observed_versions.items())
    )
    write_tsv(
        outdir / "phase_c_code_manifest.tsv",
        [
            dict(
                row,
                specification_sha256=SPECIFICATION_SHA256,
                release_id=RELEASE_ID,
            )
            for row in manifest_rows
        ],
        [
            "role",
            "relative_path",
            "bytes",
            "sha256",
            "specification_sha256",
            "release_id",
        ],
    )
    frozen_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    marker = outdir / "PHASE_C_CODE_READY"
    atomic_write_text(
        marker,
        "release_id\tstatus\tspecification_sha256\tcode_bundle_sha256\t"
        "input_bundle_sha256\tmanifest_sha256\tfrozen_at_utc\texternal_outcomes_read\n"
        f"{RELEASE_ID}\treviewed_code_frozen\t{SPECIFICATION_SHA256}\t"
        f"{code_bundle_sha256}\t{input_bundle_sha256}\t"
        f"{sha256_file(outdir / 'phase_c_code_manifest.tsv')}\t{frozen_at}\tFALSE\n",
    )
    print(
        "PHASE_C_CODE_FROZEN "
        f"specification_sha256={SPECIFICATION_SHA256} "
        f"code_bundle_sha256={code_bundle_sha256} outcomes_read=FALSE"
    )


if __name__ == "__main__":
    main()
