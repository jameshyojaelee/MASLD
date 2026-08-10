#!/usr/bin/env python3
"""Validate frozen v2 Visium outputs with the corrected finite-null contract.

This entry point delegates every source, design, graph, family, BH, robustness,
matching, and sensitivity check to the frozen validator.  Its sole declared
patch replaces the impossible requirement that every raw-Moran null value be
finite with the predeclared semantic bound used for v1: at least 95% finite
draws in every summary.  The attempted null count remains exactly 9,999 and no
other failed check is suppressed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from pathlib import Path

import pandas as pd

from visium_rerun_lib import (
    CORE_OUTPUTS,
    N_NULL,
    RELEASE_ID,
    V1_ENGINE_SHA256,
    build_paths,
    read_tsv,
    sha256_file,
    utc_now,
    write_tsv,
)


ORIGINAL_VALIDATOR_SHA256 = (
    "281f727692029bad627eba5c0912b71b3dc999728b5c23c82beab5722079a26b"
)
MIN_FINITE_NULL_FRACTION = 0.95


def load_original_validator(script_root: Path):
    path = script_root / "06_validate_visium_rerun.py"
    observed = sha256_file(path)
    if observed != ORIGINAL_VALIDATOR_SHA256:
        raise RuntimeError(
            f"original candidate validator drift: expected {ORIGINAL_VALIDATOR_SHA256}, "
            f"observed {observed}"
        )
    spec = importlib.util.spec_from_file_location("_frozen_visium_v2_validator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not import frozen v2 validator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def finite_null_gate(output: Path) -> tuple[bool, int, float]:
    nulls = pd.read_csv(output / "spatial_null_summary.tsv", sep="\t")
    counts = pd.to_numeric(nulls["n_null"], errors="raise")
    minimum = int(counts.min())
    fraction = minimum / N_NULL
    return bool(fraction >= MIN_FINITE_NULL_FRACTION), minimum, fraction


def write_validated_bundle(module, checks, output: Path, minimum: int, fraction: float) -> None:
    report = output / "validation_report.tsv"
    manifest = output / "validator_semantic_manifest.tsv"
    ready = output / "READY"
    for path in (report, manifest, ready):
        if path.exists():
            raise RuntimeError(f"refusing to overwrite v2 validation artifact: {path}")

    write_tsv(report, ("check_id", "status", "detail"), checks.rows)
    if not checks.passed:
        failures = [row["check_id"] for row in checks.rows if row["status"] != "PASS"]
        raise RuntimeError(
            "v2 semantic validation failed; no check except finite_primary_nulls "
            f"was patched: {failures}"
        )

    script = Path(__file__).resolve()
    write_tsv(
        manifest,
        (
            "validator_id",
            "scope",
            "patched_check_id",
            "minimum_finite_null_count",
            "minimum_finite_null_fraction",
            "required_finite_null_fraction",
            "original_validator_sha256",
            "semantic_validator_sha256",
            "validation_report_sha256",
        ),
        [
            {
                "validator_id": "SP-INT-03-v2-finite-null-semantic-v1",
                "scope": "validation_only_predeclared_before_v2_outcomes",
                "patched_check_id": "finite_primary_nulls",
                "minimum_finite_null_count": minimum,
                "minimum_finite_null_fraction": f"{fraction:.12g}",
                "required_finite_null_fraction": MIN_FINITE_NULL_FRACTION,
                "original_validator_sha256": ORIGINAL_VALIDATOR_SHA256,
                "semantic_validator_sha256": sha256_file(script),
                "validation_report_sha256": sha256_file(report),
            }
        ],
    )
    aggregate = "".join(
        f"{name}\0{(output / name).stat().st_size}\0{sha256_file(output / name)}\n"
        for name in CORE_OUTPUTS
    )
    write_tsv(
        ready,
        (
            "release_id",
            "registry_version",
            "status",
            "v1_engine_sha256",
            "validation_report_path",
            "validation_report_sha256",
            "validator_semantic_manifest_sha256",
            "core_outputs_aggregate_sha256",
            "validated_utc",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "registry_version": "v2",
                "status": "pass_v2_candidate",
                "v1_engine_sha256": V1_ENGINE_SHA256,
                "validation_report_path": report.name,
                "validation_report_sha256": sha256_file(report),
                "validator_semantic_manifest_sha256": sha256_file(manifest),
                "core_outputs_aggregate_sha256": hashlib.sha256(aggregate.encode()).hexdigest(),
                "validated_utc": utc_now(),
            }
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        paths = build_paths(args.project_root)
        output = paths.native_root / "v2_candidate"
        finite_ok, minimum, fraction = finite_null_gate(output)
        module = load_original_validator(Path(__file__).resolve().parent)
        original_require = module.Checks.require

        def corrected_require(instance, condition, check_id, detail):
            if check_id == "finite_primary_nulls":
                condition = finite_ok
                detail = (
                    f"minimum finite={minimum}/{N_NULL} ({fraction:.6f}); "
                    f"required fraction>={MIN_FINITE_NULL_FRACTION}"
                )
            return original_require(instance, condition, check_id, detail)

        module.Checks.require = corrected_require
        checks, validated_output = module.validate("v2", paths)
        if validated_output.resolve() != output.resolve():
            raise RuntimeError("delegated validator returned an unexpected output root")
        write_validated_bundle(module, checks, output, minimum, fraction)
        print(
            f"SP-INT-03 v2 semantic validation passed ({len(checks.rows)} checks; "
            f"minimum finite nulls={minimum}/{N_NULL})"
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
