#!/usr/bin/env python3
"""Run the frozen Plan 11 validator with one validation-only iterator repair.

The pre-outcome-frozen validator hashes to ``6ce677...c8265``.  Its
``source_abundance_fixture`` intentionally samples the first two ``v.mat``
donors, but the strict ``zip`` paired those two references with all 16 rows of
the independently certified gene-axis table.  Python therefore raised before
any validation report or READY seal was written.

This post-outcome wrapper does not alter the frozen validator, model, result
tables, multiplicity family, robustness rule, or acceptance criteria.  It
temporarily exposes the intended first two certified axis rows only while that
source-value fixture runs, then restores the original validator global.  Both
the live validator and its pre-outcome manifest row must match the pinned hash.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


FROZEN_VALIDATOR_SHA256 = (
    "6ce677fefd6f2eae262fcd8bbb5496255b7c42f8e9cdf50f87ea9420ed0c8265"
)
FAILED_VALIDATION_JOB_ID = "19637749"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_frozen_validator(script_dir: Path):
    path = script_dir / "04_validate_bundle.py"
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    spec = importlib.util.spec_from_file_location("yakubovsky_frozen_validator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load frozen validator: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if module.sha256_file(path) != FROZEN_VALIDATOR_SHA256:
        raise module.ContractError("Frozen Plan 11 validator hash drift")
    return module, path


def assert_failed_run_left_no_terminal_products(module, output: Path) -> None:
    forbidden = [
        output / "validation_report.tsv",
        output / "validation_execution_manifest.tsv",
        output / "terminal_stage_seal.tsv",
        output / "READY",
        output / "plan13_adapter" / "PLAN13_ADAPTER_READY",
    ]
    present = [str(path) for path in forbidden if path.exists()]
    if present:
        raise module.ContractError(
            "Validation hotfix refuses pre-existing terminal products: "
            + ", ".join(present)
        )


def assert_preoutcome_hash(module, output: Path, validator_path: Path) -> str:
    freeze_path = output / "analysis_freeze_manifest.tsv"
    freeze = pd.read_csv(freeze_path, sep="\t", dtype=str, keep_default_na=False)
    row = freeze[
        (freeze["record_type"] == "file_input")
        & (freeze["name"] == "producer_code::04_validate_bundle.py")
    ]
    if len(row) != 1:
        raise module.ContractError("Pre-outcome freeze lacks one validator row")
    record = row.iloc[0]
    if (
        Path(record["value"]).resolve() != validator_path.resolve()
        or record["sha256"] != FROZEN_VALIDATOR_SHA256
        or module.bool_value(record["lipid_outcome_read"])
    ):
        raise module.ContractError("Pre-outcome validator certificate drift")
    return module.sha256_file(freeze_path)


def install_fixture_repair(module) -> None:
    original_fixture = module.source_abundance_fixture
    original_axis_validator = module.validate_gene_axis_manifest

    def repaired_fixture(paths, sample, output):
        def first_two_axis_rows(*args, **kwargs):
            contract = original_axis_validator(*args, **kwargs)
            if len(contract.rows) < 2:
                raise module.ContractError(
                    "Validation fixture requires at least two certified donor axes"
                )
            return dataclasses.replace(contract, rows=contract.rows.head(2).copy())

        module.validate_gene_axis_manifest = first_two_axis_rows
        try:
            return original_fixture(paths, sample, output)
        finally:
            module.validate_gene_axis_manifest = original_axis_validator

    module.source_abundance_fixture = repaired_fixture


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    module, validator_path = load_frozen_validator(script_dir)
    output = None
    for index, argument in enumerate(sys.argv[:-1]):
        if argument == "--output-dir":
            output = Path(sys.argv[index + 1]).resolve()
            break
    if output is None:
        raise module.ContractError(
            "The validation hotfix requires an explicit --output-dir"
        )
    assert_failed_run_left_no_terminal_products(module, output)
    freeze_sha256 = assert_preoutcome_hash(module, output, validator_path)
    install_fixture_repair(module)
    module.main()

    wrapper_path = Path(__file__).resolve()
    manifest = pd.DataFrame(
        [
            {
                "release_id": module.RELEASE_ID,
                "status": "validation_only_iterator_hotfix_applied",
                "failed_validation_job_id": FAILED_VALIDATION_JOB_ID,
                "frozen_validator": str(validator_path),
                "frozen_validator_sha256": FROZEN_VALIDATOR_SHA256,
                "analysis_freeze_manifest_sha256": freeze_sha256,
                "hotfix_wrapper": str(wrapper_path),
                "hotfix_wrapper_sha256": module.sha256_file(wrapper_path),
                "repair_scope": (
                    "source_abundance_fixture_first_two_references_against_"
                    "first_two_certified_axis_rows"
                ),
                "model_or_result_recomputed": False,
                "acceptance_rule_changed": False,
                "post_outcome_validation_only": True,
                "completed_utc": utc_now(),
            }
        ]
    )
    module.atomic_write_frame(output / "validator_hotfix_manifest.tsv", manifest)


if __name__ == "__main__":
    main()
