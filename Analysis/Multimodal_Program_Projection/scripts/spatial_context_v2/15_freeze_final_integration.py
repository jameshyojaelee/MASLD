#!/usr/bin/env python3
"""Freeze/check the outcome-independent SP-INT-06/07 assembly code."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from contract_lib import ContractError, RELEASE_ID, read_tsv, sha256_file, write_tsv


MANIFEST_COLUMNS = ("release_id", "relative_path", "bytes", "sha256", "role")
FILES = (
    ("contract_lib.py", "semantic_v2_shared_contract"),
    ("02_build_synthetic_contract.py", "semantic_v2_contract_fixture_builder"),
    ("03_validate_contract.py", "semantic_v2_contract_validator"),
    ("08_build_real_adapters.py", "semantic_v2_real_adapter_builder"),
    ("09_validate_real_adapters.py", "semantic_v2_real_adapter_validator"),
    ("11_build_protein_atac_adapters.py", "semantic_v2_protein_atac_builder"),
    ("12_validate_protein_atac_adapters.py", "semantic_v2_protein_atac_validator"),
    ("final_integration_lib.py", "final_integration_library"),
    ("13_build_final_integration.py", "final_integration_builder"),
    ("14_validate_final_integration.py", "final_integration_validator"),
    ("15_freeze_final_integration.py", "code_freezer"),
    ("16_prepare_semantic_v2.py", "semantic_v2_isolated_root_preparer"),
    ("17_seal_semantic_v2.py", "semantic_v2_release_sealer"),
    ("../../../Spatial/scripts/yakubovsky2026/05_build_plan13_adapter.py", "immutable_historical_yakubovsky_adapter"),
    ("../../../Spatial/scripts/yakubovsky2026/09_build_plan13_adapter_semantic_v2.py", "semantic_v2_yakubovsky_adapter"),
    ("../../../Spatial/scripts/yakubovsky2026/tests/test_semantic_v2_adapter.py", "semantic_v2_yakubovsky_adapter_tests"),
    ("tests/test_final_integration.py", "semantic_fixture_tests"),
    ("run_final_integration.sbatch", "real_final_wrapper"),
    ("run_final_integration_fixture.sbatch", "fixture_wrapper"),
    ("run_semantic_v2_rebuild.sbatch", "semantic_v2_rebuild_wrapper"),
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def expected_rows(project_root: Path) -> list[dict[str, object]]:
    script_root = Path(__file__).resolve().parent
    rows = []
    for relative, role in FILES:
        path = (script_root / relative).resolve()
        if not path.is_file():
            raise ContractError(f"missing final-integration code-freeze input: {path}")
        rows.append(
            {
                "release_id": RELEASE_ID,
                "relative_path": path.relative_to(project_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "role": role,
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("write", "check"))
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    manifest = candidate_root / "final_integration_code_freeze.tsv"
    ready = candidate_root / "FINAL_INTEGRATION_CODE_READY"
    try:
        expected = expected_rows(project_root)
        if args.action == "write":
            if manifest.exists() or ready.exists():
                raise ContractError("refusing to overwrite final integration code freeze")
            write_tsv(manifest, MANIFEST_COLUMNS, expected)
            write_tsv(
                ready,
                ("release_id", "status", "manifest_sha256", "n_files", "real_outcomes_read_by_freezer"),
                [{"release_id": RELEASE_ID, "status": "outcome_independent_final_assembly_code_frozen", "manifest_sha256": sha256_file(manifest), "n_files": len(expected), "real_outcomes_read_by_freezer": False}],
            )
            print(f"PASS: froze {len(expected)} final-integration code files")
        else:
            _, observed = read_tsv(manifest, MANIFEST_COLUMNS)
            if observed != [{column: str(row[column]) for column in MANIFEST_COLUMNS} for row in expected]:
                raise ContractError("final integration code freeze drift")
            _, ready_rows = read_tsv(ready, ("status", "manifest_sha256", "n_files", "real_outcomes_read_by_freezer"))
            if len(ready_rows) != 1 or ready_rows[0]["status"] != "outcome_independent_final_assembly_code_frozen":
                raise ContractError("final integration code READY drift")
            if ready_rows[0]["manifest_sha256"] != sha256_file(manifest) or ready_rows[0]["n_files"] != str(len(expected)):
                raise ContractError("final integration code READY hash/count drift")
            if ready_rows[0]["real_outcomes_read_by_freezer"] != "FALSE":
                raise ContractError("code freezer claims it read real outcomes")
            print(f"PASS: {len(expected)} final-integration code hashes unchanged")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
