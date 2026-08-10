#!/usr/bin/env python3
"""Validate and seal the isolated evidence-state-v2 correction release."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from contract_lib import (
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    SEMANTIC_CONTRACT_ID,
    parse_bool,
    parse_float,
    read_tsv,
    sha256_file,
    write_tsv,
)


MANIFEST_COLUMNS = ("path_scope", "relative_path", "bytes", "sha256", "role")


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def _one_row(path: Path) -> dict[str, str]:
    _, rows = read_tsv(path)
    if len(rows) != 1:
        raise ContractError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def _verify_historical_baseline(target_root: Path) -> None:
    status = _one_row(target_root / "semantic_v2_preparation_status.tsv")
    roots = {
        "historical_plan13_root": Path(status["source_root"]),
        "historical_yakubovsky_root": Path(status["yak_native_root"]),
    }
    _, rows = read_tsv(
        target_root / "semantic_v2_source_baseline.tsv",
        ("source_scope", "relative_path", "bytes", "sha256"),
    )
    if not rows:
        raise ContractError("semantic-v2 historical baseline is empty")
    for row in rows:
        if row["source_scope"] not in roots:
            raise ContractError(f"unknown historical source scope: {row['source_scope']}")
        path = roots[row["source_scope"]] / row["relative_path"]
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256_file(path) != row["sha256"]
        ):
            raise ContractError(f"historical bundle drift after semantic remediation: {path}")


def _validate_semantics(target_root: Path) -> tuple[Counter[str], str]:
    final_root = target_root / "final_integration"
    final_ready = _one_row(final_root / "READY")
    if final_ready.get("status") != "ready_plan13_candidate_source_tables_no_pdf_no_canonical_write":
        raise ContractError("semantic-v2 final integration READY status drift")
    if parse_bool(final_ready["canonical_promotion_authorized"], "final.canonical_promotion"):
        raise ContractError("semantic-v2 final integration authorizes canonical promotion")
    _, rows = read_tsv(final_root / "integrated_program_effects.tsv", SCHEMAS["program_effects.tsv"])
    states = Counter(row["evidence_state"] for row in rows)
    for row in rows:
        if row["evidence_state"] == "tested_negative" and not row["negative_call_rule_id"].strip():
            raise ContractError("tested_negative row lacks an explicit adequate-negative rule")
        if row["evidence_state"] != "tested_negative" and row["negative_call_rule_id"].strip():
            raise ContractError("non-negative row carries a negative-call rule")

    yak = [row for row in rows if row["dataset"] == "Yakubovsky_2026"]
    if len(yak) != 2 or any(row["evidence_state"] != "indeterminate" for row in yak):
        raise ContractError("Yakubovsky non-robust effects are not both indeterminate")
    if not any(not parse_bool(row["direction_agreement"], "yak.direction_agreement") for row in yak):
        raise ContractError("Yakubovsky descriptive/inferential direction disagreement was hidden")
    for row in yak:
        if not all(
            row[column].strip()
            for column in (
                "descriptive_effect_direction",
                "inferential_test_direction",
                "heterogeneity_statistic",
                "heterogeneity_df",
                "heterogeneity_pvalue",
            )
        ):
            raise ContractError("Yakubovsky direction/heterogeneity audit is incomplete")

    atac = [row for row in rows if row["dataset"] in {"GSE281367", "GSE244832"}]
    if len(atac) != 4 or any(row["evidence_state"] != "indeterminate" for row in atac):
        raise ContractError("dynamic ATAC non-robust effects are not all indeterminate")

    native = [row for row in rows if row["dataset"] in {"GSE192741", "Vu_et_al_2025"}]
    if len(native) != 4:
        raise ContractError("native spatial semantic-v2 grid is incomplete")
    for row in native:
        matched_sd = parse_float(row["matched_null_sd"], "native.matched_null_sd")
        if matched_sd is None or matched_sd <= 0 or row["std_error"].strip():
            raise ContractError("matched-null dispersion is mislabeled as sampling uncertainty")

    unresolved = [row for row in rows if row["dataset"] in {"Vu_et_al_2025", "Govaere2026_CosMx"}]
    if not unresolved or any(
        row["biological_unit_resolution"] != "unresolved"
        or row["biological_unit"] != "unknown_public_biological_unit"
        or row["n_biological"].strip()
        or not row["n_technical"].strip()
        for row in unresolved
    ):
        raise ContractError("unresolved arrays were reported as biological replicates")
    return states, sha256_file(final_root / "READY")


def _manifest_rows(target_root: Path, yak_adapter_root: Path) -> list[dict[str, object]]:
    excluded = {"semantic_v2_release_manifest.tsv", "SEMANTIC_V2_READY"}
    rows: list[dict[str, object]] = []
    for scope, root, role in (
        ("semantic_plan13_root", target_root, "semantic_v2_candidate_artifact"),
        ("semantic_yak_adapter", yak_adapter_root, "semantic_v2_yak_adapter_artifact"),
    ):
        for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
            if path.name in excluded:
                continue
            rows.append(
                {
                    "path_scope": scope,
                    "relative_path": path.relative_to(root).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "role": role,
                }
            )
    return rows


def seal(target_root: Path, yak_adapter_root: Path) -> None:
    manifest = target_root / "semantic_v2_release_manifest.tsv"
    ready = target_root / "SEMANTIC_V2_READY"
    if manifest.exists() or ready.exists():
        raise ContractError("refusing to overwrite semantic-v2 release seal")
    _verify_historical_baseline(target_root)
    states, final_ready_sha = _validate_semantics(target_root)
    yak_ready = _one_row(yak_adapter_root / "PLAN13_ADAPTER_READY")
    if parse_bool(yak_ready["canonical_promotion_authorized"], "yak.canonical_promotion"):
        raise ContractError("semantic-v2 Yakubovsky adapter authorizes canonical promotion")
    rows = _manifest_rows(target_root, yak_adapter_root)
    write_tsv(manifest, MANIFEST_COLUMNS, rows)
    write_tsv(
        ready,
        (
            "release_id",
            "semantic_contract_id",
            "status",
            "manifest_sha256",
            "final_ready_sha256",
            "yak_adapter_ready_sha256",
            "n_files",
            "n_robust",
            "n_indeterminate",
            "n_tested_negative",
            "historical_bundles_unchanged",
            "canonical_promotion_authorized",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "semantic_contract_id": SEMANTIC_CONTRACT_ID,
                "status": "sealed_semantic_v2_candidate_no_canonical_promotion",
                "manifest_sha256": sha256_file(manifest),
                "final_ready_sha256": final_ready_sha,
                "yak_adapter_ready_sha256": sha256_file(yak_adapter_root / "PLAN13_ADAPTER_READY"),
                "n_files": len(rows),
                "n_robust": states["robust"],
                "n_indeterminate": states["indeterminate"],
                "n_tested_negative": states["tested_negative"],
                "historical_bundles_unchanged": True,
                "canonical_promotion_authorized": False,
            }
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--target-root", type=Path, default=None)
    parser.add_argument("--yak-adapter-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    release_root = (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
    )
    target_root = (
        args.target_root or (release_root / "spatial_context_semantic_v2_2026-08-08")
    ).resolve()
    yak_adapter_root = (
        args.yak_adapter_root
        or (
            project_root
            / "Analysis/Spatial/candidates"
            / RELEASE_ID
            / "semantic_v2_2026-08-08/yakubovsky_plan13_adapter"
        )
    ).resolve()
    try:
        seal(target_root, yak_adapter_root)
        print(f"PASS: sealed semantic-v2 candidate at {target_root / 'SEMANTIC_V2_READY'}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
