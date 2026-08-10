#!/usr/bin/env python3
"""Independently validate the bounded Plan 13 real-adapter seam."""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    Check,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    parse_bool,
    parse_float,
    parse_int,
    program_universe,
    read_tsv,
    sha256_file,
    validate_candidate_contract,
    validate_hotspot_seal,
    write_tsv,
)
from importlib.util import module_from_spec, spec_from_file_location


REPORT_COLUMNS = ("scope", "check_id", "status", "detail")
BUNDLE_COLUMNS = ("role", "relative_path", "bytes", "sha256")
STATUS_COLUMNS = (
    "release_id",
    "status",
    "map_selection_frozen",
    "gse287826_status",
    "govaere_geomx_status",
    "govaere_cosmx_status",
    "real_adapter_registry_sha256",
    "map_selection_manifest_sha256",
    "bundle_manifest_sha256",
    "validation_report_sha256",
    "yakubovsky_imported",
    "visium_imported",
    "proteomics_imported",
    "atac_imported",
    "plan13_complete",
    "validated_at_utc",
)
READY_COLUMNS = (
    "release_id",
    "status",
    "real_adapter_registry_sha256",
    "map_selection_manifest_sha256",
    "bundle_manifest_sha256",
    "validation_report_sha256",
    "partial_status_sha256",
    "builder_sha256",
    "validator_sha256",
    "n_real_adapters",
    "n_source_native_context_rows",
    "yakubovsky_imported",
    "visium_imported",
    "proteomics_imported",
    "atac_imported",
    "plan13_complete",
    "validated_at_utc",
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_module(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ContractError(f"cannot load validator dependency: {path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rows_by_uid(path: Path) -> dict[str, dict[str, str]]:
    _, rows = read_tsv(path, ("program_uid",))
    result = {row["program_uid"]: row for row in rows}
    if len(result) != len(rows):
        raise ContractError(f"duplicate program_uid in {path}")
    return result


def _validate_map_selection(project_root: Path, candidate_root: Path, hotspot_root: Path) -> list[Check]:
    freezer = _load_module(
        Path(__file__).resolve().parent / "07_freeze_map_selection.py",
        "spatial_map_freezer",
    )
    freezer.check_existing(project_root, candidate_root, hotspot_root)
    _, rows = read_tsv(candidate_root / "map_selection_manifest.tsv", freezer.MAP_COLUMNS)
    if any(parse_bool(row["external_outcomes_read"], "map.external_outcomes_read") for row in rows):
        raise ContractError("map selection claims external outcomes were read")
    if any(
        parse_bool(
            row["yakubovsky_program_outcomes_present_at_freeze"],
            "map.yakubovsky_program_outcomes_present_at_freeze",
        )
        for row in rows
    ):
        raise ContractError("map selection was not frozen before Yakubovsky program outcomes")
    return [
        Check(
            "map_selection_outcome_blind",
            "pass",
            f"{len(rows)} complete sealed external-test programs selected from Figure 2 only",
        )
    ]


def _validate_gse(project_root: Path, candidate_root: Path) -> list[Check]:
    adapter = candidate_root / "real_adapters/gse287826_skip"
    acquisition_gate = (
        project_root
        / "Analysis/Spatial/candidates"
        / RELEASE_ID
        / "acquisition/gse287826/gate_status.tsv"
    )
    _, gates = read_tsv(
        acquisition_gate,
        ("status", "n_biological_aois", "authoritative_donor_key_pass", "inference_authorized"),
    )
    if len(gates) != 1:
        raise ContractError("GSE287826 public source gate must contain one row")
    gate = gates[0]
    if (
        gate["status"] != "skipped_no_donor_key"
        or gate["n_biological_aois"] != "30"
        or gate["authoritative_donor_key_pass"].strip().lower() != "false"
        or gate["inference_authorized"].strip().lower() != "false"
    ):
        raise ContractError("GSE287826 skip facts no longer rederive")
    _, effects = read_tsv(adapter / "program_effects.tsv", SCHEMAS["program_effects.tsv"])
    _, scores = read_tsv(
        adapter / "per_sample_program_scores.tsv",
        SCHEMAS["per_sample_program_scores.tsv"],
    )
    if scores:
        raise ContractError("GSE287826 skip adapter contains program scores")
    for row in effects:
        if row["evidence_state"] != "skipped" or row["testable"] != "FALSE":
            raise ContractError("GSE287826 skip program row is not terminal skipped")
        if any(row[column].strip() for column in ("estimate", "std_error", "pvalue", "padj")):
            raise ContractError("GSE287826 skip program row contains an effect")
    return [
        Check(
            "gse287826_skip",
            "pass",
            f"30 AOIs; no public donor key; {len(effects)} v2 programs skipped without outcomes",
        )
    ]


def _validate_geomx(project_root: Path, candidate_root: Path) -> list[Check]:
    adapter = candidate_root / "real_adapters/govaere2026_geomx_dropped"
    decision = project_root / "Analysis/Spatial/results/govaere2026/GEOMX_DROPPED.md"
    text = decision.read_text(encoding="utf-8")
    if "NOT used as an evidence arm anywhere in the paper" not in text:
        raise ContractError("Govaere GeoMx authoritative exclusion text drifted")
    _, sources = read_tsv(adapter / "source_manifest.tsv", SCHEMAS["source_manifest.tsv"])
    prohibited_numeric_sources = {
        "geomx_de_sh_vs_pt.csv",
        "geomx_de_sh_vs_ls.csv",
        "geomx_de_sh_vs_pt_limma.csv",
        "geomx_de_sh_vs_ls_limma.csv",
    }
    if any(Path(row["relative_path"]).name in prohibited_numeric_sources for row in sources):
        raise ContractError("dropped Govaere GeoMx numeric outcomes were imported")
    _, effects = read_tsv(adapter / "program_effects.tsv", SCHEMAS["program_effects.tsv"])
    if any(row["evidence_state"] != "skipped" for row in effects):
        raise ContractError("Govaere GeoMx dropped adapter contains a non-skipped program row")
    return [
        Check(
            "govaere_geomx_exclusion",
            "pass",
            f"repository exclusion preserved; {len(effects)} program rows skipped; no numeric GeoMx result imported",
        )
    ]


def _validate_cosmx(project_root: Path, candidate_root: Path, hotspot_root: Path) -> list[Check]:
    adapter = candidate_root / "real_adapters/govaere2026_cosmx_il32_context"
    source_path = project_root / "Analysis/Spatial/results/govaere2026/il32_colocalization/il32_coloc_per_slide.csv"
    builder_path = Path(__file__).resolve().parent / "08_build_real_adapters.py"
    builder = _load_module(builder_path, "spatial_real_adapter_builder")
    source_rows = builder._read_csv(source_path)
    expected_native = builder.source_native_cosmx_rows(source_rows)
    _, observed_native = read_tsv(adapter / "source_native_context.tsv", builder.SOURCE_NATIVE_COLUMNS)
    normalized_expected = [
        {
            column: (
                "TRUE" if value is True else "FALSE" if value is False else str(value)
            )
            for column, value in ((column, row[column]) for column in builder.SOURCE_NATIVE_COLUMNS)
        }
        for row in expected_native
    ]
    if observed_native != normalized_expected:
        raise ContractError("CosMx source_native_context.tsv does not exactly rederive from per-array input")
    included = [row for row in observed_native if row["include_direction_concordance"] == "TRUE"]
    if len(included) != 3 or any(row["direction_concordant"] != "TRUE" for row in included):
        raise ContractError("CosMx accepted source context does not retain 3/3 MASH-array concordance")

    _, effects = read_tsv(adapter / "program_effects.tsv", SCHEMAS["program_effects.tsv"])
    _, scores = read_tsv(
        adapter / "per_sample_program_scores.tsv",
        SCHEMAS["per_sample_program_scores.tsv"],
    )
    if scores:
        raise ContractError("CosMx source-native adapter invented whole-program scores")
    if any(row["evidence_state"] != "not_applicable" or row["testable"] != "FALSE" for row in effects):
        raise ContractError("CosMx whole-program rows must remain not_applicable")
    if any(any(row[column].strip() for column in ("estimate", "std_error", "pvalue", "padj")) for row in effects):
        raise ContractError("CosMx whole-program rows contain an invented numeric effect")

    seal = validate_hotspot_seal(hotspot_root)
    universe = program_universe(seal, "external_test_eligible")
    membership_path = hotspot_root / "program_membership_v2.tsv"
    _, membership_rows = read_tsv(
        membership_path,
        ("program_uid", "mapped_symbol", "original_l1_weight", "membership_sha256"),
    )
    panel_path = project_root / "Analysis/Spatial/results/govaere2026/cosmx_de_Hepatocyte_MASH_vs_noMASH.csv"
    panel_rows = builder._read_csv(panel_path)
    panel = {row["gene"].strip() for row in panel_rows if row["gene"].strip()}
    expected_coverage = builder.measured_program_coverage(membership_rows, panel, universe)
    mapping = _rows_by_uid(adapter / "gene_mapping_audit.tsv")
    for uid, (n_genes, retained, _status) in expected_coverage.items():
        if parse_int(mapping[uid]["n_genes_measured"], f"CosMx.n_genes[{uid}]") != n_genes:
            raise ContractError(f"CosMx coverage gene count does not rederive for {uid}")
        observed_weight = parse_float(mapping[uid]["retained_l1_weight"], f"CosMx.weight[{uid}]")
        if observed_weight is None or abs(observed_weight - retained) > 1e-12:
            raise ContractError(f"CosMx retained L1 weight does not rederive for {uid}")
    return [
        Check(
            "govaere_cosmx_source_native",
            "pass",
            (
                f"{len(observed_native)} physical arrays preserved; 3/3 MASH directions; "
                f"{len(effects)} v2 program rows explicitly not_applicable; no program refit"
            ),
        )
    ]


def _validate_v1_preservation(candidate_root: Path) -> list[Check]:
    _, rows = read_tsv(
        candidate_root / "v1_preservation_summary.tsv",
        ("n_unchanged", "n_changed", "n_missing", "n_added", "aggregate_sha256"),
    )
    if len(rows) != 1:
        raise ContractError("v1 preservation summary must contain one row")
    row = rows[0]
    if row["n_unchanged"] != "591" or any(
        row[column] != "0" for column in ("n_changed", "n_missing", "n_added")
    ):
        raise ContractError("v1 preservation failed during real-adapter work")
    return [
        Check(
            "v1_preservation",
            "pass",
            f"591/591 unchanged; aggregate={row['aggregate_sha256']}",
        )
    ]


def run_checks(project_root: Path, candidate_root: Path, hotspot_root: Path) -> list[Check]:
    checks = _validate_map_selection(project_root, candidate_root, hotspot_root)
    checks.extend(
        validate_candidate_contract(
            candidate_root,
            candidate_root / "real_adapter_registry.tsv",
            hotspot_root,
            project_root,
        )
    )
    _, registry = read_tsv(candidate_root / "real_adapter_registry.tsv", ADAPTER_REGISTRY_COLUMNS)
    expected_ids = {
        "gse287826_skip",
        "govaere2026_geomx_dropped",
        "govaere2026_cosmx_il32_context",
    }
    if {row["adapter_id"] for row in registry} != expected_ids or len(registry) != 3:
        raise ContractError("real adapter registry is not the exact bounded three-adapter seam")
    checks.append(Check("bounded_real_registry", "pass", "exactly three prespecified real adapters"))
    checks.extend(_validate_gse(project_root, candidate_root))
    checks.extend(_validate_geomx(project_root, candidate_root))
    checks.extend(_validate_cosmx(project_root, candidate_root, hotspot_root))
    checks.extend(_validate_v1_preservation(candidate_root))
    return checks


def _bundle_manifest(candidate_root: Path) -> list[dict[str, object]]:
    paths = [candidate_root / "map_selection_manifest.tsv", candidate_root / "real_adapter_registry.tsv"]
    paths.extend(sorted((candidate_root / "real_adapters").rglob("*")))
    rows = []
    for path in paths:
        if not path.is_file():
            continue
        rows.append(
            {
                "role": "real_adapter_bundle",
                "relative_path": path.relative_to(candidate_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    if len(rows) < 2 + 3 * len(SCHEMAS):
        raise ContractError("real-adapter bundle manifest is unexpectedly small")
    return rows


def _write_outputs(project_root: Path, candidate_root: Path, checks: list[Check]) -> None:
    outputs = (
        candidate_root / "real_adapter_bundle_manifest.tsv",
        candidate_root / "real_adapter_validation.tsv",
        candidate_root / "sp_int_04_safe_seam_status.tsv",
        candidate_root / "REAL_ADAPTERS_READY",
    )
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise ContractError("refusing to overwrite real-adapter validation outputs: " + ", ".join(existing))

    bundle_rows = _bundle_manifest(candidate_root)
    write_tsv(outputs[0], BUNDLE_COLUMNS, bundle_rows)
    report_rows = [
        {"scope": "SP-INT-04-safe-seam", "check_id": row.check_id, "status": row.status, "detail": row.detail}
        for row in checks
    ]
    write_tsv(outputs[1], REPORT_COLUMNS, report_rows)
    validated = utc_now()
    registry_hash = sha256_file(candidate_root / "real_adapter_registry.tsv")
    map_hash = sha256_file(candidate_root / "map_selection_manifest.tsv")
    bundle_hash = sha256_file(outputs[0])
    report_hash = sha256_file(outputs[1])
    status_row = {
        "release_id": RELEASE_ID,
        "status": "map_gse_govaere_ready_partial_real_integration",
        "map_selection_frozen": True,
        "gse287826_status": "skipped_no_donor_key",
        "govaere_geomx_status": "skipped_repository_exclusion",
        "govaere_cosmx_status": "source_native_context_no_program_refit",
        "real_adapter_registry_sha256": registry_hash,
        "map_selection_manifest_sha256": map_hash,
        "bundle_manifest_sha256": bundle_hash,
        "validation_report_sha256": report_hash,
        "yakubovsky_imported": False,
        "visium_imported": False,
        "proteomics_imported": False,
        "atac_imported": False,
        "plan13_complete": False,
        "validated_at_utc": validated,
    }
    write_tsv(outputs[2], STATUS_COLUMNS, [status_row])
    builder_path = Path(__file__).resolve().parent / "08_build_real_adapters.py"
    ready_row = {
        "release_id": RELEASE_ID,
        "status": "ready_partial_real_adapter_seam_not_plan13_complete",
        "real_adapter_registry_sha256": registry_hash,
        "map_selection_manifest_sha256": map_hash,
        "bundle_manifest_sha256": bundle_hash,
        "validation_report_sha256": report_hash,
        "partial_status_sha256": sha256_file(outputs[2]),
        "builder_sha256": sha256_file(builder_path),
        "validator_sha256": sha256_file(Path(__file__).resolve()),
        "n_real_adapters": 3,
        "n_source_native_context_rows": 4,
        "yakubovsky_imported": False,
        "visium_imported": False,
        "proteomics_imported": False,
        "atac_imported": False,
        "plan13_complete": False,
        "validated_at_utc": validated,
    }
    write_tsv(outputs[3], READY_COLUMNS, [ready_row])


def _check_ready(candidate_root: Path) -> None:
    ready_path = candidate_root / "REAL_ADAPTERS_READY"
    _, rows = read_tsv(ready_path, READY_COLUMNS)
    if len(rows) != 1:
        raise ContractError("REAL_ADAPTERS_READY must contain one row")
    row = rows[0]
    expected = {
        "real_adapter_registry_sha256": sha256_file(candidate_root / "real_adapter_registry.tsv"),
        "map_selection_manifest_sha256": sha256_file(candidate_root / "map_selection_manifest.tsv"),
        "bundle_manifest_sha256": sha256_file(candidate_root / "real_adapter_bundle_manifest.tsv"),
        "validation_report_sha256": sha256_file(candidate_root / "real_adapter_validation.tsv"),
        "partial_status_sha256": sha256_file(candidate_root / "sp_int_04_safe_seam_status.tsv"),
        "builder_sha256": sha256_file(Path(__file__).resolve().parent / "08_build_real_adapters.py"),
        "validator_sha256": sha256_file(Path(__file__).resolve()),
    }
    if row["release_id"] != RELEASE_ID or row["status"] != "ready_partial_real_adapter_seam_not_plan13_complete":
        raise ContractError("REAL_ADAPTERS_READY status drift")
    for column, value in expected.items():
        if row[column] != value:
            raise ContractError(f"REAL_ADAPTERS_READY hash mismatch for {column}")
    if any(parse_bool(row[column], f"READY.{column}") for column in (
        "yakubovsky_imported", "visium_imported", "proteomics_imported", "atac_imported", "plan13_complete"
    )):
        raise ContractError("partial real-adapter READY overstates imported lanes or Plan 13 completion")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    )
    try:
        checks = run_checks(project_root, candidate_root, hotspot_root)
        if any(row.status != "pass" for row in checks):
            raise ContractError("one or more real-adapter checks did not pass")
        if args.check_only:
            _check_ready(candidate_root)
            print(f"PASS: rederived {len(checks)} real-adapter checks and READY hashes")
            return 0
        _write_outputs(project_root, candidate_root, checks)
        _check_ready(candidate_root)
        print(f"PASS: {len(checks)} real-adapter checks; partial seam READY")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
