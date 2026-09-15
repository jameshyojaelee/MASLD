#!/usr/bin/env python3
"""Audit the frozen GSE281367 ATAC transport object without opening outcomes."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping


class ProductionReadinessError(RuntimeError):
    """Raised when a read-only authority or fail-closed decision changes."""


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ProductionReadinessError(f"JSON authority is not an object: {path}")
    return value


def verify_file(path: Path, expected_sha256: str) -> None:
    if sha256_path(path) != expected_sha256:
        raise ProductionReadinessError(f"SHA-256 differs: {path}")


def verify_manifest_authority(
    repo_root: Path,
    spec: Mapping[str, object],
) -> tuple[Path, dict[str, Any]]:
    root = (repo_root / str(spec["path"])).resolve(strict=True)
    manifest_path = root / "ARTIFACTS.json"
    verify_file(manifest_path, str(spec["artifacts_sha256"]))
    manifest = load_json(manifest_path)
    metadata = manifest.get("metadata")
    if not isinstance(metadata, dict) or metadata.get("artifact_class") != spec["artifact_class"]:
        raise ProductionReadinessError(f"artifact class differs: {root}")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ProductionReadinessError(f"artifact inventory missing: {root}")
    paths = [row.get("path") for row in artifacts if isinstance(row, dict)]
    if len(paths) != len(artifacts) or len(set(paths)) != len(paths):
        raise ProductionReadinessError(f"artifact inventory is malformed: {root}")
    return root, manifest


def validate_contract(contract: Mapping[str, object]) -> None:
    expected_schema = "masld-bench-gse281367-atac-transport-production-readiness-v1"
    if contract.get("schema_version") != expected_schema:
        raise ProductionReadinessError("schema version differs")
    for key in (
        "automatic_activation",
        "sealed_outcomes_read",
        "biological_outcome_arrays_read",
        "condition_values_read",
        "condition_values_used",
        "source_native_input_activated",
        "external_or_sealed_evaluation",
        "champion_claim_allowed",
    ):
        if contract.get(key) is not False:
            raise ProductionReadinessError(f"{key} must remain false")
    expected_axis = contract.get("expected_axis")
    if not isinstance(expected_axis, dict):
        raise ProductionReadinessError("axis contract missing")
    if expected_axis.get("condition_column_present") is not True:
        raise ProductionReadinessError("label-bearing evaluator axis must remain explicit")
    if expected_axis.get("label_free_model_axis_frozen") is not False:
        raise ProductionReadinessError("label-free model axis is not yet frozen")
    if expected_axis.get("dataset_donor_n") != 12:
        raise ProductionReadinessError("GSE281367 donor count differs")
    if expected_axis.get("t_cell_state") != "not_applicable_in_frozen_transport_axis":
        raise ProductionReadinessError("missing T-cell state must remain not applicable")
    firewall = contract.get("production_firewall")
    if not isinstance(firewall, dict) or any(value is not False for value in firewall.values()):
        raise ProductionReadinessError("unresolved production firewalls cannot be relaxed")
    decision = contract.get("decision")
    if not isinstance(decision, dict):
        raise ProductionReadinessError("decision block missing")
    for key in ("frozen_outcome_authority_valid", "donor_grouping_contract_valid"):
        if decision.get(key) is not True:
            raise ProductionReadinessError(f"{key} must remain true")
    for key in (
        "label_agnostic_production_ready",
        "observed_atac_production_ready",
        "family_native_tournament_ready",
        "source_native_model_input_ready",
        "existing_outcome_artifact_may_be_rebuilt",
        "existing_outcome_artifact_may_be_given_to_model_jobs",
    ):
        if decision.get(key) is not False:
            raise ProductionReadinessError(f"decision {key} must remain false")
    if "gse289173" in json.dumps(contract).lower():
        raise ProductionReadinessError("sealed accession must not enter this audit")


def audit_axis_header(path: Path, expected_fields: list[str]) -> dict[str, object]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        first_line = handle.readline()
    fields = next(csv.reader([first_line], delimiter="\t"))
    if fields != expected_fields:
        raise ProductionReadinessError("donor-lineage axis fields differ")
    return {
        "axis_sha256_verified": True,
        "header_only_read": True,
        "data_rows_read": 0,
        "condition_column_present": "condition" in fields,
        "condition_values_read": False,
        "condition_values_used": False,
        "label_free_model_axis_frozen": False,
        "production_exchange_label_safe": False,
        "reason": (
            "the sole frozen row axis physically co-locates condition with "
            "model/evaluator row identifiers"
        ),
    }


def audit_windows(
    path: Path,
    fold_roles_path: Path,
    expected: Mapping[str, object],
) -> dict[str, object]:
    fold_roles = load_json(fold_roles_path)
    role_by_contig: dict[str, str] = {}
    for role in ("train", "valid", "test"):
        contigs = fold_roles.get(role)
        if not isinstance(contigs, list) or not contigs:
            raise ProductionReadinessError(f"genomic role missing: {role}")
        for contig in contigs:
            if contig in role_by_contig:
                raise ProductionReadinessError(f"contig occurs in two roles: {contig}")
            role_by_contig[str(contig)] = role
    counts: Counter[str] = Counter()
    seen_ids: set[str] = set()
    seen_indexes: set[int] = set()
    role_contigs: dict[str, set[str]] = {"valid": set(), "test": set()}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["window_index", "window_id", "role", "contig", "start", "end"]
        if list(reader.fieldnames or ()) != required:
            raise ProductionReadinessError("window fields differ")
        for row in reader:
            index = int(row["window_index"])
            start, end = int(row["start"]), int(row["end"])
            role, contig = row["role"], row["contig"]
            if role not in {"valid", "test"} or role_by_contig.get(contig) != role:
                raise ProductionReadinessError(f"window role differs at index {index}")
            if start < 0 or end - start != int(expected["window_width_bp"]):
                raise ProductionReadinessError(f"window geometry differs at index {index}")
            if index in seen_indexes or row["window_id"] in seen_ids:
                raise ProductionReadinessError(f"duplicate window at index {index}")
            seen_indexes.add(index)
            seen_ids.add(row["window_id"])
            counts[role] += 1
            role_contigs[role].add(contig)
    expected_n = int(expected["window_n"])
    if seen_indexes != set(range(expected_n)):
        raise ProductionReadinessError("window indexes are not a complete zero-based axis")
    expected_counts = {
        "valid": int(expected["valid_window_n"]),
        "test": int(expected["test_window_n"]),
    }
    if counts != expected_counts:
        raise ProductionReadinessError(f"window role counts differ: {dict(counts)}")
    if role_contigs["valid"] & role_contigs["test"]:
        raise ProductionReadinessError("valid/test contig roles overlap")
    return {
        "window_sha256_verified": True,
        "window_n": expected_n,
        "window_width_bp": int(expected["window_width_bp"]),
        "profile_bins": int(expected["profile_bins"]),
        "profile_bin_width_bp": int(expected["profile_bin_width_bp"]),
        "role_counts": dict(sorted(counts.items())),
        "train_windows_present": False,
        "valid_test_contigs_disjoint": True,
        "coordinates_read": True,
        "biological_outcome_arrays_read": False,
    }


def audit_source_locks(repo_root: Path, contract: Mapping[str, object]) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    groups = (
        (
            "builder_authority",
            (("path", "sha256"), ("test_path", "test_sha256"), ("sbatch_path", "sbatch_sha256")),
        ),
        ("membership_authority", (("builder_path", "builder_sha256"),)),
        ("task_authority", (("task_path", "task_sha256"), ("gate_path", "gate_sha256"))),
    )
    for section, pairs in groups:
        spec = contract[section]
        if not isinstance(spec, dict):
            raise ProductionReadinessError(f"source lock missing: {section}")
        for path_key, sha_key in pairs:
            path = repo_root / str(spec[path_key])
            verify_file(path, str(spec[sha_key]))
            rows.append({"path": str(path.relative_to(repo_root)), "sha256": spec[sha_key]})
    return {"all_current_sources_match_frozen_hashes": True, "sources": rows}


def audit_authorities(repo_root: Path, contract: Mapping[str, object]) -> dict[str, object]:
    roots: dict[str, Path] = {}
    for name in (
        "activation_authority",
        "outcome_authority",
        "membership_authority",
        "window_authority",
    ):
        spec = contract[name]
        if not isinstance(spec, dict):
            raise ProductionReadinessError(f"authority missing: {name}")
        roots[name], _manifest = verify_manifest_authority(repo_root, spec)
    outcome_spec = contract["outcome_authority"]
    window_spec = contract["window_authority"]
    assert isinstance(outcome_spec, dict) and isinstance(window_spec, dict)
    files = (
        (roots["outcome_authority"], "summary_path", "summary_sha256", outcome_spec),
        (roots["outcome_authority"], "axis_path", "axis_sha256", outcome_spec),
        (roots["outcome_authority"], "windows_path", "windows_sha256", outcome_spec),
        (roots["window_authority"], "summary_path", "summary_sha256", window_spec),
        (roots["window_authority"], "fold_roles_path", "fold_roles_sha256", window_spec),
    )
    for root, relative_key, sha_key, spec in files:
        verify_file(root / str(spec[relative_key]), str(spec[sha_key]))
    outcome = load_json(roots["outcome_authority"] / str(outcome_spec["summary_path"]))
    if outcome.get("artifact_role") != "evaluator_only_development_outcomes":
        raise ProductionReadinessError("outcome artifact role differs")
    if outcome.get("model_inputs_may_read_evaluator_outcomes") is not False:
        raise ProductionReadinessError("model-input firewall changed")
    if outcome.get("same_cell_pairing_claimed") is not False:
        raise ProductionReadinessError("same-cell pairing claim changed")
    if outcome.get("donors", {}).get("gse281367") != 12:
        raise ProductionReadinessError("outcome GSE281367 donor count differs")
    region_summary = load_json(roots["window_authority"] / str(window_spec["summary_path"]))
    if region_summary.get("held_atac_used_for_region_selection") is not False:
        raise ProductionReadinessError("held ATAC entered region selection")
    for key in ("test_region_source", "validation_region_source"):
        if region_summary.get(key) != "fixed_ENCODE_v4_cCRE_contract":
            raise ProductionReadinessError(f"{key} differs")
    return {
        "roots": roots,
        "outcome_summary": outcome,
        "region_summary": region_summary,
        "artifact_manifest_sha256_verified": True,
        "metadata_artifact_sha256_verified": True,
        "outcome_array_files_opened": 0,
        "outcome_array_files_rehashed": 0,
    }


def run_audit(repo_root: Path, contract: Mapping[str, object]) -> dict[str, object]:
    validate_contract(contract)
    authorities = audit_authorities(repo_root, contract)
    roots = authorities.pop("roots")
    outcome_spec = contract["outcome_authority"]
    window_spec = contract["window_authority"]
    expected_axis = contract["expected_axis"]
    expected_windows = contract["expected_windows"]
    assert isinstance(outcome_spec, dict)
    assert isinstance(window_spec, dict)
    assert isinstance(expected_axis, dict)
    assert isinstance(expected_windows, dict)
    axis = audit_axis_header(
        roots["outcome_authority"] / str(outcome_spec["axis_path"]),
        list(expected_axis["fields"]),
    )
    windows = audit_windows(
        roots["outcome_authority"] / str(outcome_spec["windows_path"]),
        roots["window_authority"] / str(window_spec["fold_roles_path"]),
        expected_windows,
    )
    source_locks = audit_source_locks(repo_root, contract)
    membership = contract["membership_authority"]
    assert isinstance(membership, dict)
    donor_safety = {
        "biological_unit": "donor",
        "dataset_donor_n": 12,
        "same_cell_pairing_claimed": False,
        "rna_state": "structurally_missing",
        "outer_fold_count": int(membership["outer_fold_count"]),
        "outer_fold_seed": int(membership["outer_fold_seed"]),
        "outer_fold_inputs": list(membership["outer_fold_inputs"]),
        "outer_fold_uses_condition": False,
        "produced_axis_data_rows_read": 0,
        "proof_basis": (
            "exact frozen membership/outcome builder hashes and successful immutable "
            "build; condition values were not reopened"
        ),
    }
    firewall = dict(contract["production_firewall"])
    blockers = [key for key, value in firewall.items() if value is False]
    decision = dict(contract["decision"])
    if not blockers or decision["label_agnostic_production_ready"] is not False:
        raise ProductionReadinessError("production decision no longer fails closed")
    return {
        "artifact_integrity": authorities,
        "donor_safety": donor_safety,
        "label_firewall": axis,
        "genomic_window_contract": windows,
        "source_locks": source_locks,
        "production_firewall": firewall,
        "unresolved_blockers": blockers,
        "decision": decision,
        "read_policy": {
            "sealed_outcomes_read": False,
            "biological_outcome_arrays_read": False,
            "condition_values_read": False,
            "condition_values_used": False,
            "window_coordinates_read": True,
            "source_and_manifest_metadata_read": True,
        },
    }


def write_outputs(output: Path, audit: Mapping[str, object]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    for name in (
        "artifact_integrity",
        "donor_safety",
        "label_firewall",
        "genomic_window_contract",
        "source_locks",
        "production_firewall",
        "read_policy",
    ):
        (output / f"{name}.json").write_bytes(canonical_json(audit[name]))
    decision = dict(audit["decision"])
    decision["unresolved_blockers"] = audit["unresolved_blockers"]
    (output / "production_decision.json").write_bytes(canonical_json(decision))
    (output / "DECISION.md").write_text(
        "# GSE281367 ATAC transport production-readiness decision\n\n"
        "**Decision: retain the immutable 32,000-window evaluator outcome authority, "
        "but do not expose it to model jobs or start production scoring.** The donor "
        "unit, 12-donor census, ATAC-only topology, exact source locks, fixed ENCODE "
        "cCRE validation/test windows, and held-chromosome roles remain valid. No "
        "biological outcome array, condition value, controlled data, or sealed outcome "
        "was opened by this audit.\n\n"
        "The current row-axis file physically includes `condition`. There is no frozen "
        "label-free prediction axis, condition-blind scorer, prediction-hash commit, "
        "or family-native evaluator command. Observed-ATAC models also lack a "
        "model-specific target-window mask and receptive-field buffer, so giving them "
        "the current query accessibility could leak the scored signal. The frozen "
        "axis has four lineages; T cell is not applicable here and must never be "
        "encoded as failure or zero.\n\n"
        "Production activation therefore waits for a new derived, immutable label-free "
        "exchange layer plus the prospective global census revision already required "
        "by the observed-multiome TaskSpec. The original outcome object is not rebuilt "
        "or modified. It remains a valid evaluator-only retrospective authority.\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    contract = tomllib.loads(args.contract.read_text(encoding="utf-8"))
    audit = run_audit(repo_root, contract)
    write_outputs(args.output, audit)
    print(json.dumps(audit["decision"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
