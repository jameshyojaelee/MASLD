"""Lock all promotion metrics and their generating artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from .config import canonical_json_bytes, write_json_exclusive
from .confirmation_matrix import load_confirmation_matrix_lock
from .contracts import ContractError, sha256_path, verify_contract_lock
from .data import verify_prepared_file
from .firewall import load_selection_lock, validate_program_firewall
from .execution import verify_execution_record
from .refit import load_refit_lock


REQUIRED_METRICS = {
    "reference_macro_f1_change_ci_low", "reference_neighborhood_jaccard_loss",
    "query_macro_f1_change", "major_lineage_f1_change",
    "alignment_improvement_vs_harmony", "alignment_improvement_vs_harmony_ci_low",
    "alignment_improvement_vs_architecture_surgery",
    "alignment_improvement_vs_architecture_surgery_ci_low", "improved_lineage_count",
    "all_lineage_alignment_pass", "maximum_protocol_stratum_worsening",
    "pooled_disease_retention", "per_study_disease_retention",
    "within_study_distance_spearman", "all_orders_reference_and_disease_pass",
    "passing_order_count", "worst_order_degradation", "order_distance_spearman",
    "order_to_seed_variability_ratio",
    "selected_seed_reference_and_disease_pass", "selected_seed_alignment_pass_count",
    "held_study_reference_and_disease_pass", "held_study_alignment_pass",
    "secondary_stress_pass", "gpu_tolerance_measured", "descriptive_projection_only",
}
ALLOWED_METRIC_SOURCE_SCHEMAS = {
    "masld-cl-benchmark-v1", "masld-cl-outcome-evaluation-v1",
    "masld-cl-confirmation-evaluation-v2", "masld-cl-order-evaluation-v2",
}


def _read_metric_file(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"metric", "scope", "value"}.issubset(rows[0]):
        raise ContractError(f"metric file is malformed: {path}")
    for row in rows:
        value = float(row["value"])
        if not (value == value and abs(value) != float("inf")):
            raise ContractError(f"metric value is non-finite: {path}")
    return [{"metric": row["metric"], "scope": row["scope"], "value": row["value"]} for row in rows]


def _verify_metric_source(
    config: dict[str, Any], selection: dict[str, Any], path: Path,
) -> tuple[dict[str, Any], Path]:
    details_path = path.with_suffix(path.suffix + ".details.json")
    if not details_path.is_file():
        raise ContractError(f"promotion metric source lacks provenance sidecar: {path}")
    with details_path.open() as handle:
        details = json.load(handle)
    schema = details.get("schema_version")
    if (
        schema not in ALLOWED_METRIC_SOURCE_SCHEMAS
        or details.get("config_sha256") != config["_config_sha256"]
        or details.get("selection_lock_sha256") != selection["lock_sha256"]
        or details.get("metrics_realpath") != str(path)
        or details.get("metrics_sha256") != sha256_path(path)
    ):
        raise ContractError(f"promotion metric provenance mismatch: {path}")
    from .confirmation_evaluation import _verify_evaluation
    if schema in {"masld-cl-benchmark-v1", "masld-cl-outcome-evaluation-v1"}:
        _verify_evaluation(path, schema, selection["lock_sha256"], config)
    elif schema == "masld-cl-confirmation-evaluation-v2":
        from .audits import verify_descriptive_projection, verify_gpu_tolerance
        for path_key, hash_key in (
            ("bundle_realpath", "bundle_sha256"),
            ("gpu_tolerance_realpath", "gpu_tolerance_sha256"),
            ("descriptive_projection_realpath", "descriptive_projection_sha256"),
        ):
            if sha256_path(details[path_key]) != details[hash_key]:
                raise ContractError("confirmation metric source changed")
        verify_gpu_tolerance(config, selection, details["gpu_tolerance_realpath"])
        verify_descriptive_projection(
            config, selection, details["descriptive_projection_realpath"]
        )
        with Path(details["bundle_realpath"]).open() as handle:
            bundle = json.load(handle)
        for entry in bundle.get("evaluations", []):
            _verify_evaluation(
                entry["benchmark_metrics"], "masld-cl-benchmark-v1",
                selection["lock_sha256"], config,
            )
            _verify_evaluation(
                entry["outcome_metrics"], "masld-cl-outcome-evaluation-v1",
                selection["lock_sha256"], config,
            )
    else:
        if sha256_path(details["bundle_realpath"]) != details["bundle_sha256"]:
            raise ContractError("order metric bundle changed")
        with Path(details["bundle_realpath"]).open() as handle:
            bundle = json.load(handle)
        evaluations = [bundle.get("joint", {}), *bundle.get("orders", [])]
        for entry in evaluations:
            _verify_evaluation(
                entry["benchmark_metrics"], "masld-cl-benchmark-v1",
                selection["lock_sha256"], config,
            )
            _verify_evaluation(
                entry["outcome_metrics"], "masld-cl-outcome-evaluation-v1",
                selection["lock_sha256"], config,
            )
    records = details.get("execution_records", [])
    if not records or any(
        not isinstance(source, dict) or set(source) != {"path", "sha256"}
        for source in records
    ):
        raise ContractError("promotion metric source lacks exact execution-record provenance")
    return details, details_path


def lock_metric_bundle(
    config: dict[str, Any], contract_lock: str | Path,
    prepared: str | Path, prepared_lock: str | Path,
    selection_lock: str | Path, refit_lock: str | Path,
    confirmation_matrix_lock: str | Path,
    metric_files: list[str | Path], execution_records: list[str | Path],
    output: str | Path,
) -> dict[str, Any]:
    contract_lock_path = Path(contract_lock).resolve()
    prepared_path = Path(prepared).resolve()
    prepared_lock_path = Path(prepared_lock).resolve()
    selection_lock_path = Path(selection_lock).resolve()
    refit_lock_path = Path(refit_lock).resolve()
    confirmation_matrix_path = Path(confirmation_matrix_lock).resolve()
    contract = verify_contract_lock(config, contract_lock_path, full_hash=True)
    prepared_identity = verify_prepared_file(prepared, prepared_lock, config, contract)
    selection = load_selection_lock(selection_lock_path, config)
    refit = load_refit_lock(refit_lock_path, config, selection)
    confirmation_matrix = load_confirmation_matrix_lock(
        confirmation_matrix_path, config, selection
    )
    firewall = validate_program_firewall(config)
    rows: list[dict[str, str]] = []
    metric_sources = []
    for value in metric_files:
        path = Path(value).resolve()
        details, details_path = _verify_metric_source(config, selection, path)
        rows.extend(_read_metric_file(path))
        metric_sources.append({
            "path": str(path), "sha256": sha256_path(path),
            "details_path": str(details_path.resolve()),
            "details_sha256": sha256_path(details_path),
            "schema_version": details["schema_version"],
        })
    observed_metrics = {row["metric"] for row in rows}
    missing = sorted(REQUIRED_METRICS - observed_metrics)
    if missing:
        raise ContractError(f"promotion metric bundle is incomplete: {missing}")
    pairs = [(row["metric"], row["scope"]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise ContractError("promotion metric bundle contains duplicate metric/scope rows")
    required_record_sources: dict[str, str] = {}
    for source in metric_sources:
        with Path(source["details_path"]).open() as handle:
            details = json.load(handle)
        for record_source in details["execution_records"]:
            path = str(Path(record_source["path"]).resolve())
            prior = required_record_sources.setdefault(path, record_source["sha256"])
            if prior != record_source["sha256"]:
                raise ContractError("metric sources disagree about an execution record")
    supplied = [str(Path(value).resolve()) for value in execution_records]
    if len(supplied) != len(set(supplied)) or set(supplied) != set(required_record_sources):
        raise ContractError("supplied execution records differ from exact metric provenance")
    manifests = []
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for value in sorted(supplied):
        path = Path(value)
        if sha256_path(path) != required_record_sources[str(path)]:
            raise ContractError("execution record changed after metric evaluation")
        record = verify_execution_record(
            path, pipeline_root, config["_config_sha256"]
        )
        with Path(record["result_identity"]["result_manifest_realpath"]).open() as handle:
            result_manifest = json.load(handle)
        if result_manifest.get("sensitivity_only") is True:
            raise ContractError("sensitivity-only BI replay cannot enter promotion")
        manifests.append({"path": str(path), "sha256": sha256_path(path)})
    if not manifests:
        raise ContractError("promotion requires explicit run manifests")
    automatic = [
        {"metric": "contract_valid", "scope": "immutable_atlas", "value": "1"},
        {"metric": "runs_complete", "scope": "all_required_runs", "value": "1"},
        {"metric": "input_hashes_valid", "scope": "all_locked_inputs", "value": "1"},
        {"metric": "frozen_program_scores_unchanged", "scope": "117_program_firewall", "value": "1"},
        {"metric": "weighted_bh_family_size", "scope": "declared_contrast", "value": str(config["program_inference_contract"]["weighted_bh_family_size"])},
        {"metric": "unweighted_bh_family_size", "scope": "declared_contrast", "value": str(config["program_inference_contract"]["unweighted_bh_family_size"])},
    ]
    rows.extend(automatic)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lock = {
        "schema_version": "masld-cl-metric-bundle-v1",
        "config_sha256": config["_config_sha256"],
        "contract_lock_sha256": contract["lock_sha256"],
        "prepared_lock_sha256": prepared_identity["lock_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "refit_lock_sha256": refit["lock_sha256"],
        "confirmation_matrix_lock_sha256": confirmation_matrix["lock_sha256"],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
        "metric_sources": metric_sources,
        "execution_records": manifests,
        "immutable_inputs": {
            "contract_lock": {"path": str(contract_lock_path), "sha256": sha256_path(contract_lock_path)},
            "prepared": {"path": str(prepared_path), "sha256": prepared_identity["prepared_file_sha256"]},
            "prepared_lock": {"path": str(prepared_lock_path), "sha256": sha256_path(prepared_lock_path)},
            "selection_lock": {"path": str(selection_lock_path), "sha256": sha256_path(selection_lock_path)},
            "refit_lock": {"path": str(refit_lock_path), "sha256": sha256_path(refit_lock_path)},
            "confirmation_matrix_lock": {
                "path": str(confirmation_matrix_path),
                "sha256": sha256_path(confirmation_matrix_path),
            },
        },
        "program_firewall": firewall,
    }
    lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
    write_json_exclusive(output.with_suffix(output.suffix + ".lock.json"), lock)
    return lock


def verify_metric_bundle(
    config: dict[str, Any], metrics: str | Path, lock_path: str | Path,
    selection: dict[str, Any],
) -> dict[str, Any]:
    with Path(lock_path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-metric-bundle-v1":
        raise ContractError("unsupported metric-bundle lock")
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if hashlib.sha256(canonical_json_bytes(payload)).hexdigest() != lock.get("lock_sha256"):
        raise ContractError("metric-bundle lock content hash mismatch")
    if lock.get("config_sha256") != config["_config_sha256"]:
        raise ContractError("metric-bundle config mismatch")
    if lock.get("selection_lock_sha256") != selection["lock_sha256"]:
        raise ContractError("metric-bundle selection mismatch")
    path = Path(metrics).resolve()
    if str(path) != lock.get("metrics_realpath") or sha256_path(path) != lock.get("metrics_sha256"):
        raise ContractError("promotion metrics changed after locking")
    immutable = lock.get("immutable_inputs", {})
    required_inputs = {
        "contract_lock", "prepared", "prepared_lock", "selection_lock",
        "refit_lock", "confirmation_matrix_lock",
    }
    if set(immutable) != required_inputs:
        raise ContractError("metric bundle lacks immutable input paths")
    for source in immutable.values():
        if sha256_path(source["path"]) != source["sha256"]:
            raise ContractError(f"locked promotion input changed: {source['path']}")
    contract = verify_contract_lock(config, immutable["contract_lock"]["path"], full_hash=True)
    verify_prepared_file(
        immutable["prepared"]["path"], immutable["prepared_lock"]["path"],
        config, contract,
    )
    stored_selection = load_selection_lock(immutable["selection_lock"]["path"], config)
    if stored_selection["lock_sha256"] != selection["lock_sha256"]:
        raise ContractError("stored metric-bundle selection differs from requested selection")
    load_refit_lock(immutable["refit_lock"]["path"], config, selection)
    load_confirmation_matrix_lock(
        immutable["confirmation_matrix_lock"]["path"], config, selection
    )
    for source in lock["metric_sources"]:
        if (
            sha256_path(source["path"]) != source["sha256"]
            or sha256_path(source["details_path"]) != source["details_sha256"]
        ):
            raise ContractError(f"locked promotion source changed: {source['path']}")
        details, details_path = _verify_metric_source(
            config, selection, Path(source["path"])
        )
        if str(details_path.resolve()) != source["details_path"] or details["schema_version"] != source["schema_version"]:
            raise ContractError("metric source identity changed")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for source in lock["execution_records"]:
        if sha256_path(source["path"]) != source["sha256"]:
            raise ContractError(f"locked execution record changed: {source['path']}")
        verify_execution_record(source["path"], pipeline_root, config["_config_sha256"])
    return lock
