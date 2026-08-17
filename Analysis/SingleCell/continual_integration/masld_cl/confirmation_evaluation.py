"""Derive seed, held-study, and stress gates from locked evaluation artifacts."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .audits import verify_descriptive_projection, verify_gpu_tolerance
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .execution import require_execution_ownership, verify_execution_record
from .firewall import load_selection_lock


def _read_metrics(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"metric", "scope", "value"}.issubset(rows[0]):
        raise ContractError(f"evaluation metrics are malformed: {path}")
    result: dict[str, list[float]] = {}
    for row in rows:
        value = float(row["value"])
        if not (value == value and abs(value) != float("inf")):
            raise ContractError(f"evaluation metric is non-finite: {path}")
        result.setdefault(row["metric"], []).append(value)
    return result


def _verify_evaluation(
    path_value: str, schema: str, selection_sha: str, config: dict[str, Any],
) -> tuple[dict[str, list[float]], dict[str, Any], Path]:
    path = Path(path_value).resolve()
    details_path = path.with_suffix(path.suffix + ".details.json")
    if not details_path.is_file():
        raise ContractError(f"evaluation lacks provenance sidecar: {path}")
    with details_path.open() as handle:
        details = json.load(handle)
    if (
        details.get("schema_version") != schema
        or details.get("config_sha256") != config["_config_sha256"]
        or details.get("selection_lock_sha256") != selection_sha
        or details.get("metrics_realpath") != str(path)
        or details.get("metrics_sha256") != sha256_path(path)
    ):
        raise ContractError(f"evaluation provenance mismatch: {path}")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    records = details.get("execution_records", [])
    if not records:
        raise ContractError("evaluation provenance lacks GPU execution records")
    verified = {}
    for source in records:
        record_path = str(Path(source["path"]).resolve())
        if sha256_path(record_path) != source["sha256"]:
            raise ContractError("evaluation GPU record changed")
        verified[record_path] = verify_execution_record(
            record_path, pipeline_root, config["_config_sha256"]
        )
    if schema == "masld-cl-benchmark-v1":
        for source in details.get("embedding_sources", {}).values():
            load_embedding(source["embedding_manifest"])
            if sha256_path(source["embedding_manifest"]) != source["embedding_manifest_sha256"]:
                raise ContractError("benchmark embedding changed after evaluation")
            if "run_manifest" in source and (
                sha256_path(source["run_manifest"]) != source["run_manifest_sha256"]
            ):
                raise ContractError("benchmark run manifest changed after evaluation")
            source_records = source.get("execution_records", [])
            if source_records:
                if len(source_records) != 1:
                    raise ContractError("benchmark source requires exactly one GPU owner")
                record_path = str(Path(source_records[0]["path"]).resolve())
                if record_path not in verified or sha256_path(record_path) != source_records[0]["sha256"]:
                    raise ContractError("benchmark source GPU ownership record changed")
                owned = [source["embedding_manifest"]]
                if "run_manifest" in source:
                    owned.append(source["run_manifest"])
                require_execution_ownership(verified[record_path], owned, role="benchmark source")
    else:
        load_embedding(details["candidate_embedding"])
        for path_key, hash_key in (
            ("candidate_embedding", "candidate_embedding_sha256"),
            ("candidate_run_manifest", "candidate_run_manifest_sha256"),
            ("raw_pca_manifest", "raw_pca_manifest_sha256"),
        ):
            if sha256_path(details[path_key]) != details[hash_key]:
                raise ContractError("outcome evaluation source changed after evaluation")
        for source in details.get("raw_pca_scope_sources", []):
            if sha256_path(source["path"]) != source["sha256"]:
                raise ContractError("raw-PCA scope changed after outcome evaluation")
        if len(verified) != 1:
            raise ContractError("outcome evaluation requires exactly one GPU owner")
        require_execution_ownership(
            next(iter(verified.values())),
            [details["candidate_embedding"], details["candidate_run_manifest"]],
            role="outcome source",
        )
    return _read_metrics(path), details, path


def _derive_passes(
    config: dict[str, Any], benchmark: dict[str, list[float]],
    outcome: dict[str, list[float]],
) -> dict[str, bool]:
    gates = config["gates"]
    required_benchmark = {
        "reference_macro_f1_change_ci_low", "reference_neighborhood_jaccard_loss",
        "query_macro_f1_change", "major_lineage_f1_change",
        "alignment_improvement_vs_harmony", "alignment_improvement_vs_harmony_ci_low",
        "alignment_improvement_vs_architecture_surgery",
        "alignment_improvement_vs_architecture_surgery_ci_low", "improved_lineage_count",
        "all_lineage_alignment_pass", "maximum_protocol_stratum_worsening",
    }
    required_outcome = {
        "pooled_disease_retention", "per_study_disease_retention",
        "within_study_distance_spearman",
    }
    if not required_benchmark.issubset(benchmark) or not required_outcome.issubset(outcome):
        raise ContractError("confirmation evaluation lacks required directly computed metrics")
    reference = (
        min(benchmark["reference_macro_f1_change_ci_low"])
        > -gates["reference_macro_f1_margin"]
        and max(benchmark["reference_neighborhood_jaccard_loss"])
        <= gates["reference_neighborhood_jaccard_loss"]
        and min(benchmark["query_macro_f1_change"])
        >= -gates["query_macro_f1_margin"]
        and min(benchmark["major_lineage_f1_change"])
        >= -gates["major_lineage_f1_margin"]
    )
    alignment = (
        min(benchmark["alignment_improvement_vs_harmony"])
        >= gates["minimum_control_alignment_improvement"]
        and min(benchmark["alignment_improvement_vs_harmony_ci_low"]) > 0
        and min(benchmark["alignment_improvement_vs_architecture_surgery"])
        >= gates["minimum_control_alignment_improvement"]
        and min(benchmark["alignment_improvement_vs_architecture_surgery_ci_low"]) > 0
        and min(benchmark["improved_lineage_count"]) >= gates["minimum_improved_lineages"]
        and min(benchmark["all_lineage_alignment_pass"]) == 1
        and max(benchmark["maximum_protocol_stratum_worsening"])
        <= gates["maximum_stratum_worsening"]
    )
    disease = (
        min(outcome["pooled_disease_retention"]) >= gates["pooled_disease_retention"]
        and min(outcome["per_study_disease_retention"])
        >= gates["per_study_disease_retention"]
        and min(outcome["within_study_distance_spearman"])
        >= gates["within_study_distance_spearman"]
    )
    return {"reference_pass": reference, "alignment_pass": alignment, "disease_pass": disease}


def _required_run_manifests(benchmark_details: dict[str, Any], outcome_details: dict[str, Any]) -> set[str]:
    required = {
        str(Path(source["run_manifest"]).resolve())
        for source in benchmark_details["embedding_sources"].values()
        if "run_manifest" in source
    }
    required.add(str(Path(outcome_details["candidate_run_manifest"]).resolve()))
    return required


def evaluate_confirmation_bundle(
    config: dict[str, Any], selection_lock: str | Path,
    bundle_path: str | Path, gpu_tolerance_path: str | Path,
    descriptive_projection_audit: str | Path, output: str | Path,
) -> list[dict[str, Any]]:
    selection = load_selection_lock(selection_lock, config)
    with Path(bundle_path).open() as handle:
        bundle = json.load(handle)
    if (
        bundle.get("schema_version") != "masld-cl-confirmation-bundle-v2"
        or bundle.get("selection_lock_sha256") != selection["lock_sha256"]
    ):
        raise ContractError("confirmation bundle schema or selection lock is invalid")
    entries = bundle.get("evaluations", [])
    expected = {
        *(('selected_seed', str(seed)) for seed in config["screen"]["confirmation_seeds"]),
        *(('held_study', study) for study in config["evaluation"]["powered_query_studies"]),
        *(('secondary_stress', study) for study in config["evaluation"]["secondary_stress_studies"]),
    }
    observed = [(row.get("evaluation_type"), row.get("evaluation_id")) for row in entries]
    if set(observed) != expected or len(observed) != len(expected):
        raise ContractError("confirmation bundle lacks the exact seed, held-study, and stress roster")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    verified_records: dict[str, dict[str, Any]] = {}
    evaluated = []
    for entry in entries:
        benchmark, benchmark_details, benchmark_path = _verify_evaluation(
            entry["benchmark_metrics"], "masld-cl-benchmark-v1", selection["lock_sha256"], config
        )
        outcome, outcome_details, outcome_path = _verify_evaluation(
            entry["outcome_metrics"], "masld-cl-outcome-evaluation-v1", selection["lock_sha256"], config
        )
        context = benchmark_details.get("candidate_context") or {}
        evaluation_type, evaluation_id = entry["evaluation_type"], entry["evaluation_id"]
        selected = selection["selected"]
        if (
            float(context.get("ewc_lambda")) != float(selected["ewc_lambda"])
            or float(context.get("replay_fraction")) != float(selected["replay_fraction"])
            or context.get("method") != "continual_learning"
        ):
            raise ContractError("confirmation candidate is not the locked selected setting")
        if evaluation_type == "selected_seed":
            if (
                context.get("seed") != int(evaluation_id)
                or context.get("held_out_datasets")
                or set(context.get("query_datasets", []))
                != set(config["evaluation"]["powered_query_studies"])
            ):
                raise ContractError("selected-seed evaluation context is incorrect")
        elif evaluation_type == "held_study":
            if (
                context.get("held_out_datasets") != [evaluation_id]
                or evaluation_id in set(context.get("query_datasets", []))
            ):
                raise ContractError("held-study evaluation did not truly hold out its study")
        elif (
            set(context.get("query_datasets", []))
            != {*config["evaluation"]["powered_query_studies"], evaluation_id}
            or context.get("held_out_datasets")
        ):
            raise ContractError("secondary stress evaluation has the wrong adaptation roster")
        records = entry.get("execution_records", [])
        if not records:
            raise ContractError("confirmation evaluation lacks GPU execution records")
        result_manifests = set()
        for value in records:
            record_path = str(Path(value).resolve())
            if record_path not in verified_records:
                verified_records[record_path] = verify_execution_record(
                    record_path, pipeline_root, config["_config_sha256"]
                )
            identity = verified_records[record_path]["result_identity"]
            result_manifests.add(identity["result_manifest_realpath"])
            output_root = Path(identity["output_realpath"])
            result_manifests.update(
                str((output_root / item["path"]).resolve())
                for item in identity["output_files"]
                if item["path"].endswith("_manifest.json")
            )
        required_runs = _required_run_manifests(benchmark_details, outcome_details)
        if not required_runs.issubset(result_manifests):
            raise ContractError("confirmation metrics are not backed by the supplied GPU records")
        required_records = {
            str(Path(source["path"]).resolve())
            for details in (benchmark_details, outcome_details)
            for source in details.get("execution_records", [])
        }
        if {str(Path(value).resolve()) for value in records} != required_records:
            raise ContractError("confirmation supplied GPU records differ from metric provenance")
        passes = _derive_passes(config, benchmark, outcome)
        evaluated.append({
            "evaluation_type": evaluation_type, "evaluation_id": evaluation_id,
            "benchmark_metrics": str(benchmark_path), "outcome_metrics": str(outcome_path),
            **passes,
        })
    gpu = verify_gpu_tolerance(config, selection, gpu_tolerance_path)
    projection = verify_descriptive_projection(
        config, selection, descriptive_projection_audit
    )
    for record_path in [
        *(source["execution_record"] for source in gpu["sources"]),
        projection["execution_record"],
    ]:
        resolved = str(Path(record_path).resolve())
        if resolved not in verified_records:
            verified_records[resolved] = verify_execution_record(
                resolved, pipeline_root, config["_config_sha256"]
            )
    projection_pass = projection.get("passed") is True
    seeds = [row for row in evaluated if row["evaluation_type"] == "selected_seed"]
    held = [row for row in evaluated if row["evaluation_type"] == "held_study"]
    stress = [row for row in evaluated if row["evaluation_type"] == "secondary_stress"]
    metrics = [
        {"metric": "selected_seed_reference_and_disease_pass", "scope": "five_seeds", "value": int(all(
            row["reference_pass"] and row["disease_pass"] for row in seeds
        ))},
        {"metric": "selected_seed_alignment_pass_count", "scope": "five_seeds", "value": sum(
            row["alignment_pass"] for row in seeds
        )},
        {"metric": "held_study_reference_and_disease_pass", "scope": "two_powered_studies", "value": int(all(
            row["reference_pass"] and row["disease_pass"] for row in held
        ))},
        {"metric": "held_study_alignment_pass", "scope": "two_powered_studies", "value": int(all(
            row["alignment_pass"] for row in held
        ))},
        {"metric": "secondary_stress_pass", "scope": "two_secondary_studies", "value": int(all(
            row["reference_pass"] and row["disease_pass"] and row["alignment_pass"]
            for row in stress
        ))},
        {"metric": "gpu_tolerance_measured", "scope": "same_seed_gpu_pair", "value": 1},
        {"metric": "descriptive_projection_only", "scope": "33_cells", "value": int(projection_pass)},
    ]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-confirmation-evaluation-v2",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "bundle_realpath": str(Path(bundle_path).resolve()),
        "bundle_sha256": sha256_path(bundle_path),
        "gpu_tolerance": gpu,
        "gpu_tolerance_realpath": str(Path(gpu_tolerance_path).resolve()),
        "gpu_tolerance_sha256": sha256_path(gpu_tolerance_path),
        "descriptive_projection": projection,
        "descriptive_projection_realpath": str(Path(descriptive_projection_audit).resolve()),
        "descriptive_projection_sha256": sha256_path(descriptive_projection_audit),
        "evaluations": evaluated,
        "execution_records": [
            {"path": path, "sha256": sha256_path(path)}
            for path in sorted(verified_records)
        ],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
    })
    return metrics
