#!/usr/bin/env python3
"""Reuse frozen GSE296875 PeakVI states for blind GSE244832 inference only."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.peakvi_cross_cohort_fit_predict import (
    LINEAGES,
    ROTATIONS,
    PeakVICrossCohortError,
    array_record,
    digest,
    encode,
    feature_rows,
    make_peakvae,
    parameter_hash,
    read_json,
)


SCHEMA = "masld-bench-peakvi-gse244832-inference-only-v1"
SEEDS = (1103, 2207, 3301)
DONORS = tuple(f"D{index:02d}" for index in range(1, 19))


def verify_bound_tree(
    root: Path,
    record: Mapping[str, Any],
    artifact_class: str,
    *,
    verify_payload: bool = True,
) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as error:
        raise PeakVICrossCohortError(f"{artifact_class} escapes benchmark root") from error
    manifest = verify_frozen_tree(path) if verify_payload else read_json(path / "ARTIFACTS.json")
    if (
        path.is_symlink()
        or digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise PeakVICrossCohortError(f"{artifact_class} authority differs")
    return path


def read_target_cells(path: Path, unit_states: Any) -> list[dict[str, str]]:
    fields = (
        "cell_index", "cell_hash", "donor_id", "lineage_id", "outer_fold",
        "selection_rank", "available_nuclei",
    )
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise PeakVICrossCohortError("GSE244832 target cell-axis fields differ")
        rows = list(reader)
    donor_index = {donor: index for index, donor in enumerate(DONORS)}
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    groups: set[tuple[int, int]] = set()
    for offset, row in enumerate(rows):
        if (
            int(row["cell_index"]) != offset
            or row["donor_id"] not in donor_index
            or row["lineage_id"] not in lineage_index
            or int(row["available_nuclei"]) < 50
        ):
            raise PeakVICrossCohortError("GSE244832 target cell axis differs")
        group = (donor_index[row["donor_id"]], lineage_index[row["lineage_id"]])
        if int(unit_states[group]) != 0:
            raise PeakVICrossCohortError("ineligible GSE244832 cell entered model input")
        groups.add(group)
    expected = set(zip(*((unit_states == 0).nonzero()), strict=True))
    if (
        len(rows) != 3049
        or len({row["cell_hash"] for row in rows}) != len(rows)
        or groups != expected
        or len(groups) != 48
    ):
        raise PeakVICrossCohortError("GSE244832 eligible cell census differs")
    return rows


def aggregate_probabilities(
    probabilities: Any,
    target_cells: Sequence[Mapping[str, str]],
    unit_states: Any,
) -> Any:
    import numpy as np

    values = np.asarray(probabilities)
    if (
        values.shape != (len(target_cells), 16_000)
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
        or np.any(values > 1)
    ):
        raise PeakVICrossCohortError("per-cell accessibility probabilities differ")
    donor_index = {donor: index for index, donor in enumerate(DONORS)}
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    sums = np.zeros((18, 4, 16_000), dtype=np.float64)
    counts = np.zeros((18, 4), dtype=np.int64)
    for offset, row in enumerate(target_cells):
        group = (donor_index[row["donor_id"]], lineage_index[row["lineage_id"]])
        if int(unit_states[group]) != 0:
            raise PeakVICrossCohortError("ineligible target cell reached aggregation")
        sums[group] += values[offset]
        counts[group] += 1
    if np.any(counts[unit_states == 0] == 0) or np.any(counts[unit_states != 0] != 0):
        raise PeakVICrossCohortError("target donor-lineage aggregation mask differs")
    output = np.full((18, 4, 16_000), np.nan, dtype=np.float32)
    output[unit_states == 0] = (
        sums[unit_states == 0] / counts[unit_states == 0][:, None]
    ).astype(np.float32)
    if (
        not np.all(np.isfinite(output[unit_states == 0]))
        or not np.all(np.isnan(output[unit_states != 0]))
    ):
        raise PeakVICrossCohortError("masked donor-lineage prediction differs")
    return output


def decode_probabilities(
    torch: Any,
    decoder: Any,
    latent: Any,
    *,
    batch_size: int,
) -> Any:
    import numpy as np

    values = torch.from_numpy(latent).cuda()
    output: list[Any] = []
    decoder.eval()
    with torch.inference_mode():
        for start in range(0, len(latent), batch_size):
            probabilities = torch.sigmoid(decoder(values[start : start + batch_size]))
            if not torch.isfinite(probabilities).all():
                raise PeakVICrossCohortError("decoder probability is non-finite")
            output.append(probabilities.cpu().numpy())
    result = np.vstack(output).astype(np.float32, copy=False)
    if result.shape != (len(latent), 16_000):
        raise PeakVICrossCohortError("decoder probability geometry differs")
    return result


def state_paths(source_result: Path, rotation_id: str, seed: int) -> tuple[Path, Path]:
    if rotation_id not in ROTATIONS or seed not in SEEDS:
        raise PeakVICrossCohortError("source state identity differs")
    seed_root = source_result / "predictions" / rotation_id / f"seed-{seed}"
    peakvi = (seed_root / "peakvi.state_dict.pt").resolve(strict=True)
    decoder = (seed_root / "decoder.state_dict.pt").resolve(strict=True)
    for path in (peakvi, decoder):
        try:
            path.relative_to(source_result)
        except ValueError as error:
            raise PeakVICrossCohortError("source state escaped frozen result") from error
        if path.is_symlink() or not path.is_file():
            raise PeakVICrossCohortError("source state path differs")
    return peakvi, decoder


def validate_config(root: Path, config_path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    config = read_json(config_path.resolve(strict=True))
    if (
        config.get("schema_version") != SCHEMA
        or config.get("source_dataset_id") != "gse296875"
        or config.get("target_dataset_id") != "gse244832"
        or config.get("source_fit_job_id") != "21101210"
        or config.get("seeds") != list(SEEDS)
        or config.get("model", {}).get("model_id") != "peakvi"
        or config.get("model", {}).get("dimensions")
        != {"hidden": 128, "latent": 16, "encoder_layers": 2, "decoder_layers": 2}
        or config.get("execution", {}).get("fit_or_refit_performed") is not False
        or config.get("execution", {}).get("target_adaptation_performed") is not False
        or config.get("execution", {}).get("scoring_authorized") is not False
    ):
        raise PeakVICrossCohortError("GSE244832 inference-only campaign differs")
    for field in (
        "condition_values_read", "phenotype_values_read", "rna_assay_read",
        "cross_assay_join_consumed", "development_outcomes_read", "sealed_data_read",
        "metrics_calculated", "target_role_atac_consumed", "missing_as_zero",
    ):
        if config.get("firewall", {}).get(field) is not False:
            raise PeakVICrossCohortError(f"GSE244832 inference firewall opened: {field}")
    paths = {
        "source_result": verify_bound_tree(
            root,
            config["source_fit_result"],
            "peakvi_cross_cohort_fit_predict_and_blind_commit",
        ),
        "source_input": verify_bound_tree(
            root,
            config["source_cell_input"],
            "peakvi_gse296875_source_cell_input",
        ),
        "target_valid": verify_bound_tree(
            root,
            config["target_inputs"]["valid"],
            "peakvi_gse244832_target_cell_context",
        ),
        "target_test": verify_bound_tree(
            root,
            config["target_inputs"]["test"],
            "peakvi_gse244832_target_cell_context",
        ),
        "axis": verify_bound_tree(
            root,
            config["exchange_axis"],
            "gse244832_reference_guarded_atac_exchange_axis",
        ),
        "registration": verify_bound_tree(
            root,
            config["execution_registration"],
            "gse244832_observed_atac_execution_registration",
        ),
        "query_valid": verify_bound_tree(
            root,
            config["query_authorities"]["valid"],
            "gse244832_label_free_query_atac",
        ),
        "query_test": verify_bound_tree(
            root,
            config["query_authorities"]["test"],
            "gse244832_label_free_query_atac",
        ),
    }
    source_receipt = read_json(paths["source_result"] / "predictions" / "receipt.json")
    source_config_path = (root / config["source_fit_campaign"]["path"]).resolve(strict=True)
    if (
        digest(source_config_path) != config["source_fit_campaign"]["sha256"]
        or source_receipt.get("campaign_sha256") != config["source_fit_campaign"]["sha256"]
        or source_receipt.get("status") != "passed"
        or source_receipt.get("metrics_calculated") is not False
    ):
        raise PeakVICrossCohortError("frozen source-fit campaign differs")
    source_config = read_json(source_config_path)
    if source_config.get("source_input") != config["source_cell_input"]:
        raise PeakVICrossCohortError("source-fit revision2 input binding differs")
    registration = read_json(paths["registration"] / "execution_contract.json")
    if (
        registration.get("allowed_model_ids", []).count("peakvi") != 1
        or registration.get("prediction_commit_required_before_outcome_open") is not True
        or registration.get("sequence_execution_authorized") is not False
        or registration.get("rna_assay_opened") is not False
    ):
        raise PeakVICrossCohortError("GSE244832 execution registration differs")
    for role in ("valid", "test"):
        target_contract = read_json(paths[f"target_{role}"] / "contract.json")
        query_contract = read_json(paths[f"query_{role}"] / "query_contract.json")
        if (
            target_contract.get("context_role") != role
            or target_contract.get("eligible_donor_lineage_units") != 48
            or target_contract.get("registered_mask_applied_before_cell_selection") is not True
            or target_contract.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
            or target_contract.get("rna_assay_read") is not False
            or query_contract.get("context_role") != role
            or query_contract.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
            or query_contract.get("metrics_calculated") is not False
        ):
            raise PeakVICrossCohortError(f"GSE244832 target context differs: {role}")
    runtime_manifest = Path(str(config["model_manifest"]["path"])).resolve(strict=True)
    if (
        runtime_manifest.is_symlink()
        or digest(runtime_manifest) != config["model_manifest"]["sha256"]
    ):
        raise PeakVICrossCohortError("PeakVI model manifest differs")
    paths["model_manifest"] = runtime_manifest
    for name, record in config.get("implementation", {}).items():
        path = (root / str(record.get("path", ""))).resolve(strict=True)
        try:
            path.relative_to(root)
        except ValueError as error:
            raise PeakVICrossCohortError(f"implementation escapes root: {name}") from error
        if path.is_symlink() or not path.is_file() or digest(path) != record.get("sha256"):
            raise PeakVICrossCohortError(f"implementation differs: {name}")
    return config, paths


def make_prediction_manifest(
    *,
    config: Mapping[str, Any],
    config_sha256: str,
    rotation_id: str,
    model_manifest: Path,
    query_authority: Path,
    predictions: Path,
    missing_state: Path,
    output: Path,
) -> dict[str, Any]:
    context_role, target_role = ROTATIONS[rotation_id]
    return {
        "schema_version": "masld-bench-atac-transport-prediction-bundle-v1",
        "bundle_id": f"peakvi-source-state-inference-{rotation_id}-{config_sha256[:16]}",
        "model_id": "peakvi",
        "input_regime": "observed_atac",
        "output_family": "masked_accessibility_count",
        "prediction_layout": "donor_lineage",
        "rotation_id": rotation_id,
        "axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
        "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse244832",
        "query_atac_authority": {
            "path": str(query_authority),
            "artifacts_sha256": config["query_authorities"][context_role]["artifacts_sha256"],
            "context_role": context_role,
        },
        "model_manifest": {
            "path": model_manifest.relative_to(output).as_posix(),
            "sha256": digest(model_manifest),
        },
        "arrays": {
            "predictions": array_record(output, predictions, "float32"),
            "missing_state": array_record(output, missing_state, "uint8"),
        },
        "native_family_contract": {
            "native_output_family": "joint_representation",
            "native_target_latents_preserved": True,
            "projected_output_family": "masked_accessibility_count",
            "projection": "frozen_source_only_logistic_fixed_window_decoder",
            "source_fit_states_reused_without_refit": True,
            "direct_native_to_profile_ranking_allowed": False,
        },
        "prediction_scale": "mean_predicted_single_nucleus_binary_accessibility_probability",
        "seeds": list(SEEDS),
        "evidence_contract": {
            "condition_labels_consumed": False,
            "phenotype_values_consumed": False,
            "rna_assay_consumed": False,
            "cross_assay_join_consumed": False,
            "sequence_features_consumed": False,
            "source_native_reference_bundle_assumed_resolved": False,
            "development_outcomes_consumed": False,
            "sealed_data_consumed": False,
            "target_role_atac_consumed": False,
            "target_role_contigs_entered_query_transform": False,
            "missing_as_zero": False,
            "scored_target_role": target_role,
            "observed_atac_context_role": context_role,
            "query_atac_consumed": True,
            "receptive_field_bp": "not_applicable_global_observed_atac_encoder",
            "query_transform_scope": "source_training_only",
            "per_donor_context_operation": "not_applicable",
        },
        "metrics_calculated": False,
        "champion_claim_allowed": False,
    }


def infer_rotation(
    *,
    torch: Any,
    config: Mapping[str, Any],
    config_sha256: str,
    paths: Mapping[str, Path],
    unit_states: Any,
    rotation_id: str,
    output: Path,
) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse

    context_role, target_role = ROTATIONS[rotation_id]
    target_root = paths[f"target_{context_role}"]
    target_cells = read_target_cells(target_root / "cell_axis.tsv", unit_states)
    source_features = feature_rows(paths["source_input"] / f"{context_role}.windows.tsv")
    target_features = feature_rows(target_root / "context.windows.tsv")
    source_target_features = feature_rows(paths["source_input"] / f"{target_role}.windows.tsv")
    if (
        source_features != target_features
        or {row["contig"] for row in source_features}
        & {row["contig"] for row in source_target_features}
    ):
        raise PeakVICrossCohortError("source/target rotation feature geometry differs")
    target_context = sparse.load_npz(target_root / "context.counts.npz").tocsr()
    if (
        target_context.shape != (3049, 16_000)
        or (target_context.data < 0).any()
        or (target_context.getnnz(axis=1) == 0).any()
    ):
        raise PeakVICrossCohortError("GSE244832 target sparse context differs")
    source_rotation = paths["source_result"] / "predictions" / rotation_id
    source_manifest = verify_frozen_tree(source_rotation)
    source_fit = read_json(source_rotation / "fit_receipt.json")
    if (
        source_manifest.get("metadata", {}).get("artifact_class")
        != "peakvi_gse281367_blind_prediction"
        or source_fit.get("rotation_id") != rotation_id
        or source_fit.get("metrics_calculated") is not False
    ):
        raise PeakVICrossCohortError("frozen source rotation differs")
    source_seeds = {int(row["seed"]): row for row in source_fit.get("seeds", [])}
    if set(source_seeds) != set(SEEDS):
        raise PeakVICrossCohortError("frozen source seed roster differs")
    output.mkdir(parents=True, mode=0o750)
    seed_predictions: list[Any] = []
    seed_receipts: list[dict[str, Any]] = []
    dimensions = config["model"]["dimensions"]
    batch_size = int(config["execution"]["batch_size"])
    for seed in SEEDS:
        peakvi_path, decoder_path = state_paths(paths["source_result"], rotation_id, seed)
        source_seed = source_seeds[seed]
        if (
            digest(peakvi_path) != source_seed["peakvi_state_file_sha256"]
            or digest(decoder_path) != source_seed["decoder_state_file_sha256"]
        ):
            raise PeakVICrossCohortError("frozen source state file hash differs")
        model = make_peakvae(16_000, dimensions)
        peakvi_state = torch.load(peakvi_path, map_location="cpu", weights_only=True)
        model.load_state_dict(peakvi_state, strict=True)
        if parameter_hash(model.state_dict()) != source_seed["peakvi_parameter_sha256"]:
            raise PeakVICrossCohortError("strict frozen PeakVI restore differs")
        target_latent = encode(
            torch,
            model,
            target_context,
            batch_size=batch_size,
            latent_dimension=int(dimensions["latent"]),
        )
        repeated = encode(
            torch,
            model,
            target_context,
            batch_size=batch_size,
            latent_dimension=int(dimensions["latent"]),
        )
        if not np.array_equal(target_latent, repeated):
            raise PeakVICrossCohortError("GSE244832 target encoding is not repeatable")
        decoder = torch.nn.Linear(int(dimensions["latent"]), 16_000).cuda()
        decoder_state = torch.load(decoder_path, map_location="cpu", weights_only=True)
        decoder.load_state_dict(decoder_state, strict=True)
        if parameter_hash(decoder.state_dict()) != source_seed["decoder_parameter_sha256"]:
            raise PeakVICrossCohortError("strict frozen decoder restore differs")
        probabilities = decode_probabilities(
            torch, decoder, target_latent, batch_size=batch_size
        )
        prediction = aggregate_probabilities(probabilities, target_cells, unit_states)
        seed_root = output / f"seed-{seed}"
        seed_root.mkdir(mode=0o750)
        np.save(seed_root / "target_latent.float32.npy", target_latent, allow_pickle=False)
        np.save(seed_root / "prediction.float32.npy", prediction, allow_pickle=False)
        seed_predictions.append(prediction)
        seed_receipts.append(
            {
                "seed": seed,
                "source_peakvi_state_path": str(peakvi_path),
                "source_peakvi_state_sha256": digest(peakvi_path),
                "source_decoder_state_path": str(decoder_path),
                "source_decoder_state_sha256": digest(decoder_path),
                "source_peakvi_parameter_sha256": source_seed["peakvi_parameter_sha256"],
                "source_decoder_parameter_sha256": source_seed["decoder_parameter_sha256"],
                "strict_state_restore_passed": True,
                "target_encoding_repeat_bit_identical": True,
                "fit_or_refit_performed": False,
            }
        )
        del model, decoder, peakvi_state, decoder_state
        del target_latent, repeated, probabilities, prediction
        torch.cuda.empty_cache()
    observed = unit_states == 0
    stacked = np.stack(seed_predictions, axis=0)
    ensemble = np.full((18, 4, 16_000), np.nan, dtype=np.float32)
    ensemble[observed] = np.mean(
        stacked[:, observed, :], axis=0, dtype=np.float64
    ).astype(np.float32)
    missing = np.broadcast_to(unit_states[:, :, None], ensemble.shape).copy()
    predictions_path = output / "predictions.float32.npy"
    missing_path = output / "missing_state.uint8.npy"
    np.save(predictions_path, ensemble, allow_pickle=False)
    np.save(missing_path, missing, allow_pickle=False)
    model_execution = {
        "schema_version": "masld-bench-peakvi-gse244832-model-execution-v1",
        "model_id": "peakvi",
        "rotation_id": rotation_id,
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse244832",
        "fit_scope": "gse296875_source_training_only",
        "source_fit_job_id": "21101210",
        "source_fit_result_path": str(paths["source_result"]),
        "source_fit_result_artifacts_sha256": config["source_fit_result"]["artifacts_sha256"],
        "source_input_path": str(paths["source_input"]),
        "source_input_artifacts_sha256": config["source_cell_input"]["artifacts_sha256"],
        "upstream_peakvi_model_manifest_path": str(paths["model_manifest"]),
        "upstream_peakvi_model_manifest_sha256": config["model_manifest"]["sha256"],
        "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
        "target_dataset_fit_or_adaptation": False,
        "fit_or_refit_performed": False,
        "sequence_features_consumed": False,
        "rna_assay_consumed": False,
        "cross_assay_join_consumed": False,
        "query_transform_scope": "source_training_only",
        "per_donor_context_operation": "not_applicable",
        "mask_applied_before_query_transform": True,
        "ineligible_query_units_read_by_transform": False,
        "fragment_cut_sites": "start_plus_4_and_end_minus_5",
        "registered_output_family": "masked_accessibility_count",
        "prediction_scale": "mean_predicted_single_nucleus_binary_accessibility_probability",
        "prediction_is_observed_fragment_count": False,
        "seeds": seed_receipts,
        "metrics_calculated": False,
    }
    model_execution_path = output / "model_execution_manifest.json"
    write_json_exclusive(model_execution_path, model_execution)
    write_json_exclusive(
        output / "native_representation_manifest.json",
        {
            "schema_version": "masld-bench-peakvi-native-representation-v1",
            "model_id": "peakvi",
            "output_family": "joint_representation",
            "rotation_id": rotation_id,
            "context_role": context_role,
            "target_cell_axis_sha256": digest(target_root / "cell_axis.tsv"),
            "target_context_features_sha256": digest(target_root / "context.windows.tsv"),
            "seed_latent_paths": [
                {"seed": seed, "target": f"seed-{seed}/target_latent.float32.npy"}
                for seed in SEEDS
            ],
            "source_latents_reused_in_frozen_source_fit_result": True,
            "source_fit_states_reused_without_refit": True,
            "direct_profile_ranking_allowed": False,
            "target_query_adaptation_performed": False,
            "target_role_atac_consumed": False,
        },
    )
    write_json_exclusive(
        output / "prediction_bundle.json",
        make_prediction_manifest(
            config=config,
            config_sha256=config_sha256,
            rotation_id=rotation_id,
            model_manifest=model_execution_path,
            query_authority=paths[f"query_{context_role}"],
            predictions=predictions_path,
            missing_state=missing_path,
            output=output,
        ),
    )
    write_json_exclusive(
        output / "inference_receipt.json",
        {
            "schema_version": "masld-bench-peakvi-gse244832-rotation-inference-receipt-v1",
            "status": "passed",
            "rotation_id": rotation_id,
            "source_fit_job_id": "21101210",
            "target_cells": len(target_cells),
            "eligible_target_units": 48,
            "structurally_missing_units": 4,
            "below_qc_units": 20,
            "seeds": seed_receipts,
            "fit_or_refit_performed": False,
            "target_adaptation_performed": False,
            "native_output_family": "joint_representation",
            "projected_output_family": "masked_accessibility_count",
            "condition_labels_read": False,
            "phenotype_values_read": False,
            "rna_assay_read": False,
            "cross_assay_join_consumed": False,
            "development_outcomes_read": False,
            "sealed_data_read": False,
            "metrics_calculated": False,
            "champion_claim_allowed": False,
        },
    )
    artifact_sha = freeze_tree(
        output,
        {
            "artifact_class": "peakvi_gse244832_blind_prediction",
            "model_id": "peakvi",
            "rotation_id": rotation_id,
            "native_output_family": "joint_representation",
            "projected_output_family": "masked_accessibility_count",
            "source_fit_job_id": "21101210",
            "fit_or_refit_performed": False,
            "rna_assay_read": False,
            "cross_assay_join_consumed": False,
            "development_outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return {
        "rotation_id": rotation_id,
        "prediction_artifacts_sha256": artifact_sha,
        "seed_count": len(SEEDS),
        "fit_or_refit_performed": False,
        "metrics_calculated": False,
    }


def run(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    import numpy as np
    import scvi
    import torch

    root = root.resolve(strict=True)
    if output.exists():
        raise PeakVICrossCohortError("refusing to overwrite GSE244832 inference output")
    config, paths = validate_config(root, config_path)
    if scvi.__version__ != "1.5.0.post1" or torch.__version__ != "2.6.0+cu124":
        raise PeakVICrossCohortError("PeakVI runtime version differs")
    if torch.cuda.device_count() != 1 or "L40S" not in torch.cuda.get_device_name(0):
        raise PeakVICrossCohortError("exactly one L40S is required")
    torch.use_deterministic_algorithms(True)
    registration = read_json(paths["registration"] / "execution_contract.json")
    unit_states = np.load(
        paths["registration"] / registration["query_lineage_missing_state_path"],
        allow_pickle=False,
    )
    if (
        unit_states.shape != (18, 4)
        or unit_states.dtype != np.uint8
        or {value: int((unit_states == value).sum()) for value in (0, 1, 3)}
        != {0: 48, 1: 4, 3: 20}
    ):
        raise PeakVICrossCohortError("GSE244832 registered missingness differs")
    output.mkdir(parents=True, mode=0o750)
    config_sha256 = digest(config_path.resolve(strict=True))
    rotations = [
        infer_rotation(
            torch=torch,
            config=config,
            config_sha256=config_sha256,
            paths=paths,
            unit_states=unit_states,
            rotation_id=rotation_id,
            output=output / rotation_id,
        )
        for rotation_id in ROTATIONS
    ]
    receipt = {
        "schema_version": "masld-bench-peakvi-gse244832-inference-only-receipt-v1",
        "status": "passed",
        "campaign_sha256": config_sha256,
        "source_fit_job_id": "21101210",
        "rotations": rotations,
        "runtime": {
            "scvi_tools": scvi.__version__,
            "torch": torch.__version__,
            "cuda_build": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(0),
        },
        "fit_or_refit_performed": False,
        "target_adaptation_performed": False,
        "native_output_family": "joint_representation",
        "projected_output_family": "masked_accessibility_count",
        "condition_labels_read": False,
        "phenotype_values_read": False,
        "rna_assay_read": False,
        "cross_assay_join_consumed": False,
        "development_outcomes_read": False,
        "sealed_data_read": False,
        "metrics_calculated": False,
        "champion_claim_allowed": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    run(arguments.root, arguments.config, arguments.output)


if __name__ == "__main__":
    main()
