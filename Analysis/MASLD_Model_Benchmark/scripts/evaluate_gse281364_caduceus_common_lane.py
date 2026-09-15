#!/usr/bin/env python3
"""Fit identical fold-safe common heads to Caduceus 131-kb embeddings."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from scripts.evaluate_gse281364_dna_lm_common_lane import (
    BOOTSTRAP_REPLICATES,
    CONTEXTS,
    SEED,
    CommonLaneError,
    _projected_features,
    _value,
    _write_gzip_tsv,
    _write_tsv,
    activity_delta,
    array_digest,
    bootstrap_metrics,
    fit_mlp_outer,
    fit_ridge_outer,
    load_outcomes,
    metric_values,
    sha256_file,
    verify_frozen_tree,
)


MODEL_ID = "caduceus"
LICENSE = "Apache-2.0"


def load_fixture(root: Path) -> tuple[list[dict[str, str]], np.ndarray]:
    path = root / "fixture/sequence_manifest.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if (
        len(rows) != 1_033
        or len({row["fixture_id"] for row in rows}) != 1_033
        or len({row["element_id"] for row in rows}) != 1_033
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != 1_033
        or {int(row["outer_fold"]) for row in rows} != set(range(5))
        or any(int(row["input_length_bp"]) != 131_072 for row in rows)
    ):
        raise CommonLaneError("Caduceus fixture census differs")
    bases = {base: index for index, base in enumerate("ACGT")}
    allele_features = np.zeros((len(rows), 16), dtype=np.float64)
    for index, row in enumerate(rows):
        ref, alt = row["ref"], row["alt"]
        if ref not in bases or alt not in bases or ref == alt:
            raise CommonLaneError("fixture allele identity differs")
        allele_features[index, 4 * bases[ref] + bases[alt]] = 1.0
    return rows, allele_features


def audit_and_load_features(
    root: Path, fixture_rows: Sequence[Mapping[str, str]]
) -> tuple[dict[int, np.ndarray], dict[str, Any]]:
    receipt = json.loads((root / "projected/receipt.json").read_text(encoding="utf-8"))
    raw_receipt = json.loads(
        (root / "raw/caduceus/receipt.json").read_text(encoding="utf-8")
    )
    if (
        receipt.get("status") != "pass_outcome_blind_fold_safe_projection"
        or receipt.get("elements") != 1_033
        or receipt.get("outer_locus_sequence_groups") != 1_033
        or receipt.get("outer_folds") != 5
        or receipt.get("projection_width") != 256
        or receipt.get("head_input_width") != 1_024
        or receipt.get("projection_fit_on_held_out_fold")
        or receipt.get("reporter_outcomes_read")
        or receipt.get("sealed_labels_read")
        or receipt.get("downstream_head_fit")
        or raw_receipt.get("status") != "pass_outcome_blind_embedding_extraction"
        or raw_receipt.get("model_id") != MODEL_ID
        or raw_receipt.get("license") != LICENSE
        or raw_receipt.get("exposure_state") != "target_label_unexposed"
        or raw_receipt.get("restricted_comparator")
        or not raw_receipt.get("open_champion_eligible_after_task_and_external_evaluation_gates")
        or raw_receipt.get("reporter_outcomes_read")
        or raw_receipt.get("sealed_labels_read")
        or raw_receipt.get("downstream_head_fit")
    ):
        raise CommonLaneError("Caduceus projection or raw receipt differs")
    expected_ids = np.asarray([row["fixture_id"] for row in fixture_rows])
    expected_groups = np.asarray(
        [row["outer_locus_sequence_group_id"] for row in fixture_rows]
    )
    expected_folds = np.asarray(
        [int(row["outer_fold"]) for row in fixture_rows], dtype=np.int64
    )
    with np.load(root / "raw/caduceus/allele_embeddings.npz", allow_pickle=False) as data:
        ids = data["fixture_ids"]
        groups = data["outer_locus_sequence_group_ids"]
        folds = data["outer_folds"].astype(np.int64)
        raw = data["embeddings"]
        alleles = tuple(data["allele_order"].tolist())
    if (
        not np.array_equal(ids, expected_ids)
        or not np.array_equal(groups, expected_groups)
        or not np.array_equal(folds, expected_folds)
        or alleles != ("REF", "ALT", "REF_RC", "ALT_RC")
        or raw.shape != (1_033, 4, 256)
        or not np.isfinite(raw).all()
    ):
        raise CommonLaneError("Caduceus raw alignment differs")
    model_receipt = receipt.get("models", {}).get(MODEL_ID, {})
    if model_receipt.get("license") != LICENSE or model_receipt.get("restricted_comparator"):
        raise CommonLaneError("Caduceus projected eligibility differs")
    features: dict[int, np.ndarray] = {}
    checks: list[dict[str, Any]] = []
    for held_fold in range(5):
        fold_root = root / f"projected/caduceus/heldout_fold{held_fold}"
        with np.load(fold_root / "projection_parameters.npz", allow_pickle=False) as data:
            parameters = {key: data[key] for key in data.files}
        with np.load(fold_root / "head_features.npz", allow_pickle=False) as data:
            fold_ids = data["fixture_ids"]
            fold_groups = data["outer_locus_sequence_group_ids"]
            fold_values = data["outer_folds"].astype(np.int64)
            block_order = tuple(data["feature_block_order"].tolist())
            saved_features = data["features"]
        train_raw = raw[folds != held_fold].reshape(-1, 256).astype(np.float64)
        expected_mean = train_raw.mean(axis=0)
        expected_std = np.maximum(train_raw.std(axis=0), 1.0e-6)
        components = np.asarray(parameters["components"], dtype=np.float64)
        if (
            int(parameters["held_out_fold"]) != held_fold
            or int(parameters["training_elements"]) != int(np.sum(folds != held_fold))
            or not np.array_equal(fold_ids, expected_ids)
            or not np.array_equal(fold_groups, expected_groups)
            or not np.array_equal(fold_values, expected_folds)
            or block_order != ("REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF")
            or saved_features.shape != (1_033, 1_024)
            or not np.allclose(parameters["mean"], expected_mean, rtol=0.0, atol=1e-12)
            or not np.allclose(parameters["standard_deviation"], expected_std, rtol=0.0, atol=1e-12)
            or components.shape != (256, 256)
            or not np.allclose(components @ components.T, np.eye(256), atol=2e-5)
            or any(row[int(np.argmax(np.abs(row)))] < 0 for row in components)
        ):
            raise CommonLaneError(f"Caduceus fold {held_fold} projection differs")
        reconstructed = _projected_features(raw, parameters)
        fold_receipt = model_receipt["folds"][held_fold]
        if (
            not np.allclose(reconstructed, saved_features, rtol=2e-6, atol=2e-6)
            or array_digest(saved_features) != fold_receipt["features_sha256"]
            or array_digest(components)
            != fold_receipt["projection_components_sha256"]
        ):
            raise CommonLaneError(f"Caduceus fold {held_fold} reconstruction differs")
        features[held_fold] = saved_features.astype(np.float64)
        checks.append(
            {
                "model_id": MODEL_ID,
                "held_out_fold": held_fold,
                "training_elements": int(np.sum(folds != held_fold)),
                "held_out_elements": int(np.sum(folds == held_fold)),
                "training_only_mean_std_rederived": True,
                "head_features_reconstructed": True,
                "orthonormal_projection": True,
            }
        )
    return features, {
        "checks": checks,
        "status": "pass_independent_projection_consistency_audit",
    }


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise CommonLaneError("evaluation output exists")
    embedding_manifest = verify_frozen_tree(arguments.embeddings, arguments.embeddings_sha256)
    outcome_manifest = verify_frozen_tree(arguments.outcomes, arguments.outcomes_sha256)
    fixture_manifest = verify_frozen_tree(arguments.fixture, arguments.fixture_sha256)
    if (
        embedding_manifest["metadata"].get("artifact_class")
        != "outcome_blind_caduceus_131k_embeddings"
        or embedding_manifest["metadata"].get("reporter_outcomes_read")
        or embedding_manifest["metadata"].get("downstream_head_fit")
        or outcome_manifest["metadata"].get("artifact_class")
        != "gse281364_replicate_safe_outcomes"
        or outcome_manifest["metadata"].get("donor_count") != 0
        or outcome_manifest["metadata"].get("outcome_role")
        != "exposed_development_MPRA_only"
        or outcome_manifest["metadata"].get("sealed_outcomes_loaded")
        or fixture_manifest["metadata"].get("artifact_class")
        != "gse281364_outcome_blind_caduceus_131k_fixture"
        or fixture_manifest["metadata"].get("outcomes_read")
    ):
        raise CommonLaneError("input role or firewall differs")
    fixture_rows, allele_features = load_fixture(arguments.fixture)
    projected, projection_audit = audit_and_load_features(arguments.embeddings, fixture_rows)
    ids = np.asarray([row["element_id"] for row in fixture_rows])
    groups = np.asarray([row["outer_locus_sequence_group_id"] for row in fixture_rows])
    folds = np.asarray([int(row["outer_fold"]) for row in fixture_rows], dtype=np.int64)
    outcomes = load_outcomes(
        arguments.outcomes / "outcomes/replicate_outcomes.tsv.gz", set(ids.tolist())
    )
    target = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in ids])
        for context in CONTEXTS
    }
    target_sd = {
        context: np.asarray([outcomes[(element, context)]["sd"] for element in ids])
        for context in CONTEXTS
    }
    arguments.output.mkdir(parents=True, mode=0o750)
    (arguments.output / "heads").mkdir()
    predictions: dict[tuple[str, str, str], np.ndarray] = {}
    selection_rows: list[dict[str, Any]] = []
    for context_index, context in enumerate(CONTEXTS):
        y = target[context]
        for baseline_id in ("zero", "training_mean"):
            oof = np.empty(len(ids), dtype=np.float64)
            for held_fold in range(5):
                test = folds == held_fold
                train = ~test
                oof[test] = 0.0 if baseline_id == "zero" else float(y[train].mean())
            predictions[("task_native_baseline", baseline_id, context)] = oof
        baseline_oof = np.empty(len(ids), dtype=np.float64)
        baseline_root = arguments.output / f"heads/task_native_baseline/{context}"
        baseline_root.mkdir(parents=True)
        for held_fold in range(5):
            test = folds == held_fold
            prediction, selection, state = fit_ridge_outer(
                allele_features, y, folds, held_fold
            )
            baseline_oof[test] = prediction
            np.savez_compressed(
                baseline_root / f"fold{held_fold}_allele_ridge.npz", **state
            )
            selection_rows.append(
                {
                    "model_id": "task_native_baseline",
                    "head_id": "allele_identity_ridge",
                    "context_id": context,
                    **selection,
                }
            )
        predictions[("task_native_baseline", "allele_identity_ridge", context)] = baseline_oof
        model_root = arguments.output / f"heads/{MODEL_ID}/{context}"
        model_root.mkdir(parents=True)
        for head_id in ("linear_ridge", "two_layer_gelu"):
            oof = np.empty(len(ids), dtype=np.float64)
            for held_fold in range(5):
                test = folds == held_fold
                values = projected[held_fold]
                if head_id == "linear_ridge":
                    prediction, selection, state = fit_ridge_outer(
                        values, y, folds, held_fold
                    )
                else:
                    prediction, selection, state = fit_mlp_outer(
                        values,
                        y,
                        folds,
                        held_fold,
                        seed=SEED + context_index * 10 + held_fold,
                    )
                oof[test] = prediction
                np.savez_compressed(model_root / f"fold{held_fold}_{head_id}.npz", **state)
                selection_rows.append(
                    {
                        "model_id": MODEL_ID,
                        "head_id": head_id,
                        "context_id": context,
                        **selection,
                    }
                )
            predictions[(MODEL_ID, head_id, context)] = oof

    prediction_rows: list[dict[str, Any]] = []
    for (model_id, head_id, context), prediction in sorted(predictions.items()):
        for index, element in enumerate(ids):
            prediction_rows.append(
                {
                    "model_id": model_id,
                    "head_id": head_id,
                    "context_id": context,
                    "element_id": element,
                    "outer_locus_sequence_group_id": groups[index],
                    "outer_fold": int(folds[index]),
                    "prediction": format(float(prediction[index]), ".17g"),
                    "observed_mean_signed_log2_activity_delta": format(
                        float(target[context][index]), ".17g"
                    ),
                    "observed_replicate_sd": format(float(target_sd[context][index]), ".17g"),
                    "experimental_replicates": 4,
                    "biological_donors": 0,
                    "outcome_role": "exposed_development_MPRA_only",
                    "champion_eligible": "false",
                }
            )
    aggregate_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []
    aggregate_lookup: dict[tuple[str, str, str], dict[str, float | None]] = {}
    for key, prediction in sorted(predictions.items()):
        model_id, head_id, context = key
        observed = target[context]
        metrics = metric_values(observed, prediction)
        aggregate_lookup[key] = metrics
        interval = bootstrap_metrics(
            observed,
            prediction,
            groups,
            seed=SEED + int(sha256("|".join(key).encode()).hexdigest()[:8], 16),
        )
        aggregate_rows.append(
            {
                "model_id": model_id,
                "head_id": head_id,
                "context_id": context,
                "elements": len(ids),
                "locus_sequence_groups": len(set(groups.tolist())),
                "pearson": _value(metrics["pearson"]),
                "pearson_ci_low": interval["pearson_ci_low"],
                "pearson_ci_high": interval["pearson_ci_high"],
                "spearman": _value(metrics["spearman"]),
                "spearman_ci_low": interval["spearman_ci_low"],
                "spearman_ci_high": interval["spearman_ci_high"],
                "rmse": _value(metrics["rmse"]),
                "rmse_ci_low": interval["rmse_ci_low"],
                "rmse_ci_high": interval["rmse_ci_high"],
                "mae": _value(metrics["mae"]),
                "r2": _value(metrics["r2"]),
                "calibration_intercept": _value(metrics["calibration_intercept"]),
                "calibration_slope": _value(metrics["calibration_slope"]),
                "license": LICENSE if model_id == MODEL_ID else "not_applicable",
                "restricted_comparator": "false",
                "open_champion_eligible_after_task_gates": str(model_id == MODEL_ID).lower(),
                "current_task_champion_eligible": "false",
                "multiple_testing": "none_descriptive_one_seed_smoke",
            }
        )
        for held_fold in range(5):
            selected = folds == held_fold
            fold_metric = metric_values(observed[selected], prediction[selected])
            fold_rows.append(
                {
                    "model_id": model_id,
                    "head_id": head_id,
                    "context_id": context,
                    "outer_fold": held_fold,
                    "elements": int(np.sum(selected)),
                    "locus_sequence_groups": int(np.sum(selected)),
                    **{name: _value(value) for name, value in fold_metric.items()},
                }
            )
    candidates = []
    for head_id in ("linear_ridge", "two_layer_gelu"):
        values = [aggregate_lookup[(MODEL_ID, head_id, context)] for context in CONTEXTS]
        candidates.append(
            (
                float(np.mean([float(value["spearman"]) for value in values])),
                -float(np.mean([float(value["rmse"]) for value in values])),
                head_id,
                values,
            )
        )
    mean_spearman, _, best_head, context_metrics = max(candidates)
    baseline_spearman = float(
        np.mean(
            [
                float(
                    aggregate_lookup[("task_native_baseline", "allele_identity_ridge", context)][
                        "spearman"
                    ]
                )
                for context in CONTEXTS
            ]
        )
    )
    positive_both = all(float(value["spearman"]) > 0 for value in context_metrics)
    mean_baseline_better = all(
        float(context_metrics[index]["rmse"])
        < float(aggregate_lookup[("task_native_baseline", "training_mean", context)]["rmse"])
        for index, context in enumerate(CONTEXTS)
    )
    recommended = (
        mean_spearman >= 0.10
        and mean_spearman - baseline_spearman >= 0.02
        and positive_both
        and mean_baseline_better
    )
    promotion_rows = [
        {
            "model_id": MODEL_ID,
            "selected_head": best_head,
            "mean_two_context_spearman": format(mean_spearman, ".17g"),
            "allele_baseline_mean_spearman": format(baseline_spearman, ".17g"),
            "spearman_gain": format(mean_spearman - baseline_spearman, ".17g"),
            "positive_spearman_both_contexts": str(positive_both).lower(),
            "rmse_better_than_training_mean_both_contexts": str(mean_baseline_better).lower(),
            "three_seed_screening_recommended": str(recommended).lower(),
            "promotion_gate": (
                "mean_spearman>=0.10;gain_over_allele_baseline>=0.02;"
                "positive_both_contexts;rmse_better_than_training_mean_both_contexts"
            ),
            "license": LICENSE,
            "restricted_comparator": "false",
            "open_champion_eligible_after_task_gates": "true",
            "current_task_champion_eligible": "false",
        }
    ]
    _write_gzip_tsv(arguments.output / "oof_predictions.tsv.gz", prediction_rows)
    _write_tsv(arguments.output / "aggregate_metrics.tsv", aggregate_rows)
    _write_tsv(arguments.output / "outer_fold_metrics.tsv", fold_rows)
    _write_tsv(arguments.output / "head_selection.tsv", selection_rows)
    _write_tsv(arguments.output / "promotion.tsv", promotion_rows)
    audit = {
        "schema_version": "masld-bench-gse281364-caduceus-input-audit-v1",
        "status": "pass_independent_frozen_and_projection_audit",
        "embedding_artifacts_sha256": arguments.embeddings_sha256,
        "outcome_artifacts_sha256": arguments.outcomes_sha256,
        "fixture_artifacts_sha256": arguments.fixture_sha256,
        "recursive_member_checksums_verified": True,
        "complete_markers_verified": True,
        "training_only_projection_statistics_rederived": True,
        "head_features_reconstructed": True,
        **projection_audit,
    }
    (arguments.output / "input_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": "masld-bench-gse281364-caduceus-common-lane-evaluation-v1",
        "status": "pass_descriptive_one_seed_development_smoke",
        "dataset_id": "gse281364",
        "models": [MODEL_ID],
        "contexts": list(CONTEXTS),
        "endpoint": "mean_signed_log2_ALT_minus_REF_RNA_over_DNA_activity",
        "elements": len(ids),
        "outer_locus_sequence_groups": len(set(groups.tolist())),
        "outer_folds": 5,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "outcome_role": "exposed_development_MPRA_only",
        "heads": ["linear_ridge", "two_layer_gelu"],
        "baselines": ["zero", "training_mean", "allele_identity_ridge"],
        "projection_fit_scope": "outer_training_locus_sequence_groups_only",
        "head_fit_scope": "outer_training_locus_sequence_groups_only",
        "hyperparameter_selection_scope": "one_prespecified_inner_outer_fold_within_outer_training_only",
        "held_out_outcomes_read_during_fit": False,
        "bootstrap_unit": "outer_locus_sequence_group",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed": SEED,
        "confirmatory_inference": False,
        "multiple_testing": "none_descriptive_one_seed_smoke",
        "sealed_outcomes_read": False,
        "external_evaluation": False,
        "champion_eligible": False,
        "clinical_claim_supported": False,
        "three_seed_screening_recommended": [MODEL_ID] if recommended else [],
        "restricted_comparator": [],
        "allowed_supervision": "locus-cross-fitted_assay-native_MPRA_activity_only",
        "cell_type_eQTL_or_ieQTL_head_fit": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--embeddings-sha256", required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--outcomes-sha256", required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--fixture-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    evaluate(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
