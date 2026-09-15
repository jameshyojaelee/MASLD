#!/usr/bin/env python3
"""Register the GSE49541 prediction requirements from the frozen source fit.

The candidate model named here is the one the source fit already chose on
GSE267145 evidence alone.  Writing the requirements from the frozen fit receipt
rather than by hand means the candidate cannot be quietly re-picked later, and
the file's own SHA-256 is what the predictor verifies before it opens any
external expression.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--source-fit", required=True, type=Path)
    parser.add_argument("--source-fit-contract", required=True, type=Path)
    parser.add_argument("--activation", required=True, type=Path)
    parser.add_argument("--external-expression", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    root = arguments.benchmark_root.resolve(strict=True)
    fit_receipt = json.loads(
        (arguments.source_fit / "fit_receipt.json").read_text(encoding="utf-8")
    )
    fit_contract = json.loads(
        arguments.source_fit_contract.read_text(encoding="utf-8")
    )
    axis_receipt = json.loads(
        (arguments.activation / "gene_axis_receipt.json").read_text(encoding="utf-8")
    )
    verdict = json.loads(
        (arguments.external_expression / "summarization_verdict.json").read_text(
            encoding="utf-8"
        )
    )
    contract = {
        "schema_version": "masld-bench-gse49541-fibrosis-transfer-prediction-v1",
        "status": "registered_prediction_only_pending",
        "prediction_id": "gse49541_fibrosis_transfer_prediction_v1",
        "task_id": "gse49541_fibrosis_transfer",
        "arm_id": fit_receipt["arm_id"],
        "gate_eligible_arm": fit_receipt["gate_eligible_arm"],
        "selected_source_model_id": fit_receipt["selected_source_model_id"],
        "selected_source_model_selection_basis": fit_receipt[
            "selected_source_model_selection_basis"
        ],
        "source_fit": {
            "path": str(arguments.source_fit.resolve().relative_to(root)),
            "artifacts_sha256": sha256_file(arguments.source_fit / "ARTIFACTS.json"),
            "required_status": fit_receipt["status"],
            "model_seeds": fit_receipt["model_seeds"],
            "fitting_participants": fit_receipt["fitting_participants"],
            "training_prevalence_advanced_f3_f4": fit_receipt[
                "training_prevalence_advanced_f3_f4"
            ],
            "contract_sha256": sha256_file(arguments.source_fit_contract),
        },
        "activation": {
            "path": str(arguments.activation.resolve().relative_to(root)),
            "artifacts_sha256": sha256_file(arguments.activation / "ARTIFACTS.json"),
            "shared_genes": axis_receipt["shared_genes"],
        },
        "external_expression": {
            "path": str(arguments.external_expression.resolve().relative_to(root)),
            "artifacts_sha256": sha256_file(
                arguments.external_expression / "ARTIFACTS.json"
            ),
            "required_status": verdict["status"],
            "participants": verdict["arrays_summarized"],
            "gene_matrix": "gse49541_gpl570_gene_matrix.tsv",
            "key_field": "ensembl_gene_id",
            "preprocessing": verdict["method"],
            "arrays_read_per_call": verdict["arrays_read_per_frma_call"],
            "labels_included": False,
        },
        "models": fit_contract["models"],
        "fitted_model_ids": fit_receipt["fitted_model_ids"],
        "unfittable_model_ids": fit_receipt["unfittable_model_ids"],
        "unfittable_models": fit_receipt["unfittable_models"],
        "external_transform": {
            "per_array_rank": "average_rank_within_one_array_over_the_shared_axis_scaled_to_the_open_unit_interval",
            "gene_median": "deposited_log2_gpl570_intensity_minus_that_same_array_median_over_the_shared_axis",
            "external_values_are_already_log2": True,
            "cross_array_pooling": False,
        },
        "prediction_aggregation": "mean_probability_across_five_source_model_seeds",
        "seed_level_results_emitted": True,
        "external_labels_read": False,
        "external_participant_table_read": False,
        "query_fit_or_calibration_allowed": False,
        "query_joint_normalization_allowed": False,
        "external_metrics_allowed": False,
        "prediction_frozen_before_evaluator_label_join": True,
        "project_sealed": False,
        "champion_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
        "source_well_posedness": {
            "training_partition_linearly_separable": fit_receipt[
                "training_partition_linearly_separable"
            ],
            "eligible_grid_configurations": fit_receipt["eligible_grid_configurations"],
            "models_with_an_isolated_surviving_configuration": fit_receipt[
                "models_with_an_isolated_surviving_configuration"
            ],
            "no_source_model_exceeds_its_own_null": fit_receipt[
                "no_source_model_exceeds_its_own_null"
            ],
            "low_positive_count_warning": fit_receipt["low_positive_count_warning"],
        },
    }
    with arguments.output.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"contract_sha256": sha256_file(arguments.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
