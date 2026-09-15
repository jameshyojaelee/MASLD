#!/usr/bin/env python3
"""Freeze the donor-safe GSE260666 activation and preprocessing selection record."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence
import tomllib


class ExternalHistologyActivationError(RuntimeError):
    """Raised when external-development activation would violate its selection record."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExternalHistologyActivationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def canonical_parameters(value: Mapping[str, Any]) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def lock_modal_recipe(parameter_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if len(parameter_rows) != 25:
        raise ExternalHistologyActivationError("source primary recipe census differs")
    counts = Counter(canonical_parameters(row) for row in parameter_rows)
    maximum = max(counts.values())
    tied = [json.loads(value) for value, count in counts.items() if count == maximum]
    selected = min(
        tied,
        key=lambda row: (
            int(row["feature_request"]),
            int(row["pca_components"]),
            float(row["c"]),
            canonical_parameters(row),
        ),
    )
    return {
        "parameters": selected,
        "selection_count": maximum,
        "distinct_recipe_count": len(counts),
        "tied_modal_recipe_count": len(tied),
        "distribution": [
            {"parameters": json.loads(value), "count": count}
            for value, count in sorted(
                counts.items(), key=lambda item: (-item[1], item[0])
            )
        ],
    }


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise ExternalHistologyActivationError(f"{label} SHA-256 differs")


def audit_activation(
    *, benchmark_root: Path, contract_path: Path, contract_sha256: str, output: Path
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise ExternalHistologyActivationError(
            f"refusing to overwrite external activation: {output}"
        )
    _check_hash(contract_path, contract_sha256, "activation contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_activation_audit_pending"
        or contract.get("evaluation_role")
        != "small_project_exposed_external_development"
        or contract.get("claim_boundary", {}).get("project_sealed") is not False
    ):
        raise ExternalHistologyActivationError("activation identity differs")

    selection_spec = contract["source_selection"]
    selection_root = benchmark_root / selection_spec["selection_path"]
    _check_hash(
        selection_root / "ARTIFACTS.json",
        selection_spec["selection_artifacts_sha256"],
        "source selection",
    )
    verify_frozen_tree(selection_root)
    selection = json.loads(
        (selection_root / "selection_receipt.json").read_text(encoding="utf-8")
    )
    if (
        selection.get("status") != selection_spec["required_status"]
        or selection.get("selected_model_id") != contract["selected_source_model_id"]
        or selection.get("secondary_endpoint_values_used_for_selection") is not False
        or selection.get("external_features_or_labels_accessed") is not False
    ):
        raise ExternalHistologyActivationError("frozen RNA selection differs")

    source = contract["source_development"]
    production = benchmark_root / source["production_path"]
    molecular = benchmark_root / source["molecular_path"]
    for root, expected, label in (
        (production, source["production_artifacts_sha256"], "source production"),
        (molecular, source["molecular_artifacts_sha256"], "source molecular"),
    ):
        _check_hash(root / "ARTIFACTS.json", expected, label)
        verify_frozen_tree(root)

    external = contract["external_fixture"]
    fixture_campaign = benchmark_root / external["fixture_campaign_path"]
    fixture = benchmark_root / external["path"]
    model_input = benchmark_root / external["model_input_path"]
    for root, expected, label in (
        (fixture_campaign, external["fixture_campaign_artifacts_sha256"], "fixture campaign"),
        (fixture, external["artifacts_sha256"], "fixture"),
        (model_input, external["model_input_artifacts_sha256"], "external model input"),
    ):
        _check_hash(root / "ARTIFACTS.json", expected, label)
    verify_frozen_tree(model_input)
    _check_hash(
        fixture_campaign / "source.sha256",
        external["fixture_source_manifest_sha256"],
        "fixture source manifest",
    )
    prospective = contract["prospective_endpoint_harmonization"]
    task_path = benchmark_root / prospective["task_spec_path"]
    _check_hash(task_path, prospective["task_spec_sha256"], "external task spec")
    source_manifest = (fixture_campaign / "source.sha256").read_text(encoding="utf-8")
    if (
        f'{prospective["task_spec_sha256"]}  {task_path.as_posix()}'
        not in source_manifest
    ):
        raise ExternalHistologyActivationError(
            "prospective task spec was not frozen in the pre-selection fixture"
        )
    with task_path.open("rb") as handle:
        task = tomllib.load(handle)
    if (
        task.get("source_to_evaluation_mapping")
        != prospective["source_to_evaluation_mapping"]
        or task.get("datasets_sealed") != []
        or "not project-sealed" not in task.get("claim_gate", "")
    ):
        raise ExternalHistologyActivationError("endpoint harmonization differs")

    _, source_features = read_tsv(molecular / "rna_feature_axis.tsv")
    _, external_features = read_tsv(model_input / "rna_feature_axis.tsv")
    source_ids = [row["stable_gene_id"] for row in source_features]
    external_ids = [row["stable_gene_id"] for row in external_features]
    if (
        len(source_ids) != len(set(source_ids))
        or len(external_ids) != len(set(external_ids))
    ):
        raise ExternalHistologyActivationError("stable-gene axis is duplicated")
    common = sorted(set(source_ids) & set(external_ids))
    preprocessing = contract["preprocessing_lock"]
    if len(common) != preprocessing["expected_common_stable_genes"]:
        raise ExternalHistologyActivationError("common stable-gene census differs")

    model_id = contract["selected_source_model_id"]
    parameter_rows: list[Mapping[str, Any]] = []
    seeds = tuple(int(value) for value in source["model_seeds"])
    for outer_fold in range(5):
        for seed in seeds:
            unit = production / f"units/outer_{outer_fold}/seed_{seed}"
            audit_path = (
                unit
                / f"fit/fit_receipts/outer_{outer_fold}/seed_{seed}"
                / f"candidate_audit--{model_id}--stage3.json"
            )
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
            selected = audit.get("selected", {})
            if (
                audit.get("model_id") != model_id
                or audit.get("endpoint_id") != "stage3"
                or "one_standard_error_threshold" not in selected
                or selected.get("outer_training_refit_valid") is not True
            ):
                raise ExternalHistologyActivationError(
                    "source stage3 one-SE recipe audit differs"
                )
            parameter_rows.append(dict(selected["parameters"]))
    modal = lock_modal_recipe(parameter_rows)
    if modal["parameters"] != preprocessing["expected_locked_parameters"]:
        raise ExternalHistologyActivationError("locked full-fit recipe differs")

    output.mkdir(parents=True)
    write_tsv(
        output / "common_stable_gene_axis.tsv",
        ("feature_index", "stable_gene_id", "observability_state"),
        [
            {
                "feature_index": index,
                "stable_gene_id": stable_id,
                "observability_state": "observed_in_both_source_and_external_assay_axes",
            }
            for index, stable_id in enumerate(common)
        ],
    )
    distribution_rows = [
        {
            "recipe_rank": index + 1,
            "count": row["count"],
            "parameters_json": canonical_parameters(row["parameters"]),
            "locked_full_fit_recipe": str(
                row["parameters"] == modal["parameters"]
            ).lower(),
        }
        for index, row in enumerate(modal["distribution"])
    ]
    write_tsv(
        output / "source_stage3_recipe_distribution.tsv",
        tuple(distribution_rows[0]),
        distribution_rows,
    )
    receipt = {
        "schema_version": "masld-bench-gse260666-external-histology-activation-audit-v1",
        "status": "passed_donor_safe_external_development_activation",
        "activation_id": contract["activation_id"],
        "evaluation_role": contract["evaluation_role"],
        "selected_source_model_id": model_id,
        "selection_artifacts_sha256": selection_spec["selection_artifacts_sha256"],
        "source_participants": source["participants"],
        "external_participants": external["participants"],
        "common_stable_genes": len(common),
        "locked_parameters": modal["parameters"],
        "locked_recipe_support_units": modal["selection_count"],
        "source_stage3_recipe_units_audited": len(parameter_rows),
        "source_feature_values_read": False,
        "external_feature_values_read": False,
        "external_labels_read": False,
        "secondary_endpoints_used": False,
        "endpoint_mapping_registered_before_source_selection": True,
        "outcome_driven_label_remapping_allowed": False,
        "nafl_nash_pooled": False,
        "full_source_fit_authorized": True,
        "prediction_only_external_label_access_allowed": False,
        "project_sealed": False,
        "champion_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
        "confirmatory_inference_allowed": False,
    }
    with (output / "activation_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = audit_activation(
        benchmark_root=arguments.benchmark_root,
        contract_path=arguments.contract,
        contract_sha256=arguments.contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
