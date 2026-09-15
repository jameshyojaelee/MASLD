#!/usr/bin/env python3
"""Fail closed on bulk-RNA and proteomic transportability authorities."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path
from typing import Any

from masld_bench.artifacts import write_json_exclusive


class BulkProteinAuditError(RuntimeError):
    """Raised when a bulk/protein transportability authority differs."""


CENSUS_RELATIVE = (
    "config/artifacts/model-authorities/bulk_protein_transportability/census.json"
)
CENSUS_SHA256 = "b4b9221ce3d134bbedd15a8c197da59663bc637dc676846e0bc824e939955b4a"

EVIDENCE_SHA256 = {
    "executions/bulk-expansion-source-audit-21065762/ARTIFACTS.json": "883f870e2a76f78c39d257c24f8e795d2e74f1bd32aba785251b52aa7cb8229e",
    "executions/bulk-expansion-matrix-audit-21065853/ARTIFACTS.json": "3f934ab9c6ee730c6bd0887b10d963f8210f274c36dc0878cb814c85026dc038",
    "executions/gse267145-authoritative-join-21064930/ARTIFACTS.json": "e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580",
    "executions/gse267145-qc-rights-21065336/ARTIFACTS.json": "ac7c92fd01b2b4240b380e4f6a58d98f1a6f853f8d0ad98d7167d7af34baaa73",
    "executions/gse268273-activation-21066074/ARTIFACTS.json": "d9a9e8b1a99739908587f0b854bb3adc472a3543b769eeeb0813544e78d57c55",
    "executions/gse260666-activation-21066072/ARTIFACTS.json": "734d45aee49221ae9b72c2b79de6b4959ae1d3d01ac85420c89135c05469e8f7",
    "executions/model-data-080-21082249-gse274114/ARTIFACTS.json": "b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300",
    "executions/model-training-067-21080776/ARTIFACTS.json": "2b95408bc92719706b728fd795ae918bb03280ff99ff8c549965318fee5d5bf3",
    "executions/model-check-070-21082520/ARTIFACTS.json": "da2e5cb0cd6da1d7f7ecb15cf9e38cd0d11a7ea5a8a0e8c989cf4f1174b6f6ef",
    "config/artifacts/models/cell_bulk_classical_baselines/checkpoints.json": "886468e16d92e8df0249b79a978d00239f278f6d59fbb60f54ad637ba1f5f938",
    "config/artifacts/models/cell_bulk_classical_baselines/development_crosswalk.json": "e5b6a5c5369dd363a386f75bd5a1b3a4b8c8f09e73eb865d1a09c39363c0b465",
    "config/artifacts/models/cell_bulk_classical_baselines/exposure_audit.json": "0d95f120b4ca1ba9db92e365e04504b8c642f2980d47ee8a2b884b0f569abffc",
    "config/artifacts/models/mofaplus/checkpoints.json": "8b7727ebe9d6d3d612b85927b91664fe39939f8107aee05ea7e1a569bd32884e",
    "config/artifacts/models/mofaplus/development_crosswalk.json": "cdc768c1e2d0b959c158eea90c8e951e706401263d39ed0540d47d6e62523163",
    "config/artifacts/models/mofaplus/exposure_audit.json": "79aef2fe01c5ccaba339afd4b007ec74dad6be4d828535a4abc732808632172b",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BulkProteinAuditError(f"JSON authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise BulkProteinAuditError(f"expected JSON object: {path}")
    return value


def load_toml(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise BulkProteinAuditError(f"TOML authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise BulkProteinAuditError(f"expected TOML table: {path}")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BulkProteinAuditError(message)


def require_hash(root: Path, relative: str, expected: str) -> Path:
    path = root / relative
    if not path.is_file() or digest(path) != expected:
        raise BulkProteinAuditError(f"authority identity differs: {relative}")
    return path


def index_rows(value: dict[str, Any], key: str, field: str) -> dict[str, dict[str, Any]]:
    rows = value.get(key)
    if not isinstance(rows, list):
        raise BulkProteinAuditError(f"{key} rows differ")
    indexed = {
        str(row.get(field)): row for row in rows if isinstance(row, dict)
    }
    if len(indexed) != len(rows):
        raise BulkProteinAuditError(f"duplicate or malformed {key} rows")
    return indexed


def validate_census(census: dict[str, Any]) -> dict[str, dict[str, Any]]:
    expected = {
        "bulk_five_cohort",
        "gse267145_znf469_human_liver",
        "gse268273_imid_masld",
        "gse260666_bulk_rna",
        "gse274114_mash_hbv",
        "gse249997",
        "antwerp_inserm_shared",
        "gse31803_gse49541_fibrosis_array",
        "pxd051911",
        "gse267031",
    }
    cohorts = index_rows(census, "cohorts", "cohort_family_id")
    require(set(cohorts) == expected, "bulk/protein cohort census differs")

    allowed = set(census.get("global_contract", {}).get("allowed_missingness_states", []))
    require(
        allowed
        == {
            "observed",
            "structurally_missing",
            "not_applicable",
            "below_qc",
            "unavailable_permission",
            "join_unresolved",
            "withheld_sealed",
            "derivable_not_processed",
        },
        "missingness state roster differs",
    )
    for cohort_id, cohort in cohorts.items():
        metadata = cohort.get("metadata")
        require(isinstance(metadata, dict) and metadata, f"{cohort_id} metadata differs")
        for field, state in metadata.items():
            require(
                isinstance(state, dict) and state.get("state") in allowed,
                f"{cohort_id}:{field} missingness state differs",
            )

    require(
        cohorts["pxd051911"].get("pairing") == "same_donor_different_tissue"
        and cohorts["pxd051911"].get("cross_cohort_rna_pairing")
        == "same_study_unpaired",
        "PXD051911 topology differs",
    )
    require(
        all(
            value.get("state") == "withheld_sealed"
            for value in cohorts["gse267031"]["metadata"].values()
        ),
        "GSE267031 metadata escaped the seal",
    )
    gse268273 = cohorts["gse268273_imid_masld"]
    require(
        gse268273["technical_runs"] == 824
        and gse268273["metadata"]["fibrosis"]["complete_participants"] == 109
        and gse268273["metadata"]["fibrosis"]["evaluator_only"] is True
        and gse268273["metadata"]["nas"]["state"] == "structurally_missing",
        "GSE268273 participant or outcome firewall differs",
    )
    gse274114 = cohorts["gse274114_mash_hbv"]
    require(
        gse274114["metadata"]["group"]["evaluator_only"] is True
        and all(
            gse274114["metadata"][field]["state"] == "structurally_missing"
            for field in ("fibrosis", "nas", "sex", "age", "bmi")
        ),
        "GSE274114 participant metadata boundary differs",
    )

    lanes = index_rows(census, "task_native_model_lanes", "lane_id")
    require(
        set(lanes)
        == {
            "bulk_rna_stage_and_fibrosis",
            "observed_bulk_rna_h3k27ac",
            "protein_transport",
        }
        and all(lane.get("runnable_new_campaign") is False for lane in lanes.values()),
        "task-native lane disposition differs",
    )
    next_campaign = census.get("next_outcome_blind_campaign", {})
    require(
        next_campaign.get("campaign_id") == "gse268273_raw_reprocessing_v1"
        and next_campaign.get("model_or_outcome_selection_allowed") is False
        and next_campaign.get("labels_or_evaluator_metadata_available_to_workers") is False
        and next_campaign.get("new_submission_authorized_by_this_reconciliation") is False
        and next_campaign.get("additional_parallel_campaign_justified") is False,
        "next outcome-blind campaign boundary differs",
    )
    scope = census.get("scope", {})
    terminal = census.get("terminal_disposition", {})
    require(
        scope.get("sealed_sources_accessed") == []
        and scope.get("models_fit") == []
        and scope.get("models_scored") == []
        and scope.get("new_jobs_submitted") == []
        and terminal.get("new_model_campaigns_authorized") == 0
        and terminal.get("new_jobs_submitted") == 0
        and terminal.get("resource_authorities_modified") is False,
        "reconciliation silently authorized work or accessed sealed data",
    )
    return cohorts


def audit(root: Path) -> dict[str, Any]:
    census = load_json(require_hash(root, CENSUS_RELATIVE, CENSUS_SHA256))
    cohorts = validate_census(census)
    evidence = {
        relative: load_json(require_hash(root, relative, expected))
        for relative, expected in EVIDENCE_SHA256.items()
    }

    datasets = {
        name: load_toml(root / f"config/datasets/{name}.toml")
        for name in ("bulk_five_cohort", "gse249997", "gse267031", "pxd051911")
    }
    gse274_dataset = load_toml(root / "config/datasets/gse274114_mash_hbv.toml")
    bulk_task = load_toml(root / "config/tasks/bulk_state_transfer.toml")
    histology_task = load_toml(
        root / "config/evaluation/gse267145_histology_state_task.toml"
    )
    gse268_task = load_toml(
        root / "config/evaluation/gse268273_fibrosis_ood_transfer_task.toml"
    )
    gse274_task = load_toml(
        root / "config/evaluation/gse274114_within_instrument_etiology_task.toml"
    )
    expansion = load_toml(root / "config/evaluation/cross_cohort_expansion.toml")
    tournament = load_toml(root / "config/evaluation/family_native_tournament.toml")
    sealed = load_toml(root / "config/evaluation/sealed_protocol.toml")
    promotion = load_toml(root / "config/evaluation/promotion_gates.toml")
    mandatory = load_toml(root / "config/models/mandatory_baselines.toml")
    multimodal = load_toml(root / "config/models/multimodal_integration.toml")
    histology_surface = load_json(
        root / "config/evaluation/gse267145_histology_benchmark_surface.json"
    )

    require(
        datasets["bulk_five_cohort"].get("expected_biological_units") == 844
        and datasets["bulk_five_cohort"].get("biological_unit") == "donor"
        and datasets["bulk_five_cohort"].get("admission_blocking") is True,
        "five-cohort bulk authority differs",
    )
    require(
        datasets["gse249997"].get("expected_biological_units") == 77
        and datasets["gse249997"].get("biological_unit") == "participant"
        and any("70 participants" in value for value in datasets["gse249997"].get("blockers", [])),
        "GSE249997 sample-participant boundary differs",
    )
    require(
        datasets["gse267031"].get("status") == "withheld_sealed"
        and datasets["gse267031"].get("outcome_url") == "WITHHELD_SEALED"
        and datasets["gse267031"].get("expected_biological_units") == 18,
        "GSE267031 seal differs",
    )
    require(
        datasets["pxd051911"].get("expected_biological_units") == 58
        and datasets["pxd051911"].get("pairing_levels") == ["same_donor_different_tissue"]
        and datasets["pxd051911"].get("admission_blocking") is True,
        "PXD051911 dataset boundary differs",
    )
    require(
        gse274_dataset.get("expected_biological_units") == 39
        and gse274_dataset.get("biological_unit") == "participant"
        and gse274_dataset.get("admission_blocking") is True,
        "GSE274114 dataset boundary differs",
    )

    require(
        bulk_task.get("datasets_sealed") == ["gse267031"]
        and "pxd051911" in bulk_task.get("datasets_development", [])
        and bulk_task.get("resampling_units") == ["donor"]
        and "prospective 80% power" in bulk_task.get("claim_gate", "")
        and "cross-sectional stage-associated remodeling" in bulk_task.get("claim_gate", ""),
        "bulk TaskSpec seal, unit, or claim gate differs",
    )
    require(
        histology_task.get("unit_of_inference") == "participant"
        and histology_task.get("datasets_sealed") == []
        and histology_task.get("evaluator_parameters", {}).get("participants") == 99
        and histology_task.get("evaluator_parameters", {}).get("fibrosis_scale")
        == [0, 1, 2, 3]
        and "No participant-level age or BMI" in " ".join(histology_task.get("admission_gates", [])),
        "GSE267145 histology TaskSpec differs",
    )
    require(
        gse268_task.get("status") == "blocked"
        and gse268_task.get("unit_of_inference") == "participant"
        and gse268_task.get("datasets_sealed") == []
        and gse268_task.get("evaluator_parameters", {}).get("participants") == 109
        and gse268_task.get("evaluator_parameters", {}).get("raw_preprocessing_complete") is False
        and gse268_task.get("evaluator_parameters", {}).get("external_outcomes_available_to_model") is False,
        "GSE268273 TaskSpec firewall differs",
    )
    require(
        gse274_task.get("unit_of_inference") == "participant"
        and gse274_task.get("datasets_sealed") == []
        and gse274_task.get("evaluator_parameters", {}).get("participants") == 39
        and gse274_task.get("evaluator_parameters", {}).get("four_class_performance_allowed") is False
        and gse274_task.get("evaluator_parameters", {}).get("participant_metadata_inputs_available") is False,
        "GSE274114 within-instrument TaskSpec differs",
    )

    expansion_rows = index_rows(expansion, "cohort_family", "family_id")
    for cohort_id in cohorts:
        if cohort_id in {"bulk_five_cohort", "gse267031", "gse249997"}:
            continue
        require(cohort_id in expansion_rows, f"{cohort_id} disappeared from expansion registry")
    require(
        expansion.get("missing_modality_encoded_as_zero") is False
        and expansion.get("same_person_all_modalities_same_outer_fold") is True
        and expansion_rows["gse268273_imid_masld"].get("technical_runs") == 824
        and expansion_rows["gse274114_mash_hbv"].get("reported_biological_units") == 39
        and expansion_rows["pxd051911"].get("topology") == "same_donor_different_tissue",
        "cross-cohort topology or missingness boundary differs",
    )

    protein_task = tournament.get("task", {}).get("protein_transport", {})
    require(
        protein_task.get("status") == "planned_registration"
        and protein_task.get("candidate_dataset_families") == ["pxd051911"]
        and protein_task.get("rna_protein_pair_assumed") is False
        and protein_task.get("champion_claim_requires_external_protein_cohort") is True,
        "protein tournament boundary differs",
    )
    require(
        sealed.get("gse267031", {}).get("activation")
        == [
            "donor_deduplication_pass",
            "stage_metadata_complete",
            "endpoint_complete",
            "prospective_power_at_least_0.80",
        ],
        "GSE267031 activation gates differ",
    )
    bulk_gate = promotion.get("tasks", {}).get("bulk_state_transfer", {})
    require(
        bulk_gate.get("sealed_holdout_dataset_ids") == ["gse267031"]
        and bulk_gate.get("delta_min") == 0.05
        and bulk_gate.get("power_minimum_effect") == 0.05,
        "bulk promotion effect or seal differs",
    )

    model_rows = index_rows(mandatory, "models", "model_id")
    for model_id in (
        "hvg_pca_nearest_centroid",
        "hvg_pca_elastic_net",
        "hvg_pca_linear_svm",
        "mean_expression_bulk",
        "assay_native_pseudobulk",
    ):
        require(model_id in model_rows, f"mandatory bulk baseline disappeared: {model_id}")
    require(
        model_rows["mean_expression_bulk"].get("admission_blocking") is True
        and "bulk_state_transfer"
        not in model_rows["hvg_pca_elastic_net"].get("execution", {}).get("executable_tasks", [])
        and "bulk_state_transfer"
        not in model_rows["hvg_pca_linear_svm"].get("execution", {}).get("executable_tasks", [])
        and "bulk_state_transfer"
        not in model_rows["assay_native_pseudobulk"].get("execution", {}).get("executable_tasks", []),
        "an unresolved baseline silently became bulk-executable",
    )
    multimodal_rows = index_rows(multimodal, "models", "model_id")
    require(
        multimodal_rows["mofaplus"].get("status") == "deferred"
        and multimodal_rows["mofaplus"].get("supported_tasks") == []
        and multimodal_rows["mofaplus"].get("admission_blocking") is True,
        "MOFA+ inductive boundary differs",
    )

    join = load_json(root / "executions/gse267145-authoritative-join-21064930/join_evidence.json")
    qc = load_json(root / "executions/gse267145-qc-rights-21065336/qc_summary.json")
    gse268_activation = load_json(root / "executions/gse268273-activation-21066074/activation_audit.json")
    gse260_activation = load_json(root / "executions/gse260666-activation-21066072/activation_audit.json")
    gse274_activation = load_json(root / "executions/model-data-080-21082249-gse274114/activation_audit.json")
    preflight = load_json(root / "executions/model-training-067-21080776/preflight_receipt.json")
    diagnostic = load_json(root / "executions/model-check-070-21082520/units/outer_1--seed_1701/fit/receipt.json")
    require(
        join.get("participants_joined") == 99
        and join.get("pairing") == "same_sample_different_aliquot"
        and qc.get("sex_counts") == {"F": 85, "M": 14}
        and qc.get("hard_qc_failures") == 0,
        "GSE267145 participant evidence differs",
    )
    require(
        gse268_activation.get("participants") == 109
        and gse268_activation.get("fibrosis_counts")
        == {"F0": 36, "F1": 41, "F2": 13, "F3": 14, "F4": 5}
        and gse268_activation.get("model_training_activated") is False,
        "GSE268273 activation evidence differs",
    )
    require(
        gse260_activation.get("participants") == 16
        and gse260_activation.get("model_training_activated") is False,
        "GSE260666 activation evidence differs",
    )
    require(
        gse274_activation.get("participants") == 39
        and gse274_activation.get("model_fit_or_scoring_run") is False
        and gse274_activation.get("four_class_performance_allowed") is False
        and gse274_activation.get("group_aggregate_metadata_descriptive_only") is True,
        "GSE274114 no-fit evidence differs",
    )
    require(
        histology_surface.get("production_fit_authorized") is False
        and histology_surface.get("inputs", {}).get("evaluator_outcomes", {}).get("model_environment_access") is False
        and preflight.get("status") == "passed"
        and preflight.get("metrics_calculated") is False
        and preflight.get("production_predictions_complete") is False
        and diagnostic.get("outer_test_metrics_calculated") is False
        and diagnostic.get("outer_test_outcomes_read") is False
        and diagnostic.get("recorded_sex_read") is False,
        "GSE267145 unscored preflight boundary differs",
    )

    baseline_crosswalk = evidence[
        "config/artifacts/models/cell_bulk_classical_baselines/development_crosswalk.json"
    ]
    require(
        baseline_crosswalk.get("open_champion_contract", {}).get("bulk", "").startswith(
            "No baseline in this bundle can currently earn an open bulk champion"
        ),
        "classical baseline bulk admission summary differs",
    )
    mofa = evidence["config/artifacts/models/mofaplus/checkpoints.json"]
    require(
        "No champion-eligible unseen-sample transform" in mofa.get("inductive_inference_contract", {}).get("ordinary_mofa", ""),
        "MOFA+ unseen-participant limitation differs",
    )

    return {
        "schema_version": "masld-bench-bulk-protein-transportability-audit-v1",
        "status": "pass_fail_closed",
        "cohorts_reconciled": sorted(cohorts),
        "cohort_count": len(cohorts),
        "task_native_lanes": [
            "bulk_rna_stage_and_fibrosis",
            "observed_bulk_rna_h3k27ac",
            "protein_transport",
        ],
        "next_outcome_blind_campaign": "gse268273_raw_reprocessing_v1",
        "additional_parallel_campaign_justified": False,
        "submission_authorized": False,
        "models_fit": [],
        "models_scored": [],
        "sealed_sources_accessed": [],
        "bulk_champion_created": False,
        "protein_champion_created": False,
        "resource_authorities_modified": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.project_root.resolve())
    args.output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(args.output / "bulk_protein_transportability_audit.json", result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
