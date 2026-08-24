from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Mapping
from unittest.mock import patch

from masld_bench.artifacts import (
    canonical_hash,
    freeze_tree,
    sha256_file,
    write_json_exclusive,
    write_text_exclusive,
)
from masld_bench.campaign import execute_run, verify_run_execution_attempt
from masld_bench.contracts import ArtifactRef, MissingState, PredictionBundle
from masld_bench.evaluators.endpoint_power import recompute_endpoint_bootstrap
from masld_bench.firewall import (
    SEALED_OUTCOME_BUNDLE_SCHEMA_VERSION,
    PowerDecision,
    PredictionCommit,
    commit_predictions,
    freeze_development_power_evidence,
    freeze_power_decision,
    resource_snapshot,
)
from masld_bench.planner import freeze_campaign, load_frozen_plan
from masld_bench.selection import freeze_selection_lock, verify_selection_lock
from masld_bench.tournament import (
    CHAMPION_GATE_THRESHOLDS,
    freeze_development_outcome_bundle,
    freeze_champion_gate_decision,
    freeze_finalist_campaign_universe,
    freeze_finalist_candidate_ledger,
    freeze_finalist_development_metric_bundle,
    freeze_scientific_run_receipt,
    freeze_terminal_evaluation_authorization,
    freeze_variant_secondary_run_receipt,
    verify_finalist_campaign_universe,
    verify_finalist_development_metric_bundle,
    verify_scientific_run_receipt,
    verify_selection_candidate_ledger,
    verify_variant_secondary_run_receipt,
)


TASK_IDS = ("cell_state_mapping", "variant_to_regulation")
SEEDS = (1103, 2909, 4721, 6673, 8111)
SELECTED_MODEL = "foundation_fixture"
BASELINE_MODEL = "baseline_fixture"
SELECTED_REGIME = "full_finetune"
BASELINE_REGIME = "fixed_baseline"
SECONDARY_COMPARATOR_MODELS = ("abc", "nearest_gene", "re2g")
SECONDARY_COMPARATOR_SEED = SEEDS[0]
SECONDARY_COMPARATOR_REGIME = "fixed_secondary_baseline"
SEALED_DATASET = "gse289173"
PREDICTION_FIRST_STRESS_DATASET = "stress_fixture"

_CELL_CLASSES = ("cholangiocyte", "hepatocyte", "stellate")
_VARIANT_LINEAGES = ("hepatocyte", "immune", "stellate")
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
_FIXTURE_RESOURCE_FIREWALL = resource_snapshot(
    PACKAGE_ROOT.parents[1],
    ("Analysis/MASLD_Model_Benchmark/pyproject.toml",),
)


def _fixture_resource_snapshot(*_args: Any, **_kwargs: Any) -> dict[str, Any]:
    return _FIXTURE_RESOURCE_FIREWALL


def _promotion_thresholds(task_id: str) -> dict[str, Any]:
    common: dict[str, Any] = {
        "promotion_mode": "sealed_confirmatory",
        "confirmatory_family_id": f"{task_id}:primary",
        "confirmatory_fwer": 0.05,
        "sealed_outcome_bundle_id": "gse289173_joint_cell_and_regulatory_bundle",
        "power": {
            "metric_id": (
                "delta_donor_class_balanced_macro_f1"
                if task_id == "cell_state_mapping"
                else "delta_mean_fisher_z_celltype_spearman"
            ),
            "minimum_effect": 0.02 if task_id == "cell_state_mapping" else 0.03,
            "alpha": 0.05,
            "target_power": 0.80,
            "n_primary_claims": 5,
            "min_units": 10,
            "two_sided": True,
        },
    }
    if task_id == "cell_state_mapping":
        common.update(
            {
                "primary_min": 0.70,
                "delta_min": 0.02,
                "delta_ci_low_strict_min": 0.0,
                "stratum_min": 0.50,
                "brier_delta_max": 0.01,
                "required_class_ids": list(_CELL_CLASSES),
            }
        )
    else:
        common.update(
            {
                "delta_min": 0.03,
                "delta_ci_low_strict_min": 0.0,
                "stratum_min": -0.03,
                "required_major_lineage_ids": list(_VARIANT_LINEAGES),
                "eqtl_auprc_noninferiority_margin": 0.01,
            }
        )
    return common


def digest(label: str) -> str:
    return canonical_hash({"fixture": label})


def _variant_capability_record(
    model_id: str, *, role: str, is_mandatory_baseline: bool
) -> dict[str, Any]:
    return {
        "model_id": model_id,
        "role": role,
        "native_outputs": ["gene_expression_delta"],
        "allowed_endpoints": [
            "allelic_direction",
            "eqtl_retrieval",
            "signed_cell_type_eqtl_effect",
        ],
        "primary_eligible": True,
        "requires_fitted_head": True,
        "requires_observed_target_context": False,
        "is_mandatory_baseline": is_mandatory_baseline,
    }


def _secondary_variant_capability_record(model_id: str) -> dict[str, Any]:
    if model_id not in SECONDARY_COMPARATOR_MODELS:
        raise ValueError(f"unsupported secondary comparator: {model_id}")
    return {
        "model_id": model_id,
        "role": "baseline" if model_id == "nearest_gene" else "link_only",
        "native_outputs": ["enhancer_gene_link_score"],
        "allowed_endpoints": ["enhancer_gene_link", "eqtl_retrieval"],
        "primary_eligible": False,
        "requires_fitted_head": False,
        "requires_observed_target_context": False,
        "is_mandatory_baseline": True,
    }


def _simple_candidate_run_identity(
    *, task_id: str, model_id: str, seed: int, regime: str
) -> dict[str, Any]:
    identity: dict[str, Any] = {
        "schema_version": "masld-bench-run-v1",
        "campaign_id": "selection-fixture",
        "wave": "selection",
        "model_id": model_id,
        "task_id": task_id,
        "fold": "outer_pending",
        "seed": seed,
        "adaptation_regime": regime,
        "resource_profile": "fixture",
        "runtime_id": "fixture-runtime",
        "runtime_registry_sha256": digest("runtime"),
        "model_registry_sha256": digest(f"model-{model_id}"),
        "task_registry_sha256": digest(f"task-{task_id}"),
        "immutable_inputs": {
            "development_fixture": digest("development-activation")
        },
    }
    if task_id == "variant_to_regulation":
        if model_id in SECONDARY_COMPARATOR_MODELS:
            capability = _secondary_variant_capability_record(model_id)
            primary_scoring_allowed = False
        else:
            capability = _variant_capability_record(
                model_id,
                role=(
                    "primary_candidate"
                    if model_id == SELECTED_MODEL
                    else "baseline"
                ),
                is_mandatory_baseline=model_id == BASELINE_MODEL,
            )
            primary_scoring_allowed = True
        identity["metadata"] = {
            "primary_evaluator_id": (
                "variant_ld_block_fisher_z_spearman_gain_v1"
            ),
            "variant_capability": capability,
            "variant_capability_sha256": canonical_hash(capability),
            "variant_capability_registry_sha256": digest(
                "variant-capability-registry"
            ),
            "primary_endpoint_scoring_allowed": primary_scoring_allowed,
        }
    return identity


def _simple_secondary_run_id(model_id: str) -> str:
    return canonical_hash(
        _simple_candidate_run_identity(
            task_id="variant_to_regulation",
            model_id=model_id,
            seed=SECONDARY_COMPARATOR_SEED,
            regime=SECONDARY_COMPARATOR_REGIME,
        )
    )


def _variant_capability_binding(
    selected_model_id: str = SELECTED_MODEL,
    baseline_model_id: str = BASELINE_MODEL,
) -> dict[str, Any]:
    selected = _variant_capability_record(
        selected_model_id,
        role="primary_candidate",
        is_mandatory_baseline=False,
    )
    baseline = _variant_capability_record(
        baseline_model_id,
        role="baseline",
        is_mandatory_baseline=True,
    )
    return {
        "schema_version": "masld-bench-locked-variant-primary-capability-v1",
        "task_id": "variant_to_regulation",
        "primary_endpoint_id": "signed_cell_type_eqtl_effect",
        "primary_evaluator_id": "variant_ld_block_fisher_z_spearman_gain_v1",
        "capability_registry_sha256": digest("variant-capability-registry"),
        "selected": {
            "model_id": selected_model_id,
            "capability": selected,
            "capability_sha256": canonical_hash(selected),
        },
        "baseline": {
            "model_id": baseline_model_id,
            "capability": baseline,
            "capability_sha256": canonical_hash(baseline),
        },
    }


@dataclass(frozen=True, slots=True)
class ScientificSelectionFixture:
    candidate_dir: Path
    selection_lock_dir: Path
    selection_candidate_ledger_dir: Path
    finalist_campaign_universe_dir: Path
    finalist_metric_bundle_dirs: dict[str, Path]
    execution_attempt_dirs: dict[str, Path]
    development_prediction_bundle_paths: dict[str, Path]
    scientific_receipt_dirs: tuple[Path, ...]
    variant_secondary_receipt_dirs: tuple[Path, ...]
    development_outcome_bundle_dirs: dict[str, Path]


_SCIENTIFIC_SELECTION_CACHE_DIR: tempfile.TemporaryDirectory[str] | None = None
_SCIENTIFIC_SELECTION_CACHE: ScientificSelectionFixture | None = None


def create_candidate(root: Path) -> tuple[Path, dict[tuple[str, str, int], str]]:
    candidate = root / "candidate"
    candidate.mkdir()
    development_registry_sha256 = digest("development-registry")
    development_activation_sha256 = digest("development-activation")
    sealed_registry_sha256 = digest("sealed-registry")
    runs: list[dict[str, Any]] = []
    run_ids: dict[tuple[str, str, int], str] = {}
    for task_id in TASK_IDS:
        for model_id, regime in (
            (SELECTED_MODEL, SELECTED_REGIME),
            (BASELINE_MODEL, BASELINE_REGIME),
        ):
            for seed in SEEDS:
                identity = _simple_candidate_run_identity(
                    task_id=task_id,
                    model_id=model_id,
                    seed=seed,
                    regime=regime,
                )
                run_id = canonical_hash(identity)
                runs.append({**identity, "run_id": run_id})
                run_ids[(task_id, model_id, seed)] = run_id
    for model_id in SECONDARY_COMPARATOR_MODELS:
        identity = _simple_candidate_run_identity(
            task_id="variant_to_regulation",
            model_id=model_id,
            seed=SECONDARY_COMPARATOR_SEED,
            regime=SECONDARY_COMPARATOR_REGIME,
        )
        run_id = canonical_hash(identity)
        runs.append({**identity, "run_id": run_id})
        run_ids[("variant_to_regulation", model_id, SECONDARY_COMPARATOR_SEED)] = (
            run_id
        )
    runs.sort(key=lambda item: item["run_id"])
    plan: dict[str, Any] = {
        "schema_version": "masld-bench-plan-v1",
        "campaign": {"campaign_id": "selection-fixture", "wave": "selection"},
        "registry_snapshot_sha256": digest("registry-snapshot"),
        "source_lock": {"source_sha256": digest("source")},
        "resource_firewall": {"snapshot_sha256": digest("resource")},
        "dataset_locks": {
            "development_fixture": {
                "registry_sha256": development_registry_sha256,
                "activation_sha256": development_activation_sha256,
                "role": "external_development",
                "status": "available",
            },
            SEALED_DATASET: {
                "registry_sha256": sealed_registry_sha256,
                "role": "withheld_sealed",
                "status": "withheld_sealed",
            },
            PREDICTION_FIRST_STRESS_DATASET: {
                "registry_sha256": digest("stress-registry"),
                "activation_sha256": digest("stress-activation"),
                "role": "intervention_stress_test",
                "status": "available",
                "prediction_first_policy": {
                    "scope": "nonchampion_stress_test",
                },
            },
        },
        "model_dispositions": [
            {
                "model_id": SELECTED_MODEL,
                "family_id": "foundation_family_fixture",
                "disposition": "eligible_open",
                "champion_eligible": True,
            },
            {
                "model_id": BASELINE_MODEL,
                "family_id": "baseline_family_fixture",
                "disposition": "eligible_open",
                "champion_eligible": True,
            },
            *(
                {
                    "model_id": model_id,
                    "family_id": f"secondary_{model_id}_fixture",
                    "disposition": "eligible_restricted",
                    "champion_eligible": False,
                }
                for model_id in SECONDARY_COMPARATOR_MODELS
            ),
        ],
        "task_dispositions": [
            {
                "task_id": task_id,
                "status": "available",
                "disposition": "eligible",
                "blockers": [],
                "datasets_train": [],
                "datasets_development": [
                    "development_fixture",
                    PREDICTION_FIRST_STRESS_DATASET,
                ],
                "datasets_sealed": [SEALED_DATASET],
                "baseline_model_ids": sorted(
                    (BASELINE_MODEL, *SECONDARY_COMPARATOR_MODELS)
                    if task_id == "variant_to_regulation"
                    else (BASELINE_MODEL,)
                ),
            }
            for task_id in TASK_IDS
        ],
        "runs": runs,
    }
    plan["plan_sha256"] = canonical_hash(plan)
    write_json_exclusive(candidate / "plan.json", plan)
    freeze_tree(candidate, {"artifact_class": "candidate_plan"})
    return candidate, run_ids


def decision_document(
    *,
    metrics_sha256: str | None = None,
    evaluator_sha256: str | None = None,
    selected_regime: str = SELECTED_REGIME,
    baseline_regime: str = BASELINE_REGIME,
    secondary_run_ids: dict[str, str] | None = None,
    selection_artifact_hashes_by_task: Mapping[str, Mapping[str, str]] | None = None,
    open_champion: bool = False,
) -> dict[str, Any]:
    decisions: list[dict[str, Any]] = []
    for task_id in TASK_IDS:
        artifact_hashes = (
            selection_artifact_hashes_by_task[task_id]
            if selection_artifact_hashes_by_task is not None
            else {
                "checkpoint_sha256": digest(f"checkpoint-{task_id}"),
                "preprocessing_sha256": digest(f"preprocessing-{task_id}"),
                "task_head_sha256": digest(f"head-{task_id}"),
                "calibration_sha256": digest(f"calibration-{task_id}"),
            }
        )
        decisions.append(
            {
                "task_id": task_id,
                "selected_model_id": SELECTED_MODEL,
                "baseline_model_id": BASELINE_MODEL,
                "selected_adaptation_regime": selected_regime,
                "baseline_adaptation_regime": baseline_regime,
                "checkpoint_sha256": artifact_hashes["checkpoint_sha256"],
                "preprocessing_sha256": artifact_hashes[
                    "preprocessing_sha256"
                ],
                "task_head_sha256": artifact_hashes["task_head_sha256"],
                "seeds": list(SEEDS),
                "calibration_sha256": artifact_hashes["calibration_sha256"],
                "thresholds": {"promotion": _promotion_thresholds(task_id)},
                "evaluator_sha256": (
                    evaluator_sha256
                    if evaluator_sha256 is not None
                    else digest(f"evaluator-{task_id}")
                ),
                "promotion_gate": CHAMPION_GATE_THRESHOLDS[task_id][
                    "promotion_gate_id"
                ],
                "open_champion": open_champion,
            }
        )
        if task_id == "variant_to_regulation":
            decisions[-1]["variant_secondary_comparators"] = [
                {
                    "model_id": model_id,
                    "run_id": (
                        secondary_run_ids[model_id]
                        if secondary_run_ids is not None
                        else _simple_secondary_run_id(model_id)
                    ),
                }
                for model_id in SECONDARY_COMPARATOR_MODELS
            ]
    return {
        "metrics_sha256": metrics_sha256 or digest("development-metrics"),
        "task_decisions": decisions,
        "conditional_model_decision": {"status": "not_evaluated"},
        "power_decisions": {"status": "pending_external_seal"},
        "sealed_results_used": False,
    }


def _write_fixture_file(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_text_exclusive(path, text.strip() + "\n")
    return path


def _artifact_toml(
    *, path: str, sha256: str, size_bytes: int, media_type: str, role: str
) -> str:
    return (
        "{ "
        f'path = "{path}", sha256 = "{sha256}", size_bytes = {size_bytes}, '
        f'media_type = "{media_type}", role = "{role}"'
        " }"
    )


def _scientific_fixture_config(root: Path) -> Path:
    config = root / "scientific_config"
    artifacts = config / "artifacts"
    environment = _write_fixture_file(
        artifacts / "fixture-environment.lock",
        (
            f"python={Path(sys.executable).resolve().as_posix()}\n"
            f"version={sys.version_info.major}.{sys.version_info.minor}."
            f"{sys.version_info.micro}"
        ),
    )

    fixture_dataset_ids = (
        "development_fixture",
        "gse244832",
        "gse281160",
        "gse281364",
        "gse296875",
    )
    fixture_dataset_topology = {
        "development_fixture": {
            "modalities": ["rna", "regulatory"],
            "pairing_levels": ["same_sample_different_aliquot"],
        },
        "gse244832": {
            "modalities": ["single_nucleus_rna", "single_nucleus_atac"],
            "pairing_levels": ["same_sample_different_aliquot"],
        },
        "gse281160": {
            "modalities": ["single_cell_rna", "crispri_guide_assignment"],
            "pairing_levels": ["same_cell"],
        },
        "gse281364": {
            "modalities": ["mpra_dna_counts", "mpra_rna_counts"],
            "pairing_levels": ["same_sample_different_aliquot"],
        },
        "gse296875": {
            "modalities": ["single_nucleus_rna", "single_nucleus_atac"],
            "pairing_levels": ["same_nucleus"],
        },
    }
    for dataset_id in fixture_dataset_ids:
        topology = fixture_dataset_topology[dataset_id]
        dataset_artifact = _write_fixture_file(
            artifacts / f"{dataset_id}.manifest.json",
            json.dumps(
                {"dataset_id": dataset_id, "units": 10},
                sort_keys=True,
                separators=(",", ":"),
            ),
        )
        dataset_authorities: list[str] = []
        for authority in (
            "rights",
            "topology",
            "donor_join",
            "labels",
            "qc",
            "reference",
            "exposure_audit",
        ):
            relative = (
                f"artifacts/dataset-authorities/{dataset_id}-{authority}.json"
            )
            source = _write_fixture_file(
                config / relative,
                json.dumps(
                    {"authority": authority, "dataset_id": dataset_id},
                    sort_keys=True,
                    separators=(",", ":"),
                ),
            )
            dataset_authorities.append(
                _artifact_toml(
                    path=relative,
                    sha256=sha256_file(source),
                    size_bytes=source.stat().st_size,
                    media_type="application/json",
                    role=f"dataset_authority:{dataset_id}:{authority}",
                )
            )
        _write_fixture_file(
            config / "datasets" / f"{dataset_id}.toml",
            f'''
schema_version = "masld-bench-dataset-v1"
dataset_id = "{dataset_id}"
title = "Synthetic development endpoint fixture"
role = "external_development"
status = "available"
species = "human"
accession = ["FIXTURE-{dataset_id.upper()}"]
access_tier = "public"
automatic_download = false
redistribution_class = "fixture_only"
license_status = "project_owned"
biological_unit = "donor_or_ld_block"
expected_biological_units = 10
native_genome_build = "GRCh38"
native_annotation_release = "GENCODE_v49_fixture"
modalities = {json.dumps(topology["modalities"])}
pairing_levels = {json.dumps(topology["pairing_levels"])}
modality_status_default = "observed"
exposure_status = "target_label_unexposed"
admission_blocking = false
blockers = []
notes = "Synthetic task-native endpoint rows only."
[provenance]
source_url = "https://example.invalid/development-fixture"
verified_date = "2026-08-21"
local_exposure = "synthetic"
[split]
group_key = "fixture_group"
outer_role = "development"
label_visibility = "development_visible"
[activation]
schema_version = "masld-bench-dataset-activation-v1"
dataset_version = "fixture-v1"
rights_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-rights.json')}"
topology_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-topology.json')}"
donor_join_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-donor_join.json')}"
labels_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-labels.json')}"
qc_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-qc.json')}"
reference_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-reference.json')}"
exposure_audit_sha256 = "{sha256_file(config / f'artifacts/dataset-authorities/{dataset_id}-exposure_audit.json')}"
evidence_artifacts = [{', '.join(dataset_authorities)}]
artifact_manifest = {_artifact_toml(path=f'artifacts/{dataset_id}.manifest.json', sha256=sha256_file(dataset_artifact), size_bytes=dataset_artifact.stat().st_size, media_type='application/json', role=f'dataset:{dataset_id}')}
ready = true
''',
        )
    _write_fixture_file(
        config / "datasets" / f"{SEALED_DATASET}.toml",
        f'''
schema_version = "masld-bench-dataset-v1"
dataset_id = "{SEALED_DATASET}"
title = "Synthetic project-sealed final fixture"
role = "withheld_sealed"
status = "withheld_sealed"
species = "human"
accession = ["GSE289173"]
access_tier = "public"
automatic_download = false
outcome_url = "WITHHELD_SEALED"
redistribution_class = "source_reference_only"
license_status = "source_terms_apply"
biological_unit = "donor_or_ld_block"
expected_biological_units = 10
native_genome_build = "GRCh38"
native_annotation_release = "GENCODE_v49_fixture"
modalities = ["rna", "regulatory"]
pairing_levels = ["same_sample_different_aliquot"]
modality_status_default = "withheld_sealed"
exposure_status = "target_label_unexposed"
admission_blocking = true
blockers = ["outcomes_withheld_until_prediction_power_commit"]
notes = "Synthetic sealed-outcome identity only."
[provenance]
source_url = "https://example.invalid/gse289173-fixture"
verified_date = "2026-08-21"
local_exposure = "withheld"
[split]
group_key = "fixture_group"
outer_role = "sealed_final"
label_visibility = "withheld_sealed"
''',
    )

    model_authorities: dict[str, list[str]] = {}
    model_authority_hashes: dict[tuple[str, str], str] = {}
    for model_id in (
        SELECTED_MODEL,
        BASELINE_MODEL,
        *SECONDARY_COMPARATOR_MODELS,
    ):
        values: list[str] = []
        for authority in (
            "preprocessing",
            "code_license",
            "weights_license",
            "derivative_weights_license",
            "declared_corpora",
            "cellxgene_uuids",
            "exposure_audit",
        ):
            relative = f"artifacts/model-authorities/{model_id}-{authority}.txt"
            license_values = {
                "code_license": "project_owned",
                "weights_license": "not_applicable_classical",
                "derivative_weights_license": "project_owned",
            }
            if authority in license_values:
                content = json.dumps(
                    {
                        "schema_version": "masld-bench-license-authority-v1",
                        "model_id": model_id,
                        "authority": authority,
                        "declared_license": license_values[authority],
                        "use_allowed": True,
                        "redistribution_allowed": True,
                        "derivative_redistribution_allowed": True,
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                )
                media_type = "application/json"
            else:
                content = f"{model_id} {authority} fixture"
                media_type = "text/plain"
            source = _write_fixture_file(
                config / relative, content
            )
            digest_value = sha256_file(source)
            model_authority_hashes[(model_id, authority)] = digest_value
            values.append(
                _artifact_toml(
                    path=relative,
                    sha256=digest_value,
                    size_bytes=source.stat().st_size,
                    media_type=media_type,
                    role=f"model_authority:{model_id}:{authority}",
                )
            )
        model_authorities[model_id] = values

    environment_ref = _artifact_toml(
        path="artifacts/fixture-environment.lock",
        sha256=sha256_file(environment),
        size_bytes=environment.stat().st_size,
        media_type="text/plain",
        role="environment:fixture_cpu",
    )
    for model_id, family_id in (
        (SELECTED_MODEL, "foundation_family_fixture"),
        (BASELINE_MODEL, "baseline_family_fixture"),
    ):
        _write_fixture_file(
            config / "models" / f"{model_id}.toml",
            f'''
schema_version = "masld-bench-model-family-v1"
family_id = "{family_id}"
modality = ["transcriptome"]
status = "candidate"
[[models]]
model_id = "{model_id}"
display_name = "{model_id}"
implementation_type = "classical_baseline"
checkpoint_revision = "fixture-v1"
checkpoint_sha256 = "{digest(f'checkpoint-{model_id}')}"
license_status = "open"
exposure_status = "clean_declared"
status = "candidate"
admission_blocking = false
blockers = []
supported_tasks = ["{TASK_IDS[0]}", "{TASK_IDS[1]}"]
[models.execution]
schema_version = "masld-bench-model-execution-v1"
upstream_repository = "fixture://{model_id}"
upstream_commit = "fixture-v1"
preprocessing_sha256 = "{model_authority_hashes[(model_id, 'preprocessing')]}"
vocabulary_or_build = "not_applicable_classical"
genome_build = "GRCh38"
code_license = "project_owned"
weights_license = "not_applicable_classical"
derivative_weights_license = "project_owned"
training_cutoff = "not_applicable_classical"
declared_corpora_sha256 = "{model_authority_hashes[(model_id, 'declared_corpora')]}"
cellxgene_uuids_sha256 = "{model_authority_hashes[(model_id, 'cellxgene_uuids')]}"
exposure_audit_sha256 = "{model_authority_hashes[(model_id, 'exposure_audit')]}"
evidence_artifacts = [{', '.join(model_authorities[model_id])}]
supported_actions = ["probe", "prepare", "fit", "predict", "export"]
supported_adaptation_regimes = ["common_lane", "full_finetune", "native_lane"]
runtime_id = "fixture_cpu"
adapter_command = ["-m", "masld_bench.adapters.stub"]
environment_artifact = {environment_ref}
capture_r_session = false
ready = true
''',
        )

    _write_fixture_file(
        config / "models" / "context_borzoi.toml",
        f'''
schema_version = "masld-bench-model-family-v1"
family_id = "conditional_context_fixture"
modality = ["sequence", "transcriptome"]
status = "deferred"
[[models]]
model_id = "context_borzoi"
display_name = "Synthetic deferred conditional context model"
implementation_type = "pretrained_adapter"
checkpoint_revision = "fixture-v1"
checkpoint_sha256 = "{digest('checkpoint-context-borzoi')}"
license_status = "open"
exposure_status = "clean_declared"
status = "deferred"
admission_blocking = true
blockers = ["synthetic_fixture_deferred"]
supported_tasks = ["variant_to_regulation"]
''',
        )

    secondary_model_blocks = []
    for model_id in SECONDARY_COMPARATOR_MODELS:
        secondary_model_blocks.append(
            "\n".join(
                (
                    "[[models]]",
                    f'model_id = "{model_id}"',
                    f'display_name = "Synthetic {model_id} retrieval baseline"',
                    'implementation_type = "classical_baseline"',
                    'checkpoint_revision = "train_from_registered_inputs"',
                    'checkpoint_sha256 = "NOT_APPLICABLE"',
                    'license_status = "open"',
                    'exposure_status = "target_label_unexposed"',
                    'status = "candidate"',
                    "admission_blocking = false",
                    "blockers = []",
                    'supported_tasks = ["variant_to_regulation"]',
                    "[models.execution]",
                    'schema_version = "masld-bench-model-execution-v1"',
                    f'upstream_repository = "fixture://{model_id}"',
                    'upstream_commit = "fixture-v1"',
                    (
                        'preprocessing_sha256 = '
                        f'"{model_authority_hashes[(model_id, "preprocessing")]}"'
                    ),
                    'vocabulary_or_build = "not_applicable_classical"',
                    'genome_build = "GRCh38"',
                    'code_license = "project_owned"',
                    'weights_license = "not_applicable_classical"',
                    'derivative_weights_license = "project_owned"',
                    'training_cutoff = "not_applicable_classical"',
                    (
                        'declared_corpora_sha256 = '
                        f'"{model_authority_hashes[(model_id, "declared_corpora")]}"'
                    ),
                    (
                        'cellxgene_uuids_sha256 = '
                        f'"{model_authority_hashes[(model_id, "cellxgene_uuids")]}"'
                    ),
                    (
                        'exposure_audit_sha256 = '
                        f'"{model_authority_hashes[(model_id, "exposure_audit")]}"'
                    ),
                    (
                        "evidence_artifacts = ["
                        f'{", ".join(model_authorities[model_id])}]'
                    ),
                    'supported_actions = ["probe", "prepare", "fit", "predict", "export"]',
                    'supported_adaptation_regimes = ["full_finetune"]',
                    'runtime_id = "fixture_cpu"',
                    'adapter_command = ["-m", "masld_bench.adapters.stub"]',
                    f"environment_artifact = {environment_ref}",
                    "capture_r_session = false",
                    "ready = true",
                )
            )
        )
    _write_fixture_file(
        config / "models" / "secondary_variant_baselines.toml",
        "\n".join(
            (
                'schema_version = "masld-bench-model-family-v1"',
                'family_id = "secondary_variant_baseline_fixture"',
                'modality = ["dna_sequence", "regulatory_graph"]',
                'status = "candidate"',
                *secondary_model_blocks,
            )
        ),
    )
    variant_records = [
        *(
            _secondary_variant_capability_record(model_id)
            for model_id in SECONDARY_COMPARATOR_MODELS
        ),
        _variant_capability_record(
            BASELINE_MODEL,
            role="baseline",
            is_mandatory_baseline=True,
        ),
        _variant_capability_record(
            "context_borzoi",
            role="conditional_only",
            is_mandatory_baseline=False,
        ),
        _variant_capability_record(
            SELECTED_MODEL,
            role="primary_candidate",
            is_mandatory_baseline=False,
        ),
    ]
    variant_records.sort(key=lambda record: record["model_id"])
    variant_registry = {
        "schema_version": "masld-bench-variant-capability-v1",
        "registry_id": "variant_to_regulation_model_capability",
        "task_id": "variant_to_regulation",
        "primary_endpoint_id": "signed_cell_type_eqtl_effect",
        "primary_evaluator_id": "variant_ld_block_fisher_z_spearman_gain_v1",
        "sealed_available_context": ["dna_sequence", "single_nucleus_rna"],
        "sealed_unavailable_context": [
            "cell_type_eqtl_summary",
            "interaction_eqtl_summary",
            "single_nucleus_atac",
        ],
        "endpoint_ids": [
            "accessibility_delta",
            "allelic_direction",
            "enhancer_gene_link",
            "eqtl_retrieval",
            "mpra_activity",
            "signed_cell_type_eqtl_effect",
        ],
        "role_ids": [
            "baseline",
            "conditional_only",
            "link_only",
            "negative_control",
            "observed_context_only",
            "primary_candidate",
            "representation_only",
            "secondary_only",
        ],
        "models": variant_records,
        "admission_blocking": True,
        "blockers": ["synthetic_fixture_capability_contract"],
    }
    _write_fixture_file(
        config / "evaluation" / "variant_to_regulation_capabilities.toml",
        "\n".join(
            (
                'schema_version = "masld-bench-variant-capability-v1"',
                'registry_id = "variant_to_regulation_model_capability"',
                'task_id = "variant_to_regulation"',
                'primary_endpoint_id = "signed_cell_type_eqtl_effect"',
                'primary_evaluator_id = "variant_ld_block_fisher_z_spearman_gain_v1"',
                'sealed_available_context = ["dna_sequence", "single_nucleus_rna"]',
                'sealed_unavailable_context = ["cell_type_eqtl_summary", "interaction_eqtl_summary", "single_nucleus_atac"]',
                'endpoint_ids = ["accessibility_delta", "allelic_direction", "enhancer_gene_link", "eqtl_retrieval", "mpra_activity", "signed_cell_type_eqtl_effect"]',
                'role_ids = ["baseline", "conditional_only", "link_only", "negative_control", "observed_context_only", "primary_candidate", "representation_only", "secondary_only"]',
                "admission_blocking = true",
                'blockers = ["synthetic_fixture_capability_contract"]',
                *(
                    "\n".join(
                        (
                            "[[models]]",
                            f'model_id = {json.dumps(record["model_id"])}',
                            f'role = {json.dumps(record["role"])}',
                            f'native_outputs = {json.dumps(record["native_outputs"])}',
                            f'allowed_endpoints = {json.dumps(record["allowed_endpoints"])}',
                            f'primary_eligible = {str(record["primary_eligible"]).lower()}',
                            f'requires_fitted_head = {str(record["requires_fitted_head"]).lower()}',
                            "requires_observed_target_context = false",
                            f'is_mandatory_baseline = {str(record["is_mandatory_baseline"]).lower()}',
                        )
                    )
                    for record in variant_registry["models"]
                ),
            )
        ),
    )
    _write_fixture_file(
        config
        / "evaluation"
        / "variant_development_proxy_selection.json",
        (
            PACKAGE_ROOT
            / "config"
            / "evaluation"
            / "variant_development_proxy_selection.json"
        ).read_text(encoding="utf-8"),
    )

    promotion_source = PACKAGE_ROOT / "config" / "evaluation" / "promotion_gates.toml"
    promotion_target = _write_fixture_file(
        config / "evaluation" / "promotion_gates.toml",
        promotion_source.read_text(encoding="utf-8"),
    )
    promotion_sha256 = sha256_file(promotion_target)
    for task_id, evaluator_id, roster_field, roster, unit in (
        (
            TASK_IDS[0],
            "cell_donor_balanced_macro_f1_v1",
            "class_roster",
            _CELL_CLASSES,
            "donor",
        ),
        (
            TASK_IDS[1],
            "variant_ld_block_fisher_z_spearman_gain_v1",
            "strata",
            _VARIANT_LINEAGES,
            "ld_block",
        ),
    ):
        roster_relative = f"artifacts/task-rosters/{task_id}.json"
        roster_path = _write_fixture_file(
            config / roster_relative,
            json.dumps(
                {
                    "schema_version": "masld-bench-evaluator-roster-v1",
                    "task_id": task_id,
                    "primary_evaluator_id": evaluator_id,
                    "roster_field": roster_field,
                    "roster": list(roster),
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
        )
        roster_ref = _artifact_toml(
            path=roster_relative,
            sha256=sha256_file(roster_path),
            size_bytes=roster_path.stat().st_size,
            media_type="application/json",
            role=f"task_evaluator_roster:{task_id}",
        )
        _write_fixture_file(
            config / "tasks" / f"{task_id}.toml",
            f'''
schema_version = "masld-bench-task-v1"
task_id = "{task_id}"
status = "candidate"
unit_of_inference = "{unit}"
input_modalities = ["rna"]
endpoint = "Synthetic registered endpoint fixture."
metrics = ["primary_metric"]
datasets_train = {json.dumps(["gse244832", "gse281364", "gse296875"] if task_id == "variant_to_regulation" else ["development_fixture"])}
datasets_development = {json.dumps(["gse281160"] if task_id == "variant_to_regulation" else [])}
datasets_sealed = ["{SEALED_DATASET}"]
split_id = "fixture_grouped"
required_pairing_levels = []
missingness_policy = "explicit_state_and_mask"
admission_gates = []
baseline_model_ids = {json.dumps(sorted((BASELINE_MODEL, *SECONDARY_COMPARATOR_MODELS)) if task_id == "variant_to_regulation" else [BASELINE_MODEL])}
uncertainty_method = "paired_{unit}_bootstrap"
resampling_units = ["{unit}"]
bootstrap_replicates = 10000
multiplicity_family = "fixture_confirmatory"
multiplicity_method = "Holm"
promotion_gate_id = "{CHAMPION_GATE_THRESHOLDS[task_id]['promotion_gate_id']}"
promotion_gate_config_sha256 = "{promotion_sha256}"
primary_evaluator_id = "{evaluator_id}"
evaluator_parameters = {{ {roster_field} = {json.dumps(list(roster))}, roster_authority = {roster_ref} }}
claim_gate = "Synthetic fixture only."
''',
        )

    _write_fixture_file(
        config / "splits" / "fixture_grouped.toml",
        '''
schema_version = "masld-bench-split-v1"
split_id = "fixture_grouped"
entity = "biological_unit"
group_keys = ["fixture_group"]
roles = ["train", "development", "withheld_sealed"]
seed = 20260821
locus_grouping = "not_applicable"
label_policy = "Synthetic biological units remain in one outer role."
exposure_policy = "Synthetic checkpoint exposure is audited separately."
admission_gates = ["Synthetic fixture grouping is frozen."]
''',
    )
    _write_fixture_file(
        config / "runtimes" / "fixture.toml",
        f'''
schema_version = "masld-bench-runtime-v1"
runtime_id = "fixture_cpu"
scheduler = "slurm"
resource_profile = "cpu_contract"
status = "ready"
environment_lock = "fixture-python-{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
container_digest = "NOT_APPLICABLE"
modules = []
command_prefix = ["{Path(sys.executable).resolve().as_posix()}"]
admission_blocking = false
blockers = []
''',
    )
    _write_fixture_file(
        config / "resources.toml",
        (PACKAGE_ROOT / "config" / "resources.toml").read_text(encoding="utf-8"),
    )
    campaign = config / "campaigns" / "scientific_fixture.toml"
    _write_fixture_file(
        campaign,
        f'''
schema_version = "masld-bench-campaign-v1"
campaign_id = "selection-fixture-scientific"
wave = "full_specialist_screen"
status = "ready"
allow_arrays = false
allow_sealed_features = false
allow_sealed_labels = false
submit_enabled = false
task_ids = ["{TASK_IDS[0]}", "{TASK_IDS[1]}"]
prerequisites = []
[selection]
model_ids = {json.dumps(sorted((SELECTED_MODEL, BASELINE_MODEL, *SECONDARY_COMPARATOR_MODELS)))}
model_statuses = ["candidate"]
include_mandatory_baselines = true
[execution]
default_resource_profile = "cpu_contract"
runtime_id = "fixture_cpu"
seeds = [{', '.join(str(seed) for seed in SEEDS)}]
folds = [0]
adaptation_regime = "{SELECTED_REGIME}"
[execution.adapter_request]
row_ids = ["fixture-row-1", "fixture-row-2"]
''',
    )
    return campaign


def _development_cell_probabilities(
    *, unit_index: int, observed: str, class_index: int, model_id: str, seed: int
) -> dict[str, float]:
    """Per-seed class probabilities whose ensemble is not the per-seed mean.

    The selected model is confidently correct on four seeds and weakly wrong on
    the fifth, rotating which seed is wrong.  Every seed's own argmax is a valid
    hard label, each seed's macro-F1 is below one, and the probability-averaged
    ensemble recovers the correct class on every row.  That is precisely the
    nonlinearity that makes mean(metric(seed)) != metric(mean(prediction)).
    """

    wrong = _CELL_CLASSES[(class_index + 1) % len(_CELL_CLASSES)]
    third = next(
        class_id
        for class_id in _CELL_CLASSES
        if class_id not in {observed, wrong}
    )
    if model_id != SELECTED_MODEL:
        # The baseline is deterministically shifted and seed-invariant.
        return {
            class_id: (0.9 if class_id == wrong else 0.05)
            for class_id in _CELL_CLASSES
        }
    seed_index = SEEDS.index(seed)
    if unit_index % len(SEEDS) == seed_index:
        return {observed: 0.35, wrong: 0.40, third: 0.25}
    return {
        class_id: (0.9 if class_id == observed else 0.05)
        for class_id in _CELL_CLASSES
    }


def _development_cell_probability_index(
    *, model_id: str, seed: int
) -> dict[str, dict[str, float]]:
    """row_hash -> class probabilities, built before the rows are sorted."""

    index: dict[str, dict[str, float]] = {}
    for unit_index in range(10):
        for class_index, observed in enumerate(_CELL_CLASSES):
            row_hash = digest(
                f"development-cell-row-{unit_index}-{observed}"
            )
            index[row_hash] = _development_cell_probabilities(
                unit_index=unit_index,
                observed=observed,
                class_index=class_index,
                model_id=model_id,
                seed=seed,
            )
    return index


def _development_rows(
    task_id: str, model_id: str | None = None, seed: int | None = None
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    if task_id == "cell_state_mapping":
        if model_id is None:
            fields = ("row_hash", "unit_hash", "observed_class")
        else:
            fields = ("row_hash", "unit_hash", "predicted_class")
        for unit_index in range(10):
            unit_hash = digest(f"development-cell-donor-{unit_index}")
            for class_index, observed in enumerate(_CELL_CLASSES):
                row: dict[str, str] = {
                    "row_hash": digest(
                        f"development-cell-row-{unit_index}-{observed}"
                    ),
                    "unit_hash": unit_hash,
                }
                if model_id is None:
                    row["observed_class"] = observed
                else:
                    if seed is None:
                        raise ValueError(
                            "cell development predictions require a seed"
                        )
                    probabilities = _development_cell_probabilities(
                        unit_index=unit_index,
                        observed=observed,
                        class_index=class_index,
                        model_id=model_id,
                        seed=seed,
                    )
                    row["predicted_class"] = min(
                        _CELL_CLASSES,
                        key=lambda item: (-probabilities[item], item),
                    )
                rows.append(row)
    elif task_id == "variant_to_regulation":
        if model_id is None:
            fields = (
                "row_hash",
                "unit_hash",
                "block_hash",
                "stratum",
                "observed",
            )
        else:
            fields = (
                "row_hash",
                "unit_hash",
                "block_hash",
                "stratum",
                "predicted",
            )
        for block_index in range(10):
            block_hash = digest(f"development-variant-block-{block_index}")
            for stratum_index, stratum in enumerate(_VARIANT_LINEAGES):
                for position in range(3):
                    observed = float(
                        (block_index + 1) * 100
                        + (stratum_index + 1) * 10
                        + position
                    )
                    row = {
                        "row_hash": digest(
                            "development-variant-row-"
                            f"{block_index}-{stratum}-{position}"
                        ),
                        "unit_hash": digest(
                            "development-variant-unit-"
                            f"{block_index}-{stratum}-{position}"
                        ),
                        "block_hash": block_hash,
                        "stratum": stratum,
                    }
                    if model_id is None:
                        row["observed"] = str(observed)
                    else:
                        row["predicted"] = str(
                            observed
                            if model_id == SELECTED_MODEL
                            else (
                                1.0 / observed
                                if model_id in SECONDARY_COMPARATOR_MODELS
                                else -observed
                            )
                        )
                    rows.append(row)
    else:
        raise ValueError(f"unsupported development fixture task: {task_id}")
    rows.sort(key=lambda row: row["row_hash"])
    return fields, rows


def _freeze_development_prediction_bundle(
    root: Path, *, run: dict[str, Any]
) -> Path:
    task_id = str(run["task_id"])
    model_id = str(run["model_id"])
    bundle_root = root / "development_prediction_bundles" / str(run["run_id"])
    bundle_root.mkdir(parents=True)
    seed = int(run["seed"])
    fields, rows = _development_rows(task_id, model_id, seed)
    table_path = bundle_root / "predictions.tsv"
    row_path = bundle_root / "row_ids.tsv"
    write_text_exclusive(table_path, _tsv(fields, rows))
    write_text_exclusive(row_path, _tsv(("row_hash", "unit_hash"), rows))
    table = ArtifactRef.from_path(
        table_path,
        relative_to=bundle_root,
        media_type="text/tab-separated-values",
        role=f"standardized_prediction_table:{task_id}",
    )
    row_ids = ArtifactRef.from_path(
        row_path,
        relative_to=bundle_root,
        media_type="text/tab-separated-values",
        role=f"prediction_row_ids:{task_id}",
    )
    namespace = f"development_fixture:{task_id}:v1"
    biological_unit = "donor" if task_id == TASK_IDS[0] else "ld_block"
    join = canonical_hash(
        {
            "task_id": task_id,
            "dataset_ids": list(run["dataset_ids"]),
            "split_id": run["split_id"],
            "row_id_field": "row_hash",
            "unit_id_field": "unit_hash",
            "unit_id_namespace": namespace,
            "biological_unit": biological_unit,
        }
    )
    artifacts = [table, row_ids]
    if task_id == "cell_state_mapping":
        probability_fields = (
            "row_hash",
            "unit_hash",
            *(f"probability::{class_id}" for class_id in _CELL_CLASSES),
        )
        probability_index = _development_cell_probability_index(
            model_id=model_id, seed=seed
        )
        probability_rows = [
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                **{
                    f"probability::{class_id}": repr(
                        probability_index[row["row_hash"]][class_id]
                    )
                    for class_id in _CELL_CLASSES
                },
            }
            for row in rows
        ]
        probability_path = bundle_root / "class_probabilities.tsv"
        write_text_exclusive(
            probability_path, _tsv(probability_fields, probability_rows)
        )
        artifacts.append(
            ArtifactRef.from_path(
                probability_path,
                relative_to=bundle_root,
                media_type="text/tab-separated-values",
                role=f"class_probabilities:{task_id}",
            )
        )
    bundle = PredictionBundle(
        schema_version="masld-bench-prediction-bundle-v1",
        bundle_id=f"development-{str(run['run_id'])[:16]}",
        run_id=str(run["run_id"]),
        task_id=task_id,
        model_id=model_id,
        dataset_ids=tuple(run["dataset_ids"]),
        split_id=str(run["split_id"]),
        artifacts=tuple(artifacts),
        standardized_table=table,
        row_ids=row_ids,
        n_predictions=len(rows),
        row_id_field="row_hash",
        unit_id_field="unit_hash",
        unit_id_namespace=namespace,
        biological_unit=biological_unit,
        table_schema_sha256=canonical_hash(
            {"format": "tsv", "fields": list(fields)}
        ),
        source_join_key_sha256=join,
        format_version="tsv-v1",
        missing_state=MissingState.OBSERVED,
        metadata=(
            {
                "variant_endpoint_id": "eqtl_retrieval",
                "variant_evaluator_id": (
                    "variant_ld_block_eqtl_retrieval_auprc_v1"
                ),
                "variant_score_transform_id": "identity_link_score_v1",
                "primary_endpoint_scoring_allowed": False,
                "variant_capability_sha256": run["metadata"][
                    "variant_capability_sha256"
                ],
            }
            if model_id in SECONDARY_COMPARATOR_MODELS
            else {}
        ),
    )
    path = bundle_root / "prediction_bundle.json"
    write_json_exclusive(path, bundle.to_dict())
    return path


@patch(
    "masld_bench.planner.resource_snapshot",
    new=_fixture_resource_snapshot,
)
def _build_scientific_selection(root: Path) -> ScientificSelectionFixture:
    """Build the full succeeded-run to five-seed finalist authority chain."""

    campaign = _scientific_fixture_config(root)
    candidate = freeze_campaign(
        campaign,
        config_root=campaign.parents[1],
        package_root=PACKAGE_ROOT,
        output_root=root / "scientific_candidates",
    )
    # Freeze the complete five-seed campaign universe BEFORE any execution,
    # prediction, outcome, or receipt exists.  That ordering is what makes the
    # artifact provably pre-scoring.
    universe = freeze_finalist_campaign_universe(
        candidate_campaign_dirs=[candidate],
        output_root=root / "finalist_campaign_universes",
    )
    universe_payload = verify_finalist_campaign_universe(universe)
    plan = load_frozen_plan(candidate)
    attempts: dict[str, Path] = {}
    predictions: dict[str, Path] = {}
    for raw_run in plan["runs"]:
        run = dict(raw_run)
        run_id = str(run["run_id"])
        attempt = execute_run(
            candidate=candidate,
            run_id=run_id,
            output_root=root / "scientific_executions",
        )
        verify_run_execution_attempt(attempt, require_succeeded=True)
        attempts[run_id] = attempt
        predictions[run_id] = _freeze_development_prediction_bundle(root, run=run)

    run_index = {
        (str(run["task_id"]), str(run["model_id"]), int(run["seed"])): dict(run)
        for run in plan["runs"]
    }
    outcome_dirs: dict[str, Path] = {}
    for task_id in TASK_IDS:
        sample_run = run_index[(task_id, SELECTED_MODEL, SEEDS[0])]
        fields, rows = _development_rows(task_id)
        raw_table = root / "development_outcome_sources" / f"{task_id}.tsv"
        raw_table.parent.mkdir(parents=True, exist_ok=True)
        write_text_exclusive(raw_table, _tsv(fields, rows))
        namespace = f"development_fixture:{task_id}:v1"
        biological_unit = "donor" if task_id == TASK_IDS[0] else "ld_block"
        outcome_dirs[task_id] = freeze_development_outcome_bundle(
            task_id=task_id,
            dataset_ids=sample_run["dataset_ids"],
            dataset_registry_sha256s={
                dataset_id: plan["dataset_locks"][dataset_id]["registry_sha256"]
                for dataset_id in sample_run["dataset_ids"]
            },
            split_id=str(sample_run["split_id"]),
            row_id_field="row_hash",
            unit_id_field="unit_hash",
            unit_id_namespace=namespace,
            biological_unit=biological_unit,
            outcome_table_path=raw_table,
            output_root=root / "development_outcome_bundles",
        )

    receipts: list[Path] = []
    for task_id in TASK_IDS:
        for seed in SEEDS:
            candidate_run = run_index[(task_id, SELECTED_MODEL, seed)]
            baseline_run = run_index[(task_id, BASELINE_MODEL, seed)]
            receipt = freeze_scientific_run_receipt(
                candidate_campaign_dir=candidate,
                candidate_execution_attempt_dir=attempts[
                    str(candidate_run["run_id"])
                ],
                candidate_prediction_bundle_path=predictions[
                    str(candidate_run["run_id"])
                ],
                baseline_campaign_dir=candidate,
                baseline_execution_attempt_dir=attempts[
                    str(baseline_run["run_id"])
                ],
                baseline_prediction_bundle_path=predictions[
                    str(baseline_run["run_id"])
                ],
                development_outcome_bundle_dir=outcome_dirs[task_id],
                output_root=root / "scientific_run_receipts",
            )
            verify_scientific_run_receipt(receipt)
            receipts.append(receipt)
    secondary_receipts: list[Path] = []
    for model_id in SECONDARY_COMPARATOR_MODELS:
        task_id = "variant_to_regulation"
        candidate_run = run_index[
            (task_id, model_id, SECONDARY_COMPARATOR_SEED)
        ]
        receipt = freeze_variant_secondary_run_receipt(
            candidate_campaign_dir=candidate,
            execution_attempt_dir=attempts[str(candidate_run["run_id"])],
            prediction_bundle_path=predictions[str(candidate_run["run_id"])],
            output_root=root / "variant_secondary_run_receipts",
        )
        verify_variant_secondary_run_receipt(receipt)
        secondary_receipts.append(receipt)
    ledger = freeze_finalist_candidate_ledger(
        finalist_campaign_universe_dir=universe,
        reviewed_universe_id=universe_payload["universe_id"],
        scientific_run_receipt_dirs=receipts,
        variant_secondary_run_receipt_dirs=secondary_receipts,
        output_root=root / "selection_candidate_ledgers",
    )
    ledger_payload = verify_selection_candidate_ledger(ledger)
    selection_artifact_hashes_by_task = {
        task_id: dict(
            next(
                candidate
                for candidate in ledger_payload["model_candidates"]
                if candidate["task_id"] == task_id
                and candidate["model_id"] == SELECTED_MODEL
                and candidate["adaptation_regime"] == SELECTED_REGIME
            )["selection_artifact_ensemble_sha256s"]
        )
        for task_id in TASK_IDS
    }
    evaluator_sha256 = sha256_file(
        Path(recompute_endpoint_bootstrap.__code__.co_filename)
    )
    decisions_path = root / "scientific_selection_decisions.json"
    write_json_exclusive(
        decisions_path,
        decision_document(
            metrics_sha256=str(ledger_payload["ledger_id"]),
            evaluator_sha256=evaluator_sha256,
            selection_artifact_hashes_by_task=(
                selection_artifact_hashes_by_task
            ),
            selected_regime=SELECTED_REGIME,
            baseline_regime=SELECTED_REGIME,
            open_champion=True,
            secondary_run_ids={
                model_id: str(
                    run_index[
                        (
                            "variant_to_regulation",
                            model_id,
                            SECONDARY_COMPARATOR_SEED,
                        )
                    ]["run_id"]
                )
                for model_id in SECONDARY_COMPARATOR_MODELS
            },
        ),
    )
    lock_dir = freeze_selection_lock(
        candidate=ledger,
        decisions_path=decisions_path,
        output_root=root / "locks",
    )
    verify_selection_lock(lock_dir)
    finalist_dirs: dict[str, Path] = {}
    for task_id in TASK_IDS:
        finalist = freeze_finalist_development_metric_bundle(
            selection_lock_dir=lock_dir,
            selection_candidate_ledger_dir=ledger,
            task_id=task_id,
            output_root=root / "finalist_metric_bundles" / task_id,
        )
        verify_finalist_development_metric_bundle(finalist, reverify_sources=True)
        finalist_dirs[task_id] = finalist
    return ScientificSelectionFixture(
        candidate_dir=candidate,
        selection_lock_dir=lock_dir,
        selection_candidate_ledger_dir=ledger,
        finalist_campaign_universe_dir=universe,
        finalist_metric_bundle_dirs=finalist_dirs,
        execution_attempt_dirs=attempts,
        development_prediction_bundle_paths=predictions,
        scientific_receipt_dirs=tuple(receipts),
        variant_secondary_receipt_dirs=tuple(secondary_receipts),
        development_outcome_bundle_dirs=outcome_dirs,
    )


def create_scientific_selection(_root: Path) -> ScientificSelectionFixture:
    """Reuse one immutable synthetic authority chain across independent tests."""

    global _SCIENTIFIC_SELECTION_CACHE_DIR, _SCIENTIFIC_SELECTION_CACHE
    if _SCIENTIFIC_SELECTION_CACHE is None:
        _SCIENTIFIC_SELECTION_CACHE_DIR = tempfile.TemporaryDirectory(
            prefix="masld-bench-selection-fixture-"
        )
        cache_root = Path(_SCIENTIFIC_SELECTION_CACHE_DIR.name) / "fixture"
        cache_root.mkdir()
        _SCIENTIFIC_SELECTION_CACHE = _build_scientific_selection(cache_root)
    return _SCIENTIFIC_SELECTION_CACHE


def create_selection(root: Path) -> Path:
    fixture = create_scientific_selection(root)
    finalist_root = root / "finalist_metric_bundles"
    if not finalist_root.exists():
        for task_id, source in fixture.finalist_metric_bundle_dirs.items():
            shutil.copytree(source.parent, finalist_root / task_id)
    return fixture.selection_lock_dir


def _sealed_prediction_rows(task_id: str, model_id: str) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    if task_id == "cell_state_mapping":
        fields = ("row_hash", "unit_hash", "predicted_class")
        for index in range(30):
            observed = _CELL_CLASSES[index % len(_CELL_CLASSES)]
            predicted = (
                observed
                if model_id == SELECTED_MODEL
                else _CELL_CLASSES[(index + 1) % len(_CELL_CLASSES)]
            )
            rows.append(
                {
                    "row_hash": digest(f"sealed-{task_id}-row-{index}"),
                    "unit_hash": digest(f"sealed-{task_id}-donor-{index}"),
                    "predicted_class": predicted,
                }
            )
    elif task_id == "variant_to_regulation":
        fields = ("row_hash", "unit_hash", "block_hash", "stratum", "predicted")
        for block_index in range(10):
            block_hash = digest(f"sealed-{task_id}-block-{block_index}")
            for stratum_index, stratum in enumerate(_VARIANT_LINEAGES):
                for position in range(3):
                    observed = float(
                        (block_index + 1) * 100
                        + (stratum_index + 1) * 10
                        + position
                    )
                    identity = f"{block_index}-{stratum}-{position}"
                    rows.append(
                        {
                            "row_hash": digest(f"sealed-{task_id}-row-{identity}"),
                            "unit_hash": digest(
                                f"sealed-{task_id}-variant-{identity}"
                            ),
                            "block_hash": block_hash,
                            "stratum": stratum,
                            "predicted": str(
                                observed
                                if model_id == SELECTED_MODEL
                                else (
                                    1.0 / observed
                                    if model_id in SECONDARY_COMPARATOR_MODELS
                                    else -observed
                                )
                            ),
                        }
                    )
    else:
        raise ValueError(f"unsupported sealed fixture task: {task_id}")
    rows.sort(key=lambda row: row["row_hash"])
    return fields, rows


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    return "\n".join(
        ["\t".join(fields)]
        + ["\t".join(row[field] for field in fields) for row in rows]
    ) + "\n"


def commit_prediction_set(
    root: Path,
    *,
    selection_lock_dir: Path,
    evaluator_state_dir: Path,
) -> dict[str, tuple[PredictionCommit, ...]]:
    lock = verify_selection_lock(selection_lock_dir)
    by_task: dict[str, tuple[PredictionCommit, ...]] = {}
    for decision in lock.task_decisions:
        task_id = str(decision["task_id"])
        task_commits: list[PredictionCommit] = []
        runs_by_id: dict[
            str, tuple[dict[str, Any], str, dict[str, Any] | None]
        ] = {}
        for runs_field, model_field in (
            ("selected_runs", "selected_model_id"),
            ("baseline_runs", "baseline_model_id"),
        ):
            for raw_run in decision[runs_field]:
                run = dict(raw_run)
                runs_by_id[str(run["run_id"])] = (
                    run,
                    str(decision[model_field]),
                    None,
                )
        secondary = decision.get("variant_secondary_evaluation")
        # SelectionLock normalizes nested objects to mappingproxy, not dict, so
        # `isinstance(secondary, dict)` silently skipped every secondary
        # comparator and produced an incomplete prediction set.
        if isinstance(secondary, Mapping):
            for raw_comparator in secondary["comparators"]:
                comparator = dict(raw_comparator)
                run_id = str(comparator["run_id"])
                runs_by_id[run_id] = (
                    {"run_id": run_id},
                    str(comparator["model_id"]),
                    comparator,
                )
        for run_id, (run, model_id, secondary_record) in sorted(
            runs_by_id.items()
        ):
            run_id = str(run["run_id"])
            bundle_root = root / "prediction_bundles" / task_id / run_id
            bundle_root.mkdir(parents=True)
            fields, rows = _sealed_prediction_rows(task_id, model_id)
            prediction_path = bundle_root / "predictions.tsv"
            row_ids_path = bundle_root / "row_ids.tsv"
            write_text_exclusive(
                prediction_path,
                _tsv(fields, rows),
            )
            row_fields = ("row_hash", "unit_hash")
            write_text_exclusive(row_ids_path, _tsv(row_fields, rows))
            standardized_table = ArtifactRef.from_path(
                prediction_path,
                relative_to=bundle_root,
                media_type="text/tab-separated-values",
                role=f"standardized_prediction_table:{task_id}",
            )
            row_ids = ArtifactRef.from_path(
                row_ids_path,
                relative_to=bundle_root,
                media_type="text/tab-separated-values",
                role=f"prediction_row_ids:{task_id}",
            )
            artifacts = [standardized_table, row_ids]
            if task_id == "cell_state_mapping":
                probability_fields = (
                    "row_hash",
                    "unit_hash",
                    *(f"probability::{class_id}" for class_id in _CELL_CLASSES),
                )
                probability_rows: list[dict[str, str]] = []
                for row in rows:
                    predicted_class = row["predicted_class"]
                    probability_rows.append(
                        {
                            "row_hash": row["row_hash"],
                            "unit_hash": row["unit_hash"],
                            **{
                                f"probability::{class_id}": (
                                    "0.9" if class_id == predicted_class else "0.05"
                                )
                                for class_id in _CELL_CLASSES
                            },
                        }
                    )
                probability_path = bundle_root / "class_probabilities.tsv"
                write_text_exclusive(
                    probability_path,
                    _tsv(probability_fields, probability_rows),
                )
                artifacts.append(
                    ArtifactRef.from_path(
                        probability_path,
                        relative_to=bundle_root,
                        media_type="text/tab-separated-values",
                        role="class_probabilities:cell_state_mapping",
                    )
                )
            unit_id_namespace = f"{SEALED_DATASET}:{task_id}:v1"
            biological_unit = "donor" if task_id == "cell_state_mapping" else "ld_block"
            source_join_key_sha256 = canonical_hash(
                {
                    "task_id": task_id,
                    "dataset_ids": [SEALED_DATASET],
                    "split_id": "joint-sealed-final",
                    "row_id_field": "row_hash",
                    "unit_id_field": "unit_hash",
                    "unit_id_namespace": unit_id_namespace,
                    "biological_unit": biological_unit,
                }
            )
            bundle = PredictionBundle(
                schema_version="masld-bench-prediction-bundle-v1",
                bundle_id=f"{task_id}-{run_id[:12]}",
                run_id=run_id,
                task_id=task_id,
                model_id=model_id,
                dataset_ids=(SEALED_DATASET,),
                split_id="joint-sealed-final",
                artifacts=tuple(artifacts),
                standardized_table=standardized_table,
                row_ids=row_ids,
                n_predictions=len(rows),
                row_id_field="row_hash",
                unit_id_field="unit_hash",
                unit_id_namespace=unit_id_namespace,
                biological_unit=biological_unit,
                table_schema_sha256=canonical_hash(
                    {"format": "tsv", "fields": list(fields)}
                ),
                source_join_key_sha256=source_join_key_sha256,
                format_version="tsv-v1",
                missing_state=MissingState.OBSERVED,
                metadata=(
                    {
                        "variant_endpoint_id": "eqtl_retrieval",
                        "variant_evaluator_id": (
                            "variant_ld_block_eqtl_retrieval_auprc_v1"
                        ),
                        "variant_score_transform_id": "identity_link_score_v1",
                        "primary_endpoint_scoring_allowed": False,
                        "variant_secondary_comparator_sha256": canonical_hash(
                            secondary_record
                        ),
                    }
                    if secondary_record is not None
                    else {}
                ),
            )
            bundle_path = bundle_root / "prediction_bundle.json"
            write_json_exclusive(bundle_path, bundle.to_dict())
            task_commits.append(
                commit_predictions(
                    selection_lock_dir=selection_lock_dir,
                    prediction_bundle_path=bundle_path,
                    evaluator_state_dir=evaluator_state_dir,
                )
            )
        by_task[task_id] = tuple(task_commits)
    return by_task


def create_development_endpoint_table(
    root: Path,
    *,
    task_id: str,
    n_units: int = 10,
) -> tuple[Path, str, dict[str, object]]:
    table_root = root / "development_endpoint_tables"
    table_root.mkdir(parents=True, exist_ok=True)
    table = table_root / f"{task_id}.tsv"
    lines: list[str] = []
    if task_id == "cell_state_mapping":
        evaluator_id = "cell_donor_balanced_macro_f1_v1"
        parameters: dict[str, object] = {"class_roster": list(_CELL_CLASSES)}
        lines.append(
            "\t".join(
                (
                    "row_hash",
                    "unit_hash",
                    "observed_class",
                    "candidate_class",
                    "baseline_class",
                )
            )
        )
        for unit_index in range(n_units):
            unit_hash = digest(f"development-cell-donor-{unit_index}")
            for class_index, observed in enumerate(_CELL_CLASSES):
                baseline = _CELL_CLASSES[(class_index + 1) % len(_CELL_CLASSES)]
                lines.append(
                    "\t".join(
                        (
                            digest(f"development-cell-row-{unit_index}-{observed}"),
                            unit_hash,
                            observed,
                            observed,
                            baseline,
                        )
                    )
                )
    elif task_id == "variant_to_regulation":
        evaluator_id = "variant_ld_block_fisher_z_spearman_gain_v1"
        parameters = {"strata": list(_VARIANT_LINEAGES)}
        lines.append(
            "\t".join(
                (
                    "row_hash",
                    "unit_hash",
                    "block_hash",
                    "stratum",
                    "observed",
                    "candidate",
                    "baseline",
                )
            )
        )
        for block_index in range(n_units):
            block_hash = digest(f"development-ld-block-{block_index}")
            for stratum_index, stratum in enumerate(_VARIANT_LINEAGES):
                for position in range(3):
                    observed = float(
                        (block_index + 1) * 100
                        + (stratum_index + 1) * 10
                        + position
                    )
                    lines.append(
                        "\t".join(
                            (
                                digest(
                                    "development-variant-row-"
                                    f"{block_index}-{stratum}-{position}"
                                ),
                                digest(
                                    "development-variant-unit-"
                                    f"{block_index}-{stratum}-{position}"
                                ),
                                block_hash,
                                stratum,
                                str(observed),
                                str(observed),
                                str(-observed),
                            )
                        )
                    )
    else:
        raise ValueError(f"unsupported sealed fixture task: {task_id}")
    write_text_exclusive(table, "\n".join(lines) + "\n")
    return table, evaluator_id, parameters


def freeze_power_set(
    *,
    evaluator_state_dir: Path,
    commits: dict[str, tuple[PredictionCommit, ...]],
    underpowered_tasks: set[str] | None = None,
) -> dict[str, PowerDecision]:
    if underpowered_tasks:
        raise ValueError(
            "production power fixtures cannot override recursively verified unit counts"
        )
    decisions: dict[str, PowerDecision] = {}
    for task_id, task_commits in commits.items():
        selection_lock_dir = Path(task_commits[0].selection_lock_dir)
        fixture_root = selection_lock_dir.parents[1]
        metric_candidates = sorted(
            path
            for path in (fixture_root / "finalist_metric_bundles" / task_id).iterdir()
            if path.is_dir()
        )
        if len(metric_candidates) != 1:
            raise ValueError(
                f"fixture lacks exactly one finalist metric bundle for {task_id}"
            )
        evidence_dir = freeze_development_power_evidence(
            finalist_development_metric_bundle_dir=metric_candidates[0],
            output_root=evaluator_state_dir,
        )
        decisions[task_id] = freeze_power_decision(
            evaluator_state_dir=evaluator_state_dir,
            prediction_commit_sha256s=(
                commit.prediction_commit_sha256 for commit in task_commits
            ),
            development_power_evidence_dir=evidence_dir,
        )
    return decisions


def freeze_terminal_authorization(
    root: Path,
    *,
    selection_lock_dir: Path,
    evaluator_state_dir: Path,
    commits: dict[str, tuple[PredictionCommit, ...]],
    powers: dict[str, PowerDecision],
    sealed_outcome_bundle_dir: Path,
    outcome_consumption_path: Path,
) -> Path:
    """Freeze task-complete terminal failures from the real consumed seal chain."""

    gates: list[Path] = []
    for task_id in TASK_IDS:
        gate = freeze_champion_gate_decision(
            selection_lock_dir=selection_lock_dir,
            task_id=task_id,
            prediction_commit_dirs=(
                evaluator_state_dir
                / "prediction_commits"
                / commit.prediction_commit_sha256
                for commit in commits[task_id]
            ),
            power_decision_dir=(
                evaluator_state_dir
                / "power_decisions"
                / powers[task_id].power_decision_sha256
            ),
            sealed_outcome_bundle_dir=sealed_outcome_bundle_dir,
            outcome_consumption_path=outcome_consumption_path,
            output_root=root / "champion_gate_decisions",
        )
        gates.append(gate)
    return freeze_terminal_evaluation_authorization(
        selection_lock_dir=selection_lock_dir,
        champion_gate_decision_dirs=gates,
        output_root=root / "terminal_authorizations",
    )


def create_joint_outcome_bundle(
    root: Path,
    *,
    unit_id_field: str = "unit_hash",
    unit_offsets: dict[str, int] | None = None,
) -> Path:
    outcome = root / "sealed_outcome"
    outcome.mkdir(parents=True)
    artifacts: list[ArtifactRef] = []
    for task_id in TASK_IDS:
        offset = (unit_offsets or {}).get(task_id, 0)
        rows: list[dict[str, str]] = []
        if task_id == "cell_state_mapping":
            fields = ("row_hash", unit_id_field, "observed_class")
            for index in range(30):
                rows.append(
                    {
                        "row_hash": digest(f"sealed-{task_id}-row-{index}"),
                        unit_id_field: digest(
                            f"sealed-{task_id}-donor-{(index + offset) % 30}"
                        ),
                        "observed_class": _CELL_CLASSES[index % len(_CELL_CLASSES)],
                    }
                )
        else:
            fields = (
                "row_hash",
                unit_id_field,
                "block_hash",
                "stratum",
                "observed",
            )
            for block_index in range(10):
                block_hash = digest(f"sealed-{task_id}-block-{block_index}")
                for stratum_index, stratum in enumerate(_VARIANT_LINEAGES):
                    for position in range(3):
                        observed = float(
                            (block_index + 1) * 100
                            + (stratum_index + 1) * 10
                            + position
                        )
                        identity = f"{block_index}-{stratum}-{position}"
                        shifted_block = (block_index + offset) % 10
                        rows.append(
                            {
                                "row_hash": digest(
                                    f"sealed-{task_id}-row-{identity}"
                                ),
                                unit_id_field: digest(
                                    "sealed-variant_to_regulation-variant-"
                                    f"{shifted_block}-{stratum}-{position}"
                                ),
                                "block_hash": block_hash,
                                "stratum": stratum,
                                "observed": str(observed),
                            }
                        )
        rows.sort(key=lambda row: row["row_hash"])
        outcome_path = outcome / f"{task_id}.tsv"
        write_text_exclusive(outcome_path, _tsv(fields, rows))
        artifacts.append(
            ArtifactRef.from_path(
                outcome_path,
                relative_to=outcome,
                media_type="text/tab-separated-values",
                role=f"outcome:{task_id}",
            )
        )
        if task_id == "variant_to_regulation":
            secondary_fields = (
                "row_hash",
                unit_id_field,
                "block_hash",
                "stratum",
                "observed_binary",
                "summary_statistics_complete",
            )
            secondary_rows = [
                {
                    "row_hash": row["row_hash"],
                    unit_id_field: row[unit_id_field],
                    "block_hash": row["block_hash"],
                    "stratum": row["stratum"],
                    "observed_binary": (
                        "1" if float(row["observed"]) >= 600.0 else "0"
                    ),
                    "summary_statistics_complete": "true",
                }
                for row in rows
            ]
            secondary_path = outcome / f"{task_id}.secondary.tsv"
            write_text_exclusive(
                secondary_path,
                _tsv(secondary_fields, secondary_rows),
            )
            artifacts.append(
                ArtifactRef.from_path(
                    secondary_path,
                    relative_to=outcome,
                    media_type="text/tab-separated-values",
                    role=f"secondary_outcome:{task_id}",
                )
            )
    identity = {
        "schema_version": SEALED_OUTCOME_BUNDLE_SCHEMA_VERSION,
        "bundle_id": "gse289173_joint_cell_and_regulatory_bundle",
        "joint_bundle": True,
        "task_bindings": [
            {"task_id": task_id, "dataset_ids": [SEALED_DATASET]}
            for task_id in TASK_IDS
        ],
        "unit_id_fields": {task_id: unit_id_field for task_id in TASK_IDS},
        "artifacts": [artifact.to_dict() for artifact in artifacts],
    }
    bundle_sha256 = canonical_hash(identity)
    write_json_exclusive(
        outcome / "sealed_outcome_bundle.json",
        {**identity, "bundle_sha256": bundle_sha256},
    )
    freeze_tree(
        outcome,
        {
            "artifact_class": "sealed_outcome_bundle",
            "bundle_sha256": bundle_sha256,
        },
    )
    return outcome
