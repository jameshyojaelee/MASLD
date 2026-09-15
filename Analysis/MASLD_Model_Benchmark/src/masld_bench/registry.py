"""Deterministic loading and cross-validation of benchmark registries."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tomllib
from types import MappingProxyType
from typing import Any, Mapping, TypeVar

from .contracts import (
    ArtifactRef,
    ContractError,
    DatasetManifest,
    DatasetViewContract,
    ExposureState,
    MissingState,
    ModelManifest,
    PairingState,
    PredictionBundle,
    RunSpec,
    SelectionLock,
    SplitSpec,
    StrictContract,
    TaskSpec,
    _as_bool,
    _as_int,
    _as_string,
    _as_string_tuple,
    _strict_mapping,
)
from .artifacts import ArtifactError, verify_frozen_tree
from .hashing import canonical_sha256, canonicalize, sha256_file
from .reference import ReferenceBundle, ReferenceError


class RegistryError(ContractError):
    """Raised when a registry document or cross-reference is invalid."""


_EXPECTED_TASK_EVALUATORS: Mapping[str, Mapping[str, Any]] = MappingProxyType(
    {
        "bulk_state_transfer": {
            "primary_evaluator_id": "bulk_paired_spearman_gain_v1",
            "evaluator_parameters": {},
        },
        "typed_evidence_graph": {
            "primary_evaluator_id": "graph_ld_block_auprc_gain_v1",
            "evaluator_parameters": {},
        },
        "perturbation_transfer": {
            "primary_evaluator_id": "not_applicable_nonpromotable_v1",
            "evaluator_parameters": {},
        },
        "gene_perturbation_response": {
            "primary_evaluator_id": "not_applicable_nonpromotable_v1",
            "evaluator_parameters": {},
        },
        "unified_model": {
            "primary_evaluator_id": "not_applicable_nonpromotable_v1",
            "evaluator_parameters": {},
        },
    }
)
_ROSTER_TASK_EVALUATORS: Mapping[str, tuple[str, str]] = MappingProxyType(
    {
        "cell_state_mapping": (
            "cell_donor_balanced_macro_f1_v1",
            "class_roster",
        ),
        "variant_to_regulation": (
            "variant_ld_block_fisher_z_spearman_gain_v1",
            "strata",
        ),
        "rna_conditioned_atac": (
            "rna_atac_two_way_deviance_reduction_v1",
            "strata",
        ),
    }
)
_LICENSE_PERMISSIONS: Mapping[str, Mapping[str, bool]] = MappingProxyType(
    {
        **{
            term: {
                "use_allowed": True,
                "redistribution_allowed": True,
                "derivative_redistribution_allowed": True,
            }
            for term in (
                "MIT",
                "BSD-2-Clause",
                "BSD-3-Clause",
                "Apache-2.0",
                "GPL-2.0-only",
                "GPL-2.0-or-later",
                "GPL-3.0-only",
                "GPL-3.0-or-later",
                "LGPL-3.0-only",
                "LGPL-3.0-or-later",
                "CC-BY-4.0",
                "CC-BY-SA-4.0",
                "CC0-1.0",
                "project_owned",
                "not_applicable_classical",
            )
        },
        **{
            term: {
                "use_allowed": True,
                "redistribution_allowed": False,
                "derivative_redistribution_allowed": False,
            }
            for term in (
                "research_only",
                "noncommercial",
                "CC-BY-NC-4.0",
                "CC-BY-NC-SA-4.0",
                "CC-BY-NC-ND-4.0",
            )
        },
        **{
            term: {
                "use_allowed": False,
                "redistribution_allowed": False,
                "derivative_redistribution_allowed": False,
            }
            for term in (
                "restricted",
                "redistribution_prohibited",
                "source_terms_apply",
            )
        },
    }
)


_ContractT = TypeVar("_ContractT", bound=StrictContract)


def _load_document(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    try:
        if source.suffix == ".toml":
            with source.open("rb") as handle:
                value = tomllib.load(handle)
        elif source.suffix == ".json":
            value = json.loads(source.read_text(encoding="utf-8"))
        else:
            raise RegistryError(f"unsupported registry extension: {source}")
    except (OSError, tomllib.TOMLDecodeError, json.JSONDecodeError) as error:
        raise RegistryError(f"cannot parse registry document {source}: {error}") from error
    if not isinstance(value, dict):
        raise RegistryError(f"registry document must contain an object/table: {source}")
    return value


def load_contract(
    path: str | Path,
    contract_type: type[_ContractT],
    *,
    table: str | None = None,
) -> _ContractT:
    value = _load_document(path)
    if table is not None:
        root = _strict_mapping(
            value,
            allowed={table},
            required={table},
            label=f"{Path(path).name} envelope",
        )
        value = root[table]
    try:
        return contract_type.from_dict(value)
    except ContractError as error:
        raise RegistryError(f"invalid {contract_type.__name__} in {path}: {error}") from error


def load_dataset_manifest(path: str | Path) -> DatasetManifest:
    return load_contract(path, DatasetManifest)


def load_dataset_view_contract(path: str | Path) -> DatasetViewContract:
    return load_contract(path, DatasetViewContract)


def load_model_manifest(path: str | Path) -> ModelManifest:
    return load_contract(path, ModelManifest)


def load_task_spec(path: str | Path) -> TaskSpec:
    return load_contract(path, TaskSpec)


def load_split_spec(path: str | Path) -> SplitSpec:
    return load_contract(path, SplitSpec)


def load_run_spec(path: str | Path) -> RunSpec:
    return load_contract(path, RunSpec)


def _validate_ready_evidence_artifact(
    artifact: ArtifactRef,
    *,
    config_root: Path,
    label: str,
) -> Path:
    configured = Path(artifact.path)
    try:
        return artifact.validate(None if configured.is_absolute() else config_root)
    except ContractError as error:
        raise RegistryError(f"{label} failed artifact verification: {error}") from error


def _validate_model_license_authority(
    *,
    path: Path,
    model_id: str,
    authority: str,
    declared_license: str,
) -> None:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RegistryError(
            f"model {model_id} {authority} authority is not valid JSON: {error}"
        ) from error
    permissions = _LICENSE_PERMISSIONS.get(declared_license)
    if permissions is None:
        raise RegistryError(
            f"model {model_id} {authority} has no recognized license policy"
        )
    expected = {
        "schema_version": "masld-bench-license-authority-v1",
        "model_id": model_id,
        "authority": authority,
        "declared_license": declared_license,
        **permissions,
    }
    if canonicalize(document) != canonicalize(expected):
        raise RegistryError(
            f"model {model_id} {authority} authority content differs from its declared terms"
        )


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RegistryError(f"{label} is not valid JSON: {error}") from error
    if not isinstance(value, dict):
        raise RegistryError(f"{label} must contain a JSON object")
    return value


def _sorted_capability_ids(
    document: Mapping[str, Any],
    field: str,
    *,
    allow_empty: bool = False,
    label: str = "rna_atac_capabilities",
) -> tuple[str, ...]:
    values = _as_string_tuple(
        document[field],
        f"{label}.{field}",
        allow_empty=allow_empty,
    )
    if values != tuple(sorted(set(values))):
        raise RegistryError(
            f"{label}.{field} must be sorted and unique"
        )
    return values


def _validate_rna_atac_capability_registry(
    document: Mapping[str, Any],
    *,
    models: Mapping[str, ModelManifest],
    task: TaskSpec,
) -> dict[str, Any]:
    fields = frozenset(
        {
            "schema_version",
            "registry_id",
            "task_id",
            "profile_output_definition",
            "latent_integration_definition",
            "profile_metrics",
            "retrieval_metrics",
            "direct_profile_candidate_models",
            "conditional_profile_candidate_models",
            "profile_baseline_models",
            "sequence_required_models",
            "latent_only_models",
            "profile_fixture_passed_models",
            "latent_only_status",
            "admission_blocking",
            "blockers",
        }
    )
    root = _strict_mapping(
        document,
        allowed=fields,
        required=fields,
        label="rna_atac_capabilities",
    )
    if root["schema_version"] != "masld-bench-atac-capability-v2":
        raise RegistryError("unsupported RNA-ATAC capability schema")
    if root["registry_id"] != "rna_conditioned_atac_model_capability":
        raise RegistryError("RNA-ATAC capability registry_id is not canonical")
    if root["task_id"] != task.task_id or task.task_id != "rna_conditioned_atac":
        raise RegistryError("RNA-ATAC capability task binding differs from TaskSpec")
    for field in ("profile_output_definition", "latent_integration_definition"):
        _as_string(root[field], f"rna_atac_capabilities.{field}")

    profile_metrics = _sorted_capability_ids(root, "profile_metrics")
    retrieval_metrics = _sorted_capability_ids(root, "retrieval_metrics")
    if set(profile_metrics) & set(retrieval_metrics):
        raise RegistryError("RNA-ATAC profile and retrieval metrics must be disjoint")
    task_metrics = set(task.metrics)
    if not set(profile_metrics).issubset(task_metrics):
        raise RegistryError("RNA-ATAC TaskSpec omits a registered profile metric")
    if set(retrieval_metrics) & task_metrics:
        raise RegistryError(
            "RNA-ATAC retrieval metrics require a separate nonpromotable TaskSpec"
        )

    categories = {
        field: _sorted_capability_ids(root, field)
        for field in (
            "direct_profile_candidate_models",
            "conditional_profile_candidate_models",
            "profile_baseline_models",
            "latent_only_models",
        )
    }
    category_names = tuple(categories)
    for index, left_name in enumerate(category_names):
        for right_name in category_names[index + 1 :]:
            overlap = set(categories[left_name]) & set(categories[right_name])
            if overlap:
                raise RegistryError(
                    "RNA-ATAC capability categories must be disjoint: "
                    + ",".join(sorted(overlap))
                )
    registered_ids = set().union(*(set(values) for values in categories.values()))
    unknown = sorted(registered_ids.difference(models))
    if unknown:
        raise RegistryError(
            "RNA-ATAC capability registry names unknown models: " + ", ".join(unknown)
        )

    profile_models = set(categories["direct_profile_candidate_models"])
    profile_models.update(categories["conditional_profile_candidate_models"])
    profile_models.update(categories["profile_baseline_models"])
    sequence_required = set(
        _sorted_capability_ids(root, "sequence_required_models")
    )
    if not sequence_required.issubset(profile_models):
        raise RegistryError(
            "RNA-ATAC sequence-required roster includes a non-profile model"
        )
    supported = {
        model.model_id
        for model in models.values()
        if task.task_id in model.supported_tasks
    }
    if supported != profile_models:
        missing = sorted(profile_models.difference(supported))
        extra = sorted(supported.difference(profile_models))
        raise RegistryError(
            "RNA-ATAC TaskSpec/model capability mismatch: "
            f"missing={missing}; extra={extra}"
        )
    if task.baseline_model_ids != categories["profile_baseline_models"]:
        raise RegistryError(
            "RNA-ATAC TaskSpec baseline roster differs from capability registry"
        )

    for model_id in categories["latent_only_models"]:
        model_tasks = set(models[model_id].supported_tasks)
        if task.task_id in model_tasks:
            raise RegistryError(
                f"latent-only model {model_id} has an invalid task registration"
            )
    if (
        root["latent_only_status"]
        != "registered_not_schedulable_until_separate_nonpromotable_task_exists"
    ):
        raise RegistryError("RNA-ATAC latent-only scheduling status is not fail-closed")

    passed = set(
        _sorted_capability_ids(root, "profile_fixture_passed_models", allow_empty=True)
    )
    if not passed.issubset(profile_models):
        raise RegistryError("RNA-ATAC fixture registry includes a non-profile model")
    unpassed = profile_models.difference(passed)
    admission_blocking = _as_bool(
        root["admission_blocking"], "rna_atac_capabilities.admission_blocking"
    )
    if unpassed and not admission_blocking:
        raise RegistryError("unpassed RNA-ATAC profile fixtures must block admission")
    for model_id in sorted(unpassed):
        if not models[model_id].admission_blocking:
            raise RegistryError(
                f"model {model_id} lacks its required RNA-ATAC fixture blocker"
            )
    _as_string_tuple(
        root["blockers"], "rna_atac_capabilities.blockers", allow_empty=False
    )
    return canonicalize(root)


def _validate_variant_capability_registry(
    document: Mapping[str, Any],
    *,
    models: Mapping[str, ModelManifest],
    task: TaskSpec,
) -> dict[str, Any]:
    label = "variant_capabilities"
    fields = frozenset(
        {
            "schema_version",
            "registry_id",
            "task_id",
            "primary_endpoint_id",
            "primary_evaluator_id",
            "sealed_available_context",
            "sealed_unavailable_context",
            "endpoint_ids",
            "role_ids",
            "models",
            "admission_blocking",
            "blockers",
        }
    )
    root = _strict_mapping(
        document,
        allowed=fields,
        required=fields,
        label=label,
    )
    if root["schema_version"] != "masld-bench-variant-capability-v1":
        raise RegistryError("unsupported variant capability schema")
    if root["registry_id"] != "variant_to_regulation_model_capability":
        raise RegistryError("variant capability registry_id is not canonical")
    if root["task_id"] != task.task_id or task.task_id != "variant_to_regulation":
        raise RegistryError("variant capability task binding differs from TaskSpec")
    if root["primary_evaluator_id"] != task.primary_evaluator_id:
        raise RegistryError("variant capability primary evaluator differs from TaskSpec")

    expected_endpoints = (
        "accessibility_delta",
        "allelic_direction",
        "enhancer_gene_link",
        "eqtl_retrieval",
        "mpra_activity",
        "signed_cell_type_eqtl_effect",
    )
    endpoint_ids = _sorted_capability_ids(root, "endpoint_ids", label=label)
    if endpoint_ids != expected_endpoints:
        raise RegistryError("variant capability endpoint roster is not canonical")
    primary_endpoint = _as_string(
        root["primary_endpoint_id"], f"{label}.primary_endpoint_id"
    )
    if primary_endpoint != "signed_cell_type_eqtl_effect":
        raise RegistryError("variant capability primary endpoint is not canonical")

    expected_roles = (
        "baseline",
        "conditional_only",
        "link_only",
        "negative_control",
        "observed_context_only",
        "primary_candidate",
        "representation_only",
        "secondary_only",
    )
    role_ids = _sorted_capability_ids(root, "role_ids", label=label)
    if role_ids != expected_roles:
        raise RegistryError("variant capability role roster is not canonical")
    available_context = _sorted_capability_ids(
        root, "sealed_available_context", label=label
    )
    unavailable_context = _sorted_capability_ids(
        root, "sealed_unavailable_context", label=label
    )
    if available_context != ("dna_sequence", "single_nucleus_rna"):
        raise RegistryError("variant sealed input context is not canonical")
    if unavailable_context != (
        "cell_type_eqtl_summary",
        "interaction_eqtl_summary",
        "single_nucleus_atac",
    ):
        raise RegistryError("variant sealed unavailable context is not canonical")
    if set(available_context) & set(unavailable_context):
        raise RegistryError("variant sealed available and unavailable context overlap")

    raw_records = root["models"]
    if not isinstance(raw_records, list) or not raw_records:
        raise RegistryError("variant capability models must be a non-empty array")
    record_fields = frozenset(
        {
            "model_id",
            "role",
            "native_outputs",
            "allowed_endpoints",
            "primary_eligible",
            "requires_fitted_head",
            "requires_observed_target_context",
            "is_mandatory_baseline",
        }
    )
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    primary_candidates: set[str] = set()
    primary_baselines: set[str] = set()
    mandatory_baselines: set[str] = set()
    conditional_models: set[str] = set()
    for index, raw_record in enumerate(raw_records):
        record = _strict_mapping(
            raw_record,
            allowed=record_fields,
            required=record_fields,
            label=f"{label}.models[{index}]",
        )
        model_id = _as_string(record["model_id"], f"{label}.models[{index}].model_id")
        if model_id in seen:
            raise RegistryError(f"duplicate variant capability model: {model_id}")
        seen.add(model_id)
        if model_id not in models:
            raise RegistryError(f"variant capability registry names unknown model: {model_id}")
        role = _as_string(record["role"], f"{label}.{model_id}.role")
        if role not in role_ids:
            raise RegistryError(f"variant model {model_id} has an unknown capability role")
        native_outputs = _sorted_capability_ids(
            record, "native_outputs", label=f"{label}.{model_id}"
        )
        allowed_endpoints = _sorted_capability_ids(
            record, "allowed_endpoints", label=f"{label}.{model_id}"
        )
        unknown_endpoints = sorted(set(allowed_endpoints).difference(endpoint_ids))
        if unknown_endpoints:
            raise RegistryError(
                f"variant model {model_id} names unknown endpoints: "
                + ", ".join(unknown_endpoints)
            )
        primary_eligible = _as_bool(
            record["primary_eligible"], f"{label}.{model_id}.primary_eligible"
        )
        requires_fitted_head = _as_bool(
            record["requires_fitted_head"],
            f"{label}.{model_id}.requires_fitted_head",
        )
        requires_observed_context = _as_bool(
            record["requires_observed_target_context"],
            f"{label}.{model_id}.requires_observed_target_context",
        )
        is_mandatory_baseline = _as_bool(
            record["is_mandatory_baseline"],
            f"{label}.{model_id}.is_mandatory_baseline",
        )
        if primary_eligible != (primary_endpoint in allowed_endpoints):
            raise RegistryError(
                f"variant model {model_id} primary eligibility and endpoint differ"
            )
        if primary_eligible:
            if role not in {"baseline", "conditional_only", "primary_candidate"}:
                raise RegistryError(
                    f"variant model {model_id} has an invalid primary capability role"
                )
            if requires_observed_context:
                raise RegistryError(
                    f"variant model {model_id} cannot use unavailable observed context on the primary endpoint"
                )
            if not set(native_outputs) & {"gene_expression_delta", "rna_coverage_delta"}:
                raise RegistryError(
                    f"variant model {model_id} lacks a signed gene-level native output"
                )
            if role == "primary_candidate":
                primary_candidates.add(model_id)
            if role == "baseline":
                primary_baselines.add(model_id)
        if (role == "observed_context_only") != requires_observed_context:
            raise RegistryError(
                f"variant model {model_id} observed-context role and requirement differ"
            )
        if role == "representation_only" and (
            not requires_fitted_head or primary_eligible
        ):
            raise RegistryError(
                f"variant representation-only model {model_id} must require a head and remain secondary"
            )
        if role == "link_only" and (
            primary_eligible
            or not set(allowed_endpoints).issubset(
                {"enhancer_gene_link", "eqtl_retrieval"}
            )
        ):
            raise RegistryError(f"variant link-only model {model_id} has invalid endpoints")
        if role == "negative_control" and primary_eligible:
            raise RegistryError(
                f"variant negative control {model_id} cannot be primary eligible"
            )
        if role == "conditional_only":
            conditional_models.add(model_id)
            if models[model_id].status != "deferred":
                raise RegistryError(
                    f"conditional variant model {model_id} must remain deferred"
                )
        if is_mandatory_baseline:
            mandatory_baselines.add(model_id)
        records.append(canonicalize(record))

    record_ids = tuple(record["model_id"] for record in records)
    if record_ids != tuple(sorted(record_ids)):
        raise RegistryError("variant capability model records must be sorted by model_id")
    supported = {
        model.model_id
        for model in models.values()
        if task.task_id in model.supported_tasks
    }
    if seen != supported:
        missing = sorted(supported.difference(seen))
        extra = sorted(seen.difference(supported))
        raise RegistryError(
            "variant TaskSpec/model capability mismatch: "
            f"missing={missing}; extra={extra}"
        )
    if mandatory_baselines != set(task.baseline_model_ids):
        missing = sorted(set(task.baseline_model_ids).difference(mandatory_baselines))
        extra = sorted(mandatory_baselines.difference(task.baseline_model_ids))
        raise RegistryError(
            "variant mandatory baseline capability mismatch: "
            f"missing={missing}; extra={extra}"
        )
    if not primary_candidates or not primary_baselines:
        raise RegistryError(
            "variant signed-effect endpoint requires primary candidates and baselines"
        )
    if conditional_models != {"context_borzoi"}:
        raise RegistryError("variant conditional roster must contain only context_borzoi")
    if not _as_bool(root["admission_blocking"], f"{label}.admission_blocking"):
        raise RegistryError("variant output compatibility must remain admission-blocking")
    _as_string_tuple(root["blockers"], f"{label}.blockers", allow_empty=False)
    normalized = canonicalize(root)
    normalized["models"] = records
    return normalized


def _validate_variant_proxy_selection_policy(
    document: Mapping[str, Any],
    *,
    datasets: Mapping[str, DatasetManifest],
    task: TaskSpec,
) -> dict[str, Any]:
    """Freeze development proxies without relabeling them as eQTL outcomes."""

    expected = {
        "schema_version": "masld-bench-variant-proxy-selection-v1",
        "task_id": "variant_to_regulation",
        "policy_id": "variant_development_proxy_pareto_v1",
        "confirmation_endpoint_id": "signed_cell_type_eqtl_effect",
        "confirmation_evaluator_id": "variant_ld_block_fisher_z_spearman_gain_v1",
        "confirmation_dataset_id": "gse289173",
        "proxy_endpoints": [
            {
                "endpoint_id": "mpra_allelic_direction",
                "dataset_ids": ["gse281364"],
                "metric_id": "locus_heldout_allelic_spearman_v1",
                "inference_unit": "merged_locus",
                "gain_scale": "absolute_correlation",
                "minimum_gain": 0.02,
                "cross_fitted": True,
            },
            {
                "endpoint_id": "rna_atac_regulatory_profile",
                "dataset_ids": ["gse244832", "gse296875"],
                "metric_id": "two_way_profile_deviance_reduction_v1",
                "inference_unit": "donor+genomic_block",
                "gain_scale": "relative_deviance",
                "minimum_gain": 0.05,
                "cross_fitted": True,
            },
        ],
        "diagnostic_only_dataset_ids": ["gse281160"],
        "eqtl_specific_head_fitting_allowed": False,
        "weighted_composite_allowed": False,
        "selection_rule": (
            "family_top_two_union_pareto_union_one_standard_error_without_"
            "weighted_composite_v1"
        ),
        "evaluated_seed_count": 5,
        "minimum_qualifying_seed_count": 4,
        "minimum_study_count": 2,
        "sealed_results_used": False,
        "claim_boundary": (
            "Development proxy results are assay-native mechanistic selection "
            "evidence, not eQTL performance. Signed cell-type eQTL and ieQTL "
            "claims require the separately locked GSE289173 evaluator."
        ),
    }
    normalized = canonicalize(document)
    if normalized != canonicalize(expected):
        raise RegistryError(
            "variant development proxy-selection policy differs from its frozen contract"
        )
    if task.task_id != expected["task_id"]:
        raise RegistryError("variant proxy-selection policy has the wrong TaskSpec")
    if task.primary_evaluator_id != expected["confirmation_evaluator_id"]:
        raise RegistryError(
            "variant proxy confirmation evaluator differs from the TaskSpec"
        )
    if expected["confirmation_dataset_id"] not in task.datasets_sealed:
        raise RegistryError(
            "variant proxy confirmation dataset is not sealed in the TaskSpec"
        )
    proxy_datasets = {
        dataset_id
        for endpoint in expected["proxy_endpoints"]
        for dataset_id in endpoint["dataset_ids"]
    }
    if not proxy_datasets.issubset(set(task.datasets_train)):
        raise RegistryError(
            "variant proxy datasets must be training/development sources"
        )
    if set(expected["diagnostic_only_dataset_ids"]) != set(
        task.datasets_development
    ):
        raise RegistryError(
            "variant diagnostic-only datasets differ from TaskSpec development sources"
        )
    named_datasets = proxy_datasets | {
        expected["confirmation_dataset_id"],
        *expected["diagnostic_only_dataset_ids"],
    }
    if not named_datasets.issubset(datasets):
        raise RegistryError("variant proxy-selection policy names unknown datasets")
    return normalized


def _validate_cell_state_capability_registry(
    document: Mapping[str, Any],
    *,
    datasets: Mapping[str, DatasetManifest],
    models: Mapping[str, ModelManifest],
    task: TaskSpec,
) -> dict[str, Any]:
    label = "cell_state_capabilities"
    fields = frozenset(
        {
            "schema_version",
            "registry_id",
            "task_id",
            "sealed_dataset_id",
            "sealed_input_modalities",
            "sealed_output_definition",
            "rna_candidate_families",
            "rna_baseline_models",
            "atac_only_development_models",
            "atac_only_status",
            "admission_blocking",
            "blockers",
        }
    )
    root = _strict_mapping(
        document,
        allowed=fields,
        required=fields,
        label=label,
    )
    if root["schema_version"] != "masld-bench-cell-state-capability-v1":
        raise RegistryError("unsupported cell-state capability schema")
    if root["registry_id"] != "cell_state_mapping_model_capability":
        raise RegistryError("cell-state capability registry_id is not canonical")
    if root["task_id"] != task.task_id or task.task_id != "cell_state_mapping":
        raise RegistryError("cell-state capability task binding differs from TaskSpec")
    sealed_dataset_id = _as_string(
        root["sealed_dataset_id"], f"{label}.sealed_dataset_id"
    )
    if (
        sealed_dataset_id not in datasets
        or sealed_dataset_id not in task.datasets_sealed
    ):
        raise RegistryError(
            "cell-state capability sealed dataset differs from the TaskSpec"
        )
    sealed_modalities = _sorted_capability_ids(
        root, "sealed_input_modalities", label=label
    )
    if sealed_modalities != ("single_nucleus_rna",):
        raise RegistryError("cell-state sealed input must remain single-nucleus RNA")
    if not set(sealed_modalities).issubset(datasets[sealed_dataset_id].modalities):
        raise RegistryError("cell-state sealed input is absent from the sealed dataset")
    if set(task.input_modalities) != {"single_cell_rna", "single_nucleus_rna"}:
        raise RegistryError("cell-state TaskSpec must remain RNA-only")
    _as_string(root["sealed_output_definition"], f"{label}.sealed_output_definition")

    rna_families = _sorted_capability_ids(
        root, "rna_candidate_families", label=label
    )
    known_families = {model.family_id for model in models.values()}
    unknown_families = sorted(set(rna_families).difference(known_families))
    if unknown_families:
        raise RegistryError(
            "cell-state capability registry names unknown model families: "
            + ", ".join(unknown_families)
        )
    baselines = _sorted_capability_ids(root, "rna_baseline_models", label=label)
    atac_only = _sorted_capability_ids(
        root, "atac_only_development_models", label=label
    )
    named_models = set(baselines) | set(atac_only)
    unknown_models = sorted(named_models.difference(models))
    if unknown_models:
        raise RegistryError(
            "cell-state capability registry names unknown models: "
            + ", ".join(unknown_models)
        )
    if set(baselines) != set(task.baseline_model_ids):
        raise RegistryError("cell-state RNA baseline roster differs from TaskSpec")
    if set(baselines) & set(atac_only):
        raise RegistryError("cell-state RNA and ATAC capability categories overlap")

    supported = {
        model.model_id
        for model in models.values()
        if task.task_id in model.supported_tasks
    }
    rna_candidates = {
        model.model_id
        for model in models.values()
        if task.task_id in model.supported_tasks
        and model.family_id in set(rna_families)
    }
    for model_id in atac_only:
        model = models[model_id]
        if task.task_id in model.supported_tasks:
            raise RegistryError(
                f"ATAC-only model {model_id} cannot support the RNA-only cell-state task"
            )
        if not model.admission_blocking or model.execution.ready:
            raise RegistryError(
                f"ATAC-only cell-state model {model_id} must remain execution-blocked"
            )
    schedulable = rna_candidates | set(baselines)
    if schedulable != supported:
        missing = sorted(supported.difference(schedulable))
        extra = sorted(schedulable.difference(supported))
        raise RegistryError(
            "cell-state TaskSpec/model capability mismatch: "
            f"missing={missing}; extra={extra}"
        )
    if root["atac_only_status"] != "registered_not_schedulable_on_sealed_rna_task":
        raise RegistryError("cell-state ATAC-only scheduling status is not fail-closed")
    if not _as_bool(root["admission_blocking"], f"{label}.admission_blocking"):
        raise RegistryError("cell-state modality boundary must remain admission-blocking")
    _as_string_tuple(root["blockers"], f"{label}.blockers", allow_empty=False)
    return canonicalize(root)


def _validate_dataset_view(
    *,
    view: DatasetViewContract,
    parent: DatasetManifest,
    parent_path: Path,
    config_root: Path,
) -> None:
    """Recursively verify one complete smoke view and its parent binding."""

    if sha256_file(parent_path) != view.parent_registry_sha256:
        raise RegistryError(
            f"dataset view {view.view_id} parent registry binding changed"
        )
    if not set(view.modalities).issubset(parent.modalities):
        raise RegistryError(
            f"dataset view {view.view_id} adds a modality absent from its parent"
        )
    if not set(view.pairing_levels).issubset(parent.pairing_levels):
        raise RegistryError(
            f"dataset view {view.view_id} adds a pairing topology absent from its parent"
        )
    if view.label_visibility != parent.split.label_visibility:
        raise RegistryError(
            f"dataset view {view.view_id} changes parent label visibility"
        )
    if (
        parent.role.value not in {"train_development", "external_development"}
        or parent.status.value != "available"
    ):
        raise RegistryError(
            f"dataset view {view.view_id} parent must be an available development dataset"
        )

    direct = {
        "manifest": view.artifact_manifest,
        "complete receipt": view.complete_receipt,
        "data": view.data_artifact,
        "selection": view.selection_artifact,
        "ontology": view.ontology_artifact,
        "derivation": view.derivation_artifact,
    }
    paths = {
        label: _validate_ready_evidence_artifact(
            artifact,
            config_root=config_root,
            label=f"dataset view {view.view_id} {label}",
        )
        for label, artifact in direct.items()
    }
    for artifact in view.authority_artifacts:
        authority_path = _validate_ready_evidence_artifact(
            artifact,
            config_root=config_root,
            label=f"dataset view {view.view_id} authority {artifact.role}",
        )
        assert artifact.role is not None
        authority = artifact.role.rsplit(":", 1)[1]
        document = _load_json_object(
            authority_path, f"dataset view {view.view_id} {authority} authority"
        )
        if set(document) != {
            "schema_version",
            "view_id",
            "parent_dataset_id",
            "authority",
            "admitted",
            "evidence",
        }:
            raise RegistryError(
                f"dataset view {view.view_id} {authority} authority has an invalid schema"
            )
        if (
            document["schema_version"] != "masld-bench-dataset-view-authority-v1"
            or document["view_id"] != view.view_id
            or document["parent_dataset_id"] != view.parent_dataset_id
            or document["authority"] != authority
            or document["admitted"] is not True
            or not isinstance(document["evidence"], Mapping)
            or not document["evidence"]
        ):
            raise RegistryError(
                f"dataset view {view.view_id} {authority} authority does not admit this view"
            )

    artifact_manifest = _load_json_object(
        paths["manifest"], f"dataset view {view.view_id} artifact manifest"
    )
    try:
        recursively_verified_manifest = verify_frozen_tree(paths["manifest"].parent)
    except ArtifactError as error:
        raise RegistryError(
            f"dataset view {view.view_id} is not an intact frozen tree: {error}"
        ) from error
    if canonicalize(recursively_verified_manifest) != canonicalize(artifact_manifest):
        raise RegistryError(
            f"dataset view {view.view_id} recursive manifest verification changed content"
        )
    expected_records = {
        artifact.path: {
            "path": Path(artifact.path).name,
            "sha256": artifact.sha256,
            "size_bytes": artifact.size_bytes,
        }
        for artifact in (
            view.data_artifact,
            view.selection_artifact,
            view.ontology_artifact,
            view.derivation_artifact,
        )
    }
    observed_records = artifact_manifest.get("artifacts")
    artifact_metadata = artifact_manifest.get("metadata", {})
    artifact_class = artifact_metadata.get("artifact_class")
    artifact_view_id = artifact_metadata.get(
        "subset_id", artifact_metadata.get("view_id")
    )
    if (
        artifact_manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or artifact_class
        not in {
            "geneformer_smoke_subset",
            "multimodal_smoke_subset",
            "unlabeled_snrna_smoke_subset",
        }
        or artifact_view_id != view.view_id
        or not isinstance(observed_records, list)
        or sorted(observed_records, key=lambda item: str(item.get("path")))
        != sorted(expected_records.values(), key=lambda item: str(item["path"]))
    ):
        raise RegistryError(
            f"dataset view {view.view_id} artifact inventory differs from its contract"
        )
    complete = _load_json_object(
        paths["complete receipt"], f"dataset view {view.view_id} complete receipt"
    )
    expected_complete = {
        "schema_version": "masld-bench-complete-v1",
        "manifest_sha256": view.artifact_manifest.sha256,
        "artifact_count": len(expected_records),
    }
    if canonicalize(complete) != canonicalize(expected_complete):
        raise RegistryError(
            f"dataset view {view.view_id} completion receipt differs from its manifest"
        )

    derivation = _load_json_object(
        paths["derivation"], f"dataset view {view.view_id} derivation receipt"
    )
    selection = derivation.get("selection")
    matrix = derivation.get("matrix")
    artifacts = derivation.get("artifacts")
    if not isinstance(matrix, Mapping) or not isinstance(artifacts, Mapping):
        raise RegistryError(
            f"dataset view {view.view_id} derivation lacks matrix or artifact records"
        )
    derivation_schema = derivation.get("schema_version")
    if derivation_schema == "masld-bench-geneformer-smoke-subset-v1":
        input_ready = derivation.get("geneformer_input_contract", {}).get("ready")
        row_shapes = [matrix.get("shape", [None])[0]]
        data_sha256 = artifacts.get("h5ad", {}).get("sha256")
        expected_artifact_class = "geneformer_smoke_subset"
    elif derivation_schema == "masld-bench-multimodal-smoke-subset-v1":
        input_ready = derivation.get("multimodal_input_contract", {}).get(
            "ready"
        )
        row_shapes = [
            matrix.get(modality, {}).get("shape", [None])[0]
            for modality in ("rna", "atac")
        ]
        data_sha256 = artifacts.get("multimodal_h5", {}).get("sha256")
        expected_artifact_class = "multimodal_smoke_subset"
    elif derivation_schema == "masld-bench-unlabeled-snrna-smoke-subset-v1":
        input_ready = derivation.get("unlabeled_input_contract", {}).get("ready")
        row_shapes = [matrix.get("shape", [None])[0]]
        data_sha256 = artifacts.get("h5ad", {}).get("sha256")
        expected_artifact_class = "unlabeled_snrna_smoke_subset"
    else:
        raise RegistryError(
            f"dataset view {view.view_id} has an unsupported derivation schema"
        )
    if (
        artifact_class != expected_artifact_class
        or derivation.get("subset_id") != view.view_id
        or derivation.get("dataset_id") != view.parent_dataset_id
        or input_ready is not True
        or not isinstance(selection, Mapping)
        or selection.get("cell_budget") != view.row_count
        or selection.get("selection_policy") != view.selection_policy
        or selection.get("selection_seed") != view.selection_seed
        or selection.get("selection_outcomes_used") is not False
        or selection.get("sealed_outcomes_used") is not False
        or canonicalize(selection.get("counts_by_class"))
        != canonicalize(view.class_counts)
        or any(row_count != view.row_count for row_count in row_shapes)
        or not isinstance(artifacts, Mapping)
        or data_sha256 != view.data_artifact.sha256
        or artifacts.get("selection", {}).get("sha256")
        != view.selection_artifact.sha256
        or artifacts.get("ontology", {}).get("sha256")
        != view.ontology_artifact.sha256
    ):
        raise RegistryError(
            f"dataset view {view.view_id} derivation receipt differs from its contract"
        )

    ontology = _load_json_object(
        paths["ontology"], f"dataset view {view.view_id} ontology"
    )
    if (
        ontology.get("schema_version") != "masld-bench-cell-ontology-v1"
        or set(ontology.get("classes", {})) != set(view.class_counts)
        or ontology.get("selection_outcomes_used") is not False
        or ontology.get("sealed_outcomes_used") is not False
    ):
        raise RegistryError(
            f"dataset view {view.view_id} ontology differs from its class contract"
        )

    try:
        with paths["selection"].open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except OSError as error:
        raise RegistryError(
            f"dataset view {view.view_id} selection is unreadable: {error}"
        ) from error
    required = {"cell_id", "donor_id", "broad_label"}
    if not rows or not required.issubset(rows[0]):
        raise RegistryError(
            f"dataset view {view.view_id} selection lacks required columns"
        )
    observed_class_counts: dict[str, int] = {}
    for row in rows:
        label = row["broad_label"]
        observed_class_counts[label] = observed_class_counts.get(label, 0) + 1
    if (
        len(rows) != view.row_count
        or len({row["cell_id"] for row in rows}) != view.row_count
        or len({row["donor_id"] for row in rows}) != view.biological_unit_count
        or canonicalize(observed_class_counts) != canonicalize(view.class_counts)
    ):
        raise RegistryError(
            f"dataset view {view.view_id} selection counts differ from its contract"
        )


def load_prediction_bundle(path: str | Path) -> PredictionBundle:
    return load_contract(path, PredictionBundle)


def load_selection_lock(path: str | Path) -> SelectionLock:
    return load_contract(path, SelectionLock)


def _registry_files(directory: Path) -> tuple[Path, ...]:
    if not directory.is_dir():
        raise RegistryError(f"required registry directory does not exist: {directory}")
    unsupported = sorted(
        path for path in directory.iterdir() if path.is_file() and path.suffix not in {".toml", ".json"}
    )
    if unsupported:
        raise RegistryError(
            "unsupported files in registry directory: "
            + ", ".join(path.name for path in unsupported)
        )
    return tuple(
        sorted(
            path
            for path in directory.iterdir()
            if path.is_file() and path.suffix in {".toml", ".json"}
        )
    )


def _insert_unique(
    target: dict[str, Any], identifier: str, value: Any, *, source: Path
) -> None:
    if identifier in target:
        raise RegistryError(f"duplicate registry identifier {identifier!r} at {source}")
    target[identifier] = value


def _load_model_family(path: Path) -> tuple[ModelManifest, ...]:
    document = _load_document(path)
    root = _strict_mapping(
        document,
        allowed=ModelManifest._FAMILY_FIELDS,
        required=ModelManifest._FAMILY_FIELDS,
        label=f"model family {path.name}",
    )
    models = root["models"]
    if not isinstance(models, list) or not models:
        raise RegistryError(f"model family has no [[models]] entries: {path}")
    try:
        result = tuple(ModelManifest.from_family_entry(root, model) for model in models)
    except ContractError as error:
        raise RegistryError(f"invalid model family {path}: {error}") from error
    return result


_RESOURCE_ROOT_FIELDS = frozenset(
    {
        "schema_version",
        "reference",
        "coordinate_contract",
        "semantics",
        "prohibited_references",
        "redistribution",
        "cluster",
        "profiles",
    }
)


def _validate_resources(document: Mapping[str, Any], path: Path) -> ReferenceBundle:
    root = _strict_mapping(
        document,
        allowed=_RESOURCE_ROOT_FIELDS,
        required=_RESOURCE_ROOT_FIELDS,
        label=f"resource registry {path.name}",
    )
    _as_string(root["schema_version"], "resources.schema_version")

    semantics_fields = frozenset(
        {
            "pairing_levels",
            "modality_statuses",
            "exposure_statuses",
            "evidence_statuses",
            "missingness_policy",
        }
    )
    semantics = _strict_mapping(
        root["semantics"],
        allowed=semantics_fields,
        required=semantics_fields,
        label="resources.semantics",
    )
    pairing_values = _as_string_tuple(
        semantics["pairing_levels"], "resources.semantics.pairing_levels", allow_empty=False
    )
    missing_values = _as_string_tuple(
        semantics["modality_statuses"],
        "resources.semantics.modality_statuses",
        allow_empty=False,
    )
    exposure_values = _as_string_tuple(
        semantics["exposure_statuses"],
        "resources.semantics.exposure_statuses",
        allow_empty=False,
    )
    for value in pairing_values:
        try:
            PairingState(value)
        except ValueError as error:
            raise RegistryError(f"unknown pairing level in resources: {value}") from error
    for value in missing_values:
        try:
            MissingState(value)
        except ValueError as error:
            raise RegistryError(f"unknown modality status in resources: {value}") from error
    for value in exposure_values:
        try:
            ExposureState(value)
        except ValueError as error:
            raise RegistryError(f"unknown exposure status in resources: {value}") from error
    for label, observed, enum_type in (
        ("pairing_levels", pairing_values, PairingState),
        ("modality_statuses", missing_values, MissingState),
        ("exposure_statuses", exposure_values, ExposureState),
    ):
        expected = {item.value for item in enum_type}
        if set(observed) != expected:
            raise RegistryError(
                f"resources.semantics.{label} must enumerate the exact contract domain"
            )
    _as_string_tuple(
        semantics["evidence_statuses"],
        "resources.semantics.evidence_statuses",
        allow_empty=False,
    )
    _as_string(semantics["missingness_policy"], "resources.semantics.missingness_policy")

    redistribution_fields = frozenset(
        {"default_data_policy", "human_raw_data_policy", "controlled_data_policy", "model_weight_policy"}
    )
    redistribution = _strict_mapping(
        root["redistribution"],
        allowed=redistribution_fields,
        required=redistribution_fields,
        label="resources.redistribution",
    )
    for key, value in redistribution.items():
        _as_string(value, f"resources.redistribution.{key}")

    cluster_fields = frozenset(
        {"default_qos", "aggregate_memory_cap_gb", "array_submission_default", "production_submit_default"}
    )
    cluster = _strict_mapping(
        root["cluster"],
        allowed=cluster_fields,
        required=cluster_fields,
        label="resources.cluster",
    )
    _as_string(cluster["default_qos"], "resources.cluster.default_qos")
    _as_int(
        cluster["aggregate_memory_cap_gb"],
        "resources.cluster.aggregate_memory_cap_gb",
        minimum=1,
    )
    _as_bool(cluster["array_submission_default"], "resources.cluster.array_submission_default")
    _as_bool(cluster["production_submit_default"], "resources.cluster.production_submit_default")

    profiles = root["profiles"]
    if not isinstance(profiles, Mapping) or not profiles:
        raise RegistryError("resources.profiles must be a non-empty table")
    profile_fields = frozenset(
        {
            "partition",
            "cpus",
            "memory_gb",
            "wall_time",
            "gpus",
            "accelerator",
            "qos",
            "admission_blocking",
            "blockers",
        }
    )
    profile_required = profile_fields - {"admission_blocking", "blockers"}
    for profile_id, profile_value in profiles.items():
        _as_string(profile_id, "resources profile id")
        profile = _strict_mapping(
            profile_value,
            allowed=profile_fields,
            required=profile_required,
            label=f"resources.profiles.{profile_id}",
        )
        for key in ("partition", "wall_time", "accelerator", "qos"):
            _as_string(profile[key], f"resources.profiles.{profile_id}.{key}")
        for key, minimum in (("cpus", 1), ("memory_gb", 1), ("gpus", 0)):
            _as_int(profile[key], f"resources.profiles.{profile_id}.{key}", minimum=minimum)
        admission_blocking = _as_bool(
            profile.get("admission_blocking", False),
            f"resources.profiles.{profile_id}.admission_blocking",
        )
        blockers = _as_string_tuple(
            profile.get("blockers", ()),
            f"resources.profiles.{profile_id}.blockers",
        )
        if profile["accelerator"] == "UNRESOLVED" and (
            not admission_blocking or not blockers
        ):
            raise RegistryError(
                f"resources.profiles.{profile_id} has an UNRESOLVED accelerator "
                "but is not explicitly admission-blocking with a blocker"
            )

    try:
        return ReferenceBundle.from_dict(root)
    except ReferenceError as error:
        raise RegistryError(f"invalid reference bundle in {path}: {error}") from error


class Registry:
    """A read-only, cross-referenced snapshot of datasets, models, and tasks."""

    __slots__ = (
        "root",
        "datasets",
        "dataset_views",
        "models",
        "tasks",
        "splits",
        "reference_bundle",
        "resources",
        "cell_state_capabilities",
        "rna_atac_capabilities",
        "variant_capabilities",
        "variant_proxy_selection",
        "source_files",
    )

    def __init__(
        self,
        *,
        root: Path,
        datasets: Mapping[str, DatasetManifest],
        dataset_views: Mapping[str, DatasetViewContract],
        models: Mapping[str, ModelManifest],
        tasks: Mapping[str, TaskSpec],
        splits: Mapping[str, SplitSpec],
        reference_bundle: ReferenceBundle,
        resources: Mapping[str, Any],
        cell_state_capabilities: Mapping[str, Any] | None,
        rna_atac_capabilities: Mapping[str, Any] | None,
        variant_capabilities: Mapping[str, Any] | None,
        variant_proxy_selection: Mapping[str, Any] | None,
        source_files: tuple[Path, ...],
    ) -> None:
        self.root = root.resolve()
        self.datasets = MappingProxyType(dict(sorted(datasets.items())))
        self.dataset_views = MappingProxyType(dict(sorted(dataset_views.items())))
        self.models = MappingProxyType(dict(sorted(models.items())))
        self.tasks = MappingProxyType(dict(sorted(tasks.items())))
        self.splits = MappingProxyType(dict(sorted(splits.items())))
        self.reference_bundle = reference_bundle
        normalized_resources = canonicalize(resources)
        self.resources = MappingProxyType(normalized_resources)
        self.cell_state_capabilities = (
            None
            if cell_state_capabilities is None
            else MappingProxyType(canonicalize(cell_state_capabilities))
        )
        self.rna_atac_capabilities = (
            None
            if rna_atac_capabilities is None
            else MappingProxyType(canonicalize(rna_atac_capabilities))
        )
        self.variant_capabilities = (
            None
            if variant_capabilities is None
            else MappingProxyType(canonicalize(variant_capabilities))
        )
        self.variant_proxy_selection = (
            None
            if variant_proxy_selection is None
            else MappingProxyType(canonicalize(variant_proxy_selection))
        )
        self.source_files = tuple(sorted(path.resolve() for path in source_files))

    @classmethod
    def load(
        cls,
        root: str | Path,
        *,
        validate_references: bool = False,
    ) -> "Registry":
        config_root = Path(root)
        if not config_root.is_dir():
            raise RegistryError(f"configuration root does not exist: {config_root}")

        datasets: dict[str, DatasetManifest] = {}
        dataset_views: dict[str, DatasetViewContract] = {}
        models: dict[str, ModelManifest] = {}
        tasks: dict[str, TaskSpec] = {}
        splits: dict[str, SplitSpec] = {}
        source_files: list[Path] = []

        for path in _registry_files(config_root / "datasets"):
            dataset = load_dataset_manifest(path)
            if path.stem != dataset.dataset_id:
                raise RegistryError(
                    f"dataset filename {path.stem!r} does not match dataset_id {dataset.dataset_id!r}"
                )
            _insert_unique(datasets, dataset.dataset_id, dataset, source=path)
            if dataset.activation.ready:
                for artifact in dataset.activation.evidence_artifacts:
                    _validate_ready_evidence_artifact(
                        artifact,
                        config_root=config_root,
                        label=f"dataset {dataset.dataset_id} authority {artifact.role}",
                    )
            source_files.append(path)

        dataset_view_root = config_root / "dataset_views"
        if dataset_view_root.exists():
            for path in _registry_files(dataset_view_root):
                view = load_dataset_view_contract(path)
                if path.stem != view.view_id:
                    raise RegistryError(
                        f"dataset view filename {path.stem!r} does not match "
                        f"view_id {view.view_id!r}"
                    )
                if view.parent_dataset_id not in datasets:
                    raise RegistryError(
                        f"dataset view {view.view_id} names unknown parent "
                        f"{view.parent_dataset_id}"
                    )
                _insert_unique(dataset_views, view.view_id, view, source=path)
                _validate_dataset_view(
                    view=view,
                    parent=datasets[view.parent_dataset_id],
                    parent_path=config_root
                    / "datasets"
                    / f"{view.parent_dataset_id}.toml",
                    config_root=config_root,
                )
                source_files.append(path)

        for path in _registry_files(config_root / "models"):
            for model in _load_model_family(path):
                _insert_unique(models, model.model_id, model, source=path)
                if model.execution.ready:
                    authority_paths: dict[str, Path] = {}
                    for artifact in model.execution.evidence_artifacts:
                        verified_path = _validate_ready_evidence_artifact(
                            artifact,
                            config_root=config_root,
                            label=f"model {model.model_id} authority {artifact.role}",
                        )
                        assert artifact.role is not None
                        authority_paths[artifact.role.rsplit(":", 1)[1]] = verified_path
                    for authority in (
                        "code_license",
                        "weights_license",
                        "derivative_weights_license",
                    ):
                        _validate_model_license_authority(
                            path=authority_paths[authority],
                            model_id=model.model_id,
                            authority=authority,
                            declared_license=getattr(model.execution, authority),
                        )
            source_files.append(path)

        for path in _registry_files(config_root / "tasks"):
            task = load_task_spec(path)
            if path.stem != task.task_id:
                raise RegistryError(
                    f"task filename {path.stem!r} does not match task_id {task.task_id!r}"
                )
            _insert_unique(tasks, task.task_id, task, source=path)
            source_files.append(path)

        for path in _registry_files(config_root / "splits"):
            split = load_split_spec(path)
            if path.stem != split.split_id:
                raise RegistryError(
                    f"split filename {path.stem!r} does not match split_id {split.split_id!r}"
                )
            _insert_unique(splits, split.split_id, split, source=path)
            source_files.append(path)

        resources_path = config_root / "resources.toml"
        resources = _load_document(resources_path)
        reference_bundle = _validate_resources(resources, resources_path)
        _set_reference_base(reference_bundle, resources_path.parent)
        source_files.append(resources_path)

        missing_datasets: list[str] = []
        for task in tasks.values():
            missing_datasets.extend(
                f"{task.task_id}:{dataset_id}"
                for dataset_id in task.dataset_ids
                if dataset_id not in datasets
            )
        if missing_datasets:
            raise RegistryError(
                "task registry references unknown datasets: " + ", ".join(missing_datasets)
            )

        missing_splits = sorted(
            f"{task.task_id}:{task.split_id}"
            for task in tasks.values()
            if task.split_id not in splits
        )
        if missing_splits:
            raise RegistryError(
                "task registry references unknown splits: " + ", ".join(missing_splits)
            )

        missing_tasks = sorted(
            {
                f"{model.model_id}:{task_id}"
                for model in models.values()
                for task_id in model.supported_tasks
                if task_id not in tasks
            }
        )
        if missing_tasks:
            raise RegistryError(
                "model registry references unknown tasks: " + ", ".join(missing_tasks)
            )

        missing_baselines = sorted(
            {
                f"{task.task_id}:{model_id}"
                for task in tasks.values()
                for model_id in task.baseline_model_ids
                if model_id not in models
                or task.task_id not in models[model_id].supported_tasks
            }
        )
        if missing_baselines:
            raise RegistryError(
                "task baseline roster references a missing or incompatible model: "
                + ", ".join(missing_baselines)
            )

        capability_path = (
            config_root / "evaluation" / "rna_conditioned_atac_capabilities.toml"
        )
        has_capability_task = "rna_conditioned_atac" in tasks
        has_capability_file = capability_path.is_file()
        if has_capability_task != has_capability_file:
            raise RegistryError(
                "RNA-ATAC TaskSpec and capability registry must be registered together"
            )
        rna_atac_capabilities: dict[str, Any] | None = None
        if has_capability_task:
            rna_atac_capabilities = _validate_rna_atac_capability_registry(
                _load_document(capability_path),
                models=models,
                task=tasks["rna_conditioned_atac"],
            )
            source_files.append(capability_path)

        cell_capability_path = (
            config_root / "evaluation" / "cell_state_mapping_capabilities.toml"
        )
        has_cell_task = "cell_state_mapping" in tasks
        has_cell_capability_file = cell_capability_path.is_file()
        if has_cell_capability_file and not has_cell_task:
            raise RegistryError(
                "cell-state capability registry requires its TaskSpec"
            )
        cell_state_capabilities: dict[str, Any] | None = None
        if has_cell_capability_file:
            cell_state_capabilities = _validate_cell_state_capability_registry(
                _load_document(cell_capability_path),
                datasets=datasets,
                models=models,
                task=tasks["cell_state_mapping"],
            )
            source_files.append(cell_capability_path)

        variant_capability_path = (
            config_root / "evaluation" / "variant_to_regulation_capabilities.toml"
        )
        has_variant_task = "variant_to_regulation" in tasks
        has_variant_capability_file = variant_capability_path.is_file()
        if has_variant_task != has_variant_capability_file:
            raise RegistryError(
                "variant TaskSpec and capability registry must be registered together"
            )
        variant_capabilities: dict[str, Any] | None = None
        variant_proxy_selection: dict[str, Any] | None = None
        if has_variant_task:
            variant_capabilities = _validate_variant_capability_registry(
                _load_document(variant_capability_path),
                models=models,
                task=tasks["variant_to_regulation"],
            )
            source_files.append(variant_capability_path)
            proxy_policy_path = (
                config_root
                / "evaluation"
                / "variant_development_proxy_selection.json"
            )
            if not proxy_policy_path.is_file():
                raise RegistryError(
                    "variant TaskSpec requires its development proxy-selection policy"
                )
            variant_proxy_selection = _validate_variant_proxy_selection_policy(
                _load_document(proxy_policy_path),
                datasets=datasets,
                task=tasks["variant_to_regulation"],
            )
            source_files.append(proxy_policy_path)

        promotion_path = config_root / "evaluation" / "promotion_gates.toml"
        promotion_document = _load_document(promotion_path)
        promotion_tasks = promotion_document.get("tasks")
        if not isinstance(promotion_tasks, Mapping):
            raise RegistryError("promotion gate registry must contain a tasks table")
        promotion_sha256 = sha256_file(promotion_path)
        for task in tasks.values():
            if task.bootstrap_replicates != 10_000:
                raise RegistryError(
                    f"task {task.task_id} must freeze exactly 10,000 bootstrap replicates"
                )
            if task.promotion_gate_config_sha256 != promotion_sha256:
                raise RegistryError(
                    f"task {task.task_id} does not bind the current promotion-gate registry"
                )
            gate = promotion_tasks.get(task.task_id)
            if not isinstance(gate, Mapping):
                aliases = [
                    candidate
                    for candidate in promotion_tasks.values()
                    if isinstance(candidate, Mapping)
                    and candidate.get("promotion_gate_id") == task.promotion_gate_id
                ]
                if (
                    task.primary_evaluator_id != "not_applicable_nonpromotable_v1"
                    or len(aliases) != 1
                    or aliases[0].get("allowed_promotion_modes") != ["exploratory_only"]
                    or aliases[0].get("sealed_holdout_dataset_ids") != []
                ):
                    raise RegistryError(
                        f"promotion gate registry has no task entry for {task.task_id}"
                    )
                gate = aliases[0]
            if gate.get("promotion_gate_id") != task.promotion_gate_id:
                raise RegistryError(
                    f"task {task.task_id} promotion_gate_id differs from its executable gate"
                )
            observed_evaluator = {
                "primary_evaluator_id": task.primary_evaluator_id,
                "evaluator_parameters": task.to_dict()["evaluator_parameters"],
            }
            roster_contract = _ROSTER_TASK_EVALUATORS.get(task.task_id)
            if roster_contract is not None:
                evaluator_id, roster_field = roster_contract
                parameters = observed_evaluator["evaluator_parameters"]
                unresolved_parameters = {
                    roster_field: [],
                    "roster_status": "UNRESOLVED",
                }
                if task.primary_evaluator_id == "UNRESOLVED":
                    if canonicalize(parameters) != canonicalize(unresolved_parameters):
                        raise RegistryError(
                            f"task {task.task_id} unresolved roster contract is ambiguous"
                        )
                else:
                    if task.primary_evaluator_id != evaluator_id or not isinstance(
                        parameters, Mapping
                    ) or set(parameters) != {roster_field, "roster_authority"}:
                        raise RegistryError(
                            f"task {task.task_id} differs from its frozen roster evaluator contract"
                        )
                    roster = parameters[roster_field]
                    if (
                        not isinstance(roster, list)
                        or not roster
                        or any(not isinstance(item, str) or not item for item in roster)
                        or roster != sorted(set(roster))
                    ):
                        raise RegistryError(
                            f"task {task.task_id} {roster_field} must be non-empty, unique, and sorted"
                        )
                    try:
                        roster_authority = ArtifactRef.from_dict(
                            parameters["roster_authority"]
                        )
                    except ContractError as error:
                        raise RegistryError(
                            f"task {task.task_id} roster authority is invalid: {error}"
                        ) from error
                    expected_role = f"task_evaluator_roster:{task.task_id}"
                    if roster_authority.role != expected_role:
                        raise RegistryError(
                            f"task {task.task_id} roster authority must use role {expected_role}"
                        )
                    roster_path = _validate_ready_evidence_artifact(
                        roster_authority,
                        config_root=config_root,
                        label=f"task {task.task_id} roster authority",
                    )
                    try:
                        roster_document = json.loads(
                            roster_path.read_text(encoding="utf-8")
                        )
                    except (OSError, json.JSONDecodeError) as error:
                        raise RegistryError(
                            f"task {task.task_id} roster authority is not valid JSON: {error}"
                        ) from error
                    expected_roster_document = {
                        "schema_version": "masld-bench-evaluator-roster-v1",
                        "task_id": task.task_id,
                        "primary_evaluator_id": task.primary_evaluator_id,
                        "roster_field": roster_field,
                        "roster": roster,
                    }
                    if canonicalize(roster_document) != canonicalize(
                        expected_roster_document
                    ):
                        raise RegistryError(
                            f"task {task.task_id} roster authority content differs from TaskSpec"
                        )
            else:
                expected_evaluator = _EXPECTED_TASK_EVALUATORS.get(
                    task.task_id,
                    {
                        "primary_evaluator_id": "not_applicable_nonpromotable_v1",
                        "evaluator_parameters": {},
                    },
                )
                if canonicalize(observed_evaluator) != canonicalize(expected_evaluator):
                    raise RegistryError(
                        f"task {task.task_id} differs from its frozen primary evaluator contract"
                    )
        source_files.append(promotion_path)

        registry = cls(
            root=config_root,
            datasets=datasets,
            dataset_views=dataset_views,
            models=models,
            tasks=tasks,
            splits=splits,
            reference_bundle=reference_bundle,
            resources=resources,
            cell_state_capabilities=cell_state_capabilities,
            rna_atac_capabilities=rna_atac_capabilities,
            variant_capabilities=variant_capabilities,
            variant_proxy_selection=variant_proxy_selection,
            source_files=tuple(source_files),
        )
        if validate_references:
            registry.reference_bundle.validate()
        return registry

    def snapshot(self) -> dict[str, Any]:
        return {
            "schema_version": "masld-bench-registry-snapshot-v1",
            "datasets": [self.datasets[key].to_dict() for key in sorted(self.datasets)],
            "dataset_views": [
                self.dataset_views[key].to_dict() for key in sorted(self.dataset_views)
            ],
            "models": [self.models[key].to_dict() for key in sorted(self.models)],
            "tasks": [self.tasks[key].to_dict() for key in sorted(self.tasks)],
            "splits": [self.splits[key].to_dict() for key in sorted(self.splits)],
            "reference_bundle": self.reference_bundle.to_dict(),
            "resources": canonicalize(self.resources),
            "cell_state_capabilities": (
                None
                if self.cell_state_capabilities is None
                else canonicalize(self.cell_state_capabilities)
            ),
            "rna_atac_capabilities": (
                None
                if self.rna_atac_capabilities is None
                else canonicalize(self.rna_atac_capabilities)
            ),
            "variant_capabilities": (
                None
                if self.variant_capabilities is None
                else canonicalize(self.variant_capabilities)
            ),
            "variant_proxy_selection": (
                None
                if self.variant_proxy_selection is None
                else canonicalize(self.variant_proxy_selection)
            ),
        }

    @property
    def snapshot_sha256(self) -> str:
        return canonical_sha256(self.snapshot())

    def get_dataset(self, dataset_id: str) -> DatasetManifest:
        try:
            return self.datasets[dataset_id]
        except KeyError as error:
            raise RegistryError(f"unknown dataset_id: {dataset_id}") from error

    def get_dataset_view(self, view_id: str) -> DatasetViewContract:
        try:
            return self.dataset_views[view_id]
        except KeyError as error:
            raise RegistryError(f"unknown dataset view_id: {view_id}") from error

    def get_model(self, model_id: str) -> ModelManifest:
        try:
            return self.models[model_id]
        except KeyError as error:
            raise RegistryError(f"unknown model_id: {model_id}") from error

    def get_task(self, task_id: str) -> TaskSpec:
        try:
            return self.tasks[task_id]
        except KeyError as error:
            raise RegistryError(f"unknown task_id: {task_id}") from error

    def get_split(self, split_id: str) -> SplitSpec:
        try:
            return self.splits[split_id]
        except KeyError as error:
            raise RegistryError(f"unknown split_id: {split_id}") from error


def _set_reference_base(bundle: ReferenceBundle, root: Path) -> None:
    object.__setattr__(bundle, "_base_dir", root.resolve())


__all__ = [
    "Registry",
    "RegistryError",
    "load_contract",
    "load_dataset_manifest",
    "load_dataset_view_contract",
    "load_model_manifest",
    "load_prediction_bundle",
    "load_run_spec",
    "load_selection_lock",
    "load_split_spec",
    "load_task_spec",
]
