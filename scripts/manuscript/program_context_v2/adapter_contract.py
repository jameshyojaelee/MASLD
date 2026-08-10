#!/usr/bin/env python3
"""Frozen, exact input/producer contract for Plan 60 presentation adapters.

The coordinator, REL01, and REL05 import this one contract.  Provenance is
therefore not an open-ended list of files that happened to be recorded: every
adapter must bind exactly the paths declared here, with no omissions, extras,
or duplicates.  File hashes and byte sizes remain run-specific and are stored
in ``adapter_provenance.tsv``.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping

from release_common import CANDIDATE_ID, ReleaseContractError


ADAPTER_PROVENANCE_VERSION = "plan60_adapter_provenance_v2"
ADAPTER_PROVENANCE_FIELDS = (
    "contract_version",
    "adapter_contract_sha256",
    "artifact_id",
    "output_relative_path",
    "output_sha256",
    "output_bytes",
    "producer_manifest_relative_path",
    "producer_manifest_sha256",
    "producer_bindings_json",
    "input_bindings_json",
    "canonical_promotion_authorized",
)
RECURSIVE_PRODUCER_FIELDS = (
    "producer_id",
    "repository_path",
    "sha256",
    "bytes",
)

_PROGRAM_CANDIDATE = "Analysis/Multimodal_Program_Projection/candidates/" + CANDIDATE_ID
_HOTSPOT = _PROGRAM_CANDIDATE + "/hotspot"
_NMF = _HOTSPOT + "/nmf_continuous_supplement"
_PASSPORT = "RNA-seq/results/evidence_passports/candidates/" + CANDIDATE_ID
_INTEGRATION = "RNA-seq/Human/Patient_Cohorts/analysis/integration"

_COMMON_PRODUCERS = (
    "scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py",
    "scripts/manuscript/program_context_v2/coordinator/prepare_real_release.py",
)
NMF_FREEZER_PRODUCER_PATH = (
    "Analysis/SingleCell/scripts/hotspot_modules/"
    "519_freeze_nmf_continuous_supplement.py"
)

# Keep tuples ordered for deterministic provenance serialization.  Equality
# checks below deliberately treat them as sets after separately rejecting
# duplicates.
ADAPTER_CONTRACT: Mapping[str, Mapping[str, object]] = {
    "cohort_overview": {
        "output_relative_path": "cohort_overview.tsv",
        "required_inputs": (
            _INTEGRATION + "/metadata/unified_metadata.csv",
            _INTEGRATION + "/qc/sample_qc_report.csv",
        ),
        "required_producers": _COMMON_PRODUCERS,
    },
    "sample_composition": {
        "output_relative_path": "sample_composition.tsv",
        "required_inputs": (
            _HOTSPOT + "/composition_sample_qc_v2.tsv",
            _HOTSPOT + "/composition_sample_qc_testability.tsv",
            _HOTSPOT + "/COMPOSITION_QC_READY",
        ),
        "required_producers": _COMMON_PRODUCERS,
    },
    "hotspot_programs": {
        "output_relative_path": "hotspot_programs.tsv",
        "required_inputs": (
            _HOTSPOT + "/program_registry_v2.tsv",
            _HOTSPOT + "/fig2_program_source.tsv",
            _HOTSPOT + "/program_registry_v2_semantic_adjudication.tsv",
            _HOTSPOT + "/stage_dataset_design_audit.tsv",
            _HOTSPOT + "/cohort_and_lodo_effects.tsv",
            _NMF + "/source_manifest.tsv",
            _NMF + "/nmf_continuous_loadings.tsv",
            _NMF + "/nmf_continuous_supplement.tsv",
            _NMF + "/validation_report.tsv",
            _NMF + "/producer_manifest.tsv",
            _NMF + "/NMF_CONTINUOUS_SUPPLEMENT_READY",
            _NMF + "/source_inputs/nmf_assignments_k4_pre_k6restore.csv",
            _NMF + "/source_inputs/nmf_assignments_k6.csv",
            _NMF + "/source_inputs/program_labels_k4_pre_k6restore.csv",
            _NMF + "/source_inputs/program_labels_k6.csv",
            _NMF + "/source_inputs/per_seed_metric_summary.csv",
            _NMF + "/source_inputs/cross_seed_stability_summary.csv",
        ),
        "required_producers": (
            *_COMMON_PRODUCERS,
            NMF_FREEZER_PRODUCER_PATH,
        ),
    },
    "genetics_evidence": {
        "output_relative_path": "genetics_evidence.tsv",
        "required_inputs": (
            _PROGRAM_CANDIDATE + "/genetics_context/phenotype_registry.tsv",
            _PROGRAM_CANDIDATE + "/genetics_context/power_stratified_interface.tsv",
            _PROGRAM_CANDIDATE + "/genetics_context/terminal_closure.tsv",
        ),
        "required_producers": _COMMON_PRODUCERS,
    },
    "spatial_evidence": {
        "output_relative_path": "spatial_evidence.tsv",
        "required_inputs": (
            _PROGRAM_CANDIDATE
            + "/spatial_context_semantic_v2_2026-08-08/final_integration/"
            "figure4_program_matrix.tsv",
        ),
        "required_producers": _COMMON_PRODUCERS,
    },
    "passport_evidence": {
        "output_relative_path": "passport_evidence.tsv",
        "required_inputs": (
            _PASSPORT + "/passport_gene_index.parquet",
            _PASSPORT + "/passport_evidence_long.parquet",
            _PASSPORT + "/passport_coverage_long.parquet",
            _PASSPORT + "/passport_program_context.parquet",
            _PASSPORT + "/passport_next_experiment.tsv",
            _PASSPORT + "/passport_source_nodes.tsv",
            _PASSPORT + "/passport_source_edges.tsv",
        ),
        "required_producers": _COMMON_PRODUCERS,
    },
    "myojin_supplement": {
        "output_relative_path": "myojin_supplement.tsv",
        "required_inputs": (
            _PROGRAM_CANDIDATE + "/myojin_hlf/class_effects.tsv",
            _PROGRAM_CANDIDATE + "/myojin_hlf/program_effects.tsv",
        ),
        "required_producers": _COMMON_PRODUCERS,
    },
}

EXPECTED_ADAPTER_IDS = frozenset(ADAPTER_CONTRACT)


def _contract_payload() -> dict[str, object]:
    return {
        "contract_version": ADAPTER_PROVENANCE_VERSION,
        "adapters": {
            artifact_id: {
                "output_relative_path": str(spec["output_relative_path"]),
                "required_inputs": list(spec["required_inputs"]),
                "required_producers": list(spec["required_producers"]),
            }
            for artifact_id, spec in sorted(ADAPTER_CONTRACT.items())
        },
    }


ADAPTER_CONTRACT_SHA256 = hashlib.sha256(
    json.dumps(_contract_payload(), sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
).hexdigest()


def adapter_input_paths(project_root: Path, artifact_id: str) -> tuple[Path, ...]:
    """Resolve the exact frozen input set for one adapter below a project root."""

    if artifact_id not in ADAPTER_CONTRACT:
        raise ReleaseContractError(f"unknown adapter contract ID: {artifact_id}")
    project = project_root.resolve()
    return tuple(
        project / str(relative)
        for relative in ADAPTER_CONTRACT[artifact_id]["required_inputs"]
    )


def require_exact_adapter_binding_paths(
    artifact_id: str,
    observed_inputs: Iterable[str],
    observed_producers: Iterable[str],
) -> None:
    """Reject any omitted, extra, or duplicate input/producer binding."""

    if artifact_id not in ADAPTER_CONTRACT:
        raise ReleaseContractError(f"unknown adapter contract ID: {artifact_id}")
    inputs = tuple(str(path) for path in observed_inputs)
    producers = tuple(str(path) for path in observed_producers)
    if len(inputs) != len(set(inputs)):
        raise ReleaseContractError(
            f"adapter contract has duplicate input bindings: {artifact_id}"
        )
    if len(producers) != len(set(producers)):
        raise ReleaseContractError(
            f"adapter contract has duplicate producer bindings: {artifact_id}"
        )
    expected_inputs = set(ADAPTER_CONTRACT[artifact_id]["required_inputs"])
    expected_producers = set(ADAPTER_CONTRACT[artifact_id]["required_producers"])
    observed_input_set = set(inputs)
    observed_producer_set = set(producers)
    if observed_input_set != expected_inputs:
        raise ReleaseContractError(
            "adapter input set differs from frozen contract: "
            f"{artifact_id}; missing={sorted(expected_inputs - observed_input_set)}; "
            f"extra={sorted(observed_input_set - expected_inputs)}"
        )
    if observed_producer_set != expected_producers:
        raise ReleaseContractError(
            "adapter producer set differs from frozen contract: "
            f"{artifact_id}; missing={sorted(expected_producers - observed_producer_set)}; "
            f"extra={sorted(observed_producer_set - expected_producers)}"
        )
