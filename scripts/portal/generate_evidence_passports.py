#!/usr/bin/env python3
"""MASLD Gene Catalog schema, rulebook, signed-selection gate, and fixtures.

Legacy ``passport_*`` identifiers remain part of the sealed v1 file contract.

PASS-01 owns the shared contract. Real PASS-02--PASS-06 assembly is delegated
to ``build_evidence_passport_bundle.py`` and remains impossible without a
coordinator-attested selection plus terminal GEN and Plan 13 gates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd
import pyarrow as pa


SCRIPT_PATH = Path(__file__).resolve()
LOGICAL_PRODUCER_ID = "scripts/portal/generate_evidence_passports.py"
PROJECT_ROOT = Path(
    os.environ.get("MASLD_PROJECT_ROOT", os.fspath(SCRIPT_PATH.parents[2]))
).resolve()
# PASSPORT_CANDIDATE_ROOT and PASSPORT_ANALYSIS_RELEASE_ID move every production
# step (adjudication, selection, build, validation) to a new candidate root and
# release together; unset, they are the 2026-08-07 v1 candidate.
CANDIDATE_REL = Path(os.environ.get(
    "PASSPORT_CANDIDATE_ROOT",
    "RNA-seq/results/evidence_passports/candidates/"
    "program-context-v2-candidate-2026-08-07",
))
DEFAULT_CANDIDATE_ROOT = PROJECT_ROOT / CANDIDATE_REL
DEFAULT_FIXTURE_ROOT = DEFAULT_CANDIDATE_ROOT / "fixtures"
ANALYSIS_RELEASE_ID = "fixture-passport-pass01-v1"
PRODUCTION_ANALYSIS_RELEASE_ID = os.environ.get(
    "PASSPORT_ANALYSIS_RELEASE_ID", "program-context-v2-candidate-2026-08-07-passports-v1"
)
FIXTURE_SIGNED_AT = "2026-08-07T00:00:00Z"
DISCLAIMER = (
    "This experiment is a discriminating research proposal, not a calibrated "
    "probability of therapeutic success or a clinical assay recommendation."
)


class PassportContractError(RuntimeError):
    """Fail-closed contract error with a stable machine-readable code."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


CALL_STATES = OrderedDict(
    [
        ("supported", "Frozen result passed its prespecified positive and robustness gate."),
        ("discordant", "Frozen upstream comparison explicitly reported conflicting accepted results."),
        ("tested_negative", "Adequately testable prespecified assay was explicitly called negative upstream."),
        ("indeterminate", "A test ran but uncertainty or design limitations precluded either directional call."),
        ("untestable", "Coverage, power, detectability, mapping, or eligibility failed."),
        ("not_applicable", "The assay cannot answer this gene/context question by design."),
    ]
)

TESTABILITY_STATES = OrderedDict(
    [
        ("testable", "The originating workstream's complete testability gate passed."),
        ("insufficient_gene_coverage", "Too few mapped genes or probes were observed."),
        ("insufficient_retained_weight", "Too little prespecified program weight was retained."),
        ("not_detected", "The feature was below the assay's prespecified detectability gate."),
        ("underpowered_source", "The source lacked adequate power for a positive or negative interpretation."),
        ("mapping_ambiguous", "Stable identity or allele mapping remained ambiguous."),
        ("insufficient_shared_posterior", "No SuSiE-COLOC signal pair kept enough posterior mass on shared SNPs to be tested."),
        ("biological_unit_unresolved", "Independent biological units could not be resolved."),
        ("assay_out_of_scope", "The assay is not designed to answer the stated question."),
        ("source_gate_failed", "A dataset-level source gate failed; no gene-level pseudo-result is authorized."),
    ]
)

PROVENANCE_STATES = OrderedDict(
    [
        ("independent", "No samples or outcome-derived construction features overlap."),
        ("partially_dependent", "Some cohorts, samples, or construction features overlap."),
        ("reused_source", "The same source or an outcome-informed derivative was reused directly."),
    ]
)

ROLE_HYPOTHESES = OrderedDict(
    [
        ("inherited_regulatory_anchor", "Hypothesis: inherited regulatory evidence anchors the candidate."),
        ("context_regulatory_candidate", "Hypothesis: context-specific regulation is the discriminating layer."),
        ("established_state_marker", "Hypothesis: the gene primarily marks established disease state."),
        ("downstream_effector_candidate", "Hypothesis: perturbation could distinguish an effector from a marker."),
        ("multi_role", "Hypothesis: more than one prespecified evidence role is represented."),
        ("unresolved", "Evidence topology does not support a more specific role hypothesis."),
    ]
)

PRIMARY_EVIDENCE_CLASSES = OrderedDict(
    [
        ("genetically_anchored", "Accepted static genetic evidence."),
        ("established_state_associated", "Accepted established-state association."),
        ("context_supported", "Accepted source-faithful context annotation or bridge."),
        ("concordant", "Accepted genetic and state evidence are concordant."),
        ("discordant", "Accepted upstream comparison explicitly classified discordance."),
        ("untested", "The decisive layer was not tested."),
        ("unresolved", "No stronger class is authorized."),
    ]
)

DIRECTIONS = OrderedDict(
    [
        ("positive", "Positive native-assay effect direction."),
        ("negative", "Negative native-assay effect direction."),
        ("discordant", "Explicitly discordant upstream direction."),
        ("no_effect", "Source workstream explicitly classified no effect in an adequate test."),
        ("not_directional", "The statistic does not encode a biological direction."),
        ("unknown", "Direction could not be resolved without semantic inflation."),
    ]
)

EVIDENCE_DOMAINS = OrderedDict(
    [
        ("genetics", "GWAS/eQTL regulatory evidence."),
        ("transcriptomics", "Bulk or single-cell transcript abundance/program evidence."),
        ("spatial", "Spatial localization or exposure evidence."),
        ("proteomics", "Protein-abundance evidence."),
        ("chromatin", "Chromatin or regulatory-element evidence."),
        ("functional", "Perturbational or functional-challenge evidence."),
        ("cross_species", "Cross-species support with explicit species boundary."),
        ("cell_context", "Cell-type or state-resolved contextual evidence."),
    ]
)

NODE_TYPES = OrderedDict(
    [
        ("publication", "Source publication."),
        ("cohort", "Biological cohort or participant pool."),
        ("dataset", "Deposited assay dataset."),
        ("frozen_release", "Immutable analysis release."),
        ("derived_object", "Accepted result derived from source nodes."),
    ]
)
EDGE_TYPES = OrderedDict(
    [
        ("derived_from", "Directed derivation relationship."),
        ("shares_samples_with", "Explicit cohort/sample overlap relationship."),
        ("tested_by", "Result was evaluated by the target assay or dataset."),
    ]
)
INCLUSION_DESTINATIONS = OrderedDict(
    [
        ("main", "Eligible for a main-text candidate view."),
        ("supplement", "Eligible for supplementary presentation."),
        ("passport_only", "Eligible only for provenance-preserving Gene Catalog display."),
        ("exclude", "Retained in audit but excluded from gene-level evidence."),
    ]
)
GATE_VERDICTS = OrderedDict(
    [
        ("include", "Accepted positive or descriptive artifact."),
        (
            "complete_nonconfirmatory",
            "Complete prespecified analysis that did not authorize a positive or informative-negative call.",
        ),
        ("skipped_source_gate", "Dataset-level source gate failed; no gene-level pseudo-results."),
    ]
)

VOCABULARIES: OrderedDict[str, OrderedDict[str, str]] = OrderedDict(
    [
        ("call_state", CALL_STATES),
        ("testability_state", TESTABILITY_STATES),
        ("provenance_state", PROVENANCE_STATES),
        ("role_hypothesis", ROLE_HYPOTHESES),
        ("primary_evidence_class", PRIMARY_EVIDENCE_CLASSES),
        ("direction", DIRECTIONS),
        ("evidence_domain", EVIDENCE_DOMAINS),
        ("source_node_type", NODE_TYPES),
        ("source_edge_type", EDGE_TYPES),
        ("inclusion_destination", INCLUSION_DESTINATIONS),
        ("gate_verdict", GATE_VERDICTS),
    ]
)


PARQUET_SCHEMAS: dict[str, OrderedDict[str, pa.DataType]] = {
    "passport_gene_index.parquet": OrderedDict(
        [
            ("analysis_release_id", pa.string()),
            ("passport_id", pa.string()),
            ("ensembl_id", pa.string()),
            ("symbol", pa.string()),
            ("symbol_collision", pa.bool_()),
            ("primary_evidence_class", pa.string()),
            ("role_hypothesis", pa.string()),
            ("role_rule_id", pa.string()),
            ("next_experiment_rule_id", pa.string()),
            ("allowed_wording", pa.string()),
            ("limitation", pa.string()),
        ]
    ),
    "passport_evidence_long.parquet": OrderedDict(
        [
            ("analysis_release_id", pa.string()),
            ("evidence_result_id", pa.string()),
            ("passport_id", pa.string()),
            ("ensembl_id", pa.string()),
            ("symbol", pa.string()),
            ("evidence_domain", pa.string()),
            ("assay", pa.string()),
            ("dataset_id", pa.string()),
            ("source_node_id", pa.string()),
            ("source_release_id", pa.string()),
            ("phenotype", pa.string()),
            ("context", pa.string()),
            ("biological_unit", pa.string()),
            ("contrast_or_exposure", pa.string()),
            ("effect_unit", pa.string()),
            ("estimate", pa.float64()),
            ("standard_error", pa.float64()),
            ("ci_lower", pa.float64()),
            ("ci_upper", pa.float64()),
            ("p_value", pa.float64()),
            ("q_value", pa.float64()),
            ("direction", pa.string()),
            ("testability_state", pa.string()),
            ("testability_reason", pa.string()),
            ("call_state", pa.string()),
            ("negative_call_rule_id", pa.string()),
            ("negative_decision_boundary", pa.float64()),
            ("negative_margin", pa.float64()),
            ("negative_call_passed", pa.bool_()),
            ("gate_id", pa.string()),
            ("n_biological_units", pa.int64()),
            ("n_technical_units", pa.int64()),
            ("provenance_state", pa.string()),
            ("source_dependent", pa.bool_()),
            ("source_artifact_sha256", pa.string()),
            ("source_input_row_id", pa.string()),
            ("source_call_state", pa.string()),
            ("source_testability_state", pa.string()),
            ("source_provenance_state", pa.string()),
            ("source_negative_call_rule_id", pa.string()),
            ("source_negative_decision_boundary", pa.float64()),
            ("source_negative_margin", pa.float64()),
            ("source_negative_call_passed", pa.bool_()),
            ("source_row_sha256", pa.string()),
            ("allowed_wording", pa.string()),
            ("limitation", pa.string()),
        ]
    ),
    "passport_coverage_long.parquet": OrderedDict(
        [
            ("analysis_release_id", pa.string()),
            ("coverage_result_id", pa.string()),
            ("passport_id", pa.string()),
            ("ensembl_id", pa.string()),
            ("symbol", pa.string()),
            ("evidence_domain", pa.string()),
            ("assay", pa.string()),
            ("dataset_id", pa.string()),
            ("source_release_id", pa.string()),
            ("testability_state", pa.string()),
            ("testability_reason", pa.string()),
            ("call_state", pa.string()),
            ("n_biological_units", pa.int64()),
            ("n_technical_units", pa.int64()),
            ("coverage_denominator", pa.string()),
            ("limitation", pa.string()),
        ]
    ),
}


PRODUCTION_ONLY_PARQUET_SCHEMAS: dict[str, OrderedDict[str, pa.DataType]] = {
    "passport_program_index.parquet": OrderedDict(
        [
            ("analysis_release_id", pa.string()),
            ("program_uid", pa.string()),
            ("program_label", pa.string()),
            ("cell_type", pa.string()),
            ("membership_sha256", pa.string()),
            ("source_release_id", pa.string()),
            ("program_call_scope", pa.string()),
            ("member_gene_call_expansion_authorized", pa.bool_()),
        ]
    ),
    "passport_program_membership.parquet": OrderedDict(
        [
            ("analysis_release_id", pa.string()),
            ("program_uid", pa.string()),
            ("ensembl_id", pa.string()),
            ("symbol", pa.string()),
            ("source_weight", pa.float64()),
            ("original_l1_weight", pa.float64()),
            ("membership_sha256", pa.string()),
            ("membership_only", pa.bool_()),
            ("gene_call_expansion_authorized", pa.bool_()),
            ("source_input_row_id", pa.string()),
            ("source_artifact_sha256", pa.string()),
            ("source_row_sha256", pa.string()),
        ]
    ),
    "passport_program_context.parquet": OrderedDict(
        [
            ("analysis_release_id", pa.string()),
            ("program_context_id", pa.string()),
            ("program_uid", pa.string()),
            ("program_label", pa.string()),
            ("cell_type", pa.string()),
            ("dataset_id", pa.string()),
            ("assay", pa.string()),
            ("source_node_id", pa.string()),
            ("source_release_id", pa.string()),
            ("context", pa.string()),
            ("biological_unit", pa.string()),
            ("contrast_or_exposure", pa.string()),
            ("effect_unit", pa.string()),
            ("estimate", pa.float64()),
            ("standard_error", pa.float64()),
            ("ci_lower", pa.float64()),
            ("ci_upper", pa.float64()),
            ("p_value", pa.float64()),
            ("q_value", pa.float64()),
            ("direction", pa.string()),
            ("testability_state", pa.string()),
            ("testability_reason", pa.string()),
            ("call_state", pa.string()),
            ("negative_call_rule_id", pa.string()),
            ("negative_decision_boundary", pa.float64()),
            ("negative_margin", pa.float64()),
            ("negative_call_passed", pa.bool_()),
            ("gate_id", pa.string()),
            ("n_biological_units", pa.int64()),
            ("n_technical_units", pa.int64()),
            ("provenance_state", pa.string()),
            ("source_dependent", pa.bool_()),
            ("source_artifact_sha256", pa.string()),
            ("source_input_row_id", pa.string()),
            ("source_row_sha256", pa.string()),
            ("allowed_wording", pa.string()),
            ("limitation", pa.string()),
            ("gene_call_expansion_authorized", pa.bool_()),
        ]
    ),
}

PRODUCTION_PARQUET_SCHEMAS = {
    **PARQUET_SCHEMAS,
    **PRODUCTION_ONLY_PARQUET_SCHEMAS,
}

PORTAL_PARQUET_SCHEMAS: dict[str, OrderedDict[str, pa.DataType]] = {
    "portal_export/passport_summary.parquet": OrderedDict(
        [
            *PARQUET_SCHEMAS["passport_gene_index.parquet"].items(),
            ("alphabetical_key", pa.string()),
        ]
    ),
    "portal_export/passport_evidence.parquet": PARQUET_SCHEMAS[
        "passport_evidence_long.parquet"
    ],
    "portal_export/passport_program_index.parquet": PRODUCTION_ONLY_PARQUET_SCHEMAS[
        "passport_program_index.parquet"
    ],
    "portal_export/passport_program_membership.parquet": PRODUCTION_ONLY_PARQUET_SCHEMAS[
        "passport_program_membership.parquet"
    ],
    "portal_export/passport_program_context.parquet": PRODUCTION_ONLY_PARQUET_SCHEMAS[
        "passport_program_context.parquet"
    ],
}


TSV_SCHEMAS: dict[str, list[str]] = {
    "passport_source_nodes.tsv": [
        "source_node_id", "node_type", "label", "source_release_id", "uri",
        "artifact_sha256", "access_status", "biological_unit", "source_datasets",
    ],
    "passport_source_edges.tsv": [
        "source_edge_id", "from_node_id", "to_node_id", "edge_type",
        "evidence_result_id", "note",
    ],
    "passport_next_experiment.tsv": [
        "analysis_release_id", "passport_id", "ensembl_id", "symbol",
        "experiment_rule_id", "role_rule_id", "biological_model", "context",
        "perturbation", "primary_readout", "falsifying_outcome", "disclaimer",
    ],
    "passport_rulebook.tsv": [
        "rule_id", "rule_type", "rule_version", "condition_json", "output_value",
        "biological_model", "context", "perturbation", "primary_readout",
        "falsifying_outcome", "prohibited_use",
    ],
    "passport_controlled_vocabularies.tsv": [
        "vocabulary_name", "value", "definition", "permitted_source",
        "prohibited_interpretation",
    ],
    "passport_data_dictionary.tsv": [
        "table_name", "column_name", "storage_type", "required",
        "allowed_values_vocabulary", "unit_semantics", "null_semantics",
        "ui_label", "description",
    ],
    "passport_gate_status.tsv": ["gate_id", "status", "promotion_allowed", "reason"],
    "passport_release_manifest.tsv": [
        "relative_path", "sha256", "bytes", "producer", "environment",
        "producer_sha256", "upstream_artifact_sha256", "release_status",
    ],
}

PRODUCTION_ONLY_TSV_SCHEMAS: dict[str, list[str]] = {
    "accepted_input_audit.tsv": [
        "input_id", "workstream_id", "artifact_role", "adapter_id",
        "artifact_grain", "artifact_path", "artifact_sha256", "artifact_bytes",
        "source_release_id", "gate_verdict", "terminal_gate_input_id",
        "provenance_state", "claim_scope", "inclusion_destination",
        "selection_sha256", "validation_status", "validation_reason",
    ],
    "passport_adapter_audit.tsv": [
        "input_id", "adapter_id", "artifact_role", "artifact_grain",
        "source_rows", "gene_rows", "evidence_rows", "coverage_rows",
        "program_rows", "membership_rows", "dataset_status_rows",
        "assay_status_rows", "quarantine_rows", "status", "reason",
    ],
    "passport_identity_quarantine.tsv": [
        "quarantine_id", "input_id", "adapter_id", "artifact_role",
        "source_input_row_id", "source_symbol", "source_ensembl_text",
        "normalized_ensembl_candidates", "reason_code", "reason",
        "source_artifact_sha256", "source_row_sha256",
    ],
    "passport_gen_identity_status.tsv": [
        "release_id", "status", "source_artifact_sha256",
        "identity_artifact_sha256", "adjudicated_artifact_sha256",
        "quarantine_artifact_sha256", "audit_artifact_sha256",
        "producer_sha256", "n_source_rows", "n_adjudicated_genes",
        "n_quarantined_source_rows", "n_irreducible_multi_id",
        "n_absent_v49", "n_identity_conflict", "n_unresolved",
        "n_invalid_token", "outcome_fields_used", "merge_discordant_calls",
        "one_row_per_ensembl", "canonical_promotion_authorized",
        "fixture_only",
    ],
    "passport_domain_coverage.tsv": [
        "analysis_release_id", "evidence_domain", "gene_evidence_status",
        "program_context_status", "assay_status", "gene_evidence_authorized",
        "member_gene_expansion_authorized", "n_gene_evidence_rows",
        "n_programs", "n_program_membership_rows", "n_program_context_rows",
        "n_assay_status_rows", "source_adapter_ids", "coverage_boundary",
        "allowed_wording", "limitation",
    ],
    "passport_dataset_status.tsv": [
        "dataset_status_id", "dataset_id", "assay", "source_node_id",
        "source_release_id", "dataset_gate", "status", "reason",
        "biological_unit", "n_biological_units", "n_technical_units",
        "provenance_state", "source_dependent", "source_artifact_sha256",
        "source_input_row_id", "source_row_sha256", "gene_expansion_authorized",
        "allowed_wording", "limitation",
    ],
    "passport_assay_status.tsv": [
        "assay_status_id", "result_grain", "entity_id", "assay", "dataset_id",
        "source_node_id", "source_release_id", "status", "call_state",
        "testability_state", "estimate", "standard_error", "ci_lower",
        "ci_upper", "p_value", "q_value", "effect_unit", "provenance_state",
        "source_dependent", "source_artifact_sha256", "source_input_row_id",
        "source_row_sha256", "gene_expansion_authorized", "allowed_wording",
        "limitation",
    ],
    "passport_build_status.tsv": [
        "stage_id", "status", "promotion_allowed", "selection_sha256", "reason",
    ],
    "passport_validation_report.tsv": [
        "check_id", "status", "detail",
    ],
    "passport_plan60_handoff.tsv": [
        "analysis_release_id", "bundle_path", "manifest_sha256",
        "selection_sha256", "automated_validation", "manual_acceptance",
        "promotion_allowed", "build_command", "environment",
        "scientific_call_recomputed", "omitted_input_inventory",
    ],
}

PRODUCTION_TSV_SCHEMAS = {
    **TSV_SCHEMAS,
    **PRODUCTION_ONLY_TSV_SCHEMAS,
}

PASS06_TERMINAL_COLUMNS = [
    "analysis_release_id", "status", "selection_sha256", "manifest_sha256",
    "validation_report_sha256", "n_genes", "n_evidence_rows", "n_programs",
    "n_program_context_rows", "automated_validation", "manual_acceptance",
    "handoff_allowed", "canonical_promotion_authorized",
    "scientific_call_recomputed", "validated_at_utc",
]
PASS06_TERMINAL_PROVENANCE_COLUMNS = [
    "analysis_release_id", "bundle_uri", "artifact_role", "relative_path",
    "sha256", "bytes", "producer", "producer_sha256",
]
MANUAL_ACCEPTANCE_COLUMNS = [
    "review_id", "analysis_release_id", "selection_sha256", "reviewer",
    "reviewed_at_utc", "ui_index_sha256", "ui_contract_sha256",
    "ui_source_manifest_sha256", "ui_review_manifest_sha256",
    "call_states_reviewed", "provenance_states_reviewed",
    "hero_genes_reviewed", "visible_boundary_pass", "visible_testability_pass",
    "visible_source_dependence_pass", "visible_falsifier_pass", "decision",
]
MANUAL_ACCEPTANCE_ATTESTATION_VERSION = "passport_manual_acceptance_v2"


ROLE_ROUTES = OrderedDict(
    [
        ("genetically_anchored", ("inherited_regulatory_anchor", "ROLE_INHERITED_V1", "EXP_ALLELE_AWARE_V1")),
        ("context_supported", ("context_regulatory_candidate", "ROLE_CONTEXT_V1", "EXP_CONTEXT_PERTURB_V1")),
        ("established_state_associated", ("established_state_marker", "ROLE_STATE_V1", "EXP_STRESS_PERTURB_V1")),
        ("concordant", ("multi_role", "ROLE_MULTI_V1", "EXP_MULTI_LAYER_V1")),
        ("discordant", ("unresolved", "ROLE_UNRESOLVED_V1", "EXP_MATCHED_RNA_PROTEIN_V1")),
        ("untested", ("unresolved", "ROLE_UNRESOLVED_V1", "EXP_DIRECT_ASSAY_V1")),
        ("unresolved", ("unresolved", "ROLE_UNRESOLVED_V1", "EXP_UNRESOLVED_V1")),
    ]
)

EXPERIMENT_RULES = {
    "EXP_ALLELE_AWARE_V1": {
        "biological_model": "isogenic human hepatocyte model",
        "context": "basal and lipid-stress conditions",
        "perturbation": "allele-aware regulatory editing",
        "primary_readout": "target transcript plus locus-linked molecular phenotype",
        "falsifying_outcome": "the edited risk allele does not change the target in either context",
    },
    "EXP_CONTEXT_PERTURB_V1": {
        "biological_model": "matched primary-like liver cell model",
        "context": "source-matched disease or regulatory state",
        "perturbation": "context-specific enhancer or gene perturbation",
        "primary_readout": "context-stratified target expression and state phenotype",
        "falsifying_outcome": "the proposed context does not modify the regulatory or phenotypic effect",
    },
    "EXP_STRESS_PERTURB_V1": {
        "biological_model": "human liver cell model with relevant stress response",
        "context": "lipotoxic or inflammatory challenge",
        "perturbation": "gene knockdown and rescue",
        "primary_readout": "cell-state program and stress phenotype",
        "falsifying_outcome": "perturbation changes neither the state program nor the stress phenotype",
    },
    "EXP_MULTI_LAYER_V1": {
        "biological_model": "isogenic liver model with orthogonal molecular assays",
        "context": "basal and disease-relevant challenge",
        "perturbation": "allele-aware edit plus gene perturbation",
        "primary_readout": "regulation, transcript, protein, and phenotype measured in the same model",
        "falsifying_outcome": "the inherited and established-state signals cannot be connected in one model",
    },
    "EXP_MATCHED_RNA_PROTEIN_V1": {
        "biological_model": "matched donor or isogenic liver model",
        "context": "the context in which discordance was observed",
        "perturbation": "targeted perturbation with paired RNA and protein sampling",
        "primary_readout": "matched transcript and protein abundance",
        "falsifying_outcome": "the reported discordance is not reproduced under paired measurement",
    },
    "EXP_DIRECT_ASSAY_V1": {
        "biological_model": "source-appropriate human tissue or cell model",
        "context": "the unresolved biological context",
        "perturbation": "acquire the missing decisive assay without outcome-informed selection",
        "primary_readout": "prespecified detectability and effect statistic",
        "falsifying_outcome": "adequate direct measurement does not support the proposed evidence layer",
    },
    "EXP_UNRESOLVED_V1": {
        "biological_model": "source-appropriate human model",
        "context": "prespecified unresolved context",
        "perturbation": "measure the single most decisive missing layer",
        "primary_readout": "prespecified assay-native endpoint",
        "falsifying_outcome": "the decisive layer remains unsupported under adequate measurement",
    },
}


SELECTION_REQUIRED_COLUMNS = [
    "input_id",
    "workstream_id",
    "artifact_role",
    "adapter_id",
    "artifact_grain",
    "artifact_path",
    "artifact_sha256",
    "artifact_bytes",
    "producer_script",
    "producer_sha256",
    "source_release_id",
    "gate_verdict",
    "terminal_gate_input_id",
    "allowed_claim_wording",
    "limitation",
    "source_datasets",
    "source_cohorts",
    "source_publication",
    "source_url",
    "biological_unit",
    "n_biological_units",
    "n_technical_units",
    "provenance_state",
    "claim_scope",
    "inclusion_destination",
    "selection_decision",
    "parent_input_ids",
    "fixture_only",
]

SELECTION_ATTESTATION_VERSION = "passport_selection_attestation_v2"
SELECTION_ARTIFACT_GRAINS = {
    "reference", "gate", "gene", "program", "program_membership",
    "dataset", "assay", "provenance",
}
SELECTION_CLAIM_SCOPES = {
    "identity_only", "gate_only", "gene_evidence", "program_context",
    "program_membership_only", "dataset_status_only", "assay_status_only",
    "source_provenance_only",
}
SELECTION_ADAPTER_IDS = {
    "synthetic_gene_evidence_v2",
    "gencode_v49_identity_v1",
    "gen_terminal_closure_v1",
    "gen_terminal_ready_v1",
    "gen_frozen_classes_v1",
    "gen_raw_identity_source_v1",
    "gen_identity_adjudication_ready_v1",
    "gen_identity_adjudication_audit_v1",
    "gen_identity_adjudication_quarantine_v1",
    "gen_adjudicated_classes_v1",
    "accepted_gene_evidence_v1",
    "hotspot_semantic_ready_v1",
    "hotspot_program_semantics_v1",
    "hotspot_program_membership_v1",
    "plan13_terminal_ready_v1",
    "plan13_final_ready_v1",
    "plan13_program_context_v1",
    "dataset_status_v1",
    "myojin_terminal_ready_v1",
    "myojin_assay_verdict_v1",
    "myojin_class_nonsupport_v1",
    "myojin_program_nonsupport_v1",
    "source_provenance_v1",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def stable_sha256(payload: Mapping[str, Any]) -> str:
    data = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(data.encode("utf-8")).hexdigest()


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def assert_candidate_root(candidate_root: Path) -> Path:
    expected = DEFAULT_CANDIDATE_ROOT.resolve()
    observed = candidate_root.resolve()
    if observed != expected:
        raise PassportContractError("OUTPUT_ROOT", f"expected {expected}; received {observed}")
    return observed


def assert_fixture_path(path: Path, candidate_root: Path) -> Path:
    resolved = path.resolve()
    fixture_root = (candidate_root / "fixtures").resolve()
    if not _is_relative_to(resolved, fixture_root):
        raise PassportContractError("FIXTURE_OUTPUT_ROOT", f"fixture path escapes {fixture_root}: {resolved}")
    return resolved


def atomic_write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    ) as handle:
        tmp = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PassportContractError("EMPTY_TABLE", str(path))
        return list(reader)


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    ) as handle:
        tmp = Path(handle.name)
        json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")
    os.replace(tmp, path)


def _pandas_dtype(arrow_type: pa.DataType) -> str:
    if pa.types.is_string(arrow_type):
        return "string"
    if pa.types.is_boolean(arrow_type):
        return "boolean"
    if pa.types.is_integer(arrow_type):
        return "Int64"
    if pa.types.is_floating(arrow_type):
        return "Float64"
    raise PassportContractError("UNSUPPORTED_TYPE", str(arrow_type))


def coerce_dataframe(rows: list[dict[str, Any]], schema: OrderedDict[str, pa.DataType]) -> pd.DataFrame:
    frame = pd.DataFrame(rows, columns=list(schema))
    for column, arrow_type in schema.items():
        frame[column] = frame[column].astype(_pandas_dtype(arrow_type))
    return frame


def atomic_write_parquet(path: Path, frame: pd.DataFrame, schema: OrderedDict[str, pa.DataType]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = frame[list(schema)].copy()
    sort_columns = [
        column for column in ("passport_id", "evidence_result_id", "coverage_result_id")
        if column in ordered.columns
    ]
    if sort_columns:
        ordered = ordered.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp",
        delete=False,
    ) as handle:
        tmp = Path(handle.name)
    try:
        ordered.to_parquet(tmp, engine="pyarrow", index=False, compression="zstd")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def validate_signed_selection(selection_path: Path, allow_fixture: bool = False) -> list[dict[str, str]]:
    if not selection_path.is_file():
        raise PassportContractError(
            "SIGNED_SELECTION_MISSING",
            f"production requires coordinator-selected {selection_path}",
        )
    signature_path = selection_path.with_name("passport_input_selection.signature.json")
    if not signature_path.is_file():
        raise PassportContractError("SELECTION_ATTESTATION_MISSING", str(signature_path))
    with selection_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PassportContractError("SELECTION_HEADER", "selection is empty")
        missing = sorted(set(SELECTION_REQUIRED_COLUMNS) - set(reader.fieldnames))
        if missing:
            raise PassportContractError("SELECTION_SCHEMA", f"missing columns: {missing}")
        rows = list(reader)
    if not rows:
        raise PassportContractError("SELECTION_EMPTY", "no selected artifacts")
    signature = json.loads(signature_path.read_text(encoding="utf-8"))
    required_signature = {
        "attestation_version", "selection_sha256", "signed_by", "signed_at_utc",
        "decision_register_id", "authority_document", "fixture_only",
        "coordinator_attested", "analysis_release_id", "gen_terminal_closure_input_id",
        "gen_terminal_ready_input_id", "plan13_terminal_ready_input_id",
    }
    if not required_signature.issubset(signature):
        raise PassportContractError(
            "SELECTION_ATTESTATION_SCHEMA",
            f"missing keys: {sorted(required_signature - set(signature))}",
        )
    if signature["selection_sha256"] != sha256_file(selection_path):
        raise PassportContractError("SELECTION_ATTESTATION_HASH", "selection checksum mismatch")
    if signature["attestation_version"] != SELECTION_ATTESTATION_VERSION:
        raise PassportContractError(
            "SELECTION_ATTESTATION_VERSION",
            f"expected {SELECTION_ATTESTATION_VERSION}; observed {signature['attestation_version']}",
        )
    for key in ("fixture_only", "coordinator_attested"):
        if type(signature[key]) is not bool:
            raise PassportContractError(
                "SELECTION_ATTESTATION_BOOLEAN", f"{key} must be a JSON Boolean"
            )
    if not str(signature["signed_by"]).strip() or not str(signature["decision_register_id"]).strip():
        raise PassportContractError("SELECTION_ATTESTATION_IDENTITY", "signer and decision ID are required")
    if not str(signature["analysis_release_id"]).strip():
        raise PassportContractError("SELECTION_RELEASE_ID", "analysis_release_id is required")
    fixture_rows = {row["fixture_only"].strip().lower() for row in rows}
    if fixture_rows - {"true", "false"} or len(fixture_rows) != 1:
        raise PassportContractError("SELECTION_FIXTURE_FLAG", "fixture_only must be one consistent Boolean")
    if (fixture_rows == {"true"}) != signature["fixture_only"]:
        raise PassportContractError(
            "SELECTION_FIXTURE_ATTESTATION_MISMATCH",
            "row fixture_only flags must equal the signed JSON fixture_only flag",
        )
    is_fixture = signature["fixture_only"]
    if is_fixture and not allow_fixture:
        raise PassportContractError("FIXTURE_SELECTION_PRODUCTION", "fixture attestation cannot unlock production")
    if not is_fixture and not signature["coordinator_attested"]:
        raise PassportContractError(
            "UNSIGNED_COORDINATOR_SELECTION",
            "production requires coordinator_attested=true in the signed attestation",
        )
    input_ids = [row["input_id"].strip() for row in rows]
    if any(not value for value in input_ids) or len(input_ids) != len(set(input_ids)):
        raise PassportContractError("SELECTION_INPUT_ID", "input_id must be nonempty and unique")
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", value) for value in input_ids):
        raise PassportContractError(
            "SELECTION_INPUT_ID",
            "input_id may contain only letters, digits, dot, underscore, and hyphen",
        )
    input_id_set = set(input_ids)
    for row in rows:
        if row["selection_decision"] != "accepted":
            raise PassportContractError("SELECTION_DECISION", "only coordinator-accepted rows may enter")
        if row["adapter_id"] not in SELECTION_ADAPTER_IDS:
            raise PassportContractError("SELECTION_ADAPTER", row["adapter_id"])
        if row["artifact_grain"] not in SELECTION_ARTIFACT_GRAINS:
            raise PassportContractError("SELECTION_ARTIFACT_GRAIN", row["artifact_grain"])
        if row["claim_scope"] not in SELECTION_CLAIM_SCOPES:
            raise PassportContractError("SELECTION_CLAIM_SCOPE", row["claim_scope"])
        if row["gate_verdict"] not in GATE_VERDICTS:
            raise PassportContractError("SELECTION_GATE_VERDICT", row["gate_verdict"])
        if row["inclusion_destination"] not in INCLUSION_DESTINATIONS:
            raise PassportContractError("SELECTION_DESTINATION", row["inclusion_destination"])
        if row["provenance_state"] not in PROVENANCE_STATES:
            raise PassportContractError("SELECTION_PROVENANCE", row["provenance_state"])
        for count_field in ("n_biological_units", "n_technical_units"):
            value = row[count_field].strip()
            if value and (not value.isdigit() or int(value) < 0):
                raise PassportContractError(
                    "SELECTION_UNIT_COUNT", f"{row['input_id']}:{count_field}={value!r}"
                )
        referenced = [
            value.strip()
            for value in [row["terminal_gate_input_id"], *row["parent_input_ids"].split(";")]
            if value.strip()
        ]
        unresolved = sorted(set(referenced) - input_id_set)
        if unresolved:
            raise PassportContractError(
                "SELECTION_INPUT_REFERENCE", f"{row['input_id']}: {unresolved}"
            )
        artifact = Path(row["artifact_path"])
        if not artifact.is_absolute():
            if not allow_fixture:
                raise PassportContractError("SELECTION_IMMUTABLE_PATH", "production artifact paths must be absolute")
            artifact = selection_path.parent / artifact
        if not artifact.is_file():
            raise PassportContractError("SELECTION_ARTIFACT_MISSING", str(artifact))
        if sha256_file(artifact) != row["artifact_sha256"]:
            raise PassportContractError("SELECTION_ARTIFACT_HASH", str(artifact))
        if artifact.stat().st_size != int(row["artifact_bytes"]):
            raise PassportContractError("SELECTION_ARTIFACT_BYTES", str(artifact))
        producer = Path(row["producer_script"])
        if not producer.is_absolute():
            producer = PROJECT_ROOT / producer
        if not producer.is_file():
            raise PassportContractError("SELECTION_PRODUCER_MISSING", str(producer))
        if sha256_file(producer) != row["producer_sha256"]:
            raise PassportContractError("SELECTION_PRODUCER_HASH", str(producer))
    for signed_key in (
        "gen_terminal_closure_input_id", "gen_terminal_ready_input_id",
        "plan13_terminal_ready_input_id",
    ):
        value = str(signature[signed_key]).strip()
        if not is_fixture and value not in input_id_set:
            raise PassportContractError(
                "SELECTION_SIGNED_GATE_REFERENCE", f"{signed_key}={value!r}"
            )
    return rows


def controlled_vocabulary_rows() -> list[dict[str, str]]:
    rows = []
    for vocabulary, values in VOCABULARIES.items():
        for value, definition in values.items():
            rows.append(
                {
                    "vocabulary_name": vocabulary,
                    "value": value,
                    "definition": definition,
                    "permitted_source": (
                        "frozen upstream verdict"
                        if vocabulary in {"call_state", "testability_state", "provenance_state"}
                        else "PASS deterministic rulebook"
                        if vocabulary == "role_hypothesis"
                        else "frozen contract or coordinator selection"
                    ),
                    "prohibited_interpretation": {
                        "call_state": "Do not convert missingness or nonsignificance into a negative call.",
                        "testability_state": "Do not upgrade testability from a different assay.",
                        "provenance_state": "Do not treat source dependence as positive or negative evidence.",
                        "role_hypothesis": "Hypothesis only; never a causal-gene or best-target label.",
                    }.get(vocabulary, "Do not convert this categorical field into a score or rank."),
                }
            )
    return rows


def rulebook_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for evidence_class, (role, role_rule, _) in ROLE_ROUTES.items():
        rows.append(
            {
                "rule_id": role_rule,
                "rule_type": "role_hypothesis",
                "rule_version": "1.0.0",
                "condition_json": json.dumps(
                    {"primary_evidence_class": {"equals": evidence_class}},
                    sort_keys=True, separators=(",", ":"),
                ),
                "output_value": role,
                "biological_model": "",
                "context": "",
                "perturbation": "",
                "primary_readout": "",
                "falsifying_outcome": "",
                "prohibited_use": "Cannot alter frozen evidence calls or imply causality.",
            }
        )
    for rule_id, fields in EXPERIMENT_RULES.items():
        rows.append(
            {
                "rule_id": rule_id,
                "rule_type": "next_experiment",
                "rule_version": "1.0.0",
                "condition_json": json.dumps(
                    {"role_route": {"experiment_rule_id": rule_id}},
                    sort_keys=True, separators=(",", ":"),
                ),
                "output_value": "one_discriminating_experiment",
                **fields,
                "prohibited_use": "Research proposal only; not a therapeutic-success or clinical recommendation.",
            }
        )
    # ROLE_UNRESOLVED is reused for three evidence classes; deduplicate exact rule rows.
    unique: OrderedDict[str, dict[str, str]] = OrderedDict()
    for row in rows:
        if row["rule_id"] not in unique:
            unique[row["rule_id"]] = row
        elif unique[row["rule_id"]] != row:
            # A shared rule is allowed only when its output is identical. Its condition
            # is widened deterministically rather than assigned a hidden precedence.
            old = unique[row["rule_id"]]
            old_condition = json.loads(old["condition_json"])
            new_condition = json.loads(row["condition_json"])
            old_values = old_condition.get("primary_evidence_class", {}).get("in", [])
            if not old_values:
                old_values = [old_condition["primary_evidence_class"]["equals"]]
            new_value = new_condition["primary_evidence_class"]["equals"]
            old["condition_json"] = json.dumps(
                {"primary_evidence_class": {"in": sorted(set(old_values + [new_value]))}},
                sort_keys=True, separators=(",", ":"),
            )
    return sorted(unique.values(), key=lambda row: row["rule_id"])


def _storage_type(value: pa.DataType | str) -> str:
    return str(value)


def data_dictionary_rows(include_production: bool = False) -> list[dict[str, str]]:
    important_metadata = {
        "passport_id": ("", "never null", "Catalog entry ID", "Stable release-plus-Ensembl key."),
        "ensembl_id": ("", "never null", "Ensembl gene ID", "Version-stripped human Ensembl primary key."),
        "symbol": ("", "nullable display alias", "Gene symbol", "Display/search alias; not a database key."),
        "call_state": ("call_state", "never null", "Assay call", "Exact frozen upstream assay verdict."),
        "negative_call_rule_id": ("", "required only for tested_negative", "Negative-call rule", "Explicit upstream equivalence, noninferiority, or informative-negative rule identifier."),
        "negative_decision_boundary": ("", "required only for tested_negative", "Negative decision boundary", "Prespecified effect boundary used by the informative-negative rule."),
        "negative_margin": ("", "required only for tested_negative", "Negative margin", "Observed nonnegative margin by which the informative-negative rule passed."),
        "negative_call_passed": ("", "true only for tested_negative; null otherwise", "Negative rule passed", "Explicit Boolean authorization; nonsignificance alone never sets this field."),
        "testability_state": ("testability_state", "never null", "Testability", "Originating workstream testability state."),
        "testability_reason": ("testability_state", "never null", "Testability reason", "Machine-readable eligibility or failure reason."),
        "provenance_state": ("provenance_state", "never null", "Source dependence", "Resolved source-overlap class."),
        "source_dependent": ("", "never null", "Source dependent", "Boolean derived only as provenance_state != independent."),
        "estimate": ("", "null when no assay-native estimate", "Estimate", "Native assay estimate; never cross-assay standardized by PASS."),
        "standard_error": ("", "nullable", "Standard error", "Native-assay uncertainty when deposited."),
        "p_value": ("", "nullable; distinct from zero", "P value", "Frozen source p value; PASS never recomputes it."),
        "q_value": ("", "nullable; distinct from zero", "Adjusted p value", "Frozen source multiplicity-adjusted value."),
        "role_hypothesis": ("role_hypothesis", "never null", "Mechanistic-role hypothesis", "Non-scoring hypothesis routed by a deterministic rule ID."),
        "primary_evidence_class": ("primary_evidence_class", "never null", "Primary evidence class", "Frozen class supplied by an accepted source."),
        "direction": ("direction", "never null", "Direction", "Native-assay direction; not compared across incompatible units."),
        "evidence_domain": ("evidence_domain", "never null", "Evidence domain", "Non-scoring assay-domain label."),
        "gene_evidence_status": ("", "never null", "Gene evidence", "Whether accepted evidence is populated at gene grain for this domain."),
        "program_context_status": ("", "never null", "Program context", "Whether accepted program definitions or context results are populated for this domain."),
        "assay_status": ("", "never null", "Assay status", "Whether accepted assay-, class-, or program-grain status is populated."),
        "coverage_boundary": ("", "never null", "Coverage boundary", "Native representation grain; never permission to expand a program or assay call to member genes."),
        "gene_evidence_authorized": ("", "never null", "Gene evidence authorized", "True only when a signed source supplied a gene-grain evidence row in the narrowed candidate scope."),
        "member_gene_expansion_authorized": ("", "always false", "Member-gene expansion", "Frozen false: program membership/context never creates member-gene evidence."),
        "node_type": ("source_node_type", "never null", "Source-node type", "Source graph node category."),
        "edge_type": ("source_edge_type", "never null", "Source-edge type", "Directed source graph relationship."),
    }
    rows: list[dict[str, str]] = []
    parquet_schemas = (
        {**PRODUCTION_PARQUET_SCHEMAS, **PORTAL_PARQUET_SCHEMAS}
        if include_production else PARQUET_SCHEMAS
    )
    tsv_schemas = PRODUCTION_TSV_SCHEMAS if include_production else TSV_SCHEMAS
    for table_name, schema in parquet_schemas.items():
        for column, arrow_type in schema.items():
            vocab, nulls, label, description = important_metadata.get(
                column,
                ("", "nullable only when source semantics permit", column.replace("_", " ").title(), "Provenance-preserving Gene Catalog field."),
            )
            unit = "native assay unit" if column in {"estimate", "standard_error", "ci_lower", "ci_upper"} else "not applicable"
            rows.append(
                {
                    "table_name": table_name,
                    "column_name": column,
                    "storage_type": _storage_type(arrow_type),
                    "required": "true",
                    "allowed_values_vocabulary": vocab,
                    "unit_semantics": unit,
                    "null_semantics": nulls,
                    "ui_label": label,
                    "description": description,
                }
            )
    for table_name, columns in tsv_schemas.items():
        if table_name == "passport_data_dictionary.tsv":
            continue
        for column in columns:
            vocab, nulls, label, description = important_metadata.get(
                column,
                ("", "empty only where contract permits", column.replace("_", " ").title(), "Provenance-preserving Gene Catalog field."),
            )
            rows.append(
                {
                    "table_name": table_name,
                    "column_name": column,
                    "storage_type": "string",
                    "required": "true",
                    "allowed_values_vocabulary": vocab,
                    "unit_semantics": "not applicable",
                    "null_semantics": nulls,
                    "ui_label": label,
                    "description": description,
                }
            )
    return rows


def source_row_hash(row: Mapping[str, Any]) -> str:
    excluded = {"source_row_sha256"}
    # Hash the exact TSV-level representation so typed producer values and
    # strings read back from the frozen artifact have one canonical digest.
    canonical = {
        key: "" if row[key] is None else str(row[key])
        for key in sorted(row)
        if key not in excluded
    }
    return stable_sha256(canonical)


def synthetic_source_rows(include_tested_negative: bool = True) -> list[dict[str, Any]]:
    base_rows = [
        {
            "source_input_row_id": "SRC001", "ensembl_id": "ENSG90000000001", "symbol": "COLLIDE",
            "primary_evidence_class": "genetically_anchored", "evidence_domain": "genetics",
            "assay": "SuSiE_COLOC", "dataset_id": "FIXTURE_DIRECT_DX", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "MASLD diagnosis", "context": "bulk liver cis regulation",
            "biological_unit": "GWAS cohort and eQTL donor", "contrast_or_exposure": "trait-locus colocalization",
            "effect_unit": "PP.H4", "estimate": 0.82, "standard_error": None, "ci_lower": None,
            "ci_upper": None, "p_value": None, "q_value": None, "direction": "not_directional",
            "source_testability_state": "testable", "source_call_state": "supported", "gate_id": "FIXTURE_GATE_01",
            "n_biological_units": 1183, "n_technical_units": 1, "source_provenance_state": "independent",
            "allowed_wording": "Genetically anchored in the fixture source.", "limitation": "Synthetic PASS-01 fixture only.",
        },
        {
            "source_input_row_id": "SRC002", "ensembl_id": "ENSG90000000002", "symbol": "COLLIDE",
            "primary_evidence_class": "established_state_associated", "evidence_domain": "functional",
            "assay": "CRISPR_palmitate", "dataset_id": "FIXTURE_FUNCTIONAL", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "palmitate-specific fitness", "context": "lipotoxic challenge",
            "biological_unit": "gene", "contrast_or_exposure": "palmitate minus vehicle",
            "effect_unit": "log2 fold change", "estimate": 0.02, "standard_error": 0.10, "ci_lower": -0.18,
            "ci_upper": 0.22, "p_value": 0.82, "q_value": 0.91, "direction": "no_effect",
            "source_testability_state": "testable", "source_call_state": "tested_negative", "gate_id": "FIXTURE_GATE_02",
            "n_biological_units": 4, "n_technical_units": 8, "source_provenance_state": "partially_dependent",
            "allowed_wording": "Explicit tested negative in an adequately testable synthetic challenge.",
            "limitation": "Synthetic PASS-01 fixture only.",
        },
        {
            "source_input_row_id": "SRC003", "ensembl_id": "ENSG90000000003", "symbol": "DISCORD",
            "primary_evidence_class": "discordant", "evidence_domain": "proteomics",
            "assay": "paired_RNA_protein", "dataset_id": "FIXTURE_PAIRED", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "MASH", "context": "matched tissue",
            "biological_unit": "donor", "contrast_or_exposure": "MASH minus control",
            "effect_unit": "standardized paired effect", "estimate": 1.2, "standard_error": 0.3, "ci_lower": 0.6,
            "ci_upper": 1.8, "p_value": 0.001, "q_value": 0.01, "direction": "discordant",
            "source_testability_state": "testable", "source_call_state": "discordant", "gate_id": "FIXTURE_GATE_03",
            "n_biological_units": 12, "n_technical_units": 24, "source_provenance_state": "reused_source",
            "allowed_wording": "Upstream fixture explicitly supplies RNA-protein discordance.",
            "limitation": "Synthetic PASS-01 fixture only.",
        },
        {
            "source_input_row_id": "SRC004", "ensembl_id": "ENSG90000000004", "symbol": "UNCERT",
            "primary_evidence_class": "context_supported", "evidence_domain": "cell_context",
            "assay": "interaction_eQTL", "dataset_id": "FIXTURE_CONTEXT", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "disease-state interaction", "context": "hepatocyte state",
            "biological_unit": "donor", "contrast_or_exposure": "genotype by state interaction",
            "effect_unit": "interaction beta", "estimate": 0.18, "standard_error": 0.15, "ci_lower": -0.11,
            "ci_upper": 0.47, "p_value": 0.24, "q_value": 0.61, "direction": "unknown",
            "source_testability_state": "underpowered_source", "source_call_state": "indeterminate", "gate_id": "FIXTURE_GATE_04",
            "n_biological_units": 7, "n_technical_units": 7, "source_provenance_state": "independent",
            "allowed_wording": "Indeterminate context signal in an underpowered synthetic source.",
            "limitation": "Synthetic PASS-01 fixture only.",
        },
        {
            "source_input_row_id": "SRC005", "ensembl_id": "ENSG90000000005", "symbol": "COVERAGE",
            "primary_evidence_class": "untested", "evidence_domain": "spatial",
            "assay": "spatial_program_test", "dataset_id": "FIXTURE_SPATIAL", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "local lipid exposure", "context": "hepatocyte spots",
            "biological_unit": "donor", "contrast_or_exposure": "lipid exposure",
            "effect_unit": "standardized donor slope", "estimate": None, "standard_error": None, "ci_lower": None,
            "ci_upper": None, "p_value": None, "q_value": None, "direction": "unknown",
            "source_testability_state": "insufficient_gene_coverage", "source_call_state": "untestable", "gate_id": "FIXTURE_GATE_05",
            "n_biological_units": 0, "n_technical_units": 0, "source_provenance_state": "independent",
            "allowed_wording": "Untestable because the synthetic assay lacks sufficient gene coverage.",
            "limitation": "No gene-level negative is authorized.",
        },
        {
            "source_input_row_id": "SRC006", "ensembl_id": "ENSG90000000006", "symbol": "OUTSCOPE",
            "primary_evidence_class": "unresolved", "evidence_domain": "functional",
            "assay": "lineage_specific_assay", "dataset_id": "FIXTURE_OUTSCOPE", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "lineage response", "context": "non-target lineage",
            "biological_unit": "not applicable", "contrast_or_exposure": "not applicable",
            "effect_unit": "not applicable", "estimate": None, "standard_error": None, "ci_lower": None,
            "ci_upper": None, "p_value": None, "q_value": None, "direction": "unknown",
            "source_testability_state": "assay_out_of_scope", "source_call_state": "not_applicable", "gate_id": "FIXTURE_GATE_06",
            "n_biological_units": 0, "n_technical_units": 0, "source_provenance_state": "independent",
            "allowed_wording": "The synthetic lineage assay is not applicable by design.",
            "limitation": "Not missing data and not a negative result.",
        },
        {
            "source_input_row_id": "SRC007", "ensembl_id": "ENSG90000000007", "symbol": "STATE",
            "primary_evidence_class": "concordant", "evidence_domain": "transcriptomics",
            "assay": "donor_pseudobulk", "dataset_id": "FIXTURE_STATE", "source_release_id": "fixture-source-v1",
            "source_node_id": "node:derived:fixture", "phenotype": "MASH", "context": "cross-sectional disease state",
            "biological_unit": "donor", "contrast_or_exposure": "MASH minus healthy",
            "effect_unit": "log2 fold change", "estimate": 0.7, "standard_error": 0.12, "ci_lower": 0.46,
            "ci_upper": 0.94, "p_value": 1e-6, "q_value": 1e-4, "direction": "positive",
            "source_testability_state": "testable", "source_call_state": "supported", "gate_id": "FIXTURE_GATE_07",
            "n_biological_units": 20, "n_technical_units": 20, "source_provenance_state": "partially_dependent",
            "allowed_wording": "Concordant established-state association in the synthetic source.",
            "limitation": "Cross-sectional association, not longitudinal progression.",
        },
    ]
    for row in base_rows:
        row["negative_call_rule_id"] = None
        row["negative_decision_boundary"] = None
        row["negative_margin"] = None
        row["negative_call_passed"] = None
    explicit_negative = next(
        row for row in base_rows if row["source_input_row_id"] == "SRC002"
    )
    explicit_negative.update(
        {
            "negative_call_rule_id": "FIXTURE_EQUIVALENCE_V1",
            "negative_decision_boundary": 0.25,
            "negative_margin": 0.03,
            "negative_call_passed": True,
        }
    )
    if not include_tested_negative:
        explicit_negative.update(
            {
                "source_call_state": "indeterminate",
                "direction": "unknown",
                "allowed_wording": "Complete synthetic test with no authorized positive or informative-negative call.",
                "negative_call_rule_id": None,
                "negative_decision_boundary": None,
                "negative_margin": None,
                "negative_call_passed": None,
            }
        )
    for row in base_rows:
        row["source_row_sha256"] = source_row_hash(row)
    return base_rows


def _source_columns() -> list[str]:
    return list(synthetic_source_rows()[0])


def write_fixture_selection(bundle: Path, source_path: Path) -> tuple[Path, Path]:
    selection_path = bundle / "passport_input_selection.tsv"
    selection_rows = [
        {
            "input_id": "FIXTURE_INPUT_001",
            "workstream_id": "PASS01_FIXTURE",
            "artifact_role": "synthetic_frozen_evidence",
            "adapter_id": "synthetic_gene_evidence_v2",
            "artifact_grain": "gene",
            "artifact_path": source_path.name,
            "artifact_sha256": sha256_file(source_path),
            "artifact_bytes": source_path.stat().st_size,
            "producer_script": LOGICAL_PRODUCER_ID,
            "producer_sha256": sha256_file(SCRIPT_PATH),
            "source_release_id": "fixture-source-v1",
            "gate_verdict": "include",
            "terminal_gate_input_id": "",
            "allowed_claim_wording": "Synthetic PASS-01 contract fixture only.",
            "limitation": "Synthetic input cannot unlock production.",
            "source_datasets": "FIXTURE_ONLY",
            "source_cohorts": "FIXTURE_ONLY",
            "source_publication": "Synthetic PASS-01 fixture",
            "source_url": "fixture://source",
            "biological_unit": "synthetic mixed units",
            "n_biological_units": "7",
            "n_technical_units": "7",
            "provenance_state": "independent",
            "claim_scope": "gene_evidence",
            "inclusion_destination": "passport_only",
            "selection_decision": "accepted",
            "parent_input_ids": "",
            "fixture_only": "true",
        }
    ]
    atomic_write_tsv(selection_path, selection_rows, SELECTION_REQUIRED_COLUMNS)
    signature_path = bundle / "passport_input_selection.signature.json"
    atomic_write_json(
        signature_path,
        {
            "attestation_version": SELECTION_ATTESTATION_VERSION,
            "selection_sha256": sha256_file(selection_path),
            "signed_by": "PASS01_SYNTHETIC_FIXTURE",
            "signed_at_utc": FIXTURE_SIGNED_AT,
            "decision_register_id": "FIXTURE_DECISION_001",
            "authority_document": "docs/PAPER.md",
            "fixture_only": True,
            "coordinator_attested": False,
            "analysis_release_id": ANALYSIS_RELEASE_ID,
            "gen_terminal_closure_input_id": "",
            "gen_terminal_ready_input_id": "",
            "plan13_terminal_ready_input_id": "",
        },
    )
    return selection_path, signature_path


def build_source_graph(
    source_artifact_sha: str, source_rows: Sequence[Mapping[str, Any]]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    nodes = [
        {
            "source_node_id": "node:publication:fixture", "node_type": "publication",
            "label": "Synthetic PASS-01 source publication", "source_release_id": "fixture-source-v1",
            "uri": "fixture://publication", "artifact_sha256": source_artifact_sha,
            "access_status": "synthetic_public_fixture", "biological_unit": "not applicable",
            "source_datasets": "FIXTURE_ONLY",
        },
        {
            "source_node_id": "node:cohort:fixture", "node_type": "cohort",
            "label": "Synthetic donor cohort", "source_release_id": "fixture-source-v1",
            "uri": "fixture://cohort", "artifact_sha256": source_artifact_sha,
            "access_status": "synthetic_public_fixture", "biological_unit": "donor",
            "source_datasets": "FIXTURE_ONLY",
        },
        {
            "source_node_id": "node:dataset:fixture", "node_type": "dataset",
            "label": "Synthetic mixed-assay dataset", "source_release_id": "fixture-source-v1",
            "uri": "fixture://dataset", "artifact_sha256": source_artifact_sha,
            "access_status": "synthetic_public_fixture", "biological_unit": "source-defined",
            "source_datasets": "FIXTURE_ONLY",
        },
        {
            "source_node_id": "node:release:fixture", "node_type": "frozen_release",
            "label": "Synthetic frozen workstream release", "source_release_id": "fixture-source-v1",
            "uri": "fixture://release", "artifact_sha256": source_artifact_sha,
            "access_status": "synthetic_public_fixture", "biological_unit": "source-defined",
            "source_datasets": "FIXTURE_ONLY",
        },
        {
            "source_node_id": "node:derived:fixture", "node_type": "derived_object",
            "label": "Synthetic accepted evidence artifact", "source_release_id": "fixture-source-v1",
            "uri": "fixture://derived", "artifact_sha256": source_artifact_sha,
            "access_status": "synthetic_public_fixture", "biological_unit": "source-defined",
            "source_datasets": "FIXTURE_ONLY",
        },
    ]
    edge_specs = [
        ("node:dataset:fixture", "node:cohort:fixture", "derived_from", "Dataset derives from cohort."),
        ("node:dataset:fixture", "node:publication:fixture", "derived_from", "Dataset is described by publication."),
        ("node:release:fixture", "node:dataset:fixture", "derived_from", "Frozen release derives from dataset."),
        ("node:derived:fixture", "node:release:fixture", "derived_from", "Accepted artifact derives from frozen release."),
        ("node:cohort:fixture", "node:dataset:fixture", "shares_samples_with", "Visible synthetic source-dependence edge."),
    ]
    edges = []
    for index, (source, target, edge_type, note) in enumerate(edge_specs, start=1):
        edges.append(
            {
                "source_edge_id": f"edge:fixture:{index:02d}",
                "from_node_id": source,
                "to_node_id": target,
                "edge_type": edge_type,
                "evidence_result_id": "",
                "note": note,
            }
        )
    for row in source_rows:
        edges.append(
            {
                "source_edge_id": f"edge:tested:{row['source_input_row_id']}",
                "from_node_id": "node:derived:fixture",
                "to_node_id": "node:dataset:fixture",
                "edge_type": "tested_by",
                "evidence_result_id": f"evidence:{row['source_input_row_id']}",
                "note": "Synthetic inferential result resolves to its accepted source dataset.",
            }
        )
    return nodes, edges


def route_hypothesis(primary_class: str) -> tuple[str, str, str]:
    try:
        return ROLE_ROUTES[primary_class]
    except KeyError as exc:
        raise PassportContractError("ROLE_ROUTE_UNDEFINED", primary_class) from exc


def build_passport_rows(source_rows: list[dict[str, Any]], source_artifact_sha: str):
    genes: OrderedDict[str, dict[str, Any]] = OrderedDict()
    evidence: list[dict[str, Any]] = []
    coverage: list[dict[str, Any]] = []
    experiments: list[dict[str, str]] = []
    symbol_counts: dict[str, int] = {}
    for row in source_rows:
        symbol_counts[row["symbol"]] = symbol_counts.get(row["symbol"], 0) + 1

    for row in source_rows:
        ensembl = row["ensembl_id"].split(".", 1)[0]
        passport_id = f"{ANALYSIS_RELEASE_ID}:{ensembl}"
        role, role_rule, experiment_rule = route_hypothesis(row["primary_evidence_class"])
        genes[ensembl] = {
            "analysis_release_id": ANALYSIS_RELEASE_ID,
            "passport_id": passport_id,
            "ensembl_id": ensembl,
            "symbol": row["symbol"],
            "symbol_collision": symbol_counts[row["symbol"]] > 1,
            "primary_evidence_class": row["primary_evidence_class"],
            "role_hypothesis": role,
            "role_rule_id": role_rule,
            "next_experiment_rule_id": experiment_rule,
            "allowed_wording": row["allowed_wording"],
            "limitation": row["limitation"],
        }
        evidence_result_id = f"evidence:{row['source_input_row_id']}"
        provenance = row["source_provenance_state"]
        call = row["source_call_state"]
        testability = row["source_testability_state"]
        evidence.append(
            {
                "analysis_release_id": ANALYSIS_RELEASE_ID,
                "evidence_result_id": evidence_result_id,
                "passport_id": passport_id,
                "ensembl_id": ensembl,
                "symbol": row["symbol"],
                "evidence_domain": row["evidence_domain"],
                "assay": row["assay"],
                "dataset_id": row["dataset_id"],
                "source_node_id": row["source_node_id"],
                "source_release_id": row["source_release_id"],
                "phenotype": row["phenotype"],
                "context": row["context"],
                "biological_unit": row["biological_unit"],
                "contrast_or_exposure": row["contrast_or_exposure"],
                "effect_unit": row["effect_unit"],
                "estimate": row["estimate"],
                "standard_error": row["standard_error"],
                "ci_lower": row["ci_lower"],
                "ci_upper": row["ci_upper"],
                "p_value": row["p_value"],
                "q_value": row["q_value"],
                "direction": row["direction"],
                "testability_state": testability,
                "testability_reason": testability,
                "call_state": call,
                "negative_call_rule_id": row["negative_call_rule_id"],
                "negative_decision_boundary": row["negative_decision_boundary"],
                "negative_margin": row["negative_margin"],
                "negative_call_passed": row["negative_call_passed"],
                "gate_id": row["gate_id"],
                "n_biological_units": row["n_biological_units"],
                "n_technical_units": row["n_technical_units"],
                "provenance_state": provenance,
                "source_dependent": provenance != "independent",
                "source_artifact_sha256": source_artifact_sha,
                "source_input_row_id": row["source_input_row_id"],
                "source_call_state": call,
                "source_testability_state": testability,
                "source_provenance_state": provenance,
                "source_negative_call_rule_id": row["negative_call_rule_id"],
                "source_negative_decision_boundary": row["negative_decision_boundary"],
                "source_negative_margin": row["negative_margin"],
                "source_negative_call_passed": row["negative_call_passed"],
                "source_row_sha256": row["source_row_sha256"],
                "allowed_wording": row["allowed_wording"],
                "limitation": row["limitation"],
            }
        )
        coverage.append(
            {
                "analysis_release_id": ANALYSIS_RELEASE_ID,
                "coverage_result_id": f"coverage:{row['source_input_row_id']}",
                "passport_id": passport_id,
                "ensembl_id": ensembl,
                "symbol": row["symbol"],
                "evidence_domain": row["evidence_domain"],
                "assay": row["assay"],
                "dataset_id": row["dataset_id"],
                "source_release_id": row["source_release_id"],
                "testability_state": testability,
                "testability_reason": testability,
                "call_state": call,
                "n_biological_units": row["n_biological_units"],
                "n_technical_units": row["n_technical_units"],
                "coverage_denominator": row["biological_unit"],
                "limitation": row["limitation"],
            }
        )

    for gene in genes.values():
        experiment = EXPERIMENT_RULES[gene["next_experiment_rule_id"]]
        experiments.append(
            {
                "analysis_release_id": ANALYSIS_RELEASE_ID,
                "passport_id": gene["passport_id"],
                "ensembl_id": gene["ensembl_id"],
                "symbol": gene["symbol"],
                "experiment_rule_id": gene["next_experiment_rule_id"],
                "role_rule_id": gene["role_rule_id"],
                **experiment,
                "disclaimer": DISCLAIMER,
            }
        )
    return list(genes.values()), evidence, coverage, experiments


def write_release_manifest(bundle: Path) -> Path:
    manifest_path = bundle / "passport_release_manifest.tsv"
    rows = []
    for path in sorted(bundle.iterdir()):
        if not path.is_file() or path == manifest_path:
            continue
        rows.append(
            {
                "relative_path": path.name,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "producer": LOGICAL_PRODUCER_ID,
                "environment": f"spatial:pandas-{pd.__version__}:pyarrow-{pa.__version__}",
                "producer_sha256": sha256_file(SCRIPT_PATH),
                "upstream_artifact_sha256": sha256_file(bundle / "fixture_source_evidence.tsv"),
                "release_status": "synthetic_fixture_not_for_promotion",
            }
        )
    atomic_write_tsv(manifest_path, rows, TSV_SCHEMAS[manifest_path.name])
    return manifest_path


def build_positive_fixture(
    bundle: Path, candidate_root: Path, include_tested_negative: bool = True
) -> None:
    bundle = assert_fixture_path(bundle, candidate_root)
    bundle.mkdir(parents=True, exist_ok=True)
    source_rows = synthetic_source_rows(include_tested_negative=include_tested_negative)
    source_path = bundle / "fixture_source_evidence.tsv"
    atomic_write_tsv(source_path, source_rows, _source_columns())
    source_sha = sha256_file(source_path)
    write_fixture_selection(bundle, source_path)
    validate_signed_selection(bundle / "passport_input_selection.tsv", allow_fixture=True)

    genes, evidence, coverage, experiments = build_passport_rows(source_rows, source_sha)
    for filename, rows in (
        ("passport_gene_index.parquet", genes),
        ("passport_evidence_long.parquet", evidence),
        ("passport_coverage_long.parquet", coverage),
    ):
        schema = PARQUET_SCHEMAS[filename]
        atomic_write_parquet(bundle / filename, coerce_dataframe(rows, schema), schema)

    nodes, edges = build_source_graph(source_sha, source_rows)
    atomic_write_tsv(bundle / "passport_source_nodes.tsv", nodes, TSV_SCHEMAS["passport_source_nodes.tsv"])
    atomic_write_tsv(bundle / "passport_source_edges.tsv", edges, TSV_SCHEMAS["passport_source_edges.tsv"])
    atomic_write_tsv(
        bundle / "passport_next_experiment.tsv", experiments,
        TSV_SCHEMAS["passport_next_experiment.tsv"],
    )
    atomic_write_tsv(
        bundle / "passport_rulebook.tsv", rulebook_rows(), TSV_SCHEMAS["passport_rulebook.tsv"]
    )
    atomic_write_tsv(
        bundle / "passport_controlled_vocabularies.tsv", controlled_vocabulary_rows(),
        TSV_SCHEMAS["passport_controlled_vocabularies.tsv"],
    )
    atomic_write_tsv(
        bundle / "passport_data_dictionary.tsv", data_dictionary_rows(),
        TSV_SCHEMAS["passport_data_dictionary.tsv"],
    )
    atomic_write_tsv(
        bundle / "passport_gate_status.tsv",
        [
            {
                "gate_id": "PASS01_SCHEMA_FIXTURE", "status": "pass",
                "promotion_allowed": "false", "reason": "Synthetic schema fixture only.",
            },
            {
                "gate_id": "PASS00_SIGNED_REAL_SELECTION", "status": "not_evaluated_fixture",
                "promotion_allowed": "false", "reason": "Fixture attestation cannot unlock production.",
            },
            {
                "gate_id": "PORTAL_PRODUCTION_READY", "status": "blocked_pass01_only",
                "promotion_allowed": "false", "reason": "No real evidence was populated and PASS-02 through PASS-06 are incomplete.",
            },
        ],
        TSV_SCHEMAS["passport_gate_status.tsv"],
    )
    write_release_manifest(bundle)


def _copy_bundle(source: Path, destination: Path, candidate_root: Path) -> None:
    assert_fixture_path(destination, candidate_root)
    shutil.copytree(source, destination, dirs_exist_ok=True)


def _rewrite_parquet(path: Path, mutate) -> None:
    frame = pd.read_parquet(path, engine="pyarrow")
    frame = mutate(frame.copy())
    schema = PARQUET_SCHEMAS[path.name]
    # Hidden-column fixture deliberately carries one extra field.
    if set(frame.columns) == set(schema):
        frame = coerce_dataframe(frame.to_dict("records"), schema)
        atomic_write_parquet(path, frame, schema)
    else:
        sort_columns = [column for column in ("passport_id", "evidence_result_id") if column in frame]
        if sort_columns:
            frame = frame.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)
        with tempfile.NamedTemporaryFile(mode="wb", dir=path.parent, delete=False) as handle:
            tmp = Path(handle.name)
        try:
            frame.to_parquet(tmp, engine="pyarrow", index=False, compression="zstd")
            os.replace(tmp, path)
        finally:
            if tmp.exists():
                tmp.unlink()


def build_negative_fixtures(positive: Path, fixture_root: Path, candidate_root: Path) -> dict[str, str]:
    expectations = OrderedDict(
        [
            ("negative_missing_as_negative", "TESTED_NEGATIVE_REQUIRES_TESTABLE"),
            ("negative_semantic_spoof", "TESTED_NEGATIVE_RULE_REQUIRED"),
            ("negative_multi_ensembl_ambiguity", "MULTI_ENSEMBL_AMBIGUITY"),
            ("negative_duplicate_id", "DUPLICATE_PASSPORT_ID"),
            ("negative_collapsed_symbol", "SOURCE_GENE_LOSS_OR_GAIN"),
            ("negative_missing_provenance", "UNKNOWN_PROVENANCE"),
            ("negative_altered_source_call", "SOURCE_CALL_MUTATION"),
            ("negative_hidden_score_rank", "FORBIDDEN_SCORE_OR_RANK"),
        ]
    )
    for name in expectations:
        destination = fixture_root / name
        _copy_bundle(positive, destination, candidate_root)

    missing_negative = fixture_root / "negative_missing_as_negative" / "passport_evidence_long.parquet"
    def mutate_missing(frame):
        mask = frame["source_input_row_id"] == "SRC005"
        frame.loc[mask, "call_state"] = "tested_negative"
        frame.loc[mask, "source_call_state"] = "tested_negative"
        return frame
    _rewrite_parquet(missing_negative, mutate_missing)

    semantic_spoof = fixture_root / "negative_semantic_spoof" / "passport_evidence_long.parquet"
    def mutate_semantic_spoof(frame):
        mask = frame["source_input_row_id"] == "SRC001"
        frame.loc[mask, "call_state"] = "tested_negative"
        frame.loc[mask, "source_call_state"] = "tested_negative"
        return frame
    _rewrite_parquet(semantic_spoof, mutate_semantic_spoof)

    ambiguous_dir = fixture_root / "negative_multi_ensembl_ambiguity"
    ambiguous_source = ambiguous_dir / "fixture_source_evidence.tsv"
    ambiguous_rows = read_tsv_rows(ambiguous_source)
    ambiguous_rows[0]["ensembl_id"] = "ENSG90000000001;ENSG90000000002"
    atomic_write_tsv(ambiguous_source, ambiguous_rows, _source_columns())
    ambiguous_selection_path = ambiguous_dir / "passport_input_selection.tsv"
    ambiguous_selection = read_tsv_rows(ambiguous_selection_path)
    ambiguous_selection[0]["artifact_sha256"] = sha256_file(ambiguous_source)
    ambiguous_selection[0]["artifact_bytes"] = str(ambiguous_source.stat().st_size)
    atomic_write_tsv(
        ambiguous_selection_path, ambiguous_selection, SELECTION_REQUIRED_COLUMNS
    )
    signature_path = ambiguous_dir / "passport_input_selection.signature.json"
    signature = json.loads(signature_path.read_text(encoding="utf-8"))
    signature["selection_sha256"] = sha256_file(ambiguous_selection_path)
    atomic_write_json(signature_path, signature)

    duplicate = fixture_root / "negative_duplicate_id" / "passport_gene_index.parquet"
    _rewrite_parquet(duplicate, lambda frame: pd.concat([frame, frame.iloc[[0]]], ignore_index=True))

    collapsed_dir = fixture_root / "negative_collapsed_symbol"
    collapsed_ens = "ENSG90000000002"
    def mutate_collapsed_index(frame):
        frame = frame[frame["ensembl_id"] != collapsed_ens].reset_index(drop=True)
        frame.loc[frame["symbol"] == "COLLIDE", "symbol_collision"] = False
        return frame
    _rewrite_parquet(
        collapsed_dir / "passport_gene_index.parquet",
        mutate_collapsed_index,
    )
    _rewrite_parquet(
        collapsed_dir / "passport_evidence_long.parquet",
        lambda frame: frame[frame["ensembl_id"] != collapsed_ens].reset_index(drop=True),
    )
    _rewrite_parquet(
        collapsed_dir / "passport_coverage_long.parquet",
        lambda frame: frame[frame["ensembl_id"] != collapsed_ens].reset_index(drop=True),
    )
    next_path = collapsed_dir / "passport_next_experiment.tsv"
    with next_path.open(newline="", encoding="utf-8") as handle:
        next_rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["ensembl_id"] != collapsed_ens]
    atomic_write_tsv(next_path, next_rows, TSV_SCHEMAS[next_path.name])

    provenance = fixture_root / "negative_missing_provenance" / "passport_evidence_long.parquet"
    def mutate_provenance(frame):
        frame.loc[frame.index[0], "provenance_state"] = pd.NA
        return frame
    _rewrite_parquet(provenance, mutate_provenance)

    altered = fixture_root / "negative_altered_source_call" / "passport_evidence_long.parquet"
    def mutate_call(frame):
        frame.loc[frame.index[0], "call_state"] = "discordant"
        return frame
    _rewrite_parquet(altered, mutate_call)

    hidden = fixture_root / "negative_hidden_score_rank" / "passport_gene_index.parquet"
    def mutate_hidden(frame):
        frame["priority_score"] = range(1, len(frame) + 1)
        return frame
    _rewrite_parquet(hidden, mutate_hidden)

    for name in expectations:
        write_release_manifest(fixture_root / name)
    atomic_write_tsv(
        fixture_root / "negative_fixture_expectations.tsv",
        [{"fixture": name, "expected_error_code": code} for name, code in expectations.items()],
        ["fixture", "expected_error_code"],
    )
    return dict(expectations)


def build_fixture_suite(candidate_root: Path = DEFAULT_CANDIDATE_ROOT) -> Path:
    candidate = assert_candidate_root(candidate_root)
    fixture_root = assert_fixture_path(candidate / "fixtures", candidate)
    fixture_root.mkdir(parents=True, exist_ok=True)
    positive = fixture_root / "positive"
    repeat = fixture_root / "positive_repeat"
    zero_negative = fixture_root / "positive_zero_tested_negative"
    build_positive_fixture(positive, candidate)
    build_positive_fixture(repeat, candidate)
    build_positive_fixture(zero_negative, candidate, include_tested_negative=False)
    build_negative_fixtures(positive, fixture_root, candidate)
    return fixture_root


def compare_fixture_determinism(first: Path, second: Path) -> list[str]:
    first_files = sorted(path.name for path in first.iterdir() if path.is_file())
    second_files = sorted(path.name for path in second.iterdir() if path.is_file())
    if first_files != second_files:
        return ["file inventory differs"]
    return [name for name in first_files if sha256_file(first / name) != sha256_file(second / name)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("production", "fixtures"), default="production")
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    parser.add_argument("--selection", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidate = assert_candidate_root(args.candidate_root)
    if args.mode == "fixtures":
        fixture_root = build_fixture_suite(candidate)
        differences = compare_fixture_determinism(
            fixture_root / "positive", fixture_root / "positive_repeat"
        )
        if differences:
            raise PassportContractError("NONDETERMINISTIC_FIXTURE", ", ".join(differences))
        print(f"PASS: deterministic PASS-01 fixtures written beneath {fixture_root}")
        return

    selection = args.selection or candidate / "passport_input_selection.tsv"
    validate_signed_selection(selection, allow_fixture=False)
    print(
        "PASS: signed PASS-00 selection satisfies PASS-01 contract; "
        "PASS-02--PASS-05 assembly may proceed through build_evidence_passport_bundle.py"
    )


if __name__ == "__main__":
    try:
        main()
    except PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
