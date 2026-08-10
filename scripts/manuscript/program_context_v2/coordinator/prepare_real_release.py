#!/usr/bin/env python3
"""Prepare, but never execute or promote, real Plan 60 REL00--05 inputs."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

PROGRAM_ROOT = Path(__file__).resolve().parents[1]
if str(PROGRAM_ROOT) not in sys.path:
    sys.path.insert(0, str(PROGRAM_ROOT))

from adapter_contract import (  # noqa: E402
    ADAPTER_CONTRACT,
    ADAPTER_CONTRACT_SHA256,
    ADAPTER_PROVENANCE_FIELDS,
    ADAPTER_PROVENANCE_VERSION,
    EXPECTED_ADAPTER_IDS,
    adapter_input_paths,
)
from fibrosis_candidate_contract import (  # noqa: E402
    FIBROSIS_BUNDLE_FILES,
    FIBROSIS_PRIMARY_RELATIVE,
    fibrosis_artifact_id,
    fibrosis_artifact_role,
    fibrosis_snapshot_relpath,
    fixed_live_bundle_root,
    validate_fibrosis_candidate_bundle,
)
from release_common import ReleaseContractError  # noqa: E402

from coordinator_contract import (  # noqa: E402
    BASE_SELECTION_FIELDS,
    CANDIDATE_ID,
    CLOSURE_FIELDS,
    COMPOSITION_READY_FIELDS,
    COMPOSITION_RESULT_FIELDS,
    COMPOSITION_TESTABILITY_FIELDS,
    CORE_PACKAGE_FIELDS,
    CORE_PRODUCER_FIELDS,
    CORE_RUNTIME_FIELDS,
    CORE_IMPORT_FIELDS,
    CORE_SYSPATH_FIELDS,
    NMF_CONTINUOUS_FIELDS,
    NMF_SOURCE_MANIFEST_FIELDS,
    PREPARATION_CONTRACT_VERSION,
    PROTECTED_BASELINE_FIELDS,
    PROTECTED_SCOPE_FIELDS,
    SOURCE_EVIDENCE_FIELDS,
    WORKSTREAM_ORDER,
    CoordinatorContractError,
    RealPaths,
    atomic_write_json,
    atomic_write_tsv,
    build_cohort_overview_rows,
    build_composition_rows,
    build_genetics_rows,
    build_hotspot_rows,
    build_myojin_rows,
    build_nmf_continuous_supplement,
    build_passport_rows,
    build_spatial_rows,
    closure_rows,
    coordinator_signature,
    core_package_export_rows,
    core_import_resolution_rows,
    core_runtime_environment_rows,
    core_sys_path_rows,
    locked_release_blueprint,
    one_row_tsv,
    project_relative,
    protected_baseline_rows,
    protected_scope_rows,
    recursive_release_producer_rows,
    read_tsv_exact,
    read_tsv_flexible,
    selection_row,
    sha256_file,
    validate_all_handoffs,
    validate_nmf_continuous_bundle,
)


PREPARATION_MANIFEST_FIELDS = (
    "artifact_id",
    "relative_path",
    "sha256",
    "bytes",
    "artifact_role",
    "canonical_promotion_authorized",
)
SOURCE_ALIAS_AUDIT_FIELDS = (
    "source_id",
    "canonical_source_id",
    "resolution_basis",
    "alias_equivalence_status",
    "used_as_discovery",
    "used_as_evaluation",
    "used_in_independent_row",
    "independence_scope",
)


def read_parquet_records(path: Path) -> list[dict[str, object]]:
    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        raise CoordinatorContractError(
            "real passport adaptation requires the existing portal/spatial pyarrow environment"
        ) from error
    table = pq.read_table(path)
    return table.to_pylist()


def assert_safe_output(project_root: Path, output_root: Path) -> Path:
    project = project_root.resolve()
    expected_parent = (
        project / "scripts/manuscript/program_context_v2/coordinator/prepared"
    ).resolve()
    output = Path(output_root.absolute())
    try:
        relative = output.relative_to(expected_parent)
    except ValueError as error:
        raise CoordinatorContractError(
            f"coordinator output must be beneath {expected_parent}: {output}"
        ) from error
    if len(relative.parts) != 1 or relative.parts[0] in {"", ".", ".."}:
        raise CoordinatorContractError(
            "coordinator output requires one decision-ID child directory"
        )
    if output.exists() or output.is_symlink():
        raise CoordinatorContractError(
            f"refusing to overwrite coordinator preparation: {output}"
        )
    return output


def write_table(path: Path, rows, fields) -> None:
    atomic_write_tsv(path, rows, fields)


def build_source_alias_audit_rows(
    evidence_rows,
    graph_ids: set[str],
) -> list[dict[str, str]]:
    """Record exact graph/accession resolution and declared-ID-only limits."""
    dependency_usage: dict[str, dict[str, bool]] = {}
    for row in evidence_rows:
        discovery = set(str(row["discovery_sources"]).split(";"))
        evaluation = set(str(row["evaluation_sources"]).split(";"))
        if row[
            "source_dependence"
        ] == "independent" and "cross-alias equivalence" not in str(
            row["independence_boundary"]
        ):
            raise CoordinatorContractError(
                "independent row lacks the explicit alias-audit limitation: "
                f"{row['record_id']}"
            )
        for source_id in discovery | evaluation:
            usage = dependency_usage.setdefault(
                source_id,
                {"discovery": False, "evaluation": False, "independent": False},
            )
            usage["discovery"] |= source_id in discovery
            usage["evaluation"] |= source_id in evaluation
            usage["independent"] |= row["source_dependence"] == "independent"
    audit = []
    for source_id, usage in sorted(dependency_usage.items()):
        if source_id in graph_ids:
            basis = "Plan50 frozen source-node/edge graph"
            alias_status = "audited_exact_graph_id"
        elif re.fullmatch(r"dataset_(?:GSE|PXD)\d+", source_id):
            basis = "stable public accession encoded in the declared source ID"
            alias_status = "audited_stable_accession_identity"
        elif re.fullmatch(r"PLAN(?:13|20|30|40|50)", source_id):
            basis = "frozen workstream identity"
            alias_status = "audited_workstream_identity"
        else:
            basis = "literal coordinator contract ID only"
            alias_status = "not_audited_cross_alias_equivalence"
        audit.append(
            {
                "source_id": source_id,
                "canonical_source_id": source_id,
                "resolution_basis": basis,
                "alias_equivalence_status": alias_status,
                "used_as_discovery": str(usage["discovery"]).lower(),
                "used_as_evaluation": str(usage["evaluation"]).lower(),
                "used_in_independent_row": str(usage["independent"]).lower(),
                "independence_scope": (
                    "independent under declared source IDs only; cross-alias equivalence is not audited"
                    if usage["independent"]
                    else "no independent claim uses this ID"
                ),
            }
        )
    return audit


def prepare(
    project_root: Path,
    output_root: Path,
    coordinator: str,
    signing_date: str,
    signed_at_utc: str,
    decision_register_id: str,
) -> dict[str, object]:
    project = project_root.resolve()
    output = assert_safe_output(project, output_root)
    if not coordinator.strip() or not decision_register_id.strip():
        raise CoordinatorContractError(
            "coordinator and decision-register identity are required"
        )

    # Complete every scientific/source preflight before the first write.
    paths = RealPaths.build(project)
    handoffs = validate_all_handoffs(paths)

    integration = project / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
    fibrosis_bundle = fixed_live_bundle_root(project)
    metadata = integration / "metadata/unified_metadata.csv"
    qc = integration / "qc/sample_qc_report.csv"
    for source in (metadata, qc):
        if not source.is_file() or source.is_symlink():
            raise CoordinatorContractError(
                f"base source is missing or unsafe: {source}"
            )
    validate_fibrosis_candidate_bundle(
        fibrosis_bundle,
        project_root=project,
        verify_external_sources=True,
    )

    _, metadata_rows = read_tsv_flexible(metadata)
    _, qc_rows = read_tsv_flexible(qc)
    pooled_datasets = {"GSE126848", "GSE135251", "GSE130970", "GSE213621", "GSE162694"}
    cohort_rows = build_cohort_overview_rows(metadata_rows, qc_rows, pooled_datasets)

    _, hotspot_registry = read_tsv_flexible(
        paths.plan20_root / "program_registry_v2.tsv"
    )
    _, hotspot_figure = read_tsv_flexible(paths.plan20_root / "fig2_program_source.tsv")
    _, hotspot_semantic = read_tsv_flexible(
        paths.plan20_root / "program_registry_v2_semantic_adjudication.tsv"
    )
    _, hotspot_design = read_tsv_flexible(
        paths.plan20_root / "stage_dataset_design_audit.tsv"
    )
    _, hotspot_lodo = read_tsv_flexible(
        paths.plan20_root / "cohort_and_lodo_effects.tsv"
    )
    hotspot_rows = build_hotspot_rows(
        hotspot_registry, hotspot_figure, hotspot_semantic, hotspot_design, hotspot_lodo
    )
    # Rebuild Figure S2 only from Plan20-owned frozen source copies.  The
    # handoff evidence remains the staged source, while this explicit rebuild
    # proves its exact rows without consulting mutable live inputs.
    frozen_nmf_loadings, frozen_nmf_evidence = validate_nmf_continuous_bundle(paths)
    nmf_root = paths.plan20_root / "nmf_continuous_supplement"
    nmf_source_manifest = read_tsv_exact(
        nmf_root / "source_manifest.tsv", NMF_SOURCE_MANIFEST_FIELDS
    )
    nmf_sources = {row["source_id"]: row for row in nmf_source_manifest}

    def nmf_source_rows(source_id: str):
        _, rows = read_tsv_flexible(
            nmf_root / nmf_sources[source_id]["bundle_relative_path"]
        )
        return rows

    rebuilt_nmf_loadings, rebuilt_nmf_evidence = build_nmf_continuous_supplement(
        nmf_source_rows("k4_loadings"),
        nmf_source_rows("k6_loadings"),
        nmf_source_rows("k4_labels"),
        nmf_source_rows("k6_labels"),
        nmf_source_rows("three_seed_metrics"),
        nmf_source_rows("three_seed_stability"),
        nmf_sources["k4_loadings"]["bundle_sha256"],
        nmf_sources["k6_loadings"]["bundle_sha256"],
    )
    if [
        {field: str(row[field]) for field in NMF_CONTINUOUS_FIELDS}
        for row in rebuilt_nmf_loadings
    ] != frozen_nmf_loadings or [
        {field: str(row[field]) for field in SOURCE_EVIDENCE_FIELDS}
        for row in rebuilt_nmf_evidence
    ] != frozen_nmf_evidence:
        raise CoordinatorContractError(
            "Plan20 continuous NMF preparation rederivation drift"
        )
    hotspot_and_nmf_rows = [*hotspot_rows, *rebuilt_nmf_evidence]
    composition_results = read_tsv_exact(
        paths.plan20_root / "composition_sample_qc_v2.tsv", COMPOSITION_RESULT_FIELDS
    )
    composition_testability = read_tsv_exact(
        paths.plan20_root / "composition_sample_qc_testability.tsv",
        COMPOSITION_TESTABILITY_FIELDS,
    )
    composition_ready = one_row_tsv(
        paths.plan20_root / "COMPOSITION_QC_READY", COMPOSITION_READY_FIELDS
    )
    composition_rows = build_composition_rows(
        composition_results, composition_testability, composition_ready
    )

    _, phenotype = read_tsv_flexible(paths.plan30_root / "phenotype_registry.tsv")
    _, power = read_tsv_flexible(paths.plan30_root / "power_stratified_interface.tsv")
    _, terminal_rows = read_tsv_flexible(paths.plan30_root / "terminal_closure.tsv")
    if len(terminal_rows) != 1:
        raise CoordinatorContractError(
            "Plan30 terminal_closure.tsv must contain one row"
        )
    genetics_rows = build_genetics_rows(phenotype, power, terminal_rows[0])

    _, spatial_matrix = read_tsv_flexible(
        paths.plan13_root / "final_integration/figure4_program_matrix.tsv"
    )
    spatial_rows = build_spatial_rows(spatial_matrix)

    _, class_effects = read_tsv_flexible(paths.plan40_root / "class_effects.tsv")
    _, program_effects = read_tsv_flexible(paths.plan40_root / "program_effects.tsv")
    myojin_rows = build_myojin_rows(class_effects, program_effects)

    gene_rows = read_parquet_records(paths.plan50_root / "passport_gene_index.parquet")
    passport_evidence = read_parquet_records(
        paths.plan50_root / "passport_evidence_long.parquet"
    )
    passport_coverage = read_parquet_records(
        paths.plan50_root / "passport_coverage_long.parquet"
    )
    program_index = read_parquet_records(
        paths.plan50_root / "passport_program_index.parquet"
    )
    program_context = read_parquet_records(
        paths.plan50_root / "passport_program_context.parquet"
    )
    _, experiments = read_tsv_flexible(
        paths.plan50_root / "passport_next_experiment.tsv"
    )
    _, passport_source_nodes = read_tsv_flexible(
        paths.plan50_root / "passport_source_nodes.tsv"
    )
    _, passport_source_edges = read_tsv_flexible(
        paths.plan50_root / "passport_source_edges.tsv"
    )
    passport_ready = one_row_tsv(
        paths.plan50_root / "PASS06_VALIDATED",
        (
            "analysis_release_id",
            "status",
            "selection_sha256",
            "manifest_sha256",
            "validation_report_sha256",
            "n_genes",
            "n_evidence_rows",
            "n_programs",
            "n_program_context_rows",
            "automated_validation",
            "manual_acceptance",
            "handoff_allowed",
            "canonical_promotion_authorized",
            "scientific_call_recomputed",
            "validated_at_utc",
        ),
    )
    if (
        int(passport_ready["n_genes"]) != len(gene_rows)
        or int(passport_ready["n_evidence_rows"]) != len(passport_evidence)
        or int(passport_ready["n_programs"]) != len(program_index)
        or len({row["program_uid"] for row in program_index}) != len(program_index)
        or int(passport_ready["n_program_context_rows"]) != len(program_context)
    ):
        raise CoordinatorContractError(
            "Plan50 PASS06 row counts disagree with frozen Parquet tables"
        )
    passport_rows = build_passport_rows(
        gene_rows,
        passport_evidence,
        passport_coverage,
        experiments,
        program_context,
        passport_source_nodes,
        passport_source_edges,
    )
    all_evidence_rows = [
        *cohort_rows,
        *composition_rows,
        *hotspot_and_nmf_rows,
        *genetics_rows,
        *spatial_rows,
        *passport_rows,
        *myojin_rows,
    ]
    graph_ids = {
        *(str(row["source_node_id"]) for row in passport_source_nodes),
        *(str(row["source_edge_id"]) for row in passport_source_edges),
    }
    source_alias_rows = build_source_alias_audit_rows(all_evidence_rows, graph_ids)

    scopes = protected_scope_rows(project)
    baseline = protected_baseline_rows(project, scopes)
    blueprint = locked_release_blueprint()
    runtime_environment = core_runtime_environment_rows()
    package_export = core_package_export_rows()
    sys_path_export = core_sys_path_rows()
    import_resolution = core_import_resolution_rows()
    recursive_producers = recursive_release_producer_rows(project)

    # Owner manifests are upstream-owned.  Closure records their exact hashes.
    closure = closure_rows(project, handoffs, coordinator, signing_date)

    # All data are now validated in memory. Materialize one no-overwrite coordinator bundle.
    output.mkdir(parents=True)
    evidence_root = output / "evidence_rows"
    evidence_root.mkdir()
    environment_root = output / "environment"
    environment_root.mkdir()
    evidence_products = {
        "cohort_overview": (evidence_root / "cohort_overview.tsv", cohort_rows),
        "sample_composition": (
            evidence_root / "sample_composition.tsv",
            composition_rows,
        ),
        "hotspot_programs": (
            evidence_root / "hotspot_programs.tsv",
            hotspot_and_nmf_rows,
        ),
        "genetics_evidence": (evidence_root / "genetics_evidence.tsv", genetics_rows),
        "spatial_evidence": (evidence_root / "spatial_evidence.tsv", spatial_rows),
        "passport_evidence": (evidence_root / "passport_evidence.tsv", passport_rows),
        "myojin_supplement": (evidence_root / "myojin_supplement.tsv", myojin_rows),
    }
    for path, rows in evidence_products.values():
        write_table(path, rows, SOURCE_EVIDENCE_FIELDS)
    source_alias_audit_path = output / "source_alias_audit.tsv"
    write_table(source_alias_audit_path, source_alias_rows, SOURCE_ALIAS_AUDIT_FIELDS)
    environment_products = {
        "core_coordinator_runtime": (
            environment_root / "core_coordinator_runtime.tsv",
            runtime_environment,
            CORE_RUNTIME_FIELDS,
        ),
        "core_coordinator_packages": (
            environment_root / "core_coordinator_packages.tsv",
            package_export,
            CORE_PACKAGE_FIELDS,
        ),
        "core_coordinator_sys_path": (
            environment_root / "core_coordinator_sys_path.tsv",
            sys_path_export,
            CORE_SYSPATH_FIELDS,
        ),
        "core_coordinator_imports": (
            environment_root / "core_coordinator_imports.tsv",
            import_resolution,
            CORE_IMPORT_FIELDS,
        ),
        "recursive_release_producers": (
            environment_root / "recursive_release_producers.tsv",
            recursive_producers,
            CORE_PRODUCER_FIELDS,
        ),
    }
    for path, rows, fields in environment_products.values():
        write_table(path, rows, fields)

    producer_manifest_path = environment_products["recursive_release_producers"][0]
    producer_manifest_sha256 = sha256_file(producer_manifest_path)
    producer_by_path = {str(row["repository_path"]): row for row in recursive_producers}
    if set(evidence_products) != EXPECTED_ADAPTER_IDS:
        raise CoordinatorContractError(
            "prepared adapter universe differs from the frozen adapter contract"
        )
    adapter_provenance_rows = []
    for artifact_id, (adapter_path, _) in evidence_products.items():
        spec = ADAPTER_CONTRACT[artifact_id]
        if adapter_path.name != spec["output_relative_path"]:
            raise CoordinatorContractError(
                f"adapter output name differs from frozen contract: {artifact_id}"
            )
        producer_paths = tuple(str(path) for path in spec["required_producers"])
        if any(path not in producer_by_path for path in producer_paths):
            raise CoordinatorContractError(
                f"adapter producer is absent from recursive manifest: {artifact_id}"
            )
        producer_bindings = [
            {
                "producer_id": producer_by_path[path]["producer_id"],
                "repository_path": path,
                "sha256": producer_by_path[path]["sha256"],
                "bytes": int(producer_by_path[path]["bytes"]),
            }
            for path in producer_paths
        ]
        inputs = []
        for source in adapter_input_paths(project, artifact_id):
            if not source.is_file() or source.is_symlink():
                raise CoordinatorContractError(
                    f"adapter source is missing or unsafe: {artifact_id}: {source}"
                )
            inputs.append(
                {
                    "source_path": project_relative(project, source),
                    "sha256": sha256_file(source),
                    "bytes": source.stat().st_size,
                }
            )
        adapter_provenance_rows.append(
            {
                "contract_version": ADAPTER_PROVENANCE_VERSION,
                "adapter_contract_sha256": ADAPTER_CONTRACT_SHA256,
                "artifact_id": artifact_id,
                # REL01 snapshots each selected adapter directly beneath
                # inputs/BASE. Bind the staged basename, not this preparation
                # bundle's private evidence_rows/ layout.
                "output_relative_path": adapter_path.name,
                "output_sha256": sha256_file(adapter_path),
                "output_bytes": adapter_path.stat().st_size,
                "producer_manifest_relative_path": "inputs/BASE/recursive_release_producers.tsv",
                "producer_manifest_sha256": producer_manifest_sha256,
                "producer_bindings_json": json.dumps(
                    producer_bindings, sort_keys=True, separators=(",", ":")
                ),
                "input_bindings_json": json.dumps(
                    inputs, sort_keys=True, separators=(",", ":")
                ),
                "canonical_promotion_authorized": "false",
            }
        )
    adapter_provenance_path = output / "adapter_provenance.tsv"
    write_table(
        adapter_provenance_path,
        adapter_provenance_rows,
        ADAPTER_PROVENANCE_FIELDS,
    )

    blueprint_path = output / "release_blueprint.json"
    atomic_write_json(blueprint_path, blueprint)
    closure_path = output / "workstream_closure.tsv"
    write_table(closure_path, closure, CLOSURE_FIELDS)
    scopes_path = output / "protected_scopes.tsv"
    baseline_path = output / "protected_release_baseline.tsv"
    write_table(scopes_path, scopes, PROTECTED_SCOPE_FIELDS)
    write_table(baseline_path, baseline, PROTECTED_BASELINE_FIELDS)

    selection = []
    for relative in FIBROSIS_BUNDLE_FILES:
        primary = relative == FIBROSIS_PRIMARY_RELATIVE
        selection.append(
            selection_row(
                project,
                fibrosis_bundle / relative,
                fibrosis_artifact_id(relative),
                fibrosis_snapshot_relpath(relative),
                fibrosis_artifact_role(relative),
                coordinator,
                signing_date,
                decision_register_id,
                (
                    "cross-sectional true-Kleiner stage-associated DEG counts"
                    if primary
                    else "validation/provenance artifact for the sealed true-Kleiner candidate bundle"
                ),
                (
                    "longitudinal progression, stale transition summary, or stand-alone unsealed coefficient table"
                    if primary
                    else "scientific evidence in isolation, canonical promotion, or substitution outside the sealed bundle"
                ),
            )
        )
    selection.extend(
        [
            selection_row(
                project,
                metadata,
                "cohort_metadata_raw",
                "inputs/BASE/unified_metadata.csv",
                "cohort_overview_raw",
                coordinator,
                signing_date,
                decision_register_id,
                "checksummed cohort-census provenance",
                "direct figure read or unadjudicated sample count",
            ),
            selection_row(
                project,
                qc,
                "cohort_qc_raw",
                "inputs/BASE/sample_qc_report.csv",
                "cohort_overview_raw",
                coordinator,
                signing_date,
                decision_register_id,
                "checksummed QC-census provenance",
                "direct figure read or unadjudicated sample count",
            ),
            selection_row(
                project,
                blueprint_path,
                "release_blueprint",
                "inputs/BASE/release_blueprint.json",
                "release_blueprint",
                coordinator,
                signing_date,
                decision_register_id,
                "locked five-figure Cell Genomics Resource branch",
                "sixth main figure, main-figure Myojin, or canonical promotion",
            ),
        ]
    )
    adapter_specs = {
        "cohort_overview": ("inputs/BASE/cohort_overview.tsv", "cohort_overview"),
        "sample_composition": (
            "inputs/BASE/sample_composition.tsv",
            "figure_source_adapter",
        ),
        "hotspot_programs": (
            "inputs/BASE/hotspot_programs.tsv",
            "figure_source_adapter",
        ),
        "genetics_evidence": (
            "inputs/BASE/genetics_evidence.tsv",
            "figure_source_adapter",
        ),
        "spatial_evidence": (
            "inputs/BASE/spatial_evidence.tsv",
            "figure_source_adapter",
        ),
        "passport_evidence": (
            "inputs/BASE/passport_evidence.tsv",
            "figure_source_adapter",
        ),
        "myojin_supplement": (
            "inputs/BASE/myojin_supplement.tsv",
            "figure_source_adapter",
        ),
    }
    for artifact_id, (snapshot, role) in adapter_specs.items():
        source = evidence_products[artifact_id][0]
        selection.append(
            selection_row(
                project,
                source,
                artifact_id,
                snapshot,
                role,
                coordinator,
                signing_date,
                decision_register_id,
                "presentation-only adapter preserving frozen upstream calls",
                "scientific recomputation, outcome-informed selection, score, rank, or canonical write",
            )
        )
    selection.append(
        selection_row(
            project,
            adapter_provenance_path,
            "adapter_provenance",
            "inputs/BASE/adapter_provenance.tsv",
            "adapter_provenance_manifest",
            coordinator,
            signing_date,
            decision_register_id,
            "exact adapter-output, input, and producer-hash bindings",
            "unbound adapter output, mutable producer substitution, or promotion authorization",
        )
    )
    selection.append(
        selection_row(
            project,
            source_alias_audit_path,
            "source_alias_audit",
            "inputs/BASE/source_alias_audit.tsv",
            "source_alias_audit",
            coordinator,
            signing_date,
            decision_register_id,
            "canonical graph/accession resolution where available and an explicit declared-ID-only limitation otherwise",
            "unqualified independence beyond declared IDs or hidden alias equivalence",
        )
    )
    environment_specs = {
        "core_coordinator_runtime": (
            "inputs/BASE/core_coordinator_runtime.tsv",
            "environment_manifest",
        ),
        "core_coordinator_packages": (
            "inputs/BASE/core_coordinator_packages.tsv",
            "environment_manifest",
        ),
        "core_coordinator_sys_path": (
            "inputs/BASE/core_coordinator_sys_path.tsv",
            "environment_manifest",
        ),
        "core_coordinator_imports": (
            "inputs/BASE/core_coordinator_imports.tsv",
            "environment_manifest",
        ),
        "recursive_release_producers": (
            "inputs/BASE/recursive_release_producers.tsv",
            "producer_manifest",
        ),
    }
    for artifact_id, (snapshot, role) in environment_specs.items():
        source = environment_products[artifact_id][0]
        selection.append(
            selection_row(
                project,
                source,
                artifact_id,
                snapshot,
                role,
                coordinator,
                signing_date,
                decision_register_id,
                "frozen runtime/package/producer provenance for candidate reproducibility",
                "scientific evidence, biological result, or promotion authorization",
            )
        )
    for workstream in WORKSTREAM_ORDER:
        signature = handoffs[workstream].signature
        selection.append(
            selection_row(
                project,
                signature,
                f"{workstream.lower()}_owner_attestation",
                f"inputs/BASE/owner_attestations/{workstream}.json",
                "owner_attestation",
                coordinator,
                signing_date,
                decision_register_id,
                "detached workstream-owner attestation",
                "scientific evidence or promotion authorization",
            )
        )
    selection_path = output / "base_input_selection.tsv"
    write_table(selection_path, selection, BASE_SELECTION_FIELDS)

    closure_sig = coordinator_signature(
        closure_path.relative_to(project),
        sha256_file(closure_path),
        coordinator,
        signed_at_utc,
        decision_register_id,
        "workstream_closure",
    )
    selection_sig = coordinator_signature(
        selection_path.relative_to(project),
        sha256_file(selection_path),
        coordinator,
        signed_at_utc,
        decision_register_id,
        "base_input_selection",
    )
    blueprint_sig = coordinator_signature(
        blueprint_path.relative_to(project),
        sha256_file(blueprint_path),
        coordinator,
        signed_at_utc,
        decision_register_id,
        "release_blueprint",
    )
    atomic_write_json(output / "workstream_closure.signature.json", closure_sig)
    atomic_write_json(output / "base_input_selection.signature.json", selection_sig)
    atomic_write_json(output / "release_blueprint.signature.json", blueprint_sig)

    manifested = []
    for path in sorted(item for item in output.rglob("*") if item.is_file()):
        if path.name in {
            "coordinator_preparation_manifest.tsv",
            "COORDINATOR_PREPARATION_READY",
        }:
            continue
        manifested.append(
            {
                "artifact_id": path.relative_to(output).as_posix().replace("/", "__"),
                "relative_path": path.relative_to(output).as_posix(),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "artifact_role": "coordinator_preparation",
                "canonical_promotion_authorized": "false",
            }
        )
    manifest_path = output / "coordinator_preparation_manifest.tsv"
    write_table(manifest_path, manifested, PREPARATION_MANIFEST_FIELDS)
    ready = {
        "contract_version": PREPARATION_CONTRACT_VERSION,
        "candidate_id": CANDIDATE_ID,
        "status": "COORDINATOR_PREPARATION_READY_NOT_EXECUTED_NOT_PROMOTED",
        "decision_register_id": decision_register_id,
        "closure_sha256": sha256_file(closure_path),
        "base_selection_sha256": sha256_file(selection_path),
        "blueprint_sha256": sha256_file(blueprint_path),
        "protected_baseline_sha256": sha256_file(baseline_path),
        "preparation_manifest_sha256": sha256_file(manifest_path),
        "rel00_05_executed": False,
        "candidate_root_created": False,
        "canonical_promotion_authorized": False,
    }
    atomic_write_json(output / "COORDINATOR_PREPARATION_READY", ready)
    return ready


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--coordinator", required=True)
    parser.add_argument("--signing-date", required=True)
    parser.add_argument("--signed-at-utc", required=True)
    parser.add_argument("--decision-register-id", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = prepare(
            args.project_root,
            args.output_root,
            args.coordinator,
            args.signing_date,
            args.signed_at_utc,
            args.decision_register_id,
        )
    except ReleaseContractError as error:
        print(f"PLAN60_COORDINATOR_BLOCKED: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
