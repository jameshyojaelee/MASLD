#!/usr/bin/env python3
"""Adversarial, write-isolated tests for the real Plan60 coordinator contract."""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


COORDINATOR_ROOT = Path(__file__).resolve().parents[1]
PROGRAM_ROOT = COORDINATOR_ROOT.parent
PROJECT_ROOT = PROGRAM_ROOT.parents[2]
TEST_FIXTURE_ROOT = PROGRAM_ROOT / "tests"
for path in (COORDINATOR_ROOT, PROGRAM_ROOT, TEST_FIXTURE_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from build_scientific_validation_registry import (  # noqa: E402
    bh_adjust,
    scientific_producer_binding,
)
from compare_clean_rebuilds import (  # noqa: E402
    CANDIDATE_RELATIVE,
    compare_clean_rebuilds,
    derive_clean_rebuild_rows,
    revalidate_clean_rebuild_report,
)
from coordinator_contract import (  # noqa: E402
    BASE_SELECTION_FIELDS,
    CANDIDATE_ID,
    COMPOSITION_READY_FIELDS,
    COMPOSITION_RESULT_FIELDS,
    COMPOSITION_TESTABILITY_FIELDS,
    CORE_PRODUCER_FIELDS,
    EXPECTED_COHORT_COUNTS,
    FIXED_HERO_GENES,
    HANDOFF_SIGNATURE_VERSION,
    Handoff,
    PLAN50_ANALYSIS_RELEASE_ID,
    RealPaths,
    SOURCE_EVIDENCE_FIELDS,
    CoordinatorContractError,
    build_cohort_overview_rows,
    build_composition_rows,
    build_hotspot_rows,
    build_myojin_rows,
    build_passport_rows,
    build_spatial_rows,
    closure_rows,
    core_package_export_rows,
    core_import_resolution_rows,
    core_runtime_environment_rows,
    core_sys_path_rows,
    locked_release_blueprint,
    one_row_tsv,
    read_tsv_flexible,
    read_tsv_exact,
    recursive_release_producer_rows,
    sha256_file,
    source_row,
    validate_owner_signature,
    validate_hotspot_registry_projection,
    validate_composition_contract,
    validate_passport_negative_semantics,
    validate_source_rows,
    validate_nmf_continuous_bundle,
)
from prepare_real_release import assert_safe_output, prepare  # noqa: E402
from rel03_render_candidate_panels import build_pdf  # noqa: E402
from release_products import EXPECTED_FIGURES, EXPECTED_FIGURE_TITLES  # noqa: E402
from release_common import (  # noqa: E402
    PLAN50_FROZEN_PRODUCER_PATHS,
    PLAN50_MANIFEST_PRODUCER_PATHS,
    ReleaseContractError,
    SCIENTIFIC_REPORT_PRODUCER_PATHS,
    atomic_write_tsv,
)
from release_products import TRANSITION_PRODUCT_FIELDS  # noqa: E402
from fixture_factory import SyntheticProject  # noqa: E402
from rel01_snapshot_candidate import (  # noqa: E402
    producer_script_inventory,
    snapshot_candidate,
    validate_plan50_frozen_producer_binding,
    validate_prepared_adapter_provenance,
)
from rel02_stage_base_inputs import stage_base_inputs  # noqa: E402
from rel02_build_candidate_tables import build_candidate_tables  # noqa: E402
from rel03_render_candidate_panels import render_candidate_panels  # noqa: E402
from rel04_build_candidate_manuscript import build_candidate_manuscript  # noqa: E402
from rel05_validate_candidate import validate_snapshot_base  # noqa: E402
import publish_plan60_workstream_handoff as handoff_publisher  # noqa: E402
import prepare_real_release as preparation_module  # noqa: E402
import materialize_clean_rebuild_project as materializer  # noqa: E402
import coordinator_contract as coordinator_module  # noqa: E402
import rel01_snapshot_candidate as rel01_module  # noqa: E402
import rel05_validate_candidate as rel05_module  # noqa: E402
import build_scientific_validation_registry as scientific_registry_module  # noqa: E402
from adapter_contract import (  # noqa: E402
    ADAPTER_CONTRACT,
    ADAPTER_CONTRACT_SHA256,
    NMF_FREEZER_PRODUCER_PATH,
    require_exact_adapter_binding_paths,
)
from fibrosis_candidate_contract import (  # noqa: E402
    FIBROSIS_ACCEPTED_PRIMARY_RESULT_SHA256,
    FIBROSIS_ACCEPTED_READY_SHA256,
    FIBROSIS_ACCEPTED_VALIDATED_MANIFEST_SHA256,
    FIBROSIS_BUNDLE_FILES,
    FIBROSIS_PRIMARY_ARTIFACT_ID,
    FIBROSIS_VALIDATION_ROLE,
    fibrosis_artifact_id,
)


def hotspot_fixture():
    registries = []
    figures = []
    semantics = []
    designs = []
    lodo = []
    for index in range(117):
        uid = f"hotspot_hepatocytes_{index:03d}"
        robust = index in {7, 19}
        figure = {
            "release_id": CANDIDATE_ID,
            "program_uid": uid,
            "membership_sha256": f"{index:064x}",
            "cell_type": "hepatocytes",
            "module": str(index + 1),
            "module_name": f"Program {index + 1}",
            "module_hallmark": "Hallmark",
            "module_top_pathway": "Pathway",
            "beta": str(index / 100 + 0.01),
            "se": "0.125",
            "statistic": "1.5",
            "pvalue": "0.001" if robust else "0.5",
            "qvalue": "0.01" if robust else "0.8",
            "hc3_se": "0.13",
            "hc3_statistic": "1.4",
            "hc3_pvalue": "0.16",
            "hc3_qvalue": "0.3",
            "max_leverage": "0.2",
            "direction": "positive",
            "stability_score": "0.75",
            "stability_mean": "0.74",
            "stability_median": "0.76",
            "primary_selected": "TRUE" if robust else "FALSE",
            "hc3_supported": "TRUE" if robust else "FALSE",
            "scoring_direction_agree": "TRUE",
            "selected_hc3_fragile": "FALSE",
            "robust_display": "TRUE" if robust else "FALSE",
            "source_stability": "multi_dataset_supported",
            "source_stability_reason": "fixture",
            "n_donors": "45",
            "n_datasets": "4",
            "n_healthy": "10",
            "n_steatosis": "9",
            "n_steatohepatitis": "26",
        }
        figures.append(figure)
        registries.append(
            {
                "release_id": figure["release_id"],
                "program_uid": uid,
                "membership_sha256": figure["membership_sha256"],
                "cell_type": figure["cell_type"],
                "module": figure["module"],
                "module_name": figure["module_name"],
                "module_hallmark": figure["module_hallmark"],
                "module_top_pathway": figure["module_top_pathway"],
                "primary_beta": figure["beta"],
                "primary_se": figure["se"],
                "primary_statistic": figure["statistic"],
                "primary_pvalue": figure["pvalue"],
                "primary_qvalue": figure["qvalue"],
                "primary_hc3_se": figure["hc3_se"],
                "primary_hc3_statistic": figure["hc3_statistic"],
                "primary_hc3_pvalue": figure["hc3_pvalue"],
                "primary_hc3_qvalue": figure["hc3_qvalue"],
                "primary_max_leverage": figure["max_leverage"],
                "primary_direction": figure["direction"],
                "stability_score": figure["stability_score"],
                "stability_mean": figure["stability_mean"],
                "stability_median": figure["stability_median"],
                "primary_selected": figure["primary_selected"],
                "hc3_supported": figure["hc3_supported"],
                "scoring_direction_agree": figure["scoring_direction_agree"],
                "selected_hc3_fragile": figure["selected_hc3_fragile"],
                "robust_display": figure["robust_display"],
                "source_stability": figure["source_stability"],
                "source_stability_reason": figure["source_stability_reason"],
                "primary_n_donors": figure["n_donors"],
                "primary_n_datasets": figure["n_datasets"],
                "primary_n_stage0": figure["n_healthy"],
                "primary_n_stage1": figure["n_steatosis"],
                "primary_n_stage2": figure["n_steatohepatitis"],
            }
        )
        semantics.append(
            {
                "program_uid": uid,
                "module": str(index + 1),
                "module_name": f"Program {index + 1}",
                "cell_type": "hepatocytes",
                "robust_display": "TRUE" if robust else "FALSE",
                "adjudicated_state": (
                    "supported_internal_stage_association"
                    if robust
                    else "indeterminate_nonconfirmatory"
                ),
                "tested_negative_authorized": "FALSE",
                "legacy_tested_negative": "TRUE",
            }
        )
        if robust:
            for stage in ("0", "1", "2"):
                designs.append(
                    {
                        "program_uid": uid,
                        "dataset": "GSE244832",
                        "stage_ordinal": stage,
                        "stage_dataset_n_donors": "3",
                    }
                )
            for stage in ("1", "2"):
                designs.append(
                    {
                        "program_uid": uid,
                        "dataset": "GSE202379",
                        "stage_ordinal": stage,
                        "stage_dataset_n_donors": "4",
                    }
                )
            lodo.append(
                {
                    "program_uid": uid,
                    "analysis_type": "lodo",
                    "held_out_dataset": "GSE244832",
                    "estimable": "TRUE",
                    "direction": "positive",
                }
            )
    return registries, figures, semantics, designs, lodo


def passport_fixture():
    genes = []
    evidence = []
    coverage = []
    experiments = []
    source_nodes = [
        {
            "source_node_id": "node:release:fixture",
            "source_release_id": "fixture-source-v1",
        },
        {
            "source_node_id": "node:dataset:fixture",
            "source_release_id": "fixture-source-v1",
        },
    ]
    source_edges = []
    for index, symbol in enumerate(FIXED_HERO_GENES, 1):
        passport = f"passport-{index}"
        ensembl = f"ENSG{index:011d}"
        genes.append(
            {
                "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
                "passport_id": passport,
                "ensembl_id": ensembl,
                "symbol": symbol,
                "symbol_collision": False,
                "primary_evidence_class": "unresolved",
            }
        )
        evidence.append(
            {
                "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
                "evidence_result_id": f"ev-{index}",
                "passport_id": passport,
                "ensembl_id": ensembl,
                "symbol": symbol,
                "evidence_domain": "genetics",
                "assay": "COLOC",
                "dataset_id": "frozen",
                "biological_unit": "GWAS locus",
                "contrast_or_exposure": "trait association",
                "effect_unit": "PP.H4",
                "estimate": 0.8,
                "p_value": None,
                "q_value": None,
                "testability_state": "testable",
                "call_state": "supported",
                "negative_call_rule_id": None,
                "negative_decision_boundary": None,
                "negative_margin": None,
                "negative_call_passed": None,
                "n_biological_units": 1,
                "provenance_state": "independent",
                "source_node_id": "node:release:fixture",
                "source_release_id": "fixture-source-v1",
            }
        )
        coverage.append(
            {
                "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
                "coverage_result_id": f"cov-{index}",
                "passport_id": passport,
                "ensembl_id": ensembl,
                "symbol": symbol,
                "evidence_domain": "proteomics",
                "assay": "DIA-MS",
                "dataset_id": "frozen",
                "call_state": "untestable",
                "testability_state": "not_detected",
                "n_biological_units": None,
                "coverage_denominator": "frozen protein panel",
                "source_release_id": "fixture-source-v1",
            }
        )
        experiments.append(
            {
                "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
                "passport_id": passport,
                "ensembl_id": ensembl,
                "biological_model": "primary-like hepatocyte",
                "context": "lipid stress",
                "perturbation": "allele-aware perturbation",
                "primary_readout": "target RNA",
                "falsifying_outcome": "no allele-specific response",
            }
        )
        source_edges.append(
            {
                "source_edge_id": f"edge:tested:{index}",
                "from_node_id": "node:release:fixture",
                "to_node_id": "node:dataset:fixture",
                "edge_type": "tested_by",
                "evidence_result_id": f"ev-{index}",
            }
        )
    contexts = [
        {
            "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
            "program_context_id": "ctx-1",
            "call_state": "indeterminate",
            "testability_state": "testable",
            "negative_call_rule_id": None,
            "negative_decision_boundary": None,
            "negative_margin": None,
            "negative_call_passed": None,
            "gene_call_expansion_authorized": False,
        }
    ]
    return genes, evidence, coverage, experiments, contexts, source_nodes, source_edges


class CoordinatorContractTests(unittest.TestCase):
    def test_clean_rebuild_materializer_treats_authority_diff_as_provenance(self):
        spec = {
            "repository": {
                "repository_commit": "1" * 40,
                "tracked_diff_sha256": "2" * 64,
            }
        }
        rebuild_spec = {"repository_commit": "1" * 40}
        observed = materializer.authority_repository_provenance(spec, rebuild_spec)
        self.assertEqual(observed["authority_repository_commit"], "1" * 40)
        self.assertEqual(observed["authority_tracked_diff_sha256"], "2" * 64)

    def test_sparse_materializer_does_not_copy_unrelated_authority_diff(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            authority = root / "authority"
            destination = root / "retained"
            authority.mkdir()
            subprocess.run(
                ["git", "init", str(authority)], check=True, capture_output=True
            )
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(authority),
                    "config",
                    "user.email",
                    "test@example.org",
                ],
                check=True,
            )
            subprocess.run(
                ["git", "-C", str(authority), "config", "user.name", "Plan60 test"],
                check=True,
            )
            tracked = authority / "tracked.txt"
            tracked.write_text("committed\n", encoding="utf-8")
            subprocess.run(
                ["git", "-C", str(authority), "add", "tracked.txt"], check=True
            )
            subprocess.run(
                ["git", "-C", str(authority), "commit", "-m", "fixture"],
                check=True,
                capture_output=True,
            )
            commit = subprocess.run(
                ["git", "-C", str(authority), "rev-parse", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            tracked.write_text("unrelated authority modification\n", encoding="utf-8")
            materializer.create_sparse_worktree(
                authority, destination, commit, {"tracked.txt"}
            )
            self.assertEqual(
                (destination / "tracked.txt").read_text(encoding="utf-8"),
                "committed\n",
            )
            state = materializer.validate_isolated_repository_state(
                destination, commit, set()
            )
            self.assertEqual(state["isolated_tracked_overlay_paths"], [])

    def test_isolated_rebuild_rejects_dirty_gitlink(self):
        patch = (
            b"diff --git a/GWAS/finemapping/src/chrombpnet_variant_scorer "
            b"b/GWAS/finemapping/src/chrombpnet_variant_scorer\n"
            b"-Subproject commit 0e1e34199e63112aa618748bb79a206fc491300a\n"
            b"+Subproject commit 0e1e34199e63112aa618748bb79a206fc491300a-dirty\n"
        )

        def fake_git(_project, *arguments, input_bytes=None):
            del input_bytes
            if arguments[:2] == ("rev-parse", "HEAD"):
                return b"1" * 40 + b"\n"
            if arguments[:3] == ("diff", "--binary", "HEAD"):
                return patch
            if arguments[:3] == ("diff", "--name-only", "-z"):
                return b"GWAS/finemapping/src/chrombpnet_variant_scorer\0"
            if arguments[:2] == ("ls-files", "--others"):
                return b""
            raise AssertionError(arguments)

        with mock.patch.object(materializer, "git_output", side_effect=fake_git):
            with self.assertRaisesRegex(
                materializer.MaterializationError, "dirty gitlink"
            ):
                materializer.validate_isolated_repository_state(
                    Path("/retained"), "1" * 40, set()
                )

    def test_materializer_modified_frozen_source_blocks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.tsv"
            destination = root / "destination.tsv"
            source.write_text("frozen\n", encoding="utf-8")
            digest = sha256_file(source)
            size = source.stat().st_size
            source.write_text("modified\n", encoding="utf-8")
            with self.assertRaisesRegex(
                materializer.MaterializationError, "frozen source hash/byte drift"
            ):
                materializer.copy_frozen_file(source, destination, digest, size)

    def test_rebuild_source_modified_frozen_producer_blocks(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = SyntheticProject(Path(temporary))
            snapshot_candidate(
                project.root,
                project.closure_path,
                project.protected_scopes_path,
                project.protected_baseline_path,
                fixture_mode=True,
            )
            producer_rows = read_tsv_exact(
                project.candidate_root / "manifests/producer_script_manifest.tsv",
                (
                    "producer_id",
                    "repository_path",
                    "sha256",
                    "bytes",
                    "snapshot_path",
                ),
            )
            producer = project.candidate_root / producer_rows[0]["snapshot_path"]
            producer.write_bytes(producer.read_bytes() + b"\n# tampered\n")
            with self.assertRaisesRegex(
                ReleaseContractError, "frozen producer copy drift"
            ):
                validate_snapshot_base(project.root, fixture_mode=True)

    def test_rel01_snapshots_and_rehydrates_external_nmf_freezer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = SyntheticProject(root / "authority")
            snapshot_candidate(
                project.root,
                project.closure_path,
                project.protected_scopes_path,
                project.protected_baseline_path,
                fixture_mode=True,
            )
            producer_rows = read_tsv_exact(
                project.candidate_root / "manifests/producer_script_manifest.tsv",
                (
                    "producer_id",
                    "repository_path",
                    "sha256",
                    "bytes",
                    "snapshot_path",
                ),
            )
            external = next(
                row
                for row in producer_rows
                if row["repository_path"] == NMF_FREEZER_PRODUCER_PATH
            )
            frozen_copy = project.candidate_root / external["snapshot_path"]
            self.assertEqual(sha256_file(frozen_copy), external["sha256"])
            destination = root / "retained_project"
            restored_count = materializer.restore_frozen_producers(
                project.candidate_root, destination
            )
            self.assertEqual(restored_count, len(producer_rows))
            restored = destination / NMF_FREEZER_PRODUCER_PATH
            self.assertEqual(sha256_file(restored), external["sha256"])
            self.assertEqual(restored.stat().st_size, int(external["bytes"]))

    def test_adapter_contract_requires_exact_input_and_producer_sets(self):
        self.assertEqual(len(ADAPTER_CONTRACT_SHA256), 64)
        for artifact_id, spec in ADAPTER_CONTRACT.items():
            inputs = tuple(spec["required_inputs"])
            producers = tuple(spec["required_producers"])
            require_exact_adapter_binding_paths(artifact_id, inputs, producers)
            with self.assertRaisesRegex(
                ReleaseContractError, "input set differs from frozen contract"
            ):
                require_exact_adapter_binding_paths(artifact_id, inputs[:-1], producers)
            with self.assertRaisesRegex(
                ReleaseContractError, "producer set differs from frozen contract"
            ):
                require_exact_adapter_binding_paths(
                    artifact_id,
                    inputs,
                    (*producers, "scripts/tampered_adapter.py"),
                )

    def test_adapter_contract_rejects_duplicate_binding_tamper(self):
        artifact_id = "cohort_overview"
        spec = ADAPTER_CONTRACT[artifact_id]
        inputs = tuple(spec["required_inputs"])
        producers = tuple(spec["required_producers"])
        with self.assertRaisesRegex(ReleaseContractError, "duplicate input"):
            require_exact_adapter_binding_paths(
                artifact_id, (*inputs, inputs[0]), producers
            )
        with self.assertRaisesRegex(ReleaseContractError, "duplicate producer"):
            require_exact_adapter_binding_paths(
                artifact_id, inputs, (*producers, producers[0])
            )

    def test_real_nmf_freezer_rederives_and_renders_both_panels(self):
        generator_path = (
            PROJECT_ROOT
            / "Analysis/SingleCell/scripts/hotspot_modules"
            / "519_freeze_nmf_continuous_supplement.py"
        )
        specification = importlib.util.spec_from_file_location(
            "plan20_nmf_freezer_test", generator_path
        )
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            shutil.copytree(
                PROGRAM_ROOT,
                project / "scripts/manuscript/program_context_v2",
            )
            source_paths = [
                "RNA-seq/results/subtypes/nmf_assignments_k4_pre_k6restore.csv",
                "RNA-seq/results/subtypes/nmf_assignments.csv",
                "RNA-seq/results/subtypes/program_labels_k4_pre_k6restore.csv",
                "RNA-seq/results/subtypes/program_labels.csv",
                "RNA-seq/results/subtypes/nmf_3seed_canonical/per_seed_metric_summary.csv",
                "RNA-seq/results/subtypes/nmf_3seed_canonical/cross_seed_stability_summary.csv",
                "Analysis/SingleCell/scripts/hotspot_modules/519_freeze_nmf_continuous_supplement.py",
                *PLAN50_FROZEN_PRODUCER_PATHS,
            ]
            for relative in source_paths:
                destination = project / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(PROJECT_ROOT / relative, destination)
            result = module.freeze(project, "2026-08-08")
            self.assertEqual(result["n_loading_rows"], 11040)
            paths = RealPaths.build(project)
            loadings, evidence = validate_nmf_continuous_bundle(paths)
            self.assertEqual((len(loadings), len(evidence)), (11040, 18))
            s2a = [row for row in evidence if row["panel_id"] == "S2A"]
            s2b = [row for row in evidence if row["panel_id"] == "S2B"]
            self.assertEqual((len(s2a), len(s2b)), (10, 8))
            self.assertTrue(all(row["plot_role"] == "mark" for row in s2a + s2b))
            self.assertIn(
                b"median=", build_pdf("Continuous k4/k6 NMF loading distributions", s2a)
            )
            self.assertIn(
                b"mean Hungaria",
                build_pdf(
                    "Three-seed factor stability and hard-partition boundary", s2b
                ),
            )

            nmf_root = paths.plan20_root / "nmf_continuous_supplement"
            historical_manifest = nmf_root / "producer_manifest.tsv"
            original_historical_bytes = historical_manifest.read_bytes()
            historical_manifest.write_bytes(original_historical_bytes + b"\n")
            with self.assertRaisesRegex(
                CoordinatorContractError,
                "producer_manifest_sha256 drifted after READY",
            ):
                validate_nmf_continuous_bundle(paths)
            historical_manifest.write_bytes(original_historical_bytes)

            current_rows = recursive_release_producer_rows(project)
            tampered_current_rows = copy.deepcopy(current_rows)
            freezer_relative = (
                "Analysis/SingleCell/scripts/hotspot_modules/"
                "519_freeze_nmf_continuous_supplement.py"
            )
            next(
                row
                for row in tampered_current_rows
                if row["repository_path"] == freezer_relative
            )["sha256"] = "0" * 64
            with (
                mock.patch.object(
                    coordinator_module,
                    "recursive_release_producer_rows",
                    return_value=tampered_current_rows,
                ),
                self.assertRaisesRegex(
                    CoordinatorContractError, "current producer drift"
                ),
            ):
                validate_nmf_continuous_bundle(paths)

            real_builder = coordinator_module.build_nmf_continuous_supplement

            def tampered_builder(*args, **kwargs):
                rebuilt_loadings, rebuilt_evidence = real_builder(*args, **kwargs)
                altered_evidence = copy.deepcopy(rebuilt_evidence)
                altered_evidence[0]["value"] = "999"
                return rebuilt_loadings, altered_evidence

            with (
                mock.patch.object(
                    coordinator_module,
                    "build_nmf_continuous_supplement",
                    side_effect=tampered_builder,
                ),
                self.assertRaisesRegex(
                    CoordinatorContractError,
                    "frozen products do not rederive from source copies",
                ),
            ):
                validate_nmf_continuous_bundle(paths)

    def test_frozen_nmf_historical_producer_hash_is_not_live_authority(self):
        paths = RealPaths.build(PROJECT_ROOT)
        historical = read_tsv_exact(
            paths.plan20_root / "nmf_continuous_supplement/producer_manifest.tsv",
            CORE_PRODUCER_FIELDS,
        )
        current = recursive_release_producer_rows(PROJECT_ROOT)
        historical_by_path = {row["repository_path"]: row for row in historical}
        current_by_path = {row["repository_path"]: row for row in current}
        contract_relative = (
            "scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py"
        )
        self.assertNotEqual(
            historical_by_path[contract_relative]["sha256"],
            current_by_path[contract_relative]["sha256"],
        )
        loadings, evidence = validate_nmf_continuous_bundle(paths)
        self.assertEqual((len(loadings), len(evidence)), (11040, 18))

    def test_real_cohort_census_rederives_corrected_denominator(self):
        integration = (
            PROJECT_ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
        )
        _, metadata = read_tsv_flexible(integration / "metadata/unified_metadata.csv")
        _, qc = read_tsv_flexible(integration / "qc/sample_qc_report.csv")
        rows = build_cohort_overview_rows(
            metadata,
            qc,
            {"GSE126848", "GSE135251", "GSE130970", "GSE213621", "GSE162694"},
        )
        observed = {row["number_role"]: int(row["value"]) for row in rows[:4]}
        self.assertEqual(observed, EXPECTED_COHORT_COUNTS)
        self.assertEqual(observed["qc_samples"], 1260)

    def test_recursive_producer_and_environment_inventory_is_complete(self):
        runtime = core_runtime_environment_rows()
        packages = core_package_export_rows()
        sys_path = core_sys_path_rows()
        imports = core_import_resolution_rows()
        producers = recursive_release_producer_rows(PROJECT_ROOT)
        self.assertEqual(len(runtime), 1)
        self.assertRegex(str(runtime[0]["executable_sha256"]), r"^[0-9a-f]{64}$")
        self.assertTrue(packages)
        self.assertTrue(sys_path)
        self.assertEqual([row["module"] for row in imports], ["pyarrow"])
        self.assertRegex(str(imports[0]["imported_file_sha256"]), r"^[0-9a-f]{64}$")
        jinja = [row for row in packages if row["package"] == "jinja2"]
        if len({row["version"] for row in jinja}) > 1:
            self.assertTrue(
                all(
                    row["resolution_status"] == "ambiguous_multiple_versions"
                    for row in jinja
                )
            )
        paths = {row["repository_path"] for row in producers}
        self.assertIn(
            "scripts/manuscript/program_context_v2/publish_plan60_workstream_handoff.py",
            paths,
        )
        self.assertIn(
            "scripts/manuscript/program_context_v2/coordinator/prepare_real_release.py",
            paths,
        )
        self.assertIn(
            "scripts/manuscript/program_context_v2/coordinator/tests/test_coordinator_contract.py",
            paths,
        )

    def test_publisher_contract_has_exact_unique_expanded_artifact_sets(self):
        contracts = handoff_publisher.workstream_contract(PROJECT_ROOT)
        self.assertEqual(
            {
                workstream: len(specifications)
                for workstream, (_, specifications) in contracts.items()
                if workstream != "PLAN50"
            },
            {"PLAN13": 13, "PLAN20": 37, "PLAN30": 15, "PLAN40": 21},
        )
        for workstream, (_, specifications) in contracts.items():
            artifact_ids = [row[0] for row in specifications]
            snapshots = [row[3] for row in specifications]
            self.assertEqual(len(artifact_ids), len(set(artifact_ids)), workstream)
            self.assertEqual(len(snapshots), len(set(snapshots)), workstream)
        plan20 = {row[0]: row for row in contracts["PLAN20"][1]}
        self.assertEqual(
            plan20["nmf_continuous_supplement"][3],
            "nmf_continuous_supplement.tsv",
        )
        plan50_root, plan50_specifications = contracts["PLAN50"]
        payload_rows = read_tsv_exact(
            plan50_root / "passport_release_manifest.tsv",
            handoff_publisher.PLAN50_PAYLOAD_MANIFEST_FIELDS,
        )
        payload_paths = {row["relative_path"] for row in payload_rows}
        terminal_paths = {
            relative for _, relative, _ in handoff_publisher.PLAN50_TERMINAL_CHAIN
        }
        self.assertEqual(
            terminal_paths,
            {
                "PASS06_VALIDATED",
                "passport_release_manifest.tsv",
                "passport_validation_report.tsv",
                "passport_terminal_provenance.tsv",
                "passport_plan60_handoff.tsv",
            },
        )
        self.assertEqual(
            len(plan50_specifications),
            len(payload_rows) + len(handoff_publisher.PLAN50_TERMINAL_CHAIN),
        )
        self.assertEqual(
            {row[3] for row in plan50_specifications},
            payload_paths | terminal_paths,
        )

    def test_prep_plan50_stale_manifest_producer_hash_blocks(self):
        producer = PLAN50_MANIFEST_PRODUCER_PATHS[0]
        producer_sha256 = sha256_file(PROJECT_ROOT / producer)
        frozen = {
            producer: {
                "sha256": producer_sha256,
            }
        }
        row = {"producer": producer, "producer_sha256": producer_sha256}
        coordinator_module._validate_plan50_producer_binding(
            row, frozen, "Plan50 payload manifest fixture.tsv"
        )
        row["producer_sha256"] = "0" * 64
        with self.assertRaisesRegex(
            CoordinatorContractError,
            "producer hash differs from the frozen REL01 producer",
        ):
            coordinator_module._validate_plan50_producer_binding(
                row, frozen, "Plan50 payload manifest fixture.tsv"
            )

    def test_rel01_plan50_stale_manifest_producer_hash_blocks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "passport_release_manifest.tsv"
            terminal_provenance = root / "passport_terminal_provenance.tsv"
            frozen_producers = [
                {
                    "repository_path": producer,
                    "sha256": sha256_file(PROJECT_ROOT / producer),
                }
                for producer in PLAN50_FROZEN_PRODUCER_PATHS
            ]
            frozen_hashes = {
                row["repository_path"]: row["sha256"] for row in frozen_producers
            }
            payload_rows = [
                {
                    "relative_path": f"payload_{index}.tsv",
                    "sha256": f"{index + 1:064x}",
                    "bytes": "1",
                    "producer": producer,
                    "environment": "synthetic-test",
                    "producer_sha256": frozen_hashes[producer],
                    "upstream_artifact_sha256": f"{index + 2:064x}",
                    "release_status": "validated_candidate_handoff",
                }
                for index, producer in enumerate(PLAN50_MANIFEST_PRODUCER_PATHS)
            ]
            atomic_write_tsv(
                manifest,
                payload_rows,
                handoff_publisher.PLAN50_PAYLOAD_MANIFEST_FIELDS,
            )
            terminal_rows = [
                {
                    "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
                    "bundle_uri": f"release://{PLAN50_ANALYSIS_RELEASE_ID}/",
                    "artifact_role": role,
                    "relative_path": relative,
                    "sha256": f"{index + 20:064x}",
                    "bytes": "1",
                    "producer": PLAN50_MANIFEST_PRODUCER_PATHS[-1],
                    "producer_sha256": frozen_hashes[
                        PLAN50_MANIFEST_PRODUCER_PATHS[-1]
                    ],
                }
                for index, (role, relative) in enumerate(
                    (
                        ("payload_manifest", "passport_release_manifest.tsv"),
                        ("plan60_handoff", "passport_plan60_handoff.tsv"),
                    )
                )
            ]
            atomic_write_tsv(
                terminal_provenance,
                terminal_rows,
                rel01_module.PLAN50_TERMINAL_PROVENANCE_FIELDS,
            )
            closure = mock.Mock(
                workstreams=[
                    mock.Mock(
                        row={"workstream_id": "PLAN50"},
                        artifacts=[
                            mock.Mock(source_path=manifest),
                            mock.Mock(source_path=terminal_provenance),
                        ],
                    )
                ]
            )
            result = validate_plan50_frozen_producer_binding(closure, frozen_producers)
            self.assertEqual(
                result["manifest_producer_count"],
                len(PLAN50_MANIFEST_PRODUCER_PATHS),
            )
            payload_rows[0]["producer_sha256"] = "0" * 64
            atomic_write_tsv(
                manifest,
                payload_rows,
                handoff_publisher.PLAN50_PAYLOAD_MANIFEST_FIELDS,
            )
            with self.assertRaisesRegex(
                ReleaseContractError,
                "producer hash differs from REL01 frozen bytes",
            ):
                validate_plan50_frozen_producer_binding(closure, frozen_producers)

    def test_stale_cohort_count_blocks(self):
        metadata = [{"sample_id": f"S{i}", "dataset": "GSE126848"} for i in range(1259)]
        qc = [
            {"sample_id": row["sample_id"], "pass_technical": "TRUE"}
            for row in metadata
        ]
        with self.assertRaisesRegex(
            CoordinatorContractError, "canonical numerical authority"
        ):
            build_cohort_overview_rows(metadata, qc, {"GSE126848"})

    def test_real_composition_adapter_preserves_22_column_universe(self):
        root = (
            PROJECT_ROOT
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "hotspot"
        )
        result_fields, results = read_tsv_flexible(
            root / "composition_sample_qc_v2.tsv"
        )
        testability_fields, testability = read_tsv_flexible(
            root / "composition_sample_qc_testability.tsv"
        )
        ready = one_row_tsv(root / "COMPOSITION_QC_READY", COMPOSITION_READY_FIELDS)
        self.assertEqual(result_fields, COMPOSITION_RESULT_FIELDS)
        self.assertEqual(testability_fields, COMPOSITION_TESTABILITY_FIELDS)
        rows = build_composition_rows(results, testability, ready)
        self.assertEqual(len(rows), 23)
        self.assertEqual(sum(row["plot_role"] == "mark" for row in rows), 22)
        self.assertEqual(
            sum(row["number_role"] == "composition_effect_arcsin_sqrt" for row in rows),
            16,
        )
        untestable = [row for row in rows if row["evidence_status"] == "untestable"]
        self.assertEqual(len(untestable), 6)
        self.assertTrue(
            all(row["value"] == "1" and not row["p_value"] for row in untestable)
        )
        self.assertTrue(all("not an effect" in row["unit"] for row in untestable))
        self.assertEqual(
            {row["panel_title"] for row in rows},
            {"QC-passing sample-level cell composition"},
        )
        self.assertTrue(
            all(
                row["denominator"] == "1221"
                for row in rows
                if row["number_role"] == "composition_effect_arcsin_sqrt"
            )
        )
        summary = next(
            row
            for row in rows
            if row["number_role"] == "bh_significant_testable_celltype_count"
        )
        self.assertEqual(summary["value"], "11")
        self.assertEqual(summary["denominator"], "16")
        self.assertEqual(summary["manuscript_included"], "true")

    def test_composition_bh_tamper_blocks(self):
        root = (
            PROJECT_ROOT
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "hotspot"
        )
        _, results = read_tsv_flexible(root / "composition_sample_qc_v2.tsv")
        _, testability = read_tsv_flexible(
            root / "composition_sample_qc_testability.tsv"
        )
        ready = one_row_tsv(root / "COMPOSITION_QC_READY", COMPOSITION_READY_FIELDS)
        broken = copy.deepcopy(results)
        broken[0]["padj"] = "0.999"
        with self.assertRaisesRegex(CoordinatorContractError, "BH"):
            validate_composition_contract(broken, testability, ready)

    def test_hotspot_adapter_uses_semantic_adjudication(self):
        registries, figures, semantics, designs, lodo = hotspot_fixture()
        rows = build_hotspot_rows(registries, figures, semantics, designs, lodo)
        self.assertEqual(len(rows), 7)
        self.assertEqual(
            {row["number_role"] for row in rows},
            {
                "stage_ordinal_beta",
                "stage_ordinal_se",
                "robust_display_flag",
                "leave_gse244832_direction_agreement_count",
            },
        )
        boundary = next(
            row
            for row in rows
            if row["number_role"] == "leave_gse244832_direction_agreement_count"
        )
        self.assertEqual(boundary["manuscript_included"], "true")
        self.assertIn("Only GSE244832 spans", boundary["claim_text"])
        self.assertFalse(
            any(row["evidence_status"] == "tested_negative" for row in rows)
        )

    def test_hotspot_tested_negative_authorization_blocks(self):
        registries, figures, semantics, designs, lodo = hotspot_fixture()
        semantics[7]["tested_negative_authorized"] = "TRUE"
        with self.assertRaisesRegex(CoordinatorContractError, "tested-negative"):
            build_hotspot_rows(registries, figures, semantics, designs, lodo)

    def test_hotspot_registry_projection_blocks_tampered_beta_q_or_n(self):
        registries, figures, _, _, _ = hotspot_fixture()
        for field in ("beta", "qvalue", "n_donors"):
            broken = copy.deepcopy(figures)
            broken[7][field] = "999"
            with self.assertRaisesRegex(CoordinatorContractError, "drift"):
                validate_hotspot_registry_projection(registries, broken)

    def test_negative_rule_requires_complete_passing_margin(self):
        row = {
            "evidence_result_id": "negative",
            "call_state": "tested_negative",
            "testability_state": "testable",
            "negative_call_rule_id": "EQ_V1",
            "negative_decision_boundary": 0.2,
            "negative_margin": 0.05,
            "negative_call_passed": True,
        }
        validate_passport_negative_semantics([row])
        for field, value in (
            ("negative_call_rule_id", None),
            ("negative_margin", -0.1),
            ("negative_call_passed", False),
        ):
            broken = dict(row)
            broken[field] = value
            with self.assertRaises(CoordinatorContractError):
                validate_passport_negative_semantics([broken])

    def test_nonnegative_row_rejects_reserved_negative_metadata(self):
        row = {
            "evidence_result_id": "indeterminate",
            "call_state": "indeterminate",
            "testability_state": "testable",
            "negative_call_rule_id": None,
            "negative_decision_boundary": 0.2,
            "negative_margin": None,
            "negative_call_passed": None,
        }
        with self.assertRaisesRegex(
            CoordinatorContractError, "reserved negative metadata"
        ):
            validate_passport_negative_semantics([row])

    def test_passport_adapter_is_ensembl_first_and_complete(self):
        fixture = passport_fixture()
        rows = build_passport_rows(*fixture)
        self.assertEqual(len(rows), 16)
        self.assertEqual({row["figure_id"] for row in rows}, {"Figure5"})
        self.assertEqual(
            {row["source_dependence"] for row in rows}, {"source_dependent"}
        )
        self.assertTrue(
            all(
                "figure5_vignette_selection_posthoc" in row["discovery_sources"]
                for row in rows
            )
        )
        self.assertFalse(
            any("program-to-gene" in row["allowed_wording"] for row in rows)
        )

    def test_passport_program_to_gene_expansion_blocks(self):
        genes, evidence, coverage, experiments, contexts, nodes, edges = (
            passport_fixture()
        )
        contexts[0]["gene_call_expansion_authorized"] = True
        with self.assertRaisesRegex(CoordinatorContractError, "program-to-gene"):
            build_passport_rows(
                genes, evidence, coverage, experiments, contexts, nodes, edges
            )

    def test_passport_duplicate_symbol_blocks(self):
        genes, evidence, coverage, experiments, contexts, nodes, edges = (
            passport_fixture()
        )
        genes.append(dict(genes[0]))
        with self.assertRaisesRegex(CoordinatorContractError, "maps to 2"):
            build_passport_rows(
                genes, evidence, coverage, experiments, contexts, nodes, edges
            )

    def test_real_spatial_adapter_preserves_native_uncertainty_and_directions(self):
        path = (
            PROJECT_ROOT
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "spatial_context_semantic_v2_2026-08-08/final_integration/figure4_program_matrix.tsv"
        )
        _, matrix = read_tsv_flexible(path)
        rows = build_spatial_rows(matrix)
        self.assertEqual(sum(row["plot_role"] == "mark" for row in rows), 18)
        self.assertEqual(
            sum(
                row["number_role"] == "matched_null_standard_deviation" for row in rows
            ),
            4,
        )
        self.assertEqual(
            {
                row["number_role"]
                for row in rows
                if row["record_id"].startswith("yak_module8_")
            },
            {
                "descriptive_median_slope_direction",
                "inferential_signed_stouffer_direction",
            },
        )
        self.assertTrue(
            all(
                row["group_id"] == "yak_module8"
                for row in rows
                if row["record_id"].startswith("yak_module8_")
            )
        )
        panel4b = [
            row
            for row in rows
            if row["figure_id"] == "Figure4" and row["panel_id"] == "4B"
        ]
        payload = build_pdf("Zonation-adjusted early lipid context", panel4b)
        self.assertIn(b"descriptive_median_slope_direction=positive", payload)
        self.assertIn(b"inferential_signed_stouffer_direction=negative", payload)
        self.assertFalse(
            any(row["evidence_status"] == "tested_negative" for row in rows)
        )

    def test_spatial_unresolved_biological_n_blocks(self):
        path = (
            PROJECT_ROOT
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "spatial_context_semantic_v2_2026-08-08/final_integration/figure4_program_matrix.tsv"
        )
        _, matrix = read_tsv_flexible(path)
        broken = copy.deepcopy(matrix)
        unresolved = next(
            row for row in broken if row["biological_unit_resolution"] == "unresolved"
        )
        unresolved["n_biological"] = "30"
        with self.assertRaisesRegex(
            CoordinatorContractError, "unresolved spatial unit"
        ):
            build_spatial_rows(broken)

    def test_real_myojin_adapter_uses_categorical_untestable_indicator(self):
        root = (
            PROJECT_ROOT
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "myojin_hlf"
        )
        _, classes = read_tsv_flexible(root / "class_effects.tsv")
        _, programs = read_tsv_flexible(root / "program_effects.tsv")
        rows = build_myojin_rows(classes, programs)
        self.assertEqual(len(rows), len(classes) + len(programs))
        class_rows = [row for row in rows if row["panel_id"] == "S1A"]
        omnibus = [
            row for row in class_rows if row["number_role"] == "class_omnibus_partial_f"
        ]
        pairwise = [
            row
            for row in class_rows
            if row["number_role"] == "class_pairwise_standardized_effect"
        ]
        self.assertEqual((len(omnibus), len(pairwise)), (7, 21))
        self.assertTrue(all(row["plot_role"] == "annotation" for row in omnibus))
        self.assertTrue(all(row["plot_role"] == "mark" for row in pairwise))
        self.assertTrue(
            all(
                "not comparable" in row["unit"] and not row["q_value"]
                for row in omnibus
            )
        )
        untestable = [row for row in rows if row["evidence_status"] == "untestable"]
        self.assertEqual(len(untestable), 1)
        self.assertEqual(untestable[0]["number_role"], "program_testability_indicator")
        self.assertEqual(untestable[0]["value"], "1")
        self.assertIn("not an effect", untestable[0]["unit"])
        self.assertEqual(untestable[0]["p_value"], "")

    def test_locked_blueprint_has_five_main_and_four_supplement_panels(self):
        blueprint = locked_release_blueprint()
        self.assertEqual(blueprint["figure_count"], 5)
        self.assertEqual(
            [row["figure_id"] for row in blueprint["figures"]], list(EXPECTED_FIGURES)
        )
        self.assertEqual(
            blueprint["figures"][3]["title"], EXPECTED_FIGURE_TITLES["Figure4"]
        )
        self.assertEqual(
            blueprint["figures"][1]["panels"][1]["title"],
            "QC-passing sample-level cell composition",
        )
        serialized = json.dumps(blueprint, sort_keys=True)
        self.assertNotIn("composition_shifts.csv", serialized)
        self.assertNotIn("cell_composition_raw", serialized)
        self.assertEqual(blueprint["myojin_role"], "supplement_only")
        self.assertEqual(
            [
                panel["panel_id"]
                for panel in blueprint["supplementary_figures"][0]["panels"]
            ],
            ["S1A", "S1B"],
        )
        self.assertEqual(
            [
                panel["panel_id"]
                for panel in blueprint["supplementary_figures"][1]["panels"]
            ],
            ["S2A", "S2B"],
        )
        nmf_source = next(
            source
            for source in blueprint["sources"]
            if source["source_key"] == "PLAN20:nmf_continuous_supplement"
        )
        self.assertEqual(
            nmf_source["snapshot_path"],
            "inputs/HS-V2/nmf_continuous_supplement.tsv",
        )

    def test_closure_has_exact_five_roles_and_myojin_supplement(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            handoffs = {}
            for workstream in ("PLAN13", "PLAN20", "PLAN30", "PLAN40", "PLAN50"):
                root = project / workstream
                root.mkdir()
                manifest = root / "plan60_terminal_artifacts.tsv"
                manifest.write_text("fixture\n", encoding="utf-8")
                handoffs[workstream] = Handoff(
                    workstream,
                    root,
                    manifest,
                    root / "signature.json",
                    sha256_file(manifest),
                    tuple(),
                )
            rows = closure_rows(project, handoffs, "coordinator", "2026-08-08")
            self.assertEqual([row["workstream_id"] for row in rows], list(handoffs))
            plan40 = next(row for row in rows if row["workstream_id"] == "PLAN40")
            self.assertEqual(plan40["terminal_state"], "accepted_supplement")
            self.assertEqual(plan40["include_main"], "false")
            self.assertEqual(plan40["include_supplement"], "true")

    def test_owner_signature_is_exact_and_nonpromoting(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            manifest.write_text("x\n", encoding="utf-8")
            digest = sha256_file(manifest)
            signature = root / "signature.json"
            signature.write_text(
                json.dumps(
                    {
                        "attestation_version": HANDOFF_SIGNATURE_VERSION,
                        "candidate_id": CANDIDATE_ID,
                        "workstream_id": "PLAN13",
                        "manifest_sha256": digest,
                        "signed_by": "owner",
                        "signed_at_utc": "2026-08-08T00:00:00Z",
                        "canonical_promotion_authorized": False,
                    }
                ),
                encoding="utf-8",
            )
            validate_owner_signature(signature, "PLAN13", digest)
            payload = json.loads(signature.read_text())
            payload["canonical_promotion_authorized"] = True
            signature.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(CoordinatorContractError, "promotion"):
                validate_owner_signature(signature, "PLAN13", digest)

    def test_handoff_publisher_validation_failure_is_fail_atomic(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            root = (
                project
                / "Analysis/Multimodal_Program_Projection/candidates"
                / CANDIDATE_ID
                / "spatial_context"
            )
            root.mkdir(parents=True)
            source = root / "READY"
            source.write_text("fixture\n", encoding="utf-8")
            broken_contract = {
                "PLAN13": (
                    root,
                    (("terminal", "READY", "", "READY"),),
                )
            }
            with mock.patch.object(
                handoff_publisher, "workstream_contract", return_value=broken_contract
            ):
                with self.assertRaisesRegex(
                    ReleaseContractError, "blank artifact role"
                ):
                    handoff_publisher.publish(
                        project, "PLAN13", "owner", "2026-08-08T00:00:00Z"
                    )
            self.assertFalse((root / "plan60_terminal_artifacts.tsv").exists())
            self.assertFalse(
                (root / "plan60_terminal_artifacts.signature.json").exists()
            )

    def test_handoff_publisher_post_link_failure_rolls_back_pair(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            root = (
                project
                / "Analysis/Multimodal_Program_Projection/candidates"
                / CANDIDATE_ID
                / "spatial_context"
            )
            root.mkdir(parents=True)
            source = root / "READY"
            source.write_text("fixture\n", encoding="utf-8")
            contract = {
                "PLAN13": (
                    root,
                    (("terminal", "READY", "terminal_verdict", "READY"),),
                )
            }
            with (
                mock.patch.object(
                    handoff_publisher, "workstream_contract", return_value=contract
                ),
                mock.patch.object(
                    handoff_publisher,
                    "validate_artifact_manifest",
                    side_effect=(
                        (object(),),
                        ReleaseContractError("post-link fixture"),
                    ),
                ),
            ):
                with self.assertRaisesRegex(ReleaseContractError, "post-link fixture"):
                    handoff_publisher.publish(
                        project, "PLAN13", "owner", "2026-08-08T00:00:00Z"
                    )
            self.assertFalse((root / "plan60_terminal_artifacts.tsv").exists())
            self.assertFalse(
                (root / "plan60_terminal_artifacts.signature.json").exists()
            )

    def test_two_clean_rebuild_comparison_requires_distinct_identical_builds(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            authority = SyntheticProject(base)

            def build(
                project: SyntheticProject, task_id: str, tracked_diff: str
            ) -> str:
                environment = {
                    "SLURM_JOB_ID": "900001",
                    "SLURM_ARRAY_TASK_ID": task_id,
                    "SLURM_STEP_ID": "batch",
                }
                with (
                    mock.patch.dict(os.environ, environment, clear=False),
                    mock.patch.object(
                        sys,
                        "argv",
                        [
                            "rel01_snapshot_candidate.py",
                            "--project-root",
                            str(project.root),
                        ],
                    ),
                    mock.patch.object(
                        rel01_module,
                        "repository_state",
                        return_value={
                            "repository_commit": "synthetic-fixture",
                            "tracked_diff_sha256": tracked_diff,
                        },
                    ),
                ):
                    snapshot_candidate(
                        project.root,
                        project.closure_path,
                        project.protected_scopes_path,
                        project.protected_baseline_path,
                        fixture_mode=True,
                    )
                stage_base_inputs(
                    project.root,
                    project.base_selection_path,
                    sha256_file(project.base_selection_path),
                    fixture_mode=True,
                )
                build_candidate_tables(project.root, fixture_mode=True)
                render_candidate_panels(project.root, fixture_mode=True)
                build_candidate_manuscript(project.root, fixture_mode=True)
                self.assertFalse(
                    (
                        project.candidate_root
                        / "manifests/scientific_validation_registry.tsv"
                    ).exists()
                )
                return f"slurm:900001:array:{task_id}:step:batch"

            authority_identity = build(authority, "0", "a" * 64)
            comparison = (
                authority.candidate_root / "retained_clean_rebuilds/comparison_1"
            )
            build_a = SyntheticProject(comparison / "build_a/project")
            build_b = SyntheticProject(comparison / "build_b/project")
            build_a_identity = build(build_a, "1", "b" * 64)
            build_b_identity = build(build_b, "2", "b" * 64)
            self.assertNotEqual(authority_identity, build_a_identity)

            with self.assertRaisesRegex(ReleaseContractError, "CLI process IDs"):
                derive_clean_rebuild_rows(
                    base,
                    "build_a",
                    "slurm:operator:array:1:step:batch",
                    "build_b",
                    build_b_identity,
                    "comparison_1",
                    fixture_mode=True,
                )
            rows = compare_clean_rebuilds(
                base,
                "build_a",
                build_a_identity,
                "build_b",
                build_b_identity,
                "comparison_1",
                fixture_mode=True,
            )
            self.assertGreater(len(rows), 1)
            self.assertEqual(rows[0]["equal"], "true")
            self.assertEqual(
                len(
                    {
                        rows[0]["authority_snapshot_spec_sha256"],
                        rows[0]["build_a_snapshot_spec_sha256"],
                    }
                ),
                2,
            )
            self.assertEqual(
                rows[0]["rebuild_source_spec_sha256"],
                json.loads(
                    (
                        authority.candidate_root / "manifests/rebuild_source_spec.json"
                    ).read_text(encoding="utf-8")
                )["rebuild_source_spec_sha256"],
            )
            self.assertEqual(
                rows[0]["authority_candidate_root"], str(authority.candidate_root)
            )
            self.assertEqual(
                rows[0]["comparison_root"],
                str(comparison),
            )
            self.assertEqual(
                len(
                    {
                        rows[0]["authority_execution_context_sha256"],
                        rows[0]["build_a_execution_context_sha256"],
                        rows[0]["build_b_execution_context_sha256"],
                    }
                ),
                3,
            )
            self.assertTrue(
                (
                    authority.candidate_root
                    / "manifests/two_clean_rebuild_comparison.tsv"
                ).is_file()
            )
            report = (
                authority.candidate_root / "manifests/two_clean_rebuild_comparison.tsv"
            )
            self.assertEqual(
                len(revalidate_clean_rebuild_report(base, report, fixture_mode=True)),
                len(rows),
            )
            product_rows = read_tsv_exact(
                build_b.candidate_root / "manifests/rel02_04_product_manifest.tsv",
                TRANSITION_PRODUCT_FIELDS,
            )
            tampered = build_b.root / product_rows[0]["project_relative_path"]
            tampered.write_bytes(tampered.read_bytes() + b"tamper")
            with self.assertRaisesRegex(
                ReleaseContractError, "transition product drift"
            ):
                revalidate_clean_rebuild_report(base, report, fixture_mode=True)

    def test_two_clean_rebuild_comparator_rejects_skeletal_candidates(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            authority = (
                base / "RNA-seq/results/manuscript_release/candidates" / CANDIDATE_ID
            )
            comparison = authority / "retained_clean_rebuilds/comparison_1"
            candidates = (
                authority,
                comparison / "build_a/project" / CANDIDATE_RELATIVE,
                comparison / "build_b/project" / CANDIDATE_RELATIVE,
            )
            for candidate in candidates:
                (candidate / "manifests").mkdir(parents=True)
            with self.assertRaisesRegex(
                ReleaseContractError, "required candidate directory"
            ):
                derive_clean_rebuild_rows(
                    base,
                    "build_a",
                    "local:a:1:20260808T000000000000Z",
                    "build_b",
                    "local:b:2:20260808T000001000000Z",
                    "comparison_1",
                    fixture_mode=True,
                )

    def test_source_rows_require_negative_adequacy_rule(self):
        row = source_row(**{field: "x" for field in SOURCE_EVIDENCE_FIELDS})
        row.update(
            record_id="negative",
            value="0",
            p_value="",
            q_value="",
            evidence_status="tested_negative",
            negative_adequacy_criterion="",
            source_dependence="reused_source",
            discovery_sources="fixture_source",
            evaluation_sources="fixture_source",
            reuse_detail="fixture reuse",
            independence_boundary="fixture is not independent",
            manuscript_included="false",
            plot_order="1",
            is_control="false",
        )
        with self.assertRaisesRegex(CoordinatorContractError, "adequate-negative rule"):
            validate_source_rows([row], "fixture")

    def test_bh_rederivation(self):
        adjusted = bh_adjust([0.01, 0.04, 0.03, 0.002])
        for observed, expected in zip(adjusted, [0.02, 0.04, 0.04, 0.008], strict=True):
            self.assertAlmostEqual(observed, expected)

    def test_output_guard_refuses_outside_and_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            safe = (
                project
                / "scripts/manuscript/program_context_v2/coordinator/prepared/DECISION"
            )
            self.assertEqual(assert_safe_output(project, safe), safe.absolute())
            with self.assertRaisesRegex(CoordinatorContractError, "beneath"):
                assert_safe_output(project, project / "elsewhere")
            safe.mkdir(parents=True)
            with self.assertRaisesRegex(CoordinatorContractError, "overwrite"):
                assert_safe_output(project, safe)

    def test_prepare_missing_handoff_blocks_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            output = (
                project
                / "scripts/manuscript/program_context_v2/coordinator/prepared/DECISION"
            )
            with self.assertRaisesRegex(
                CoordinatorContractError, "manifest is missing"
            ):
                prepare(
                    project,
                    output,
                    "coordinator",
                    "2026-08-08",
                    "2026-08-08T00:00:00Z",
                    "DECISION",
                )
            self.assertFalse(output.exists())

    def test_real_preparation_adapters_are_hash_bound_and_nonpromoting(self):
        """Exercise every real adapter without manufacturing owner approval."""

        paths = RealPaths.build(PROJECT_ROOT)
        terminal_by_workstream = {
            "PLAN13": paths.plan13_root / "SEMANTIC_V2_READY",
            "PLAN20": paths.plan20_root / "SEMANTIC_ADJUDICATION_READY",
            "PLAN30": paths.plan30_root / "GEN_TERMINAL_CLOSURE_READY",
            "PLAN40": paths.plan40_root / "PHASE_C_VALIDATED",
            # Plan50 is being re-frozen independently.  Its immutable scientific
            # products are exercised here, while the terminal-count seal below
            # is synthesized in memory and never published as an owner handoff.
            "PLAN50": paths.plan50_root / "passport_build_status.tsv",
        }
        handoffs = {}
        for workstream, terminal in terminal_by_workstream.items():
            self.assertTrue(terminal.is_file(), terminal)
            digest = sha256_file(terminal)
            handoffs[workstream] = Handoff(
                workstream, terminal.parent, terminal, terminal, digest, tuple()
            )

        gene_rows = preparation_module.read_parquet_records(
            paths.plan50_root / "passport_gene_index.parquet"
        )
        evidence_rows = preparation_module.read_parquet_records(
            paths.plan50_root / "passport_evidence_long.parquet"
        )
        program_index_rows = preparation_module.read_parquet_records(
            paths.plan50_root / "passport_program_index.parquet"
        )
        program_context_rows = preparation_module.read_parquet_records(
            paths.plan50_root / "passport_program_context.parquet"
        )
        synthetic_plan50_ready = {
            "analysis_release_id": PLAN50_ANALYSIS_RELEASE_ID,
            "status": "synthetic_test_count_seal_only",
            "selection_sha256": "0" * 64,
            "manifest_sha256": "0" * 64,
            "validation_report_sha256": "0" * 64,
            "n_genes": str(len(gene_rows)),
            "n_evidence_rows": str(len(evidence_rows)),
            "n_programs": str(len(program_index_rows)),
            "n_program_context_rows": str(len(program_context_rows)),
            "automated_validation": "test_only",
            "manual_acceptance": "false",
            "handoff_allowed": "false",
            "canonical_promotion_authorized": "false",
            "scientific_call_recomputed": "false",
            "validated_at_utc": "test_only",
        }

        def terminal_count_seal(path, fields):
            if Path(path).name == "PASS06_VALIDATED":
                return {field: synthetic_plan50_ready[field] for field in fields}
            return one_row_tsv(path, fields)

        prepared_parent = (
            PROJECT_ROOT / "scripts/manuscript/program_context_v2/coordinator/prepared"
        )
        prepared_parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix="real-adapter-test-", dir=prepared_parent
        ) as temporary:
            output = Path(temporary) / "decision"
            protected_scopes = [
                {
                    "scope_id": "test_scope",
                    "protection_class": "test_only",
                    "root_path": "README.md",
                }
            ]
            protected_baseline = [
                {
                    "scope_id": "test_scope",
                    "relative_path": "README.md",
                    "sha256": sha256_file(PROJECT_ROOT / "README.md"),
                    "bytes": (PROJECT_ROOT / "README.md").stat().st_size,
                }
            ]
            with (
                mock.patch.object(
                    preparation_module, "validate_all_handoffs", return_value=handoffs
                ),
                mock.patch.object(
                    preparation_module, "assert_safe_output", return_value=output
                ),
                mock.patch.object(
                    preparation_module,
                    "protected_scope_rows",
                    return_value=protected_scopes,
                ),
                mock.patch.object(
                    preparation_module,
                    "protected_baseline_rows",
                    return_value=protected_baseline,
                ),
                mock.patch.object(
                    preparation_module, "one_row_tsv", side_effect=terminal_count_seal
                ),
            ):
                ready = preparation_module.prepare(
                    PROJECT_ROOT,
                    output,
                    "test-coordinator",
                    "2026-08-08",
                    "2026-08-08T00:00:00Z",
                    "real-adapter-test",
                )

            self.assertFalse(ready["rel00_05_executed"])
            self.assertFalse(ready["candidate_root_created"])
            self.assertFalse(ready["canonical_promotion_authorized"])

            selection_rows = read_tsv_exact(
                output / "base_input_selection.tsv", BASE_SELECTION_FIELDS
            )
            fibrosis_rows = [
                row
                for row in selection_rows
                if row["artifact_role"]
                in {"fibrosis_transition_raw", FIBROSIS_VALIDATION_ROLE}
            ]
            self.assertEqual(len(fibrosis_rows), len(FIBROSIS_BUNDLE_FILES))
            self.assertEqual(
                sum(
                    row["artifact_id"] == FIBROSIS_PRIMARY_ARTIFACT_ID
                    and row["artifact_role"] == "fibrosis_transition_raw"
                    for row in fibrosis_rows
                ),
                1,
            )
            self.assertEqual(
                sum(
                    row["artifact_role"] == FIBROSIS_VALIDATION_ROLE
                    for row in fibrosis_rows
                ),
                len(FIBROSIS_BUNDLE_FILES) - 1,
            )
            fibrosis_by_id = {row["artifact_id"]: row for row in fibrosis_rows}
            self.assertEqual(
                fibrosis_by_id[fibrosis_artifact_id("READY")]["source_sha256"],
                FIBROSIS_ACCEPTED_READY_SHA256,
            )
            self.assertEqual(
                fibrosis_by_id[
                    fibrosis_artifact_id("manifests/validated_artifact_manifest.tsv")
                ]["source_sha256"],
                FIBROSIS_ACCEPTED_VALIDATED_MANIFEST_SHA256,
            )
            self.assertEqual(
                fibrosis_by_id[FIBROSIS_PRIMARY_ARTIFACT_ID]["source_sha256"],
                FIBROSIS_ACCEPTED_PRIMARY_RESULT_SHA256,
            )

            adapter_rows = read_tsv_exact(
                output / "adapter_provenance.tsv",
                preparation_module.ADAPTER_PROVENANCE_FIELDS,
            )
            self.assertEqual(len(adapter_rows), 7)
            producer_rows = read_tsv_exact(
                output / "environment/recursive_release_producers.tsv",
                CORE_PRODUCER_FIELDS,
            )
            producer_by_path = {row["repository_path"]: row for row in producer_rows}
            self.assertIn(NMF_FREEZER_PRODUCER_PATH, producer_by_path)

            # REL01 must freeze the complete prepared recursive producer
            # universe, not only the producer subset named by adapter rows.
            frozen_rel01 = producer_script_inventory()
            frozen_by_path = {str(row["repository_path"]): row for row in frozen_rel01}
            self.assertTrue(set(producer_by_path).issubset(frozen_by_path))
            for repository_path, recursive in producer_by_path.items():
                frozen = frozen_by_path[repository_path]
                self.assertEqual(frozen["producer_id"], recursive["producer_id"])
                self.assertEqual(frozen["sha256"], recursive["sha256"])
                self.assertEqual(int(frozen["bytes"]), int(recursive["bytes"]))
            file_set_before = {
                path.relative_to(output).as_posix()
                for path in output.rglob("*")
                if path.is_file()
            }
            validated = validate_prepared_adapter_provenance(
                PROJECT_ROOT,
                output / "base_input_selection.tsv",
                frozen_rel01,
            )
            self.assertEqual(validated["n_adapters"], 7)

            omitted_nmf = [
                row
                for row in frozen_rel01
                if row["repository_path"] != NMF_FREEZER_PRODUCER_PATH
            ]
            with self.assertRaisesRegex(
                ReleaseContractError,
                "recursive producer absent from REL01 freeze.*519_freeze_nmf",
            ):
                validate_prepared_adapter_provenance(
                    PROJECT_ROOT,
                    output / "base_input_selection.tsv",
                    omitted_nmf,
                )

            # Simulate release_products.py changing after preparation.  The
            # recursive/frozen temporal gate must block before REL01 creates a
            # candidate root or writes any other artifact.
            post_preparation_mutation = copy.deepcopy(frozen_rel01)
            release_products_row = next(
                row
                for row in post_preparation_mutation
                if row["repository_path"]
                == "scripts/manuscript/program_context_v2/release_products.py"
            )
            release_products_row["sha256"] = "0" * 64
            with self.assertRaisesRegex(
                ReleaseContractError,
                "recursive producer differs from REL01 freeze.*release_products.py",
            ):
                validate_prepared_adapter_provenance(
                    PROJECT_ROOT,
                    output / "base_input_selection.tsv",
                    post_preparation_mutation,
                )
            self.assertEqual(
                file_set_before,
                {
                    path.relative_to(output).as_posix()
                    for path in output.rglob("*")
                    if path.is_file()
                },
            )

            # Builder and REL05 import one foundational producer universe, and
            # the actual non-fixture preparation payload can satisfy it.
            self.assertIs(
                scientific_registry_module.SCIENTIFIC_REPORT_PRODUCER_PATHS,
                SCIENTIFIC_REPORT_PRODUCER_PATHS,
            )
            self.assertIs(
                rel05_module.SCIENTIFIC_REPORT_PRODUCER_PATHS,
                SCIENTIFIC_REPORT_PRODUCER_PATHS,
            )
            producer_manifest = output / "environment/recursive_release_producers.tsv"
            source_index = {
                "BASE:recursive_release_producers": {
                    "workstream_id": "BASE",
                    "artifact_role": "producer_manifest",
                    "snapshot_path": producer_manifest.relative_to(output).as_posix(),
                    "sha256": sha256_file(producer_manifest),
                }
            }
            report_binding = scientific_producer_binding(
                PROJECT_ROOT,
                output,
                source_index,
            )
            self.assertEqual(
                tuple(
                    row["repository_path"]
                    for row in report_binding["producer_bindings"]
                ),
                SCIENTIFIC_REPORT_PRODUCER_PATHS,
            )
            self.assertIn(
                "scripts/manuscript/program_context_v2/rebuild_source_spec.py",
                SCIENTIFIC_REPORT_PRODUCER_PATHS,
            )
            for row in adapter_rows:
                self.assertEqual(
                    len(Path(row["output_relative_path"]).parts),
                    1,
                    row["artifact_id"],
                )
                self.assertEqual(
                    row["producer_manifest_relative_path"],
                    "inputs/BASE/recursive_release_producers.tsv",
                )
                self.assertEqual(row["canonical_promotion_authorized"], "false")
                bindings = json.loads(row["producer_bindings_json"])
                expected_producers = set(
                    ADAPTER_CONTRACT[row["artifact_id"]]["required_producers"]
                )
                self.assertEqual(
                    {binding["repository_path"] for binding in bindings},
                    expected_producers,
                )
                for binding in bindings:
                    frozen = producer_by_path[binding["repository_path"]]
                    self.assertEqual(binding["sha256"], frozen["sha256"])
                    self.assertEqual(binding["bytes"], int(frozen["bytes"]))
                for binding in json.loads(row["input_bindings_json"]):
                    source = PROJECT_ROOT / binding["source_path"]
                    self.assertEqual(binding["sha256"], sha256_file(source))
                    self.assertEqual(binding["bytes"], source.stat().st_size)

            hotspot_rows = read_tsv_exact(
                output / "evidence_rows/hotspot_programs.tsv",
                SOURCE_EVIDENCE_FIELDS,
            )
            supplement_rows = [
                row for row in hotspot_rows if row["figure_id"] == "FigureS2"
            ]
            self.assertEqual(len(supplement_rows), 18)
            self.assertEqual(
                {row["panel_id"] for row in supplement_rows}, {"S2A", "S2B"}
            )
            hotspot_provenance = next(
                row for row in adapter_rows if row["artifact_id"] == "hotspot_programs"
            )
            self.assertIn(
                "Analysis/SingleCell/scripts/hotspot_modules/"
                "519_freeze_nmf_continuous_supplement.py",
                {
                    binding["repository_path"]
                    for binding in json.loads(
                        hotspot_provenance["producer_bindings_json"]
                    )
                },
            )
            hotspot_inputs = {
                binding["source_path"]
                for binding in json.loads(hotspot_provenance["input_bindings_json"])
            }
            self.assertIn(
                "Analysis/Multimodal_Program_Projection/candidates/"
                f"{CANDIDATE_ID}/hotspot/nmf_continuous_supplement/source_manifest.tsv",
                hotspot_inputs,
            )

            alias_rows = read_tsv_exact(
                output / "source_alias_audit.tsv",
                preparation_module.SOURCE_ALIAS_AUDIT_FIELDS,
            )
            self.assertTrue(
                any(
                    row["alias_equivalence_status"]
                    == "not_audited_cross_alias_equivalence"
                    for row in alias_rows
                )
            )
            blueprint = json.loads(
                (output / "release_blueprint.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                sum(len(figure["panels"]) for figure in blueprint["figures"]), 11
            )
            self.assertEqual(
                sum(
                    len(figure["panels"])
                    for figure in blueprint["supplementary_figures"]
                ),
                4,
            )

    def test_plan50_nested_release_identity_is_explicit(self):
        self.assertEqual(
            PLAN50_ANALYSIS_RELEASE_ID,
            "program-context-v2-candidate-2026-08-07-passports-v1",
        )
        self.assertNotEqual(PLAN50_ANALYSIS_RELEASE_ID, CANDIDATE_ID)


if __name__ == "__main__":
    unittest.main()
