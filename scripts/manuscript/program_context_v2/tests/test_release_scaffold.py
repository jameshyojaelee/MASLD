#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


TEST_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TEST_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(TEST_DIR))

from fixture_factory import SyntheticProject  # noqa: E402
from fibrosis_candidate_contract import (  # noqa: E402
    FIBROSIS_ACCEPTED_PRIMARY_RESULT_SHA256,
    FIBROSIS_EXPECTED_N_GENES,
    FIBROSIS_PRIMARY_ARTIFACT_ID,
    validate_fibrosis_candidate_bundle,
)
from release_common import (  # noqa: E402
    CANDIDATE_ID,
    PLAN50_FROZEN_PRODUCER_PATHS,
    ReleaseContractError,
    allowed_workstream_roots,
    atomic_write_tsv,
    assert_candidate_root,
    candidate_root,
    read_tsv_exact,
    sha256_file,
    validate_closure,
)
from rel01_snapshot_candidate import (  # noqa: E402
    PRODUCER_SCRIPT_FIELDS,
    snapshot_candidate,
)
from rel02_build_candidate_tables import build_candidate_tables  # noqa: E402
from rel02_stage_base_inputs import stage_base_inputs  # noqa: E402
from rel03_render_candidate_panels import build_pdf, render_candidate_panels  # noqa: E402
from rel04_build_candidate_manuscript import build_candidate_manuscript  # noqa: E402
from rel05_validate_candidate import validate_candidate  # noqa: E402
from release_products import (  # noqa: E402
    BASE_SELECTION_FIELDS,
    FIGURE_SOURCE_MANIFEST_FIELDS,
    NUMBERS_FIELDS,
    SOURCE_ROW_FIELDS,
    contains_nmf_discrete_class_claim,
    validate_nmf_source_language,
)


class ReleaseScaffoldTests(unittest.TestCase):
    def test_nmf_subtype_predicate_is_clause_local_and_paraphrase_aware(self):
        for wording in (
            "six subtypes were recovered",
            "6 patient clusters were identified",
            "the model assigned 4 groups",
            "six clusters were not stable",
            "six clusters, not subtypes, were retained",
            "hard clusters separate the cohort",
            "the model separated patients into groups",
            "not reproducible, but stable patient clusters emerged",
        ):
            self.assertTrue(contains_nmf_discrete_class_claim(wording), wording)
        for wording in (
            "not used as a patient stratification system",
            "no patient groups are inferred",
            "patient groups are not inferred",
            "these are not discrete subtypes",
            "factor stability does not establish stable patient classes",
            "factor stability does not imply stable patient partitions",
        ):
            self.assertFalse(contains_nmf_discrete_class_claim(wording), wording)

    def test_nmf_generated_language_is_an_exact_closed_template(self):
        project = self.make_project()
        nmf = next(
            path
            for path in project.source_paths["PLAN20"]
            if path.name == "nmf_continuous_supplement.tsv"
        )
        rows = read_tsv_exact(nmf, SOURCE_ROW_FIELDS[:-2])
        for row in rows:
            validate_nmf_source_language(row)
        omitted_boundary = dict(rows[0])
        omitted_boundary["claim_text"] = "k=4 P1 is a continuous axis."
        with self.assertRaisesRegex(ReleaseContractError, "approved NMF language"):
            validate_nmf_source_language(omitted_boundary)
        numeric_assignment = dict(rows[0])
        numeric_assignment["claim_text"] = "The model assigned 4 patient clusters."
        with self.assertRaisesRegex(ReleaseContractError, "approved NMF language"):
            validate_nmf_source_language(numeric_assignment)

    def test_plan13_allows_only_the_locked_semantic_v2_sibling_in_addition_to_legacy_roots(
        self,
    ):
        project = self.make_project()
        roots = allowed_workstream_roots(project.root, "PLAN13")
        self.assertEqual(len(roots), 3)
        self.assertIn(
            (
                project.root
                / "Analysis/Multimodal_Program_Projection/candidates"
                / CANDIDATE_ID
                / "spatial_context_semantic_v2_2026-08-08"
            ).resolve(),
            roots,
        )

    def make_project(self):
        temporary = tempfile.TemporaryDirectory()
        project = SyntheticProject(Path(temporary.name))
        self.addCleanup(temporary.cleanup)
        return project

    def snapshot(self, project: SyntheticProject):
        return snapshot_candidate(
            project.root,
            project.closure_path,
            project.protected_scopes_path,
            project.protected_baseline_path,
            fixture_mode=True,
        )

    def build_products(self, project: SyntheticProject):
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        build_candidate_tables(project.root, fixture_mode=True)
        render_candidate_panels(project.root, fixture_mode=True)
        return build_candidate_manuscript(project.root, fixture_mode=True)

    def mutate_evidence(self, project, workstream, filename, transform):
        source = next(
            path for path in project.source_paths[workstream] if path.name == filename
        )
        rows = read_tsv_exact(source, SOURCE_ROW_FIELDS[:-2])
        atomic_write_tsv(source, transform(rows), SOURCE_ROW_FIELDS[:-2])

        def update_manifest(manifest_rows):
            for row in manifest_rows:
                if row["source_path"].endswith(filename):
                    row["source_sha256"] = sha256_file(source)
                    row["source_bytes"] = str(source.stat().st_size)
            return manifest_rows

        project.rewrite_artifact_manifest(workstream, update_manifest)
        return source

    def test_nonconfirmatory_supplement_closure_is_idempotent_and_validates(self):
        project = self.make_project()
        bundle = validate_closure(project.root, project.closure_path)
        self.assertEqual(len(bundle.workstreams), 5)
        self.assertEqual(
            next(
                item.row["terminal_state"]
                for item in bundle.workstreams
                if item.row["workstream_id"] == "PLAN40"
            ),
            "accepted_supplement",
        )
        first = self.snapshot(project)
        tracked = [
            project.candidate_root / "manifests/snapshot_spec.json",
            project.candidate_root / "manifests/input_snapshot_manifest.tsv",
            project.candidate_root / "manifests/downstream_inputs.tsv",
            project.candidate_root / "manifests/candidate_state.json",
        ]
        first_hashes = {path.name: sha256_file(path) for path in tracked}
        producers = read_tsv_exact(
            project.candidate_root / "manifests/producer_script_manifest.tsv",
            PRODUCER_SCRIPT_FIELDS,
        )
        frozen_repository_paths = {row["repository_path"] for row in producers}
        self.assertTrue(
            set(PLAN50_FROZEN_PRODUCER_PATHS).issubset(frozen_repository_paths)
        )
        second = self.snapshot(project)
        second_hashes = {path.name: sha256_file(path) for path in tracked}
        self.assertEqual(first["spec_sha256"], second["spec_sha256"])
        self.assertEqual(first_hashes, second_hashes)
        for path in (project.candidate_root / "inputs").rglob("*"):
            if path.is_file():
                self.assertFalse(path.is_symlink())
        product_result = self.build_products(project)
        self.assertEqual(product_result["figure_count"], 5)
        self.assertEqual(product_result["myojin_role"], "supplement_only")
        figure_manifest = read_tsv_exact(
            project.candidate_root / "manifests/figure_source_manifest.tsv",
            FIGURE_SOURCE_MANIFEST_FIELDS,
        )
        self.assertEqual(len(figure_manifest), 15)
        self.assertEqual(
            {
                row["figure_id"]
                for row in figure_manifest
                if not row["figure_id"].startswith("FigureS")
            },
            {"Figure1", "Figure2", "Figure3", "Figure4", "Figure5"},
        )
        self.assertEqual(
            {
                (row["figure_id"], row["panel_id"])
                for row in figure_manifest
                if row["figure_id"].startswith("FigureS")
            },
            {
                ("FigureS1", "S1A"),
                ("FigureS1", "S1B"),
                ("FigureS2", "S2A"),
                ("FigureS2", "S2B"),
            },
        )
        for row in figure_manifest:
            panel = project.root / row["pdf_path"]
            self.assertEqual(panel.suffix, ".pdf")
            payload = panel.read_bytes()
            self.assertIn(b"/BaseFont /Helvetica", payload)
            self.assertIn(b"/F1 6 Tf", payload)
        numbers = read_tsv_exact(
            project.candidate_root / "tables/numbers_ledger.tsv", NUMBERS_FIELDS
        )
        observed_stage = {
            row["model_contrast"]: row["raw_value"]
            for row in numbers
            if row["number_role"] == "deg_count"
        }
        self.assertEqual(
            observed_stage,
            {
                "F1_vs_F0": "0",
                "F2_vs_F1": "0",
                "F3_vs_F2": "0",
                "F4_vs_F3": "0",
            },
        )
        manuscript = (
            project.root / "docs/manuscript/candidates" / CANDIDATE_ID / "manuscript.md"
        ).read_text(encoding="utf-8")
        self.assertNotIn("Figure 6", manuscript)
        abstract = manuscript.split("## Abstract", 1)[1].split("## Introduction", 1)[0]
        self.assertNotIn("Myojin", abstract)
        self.assertIn("supplementary", manuscript.lower())
        deterministic_products = [
            project.candidate_root / "manifests/base_input_transition.json",
            project.candidate_root / "manifests/rel02_state.json",
            project.candidate_root / "manifests/figure_source_manifest.tsv",
            project.candidate_root / "manifests/rel03_state.json",
            project.candidate_root / "manifests/rel02_04_transition.json",
            project.root
            / "docs/manuscript/candidates"
            / CANDIDATE_ID
            / "manuscript.md",
        ]
        first_product_hashes = {
            path.relative_to(project.root).as_posix(): sha256_file(path)
            for path in deterministic_products
        }
        self.build_products(project)
        second_product_hashes = {
            path.relative_to(project.root).as_posix(): sha256_file(path)
            for path in deterministic_products
        }
        self.assertEqual(first_product_hashes, second_product_hashes)
        project.install_scientific_registry()
        report = validate_candidate(project.root, fixture_mode=True, write_report=True)
        self.assertEqual(report["status"], "SYNTHETIC_PASS")
        self.assertFalse(report["full_release_pass"])
        self.assertFalse(report["canonical_promotion_authorized"])

    def test_post_snapshot_build_does_not_require_live_upstream_sources(self):
        project = self.make_project()
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        upstream_paths = {
            project.closure_path,
            project.base_selection_path,
            *project.base_source_paths.values(),
            *(path for paths in project.source_paths.values() for path in paths),
            *project.manifest_paths.values(),
        }
        for path in upstream_paths:
            path.unlink(missing_ok=True)
        build_candidate_tables(project.root, fixture_mode=True)
        render_candidate_panels(project.root, fixture_mode=True)
        build_candidate_manuscript(project.root, fixture_mode=True)
        project.install_scientific_registry()
        report = validate_candidate(project.root, fixture_mode=True, write_report=False)
        self.assertEqual(report["status"], "SYNTHETIC_PASS")

    def test_fixture_pdf_renderer_accepts_an_all_null_panel(self):
        payload = build_pdf(
            "All-null compatibility fixture",
            [
                {
                    "plot_role": "mark",
                    "group_id": "null_a",
                    "value": "0",
                    "is_control": "false",
                    "label": "program A",
                    "display_value": "0.000",
                },
                {
                    "plot_role": "mark",
                    "group_id": "null_b",
                    "value": "0",
                    "is_control": "true",
                    "label": "program B",
                    "display_value": "0.000",
                },
            ],
        )
        self.assertTrue(payload.startswith(b"%PDF-1.4"))
        self.assertTrue(payload.rstrip().endswith(b"%%EOF"))
        self.assertIn(b"/BaseFont /Helvetica", payload)
        self.assertIn(b"/F1 6 Tf", payload)
        self.assertIn(b"All prespecified rows retained", payload)

    def test_fixture_pdf_renderer_marks_untestable_without_a_zero_effect_bar(self):
        payload = build_pdf(
            "Untestable compatibility fixture",
            [
                {
                    "plot_role": "mark",
                    "group_id": "untestable_a",
                    "value": "1",
                    "is_control": "false",
                    "label": "program A",
                    "display_value": "untestable",
                    "evidence_status": "untestable",
                }
            ],
        )
        self.assertIn(b"X=untestable", payload)
        self.assertIn(b"untestable", payload)

    def test_fixture_pdf_renderer_uses_neutral_non_effect_gate_glyphs(self):
        payload = build_pdf(
            "Non-effect gate compatibility fixture",
            [
                {
                    "plot_role": "mark",
                    "group_id": "skip_a",
                    "value": "1",
                    "is_control": "false",
                    "label": "dataset A",
                    "display_value": "skipped",
                    "evidence_status": "skipped",
                },
                {
                    "plot_role": "mark",
                    "group_id": "na_b",
                    "value": "1",
                    "is_control": "false",
                    "label": "dataset B",
                    "display_value": "not applicable",
                    "evidence_status": "not_applicable",
                },
            ],
        )
        self.assertIn(b"hollow=skipped", payload)
        self.assertIn(b"/=not applicable", payload)
        self.assertIn(b"none are zero effects", payload)

    def test_fixture_pdf_renderer_keeps_attached_and_orphan_direction_annotations(self):
        payload = build_pdf(
            "Direction annotation compatibility fixture",
            [
                {
                    "plot_role": "mark",
                    "group_id": "program_a",
                    "value": "0.2",
                    "is_control": "false",
                    "label": "program A",
                    "display_value": "0.200",
                    "evidence_status": "indeterminate",
                    "plot_order": "1",
                },
                {
                    "plot_role": "annotation",
                    "group_id": "program_a",
                    "value": "0.1",
                    "is_control": "false",
                    "label": "Stromal ECM",
                    "display_value": "[-0.1, 0.3]",
                    "number_role": "confidence_interval",
                    "plot_order": "2",
                },
                {
                    "plot_role": "annotation",
                    "group_id": "program_a",
                    "value": "-1",
                    "is_control": "false",
                    "label": "Stromal ECM",
                    "display_value": "negative",
                    "number_role": "inferential_signed_stouffer_direction",
                    "plot_order": "3",
                },
                {
                    "plot_role": "annotation",
                    "group_id": "orphan_direction",
                    "value": "1",
                    "is_control": "false",
                    "label": "Stromal ECM",
                    "display_value": "positive",
                    "number_role": "descriptive_median_slope_direction",
                    "plot_order": "4",
                },
            ],
        )
        self.assertIn(b"inferential_signed_stouffer_direction=negative", payload)
        self.assertIn(b"descriptive_median_slope_direction=positive", payload)

    def test_closure_requires_each_workstream_exactly_once(self):
        project = self.make_project()

        def remove_plan50(rows):
            return [row for row in rows if row["workstream_id"] != "PLAN50"]

        project.rewrite_closure(remove_plan50)
        with self.assertRaisesRegex(ReleaseContractError, "exactly one row"):
            validate_closure(project.root, project.closure_path)

    def test_coverage_limited_is_a_gate_verdict_not_a_terminal_state(self):
        project = self.make_project()

        def drift_plan30_state(rows):
            for row in rows:
                if row["workstream_id"] == "PLAN30":
                    row["terminal_state"] = "coverage_limited"
                    row["gate_verdict"] = "coverage_limited_terminal"
            return rows

        project.rewrite_closure(drift_plan30_state)
        with self.assertRaisesRegex(ReleaseContractError, "unsupported terminal state"):
            validate_closure(project.root, project.closure_path)

    def test_source_gate_skip_cannot_be_recast_as_main_inference(self):
        project = self.make_project()

        def skip_plan13(rows):
            for row in rows:
                if row["workstream_id"] == "PLAN13":
                    row["terminal_state"] = "skipped_source_gate"
                    row["include_main"] = "false"
                    row["exclusion_reason"] = "authoritative source gate failed"
            return rows

        project.rewrite_closure(skip_plan13)
        validate_closure(project.root, project.closure_path)

        def promote_skip(rows):
            for row in rows:
                if row["workstream_id"] == "PLAN13":
                    row["include_main"] = "true"
            return rows

        project.rewrite_closure(promote_skip)
        with self.assertRaisesRegex(ReleaseContractError, "cannot enter main"):
            validate_closure(project.root, project.closure_path)

    def test_source_hash_mismatch_blocks_before_candidate_creation(self):
        project = self.make_project()
        project.source_paths["PLAN13"][0].write_text(
            "mutated after terminal manifest\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ReleaseContractError, "hash/byte mismatch"):
            self.snapshot(project)
        self.assertFalse(project.candidate_root.exists())

    def test_base_selection_requires_exact_signed_hash(self):
        project = self.make_project()
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "selection SHA256 mismatch"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                "0" * 64,
                fixture_mode=True,
            )
        self.assertEqual(list((project.candidate_root / "inputs/BASE").iterdir()), [])

    def test_base_sources_are_all_preflighted_before_first_candidate_write(self):
        project = self.make_project()
        self.snapshot(project)
        project.base_source_paths["sample_composition"].write_text(
            "mutated after coordinator selection\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ReleaseContractError, "pre-write hash/byte"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )
        self.assertFalse(
            (project.candidate_root / "manifests/base_input_selection.tsv").exists()
        )
        self.assertEqual(list((project.candidate_root / "inputs/BASE").iterdir()), [])

    def test_old_source_renamed_to_primary_basename_cannot_replace_sealed_bundle(self):
        project = self.make_project()
        self.snapshot(project)
        stale = (
            project.root
            / "fixture_base_sources/stale/fibrosis_consecutive_lvqw_true_kleiner.csv"
        )
        stale.parent.mkdir()
        stale.write_bytes(
            project.base_source_paths["fibrosis_transitions"].read_bytes()
        )
        rows = read_tsv_exact(project.base_selection_path, BASE_SELECTION_FIELDS)
        for row in rows:
            if row["artifact_id"] == FIBROSIS_PRIMARY_ARTIFACT_ID:
                row["source_path"] = stale.relative_to(project.root).as_posix()
                row["source_sha256"] = sha256_file(stale)
                row["source_bytes"] = str(stale.stat().st_size)
        atomic_write_tsv(project.base_selection_path, rows, BASE_SELECTION_FIELDS)
        with self.assertRaisesRegex(
            ReleaseContractError, "one exact bundle|cannot be reconstructed"
        ):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_fibrosis_artifact_changed_after_validation_blocks_before_copy(self):
        project = self.make_project()
        fibrosis = project.base_source_paths["fibrosis_transitions"]
        with fibrosis.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        rows[0]["padj"] = "0.5"
        with fibrosis.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        selection = read_tsv_exact(project.base_selection_path, BASE_SELECTION_FIELDS)
        for row in selection:
            if row["artifact_id"] == FIBROSIS_PRIMARY_ARTIFACT_ID:
                row["source_sha256"] = sha256_file(fibrosis)
                row["source_bytes"] = str(fibrosis.stat().st_size)
        atomic_write_tsv(project.base_selection_path, selection, BASE_SELECTION_FIELDS)
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "changed after validation"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )
        self.assertFalse(
            (project.candidate_root / "manifests/base_input_selection.tsv").exists()
        )

    def test_fibrosis_missing_ready_blocks_before_copy(self):
        project = self.make_project()
        (project.fibrosis_bundle_root / "READY").unlink()
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "file universe is not exact"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_fibrosis_validation_report_hash_drift_blocks_before_copy(self):
        project = self.make_project()
        report = project.fibrosis_bundle_root / "manifests/validation_report.json"
        report.write_bytes(report.read_bytes() + b" ")
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "changed after validation"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_fibrosis_wrong_candidate_identity_blocks_even_when_selection_hashes_match(
        self,
    ):
        project = self.make_project()
        ready = project.fibrosis_bundle_root / "READY"
        payload = json.loads(ready.read_text(encoding="utf-8"))
        payload["candidate_id"] = "stale-candidate"
        ready.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
        project.reseal_fibrosis_bundle()
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "READY candidate_id drift"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_fibrosis_wrong_census_blocks_after_full_reseal(self):
        project = self.make_project()
        manifest = project.fibrosis_bundle_root / "audits/sample_manifest.tsv"
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows.pop()
        atomic_write_tsv(manifest, rows, tuple(rows[0]))
        project.reseal_fibrosis_bundle()
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "row census drift"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_fibrosis_duplicate_biological_donor_blocks_after_full_reseal(self):
        project = self.make_project()
        manifest = project.fibrosis_bundle_root / "audits/sample_manifest.tsv"
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows[1]["donor_id"] = rows[0]["donor_id"]
        atomic_write_tsv(manifest, rows, tuple(rows[0]))
        project.reseal_fibrosis_bundle()
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "biological donor"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_fibrosis_second_biopsy_blocks_after_full_reseal(self):
        project = self.make_project()
        manifest = project.fibrosis_bundle_root / "audits/sample_manifest.tsv"
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        gse = next(row for row in rows if row["dataset"] == "GSE193066")
        gse["biopsy"] = "2nd biopsy"
        atomic_write_tsv(manifest, rows, tuple(rows[0]))
        project.reseal_fibrosis_bundle()
        self.snapshot(project)
        with self.assertRaisesRegex(ReleaseContractError, "eligible first biopsy"):
            stage_base_inputs(
                project.root,
                project.base_selection_path,
                sha256_file(project.base_selection_path),
                fixture_mode=True,
            )

    def test_valid_zero_significant_fibrosis_family_builds_without_positivity_gate(
        self,
    ):
        project = self.make_project()
        self.snapshot(project)
        self.build_products(project)
        numbers = read_tsv_exact(
            project.candidate_root / "tables/numbers_ledger.tsv", NUMBERS_FIELDS
        )
        stage_counts = [
            row["raw_value"] for row in numbers if row["number_role"] == "deg_count"
        ]
        self.assertEqual(stage_counts, ["0", "0", "0", "0"])

    def test_real_fibrosis_bundle_is_exact_27638_by_four_and_hash_pinned(self):
        project_root = SCRIPT_DIR.parents[2]
        bundle = (
            project_root
            / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates"
            / CANDIDATE_ID
            / "fibrosis-adjacent-true-kleiner-lvqw-v1"
        )
        validated = validate_fibrosis_candidate_bundle(
            bundle,
            expected_n_genes=FIBROSIS_EXPECTED_N_GENES,
            project_root=project_root,
            verify_external_sources=True,
        )
        self.assertEqual(validated["n_genes_per_contrast"], 27_638)
        self.assertEqual(
            validated["primary_result_sha256"],
            FIBROSIS_ACCEPTED_PRIMARY_RESULT_SHA256,
        )
        with tempfile.TemporaryDirectory() as temporary:
            replacement = Path(temporary) / "self_consistent_replacement"
            shutil.copytree(bundle, replacement)
            ready = replacement / "READY"
            ready.write_bytes(ready.read_bytes() + b" ")
            with self.assertRaisesRegex(
                ReleaseContractError, "differs from accepted upstream release: READY"
            ):
                validate_fibrosis_candidate_bundle(
                    replacement,
                    expected_n_genes=FIBROSIS_EXPECTED_N_GENES,
                )

    def test_nonconfirmatory_result_cannot_be_recast_tested_negative(self):
        project = self.make_project()
        source = next(
            path
            for path in project.source_paths["PLAN13"]
            if path.name == "spatial_evidence.tsv"
        )
        rows = read_tsv_exact(source, SOURCE_ROW_FIELDS[:-2])
        rows[0]["evidence_status"] = "tested_negative"
        rows[0]["allowed_wording"] = "nonconfirmatory result"
        atomic_write_tsv(source, rows, SOURCE_ROW_FIELDS[:-2])

        def update_manifest(manifest_rows):
            for row in manifest_rows:
                if row["artifact_id"] == "spatial_evidence":
                    row["source_sha256"] = sha256_file(source)
                    row["source_bytes"] = str(source.stat().st_size)
            return manifest_rows

        project.rewrite_artifact_manifest("PLAN13", update_manifest)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(ReleaseContractError, "adequate-negative"):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_overlapping_sources_cannot_be_labeled_independent(self):
        project = self.make_project()

        def mutate(rows):
            rows[0]["source_dependence"] = "independent"
            rows[0]["discovery_sources"] = "shared_source"
            rows[0]["evaluation_sources"] = "shared_source"
            rows[0]["reuse_detail"] = "The same source appears on both sides."
            rows[0]["independence_boundary"] = "No valid independence boundary."
            return rows

        self.mutate_evidence(project, "PLAN13", "spatial_evidence.tsv", mutate)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(
            ReleaseContractError, "overlapping sources independent"
        ):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_yak_module8_directions_cannot_be_collapsed(self):
        project = self.make_project()
        source = next(
            path
            for path in project.source_paths["PLAN13"]
            if path.name == "spatial_evidence.tsv"
        )
        rows = read_tsv_exact(source, SOURCE_ROW_FIELDS[:-2])
        rows = [
            row
            for row in rows
            if row["number_role"] != "inferential_signed_stouffer_direction"
        ]
        atomic_write_tsv(source, rows, SOURCE_ROW_FIELDS[:-2])

        def update_manifest(manifest_rows):
            for row in manifest_rows:
                if row["artifact_id"] == "spatial_evidence":
                    row["source_sha256"] = sha256_file(source)
                    row["source_bytes"] = str(source.stat().st_size)
            return manifest_rows

        project.rewrite_artifact_manifest("PLAN13", update_manifest)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(ReleaseContractError, "keep descriptive"):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_array_count_cannot_masquerade_as_biological_donor_n(self):
        project = self.make_project()

        def mutate(rows):
            for row in rows:
                if row["record_id"] == "vu_array_module1":
                    row["unit"] = "biological donors"
                    row["biological_unit"] = "10 donors"
            return rows

        self.mutate_evidence(project, "PLAN13", "spatial_evidence.tsv", mutate)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(ReleaseContractError, "technical units"):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_moran_null_dispersion_cannot_be_sampling_se(self):
        project = self.make_project()

        def mutate(rows):
            for row in rows:
                if row["record_id"] == "spatial_moran_module1":
                    row["effect_unit"] = "sampling SE"
            return rows

        self.mutate_evidence(project, "PLAN13", "spatial_evidence.tsv", mutate)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(ReleaseContractError, "matched-null Moran"):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_myojin_precise_null_language_is_rejected(self):
        project = self.make_project()

        def mutate(rows):
            rows[0]["allowed_wording"] = "assay-specific non-support; precise null"
            return rows

        self.mutate_evidence(project, "PLAN40", "myojin_supplement.tsv", mutate)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(ReleaseContractError, "overstates the Myojin"):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_nmf_stable_molecular_subtyping_paraphrase_is_rejected(self):
        project = self.make_project()

        def mutate(rows):
            for row in rows:
                if row["section_id"] == "supplementary_nmf":
                    row["allowed_wording"] = (
                        "continuous NMF axes define stable molecular subtyping"
                    )
                    break
            return rows

        self.mutate_evidence(project, "PLAN20", "nmf_continuous_supplement.tsv", mutate)
        self.snapshot(project)
        stage_base_inputs(
            project.root,
            project.base_selection_path,
            sha256_file(project.base_selection_path),
            fixture_mode=True,
        )
        with self.assertRaisesRegex(ReleaseContractError, "approved NMF language"):
            build_candidate_tables(project.root, fixture_mode=True)

    def test_symlinked_source_is_rejected_even_when_content_matches(self):
        project = self.make_project()
        source = project.source_paths["PLAN13"][0]
        target = source.with_name("same_content_target.tsv")
        target.write_bytes(source.read_bytes())
        source.unlink()
        source.symlink_to(target.name)
        with self.assertRaisesRegex(ReleaseContractError, "non-symlink"):
            validate_closure(project.root, project.closure_path)

    def test_mutable_downstream_path_is_rejected_at_closure(self):
        project = self.make_project()

        def mutate(rows):
            rows[0]["downstream_read_path"] = "figures/main/live_panel.tsv"
            return rows

        project.rewrite_artifact_manifest("PLAN13", mutate)
        with self.assertRaisesRegex(ReleaseContractError, "copied snapshot"):
            validate_closure(project.root, project.closure_path)

    def test_existing_candidate_with_different_spec_is_refused(self):
        project = self.make_project()
        self.snapshot(project)

        def mutate(rows):
            rows[0]["allowed_wording"] = "different signed wording"
            return rows

        project.rewrite_closure(mutate)
        with self.assertRaisesRegex(ReleaseContractError, "spec differs"):
            self.snapshot(project)

    def test_existing_generated_manifest_is_never_repaired_in_place(self):
        project = self.make_project()
        self.snapshot(project)
        manifest = project.candidate_root / "manifests/input_snapshot_manifest.tsv"
        manifest.write_text("mutated generated manifest\n", encoding="utf-8")
        with self.assertRaisesRegex(ReleaseContractError, "refusing replacement"):
            self.snapshot(project)

    def test_corrupt_existing_copy_is_never_overwritten_on_resume(self):
        project = self.make_project()
        self.snapshot(project)
        frozen = project.candidate_root / "inputs/SP-INT/plan13_terminal.tsv"
        frozen.write_text("corrupted candidate copy\n", encoding="utf-8")
        with self.assertRaisesRegex(ReleaseContractError, "existing snapshot differs"):
            self.snapshot(project)

    def test_protected_named_release_mutation_blocks_consistency(self):
        project = self.make_project()
        self.snapshot(project)
        self.build_products(project)
        project.install_scientific_registry()
        protected = (
            project.root
            / "RNA-seq/results/manuscript_release/2026-01-01-r1/manifest.tsv"
        )
        protected.write_text("release\tvalue\nr1\tmutated\n", encoding="utf-8")
        with self.assertRaisesRegex(ReleaseContractError, "protected release drift"):
            validate_candidate(project.root, fixture_mode=True, write_report=False)

    def test_new_unregistered_named_release_blocks_snapshot(self):
        project = self.make_project()
        unexpected = project.root / "RNA-seq/results/manuscript_release/2026-03-03-r3"
        unexpected.mkdir()
        (unexpected / "manifest.tsv").write_text(
            "release\tvalue\nr3\tunregistered\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ReleaseContractError, "every pre-existing named"):
            self.snapshot(project)

    def test_missing_scientific_registry_never_emits_pass(self):
        project = self.make_project()
        self.snapshot(project)
        self.build_products(project)
        with self.assertRaisesRegex(ReleaseContractError, "missing or symlinked"):
            validate_candidate(project.root, fixture_mode=True, write_report=False)

    def test_scientific_registry_requires_every_check_once(self):
        project = self.make_project()
        self.snapshot(project)
        self.build_products(project)
        project.install_scientific_registry(omit_check="null_compatibility")
        with self.assertRaisesRegex(ReleaseContractError, "every required check"):
            validate_candidate(project.root, fixture_mode=True, write_report=False)

    def test_warning_without_written_adjudication_blocks(self):
        project = self.make_project()
        self.snapshot(project)
        self.build_products(project)
        project.install_scientific_registry(
            warning_check="numbers_claim_ledger",
            adjudicate_warning=False,
        )
        with self.assertRaisesRegex(ReleaseContractError, "lack written adjudication"):
            validate_candidate(project.root, fixture_mode=True, write_report=False)

    def test_candidate_root_guard_rejects_canonical_paths(self):
        project = self.make_project()
        with self.assertRaisesRegex(ReleaseContractError, "must be exactly"):
            assert_candidate_root(project.root, project.root / "figures/main")
        self.assertEqual(
            assert_candidate_root(project.root, candidate_root(project.root)),
            project.candidate_root,
        )
        self.assertEqual(project.candidate_root.name, CANDIDATE_ID)

    def test_candidate_root_symlink_is_rejected_before_any_copy(self):
        project = self.make_project()
        project.candidate_root.parent.mkdir(parents=True, exist_ok=True)
        detour = project.root / "candidate_detour"
        detour.mkdir()
        project.candidate_root.symlink_to(detour)
        with self.assertRaisesRegex(ReleaseContractError, "non-symlink"):
            self.snapshot(project)
        self.assertEqual(list(detour.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
