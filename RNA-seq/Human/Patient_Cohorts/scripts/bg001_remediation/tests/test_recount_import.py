#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import importlib.util
import inspect
import io
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT_DIR = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


importer = load_module("bg001_import_recount", SCRIPT_DIR / "import_recount.py")
finalizer = load_module("bg001_finalize_recount", SCRIPT_DIR / "finalize_recount.py")
verifier = load_module("bg001_verify_recount", SCRIPT_DIR / "verify_recount_artifacts.py")


def make_dataset_validation_fixture(root: Path, dataset: str = "GSE130970"):
    root = root.resolve()
    dataset_root = root / "counts" / dataset
    dataset_root.mkdir(parents=True)
    command = (
        '# Program:featureCounts v2.1.1; Command:"featureCounts" "-T" "8" '
        '"-p" "--countReadPairs" "-B" "-s" "2" "-a" "/gencode.v49.gtf.gz"'
    )
    paths = {
        "count_path": dataset_root / "gene_counts.txt",
        "summary_path": dataset_root / "gene_counts.txt.summary",
        "log_path": dataset_root / "featureCounts.log",
        "environment_path": dataset_root / "environment.txt",
    }
    paths["count_path"].write_text(command + "\nGeneid\n")
    paths["summary_path"].write_text("Status\n")
    paths["log_path"].write_text("paired fragments\n")
    paths["environment_path"].write_text("environment\n")
    samples = [f"S{index:03d}" for index in range(verifier.EXPECTED_SAMPLES[dataset])]
    payload = {
        "status": "PASS",
        "dataset": dataset,
        "n_samples": len(samples),
        "n_genes": 86_369,
        "featurecounts_version": "2.1.1",
        "gtf_sha256": verifier.EXPECTED_GTF_SHA256,
        "samples": samples,
        "command": command,
        "counts_sha256": verifier.sha256(paths["count_path"]),
        "summary_sha256": verifier.sha256(paths["summary_path"]),
        "log_sha256": verifier.sha256(paths["log_path"]),
        "environment_sha256": verifier.sha256(paths["environment_path"]),
    }
    validation = dataset_root / "validation.json"
    validation.write_text(json.dumps(payload, sort_keys=True) + "\n")
    row = {field: str(path) for field, path in paths.items()}
    row.update(
        {
            "count_sha256": payload["counts_sha256"],
            "summary_sha256": payload["summary_sha256"],
            "log_sha256": payload["log_sha256"],
            "environment_sha256": payload["environment_sha256"],
            "validation_sha256": verifier.sha256(validation),
        }
    )
    return row, [command], samples, validation


def make_task_checksum_fixture(root: Path):
    manifest_rows = {}
    for dataset, count in verifier.EXPECTED_SAMPLES.items():
        manifest_rows[dataset] = [
            {
                "dataset": dataset,
                "sample_id": f"{dataset}_S{index:03d}",
                "bam_path": f"/bam/{dataset}_S{index:03d}.bam",
                "layout": "paired",
                "strandedness": "reverse",
                "expected_status": "included",
                "exclusion_reason": "",
                "bam_size": str(1000 + index),
                "bam_mtime": str(2000 + index),
                "bam_sha256": f"{index % 16:x}" * 64,
            }
            for index in range(count)
        ]

    contracts = finalizer.task_checksum_contracts(root, manifest_rows)
    for task_dir, expected_rows in contracts:
        task_dir.mkdir(parents=True, exist_ok=True)
        with (task_dir / "bam_checksums.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=finalizer.CHECKSUM_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(expected_rows)
    return manifest_rows


class ProducerCompatibilityTests(unittest.TestCase):
    def test_reviewed_legacy_to_guarded_checker_transition_reproduces_exactly(self):
        project = SCRIPT_DIR.parents[4]
        source_root = (
            project / "results/remediation/bg001/"
            "bg001-fragment-v211-gencode49-20260806T055119Z"
        )
        source_checker = source_root / "source_snapshot" / importer.REGRESSION_CHECKER_REL
        destination_checker = project / importer.REGRESSION_CHECKER_REL
        if not source_checker.is_file():
            self.skipTest("sealed execution-source checker is unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            destination_root = Path(tmp) / "destination"
            destination_root.mkdir()
            (destination_root / "source_snapshot").symlink_to(project, target_is_directory=True)
            source_manifest = {
                importer.REGRESSION_CHECKER_REL: {
                    "size": source_checker.stat().st_size,
                    "sha256": importer.sha256(source_checker),
                }
            }
            destination_manifest = {
                importer.REGRESSION_CHECKER_REL: {
                    "size": destination_checker.stat().st_size,
                    "sha256": importer.sha256(destination_checker),
                }
            }
            transition = importer.compare_regression_checker(
                source_root, destination_root, source_manifest, destination_manifest
            )
        self.assertEqual(transition["execution_generation"], "legacy_055119Z")
        self.assertEqual(transition["destination_generation"], "guarded_routes_v1")
        self.assertEqual(
            {
                key: transition["destination_checker_on_execution_snapshot"][key]
                for key in importer.GUARDED_CHECKER_ON_LEGACY_EXPECTED
            },
            importer.GUARDED_CHECKER_ON_LEGACY_EXPECTED,
        )

    def test_destination_guard_mitigation_bundle_is_exactly_source_manifest_bound(self):
        project = SCRIPT_DIR.parents[4]
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp)
            manifest = {}
            for relative in importer.DESTINATION_GUARD_MITIGATION_PATHS:
                source = project / relative
                target = destination / "source_snapshot" / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                manifest[relative] = {
                    "size": target.stat().st_size,
                    "sha256": importer.sha256(target),
                }
            rows = importer.validate_destination_guard_mitigations(destination, manifest)
            self.assertEqual(len(rows), 4)
            changed = destination / "source_snapshot" / importer.DESTINATION_GUARD_MITIGATION_PATHS[0]
            changed.write_bytes(changed.read_bytes() + b"\n")
            with self.assertRaisesRegex(SystemExit, "differs from source manifest"):
                importer.validate_destination_guard_mitigations(destination, manifest)

    def test_only_exact_producer_bundle_is_reusable(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source = base / "source"
            destination = base / "destination"
            source_manifest = {}
            destination_manifest = {}
            for root in (source, destination):
                (root / "source_snapshot").mkdir(parents=True)
            for relative in importer.PRODUCER_PATHS:
                payload = f"producer:{relative}\n".encode()
                digest = hashlib.sha256(payload).hexdigest()
                for root in (source, destination):
                    path = root / "source_snapshot" / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(payload)
                record = {"size": len(payload), "sha256": digest}
                source_manifest[relative] = dict(record)
                destination_manifest[relative] = dict(record)
            rows = importer.compare_producer_manifests(
                source, destination, source_manifest, destination_manifest
            )
            self.assertEqual(len(rows), len(importer.PRODUCER_PATHS))

            changed = importer.PRODUCER_PATHS[0]
            destination_manifest[changed]["sha256"] = "0" * 64
            with self.assertRaisesRegex(SystemExit, "producer contract differs"):
                importer.compare_producer_manifests(
                    source, destination, source_manifest, destination_manifest
                )

    def test_nonproducer_source_differences_are_recorded(self):
        source = {"protected.tsv": {"size": 1, "sha256": "a" * 64}}
        destination = {"protected.tsv": {"size": 2, "sha256": "b" * 64}}
        rows = importer.source_manifest_differences(source, destination)
        self.assertEqual(rows[0]["status"], "changed")
        self.assertEqual(rows[0]["execution_source_sha256"], "a" * 64)
        self.assertEqual(rows[0]["destination_analysis_sha256"], "b" * 64)


class DigestSafetyTests(unittest.TestCase):
    def test_digest_path_escape_is_rejected_before_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "seal.sha256"
            manifest.write_text(f"{'0' * 64}\t../escape\n")
            with self.assertRaisesRegex(SystemExit, "Escaping or duplicate"):
                importer.parse_digest_manifest(root, manifest)

    def test_duplicate_archived_digest_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "seal.sha256"
            manifest.write_text(f"{'a' * 64}\tcounts/x\n{'b' * 64}\tcounts/x\n")
            with self.assertRaisesRegex(SystemExit, "Escaping or duplicate"):
                importer.read_digest_records(manifest)

    def test_hardlinked_copy_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            target = root / "target"
            source.write_text("sealed bytes\n")

            def hardlink_copy(source_path, target_path, **_kwargs):
                os.link(source_path, target_path)
                return str(target_path)

            with mock.patch.object(importer.shutil, "copy2", side_effect=hardlink_copy):
                with self.assertRaisesRegex(SystemExit, "hard-linked"):
                    importer.copy_independent_regular(
                        source,
                        target,
                        importer.sha256(source),
                        "test source",
                    )


class AtomicPublicationSafetyTests(unittest.TestCase):
    def test_finalize_rejects_legacy_temp_symlink_without_following_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp).resolve()
            root = base / "candidate"
            manifests = root / "manifests"
            manifests.mkdir(parents=True)
            target = manifests / "count_manifest.tsv"
            legacy_temp = target.with_suffix(target.suffix + ".tmp")
            outside = base / "outside.tsv"
            outside.write_text("outside must remain unchanged\n")
            legacy_temp.symlink_to(outside)

            with self.assertRaisesRegex(SystemExit, "existing output or temporary"):
                finalizer.atomic_write_tsv(
                    target,
                    ("dataset", "count_path"),
                    ({"dataset": "GSE1", "count_path": "/candidate/counts"},),
                    root=root,
                )
            self.assertEqual(outside.read_text(), "outside must remain unchanged\n")
            self.assertFalse(target.exists() or target.is_symlink())

    def test_recount_seal_rejects_legacy_temp_symlink_without_following_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp).resolve()
            root = base / "candidate"
            comparisons = root / "comparisons"
            comparisons.mkdir(parents=True)
            target = comparisons / "structural_validation.json"
            legacy_temp = target.with_suffix(".json.tmp")
            outside = base / "outside.json"
            outside.write_text("outside must remain unchanged\n")
            legacy_temp.symlink_to(outside)

            with self.assertRaisesRegex(SystemExit, "existing output or temporary"):
                verifier.atomic_write_text(
                    target,
                    '{"status":"PASS"}\n',
                    root=root,
                    legacy_temporary=legacy_temp,
                )
            self.assertEqual(outside.read_text(), "outside must remain unchanged\n")
            self.assertFalse(target.exists() or target.is_symlink())

    def test_publishers_create_regular_single_link_outputs_and_remove_uuid_temps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            manifests = root / "manifests"
            comparisons = root / "comparisons"
            manifests.mkdir()
            comparisons.mkdir()

            table = manifests / "count_manifest.tsv"
            finalizer.atomic_write_tsv(
                table,
                ("dataset", "count_path"),
                ({"dataset": "GSE1", "count_path": "/candidate/counts"},),
                root=root,
            )
            report = comparisons / "structural_validation.json"
            verifier.atomic_write_text(
                report,
                '{"status":"PASS"}\n',
                root=root,
                legacy_temporary=report.with_suffix(".json.tmp"),
            )

            for path in (table, report):
                self.assertTrue(path.is_file())
                self.assertFalse(path.is_symlink())
                self.assertEqual(path.stat().st_nlink, 1)
            self.assertFalse(list(root.rglob(".*.tmp.*")))


class TaskChecksumAttributionTests(unittest.TestCase):
    def test_exact_task_slices_pass_both_finalizer_and_independent_verifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            manifest_rows = make_task_checksum_fixture(root)
            observed = finalizer.validate_task_checksum_tables(root, manifest_rows)
            self.assertEqual(len(observed), sum(verifier.EXPECTED_SAMPLES.values()))
            verifier.validate_task_checksum_tables(root, manifest_rows)

    def test_cross_chunk_checksum_swap_is_rejected_by_both_layers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            manifest_rows = make_task_checksum_fixture(root)
            chunk0 = root / "counts/GSE213621/chunks/chunk0/bam_checksums.tsv"
            chunk1 = root / "counts/GSE213621/chunks/chunk1/bam_checksums.tsv"
            chunk0_text, chunk1_text = chunk0.read_text(), chunk1.read_text()
            chunk0.write_text(chunk1_text)
            chunk1.write_text(chunk0_text)

            for module in (finalizer, verifier):
                with self.subTest(module=module.__name__):
                    with self.assertRaisesRegex(SystemExit, "frozen manifest slice"):
                        module.validate_task_checksum_tables(root, manifest_rows)

    def test_checksum_schema_drift_is_rejected_by_both_layers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            manifest_rows = make_task_checksum_fixture(root)
            path = root / "counts/GSE130970/bam_checksums.tsv"
            lines = path.read_text().splitlines()
            lines[0] += "\textra"
            path.write_text("\n".join(lines) + "\n")

            for module in (finalizer, verifier):
                with self.subTest(module=module.__name__):
                    with self.assertRaisesRegex(SystemExit, "checksum schema"):
                        module.validate_task_checksum_tables(root, manifest_rows)


class MetadataRebaseTests(unittest.TestCase):
    def test_provenance_rebase_preserves_payload_hashes_and_uses_destination_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / "stage"
            stage.mkdir()
            for name in importer.PROVENANCE_PAYLOADS:
                (stage / name).write_text(f"payload:{name}\n")
            final = Path(tmp) / "destination/counts/GSE130970"
            importer.write_rebased_provenance(stage, final)
            lines = (stage / "provenance.sha256").read_text().splitlines()
            self.assertEqual(len(lines), len(importer.PROVENANCE_PAYLOADS))
            for line, name in zip(lines, importer.PROVENANCE_PAYLOADS):
                digest, path = line.split(maxsplit=1)
                self.assertEqual(digest, importer.sha256(stage / name))
                self.assertEqual(path, str(final / name))

    def test_merge_publication_rebase_changes_only_run_bound_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)
            for name in importer.MERGE_PAYLOADS:
                (stage / name).write_text(f"merged:{name}\n")
            run_id = "bg001-fragment-v211-gencode49-20260807T000000Z"
            importer.write_rebased_merge_publication(stage, run_id)
            publication = json.loads((stage / "merge_publication.json").read_text())
            self.assertEqual(publication["run_id"], run_id)
            self.assertEqual(publication["publish_order"], list(importer.MERGE_PUBLISH_ORDER))
            for name in importer.MERGE_PAYLOADS:
                self.assertEqual(publication["hashes"][name], importer.sha256(stage / name))
            self.assertEqual(
                (stage / "MERGE_COMPLETE").read_text().strip(),
                importer.sha256(stage / "merge_publication.json"),
            )

    def test_copy_archives_original_path_metadata_and_does_not_mutate_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source = base / "source"
            destination = base / "destination"
            staging = destination / ".stage"
            source.mkdir()
            destination.mkdir()
            staging.mkdir()
            task_relatives = [
                f"counts/{dataset}"
                for dataset in ("GSE130970", "GSE135251", "GSE174478", "GSE240729")
            ] + [f"counts/GSE213621/chunks/chunk{index}" for index in range(4)]
            for relative in task_relatives:
                task = source / relative
                task.mkdir(parents=True)
                for name in importer.PROVENANCE_PAYLOADS:
                    (task / name).write_text(f"{relative}:{name}\n")
                (task / "provenance.sha256").write_text(f"source-only:{relative}\n")
                (task / "COMPLETE").touch()
            merged = source / "counts/GSE213621"
            for name in importer.MERGE_PAYLOADS:
                (merged / name).write_text(f"merged:{name}\n")
            (merged / "merge_publication.json").write_text(
                json.dumps({"run_id": "source", "hashes": {}}, sort_keys=True) + "\n"
            )
            (merged / "MERGE_COMPLETE").write_text("source-marker\n")
            source_manifest = source / "manifests/bam_manifest.tsv"
            source_manifest.parent.mkdir(parents=True)
            source_manifest.write_text("source-manifest\n")
            sealed = {
                str(path.relative_to(source)): importer.sha256(path)
                for path in source.rglob("*")
                if path.is_file()
            }
            source_hashes_before = dict(sealed)
            _, lineage, inventory = importer.copy_and_rebase_counts(
                source, destination, sealed, staging
            )
            self.assertEqual(
                {
                    str(path.relative_to(source)): importer.sha256(path)
                    for path in source.rglob("*")
                    if path.is_file()
                },
                source_hashes_before,
            )
            by_source = {row["source_relative_path"]: row for row in inventory}
            for relative in importer.PATH_BOUND_COUNT_FILES:
                archived = lineage / "source_originals" / relative
                self.assertEqual(importer.sha256(archived), sealed[relative])
                self.assertEqual(
                    by_source[relative]["disposition"], "archived_original_and_rebased"
                )
                self.assertNotEqual(
                    by_source[relative]["effective_sha256"], sealed[relative]
                )
            manifest_row = by_source["manifests/bam_manifest.tsv"]
            self.assertEqual(manifest_row["disposition"], "archived_original")
            self.assertEqual(
                manifest_row["preserved_destination_relative_path"],
                "manifests/recount_import/source_originals/manifests/bam_manifest.tsv",
            )


class ProtectedExpansionTests(unittest.TestCase):
    def test_protected_gene_reasons_and_sources_may_only_expand(self):
        old = {
            "ENSG1": {
                "symbol": "GENE1",
                "reasons": {"canonical_treat"},
                "sources": {"canonical.csv"},
            }
        }
        expanded = {
            "ENSG1": {
                "symbol": "GENE1",
                "reasons": {"canonical_treat", "active_curated_figure"},
                "sources": {"canonical.csv", "labels.tsv"},
            },
            "ENSG2": {
                "symbol": "GENE2",
                "reasons": {"active_curated_figure"},
                "sources": {"labels.tsv"},
            },
        }
        importer.require_monotonic_protected_genes(old, expanded)
        narrowed = {"ENSG1": {"symbol": "GENE1", "reasons": set(), "sources": set()}}
        with self.assertRaisesRegex(SystemExit, "not monotonic"):
            importer.require_monotonic_protected_genes(old, narrowed)
        with self.assertRaisesRegex(SystemExit, "removed protected genes"):
            importer.require_monotonic_protected_genes(old, {})

    def test_read_count_manifests_are_compared_after_run_path_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            normalized = []
            for run_name in ("source", "destination"):
                root = base / run_name
                manifests = root / "manifests"
                manifests.mkdir(parents=True)
                source_lines = [
                    "dataset\tsource_count_path\tfrozen_count_path\t"
                    "source_summary_path\tfrozen_summary_path\n"
                ]
                override_lines = ["dataset\tcount_path\n"]
                for index in range(9):
                    dataset = f"GSE{index}"
                    frozen = root / "frozen_sets/read_counts" / dataset
                    frozen.mkdir(parents=True)
                    count = frozen / "gene_counts.txt"
                    summary = frozen / "gene_counts.txt.summary"
                    count.write_text("count\n")
                    summary.write_text("summary\n")
                    source_lines.append(
                        f"{dataset}\t/live/{dataset}/gene_counts.txt\t{count}\t"
                        f"/live/{dataset}/gene_counts.txt.summary\t{summary}\n"
                    )
                    override_lines.append(f"{dataset}\t{count}\n")
                (manifests / "read_count_sources.tsv").write_text("".join(source_lines))
                (manifests / "read_count_overrides.tsv").write_text("".join(override_lines))
                normalized.append(importer.normalized_read_count_manifests(root))
            self.assertEqual(normalized[0], normalized[1])


class GuardedBaselineTransitionTests(unittest.TestCase):
    def test_legacy_source_guard_transition_is_explicit_and_never_synthesized(self):
        destination_records = {
            relative: hashlib.sha256(relative.encode()).hexdigest()
            for relative in (
                importer.BAM_GUARD_BASELINE_INPUTS | importer.RUNTIME_BASELINE_INPUTS
            )
        }
        with (
            mock.patch.object(
                importer,
                "verify_baseline",
                side_effect=[({}, "legacy_zero_marker"), (destination_records, "structured_v1")],
            ),
            mock.patch.object(importer, "normalized_read_count_manifests", return_value={}),
            mock.patch.object(importer, "normalized_protected_sources", return_value={}),
            mock.patch.object(importer, "protected_gene_rows", return_value={}),
            mock.patch.object(importer, "not_testable_rows", return_value=set()),
        ):
            rows = importer.compare_baselines(Path("/source"), Path("/destination"))
        guard_rows = [
            row for row in rows
            if row["relative_path"] in importer.BAM_GUARD_BASELINE_INPUTS
        ]
        self.assertEqual(
            {row["change_class"] for row in guard_rows},
            {"legacy_source_no_guard_to_guarded_destination"},
        )
        self.assertEqual(
            {row["relative_path"] for row in guard_rows},
            importer.BAM_GUARD_BASELINE_INPUTS,
        )

    def test_legacy_baseline_rejects_unsealed_guard_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "BASELINE_FROZEN").touch()
            guard_path = root / next(iter(importer.BAM_GUARD_BASELINE_INPUTS))
            guard_path.parent.mkdir()
            guard_path.write_text('{}\n')
            with self.assertRaisesRegex(SystemExit, "unsealed/synthesized BAM guard"):
                importer.baseline_generation(root)

    def test_recount_authoritative_seal_includes_both_guard_records(self):
        source = inspect.getsource(verifier.authoritative_manifest_files)
        for relative in importer.BAM_GUARD_BASELINE_INPUTS:
            self.assertIn(Path(relative).name, source)


class ImportAwareVerifierTests(unittest.TestCase):
    def test_manifest_schema_separates_execution_and_analysis_snapshots(self):
        fields = verifier.COUNT_MANIFEST_FIELDS
        self.assertIn("artifact_origin", fields)
        self.assertIn("execution_source_run_id", fields)
        self.assertIn("execution_source_manifest_sha256", fields)
        self.assertIn("destination_analysis_source_manifest_sha256", fields)
        self.assertLess(fields.index("source_manifest_sha256"), fields.index("artifact_origin"))

    def test_analysis_runner_blocks_incomplete_import_and_reverifies_lineage(self):
        script = (SCRIPT_DIR / "run_four_arms.sh").read_text()
        self.assertIn("RECOUNT_IMPORT_IN_PROGRESS", script)
        self.assertIn("RECOUNT_IMPORT_FAILED.json", script)
        self.assertIn('import_recount.py" --verify-import', script)

    def test_recount_symlink_is_rejected_before_artifact_enumeration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "counts").mkdir()
            (root / "manifests").mkdir()
            target = root / "outside-count"
            target.write_text("count\n")
            (root / "counts/gene_counts.txt").symlink_to(target)
            with self.assertRaisesRegex(SystemExit, "contains symlink"):
                verifier.reject_recount_symlinks(root)

    def test_swapped_cohort_artifact_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            row, commands, samples, _ = make_dataset_validation_fixture(root)
            row["count_path"] = str(root / "counts/GSE135251/gene_counts.txt")
            with self.assertRaisesRegex(SystemExit, "exact GSE130970 destination"):
                verifier.validate_dataset_artifacts(
                    root, "GSE130970", row, commands, samples
                )

    def test_swapped_validation_payload_is_rejected_even_when_rehashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            row, commands, samples, validation = make_dataset_validation_fixture(root)
            payload = json.loads(validation.read_text())
            payload["dataset"] = "GSE135251"
            payload["n_samples"] = 216
            validation.write_text(json.dumps(payload, sort_keys=True) + "\n")
            row["validation_sha256"] = verifier.sha256(validation)
            with self.assertRaisesRegex(SystemExit, "validation.json: dataset"):
                verifier.validate_dataset_artifacts(
                    root, "GSE130970", row, commands, samples
                )

    def test_check_only_returns_before_any_import_write(self):
        prepared = (
            Path("/source"),
            Path("/destination"),
            {"source_manifest_sha256": "a" * 64},
            {"source_manifest_sha256": "b" * 64},
            {},
            [],
            {
                "execution_generation": "legacy_055119Z",
                "destination_generation": "slurm_directive_aware",
            },
            [],
            [],
            [],
            {"bam_manifest_sha256": "c" * 64},
            {
                "raw_sacct_psv": "scheduler evidence\n",
                "validation": {
                    "status": "PASS",
                    "schema_version": "bg001-source-slurm-attestation-v1",
                    "ledger_sha256": "d" * 64,
                    "raw_sacct_sha256": "e" * 64,
                    "logical_jobs": 20,
                    "completed_jobs": 19,
                    "cancelled_unstarted_jobs": 1,
                },
            },
            {},
            {},
        )
        with (
            mock.patch.object(importer, "prepare_import", return_value=prepared),
            mock.patch.object(importer, "baseline_generation", return_value="legacy_zero_marker"),
            mock.patch.object(
                importer,
                "atomic_write_json",
                side_effect=AssertionError("check-only attempted a write"),
            ),
        ):
            result = importer.perform_import(
                Path("/source"), Path("/destination"), check_only=True
            )
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["execution_source_run_id"], "source")

    def test_scheduler_attestation_is_archived_independently_and_hashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            lineage = root / "lineage"
            (source / "contract").mkdir(parents=True)
            lineage.mkdir()
            ledger = source / "contract/slurm_submission.tsv"
            ledger.write_text("stage\tjob_id\nrecount\t123\n")
            raw = "123|featureCounts|COMPLETED|0:0\n"
            attestation = {
                "raw_sacct_psv": raw,
                "validation": {
                    "status": "PASS",
                    "schema_version": "bg001-source-slurm-attestation-v1",
                    "ledger_sha256": importer.sha256(ledger),
                    "raw_sacct_sha256": importer.sha256_bytes(raw.encode()),
                    "logical_jobs": 20,
                    "completed_jobs": 19,
                    "cancelled_unstarted_jobs": 1,
                },
            }
            result = importer.archive_scheduler_attestation(
                source, lineage, attestation
            )
            scheduler = lineage / "slurm"
            self.assertEqual(
                (scheduler / "source_submission_ledger.tsv").read_bytes(),
                ledger.read_bytes(),
            )
            self.assertEqual((scheduler / "sacct.psv").read_text(), raw)
            self.assertEqual(result["logical_jobs"], 20)
            self.assertEqual(result["completed_jobs"], 19)
            self.assertEqual(result["cancelled_unstarted_jobs"], 1)
            self.assertEqual(
                result["validation_sha256"],
                importer.sha256(scheduler / "validation.json"),
            )

    def test_scheduler_attestation_rejects_ledger_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            lineage = root / "lineage"
            (source / "contract").mkdir(parents=True)
            lineage.mkdir()
            ledger = source / "contract/slurm_submission.tsv"
            ledger.write_text("original\n")
            raw = "scheduler evidence\n"
            attestation = {
                "raw_sacct_psv": raw,
                "validation": {
                    "status": "PASS",
                    "ledger_sha256": importer.sha256(ledger),
                    "raw_sacct_sha256": importer.sha256_bytes(raw.encode()),
                },
            }
            ledger.write_text("changed\n")
            with self.assertRaisesRegex(SystemExit, "changed after attestation"):
                importer.archive_scheduler_attestation(source, lineage, attestation)

    def test_archived_scheduler_attestation_is_reverified_from_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            lineage = Path(tmp) / "lineage"
            scheduler = lineage / "slurm"
            scheduler.mkdir(parents=True)
            ledger = scheduler / "source_submission_ledger.tsv"
            raw = scheduler / "sacct.psv"
            validation_path = scheduler / "validation.json"
            ledger.write_text("submission ledger\n")
            scheduler_module = importer.load_scheduler_attestation_module()
            normalized_records = []
            raw_buffer = io.StringIO()
            raw_writer = csv.writer(raw_buffer, delimiter="|", lineterminator="\n")
            for index in range(20):
                record = {
                    field: f"{field}-{index}"
                    for field in scheduler_module.SACCT_FIELDS
                }
                record["JobID"] = str(index)
                normalized_records.append({"stage": "fixture", **record})
                raw_writer.writerow(
                    [record[field] for field in scheduler_module.SACCT_FIELDS]
                )
            raw.write_text(raw_buffer.getvalue())
            concurrency = {"bam_manifest_hash": 4, "fragment_recount": 4}
            validation = {
                "status": "PASS",
                "schema_version": importer.EXPECTED_SCHEDULER_ATTESTATION_SCHEMA,
                "source_run_id": "source-run",
                "ledger_sha256": importer.sha256(ledger),
                "raw_sacct_sha256": importer.sha256(raw),
                "logical_jobs": 20,
                "completed_jobs": 19,
                "cancelled_unstarted_jobs": 1,
                "normalized_records": normalized_records,
                "sacct_fields": list(scheduler_module.SACCT_FIELDS),
                "observed_array_max_concurrency": concurrency,
            }
            validation_path.write_text(json.dumps(validation, sort_keys=True) + "\n")
            summary = {
                "source_submission_ledger_sha256": importer.sha256(ledger),
                "raw_sacct_sha256": importer.sha256(raw),
                "validation_sha256": importer.sha256(validation_path),
                "logical_jobs": 20,
                "completed_jobs": 19,
                "cancelled_unstarted_jobs": 1,
                "observed_array_max_concurrency": concurrency,
                "schema_version": importer.EXPECTED_SCHEDULER_ATTESTATION_SCHEMA,
            }
            lineage_json = {
                "execution_source_run_id": "source-run",
                "scheduler_attestation": summary,
                "source_scheduler_attestation_schema": summary["schema_version"],
                "source_scheduler_ledger_sha256": summary[
                    "source_submission_ledger_sha256"
                ],
                "source_scheduler_raw_sacct_sha256": summary["raw_sacct_sha256"],
                "source_scheduler_logical_jobs": 20,
                "source_scheduler_completed_jobs": 19,
                "source_scheduler_cancelled_unstarted_jobs": 1,
                "source_scheduler_observed_array_max_concurrency": concurrency,
            }
            importer.verify_archived_scheduler_attestation(lineage, lineage_json)
            first_row = raw.read_text().splitlines()[0]
            with raw.open("a") as handle:
                handle.write(first_row + "\n")
            validation["raw_sacct_sha256"] = importer.sha256(raw)
            validation_path.write_text(json.dumps(validation, sort_keys=True) + "\n")
            summary["raw_sacct_sha256"] = importer.sha256(raw)
            summary["validation_sha256"] = importer.sha256(validation_path)
            lineage_json["source_scheduler_raw_sacct_sha256"] = importer.sha256(raw)
            with self.assertRaisesRegex(SystemExit, "internally inconsistent"):
                importer.verify_archived_scheduler_attestation(lineage, lineage_json)

    def test_transaction_marker_brackets_every_publish_step(self):
        source = inspect.getsource(importer.perform_import)
        self.assertLess(source.index("if check_only:"), source.index("in_progress ="))
        self.assertLess(source.index("atomic_write_json(\n        in_progress"), source.index("publish_staging("))
        self.assertLess(source.index("publish_staging("), source.index('destination_root / "RECOUNT_COMPLETE"'))
        self.assertLess(source.index('destination_root / "RECOUNT_COMPLETE"'), source.index("in_progress.unlink()"))


if __name__ == "__main__":
    unittest.main()
