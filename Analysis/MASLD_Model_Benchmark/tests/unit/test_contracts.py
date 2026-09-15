from __future__ import annotations

from dataclasses import replace
import gzip
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from masld_bench.contracts import (
    AccessState,
    ArtifactRef,
    ContractError,
    DatasetManifest,
    ExposureState,
    MissingState,
    ModelManifest,
    PairingState,
    ReleaseState,
    RunSpec,
    SelectionLock,
)
from masld_bench.hashing import (
    HashingError,
    canonical_json,
    canonical_sha256,
    sha256_file,
)
from masld_bench.reference import (
    AnnotationReference,
    CoordinateContract,
    GenomeReference,
    KNOWN_TRUNCATED_GENCODE_V49_PATHS,
    ReferenceBundle,
    ReferenceError,
)
from masld_bench.registry import Registry, RegistryError


def dataset_document() -> dict[str, object]:
    return {
        "schema_version": "masld-bench-dataset-v1",
        "dataset_id": "fixture",
        "title": "Fixture dataset",
        "role": "train_development",
        "status": "available",
        "species": "human",
        "accession": ["GSE000000"],
        "access_tier": "public",
        "automatic_download": False,
        "outcome_url": "NOT_REGISTERED",
        "redistribution_class": "source_reference_only",
        "license_status": "source_terms_apply",
        "biological_unit": "donor",
        "expected_biological_units": 4,
        "native_genome_build": "GRCh38",
        "native_annotation_release": "GENCODE_v49_fixture",
        "modalities": ["single_nucleus_rna"],
        "pairing_levels": ["same_study_unpaired"],
        "modality_status_default": "observed",
        "exposure_status": "target_label_unexposed",
        "admission_blocking": False,
        "blockers": [],
        "notes": "Synthetic contract fixture.",
        "provenance": {
            "source_url": "https://example.org/fixture",
            "verified_date": "2026-08-21",
            "local_exposure": "No target labels used for model selection.",
        },
        "split": {
            "group_key": "person_key",
            "outer_role": "train",
            "label_visibility": "development_visible",
        },
        "activation": {
            "schema_version": "masld-bench-dataset-activation-v1",
            "dataset_version": "fixture-v1",
            "rights_sha256": "1" * 64,
            "topology_sha256": "2" * 64,
            "donor_join_sha256": "3" * 64,
            "labels_sha256": "4" * 64,
            "qc_sha256": "5" * 64,
            "reference_sha256": "6" * 64,
            "exposure_audit_sha256": "7" * 64,
            "evidence_artifacts": [
                {
                    "path": f"authority-{authority}.json",
                    "sha256": str(index) * 64,
                    "size_bytes": index,
                    "media_type": "application/json",
                    "role": f"dataset_authority:fixture:{authority}",
                }
                for index, authority in enumerate(
                    (
                        "rights",
                        "topology",
                        "donor_join",
                        "labels",
                        "qc",
                        "reference",
                        "exposure_audit",
                    ),
                    start=1,
                )
            ],
            "artifact_manifest": {
                "path": "fixture-dataset-manifest.json",
                "sha256": "8" * 64,
                "size_bytes": 1,
                "media_type": "application/json",
                "role": "dataset:fixture",
            },
            "ready": True,
        },
    }


class EnumContractTests(unittest.TestCase):
    def test_fail_closed_enum_domains_are_exact(self) -> None:
        self.assertEqual(
            {item.value for item in PairingState},
            {
                "same_molecule",
                "same_cell",
                "same_nucleus",
                "same_section",
                "adjacent_section",
                "same_sample_different_aliquot",
                "same_donor_different_tissue",
                "same_study_unpaired",
            },
        )
        self.assertEqual(len(MissingState), 8)
        self.assertEqual(
            {item.value for item in ExposureState},
            {
                "clean_declared",
                "target_label_unexposed",
                "encoder_seen",
                "continual_seen",
                "reference_only",
                "downstream_demo",
                "unknown",
            },
        )
        self.assertEqual(
            {item.value for item in AccessState}, {"public", "controlled", "conditional"}
        )

    def test_unknown_and_nested_unknown_fields_are_rejected(self) -> None:
        top_level = dataset_document()
        top_level["typo"] = True
        with self.assertRaisesRegex(ContractError, "unknown field"):
            DatasetManifest.from_dict(top_level)

        nested = dataset_document()
        nested["split"] = {**nested["split"], "typo": True}  # type: ignore[arg-type]
        with self.assertRaisesRegex(ContractError, "unknown field"):
            DatasetManifest.from_dict(nested)

    def test_json_loader_round_trip(self) -> None:
        manifest = DatasetManifest.from_json(json.dumps(dataset_document()))
        self.assertEqual(manifest.dataset_id, "fixture")
        self.assertEqual(manifest.pairing_levels, (PairingState.SAME_STUDY_UNPAIRED,))
        self.assertEqual(DatasetManifest.from_json(manifest.to_json()), manifest)

    def test_unknown_exposure_must_block_admission(self) -> None:
        document = dataset_document()
        document["exposure_status"] = "unknown"
        with self.assertRaisesRegex(ContractError, "admission-blocking"):
            DatasetManifest.from_dict(document)

    def test_unresolved_native_reference_identity_must_block_admission(self) -> None:
        for field, value in (
            ("native_genome_build", "mixed_source_UNRESOLVED"),
            ("native_annotation_release", "mixed_source_UNRESOLVED"),
        ):
            document = dataset_document()
            document[field] = value
            with self.subTest(field=field), self.assertRaisesRegex(
                ContractError, "admission-blocking"
            ):
                DatasetManifest.from_dict(document)

    def test_only_deferred_models_may_have_no_supported_tasks(self) -> None:
        family = {
            "schema_version": "masld-bench-model-v1",
            "family_id": "fixture_family",
            "modality": ["single_nucleus_atac"],
            "status": "candidate",
            "models": [],
        }
        model = {
            "model_id": "transductive_fixture",
            "display_name": "Transductive fixture",
            "implementation_type": "train_locally",
            "checkpoint_revision": "fixture/repository@" + "a" * 40,
            "checkpoint_sha256": "NOT_APPLICABLE",
            "license_status": "MIT",
            "exposure_status": "target_label_unexposed",
            "status": "deferred",
            "admission_blocking": True,
            "blockers": ["No inductive held-query transform."],
            "supported_tasks": [],
        }
        manifest = ModelManifest.from_family_entry(family, model)
        self.assertEqual(manifest.supported_tasks, ())

        model["status"] = "candidate"
        with self.assertRaisesRegex(ContractError, "only deferred models"):
            ModelManifest.from_family_entry(family, model)

    def test_embedded_unresolved_model_fields_must_block_admission(self) -> None:
        family = {
            "schema_version": "masld-bench-model-v1",
            "family_id": "fixture_family",
            "modality": ["single_nucleus_rna"],
            "status": "candidate",
            "models": [],
        }
        base = {
            "model_id": "fixture_model",
            "display_name": "Fixture model",
            "implementation_type": "train_locally",
            "checkpoint_revision": "fixture/repository@" + "a" * 40,
            "checkpoint_sha256": "NOT_APPLICABLE",
            "license_status": "MIT",
            "exposure_status": "target_label_unexposed",
            "status": "candidate",
            "admission_blocking": False,
            "blockers": [],
            "supported_tasks": ["cell_state_mapping"],
        }
        for field, value in (
            ("checkpoint_revision", "fixture/repository@UNRESOLVED"),
            ("license_status", "code_MIT_weights_UNRESOLVED"),
        ):
            model = {**base, field: value}
            with self.subTest(field=field), self.assertRaisesRegex(
                ContractError, "UNRESOLVED model fields"
            ):
                ModelManifest.from_family_entry(family, model)

    def test_ready_execution_requires_exact_adaptation_capabilities(self) -> None:
        config_root = Path(__file__).resolve().parents[2] / "config"
        execution = Registry.load(config_root).models[
            "hvg_pca_logistic"
        ].execution
        self.assertEqual(execution.supported_adaptation_regimes, ("native_lane",))
        with self.assertRaisesRegex(ContractError, "requires exact repository"):
            replace(execution, supported_adaptation_regimes=())
        with self.assertRaisesRegex(ContractError, "lowercase identifiers"):
            replace(execution, supported_adaptation_regimes=("Native Lane",))
        for value in ("source@UNRESOLVED", "source@unresolved"):
            with self.subTest(value=value), self.assertRaisesRegex(
                ContractError, "requires exact repository"
            ):
                replace(execution, upstream_commit=value)
        copyleft = replace(
            execution,
            code_license="GPL-2.0-only",
            weights_license="CC-BY-SA-4.0",
            derivative_weights_license="CC-BY-SA-4.0",
        )
        self.assertEqual(copyleft.code_license, "GPL-2.0-only")
        self.assertEqual(copyleft.weights_license, "CC-BY-SA-4.0")

    def test_ready_activation_requires_exact_authority_artifacts(self) -> None:
        missing = dataset_document()
        del missing["activation"]["evidence_artifacts"]  # type: ignore[index]
        with self.assertRaisesRegex(ContractError, "missing field.*evidence_artifacts"):
            DatasetManifest.from_dict(missing)

        mismatched = dataset_document()
        mismatched_evidence = mismatched["activation"]["evidence_artifacts"]  # type: ignore[index]
        mismatched_evidence[0]["sha256"] = "f" * 64  # type: ignore[index]
        with self.assertRaisesRegex(ContractError, "authority artifact SHA-256"):
            DatasetManifest.from_dict(mismatched)

        wrong_namespace = dataset_document()
        namespace_evidence = wrong_namespace["activation"]["evidence_artifacts"]  # type: ignore[index]
        namespace_evidence[0]["role"] = "dataset_authority:other:rights"  # type: ignore[index]
        with self.assertRaisesRegex(ContractError, "bind the DatasetManifest.dataset_id"):
            DatasetManifest.from_dict(wrong_namespace)


class HashingAndRunIdentityTests(unittest.TestCase):
    def test_canonical_hash_ignores_mapping_insertion_order(self) -> None:
        left = {"b": 2, "a": [1, {"d": 4, "c": 3}]}
        right = {"a": [1, {"c": 3, "d": 4}], "b": 2}
        self.assertEqual(canonical_json(left), canonical_json(right))
        self.assertEqual(canonical_sha256(left), canonical_sha256(right))
        with self.assertRaises(HashingError):
            canonical_json({"unordered": {1, 2}})
        with self.assertRaises(HashingError):
            canonical_json(float("nan"))

    def test_run_id_excludes_only_output_location_and_attempt(self) -> None:
        run = RunSpec(
            schema_version="masld-bench-run-v1",
            campaign_id="fixture-campaign",
            task_id="cell_state_mapping",
            model_id="fixture-model",
            dataset_ids=("fixture",),
            split_id="donor_outer",
            seed=20260821,
            action=("predict",),
            adapter_command=("python", "adapter.py"),
            immutable_inputs={"fixture": "2" * 64},
            runtime_id="runtime-lock-v1",
            runtime_registry_sha256="3" * 64,
            resource_profile="cpu_contract",
            code_lock_sha256="1" * 64,
            output_dir="attempt-0",
            attempt=0,
        )
        retried = replace(run, output_dir="attempt-7", attempt=7)
        changed_seed = replace(run, seed=20260822)
        self.assertEqual(run.run_id, retried.run_id)
        self.assertNotEqual(run.run_id, changed_seed.run_id)


class ArtifactAndReferenceTests(unittest.TestCase):
    def test_artifact_size_hash_and_root_containment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            root = parent / "bundle"
            root.mkdir()
            artifact_path = root / "predictions.jsonl"
            artifact_path.write_bytes(b"first")
            artifact = ArtifactRef.from_path(artifact_path, relative_to=root)
            self.assertEqual(artifact.validate(root), artifact_path.resolve())

            artifact_path.write_bytes(b"other")
            with self.assertRaisesRegex(ContractError, "SHA-256 mismatch"):
                artifact.validate(root)

            outside = parent / "outside.txt"
            outside.write_bytes(b"outside")
            escaping = ArtifactRef.from_path(outside)
            escaping = replace(escaping, path="../outside.txt")
            with self.assertRaisesRegex(ContractError, "escapes declared root"):
                escaping.validate(root, require_relative=True)

    def test_artifact_rejects_symlink_file_and_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            root = parent / "bundle"
            root.mkdir()
            target = root / "target.txt"
            target.write_text("target\n", encoding="utf-8")
            link = root / "link.txt"
            link.symlink_to(target)
            with self.assertRaisesRegex(ContractError, "symlink"):
                ArtifactRef.from_path(link)

            root_link = parent / "bundle-link"
            root_link.symlink_to(root, target_is_directory=True)
            with self.assertRaisesRegex(ContractError, "root.*symlink"):
                ArtifactRef.from_path(target, relative_to=root_link)

            artifact = ArtifactRef.from_path(target, relative_to=root)
            with self.assertRaisesRegex(ContractError, "root.*symlink"):
                artifact.validate(root_link, require_relative=True)

            outer = parent / "outer"
            outer.mkdir()
            intermediate = outer / "intermediate"
            intermediate.symlink_to(root, target_is_directory=True)
            with self.assertRaisesRegex(ContractError, "traverse a symlink"):
                ArtifactRef.from_path(intermediate / "target.txt")
            with self.assertRaisesRegex(ContractError, "traverse a symlink"):
                artifact.validate(intermediate, require_relative=True)

    def test_artifact_rejects_change_during_hashing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "authority.json"
            path.write_text("before\n", encoding="utf-8")
            artifact = ArtifactRef.from_path(path)

            def changing_hash(source: str | Path) -> str:
                digest = sha256_file(source)
                Path(source).write_text("changed after digest\n", encoding="utf-8")
                return digest

            with patch("masld_bench.contracts.sha256_file", side_effect=changing_hash):
                with self.assertRaisesRegex(ContractError, "changed while hashing"):
                    artifact.validate()

    def _bundle(self, root: Path) -> ReferenceBundle:
        genome_path = root / "genome.fa.gz"
        annotation_path = root / "annotation.gtf.gz"
        with gzip.open(genome_path, "wb") as handle:
            handle.write(b">chr1\nACGTACGT\n")
        with gzip.open(annotation_path, "wb") as handle:
            handle.write(b"chr1\ttest\tgene\t1\t8\t.\t+\t.\tgene_id \"G1\";\n")
        return ReferenceBundle(
            bundle_id="fixture-reference",
            assembly="fixture",
            assembly_patch="p0",
            annotation_release="fixture-v1",
            genome=GenomeReference(
                path=genome_path.name,
                sha256=sha256_file(genome_path),
                compression="gzip",
                integrity="gzip_and_sha256_verified",
            ),
            annotation=AnnotationReference(
                path=annotation_path.name,
                sha256=sha256_file(annotation_path),
                format="gtf_gzip",
                integrity="gzip_and_sha256_verified",
            ),
            sequence_extraction_ready=True,
            sequence_extraction_blockers=(),
            coordinate_contract=CoordinateContract(
                bed_system="zero_based_half_open",
                contig_policy="fixture",
                analysis_contigs=("chr1",),
                excluded_contig_policy="reject",
                variant_normalization="fixture",
                liftover_policy="fixture",
                native_ld_policy="fixture",
            ),
            prohibited_paths=tuple(sorted(KNOWN_TRUNCATED_GENCODE_V49_PATHS)),
            prohibited_reason="Verified truncated fixtures are prohibited.",
        )

    def test_reference_streams_gzip_and_checks_sha(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self._bundle(root)
            validated = bundle.validate(root)
            self.assertIsInstance(validated["genome"], str)

            genome_path = root / "genome.fa.gz"
            payload = genome_path.read_bytes()
            genome_path.write_bytes(payload[:-4])
            truncated = replace(
                bundle,
                genome=replace(bundle.genome, sha256=sha256_file(genome_path)),
            )
            with self.assertRaisesRegex(ReferenceError, "gzip integrity check failed"):
                truncated.validate(root)

    def test_known_truncated_v49_paths_are_denylisted(self) -> None:
        expected = {
            "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/GRCh38.p14.genome.fa.gz",
            "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.transcripts.fa.gz",
        }
        self.assertEqual(KNOWN_TRUNCATED_GENCODE_V49_PATHS, expected)

    def test_reference_rejects_annotation_contig_absent_from_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self._bundle(root)
            annotation_path = root / "annotation.gtf.gz"
            with gzip.open(annotation_path, "wb") as handle:
                handle.write(
                    b"chr2\ttest\tgene\t1\t8\t.\t+\t.\tgene_id \"G2\";\n"
                )
            mismatch = replace(
                bundle,
                annotation=replace(
                    bundle.annotation,
                    sha256=sha256_file(annotation_path),
                ),
            )
            with self.assertRaisesRegex(ReferenceError, "absent from the genome"):
                mismatch.validate(root)

    def test_reference_rejects_analysis_contig_absent_from_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self._bundle(root)
            mismatch = replace(
                bundle,
                coordinate_contract=replace(
                    bundle.coordinate_contract,
                    analysis_contigs=("chr1", "chr2"),
                ),
            )
            with self.assertRaisesRegex(ReferenceError, "analysis coordinate contract"):
                mismatch.validate(root)

    def test_coordinate_contract_rejects_alt_patch_and_scaffold_contigs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            contract = self._bundle(Path(temporary)).coordinate_contract
            self.assertEqual(contract.require_analysis_contig("chr1"), "chr1")
            with self.assertRaisesRegex(ReferenceError, "excluded"):
                contract.require_analysis_contig("chr1_KI270706v1_random")

    def test_sequence_jobs_require_an_admitted_indexed_copy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            unindexed = replace(
                self._bundle(root),
                sequence_extraction_ready=False,
                sequence_extraction_blockers=("derived indexes are absent",),
            )
            unindexed.validate(root)
            with self.assertRaisesRegex(ReferenceError, "indexed sequence extraction"):
                unindexed.validate(root, require_indexed_sequence=True)


class SelectionAndRegistryTests(unittest.TestCase):
    def test_selection_lock_is_fail_closed(self) -> None:
        candidate = "f" * 64
        lock = SelectionLock(
            schema_version="masld-bench-selection-lock-v1",
            lock_id="a" * 64,
            campaign_id="fixture-campaign",
            plan_sha256="b" * 64,
            candidate_manifest_sha256="c" * 64,
            registry_sha256="d" * 64,
            metrics_sha256="e" * 64,
            locked=True,
            outcomes_unlocked=False,
            release_state=ReleaseState.CANDIDATE,
            candidate_run_ids=(candidate,),
            selected_run_ids=(candidate,),
            task_decisions=({"task_id": "cell_state_mapping"},),
            conditional_model_decision={"decision": "not_triggered"},
            power_decisions={},
            multiplicity_plan="holm_confirmatory_bh_secondary",
            terminal_policy="no_reselection_recalibration_threshold_change_or_repair",
        )
        self.assertIs(lock.require_locked(for_external_scoring=True), lock)
        with self.assertRaisesRegex(ContractError, "not locked"):
            replace(lock, locked=False).require_locked()
        with self.assertRaisesRegex(ContractError, "after outcomes were unlocked"):
            replace(lock, outcomes_unlocked=True).require_locked(
                for_external_scoring=True
            )

    def test_checked_in_registry_loads_deterministically(self) -> None:
        config_root = Path(__file__).resolve().parents[2] / "config"
        first = Registry.load(config_root)
        second = Registry.load(config_root)
        self.assertGreaterEqual(len(first.datasets), 1)
        self.assertGreaterEqual(len(first.models), 1)
        self.assertGreaterEqual(len(first.tasks), 1)
        self.assertIn(
            "resource_atlas_geneformer_smoke_1000_v1", first.dataset_views
        )
        view = first.get_dataset_view(
            "resource_atlas_geneformer_smoke_1000_v1"
        )
        self.assertEqual(view.parent_dataset_id, "resource_atlas_current")
        self.assertEqual(view.allowed_waves, ("smoke",))
        self.assertEqual(view.row_count, 1000)
        self.assertFalse(view.selection_outcomes_used)
        self.assertFalse(view.sealed_outcomes_used)
        self.assertIsNotNone(first.variant_capabilities)
        self.assertIsNotNone(first.variant_proxy_selection)
        self.assertFalse(
            first.variant_proxy_selection["eqtl_specific_head_fitting_allowed"]
        )
        self.assertFalse(
            first.variant_proxy_selection["weighted_composite_allowed"]
        )
        self.assertEqual(
            len(first.variant_capabilities["models"]),
            len(
                [
                    model
                    for model in first.models.values()
                    if "variant_to_regulation" in model.supported_tasks
                ]
            ),
        )
        self.assertEqual(first.snapshot_sha256, second.snapshot_sha256)

    def test_registry_rejects_variant_proxy_policy_relabeling(self) -> None:
        config_root = Path(__file__).resolve().parents[2] / "config"
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "config"
            shutil.copytree(config_root, copied)
            policy = (
                copied
                / "evaluation"
                / "variant_development_proxy_selection.json"
            )
            document = json.loads(policy.read_text(encoding="utf-8"))
            document["eqtl_specific_head_fitting_allowed"] = True
            policy.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(RegistryError, "proxy-selection policy"):
                Registry.load(copied)

    def test_registry_rejects_overlapping_rna_atac_capabilities(self) -> None:
        config_root = Path(__file__).resolve().parents[2] / "config"
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "config"
            shutil.copytree(config_root, copied)
            capability = copied / "evaluation" / "rna_conditioned_atac_capabilities.toml"
            text = capability.read_text(encoding="utf-8")
            text = text.replace(
                '"mofaplus", "pca_procrustes"',
                '"mofaplus", "multivi", "pca_procrustes"',
            )
            capability.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(RegistryError, "categories must be disjoint"):
                Registry.load(copied)

    def test_registry_rejects_atac_only_model_on_rna_cell_task(self) -> None:
        config_root = Path(__file__).resolve().parents[2] / "config"
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "config"
            shutil.copytree(config_root, copied)
            models = copied / "models" / "chromatin.toml"
            text = models.read_text(encoding="utf-8")
            prefix, epiagent = text.split('model_id = "epiagent"', maxsplit=1)
            epiagent = epiagent.replace(
                "supported_tasks = []",
                'supported_tasks = ["cell_state_mapping"]',
                1,
            )
            models.write_text(
                prefix + 'model_id = "epiagent"' + epiagent,
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                RegistryError,
                "ATAC-only model epiagent cannot support the RNA-only cell-state task",
            ):
                Registry.load(copied)

    def test_registry_rejects_observed_context_variant_as_primary(self) -> None:
        config_root = Path(__file__).resolve().parents[2] / "config"
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "config"
            shutil.copytree(config_root, copied)
            capability = (
                copied
                / "evaluation"
                / "variant_to_regulation_capabilities.toml"
            )
            text = capability.read_text(encoding="utf-8")
            marker = 'model_id = "epibert"'
            prefix, record = text.split(marker, maxsplit=1)
            record = record.replace(
                'allowed_endpoints = ["accessibility_delta", "allelic_direction", "eqtl_retrieval"]',
                'allowed_endpoints = ["accessibility_delta", "allelic_direction", "eqtl_retrieval", "signed_cell_type_eqtl_effect"]',
                1,
            ).replace("primary_eligible = false", "primary_eligible = true", 1)
            capability.write_text(prefix + marker + record, encoding="utf-8")
            with self.assertRaisesRegex(
                RegistryError,
                "invalid primary capability role",
            ):
                Registry.load(copied)

if __name__ == "__main__":
    unittest.main()
