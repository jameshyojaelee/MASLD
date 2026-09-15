from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest.mock import patch

from masld_bench.artifacts import canonical_hash, verify_frozen_tree
from masld_bench.campaign import (
    CampaignError,
    _assert_config_snapshot,
    _observed_source_lock,
    execute_admission_run,
    execute_run,
    submission_commands,
    verify_run_execution_attempt,
)
from masld_bench.hashing import sha256_file
from masld_bench.firewall import resource_snapshot
from masld_bench.planner import (
    PlanningError,
    _contains_unresolved,
    _model_has_runnable_task_contract,
    _resolved_model_hyperparameters,
    build_plan,
    freeze_campaign,
    load_frozen_plan,
)
from masld_bench.registry import Registry


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REPOSITORY_ROOT = PACKAGE_ROOT.parents[1]
_FIXTURE_RESOURCE_FIREWALL = resource_snapshot(
    REPOSITORY_ROOT,
    ("Analysis/MASLD_Model_Benchmark/pyproject.toml",),
)


def _test_source_parent() -> Path:
    path = PACKAGE_ROOT / "executions" / "tests"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content).strip() + "\n", encoding="utf-8")


def _fixture_config(root: Path) -> Path:
    config = root / "config"
    dataset_artifact = config / "artifacts" / "fixture_development.manifest.json"
    _write(
        dataset_artifact,
        '{"dataset_id":"fixture_development","donor_ids":["donor-1","donor-2"]}',
    )
    dataset_artifact_sha256 = sha256_file(dataset_artifact)
    dataset_artifact_size = dataset_artifact.stat().st_size
    environment_artifact = config / "artifacts" / "fixture-environment.lock"
    _write(
        environment_artifact,
        f"python={Path(sys.executable).resolve().as_posix()}\nversion={sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
    )
    environment_artifact_sha256 = sha256_file(environment_artifact)
    environment_artifact_size = environment_artifact.stat().st_size
    dataset_authorities: dict[str, tuple[Path, str, int]] = {}
    for authority in (
        "rights",
        "topology",
        "donor_join",
        "labels",
        "qc",
        "reference",
        "exposure_audit",
    ):
        path = config / "artifacts" / "dataset-authorities" / f"{authority}.json"
        _write(path, f'{{"authority":"{authority}","dataset_id":"fixture_development"}}')
        dataset_authorities[authority] = (
            path,
            sha256_file(path),
            path.stat().st_size,
        )
    dataset_evidence_toml = ", ".join(
        (
            f'{{ path = "artifacts/dataset-authorities/{authority}.json", '
            f'sha256 = "{digest}", size_bytes = {size}, '
            f'media_type = "application/json", '
            f'role = "dataset_authority:fixture_development:{authority}" }}'
        )
        for authority, (_, digest, size) in dataset_authorities.items()
    )
    model_authorities: dict[str, tuple[Path, str, int]] = {}
    license_values = {
        "code_license": "project_owned",
        "weights_license": "not_applicable_classical",
        "derivative_weights_license": "project_owned",
    }
    for authority in (
        "preprocessing",
        "code_license",
        "weights_license",
        "derivative_weights_license",
        "declared_corpora",
        "cellxgene_uuids",
        "exposure_audit",
    ):
        path = config / "artifacts" / "model-authorities" / f"{authority}.txt"
        if authority in license_values:
            _write(
                path,
                (
                    '{"authority":"%s","declared_license":"%s",'
                    '"derivative_redistribution_allowed":true,'
                    '"model_id":"fixture_open","redistribution_allowed":true,'
                    '"schema_version":"masld-bench-license-authority-v1",'
                    '"use_allowed":true}'
                )
                % (authority, license_values[authority]),
            )
        else:
            _write(path, f"fixture_open {authority}")
        model_authorities[authority] = (
            path,
            sha256_file(path),
            path.stat().st_size,
        )
    model_evidence_toml = ", ".join(
        (
            f'{{ path = "artifacts/model-authorities/{authority}.txt", '
            f'sha256 = "{digest}", size_bytes = {size}, media_type = "text/plain", '
            f'role = "model_authority:fixture_open:{authority}" }}'
        )
        for authority, (_, digest, size) in model_authorities.items()
    )
    _write(
        config / "datasets" / "fixture_development.toml",
        f"""
        schema_version = "masld-bench-dataset-v1"
        dataset_id = "fixture_development"
        title = "Fixture"
        role = "train_development"
        status = "available"
        species = "human"
        accession = ["FIXTURE-DEVELOPMENT"]
        access_tier = "public"
        automatic_download = false
        redistribution_class = "fixture_only"
        license_status = "project_owned"
        biological_unit = "donor"
        expected_biological_units = 2
        native_genome_build = "GRCh38"
        native_annotation_release = "GENCODE_v49_fixture"
        modalities = ["single_cell_rna"]
        pairing_levels = ["same_donor_different_tissue"]
        modality_status_default = "observed"
        exposure_status = "target_label_unexposed"
        admission_blocking = false
        blockers = []
        notes = "Synthetic control-plane fixture only."
        [provenance]
        source_url = "https://example.invalid/fixture"
        verified_date = "2026-08-21"
        local_exposure = "synthetic"
        [split]
        group_key = "cohort_family_id+donor_id"
        outer_role = "development"
        label_visibility = "development"
        [activation]
        schema_version = "masld-bench-dataset-activation-v1"
        dataset_version = "fixture-v1"
        rights_sha256 = "{dataset_authorities['rights'][1]}"
        topology_sha256 = "{dataset_authorities['topology'][1]}"
        donor_join_sha256 = "{dataset_authorities['donor_join'][1]}"
        labels_sha256 = "{dataset_authorities['labels'][1]}"
        qc_sha256 = "{dataset_authorities['qc'][1]}"
        reference_sha256 = "{dataset_authorities['reference'][1]}"
        exposure_audit_sha256 = "{dataset_authorities['exposure_audit'][1]}"
        evidence_artifacts = [{dataset_evidence_toml}]
        artifact_manifest = {{ path = "artifacts/fixture_development.manifest.json", sha256 = "{dataset_artifact_sha256}", size_bytes = {dataset_artifact_size}, media_type = "application/json", role = "dataset:fixture_development" }}
        ready = true
        """,
    )
    _write(
        config / "models" / "fixture_family.toml",
        f"""
        schema_version = "masld-bench-model-family-v1"
        family_id = "fixture_family"
        modality = ["transcriptome"]
        status = "candidate"
        [[models]]
        model_id = "fixture_open"
        display_name = "Fixture"
        implementation_type = "classical_baseline"
        checkpoint_revision = "fixture-v1"
        checkpoint_sha256 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        license_status = "open"
        exposure_status = "clean_declared"
        status = "candidate"
        admission_blocking = false
        blockers = []
        supported_tasks = ["fixture_task"]
        [models.execution]
        schema_version = "masld-bench-model-execution-v1"
        upstream_repository = "fixture://classical-baseline"
        upstream_commit = "fixture-v1"
        preprocessing_sha256 = "{model_authorities['preprocessing'][1]}"
        vocabulary_or_build = "not_applicable_classical"
        genome_build = "GRCh38"
        code_license = "project_owned"
        weights_license = "not_applicable_classical"
        derivative_weights_license = "project_owned"
        training_cutoff = "not_applicable_classical"
        declared_corpora_sha256 = "{model_authorities['declared_corpora'][1]}"
        cellxgene_uuids_sha256 = "{model_authorities['cellxgene_uuids'][1]}"
        exposure_audit_sha256 = "{model_authorities['exposure_audit'][1]}"
        evidence_artifacts = [{model_evidence_toml}]
        supported_actions = ["probe", "prepare", "fit", "predict", "export"]
        supported_adaptation_regimes = ["common_lane", "native_lane"]
        runtime_id = "fixture_cpu"
        adapter_command = ["-m", "masld_bench.adapters.stub"]
        environment_artifact = {{ path = "artifacts/fixture-environment.lock", sha256 = "{environment_artifact_sha256}", size_bytes = {environment_artifact_size}, media_type = "text/plain", role = "environment:fixture_cpu" }}
        capture_r_session = false
        ready = true
        [[models]]
        model_id = "blocked_unknown"
        display_name = "Blocked"
        implementation_type = "classical_baseline"
        checkpoint_revision = "UNRESOLVED"
        checkpoint_sha256 = "UNRESOLVED"
        license_status = "UNRESOLVED"
        exposure_status = "unknown"
        status = "blocked_terms"
        admission_blocking = true
        blockers = ["terms_unresolved"]
        supported_tasks = ["fixture_task"]
        """,
    )
    promotion = config / "evaluation" / "promotion_gates.toml"
    _write(
        promotion,
        """
        schema_version = "masld-bench-promotion-gates-v1"
        confirmatory_method = "Holm"
        confirmatory_fwer = 0.05
        [tasks.fixture_task]
        promotion_gate_id = "fixture_promotion_v1"
        allowed_promotion_modes = ["development_only"]
        sealed_holdout_dataset_ids = []
        """,
    )
    promotion_sha256 = sha256_file(promotion)
    _write(
        config / "tasks" / "fixture_task.toml",
        f"""
        schema_version = "masld-bench-task-v1"
        task_id = "fixture_task"
        status = "candidate"
        unit_of_inference = "donor"
        input_modalities = ["rna"]
        endpoint = "fixture"
        metrics = ["macro_f1"]
        datasets_train = ["fixture_development"]
        datasets_development = []
        datasets_sealed = []
        split_id = "donor_grouped"
        required_pairing_levels = ["same_donor_different_tissue"]
        missingness_policy = "explicit_state_and_mask"
        admission_gates = ["fixture"]
        baseline_model_ids = ["fixture_open"]
        uncertainty_method = "paired_donor_cluster_bootstrap"
        resampling_units = ["donor"]
        bootstrap_replicates = 10000
        multiplicity_family = "fixture_confirmatory"
        multiplicity_method = "Holm"
        promotion_gate_id = "fixture_promotion_v1"
        promotion_gate_config_sha256 = "{promotion_sha256}"
        primary_evaluator_id = "not_applicable_nonpromotable_v1"
        evaluator_parameters = {{}}
        claim_gate = "none"
        """,
    )
    _write(
        config / "splits" / "donor_grouped.toml",
        """
        schema_version = "masld-bench-split-v1"
        split_id = "donor_grouped"
        entity = "donor"
        group_keys = ["cohort_family_id", "donor_id"]
        roles = ["train", "development", "withheld_sealed"]
        seed = 20260821
        locus_grouping = "not_applicable"
        label_policy = "All donor-linked observations remain in one outer role."
        exposure_policy = "Checkpoint exposure is audited separately."
        admission_gates = ["Fixture donor grouping is frozen."]
        """,
    )
    _write(
        config / "runtimes" / "fixture.toml",
        f"""
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
        """,
    )
    _write(config / "resources.toml", (PACKAGE_ROOT / "config" / "resources.toml").read_text(encoding="utf-8"))
    campaign = config / "campaigns" / "fixture.toml"
    _write(
        campaign,
        """
        schema_version = "masld-bench-campaign-v1"
        campaign_id = "fixture-admission"
        wave = "admission"
        allow_arrays = false
        allow_sealed_features = false
        allow_sealed_labels = false
        submit_enabled = false
        task_ids = ["fixture_task"]
        [selection]
        model_statuses = ["candidate", "blocked_terms"]
        [execution]
        default_resource_profile = "cpu_contract"
        runtime_id = "fixture_cpu"
        seeds = [20260821]
        """,
    )
    return campaign


def _fixture_scientific_campaign(root: Path, prerequisite_attempt: Path) -> Path:
    config = root / "config"
    prerequisite = prerequisite_attempt / "run_execution_receipt.json"
    campaign = config / "campaigns" / "fixture_smoke.toml"
    _write(
        campaign,
        f"""
        schema_version = "masld-bench-campaign-v1"
        campaign_id = "fixture-smoke"
        wave = "smoke"
        status = "ready"
        allow_arrays = false
        allow_sealed_features = false
        allow_sealed_labels = false
        submit_enabled = false
        task_ids = ["fixture_task"]
        prerequisites = ["fixture_admission"]
        [prerequisite_receipts.fixture_admission]
        path = "{prerequisite.resolve().as_posix()}"
        sha256 = "{sha256_file(prerequisite)}"
        size_bytes = {prerequisite.stat().st_size}
        media_type = "application/json"
        role = "prerequisite:fixture_admission"
        [selection]
        model_statuses = ["candidate"]
        include_mandatory_baselines = true
        [execution]
        default_resource_profile = "cpu_contract"
        runtime_id = "fixture_cpu"
        seeds = [20260821]
        folds = [0]
        adaptation_regimes = ["common_lane"]
        [execution.adapter_request]
        row_ids = ["donor-1", "donor-2"]
        """,
    )
    return campaign


def _add_fixture_dataset_view(root: Path) -> str:
    """Add a complete two-row smoke view over an unresolved fixture parent."""

    config = root / "config"
    parent = config / "datasets" / "fixture_development.toml"
    document = parent.read_text(encoding="utf-8")
    document = document[: document.index("[activation]")]
    document = document.replace(
        "admission_blocking = false", "admission_blocking = true", 1
    ).replace(
        "blockers = []", 'blockers = ["full_parent_activation_unresolved"]', 1
    ).replace(
        'label_visibility = "development"',
        'label_visibility = "development_visible"',
        1,
    )
    parent.write_text(document, encoding="utf-8")

    view_id = "fixture_smoke_view_v1"
    artifact_root = config / "artifacts" / "dataset-views" / view_id
    data = artifact_root / "data.json"
    selection = artifact_root / "selection.tsv"
    ontology = artifact_root / "ontology.json"
    derivation = artifact_root / "derivation.json"
    manifest = artifact_root / "ARTIFACTS.json"
    complete = artifact_root / "COMPLETE"
    _write(data, '{"rows":2}')
    _write(
        selection,
        """
        cell_id\tdonor_id\tbroad_label
        cell-1\tdonor-1\tclass_a
        cell-2\tdonor-2\tclass_b
        """,
    )
    _write(
        ontology,
        json.dumps(
            {
                "schema_version": "masld-bench-cell-ontology-v1",
                "ontology_id": "fixture_two_class_v1",
                "classes": {"class_a": ["class_a"], "class_b": ["class_b"]},
                "selection_outcomes_used": False,
                "sealed_outcomes_used": False,
            },
            sort_keys=True,
        ),
    )
    _write(
        derivation,
        json.dumps(
            {
                "schema_version": "masld-bench-geneformer-smoke-subset-v1",
                "subset_id": view_id,
                "dataset_id": "fixture_development",
                "manifest_sha256": "a" * 64,
                "geneformer_input_contract": {"ready": True},
                "selection": {
                    "cell_budget": 2,
                    "selection_policy": "fixture_donor_round_robin",
                    "selection_seed": 20260821,
                    "selection_outcomes_used": False,
                    "sealed_outcomes_used": False,
                    "counts_by_class": {"class_a": 1, "class_b": 1},
                },
                "matrix": {"shape": [2, 1]},
                "artifacts": {
                    "h5ad": {"sha256": sha256_file(data)},
                    "selection": {"sha256": sha256_file(selection)},
                    "ontology": {"sha256": sha256_file(ontology)},
                },
            },
            sort_keys=True,
        ),
    )
    records = [
        {
            "path": path.name,
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in (data, selection, ontology, derivation)
    ]
    _write(
        manifest,
        json.dumps(
            {
                "schema_version": "masld-bench-artifacts-v1",
                "metadata": {
                    "artifact_class": "geneformer_smoke_subset",
                    "subset_id": view_id,
                },
                "artifacts": sorted(records, key=lambda item: item["path"]),
            },
            sort_keys=True,
        ),
    )
    _write(
        complete,
        json.dumps(
            {
                "schema_version": "masld-bench-complete-v1",
                "manifest_sha256": sha256_file(manifest),
                "artifact_count": 4,
            },
            sort_keys=True,
        ),
    )

    authorities: dict[str, Path] = {}
    authority_root = config / "artifacts" / "dataset-view-authorities" / view_id
    for authority in (
        "rights",
        "topology",
        "donor_join",
        "labels",
        "qc",
        "reference",
        "exposure_audit",
    ):
        path = authority_root / f"{authority}.json"
        _write(
            path,
            json.dumps(
                {
                    "schema_version": "masld-bench-dataset-view-authority-v1",
                    "view_id": view_id,
                    "parent_dataset_id": "fixture_development",
                    "authority": authority,
                    "admitted": True,
                    "evidence": {"fixture": True},
                },
                sort_keys=True,
            ),
        )
        authorities[authority] = path

    def artifact_ref(path: Path, role: str, media_type: str) -> str:
        relative = path.relative_to(config).as_posix()
        return (
            f'{{ path = "{relative}", sha256 = "{sha256_file(path)}", '
            f'size_bytes = {path.stat().st_size}, media_type = "{media_type}", '
            f'role = "{role}" }}'
        )

    authority_toml = ",\n          ".join(
        artifact_ref(
            path,
            f"dataset_view_authority:{view_id}:{authority}",
            "application/json",
        )
        for authority, path in authorities.items()
    )
    _write(
        config / "dataset_views" / f"{view_id}.toml",
        f"""
        schema_version = "masld-bench-dataset-view-v1"
        view_id = "{view_id}"
        parent_dataset_id = "fixture_development"
        parent_registry_sha256 = "{sha256_file(parent)}"
        purpose = "compatibility_smoke"
        allowed_waves = ["smoke"]
        row_count = 2
        biological_unit_count = 2
        modalities = ["single_cell_rna"]
        pairing_levels = ["same_donor_different_tissue"]
        label_visibility = "development_visible"
        selection_policy = "fixture_donor_round_robin"
        selection_seed = 20260821
        selection_outcomes_used = false
        sealed_outcomes_used = false
        class_counts = {{ class_a = 1, class_b = 1 }}
        artifact_manifest = {artifact_ref(manifest, f"dataset_view_manifest:{view_id}", "application/json")}
        complete_receipt = {artifact_ref(complete, f"dataset_view_complete:{view_id}", "application/json")}
        data_artifact = {artifact_ref(data, f"dataset_view_data:{view_id}", "application/json")}
        selection_artifact = {artifact_ref(selection, f"dataset_view_selection:{view_id}", "text/tab-separated-values")}
        ontology_artifact = {artifact_ref(ontology, f"dataset_view_ontology:{view_id}", "application/json")}
        derivation_artifact = {artifact_ref(derivation, f"dataset_view_derivation:{view_id}", "application/json")}
        authority_artifacts = [
          {authority_toml}
        ]
        ready = true
        """,
    )
    return view_id


def _rewrite_fixture_license_authority(
    root: Path,
    *,
    authority: str,
    old_license: str,
    new_license: str,
    use_allowed: bool,
    redistribution_allowed: bool,
    derivative_redistribution_allowed: bool,
) -> None:
    path = root / "config" / "artifacts" / "model-authorities" / f"{authority}.txt"
    old_sha256 = sha256_file(path)
    old_size = path.stat().st_size
    _write(
        path,
        json.dumps(
            {
                "schema_version": "masld-bench-license-authority-v1",
                "model_id": "fixture_open",
                "authority": authority,
                "declared_license": new_license,
                "use_allowed": use_allowed,
                "redistribution_allowed": redistribution_allowed,
                "derivative_redistribution_allowed": (
                    derivative_redistribution_allowed
                ),
            },
            sort_keys=True,
            separators=(",", ":"),
        ),
    )
    new_sha256 = sha256_file(path)
    new_size = path.stat().st_size
    model = root / "config" / "models" / "fixture_family.toml"
    document = model.read_text(encoding="utf-8")
    document = document.replace(
        f'{authority} = "{old_license}"',
        f'{authority} = "{new_license}"',
        1,
    )
    old_fragment = (
        f'path = "artifacts/model-authorities/{authority}.txt", '
        f'sha256 = "{old_sha256}", size_bytes = {old_size}'
    )
    new_fragment = (
        f'path = "artifacts/model-authorities/{authority}.txt", '
        f'sha256 = "{new_sha256}", size_bytes = {new_size}'
    )
    if old_fragment not in document:
        raise AssertionError(f"fixture authority fragment is missing: {authority}")
    model.write_text(
        document.replace(old_fragment, new_fragment, 1), encoding="utf-8"
    )


def _execute_fixture_admission(root: Path) -> Path:
    campaign = _fixture_config(root)
    candidate = freeze_campaign(
        campaign,
        config_root=root / "config",
        package_root=PACKAGE_ROOT,
        output_root=root / "candidates",
    )
    frozen = load_frozen_plan(candidate)
    return execute_admission_run(
        candidate=candidate,
        run_id=frozen["runs"][0]["run_id"],
        output_root=root / "executions" / "admission",
    )


class CampaignIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        # Integration cases exercise plan/run firewall bindings against one real,
        # immutable repository authority. The production candidate freeze after
        # the suite still inventories the complete Resource authority surface.
        self._resource_firewall_patch = patch(
            "masld_bench.planner.resource_snapshot",
            return_value=_FIXTURE_RESOURCE_FIREWALL,
        )
        self._resource_firewall_patch.start()

    def tearDown(self) -> None:
        self._resource_firewall_patch.stop()

    def test_rna_atac_reference_and_baseline_gates_use_runnable_models(self) -> None:
        registry = Registry.load(PACKAGE_ROOT / "config", validate_references=False)

        def runnable(model_id: str) -> bool:
            return _model_has_runnable_task_contract(
                registry.models[model_id],
                task_id="rna_conditioned_atac",
                actions=("prepare", "fit", "predict"),
                adaptation_regimes=("native_lane",),
            )

        self.assertTrue(runnable("scpair"))
        self.assertTrue(runnable("shrunken_pseudobulk"))
        self.assertFalse(runnable("chrombpnet"))
        self.assertFalse(runnable("scbasset"))

    def test_model_hyperparameters_replace_instead_of_merge(self) -> None:
        resolved = _resolved_model_hyperparameters(
            "baseline",
            default={"scpair_only": 1, "shared": "wrong"},
            by_model={"baseline": {"baseline_only": 2, "shared": "right"}},
        )
        self.assertEqual(resolved, {"baseline_only": 2, "shared": "right"})
        self.assertNotIn("scpair_only", resolved)

    def test_unresolved_control_markers_do_not_hide_in_compound_values(self) -> None:
        self.assertTrue(
            _contains_unresolved({"environment_lock": "fixture_unresolved_lock"})
        )
        self.assertFalse(
            _contains_unresolved(
                {"description": "This prose describes an unresolved audit."}
            )
        )

    def test_freeze_admit_and_submission_guards(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            plan = build_plan(campaign, config_root=root / "config", package_root=PACKAGE_ROOT)
            self.assertEqual(len(plan["runs"]), 1)
            self.assertEqual(
                plan["model_census_sha256"],
                canonical_hash(sorted(item["model_id"] for item in plan["model_dispositions"])),
            )
            task_disposition = plan["task_dispositions"][0]
            self.assertEqual(task_disposition["split_id"], "donor_grouped")
            self.assertEqual(
                task_disposition["split_contract_sha256"],
                canonical_hash(task_disposition["split_contract"]),
            )
            self.assertEqual(
                plan["runs"][0]["metadata"]["split_contract_sha256"],
                task_disposition["split_contract_sha256"],
            )
            self.assertEqual(
                plan["runs"][0]["metadata"]["fit_dataset_ids"],
                ["fixture_development"],
            )
            self.assertEqual(
                plan["runs"][0]["metadata"]["development_prediction_dataset_ids"],
                [],
            )
            control_plane = plan["control_plane_lock"]
            self.assertEqual(
                control_plane["python_executable_sha256"],
                sha256_file(Path(control_plane["python_executable"])),
            )
            self.assertEqual(control_plane["python_version"], sys.version)
            blocked = {item["model_id"]: item for item in plan["model_dispositions"]}
            self.assertEqual(blocked["blocked_unknown"]["disposition"], "blocked")
            self.assertTrue(blocked["fixture_open"]["open_champion_eligible"])
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "candidates",
            )
            verify_frozen_tree(candidate)
            frozen = load_frozen_plan(candidate)
            commands = submission_commands(
                candidate=candidate,
                approved_campaign_sha256=sha256_file(candidate / "ARTIFACTS.json"),
                execute=False,
            )
            self.assertEqual(len(commands), 1)
            with self.assertRaises(CampaignError):
                submission_commands(
                    candidate=candidate,
                    approved_campaign_sha256="0" * 64,
                    execute=False,
                )
            real_execution = root / "real-executions"
            real_execution.mkdir()
            execution_link = root / "execution-link"
            execution_link.symlink_to(real_execution, target_is_directory=True)
            with self.assertRaisesRegex(CampaignError, "symlink"):
                execute_admission_run(
                    candidate=candidate,
                    run_id=frozen["runs"][0]["run_id"],
                    output_root=execution_link,
                )
            real_nested_output = root / "real-nested-output"
            real_nested_output.mkdir()
            nested_output_parent_link = root / "nested-output-parent-link"
            nested_output_parent_link.symlink_to(
                real_nested_output, target_is_directory=True
            )
            with self.assertRaisesRegex(CampaignError, "symlink"):
                execute_admission_run(
                    candidate=candidate,
                    run_id=frozen["runs"][0]["run_id"],
                    output_root=nested_output_parent_link / "executions",
                )
            attempt = execute_admission_run(
                candidate=candidate,
                run_id=frozen["runs"][0]["run_id"],
                output_root=root / "executions",
            )
            verify_frozen_tree(attempt)
            receipt = verify_run_execution_attempt(attempt)
            self.assertEqual(receipt["status"], "succeeded")
            self.assertEqual(
                [item["action"] for item in receipt["action_receipts"]],
                ["probe"],
            )
            candidate_link = root / "candidate-link"
            candidate_link.symlink_to(candidate, target_is_directory=True)
            with self.assertRaisesRegex(CampaignError, "symlink"):
                submission_commands(
                    candidate=candidate_link,
                    approved_campaign_sha256=sha256_file(candidate / "ARTIFACTS.json"),
                    execute=False,
                )
            candidate_parent_link = root / "candidate-parent-link"
            candidate_parent_link.symlink_to(
                candidate.parent, target_is_directory=True
            )
            with self.assertRaisesRegex(CampaignError, "symlink"):
                submission_commands(
                    candidate=candidate_parent_link / candidate.name,
                    approved_campaign_sha256=sha256_file(candidate / "ARTIFACTS.json"),
                    execute=False,
                )
            attempt_link = root / "attempt-link"
            attempt_link.symlink_to(attempt, target_is_directory=True)
            with self.assertRaisesRegex(CampaignError, "symlink"):
                verify_run_execution_attempt(attempt_link)
            attempt_parent_link = root / "attempt-parent-link"
            attempt_parent_link.symlink_to(attempt.parent, target_is_directory=True)
            with self.assertRaisesRegex(CampaignError, "symlink"):
                verify_run_execution_attempt(attempt_parent_link / attempt.name)

    def test_review_only_admission_has_complete_census_and_no_jobs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            campaign.write_text(
                campaign.read_text(encoding="utf-8").replace(
                    'wave = "admission"',
                    'wave = "admission"\nreview_only = true',
                    1,
                ),
                encoding="utf-8",
            )
            plan = build_plan(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
            )
            self.assertEqual(plan["runs"], [])
            self.assertEqual(plan["dag"]["run_nodes"], [])
            self.assertEqual(plan["dag"]["edges"], [])
            self.assertEqual(plan["resource_totals"]["jobs"], 0)
            review = plan["admission_review"]
            self.assertIs(review["review_only"], True)
            self.assertEqual(review["model_count"], len(plan["model_dispositions"]))
            self.assertEqual(review["task_count"], len(plan["task_dispositions"]))
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "review-candidates",
            )
            self.assertFalse((candidate / "jobs").exists())
            self.assertEqual((candidate / "sbatch_headers.txt").read_text(), "")

    def test_stale_model_census_hash_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            campaign.write_text(
                campaign.read_text(encoding="utf-8").replace(
                    "[selection]",
                    '[selection]\ncensus_model_ids_sha256 = "' + "0" * 64 + '"',
                    1,
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(PlanningError, "frozen model census"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_production_submission_is_single_use_and_records_job_id(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            campaign.write_text(
                campaign.read_text(encoding="utf-8").replace(
                    "submit_enabled = false", "submit_enabled = true"
                ),
                encoding="utf-8",
            )
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "candidates",
            )
            approved = sha256_file(candidate / "ARTIFACTS.json")
            ledger_root = root / "submission-ledgers"
            with patch(
                "masld_bench.campaign.subprocess.run",
                return_value=subprocess.CompletedProcess(
                    args=[], returncode=0, stdout="123456;cluster\n", stderr=""
                ),
            ):
                commands = submission_commands(
                    candidate=candidate,
                    approved_campaign_sha256=approved,
                    execute=True,
                    submission_ledger_root=ledger_root,
                )
            self.assertEqual(len(commands), 1)
            ledger = next(path for path in ledger_root.iterdir() if path.is_dir())
            verify_frozen_tree(ledger)
            self.assertIn(
                '"job_id":"123456"',
                (ledger / "records" / "00001.json").read_text(encoding="utf-8"),
            )
            with self.assertRaisesRegex(CampaignError, "already has"):
                submission_commands(
                    candidate=candidate,
                    approved_campaign_sha256=approved,
                    execute=True,
                    submission_ledger_root=ledger_root,
                )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            campaign.write_text(
                campaign.read_text(encoding="utf-8").replace(
                    "submit_enabled = false", "submit_enabled = true"
                ),
                encoding="utf-8",
            )
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "candidates",
            )
            ledger_parent = root / "real-ledger-parent"
            ledger_parent.mkdir()
            ledger_parent_link = root / "ledger-parent-link"
            ledger_parent_link.symlink_to(ledger_parent, target_is_directory=True)
            with self.assertRaisesRegex(CampaignError, "symlink"):
                submission_commands(
                    candidate=candidate,
                    approved_campaign_sha256=sha256_file(
                        candidate / "ARTIFACTS.json"
                    ),
                    execute=True,
                    submission_ledger_root=ledger_parent_link / "ledgers",
                )

    def test_gpu_campaign_submission_requires_central_dispatcher(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            runtime = root / "config" / "runtimes" / "fixture.toml"
            runtime.write_text(
                runtime.read_text(encoding="utf-8").replace(
                    'resource_profile = "cpu_contract"',
                    'resource_profile = "gpu_single"',
                ),
                encoding="utf-8",
            )
            campaign.write_text(
                campaign.read_text(encoding="utf-8")
                .replace("submit_enabled = false", "submit_enabled = true")
                .replace(
                    'default_resource_profile = "cpu_contract"',
                    'default_resource_profile = "gpu_single"',
                ),
                encoding="utf-8",
            )
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "candidates",
            )
            approved = sha256_file(candidate / "ARTIFACTS.json")
            with patch("masld_bench.campaign.subprocess.run") as submit:
                with self.assertRaisesRegex(CampaignError, "dispatcher-only"):
                    submission_commands(
                        candidate=candidate,
                        approved_campaign_sha256=approved,
                        execute=True,
                        submission_ledger_root=root / "submission-ledgers",
                    )
            submit.assert_not_called()

    def test_scientific_run_requires_and_executes_frozen_actions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "candidates",
            )
            frozen = load_frozen_plan(candidate)
            self.assertEqual(len(frozen["runs"]), 1)
            self.assertEqual(
                frozen["runs"][0]["action"], ["prepare", "fit", "predict"]
            )
            attempt = execute_run(
                candidate=candidate,
                run_id=frozen["runs"][0]["run_id"],
                output_root=root / "executions",
            )
            receipt = verify_run_execution_attempt(attempt)
            self.assertEqual(receipt["status"], "succeeded")
            self.assertEqual(
                [item["action"] for item in receipt["action_receipts"]],
                ["prepare", "fit", "predict"],
            )
            requests = sorted((attempt / "requests").glob("*.json"))
            self.assertEqual(len(requests), 3)
            payloads = [
                json.loads(path.read_text(encoding="utf-8")) for path in requests
            ]
            for payload in payloads:
                scoped_roles = {
                    item["role"] for item in payload["run_spec"]["inputs"]
                }
                owners = payload["run_spec"]["metadata"]["input_owner_by_role"]
                allowed = set(payload["action_dataset_ids"])
                self.assertEqual(
                    scoped_roles,
                    {
                        role
                        for role, owner in owners.items()
                        if owner is None or owner in allowed
                    },
                )
            self.assertEqual(
                [len(item["prior_action_outputs"]) for item in payloads],
                [0, 1, 2],
            )
            for action_index, payload in enumerate(payloads, start=1):
                self.assertEqual(
                    [item["action"] for item in payload["prior_action_outputs"]],
                    frozen["runs"][0]["action"][: action_index - 1],
                )
            fit_payload = next(item for item in payloads if item["action"] == "fit")
            self.assertEqual(
                fit_payload["action_dataset_ids"],
                frozen["runs"][0]["metadata"]["fit_dataset_ids"],
            )
            self.assertIn(
                "dataset:fixture_development",
                {item["role"] for item in fit_payload["run_spec"]["inputs"]},
            )
            with self.assertRaisesRegex(CampaignError, "prior attempt"):
                execute_run(
                    candidate=candidate,
                    run_id=frozen["runs"][0]["run_id"],
                    output_root=root / "executions",
                )

    def test_scientific_plan_omits_unsupported_adaptation_lanes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            campaign_document = campaign.read_text(encoding="utf-8")
            campaign.write_text(
                campaign_document.replace(
                    'adaptation_regimes = ["common_lane"]',
                    'adaptation_regime = "common_and_native"',
                    1,
                ),
                encoding="utf-8",
            )
            model = root / "config" / "models" / "fixture_family.toml"
            model_document = model.read_text(encoding="utf-8")
            model.write_text(
                model_document.replace(
                    'supported_adaptation_regimes = ["common_lane", "native_lane"]',
                    'supported_adaptation_regimes = ["common_lane"]',
                    1,
                ),
                encoding="utf-8",
            )
            plan = build_plan(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
            )
            self.assertEqual(
                {run["adaptation_regime"] for run in plan["runs"]},
                {"common_lane"},
            )
            disposition = next(
                item
                for item in plan["model_dispositions"]
                if item["model_id"] == "fixture_open"
            )
            self.assertEqual(
                disposition["omitted_adaptation_regimes"], ["native_lane"]
            )
            self.assertEqual(
                disposition["executable_adaptation_regimes"], ["common_lane"]
            )

    def test_common_heads_are_separate_runs_and_do_not_expand_native_lane(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            campaign.write_text(
                campaign.read_text(encoding="utf-8")
                .replace(
                    'adaptation_regimes = ["common_lane"]',
                    'adaptation_regime = "common_and_native"\n'
                    'common_head_ids = ["linear", "two_layer_mlp"]',
                    1,
                ),
                encoding="utf-8",
            )

            plan = build_plan(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
            )

            self.assertEqual(len(plan["runs"]), 3)
            common_runs = [
                run
                for run in plan["runs"]
                if run["adaptation_regime"] == "common_lane"
            ]
            native_runs = [
                run
                for run in plan["runs"]
                if run["adaptation_regime"] == "native_lane"
            ]
            self.assertEqual(
                {run["hyperparameters"]["common_head_id"] for run in common_runs},
                {"linear", "two_layer_mlp"},
            )
            self.assertEqual(
                {run["metadata"]["common_head_id"] for run in common_runs},
                {"linear", "two_layer_mlp"},
            )
            self.assertEqual(len({run["run_id"] for run in common_runs}), 2)
            self.assertEqual(len(native_runs), 1)
            self.assertNotIn("common_head_id", native_runs[0]["hyperparameters"])
            self.assertIsNone(native_runs[0]["metadata"]["common_head_id"])

    def test_common_head_configuration_rejects_ambiguous_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            campaign.write_text(
                campaign.read_text(encoding="utf-8")
                .replace(
                    'adaptation_regimes = ["common_lane"]',
                    'adaptation_regimes = ["common_lane"]\n'
                    'common_head_ids = ["linear"]',
                    1,
                )
                .replace(
                    "[execution.adapter_request]",
                    '[execution.hyperparameters]\ncommon_head_id = "linear"\n'
                    "[execution.adapter_request]",
                    1,
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(PlanningError, "conflicts"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_unfreezable_adapter_output_gets_separate_terminal_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "candidates",
            )
            frozen = load_frozen_plan(candidate)
            run_id = frozen["runs"][0]["run_id"]

            def leave_unfreezable_output(_action, _request_path, output):
                output.mkdir(parents=True, exist_ok=False)
                (output / "outside-link").symlink_to(root, target_is_directory=True)
                raise RuntimeError("fixture adapter stopped after unsafe output")

            with patch(
                "masld_bench.campaign.SubprocessAdapter.invoke",
                side_effect=leave_unfreezable_output,
            ):
                with self.assertRaisesRegex(
                    CampaignError, "minimal terminal receipt"
                ):
                    execute_admission_run(
                        candidate=candidate,
                        run_id=run_id,
                        output_root=root / "executions",
                    )

            attempt = root / "executions" / "runs" / run_id / "attempt-001"
            receipt = verify_run_execution_attempt(
                attempt, require_succeeded=False
            )
            self.assertEqual(receipt["status"], "invalid_output")
            self.assertEqual(receipt["action_receipts"], [])
            context = json.loads(
                (attempt / "failure_context.json").read_text(encoding="utf-8")
            )
            quarantine = Path(context["untrusted_quarantine_path"])
            self.assertTrue(context["untrusted_quarantine_published"])
            self.assertTrue(quarantine.is_dir())
            self.assertTrue(
                (quarantine / "adapter_actions" / "001-probe" / "outside-link").is_symlink()
            )

    def test_scientific_plan_fails_closed_without_dataset_activation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            dataset = root / "config" / "datasets" / "fixture_development.toml"
            document = dataset.read_text(encoding="utf-8")
            document = document[: document.index("[activation]")]
            document = document.replace(
                "admission_blocking = false",
                "admission_blocking = true",
                1,
            ).replace("blockers = []", 'blockers = ["activation_unresolved"]', 1)
            dataset.write_text(document, encoding="utf-8")
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            with self.assertRaisesRegex(PlanningError, "activation"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_smoke_view_activates_only_the_frozen_parent_view(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            admission_campaign = _fixture_config(root)
            view_id = _add_fixture_dataset_view(root)
            admission_candidate = freeze_campaign(
                admission_campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "admission-candidates",
            )
            admission_plan = load_frozen_plan(admission_candidate)
            prerequisite_attempt = execute_admission_run(
                candidate=admission_candidate,
                run_id=admission_plan["runs"][0]["run_id"],
                output_root=root / "admission-executions",
            )
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            campaign.write_text(
                campaign.read_text(encoding="utf-8")
                .replace(
                    'task_ids = ["fixture_task"]',
                    f'task_ids = ["fixture_task"]\ndataset_view_ids = ["{view_id}"]',
                    1,
                )
                .replace("folds = [0]", "folds = [0]\ncell_budget = 2", 1),
                encoding="utf-8",
            )
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "scientific-candidates",
            )
            plan = load_frozen_plan(candidate)
            lock = plan["dataset_locks"]["fixture_development"]
            self.assertFalse(lock["parent_activation_ready"])
            self.assertTrue(lock["activation_ready"])
            self.assertEqual(lock["activation_source"], "dataset_view")
            self.assertEqual(lock["dataset_view"]["view_id"], view_id)
            run = plan["runs"][0]
            self.assertEqual(run["dataset_ids"], ["fixture_development"])
            self.assertEqual(
                run["metadata"]["dataset_routing"]["fixture_development"][
                    "dataset_view_id"
                ],
                view_id,
            )
            input_roles = {item["role"] for item in run["inputs"]}
            self.assertIn(f"dataset_view_data:{view_id}", input_roles)
            self.assertIn(f"dataset_view_manifest:{view_id}", input_roles)
            attempt = execute_run(
                candidate=candidate,
                run_id=run["run_id"],
                output_root=root / "scientific-executions",
            )
            self.assertEqual(
                verify_run_execution_attempt(attempt)["status"], "succeeded"
            )

            data = (
                root
                / "config"
                / "artifacts"
                / "dataset-views"
                / view_id
                / "data.json"
            )
            data.write_text('{"tampered":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(PlanningError, "artifact verification"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_incremental_smoke_can_run_one_explicit_mandatory_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            campaign.write_text(
                campaign.read_text(encoding="utf-8")
                .replace(
                    'model_statuses = ["candidate"]',
                    'model_statuses = ["candidate"]\nmodel_ids = ["fixture_open"]',
                    1,
                )
                .replace(
                    "include_mandatory_baselines = true",
                    "include_mandatory_baselines = false",
                    1,
                ),
                encoding="utf-8",
            )

            plan = build_plan(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
            )

            self.assertEqual(len(plan["runs"]), 1)
            self.assertEqual(plan["runs"][0]["model_id"], "fixture_open")
            self.assertEqual(plan["runs"][0]["task_id"], "fixture_task")

    def test_incremental_smoke_rejects_nonbaseline_only_request(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            admission_campaign = _fixture_config(root)
            task = root / "config" / "tasks" / "fixture_task.toml"
            task.write_text(
                task.read_text(encoding="utf-8").replace(
                    'baseline_model_ids = ["fixture_open"]',
                    'baseline_model_ids = ["blocked_unknown"]',
                    1,
                ),
                encoding="utf-8",
            )
            admission_candidate = freeze_campaign(
                admission_campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
                output_root=root / "admission-candidates",
            )
            admission_plan = load_frozen_plan(admission_candidate)
            prerequisite_attempt = execute_admission_run(
                candidate=admission_candidate,
                run_id=admission_plan["runs"][0]["run_id"],
                output_root=root / "admission-executions",
            )
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            campaign.write_text(
                campaign.read_text(encoding="utf-8")
                .replace(
                    'model_statuses = ["candidate"]',
                    'model_statuses = ["candidate"]\nmodel_ids = ["fixture_open"]',
                    1,
                )
                .replace(
                    "include_mandatory_baselines = true",
                    "include_mandatory_baselines = false",
                    1,
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                PlanningError, "cannot drop mandatory baseline runs"
            ):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_scientific_plan_rejects_tampered_dataset_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            artifact = root / "config" / "artifacts" / "fixture_development.manifest.json"
            artifact.write_text('{"tampered":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(PlanningError, "artifact verification"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_unavailable_optional_development_source_does_not_block_task(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prerequisite_attempt = _execute_fixture_admission(root)
            _write(
                root / "config" / "datasets" / "fixture_optional.toml",
                """
                schema_version = "masld-bench-dataset-v1"
                dataset_id = "fixture_optional"
                title = "Unavailable optional development source"
                role = "blocked"
                status = "blocked"
                species = "human"
                accession = ["FIXTURE-OPTIONAL"]
                access_tier = "conditional"
                automatic_download = false
                redistribution_class = "fixture_only"
                license_status = "UNRESOLVED"
                biological_unit = "donor"
                expected_biological_units = "UNRESOLVED"
                native_genome_build = "UNRESOLVED"
                native_annotation_release = "UNRESOLVED"
                modalities = ["single_cell_rna"]
                pairing_levels = ["same_study_unpaired"]
                modality_status_default = "unavailable_permission"
                exposure_status = "unknown"
                admission_blocking = true
                blockers = ["optional_source_not_activated"]
                notes = "Synthetic blocked optional source."
                [provenance]
                source_url = "https://example.invalid/optional"
                verified_date = "2026-08-21"
                local_exposure = "none"
                [split]
                group_key = "cohort_family_id+donor_id"
                outer_role = "development"
                label_visibility = "unavailable_permission"
                """,
            )
            task = root / "config" / "tasks" / "fixture_task.toml"
            task.write_text(
                task.read_text(encoding="utf-8").replace(
                    "datasets_development = []",
                    'datasets_development = ["fixture_optional"]',
                ),
                encoding="utf-8",
            )
            campaign = _fixture_scientific_campaign(root, prerequisite_attempt)
            plan = build_plan(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
            )
            disposition = plan["task_dispositions"][0]
            self.assertEqual(disposition["disposition"], "eligible_partial")
            self.assertEqual(disposition["run_dataset_ids"], ["fixture_development"])
            self.assertEqual(
                disposition["omitted_datasets"],
                [
                    {
                        "dataset_id": "fixture_optional",
                        "reason": "dataset_activation_not_ready:fixture_optional",
                    }
                ],
            )
            self.assertTrue(plan["runs"])

    def test_ready_dataset_and_model_require_concrete_authority_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            dataset = root / "config" / "datasets" / "fixture_development.toml"
            dataset.write_text(
                "\n".join(
                    line
                    for line in dataset.read_text(encoding="utf-8").splitlines()
                    if not line.strip().startswith("evidence_artifacts =")
                )
                + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(PlanningError, "evidence_artifacts"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            model = root / "config" / "models" / "fixture_family.toml"
            model.write_text(
                "\n".join(
                    line
                    for line in model.read_text(encoding="utf-8").splitlines()
                    if not line.strip().startswith("evidence_artifacts =")
                )
                + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(PlanningError, "model execution"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_open_champion_decision_keeps_license_dimensions_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            _rewrite_fixture_license_authority(
                root,
                authority="weights_license",
                old_license="not_applicable_classical",
                new_license="research_only",
                use_allowed=True,
                redistribution_allowed=False,
                derivative_redistribution_allowed=False,
            )
            plan = build_plan(
                campaign,
                config_root=root / "config",
                package_root=PACKAGE_ROOT,
            )
            disposition = next(
                item
                for item in plan["model_dispositions"]
                if item["model_id"] == "fixture_open"
            )
            decision = disposition["open_champion_license_decision"]
            self.assertEqual(disposition["disposition"], "comparator_only")
            self.assertFalse(disposition["open_champion_eligible"])
            self.assertTrue(decision["code_use_open"])
            self.assertFalse(decision["weights_redistributable"])
            self.assertTrue(decision["derivative_weights_redistributable"])

    def test_ready_model_rejects_unrecognized_license_term(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            model = root / "config" / "models" / "fixture_family.toml"
            model.write_text(
                model.read_text(encoding="utf-8").replace(
                    'weights_license = "not_applicable_classical"',
                    'weights_license = "mystery-license"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(PlanningError, "unrecognized weights_license"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_ready_model_rejects_arbitrary_hashed_license_authority(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign = _fixture_config(root)
            authority = (
                root
                / "config"
                / "artifacts"
                / "model-authorities"
                / "weights_license.txt"
            )
            old_sha256 = sha256_file(authority)
            old_size = authority.stat().st_size
            _write(authority, '{"arbitrary":true}')
            new_sha256 = sha256_file(authority)
            new_size = authority.stat().st_size
            model = root / "config" / "models" / "fixture_family.toml"
            document = model.read_text(encoding="utf-8").replace(
                (
                    'path = "artifacts/model-authorities/weights_license.txt", '
                    f'sha256 = "{old_sha256}", size_bytes = {old_size}'
                ),
                (
                    'path = "artifacts/model-authorities/weights_license.txt", '
                    f'sha256 = "{new_sha256}", size_bytes = {new_size}'
                ),
                1,
            )
            model.write_text(document, encoding="utf-8")
            with self.assertRaisesRegex(PlanningError, "authority content differs"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=PACKAGE_ROOT,
                )

    def test_supplementary_toml_changes_registry_run_and_plan_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory(
            prefix=".planner-source-", dir=_test_source_parent()
        ) as source_temporary:
            root = Path(temporary)
            source_root = Path(source_temporary)
            _write(source_root / "src" / "fixture.py", "VALUE = 1")
            campaign = _fixture_config(root)
            evaluation = root / "config" / "evaluation" / "supplementary.toml"
            _write(
                evaluation,
                """
                schema_version = "fixture-evaluation-v1"
                threshold = 0.1
                """,
            )
            before = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )
            manifest = {
                record["path"]: record
                for record in before["registry_snapshot"]["files"]
            }
            self.assertIn("evaluation/supplementary.toml", manifest)
            self.assertEqual(
                before["registry_snapshot_sha256"],
                canonical_hash(before["registry_snapshot"]),
            )
            run_before = before["runs"][0]
            self.assertEqual(
                run_before["metadata"]["registry_snapshot_sha256"],
                before["registry_snapshot_sha256"],
            )
            self.assertEqual(
                run_before["metadata"]["campaign_file_sha256"],
                before["campaign_file"]["sha256"],
            )
            self.assertEqual(
                run_before["code_lock_sha256"], before["source_lock_sha256"]
            )

            _write(
                evaluation,
                """
                schema_version = "fixture-evaluation-v1"
                threshold = 0.2
                """,
            )
            after = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )

            self.assertEqual(
                before["core_registry_contract_sha256"],
                after["core_registry_contract_sha256"],
            )
            self.assertEqual(
                before["campaign_file"]["sha256"], after["campaign_file"]["sha256"]
            )
            self.assertEqual(before["source_lock_sha256"], after["source_lock_sha256"])
            self.assertNotEqual(
                before["registry_snapshot_sha256"], after["registry_snapshot_sha256"]
            )
            self.assertNotEqual(before["runs"][0]["run_id"], after["runs"][0]["run_id"])
            self.assertNotEqual(before["plan_sha256"], after["plan_sha256"])

    def test_supplementary_json_changes_registry_run_and_plan_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory(
            prefix=".planner-source-", dir=_test_source_parent()
        ) as source_temporary:
            root = Path(temporary)
            source_root = Path(source_temporary)
            _write(source_root / "src" / "fixture.py", "VALUE = 1")
            campaign = _fixture_config(root)
            audit = root / "config" / "artifacts" / "models" / "fixture" / "audit.json"
            _write(audit, '{"schema_version":"fixture-audit-v1","value":1}')
            before = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )
            manifest = {
                record["path"]: record
                for record in before["registry_snapshot"]["files"]
            }
            self.assertIn("artifacts/models/fixture/audit.json", manifest)

            _write(audit, '{"schema_version":"fixture-audit-v1","value":2}')
            after = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )
            self.assertNotEqual(
                before["registry_snapshot_sha256"],
                after["registry_snapshot_sha256"],
            )
            self.assertNotEqual(before["runs"][0]["run_id"], after["runs"][0]["run_id"])
            self.assertNotEqual(before["plan_sha256"], after["plan_sha256"])

    def test_source_change_changes_run_identity_but_outputs_do_not(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory(
            prefix=".planner-source-", dir=_test_source_parent()
        ) as source_temporary:
            root = Path(temporary)
            source_root = Path(source_temporary)
            code = source_root / "src" / "fixture.py"
            _write(code, "VALUE = 1")
            campaign = _fixture_config(root)
            before = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )

            source_payload = dict(before["source_lock"])
            claimed_source_sha256 = source_payload.pop("source_lock_sha256")
            self.assertEqual(claimed_source_sha256, canonical_hash(source_payload))
            self.assertEqual(before["source_lock_sha256"], claimed_source_sha256)
            self.assertIn("git_commit", source_payload)
            self.assertIn("dirty_patch_sha256", source_payload)
            self.assertIn("tree_manifest", source_payload)

            _write(code, "VALUE = 2")
            changed = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )
            self.assertEqual(
                before["registry_snapshot_sha256"],
                changed["registry_snapshot_sha256"],
            )
            self.assertEqual(
                before["campaign_file"]["sha256"],
                changed["campaign_file"]["sha256"],
            )
            self.assertEqual(
                before["source_lock"]["git_commit"],
                changed["source_lock"]["git_commit"],
            )
            self.assertNotEqual(
                before["source_lock"]["tree_sha256"],
                changed["source_lock"]["tree_sha256"],
            )
            self.assertNotEqual(before["source_lock_sha256"], changed["source_lock_sha256"])
            self.assertNotEqual(before["runs"][0]["run_id"], changed["runs"][0]["run_id"])
            self.assertNotEqual(before["plan_sha256"], changed["plan_sha256"])

            with self.assertRaisesRegex(PlanningError, "provenance-excluded"):
                freeze_campaign(
                    campaign,
                    config_root=root / "config",
                    package_root=source_root,
                    output_root=source_root / "untracked-output",
                )
            candidate = freeze_campaign(
                campaign,
                config_root=root / "config",
                package_root=source_root,
                output_root=root / "candidates",
            )
            frozen = load_frozen_plan(candidate)
            _write(source_root / "selections" / "fixture.json", "{}")
            after_output = build_plan(
                campaign,
                config_root=root / "config",
                package_root=source_root,
            )
            self.assertEqual(frozen["source_lock"], after_output["source_lock"])
            self.assertEqual(frozen["runs"][0]["run_id"], after_output["runs"][0]["run_id"])
            self.assertEqual(frozen["plan_sha256"], after_output["plan_sha256"])

    def test_scoped_cobolt_lock_ignores_queue_but_detects_bound_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory(
            prefix=".planner-source-", dir=_test_source_parent()
        ) as source_temporary:
            root = Path(temporary)
            source_root = Path(source_temporary)
            bound_source = source_root / "src" / "cobolt_bound.py"
            _write(bound_source, "VALUE = 1")
            campaign = _fixture_config(root)
            config = root / "config"
            registry_paths = sorted(
                path.relative_to(config).as_posix()
                for path in config.rglob("*")
                if path.is_file()
            )
            document = campaign.read_text(encoding="utf-8").replace(
                'task_ids = ["fixture_task"]',
                '\n'.join(
                    (
                        'task_ids = ["fixture_task"]',
                        'source_lock_paths = ["src/cobolt_bound.py"]',
                        "registry_snapshot_paths = "
                        + json.dumps(registry_paths, separators=(",", ":")),
                    )
                ),
                1,
            )
            campaign.write_text(document, encoding="utf-8")

            plan = build_plan(
                campaign,
                config_root=config,
                package_root=source_root,
            )
            self.assertEqual(
                plan["source_lock"]["schema_version"],
                "masld-bench-source-lock-v3",
            )
            self.assertEqual(
                plan["registry_snapshot"]["schema_version"],
                "masld-bench-registry-file-snapshot-v2",
            )

            _write(
                config / "campaigns" / "gpu_bundle_queue" / "unrelated.json",
                '{"model_id":"unrelated"}',
            )
            _write(
                source_root
                / "config"
                / "campaigns"
                / "gpu_bundle_queue"
                / "unrelated.json",
                '{"model_id":"unrelated"}',
            )
            _assert_config_snapshot(plan)
            self.assertEqual(
                _observed_source_lock(plan, source_root), plan["source_lock"]
            )

            _write(bound_source, "VALUE = 2")
            self.assertNotEqual(
                _observed_source_lock(plan, source_root), plan["source_lock"]
            )

    def test_configuration_snapshot_rejects_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory(
            prefix=".planner-source-", dir=_test_source_parent()
        ) as source_temporary:
            root = Path(temporary)
            source_root = Path(source_temporary)
            _write(source_root / "src" / "fixture.py", "VALUE = 1")
            campaign = _fixture_config(root)
            target = root / "outside.toml"
            _write(target, 'schema_version = "outside-v1"')
            link = root / "config" / "evaluation" / "linked.toml"
            link.parent.mkdir(parents=True, exist_ok=True)
            link.symlink_to(target)

            with self.assertRaisesRegex(PlanningError, "symlink"):
                build_plan(
                    campaign,
                    config_root=root / "config",
                    package_root=source_root,
                )


if __name__ == "__main__":
    unittest.main()
