"""Contract tests for the Cobolt v10 recovery revision.

GPU job 21100294 ran all 25 Cobolt v9 runs and failed every one at run
admission, before any adapter action, because concurrent Figure 4
Resource-authority promotions changed protected files while the job held the
frozen v9 firewall baseline.  The frozen continuation policy forbids mid-run
resume and demands a new campaign revision, a distinct immutable bundle, and a
distinct queue item.

These tests pin the three properties that make the recovery trustworthy:
the revision changed the firewall baseline and nothing else scientific, the
new identities are disjoint from the failed ones, and the firewall itself was
re-baselined rather than weakened.
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import unittest

from scripts.gpu_bundle_dispatcher import load_item
from scripts.gpu_dependency_dispatcher import validate_item_route


ROOT = Path(__file__).resolve().parents[2]
BUNDLE_ID = "cobolt-5seed-25run-v10-generic7"
SUPERSEDED_BUNDLE_ID = "cobolt-5seed-25run-v9-generic6"
JOB_NAME = "model-work-354"
PRIORITY = 102

WRAPPER = ROOT / "slurm/run_cobolt_5seed_25run_v10_generic7_bundle.sbatch"
STAGED_ITEM = (
    ROOT / "config/campaigns/gpu_bundle_queue_staging" / f"{PRIORITY}-{BUNDLE_ID}.json"
)
LIVE_ITEM = ROOT / "config/campaigns/gpu_bundle_queue" / f"{PRIORITY}-{BUNDLE_ID}.json"
BUNDLE = ROOT / "executions/dispatch-bundles" / BUNDLE_ID
SUPERSEDED_BUNDLE = ROOT / "executions/dispatch-bundles" / SUPERSEDED_BUNDLE_ID
SUPERSEDED_CANDIDATE = (
    ROOT / "candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58"
)
REVISION = (
    ROOT / "config/campaigns/v1_rna_atac_cobolt_5seed_smoke_revision10_20260825.json"
)
INCIDENT = (
    ROOT
    / "config/artifacts/incidents/cobolt_v9_resource_authority_firewall_21100294.json"
)
DISPATCH_POLICY = (
    ROOT / "config/campaigns/gpu_dependency_dispatcher_policy_20260825.json"
)
FAILED_JOB_ID = 21100294
SEEDS = (1103, 2107, 3109, 4111, 5113)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=4)
def load_plan(candidate: Path) -> dict:
    """Frozen plans carry the whole 150k-file firewall snapshot, so parse once."""
    return read_json(candidate / "plan.json")


class CoboltV10QueueItemTest(unittest.TestCase):
    def test_staged_item_carries_exactly_the_queue_schema_keys(self) -> None:
        value = read_json(STAGED_ITEM)
        self.assertEqual(
            set(value),
            {
                "schema_version",
                "bundle_id",
                "priority",
                "enabled",
                "wrapper_path",
                "wrapper_sha256",
                "exports",
                "required_paths",
                "logical_tasks",
                "family",
            },
        )
        self.assertEqual(value["schema_version"], "masld-bench-gpu-bundle-queue-item-v1")
        self.assertEqual(value["bundle_id"], BUNDLE_ID)
        self.assertEqual(value["priority"], PRIORITY)
        self.assertIs(value["enabled"], True)
        self.assertEqual(value["logical_tasks"], 25)
        self.assertEqual(value["family"], "observed_multiome")
        self.assertEqual(value["exports"], {})
        self.assertEqual(
            value["wrapper_path"],
            "slurm/run_cobolt_5seed_25run_v10_generic7_bundle.sbatch",
        )
        self.assertEqual(value["wrapper_sha256"], digest(WRAPPER))

    def test_staged_item_is_ready_and_routes_to_a_generic_l40s_lane(self) -> None:
        item = validate_item_route(load_item(STAGED_ITEM, ROOT.resolve(strict=True)))
        self.assertIs(item["_ready"], True)
        self.assertEqual(item["_job_name"], JOB_NAME)
        self.assertEqual(item["_hours"], 6.0)

    def test_required_paths_bind_the_new_bundle_candidate_and_incident(self) -> None:
        value = read_json(STAGED_ITEM)
        required = value["required_paths"]
        self.assertIn(f"executions/dispatch-bundles/{BUNDLE_ID}/ARTIFACTS.json", required)
        self.assertIn(
            "config/artifacts/incidents/"
            "cobolt_v9_resource_authority_firewall_21100294.json",
            required,
        )
        for relative in required:
            self.assertTrue((ROOT / relative).is_file(), relative)

    def test_queue_identity_is_distinct_from_the_failed_bundle(self) -> None:
        self.assertNotEqual(BUNDLE_ID, SUPERSEDED_BUNDLE_ID)
        self.assertNotEqual(STAGED_ITEM.name, f"057-{SUPERSEDED_BUNDLE_ID}.json")
        superseded_item = (
            ROOT / "config/campaigns/gpu_bundle_queue/057-cobolt-5seed-25run-v9-generic6.json"
        )
        self.assertNotEqual(
            read_json(superseded_item)["wrapper_sha256"], digest(WRAPPER)
        )


class CoboltV10WrapperTest(unittest.TestCase):
    def setUp(self) -> None:
        self.text = WRAPPER.read_text(encoding="utf-8")

    def test_wrapper_is_exact_nslab_l40s_allocation_with_a_generic_name(self) -> None:
        for expected in (
            f"#SBATCH --job-name={JOB_NAME}",
            "#SBATCH --account=nslab",
            "#SBATCH --partition=gpu",
            "#SBATCH --qos=nslab",
            "#SBATCH --gres=gpu:l40s:1",
            "#SBATCH --cpus-per-task=4",
            "#SBATCH --mem=64G",
            "#SBATCH --time=06:00:00",
        ):
            self.assertIn(expected, self.text)
        self.assertNotIn("#SBATCH --array", self.text)
        self.assertNotIn("innovation", self.text.lower())
        name_line = next(
            line for line in self.text.splitlines() if "--job-name=" in line
        )
        self.assertNotIn("masld", name_line.lower())
        self.assertNotIn("cobolt", name_line.lower())

    def test_wrapper_never_submits_another_job(self) -> None:
        body = "\n".join(
            line for line in self.text.splitlines() if not line.lstrip().startswith("#")
        )
        self.assertNotIn("sbatch", body)

    def test_wrapper_reverifies_the_frozen_sources_before_any_run(self) -> None:
        self.assertIn("source_lock(root, included_paths=included_paths)", self.text)
        self.assertIn("zero runs started", self.text)
        self.assertIn("sha256sum --check --strict", self.text)
        self.assertIn("--execution-id", self.text)
        self.assertIn(BUNDLE_ID, self.text)

    def test_wrapper_pins_the_frozen_manifest_and_runner_digests(self) -> None:
        manifest = BUNDLE / "bundle_manifest.json"
        self.assertIn(digest(manifest), self.text)
        self.assertIn(read_json(manifest)["runner_sha256"], self.text)


class CoboltV10BundleTest(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = read_json(BUNDLE / "bundle_manifest.json")
        self.candidate = Path(self.manifest["candidate"])
        self.plan = load_plan(self.candidate)
        self.superseded_plan = load_plan(SUPERSEDED_CANDIDATE)

    def test_bundle_carries_the_unchanged_five_seed_by_five_fold_grid(self) -> None:
        observed = [(int(r["seed"]), int(r["fold"])) for r in self.manifest["runs"]]
        self.assertEqual(
            observed, [(seed, fold) for seed in SEEDS for fold in range(5)]
        )
        self.assertEqual(len(self.manifest["runs"]), 25)
        self.assertEqual(
            [r["index"] for r in self.manifest["runs"]], list(range(1, 26))
        )

    def test_continuation_policy_still_forbids_mid_run_resume(self) -> None:
        policy = self.manifest["continuation_policy"]
        self.assertIs(policy["mid_run_resume"], False)
        self.assertEqual(
            policy["interrupted_run_recovery"], "new_campaign_revision_required"
        )
        self.assertEqual(
            policy["bundle_retry_identity"],
            "distinct_immutable_bundle_and_queue_item_required",
        )

    def test_run_identities_and_execution_root_are_disjoint_from_the_failed_attempt(
        self,
    ) -> None:
        superseded_ids = {str(r["run_id"]) for r in self.superseded_plan["runs"]}
        new_ids = {str(r["run_id"]) for r in self.manifest["runs"]}
        self.assertEqual(len(new_ids), 25)
        self.assertEqual(superseded_ids & new_ids, set())
        superseded_root = (
            SUPERSEDED_CANDIDATE.parent / "executions" / SUPERSEDED_CANDIDATE.name
        )
        self.assertNotEqual(
            Path(self.manifest["execution_output_root"]), superseded_root
        )
        self.assertTrue(superseded_root.is_dir(), "failed run roots must be preserved")

    def test_resources_are_unchanged_from_the_failed_bundle(self) -> None:
        superseded = read_json(SUPERSEDED_BUNDLE / "bundle_manifest.json")
        self.assertEqual(self.manifest["resources"], superseded["resources"])
        self.assertEqual(self.manifest["runner_sha256"], superseded["runner_sha256"])
        for key in (
            "expected_resource_profile",
            "expected_runtime_id",
            "expected_model_id",
            "expected_task_id",
            "expected_actions",
        ):
            self.assertEqual(self.manifest[key], superseded[key])


class CoboltV10RevisionScopeTest(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = read_json(BUNDLE / "bundle_manifest.json")
        self.plan = load_plan(Path(self.manifest["candidate"]))
        self.superseded_plan = load_plan(SUPERSEDED_CANDIDATE)

    def test_revision_is_a_distinct_campaign(self) -> None:
        self.assertEqual(
            self.plan["campaign"]["campaign_id"], "v1-rna-atac-cobolt-5seed-smoke-v10"
        )
        self.assertNotEqual(
            self.plan["plan_sha256"], self.superseded_plan["plan_sha256"]
        )

    def test_scoped_source_lock_is_byte_identical(self) -> None:
        self.assertEqual(self.plan["source_lock"], self.superseded_plan["source_lock"])

    def test_seeds_folds_hyperparameters_and_bindings_are_unchanged(self) -> None:
        def shape(plan: dict) -> list[str]:
            rows = []
            for run in plan["runs"]:
                row = {k: v for k, v in run.items() if k not in ("run_id", "campaign_id")}
                metadata = dict(row["metadata"])
                metadata.pop("campaign_file_sha256", None)
                row["metadata"] = metadata
                rows.append(json.dumps(row, sort_keys=True))
            return sorted(rows)

        self.assertEqual(shape(self.plan), shape(self.superseded_plan))
        self.assertEqual(
            self.plan["resource_totals"], self.superseded_plan["resource_totals"]
        )
        self.assertEqual(self.plan["safety"], self.superseded_plan["safety"])

    def test_firewall_was_rebaselined_and_not_weakened(self) -> None:
        new = self.plan["resource_firewall"]
        old = self.superseded_plan["resource_firewall"]
        self.assertEqual(new["schema_version"], "masld-resource-firewall-v2")
        self.assertNotEqual(new["snapshot_sha256"], old["snapshot_sha256"])
        self.assertEqual(new["authority_roots"], old["authority_roots"])
        self.assertGreaterEqual(len(new["files"]), 1)
        old_paths = {item["path"] for item in old["files"]}
        new_paths = {item["path"] for item in new["files"]}
        # Re-baselining may add promoted files; it must never stop watching a
        # directory that was protected before.
        self.assertEqual(
            {Path(p).parts[0] for p in old_paths} - {Path(p).parts[0] for p in new_paths},
            set(),
        )

    def test_revision_record_names_the_failed_job_and_its_scope(self) -> None:
        value = read_json(REVISION)
        self.assertEqual(value["revision"], 10)
        self.assertEqual(value["supersedes_failed_job_ids"], [FAILED_JOB_ID])
        self.assertIn("firewall", value["revision_scope"].lower())
        self.assertIs(value["manual_gpu_sbatch_allowed"], False)
        self.assertIs(value["firewall"]["watched_authority_roots_reduced"], False)
        self.assertIs(value["firewall"]["mid_run_resume"], False)


class CoboltV10IncidentTest(unittest.TestCase):
    def test_incident_records_a_resource_side_delta_and_zero_model_work(self) -> None:
        value = read_json(INCIDENT)
        self.assertEqual(value["failed_job_id"], str(FAILED_JOB_ID))
        self.assertEqual(value["observed_failure"]["runs_failed"], 25)
        self.assertEqual(value["observed_failure"]["runs_succeeded"], 0)
        self.assertEqual(value["observed_failure"]["adapter_actions_executed"], 0)
        self.assertIs(value["attribution"]["delta_is_resource_side"], True)
        self.assertIs(value["attribution"]["delta_caused_by_any_model_job"], False)
        self.assertEqual(value["attribution"]["changed_files_under_benchmark_root"], 0)
        self.assertIs(value["model_fit_started"], False)
        self.assertIs(value["metrics_calculated"], False)
        self.assertIs(value["development_outcomes_read"], False)
        self.assertIs(value["sealed_outcomes_read"], False)
        self.assertIs(value["correction"]["firewall_weakened"], False)
        self.assertIs(
            value["continuation_policy_applied"]["failed_run_roots_reused"], False
        )
        self.assertIs(
            value["continuation_policy_applied"]["failed_run_roots_deleted"], False
        )

    def test_no_changed_protected_file_lives_under_the_benchmark_root(self) -> None:
        value = read_json(INCIDENT)
        for relative in value["changed_protected_files"]["paths"]:
            self.assertFalse(relative.startswith("Analysis/MASLD_Model_Benchmark"))


class CoboltV10DispatchRouteTest(unittest.TestCase):
    def test_dispatcher_source_bindings_are_untouched_by_this_revision(self) -> None:
        policy = read_json(DISPATCH_POLICY)
        for label, record in policy["source_bindings"].items():
            path = ROOT / record["path"]
            self.assertEqual(digest(path), record["sha256"], label)

    def test_publication_is_a_hard_link_from_staging_not_a_manual_gpu_sbatch(
        self,
    ) -> None:
        validator = ROOT / "slurm/validate_cobolt_v10_generic7_queue_cpu.sbatch"
        text = validator.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=cpu", text)
        self.assertIn("os.link(staged, live", text)
        body = "\n".join(
            line for line in text.splitlines() if not line.lstrip().startswith("#")
        )
        self.assertNotIn("sbatch", body)


if __name__ == "__main__":
    unittest.main()
