"""The NO_TASKSPEC blocker was re-grounded, not lifted.

On 2026-08-24 the observed_multiome_integration record asserted that no
observed-multiome TaskSpec existed.  A standalone TaskSpec was written to
config/evaluation on 2026-08-25 at 01:58:03, after the record, so that clause
became false.  Correcting it is only safe if the replacement check is at least
as strict as the one it replaces.

The old check was "neither TaskSpec path exists".  The new check is the pair
"the canonical registration path is still absent AND the standalone spec's own
gates are still shut".  These tests drive the real ``audit_disposition`` against
the real tree with a mutated copy of the authority record, and prove that each
way the correction could have become an exemption still raises.

The audit cannot be run against a mirrored root: its frozen receipts pin source
members by absolute path, and it refuses both a source member and an authority
outside its root.  So the gate files themselves are never mutated here -- writing
to the live tree is not acceptable.  They are guarded in two layers instead: a recorded digest that
catches any edit, and content assertions inside the audit that catch an edit
made together with a re-pin.  ``test_the_live_gates_are_still_shut`` asserts the
gate contents directly, so opening a gate fails here even before the audit runs.
"""
import json
import sys
import tomllib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.audit_observed_multiome_integration_disposition import (  # noqa: E402
    ObservedMultiomeDispositionError,
    audit_disposition,
)


ROOT = Path(__file__).resolve().parents[2]
RECORD = (
    ROOT
    / "config/artifacts/models/observed_multiome_integration/disposition_20260824.json"
)
CANONICAL = "config/tasks/observed_multiome.toml"
STANDALONE = "config/evaluation/observed_multiome_task.toml"
DEVGATE = "config/evaluation/observed_multiome_development_gate.toml"
SUBJECTS = ("multivi", "scglue", "mofaplus", "seurat_wnn")


class ObservedMultiomeTaskSpecRegroundingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def audit_with(self, mutate=None) -> dict:
        """Run the real audit against the real root, optionally mutating the record."""
        if mutate is None:
            return audit_disposition(ROOT, RECORD)
        body = json.loads(RECORD.read_text(encoding="utf-8"))
        mutate(body)
        # The audit refuses an authority outside its own root, so the mutated
        # copy has to live under ROOT.  It is created and removed per call and
        # nothing in the tree is modified.
        with TemporaryDirectory(dir=ROOT) as tmp:
            path = Path(tmp) / "disposition.json"
            path.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
            return audit_disposition(ROOT, path)

    def assert_refuses(self, mutate, fragment: str | None = None) -> None:
        with self.assertRaises(ObservedMultiomeDispositionError) as caught:
            self.audit_with(mutate)
        if fragment is not None:
            self.assertIn(fragment, str(caught.exception))

    # -- the corrected audit passes because it was re-derived ---------------

    def test_the_corrected_audit_passes_on_the_live_tree(self) -> None:
        receipt = self.audit_with()
        self.assertEqual(
            receipt["status"], "pass_metadata_audit_all_model_lanes_blocked"
        )
        self.assertEqual(receipt["verified_metadata_files"], 34)
        self.assertEqual(receipt["admitted_model_count"], 0)
        self.assertFalse(receipt["canonical_task_spec_exists"])
        self.assertTrue(receipt["standalone_task_spec_exists"])
        self.assertFalse(receipt["standalone_task_spec_execution_authorized"])
        self.assertFalse(receipt["external_seal_bound"])
        self.assertFalse(receipt["scoring_authorized"])

    # -- the live gates are shut, asserted directly ------------------------

    def test_the_live_gates_are_still_shut(self) -> None:
        """If anyone opens a gate, this fails before the audit is even consulted."""
        with (ROOT / DEVGATE).open("rb") as handle:
            gate = tomllib.load(handle)
        self.assertIs(gate["global_census_revision_required_before_execution"], True)
        for key in (
            "model_promotion_allowed",
            "champion_claim_allowed",
            "selection_lock_allowed",
            "external_seal_bound",
            "rna_conditioned_atac_claim_allowed",
        ):
            self.assertIs(gate[key], False, key)
        with (ROOT / STANDALONE).open("rb") as handle:
            spec = tomllib.load(handle)
        self.assertEqual(spec["task_id"], "observed_multiome")
        self.assertEqual(spec["datasets_sealed"], [])
        self.assertIn(
            "execution remains blocked", " ".join(spec["admission_gates"])
        )

    # -- each way the correction could become an exemption still raises ----

    def test_claiming_the_canonical_path_exists_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["canonical_task_spec_exists"] = True
        self.assert_refuses(mutate, "task gate opened")

    def test_denying_the_standalone_spec_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["standalone_task_spec_exists"] = False
        self.assert_refuses(mutate, "task gate opened")

    def test_a_wrong_path_map_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["task_spec_paths"] = {CANONICAL: False}
        self.assert_refuses(mutate, "TaskSpec path check differs")

    def test_inverting_the_path_map_is_refused(self) -> None:
        """The old world -- both absent -- must no longer validate."""
        def mutate(body):
            body["task_gate"]["task_spec_paths"] = {CANONICAL: False, STANDALONE: False}
        self.assert_refuses(mutate)

    def test_a_repinned_standalone_spec_digest_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["standalone_task_spec"]["sha256"] = "0" * 64
        self.assert_refuses(mutate, "standalone TaskSpec binding differs")

    def test_a_repinned_development_gate_digest_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["standalone_development_gate"]["sha256"] = "0" * 64
        self.assert_refuses(mutate, "development gate binding differs")

    def test_claiming_the_spec_is_registered_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["standalone_task_spec"]["registered_in_config_tasks"] = True
        self.assert_refuses(mutate, "standalone TaskSpec binding differs")

    def test_claiming_the_spec_has_a_sealed_dataset_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["standalone_task_spec"]["datasets_sealed"] = ["gse289173"]
        self.assert_refuses(mutate, "standalone TaskSpec binding differs")

    def test_misreporting_the_spec_document_status_is_refused(self) -> None:
        def mutate(body):
            body["task_gate"]["standalone_task_spec"]["document_status"] = "active"
        self.assert_refuses(mutate, "standalone TaskSpec drifted")

    def test_any_model_becoming_admitted_is_refused(self) -> None:
        for subject in SUBJECTS:
            with self.subTest(subject=subject):
                def mutate(body, subject=subject):
                    body["model_dispositions"][subject]["admitted"] = True
                self.assert_refuses(mutate, "model admission gate opened")

    def test_reverting_to_the_v1_schema_or_status_is_refused(self) -> None:
        def revert_schema(body):
            body["schema_version"] = (
                "masld-bench-observed-multiome-integration-disposition-v1"
            )
        self.assert_refuses(revert_schema, "disposition schema differs")

        def revert_status(body):
            body["status"] = "FAIL_CLOSED_NO_TASKSPEC_OR_ADMITTED_MODEL_LANE"
        self.assert_refuses(revert_status, "disposition is not fail-closed")

    def test_reopening_an_execution_authorization_is_refused(self) -> None:
        def mutate(body):
            body["execution_disposition"]["training_authorized"] = True
        self.assert_refuses(mutate, "execution gate opened")

    # -- the correction did not weaken anything ----------------------------

    def test_no_disposition_moved(self) -> None:
        for subject in SUBJECTS:
            self.assertIs(
                self.record["model_dispositions"][subject]["admitted"], False, subject
            )

    def test_every_gate_kept_its_conjunct_count(self) -> None:
        """NO_TASKSPEC became TASKSPEC_EXECUTION_BLOCKED; nothing else moved."""
        for subject in SUBJECTS:
            entry = self.record["model_dispositions"][subject]
            before = entry["gate_superseded_20260825"]
            after = entry["gate"]
            self.assertTrue(before.startswith("BLOCKED_NO_TASKSPEC_"), subject)
            self.assertTrue(
                after.startswith("BLOCKED_TASKSPEC_EXECUTION_BLOCKED_"), subject
            )
            self.assertEqual(
                before[len("BLOCKED_NO_TASKSPEC_"):],
                after[len("BLOCKED_TASKSPEC_EXECUTION_BLOCKED_"):],
                f"{subject}: the surviving conjuncts must be untouched",
            )

    def test_the_two_status_fields_describe_different_axes(self) -> None:
        """There was no conflict, so nothing was made to yield."""
        spec_statuses = set()
        for path in sorted((ROOT / "config" / "tasks").glob("*.toml")):
            with path.open("rb") as handle:
                spec_statuses.add(tomllib.load(handle)["status"])
        with (ROOT / STANDALONE).open("rb") as handle:
            spec_statuses.add(tomllib.load(handle)["status"])
        with (ROOT / "config/evaluation/family_native_tournament.toml").open("rb") as h:
            tournament = tomllib.load(h)
        task_statuses = {v["status"] for v in tournament["task"].values()}

        self.assertEqual(spec_statuses & task_statuses, set())
        self.assertIn("candidate", spec_statuses)
        self.assertIn("planned_registration", task_statuses)
        for task_id in ("rna_conditioned_atac", "variant_to_regulation"):
            with (ROOT / "config" / "tasks" / f"{task_id}.toml").open("rb") as handle:
                self.assertEqual(tomllib.load(handle)["status"], "candidate", task_id)
            self.assertEqual(tournament["task"][task_id]["status"], "active", task_id)
        self.assertEqual(
            self.record["task_gate"]["registry_status_observed"], "planned_registration"
        )

    def test_the_two_path_reality_is_recorded_as_two_values(self) -> None:
        gate = self.record["task_gate"]
        self.assertIs(gate["canonical_task_spec_exists"], False)
        self.assertIs(gate["standalone_task_spec_exists"], True)
        self.assertEqual(gate["task_spec_paths"], {CANONICAL: False, STANDALONE: True})
        self.assertFalse((ROOT / CANONICAL).exists())
        self.assertTrue((ROOT / STANDALONE).is_file())
        self.assertNotIn("task_spec_exists", gate)

    def test_half_the_subjects_are_absent_from_every_authority(self) -> None:
        """Why this is a correction and not a reassessment."""
        for subject in ("scglue", "mofaplus"):
            presence = self.record["model_dispositions"][subject][
                "observed_multiome_authority_presence"
            ]
            for key in (
                "observed_multiome_task_spec",
                "capabilities",
                "family_native_tournament",
                "additive_census_model_check_288",
            ):
                self.assertEqual(presence[key], "absent", f"{subject}:{key}")
        spec_baselines = tomllib.loads(
            (ROOT / STANDALONE).read_text(encoding="utf-8")
        )["baseline_model_ids"]
        self.assertIn("seurat_wnn", spec_baselines)
        self.assertNotIn("multivi", spec_baselines)


if __name__ == "__main__":
    unittest.main()
