"""Scoping the audit root by species must not open a laundering route.

config/artifacts/models/ holds two species: census bundles, and cross-cutting
task-scoped campaign disposition records whose subject models each already own
a census bundle.  The second kind cannot be registered as a group binding --
binding their subjects re-resolves them away from their own directories and
orphans those instead, taking the audit root from 104 bundles / 2 orphans to
89 / 17 -- so they are declared explicitly instead.

A declaration is a hole unless something guards it.  These tests exercise the
real ``assert_non_census_record_dir`` from the contract suite against synthetic
directories, and prove that each way of abusing the declaration is flagged
rather than skipped.  Same discipline as the campaign['selection'] scoping fix.
"""
import json
import tomllib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from tests.contract.test_registry_files import (
    assert_non_census_record_dir,
    audit_root_species,
)


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY_PATH = ROOT / "config" / "evaluation" / "model_audit_bindings.toml"
AUDIT_ROOT = ROOT / "config" / "artifacts" / "models"
REQUIRED = {"checkpoints.json", "development_crosswalk.json", "exposure_audit.json"}


def load_authority() -> dict:
    with AUTHORITY_PATH.open("rb") as handle:
        return tomllib.load(handle)


def registry_models() -> dict[str, dict]:
    models: dict[str, dict] = {}
    for path in sorted((ROOT / "config" / "models").glob("*.toml")):
        with path.open("rb") as handle:
            for model in tomllib.load(handle).get("models", []):
                models[model["model_id"]] = model
    return models


class NonCensusRecordDirectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authority = load_authority()
        cls.models = registry_models()
        grouped = {
            model_id: binding["bundle_id"]
            for binding in cls.authority["group_bindings"]
            for model_id in binding["model_ids"]
        }
        cls.resolved = {m: grouped.get(m, m) for m in cls.models}
        cls.bundle_ids = set(cls.resolved.values())
        cls.declared, cls.rules = audit_root_species(cls.authority)

    def check(self, bundle: Path, *, audit_root: Path | None = None) -> None:
        assert_non_census_record_dir(
            self,
            bundle,
            required=REQUIRED,
            rules=self.rules,
            resolved=self.resolved,
            bundle_ids=self.bundle_ids,
            audit_root=audit_root if audit_root is not None else AUDIT_ROOT,
        )

    # -- the declaration is what we think it is ----------------------------

    def test_declaration_names_exactly_the_two_known_record_directories(self) -> None:
        self.assertEqual(
            sorted(self.declared),
            ["observed_multiome_integration", "remaining_multimodal_families"],
        )
        for name in self.declared:
            self.assertTrue((AUDIT_ROOT / name).is_dir(), name)

    def test_every_live_declared_directory_satisfies_the_contract(self) -> None:
        for name in sorted(self.declared):
            self.check(AUDIT_ROOT / name)

    def test_the_aggregate_rule_has_no_exceptions(self) -> None:
        """A group binding exists because its members own no directory.

        This is what separates a real aggregate family from a cross-cutting
        record, so it is the property the declaration must never blur.
        """
        for binding in self.authority["group_bindings"]:
            for model_id in binding["model_ids"]:
                self.assertFalse(
                    (AUDIT_ROOT / model_id).is_dir(),
                    f"{binding['bundle_id']}:{model_id}",
                )

    def test_declared_directories_never_shadow_a_census_bundle(self) -> None:
        self.assertEqual(self.declared & self.bundle_ids, set())

    # -- each abuse of the declaration is caught ---------------------------

    def _record_dir(self, tmp: str, name: str, payload: dict) -> Path:
        bundle = Path(tmp) / name
        bundle.mkdir()
        (bundle / "disposition_20260824.json").write_text(
            json.dumps(payload), encoding="utf-8"
        )
        return bundle

    def test_a_record_dir_that_grows_a_member_file_is_rejected(self) -> None:
        """A bundle pretending to be a record must go back into the census."""
        with TemporaryDirectory() as tmp:
            bundle = self._record_dir(
                tmp, "pretend_record", {"model_dispositions": {"stabmap": {}}}
            )
            self.check(bundle)  # clean so far
            (bundle / "checkpoints.json").write_text("{}", encoding="utf-8")
            with self.assertRaises(AssertionError):
                self.check(bundle)

    def test_a_record_dir_with_no_disposition_record_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            bundle = Path(tmp) / "empty_record"
            bundle.mkdir()
            (bundle / "notes.json").write_text(
                json.dumps({"model_dispositions": {"stabmap": {}}}), encoding="utf-8"
            )
            with self.assertRaises(AssertionError):
                self.check(bundle)

    def test_a_record_naming_a_grouped_model_is_rejected(self) -> None:
        """This is the laundering route: hiding a real family as a record.

        A grouped model_id has no directory of its own, so a record naming one
        is describing a family, not cross-cutting work over existing bundles.
        """
        grouped_model = next(
            model_id
            for model_id, bundle_id in self.resolved.items()
            if bundle_id != model_id
        )
        with TemporaryDirectory() as tmp:
            bundle = self._record_dir(
                tmp, "family_in_disguise", {"model_dispositions": {grouped_model: {}}}
            )
            with self.assertRaises(AssertionError):
                self.check(bundle)

    def test_a_record_naming_no_registry_model_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            bundle = self._record_dir(tmp, "names_nothing", {"note": "unrelated"})
            with self.assertRaises(AssertionError):
                self.check(bundle)

    def test_a_declaration_that_shadows_a_census_bundle_is_rejected(self) -> None:
        census_bundle = sorted(self.bundle_ids)[0]
        with TemporaryDirectory() as tmp:
            bundle = self._record_dir(
                tmp, census_bundle, {"model_dispositions": {"stabmap": {}}}
            )
            with self.assertRaises(AssertionError):
                self.check(bundle)

    def test_a_named_model_without_its_own_directory_is_rejected(self) -> None:
        """Word-boundary naming plus an own-directory requirement, together."""
        with TemporaryDirectory() as tmp:
            fake_root = Path(tmp) / "audit"
            fake_root.mkdir()
            bundle = self._record_dir(
                tmp, "missing_subject_dir", {"model_dispositions": {"stabmap": {}}}
            )
            with self.assertRaises(AssertionError):
                self.check(bundle, audit_root=fake_root)

    # -- an undeclared directory still fails, exactly as before -------------

    def test_an_undeclared_extra_directory_is_still_a_failure(self) -> None:
        observed = {p.name for p in AUDIT_ROOT.iterdir() if p.is_dir()}
        self.assertEqual(observed, self.bundle_ids | self.declared)
        intruder = observed | {"some_new_untracked_directory"}
        self.assertNotEqual(intruder, self.bundle_ids | self.declared)

    def test_removing_the_declaration_makes_the_dirs_fail_again(self) -> None:
        """The scoping narrows nothing that was previously covered."""
        empty_declaration, _ = audit_root_species({})
        self.assertEqual(empty_declaration, set())
        observed = {p.name for p in AUDIT_ROOT.iterdir() if p.is_dir()}
        self.assertNotEqual(observed, self.bundle_ids | empty_declaration)


if __name__ == "__main__":
    unittest.main()
