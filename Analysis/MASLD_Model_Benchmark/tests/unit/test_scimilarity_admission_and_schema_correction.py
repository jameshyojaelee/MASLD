"""Two scimilarity corrections, and the invariants that keep them honest.

Item 1 was sanctioned schema evolution recorded as drift: scimilarity is the
only bundle on -v2, and all three of its v2 files are strict supersets of the
v1 shape.  The contract test's accepted sets were widened rather than the
bundle migrated back, because migrating back would delete real structure to
satisfy a frozen list.

Item 2 was a test asserting a state the registry refuses.  Two lanes edited
scimilarity the same day in opposite directions -- assertions were added
requiring an admitted model while the registry was rewritten to keep it
blocked -- and the registry won because ``admission_blocking == bool(blockers)``
holds for every model in the census.

Both corrections widen or relax something, so both need a guard.  These tests
are those guards.
"""
import json
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "config" / "artifacts" / "models"
SCIMILARITY = MODELS / "scimilarity"


def registry_models() -> dict[str, dict]:
    models: dict[str, dict] = {}
    for path in sorted((ROOT / "config" / "models").glob("*.toml")):
        with path.open("rb") as handle:
            for model in tomllib.load(handle).get("models", []):
                models[model["model_id"]] = model
    return models


class AdmissionBlockingInvariantTests(unittest.TestCase):
    """admission_blocking means 'has unresolved blockers that block admission'."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.models = registry_models()

    def test_admission_blocking_equals_having_blockers_for_every_model(self) -> None:
        """The invariant the scimilarity correction was decided on.

        This is the guard: nobody can clear admission_blocking on any model
        while leaving its blockers in place, or vice versa, without failing
        here.  Admitting a model means clearing its blockers first.
        """
        violations = [
            (model_id, model["admission_blocking"], len(model["blockers"]))
            for model_id, model in self.models.items()
            if model["admission_blocking"] != bool(model["blockers"])
        ]
        self.assertEqual(violations, [], "admission_blocking must track blockers")
        self.assertEqual(len(self.models), 136)

    def test_scimilarity_is_admission_blocked_with_its_stated_reason(self) -> None:
        model = self.models["scimilarity_v1_1"]
        self.assertIs(model["admission_blocking"], True)
        self.assertEqual(len(model["blockers"]), 1)
        self.assertIn("Keep this model admission-blocked", model["blockers"][0])

    def test_the_admission_blocker_is_the_contract_not_the_exposure_class(self) -> None:
        """Two different gates; clearing one would not clear the other.

        A reader conflating them would think admitting scimilarity requires
        resolving its exposure class, which is a different and much harder
        problem than materialising a ModelExecutionContract.
        """
        blocker = self.models["scimilarity_v1_1"]["blockers"][0]
        self.assertIn("ModelExecutionContract", blocker)
        audit = json.loads(
            (SCIMILARITY / "exposure_audit.json").read_text(encoding="utf-8")
        )
        # encoder_seen drives champion ineligibility, a separate gate.
        self.assertEqual(
            self.models["scimilarity_v1_1"]["exposure_status"], "encoder_seen"
        )
        self.assertFalse(audit["sealed_champion_eligible"])
        self.assertEqual(
            audit["sealed_champion_eligibility"],
            "ineligible_aggregate_development_encoder_seen",
        )

    def test_scimilarity_is_still_the_only_encoder_seen_model(self) -> None:
        seen = sorted(
            model_id
            for model_id, model in self.models.items()
            if model["exposure_status"] == "encoder_seen"
        )
        self.assertEqual(seen, ["scimilarity_v1_1"])


class VersionedSchemaSupersetTests(unittest.TestCase):
    """A schema bump may add structure; it may never shed a requirement."""

    FAMILIES = {
        "checkpoints.json": "masld-bench-upstream-checkpoint-preflight-v",
        "development_crosswalk.json": "masld-bench-development-crosswalk-v",
        "exposure_audit.json": "masld-bench-checkpoint-exposure-audit-v",
    }

    @staticmethod
    def shape(node, prefix: str = "") -> set[str]:
        paths: set[str] = set()
        if isinstance(node, dict):
            for key, value in node.items():
                paths.add(prefix + key)
                paths |= VersionedSchemaSupersetTests.shape(value, f"{prefix}{key}.")
        elif isinstance(node, list) and node and isinstance(node[0], dict):
            paths |= VersionedSchemaSupersetTests.shape(node[0], f"{prefix}[].")
        return paths

    def test_scimilarity_is_the_only_versioned_successor_bundle(self) -> None:
        successors = []
        for filename, stem in self.FAMILIES.items():
            for bundle in sorted(p for p in MODELS.iterdir() if p.is_dir()):
                member = bundle / filename
                if not member.is_file():
                    continue
                version = json.loads(member.read_text(encoding="utf-8")).get(
                    "schema_version", ""
                )
                if version.startswith(stem) and version != f"{stem}1":
                    successors.append((bundle.name, filename, version))
        self.assertEqual(
            sorted({name for name, _, _ in successors}),
            ["scimilarity"],
            "a new versioned bundle appeared; confirm it is real evolution",
        )
        self.assertEqual(len(successors), 3)

    def test_each_v2_file_is_a_strict_superset_of_the_v1_shape(self) -> None:
        for filename, stem in self.FAMILIES.items():
            v1_shapes = []
            for bundle in sorted(p for p in MODELS.iterdir() if p.is_dir()):
                member = bundle / filename
                if not member.is_file():
                    continue
                body = json.loads(member.read_text(encoding="utf-8"))
                if body.get("schema_version") == f"{stem}1":
                    v1_shapes.append(self.shape(body))
            self.assertTrue(v1_shapes, filename)
            universal = set.intersection(*v1_shapes)
            v2 = self.shape(
                json.loads((SCIMILARITY / filename).read_text(encoding="utf-8"))
            )
            self.assertEqual(sorted(universal - v2), [], f"scimilarity/{filename}")

    def test_the_v2_files_carry_the_structure_that_justified_the_bump(self) -> None:
        """If these disappear, the bump was not evolution after all."""
        checkpoints = json.loads(
            (SCIMILARITY / "checkpoints.json").read_text(encoding="utf-8")
        )
        self.assertIn("staged_model_contract", checkpoints)
        exposure = json.loads(
            (SCIMILARITY / "exposure_audit.json").read_text(encoding="utf-8")
        )
        self.assertIn("frozen_fixture_contract", exposure)
        self.assertIn("sealed_champion_eligible", exposure)
        crosswalk = json.loads(
            (SCIMILARITY / "development_crosswalk.json").read_text(encoding="utf-8")
        )
        for key in ("encoder_training", "reference_only", "sealed_timing"):
            self.assertIn(key, crosswalk["rules"], key)
        for key in ("annotated_training_cells", "annotated_test_cells"):
            self.assertIn(key, crosswalk["corpus"], key)


if __name__ == "__main__":
    unittest.main()
