from __future__ import annotations

import json
from pathlib import Path
import tomllib
import unittest

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.contracts import ModelManifest


ROOT = Path(__file__).resolve().parents[2]
EXECUTIONS = ROOT / "executions"


class DNALanguageRuntimeStatusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (ROOT / "config/models/dna_language.toml").open("rb") as handle:
            cls.family = tomllib.load(handle)
        cls.models = {item["model_id"]: item for item in cls.family["models"]}

    def test_family_still_loads_as_strict_model_manifests(self) -> None:
        manifests = tuple(
            ModelManifest.from_family_entry(self.family, entry)
            for entry in self.family["models"]
        )
        self.assertEqual(
            {model.model_id for model in manifests},
            {"dnabert2", "nucleotide_transformer", "hyenadna", "caduceus", "evo"},
        )

    def test_three_common_runtime_artifacts_are_frozen_and_bound(self) -> None:
        expected = {
            "dnabert2": (
                "dnabert2-common-runtime-probe-21065134",
                "a0324f4b774c0af4d8bfbdc9bb4ce4514d62e704e4f0ae34bf8cd8d2c67f9545",
            ),
            "nucleotide_transformer": (
                "nucleotide-transformer-common-runtime-probe-21065102",
                "f0c65ebdb5880b9392dfd131e43b3f320df88cd8eccb7accf1f7cc9dd358e325",
            ),
            "hyenadna": (
                "hyenadna-common-runtime-probe-21064769",
                "9deb82d5faae1636cb6129ce6d30b70b8264c6cc30313caa8fb7506293aa12b6",
            ),
        }
        for model_id, (directory, digest) in expected.items():
            with self.subTest(model_id=model_id):
                root = EXECUTIONS / directory
                verify_frozen_tree(root)
                import hashlib

                observed = hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
                self.assertEqual(observed, digest)
                self.assertIn(digest, " ".join(self.models[model_id]["blockers"]))

    def test_nucleotide_transformer_counts_remain_distinctly_labeled(self) -> None:
        root = EXECUTIONS / "nucleotide-transformer-common-runtime-probe-21065102"
        receipt = json.loads(
            (root / "predictions/runtime_probe_receipt.json").read_text(encoding="utf-8")
        )
        self.assertEqual(receipt["model_parameter_count"], 485_699_306)
        blockers = " ".join(self.models["nucleotide_transformer"]["blockers"])
        self.assertIn("485,699,306 executable learned parameters", blockers)
        self.assertIn("485,700,308 tensor elements", blockers)
        self.assertIn("485,729,545", blockers)
        self.assertIn("30,239", blockers)

    def test_runtime_fixtures_do_not_claim_task_head_or_outcome_use(self) -> None:
        for directory in (
            "dnabert2-common-runtime-probe-21065134",
            "nucleotide-transformer-common-runtime-probe-21065102",
            "hyenadna-common-runtime-probe-21064769",
        ):
            receipt = json.loads(
                (
                    EXECUTIONS
                    / directory
                    / "predictions/runtime_probe_receipt.json"
                ).read_text(encoding="utf-8")
            )
            self.assertFalse(receipt["head_fit"])
            self.assertFalse(receipt["observed_outcomes_loaded"])
            self.assertFalse(receipt["sealed_outcomes_loaded"])


if __name__ == "__main__":
    unittest.main()
