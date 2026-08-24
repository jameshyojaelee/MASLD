"""Lock the distinction between histopathology and adjudicated MASLD status."""

from __future__ import annotations

from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]


def _load(relative: str) -> dict:
    with (ROOT / relative).open("rb") as handle:
        return tomllib.load(handle)


class GSE296875PhenotypeContractTests(unittest.TestCase):
    def test_histopathology_is_present_but_not_a_masld_case_label(self) -> None:
        dataset = _load("config/datasets/gse296875.toml")
        evaluation = _load("config/evaluation/gse296875_histopathology.toml")

        self.assertIn("pathology_scores", dataset["modalities"])
        self.assertEqual(dataset["expected_biological_units"], 39)
        self.assertEqual(evaluation["steatosis"]["observed_donors"], 38)
        self.assertEqual(evaluation["fibrosis"]["observed_donors"], 37)
        self.assertIs(evaluation["external_or_sealed"], False)
        self.assertIs(evaluation["champion_eligible"], False)
        self.assertEqual(evaluation["role"], "secondary_development_only")
        self.assertIn("MASLD diagnosis", evaluation["claims"]["forbidden"])
        self.assertIn("MASH diagnosis", evaluation["claims"]["forbidden"])
        self.assertIn("standard fibrosis stage", evaluation["claims"]["forbidden"])
        self.assertIs(dataset["admission_blocking"], True)

    def test_missing_pathology_cannot_be_encoded_as_absence(self) -> None:
        evaluation = _load("config/evaluation/gse296875_histopathology.toml")
        self.assertEqual(evaluation["steatosis"]["missing_donors"], 1)
        self.assertEqual(evaluation["fibrosis"]["missing_donors"], 2)
        self.assertEqual(
            evaluation["inference"]["missingness"],
            "explicit_endpoint_mask_never_absence",
        )
        self.assertEqual(evaluation["inference"]["bootstrap_replicates"], 10_000)
        self.assertEqual(evaluation["inference"]["uncertainty"], "paired_donor_cluster_bootstrap")


if __name__ == "__main__":
    unittest.main()
