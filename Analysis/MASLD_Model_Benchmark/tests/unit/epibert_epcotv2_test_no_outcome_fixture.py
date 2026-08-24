from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE = Path(__file__).parents[2] / "scripts" / "epibert_epcotv2_no_outcome_fixture.py"
SPEC = importlib.util.spec_from_file_location("observed_atac_fixture_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)


class ObservedATACNoOutcomeFixtureTests(unittest.TestCase):
    def test_epibert_masks_scored_atac_and_is_not_rna_conditioned(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            receipt = fixture.build("epibert", Path(directory) / "epibert", "numpy")
            self.assertFalse(receipt["outcomes_read"])
            self.assertTrue(
                receipt["query_semantics"]["scored_atac_span_removed_from_input"]
            )
            self.assertFalse(receipt["query_semantics"]["rna_conditioning"])
            self.assertFalse(receipt["rna_conditioned_atac_eligible"])

    def test_epcotv2_uses_observed_atac_but_has_no_atac_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            receipt = fixture.build("epcotv2", Path(directory) / "epcotv2", "numpy")
            self.assertFalse(receipt["outcomes_read"])
            self.assertTrue(receipt["query_semantics"]["observed_atac_required_at_query"])
            self.assertFalse(receipt["query_semantics"]["atac_is_output_head"])
            self.assertFalse(receipt["query_semantics"]["rna_conditioning"])

    def test_outcome_like_keys_fail_closed(self) -> None:
        with self.assertRaisesRegex(fixture.NoOutcomeFixtureError, "outcome-like"):
            fixture._validate_keys(["sequence", "observed_y"])


if __name__ == "__main__":
    unittest.main()
