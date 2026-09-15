from __future__ import annotations

import unittest

from scripts.diagnose_microarray_cross_accession_timepoints import token_disposition


class CrossAccessionTimepointDiagnosticTests(unittest.TestCase):
    def test_exact_token_is_retained(self) -> None:
        self.assertEqual(
            token_disposition("baseline", "baseline"),
            "exact_same_source_token",
        )

    def test_hyphenation_difference_is_encoding_semantics(self) -> None:
        self.assertEqual(
            token_disposition("followup", "follow-up"),
            "encoding_semantics_punctuation_only",
        )

    def test_different_visit_is_not_coerced(self) -> None:
        self.assertEqual(
            token_disposition("baseline", "follow-up"),
            "semantic_mismatch_review_required",
        )


if __name__ == "__main__":
    unittest.main()
