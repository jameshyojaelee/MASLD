from __future__ import annotations

import unittest

import numpy as np

from scripts.audit_microarray_cross_accession_overlap import (
    _correlations,
    _evidence_digest,
)


class CrossAccessionOverlapAuditTests(unittest.TestCase):
    def test_correlation_identifies_exact_and_reversed_vectors(self) -> None:
        left = np.asarray([[1.0, 4.0], [2.0, 3.0], [3.0, 2.0], [4.0, 1.0]])
        right = np.asarray([[10.0, 1.0], [20.0, 2.0], [30.0, 3.0], [40.0, 4.0]])
        observed = _correlations(left, right)
        self.assertAlmostEqual(observed[0, 0], 1.0)
        self.assertAlmostEqual(observed[1, 0], -1.0)

    def test_nonreversible_evidence_digest_is_domain_bound(self) -> None:
        first = _evidence_digest("a" * 64, "b" * 64, "c" * 64, "GSM1", "GSM2", 0.99999)
        second = _evidence_digest("a" * 64, "b" * 64, "c" * 64, "GSM1", "GSM3", 0.99999)
        self.assertEqual(len(first), 64)
        self.assertNotEqual(first, second)


if __name__ == "__main__":
    unittest.main()
