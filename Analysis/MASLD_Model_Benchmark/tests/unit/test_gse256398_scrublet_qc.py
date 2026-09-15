from __future__ import annotations

import unittest

import numpy as np
from scipy import sparse

from scripts.build_gse256398_qc_membership import (
    GSE256398MembershipError,
    source_identity,
)
from scripts.probe_gse256398_scrublet_qc import basic_qc


class GSE256398ScrubletQCTests(unittest.TestCase):
    def test_source_identity_is_donor_scoped(self) -> None:
        gsm, sample = source_identity(
            __import__("pathlib").Path(
                "GSM8097071_S1_CB_raw_feature_bc_matrix_filtered.h5"
            )
        )
        self.assertEqual((gsm, sample), ("GSM8097071", "S1"))
        with self.assertRaises(GSE256398MembershipError):
            source_identity(__import__("pathlib").Path("unknown.h5"))

    def test_basic_qc_uses_strict_author_thresholds(self) -> None:
        names = [f"G{i}" for i in range(201)] + ["MT-X"]
        rows = [
            np.ones(202, dtype=np.int32),
            np.concatenate([np.ones(201, dtype=np.int32), np.array([100])]),
            np.full(202, 200, dtype=np.int32),
        ]
        metrics = basic_qc(sparse.csr_matrix(np.stack(rows)), names)
        self.assertEqual(metrics["basic_qc_retained"].tolist(), [True, False, False])


if __name__ == "__main__":
    unittest.main()
