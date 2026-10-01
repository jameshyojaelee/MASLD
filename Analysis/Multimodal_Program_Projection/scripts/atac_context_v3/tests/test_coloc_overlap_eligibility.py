#!/usr/bin/env python3
"""Runs the synthetic coloc.susie shared-posterior eligibility fixture in R."""

from __future__ import annotations

import os
import shutil
import subprocess
import unittest
from pathlib import Path

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parents[4]
HELPER = ROOT / "GWAS/finemapping/src/coloc_signal_eligibility.R"
FIXTURE = TESTS / "smoke_coloc_overlap_eligibility.R"
RSCRIPT = os.environ.get("COLOC_RSCRIPT", "Rscript")


def coloc_available() -> bool:
    if shutil.which(RSCRIPT) is None:
        return False
    probe = subprocess.run(
        [RSCRIPT, "-e", "quit(status = !requireNamespace('coloc', quietly = TRUE))"],
        capture_output=True,
        check=False,
    )
    return probe.returncode == 0


@unittest.skipUnless(coloc_available(), "Rscript with coloc not on PATH (use the rnaseq env)")
class ColocOverlapEligibilityTests(unittest.TestCase):
    def test_fixture(self) -> None:
        result = subprocess.run(
            [RSCRIPT, str(FIXTURE), str(HELPER)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("COLOC_OVERLAP_ELIGIBILITY_FIXTURE\tPASS", result.stdout)


if __name__ == "__main__":
    unittest.main()
