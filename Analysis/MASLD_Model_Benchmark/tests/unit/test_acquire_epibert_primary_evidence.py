#!/usr/bin/env python3
"""Unit checks for EpiBERT primary-evidence acquisition guards."""

from __future__ import annotations

import importlib.util
from io import BytesIO
from pathlib import Path
import tarfile
import tempfile
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/acquire_epibert_primary_evidence.py"
SPEC = importlib.util.spec_from_file_location("epibert_evidence", MODULE)
assert SPEC and SPEC.loader
evidence = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evidence)


class EpiBERTPrimaryEvidenceTest(unittest.TestCase):
    def test_oa_package_url_accepts_only_exact_ncbi_tgz(self) -> None:
        payload = b'<OA><records><record id="PMC1"><link format="tgz" href="ftp://ftp.ncbi.nlm.nih.gov/pub/a.tar.gz" /></record></records></OA>'
        self.assertEqual(
            evidence.oa_package_url(payload, "PMC1"),
            "ftp://ftp.ncbi.nlm.nih.gov/pub/a.tar.gz",
        )
        with self.assertRaises(evidence.EvidenceAcquisitionError):
            evidence.oa_package_url(payload, "PMC2")

    def test_safe_member_rejects_traversal_and_links(self) -> None:
        with self.assertRaises(evidence.EvidenceAcquisitionError):
            evidence.safe_member_name("../escape")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "bad.tar.gz"
            with tarfile.open(archive, "w:gz") as handle:
                member = tarfile.TarInfo("paper/link")
                member.type = tarfile.SYMTYPE
                member.linkname = "elsewhere"
                handle.addfile(member)
            with self.assertRaises(evidence.EvidenceAcquisitionError):
                evidence.extract_pmc_package(archive, root / "output")

    def test_jats_supplement_hrefs_require_safe_basenames(self) -> None:
        payload = b'''<article xmlns:xlink="http://www.w3.org/1999/xlink"><supplementary-material id="mmc1"><caption><title>Document S1</title></caption><media xlink:href="mmc1.pdf"/></supplementary-material><supplementary-material id="mmc2"><caption><title>Table S1</title></caption><media xlink:href="mmc2.xlsx"/></supplementary-material></article>'''
        self.assertEqual(evidence.jats_supplement_hrefs(payload), ["mmc1.pdf", "mmc2.xlsx"])
        self.assertEqual(
            evidence.jats_supplement_inventory(payload)[0],
            {"id": "mmc1", "title": "Document S1", "href": "mmc1.pdf"},
        )
        unsafe = b'''<article xmlns:xlink="http://www.w3.org/1999/xlink"><supplementary-material><media xlink:href="../mmc1.pdf"/></supplementary-material></article>'''
        with self.assertRaises(evidence.EvidenceAcquisitionError):
            evidence.jats_supplement_hrefs(unsafe)

    def test_regular_pmc_member_is_hashed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "paper.tar.gz"
            with tarfile.open(archive, "w:gz") as handle:
                payload = b"primary source"
                member = tarfile.TarInfo("paper/article.nxml")
                member.size = len(payload)
                handle.addfile(member, BytesIO(payload))
            rows = evidence.extract_pmc_package(archive, root / "output")
            self.assertEqual(rows[0]["path"], "paper/article.nxml")
            self.assertEqual(rows[0]["size_bytes"], len(payload))


if __name__ == "__main__":
    unittest.main()
