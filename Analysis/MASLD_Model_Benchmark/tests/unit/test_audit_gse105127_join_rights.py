from __future__ import annotations

import unittest

from scripts.audit_gse105127_join_rights import (
    audit_article,
    fold_mapping,
    parse_group_mapping,
    validate_supplementary_files,
)


class GSE105127JoinRightsTests(unittest.TestCase):
    def test_group_mapping_is_participant_unique(self) -> None:
        design = (
            "normal control (NC = 1, 2, 3, 4), healthy obese "
            "(HO = 5, 6, 7, 8, 9), bland steatosis "
            "(STEA = 10, 11, 12, 13, 14) and early NASH "
            "(EARLY = 15, 16, 17, 18, 19)."
        )
        observed = parse_group_mapping(design)
        self.assertEqual(len(observed), 19)
        self.assertEqual(observed["1"], "NC")
        self.assertEqual(observed["19"], "EARLY")

    def test_stage_balanced_folds_keep_participants_intact(self) -> None:
        mapping = {
            **{str(value): "NC" for value in range(1, 5)},
            **{str(value): "HO" for value in range(5, 10)},
            **{str(value): "STEA" for value in range(10, 15)},
            **{str(value): "EARLY" for value in range(15, 20)},
        }
        folds = fold_mapping(mapping)
        self.assertEqual(sorted(set(folds.values())), [0, 1, 2, 3, 4])
        self.assertEqual(sum(value == 4 for value in folds.values()), 3)

    def test_rrbs_requires_methylation_and_coverage_files(self) -> None:
        files = [
            "ftp://example/sample.CG.bed.gz",
            "ftp://example/sample.CG.bw",
            "ftp://example/sample.CG.ct_coverage.bw",
        ]
        self.assertEqual(
            validate_supplementary_files("RRBS", files), files[0]
        )

    def test_article_proves_doi_and_adjacent_sections(self) -> None:
        xml = b"""<article xmlns:xlink="http://www.w3.org/1999/xlink">
          <front><article-meta>
            <article-id pub-id-type="doi">10.1038/s41467-018-06611-5</article-id>
            <permissions><license><license-ref xlink:href="https://example/license"/></license></permissions>
          </article-meta></front>
          <body><p>Adjacent cryosections were used to generate matching reduced representation bisulfite sequencing data.</p></body>
        </article>"""
        observed = audit_article(xml)
        self.assertTrue(observed["adjacent_section_evidence_present"])
        self.assertEqual(observed["doi"], "10.1038/s41467-018-06611-5")


if __name__ == "__main__":
    unittest.main()
