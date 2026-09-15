from __future__ import annotations

import gzip
import json
from pathlib import Path
import tempfile
import unittest

from scripts.audit_microarray_transfer_platforms import (
    MicroarrayPlatformAuditError,
    audit_overlap,
    parse_platform_soft,
)


class MicroarrayTransferPlatformAuditTests(unittest.TestCase):
    def test_platform_parser_retains_exact_tabular_axis(self) -> None:
        value = (
            "^PLATFORM = GPL1\n"
            "!Platform_title = test\n"
            "!Platform_data_row_count = 2\n"
            "!platform_table_begin\n"
            "ID\tGB_ACC\n"
            "a\tNM_1\n"
            "b\tNM_2\n"
            "!platform_table_end\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "family.soft.gz"
            path.write_bytes(gzip.compress(value.encode("utf-8")))
            metadata, header, rows = parse_platform_soft(path, "GPL1")
        self.assertEqual(metadata["!Platform_title"], ["test"])
        self.assertEqual(header, ("ID", "GB_ACC"))
        self.assertEqual([row["ID"] for row in rows], ["a", "b"])

    def test_platform_parser_rejects_duplicate_columns(self) -> None:
        value = (
            "^PLATFORM = GPL1\n"
            "!platform_table_begin\n"
            "ID\tID\n"
            "a\tb\n"
            "!platform_table_end\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "family.soft.gz"
            path.write_bytes(gzip.compress(value.encode("utf-8")))
            with self.assertRaisesRegex(MicroarrayPlatformAuditError, "header"):
                parse_platform_soft(path, "GPL1")

    def test_overlap_requires_gse31803_parent_alias(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "source_audit.json").write_text(
                json.dumps(
                    {
                        "cohorts": {
                            "GSE49541": {"relations": ["BioProject: PRJNA214553"]},
                            "GSE83452": {"relations": []},
                        }
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MicroarrayPlatformAuditError, "parent"):
                audit_overlap(root)

    def test_overlap_collapses_antwerp_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "source_audit.json").write_text(
                json.dumps(
                    {
                        "cohorts": {
                            "GSE49541": {
                                "relations": ["SubSeries of: GSE31803"]
                            },
                            "GSE83452": {"relations": []},
                        }
                    }
                ),
                encoding="utf-8",
            )
            result = audit_overlap(root)
        self.assertEqual(
            result["GSE83452"]["cohort_family_id"],
            "antwerp_inserm_shared",
        )
        self.assertEqual(
            set(result["GSE83452"]["accession_aliases"]),
            {"GSE106737", "GSE83452"},
        )


if __name__ == "__main__":
    unittest.main()
