from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from scripts import build_gse296875_fragment_membership as membership


class GSE296875FragmentMembershipTests(unittest.TestCase):
    def test_fold_index_matches_registered_donor_hash(self) -> None:
        expected = int.from_bytes(
            sha256(b"20260821\x00331").digest()[:8], "big"
        ) % 5
        self.assertEqual(membership.fold_index("331"), expected)
        with self.assertRaises(membership.MembershipError):
            membership.fold_index("")

    def test_parse_fragment_roster_rejects_missing_well(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "filelist.txt"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["#Archive/File", "Name", "Time", "Size", "Type"])
                writer.writerow(
                    [
                        "File",
                        "GSM1_well1_atac_fragments.tsv.gz",
                        "date",
                        "1",
                        "TSV",
                    ]
                )
            with mock.patch.object(
                membership, "sha256_file", return_value=membership.FILELIST_SHA256
            ):
                with self.assertRaises(membership.MembershipError):
                    membership.parse_fragment_roster(path)

    def test_registered_lineage_roles_are_disjoint_and_complete(self) -> None:
        primary = set(membership.PRIMARY_LINEAGES)
        secondary = set(membership.SECONDARY_LINEAGES)
        self.assertFalse(primary & secondary)
        self.assertEqual(
            primary | secondary,
            {value[0] for value in membership.LABEL_CONTRACT.values()},
        )
        self.assertEqual(sum(value[0] for value in membership.EXPECTED_LABEL_CENSUS.values()), 68_398)


if __name__ == "__main__":
    unittest.main()
