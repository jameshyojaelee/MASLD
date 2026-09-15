from __future__ import annotations

import unittest
from pathlib import Path
import tempfile
from unittest import mock

from scripts.freeze_gse105127_reference_bundle import GSE105127ReferenceError
from scripts.freeze_gse105127_reference_bundle_legacy_exact import (
    SOURCE_BYTES,
    SOURCE_BYTES_UNCOMPRESSED,
    SOURCE_MD5_UNCOMPRESSED,
    SOURCE_SHA256_COMPRESSED,
    SOURCE_SHA256_UNCOMPRESSED,
    SOURCE_CONTIGS_TOTAL,
    LEGACY_ALPHABET_COUNTS,
    fasta_inventory_legacy_exact,
    validate_alphabet_diagnostic,
    validate_legacy_diagnostic,
    validate_source_cache,
    validate_source_inventory_names,
)


class GSE105127LegacyReferenceTests(unittest.TestCase):
    def test_diagnostic_admission_is_exact_and_fail_closed(self) -> None:
        receipt = {
            "schema_version": "masld-bench-gse105127-legacy-reference-stream-diagnostic-v1",
            "status": "pass_exact_legacy_stream",
            "compressed_bytes": SOURCE_BYTES,
            "compressed_sha256": SOURCE_SHA256_COMPRESSED,
            "gzip_exit_status": 2,
            "legacy_warning": "trailing garbage ignored",
            "uncompressed_bytes": SOURCE_BYTES_UNCOMPRESSED,
            "uncompressed_md5": SOURCE_MD5_UNCOMPRESSED,
            "uncompressed_sha256": SOURCE_SHA256_UNCOMPRESSED,
            "fatal_gzip_diagnostic": False,
            "fit_or_score_performed": False,
            "labels_accessed": False,
        }
        validate_legacy_diagnostic(receipt)
        for changed in (
            dict(receipt, gzip_exit_status=0),
            dict(receipt, fatal_gzip_diagnostic=True),
            dict(receipt, extra_field=True),
        ):
            with self.assertRaisesRegex(GSE105127ReferenceError, "diagnostic differs"):
                validate_legacy_diagnostic(changed)

    def test_alphabet_diagnostic_admission_is_exact(self) -> None:
        receipt = {
            "schema_version": "masld-bench-gse105127-reference-alphabet-diagnostic-v1",
            "status": "passed_outcome_free_census",
            "source_bytes": SOURCE_BYTES_UNCOMPRESSED,
            "source_sha256": SOURCE_SHA256_UNCOMPRESSED,
            "header_rows": 84,
            "sequence_rows": 51_696_747,
            "alphabet_counts": LEGACY_ALPHABET_COUNTS,
            "alphabet": sorted(LEGACY_ALPHABET_COUNTS),
            "labels_accessed": False,
            "fit_or_score_performed": False,
        }
        validate_alphabet_diagnostic(receipt)
        with self.assertRaisesRegex(GSE105127ReferenceError, "alphabet diagnostic differs"):
            validate_alphabet_diagnostic(dict(receipt, alphabet=sorted(LEGACY_ALPHABET_COUNTS) + ["Y"]))

    def test_inventory_rejects_any_unobserved_symbol(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.fa"
            path.write_bytes(b">chr1\nACGMRNTY\n")
            with self.assertRaisesRegex(GSE105127ReferenceError, "alphabet differs"):
                fasta_inventory_legacy_exact(path)

    def test_source_inventory_accepts_primary_prefix_and_exact_total(self) -> None:
        primary = [str(value) for value in range(1, 23)] + ["X", "Y", "MT"]
        inventory = [
            {"contig": name, "length": 1}
            for name in primary
            + [f"GL{index:06d}.1" for index in range(SOURCE_CONTIGS_TOTAL - len(primary))]
        ]
        self.assertEqual(validate_source_inventory_names(inventory), [row["contig"] for row in inventory])
        with self.assertRaisesRegex(GSE105127ReferenceError, "contig count differs"):
            validate_source_inventory_names(inventory[:-1])
        swapped = list(inventory)
        swapped[0], swapped[1] = swapped[1], swapped[0]
        with self.assertRaisesRegex(GSE105127ReferenceError, "primary contig prefix differs"):
            validate_source_inventory_names(swapped)

    def test_source_cache_manifest_is_required(self) -> None:
        with tempfile.TemporaryDirectory() as directory, mock.patch(
            "scripts.freeze_gse105127_reference_bundle_legacy_exact.verify_frozen_tree"
        ), mock.patch(
            "scripts.freeze_gse105127_reference_bundle_legacy_exact.digest_file",
            return_value="not-the-bound-manifest",
        ):
            with self.assertRaisesRegex(GSE105127ReferenceError, "source cache ARTIFACTS"):
                validate_source_cache(Path(directory))


if __name__ == "__main__":
    unittest.main()
