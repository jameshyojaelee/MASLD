from __future__ import annotations

import csv
import gzip
from hashlib import sha256
import io
from pathlib import Path
import tarfile
import tempfile
import unittest

from scripts.stage_microarray_cel_arrays import (
    CELStagingError,
    read_manifest,
    stage_arrays,
)

COLUMNS = [
    "series",
    "sample_accession",
    "member_name",
    "compressed_bytes",
    "decompressed_bytes",
    "decompressed_sha256",
    "header_hex_32",
]


def payload(accession: str) -> bytes:
    return f"CEL-BODY-{accession}".encode("ascii") * 8


def build_source(
    base: Path, *, series: str = "GSE83452", records: int = 3, other_series: int = 0
) -> tuple[Path, Path, list[str]]:
    accessions = [f"GSM{index:06d}" for index in range(records)]
    tar_path = base / "RAW.tar"
    rows = []
    with tarfile.open(tar_path, "w") as archive:
        for accession in accessions:
            body = payload(accession)
            compressed = gzip.compress(body, mtime=0)
            info = tarfile.TarInfo(f"{accession}_01_HuGene-2_0-st_.CEL.gz")
            info.size = len(compressed)
            archive.addfile(info, io.BytesIO(compressed))
            rows.append(
                {
                    "series": series,
                    "sample_accession": accession,
                    "member_name": info.name,
                    "compressed_bytes": str(len(compressed)),
                    "decompressed_bytes": str(len(body)),
                    "decompressed_sha256": sha256(body).hexdigest(),
                    "header_hex_32": "00",
                }
            )
    # Rows from a different series in the same manifest must be ignored.
    for index in range(other_series):
        rows.append(
            {
                "series": "GSE49541",
                "sample_accession": f"GSMOTHER{index}",
                "member_name": f"other{index}.CEL.gz",
                "compressed_bytes": "1",
                "decompressed_bytes": "1",
                "decompressed_sha256": "0" * 64,
                "header_hex_32": "00",
            }
        )
    manifest_path = base / "manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return tar_path, manifest_path, accessions


class ManifestTests(unittest.TestCase):
    def test_foreign_series_rows_are_filtered_out(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            _, manifest_path, accessions = build_source(base, records=3, other_series=5)
            rows = read_manifest(manifest_path, series="GSE83452", expected_records=3)
            self.assertEqual([row["sample_accession"] for row in rows], sorted(accessions))

    def test_wrong_record_count_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            _, manifest_path, _ = build_source(base, records=3)
            with self.assertRaisesRegex(CELStagingError, "not 231 records"):
                read_manifest(manifest_path, series="GSE83452", expected_records=231)

    def test_unknown_series_yields_no_records(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            _, manifest_path, _ = build_source(base, records=3)
            with self.assertRaisesRegex(CELStagingError, "not 3 records"):
                read_manifest(manifest_path, series="GSE00000", expected_records=3)


class StageArraysTests(unittest.TestCase):
    def test_stages_every_array_and_verifies_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, accessions = build_source(base, records=3, other_series=2)
            output = base / "cel"
            receipt = stage_arrays(
                tar_path=tar_path,
                manifest_path=manifest_path,
                output=output,
                series="GSE83452",
                expected_records=3,
                platform_id="GPL16686",
                accessions=None,
            )
            self.assertEqual(receipt["staged_arrays"], 3)
            self.assertEqual(receipt["platform_id"], "GPL16686")
            self.assertIs(receipt["every_array_matches_frozen_sha256"], True)
            self.assertIs(receipt["raw_CEL_written_into_published_artifact"], False)
            for accession in accessions:
                self.assertEqual((output / f"{accession}.CEL").read_bytes(), payload(accession))

    def test_subset_staging_writes_only_requested_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, accessions = build_source(base, records=3)
            output = base / "cel"
            receipt = stage_arrays(
                tar_path=tar_path,
                manifest_path=manifest_path,
                output=output,
                series="GSE83452",
                expected_records=3,
                platform_id="GPL16686",
                accessions=[accessions[2]],
            )
            self.assertEqual(receipt["staged_arrays"], 1)
            self.assertEqual(
                [path.name for path in output.iterdir()], [f"{accessions[2]}.CEL"]
            )

    def test_corrupted_payload_fails_the_frozen_hash(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, accessions = build_source(base, records=3)
            rows = read_manifest(manifest_path, series="GSE83452", expected_records=3)
            rows[0]["decompressed_sha256"] = "0" * 64
            with manifest_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(CELStagingError, "SHA-256 differs"):
                stage_arrays(
                    tar_path=tar_path,
                    manifest_path=manifest_path,
                    output=base / "cel",
                    series="GSE83452",
                    expected_records=3,
                    platform_id="GPL16686",
                    accessions=[accessions[0]],
                )

    def test_unknown_accession_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, _ = build_source(base, records=3)
            with self.assertRaisesRegex(CELStagingError, "not in the frozen manifest"):
                stage_arrays(
                    tar_path=tar_path,
                    manifest_path=manifest_path,
                    output=base / "cel",
                    series="GSE83452",
                    expected_records=3,
                    platform_id="GPL16686",
                    accessions=["GSM999999"],
                )


if __name__ == "__main__":
    unittest.main()
