from __future__ import annotations

import csv
import gzip
from hashlib import sha256
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
import unittest.mock

from scripts.stage_gse49541_cel_arrays import (
    CELStagingError,
    read_manifest,
    stage_arrays,
)


def payload(accession: str) -> bytes:
    return f"CEL-BODY-{accession}".encode("ascii") * 8


def build_source(base: Path, records: int = 3) -> tuple[Path, Path, list[str]]:
    accessions = [f"GSM{index:06d}" for index in range(records)]
    tar_path = base / "GSE49541_RAW.tar"
    rows = []
    with tarfile.open(tar_path, "w") as archive:
        for accession in accessions:
            body = payload(accession)
            compressed = gzip.compress(body, mtime=0)
            info = tarfile.TarInfo(f"{accession}.CEL.gz")
            info.size = len(compressed)
            archive.addfile(info, io.BytesIO(compressed))
            rows.append(
                {
                    "series": "GSE49541",
                    "sample_accession": accession,
                    "member_name": info.name,
                    "compressed_bytes": str(len(compressed)),
                    "decompressed_bytes": str(len(body)),
                    "decompressed_sha256": sha256(body).hexdigest(),
                    "header_hex_32": "00",
                }
            )
    manifest_path = base / "manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return tar_path, manifest_path, accessions


class StageArraysTests(unittest.TestCase):
    def setUp(self) -> None:
        patcher = unittest.mock.patch(
            "scripts.stage_gse49541_cel_arrays.EXPECTED_RECORDS", 3
        )
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_stages_every_array_and_verifies_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, accessions = build_source(base)
            output = base / "cel"
            receipt = stage_arrays(
                tar_path=tar_path, manifest_path=manifest_path, output=output, accessions=None
            )
            self.assertEqual(receipt["staged_arrays"], 3)
            self.assertIs(receipt["every_array_matches_frozen_sha256"], True)
            self.assertIs(receipt["raw_CEL_written_into_published_artifact"], False)
            for accession in accessions:
                self.assertEqual((output / f"{accession}.CEL").read_bytes(), payload(accession))

    def test_subset_staging_writes_only_requested_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, accessions = build_source(base)
            output = base / "cel"
            receipt = stage_arrays(
                tar_path=tar_path,
                manifest_path=manifest_path,
                output=output,
                accessions=[accessions[2], accessions[0]],
            )
            self.assertEqual(receipt["staged_arrays"], 2)
            self.assertEqual(
                sorted(path.name for path in output.iterdir()),
                [f"{accessions[0]}.CEL", f"{accessions[2]}.CEL"],
            )

    def test_unknown_accession_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, _ = build_source(base)
            with self.assertRaisesRegex(CELStagingError, "not in the frozen manifest"):
                stage_arrays(
                    tar_path=tar_path,
                    manifest_path=manifest_path,
                    output=base / "cel",
                    accessions=["GSM999999"],
                )

    def test_corrupted_payload_fails_the_frozen_hash(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            tar_path, manifest_path, accessions = build_source(base)
            rows = read_manifest(manifest_path)
            rows[0]["decompressed_sha256"] = "0" * 64
            with manifest_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(CELStagingError, "SHA-256 differs"):
                stage_arrays(
                    tar_path=tar_path,
                    manifest_path=manifest_path,
                    output=base / "cel",
                    accessions=[accessions[0]],
                )

    def test_manifest_of_wrong_size_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            _, manifest_path, _ = build_source(base, records=2)
            with self.assertRaisesRegex(CELStagingError, "not 3 records"):
                read_manifest(manifest_path)

    def test_manifest_is_returned_in_sorted_accession_order(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            _, manifest_path, accessions = build_source(base)
            rows = read_manifest(manifest_path)
            self.assertEqual([row["sample_accession"] for row in rows], sorted(accessions))


if __name__ == "__main__":
    unittest.main()
