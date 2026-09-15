from __future__ import annotations

import gzip
import io
from pathlib import Path
import tarfile
import tempfile
import unittest

from scripts.acquire_audit_microarray_transfer_raw import (
    MicroarrayRawAuditError,
    audit_tar,
    identify_cel_format,
)


def write_tar(path: Path, name: str, value: bytes) -> None:
    compressed = gzip.compress(value)
    with tarfile.open(path, mode="w") as archive:
        info = tarfile.TarInfo(name)
        info.size = len(compressed)
        archive.addfile(info, io.BytesIO(compressed))


class MicroarrayTransferRawAuditTests(unittest.TestCase):
    def test_audit_streams_gzip_without_extracting(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "raw.tar"
            write_tar(path, "GSM123.CEL.gz", b"[CEL]\nVersion=3\n")
            values = audit_tar(path, {"GSM123"})
        self.assertEqual(values[0]["decompressed_bytes"], 16)
        self.assertEqual(values[0]["cel_format"], "text_v3")
        self.assertTrue(str(values[0]["decompressed_sha256"]))

    def test_cel_signatures_cover_text_xda_and_calvin(self) -> None:
        self.assertEqual(identify_cel_format(b"[CEL]\nVersion=3\n"), "text_v3")
        self.assertEqual(
            identify_cel_format(b"\x40\x00\x00\x00\x04\x00\x00\x00"),
            "xda_v4",
        )
        self.assertEqual(identify_cel_format(b"\x3b\x01\x00\x00\x00\x01"), "calvin_v1")

    def test_cel_signature_rejects_arbitrary_nonempty_payload(self) -> None:
        with self.assertRaisesRegex(MicroarrayRawAuditError, "signature"):
            identify_cel_format(b"not a CEL file")

    def test_audit_rejects_sample_axis_difference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "raw.tar"
            write_tar(path, "GSM123.CEL.gz", b"[CEL]\nVersion=3\n")
            with self.assertRaisesRegex(MicroarrayRawAuditError, "axes"):
                audit_tar(path, {"GSM999"})


if __name__ == "__main__":
    unittest.main()
