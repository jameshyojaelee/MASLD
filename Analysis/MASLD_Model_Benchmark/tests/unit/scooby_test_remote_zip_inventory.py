from __future__ import annotations

import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "scooby_remote_zip_inventory.py"
SPEC = importlib.util.spec_from_file_location("scooby_remote_zip_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
remote_zip = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(remote_zip)


def _zip_bytes(names: list[str]) -> bytes:
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for index, name in enumerate(names):
            handle.writestr(name, f"payload-{index}".encode("ascii"))
    return stream.getvalue()


class RemoteZipInventoryTests(unittest.TestCase):
    def test_bounded_inventory_finds_context_candidates(self) -> None:
        payload = _zip_bytes(
            [
                "training/README.txt",
                "training/onek1k_scPoli_reference.pt",
                "training/cell_embeddings.parquet",
            ]
        )
        reader = remote_zip.BytesRangeReader(payload)
        receipt = remote_zip.inventory_reader(reader)
        self.assertEqual(receipt["entry_count"], 3)
        self.assertEqual(receipt["context_candidate_count"], 2)
        self.assertFalse(receipt["archive_payload_downloaded"])
        self.assertLessEqual(reader.request_count, 3)
        self.assertLessEqual(reader.bytes_downloaded, len(payload) * 2)

    def test_unsafe_member_path_fails_closed(self) -> None:
        payload = _zip_bytes(["../escape.pt"])
        with self.assertRaisesRegex(
            remote_zip.RemoteZipInventoryError, "unsafe ZIP member path"
        ):
            remote_zip.inventory_reader(remote_zip.BytesRangeReader(payload))

    def test_declared_md5_must_be_frozen(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                remote_zip.RemoteZipInventoryError, "MD5 contract"
            ):
                remote_zip.inventory_remote(
                    "https://example.invalid/archive.zip",
                    100,
                    "not-an-md5",
                    Path(directory) / "receipt.json",
                )


if __name__ == "__main__":
    unittest.main()
