from __future__ import annotations

import importlib.util
from io import BytesIO
from pathlib import Path
import tarfile
import tempfile
import unittest


MODULE = Path(__file__).parents[2] / "scripts" / "epibert_safe_archive_inventory.py"
SPEC = importlib.util.spec_from_file_location("epibert_safe_inventory_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


class EpiBERTSafeArchiveInventoryTests(unittest.TestCase):
    def _archive(self, path: Path, name: str, *, linked: bool = False) -> None:
        with tarfile.open(path, "w:gz") as handle:
            member = tarfile.TarInfo(name)
            if linked:
                member.type = tarfile.SYMTYPE
                member.linkname = "elsewhere"
                handle.addfile(member)
            else:
                payload = b"checkpoint fixture"
                member.size = len(payload)
                handle.addfile(member, BytesIO(payload))

    def test_regular_archive_is_hashed_without_extraction(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "fixture.tar.gz"
            output = root / "inventory.json"
            self._archive(archive, "release/model.index")
            result = inventory.inventory(
                archive,
                output,
                expected_size=archive.stat().st_size,
                expected_md5=None,
                maximum_members=10,
                maximum_uncompressed_bytes=1024,
            )
            self.assertTrue(result["safe_inventory_passed"])
            self.assertFalse(result["archive_extracted"])
            self.assertEqual(result["checkpoint_index_members"], ["release/model.index"])
            self.assertFalse((root / "release").exists())

    def test_path_traversal_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "traversal.tar.gz"
            self._archive(archive, "../escape")
            with self.assertRaisesRegex(inventory.ArchiveInventoryError, "traversal"):
                inventory.inventory(
                    archive,
                    root / "inventory.json",
                    expected_size=None,
                    expected_md5=None,
                    maximum_members=10,
                    maximum_uncompressed_bytes=1024,
                )

    def test_link_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "link.tar.gz"
            self._archive(archive, "release/link", linked=True)
            with self.assertRaisesRegex(inventory.ArchiveInventoryError, "linked"):
                inventory.inventory(
                    archive,
                    root / "inventory.json",
                    expected_size=None,
                    expected_md5=None,
                    maximum_members=10,
                    maximum_uncompressed_bytes=1024,
                )


if __name__ == "__main__":
    unittest.main()
