from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from masld_bench.artifacts import (
    ArtifactError,
    freeze_tree,
    publish_directory_noreplace,
    verify_frozen_tree,
    write_text_exclusive,
)
from masld_bench.hashing import sha256_file


class ArtifactTests(unittest.TestCase):
    def test_directory_publication_never_replaces_an_empty_target(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "staging"
            source.mkdir()
            write_text_exclusive(source / "value.txt", "alpha\n")
            occupied = root / "occupied"
            occupied.mkdir()
            with self.assertRaisesRegex(ArtifactError, "refusing to replace"):
                publish_directory_noreplace(source, occupied)
            self.assertTrue(source.is_dir())
            self.assertEqual(list(occupied.iterdir()), [])

            target = root / "published"
            publish_directory_noreplace(source, target)
            self.assertFalse(source.exists())
            self.assertEqual((target / "value.txt").read_text(), "alpha\n")

    def test_reserved_directory_fallback_publishes_complete_last(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "staging"
            source.mkdir()
            write_text_exclusive(source / "value.txt", "alpha\n")
            freeze_tree(source, {"fixture": True})
            target = root / "published"

            with patch(
                "masld_bench.artifacts._renameat2_noreplace", return_value=False
            ):
                publish_directory_noreplace(source, target)

            self.assertFalse(source.exists())
            self.assertEqual(
                verify_frozen_tree(target)["metadata"], {"fixture": True}
            )

    def test_reserved_directory_fallback_never_replaces_target(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "staging"
            source.mkdir()
            write_text_exclusive(source / "value.txt", "alpha\n")
            target = root / "published"
            target.mkdir()

            with patch(
                "masld_bench.artifacts._renameat2_noreplace", return_value=False
            ), self.assertRaisesRegex(ArtifactError, "refusing to replace"):
                publish_directory_noreplace(source, target)

            self.assertTrue(source.is_dir())
            self.assertEqual(list(target.iterdir()), [])

    def test_reserved_directory_fallback_preserves_read_only_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "staging"
            source.mkdir()
            write_text_exclusive(source / "value.txt", "alpha\n")
            freeze_tree(source, {"fixture": True})
            source.chmod(0o555)
            target = root / "published"

            try:
                with patch(
                    "masld_bench.artifacts._renameat2_noreplace", return_value=False
                ):
                    publish_directory_noreplace(source, target)
                self.assertEqual(target.stat().st_mode & 0o777, 0o555)
                self.assertEqual(
                    verify_frozen_tree(target)["metadata"], {"fixture": True}
                )
            finally:
                if target.exists():
                    target.chmod(0o700)

    def test_frozen_tree_is_strict_and_exclusive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_text_exclusive(root / "value.txt", "alpha\n")
            with self.assertRaises(ArtifactError):
                write_text_exclusive(root / "value.txt", "beta\n")
            freeze_tree(root, {"fixture": True})
            manifest = verify_frozen_tree(root)
            self.assertEqual([item["path"] for item in manifest["artifacts"]], ["value.txt"])
            (root / "unexpected.txt").write_text("not frozen\n", encoding="utf-8")
            with self.assertRaisesRegex(ArtifactError, "unexpected"):
                verify_frozen_tree(root)

    def test_tamper_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_text_exclusive(root / "value.txt", "alpha\n")
            freeze_tree(root)
            (root / "value.txt").write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ArtifactError, "size mismatch|checksum mismatch"):
                verify_frozen_tree(root)

    def test_symlinks_are_rejected_from_frozen_trees(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            root = parent / "tree"
            root.mkdir()
            outside = parent / "outside.txt"
            outside.write_text("mutable\n", encoding="utf-8")
            os.symlink(outside, root / "linked.txt")
            with self.assertRaisesRegex(ArtifactError, "symlinks are forbidden"):
                freeze_tree(root)

    def test_root_symlinks_are_rejected_before_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            actual = parent / "actual"
            actual.mkdir()
            write_text_exclusive(actual / "value.txt", "alpha\n")
            linked = parent / "linked"
            os.symlink(actual, linked)
            with self.assertRaisesRegex(ArtifactError, "root.*symlink"):
                freeze_tree(linked)
            freeze_tree(actual)
            with self.assertRaisesRegex(ArtifactError, "root.*symlink"):
                verify_frozen_tree(linked)

            outer = parent / "outer"
            outer.mkdir()
            intermediate = outer / "intermediate"
            intermediate.symlink_to(parent, target_is_directory=True)
            with self.assertRaisesRegex(ArtifactError, "traverse a symlink"):
                verify_frozen_tree(intermediate / "actual")

    def test_exclusive_write_rejects_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            actual = parent / "actual"
            actual.mkdir()
            linked = parent / "linked"
            os.symlink(actual, linked)
            with self.assertRaisesRegex(ArtifactError, "traverse a symlink"):
                write_text_exclusive(linked / "escaped.txt", "unsafe\n")
            self.assertFalse((actual / "escaped.txt").exists())

    def test_control_documents_have_strict_schemas(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_text_exclusive(root / "value.txt", "alpha\n")
            freeze_tree(root)
            manifest_path = root / "ARTIFACTS.json"
            completion_path = root / "COMPLETE"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["unregistered"] = True
            manifest_path.write_text(
                json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            completion["manifest_sha256"] = sha256_file(manifest_path)
            completion_path.write_text(
                json.dumps(completion, sort_keys=True, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ArtifactError, "exactly its registered fields"):
                verify_frozen_tree(root)


if __name__ == "__main__":
    unittest.main()
