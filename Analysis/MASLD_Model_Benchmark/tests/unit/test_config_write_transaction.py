from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.config_write_transaction import (
    TransactionError,
    digest,
    restore,
    snapshot,
)


class TransactionTests(unittest.TestCase):
    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.root = Path(holder.name)
        (self.root / "config").mkdir()
        self.a = self.root / "config/a.json"
        self.b = self.root / "config/b.json"
        self.a.write_text('{"v": 1}', encoding="utf-8")
        self.b.write_text('{"v": 2}', encoding="utf-8")
        self.into = self.root / "staging/rollback"

    def test_snapshot_records_pre_write_digests(self) -> None:
        m = snapshot(self.root, ["config/a.json", "config/b.json"], self.into)
        self.assertEqual(len(m["files"]), 2)
        self.assertEqual(m["files"][0]["sha256_before"], digest(self.a))

    def test_restore_undoes_a_partial_write_and_verifies_it(self) -> None:
        """The failure this exists for: a run wrote, then died."""

        before = digest(self.a)
        snapshot(self.root, ["config/a.json", "config/b.json"], self.into)
        self.a.write_text('{"v": 999}', encoding="utf-8")
        result = restore(self.root, self.into)
        self.assertEqual(digest(self.a), before)
        self.assertEqual(result["files_the_run_had_modified"], 1)
        self.assertTrue(result["all_verified"])

    def test_restore_reports_which_files_the_run_actually_moved(self) -> None:
        snapshot(self.root, ["config/a.json", "config/b.json"], self.into)
        self.a.write_text('{"v": 3}', encoding="utf-8")
        result = restore(self.root, self.into)
        moved = {r["path"]: r["was_modified_by_the_run"] for r in result["restore_results"]}
        self.assertTrue(moved["config/a.json"])
        self.assertFalse(moved["config/b.json"])

    def test_untouched_run_restores_to_a_no_op(self) -> None:
        snapshot(self.root, ["config/a.json"], self.into)
        result = restore(self.root, self.into)
        self.assertEqual(result["files_the_run_had_modified"], 0)

    def test_snapshot_refuses_a_missing_target(self) -> None:
        with self.assertRaises(TransactionError):
            snapshot(self.root, ["config/absent.json"], self.into)

    def test_snapshot_refuses_to_overwrite_an_existing_snapshot(self) -> None:
        snapshot(self.root, ["config/a.json"], self.into)
        with self.assertRaises(TransactionError):
            snapshot(self.root, ["config/a.json"], self.into)

    def test_restore_without_a_manifest_is_refused(self) -> None:
        (self.root / "empty").mkdir()
        with self.assertRaises(TransactionError):
            restore(self.root, self.root / "empty")

    def test_corrupt_snapshot_copy_fails_loudly_rather_than_silently(self) -> None:
        snapshot(self.root, ["config/a.json"], self.into)
        stored = next(p for p in self.into.iterdir() if p.name != "MANIFEST.json")
        stored.write_text('{"v": "corrupted"}', encoding="utf-8")
        with self.assertRaises(TransactionError) as caught:
            restore(self.root, self.into)
        self.assertIn("restored to", str(caught.exception))

    def test_manifest_records_the_restore_for_the_receipt(self) -> None:
        snapshot(self.root, ["config/a.json"], self.into)
        self.a.write_text("x", encoding="utf-8")
        restore(self.root, self.into)
        m = json.loads((self.into / "MANIFEST.json").read_text())
        self.assertTrue(m["restored"])
        self.assertTrue(m["all_verified"])


if __name__ == "__main__":
    unittest.main()
