from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.refresh_pin_companion_metadata import apply, stale_companions


class CompanionRefreshTests(unittest.TestCase):
    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.root = Path(holder.name)
        (self.root / "config").mkdir()
        self.target = self.root / "config/t.toml"
        self.target.write_text("x = 1\n", encoding="utf-8")
        self.record = self.root / "r.json"

    def write(self, node: dict) -> None:
        self.record.write_text(json.dumps({"source_records": [node]}), encoding="utf-8")

    def test_stale_size_is_detected_and_corrected(self) -> None:
        self.write({"path": "config/t.toml", "sha256": "a" * 64, "size_bytes": 1})
        found = stale_companions(self.root, self.record)
        self.assertEqual(found[0]["to"], self.target.stat().st_size)
        apply(self.root, self.record)
        doc = json.loads(self.record.read_text())
        self.assertEqual(doc["source_records"][0]["size_bytes"],
                         self.target.stat().st_size)

    def test_correct_size_is_left_alone(self) -> None:
        self.write({"path": "config/t.toml", "sha256": "a" * 64,
                    "size_bytes": self.target.stat().st_size})
        self.assertEqual(stale_companions(self.root, self.record), [])
        self.assertEqual(apply(self.root, self.record), [])

    def test_unresolvable_path_is_skipped_not_guessed(self) -> None:
        self.write({"path": "config/absent.toml", "size_bytes": 99})
        self.assertEqual(stale_companions(self.root, self.record), [])

    def test_digest_currency_is_reported_so_provenance_stays_visible(self) -> None:
        """Flags whether the companion belongs to a digest already current."""

        self.write({"path": "config/t.toml", "sha256": "a" * 64, "size_bytes": 1})
        self.assertFalse(stale_companions(self.root, self.record)[0]["digest_current"])

    def test_only_size_fields_move(self) -> None:
        self.write({"path": "config/t.toml", "sha256": "a" * 64, "size_bytes": 1,
                    "note": "keep me"})
        apply(self.root, self.record)
        doc = json.loads(self.record.read_text())["source_records"][0]
        self.assertEqual(doc["sha256"], "a" * 64)
        self.assertEqual(doc["note"], "keep me")


if __name__ == "__main__":
    unittest.main()
