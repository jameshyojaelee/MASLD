from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from scripts.refresh_authorized_pins import (
    PinRefreshError,
    apply_refresh,
    navigate,
    plan_refresh,
)


STALE = "a" * 64


def build_tree(directory: Path, recorded: str, target_text: str) -> str:
    (directory / "config/models").mkdir(parents=True)
    (directory / "config/artifacts/models/demo").mkdir(parents=True)
    target = directory / "config/models/mandatory_baselines.toml"
    target.write_text(target_text, encoding="utf-8")
    bundle = directory / "config/artifacts/models/demo/checkpoints.json"
    bundle.write_text(
        json.dumps({"source_records": [{"identity": "x", "sha256": recorded}]}),
        encoding="utf-8",
    )
    from hashlib import sha256

    return sha256(target.read_bytes()).hexdigest()


def authorization(live: str, stale: str = STALE, pinned_by=None) -> dict:
    return {
        "authorization_id": "test",
        "authorized_by": "test",
        "authorized_pins": [
            {
                "bundle_path": "config/artifacts/models/demo/checkpoints.json",
                "json_path": ".source_records[0]",
                "key": "sha256",
                "target_path": "config/models/mandatory_baselines.toml",
                "expected_stale_sha256": stale,
                "expected_live_sha256": live,
                "bundle_is_pinned_by": pinned_by or [],
            }
        ],
    }


class NavigateTests(unittest.TestCase):
    def test_indexed_and_nested_paths(self) -> None:
        doc = {"a": [{"b": 1}], "c": {"d": {"e": 2}}}
        self.assertEqual(navigate(doc, ".a[0]"), {"b": 1})
        self.assertEqual(navigate(doc, ".c.d"), {"e": 2})


class PlanRefreshTests(unittest.TestCase):
    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.root = Path(holder.name)

    def test_valid_authorization_plans_one_action(self) -> None:
        live = build_tree(self.root, STALE, "baselines = 1\n")
        actions = plan_refresh(self.root, authorization(live))
        self.assertEqual(actions[0]["from"], STALE)
        self.assertEqual(actions[0]["to"], live)

    def test_unexpected_recorded_value_aborts(self) -> None:
        live = build_tree(self.root, "b" * 64, "baselines = 1\n")
        with self.assertRaises(PinRefreshError) as caught:
            plan_refresh(self.root, authorization(live))
        self.assertIn("not the authorized stale value", str(caught.exception))

    def test_source_moving_again_aborts_rather_than_writing(self) -> None:
        """If the target moved since authorization, the authorization is stale."""

        build_tree(self.root, STALE, "baselines = 1\n")
        with self.assertRaises(PinRefreshError) as caught:
            plan_refresh(self.root, authorization("c" * 64))
        self.assertIn("moved again", str(caught.exception))

    def test_nothing_is_written_when_planning_fails(self) -> None:
        build_tree(self.root, STALE, "baselines = 1\n")
        bundle = self.root / "config/artifacts/models/demo/checkpoints.json"
        before = bundle.read_text()
        with self.assertRaises(PinRefreshError):
            plan_refresh(self.root, authorization("c" * 64))
        self.assertEqual(bundle.read_text(), before)

    def test_apply_rewrites_only_the_named_pin(self) -> None:
        live = build_tree(self.root, STALE, "baselines = 1\n")
        actions = apply_refresh(self.root, plan_refresh(self.root, authorization(live)))
        bundle = json.loads(
            (self.root / "config/artifacts/models/demo/checkpoints.json").read_text()
        )
        self.assertEqual(bundle["source_records"][0]["sha256"], live)
        self.assertEqual(bundle["source_records"][0]["identity"], "x")
        self.assertNotEqual(
            actions[0]["bundle_sha256_before"], actions[0]["bundle_sha256_after"]
        )


if __name__ == "__main__":
    unittest.main()


class BranchClosureGuardTests(unittest.TestCase):
    """A pinned bundle may only move inside a declared, complete branch."""

    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.root = Path(holder.name)
        self.live = build_tree(self.root, STALE, "baselines = 1\n")

    def _auth(self, **over):
        a = authorization(self.live, pinned_by=["config/root.json"])
        a.update(over)
        return a

    def test_pinned_bundle_without_a_declared_closure_is_refused(self) -> None:
        with self.assertRaises(PinRefreshError) as c:
            plan_refresh(self.root, self._auth())
        self.assertIn("no branch_closure is declared", str(c.exception))

    def test_pinner_outside_the_declared_closure_is_refused(self) -> None:
        with self.assertRaises(PinRefreshError) as c:
            plan_refresh(self.root, self._auth(branch_closure=["config/other.json"]))
        self.assertIn("declared branch_closure does not cover", str(c.exception))

    def test_declared_pinner_that_this_run_does_not_refresh_is_refused(self) -> None:
        """Naming the closure is not enough; the run must actually update it."""

        with self.assertRaises(PinRefreshError) as c:
            plan_refresh(self.root, self._auth(branch_closure=["config/root.json"]))
        self.assertIn("does not re-freeze", str(c.exception))

    def test_unpinned_bundle_still_needs_no_closure(self) -> None:
        actions = plan_refresh(self.root, authorization(self.live))
        self.assertEqual(len(actions), 1)


class DeferredTargetTests(unittest.TestCase):
    """A pin whose target moves earlier in the same run defers its value."""

    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.root = Path(holder.name)
        self.live = build_tree(self.root, STALE, "baselines = 1\n")
        root_doc = self.root / "config/root.json"
        root_doc.write_text(json.dumps({"authorities": {"demo": {
            "path": "config/artifacts/models/demo/checkpoints.json",
            "sha256": "f" * 64}}}))

    def _auth(self, defer_target: str):
        seed = authorization(self.live, pinned_by=["config/root.json"])
        seed["branch_closure"] = ["config/root.json",
                                  "config/artifacts/models/demo/checkpoints.json"]
        seed["authorized_pins"].append({
            "bundle_path": "config/root.json",
            "json_path": ".authorities.demo",
            "key": "sha256",
            "target_path": defer_target,
            "expected_stale_sha256": "f" * 64,
            "expected_live_sha256": "RECOMPUTED_AFTER_SEEDS",
            "bundle_is_pinned_by": []})
        return seed

    def test_deferred_pin_resolves_to_the_post_refresh_digest(self) -> None:
        target = "config/artifacts/models/demo/checkpoints.json"
        actions = plan_refresh(self.root, self._auth(target))
        apply_refresh(self.root, actions)
        written = json.loads((self.root / "config/root.json").read_text())
        actual = hashlib.sha256((self.root / target).read_bytes()).hexdigest()
        self.assertEqual(written["authorities"]["demo"]["sha256"], actual)
        self.assertNotEqual(actual, "f" * 64)

    def test_deferring_to_a_target_this_run_does_not_refresh_is_refused(self) -> None:
        with self.assertRaises(PinRefreshError) as c:
            plan_refresh(self.root, self._auth("config/models/mandatory_baselines.toml"))
        self.assertIn("deferral requires the target to move first", str(c.exception))


class TomlPinTests(unittest.TestCase):
    """TOML pins are edited surgically, never round-tripped."""

    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.path = Path(holder.name) / "a.toml"
        self.path.write_text(
            '# a comment that must survive\n'
            'schema_version = "v1"\n\n'
            '[frozen_file_authorities.other]\n'
            'path = "config/other.json"\n'
            'sha256 = "' + "b" * 64 + '"\n\n'
            '[frozen_file_authorities.target]\n'
            'path = "config/target.json"\n'
            'sha256 = "' + "a" * 64 + '"\n',
            encoding="utf-8")

    def test_reads_the_right_table(self) -> None:
        from scripts.refresh_authorized_pins import read_toml_pin
        self.assertEqual(
            read_toml_pin(self.path, ".frozen_file_authorities.target", "sha256"),
            "a" * 64)

    def test_write_touches_one_line_and_keeps_comments(self) -> None:
        from scripts.refresh_authorized_pins import read_toml_pin, write_toml_pin
        before = self.path.read_text().splitlines()
        write_toml_pin(self.path, ".frozen_file_authorities.target", "sha256", "c" * 64)
        after = self.path.read_text().splitlines()
        self.assertEqual(len(before), len(after))
        self.assertEqual(sum(1 for x, y in zip(before, after) if x != y), 1)
        self.assertIn("# a comment that must survive", self.path.read_text())
        self.assertEqual(
            read_toml_pin(self.path, ".frozen_file_authorities.other", "sha256"),
            "b" * 64, "the sibling table must not move")
        self.assertEqual(
            read_toml_pin(self.path, ".frozen_file_authorities.target", "sha256"),
            "c" * 64)

    def test_missing_table_is_refused(self) -> None:
        from scripts.refresh_authorized_pins import write_toml_pin
        with self.assertRaises(PinRefreshError):
            write_toml_pin(self.path, ".frozen_file_authorities.absent", "sha256", "d" * 64)
