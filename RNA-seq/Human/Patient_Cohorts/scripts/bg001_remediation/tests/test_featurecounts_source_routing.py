#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SCRIPT_DIR.parents[4]


def load_validator():
    path = SCRIPT_DIR / "validate_featurecounts_sources.py"
    spec = importlib.util.spec_from_file_location("bg001_source_routes", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


validator = load_validator()


class ActiveFeatureCountsRouteTests(unittest.TestCase):
    @staticmethod
    def copy_contract_tree(destination: Path) -> None:
        paths = set(
            validator.GUARDED_CUSTOM_SNAKEFILES
            + validator.DYNAMIC_REQUIRED
            + validator.SHARED_AFFECTED_CONFIGS
            + tuple(route[1] for route in validator.ACTIVE_COUNT_ROUTES)
        )
        for relative in paths:
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(PROJECT_ROOT / relative, target)

    def test_authoritative_launcher_workflow_config_composition_is_fail_closed(self):
        self.assertEqual(validator.validate_actual_count_routes(PROJECT_ROOT), [])

    def test_each_custom_count_rule_guard_is_mandatory(self):
        for relative in validator.GUARDED_CUSTOM_SNAKEFILES:
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                self.copy_contract_tree(root)
                path = root / relative
                path.write_text(path.read_text().replace("exit 64", "exit 0", 1))
                failures = validator.validate_actual_count_routes(root)
                self.assertTrue(
                    any(f"custom live count rule is not fail-closed: {relative}" == item for item in failures),
                    failures,
                )

    def test_shared_count_merge_guards_and_exact_strand_are_mandatory(self):
        mutations = (
            (
                validator.DYNAMIC_REQUIRED[0],
                "\"echo '{BG001_RECOUNT_REFUSAL}' >&2; exit 64\"",
                "\"echo '{BG001_RECOUNT_REFUSAL}' >&2; exit 0\"",
                "shared affected-dataset count/merge guards are incomplete",
            ),
            (
                validator.DYNAMIC_REQUIRED[0],
                "if BG001_RECOUNT_REQUIRED:",
                "if False:",
                "shared affected-dataset count/merge guards are incomplete",
            ),
            (
                validator.DYNAMIC_REQUIRED[0],
                "BG001_MERGE_GUARD_OUTPUTS",
                "BG001_BROKEN_GUARD_OUTPUTS",
                "shared affected-dataset count/merge guards are incomplete",
            ),
            (
                validator.SHARED_AFFECTED_CONFIGS[0],
                "strandedness: 2",
                "strandedness: 0",
                "affected shared-workflow config is not exact -s 2",
            ),
        )
        for relative, old, new, expected in mutations:
            with self.subTest(relative=relative, old=old), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                self.copy_contract_tree(root)
                path = root / relative
                text = path.read_text()
                self.assertIn(old, text)
                path.write_text(text.replace(old, new, 1))
                failures = validator.validate_actual_count_routes(root)
                self.assertTrue(any(expected in item for item in failures), failures)

    def test_active_shared_launcher_cannot_switch_to_retroactive_direct_workflow(self):
        relative = next(
            launcher for dataset, launcher, kind in validator.ACTIVE_COUNT_ROUTES
            if dataset == "GSE174478" and kind == "shared"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.copy_contract_tree(root)
            launcher = root / relative
            text = launcher.read_text()
            self.assertIn("pipelines/shared/Snakefile", text)
            launcher.write_text(text.replace("pipelines/shared/Snakefile", "pipelines/custom/GSE174478/workflow/Snakefile"))
            failures = validator.validate_actual_count_routes(root)
            self.assertTrue(any("does not select shared Snakefile" in item for item in failures), failures)


if __name__ == "__main__":
    unittest.main()
