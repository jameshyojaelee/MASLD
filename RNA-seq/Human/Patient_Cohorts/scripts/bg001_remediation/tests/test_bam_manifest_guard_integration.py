#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.util
import tempfile
import unittest
from fnmatch import fnmatch
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]


def load_prepare_run():
    spec = importlib.util.spec_from_file_location(
        "bg001_prepare_run_guard_test", SCRIPT_DIR / "prepare_run.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


prepare_run = load_prepare_run()


def load_source_validator():
    spec = importlib.util.spec_from_file_location(
        "bg001_source_validator_guard_test",
        SCRIPT_DIR / "validate_featurecounts_sources.py",
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


source_validator = load_source_validator()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BamManifestGuardIntegrationTests(unittest.TestCase):
    def test_guarded_hash_wrapper_preserves_resources_and_exact_inner_producer(self):
        wrapper = (SCRIPT_DIR / "bam_manifest_guarded_hash_array.sbatch").read_text()
        for directive in (
            "#SBATCH --job-name=samtools",
            "#SBATCH --partition=io",
            "#SBATCH --qos=interactive",
            "#SBATCH --array=0-7%4",
            "#SBATCH --cpus-per-task=1",
            "#SBATCH --mem=8G",
            "#SBATCH --time=48:00:00",
        ):
            self.assertIn(directive, wrapper)
        self.assertLess(wrapper.index("umask 077"), wrapper.index("exec bash"))
        self.assertIn('exec bash "$SCRIPT_DIR/bam_manifest_hash_array.sbatch"', wrapper)
        self.assertEqual(
            sha256(SCRIPT_DIR / "bam_manifest_hash_array.sbatch"),
            "8e30bb390d4a4eba3ccf4aa74f38d59ce4070a952f64e99d06a4da91dad90fd9",
        )

    def test_guarded_finalize_orders_pre_exact_producer_post_and_verify(self):
        wrapper = (SCRIPT_DIR / "bam_manifest_guarded_finalize.sbatch").read_text()
        for directive in (
            "#SBATCH --job-name=samtools",
            "#SBATCH --partition=io",
            "#SBATCH --qos=interactive",
            "#SBATCH --cpus-per-task=1",
            "#SBATCH --mem=8G",
            "#SBATCH --time=48:00:00",
        ):
            self.assertIn(directive, wrapper)
        pre = wrapper.index('bam_manifest_finalization_guard.py" pre')
        producer = wrapper.index('bash "$SCRIPT_DIR/bam_manifest_finalize.sbatch"')
        post = wrapper.index('bam_manifest_finalization_guard.py" post')
        verify = wrapper.index('bam_manifest_finalization_guard.py" verify')
        self.assertLess(wrapper.index("umask 077"), pre)
        self.assertLess(wrapper.index("micromamba activate rnaseq"), pre)
        self.assertLess(pre, producer)
        self.assertLess(producer, post)
        self.assertLess(post, verify)
        self.assertEqual(
            sha256(SCRIPT_DIR / "bam_manifest_finalize.sbatch"),
            "cfafb5b3ec934f27167b262e46ac0d93be0f642e34b588a6c3f29070a7e2f5c1",
        )

    def test_prepare_security_check_rejects_group_writable_private_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "candidate"
            path.mkdir(mode=0o700)
            path.chmod(0o770)
            with self.assertRaisesRegex(SystemExit, "owner-only"):
                prepare_run.require_secure_directory(
                    path, "candidate fixture", owner_private=True
                )

    def test_prepare_security_check_rejects_symlink(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            target = base / "target"
            target.mkdir(mode=0o700)
            link = base / "link"
            link.symlink_to(target)
            with self.assertRaisesRegex(SystemExit, "non-symlink"):
                prepare_run.require_secure_directory(
                    link, "candidate fixture", owner_private=True
                )

    def test_prepare_source_requests_private_modes_and_umask(self):
        source = (SCRIPT_DIR / "prepare_run.py").read_text()
        self.assertIn("os.umask(0o077)", source)
        self.assertIn("run_root.mkdir(mode=0o700)", source)
        self.assertIn("directory.mkdir(mode=0o700)", source)
        self.assertIn("os.chmod(draft_manifest, 0o600)", source)
        self.assertLess(source.index("live_source_gate ="), source.index("shutil.copytree"))
        self.assertIn('"LIVE WORKTREE VALIDATION\\n"', source)
        self.assertIn('"FROZEN SNAPSHOT VALIDATION\\n"', source)
        self.assertIn(
            '"RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/workflow/config.yaml"',
            source,
        )
        self.assertIn(
            '"RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/workflow/config.yaml"',
            source,
        )
        self.assertIn(
            '"RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE130970/scripts/submit_pipeline.sh"',
            source,
        )
        self.assertIn(
            '"RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE135251/scripts/submit_pipeline.sh"',
            source,
        )

    def test_prepare_snapshot_closes_over_validator_required_paths(self):
        required = set(
            source_validator.HUMAN_REQUIRED
            + source_validator.DYNAMIC_REQUIRED
            + source_validator.GUARDED_CUSTOM_SNAKEFILES
            + source_validator.SHARED_AFFECTED_CONFIGS
            + source_validator.MOUSE_REQUIRED
            + source_validator.RETIRED_REQUIRED
            + source_validator.DISABLED_LIVE_PRODUCERS
        )
        required.update(
            launcher for _, launcher, _ in source_validator.ACTIVE_COUNT_ROUTES
        )
        explicit = set(prepare_run.SOURCE_FILES)
        missing = {
            relative
            for relative in required
            if not relative.startswith(
                "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/"
            )
            and relative not in explicit
            and not any(fnmatch(relative, pattern) for pattern in prepare_run.SOURCE_GLOBS)
        }
        self.assertEqual(missing, set())


if __name__ == "__main__":
    unittest.main()
