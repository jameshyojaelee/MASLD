"""Regression guard for the GSE267145 production child module search path.

SLURM jobs 21088904 and 21088907 wrote no predictions because the production
runner hands every child ``PYTHONPATH=<root>/src`` while the fitter chain
imports the repository-root ``scripts`` package and only the aggregator
imports ``masld_bench``.  These tests capture the environment the runner
actually builds and then execute both real entry points inside it.
"""

from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import scripts.run_gse267145_histology_production_bundle as v2


ROOT = Path(__file__).resolve().parents[2]
FITTER = ROOT / "scripts/fit_gse267145_histology_baselines_v4.py"
V3_FITTER = ROOT / "scripts/fit_gse267145_histology_baselines_v3.py"
AGGREGATOR = ROOT / "scripts/aggregate_gse267145_histology_production_predictions.py"


class _ChildEnvironmentCaptured(Exception):
    """Raised in place of the real child process once its environment is known."""


def _capture_child_environment(build: str) -> dict[str, str]:
    """Return the environment the runner passes to a fit or aggregate child."""

    captured: dict[str, str] = {}

    def _record(command, **keywords):
        captured.update(keywords["env"])
        raise _ChildEnvironmentCaptured

    with tempfile.TemporaryDirectory() as directory:
        arguments = Namespace(
            root=ROOT,
            campaign=Path(directory) / "campaign",
            python=Path(sys.executable),
            fitter=FITTER,
            aggregator=AGGREGATOR,
            molecular=ROOT / "molecular",
            folds=ROOT / "folds",
            fit_views=ROOT / "fit_views",
            surface=ROOT / "surface.json",
            blas_threads=3,
        )
        with patch.object(v2.subprocess, "run", _record):
            try:
                if build == "fit":
                    v2._run_unit(
                        arguments,
                        outer_fold=1,
                        seed=1721,
                        job_id="child-environment-test",
                        campaign_spec_sha256="0" * 64,
                    )
                else:
                    v2._aggregate(arguments, "child-environment-test")
            except _ChildEnvironmentCaptured:
                pass
            else:  # pragma: no cover - the runner must spawn a child
                raise AssertionError(f"runner did not spawn a {build} child process")
    return captured


class GSE267145ProductionChildEnvironmentTests(unittest.TestCase):
    def test_fit_child_environment_imports_scripts_package(self) -> None:
        environment = _capture_child_environment("fit")
        self.assertEqual(environment["PYTHONPATH"], (ROOT / "src").as_posix())
        completed = subprocess.run(
            [sys.executable, FITTER.as_posix(), "--software-self-test"],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertNotIn("ModuleNotFoundError", completed.stderr)
        self_test = json.loads(completed.stdout)
        self.assertEqual(self_test["status"], "passed_child_import_contract")
        self.assertTrue(self_test["repository_root_on_sys_path"])
        self.assertFalse(self_test["outcomes_read"])
        self.assertFalse(self_test["metrics_calculated"])
        self.assertFalse(self_test["scorer_called"])

    def test_aggregate_child_environment_imports_masld_bench(self) -> None:
        environment = _capture_child_environment("aggregate")
        self.assertEqual(environment["PYTHONPATH"], (ROOT / "src").as_posix())
        completed = subprocess.run(
            [sys.executable, AGGREGATOR.as_posix(), "--help"],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertNotIn("ModuleNotFoundError", completed.stderr)

    def test_unshimmed_v3_fitter_still_fails_in_the_child_environment(self) -> None:
        """The v4 shim, not the child environment, is what supplies ``scripts``."""

        environment = _capture_child_environment("fit")
        completed = subprocess.run(
            [sys.executable, V3_FITTER.as_posix(), "--software-self-test"],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("No module named 'scripts'", completed.stderr)


if __name__ == "__main__":
    unittest.main()
