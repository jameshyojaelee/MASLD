"""Structural proof that the protein_transport power floors were computed blind.

A power calculation informed by the effect it is meant to bound is circular. Promising
non-circularity in prose is unverifiable, so this walks the AST of the power module and
its whole import closure and asserts that none of them can perform I/O at all. A module
that cannot open a file cannot have consulted a result.
"""

from __future__ import annotations

import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# The power module plus every first-party module it imports. All three must be blind;
# checking only the top module would leave a hole one import deep.
BLIND_MODULES = (
    ROOT / "src" / "masld_bench" / "protein_transport_power.py",
    ROOT / "src" / "masld_bench" / "evaluators" / "auprc_reference.py",
    ROOT / "src" / "masld_bench" / "evaluators" / "metrics.py",
)

# Names that are I/O whatever the receiver. `get`, `run`, `post` and `connect` are
# deliberately NOT here: they are ordinary method names (`row.get(label, 0.0)` in
# metrics.py is a dict lookup), and flagging them by name alone produces false
# positives that would train a reader to ignore this test.
UNAMBIGUOUS_CALLS = frozenset(
    {
        "open",
        "read_text",
        "read_bytes",
        "write_text",
        "write_bytes",
        "urlopen",
        "urlretrieve",
        "Popen",
        "check_output",
        "check_call",
        "system",
        "popen",
        "input",
        "eval",
        "exec",
        "compile",
        "__import__",
        "getenv",
        "read_csv",
        "read_parquet",
        "loadtxt",
        "fromfile",
    }
)

# Receivers whose every call is suspect, so `json.load`, `subprocess.run` and
# `requests.get` are caught by qualification rather than by bare name.
IO_RECEIVERS = frozenset(
    {
        "os",
        "io",
        "sys",
        "json",
        "csv",
        "pathlib",
        "Path",
        "subprocess",
        "socket",
        "shutil",
        "sqlite3",
        "pickle",
        "urllib",
        "requests",
        "http",
        "tomllib",
        "np",
        "numpy",
        "pd",
        "pandas",
        "h5py",
    }
)

FORBIDDEN_IMPORTS = frozenset(
    {
        "os",
        "io",
        "sys",
        "json",
        "csv",
        "pathlib",
        "subprocess",
        "socket",
        "shutil",
        "sqlite3",
        "pickle",
        "urllib",
        "urllib.request",
        "requests",
        "http",
        "http.client",
        "tomllib",
        "numpy",
        "pandas",
        "h5py",
        "anndata",
        "glob",
        "tempfile",
        "importlib",
    }
)

ALLOWED_FIRST_PARTY_IMPORTS = frozenset(
    {
        "masld_bench.evaluators.metrics",
        ".evaluators.metrics",
        ".metrics",
    }
)


def _called_name(node: ast.Call) -> str:
    function = node.func
    if isinstance(function, ast.Name):
        return function.id
    if isinstance(function, ast.Attribute):
        return function.attr
    return ""


def _io_offence(node: ast.Call) -> str:
    """Return a description if this call can reach outside the process, else ""."""
    function = node.func
    if isinstance(function, ast.Name) and function.id in UNAMBIGUOUS_CALLS:
        return function.id
    if isinstance(function, ast.Attribute):
        if function.attr in UNAMBIGUOUS_CALLS:
            return function.attr
        receiver = function.value
        if isinstance(receiver, ast.Name) and receiver.id in IO_RECEIVERS:
            return f"{receiver.id}.{function.attr}"
    return ""


class PowerModuleIsBlind(unittest.TestCase):
    def test_every_blind_module_exists(self) -> None:
        for path in BLIND_MODULES:
            self.assertTrue(path.is_file(), path)

    def test_no_module_performs_input_or_output(self) -> None:
        for path in BLIND_MODULES:
            tree = ast.parse(path.read_text(), filename=str(path))
            offenders = sorted(
                {
                    f"{_io_offence(node)} at line {node.lineno}"
                    for node in ast.walk(tree)
                    if isinstance(node, ast.Call) and _io_offence(node)
                }
            )
            self.assertEqual(
                offenders,
                [],
                f"{path.name} can reach outside itself, so its output cannot be called blind: {offenders}",
            )

    def test_no_module_imports_an_io_capable_package(self) -> None:
        for path in BLIND_MODULES:
            tree = ast.parse(path.read_text(), filename=str(path))
            imported: set[str] = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        imported.add(alias.name)
                        imported.add(alias.name.split(".")[0])
                elif isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                    if node.level == 0:
                        imported.add(module)
                        imported.add(module.split(".")[0])
            self.assertEqual(
                sorted(imported & FORBIDDEN_IMPORTS),
                [],
                f"{path.name} imports an I/O capable package",
            )

    def test_no_module_reads_a_dunder_file_or_environment(self) -> None:
        for path in BLIND_MODULES:
            source = path.read_text()
            for token in ("__file__", "environ", "sys.argv", "stdin"):
                self.assertNotIn(token, source, f"{path.name} references {token}")

    def test_the_power_module_has_no_module_level_mutable_state(self) -> None:
        """A cached result smuggled in at import time would defeat the AST check."""
        path = BLIND_MODULES[0]
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        self.assertIsInstance(
                            node.value,
                            (ast.Constant, ast.Tuple, ast.Call),
                            f"module level {target.id} is not an immutable constant",
                        )
                        if isinstance(node.value, ast.Call):
                            self.assertEqual(_called_name(node.value), "frozenset")


class PowerFloorsAreFunctionsOfTheDesignAlone(unittest.TestCase):
    """The floors must move with n and alpha and with nothing else."""

    def setUp(self) -> None:
        from masld_bench import protein_transport_power as power

        self.power = power

    def test_spearman_floor_falls_as_n_rises(self) -> None:
        floors = [self.power.detectable_spearman(n) for n in (25, 56, 58, 200)]
        self.assertEqual(floors, sorted(floors, reverse=True))

    def test_spearman_floor_rises_as_alpha_tightens(self) -> None:
        loose = self.power.detectable_spearman(58, alpha=0.05)
        tight = self.power.detectable_spearman(58, alpha=0.05 / 6)
        self.assertGreater(tight, loose)

    def test_auroc_floor_is_worst_at_the_most_imbalanced_split(self) -> None:
        balanced = self.power.detectable_auroc(29, 29)
        imbalanced = self.power.detectable_auroc(14, 44)
        self.assertGreater(imbalanced, balanced)

    def test_bh_single_discovery_alpha_divides_by_the_family(self) -> None:
        self.assertAlmostEqual(self.power.bh_single_discovery_alpha(0.05, 6), 0.05 / 6)
        with self.assertRaises(self.power.PowerError):
            self.power.bh_single_discovery_alpha(0.05, 0)

    def test_average_precision_null_sits_well_above_prevalence(self) -> None:
        """The error this whole contract exists to prevent."""
        null_mean = self.power.average_precision_null_mean(14, 44, n_draws=4000)
        self.assertGreater(null_mean, self.power.prevalence(14, 44))

    def test_a_rare_level_cannot_reach_every_fold(self) -> None:
        table = self.power.smallest_estimable_stratum(
            {"F0": 13, "F1": 31, "F2": 11, "F3": 3}, outer_folds=5
        )
        self.assertFalse(table["F3"]["can_appear_in_every_fold"])
        self.assertEqual(table["F3"]["maximum_folds_reachable"], 3)
        self.assertTrue(table["F1"]["can_appear_in_every_fold"])

    def test_rank_estimability_needs_only_a_non_constant_fold(self) -> None:
        self.assertTrue(self.power.rank_endpoint_is_estimable([5, 6, 5, 9, 5]))
        self.assertFalse(self.power.rank_endpoint_is_estimable([5, 1, 5, 9, 5]))

    def test_classification_estimability_flags_the_vanishing_level(self) -> None:
        verdict = self.power.classification_endpoint_per_fold_estimability(
            {"MASH": [2, 4, 4, 6, 1], "MASL": [5, 9, 6, 5, 4], "No_MASLD": [3, 4, 0, 3, 2]}
        )
        self.assertEqual(verdict["levels_absent_from_some_fold"], ["No_MASLD"])
        self.assertFalse(verdict["per_fold_metric_estimable"])
        self.assertTrue(verdict["pooled_out_of_fold_estimable"])


if __name__ == "__main__":
    unittest.main()


class TheBlindnessCheckItselfWorks(unittest.TestCase):
    """A guard that cannot fail is not a guard. These are its positive controls."""

    def test_it_catches_a_bare_open(self) -> None:
        tree = ast.parse("def f():\n    return open('x').read()\n")
        self.assertTrue(any(_io_offence(n) for n in ast.walk(tree) if isinstance(n, ast.Call)))

    def test_it_catches_a_qualified_module_call(self) -> None:
        for source in ("import json\njson.load(f)\n", "requests.get(u)\n", "subprocess.run(c)\n"):
            tree = ast.parse(source)
            self.assertTrue(
                any(_io_offence(n) for n in ast.walk(tree) if isinstance(n, ast.Call)), source
            )

    def test_it_does_not_flag_an_ordinary_dict_get(self) -> None:
        tree = ast.parse("row.get(label, 0.0)\n")
        self.assertFalse(any(_io_offence(n) for n in ast.walk(tree) if isinstance(n, ast.Call)))

    def test_it_does_not_flag_an_ordinary_method_named_run(self) -> None:
        tree = ast.parse("simulation.run(draws)\n")
        self.assertFalse(any(_io_offence(n) for n in ast.walk(tree) if isinstance(n, ast.Call)))
