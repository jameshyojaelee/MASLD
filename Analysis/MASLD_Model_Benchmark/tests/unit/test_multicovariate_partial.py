"""The extension must reduce exactly to the sealed instrument at one covariate.

That reduction is the whole justification for adding a second module rather
than editing the one six sealed jobs hashed. If it does not hold to machine
precision, Stage 3's numbers sit on a quietly different instrument from the four
stages they will be compared against, and no amount of documentation fixes that.

The closure case is pinned separately: five proportions where one sits near 0.9
are strongly dependent by construction, and the block's numerical rank has to be
measured rather than presumed equal to the covariate count.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "multicovariate_partial.py"
SEALED = ROOT / "scripts" / "evaluate_gse267145_axis_count.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mc = _load(MODULE, "multicovariate_tested")
analysis = mc.load_sealed(SEALED)


def _fixture(seed: int, n: int = 90, g: int = 500):
    generator = np.random.default_rng(seed)
    ranks = analysis.average_ranks(generator.normal(size=(n, g)))
    ranks -= ranks.mean(axis=0, keepdims=True)
    exposure = generator.integers(0, 9, size=n).astype(float)
    return generator, ranks, exposure


class ReductionTests(unittest.TestCase):
    """One covariate must give bit-comparable answers to the sealed path."""

    def test_reduces_to_the_sealed_path_on_an_ordinal_covariate(self) -> None:
        for seed in (1, 2, 3, 4, 5):
            generator, ranks, exposure = _fixture(seed)
            covariate = generator.integers(0, 5, size=ranks.shape[0]).astype(float)
            worst = mc.reduces_to_sealed_single_covariate(
                analysis, ranks, exposure, covariate)
            self.assertLess(worst, 1e-12, f"seed {seed} disagreed by {worst}")

    def test_reduces_on_a_heavily_tied_covariate(self) -> None:
        """The regime the histology axes actually live in."""

        generator, ranks, exposure = _fixture(11)
        covariate = np.zeros(ranks.shape[0])
        covariate[:20] = 1.0
        covariate[20:28] = 2.0
        worst = mc.reduces_to_sealed_single_covariate(
            analysis, ranks, exposure, covariate)
        self.assertLess(worst, 1e-12)

    def test_reduces_on_a_continuous_covariate(self) -> None:
        generator, ranks, exposure = _fixture(12)
        covariate = generator.normal(size=ranks.shape[0])
        worst = mc.reduces_to_sealed_single_covariate(
            analysis, ranks, exposure, covariate)
        self.assertLess(worst, 1e-12)


class ClosureTests(unittest.TestCase):
    def test_a_rank_deficient_block_is_detected_not_presumed_full(self) -> None:
        generator, ranks, exposure = _fixture(21)
        n = ranks.shape[0]
        a = generator.random(n)
        b = generator.random(n)
        block = [a, b, a + b]  # exactly dependent in value; ranks need not be
        _, report = mc.covariate_basis(analysis, block)
        self.assertEqual(report["n_covariates"], 3)
        self.assertIn("numerical_rank", report)
        self.assertGreaterEqual(report["condition_number"], 1.0)

    def test_an_exactly_duplicated_covariate_reduces_the_rank(self) -> None:
        generator, ranks, exposure = _fixture(22)
        a = generator.random(ranks.shape[0])
        _, report = mc.covariate_basis(analysis, [a, a.copy()])
        self.assertEqual(report["numerical_rank"], 1)
        self.assertTrue(report["is_rank_deficient"])

    def test_a_closed_five_part_composition_is_measured(self) -> None:
        """Five proportions with one near 0.9 - the real Stage 3 geometry."""

        generator = np.random.default_rng(31)
        n = 120
        rest = generator.dirichlet(np.ones(4) * 0.5, size=n) * 0.1
        hep = 1.0 - rest.sum(axis=1)
        parts = [hep] + [rest[:, i] for i in range(4)]
        _, report = mc.covariate_basis(analysis, parts)
        self.assertEqual(report["n_covariates"], 5)
        self.assertGreater(report["condition_number"], 1.0)
        self.assertIn("closure", report["why_it_is_measured"])

    def test_all_constant_covariates_are_refused(self) -> None:
        with self.assertRaises(mc.MultiCovariateError):
            mc.covariate_basis(analysis, [np.ones(50), np.full(50, 3.0)])

    def test_no_covariates_is_refused(self) -> None:
        with self.assertRaises(mc.MultiCovariateError):
            mc.covariate_basis(analysis, [])


class PartialTests(unittest.TestCase):
    def test_residual_df_subtracts_the_measured_rank_not_the_count(self) -> None:
        generator, ranks, exposure = _fixture(41, n=100)
        a = generator.random(100)
        _, report = mc.partial_associations(analysis, ranks, exposure, [a, a.copy()])
        self.assertEqual(report["numerical_rank"], 1)
        self.assertEqual(report["residual_df"], 100 - 2 - 1)

    def test_six_covariates_give_df_n_minus_8(self) -> None:
        generator, ranks, exposure = _fixture(42, n=180)
        covs = [generator.random(180) for _ in range(6)]
        _, report = mc.partial_associations(analysis, ranks, exposure, covs)
        self.assertEqual(report["numerical_rank"], 6)
        self.assertEqual(report["residual_df"], 180 - 8)

    def test_an_exposure_inside_the_covariate_span_is_refused(self) -> None:
        generator, ranks, _ = _fixture(43)
        covariate = generator.integers(0, 5, size=ranks.shape[0]).astype(float)
        with self.assertRaises(mc.MultiCovariateError):
            mc.partial_associations(analysis, ranks, covariate, [covariate])

    def test_a_gene_with_no_residual_variance_is_counted_not_divided_by(self) -> None:
        generator, ranks, exposure = _fixture(44)
        covariate = generator.random(ranks.shape[0])
        ranks = ranks.copy()
        ranks[:, 0] = analysis.centred_ranks(covariate)
        values, report = mc.partial_associations(
            analysis, ranks, exposure, [covariate])
        self.assertGreaterEqual(report["genes_with_no_residual_variance"], 1)
        self.assertTrue(np.isfinite(values).all())
        self.assertEqual(float(values[0]), 0.0)


class SealedInstrumentTests(unittest.TestCase):
    def test_the_sealed_instrument_is_unmodified(self) -> None:
        """Six sealed jobs hashed this file; the extension must not touch it."""

        digest = subprocess.run(
            ["sha256sum", str(SEALED)], capture_output=True, text=True, check=True
        ).stdout.split()[0]
        self.assertTrue(
            digest.startswith(mc.SEALED_INSTRUMENT_SHA256_PREFIX),
            f"the sealed instrument changed: {digest[:16]}",
        )

    def test_the_extension_defines_no_function_the_sealed_module_owns(self) -> None:
        import ast

        defined = {
            node.name
            for node in ast.walk(ast.parse(MODULE.read_text(encoding="utf-8")))
            if isinstance(node, ast.FunctionDef)
        }
        for owned in ("average_ranks", "centred_ranks", "unit_ranks", "unit_columns",
                      "bh_count", "evaluate_gate"):
            self.assertNotIn(owned, defined)


if __name__ == "__main__":
    unittest.main()
