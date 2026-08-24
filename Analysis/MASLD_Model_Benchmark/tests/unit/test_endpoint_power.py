"""Tests for registered endpoint-recomputed prospective-power inputs."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import masld_bench.evaluators.endpoint_power as endpoint_power
from masld_bench.evaluators.endpoint_power import (
    EndpointPowerError,
    recompute_endpoint_bootstrap,
    registered_endpoint_evaluators,
)


def _hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _write_table(
    root: Path, name: str, header: tuple[str, ...], rows: list[tuple[object, ...]]
) -> Path:
    path = root / name
    lines = ["\t".join(header)]
    lines.extend("\t".join(str(value) for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class EndpointPowerTests(unittest.TestCase):
    def tearDown(self) -> None:
        endpoint_power._clear_endpoint_cache_for_testing()

    def test_registry_is_exact_and_unknown_evaluator_fails_closed(self) -> None:
        self.assertEqual(
            registered_endpoint_evaluators(),
            (
                "bulk_paired_spearman_gain_v1",
                "cell_donor_balanced_macro_f1_v1",
                "graph_ld_block_auprc_gain_v1",
                "locus_heldout_allelic_spearman_v1",
                "rna_atac_two_way_deviance_reduction_v1",
                "variant_ld_block_fisher_z_spearman_gain_v1",
            ),
        )
        with self.assertRaisesRegex(EndpointPowerError, "unregistered"):
            recompute_endpoint_bootstrap(
                evaluator_id="caller_supplied_summary",
                table_path="unused.tsv",
                n_resamples=100,
                seed=1,
                parameters={},
            )

    def test_mpra_proxy_is_locus_bootstrapped_and_not_an_eqtl_endpoint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            for locus in range(20):
                block_hash = _hash(f"locus-{locus}")
                for allele in range(2):
                    observed = float(2 * locus + allele)
                    rows.append(
                        (
                            _hash(f"mpra-row-{locus}-{allele}"),
                            _hash(f"mpra-unit-{locus}-{allele}"),
                            block_hash,
                            observed,
                            observed,
                            -observed,
                        )
                    )
            table = _write_table(
                root,
                "mpra.tsv",
                (
                    "row_hash",
                    "unit_hash",
                    "block_hash",
                    "observed",
                    "candidate",
                    "baseline",
                ),
                rows,
            )
            result = recompute_endpoint_bootstrap(
                evaluator_id="locus_heldout_allelic_spearman_v1",
                table_path=table,
                n_resamples=100,
                seed=43,
                parameters={},
            )
            self.assertEqual(result.independent_unit_counts, {"merged_locus": 20})
            self.assertGreater(result.observed_effect, 0.0)
            self.assertNotIn("eqtl", result.evaluator_id)

    def test_cell_endpoint_recomputes_macro_f1_at_donor_level(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows: list[tuple[object, ...]] = []
            observations = (
                ("d1", "a", "a", "a"),
                ("d1", "b", "b", "a"),
                ("d2", "a", "a", "a"),
                ("d2", "b", "b", "a"),
                ("d3", "a", "a", "a"),
                ("d3", "b", "b", "a"),
            )
            for index, (donor, truth, candidate, baseline) in enumerate(observations):
                rows.append(
                    (_hash(f"row-{index}"), _hash(donor), truth, candidate, baseline)
                )
            table = _write_table(
                root,
                "cell.tsv",
                (
                    "row_hash",
                    "unit_hash",
                    "observed_class",
                    "candidate_class",
                    "baseline_class",
                ),
                rows,
            )
            first = recompute_endpoint_bootstrap(
                evaluator_id="cell_donor_balanced_macro_f1_v1",
                table_path=table,
                n_resamples=100,
                seed=41,
                parameters={"class_roster": ["a", "b"]},
            )
            second = recompute_endpoint_bootstrap(
                evaluator_id="cell_donor_balanced_macro_f1_v1",
                table_path=table,
                n_resamples=100,
                seed=41,
                parameters={"class_roster": ["a", "b"]},
            )
            self.assertEqual(first, second)
            self.assertGreater(first.observed_effect, 0.0)
            self.assertAlmostEqual(
                first.observed_effect,
                first.candidate_primary_metric - first.baseline_primary_metric,
            )
            self.assertTrue(
                all(
                    abs((candidate - baseline) - difference) < 1e-12
                    for candidate, baseline, difference in zip(
                        first.candidate_primary_bootstrap,
                        first.baseline_primary_bootstrap,
                        first.paired_difference_bootstrap,
                        strict=True,
                    )
                )
            )
            self.assertEqual(first.independent_unit_counts, {"donor": 3})
            self.assertEqual(first.n_rows, 6)

    def test_bulk_endpoint_resamples_donors_and_rejects_duplicate_units(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            header = (
                "row_hash",
                "unit_hash",
                "observed",
                "candidate",
                "baseline",
            )
            rows = [
                (_hash(f"row-{index}"), _hash(f"d{index}"), value, value, -value)
                for index, value in enumerate(
                    (float(2**power) for power in range(20)), start=1
                )
            ]
            table = _write_table(root, "bulk.tsv", header, rows)
            result = recompute_endpoint_bootstrap(
                evaluator_id="bulk_paired_spearman_gain_v1",
                table_path=table,
                n_resamples=100,
                seed=17,
                parameters={},
            )
            self.assertAlmostEqual(result.observed_effect, 2.0)
            self.assertAlmostEqual(result.candidate_primary_metric, 1.0)
            self.assertAlmostEqual(result.baseline_primary_metric, -1.0)
            self.assertEqual(result.n_units, 20)

            duplicate = list(rows)
            duplicate[-1] = (
                duplicate[-1][0],
                duplicate[0][1],
                duplicate[-1][2],
                duplicate[-1][3],
                duplicate[-1][4],
            )
            duplicate_table = _write_table(root, "duplicate.tsv", header, duplicate)
            with self.assertRaisesRegex(EndpointPowerError, "one row per donor"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=duplicate_table,
                    n_resamples=100,
                    seed=17,
                    parameters={},
                )

    def test_variant_and_graph_endpoints_resample_blocks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            variant_rows: list[tuple[object, ...]] = []
            for stratum in ("hep", "stellate"):
                for index, observed in enumerate(
                    (float(value * value + 1) for value in range(1, 21)), start=1
                ):
                    variant_rows.append(
                        (
                            _hash(f"variant-{stratum}-{index}"),
                            _hash(f"variant-unit-{stratum}-{index}"),
                            _hash(f"block-{index}"),
                            stratum,
                            observed,
                            observed + 0.01 * index,
                            -observed,
                        )
                    )
            variant = _write_table(
                root,
                "variant.tsv",
                (
                    "row_hash",
                    "unit_hash",
                    "block_hash",
                    "stratum",
                    "observed",
                    "candidate",
                    "baseline",
                ),
                variant_rows,
            )
            variant_result = recompute_endpoint_bootstrap(
                evaluator_id="variant_ld_block_fisher_z_spearman_gain_v1",
                table_path=variant,
                n_resamples=100,
                seed=29,
                parameters={"strata": ["hep", "stellate"]},
            )
            self.assertGreater(variant_result.observed_effect, 0.0)
            self.assertEqual(
                variant_result.independent_unit_counts, {"ld_block": 20}
            )

            graph_rows: list[tuple[object, ...]] = []
            row_index = 0
            for block in range(1, 5):
                for truth in (0, 1):
                    row_index += 1
                    graph_rows.append(
                        (
                            _hash(f"graph-row-{row_index}"),
                            _hash(f"graph-unit-{row_index}"),
                            _hash(f"graph-block-{block}"),
                            truth,
                            0.9 if truth else 0.1,
                            0.1 if truth else 0.9,
                        )
                    )
            graph = _write_table(
                root,
                "graph.tsv",
                (
                    "row_hash",
                    "unit_hash",
                    "block_hash",
                    "observed_binary",
                    "candidate",
                    "baseline",
                ),
                graph_rows,
            )
            graph_result = recompute_endpoint_bootstrap(
                evaluator_id="graph_ld_block_auprc_gain_v1",
                table_path=graph,
                n_resamples=100,
                seed=31,
                parameters={},
            )
            self.assertGreater(graph_result.observed_effect, 0.0)
            self.assertEqual(
                graph_result.independent_unit_counts, {"held_block": 4}
            )

    def test_rna_atac_endpoint_supports_sparse_donor_block_topology(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows: list[tuple[object, ...]] = []
            row_index = 0
            # The d10-b1 cell is deliberately absent. A two-way biological-unit
            # bootstrap must not require a fabricated complete Cartesian grid.
            donor_block_cells = [
                (f"d{donor}", f"b{block}")
                for donor in range(1, 11)
                for block in range(1, 11)
                if (donor, block) != (10, 1)
            ]
            for cell_index, (donor, block) in enumerate(
                donor_block_cells, start=1
            ):
                for stratum, offset in (("hep", 5.0), ("stellate", 1.0)):
                    # Two bins are required to define a profile deviance.
                    # The candidate reproduces their relative abundance while
                    # the cell-type-mean fixture deliberately does not.
                    for bin_index, multiplier in ((1, 1.0), (2, 2.0)):
                        row_index += 1
                        observed = (offset + cell_index) * multiplier
                        rows.append(
                            (
                                _hash(f"rna-row-{row_index}"),
                                _hash(donor),
                                _hash(block),
                                stratum,
                                observed,
                                observed,
                                5.0,
                            )
                        )
            table = _write_table(
                root,
                "rna-atac.tsv",
                (
                    "row_hash",
                    "donor_hash",
                    "block_hash",
                    "stratum",
                    "observed",
                    "candidate",
                    "baseline",
                ),
                rows,
            )
            result = recompute_endpoint_bootstrap(
                evaluator_id="rna_atac_two_way_deviance_reduction_v1",
                table_path=table,
                n_resamples=100,
                seed=37,
                parameters={"strata": ["hep", "stellate"]},
            )
            self.assertGreater(result.observed_effect, 0.0)
            self.assertEqual(
                result.independent_unit_counts,
                {"donor": 10, "genomic_block": 10},
            )

    def test_undefined_fixed_bootstrap_draw_is_not_retried(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = _write_table(
                root,
                "degenerate.tsv",
                ("row_hash", "unit_hash", "observed", "candidate", "baseline"),
                [
                    (_hash(f"r{i}"), _hash(f"d{i}"), i, i, -i)
                    for i in range(1, 4)
                ],
            )
            with self.assertRaisesRegex(EndpointPowerError, "may not be retried"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=table,
                    n_resamples=100,
                    seed=17,
                    parameters={},
                )

    def test_schema_and_parameter_drift_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = _write_table(
                root,
                "bad.tsv",
                ("unit_hash", "row_hash", "observed", "candidate", "baseline"),
                [(_hash("d1"), _hash("r1"), 1, 1, 0)],
            )
            with self.assertRaisesRegex(EndpointPowerError, "header must be exactly"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=table,
                    n_resamples=100,
                    seed=1,
                    parameters={},
                )
            extra = _write_table(
                root,
                "extra.tsv",
                ("row_hash", "unit_hash", "observed", "candidate", "baseline"),
                [(_hash("r1"), _hash("d1"), 1, 1, 0, "undeclared")],
            )
            with self.assertRaisesRegex(EndpointPowerError, "exactly the registered"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=extra,
                    n_resamples=100,
                    seed=1,
                    parameters={},
                )
            with self.assertRaisesRegex(EndpointPowerError, "parameters"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=table,
                    n_resamples=100,
                    seed=1,
                    parameters={"caller_summary": True},
                )

    def test_cache_is_content_addressed_copy_safe_and_invalidates_on_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            header = ("row_hash", "unit_hash", "observed", "candidate", "baseline")
            rows = [
                (_hash(f"cache-r{i}"), _hash(f"cache-d{i}"), i * i, i * i, -i)
                for i in range(1, 21)
            ]
            table = _write_table(root, "cache.tsv", header, rows)
            evaluator_id = "bulk_paired_spearman_gain_v1"
            original = endpoint_power._EVALUATORS[evaluator_id]
            calls = 0

            def counted(path, n_resamples, seed, parameters):
                nonlocal calls
                calls += 1
                return original(path, n_resamples, seed, parameters)

            endpoint_power._clear_endpoint_cache_for_testing()
            with patch.dict(
                endpoint_power._EVALUATORS, {evaluator_id: counted}, clear=False
            ):
                first = recompute_endpoint_bootstrap(
                    evaluator_id=evaluator_id,
                    table_path=table,
                    n_resamples=100,
                    seed=53,
                    parameters={},
                )
                second = recompute_endpoint_bootstrap(
                    evaluator_id=evaluator_id,
                    table_path=table,
                    n_resamples=100,
                    seed=53,
                    parameters={},
                )
                self.assertEqual(calls, 1)
                self.assertEqual(first, second)
                first.independent_unit_counts["tampered"] = 999  # type: ignore[index]
                third = recompute_endpoint_bootstrap(
                    evaluator_id=evaluator_id,
                    table_path=table,
                    n_resamples=100,
                    seed=53,
                    parameters={},
                )
                self.assertEqual(calls, 1)
                self.assertNotIn("tampered", third.independent_unit_counts)

                changed_rows = list(rows)
                changed_rows[-1] = (*changed_rows[-1][:-2], 9999, -20)
                _write_table(root, "cache.tsv", header, changed_rows)
                recompute_endpoint_bootstrap(
                    evaluator_id=evaluator_id,
                    table_path=table,
                    n_resamples=100,
                    seed=53,
                    parameters={},
                )
                self.assertEqual(calls, 2)

    def test_cache_boundary_rejects_symlinks_and_evaluator_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            header = ("row_hash", "unit_hash", "observed", "candidate", "baseline")
            rows = [
                (_hash(f"boundary-r{i}"), _hash(f"boundary-d{i}"), i * i, i * i, -i)
                for i in range(1, 21)
            ]
            table = _write_table(root, "table.tsv", header, rows)
            direct_link = root / "table-link.tsv"
            os.symlink(table, direct_link)
            with self.assertRaisesRegex(EndpointPowerError, "symlink"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=direct_link,
                    n_resamples=100,
                    seed=59,
                    parameters={},
                )

            real_parent = root / "real-parent"
            real_parent.mkdir()
            nested = _write_table(real_parent, "nested.tsv", header, rows)
            parent_link = root / "parent-link"
            os.symlink(real_parent, parent_link)
            with self.assertRaisesRegex(EndpointPowerError, "symlink"):
                recompute_endpoint_bootstrap(
                    evaluator_id="bulk_paired_spearman_gain_v1",
                    table_path=parent_link / nested.name,
                    n_resamples=100,
                    seed=59,
                    parameters={},
                )

            evaluator_id = "bulk_paired_spearman_gain_v1"
            original = endpoint_power._EVALUATORS[evaluator_id]

            def mutating(path, n_resamples, seed, parameters):
                result = original(path, n_resamples, seed, parameters)
                path.write_text(path.read_text(encoding="utf-8") + "# changed\n")
                return result

            endpoint_power._clear_endpoint_cache_for_testing()
            with patch.dict(
                endpoint_power._EVALUATORS, {evaluator_id: mutating}, clear=False
            ):
                with self.assertRaisesRegex(EndpointPowerError, "evaluator was running"):
                    recompute_endpoint_bootstrap(
                        evaluator_id=evaluator_id,
                        table_path=table,
                        n_resamples=100,
                        seed=61,
                        parameters={},
                    )


if __name__ == "__main__":
    unittest.main()
