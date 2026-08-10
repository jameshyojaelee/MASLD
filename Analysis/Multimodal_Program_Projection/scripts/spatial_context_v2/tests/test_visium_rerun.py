#!/usr/bin/env python3
"""Lightweight semantic fixtures for the compatible-Visium candidate lane."""

from __future__ import annotations

import hashlib
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from visium_rerun_lib import (  # noqa: E402
    N_NULL,
    V1_ENGINE_SHA256,
    build_paths,
    load_legacy_engine,
    prepare_engine_inputs,
    score_variants,
    sha256_file,
    universe_hash,
    v1_universe,
    v2_universe,
)


class VisiumRerunSemantics(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.paths = build_paths()
        cls.engine = load_legacy_engine(cls.paths)

    def test_pinned_engine_and_frozen_families(self) -> None:
        self.assertEqual(sha256_file(self.paths.v1_engine), V1_ENGINE_SHA256)
        v1 = v1_universe(self.paths)
        v2 = v2_universe(self.paths)
        self.assertEqual(len(v1), 22)
        self.assertEqual(len(v2), 2)
        self.assertNotEqual(universe_hash(v1), universe_hash(v2))
        self.assertEqual(
            {row["legacy_program_id"] for row in v2},
            {"hepatocytes::8", "hepatocytes::20"},
        )

    def test_v2_registry_translation_is_deterministic(self) -> None:
        registry_a, membership_a, universe_a = prepare_engine_inputs(self.paths, "v2")
        registry_b, membership_b, universe_b = prepare_engine_inputs(self.paths, "v2")
        pd.testing.assert_frame_equal(registry_a, registry_b)
        pd.testing.assert_frame_equal(membership_a, membership_b)
        self.assertEqual(universe_hash(universe_a), universe_hash(universe_b))
        self.assertEqual(registry_a["program_id"].tolist(), [row["program_id"] for row in universe_a])
        # Four source genes are retained in the immutable module definitions
        # but are not uniquely mappable in GENCODE v49, so 67/71 rows enter an
        # external assay. Their original L1 mass remains in the denominator.
        self.assertEqual(len(membership_a), 67)
        self.assertTrue(membership_a["mapped_symbol"].all())

    def test_tissue_islands_exclude_small_components(self) -> None:
        grid = np.array([(x, y) for x in range(3) for y in range(3)], dtype=float)
        coords = np.vstack([grid, grid + np.array([100.0, 0.0]), np.array([[200, 0], [201, 0], [202, 0]])])
        obs = pd.DataFrame({"sample_id": ["array"] * len(coords)})
        graphs, audit = self.engine.section_graphs(obs, coords, k=6)
        self.assertEqual(len(graphs["array"]), 2)
        self.assertEqual(int(audit.loc[0, "n_tissue_islands"]), 2)
        self.assertEqual(int(audit.loc[0, "n_spots_in_graph"]), 18)
        self.assertEqual(int(audit.loc[0, "n_spots_excluded_small_islands"]), 3)
        self.assertLessEqual(
            float(audit.loc[0, "max_observed_edge_distance"]),
            float(audit.loc[0, "max_allowed_edge_distance"]),
        )

    def test_moran_matches_hand_calculation(self) -> None:
        graphs = {
            "array": [
                (
                    np.array([0, 1, 2, 3]),
                    np.array([0, 1, 2]),
                    np.array([1, 2, 3]),
                )
            ]
        }
        observed = self.engine.moran_columns(np.array([0.0, 1.0, 2.0, 3.0]), graphs)["array"][0]
        centered = np.array([-1.5, -0.5, 0.5, 1.5])
        expected = 4 * np.sum(centered[:-1] * centered[1:]) / (3 * np.sum(centered**2))
        self.assertAlmostEqual(float(observed), float(expected), places=7)

    def test_donor_collapse_is_not_section_pseudoreplication(self) -> None:
        values = {"H35a": np.array([1.0]), "H35b": np.array([3.0]), "D2": np.array([8.0])}
        obs = pd.DataFrame(
            {
                "sample_id": ["H35a", "H35b", "D2"],
                "individual": ["H35", "H35", "D2"],
            }
        )
        observed = self.engine.donor_collapse(values, obs, "GSE192741")[0]
        self.assertEqual(observed, 5.0)  # mean(mean(1,3), 8), not mean(1,3,8)
        vu_obs = pd.DataFrame({"sample_id": ["A", "B"], "individual": ["A", "B"]})
        vu = self.engine.donor_collapse({"A": np.array([2.0]), "B": np.array([6.0])}, vu_obs, "Vu_et_al_2025")[0]
        self.assertEqual(vu, 4.0)

    def test_matched_controls_exclude_complete_program_union(self) -> None:
        genes = [f"G{i:02d}" for i in range(40)]
        pool = pd.DataFrame(
            {
                "gene": genes,
                "biotype": ["protein_coding"] * 40,
                "mt": [False] * 40,
                "ribo": [False] * 40,
                "mean_bin": [5] * 40,
                "detect_bin": [5] * 40,
                "lineage_corr_bin": [5] * 40,
                "lineage_corr": np.linspace(-0.2, 0.2, 40),
            }
        )
        measured = pd.DataFrame({"gene_symbol": ["G00", "G01"]})
        union = {"G00", "G01", "G02", "G03"}
        first, expr_relax, lineage_relax = self.engine.build_matched_sets(
            measured, pool, union, 25, np.random.default_rng(42)
        )
        second, _, _ = self.engine.build_matched_sets(
            measured, pool, union, 25, np.random.default_rng(42)
        )
        np.testing.assert_array_equal(first, second)
        selected = set(pool.iloc[np.unique(first)]["gene"])
        self.assertFalse(selected.intersection(union))
        self.assertEqual(expr_relax, [0, 0])
        self.assertEqual(lineage_relax, [0, 0])
        self.assertEqual(hashlib.sha256(first.tobytes()).hexdigest(), hashlib.sha256(second.tobytes()).hexdigest())

    def test_weighted_equal_and_leave_top_score_fixtures(self) -> None:
        z = np.array([[1.0, 2.0, 4.0], [4.0, 1.0, -1.0]])
        primary, equal, leave, top, leave_weights = score_variants(z, np.array([0.6, 0.3, 0.1]))
        np.testing.assert_allclose(primary, np.array([1.6, 2.6]))
        np.testing.assert_allclose(equal, np.mean(z, axis=1))
        self.assertEqual(top, 0)
        np.testing.assert_allclose(leave_weights, np.array([0.75, 0.25]))
        np.testing.assert_allclose(leave, z[:, 1:] @ np.array([0.75, 0.25]))

    def test_primary_production_null_is_fixed(self) -> None:
        self.assertEqual(N_NULL, 9_999)
        self.assertFalse(self.engine.SMOKE)
        self.assertEqual(self.engine.N_NULL, 9_999)
        self.assertEqual(self.engine.N_SENSITIVITY_NULL, 999)


if __name__ == "__main__":
    unittest.main(verbosity=2)
