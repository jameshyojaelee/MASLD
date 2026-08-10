#!/usr/bin/env python3
"""Synthetic, outcome-free tests for the Plan 11 statistical machinery."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import anndata as ad
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from yakubovsky_common import (  # noqa: E402
    RELEASE_ID,
    ContractError,
    bh_adjust,
    block_bootstrap_lipid_slope,
    bulk_decode_matlab_cellstr,
    build_model_design,
    build_spatial_blocks,
    gene_axis_sha256,
    ols_lipid_fit,
    read_h5ad_dataframe_index,
    score_frozen_program,
    signed_stouffer,
    sha256_file,
    validate_source_gate_summary,
)
from importlib import import_module  # noqa: E402

model_module = import_module("03_fit_binary_lipid")
scoring_module = import_module("02_score_programs")
adapter_module = import_module("01_build_adapter")
freeze_module = import_module("00_freeze_plan11")


class Plan11CoreTests(unittest.TestCase):

    def test_named_anndata_dataframe_index_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "named_index.h5ad"
            adata = ad.AnnData(
                X=np.ones((2, 2), dtype=float),
                obs=pd.DataFrame(index=pd.Index(["S1", "S2"], name="spot_id")),
                var=pd.DataFrame(index=pd.Index(["G1", "β2"], name="gene_symbol")),
            )
            adata.write_h5ad(path)
            with h5py.File(path, "r") as handle:
                self.assertEqual(handle["var"].attrs["_index"], "gene_symbol")
                self.assertNotIn("_index", handle["var"])
                np.testing.assert_array_equal(
                    read_h5ad_dataframe_index(handle["var"]),
                    np.asarray(["G1", "β2"], dtype=object),
                )

    @staticmethod
    def _write_cellstr(
        handle: h5py.File,
        axis_name: str,
        strings: list[str],
        *,
        code_units: list[list[int]] | None = None,
    ) -> h5py.Dataset:
        references = np.empty((1, len(strings)), dtype=h5py.ref_dtype)
        for index, value in enumerate(strings):
            units = code_units[index] if code_units is not None else [ord(x) for x in value]
            target = handle.create_dataset(
                f"chars/{axis_name}_{index}",
                data=np.asarray(units, dtype=np.uint16).reshape(-1, 1),
            )
            references[0, index] = target.ref
        axis = handle.create_dataset(axis_name, shape=references.shape, dtype=h5py.ref_dtype)
        axis[...] = references
        return axis

    def test_lipid_applicability_vocabulary_serializes_to_h5ad(self) -> None:
        import anndata as ad
        from scipy import sparse

        values = [
            adapter_module.h5ad_lipid_zone_value("lipid_zone", True),
            adapter_module.h5ad_lipid_zone_value("non_lipid_zone", True),
            adapter_module.h5ad_lipid_zone_value("not_applicable", False),
        ]
        self.assertEqual(values, ["true", "false", "not_applicable"])
        self.assertTrue(all(isinstance(value, str) for value in values))
        obs = pd.DataFrame(
            {"lipid_zone": values}, index=pd.Index(["S1", "S2", "S3"])
        )
        var = pd.DataFrame(index=pd.Index(["G1"]))
        fixture = ad.AnnData(X=None, obs=obs, var=var, shape=(3, 1))
        fixture.layers["source_background_corrected_abundance"] = sparse.csr_matrix(
            np.asarray([[1.0], [2.0], [0.0]])
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "three_state_obs.h5ad"
            fixture.write_h5ad(path)
            observed = ad.read_h5ad(path).obs["lipid_zone"].astype(str).tolist()
        self.assertEqual(observed, values)

    def test_p17_like_reference_only_spot_is_retained_but_never_eligible(self) -> None:
        self.assertFalse(
            adapter_module.source_spot_analysis_eligible(False, 125.0)
        )
        self.assertFalse(
            adapter_module.source_spot_analysis_eligible(True, 0.0)
        )
        self.assertTrue(
            adapter_module.source_spot_analysis_eligible(True, 125.0)
        )
        self.assertEqual(
            adapter_module.source_spot_binary_lipid_state(
                True, False, False, ""
            ),
            (False, "not_applicable"),
        )
        self.assertEqual(
            adapter_module.source_spot_binary_lipid_state(
                True, True, True, "lipid_zone"
            ),
            (True, "lipid_zone"),
        )
        with self.assertRaises(ContractError):
            adapter_module.source_spot_binary_lipid_state(
                True, True, True, "invalid"
            )

    def test_binary_source_gate_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            gate = Path(directory) / "gate.tsv"
            pd.DataFrame(
                [
                    {
                        "release_id": RELEASE_ID,
                        "dataset": "yakubovsky2026",
                        "status": "pass_ordinal_lipid",
                        "terminal_verdict": "pass_ordinal_lipid",
                        "source_gate_pass": True,
                        "integrity_gate_pass": True,
                        "zonation_gate_pass": True,
                        "continuous_lipid_gate_pass": False,
                        "ordinal_lipid_gate_pass": True,
                        "authoritative_barcode_join_verified": True,
                        "authoritative_ordinal_lipid_join_verified": True,
                        # This is the distinct continuous-lipid donor gate and
                        # must not invalidate a source-authorized ordinal branch.
                        "minimum_donor_thresholds_pass": False,
                        "n_passing_ordinal_donors": 3,
                        "passing_ordinal_donors": "D1;D2;D3",
                        "ordinal_donor_spot_counts_lipid_nonlipid_total": (
                            "D1:40/60/100;D2:55/65/120;D3:70/80/150"
                        ),
                        "ordinal_effect_interpretation": "binary_lipid_zone_vs_non_lipid_zone",
                        "v2_registry_read": False,
                    }
                ]
            ).to_csv(gate, sep="\t", index=False)
            contract = validate_source_gate_summary(gate)
            self.assertEqual(contract.passing_donors, ("D1", "D2", "D3"))

    def test_bulk_matlab_cellstr_resolver_and_axis_digest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cellstr.h5"
            with h5py.File(path, "w") as handle:
                first = self._write_cellstr(handle, "axis_a", ["G1", "G2", "G3"])
                second = self._write_cellstr(handle, "axis_b", ["G1", "G2", "G3"])
                zero_normalized = self._write_cellstr(
                    handle,
                    "axis_zero",
                    ["G1", "G2", "G3"],
                    code_units=[[71, 49, 0], [71, 0, 50], [71, 51]],
                )
                resolved = bulk_decode_matlab_cellstr(
                    {"a": first, "b": second, "z": zero_normalized}, handle
                )
                self.assertEqual(resolved.values["a"], ["G1", "G2", "G3"])
                self.assertEqual(resolved.values["a"], resolved.values["b"])
                self.assertEqual(resolved.values["a"], resolved.values["z"])
                self.assertEqual(resolved.n_references, 9)
                self.assertEqual(resolved.n_unique_objects, 9)
                base_hash = gene_axis_sha256(resolved.values["a"])
                self.assertEqual(base_hash, gene_axis_sha256(resolved.values["b"]))
                self.assertNotEqual(base_hash, gene_axis_sha256(["G1", "GX", "G3"]))
                self.assertNotEqual(base_hash, gene_axis_sha256(["G2", "G1", "G3"]))
                self.assertNotEqual(base_hash, gene_axis_sha256(["G1", "G2"]))

            with h5py.File(path, "a") as handle:
                null_axis = handle.create_dataset(
                    "null_axis", shape=(1, 1), dtype=h5py.ref_dtype
                )
                null_axis[0, 0] = h5py.Reference()
                with self.assertRaises(ContractError):
                    bulk_decode_matlab_cellstr({"null": null_axis}, handle)

    def test_gene_axis_verification_artifact_and_wrapper_order(self) -> None:
        verifier = SCRIPT_DIR / "07_verify_gene_axis.py"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gene_axis_verification.tsv"
            pd.DataFrame(
                [
                    {
                        "release_id": RELEASE_ID,
                        "status": "exact_all_donor_gene_axes_verified",
                        "n_donors": 16,
                        "n_genes": 36601,
                        "common_gene_axis_sha256": "a" * 64,
                        "v_mat_sha256": "b" * 64,
                        "manifest_sha256": "c" * 64,
                        "producer": str(verifier.resolve()),
                        "producer_sha256": sha256_file(verifier),
                        "matrix_values_read": False,
                        "lipid_fields_read": False,
                    }
                ]
            ).to_csv(path, sep="\t", index=False)
            observed = freeze_module.validate_gene_axis_verification(
                path,
                script_dir=SCRIPT_DIR,
                v_mat_sha256="b" * 64,
                manifest_sha256="c" * 64,
                common_sha256="a" * 64,
                n_donors=16,
                n_genes=36601,
            )
            self.assertEqual(observed, sha256_file(path))
            tampered = pd.read_csv(path, sep="\t")
            tampered.loc[0, "matrix_values_read"] = True
            tampered.to_csv(path, sep="\t", index=False)
            with self.assertRaises(ContractError):
                freeze_module.validate_gene_axis_verification(
                    path,
                    script_dir=SCRIPT_DIR,
                    v_mat_sha256="b" * 64,
                    manifest_sha256="c" * 64,
                    common_sha256="a" * 64,
                    n_donors=16,
                    n_genes=36601,
                )
        wrapper = (SCRIPT_DIR / "run_adapter.sbatch").read_text()
        self.assertLess(
            wrapper.index("07_verify_gene_axis.py"),
            wrapper.index("00_freeze_plan11.py"),
        )
        self.assertIn("--write --output", wrapper)
        self.assertIn("--check --output", wrapper)

    def test_chunk_aligned_reader_identity_both_orientations(self) -> None:
        values = np.asarray(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [1.5, 2.0, 3.25, 4.0, 5.0],
                [6.0, 0.0, 7.0, 8.5, 9.0],
                [10.0, 11.0, 0.0, 12.0, 13.0],
                [14.0, 15.0, 16.0, 0.0, 17.0],
                [18.0, 19.0, 20.0, 21.0, 0.0],
                [22.0, 23.0, 24.0, 25.0, 26.0],
            ],
            dtype=float,
        )
        genes = ["A", "B", "C", "D", "E"]
        mapping = {"A": "G1", "B": "G2", "C": "G1", "D": None, "E": "G3"}
        global_index = {"G1": 0, "G2": 1, "G3": 2}
        expected = np.column_stack(
            [values[:, 0] + values[:, 2], values[:, 1], values[:, 4]]
        )
        observed = []
        with tempfile.TemporaryDirectory() as directory:
            for label, payload, chunks in (
                ("spot_gene", values, (3, 2)),
                ("gene_spot", values.T, (2, 3)),
            ):
                path = Path(directory) / f"{label}.h5"
                with h5py.File(path, "w") as handle:
                    dataset = handle.create_dataset(
                        "mat", data=payload, chunks=chunks, compression="gzip"
                    )
                    result = adapter_module.stream_donor_abundance(
                        dataset,
                        len(values),
                        genes,
                        mapping,
                        global_index,
                        np.asarray([1, 4]),
                        dense_tile_mib=0.0001,
                    )
                abundance, landmarks, source_sum, detected, audit = result
                np.testing.assert_allclose(abundance.toarray(), expected, atol=0, rtol=0)
                np.testing.assert_allclose(
                    source_sum,
                    np.sum(values, axis=1, dtype=np.longdouble).astype(float),
                    atol=1e-12,
                    rtol=0,
                )
                np.testing.assert_array_equal(detected, np.count_nonzero(values, axis=1))
                np.testing.assert_allclose(landmarks, values[:, 1] + values[:, 4])
                self.assertEqual(audit["physical_chunk_coverage_min"], 1)
                self.assertEqual(audit["physical_chunk_coverage_max"], 1)
                self.assertEqual(
                    audit["optimized_estimated_chunk_touches"],
                    audit["n_physical_hdf5_chunks"],
                )
                observed.append(abundance)
            self.assertEqual(observed[0].shape, observed[1].shape)
            self.assertEqual((observed[0] != observed[1]).nnz, 0)

    def test_chunk_reader_edge_contracts_and_contiguous_fallback(self) -> None:
        values = np.asarray([[0.0, 1.0, 2.0], [3.0, 0.0, 4.5]], dtype=float)
        genes = ["A", "B", "C"]
        with tempfile.TemporaryDirectory() as directory:
            contiguous_path = Path(directory) / "contiguous.h5"
            with h5py.File(contiguous_path, "w") as handle:
                dataset = handle.create_dataset("mat", data=values)
                abundance, _, source_sum, detected, audit = (
                    adapter_module.stream_donor_abundance(
                        dataset,
                        2,
                        genes,
                        {gene: None for gene in genes},
                        {},
                        np.asarray([], dtype=int),
                        dense_tile_mib=0.00001,
                    )
                )
                self.assertEqual(abundance.shape, (2, 0))
                np.testing.assert_allclose(source_sum, values.sum(axis=1))
                np.testing.assert_array_equal(detected, [2, 2])
                self.assertEqual(audit["matrix_traversal"], "contiguous_bounded_fallback")

            square_path = Path(directory) / "square.h5"
            with h5py.File(square_path, "w") as handle:
                square = handle.create_dataset("mat", data=np.eye(3), chunks=(2, 2))
                with self.assertRaises(ContractError):
                    adapter_module.stream_donor_abundance(
                        square,
                        3,
                        genes,
                        {gene: gene for gene in genes},
                        {gene: index for index, gene in enumerate(genes)},
                        np.asarray([], dtype=int),
                        dense_tile_mib=1,
                    )

            for name, bad in (
                ("nan", np.asarray([[1.0, np.nan, 0.0], [0.0, 1.0, 2.0]])),
                ("inf", np.asarray([[1.0, np.inf, 0.0], [0.0, 1.0, 2.0]])),
                ("negative", np.asarray([[1.0, -0.1, 0.0], [0.0, 1.0, 2.0]])),
            ):
                bad_path = Path(directory) / f"{name}.h5"
                with h5py.File(bad_path, "w") as handle:
                    dataset = handle.create_dataset("mat", data=bad, chunks=(2, 2))
                    with self.assertRaises(ContractError):
                        adapter_module.stream_donor_abundance(
                            dataset,
                            2,
                            genes,
                            {gene: gene for gene in genes},
                            {gene: index for index, gene in enumerate(genes)},
                            np.asarray([], dtype=int),
                            dense_tile_mib=1e-9,
                        )

    def test_recorded_selections_cover_each_hdf5_chunk_once(self) -> None:
        class RecordingDataset:
            def __init__(self, values: np.ndarray, chunks: tuple[int, int]):
                self.values = values
                self.shape = values.shape
                self.dtype = values.dtype
                self.chunks = chunks
                self.selections: list[tuple[slice, slice]] = []

            def read_direct(self, destination, source_sel=None):
                self.selections.append(source_sel)
                destination[...] = self.values[source_sel]

        values = np.arange(11 * 7, dtype=float).reshape(11, 7)
        dataset = RecordingDataset(values, (4, 3))
        genes = [f"G{index}" for index in range(7)]
        result = adapter_module.stream_donor_abundance(
            dataset,
            11,
            genes,
            {gene: gene for gene in genes},
            {gene: index for index, gene in enumerate(genes)},
            np.asarray([0, 6]),
            dense_tile_mib=1e-9,
        )
        np.testing.assert_allclose(result[0].toarray(), values)
        coverage = np.zeros((3, 3), dtype=int)
        for spot_slice, gene_slice in dataset.selections:
            self.assertEqual(spot_slice.start % 4, 0)
            self.assertEqual(gene_slice.start % 3, 0)
            spot_stop_chunk = int(np.ceil(spot_slice.stop / 4))
            gene_stop_chunk = int(np.ceil(gene_slice.stop / 3))
            coverage[
                spot_slice.start // 4 : spot_stop_chunk,
                gene_slice.start // 3 : gene_stop_chunk,
            ] += 1
        np.testing.assert_array_equal(coverage, np.ones_like(coverage))
        self.assertTrue(result[4]["minimum_atomic_chunk_promoted_over_budget"])

    def test_three_program_two_donor_score_fixture(self) -> None:
        rng = np.random.default_rng(11)
        weights = [
            np.arange(1, 11, dtype=float),
            np.arange(10, 0, -1, dtype=float),
            np.repeat(1.0, 10),
        ]
        comparisons = 0
        for _donor in ("D1", "D2"):
            values = rng.normal(size=(120, 10))
            values = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
            for program_weights in weights:
                primary, equal, leave, top = score_frozen_program(values, program_weights)
                np.testing.assert_allclose(
                    primary,
                    values @ (program_weights / np.sum(program_weights)),
                    rtol=0,
                    atol=1e-12,
                )
                np.testing.assert_allclose(equal, values.mean(axis=1), rtol=0, atol=1e-12)
                keep = np.arange(10) != top
                np.testing.assert_allclose(
                    leave,
                    values[:, keep]
                    @ (program_weights[keep] / np.sum(program_weights[keep])),
                    rtol=0,
                    atol=1e-12,
                )
                comparisons += 3
        self.assertEqual(comparisons, 18)

    def test_spatial_blocks_and_bootstrap(self) -> None:
        grid_x, grid_y = np.meshgrid(np.arange(25), np.arange(25))
        x = grid_x.ravel().astype(float)
        y = grid_y.ravel().astype(float)
        section = np.repeat("S1", len(x))
        built = build_spatial_blocks("D1", section, x, y)
        self.assertGreaterEqual(built.n_blocks, 8)
        active = built.audit[built.audit["block_id"] >= 0]
        self.assertTrue((active["n_spots"] >= 5).all())
        self.assertEqual(active.groupby("block_id")["section"].nunique().max(), 1)

        zonation = x / x.max()
        lipid = ((grid_x.ravel() // 5 + grid_y.ravel() // 5) % 2).astype(float)
        library = 10_000.0 + 7 * x + 11 * y
        detected = 1_000 + y
        design, _ = build_model_design(lipid, zonation, library, detected, 4)
        coefficient = np.linspace(-0.2, 0.3, design.shape[1])
        coefficient[1] = 0.8
        outcome = design @ coefficient + np.random.default_rng(7).normal(0, 0.05, len(x))
        fit = ols_lipid_fit(design, outcome)
        self.assertAlmostEqual(fit.beta, 0.8, delta=0.02)
        bootstrap = block_bootstrap_lipid_slope(
            design,
            outcome,
            built.labels,
            section,
            n_bootstrap=99,
            seed=123,
        )
        self.assertLessEqual(bootstrap.n_failed, 4)
        self.assertGreater(bootstrap.ci_low, 0)

    def test_block_level_label_shuffle_is_centered(self) -> None:
        grid_x, grid_y = np.meshgrid(np.arange(25), np.arange(25))
        x = grid_x.ravel().astype(float)
        y = grid_y.ravel().astype(float)
        section = np.repeat("S1", len(x))
        built = build_spatial_blocks("D1", section, x, y)
        zonation = x / x.max()
        library = 9_000.0 + 5 * x + 13 * y
        detected = 1_100.0 + ((3 * x + 7 * y) % 29)
        rng = np.random.default_rng(99)
        outcome = 1.5 * zonation + rng.normal(0, 0.4, len(x))
        block_ids = np.unique(built.labels[built.labels >= 0])
        base_labels = np.asarray([index % 2 for index in range(len(block_ids))])
        estimates = []
        for _ in range(200):
            shuffled = rng.permutation(base_labels)
            label_map = dict(zip(block_ids, shuffled, strict=True))
            lipid = np.asarray([label_map[block] for block in built.labels], dtype=float)
            design, _ = build_model_design(lipid, zonation, library, detected, 4)
            estimates.append(ols_lipid_fit(design, outcome).beta)
        self.assertLess(abs(float(np.mean(estimates))), 0.08)

    def test_bootstrap_tail_probability_uses_only_valid_draws(self) -> None:
        # The exposure is carried by one of eight blocks, so draws that omit
        # that block are singular by construction.  Those failed draws must be
        # reported but cannot silently enlarge the tail-probability denominator.
        blocks = np.repeat(np.arange(8), 10)
        exposure = (blocks == 0).astype(float)
        design = np.column_stack([np.ones(len(blocks)), exposure])
        outcome = 1.2 * exposure + np.random.default_rng(13).normal(0, 0.1, len(blocks))
        result = block_bootstrap_lipid_slope(
            design,
            outcome,
            blocks,
            np.repeat("S1", len(blocks)),
            n_bootstrap=999,
            seed=314,
        )
        valid = result.estimates[np.isfinite(result.estimates)]
        self.assertGreater(result.n_failed, 0)
        expected = min(
            1.0,
            2
            * (1 + min(int(np.sum(valid <= 0)), int(np.sum(valid >= 0))))
            / (len(valid) + 1),
        )
        self.assertAlmostEqual(result.pvalue, expected, places=15)

    def test_stouffer_and_bh(self) -> None:
        slopes = np.asarray([0.3, 0.2, 0.4])
        pvalues = np.asarray([0.05, 0.1, 0.01])
        z, p = signed_stouffer(slopes, pvalues)
        expected_z = np.sum(__import__("scipy").stats.norm.isf(pvalues / 2)) / np.sqrt(3)
        self.assertAlmostEqual(z, expected_z, places=12)
        self.assertAlmostEqual(p, 2 * __import__("scipy").stats.norm.sf(abs(expected_z)), places=12)
        np.testing.assert_allclose(
            bh_adjust([0.01, 0.04, 0.03]),
            [0.03, 0.04, 0.04],
            atol=1e-12,
        )

    def test_all_null_is_valid(self) -> None:
        adjusted = bh_adjust([0.7, 0.9])
        self.assertTrue((adjusted >= 0.7).all())
        self.assertEqual(int((adjusted < 0.05).sum()), 0)

    def test_zero_valid_donor_tables_are_schema_safe(self) -> None:
        empty = pd.DataFrame(columns=model_module.DONOR_EFFECT_COLUMNS)
        empty["estimable"] = empty["estimable"].astype(bool)
        summary = model_module.combine_variant(empty)
        self.assertEqual(summary["n_valid_donors"], 0)
        self.assertTrue(np.isnan(summary["combined_z"]))
        q_value, q_df, q_p = model_module.heterogeneity(empty)
        self.assertTrue(np.isnan(q_value))
        self.assertEqual(q_df, 0)
        self.assertTrue(np.isnan(q_p))

    def test_zero_testable_and_zero_eligible_scoring_tables_are_schema_safe(
        self,
    ) -> None:
        for branch in ("zero_testable_programs", "zero_eligible_spots"):
            with self.subTest(branch=branch), tempfile.TemporaryDirectory() as directory:
                frames = scoring_module.assemble_scoring_output_frames([], [], [], [])
                expected_columns = (
                    scoring_module.SCORE_COLUMNS,
                    scoring_module.DONOR_SCORE_COLUMNS,
                    scoring_module.ZONATION_REFERENCE_COLUMNS,
                    scoring_module.SCORING_AUDIT_COLUMNS,
                )
                filenames = (
                    "per_sample_program_scores.tsv",
                    "per_donor_program_scores.tsv",
                    "zonation_reference.tsv",
                    "scoring_audit.tsv",
                )
                for frame, columns, filename in zip(
                    frames, expected_columns, filenames, strict=True
                ):
                    self.assertTrue(frame.empty)
                    self.assertEqual(tuple(frame.columns), columns)
                    path = Path(directory) / filename
                    scoring_module.atomic_write_frame(path, frame)
                    reread = pd.read_csv(
                        path, sep="\t", dtype=str, keep_default_na=False
                    )
                    self.assertTrue(reread.empty)
                    self.assertEqual(tuple(reread.columns), columns)
                self.assertFalse(
                    frames[0].duplicated(["spot_id", "program_uid"]).any()
                )
                model_audit = pd.DataFrame(
                    [], columns=model_module.MODEL_DESIGN_AUDIT_COLUMNS
                )
                model_audit_path = Path(directory) / "model_design_audit.tsv"
                scoring_module.atomic_write_frame(model_audit_path, model_audit)
                reread_model_audit = pd.read_csv(
                    model_audit_path, sep="\t", dtype=str, keep_default_na=False
                )
                self.assertTrue(reread_model_audit.empty)
                self.assertEqual(
                    tuple(reread_model_audit.columns),
                    model_module.MODEL_DESIGN_AUDIT_COLUMNS,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
