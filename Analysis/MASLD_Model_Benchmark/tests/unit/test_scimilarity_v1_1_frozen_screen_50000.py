from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts import build_scimilarity_v1_1_frozen_screen_50000 as builder
from scripts import extract_scimilarity_v1_1_frozen_screen_50000 as extractor


class SCimilarityV11FrozenScreen50000Test(unittest.TestCase):
    def _inputs(self, root: Path, *, forbidden: bool = False, bad_join: bool = False):
        import anndata
        import numpy as np
        import pandas as pd
        from scipy import sparse

        row_ids = [f"cell-{index}" for index in range(5)]
        studies = list(builder.EXPOSURE_BY_STUDY)[:5]
        obs_values: dict[str, object] = {
            "row_id": row_ids,
            "donor_id": [f"donor-{index}" for index in range(5)],
            "dataset_id": studies,
            "outer_fold": list(range(5)),
            "assay": ["unknown"] * 5,
        }
        if forbidden:
            obs_values["broad_label"] = ["hepatocyte"] * 5
        obs = pd.DataFrame(obs_values, index=pd.Index(row_ids, name="row_id_index"))
        ensembl = ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]
        matrix = sparse.csr_matrix(
            np.asarray(
                [[10, 5, 0, 5], [0, 2, 8, 0], [4, 0, 1, 0], [1, 1, 1, 1], [0, 3, 2, 5]],
                dtype=np.float32,
            )
        )
        source = root / "unlabeled.h5ad"
        anndata.AnnData(
            X=matrix,
            obs=obs,
            var=pd.DataFrame(index=pd.Index(ensembl, name="ensembl_id")),
        ).write_h5ad(source)
        symbol_source = root / "symbol_source.h5ad"
        symbol_obs = pd.DataFrame(
            {"broad_label": ["must-not-be-read"] * 5}, index=pd.Index(row_ids)
        )
        symbol_var = pd.DataFrame(
            {
                "ensembl_id": ensembl,
                "source_feature_id": ["A", "B", "A", "C"],
            },
            index=pd.Index(ensembl, name="ensembl_id_index"),
        )
        anndata.AnnData(X=matrix, obs=symbol_obs, var=symbol_var).write_h5ad(symbol_source)
        split = root / "split.tsv"
        with split.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["row_id", "donor_id", "dataset", "outer_fold"])
            for index, (row_id, study) in enumerate(zip(row_ids, studies, strict=True)):
                donor = "wrong" if bad_join and index == 3 else f"donor-{index}"
                writer.writerow([row_id, donor, study, index])
        gene_order = root / "gene_order.tsv"
        gene_order.write_text("A\nB\nC\nD\n", encoding="utf-8")
        library_root = root / "libraries"
        library_root.mkdir()
        library_manifest = root / "library_manifest.tsv"
        with library_manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["dataset", "library_id", "analysis_eligible"])
            for index, study in enumerate(studies):
                library_id = f"library-{index}"
                writer.writerow([study, library_id, "True"])
                anndata.AnnData(
                    X=sparse.csr_matrix((1, 4), dtype=np.float32),
                    var=pd.DataFrame(index=pd.Index(["A", "B", "C", "D"])),
                ).write_h5ad(library_root / f"{library_id}.h5ad")
        return source, symbol_source, split, gene_order, library_manifest, library_root

    def _build(self, root: Path, **kwargs):
        source, symbol_source, split, gene_order, library_manifest, library_root = self._inputs(root, **kwargs)
        output = root / "fixture"
        receipt = builder.build(
            source,
            symbol_source,
            split,
            gene_order,
            library_manifest,
            library_root,
            output,
            expected_source_sha256=builder.sha256_file(source),
            expected_symbol_source_sha256=builder.sha256_file(symbol_source),
            expected_split_sha256=builder.sha256_file(split),
            expected_gene_order_sha256=builder.sha256_file(gene_order),
            expected_library_manifest_sha256=builder.sha256_file(library_manifest),
            expected_rows=5,
            expected_donors=5,
            expected_studies=5,
            expected_dimension=4,
            minimum_overlap=3,
        )
        return output, receipt

    def test_duplicate_symbols_are_summed_then_tp10k_log1p(self) -> None:
        import numpy as np
        from scipy import sparse

        with tempfile.TemporaryDirectory() as temporary:
            output, receipt = self._build(Path(temporary))
            self.assertEqual(receipt["status"], "pass_outcome_blind_fixture")
            self.assertEqual(receipt["gene_overlap"], 3)
            self.assertFalse(receipt["feature_symbol_source_observations_read"])
            self.assertEqual(receipt["aggregate_exposure_status"], "encoder_seen")
            self.assertTrue(receipt["study_gene_observability_is_global"])
            self.assertFalse(receipt["native_model_accepts_observability_mask"])
            self.assertEqual(receipt["evaluation_label_columns_read"], [])
            observed = sparse.load_npz(output / "normalized_counts.npz").toarray()
            expected_raw = np.asarray([10.0, 5.0, 5.0, 0.0], dtype=np.float32)
            expected = np.log1p(expected_raw / expected_raw.sum() * 10_000.0)
            np.testing.assert_allclose(observed[0], expected, rtol=1e-6, atol=1e-6)
            loaded = extractor.load_fixture(output, expected_rows=5, expected_dimension=4)
            self.assertEqual(loaded["rows"], [f"cell-{index}" for index in range(5)])

    def test_rejects_labels_in_the_model_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(builder.SCimilarityFixtureError, "firewall"):
                self._build(Path(temporary), forbidden=True)

    def test_rejects_a_study_split_join_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(builder.SCimilarityFixtureError, "join"):
                self._build(Path(temporary), bad_join=True)

    def test_exposure_roster_and_reference_firewall_are_frozen(self) -> None:
        self.assertEqual(builder.EXPOSURE_BY_STUDY["GSE185477"], "encoder_seen")
        self.assertEqual(builder.EXPOSURE_BY_STUDY["GSE136103"], "reference_only")
        self.assertEqual(
            sum(value == "clean_declared" for value in builder.EXPOSURE_BY_STUDY.values()),
            5,
        )
        source = Path(extractor.__file__).read_text(encoding="utf-8")
        self.assertNotIn("CellAnnotation", source)
        self.assertNotIn("CellSearchKNN", source)
        self.assertIn('"released_reference_index_used": False', source)
        self.assertIn("weights_only=True", source)
        self.assertIn('"exposure_status": "encoder_seen"', source)
        self.assertIn('representation = output / "common_embeddings.npz"', source)
        self.assertIn('set(representation_data.files) != {"embeddings", "outer_folds", "row_ids"}', source)

    def test_fixture_receipt_rejects_reference_index_use(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output, _ = self._build(Path(temporary))
            path = output / "fixture_receipt.json"
            receipt = json.loads(path.read_text(encoding="utf-8"))
            receipt["released_reference_index_used"] = True
            path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(extractor.SCimilarityExtractionError, "receipt"):
                extractor.load_fixture(output, expected_rows=5, expected_dimension=4)


if __name__ == "__main__":
    unittest.main()
