"""Input-orientation regression checks using the unchanged v1.2 model assets."""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts import score_severity_v12_corrected as corrected


class CorrectedSeverityReaderTests(unittest.TestCase):
    def test_transposed_npz_generates_ids_for_samples(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "counts.npz")
            counts = np.asarray([[1, 2, 3], [4, 5, 6]])
            np.savez(path, counts=counts, genes=["g1", "g2", "g3"])
            mat, genes, samples = corrected.read_counts(
                path, orientation="samples_x_genes",
            )
            np.testing.assert_array_equal(mat, counts.T)
            np.testing.assert_array_equal(genes, ["g1", "g2", "g3"])
            np.testing.assert_array_equal(samples, ["sample_1", "sample_2"])

    def test_previously_accepted_npz_inputs_are_identical(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for orientation in ("samples_x_genes", "genes_x_samples"):
                for explicit_samples in (False, True):
                    counts = np.asarray([[1, 4], [2, 5], [3, 6]])
                    if orientation == "samples_x_genes":
                        counts = counts.T
                    if orientation == "samples_x_genes" and not explicit_samples:
                        continue  # The corrected, previously refused case above.
                    path = str(Path(directory) / "counts.npz")
                    kwargs = {"samples": ["a", "b"]} if explicit_samples else {}
                    np.savez(path, counts=counts, genes=["g1", "g2", "g3"], **kwargs)
                    expected = corrected._original_read_counts(path, orientation=orientation)
                    actual = corrected.read_counts(path, orientation=orientation)
                    for before, after in zip(expected, actual):
                        np.testing.assert_array_equal(before, after)

    def test_corrected_orientation_retains_count_validation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "counts.npz")
            for value in (float("nan"), float("inf"), -1.0):
                np.savez(path, counts=[[1.0, value]], genes=["g1", "g2"])
                with self.assertRaises(corrected.source.ScoreError):
                    corrected.read_counts(path, orientation="samples_x_genes")

    def test_transposed_inputs_give_identical_real_model_scores(self) -> None:
        with np.load(corrected.source.DEFAULT_WEIGHTS, allow_pickle=True) as weights:
            genes = weights["gene_axis"].astype(str)
        counts = np.arange(1, len(genes) + 1, dtype=float) * 10
        counts = np.column_stack((counts, counts[::-1] * 2))
        with tempfile.TemporaryDirectory() as directory:
            original_path = str(Path(directory) / "original.npz")
            corrected_path = str(Path(directory) / "transposed.npz")
            np.savez(original_path, counts=counts, genes=genes)
            np.savez(corrected_path, counts=counts.T, genes=genes)
            original = corrected._original_read_counts(original_path)
            repaired = corrected.read_counts(corrected_path, orientation="samples_x_genes")
            for before, after in zip(original, repaired):
                np.testing.assert_array_equal(before, after)
            expected = corrected.source.score(*original[:2])
            actual = corrected.source.score(*repaired[:2])

        def check(before, after):
            if isinstance(before, dict):
                self.assertEqual(before.keys(), after.keys())
                for key in before:
                    check(before[key], after[key])
            elif isinstance(before, np.ndarray):
                np.testing.assert_array_equal(before, after)
            else:
                self.assertEqual(before, after)

        check(expected, actual)


if __name__ == "__main__":
    unittest.main()
