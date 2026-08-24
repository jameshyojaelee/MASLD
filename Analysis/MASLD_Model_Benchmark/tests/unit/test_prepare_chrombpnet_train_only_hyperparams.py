from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts import prepare_chrombpnet_train_only_hyperparams as preparation


class TrainOnlyHyperparameterPreparationTests(unittest.TestCase):
    def test_validation_is_moved_to_outcome_blind_hyperparameter_test_role(self) -> None:
        actual = {
            "train": ["chr3", "chr4"],
            "valid": ["chr2"],
            "test": ["chr1"],
        }
        self.assertEqual(
            preparation.build_hyperparameter_fold(actual),
            {
                "train": ["chr3", "chr4"],
                "valid": [],
                "test": ["chr2", "chr1"],
            },
        )

    def test_only_fold_path_is_rewritten_in_final_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.tsv"
            rows = [
                ("counts_loss_weight", "10"),
                ("filters", "512"),
                ("n_dil_layers", "8"),
                ("bias_model_path", "/tmp/bias.h5"),
                ("inputlen", "2114"),
                ("outputlen", "1000"),
                ("max_jitter", "500"),
                ("chr_fold_path", "/tmp/hyperparameter.json"),
                ("negative_sampling_ratio", "0.1"),
            ]
            source.write_text(
                "".join(f"{key}\t{value}\n" for key, value in rows),
                encoding="utf-8",
            )
            output = root / "final.tsv"
            rewritten = preparation.write_actual_fold_parameters(
                source,
                "actual_fold.json",
                "train_only_bias_model_scaled.h5",
                output,
            )
            self.assertEqual(rewritten["chr_fold_path"], "actual_fold.json")
            self.assertEqual(rewritten["counts_loss_weight"], "10")
            self.assertEqual(
                rewritten["bias_model_path"], "train_only_bias_model_scaled.h5"
            )

    def test_held_rows_treat_numeric_rendering_as_semantically_equal(self) -> None:
        actual = {
            "train": ["chr3"],
            "valid": ["chr2"],
            "test": ["chr1"],
        }
        source = [
            ("chr1", "100", "200", "peak_a", "0", ".", "0", "-1", "-1", "50"),
            ("chr2", "300", "400", "peak_b", "0", "+", "4", "5", "6", "25"),
        ]
        rewritten = [
            ("chr1", "100", "200", "peak_a", "0", ".", "0.0", "-1.0", "-1.0", "50"),
            ("chr2", "300", "400", "peak_b", "0", "+", "4.0", "5.00", "6e0", "25"),
        ]
        self.assertEqual(
            preparation.held_rows(source, actual),
            preparation.held_rows(rewritten, actual),
        )

    def test_held_rows_preserve_explicit_missing_numeric_sentinels(self) -> None:
        actual = {"train": ["chr3"], "valid": ["chr2"], "test": ["chr1"]}
        missing = [
            ("chr1", "100", "200", ".", ".", ".", ".", ".", ".", "50")
        ]
        numeric_zero = [
            ("chr1", "100", "200", ".", ".", ".", "0", "0", "0", "50")
        ]
        self.assertEqual(
            preparation.held_rows(missing, actual),
            preparation.held_rows(list(missing), actual),
        )
        self.assertNotEqual(
            preparation.held_rows(missing, actual),
            preparation.held_rows(numeric_zero, actual),
        )

    def test_held_rows_preserve_id_strand_and_summit_sensitivity(self) -> None:
        actual = {"train": ["chr3"], "valid": ["chr2"], "test": ["chr1"]}
        source = [
            ("chr1", "100", "200", "peak_a", "0", ".", "0", "-1", "-1", "50")
        ]
        variants = [
            ("chr1", "100", "200", "peak_b", "0", ".", "0.0", "-1.0", "-1.0", "50"),
            ("chr1", "100", "200", "peak_a", "0", "+", "0.0", "-1.0", "-1.0", "50"),
            ("chr1", "100", "200", "peak_a", "0", ".", "0.0", "-1.0", "-1.0", "51"),
        ]
        for variant in variants:
            with self.subTest(variant=variant):
                self.assertNotEqual(
                    preparation.held_rows(source, actual),
                    preparation.held_rows([variant], actual),
                )

    def test_narrowpeak_semantics_reject_nonfinite_payloads(self) -> None:
        row = ("chr1", "100", "200", "peak_a", "0", ".", "nan", "-1", "-1", "50")
        with self.assertRaisesRegex(
            preparation.TrainOnlyHyperparameterError, "non-finite"
        ):
            preparation.narrowpeak_semantic_row(row)


if __name__ == "__main__":
    unittest.main()
