from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts import bpnet_contract as contract


class BPNetContractTests(unittest.TestCase):
    def test_architecture_is_exact_bias_free_control(self) -> None:
        parameters = contract.architecture_parameters(42.5)
        self.assertEqual(parameters["input_len"], 2114)
        self.assertEqual(parameters["output_profile_len"], 1000)
        self.assertEqual(parameters["motif_module_params"]["filters"], [64])
        self.assertEqual(
            parameters["syntax_module_params"]["num_dilation_layers"], 8
        )
        self.assertEqual(parameters["syntax_module_params"]["filters"], 64)
        self.assertEqual(parameters["loss_weights"], [1.0, 42.5])

    def test_invalid_counts_weight_fails_closed(self) -> None:
        for value in (0.0, -1.0, float("nan"), float("inf")):
            with self.assertRaises(contract.BPNetContractError):
                contract.architecture_parameters(value)

    def test_fold_roles_must_be_disjoint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fold.json"
            path.write_text(
                json.dumps(
                    {
                        "train": ["chr1", "chr2"],
                        "valid": ["chr2"],
                        "test": ["chr3"],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaises(contract.BPNetContractError):
                contract.read_fold(path, require_primary_partition=False)

    def test_narrowpeak_normalization_preserves_summit(self) -> None:
        row = contract.normalize_narrowpeak_fields(
            ["chr1", "100", "250", "peak", "12", ".", "8", "4", "3", "70"],
            row_number=1,
        )
        self.assertEqual(row, ("chr1", "100", "250", "peak", "12", ".", "8", "4", "3", "70"))
        short = contract.normalize_narrowpeak_fields(
            ["chr2", "100", "300", "negative"], row_number=2
        )
        self.assertEqual(short[9], "100")

    def test_prepare_narrowpeak_filters_disallowed_contigs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "regions.bed"
            output = root / "bpnet.narrowPeak"
            source.write_text(
                "chr1\t100\t300\tone\nchr2\t200\t400\ttwo\n", encoding="utf-8"
            )
            count = contract.prepare_narrowpeak(
                source, output, allowed_contigs={"chr2"}
            )
            self.assertEqual(count, 1)
            self.assertTrue(output.read_text(encoding="utf-8").startswith("chr2\t"))

    def test_prepare_narrowpeak_requires_full_jittered_window(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "regions.bed"
            output = root / "bpnet.narrowPeak"
            source.write_text(
                "chr1\t90\t110\tedge\nchr1\t4990\t5010\tsafe\n",
                encoding="utf-8",
            )
            count = contract.prepare_narrowpeak(
                source,
                output,
                allowed_contigs={"chr1"},
                chrom_sizes={"chr1": 10_000},
                required_flank=1_185,
            )
            self.assertEqual(count, 1)
            self.assertIn("\tsafe\t", output.read_text(encoding="utf-8"))

    def test_prepare_narrowpeak_uses_role_specific_flanks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "regions.bed"
            output = root / "bpnet.narrowPeak"
            source.write_text(
                "chrTrain\t1090\t1110\ttrain_edge\n"
                "chrValid\t1090\t1110\tvalid_safe\n",
                encoding="utf-8",
            )
            count = contract.prepare_narrowpeak(
                source,
                output,
                allowed_contigs={"chrTrain", "chrValid"},
                chrom_sizes={"chrTrain": 10_000, "chrValid": 10_000},
                required_flank_by_contig={"chrTrain": 1_185, "chrValid": 1_057},
            )
            self.assertEqual(count, 1)
            self.assertIn("valid_safe", output.read_text(encoding="utf-8"))

    def test_read_chrom_sizes_rejects_duplicate_contig(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "chrom.sizes"
            path.write_text("chr1\t1000\nchr1\t2000\n", encoding="utf-8")
            with self.assertRaises(contract.BPNetContractError):
                contract.read_chrom_sizes(path)

    def test_outcome_named_input_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evaluator_outcomes.tsv"
            path.write_text("x\n", encoding="utf-8")
            with self.assertRaises(contract.BPNetContractError):
                contract.validate_model_inputs((path,))

    def test_symlink_input_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target.tsv"
            link = root / "link.tsv"
            target.write_text("x\n", encoding="utf-8")
            link.symlink_to(target)
            with self.assertRaises(contract.BPNetContractError):
                contract.validate_model_inputs((link,))

    def test_one_hot_reverse_complement_and_softmax(self) -> None:
        encoded = contract.one_hot_dna(["ACGT", "TGCA"])
        restored = contract.reverse_complement_one_hot(
            contract.reverse_complement_one_hot(encoded)
        )
        np.testing.assert_array_equal(encoded, restored)
        probabilities = contract.softmax(np.asarray([[0.0, 1.0], [2.0, 2.0]]))
        np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)
        with self.assertRaises(contract.BPNetContractError):
            contract.one_hot_dna(["ACGN"])


if __name__ == "__main__":
    unittest.main()
