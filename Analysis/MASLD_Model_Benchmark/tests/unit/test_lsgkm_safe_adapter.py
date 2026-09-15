from __future__ import annotations

import gzip
import math
from pathlib import Path
import tempfile
import unittest

from scripts.lsgkm_safe_adapter import (
    LSGKMSafetyError,
    canonical_kmer,
    canonical_score_table,
    deltasvm_scores,
    direct_gkmsvm_alt_minus_ref,
    parse_predictions,
    reverse_complement,
    validate_fasta,
    validate_project_model,
)


def model_text(*, rho: str = "0.125", total_sv: int = 2) -> str:
    support = [f"0.5 {'A' * 300}", f"-0.5 {'C' * 300}"]
    return "\n".join(
        [
            "svm_type c_svc",
            "kernel_type gkm_esttrunc",
            "L 11",
            "k 7",
            "d 3",
            "norc 0",
            "nr_class 2",
            f"total_sv {total_sv}",
            f"rho {rho}",
            "label 1 -1",
            "nr_sv 1 1",
            "SV",
            *support,
            "",
        ]
    )


class LSGKMSafeAdapterTests(unittest.TestCase):
    def test_valid_fasta_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "valid.fa"
            path.write_text(f">first\n{'ACGT' * 75}\n>second\n{'TGCA' * 75}\n")
            records = validate_fasta(path, expected_length=300)
            self.assertEqual([name for name, _ in records], ["first", "second"])

    def test_fasta_rejects_ambiguity_lowercase_duplicates_and_mixed_lengths(self) -> None:
        invalid = {
            "ambiguity": ">x\nACNT\n",
            "lowercase": ">x\nacgt\n",
            "duplicate": ">x\nACGT\n>x\nACGT\n",
            "mixed": ">x\nACGT\n>y\nACGTA\n",
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, content in invalid.items():
                with self.subTest(name=name):
                    path = Path(directory) / f"{name}.fa"
                    path.write_text(content)
                    with self.assertRaises(LSGKMSafetyError):
                        validate_fasta(path)

    def test_fasta_rejects_overlength_before_native_code(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "overlength.fa"
            path.write_text(f">x\n{'A' * 2048}\n")
            with self.assertRaises(LSGKMSafetyError):
                validate_fasta(path)

    def test_valid_gzip_model_is_bounded_and_parsed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.txt.gz"
            with gzip.GzipFile(filename=path, mode="wb", mtime=0) as handle:
                handle.write(model_text().encode("ascii"))
            receipt = validate_project_model(path)
            self.assertTrue(receipt["gzip"])
            self.assertEqual(receipt["total_sv"], 2)
            self.assertEqual(receipt["label"], [1, -1])

    def test_model_rejects_header_drift_nonfinite_and_count_drift(self) -> None:
        variants = {
            "header": model_text().replace("kernel_type gkm_esttrunc", "kernel_type gkm_cnt"),
            "rho": model_text(rho="nan"),
            "count": model_text(total_sv=3),
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, content in variants.items():
                with self.subTest(name=name):
                    path = Path(directory) / f"{name}.model.txt"
                    path.write_text(content)
                    with self.assertRaises(LSGKMSafetyError):
                        validate_project_model(path)

    def test_prediction_parser_requires_exact_order_and_finite_scores(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scores.tsv"
            path.write_text("ref\t-0.25\nalt\t0.75\n")
            self.assertEqual(
                parse_predictions(path, expected_identifiers=["ref", "alt"]),
                {"ref": -0.25, "alt": 0.75},
            )
            path.unlink()
            path.write_text("ref\tnan\n")
            with self.assertRaises(LSGKMSafetyError):
                parse_predictions(path)

    def test_direct_gkmsvm_score_is_alt_minus_ref(self) -> None:
        self.assertEqual(
            direct_gkmsvm_alt_minus_ref(reference_score=-0.25, alternative_score=0.75),
            1.0,
        )

    def test_deltasvm_matches_hand_calculated_sign(self) -> None:
        reference = "AACCGGTTAACCGGTTAACCG"
        alternative = reference[:10] + "G" + reference[11:]
        keys = {
            canonical_kmer(context[start : start + 11])
            for context in (reference, alternative)
            for start in range(11)
        }
        weights = {key: float(index + 1) for index, key in enumerate(sorted(keys))}
        expected_ref = math.fsum(
            weights[canonical_kmer(reference[start : start + 11])]
            for start in range(11)
        )
        expected_alt = math.fsum(
            weights[canonical_kmer(alternative[start : start + 11])]
            for start in range(11)
        )
        observed = deltasvm_scores(reference, alternative, weights)
        self.assertEqual(observed["native_ref_minus_alt"], expected_ref - expected_alt)
        self.assertEqual(
            observed["canonical_alt_minus_ref"], expected_alt - expected_ref
        )

    def test_deltasvm_is_reverse_complement_invariant(self) -> None:
        reference = "AACCGGTTAACCGGTTAACCG"
        alternative = reference[:10] + "G" + reference[11:]
        keys = {
            canonical_kmer(context[start : start + 11])
            for context in (reference, alternative)
            for start in range(11)
        }
        weights = {key: float(index) / 7.0 for index, key in enumerate(sorted(keys))}
        forward = deltasvm_scores(reference, alternative, weights)
        reverse = deltasvm_scores(
            reverse_complement(reference), reverse_complement(alternative), weights
        )
        self.assertEqual(forward, reverse)

    def test_deltasvm_rejects_indels_off_center_and_missing_weights(self) -> None:
        reference = "AACCGGTTAACCGGTTAACCG"
        off_center = "G" + reference[1:]
        with self.assertRaises(LSGKMSafetyError):
            deltasvm_scores(reference, off_center, {})
        centered = reference[:10] + "G" + reference[11:]
        with self.assertRaises(LSGKMSafetyError):
            deltasvm_scores(reference, centered, {})

    def test_canonical_score_table_does_not_average_rc_drift(self) -> None:
        self.assertEqual(
            canonical_score_table([("AAC", 1.25), ("GTT", 1.25)]),
            {canonical_kmer("AAC"): 1.25},
        )
        with self.assertRaises(LSGKMSafetyError):
            canonical_score_table([("AAC", 1.25), ("GTT", 1.5)])


if __name__ == "__main__":
    unittest.main()
