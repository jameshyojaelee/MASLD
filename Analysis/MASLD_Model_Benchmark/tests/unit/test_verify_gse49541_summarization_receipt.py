from __future__ import annotations

import copy
import csv
from pathlib import Path
import tempfile
import unittest

from scripts.verify_gse49541_summarization_receipt import (
    EXPECTED_PACKAGE_VERSIONS,
    REQUIRED_FALSE_FLAGS,
    SummarizationError,
    read_matrix_axes,
    verify_summarization,
)

ARRAYS = [f"GSM{index:06d}" for index in range(72)]
FEATURES = [f"{index}_at" for index in range(54_675)]
GENES = [f"ENSG{index:011d}" for index in range(1_000)]


def write_matrix(path: Path, id_column: str, rows: list[str], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([id_column, *columns])
        for row in rows:
            writer.writerow([row, *["1.0"] * len(columns)])


def per_array_record(accession: str) -> dict[str, object]:
    return {
        "sample_accession": accession,
        "arrays_read_in_this_call": 1,
        "cdf_name": "HG-U133_Plus_2",
        "array_type_matches_GPL570": True,
        "summary_finite": 54_675,
        "excluded": False,
        "descriptive_flag_only": True,
    }


def receipt() -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": "masld-bench-gse49541-gpl570-single-array-summary-v1",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "method": "frma_with_exact_hgu133plus2frmavecs",
        "arrays_read_per_frma_call": 1,
        "arrays_summarized": 72,
        "platform_features": 54_675,
        "gene_features": 1_000,
        "eligible_platform_features": 39_385,
        "label_blind": True,
        "samples_excluded": 0,
        "package_versions": dict(EXPECTED_PACKAGE_VERSIONS),
        "preprocessing_object_sha256": "a" * 64,
        "platform_feature_matrix_sha256": "b" * 64,
        "gene_matrix_sha256": "c" * 64,
        "hard_failures": {
            "CEL_parse_failure": 0,
            "array_type_mismatch": 0,
            "sample_axis_mismatch": 0,
            "nonfinite_summary": 0,
            "duplicate_feature_or_sample_ID": 0,
            "no_positive_intensity": 0,
        },
        "per_array": [per_array_record(accession) for accession in ARRAYS],
    }
    for flag in REQUIRED_FALSE_FLAGS:
        payload[flag] = False
    return payload


class VerifySummarizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        base = Path(self.directory.name)
        self.platform = base / "platform.tsv"
        self.gene = base / "gene.tsv"
        write_matrix(self.platform, "platform_feature_id", FEATURES, ARRAYS)
        write_matrix(self.gene, "ensembl_gene_id", GENES, ARRAYS)

    def verify(self, payload: dict[str, object]) -> dict[str, object]:
        return verify_summarization(
            payload, platform_matrix=self.platform, gene_matrix=self.gene
        )

    def test_clean_summarization_passes(self) -> None:
        verdict = self.verify(receipt())
        self.assertEqual(verdict["status"], "pass_label_blind_single_array_summarization")
        self.assertEqual(verdict["arrays_summarized"], 72)
        self.assertEqual(verdict["platform_features"], 54_675)
        self.assertEqual(verdict["gene_features"], 1_000)
        self.assertEqual(verdict["arrays_read_per_frma_call"], 1)

    def test_batched_frma_call_is_rejected(self) -> None:
        payload = receipt()
        payload["arrays_read_per_frma_call"] = 72
        with self.assertRaisesRegex(SummarizationError, "more than one array"):
            self.verify(payload)

    def test_per_array_batching_is_rejected(self) -> None:
        payload = receipt()
        payload["per_array"][5]["arrays_read_in_this_call"] = 6  # type: ignore[index]
        with self.assertRaisesRegex(SummarizationError, "alongside another array"):
            self.verify(payload)

    def test_hard_qc_failure_is_rejected(self) -> None:
        payload = receipt()
        payload["hard_failures"]["array_type_mismatch"] = 1  # type: ignore[index]
        with self.assertRaisesRegex(SummarizationError, "hard QC failure"):
            self.verify(payload)

    def test_sample_exclusion_is_rejected(self) -> None:
        payload = receipt()
        payload["samples_excluded"] = 1
        with self.assertRaisesRegex(SummarizationError, "excluded from a label-blind stage"):
            self.verify(payload)

    def test_per_array_exclusion_is_rejected(self) -> None:
        payload = receipt()
        payload["per_array"][0]["excluded"] = True  # type: ignore[index]
        with self.assertRaisesRegex(SummarizationError, "was excluded"):
            self.verify(payload)

    def test_every_forbidden_flag_must_be_false(self) -> None:
        for flag in REQUIRED_FALSE_FLAGS:
            payload = copy.deepcopy(receipt())
            payload[flag] = True
            with self.assertRaisesRegex(SummarizationError, flag):
                self.verify(payload)

    def test_receipt_cannot_claim_a_shape_the_matrix_lacks(self) -> None:
        payload = receipt()
        payload["gene_features"] = 18_542
        with self.assertRaisesRegex(SummarizationError, "gene axis differs"):
            self.verify(payload)

    def test_disagreeing_sample_axes_are_rejected(self) -> None:
        write_matrix(self.gene, "ensembl_gene_id", GENES, list(reversed(ARRAYS)))
        with self.assertRaisesRegex(SummarizationError, "disagree on the sample axis"):
            self.verify(receipt())

    def test_per_array_qc_must_cover_the_matrix(self) -> None:
        payload = receipt()
        payload["per_array"] = payload["per_array"][:71]  # type: ignore[index]
        with self.assertRaisesRegex(SummarizationError, "does not cover 72 arrays"):
            self.verify(payload)

    def test_foreign_array_type_is_rejected(self) -> None:
        payload = receipt()
        payload["per_array"][3]["cdf_name"] = "HuGene-2_0-st"  # type: ignore[index]
        with self.assertRaisesRegex(SummarizationError, "is not a HG-U133_Plus_2 array"):
            self.verify(payload)

    def test_nonfinite_summary_is_rejected(self) -> None:
        payload = receipt()
        payload["per_array"][9]["summary_finite"] = 54_674  # type: ignore[index]
        with self.assertRaisesRegex(SummarizationError, "nonfinite value"):
            self.verify(payload)

    def test_blocked_method_is_rejected(self) -> None:
        payload = receipt()
        payload["method"] = "all_sample_RMA"
        with self.assertRaisesRegex(SummarizationError, "contract does not admit"):
            self.verify(payload)

    def test_r_length_one_lists_are_unwrapped(self) -> None:
        payload = receipt()
        payload["arrays_read_per_frma_call"] = [1]
        payload["label_blind"] = [True]
        payload["samples_excluded"] = [0]
        payload["platform_features"] = [54_675]
        payload["package_versions"] = {
            name: [value] for name, value in EXPECTED_PACKAGE_VERSIONS.items()
        }
        payload["per_array"][0]["cdf_name"] = ["HG-U133_Plus_2"]  # type: ignore[index]
        self.assertEqual(self.verify(payload)["arrays_summarized"], 72)

    def test_ragged_matrix_is_rejected(self) -> None:
        with self.platform.open("a", encoding="utf-8") as handle:
            handle.write("extra_at\t1.0\n")
        with self.assertRaisesRegex(SummarizationError, "ragged"):
            self.verify(receipt())

    def test_matrix_axes_reject_repeated_identifiers(self) -> None:
        write_matrix(self.gene, "ensembl_gene_id", ["ENSG1", "ENSG1"], ARRAYS)
        with self.assertRaisesRegex(SummarizationError, "repeats a feature identifier"):
            read_matrix_axes(self.gene, "ensembl_gene_id")


if __name__ == "__main__":
    unittest.main()
