from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.predict_gse83452_baseline_nash_transfer import (
    FORBIDDEN_MASK_KEYS,
    OUTCOME_TOKENS,
    STRUCTURALLY_WITHHELD_KEYS,
    ExternalPredictionError,
    load_outcome_blind_masks,
    read_external_matrix,
)


def make_mask_payload(records: int = 3) -> dict:
    return {
        "status": "pass_label_blind_participant_masks",
        "labels_read": False,
        "outcome_columns_read": False,
        "withheld_fields_imputed": False,
        "selection_is_label_blind": True,
        "arrays_per_participant_after_selection": 1,
        "mask_records": [
            {
                "participant_id": f"gse83452::GSM{2203254 + index}",
                "sample_accession": f"GSM{2203254 + index}",
                "timepoint": "baseline",
                "repeat_topology": "baseline_only_or_pair_anchor",
                "outer_group_key": f"antwerp_inserm_shared::gse83452::GSM{2203254 + index}",
                "prediction_eligible": True,
                "baseline_selected": True,
                "nash_status": None,
                "nash_status_observed": False,
                "age": None,
                "age_observed": False,
                "sex": None,
                "sex_observed": False,
                "intervention": None,
                "intervention_observed": False,
            }
            for index in range(records)
        ],
    }


def write_matrix(path: Path, gene_ids: list[str], accessions: list[str]) -> None:
    lines = ["\t".join(["ensembl_gene_id"] + accessions)]
    for offset, gene_id in enumerate(gene_ids):
        lines.append(
            "\t".join(
                [gene_id]
                + [f"{6.0 + offset + index * 0.5:.6f}" for index in range(len(accessions))]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class MaskGuardTests(unittest.TestCase):
    def _load(self, payload: dict) -> list[dict]:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "masks.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            return load_outcome_blind_masks(path)

    def test_a_clean_baseline_mask_loads_sorted_by_participant(self) -> None:
        records = self._load(make_mask_payload())
        self.assertEqual(len(records), 3)
        self.assertEqual(
            [record["participant_id"] for record in records],
            sorted(record["participant_id"] for record in records),
        )

    def test_a_filled_outcome_is_rejected(self) -> None:
        payload = make_mask_payload()
        payload["mask_records"][0]["nash_status"] = "nash"
        with self.assertRaises(ExternalPredictionError):
            self._load(payload)

    def test_a_filled_covariate_is_rejected(self) -> None:
        payload = make_mask_payload()
        payload["mask_records"][1]["age"] = "55"
        with self.assertRaises(ExternalPredictionError):
            self._load(payload)

    def test_a_forbidden_column_is_rejected(self) -> None:
        payload = make_mask_payload()
        payload["mask_records"][0]["primary_transfer_evaluable"] = "True"
        with self.assertRaises(ExternalPredictionError):
            self._load(payload)

    def test_an_outcome_token_smuggled_as_a_value_is_rejected(self) -> None:
        payload = make_mask_payload()
        payload["mask_records"][2]["repeat_topology"] = "nash"
        with self.assertRaises(ExternalPredictionError):
            self._load(payload)

    def test_a_non_baseline_record_is_rejected(self) -> None:
        payload = make_mask_payload()
        payload["mask_records"][0]["baseline_selected"] = False
        with self.assertRaises(ExternalPredictionError):
            self._load(payload)

    def test_a_mask_that_pools_repeat_arrays_is_rejected(self) -> None:
        payload = make_mask_payload()
        payload["arrays_per_participant_after_selection"] = 2
        with self.assertRaises(ExternalPredictionError):
            self._load(payload)

    def test_the_guard_rosters_are_disjoint(self) -> None:
        self.assertEqual(
            set(FORBIDDEN_MASK_KEYS) & set(STRUCTURALLY_WITHHELD_KEYS), set()
        )
        self.assertIn("nash", OUTCOME_TOKENS)
        self.assertIn("no_nash", OUTCOME_TOKENS)


class MatrixReaderTests(unittest.TestCase):
    def test_only_the_declared_baseline_columns_are_decoded(self) -> None:
        """Follow-up columns sit in the same file and must never be read."""

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            genes = ["ENSG00000000003", "ENSG00000000005"]
            accessions = ["GSM1", "GSM2", "GSM3"]
            write_matrix(path, genes, accessions)
            matrix = read_external_matrix(
                path,
                key_field="ensembl_gene_id",
                gene_ids=genes,
                accessions=["GSM1", "GSM3"],
            )
            self.assertEqual(matrix.shape, (2, 2))
            # Column GSM2 carries 6.5 and 7.5; neither may appear.
            self.assertNotIn(6.5, set(matrix.ravel().tolist()))

    def test_the_matrix_is_returned_in_shared_axis_order(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            write_matrix(path, ["ENSG00000000005", "ENSG00000000003"], ["GSM1"])
            matrix = read_external_matrix(
                path,
                key_field="ensembl_gene_id",
                gene_ids=["ENSG00000000003", "ENSG00000000005"],
                accessions=["GSM1"],
            )
            # ENSG...003 was deposited second and carries 7.0.
            self.assertAlmostEqual(float(matrix[0, 0]), 7.0)
            self.assertAlmostEqual(float(matrix[0, 1]), 6.0)

    def test_an_incomplete_shared_axis_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            write_matrix(path, ["ENSG00000000003"], ["GSM1"])
            with self.assertRaises(ExternalPredictionError):
                read_external_matrix(
                    path,
                    key_field="ensembl_gene_id",
                    gene_ids=["ENSG00000000003", "ENSG00000000005"],
                    accessions=["GSM1"],
                )

    def test_a_missing_declared_array_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            write_matrix(path, ["ENSG00000000003"], ["GSM1"])
            with self.assertRaises(ExternalPredictionError):
                read_external_matrix(
                    path,
                    key_field="ensembl_gene_id",
                    gene_ids=["ENSG00000000003"],
                    accessions=["GSM9"],
                )

    def test_an_outcome_column_in_the_matrix_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            write_matrix(path, ["ENSG00000000003"], ["GSM1", "nash_status"])
            with self.assertRaises(ExternalPredictionError):
                read_external_matrix(
                    path,
                    key_field="ensembl_gene_id",
                    gene_ids=["ENSG00000000003"],
                    accessions=["GSM1"],
                )


class RepresentationParityTests(unittest.TestCase):
    def test_the_predictor_imports_the_source_representation_code(self) -> None:
        """Source and target must not drift onto two definitions."""

        import scripts.fit_gse135251_nash_transfer_source_models as fit
        import scripts.predict_gse83452_baseline_nash_transfer as predict

        self.assertIs(predict.build_representation, fit.build_representation)
        self.assertIs(predict.apply_pipeline, fit.apply_pipeline)
        self.assertEqual(predict.CLASSES, fit.CLASSES)

    def test_the_per_array_transform_is_identical_on_both_sides(self) -> None:
        import scripts.fit_gse135251_nash_transfer_source_models as fit

        values = np.asarray([[1.0, 4.0, 2.0, 8.0], [3.0, 3.0, 9.0, 1.0]])
        for name in ("per_array_rank", "gene_median"):
            built = fit.build_representation(name=name, log2_shared=values)
            for index in range(values.shape[0]):
                alone = fit.REPRESENTATION_FUNCTIONS[name](values[index : index + 1])
                self.assertTrue(np.array_equal(alone[0], built[index]))


if __name__ == "__main__":
    unittest.main()
