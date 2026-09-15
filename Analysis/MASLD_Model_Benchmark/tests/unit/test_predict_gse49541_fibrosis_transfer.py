from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.fit_gse267145_fibrosis_transfer_source_models import (
    apply_pipeline,
    fit_pipeline,
)
from scripts.predict_gse49541_fibrosis_transfer import (
    ExternalPredictionError,
    load_outcome_blind_masks,
    read_external_matrix,
)


def mask_payload(records):
    return {
        "schema_version": "masld-bench-gse49541-participant-mask-v1",
        "status": "pass_label_blind_participant_masks",
        "labels_read": False,
        "outcome_columns_read": False,
        "structurally_missing_imputed": False,
        "records": records,
    }


def clean_record(index: int) -> dict:
    record = {
        "cohort_family_id": "gse31803_gse49541_fibrosis_array",
        "series": "GSE49541",
        "participant_id": f"gse49541::{index:05d}",
        "sample_accession": f"GSM{789110 + index}",
        "timepoint": "baseline_cross_sectional",
        "repeat_topology": "one_record_per_participant",
        "prediction_eligible": True,
    }
    for field in ("exact_fibrosis_stage", "age", "sex"):
        record[field] = None
        record[f"{field}_observed"] = False
    return record


def write_mask(path: Path, records) -> None:
    path.write_text(json.dumps(mask_payload(records)) + "\n", encoding="utf-8")


def write_matrix(path: Path, gene_ids, accessions, values) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(["ensembl_gene_id"] + list(accessions)) + "\n")
        for row, gene_id in enumerate(gene_ids):
            handle.write(
                "\t".join(
                    [gene_id] + [format(float(values[column][row]), ".17g") for column in range(len(accessions))]
                )
                + "\n"
            )


class MaskFirewallTests(unittest.TestCase):
    def test_clean_mask_loads_and_sorts_by_participant(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "masks.json"
            write_mask(path, [clean_record(2), clean_record(1)])
            records = load_outcome_blind_masks(path)
            self.assertEqual(
                [record["participant_id"] for record in records],
                ["gse49541::00001", "gse49541::00002"],
            )

    def test_structurally_missing_nulls_are_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "masks.json"
            write_mask(path, [clean_record(1)])
            records = load_outcome_blind_masks(path)
            self.assertIsNone(records[0]["exact_fibrosis_stage"])

    def test_outcome_column_in_a_record_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "masks.json"
            leaked = clean_record(1)
            leaked["fibrosis_stage_group"] = "advanced_f3_f4"
            write_mask(path, [leaked])
            with self.assertRaises(ExternalPredictionError):
                load_outcome_blind_masks(path)

    def test_outcome_value_hidden_under_a_benign_key_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "masks.json"
            leaked = clean_record(1)
            leaked["timepoint"] = "advanced_f3_f4"
            write_mask(path, [leaked])
            with self.assertRaises(ExternalPredictionError):
                load_outcome_blind_masks(path)

    def test_filled_structurally_missing_field_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "masks.json"
            leaked = clean_record(1)
            leaked["exact_fibrosis_stage"] = "F3"
            leaked["exact_fibrosis_stage_observed"] = True
            write_mask(path, [leaked])
            with self.assertRaises(ExternalPredictionError):
                load_outcome_blind_masks(path)

    def test_mask_claiming_label_access_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "masks.json"
            payload = mask_payload([clean_record(1)])
            payload["labels_read"] = True
            path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
            with self.assertRaises(ExternalPredictionError):
                load_outcome_blind_masks(path)


class MatrixReadingTests(unittest.TestCase):
    def test_matrix_is_returned_in_shared_axis_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv"
            deposited = ["ENSG3", "ENSG1", "ENSG2"]
            accessions = ["GSM2", "GSM1"]
            values = {0: [30.0, 10.0, 20.0], 1: [31.0, 11.0, 21.0]}
            write_matrix(path, deposited, accessions, values)
            matrix = read_external_matrix(
                path,
                key_field="ensembl_gene_id",
                gene_ids=["ENSG1", "ENSG2"],
                accessions=["GSM1", "GSM2"],
            )
            # Row order follows the requested accessions, column order the axis.
            self.assertTrue(np.allclose(matrix, [[11.0, 21.0], [10.0, 20.0]]))

    def test_missing_shared_gene_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv"
            write_matrix(path, ["ENSG1"], ["GSM1"], {0: [1.0]})
            with self.assertRaises(ExternalPredictionError):
                read_external_matrix(
                    path,
                    key_field="ensembl_gene_id",
                    gene_ids=["ENSG1", "ENSG2"],
                    accessions=["GSM1"],
                )

    def test_undeclared_array_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv"
            write_matrix(path, ["ENSG1"], ["GSM1"], {0: [1.0]})
            with self.assertRaises(ExternalPredictionError):
                read_external_matrix(
                    path,
                    key_field="ensembl_gene_id",
                    gene_ids=["ENSG1"],
                    accessions=["GSM9"],
                )

    def test_outcome_column_in_the_matrix_header_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv"
            write_matrix(path, ["ENSG1"], ["GSM1", "fibrosis_stage_group"], {0: [1.0], 1: [1.0]})
            with self.assertRaises(ExternalPredictionError):
                read_external_matrix(
                    path,
                    key_field="ensembl_gene_id",
                    gene_ids=["ENSG1"],
                    accessions=["GSM1"],
                )


class TransferApplicationTests(unittest.TestCase):
    def test_a_source_state_applies_to_a_disjoint_participant_block(self) -> None:
        """The predictor's only job: apply, never refit.

        Applying a frozen state to a held block must not depend on which other
        participants are in that block.
        """

        rng = np.random.default_rng(31)
        labels = np.asarray([0] * 20 + [1] * 6, dtype=np.int64)
        source = rng.normal(size=(26, 30))
        source[labels == 1] += 1.5
        config = {
            "representation": "gene_median",
            "reducer": "none",
            "classifier": "elastic_net",
        }
        state = fit_pipeline(
            model_id="gene_median_elastic_net",
            config=config,
            representation=source,
            eligibility=np.ones_like(source),
            gene_ids=[f"ENSG{index:05d}" for index in range(30)],
            fitting=np.arange(26, dtype=np.int64),
            labels=labels,
            feature_count=10,
            pca_components=3,
            hyperparameters={"c": 1.0, "l1_ratio": 0.5},
            seed=1701,
        )
        state.pop("_converged")
        state.pop("_design_columns")
        external = rng.normal(size=(8, 30))
        whole = apply_pipeline(config=config, state=state, representation=external)
        for index in range(8):
            alone = apply_pipeline(
                config=config, state=state, representation=external[index : index + 1]
            )
            self.assertAlmostEqual(float(alone[0]), float(whole[index]), places=12)


if __name__ == "__main__":
    unittest.main()
