from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.build_pxd051911_activation_readiness import (
    PXD051911ActivationError,
    build_audit,
    observed,
    quant_sample_columns,
    tissue_join,
)


HEADER = "ProteinAccessions\tGenes\tProteinDescriptions"

META_FIELDS = [
    "unique_identifier",
    "patient_name",
    "sample_group",
    "gender",
    "alder",
    "bmi",
    "fat_pct",
    "kleiner_fibrosis_grade",
    "steatosis_score",
    "hepatocellular_ballooning_score",
    "lobular_inflammation_score",
    "nafld_activity_score",
    "saf_diagnosis",
    "liver_proteomics_filename",
    "oWAT_proteomics_filename",
    "scWAT_proteomics_filename",
]


def meta_row(patient: str, **overrides: str) -> str:
    row = {
        "unique_identifier": f"V1_{patient}",
        "patient_name": patient,
        "sample_group": "initial_sample",
        "gender": "Female",
        "alder": "50",
        "bmi": "40",
        "fat_pct": "45",
        "kleiner_fibrosis_grade": "F1",
        "steatosis_score": "1",
        "hepatocellular_ballooning_score": "0",
        "lobular_inflammation_score": "1",
        "nafld_activity_score": "2",
        "saf_diagnosis": "MASL",
        "liver_proteomics_filename": f"liver_{patient}.raw",
        "oWAT_proteomics_filename": f"owat_{patient}.raw",
        "scWAT_proteomics_filename": f"scwat_{patient}.raw",
    }
    row.update(overrides)
    return "\t".join(row[f] for f in META_FIELDS)


def write_source(directory: Path, rows: list[str], tissue_columns: dict) -> None:
    (directory / "meta_data.txt").write_text(
        "\t".join(META_FIELDS) + "\n" + "\n".join(rows) + "\n", encoding="utf-8"
    )
    for filename, columns in tissue_columns.items():
        body = HEADER + "".join(f"\t{c}" for c in columns) + "\n"
        body += "P00001\tALB\tAlbumin" + "\t1.0" * len(columns) + "\n"
        (directory / filename).write_text(body, encoding="utf-8")


class PXD051911JoinTests(unittest.TestCase):
    def test_missing_sentinels_are_not_treated_as_observed(self) -> None:
        self.assertFalse(observed("NA"))
        self.assertFalse(observed(""))
        self.assertFalse(observed(None))
        self.assertTrue(observed("0"))

    def test_annotation_columns_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "q.txt"
            path.write_text("Wrong\tHeader\tShape\ts1.raw\n", encoding="utf-8")
            with self.assertRaises(PXD051911ActivationError):
                quant_sample_columns(path)

    def test_join_reports_orphans_on_both_sides(self) -> None:
        rows = [
            {"patient_name": "ID1", "liver_proteomics_filename": "a.raw"},
            {"patient_name": "ID2", "liver_proteomics_filename": "NA"},
        ]
        join = tissue_join(rows, "liver_proteomics_filename", ["a.raw", "ghost.raw"])
        self.assertEqual(join["exact_filename_join"], 1)
        self.assertEqual(join["quant_columns_without_metadata_row"], ["ghost.raw"])
        self.assertEqual(join["participants"], ["ID1"])


class PXD051911AuditTests(unittest.TestCase):
    def build(self, rows: list[str], tissue_columns: dict) -> dict:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        directory = Path(holder.name)
        write_source(directory, rows, tissue_columns)
        return build_audit(directory, directory / "meta_data.txt")

    def test_matched_and_partial_adipose_coverage_stay_separable(self) -> None:
        rows = [
            meta_row("ID1"),
            meta_row("ID2", oWAT_proteomics_filename="NA"),
        ]
        audit = self.build(
            rows,
            {
                "liver_protein_quant.txt": ["liver_ID1.raw", "liver_ID2.raw"],
                "scwat_protein_quant.txt": ["scwat_ID1.raw", "scwat_ID2.raw"],
                "owat_protein_quant.txt": ["owat_ID1.raw"],
            },
        )
        self.assertEqual(audit["liver_participants"], 2)
        self.assertEqual(audit["cross_tissue"]["liver_and_scwat"], 2)
        self.assertEqual(audit["cross_tissue"]["liver_and_owat"], 1)
        self.assertEqual(audit["cross_tissue"]["liver_and_both_adipose"], 1)

    def test_repeat_visits_block_a_donor_safe_split(self) -> None:
        rows = [
            meta_row("ID1"),
            meta_row(
                "ID1",
                unique_identifier="V2_ID1",
                sample_group="follow_up_sample",
                liver_proteomics_filename="liver_ID1b.raw",
                oWAT_proteomics_filename="NA",
                scWAT_proteomics_filename="NA",
            ),
        ]
        with self.assertRaises(PXD051911ActivationError):
            self.build(
                rows,
                {
                    "liver_protein_quant.txt": ["liver_ID1.raw", "liver_ID1b.raw"],
                    "scwat_protein_quant.txt": ["scwat_ID1.raw"],
                    "owat_protein_quant.txt": ["owat_ID1.raw"],
                },
            )

    def test_incomplete_histology_blocks_activation(self) -> None:
        rows = [meta_row("ID1", kleiner_fibrosis_grade="NA")]
        with self.assertRaises(PXD051911ActivationError):
            self.build(
                rows,
                {
                    "liver_protein_quant.txt": ["liver_ID1.raw"],
                    "scwat_protein_quant.txt": ["scwat_ID1.raw"],
                    "owat_protein_quant.txt": ["owat_ID1.raw"],
                },
            )

    def test_partially_observed_field_is_masked_not_zeroed(self) -> None:
        rows = [meta_row("ID1", fat_pct="NA")]
        audit = self.build(
            rows,
            {
                "liver_protein_quant.txt": ["liver_ID1.raw"],
                "scwat_protein_quant.txt": ["scwat_ID1.raw"],
                "owat_protein_quant.txt": ["owat_ID1.raw"],
            },
        )
        self.assertEqual(
            audit["labels"]["fat_pct"]["state"], "partially_observed_typed_mask"
        )
        self.assertEqual(audit["labels"]["fat_pct"]["missing"], 1)
        self.assertFalse(audit["labels"]["fat_pct"]["encoded_as_zero"])

    def test_quant_column_without_a_participant_row_is_rejected(self) -> None:
        rows = [meta_row("ID1")]
        with self.assertRaises(PXD051911ActivationError):
            self.build(
                rows,
                {
                    "liver_protein_quant.txt": ["liver_ID1.raw", "unjoined.raw"],
                    "scwat_protein_quant.txt": ["scwat_ID1.raw"],
                    "owat_protein_quant.txt": ["owat_ID1.raw"],
                },
            )


if __name__ == "__main__":
    unittest.main()
