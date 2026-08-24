from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts import assemble_chrombpnet_fold_regions as assembler


class AssembleChromBPNetFoldRegionsTests(unittest.TestCase):
    def test_held_regions_are_fixed_ccres_and_training_peak_is_train_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            contract = root / "contract"
            contract.mkdir()
            with (contract / "crossed_outer_splits.tsv").open(
                "x", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=assembler.SPLIT_FIELDS, delimiter="\t"
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "split_id": "donor0_genomic0",
                        "donor_train_folds": "2,3,4",
                        "donor_valid_fold": "1",
                        "donor_test_fold": "0",
                        "genomic_train_folds": "2,3,4",
                        "genomic_valid_fold": "1",
                        "genomic_test_fold": "0",
                    }
                )
            fold_contigs = {
                0: "chr1,chr9,chr14,chr18,chr22",
                1: "chr2,chr8,chr13,chr19,chr20",
                2: "chr3,chr12,chr17,chrX,chrY",
                3: "chr4,chr7,chr11,chr15",
                4: "chr5,chr6,chr10,chr16,chr21",
            }
            with (contract / "genomic_folds.tsv").open(
                "x", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=assembler.GENOMIC_FOLD_FIELDS, delimiter="\t"
                )
                writer.writeheader()
                for fold, contigs in fold_contigs.items():
                    writer.writerow(
                        {"genomic_fold": fold, "contigs": contigs, "total_bp": 1}
                    )
            with (contract / "ccre_evaluation_windows.tsv").open(
                "x", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=assembler.CCRE_FIELDS, delimiter="\t"
                )
                writer.writeheader()
                for fold in range(5):
                    writer.writerow(
                        {
                            "contig": f"chr{fold + 1}",
                            "output_start": 1000,
                            "output_end": 2000,
                            "window_id": f"window_{fold}",
                            "genomic_fold": fold,
                            "window_class": "encode_ccre",
                            "ccre_class": "PLS",
                            "ccre_id": f"ccre_{fold}",
                            "ccre_start": 1200,
                            "ccre_end": 1800,
                            "input_start": 443,
                            "input_end": 2557,
                            "selection_hash": str(fold),
                        }
                    )
            training = root / "training.narrowPeak"
            training.write_text(
                "chr3\t2490\t2510\tpeak_train\t100\t.\t10\t5\t5\t10\n",
                encoding="utf-8",
            )
            output = root / "regions.bed"
            fold_json = root / "fold.json"
            summary_json = root / "summary.json"
            summary = assembler.assemble(
                training_peaks=training,
                split_contract=contract,
                split_id="donor0_genomic0",
                output_peaks=output,
                output_fold=fold_json,
                output_summary=summary_json,
                expected_ccre_per_fold=1,
            )
            with output.open(encoding="utf-8") as handle:
                rows = [line.rstrip("\n").split("\t") for line in handle]
            self.assertEqual([row[0] for row in rows], ["chr1", "chr2", "chr3"])
            self.assertTrue(rows[0][3].startswith("fixed_test|"))
            self.assertTrue(rows[1][3].startswith("fixed_valid|"))
            self.assertEqual(rows[2][3], "train_macs2|peak_train")
            self.assertEqual(summary["role_counts"], {"train": 1, "valid": 1, "test": 1})
            self.assertFalse(summary["held_atac_used_for_region_selection"])
            folds = json.loads(fold_json.read_text())
            self.assertIn("chr3", folds["train"])
            self.assertEqual(folds["valid"], fold_contigs[1].split(","))
            self.assertEqual(folds["test"], fold_contigs[0].split(","))


if __name__ == "__main__":
    unittest.main()
