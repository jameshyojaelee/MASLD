from __future__ import annotations

import json
from pathlib import Path
import tempfile
import tomllib
import unittest

from scripts.build_gse281367_label_free_atac_exchange import (
    LabelFreeATACExchangeError,
    fold_index,
    read_label_free_obs,
    read_windows,
    resolve_raw_barcodes,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/gse281367_label_free_atac_exchange.toml"


class GSE281367LabelFreeATACExchangeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_real_contract_binds_nonchampion_label_free_exchange(self) -> None:
        validate_contract(ROOT, self.contract)
        self.assertFalse(self.contract["condition_labels_used_for_axis_selection"])
        self.assertTrue(self.contract["claim_boundary"]["gse281367_is_non_champion"])
        self.assertFalse(self.contract["claim_boundary"]["external_confirmation_claim_allowed"])

    def test_donor_fold_is_condition_free_and_deterministic(self) -> None:
        self.assertEqual(
            fold_index("gse281367", "Z01"),
            fold_index("gse281367", "Z01"),
        )
        self.assertIn(fold_index("gse281367", "Z12"), range(5))

    def test_condition_like_h5_field_is_rejected_before_open(self) -> None:
        with self.assertRaises(LabelFreeATACExchangeError):
            read_label_free_obs(Path("does-not-exist.h5ad"), ("condition",))

    def test_barcode_resolution_is_bijective(self) -> None:
        self.assertEqual(
            resolve_raw_barcodes(["AAAC-1-2", "TTTT-1-7"], ["AAAC-1", "TTTT-1"]),
            ["AAAC-1", "TTTT-1"],
        )
        with self.assertRaises(LabelFreeATACExchangeError):
            resolve_raw_barcodes(["AAAC-1"], ["AAAC-1", "TTTT-1"])

    def test_window_roles_must_be_whole_contig_disjoint(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "windows.tsv"
            rows = [
                f"{index}\tw{index}\t{'valid' if index < 16000 else 'test'}\t"
                f"{'chr1' if index < 16000 else 'chr2'}\t{index * 1000}\t{index * 1000 + 1000}\n"
                for index in range(32000)
            ]
            path.write_text(
                "window_index\twindow_id\trole\tcontig\tstart\tend\n" + "".join(rows),
                encoding="utf-8",
            )
            self.assertEqual(len(read_windows(path)), 32000)
            rows[-1] = rows[-1].replace("chr2", "chr1")
            path.write_text(
                "window_index\twindow_id\trole\tcontig\tstart\tend\n" + "".join(rows),
                encoding="utf-8",
            )
            with self.assertRaises(LabelFreeATACExchangeError):
                read_windows(path)

    def test_observed_and_sequence_families_remain_distinct(self) -> None:
        by_id = {row["model_id"]: row for row in self.contract["family_eligibility"]}
        self.assertEqual(by_id["scbasset"]["input_regime"], "sequence_only")
        self.assertFalse(by_id["scbasset"]["query_atac_allowed"])
        self.assertEqual(by_id["epibert"]["input_regime"], "observed_atac")
        self.assertTrue(by_id["epibert"]["query_atac_allowed"])


if __name__ == "__main__":
    unittest.main()
