from __future__ import annotations

from pathlib import Path
import tomllib
import unittest

from scripts.build_gse244832_label_free_atac_exchange import (
    DONORS,
    LABEL_CONTRACT,
    LabelFreeATACExchangeError,
    fold_index,
    read_label_free_obs,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/gse244832_label_free_atac_exchange.toml"


class GSE244832LabelFreeATACExchangeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_real_contract_binds_atac_only_nonchampion_exchange(self) -> None:
        validate_contract(ROOT, self.contract)
        self.assertFalse(self.contract["condition_labels_used_for_axis_selection"])
        self.assertFalse(self.contract["rna_assay_available_to_atac_model_jobs"])
        self.assertTrue(self.contract["claim_boundary"]["cross_assay_join_unresolved"])
        self.assertTrue(self.contract["claim_boundary"]["gse244832_is_non_champion"])

    def test_donor_fold_is_condition_free_and_deterministic(self) -> None:
        self.assertEqual(DONORS, tuple(f"D{index:02d}" for index in range(1, 19)))
        self.assertEqual(fold_index("gse244832", "D01"), fold_index("gse244832", "D01"))
        self.assertIn(fold_index("gse244832", "D18"), range(5))

    def test_condition_like_h5_field_is_rejected_before_open(self) -> None:
        with self.assertRaises(LabelFreeATACExchangeError):
            read_label_free_obs(Path("does-not-exist.h5ad"), ("condition",))

    def test_source_labels_have_explicit_roles(self) -> None:
        self.assertEqual(LABEL_CONTRACT["Fibroblasts"], ("fibroblast", "primary"))
        self.assertEqual(LABEL_CONTRACT["Low_confidence"], ("low_confidence", "excluded"))
        self.assertEqual(LABEL_CONTRACT["T_cells"], ("t_nk_cell", "secondary"))

    def test_fragment_geometry_is_dataset_native(self) -> None:
        exchange = self.contract["window_exchange"]
        self.assertEqual(exchange["query_fragment_cut_sites"], "start_plus_4_and_end_minus_5")
        self.assertFalse(exchange["cellranger_tn5_adjusted_coordinates"])


if __name__ == "__main__":
    unittest.main()
