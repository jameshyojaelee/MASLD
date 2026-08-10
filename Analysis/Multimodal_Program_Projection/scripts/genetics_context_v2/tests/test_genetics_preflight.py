#!/usr/bin/env python3

import csv
import importlib.util
import sys
import unittest
from collections import Counter
from pathlib import Path


TEST_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TEST_DIR.parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(SCRIPT_DIR))

from genetics_common import (  # noqa: E402
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    assert_candidate_root,
    clean,
    ensembl_base,
    parse_float,
)


def load_numbered(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


freeze = load_numbered("freeze_stage", "01_freeze_and_rederive.py")
hong = load_numbered("hong_source_audit", "07_audit_hong_public_contract.py")


def read_csv(name: str):
    with (TEST_DIR / "fixtures" / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


class GeneticsPreflightTests(unittest.TestCase):
    def test_missing_and_ensembl_normalization(self):
        self.assertEqual(clean("NA"), "")
        self.assertEqual(ensembl_base("ENSG000001.17"), "ENSG000001")
        self.assertEqual(parse_float("0"), 0.0)

    def test_candidate_path_is_exact_and_fail_closed(self):
        self.assertEqual(
            assert_candidate_root(PROJECT_ROOT, DEFAULT_CANDIDATE_ROOT),
            DEFAULT_CANDIDATE_ROOT.resolve(),
        )
        with self.assertRaises(ContractError):
            assert_candidate_root(PROJECT_ROOT, PROJECT_ROOT / "GWAS/finemapping/results")

    def test_bulk_zero_fdr_is_not_lost(self):
        collapsed, n_rows = freeze.collapse_bulk(read_csv("bulk_fixture.csv"))
        self.assertEqual(n_rows, 2)  # GENEA plus blank-symbol technical row
        self.assertTrue(collapsed["GENEA"]["established_state_associated"])
        self.assertFalse(collapsed["GENEB"]["established_state_associated"])
        self.assertEqual(collapsed["GENEA"]["bulk_treat_fdr"], "0")

    def test_susie_primary_does_not_fall_back_to_abf(self):
        collapsed = freeze.collapse_genetics(
            read_csv("tier_fixture.csv"), read_csv("full_fixture.csv")
        )
        self.assertTrue(collapsed["GENEA"]["primary_genetic"])
        self.assertFalse(collapsed["GENEB"]["primary_genetic"])
        self.assertFalse(collapsed["GENEC"]["primary_genetic"])

    def test_phenotype_catalog_is_one_row_per_primary_study(self):
        with (SCRIPT_DIR / "config" / "phenotype_source_catalog_v1.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 35)
        self.assertEqual(len({row["study_name"] for row in rows}), 35)
        self.assertEqual(
            Counter(row["phenotype_stratum"] for row in rows),
            Counter(
                {
                    "direct_masld_mash_diagnosis": 12,
                    "mri_pdff_or_histologic_steatosis": 3,
                    "alt_ast_or_ggt": 20,
                }
            ),
        )

    def test_broadaway_deposit_uses_joint_signal_p_for_reproduction(self):
        path = PROJECT_ROOT / "data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv"
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 9013)
        self.assertEqual(len({row["Ensembl"] for row in rows}), 6564)
        self.assertEqual(sum(row["Signal"].startswith("1:") for row in rows), 6564)
        self.assertEqual(sum(float(row["Joint_Pvalue"]) <= 1e-5 for row in rows), 9013)
        # Conditional refitting changes some Stepwise_Pvalue values; these rows
        # remain source-deposited signals because their joint p passes.
        self.assertEqual(sum(float(row["Stepwise_Pvalue"]) > 1e-5 for row in rows), 17)

    def test_quarantined_local_ieqtl_outputs_are_not_contract_inputs(self):
        text = (SCRIPT_DIR / "config" / "input_contract_v1.tsv").read_text()
        self.assertNotIn("ieqtl_disease_genes.csv", text)
        self.assertNotIn("ieqtl_summary.csv", text)

    def test_hong_variant_normalization_is_source_format_bounded(self):
        for row in read_csv("hong_variant_formats.csv"):
            if row["expected_status"] == "accepted":
                self.assertEqual(
                    hong.normalized_variant(row["source_value"]),
                    row["expected_normalized"],
                )
            else:
                self.assertEqual(row["expected_status"], "rejected")
                with self.assertRaises(ContractError):
                    hong.normalized_variant(row["source_value"])


if __name__ == "__main__":
    unittest.main()
