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

    def test_canonical_rule_uses_padj_and_logfc(self):
        rows = [
            {"gene": "ENSG000001.1", "symbol": "GENEA", "logFC": "0.6", "t": "3", "padj": "0.06", "treat_fdr": "0.01"},
            {"gene": "ENSG000001.2", "symbol": "GENEA", "logFC": "0.7", "t": "2", "padj": "0.04", "treat_fdr": "0.5"},
            {"gene": "ENSG000002.1", "symbol": "GENEB", "logFC": "0.4", "t": "5", "padj": "0.001", "treat_fdr": "0.01"},
        ]
        canonical, n_rows = freeze.collapse_bulk(rows, "canonical")
        self.assertEqual(n_rows, 1)
        self.assertTrue(canonical["GENEA"]["established_state_associated"])
        self.assertEqual(canonical["GENEA"]["bulk_padj"], "0.04")  # displayed row = lowest padj
        self.assertFalse(canonical["GENEB"]["established_state_associated"])  # |logFC| <= 0.5
        treat, n_treat = freeze.collapse_bulk(rows, "treat")
        self.assertEqual(n_treat, 2)
        self.assertTrue(treat["GENEB"]["established_state_associated"])
        self.assertEqual(treat["GENEA"]["bulk_treat_fdr"], "0.01")

    def test_primary_assembly_row_wins_over_alt_contig(self):
        # HLA-DRA-like: the alternate-haplotype row has the lower padj, the
        # primary chr6 row is kept and alone gives the call and the identity.
        rows = [
            {"gene": "ENSG000000PRI.1", "symbol": "GENEA", "logFC": "0.6", "t": "2", "padj": "0.060", "treat_fdr": "0.2"},
            {"gene": "ENSG000000ALT.1", "symbol": "GENEA", "logFC": "0.6", "t": "3", "padj": "0.048", "treat_fdr": "0.1"},
            {"gene": "ENSG000000AL1.1", "symbol": "GENEB", "logFC": "0.9", "t": "4", "padj": "0.01", "treat_fdr": "0.3"},
            {"gene": "ENSG000000AL2.1", "symbol": "GENEB", "logFC": "0.9", "t": "1", "padj": "0.20", "treat_fdr": "0.3"},
        ]
        primary = {"ENSG000000PRI"}
        collapsed, _ = freeze.collapse_bulk(rows, "canonical", primary)
        self.assertFalse(collapsed["GENEA"]["established_state_associated"])
        self.assertEqual(collapsed["GENEA"]["bulk_padj"], "0.060")
        self.assertEqual(collapsed["GENEA"]["ensembl_bulk"], "ENSG000000PRI")
        # No primary row: the lowest-q fallback, with every ID kept for identity.
        self.assertTrue(collapsed["GENEB"]["established_state_associated"])
        self.assertEqual(collapsed["GENEB"]["ensembl_bulk"], "ENSG000000AL1;ENSG000000AL2")
        # SORD2P-like: the alt row has the lower padj but fails; the primary row is a DEG.
        sord = [
            {"gene": "ENSG000000PRI.2", "symbol": "GENEC", "logFC": "0.8", "t": "3", "padj": "0.01", "treat_fdr": "0.2"},
            {"gene": "ENSG000000ALT.2", "symbol": "GENEC", "logFC": "0.1", "t": "5", "padj": "0.001", "treat_fdr": "0.2"},
        ]
        collapsed, _ = freeze.collapse_bulk(sord, "canonical", {"ENSG000000PRI"})
        self.assertTrue(collapsed["GENEC"]["established_state_associated"])

    def test_renamed_symbol_joins_by_ensembl(self):
        # COLOC still calls ENSG000001 by its old symbol; the bulk table uses the
        # GENCODE v49 name. The gene is keyed on the bulk symbol, not split in two.
        tier = [{"ensembl": "ENSG000001", "coloc_best_susie_pp4": "0.8", "coloc_best_abf_pp4": "0.7",
                 "driving_gwas": "DIRECT", "driving_trait": "MASLD"}]
        full = [{"gene": "C11orf80", "ensembl": "ENSG000001.3", "coloc_best_susie_pp4": "0.8"}]
        renamed = freeze.collapse_genetics(tier, full, None, {"ENSG000001": "TOP6BL"})
        self.assertEqual(set(renamed), {"TOP6BL"})
        self.assertEqual(set(freeze.collapse_genetics(tier, full)), {"C11orf80"})

    def test_untestable_susie_is_its_own_state_not_a_negative(self):
        # GENEB: SuSiE 0 elsewhere but an untestable pair; GENEC: no SuSiE PP4, ABF 0.9,
        # untestable pair; GENEA: SuSiE-positive, so an untestable pair does not matter.
        untestable = {"ENSG000001", "ENSG000002", "ENSG000003"}
        collapsed = freeze.collapse_genetics(
            read_csv("tier_fixture.csv"), read_csv("full_fixture.csv"), untestable
        )
        self.assertFalse(collapsed["GENEA"]["genetic_untestable"])
        self.assertTrue(collapsed["GENEB"]["genetic_untestable"])
        self.assertTrue(collapsed["GENEC"]["genetic_untestable"])
        legacy = freeze.collapse_genetics(read_csv("tier_fixture.csv"), read_csv("full_fixture.csv"))
        self.assertNotIn("genetic_untestable", legacy["GENEB"])

    def test_untestable_pairs_are_read_from_tier12_studies_only(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            tier = Path(tmp) / "tier.tsv"
            tier.write_text(
                "study_name\ttier\nT1\t1\nT3\t3\n", encoding="utf-8"
            )
            master = Path(tmp) / "master.csv"
            master.write_text(
                "gwas_name,gene,ensembl,method\n"
                f"T1,GENEB,ENSG000002.4,{freeze.UNTESTABLE_METHOD}\n"
                f"T3,GENED,ENSG000009.1,{freeze.UNTESTABLE_METHOD}\n"
                "T1,GENEA,ENSG000001.1,susie\n",
                encoding="utf-8",
            )
            self.assertEqual(freeze.untestable_ensembl(master, tier), {"ENSG000002"})

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
