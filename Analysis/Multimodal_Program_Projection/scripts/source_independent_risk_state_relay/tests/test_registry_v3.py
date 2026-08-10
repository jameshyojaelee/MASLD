#!/usr/bin/env python3
"""Small deterministic fixtures for the Plan 45 registry-v3 logic."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))
SPEC = importlib.util.spec_from_file_location(
    "plan45_registry_v3", SCRIPT_ROOT / "02_build_corrected_registry_v3.py"
)
assert SPEC is not None and SPEC.loader is not None
REGISTRY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REGISTRY)
ORIENTATION_SPEC = importlib.util.spec_from_file_location(
    "plan45_orientation_worklist", SCRIPT_ROOT / "04_prepare_orientation_worklist.py"
)
assert ORIENTATION_SPEC is not None and ORIENTATION_SPEC.loader is not None
ORIENTATION = importlib.util.module_from_spec(ORIENTATION_SPEC)
ORIENTATION_SPEC.loader.exec_module(ORIENTATION)
FREEZE_SPEC = importlib.util.spec_from_file_location(
    "plan45_orientation_freeze", SCRIPT_ROOT / "07_validate_and_freeze_orientations.py"
)
assert FREEZE_SPEC is not None and FREEZE_SPEC.loader is not None
FREEZE = importlib.util.module_from_spec(FREEZE_SPEC)
FREEZE_SPEC.loader.exec_module(FREEZE)


class RegistryV3Tests(unittest.TestCase):
    def test_class_truth_table(self) -> None:
        bulk = {
            "CONV": {"established_state_associated": True},
            "GENETIC": {"established_state_associated": False},
            "STATE": {"established_state_associated": True},
            "NEITHER": {"established_state_associated": False},
            "BULKONLY": {"established_state_associated": False},
        }
        genetics = {
            "CONV": {"primary_genetic": True},
            "GENETIC": {"primary_genetic": True},
            "GENONLYMAP": {"primary_genetic": True},
            "STATE": {"primary_genetic": False},
            "NEITHER": {"primary_genetic": False},
        }
        observed = {
            row["gene_symbol"]: row["static_class"]
            for row in REGISTRY.build_classes(bulk, genetics)
        }
        self.assertEqual(observed["CONV"], "convergent")
        self.assertEqual(observed["GENETIC"], "genetic_only")
        self.assertEqual(observed["GENONLYMAP"], "indeterminate_not_jointly_testable")
        self.assertEqual(observed["STATE"], "disease_state_only")
        self.assertEqual(observed["NEITHER"], "neither")
        self.assertEqual(observed["BULKONLY"], "indeterminate_not_jointly_testable")

    def test_cross_study_physical_locus_and_representative(self) -> None:
        base = {
            "primary_genetic": "true",
            "static_class": "genetic_only",
            "driving_trait": "NAFLD",
            "driving_phenotype_stratum": "direct_masld_mash_diagnosis",
            "driving_ancestry": "EUR",
            "driving_regulatory_ancestry_status": "ancestry_matched_eur",
            "coloc_best_abf_pp4": "0.1",
            "driving_method": "susie",
        }
        rows = [
            {
                **base,
                "gene_symbol": "A",
                "ensembl_genetic": "ENSGA",
                "driving_ensembl": "ENSGA",
                "driving_gwas": "GWAS1",
                "driving_top_snp": "1:1000000",
                "coloc_best_susie_pp4": "0.7",
            },
            {
                **base,
                "gene_symbol": "B",
                "ensembl_genetic": "ENSGB",
                "driving_ensembl": "ENSGB",
                "driving_gwas": "GWAS2",
                "driving_top_snp": "1:1800000",
                "coloc_best_susie_pp4": "0.9",
            },
            {
                **base,
                "gene_symbol": "C",
                "ensembl_genetic": "ENSGC",
                "driving_ensembl": "ENSGC",
                "driving_gwas": "GWAS1",
                "driving_top_snp": "1:3000001",
                "coloc_best_susie_pp4": "0.8",
            },
        ]
        loci = REGISTRY.build_locus_registry(rows)
        by_gene = {row["gene_symbol"]: row for row in loci}
        self.assertEqual(by_gene["A"]["coarse_locus_uid"], by_gene["B"]["coarse_locus_uid"])
        self.assertNotEqual(by_gene["A"]["coarse_locus_uid"], by_gene["C"]["coarse_locus_uid"])
        self.assertEqual(by_gene["A"]["representative_gene"], "B")
        self.assertEqual(by_gene["B"]["is_representative"], "true")
        self.assertEqual(by_gene["A"]["orientation_status"], "requires_targeted_credible_set_pair_export")

    def test_top_snp_parser(self) -> None:
        self.assertEqual(REGISTRY.parse_top_snp("chr10:123:A:G"), ("10", 123))
        self.assertEqual(REGISTRY.parse_top_snp("10_456"), ("10", 456))
        self.assertIsNone(REGISTRY.parse_top_snp(""))

    def test_orientation_eligibility_is_fail_closed(self) -> None:
        row = {
            "pair_primary_genetic": "true",
            "is_representative_gene_in_locus": "true",
            "static_class": "genetic_only",
            "tier": "1",
            "phenotype_stratum": "direct_masld_mash_diagnosis",
            "regulatory_ancestry_status": "ancestry_matched_eur",
            "method": "susie",
            "susie_pp4": "0.8",
            "chromosome": "1",
        }
        self.assertTrue(ORIENTATION.eligible(row))
        for field, disallowed in [
            ("pair_primary_genetic", "false"),
            ("is_representative_gene_in_locus", "false"),
            ("static_class", "convergent"),
            ("tier", "2"),
            ("phenotype_stratum", "alt_ast_or_ggt"),
            ("regulatory_ancestry_status", "cross_ancestry_eqtl_limited"),
            ("method", "abf_fallback"),
            ("susie_pp4", "0.5"),
            ("chromosome", ""),
        ]:
            candidate = dict(row)
            candidate[field] = disallowed
            self.assertFalse(ORIENTATION.eligible(candidate), field)

    def test_pair_level_locus_registry_retains_all_gwas_for_orientation(self) -> None:
        classes = [
            {"gene_symbol": "A", "static_class": "genetic_only"},
            {"gene_symbol": "B", "static_class": "genetic_only"},
        ]
        base = {
            "trait": "NAFLD",
            "tier": "1",
            "phenotype_stratum": "direct_masld_mash_diagnosis",
            "ancestry": "EUR",
            "regulatory_ancestry_status": "ancestry_matched_eur",
            "source_dependence": "shared_eqtl_only",
            "abf_pp4": "0.1",
            "method": "susie",
            "top_snp_pp": "0.2",
            "pair_primary_genetic": "true",
            "pair_row_count": 1,
            "coarse_locus_uid": "",
            "static_class": "",
            "is_representative_gene_in_locus": "false",
            "representative_gene": "",
            "orientation_status": "requires_targeted_credible_set_pair_export",
        }
        pairs = [
            {
                **base,
                "gene_symbol": "A", "ensembl_id": "ENSGA", "gwas_name": "GWAS1",
                "susie_pp4": "0.9", "top_snp": "1:1000000", "chromosome": "1",
                "position": 1000000,
            },
            {
                **base,
                "gene_symbol": "A", "ensembl_id": "ENSGA", "gwas_name": "GWAS2",
                "susie_pp4": "0.8", "top_snp": "1:1200000", "chromosome": "1",
                "position": 1200000,
            },
            {
                **base,
                "gene_symbol": "B", "ensembl_id": "ENSGB", "gwas_name": "GWAS3",
                "susie_pp4": "0.7", "top_snp": "1:1300000", "chromosome": "1",
                "position": 1300000,
            },
        ]
        loci, retained = REGISTRY.build_pair_level_locus_registries(classes, pairs)
        self.assertEqual(len(retained), 3)
        self.assertEqual(len(loci), 2)
        a_locus = next(row for row in loci if row["gene_symbol"] == "A")
        self.assertEqual(a_locus["n_primary_gwas_for_gene_in_locus"], 2)
        self.assertEqual(a_locus["is_representative"], "true")
        self.assertEqual(
            sum(row["is_representative_gene_in_locus"] == "true" for row in retained),
            2,
        )

    def test_pair_family_requires_every_gwas_and_direction_concordance(self) -> None:
        concordant = [
            {
                "orientation_gate_pass": "true",
                "oriented_risk_effect": "risk_increases_expression",
            },
            {
                "orientation_gate_pass": "true",
                "oriented_risk_effect": "risk_increases_expression",
            },
        ]
        self.assertEqual(
            FREEZE.pair_family_verdict(concordant),
            (
                True,
                "all_eligible_pairs_oriented_and_concordant",
                "risk_increases_expression",
            ),
        )
        conflicting = [dict(row) for row in concordant]
        conflicting[1]["oriented_risk_effect"] = "risk_decreases_expression"
        self.assertEqual(FREEZE.pair_family_verdict(conflicting)[0], False)
        unresolved = [dict(row) for row in concordant]
        unresolved[1]["orientation_gate_pass"] = "false"
        self.assertEqual(
            FREEZE.pair_family_verdict(unresolved)[1],
            "one_or_more_eligible_pairs_failed_orientation",
        )

    def test_known_ukbb_pdff_overlap_collapses_to_one_evidence_family(self) -> None:
        first = REGISTRY.gwas_evidence_family(
            "PDFF_A",
            {
                "phenotype_stratum": "mri_pdff_or_histologic_steatosis",
                "source_note": "UK Biobank imaging cohort; overlaps other PDFF studies.",
            },
        )
        second = REGISTRY.gwas_evidence_family(
            "PDFF_B",
            {
                "phenotype_stratum": "mri_pdff_or_histologic_steatosis",
                "source_note": "UK Biobank imaging cohort; overlaps other PDFF studies.",
            },
        )
        self.assertEqual(first[0], second[0])
        self.assertTrue(first[2])
        meta = REGISTRY.gwas_evidence_family(
            "META",
            {
                "phenotype_stratum": "direct_masld_mash_diagnosis",
                "source_note": "Meta-analysis includes participant pools reused by other studies.",
            },
        )
        self.assertFalse(meta[2])


if __name__ == "__main__":
    unittest.main()
