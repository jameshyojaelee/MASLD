#!/usr/bin/env python3
"""Focused tests for candidate-only Gene Catalog observability matching."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

import gene_catalog_observability as obs


def gene_fixture(n_per_group: int = 40) -> pd.DataFrame:
    rows = []
    for index in range(n_per_group):
        abundance = 1.0 + index / 10
        variability = 0.5 + index / 20
        gene_length = 1000 + 100 * index
        detection = 1 + index % 5
        rows.append(
            {
                "ensembl_id": f"ENSG_L{index:05d}",
                "gene_biotype": "lncRNA",
                "bulk_abundance": abundance,
                "expression_variability": variability,
                "gene_length_bp": gene_length,
                "cohort_detection_count": detection,
            }
        )
        rows.append(
            {
                "ensembl_id": f"ENSG_P{index:05d}",
                "gene_biotype": "protein_coding",
                "bulk_abundance": abundance,
                "expression_variability": variability,
                "gene_length_bp": gene_length,
                "cohort_detection_count": detection,
            }
        )
    return pd.DataFrame(rows)


def coverage_fixture(genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for position, row in enumerate(
        genes.sort_values("ensembl_id").itertuples(index=False)
    ):
        is_lnc = row.gene_biotype == "lncRNA"
        rows.extend(
            [
                {
                    "ensembl_id": row.ensembl_id,
                    "assay_id": "bulk-five-cohort",
                    "assay_family": "bulk_rna",
                    "molecular_scope": "rna_molecule",
                    "applicability_state": "applicable",
                    "observable": position % 3 != 0 if is_lnc else position % 4 != 0,
                },
                {
                    "ensembl_id": row.ensembl_id,
                    "assay_id": "liver-promoter-atac",
                    "assay_family": "chromatin",
                    "molecular_scope": "regulatory_locus",
                    "applicability_state": "applicable",
                    "observable": position % 2 == 0,
                },
                {
                    "ensembl_id": row.ensembl_id,
                    "assay_id": "liver-proteome",
                    "assay_family": "proteomics",
                    "molecular_scope": "protein_molecule",
                    "applicability_state": "not_applicable" if is_lnc else "applicable",
                    "observable": pd.NA if is_lnc else position % 2 == 0,
                },
            ]
        )
    return pd.DataFrame(rows)


class ObservabilityTests(unittest.TestCase):
    def test_matching_is_deterministic_and_without_replacement(self):
        genes = gene_fixture()
        first = obs.match_genes(genes)
        second = obs.match_genes(genes.sample(frac=1, random_state=99))
        pd.testing.assert_frame_equal(first.pairs, second.pairs)
        self.assertFalse(first.pairs["lncrna_ensembl_id"].duplicated().any())
        self.assertFalse(first.pairs["protein_coding_ensembl_id"].duplicated().any())
        self.assertEqual(first.gate_status, "pass")
        self.assertEqual(first.match_rate, 1.0)
        self.assertLessEqual(first.max_abs_matched_smd, 0.10)

    def test_caliper_is_exactly_point_two_logit_sd(self):
        result = obs.match_genes(gene_fixture())
        logits = result.eligible_genes["propensity_logit"].to_numpy()
        self.assertAlmostEqual(result.caliper, 0.2 * np.std(logits, ddof=1), places=12)

    def test_balance_reports_all_four_covariates(self):
        result = obs.match_genes(gene_fixture())
        self.assertEqual(set(result.balance["stage"]), {"before", "matched"})
        self.assertEqual(set(result.balance["covariate"]), set(obs.MATCH_COVARIATES))

    def test_raw_matched_differences_and_bootstrap_are_deterministic(self):
        genes = gene_fixture()
        coverage = coverage_fixture(genes)
        _, first, counts = obs.compare_observability(
            genes, coverage, bootstrap_replicates=500
        )
        _, second, _ = obs.compare_observability(
            genes, coverage.sample(frac=1, random_state=2), bootstrap_replicates=500
        )
        pd.testing.assert_frame_equal(first, second)
        bulk = first.set_index("assay_id").loc["bulk-five-cohort"]
        self.assertEqual(bulk["comparison_status"], "matched_comparison")
        self.assertTrue(np.isfinite(bulk["raw_difference_percentage_points"]))
        self.assertTrue(np.isfinite(bulk["matched_difference_percentage_points"]))
        self.assertTrue(np.isfinite(bulk["bootstrap_ci_lower_percentage_points"]))
        self.assertFalse(counts.empty)

    def test_proteomics_is_not_applicable_for_lncrna(self):
        genes = gene_fixture()
        _, comparisons, counts = obs.compare_observability(
            genes, coverage_fixture(genes), bootstrap_replicates=500
        )
        protein = comparisons.set_index("assay_id").loc["liver-proteome"]
        self.assertEqual(protein["comparison_status"], "not_applicable")
        self.assertTrue(np.isnan(protein["matched_difference_percentage_points"]))
        self.assertTrue(np.isfinite(protein["raw_protein_coding_percent"]))
        lnc_status = counts[
            (counts["assay_id"] == "liver-proteome")
            & (counts["gene_biotype"] == "lncRNA")
        ]
        self.assertEqual(set(lnc_status["applicability_state"]), {"not_applicable"})

    def test_chromatin_is_locus_observability(self):
        genes = gene_fixture()
        _, comparisons, _ = obs.compare_observability(
            genes, coverage_fixture(genes), bootstrap_replicates=500
        )
        chromatin = comparisons.set_index("assay_id").loc["liver-promoter-atac"]
        self.assertEqual(chromatin["molecular_scope"], "regulatory_locus")
        self.assertEqual(
            chromatin["interpretation"], "locus_observability_not_lncrna_molecule"
        )

    def test_wrong_chromatin_scope_fails(self):
        coverage = coverage_fixture(gene_fixture())
        coverage.loc[coverage["assay_family"] == "chromatin", "molecular_scope"] = (
            "rna_molecule"
        )
        with self.assertRaisesRegex(obs.ObservabilityContractError, "CHROMATIN_SCOPE"):
            obs.validate_coverage(coverage)

    def test_lncrna_proteomics_applicable_fails(self):
        genes = gene_fixture()
        coverage = coverage_fixture(genes)
        lnc_id = genes.loc[genes["gene_biotype"] == "lncRNA", "ensembl_id"].iloc[0]
        target = (coverage["ensembl_id"] == lnc_id) & (
            coverage["assay_family"] == "proteomics"
        )
        coverage.loc[target, "applicability_state"] = "applicable"
        coverage.loc[target, "observable"] = True
        with self.assertRaisesRegex(
            obs.ObservabilityContractError, "LNCRNA_PROTEOMICS"
        ):
            obs.compare_observability(genes, coverage, bootstrap_replicates=500)

    def test_unresolved_status_does_not_become_false_observability(self):
        coverage = coverage_fixture(gene_fixture())
        target = coverage["assay_id"] == "bulk-five-cohort"
        first = coverage.index[target][0]
        coverage.loc[first, "applicability_state"] = "indeterminate"
        coverage.loc[first, "observable"] = False
        with self.assertRaisesRegex(
            obs.ObservabilityContractError, "OBSERVABILITY_WITHOUT_APPLICABILITY"
        ):
            obs.validate_coverage(coverage)

    def test_incomplete_coverage_matrix_fails(self):
        genes = gene_fixture()
        coverage = coverage_fixture(genes).iloc[1:].copy()
        with self.assertRaisesRegex(
            obs.ObservabilityContractError, "INCOMPLETE_COVERAGE_MATRIX"
        ):
            obs.compare_observability(genes, coverage, bootstrap_replicates=500)

    def test_raw_universe_keeps_gene_with_missing_matching_covariate(self):
        genes = gene_fixture()
        coverage = coverage_fixture(genes)
        genes.loc[0, "bulk_abundance"] = np.nan
        result, comparisons, _ = obs.compare_observability(
            genes, coverage, bootstrap_replicates=500
        )
        bulk = comparisons.set_index("assay_id").loc["bulk-five-cohort"]
        self.assertEqual(
            bulk["n_raw_lncrna"],
            int((genes["gene_biotype"] == "lncRNA").sum()),
        )
        self.assertEqual(
            len(result.eligible_genes),
            len(genes) - 1,
        )

    def test_matching_gate_row_tracks_result(self):
        result = obs.match_genes(gene_fixture())
        gate = obs.matching_gate_row(result)
        self.assertEqual(gate["n_matched_pairs"], len(result.pairs))
        self.assertEqual(gate["promotion_allowed"], result.gate_status == "pass")

    def test_no_caliper_matches_is_valid_failed_gate(self):
        genes = gene_fixture()
        lnc = genes["gene_biotype"] == "lncRNA"
        genes.loc[lnc, "bulk_abundance"] += 100.0
        result = obs.match_genes(genes)
        self.assertEqual(result.gate_status, "fail")
        self.assertEqual(len(result.pairs), 0)
        _, comparisons, _ = obs.compare_observability(
            genes, coverage_fixture(genes), bootstrap_replicates=500
        )
        bulk = comparisons.set_index("assay_id").loc["bulk-five-cohort"]
        self.assertEqual(bulk["comparison_status"], "matching_gate_failed")
        self.assertTrue(np.isnan(bulk["matched_difference_percentage_points"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
