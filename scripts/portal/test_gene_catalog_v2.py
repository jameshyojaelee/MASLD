#!/usr/bin/env python3
"""Focused unit tests for the noncoding-aware MASLD Gene Catalog v2 contract."""

from __future__ import annotations

import json
import unittest

import gene_catalog_v2 as catalog


def base_v2(**updates):
    record = {
        "ensembl_id": "ENSG00000000001",
        "symbol": "GENE1",
        "primary_evidence_class": "established_state_associated",
        "gene_biotype": "protein_coding",
        "candidate_object": "protein_coding_gene",
        "lncrna_genomic_class": "not_applicable",
        "annotation_release": "GENCODE v49",
        "strand_audit_status": "passed_reverse_stranded",
        "mapping_status": "unique_stable_ensembl",
        "noncoding_dna_context": "not_applicable",
        "credible_set_id": "",
        "credible_set_coding_pip_mass": None,
        "credible_set_noncoding_pip_mass": None,
        "credible_set_unresolved_pip_mass": None,
        "credible_set_to_gene_link_status": "not_applicable",
        "target_link_basis": "none",
        "assay_applicability": {"bulk_rna": "applicable", "proteomics": "applicable"},
        "open_mechanistic_question": "protein_coding_state_function",
    }
    record.update(updates)
    return record


class GeneCatalogV2Tests(unittest.TestCase):
    def test_legacy_record_retains_v1_route(self):
        adapted = catalog.adapt_record(
            {
                "ensembl_id": "ENSG00000000001",
                "symbol": "GENE1",
                "primary_evidence_class": "genetically_anchored",
                "next_experiment_rule_id": "EXP_ALLELE_AWARE_V1",
            }
        )
        self.assertEqual(
            adapted["catalog_schema_version"], catalog.LEGACY_SCHEMA_VERSION
        )
        self.assertEqual(
            adapted["recommended_experiment_rule_id"], "EXP_ALLELE_AWARE_V1"
        )
        self.assertTrue(set(catalog.V2_FIELDS).issubset(adapted))

    def test_regulatory_dna_route(self):
        record = base_v2(
            gene_biotype="not_applicable",
            candidate_object="noncoding_regulatory_element",
            strand_audit_status="not_applicable",
            noncoding_dna_context="liver_lineage_accessible_chromatin;liver_abc_enhancer",
            credible_set_id="CS1",
            credible_set_coding_pip_mass=0.05,
            credible_set_noncoding_pip_mass=0.90,
            credible_set_unresolved_pip_mass=0.05,
            credible_set_to_gene_link_status="source_qualified_context",
            target_link_basis="abc_enhancer_gene",
            assay_applicability={
                "chromatin": "applicable",
                "proteomics": "not_applicable",
            },
            open_mechanistic_question="regulatory_element_to_target",
        )
        adapted = catalog.adapt_record(record)
        self.assertEqual(
            adapted["recommended_experiment_rule_id"], "EXP_REGULATORY_DNA_V2"
        )
        self.assertEqual(
            adapted["noncoding_dna_context"],
            "liver_abc_enhancer;liver_lineage_accessible_chromatin",
        )

    def test_disease_state_lncrna_routes_to_rna_depletion(self):
        record = base_v2(
            gene_biotype="lncRNA",
            candidate_object="lncrna_transcript",
            lncrna_genomic_class="intergenic",
            noncoding_dna_context="not_applicable",
            credible_set_to_gene_link_status="not_applicable",
            assay_applicability={
                "bulk_rna": "applicable",
                "proteomics": "not_applicable",
            },
            open_mechanistic_question="lncrna_rna_product",
        )
        card = catalog.experiment_for_record(record)
        self.assertEqual(card["experiment_rule_id"], "EXP_LNCRNA_RNA_PRODUCT_V2")
        self.assertIn("RNA depletion", card["perturbation"])

    def test_colocalized_lncrna_routes_to_disambiguation(self):
        record = base_v2(
            primary_evidence_class="genetically_anchored",
            gene_biotype="lncRNA",
            candidate_object="lncrna_transcript",
            lncrna_genomic_class="opposite_strand_promoter_proximal",
            noncoding_dna_context="promoter_proximal",
            credible_set_to_gene_link_status="colocalized_egene",
            target_link_basis="coloc_susie",
            assay_applicability={
                "bulk_rna": "applicable",
                "proteomics": "not_applicable",
            },
            open_mechanistic_question="lncrna_dna_transcription_or_rna",
        )
        card = catalog.experiment_for_record(record)
        self.assertEqual(
            card["experiment_rule_id"], "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2"
        )
        self.assertIn("paired RNA depletion", card["perturbation"])
        self.assertIn("locus-level CRISPRi", card["perturbation"])

    def test_colocalized_lncrna_cannot_route_as_rna_product_only(self):
        with self.assertRaisesRegex(catalog.GeneCatalogV2Error, "LNCRNA_MECHANISM"):
            catalog.adapt_record(
                base_v2(
                    primary_evidence_class="genetic_only",
                    gene_biotype="lncRNA",
                    candidate_object="lncrna_transcript",
                    lncrna_genomic_class="intergenic",
                    credible_set_to_gene_link_status="colocalized_egene",
                    target_link_basis="coloc_susie",
                    assay_applicability={"proteomics": "not_applicable"},
                    open_mechanistic_question="lncrna_rna_product",
                )
            )

    def test_ambiguous_lncrna_locus_has_same_disambiguation_route(self):
        adapted = catalog.adapt_record(
            base_v2(
                gene_biotype="lncRNA",
                candidate_object="ambiguous_lncrna_locus",
                lncrna_genomic_class="complex",
                mapping_status="mapping_ambiguous",
                noncoding_dna_context="unresolved",
                credible_set_to_gene_link_status="unknown",
                assay_applicability={
                    "bulk_rna": "untestable",
                    "proteomics": "not_applicable",
                },
                open_mechanistic_question="lncrna_dna_transcription_or_rna",
            )
        )
        self.assertEqual(
            adapted["recommended_experiment_rule_id"],
            "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2",
        )

    def test_protein_coding_state_route(self):
        adapted = catalog.adapt_record(base_v2())
        self.assertEqual(
            adapted["recommended_experiment_rule_id"], "EXP_PROTEIN_STATE_CONTEXT_V2"
        )

    def test_untestable_route_has_priority(self):
        adapted = catalog.adapt_record(
            base_v2(
                candidate_object="untestable_nomination",
                gene_biotype="unknown",
                strand_audit_status="not_audited",
                mapping_status="mapping_ambiguous",
                assay_applicability={"bulk_rna": "untestable"},
                open_mechanistic_question="missing_decisive_assay",
            )
        )
        self.assertEqual(
            adapted["recommended_experiment_rule_id"], "EXP_MEASURE_MISSING_ASSAY_V2"
        )

    def test_lncrna_proteomics_is_not_applicable(self):
        with self.assertRaisesRegex(catalog.GeneCatalogV2Error, "LNCRNA_PROTEOMICS"):
            catalog.adapt_record(
                base_v2(
                    gene_biotype="lncRNA",
                    candidate_object="lncrna_transcript",
                    lncrna_genomic_class="intergenic",
                    assay_applicability={"proteomics": "untestable"},
                    open_mechanistic_question="lncrna_rna_product",
                )
            )

    def test_assay_applicability_is_canonical_json(self):
        adapted = catalog.adapt_record(
            base_v2(
                assay_applicability={
                    "proteomics": "applicable",
                    "bulk_rna": "applicable",
                }
            )
        )
        self.assertEqual(
            adapted["assay_applicability"],
            json.dumps(
                {"bulk_rna": "applicable", "proteomics": "applicable"},
                sort_keys=True,
                separators=(",", ":"),
            ),
        )

    def test_pip_mass_contract(self):
        with self.assertRaisesRegex(catalog.GeneCatalogV2Error, "PIP_MASS_SUM"):
            catalog.adapt_record(
                base_v2(
                    credible_set_id="CS1",
                    credible_set_coding_pip_mass=0.2,
                    credible_set_noncoding_pip_mass=0.7,
                    credible_set_unresolved_pip_mass=0.2,
                )
            )

    def test_supplied_route_cannot_override_rule(self):
        with self.assertRaisesRegex(catalog.GeneCatalogV2Error, "ROUTE_MUTATION"):
            catalog.adapt_record(
                base_v2(recommended_experiment_rule_id="EXP_REGULATORY_DNA_V2")
            )

    def test_object_question_mismatch_is_rejected(self):
        with self.assertRaisesRegex(catalog.GeneCatalogV2Error, "OBJECT_QUESTION"):
            catalog.adapt_record(
                base_v2(open_mechanistic_question="regulatory_element_to_target")
            )

    def test_cas13_screen_fields_are_rejected(self):
        with self.assertRaisesRegex(
            catalog.GeneCatalogV2Error, "CAS13_SCREEN_FIREWALL"
        ):
            catalog.adapt_record(base_v2(cas13_screen_hit=True))
        with self.assertRaisesRegex(
            catalog.GeneCatalogV2Error, "CAS13_SCREEN_FIREWALL"
        ):
            catalog.adapt_record(
                base_v2(limitation="selected from Cas13 screen rank 3")
            )
        with self.assertRaisesRegex(
            catalog.GeneCatalogV2Error, "CAS13_SCREEN_FIREWALL"
        ):
            catalog.adapt_record(base_v2(metadata={"cas13_screen_result": "positive"}))

    def test_generic_cas13_method_is_allowed(self):
        record = base_v2(
            gene_biotype="lncRNA",
            candidate_object="lncrna_transcript",
            lncrna_genomic_class="intergenic",
            assay_applicability={"proteomics": "not_applicable"},
            open_mechanistic_question="lncrna_rna_product",
        )
        self.assertIn("Cas13", catalog.experiment_for_record(record)["perturbation"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
