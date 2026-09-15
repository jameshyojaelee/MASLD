from __future__ import annotations

import copy
import csv
import gzip
import io
import json
from pathlib import Path
import tomllib
import unittest

from scripts.audit_gse281367_activation_readiness import (
    ActivationReadinessError,
    parse_bioc,
    parse_crossref,
    parse_geo_soft,
    parse_sra_runinfo,
    validate_contract,
    validate_project_artifacts,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/gse281367_activation_readiness.toml"


def geo_fixture() -> bytes:
    lines = [
        "^SERIES = GSE281367",
        "!Series_supplementary_file = ftp://example/GSE281367_RAW.tar",
        "!Series_supplementary_file = ftp://example/GSE281367_seurat_clustered.rds.gz",
    ]
    gsm = 8619363
    for condition in ("MASH", "Normal"):
        for replicate in range(1, 7):
            prefix = f"GSM{gsm}_{condition}_rep{replicate}"
            lines.extend(
                [
                    f"^SAMPLE = GSM{gsm}",
                    f"!Sample_title = {condition}, replicate{replicate}, snATAC-seq",
                    "!Sample_library_strategy = ATAC-seq",
                    f"!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN{gsm}",
                    f"!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX{gsm}",
                    f"!Sample_supplementary_file_1 = ftp://example/{prefix}_barcodes.tsv.gz",
                    f"!Sample_supplementary_file_2 = ftp://example/{prefix}_fragments.tsv.gz",
                    f"!Sample_supplementary_file_3 = ftp://example/{prefix}_matrix.mtx.gz",
                    f"!Sample_supplementary_file_4 = ftp://example/{prefix}_peaks.bed.gz",
                ]
            )
            gsm += 1
    return gzip.compress(("\n".join(lines) + "\n").encode())


def sra_fixture(records: list[dict[str, object]]) -> str:
    fields = [
        "Run", "Experiment", "BioSample", "LibraryStrategy", "LibraryLayout",
        "Consent", "Sex", "Disease", "Histological_Type", "size_MB", "spots", "bases",
    ]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for index, record in enumerate(records, start=1):
        writer.writerow(
            {
                "Run": f"SRR{index:08d}",
                "Experiment": record["experiment"],
                "BioSample": record["biosample"],
                "LibraryStrategy": "ATAC-seq",
                "LibraryLayout": "PAIRED",
                "Consent": "public",
                "Sex": "",
                "Disease": "",
                "Histological_Type": "",
                "size_MB": "10",
                "spots": "100",
                "bases": "200",
            }
        )
    return stream.getvalue()


def bioc_fixture() -> bytes:
    text = " ".join(
        [
            "liver nuclei from 12 individuals (6 Normal and 6 MASLD)",
            "69,595 high-quality nuclei",
            "Samples from different patients were processed separately",
            "Normal and MASLD human livers for snATAC-seq were provided by the Liver Tissue Cell Distribution System",
            "Cell Ranger ATAC v2.1.0",
            "downloaded from GSE189600",
            "GSE281367, GSE281364, and GSE281160",
        ]
    )
    value = [
        {
            "documents": [
                {
                    "id": "PMC12633503",
                    "infons": {"license": "CC BY"},
                    "passages": [{"text": text}],
                }
            ]
        }
    ]
    return json.dumps(value).encode()


class GSE281367ActivationReadinessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_participant_and_atac_only_topology_are_binding(self) -> None:
        validate_contract(self.contract)
        self.assertEqual(self.contract["declared_participant_n"], 12)
        self.assertEqual(self.contract["pairing_topology"], "same_study_unpaired_atac_only")
        self.assertFalse(self.contract["rna_observed"])
        self.assertFalse(self.contract["rna_invention_allowed"])
        self.assertTrue(self.contract["assay_native_donor_grouped_split_available"])

    def test_geo_mash_and_article_masld_conflict_cannot_be_erased(self) -> None:
        phenotypes, _, _ = validate_contract(self.contract)
        by_id = {row["field_id"]: row for row in phenotypes}
        self.assertEqual(
            self.contract["source_label_semantics"],
            "CONFLICT_GEO_MASH_ARTICLE_MASLD",
        )
        self.assertEqual(
            by_id["source_condition"]["missingness_state"],
            "observed_conflicting_semantics",
        )
        self.assertFalse(self.contract["geo_mash_label_is_histology_adjudicated"])
        self.assertFalse(self.contract["normal_label_is_adjudicated_masld_negative"])

    def test_histology_and_clinical_metadata_remain_structurally_missing(self) -> None:
        phenotypes, _, _ = validate_contract(self.contract)
        by_id = {row["field_id"]: row for row in phenotypes}
        for field_id in (
            "mash_adjudication", "masld_diagnostic_criteria", "fibrosis_stage",
            "nas_total", "histology_components", "sex", "age", "bmi",
            "alcohol_and_other_etiology_exclusions",
        ):
            self.assertEqual(by_id[field_id]["missingness_state"], "structurally_missing")

    def test_family_native_tasks_exclude_rna_histology_and_champion_claims(self) -> None:
        _, tasks, _ = validate_contract(self.contract)
        by_id = {row["task_id"]: row for row in tasks}
        self.assertEqual(by_id["atac_fixed_window_donor_lineage_transport"]["eligibility"], "active_development")
        self.assertEqual(by_id["observed_atac_profile_and_count_tournament"]["eligibility"], "conditional")
        for task_id in (
            "rna_conditioned_atac_or_trans_expression_context",
            "paired_rna_atac_integration_or_translation",
            "histology_or_clinical_metadata_conditioning",
            "external_sealed_universal_or_masld_champion",
        ):
            self.assertEqual(by_id[task_id]["eligibility"], "prohibited")

    def test_public_geo_and_sra_fixtures_resolve_12_individual_atac_records(self) -> None:
        records, series_files = parse_geo_soft(geo_fixture())
        self.assertEqual(len(records), 12)
        self.assertTrue(all(row["rna_state"] == "structurally_missing" for row in records))
        self.assertEqual(sum(len(row["processed_files"]) for row in records), 48)
        self.assertEqual(len(series_files), 2)
        rows, summary = parse_sra_runinfo(sra_fixture(records), records)
        self.assertEqual(len(rows), 12)
        self.assertTrue(summary["sex_disease_histology_fields_blank"])
        self.assertFalse(summary["raw_reads_downloaded"])

    def test_bioc_and_crossref_preserve_publication_dependencies(self) -> None:
        article = parse_bioc(bioc_fixture())
        self.assertEqual(article["participant_n"], 12)
        self.assertEqual(article["source_reported_high_quality_nucleus_n"], 69595)
        self.assertEqual(article["rna_annotation_reference"], "GSE189600")
        crossref = parse_crossref(
            json.dumps(
                {
                    "message": {
                        "relation": {
                            "has-preprint": [
                                {"id": "10.21203/rs.3.rs-6984670/v1"}
                            ]
                        },
                        "license": [],
                    }
                }
            ).encode()
        )
        self.assertTrue(crossref["preprint_relation_frozen"])

    def test_cohort_overlap_and_annotation_dependency_are_not_independent_confirmation(self) -> None:
        _, _, overlaps = validate_contract(self.contract)
        by_id = {row["overlap_id"]: row for row in overlaps}
        self.assertFalse(
            by_id["same_publication_gse281364_gse281160"]["independent_confirmation_allowed"]
        )
        self.assertFalse(
            by_id["gse189600_annotation_reference"]["independent_cell_label_evaluation_allowed"]
        )
        self.assertFalse(
            by_id["other_project_liver_cohorts"]["independent_donor_claim_allowed"]
        )

    def test_firewalls_cannot_be_relaxed(self) -> None:
        for key in (
            "automatic_activation", "sealed_outcomes_read", "biological_matrices_downloaded",
            "rna_observed", "rna_invention_allowed", "same_cell_pairing_allowed",
            "geo_mash_label_is_histology_adjudicated", "external_masld_champion_eligible",
        ):
            changed = copy.deepcopy(self.contract)
            changed[key] = True
            with self.assertRaises(ActivationReadinessError):
                validate_contract(changed)

    def test_frozen_project_authorities_are_hash_bound_without_matrix_read(self) -> None:
        authority = validate_project_artifacts(ROOT, self.contract)
        self.assertEqual(authority["project_participant_n"], 12)
        self.assertEqual(authority["project_audited_atac_cell_n"], 226224)
        self.assertEqual(authority["fixed_window_n"], 32000)
        self.assertFalse(authority["biological_matrix_read_by_this_audit"])


if __name__ == "__main__":
    unittest.main()
