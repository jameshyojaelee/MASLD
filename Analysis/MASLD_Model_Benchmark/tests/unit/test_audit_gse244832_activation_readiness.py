from __future__ import annotations

import copy
import csv
import gzip
import io
from pathlib import Path
import tomllib
import unittest

from scripts.audit_gse244832_activation_readiness import (
    ActivationReadinessError,
    parse_geo_soft,
    parse_sra_runinfo,
    validate_contract,
    validate_project_artifacts,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/gse244832_activation_readiness.toml"


def geo_fixture() -> bytes:
    lines = [
        "^SERIES = GSE244832",
        "!Series_supplementary_file = ftp://example/GSE244832_RAW.tar",
        "!Series_supplementary_file = ftp://example/GSE244832_hLIVER_processed_files.tar.gz",
    ]
    conditions = ["NORMAL"] * 5 + ["MASL"] * 4 + ["MASH"] * 9
    gsm = 7830541
    for modality in ("ATAC", "RNA"):
        for index, condition in enumerate(conditions, start=1):
            title = f"MM_{400 + index}, snATACseq" if modality == "ATAC" else f"JB_{200 + index}, snRNAseq"
            strategy = "ATAC-seq" if modality == "ATAC" else "RNA-Seq"
            lines.extend(
                [
                    f"^SAMPLE = GSM{gsm}",
                    f"!Sample_title = {title}",
                    f"!Sample_description = {condition}",
                    f"!Sample_library_strategy = {strategy}",
                    f"!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN{gsm}",
                    f"!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX{gsm}",
                    (
                        f"!Sample_supplementary_file_1 = ftp://example/GSM{gsm}_{title.split(',')[0]}.sorted.bed.gz"
                        if modality == "ATAC"
                        else "!Sample_supplementary_file_1 = NONE"
                    ),
                ]
            )
            gsm += 1
    return gzip.compress(("\n".join(lines) + "\n").encode())


def sra_fixture(geo_records: list[dict[str, object]]) -> str:
    fields = [
        "Run",
        "Experiment",
        "BioSample",
        "LibraryStrategy",
        "LibraryLayout",
        "Consent",
        "size_MB",
    ]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    run_number = 1
    rna_index = 0
    for record in geo_records:
        if record["modality"] == "ATAC":
            run_count = 1
            strategy = "ATAC-seq"
        else:
            run_count = 8 if rna_index < 9 else 5
            strategy = "RNA-Seq"
            rna_index += 1
        for _ in range(run_count):
            writer.writerow(
                {
                    "Run": f"SRR{run_number:08d}",
                    "Experiment": record["experiment"],
                    "BioSample": record["biosample"],
                    "LibraryStrategy": strategy,
                    "LibraryLayout": "PAIRED",
                    "Consent": "public",
                    "size_MB": "10",
                }
            )
            run_number += 1
    return stream.getvalue()


class GSE244832ActivationReadinessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_counts_and_different_aliquot_topology_are_binding(self) -> None:
        validate_contract(self.contract)
        self.assertEqual(self.contract["declared_participant_n"], 18)
        self.assertEqual(self.contract["rna_assay_record_n"], 18)
        self.assertEqual(self.contract["atac_assay_record_n"], 18)
        self.assertEqual(self.contract["pairing_topology"], "same_sample_different_aliquot")
        self.assertFalse(self.contract["same_cell_pairing_allowed"])
        self.assertFalse(self.contract["same_nucleus_pairing_allowed"])
        self.assertFalse(self.contract["infer_cross_assay_join_by_record_order"])

    def test_only_disease_group_is_complete_model_metadata(self) -> None:
        phenotypes, _ = validate_contract(self.contract)
        by_id = {row["field_id"]: row for row in phenotypes}
        self.assertEqual(by_id["disease_group"]["missingness_state"], "observed")
        self.assertEqual(by_id["fibrosis_stage"]["missingness_state"], "join_unresolved")
        for field_id in ("nas_total", "histology_components", "sex", "age", "bmi_or_body_composition"):
            self.assertEqual(by_id[field_id]["missingness_state"], "structurally_missing")

    def test_assay_native_development_is_separate_from_cross_assay_and_external(self) -> None:
        _, tasks = validate_contract(self.contract)
        by_id = {row["task_id"]: row for row in tasks}
        self.assertEqual(
            by_id["atac_fixed_window_donor_lineage_transport"]["eligibility"],
            "active_development",
        )
        self.assertEqual(
            by_id["cross_assay_donor_pseudobulk_late_fusion"]["eligibility"],
            "conditional",
        )
        self.assertEqual(
            by_id["same_cell_or_same_nucleus_rna_atac_pairing"]["eligibility"],
            "prohibited",
        )
        self.assertEqual(
            by_id["external_or_sealed_masld_champion"]["eligibility"],
            "prohibited",
        )

    def test_firewalls_cannot_be_relaxed(self) -> None:
        for key in (
            "automatic_activation",
            "sealed_outcomes_read",
            "biological_matrices_downloaded",
            "same_cell_pairing_allowed",
            "same_nucleus_pairing_allowed",
            "infer_cross_assay_join_by_record_order",
            "external_masld_champion_eligible",
        ):
            changed = copy.deepcopy(self.contract)
            changed[key] = True
            with self.assertRaises(ActivationReadinessError):
                validate_contract(changed)

    def test_geo_parser_preserves_unjoined_assay_namespaces(self) -> None:
        records, files = parse_geo_soft(geo_fixture())
        self.assertEqual(len(records), 36)
        self.assertEqual(sum(row["modality"] == "ATAC" for row in records), 18)
        self.assertEqual(sum(row["modality"] == "RNA" for row in records), 18)
        self.assertTrue(all(row["cross_assay_participant_id"] is None for row in records))
        self.assertTrue(all("not_inferred_from_record_order" in row["cross_assay_join_state"] for row in records))
        self.assertEqual(len(files), 2)

    def test_sra_parser_resolves_raw_runs_without_downloading_reads(self) -> None:
        records, _ = parse_geo_soft(geo_fixture())
        rows, summary = parse_sra_runinfo(sra_fixture(records), records)
        self.assertEqual(len(rows), 135)
        self.assertEqual(summary["experiment_n"], 36)
        self.assertEqual(
            summary["run_n_by_library_strategy"], {"ATAC-seq": 18, "RNA-Seq": 117}
        )
        self.assertFalse(summary["raw_reads_downloaded"])
        bad = sra_fixture(records).replace("public", "controlled", 1)
        with self.assertRaises(ActivationReadinessError):
            parse_sra_runinfo(bad, records)

    def test_frozen_atac_authority_is_hash_bound_without_matrix_read(self) -> None:
        authority = validate_project_artifacts(ROOT, self.contract)
        self.assertEqual(authority["atac_participant_n"], 18)
        self.assertEqual(authority["atac_cell_n"], 88814)
        self.assertEqual(authority["pairing_topology"], "same_sample_different_aliquot")
        self.assertFalse(authority["biological_matrix_read_by_this_audit"])


if __name__ == "__main__":
    unittest.main()
