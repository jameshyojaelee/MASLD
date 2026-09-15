from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import tomllib
import unittest
import zipfile

import numpy as np

from masld_bench.contracts import TaskSpec
from scripts.build_gse268273_fibrosis_ood_fixture import (
    EXPECTED_PROMOTION_SHA256,
    GSE268273FixtureError,
    _iter_xlsx_rows,
    _numeric,
    _rank_rows,
    characteristics,
    parse_article_xml,
    parse_ena_tsv,
)


class GSE268273FibrosisOODFixtureTests(unittest.TestCase):
    def test_ena_runs_are_technical_partitions_of_unique_experiments(self) -> None:
        header = (
            "run_accession\tstudy_accession\tsample_accession\texperiment_accession\t"
            "scientific_name\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n"
        )
        rows = []
        for sample, experiment in (("SAMN1", "SRX1"), ("SAMN2", "SRX2")):
            for index in range(4):
                rows.append(
                    f"SRR{sample[-1]}{index}\tPRJNA1116068\t{sample}\t{experiment}\t"
                    f"Homo sapiens\tSINGLE\tftp.example/{sample}.{index}.fastq.gz\t100\t"
                    f"{'0' * 32}\n"
                )
        parsed, receipt = parse_ena_tsv(
            (header + "".join(rows)).encode(),
            {"SAMN1", "SAMN2"},
            {"SRX1", "SRX2"},
            expected_runs=8,
        )
        self.assertEqual(len(parsed), 8)
        self.assertEqual(receipt["biosamples"], 2)
        self.assertEqual(receipt["experiments"], 2)
        self.assertEqual(receipt["runs_per_experiment"], {"4": 2})
        self.assertEqual(receipt["run_role"], "technical_partition_not_biological_replicate")

    def test_ena_duplicate_run_fails_closed(self) -> None:
        payload = (
            "run_accession\tstudy_accession\tsample_accession\texperiment_accession\t"
            "scientific_name\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n"
            f"SRR1\tPRJNA1116068\tSAMN1\tSRX1\tHomo sapiens\tSINGLE\ta\t1\t{'0' * 32}\n"
            f"SRR1\tPRJNA1116068\tSAMN1\tSRX1\tHomo sapiens\tSINGLE\tb\t1\t{'1' * 32}\n"
        ).encode()
        with self.assertRaises(GSE268273FixtureError):
            parse_ena_tsv(payload, {"SAMN1"}, {"SRX1"}, expected_runs=2)

    def test_xlsx_reader_preserves_shared_strings_and_sparse_columns(self) -> None:
        main = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        document_rel = (
            "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
        )
        package_rel = (
            "http://schemas.openxmlformats.org/package/2006/relationships"
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.xlsx"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(
                    "xl/workbook.xml",
                    f'<workbook xmlns="{main}" xmlns:r="{document_rel}"><sheets>'
                    '<sheet name="RNAseq_expression_values_limma_" sheetId="1" '
                    'r:id="rId1"/></sheets></workbook>',
                )
                archive.writestr(
                    "xl/_rels/workbook.xml.rels",
                    f'<Relationships xmlns="{package_rel}"><Relationship Id="rId1" '
                    'Target="worksheets/sheet1.xml"/></Relationships>',
                )
                archive.writestr(
                    "xl/sharedStrings.xml",
                    f'<sst xmlns="{main}"><si><t>gene</t></si><si><t>x</t></si>'
                    '<si><t>y</t></si><si><t>P1</t></si><si><t>ENSG1</t></si></sst>',
                )
                archive.writestr(
                    "xl/worksheets/sheet1.xml",
                    f'<worksheet xmlns="{main}"><sheetData><row r="1">'
                    '<c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c>'
                    '<c r="C1" t="s"><v>2</v></c><c r="D1" t="s"><v>3</v></c>'
                    '</row><row r="2"><c r="A2" t="s"><v>4</v></c>'
                    '<c r="D2"><v>2.5</v></c></row></sheetData></worksheet>',
                )
            rows = list(_iter_xlsx_rows(path))
        self.assertEqual(rows[0], ["gene", "x", "y", "P1"])
        self.assertEqual(rows[1], ["ENSG1", None, None, 2.5])

    def test_article_contract_distinguishes_article_and_dataset_rights(self) -> None:
        required = " ".join(
            (
                "PMC11474426",
                "10.1016/j.jhepr.2024.101167",
                "cross-sectional, case-control study",
                "TSC of 109 histologically characterized cases",
                "Total RNA was extracted from freshly frozen biopsies",
                "human reference genome (GRCh38)",
                "gencode release 36",
                "RSEM (version 1.3.0)",
                "limma-voom",
                "adjusting for BMI, sex, and fibrosis severity",
            )
        )
        payload = (
            "<article><license href=\"https://creativecommons.org/licenses/by-nc-nd/4.0/\">"
            "Creative Commons Attribution-NonCommercial-NoDerivatives CC BY-NC-ND 4.0"
            "</license>"
            f"<body><p>{required}</p></body></article>"
        ).encode()
        receipt = parse_article_xml(payload)
        self.assertEqual(receipt["transcriptomic_participant_wording"], "109_histologically_characterized_cases")
        self.assertEqual(receipt["article_license"], "CC_BY_NC_ND_4_0")
        self.assertTrue(receipt["article_license_does_not_assign_dataset_terms"])

    def test_article_contract_rejects_plain_cc_by_license(self) -> None:
        required = " ".join(
            (
                "PMC11474426",
                "10.1016/j.jhepr.2024.101167",
                "cross-sectional, case-control study",
                "TSC of 109 histologically characterized cases",
                "Total RNA was extracted from freshly frozen biopsies",
                "human reference genome (GRCh38)",
                "gencode release 36",
                "RSEM (version 1.3.0)",
                "limma-voom",
                "adjusting for BMI, sex, and fibrosis severity",
            )
        )
        payload = (
            "<article><license>Creative Commons Attribution CC BY 4.0</license>"
            f"<body><p>{required}</p></body></article>"
        ).encode()
        with self.assertRaises(GSE268273FixtureError):
            parse_article_xml(payload)

    def test_metadata_parsing_preserves_missingness_and_locale_decimal(self) -> None:
        row = {
            "characteristics_json": (
                '["Sex: Female", "case/control: IMID-NAFLD", '
                '"hba1c: 5,9"]'
            )
        }
        observed = characteristics(row)
        self.assertEqual(observed["case/control"], "IMID-NAFLD")
        self.assertEqual(_numeric(observed["hba1c"]), ("5.9000000000000004", "observed"))
        self.assertEqual(_numeric(None), ("", "structurally_missing"))

    def test_tied_rank_uses_average_ranks(self) -> None:
        observed = _rank_rows(np.asarray([[5.0, 1.0, 1.0, 9.0]]))
        np.testing.assert_array_equal(observed, np.asarray([[3.0, 1.5, 1.5, 4.0]]))

    def test_task_is_blocked_native_fibrosis_transport_without_stage3_mapping(self) -> None:
        root = Path(__file__).parents[2]
        task_path = root / "config" / "evaluation" / "gse268273_fibrosis_ood_transfer_task.toml"
        gate_path = root / "config" / "evaluation" / "gse268273_fibrosis_ood_promotion_gate.json"
        task = TaskSpec.from_toml(task_path.read_text(encoding="utf-8"))
        self.assertEqual(task.status.value, "blocked")
        self.assertEqual(task.primary_metric, "participant_spearman_fibrosis_f0_f4")
        self.assertEqual(task.promotion_gate_config_sha256, EXPECTED_PROMOTION_SHA256)
        self.assertEqual(hashlib.sha256(gate_path.read_bytes()).hexdigest(), EXPECTED_PROMOTION_SHA256)
        self.assertNotIn("NOR", task.endpoint)
        self.assertIn("All 109 participants have biopsy-proven MASLD", "\n".join(task.admission_gates))
        with task_path.open("rb") as handle:
            raw = tomllib.load(handle)
        self.assertFalse(raw["evaluator_parameters"]["deposited_voom_scored_input_allowed"])
        self.assertFalse(raw["evaluator_parameters"]["imid_status_allowed_as_primary_model_input"])


if __name__ == "__main__":
    unittest.main()
