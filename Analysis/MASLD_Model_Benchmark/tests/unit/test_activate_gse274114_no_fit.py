from __future__ import annotations

import unittest

from scripts.activate_gse274114_no_fit import (
    GSE274114ActivationError,
    aggregate_quant_records,
    audit_exact_identifier_overlap,
    assign_folds,
    build_crosswalk,
    stable_gene_id,
)


class GSE274114NoFitActivationTests(unittest.TestCase):
    def test_stable_gene_parser_retains_par_y_state(self) -> None:
        self.assertEqual(stable_gene_id("ENSG00000185203.12_PAR_Y"), ("ENSG00000185203", True))
        self.assertEqual(stable_gene_id("ENSG00000185203.13"), ("ENSG00000185203", False))
        with self.assertRaises(GSE274114ActivationError):
            stable_gene_id("WASIR1")

    def test_crosswalk_marks_par_pair_for_raw_scale_aggregation(self) -> None:
        genes = {
            "ENSG00000185203": {
                "gencode_v49_gene_id": "ENSG00000185203.13",
                "gencode_v49_gene_name": "WASIR1",
                "gencode_v49_gene_type": "lncRNA",
            }
        }
        observed = build_crosswalk(
            ["ENSG00000185203.12", "ENSG00000185203.12_PAR_Y"], genes
        )
        self.assertEqual({row["target_source_row_count"] for row in observed}, {"2"})
        self.assertTrue(all(row["aggregation_rule"].startswith("sum_NumReads") for row in observed))

    def test_par_pair_is_materialized_on_one_v49_row_without_rounding(self) -> None:
        genes = {
            "ENSG00000185203": {
                "gencode_v49_gene_id": "ENSG00000185203.13",
                "gencode_v49_gene_name": "WASIR1",
                "gencode_v49_gene_type": "lncRNA",
            }
        }
        crosswalk = build_crosswalk(
            ["ENSG00000185203.12", "ENSG00000185203.12_PAR_Y"], genes
        )
        observed, audit = aggregate_quant_records(
            [
                {"Name": "ENSG00000185203.12", "Length": "100", "EffectiveLength": "80", "TPM": "1", "NumReads": "2.25"},
                {"Name": "ENSG00000185203.12_PAR_Y", "Length": "120", "EffectiveLength": "100", "TPM": "3", "NumReads": "4.75"},
            ],
            crosswalk,
        )
        self.assertEqual(len(observed), 1)
        self.assertEqual(observed[0]["Name"], "ENSG00000185203.13")
        self.assertEqual(observed[0]["NumReads"], 7.0)
        self.assertEqual(observed[0]["TPM"], 4.0)
        self.assertEqual(observed[0]["Length"], 115.0)
        self.assertEqual(observed[0]["EffectiveLength"], 95.0)
        self.assertEqual(audit["output_numreads"], audit["mapped_numreads"])

    def test_zero_tpm_par_pair_uses_explicit_no_abundance_fallback(self) -> None:
        genes = {
            "ENSG00000185203": {
                "gencode_v49_gene_id": "ENSG00000185203.13",
                "gencode_v49_gene_name": "WASIR1",
                "gencode_v49_gene_type": "lncRNA",
            }
        }
        crosswalk = build_crosswalk(
            ["ENSG00000185203.12", "ENSG00000185203.12_PAR_Y"], genes
        )
        observed, audit = aggregate_quant_records(
            [
                {"Name": "ENSG00000185203.12", "Length": "100", "EffectiveLength": "80", "TPM": "0", "NumReads": "0"},
                {"Name": "ENSG00000185203.12_PAR_Y", "Length": "120", "EffectiveLength": "100", "TPM": "0", "NumReads": "0"},
            ],
            crosswalk,
        )
        self.assertEqual(observed[0]["Length"], 110.0)
        self.assertEqual(observed[0]["EffectiveLength"], 90.0)
        self.assertEqual(audit["zero_tpm_multirow_targets"], 1)

    def test_folds_preserve_all_four_source_groups(self) -> None:
        labels = []
        for group, count in {"CTRL": 9, "ENEG": 11, "NASH": 10, "ENEG_NASH": 9}.items():
            labels.extend({"row_id": f"{group}_{index}", "source_group": group} for index in range(count))
        observed = assign_folds(labels)
        self.assertEqual(len(observed), 39)
        self.assertEqual(set(observed.values()), {0, 1, 2, 3, 4})

    def test_overlap_audit_fails_closed_on_accession(self) -> None:
        import tempfile
        from pathlib import Path

        topology = [{
            "participant_id": "P1", "sample_accession": "GSM1",
            "biosample_accession": "SAMN1", "sra_experiment": "SRX1",
        }]
        runs = [{"sra_run": "SRR1"}]
        with tempfile.TemporaryDirectory() as temporary:
            evidence = Path(temporary) / "other.toml"
            evidence.write_text('accession = ["GSE274114"]\n', encoding="utf-8")
            with self.assertRaises(GSE274114ActivationError):
                audit_exact_identifier_overlap(topology, runs, [evidence])


if __name__ == "__main__":
    unittest.main()
