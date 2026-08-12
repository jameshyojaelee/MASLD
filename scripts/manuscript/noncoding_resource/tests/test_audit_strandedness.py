import importlib.util
import csv
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "audit_strandedness.py"
SPEC = importlib.util.spec_from_file_location("audit_strandedness", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_parse_star_uses_gene_rows_only(tmp_path):
    path = tmp_path / "sample.ReadsPerGene.out.tab"
    path.write_text("N_unmapped\t1\t2\t3\nGENE1\t10\t4\t20\nGENE2\t5\t1\t10\n")
    assert MODULE.parse_star(path) == (15, 5, 30)


def test_parse_featurecounts_assigned(tmp_path):
    path = tmp_path / "counts.summary"
    path.write_text(
        "Status\t/a/S1.Aligned.bam\t/a/S2.Aligned.bam\nAssigned\t10\t20\nUnassigned_NoFeatures\t1\t2\n"
    )
    assert MODULE.parse_featurecounts_assigned(path) == {"S1": 10, "S2": 20}


def test_amended_cohort_rule_keeps_all_samples():
    row = {
        "n_audited": 100,
        "n_below_fourfold": 6,
        "median_reverse_same_ratio": 10.0,
        "minimum_ratio": 1.2,
    }
    assert MODULE.amended_cohort_pass(row)
    row["minimum_ratio"] = 1.0
    assert not MODULE.amended_cohort_pass(row)


def test_stage_concentration_gate(tmp_path):
    metadata = tmp_path / "metadata.csv"
    with metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("sample_id", "dataset", "condition")
        )
        writer.writeheader()
        rows = []
        sample_rows = []
        for group in ("Control", "Fibrosis_F0F1", "Fibrosis_F2", "Fibrosis_F3F4"):
            for index in range(20):
                sample_id = f"{group}_{index}"
                rows.append(
                    {"sample_id": sample_id, "dataset": "GSE213621", "condition": group}
                )
                sample_rows.append(
                    {
                        "sample_id": sample_id,
                        "dataset": "GSE213621",
                        "pass_fourfold": index != 0,
                    }
                )
        writer.writerows(rows)
    result = MODULE.stage_concentration_audit(sample_rows, metadata)
    assert len(result) == 4
    assert all(row["stage_concentration_gate_passed"] for row in result)
