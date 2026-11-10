"""Check fixed-family metadata and unchanged numerical inference on real results."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    with path.open("x") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--real-reference", type=Path, required=True)
    parser.add_argument("--real-output", type=Path, required=True)
    parser.add_argument("--da", type=Path, required=True)
    parser.add_argument("--network", type=Path, action="append", required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run in a compute allocation")
    args.output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[1] / "scripts/summarize_candidate.py"
    spec = importlib.util.spec_from_file_location("pooled_family_fixture", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    fixture_da = args.output / "da"
    fixture_da.mkdir()
    (fixture_da / "run_summary.json").write_text("{}")
    for filename in ["cohort_summary.tsv", "baseline_reproduction.tsv"]:
        (fixture_da / filename).write_text("fixture\n")
    keys = ["cohort", "lineage", "endpoint", "method", "baseline"]
    planned = list(module.planned_comparisons())
    fixtures = []
    for label in ["complete", "missing_test", "untestable_test"]:
        root = args.output / label
        root.mkdir()
        (root / "completion.json").write_text('{"status":"completed"}')
        (root / "run_settings.json").write_text('{"pilot":false,"analysis_scope":"genome_wide"}')
        rows = [dict(**dict(zip(keys, key)), pvalue=str((i+1)/1000), status="evaluated", ci_low="0.2",
                     family="source_provisional", family_complete="False") for i, key in enumerate(planned)]
        if label == "missing_test":
            rows.pop(0)
        elif label == "untestable_test":
            rows[0]["status"] = "untestable"
        write(root / "comparisons.tsv", rows)
        destination = args.output / (label + "_pooled")
        subprocess.run([sys.executable, str(source), "--da", str(fixture_da), "--network", str(root), "--output", str(destination)], check=True)
        pooled = read(destination / "network_comparisons.tsv")
        assert len(pooled) == 48
        assert all(row["family_complete"] == str(label == "complete") for row in pooled)
        assert all(row["family"] == "prespecified_48_cohort_lineage_endpoint_method_baseline_comparisons" for row in pooled)
        assert all(row["family_size"] == "48" for row in pooled)
        assert all(row["source_family"] == "source_provisional" for row in pooled if row["source"])
        assert all(row["source_family_complete"] == "False" for row in pooled if row["source"])
        if label == "missing_test":
            assert pooled[0]["status"] == "not_run" and pooled[0]["pvalue"] == "1.0"
        fixtures.append(label)
    command = [sys.executable, str(source), "--da", str(args.da), "--output", str(args.real_output)]
    for network in args.network:
        command += ["--network", str(network)]
    subprocess.run(command, check=True)
    before = read(args.real_reference / "network_comparisons.tsv")
    after = read(args.real_output / "network_comparisons.tsv")
    assert len(before) == len(after) == 48
    permitted = {"family", "family_complete", "source_family", "source_family_complete"}
    unchanged = sorted(set(before[0]) - permitted)
    for old, new in zip(before, after):
        assert all(old[field] == new[field] for field in unchanged)
        assert new["source_family"] == old["family"]
        assert new["source_family_complete"] == old["family_complete"]
        assert new["family_complete"] == "True"
    old_summary = json.loads((args.real_reference / "summary.json").read_text())
    new_summary = json.loads((args.real_output / "summary.json").read_text())
    assert new_summary.pop("network_family_complete") is True
    assert old_summary == new_summary
    receipt = dict(status="passed", fixtures=fixtures, real_comparisons=48, family_complete=True,
        unchanged_columns=unchanged, numerical_values_keys_and_supported_calls_identical=True,
        changed_metadata=sorted(permitted), supported_calls=sum(r["predictive_improvement"] == "True" for r in after),
        original=str(args.real_reference.resolve()), corrected=str(args.real_output.resolve()))
    (args.output / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
