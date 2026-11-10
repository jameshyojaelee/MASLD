"""Recheck final merger on finished small fixtures, without repeating model fits."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import merge_network_folds as merger
    args.output.mkdir(parents=True, exist_ok=False)
    shards = [args.fixtures / f"shard{fold}" for fold in range(5)]
    monolithic = args.fixtures / "monolithic"
    group = "GSE296875_hepatocyte"
    scientific = ["predictions.tsv", "links.tsv", "coverage.tsv", "donor_errors.tsv", "comparisons.tsv", "fold_assignments.tsv", "method_status.tsv", "method_execution_summary.tsv", "failures.tsv"]
    for label, inputs, snapshots in [("all_shards", shards, []), ("legacy_snapshot", shards[1:], [(monolithic, 0)])]:
        destination = args.output / label
        merger.merge(inputs, snapshots, destination)
        for filename in scientific:
            if (monolithic / group / filename).read_bytes() != (destination / group / filename).read_bytes():
                raise AssertionError(f"Nonidentical {label}: {filename}")
    rejected = []
    def reject(label, inputs=None, snapshots=None):
        try:
            merger.merge(shards if inputs is None else inputs, [] if snapshots is None else snapshots, args.output / label)
        except (ValueError, FileNotFoundError):
            rejected.append(label)
        else:
            raise AssertionError(f"Invalid input accepted: {label}")
    reject("missing_fold", shards[1:])
    reject("duplicate_fold", shards + [shards[0]])
    reject("unflushed_legacy_fold", shards[:4], [(monolithic, 4)])
    def mutate_table(label, filename, mutation):
        path = shards[0] / group / filename
        saved = path.read_bytes()
        try:
            reader = csv.DictReader(io.StringIO(saved.decode()), delimiter="\t")
            columns, rows = reader.fieldnames, list(reader)
            mutation(rows)
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, columns, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            reject(label)
        finally:
            path.write_bytes(saved)
    mutate_table("missing_prediction", "predictions.tsv", lambda rows: rows.pop(0))
    mutate_table("duplicate_prediction", "predictions.tsv", lambda rows: rows.insert(0, dict(rows[0])))
    mutate_table("wrong_held_out_donor", "predictions.tsv", lambda rows: rows[0].update(donor_id="not_held_out"))
    mutate_table("missing_gene", "coverage.tsv", lambda rows: rows.pop(0))
    def swap(rows):
        rows[0], rows[1] = rows[1], rows[0]
    mutate_table("reordered_gene", "coverage.tsv", swap)
    def duplicate_peak(rows):
        rows[1]["peak_id"] = rows[0]["peak_id"]
    mutate_table("duplicate_peak_different_rank", "links.tsv", duplicate_peak)
    path = shards[0] / "INVALIDATED_test.json"
    try:
        path.write_text("{}")
        reject("invalidated_source")
    finally:
        path.unlink()
    path = shards[0] / group / "shard_completion.json"
    saved = path.read_bytes()
    try:
        marker = json.loads(saved)
        marker["status"] = "failed"
        path.write_text(json.dumps(marker))
        reject("failed_shard")
    finally:
        path.write_bytes(saved)
    (args.output / "verification.json").write_text(json.dumps(dict(status="passed", exact_tables=scientific,
        rejected=rejected, merger_sha256=merger.sha(Path(merger.__file__))), indent=2) + "\n")


if __name__ == "__main__":
    main()
