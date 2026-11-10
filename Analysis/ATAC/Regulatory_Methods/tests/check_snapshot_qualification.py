"""Check finished-v3 qualification independently of reference-shard completion."""
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
    original, reference = args.fixtures / "monolithic", args.fixtures / "shard1"
    group = "GSE296875_hepatocyte"
    marker = reference / group / "shard_completion.json"
    saved = marker.read_bytes()
    try:
        marker.write_text('{"status":"still_running"}')
        qualification = merger.qualify_snapshot(original, 0, reference, args.output / "valid.json")
    finally:
        marker.write_bytes(saved)
    rejected = []
    def reject(label):
        try:
            merger.qualify_snapshot(original, 0, reference, args.output / (label + ".json"))
        except (ValueError, FileNotFoundError):
            rejected.append(label)
        else:
            raise AssertionError(f"Invalid snapshot accepted: {label}")
    def mutate(label, relative, operation):
        path = original / relative
        saved = path.read_bytes()
        try:
            path.write_bytes(operation(saved))
            reject(label)
        finally:
            path.write_bytes(saved)
    mutate("failed_cicero", group + "/cicero_fold0_outer/run.log", lambda saved: b"Models with errors: 1\nwindow_models_complete\n")
    mutate("incomplete_coverage", group + "/coverage.tsv", lambda saved: b"".join(saved.splitlines(keepends=True)[:1] + saved.splitlines(keepends=True)[2:]))
    mutate("incomplete_predictions", group + "/predictions.tsv", lambda saved: b"".join(saved.splitlines(keepends=True)[:1] + saved.splitlines(keepends=True)[2:]))
    def wrong_donor(saved):
        reader = csv.DictReader(io.StringIO(saved.decode()), delimiter="\t")
        rows, columns = list(reader), reader.fieldnames
        rows[0]["donor_id"] = "not_held_out"
        target = io.StringIO()
        writer = csv.DictWriter(target, columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        return target.getvalue().encode()
    mutate("wrong_held_out_donor", group + "/predictions.tsv", wrong_donor)
    mutate("different_membership", group + "_metacell_membership.tsv", lambda saved: saved + b"different\n")
    mutate("unreviewed_source", "source_snapshot/network_compare.py", lambda saved: saved + b"\n")
    (args.output / "verification.json").write_text(json.dumps(dict(status="passed",
        prepared_reference_sufficient=True, n_qualified_genes=qualification["n_genes"], rejected=rejected,
        merger_sha256=merger.sha(Path(merger.__file__))), indent=2) + "\n")


if __name__ == "__main__":
    main()
