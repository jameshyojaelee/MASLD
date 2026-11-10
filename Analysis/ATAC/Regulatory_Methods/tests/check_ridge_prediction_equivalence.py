#!/usr/bin/env python3
"""Compare completed ridge predictions with a bounded original-file prefix.

Original files may be growing. Their initial byte limits are captured before
reading and every captured byte is hashed. Only an incomplete final gene may be
omitted, and its written rows must equal a prefix of the complete candidate gene.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path


class PrefixRows:
    """Line-bounded snapshot without copying or writing to the original file."""
    def __init__(self, path):
        self.path = Path(path)
        self.limit = self.path.stat().st_size
        self.hash = hashlib.sha256()
        self.bytes_read = 0
        self.partial = b""
        self.header = None
        self.rows = 0

    def __iter__(self):
        with self.path.open("rb") as handle:
            while self.bytes_read < self.limit:
                raw = handle.readline(self.limit-self.bytes_read)
                if not raw:
                    raise ValueError(f"Original file shrank during bounded read: {self.path}")
                self.hash.update(raw)
                self.bytes_read += len(raw)
                if not raw.endswith(b"\n"):
                    self.partial = raw
                    break
                fields = next(csv.reader([raw.decode("utf-8")], delimiter="\t"))
                if self.header is None:
                    self.header = fields
                    if not fields or len(fields) != len(set(fields)) or "gene_id" not in fields:
                        raise ValueError("Malformed prediction or tuning header")
                    continue
                if len(fields) != len(self.header):
                    raise ValueError("Malformed complete row inside bounded source prefix")
                self.rows += 1
                yield dict(zip(self.header, fields)), raw
        if self.header is None:
            raise ValueError("Bounded source lacks a complete header")

    def receipt(self):
        if self.bytes_read != self.limit:
            raise ValueError("Source prefix was not fully read and hashed")
        return dict(path=str(self.path.resolve()), captured_bytes=self.limit,
                    captured_sha256=self.hash.hexdigest(), complete_rows=self.rows,
                    truncated_final_row_bytes=len(self.partial),
                    bytes_at_finish=self.path.stat().st_size)


def gene_groups(rows):
    previous = None
    for gene, group in itertools.groupby(rows, key=lambda record: record[0]["gene_id"]):
        if previous is not None and gene <= previous:
            raise ValueError("Prediction or tuning gene groups are not strictly ordered")
        previous = gene
        yield gene, list(group)


def with_last(iterator):
    iterator = iter(iterator)
    current = next(iterator, None)
    while current is not None:
        following = next(iterator, None)
        yield current, following is None
        current = following


def verify_rows(records, kind):
    keys = set()
    for row, raw in records:
        if row["grid"] not in {"original_replay", "wide"}:
            raise ValueError("Unexpected ridge grid in compared rows")
        if kind == "predictions":
            if row["endpoint"] != "rna":
                raise ValueError("Unexpected non-RNA endpoint")
            key = row["grid"], row["method"], row["donor_id"]
            value = float(row["standardized_squared_error"])
            if not math.isfinite(value) or value < 0:
                raise ValueError("Invalid held-donor prediction loss")
        else:
            key = row["grid"], row["method"]
        if key in keys:
            raise ValueError("Duplicate grid/method/donor key within gene")
        keys.add(key)
    for grid in ("original_replay", "wide"):
        other = {key[1:] for key in keys if key[0] == grid}
        paired = {key[1:] for key in keys if key[0] != grid}
        if other != paired:
            raise ValueError("Complete candidate gene lacks matched original/wide keys")


def compare_file(candidate_path, original_path, pilot=False, original_complete=False, kind="predictions"):
    candidate_reader, original_reader = PrefixRows(candidate_path), PrefixRows(original_path)
    candidate_groups = iter(gene_groups(candidate_reader))
    candidate = next(candidate_groups, None)
    if candidate is None:
        raise ValueError("Candidate has no complete evaluable gene")
    matched_genes, matched_rows, skipped_original_rows = 0, 0, 0
    incomplete_gene = None
    expected_partial = None
    candidate_consumed = False
    for (gene, original), last in with_last(gene_groups(original_reader)):
        if candidate is None:
            if not pilot:
                raise ValueError("Original has a gene absent from complete candidate")
            skipped_original_rows += len(original)
            continue
        candidate_gene, proposed = candidate
        verify_rows(proposed, kind)
        if candidate_gene != gene:
            raise ValueError(f"Missing or reordered middle gene: original={gene}, candidate={candidate_gene}")
        if candidate_reader.header != original_reader.header:
            raise ValueError("Original and candidate columns differ")
        if len(original) > len(proposed):
            raise ValueError(f"Original gene has extra rows: {gene}")
        for (old, old_raw), (new, new_raw) in zip(original, proposed):
            if old != new:
                changes = {key: [old.get(key), new.get(key)] for key in set(old) | set(new) if old.get(key) != new.get(key)}
                raise ValueError(f"Changed {kind} gene={gene}: {changes}")
        if len(original) < len(proposed):
            if not last or original_complete:
                raise ValueError(f"Incomplete middle or completed original gene: {gene}")
            incomplete_gene = gene
            expected_partial = proposed[len(original)][1]
            # A partial gene is excluded in its entirety from equivalence counts.
            if pilot:
                raise ValueError("Original prefix does not yet cover every candidate pilot gene")
        else:
            verify_rows(original, kind)
            matched_genes += 1
            matched_rows += len(original)
            candidate = next(candidate_groups, None)
            if candidate is None:
                candidate_consumed = True
            elif last:
                expected_partial = candidate[1][0][1]
    if not matched_genes:
        raise ValueError("No complete gene could be compared")
    if original_reader.partial and not (pilot and candidate_consumed):
        if original_complete or expected_partial is None or not expected_partial.startswith(original_reader.partial):
            raise ValueError("Truncated EOF row is not an exact candidate-row prefix")
        if incomplete_gene is None:
            incomplete_gene = candidate[0]
    if pilot and not candidate_consumed:
        raise ValueError("Not every completed pilot gene is present in original prefix")
    if original_complete and not pilot and candidate is not None:
        raise ValueError("Completed original and candidate gene coverage differs")
    # Complete candidate data must remain structurally valid beyond the prefix.
    if candidate is not None:
        verify_rows(candidate[1], kind)
    for gene, proposed in candidate_groups:
        verify_rows(proposed, kind)
    if candidate_reader.partial:
        raise ValueError("Completed candidate has a truncated final row")
    return dict(kind=kind, exact_string_equality=True, matched_complete_genes=matched_genes,
                matched_rows=matched_rows, excluded_incomplete_final_original_gene=incomplete_gene,
                untested_original_rows_after_pilot=skipped_original_rows,
                candidate=candidate_reader.receipt(), original=original_reader.receipt())


def checked_settings(candidate, original, pilot):
    if any(candidate.glob("INVALIDATED*.json")) or any(original.glob("INVALIDATED*.json")):
        raise ValueError("Invalidated candidate or original fold")
    complete = json.loads((candidate / "completion.json").read_text())
    if complete.get("status") != "completed":
        raise ValueError("Candidate must be completed")
    if bool(complete.get("technical_check")) != pilot or (complete.get("max_genes", 0) > 0) != pilot:
        raise ValueError("Pilot flag does not match completed candidate scope")
    a = json.loads((candidate / "run_settings.json").read_text())
    b = json.loads((original / "run_settings.json").read_text())
    fields = ["prepared_root", "baseline_root", "group", "lineage", "outer_fold", "original_grid", "wide_grid", "endpoint", "seed"]
    for key in fields:
        if key not in a or key not in b or a[key] != b[key] or complete.get(key) != a[key]:
            raise ValueError(f"Candidate/original scientific fold setting differs: {key}")
    candidate_inner = (candidate / "inner_fold_assignments.tsv").read_bytes()
    original_inner = (original / "inner_fold_assignments.tsv").read_bytes()
    if candidate_inner != original_inner:
        raise ValueError("Inner donor/metacell fold assignments differ")
    done = original / "completion.json"
    original_complete = done.exists() and json.loads(done.read_text()).get("status") == "completed"
    identity = {field: a[field] for field in fields}
    identity["inner_fold_assignments_sha256"] = hashlib.sha256(candidate_inner).hexdigest()
    return complete, original_complete, identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--original", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    complete, original_complete, settings = checked_settings(args.candidate, args.original, args.pilot)
    comparisons = [compare_file(args.candidate / filename, args.original / filename,
                                args.pilot, original_complete, kind)
                   for filename, kind in [("predictions.tsv", "predictions"), ("tuning.tsv", "tuning")]]
    n_prediction_rows = comparisons[0]["candidate"]["complete_rows"]
    if complete.get("n_prediction_rows") != n_prediction_rows or complete.get("n_original_predictions_verified", -1) * 2 != n_prediction_rows:
        raise ValueError("Completed candidate row/anchor counts differ from parsed data")
    result = dict(status="passed", pilot=args.pilot, original_complete_at_start=original_complete,
                  matched_complete_prediction_genes=comparisons[0]["matched_complete_genes"],
                  matched_prediction_rows=comparisons[0]["matched_rows"],
                  exact_both_grid_predictions=True, exact_tuning=True,
                  settings=settings, comparisons=comparisons,
                  interpretation="exact_implementation_equivalence_on_captured_complete_original_genes_not_new_donor_validation")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps({key: result[key] for key in ["status", "pilot", "matched_complete_prediction_genes", "matched_prediction_rows"]}))


if __name__ == "__main__":
    main()
