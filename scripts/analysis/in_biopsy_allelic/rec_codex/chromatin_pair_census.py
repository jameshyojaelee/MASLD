#!/usr/bin/env python3
"""Metadata-only availability of exact histology/sex donor contrasts."""
import csv
import hashlib
import json
import os
import platform
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
JOIN = ROOT / "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"
FOLDS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv"
FIELDS = ("steatosis", "ballooning", "lobular_inflammation", "fibrosis", "sex")


def census(rows):
    strata = Counter(tuple(row[k] for k in FIELDS) for row in rows)
    return {
        "participants": len(rows),
        "strata": len(strata),
        "participants_with_exact_partner": sum(n for n in strata.values() if n >= 2),
        "singleton_participants": sum(n for n in strata.values() if n == 1),
        "disjoint_pairs": sum(n // 2 for n in strata.values()),
        "complete_graph_edges": sum(n * (n - 1) // 2 for n in strata.values()),
        "strata_rows": [dict(zip(FIELDS, key), participants=n,
                             disjoint_pairs=n // 2) for key, n in sorted(strata.items())],
    }


def main():
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        raise RuntimeError("CPU SLURM allocation required")
    with JOIN.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with FOLDS.open() as stream:
        fold_rows = list(csv.DictReader(stream, delimiter="\t"))
    by_id = {row["participant_id"]: row for row in rows}
    folded = {row["participant_id"]: row["outer_fold"] for row in fold_rows}
    if len(rows) != 99 or len(by_id) != 99 or len(folded) != len(fold_rows):
        raise ValueError("Participant uniqueness/count mismatch")
    if by_id.keys() != folded.keys():
        raise ValueError("Participant axes differ")
    for row in rows:
        if row["outer_fold"] != folded[row["participant_id"]]:
            raise ValueError("Outer fold mismatch")
        if any(not row[field] or row[field] in ("NA", "nan") for field in FIELDS):
            raise ValueError("Missing recorded matching field")
    folds = sorted(set(folded.values()))
    if folds != list("01234"):
        raise ValueError("Expected five existing folds")
    comparisons = []
    for outer in folds:
        train = [row for row in rows if row["outer_fold"] != outer]
        test = [row for row in rows if row["outer_fold"] == outer]
        comparisons.append({"outer_fold": outer, "role": "outer_train", **census(train)})
        comparisons.append({"outer_fold": outer, "role": "outer_test", **census(test)})
        for inner in folds:
            if inner == outer:
                continue
            inner_train = [row for row in train if row["outer_fold"] != inner]
            inner_valid = [row for row in train if row["outer_fold"] == inner]
            train_ids = {row["participant_id"] for row in inner_train}
            valid_ids = {row["participant_id"] for row in inner_valid}
            test_ids = {row["participant_id"] for row in test}
            if train_ids & valid_ids or train_ids & test_ids or valid_ids & test_ids:
                raise ValueError("Donor split leakage")
            for role, subset in (("inner_train", inner_train), ("inner_valid", inner_valid)):
                comparisons.append({"outer_fold": outer, "inner_fold": inner,
                                    "role": role, **census(subset)})
    tests = [case for case in comparisons if case["role"] == "outer_test"]
    if sum(case["disjoint_pairs"] for case in tests) != 22:
        raise ValueError("Existing 22-pair evaluation census does not reproduce")
    summary = {
        "status": "metadata_only_no_model_fit_or_outcome_read",
        "matching_fields": FIELDS,
        "inner_split_rule": "Use each of four remaining existing participant folds as inner validation; pairs formed after donor split",
        "biological_unit": "participant; graph edges are dependent contrasts",
        "population": census(rows),
        "cases": comparisons,
        "checks": {"unique_99_participants": True, "fold_axes_agree": True,
                   "inner_outer_donors_disjoint": True, "existing_22_pairs_reproduced": True},
        "input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (JOIN, FOLDS)},
        "python": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "slurm_job_id": job,
        "no_stochastic_sampling": True,
    }
    out = Path(os.environ["CODEX_REC_OUTPUT"]) / f"chromatin_pair_census_{job}"
    out.mkdir(exist_ok=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output": str(out), "test_disjoint_pairs": 22,
                      "min_inner_train_disjoint_pairs": min(case["disjoint_pairs"] for case in comparisons if case["role"] == "inner_train")}))


if __name__ == "__main__":
    main()
