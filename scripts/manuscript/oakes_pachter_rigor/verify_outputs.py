#!/usr/bin/env python3
"""Independent checks of the focused rigor outputs; no analysis refitting."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()
root = args.root.resolve()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=False)
base = root / "docs/technical/oakes_pachter_rigor_20260908"
checks = []


def check(name, passed, detail):
    checks.append(dict(check=name, passed=bool(passed), detail=str(detail)))


def read(path):
    return pd.read_csv(path, sep="\t")


tab = base / "tables-v4"
den = read(tab / "supplementary_analysis_denominators.tsv")
check("denominator_keys", not den.duplicated(["cohort", "assay", "analysis"]).any(), len(den))
bulk = den[den.assay.eq("bulk_RNA")]
check("sequential_count_losses", (bulk.source_records - bulk.metadata_without_counts == bulk.count_matched_records).all(), "Metadata scope exclusions remain explicitly labeled")
check("sequential_QC_losses", (bulk.count_matched_records - bulk.technical_failures == bulk.technical_qc_records).all(), "QC union, not sum of flags")
check("sequential_analysis_losses", (bulk.technical_qc_records - bulk.technical_pass_not_in_analysis == bulk.final_analysis_records).all(), "Analysis eligibility after QC")
members = read(tab / "bulk_analysis_membership.tsv")
observed = members.groupby(["dataset", "analysis"]).size()
recount = [observed.get((r.cohort, r.analysis), 0) for r in bulk.itertuples()]
check("bulk_denominators_from_members", np.array_equal(recount, bulk.final_analysis_records), len(bulk))
check("bulk_members_unique", not members.duplicated(["dataset", "sample_id", "analysis"]).any(), len(members))
families = read(tab / "statistical_family_checks.tsv")
check("BH_families", len(families) == 67 and families["pass"].all(), families.max_abs_BH_error.max())
ids = read(tab / "identifier_and_effect_checks.tsv")
failed = ids.loc[~ids["pass"], "check"].tolist()
check("only_expected_original_interval_failure", failed == ["stage_intervals_finite"], failed)
model = read(base / "models-v2/model_checks.tsv")
check("model_checks", len(model) == 23 and model.passed.all(), len(model))
eval_members = read(base / "models-v2/evaluation_membership.tsv")
check("model_members_unique", not eval_members.duplicated(["cohort", "analysis", "participant_id"]).any(), len(eval_members))
summary = read(base / "models-v2/model_evaluation_summary.tsv")
check("model_n_from_members", sorted(eval_members.groupby(["cohort", "analysis"]).size()) == sorted(summary.iloc[:4]["n"].astype(int)), "Four reproduced evaluations")
stage = root / "figures/candidates/pi-figure-redesign-2026-08-13-v8/analysis/stage_extensions/stage_extension_all_gene_results.tsv"
original = read(stage).sort_values(["contrast", "gene_id_versioned"]).reset_index(drop=True)
repair = read(base / "stage-interval-repair/stage_extension_all_gene_results.tsv.gz").sort_values(["contrast", "gene_id_versioned"]).reset_index(drop=True)
check("repaired_stage_keys", len(repair) == 257070 and original[["contrast", "gene_id_versioned"]].equals(repair[["contrast", "gene_id_versioned"]]), len(repair))
for name in ("logFC", "SE", "t", "P.Value", "FDR"):
    error = np.max(np.abs(original[name] - repair[name]))
    check("unchanged_" + name, error < 1e-12, error)
check("finite_repaired_intervals", np.isfinite(repair[["CI_low", "CI_high"]]).all().all(), len(repair))
check("intervals_contain_effect", ((repair.CI_low <= repair.logFC) & (repair.CI_high >= repair.logFC)).all(), "Pointwise intervals")
for folder in ("tables-v4", "models-v2", "stage-interval-repair"):
    bad = []
    for record in read(base / folder / "input_hashes.tsv").itertuples():
        h = hashlib.sha256()
        with open(record.source, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                h.update(block)
        if h.hexdigest() != record.sha256:
            bad.append(record.source)
    check("source_hashes_" + folder, not bad, bad)
assessment = root / "docs/technical/OAKES_PACHTER_RIGOR_INTEGRATION_ASSESSMENT_2026-09-08.md"
missing = []
for target in re.findall(r"\]\(([^)]+)\)", assessment.read_text()):
    if not target.startswith(("https://", "http://", "#")) and not (assessment.parent / target.split("#")[0]).exists():
        missing.append(target)
check("assessment_local_links", not missing, missing)
scope = subprocess.run(["python3", "scripts/manuscript/validate_resource_scope.py"], cwd=root, capture_output=True, text=True)
(out / "resource_scope.txt").write_text(scope.stdout + scope.stderr)
check("resource_scope", scope.returncode == 0, scope.stdout.strip())
pd.DataFrame(checks).to_csv(out / "verification_checks.tsv", sep="\t", index=False)
(out / "summary.json").write_text(json.dumps(dict(passed=sum(c["passed"] for c in checks), total=len(checks)), indent=2) + "\n")
for c in checks:
    print(c)
if not all(c["passed"] for c in checks):
    raise SystemExit(1)
