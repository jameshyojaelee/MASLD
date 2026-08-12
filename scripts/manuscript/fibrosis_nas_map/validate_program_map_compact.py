#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
import re
import subprocess
import sys
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


if len(sys.argv) != 3:
    fail("Usage: validate_program_map_compact.py WORKSTREAM_ROOT OUTPUT_ROOT")

workstream = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2]).resolve(strict=True)
input_table = workstream / "discovery" / "program_map.tsv"
source_table = output / "source_table.tsv"
map_pdf = output / "fig3_program_effect_map_impact.pdf"
method_pdf = output / "fig3_program_score_method.pdf"

expected_files = {
    "artifact_manifest.tsv",
    "fig3_program_effect_map_impact.pdf",
    "fig3_program_score_method.pdf",
    "figure_caption.txt",
    "figure_summary.tsv",
    "render_spec.tsv",
    "sessionInfo.txt",
    "source_table.tsv",
}
observed_files = {path.name for path in output.iterdir() if path.is_file()}
if observed_files != expected_files:
    fail(f"Unexpected compact release files: {sorted(observed_files ^ expected_files)}")

with input_table.open(newline="") as handle:
    input_rows = list(csv.DictReader(handle, delimiter="\t"))
with source_table.open(newline="") as handle:
    source_rows = list(csv.DictReader(handle, delimiter="\t"))

if len(input_rows) != 117 or len(source_rows) != 117:
    fail("Input and source tables must each retain all 117 frozen programs")
if len({row["feature_id"] for row in source_rows}) != 117:
    fail("source_table.tsv contains duplicate feature IDs")

input_by_id = {row["feature_id"]: row for row in input_rows}
source_by_id = {row["feature_id"]: row for row in source_rows}
if input_by_id.keys() != source_by_id.keys():
    fail("Source-table feature registry differs from the sealed input")

exact_fields = (
    "module_name",
    "testable",
    "evidence_state",
    "robust_display",
    "linearity_reversal_fibrosis",
    "linearity_reversal_nas",
)
numeric_fields = (
    "beta_meta_fibrosis",
    "beta_meta_nas",
    "q_value_fibrosis",
    "q_value_nas",
)


def missing(value: str) -> bool:
    return value.strip() in {"", "NA", "NaN"}


for feature_id, source_row in source_by_id.items():
    input_row = input_by_id[feature_id]
    for field in exact_fields:
        left = source_row[field]
        right = input_row[field]
        if missing(left) and missing(right):
            continue
        if left != right:
            fail(f"Source table changed {field} for {feature_id}")
    for field in numeric_fields:
        left_raw = source_row[field]
        right_raw = input_row[field]
        if missing(left_raw) and missing(right_raw):
            continue
        if missing(left_raw) != missing(right_raw):
            fail(f"Source table changed {field} missingness for {feature_id}")
        left = float(left_raw)
        right = float(right_raw)
        if not math.isclose(left, right, rel_tol=0, abs_tol=1e-12):
            fail(f"Source table changed {field} for {feature_id}")

testable = [row for row in source_rows if row["testable"] == "TRUE"]
if len(testable) != 113:
    fail(f"Expected 113 testable programs in this sealed candidate; found {len(testable)}")

group_counts: dict[str, int] = {}
for row in testable:
    group_counts[row["display_group"]] = group_counts.get(row["display_group"], 0) + 1
if group_counts != {"Fibrosis": 27, "NAS": 4, "Unsupported": 82}:
    fail(f"Unexpected display-group counts: {group_counts}")

named = [row for row in source_rows if row["robust_display"] == "TRUE"]
if {row["display_label"] for row in named} != {"Ductular injury", "Stromal ECM"}:
    fail("The two frozen display labels changed")

with (output / "figure_summary.tsv").open(newline="") as handle:
    summary_rows = list(csv.DictReader(handle, delimiter="\t"))
summary = {row["metric"]: float(row["value"]) for row in summary_rows}
expected_summary = {
    "frozen_programs": 117,
    "testable_programs": 113,
    "untestable_programs": 4,
    "fibrosis_supported_only": 27,
    "nas_supported_only": 4,
    "both_axes_supported": 0,
    "unsupported": 82,
    "discovery_cohorts": 4,
    "discovery_biological_n": 469,
    "fibrosis_positive_loco_folds": 4,
    "nas_positive_loco_folds": 4,
    "paired_participants": 46,
}
for metric, expected in expected_summary.items():
    if summary.get(metric) != expected:
        fail(f"Unexpected figure summary {metric}: {summary.get(metric)}")
if not math.isclose(summary.get("fibrosis_loco_holm_p", math.nan), 0.0001999800019998,
                    rel_tol=0, abs_tol=1e-15):
    fail("Fibrosis LOCO Holm P changed")
if not math.isclose(summary.get("nas_loco_holm_p", math.nan), 0.0041995800419958,
                    rel_tol=0, abs_tol=1e-15):
    fail("NAS LOCO Holm P changed")
paired_expected = {
    "paired_median_cosine": 0.163569724713045,
    "paired_ci_lower": -0.0326729944275066,
    "paired_ci_upper": 0.243904407298705,
    "paired_empirical_p": 0.217578242175782,
}
for metric, expected in paired_expected.items():
    if not math.isclose(summary.get(metric, math.nan), expected,
                        rel_tol=0, abs_tol=1e-15):
        fail(f"Paired-validation summary changed: {metric}")

with (output / "artifact_manifest.tsv").open(newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))
for row in manifest:
    path = output / row["artifact"]
    if not path.is_file() or sha256(path) != row["sha256"]:
        fail(f"Artifact hash mismatch: {path.name}")
    if path.stat().st_size != int(row["size_bytes"]):
        fail(f"Artifact size mismatch: {path.name}")

def inspect_pdf(path: Path, expected_width: tuple[float, float],
                expected_height: tuple[float, float]) -> str:
    info = subprocess.run(
        ["pdfinfo", str(path)], check=True, text=True, capture_output=True
    ).stdout
    if not re.search(r"^Pages:\s+1$", info, re.MULTILINE):
        fail(f"{path.name} must be a one-page PDF")
    size_match = re.search(
        r"^Page size:\s+([0-9.]+) x ([0-9.]+) pts", info, re.MULTILINE
    )
    if not size_match:
        fail(f"Could not read page size for {path.name}")
    width, height = map(float, size_match.groups())
    if not (expected_width[0] <= width <= expected_width[1] and
            expected_height[0] <= height <= expected_height[1]):
        fail(f"Unexpected canvas for {path.name}: {width} x {height} pt")
    if re.search(rb"/Subtype\s*/Type3\b", path.read_bytes()):
        fail(f"{path.name} contains Type 3 fonts")
    return subprocess.run(
        ["pdftotext", "-layout", str(path), "-"], check=True, text=True,
        capture_output=True
    ).stdout


pdf_text = inspect_pdf(map_pdf, (258.5, 259.8), (215.4, 216.6))
required_text = (
    "Frozen programs reveal broader fibrosis-associated remodeling",
    "No testable program supported on both",
    "both axes reproduced in 4/4 held-out cohorts",
    "Fibrosis association",
    "adjusted for NAS",
    "NAS association",
    "adjusted for fibrosis",
    "Fibrosis only (27)",
    "NAS only (4)",
    "Unsupported (82)",
    "Ductular injury",
    "Stromal ECM",
)
for token in required_text:
    if token not in pdf_text:
        fail(f"Missing required panel text: {token}")
for token in ("Lineage", "Evidence state", "frozen programs", "untestable", "Cross:", "program-score SD"):
    if token in pdf_text:
        fail(f"Technical text leaked back onto compact panel: {token}")

method_text = inspect_pdf(method_pdf, (373.5, 375.0), (82.0, 83.5))
method_required = (
    "Outcome-locked program scoring",
    "no histology-based feature selection or score fitting",
    "117 frozen scRNA programs",
    "fixed genes + L1 weights",
    "L1-weighted bulk score",
    "within-cohort gene z-scores",
    "Mutually adjusted effects",
    "fibrosis | NAS; NAS | fibrosis",
    "Random-effects meta-analysis",
    "4 cohorts; n = 469",
)
for token in method_required:
    if token not in method_text:
        fail(f"Missing required method-panel text: {token}")

print("PASS: impact map, scoring schematic, registry, counts, hashes, PDF mechanics, and text scope")
