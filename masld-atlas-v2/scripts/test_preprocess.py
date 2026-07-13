#!/usr/bin/env python3
"""
test_preprocess.py
==================
Validation tests for the preprocessed MASLD Atlas data files.

Checks:
  1. atlas.parquet: exists, >30K rows, >=80 cols, contains THRB
  2. gene_index.json: sorted list, >30K entries, required keys + evidence s1-s7
  3. atlas_summary.json: total_genes >30K, total_cohorts==10, total_samples==1444
  4. featured_genes.json: >=4 genes, THRB featured, each has evidence + tagline

Usage:
  python test_preprocess.py --data-dir ../public/data
"""

import argparse
import json
import os
import sys
from pathlib import Path


class TestResult:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.errors = []

    def ok(self, name: str, detail: str = ""):
        self.passed += 1
        msg = f"  PASS: {name}"
        if detail:
            msg += f" ({detail})"
        print(msg)

    def fail(self, name: str, detail: str = ""):
        self.failed += 1
        msg = f"  FAIL: {name}"
        if detail:
            msg += f" -- {detail}"
        print(msg)
        self.errors.append(msg)

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'='*60}")
        print(f"Results: {self.passed}/{total} passed, {self.failed} failed")
        if self.errors:
            print("Failures:")
            for e in self.errors:
                print(f"  {e}")
        print(f"{'='*60}")
        return self.failed == 0


def test_atlas_parquet(data_dir: Path, result: TestResult):
    """Test atlas.parquet file."""
    print("\n--- atlas.parquet ---")
    fpath = data_dir / "atlas.parquet"

    if not fpath.exists():
        result.fail("atlas.parquet exists", "File not found")
        return
    result.ok("atlas.parquet exists", f"{fpath.stat().st_size / 1024 / 1024:.1f} MB")

    try:
        import pyarrow.parquet as pq

        table = pq.read_table(fpath)
        nrows = table.num_rows
        ncols = table.num_columns
        col_names = table.column_names
    except ImportError:
        # Fallback: try pandas
        import pandas as pd

        df = pd.read_parquet(fpath)
        nrows = len(df)
        ncols = len(df.columns)
        col_names = list(df.columns)

    if nrows > 30000:
        result.ok("rows > 30K", f"{nrows:,}")
    else:
        result.fail("rows > 30K", f"got {nrows:,}")

    if ncols >= 80:
        result.ok("cols >= 80", f"{ncols}")
    else:
        result.fail("cols >= 80", f"got {ncols}")

    if "human_symbol" in col_names:
        # Check for THRB
        try:
            import pyarrow.parquet as pq

            df = pq.read_table(fpath, columns=["human_symbol"]).to_pandas()
        except ImportError:
            import pandas as pd

            df = pd.read_parquet(fpath, columns=["human_symbol"])

        has_thrb = (df["human_symbol"] == "THRB").any()
        if has_thrb:
            result.ok("contains THRB")
        else:
            result.fail("contains THRB", "THRB not found in human_symbol column")
    else:
        result.fail("contains THRB", "human_symbol column missing")


def test_gene_index(data_dir: Path, result: TestResult):
    """Test gene_index.json file."""
    print("\n--- gene_index.json ---")
    fpath = data_dir / "gene_index.json"

    if not fpath.exists():
        result.fail("gene_index.json exists", "File not found")
        return
    result.ok("gene_index.json exists", f"{fpath.stat().st_size / 1024 / 1024:.2f} MB")

    with open(fpath, "r") as f:
        data = json.load(f)

    if not isinstance(data, list):
        result.fail("is a list", f"got {type(data).__name__}")
        return
    result.ok("is a list")

    if len(data) > 30000:
        result.ok("entries > 30K", f"{len(data):,}")
    else:
        result.fail("entries > 30K", f"got {len(data):,}")

    # Check sorting
    symbols = [d["symbol"] for d in data]
    is_sorted = all(symbols[i] <= symbols[i + 1] for i in range(len(symbols) - 1))
    if is_sorted:
        result.ok("sorted by symbol")
    else:
        result.fail("sorted by symbol")

    # Check required keys on first entry
    # Core fields always present; optional fields (sex_class, zonation_class,
    # ferroptosis_class, is_deg, is_conserved, dgidb_druggable,
    # layers_active) may be omitted when null/false/zero for compactness
    always_present_keys = {
        "symbol",
        "ensembl_id",
        "biotype",
        "bulk_logfc",
        "bulk_padj",
        "evidence",
    }
    # Optional keys that may appear depending on the gene
    optional_keys = {
        "is_deg",
        "is_conserved",
        "sex_class",
        "zonation_class",
        "ferroptosis_class",
        "dgidb_druggable",
        "layers_active",
    }
    evidence_keys = {
        "s1_human",
        "s2_genetic",
        "s3_essential",
        "s4_epigenomic",
        "s5_spatial",
        "s6_singlecell",
        "s7_mouse",
    }

    sample_entry = data[0]
    missing_core = always_present_keys - set(sample_entry.keys())
    if not missing_core:
        result.ok("required core keys present")
    else:
        result.fail("required core keys present", f"missing: {missing_core}")

    # Check optional keys appear on at least some entries
    all_keys_seen = set()
    for entry in data[:5000]:
        all_keys_seen.update(entry.keys())
    expected_optional_seen = optional_keys & all_keys_seen
    if len(expected_optional_seen) >= 3:
        result.ok("optional keys appear on some entries", f"seen: {expected_optional_seen}")
    else:
        result.fail("optional keys appear on some entries", f"only seen: {expected_optional_seen}")

    if "evidence" in sample_entry and isinstance(sample_entry["evidence"], dict):
        # Evidence keys may be omitted when 0; check that valid keys are a subset
        ev_keys_present = set(sample_entry["evidence"].keys())
        invalid_ev = ev_keys_present - evidence_keys
        if not invalid_ev:
            result.ok("evidence keys valid (subset of s1-s7)")
        else:
            result.fail("evidence keys valid", f"unexpected keys: {invalid_ev}")

        # Check that at least some entries have non-empty evidence
        entries_with_evidence = sum(1 for d in data[:5000] if d.get("evidence"))
        if entries_with_evidence > 100:
            result.ok("many entries have non-empty evidence", f"{entries_with_evidence}/5000 sample")
        else:
            result.fail("many entries have non-empty evidence", f"only {entries_with_evidence}/5000")
    else:
        result.fail("evidence field present", "evidence field missing or not a dict")

    # Check no NaN strings in JSON (should be null)
    raw_text = fpath.read_text()
    nan_occurrences = raw_text.count(":NaN") + raw_text.count(": NaN")
    if nan_occurrences == 0:
        result.ok("no NaN values in JSON")
    else:
        result.fail("no NaN values in JSON", f"found {nan_occurrences} NaN occurrences")

    # Check evidence values are in [0, 1] range
    bad_range = 0
    for entry in data[:1000]:  # spot-check first 1000
        if "evidence" in entry:
            for k, v in entry["evidence"].items():
                if v is not None and (v < 0 or v > 1):
                    bad_range += 1
    if bad_range == 0:
        result.ok("evidence values in [0,1] range (sample check)")
    else:
        result.fail(
            "evidence values in [0,1] range",
            f"{bad_range} out-of-range values in first 1000 entries",
        )


def test_atlas_summary(data_dir: Path, result: TestResult):
    """Test atlas_summary.json file."""
    print("\n--- atlas_summary.json ---")
    fpath = data_dir / "atlas_summary.json"

    if not fpath.exists():
        result.fail("atlas_summary.json exists", "File not found")
        return
    result.ok(
        "atlas_summary.json exists", f"{fpath.stat().st_size / 1024:.1f} KB"
    )

    with open(fpath, "r") as f:
        data = json.load(f)

    if data.get("total_genes", 0) > 30000:
        result.ok("total_genes > 30K", f"{data['total_genes']:,}")
    else:
        result.fail("total_genes > 30K", f"got {data.get('total_genes')}")

    if data.get("total_cohorts") == 10:
        result.ok("total_cohorts == 10")
    else:
        result.fail("total_cohorts == 10", f"got {data.get('total_cohorts')}")

    if data.get("total_samples") == 1444:
        result.ok("total_samples == 1444")
    else:
        result.fail("total_samples == 1444", f"got {data.get('total_samples')}")

    if data.get("total_degs", 0) > 5000:
        result.ok("total_degs > 5K", f"{data.get('total_degs'):,}")
    else:
        result.fail("total_degs > 5K", f"got {data.get('total_degs')}")

    if data.get("evidence_sources") == 7:
        result.ok("evidence_sources == 7")
    else:
        result.fail("evidence_sources == 7", f"got {data.get('evidence_sources')}")

    # Check all expected keys
    expected_keys = {
        "total_genes",
        "total_degs",
        "total_cohorts",
        "total_samples",
        "conserved_count",
        "coloc_genes",
        "drug_targets",
        "mouse_datasets",
        "evidence_sources",
    }
    missing = expected_keys - set(data.keys())
    if not missing:
        result.ok("all expected keys present")
    else:
        result.fail("all expected keys present", f"missing: {missing}")


def test_featured_genes(data_dir: Path, result: TestResult):
    """Test featured_genes.json file."""
    print("\n--- featured_genes.json ---")
    fpath = data_dir / "featured_genes.json"

    if not fpath.exists():
        result.fail("featured_genes.json exists", "File not found")
        return
    result.ok(
        "featured_genes.json exists", f"{fpath.stat().st_size / 1024:.1f} KB"
    )

    with open(fpath, "r") as f:
        data = json.load(f)

    if not isinstance(data, list):
        result.fail("is a list", f"got {type(data).__name__}")
        return

    if len(data) >= 4:
        result.ok("has >= 4 genes", f"{len(data)}")
    else:
        result.fail("has >= 4 genes", f"got {len(data)}")

    # Check THRB is featured
    symbols = {d.get("symbol") for d in data}
    if "THRB" in symbols:
        result.ok("THRB is featured")
    else:
        result.fail("THRB is featured", f"symbols: {symbols}")

    # Check each has evidence and tagline
    all_have_evidence = all("evidence" in d for d in data)
    all_have_tagline = all("tagline" in d for d in data)

    if all_have_evidence:
        result.ok("all entries have evidence")
    else:
        result.fail("all entries have evidence")

    if all_have_tagline:
        result.ok("all entries have tagline")
    else:
        result.fail("all entries have tagline")

    # Check evidence keys are valid (subset of s1-s7; zero-valued keys may be omitted)
    ev_keys = {"s1_human", "s2_genetic", "s3_essential", "s4_epigenomic", "s5_spatial", "s6_singlecell", "s7_mouse"}
    for entry in data:
        if "evidence" in entry:
            present = set(entry["evidence"].keys())
            invalid = present - ev_keys
            if invalid:
                result.fail(
                    f"{entry.get('symbol', '?')} evidence keys valid",
                    f"unexpected: {invalid}",
                )
            else:
                result.ok(f"{entry.get('symbol', '?')} evidence keys valid ({len(present)} active)")


def main():
    parser = argparse.ArgumentParser(description="Validate preprocessed MASLD data.")
    parser.add_argument(
        "--data-dir",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "..", "public", "data"),
        help="Directory containing preprocessed data files (default: ../public/data)",
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    print(f"Testing data in: {data_dir}")

    result = TestResult()

    test_atlas_parquet(data_dir, result)
    test_gene_index(data_dir, result)
    test_atlas_summary(data_dir, result)
    test_featured_genes(data_dir, result)

    success = result.summary()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
