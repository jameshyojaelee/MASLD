"""Per-stage output schema validator. Called at end of each script.

Usage:
    python validate.py <stage> <cell_type>
    e.g. python validate.py hotspot cholangiocytes
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import hotspot_outdir, RESULTS  # noqa: E402

import pandas as pd  # noqa: E402

REQUIRED_COLS = {
    "module_genes":   {"gene", "module", "weight"},
    "cell_scores":    {"cell_id", "module", "score"},   # parquet
    "donor_scores":   {"sample", "module", "score"},
    "autocorr":       {"gene", "Z", "FDR"},
}


def check_file(path: Path, expected_cols: set, sep: str = "\t") -> None:
    assert path.exists(), f"Missing: {path}"
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep=sep)
    missing = expected_cols - set(df.columns)
    assert not missing, f"{path.name} missing columns: {missing}"
    assert len(df) > 0, f"{path.name} is empty"


def validate_hotspot_stage(cell_type: str) -> None:
    outdir = hotspot_outdir(cell_type)
    for name, cols in REQUIRED_COLS.items():
        ext = ".parquet" if name == "cell_scores" else ".tsv"
        check_file(outdir / f"{name}{ext}", cols)
    print(f"[OK] hotspot outputs for {cell_type}")


STAGE_DISPATCH = {
    "hotspot": validate_hotspot_stage,
}

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Validate Hotspot pipeline outputs")
    ap.add_argument("stage", choices=list(STAGE_DISPATCH))
    ap.add_argument("cell_type")
    args = ap.parse_args()
    STAGE_DISPATCH[args.stage](args.cell_type)
