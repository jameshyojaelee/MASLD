#!/usr/bin/env python3
"""Read-only integration checks for the continuum manuscript candidate."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
CORE = ROOT / (
    "RNA-seq/results/histology_anchored_continuum/candidates/"
    "hac-continuum-20260818T024923Z"
)
MOLECULAR = ROOT / (
    "RNA-seq/results/histology_anchored_continuum/molecular_layers/"
    "hac-molecular-layers-20260818T173348Z"
)
COMPARISON = ROOT / (
    "figures/candidates/"
    "histology-vs-continuum-comparison-20260818T192630Z"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    require(path.is_file() and path.stat().st_size > 0, f"Missing TSV: {path}")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    core_summary = json.loads(
        (CORE / "validation/validation_summary.json").read_text()
    )
    require(core_summary.get("status") == "PASS", "Core candidate did not pass")

    molecular_summary = json.loads(
        (MOLECULAR / "validation/validation_summary.json").read_text()
    )
    require(
        molecular_summary.get("validation_pass") is True
        and molecular_summary.get("n_failed") == 0,
        "Molecular-layer candidate did not pass",
    )

    index_rows = read_tsv(FIG4 / "CANDIDATE_PANEL_INDEX.tsv")
    row4f = [row for row in index_rows if row["callout"] == "4F"]
    require(len(row4f) == 1, "Figure 4 index must contain exactly one 4F row")
    indexed_4f = (FIG4 / row4f[0]["candidate_relative_path"]).resolve()
    expected_4f = (
        MOLECULAR / "figures/panels/fig4f_continuum_program_trajectories.pdf"
    ).resolve()
    require(indexed_4f == expected_4f, "Figure 4F index points to the wrong candidate")
    require(
        indexed_4f.is_file()
        and indexed_4f.stat().st_size > 10_000
        and indexed_4f.read_bytes()[:4] == b"%PDF",
        "Figure 4F is missing or not a nonempty PDF",
    )

    molecular_figures = read_tsv(MOLECULAR / "figures/figure_manifest.tsv")
    require(molecular_figures, "Molecular figure manifest is empty")
    require(
        all(Path(row["path"]).is_file() for row in molecular_figures),
        "A molecular figure is missing",
    )
    require(
        all(row["current_figure3_modified"] == "FALSE" for row in molecular_figures),
        "A continuum figure claims to modify Figure 3",
    )

    comparison_figures = read_tsv(COMPARISON / "figure_manifest.tsv")
    require(len(comparison_figures) == 11, "Comparison figure family drifted")
    for row in comparison_figures:
        path = Path(row["path"])
        require(path.is_file() and path.stat().st_size > 0, f"Missing figure: {path}")
        require(sha256(path) == row["sha256"], f"Figure hash drift: {path}")
        require(
            row["claim_boundary"]
            == "same-substrate comparison; continuum windows descriptive only",
            f"Claim-boundary drift: {path}",
        )

    print("CONTINUUM_INTEGRATION\tPASS")
    print(f"CORE_VALIDATION\t{core_summary['status']}")
    print(f"MOLECULAR_VALIDATION\t{molecular_summary['validation_pass']}")
    print(f"FIGURE4F\t{indexed_4f}")
    print(f"MOLECULAR_FIGURES\t{len(molecular_figures)}")
    print(f"COMPARISON_FIGURES_HASHED\t{len(comparison_figures)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
