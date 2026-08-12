#!/usr/bin/env python3
"""Assemble and independently validate all leave-all-lncRNA-out program lanes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from validate_program_without_lncrna import validate_program_sensitivity


CELL_TYPES = ("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")
EXPECTED = {
    "hepatocytes": 18,
    "fibroblasts": 11,
    "macrophages": 2,
    "cholangiocytes": 14,
    "tcells": 4,
}
MINIMUM_BASELINE_REPRODUCTION_R = 0.995


def require(value: bool, message: str) -> None:
    if not value:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_manifest(root: Path) -> None:
    manifest = pd.read_csv(root / "output_manifest.tsv", sep="\t")
    require(len(manifest) == 4, f"manifest cardinality drift: {root}")
    for row in manifest.itertuples(index=False):
        path = root / row.relative_path
        require(path.is_file() and not path.is_symlink(), f"invalid artifact: {path}")
        require(path.stat().st_size == int(row.size_bytes), f"size drift: {path}")
        require(sha256(path) == row.sha256, f"checksum drift: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lane-root", required=True, type=Path)
    parser.add_argument("--program-content", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists")
    require(
        args.program_content.is_file() and not args.program_content.is_symlink(),
        "invalid content table",
    )

    summaries = []
    scores = []
    source_rows = []
    for cell_type in CELL_TYPES:
        root = args.lane_root / cell_type
        require(root.is_dir() and not root.is_symlink(), f"missing lane: {cell_type}")
        validate_manifest(root)
        summary = pd.read_csv(root / "program_without_lncrna_sensitivity.tsv", sep="\t")
        donor = pd.read_csv(root / "donor_scores.tsv", sep="\t")
        require(
            len(summary) == EXPECTED[cell_type], f"trigger count drift: {cell_type}"
        )
        require(set(summary["cell_type"]) == {cell_type}, f"lineage drift: {cell_type}")
        summaries.append(summary)
        scores.append(
            donor.rename(columns={"donor": "sample_id"})[
                ["program_uid", "sample_id", "original_score", "without_lncrna_score"]
            ]
        )
        source_rows.append(
            {
                "cell_type": cell_type,
                "path": str(root.resolve()),
                "summary_sha256": sha256(
                    root / "program_without_lncrna_sensitivity.tsv"
                ),
                "donor_scores_sha256": sha256(root / "donor_scores.tsv"),
            }
        )

    summary = pd.concat(summaries, ignore_index=True)
    score = pd.concat(scores, ignore_index=True)
    require(
        len(summary) == 49 and summary["program_uid"].nunique() == 49,
        "complete trigger family drift",
    )
    content = pd.read_csv(args.program_content, sep="\t")
    stage = summary.rename(
        columns={
            "stored_stage_beta": "original_stage_beta",
            "without_lncrna_stage_beta": "without_lncrna_stage_beta",
        }
    )[["program_uid", "original_stage_beta", "without_lncrna_stage_beta"]]
    validated = validate_program_sensitivity(content, score, stage)
    require(len(validated) == 49, "validator family drift")

    joined = summary.merge(
        validated, on=["program_uid", "lncrna_l1_fraction"], validate="one_to_one"
    )
    require(
        (
            joined["sensitivity_pass"].astype(bool)
            == joined["sensitivity_passed"].astype(bool)
        ).all(),
        "producer and independent validator disagree",
    )
    require(
        (joined["without_lncrna_score_r"] >= 0.90).all(),
        "cell-level leave-lncRNA-out correlation gate failed",
    )
    require(
        (joined["pearson_correlation"] >= 0.90).all()
        and (joined["spearman_correlation"] >= 0.90).all(),
        "donor-level leave-lncRNA-out correlation gate failed",
    )
    require(
        (joined["score_reproduction_r"] >= MINIMUM_BASELINE_REPRODUCTION_R).all(),
        "baseline reconstruction gate failed",
    )
    require(
        (joined["reconstructed_stage_beta"] - joined["stored_stage_beta"]).abs().max()
        < 1e-8,
        "stored stage coefficient not reproduced",
    )

    args.output.mkdir(parents=True)
    joined.sort_values("program_uid").to_csv(
        args.output / "program_without_lncrna_sensitivity.tsv", sep="\t", index=False
    )
    pd.DataFrame(source_rows).to_csv(
        args.output / "source_manifest.tsv", sep="\t", index=False
    )
    payload = {
        "status": "pass",
        "n_triggered_programs": 49,
        "n_sensitivity_passed": int(joined["sensitivity_passed"].sum()),
        "n_sensitivity_failed": int((~joined["sensitivity_passed"].astype(bool)).sum()),
        "minimum_baseline_reproduction_r": float(joined["score_reproduction_r"].min()),
        "minimum_without_lncrna_cell_pearson_r": float(
            joined["without_lncrna_score_r"].min()
        ),
        "minimum_without_lncrna_donor_pearson_r": float(
            joined["pearson_correlation"].min()
        ),
        "minimum_without_lncrna_donor_spearman_r": float(
            joined["spearman_correlation"].min()
        ),
        "n_stage_direction_preserved": int(joined["stage_direction_preserved_y"].sum()),
        "biological_unit": "biological_donor",
        "module_discovery": False,
        "scoring": "exact_stored_score_plus_reconstructed_gene_removal_delta",
        "baseline_gate_status": "technical_threshold_finalized_after_source_reconstruction_audit",
    }
    (args.output / "VALIDATED.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
