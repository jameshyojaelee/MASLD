#!/usr/bin/env python3
"""Freeze outcome-independent GSE238219 leave-one-gene-target-out folds."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from masld_bench.artifacts import verify_frozen_tree


SAMPLES = ("GSM7660623", "GSM7660624", "GSM7660625", "GSM7660626", "GSM7660627")
EXPECTED_TARGETS = (
    "C6orf106", "GPAM", "LYPLAL1", "NCKIPSD", "PNPLA3", "PPP1R3B",
    "RBM6", "TNKS", "TRIB1", "VKORC1", "WDR6",
)


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--materialization-array-job-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("refusing to overwrite target split manifest")
    target_families: dict[str, set[str]] = {target: set() for target in EXPECTED_TARGETS}
    receipts = []
    for index, sample in enumerate(SAMPLES):
        source = args.execution_root / (
            f"gse238219-pseudobulk-{sample}-{args.materialization_array_job_id}_{index}"
        )
        verify_frozen_tree(source)
        receipts.append(
            {
                "sample": sample,
                "artifact_manifest_sha256": digest(source / "ARTIFACTS.json"),
            }
        )
        with (source / "guide_family_summary.csv").open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["classification"] == "single_gene":
                    target = row["target_gene"]
                    if target not in target_families:
                        raise RuntimeError(f"unexpected single-gene target: {target}")
                    target_families[target].add(row["guide_family"])
    if any(not families for families in target_families.values()):
        raise RuntimeError("one or more target genes have no guide family")
    folds = []
    for index, held_target in enumerate(EXPECTED_TARGETS):
        folds.append(
            {
                "fold_id": f"gse238219_unseen_target_{index:02d}_{held_target}",
                "train_target_genes": [target for target in EXPECTED_TARGETS if target != held_target],
                "development_target_genes": [held_target],
                "held_guide_families": sorted(target_families[held_target]),
                "control_role": "observed_context_input_only_not_response_label",
                "excluded_assignment_classes": ["mixed_control_target", "multi_gene", "unassigned"],
            }
        )
    result = {
        "schema_version": "masld-bench-gse238219-unseen-target-folds-v1",
        "dataset_id": "gse238219",
        "split_id": "gene_perturbation_outer",
        "source_materializations": receipts,
        "target_genes": list(EXPECTED_TARGETS),
        "guide_families_by_target": {
            target: sorted(families) for target, families in target_families.items()
        },
        "folds": folds,
        "fold_count": len(folds),
        "split_construction_inputs": "guide-to-target identity and deposited sample identifiers only",
        "source_perturbation_effect_tables_read": False,
        "expression_values_used_to_assign_roles": False,
        "biological_unit_status": "UNRESOLVED",
        "execution_status": "blocked_until_independent_culture_ancestry_and_dataset_activation_are_frozen",
        "permitted_models_after_activation": [
            "bilinear_perturbation", "control_mean", "gears", "genepert",
            "perturbed_mean", "ridge_perturbation", "scgpt_perturbation_head", "systema",
        ],
        "not_applicable_models_for_unseen_targets": ["cellot", "cpa", "perturblib_lpm", "scgen"],
        "gse281160_nearest_gene_mapping_performed": False,
        "gse313774_accessed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "target_split_manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
