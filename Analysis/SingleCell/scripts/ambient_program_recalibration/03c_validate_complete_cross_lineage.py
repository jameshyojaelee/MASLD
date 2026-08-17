#!/usr/bin/env python3
"""Independently validate complete-atlas six-lineage localization outputs."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from common import MEMBERSHIP, READY, REGISTRY, bh_adjust, refuse_existing, require, sha256


SCRIPT_ROOT = Path(__file__).resolve().parent
CANDIDATE = Path(os.environ["CROSS_CAND_ROOT"])
RESULTS = CANDIDATE / "results"


def load_validator_helpers():
    path = SCRIPT_ROOT / "04_validate.py"
    specification = importlib.util.spec_from_file_location("ambient_validator_helpers", path)
    require(
        specification is not None and specification.loader is not None,
        "cannot load independent validation helpers",
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


helpers = load_validator_helpers()


def check_close(observed: float, expected: float, label: str) -> None:
    if np.isnan(observed) and np.isnan(expected):
        return
    require(
        np.isclose(observed, expected, atol=1e-8, rtol=1e-8),
        f"numeric validation failed: {label}: {observed} != {expected}",
    )


def main() -> None:
    validation_path = RESULTS / "validation_report.tsv"
    manifest_path = RESULTS / "release_manifest.tsv"
    validated_path = CANDIDATE / "VALIDATED"
    for path in (validation_path, manifest_path, validated_path):
        refuse_existing(path)
    rows = []

    def add(check: str, passed: bool, detail: object) -> None:
        rows.append({"check": check, "passed": bool(passed), "detail": str(detail)})
        require(bool(passed), f"validation failed: {check}: {detail}")

    ready = pd.read_csv(READY, sep="\t", dtype=str)
    add("registry_hash_matches_ready", sha256(REGISTRY) == ready.loc[0, "registry_sha256"], sha256(REGISTRY))
    add("membership_hash_matches_ready", sha256(MEMBERSHIP) == ready.loc[0, "membership_table_sha256"], sha256(MEMBERSHIP))
    scores = pd.read_csv(RESULTS / "hero_lineage_scores.tsv.gz", sep="\t")
    contrasts = pd.read_csv(RESULTS / "hero_lineage_contrasts.tsv", sep="\t")
    lodo = pd.read_csv(RESULTS / "hero_lineage_lodo.tsv", sep="\t")
    promotion = pd.read_csv(RESULTS / "hero_lineage_promotion.tsv", sep="\t")
    census = pd.read_csv(RESULTS / "hero_lineage_census.tsv", sep="\t")
    coverage = pd.read_csv(RESULTS / "hero_lineage_weight_coverage.tsv", sep="\t")
    provenance = json.loads((RESULTS / "hero_lineage_provenance.json").read_text())
    input_provenance = pd.read_csv(
        RESULTS / "complete_atlas_input_provenance.tsv", sep="\t"
    )
    complete_atlas = Path(
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
    )
    add(
        "complete_atlas_is_annotation_authority",
        Path(provenance["atlas"]) == complete_atlas.resolve()
        and provenance["atlas_sha256"] == sha256(complete_atlas),
        provenance["atlas"],
    )
    provenance_hashes_exact = True
    for row in input_provenance.itertuples(index=False):
        path = Path(row.path)
        provenance_hashes_exact &= path.exists() and sha256(path) == row.sha256
    add(
        "complete_input_provenance_hashes_exact",
        provenance_hashes_exact,
        input_provenance["role"].tolist(),
    )
    add("contrast_rows_40", len(contrasts) == 40, len(contrasts))
    add("promotion_rows_2", len(promotion) == 2, len(promotion))
    add(
        "six_lineages_present",
        set(census["lineage"])
        == {
            "Hepatocytes",
            "Fibroblasts",
            "Cholangiocytes",
            "Endothelial cells",
            "Macrophages",
            "T cells",
        },
        sorted(census["lineage"].unique()),
    )
    primary = contrasts[
        (contrasts["annotation_filter"] == "all_annotated_cells")
        & (contrasts["minimum_cells"] == 50)
    ]
    add("primary_family_10", len(primary) == 10, len(primary))
    tcell = primary[primary["comparison_lineage"] == "T cells"]
    add(
        "tcell_comparisons_estimable",
        len(tcell) == 2 and tcell["beta"].notna().all(),
        tcell[["program_name", "beta", "n_paired_donors", "n_datasets"]].to_dict("records"),
    )
    for family, group in contrasts.groupby(["annotation_filter", "minimum_cells"]):
        add(f"family_10:{family}", len(group) == 10, len(group))
        add(
            f"BH_exact:{family}",
            np.allclose(
                bh_adjust(group["pvalue"]), group["qvalue"], equal_nan=True, atol=1e-12
            ),
            10,
        )
        add(
            f"HC3_BH_exact:{family}",
            np.allclose(
                bh_adjust(group["hc3_pvalue"]),
                group["hc3_qvalue"],
                equal_nan=True,
                atol=1e-12,
            ),
            10,
        )
    for result in contrasts.itertuples(index=False):
        subset = scores[
            (scores["annotation_filter"] == result.annotation_filter)
            & (scores["minimum_cells"] == result.minimum_cells)
        ]
        try:
            fit = helpers.independent_paired_fit(
                subset, result.program_uid, result.comparison_lineage
            )
        except RuntimeError:
            add(
                f"unestimable_matches:{result.program_uid}:{result.annotation_filter}:{result.minimum_cells}:{result.comparison_lineage}",
                pd.isna(result.beta),
                result.beta,
            )
            continue
        for metric in ("beta", "se", "pvalue", "hc3_se", "hc3_pvalue"):
            check_close(
                float(getattr(result, metric)),
                float(fit[metric]),
                f"{result.program_uid}:{result.annotation_filter}:{result.minimum_cells}:{result.comparison_lineage}:{metric}",
            )
    add("all_40_contrasts_recomputed", True, 40)
    add(
        "coverage_gates_complete",
        coverage["n_observed_genes"].ge(8).all()
        and coverage["retained_l1_weight"].ge(0.80).all(),
        len(coverage),
    )
    family_eligible = bool(promotion["program_pass"].all())
    add(
        "promotion_destination_consistent",
        set(promotion["main_figure_family_eligible"].astype(bool)) == {family_eligible}
        and set(promotion["destination"])
        == ({"Figure_4E"} if family_eligible else {"Figure_S4"}),
        promotion[["program_name", "program_pass", "destination"]].to_dict("records"),
    )
    for result in promotion.itertuples(index=False):
        target_lodo = lodo[
            (lodo["program_uid"] == result.program_uid)
            & (lodo["comparison_lineage"] == result.target_lineage)
            & (lodo["annotation_filter"] == "all_annotated_cells")
            & (lodo["minimum_cells"] == 50)
            & lodo["beta"].notna()
        ]
        add(
            f"target_lodo_positive:{result.program_uid}",
            len(target_lodo) >= 3 and (target_lodo["beta"].astype(float) > 0).all(),
            target_lodo[["held_out_dataset", "beta"]].to_dict("records"),
        )

    validation = pd.DataFrame(rows)
    validation.to_csv(validation_path, sep="\t", index=False)
    candidate_files = sorted(
        path
        for path in RESULTS.rglob("*")
        if path.is_file() and path != manifest_path
    )
    release = pd.DataFrame(
        [
            {
                "relative_path": str(path.relative_to(CANDIDATE)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in candidate_files
        ]
    )
    release.to_csv(manifest_path, sep="\t", index=False)
    validated_path.write_text(
        "status\tvalidated_complete_six_lineage_candidate\n"
        f"validation_report_sha256\t{sha256(validation_path)}\n"
        f"release_manifest_sha256\t{sha256(manifest_path)}\n",
        encoding="utf-8",
    )
    print(f"PASS {len(rows)} complete-atlas cross-lineage checks")


if __name__ == "__main__":
    main()
