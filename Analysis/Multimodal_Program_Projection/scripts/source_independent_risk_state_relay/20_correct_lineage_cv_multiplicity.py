#!/usr/bin/env python3
"""Non-overwriting multiplicity correction for the sealed Plan 45 lineage CV."""

from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from pathlib import Path

from relay_common import PROJECT_ROOT, sha256_file


UPSTREAM_ID = "source-independent-risk-state-relay-lineage-disease-reference-cv-2026-08-10"
CANDIDATE_ID = os.environ.get(
    "PLAN45_DISEASE_CV_CORRECTION_ID",
    "source-independent-risk-state-relay-lineage-disease-reference-cv-fdr-correction-2026-08-10",
).strip()
CANDIDATE_PARENT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
UPSTREAM = CANDIDATE_PARENT / UPSTREAM_ID
ROOT = CANDIDATE_PARENT / CANDIDATE_ID


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, data: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def bh(pvalues: list[float]) -> list[float]:
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    m = len(pvalues)
    for reverse_index in range(m - 1, -1, -1):
        original = order[reverse_index]
        running = min(running, pvalues[original] * m / (reverse_index + 1))
        adjusted[original] = min(1.0, running)
    return adjusted


def main() -> None:
    if ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite correction release: {ROOT}")
    upstream_seal_path = UPSTREAM / "LINEAGE_DISEASE_CV_SEALED.json"
    upstream_effects_path = UPSTREAM / "cv_fold_effects.tsv"
    upstream_gate_path = UPSTREAM / "cv_gate_status.tsv"
    seal = json.loads(upstream_seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_cross_cohort_lineage_disease_reference_validation":
        raise RuntimeError("Upstream CV release is not sealed")
    for key in (
        "experimental_outcomes_inspected", "experimental_targets_frozen",
        "failed_plan43_or_plan44_score_reused",
    ):
        if seal.get(key) is not False:
            raise RuntimeError(f"Upstream firewall violation: {key}")

    effects = read_rows(upstream_effects_path)
    primary = [
        row for row in effects
        if row["heldout_dataset"] == "GSE244832"
        and row["weighting_rule"] == "continuous_all"
    ]
    if len(primary) != 5 or len({row["cell_type"] for row in primary}) != 5:
        raise RuntimeError("Primary five-lineage multiplicity family drift")
    qvalues = bh([float(row["p_one_sided_positive"]) for row in primary])
    q_by_lineage = {row["cell_type"]: q for row, q in zip(primary, qvalues)}

    corrected_effects = []
    for row in effects:
        new = dict(row)
        new["primary_gse244832_lineage_family_q"] = (
            f"{q_by_lineage[row['cell_type']]:.17g}"
            if row["heldout_dataset"] == "GSE244832"
            and row["weighting_rule"] == "continuous_all"
            else ""
        )
        corrected_effects.append(new)

    corrected_gates = []
    for row in read_rows(upstream_gate_path):
        new = dict(row)
        qvalue = q_by_lineage[row["cell_type"]]
        new["gse244832_lineage_family_q"] = f"{qvalue:.17g}"
        corrected_pass = (
            float(row["gse244832_effect"]) > 0
            and qvalue < 0.05
            and float(row["gse174748_effect"]) > 0
            and row["all_primary_leave_one_positive"].lower() == "true"
            and row["all_weighting_sensitivities_positive"].lower() == "true"
            and int(row["minimum_nonzero_loadings"]) >= 500
        )
        new["multiplicity_corrected_cv_gate_pass"] = str(corrected_pass).upper()
        new["multiplicity_correction"] = (
            "BH across five primary continuous-all GSE244832 lineage tests"
        )
        new["interpretation"] = (
            "cross_cohort_transport_supported_before_experimental_outcomes"
            if corrected_pass
            else "lineage_disease_projection_not_validated_for_primary_relay_endpoint"
        )
        corrected_gates.append(new)

    CANDIDATE_PARENT.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(prefix=f".{CANDIDATE_ID}.", dir=CANDIDATE_PARENT))
    try:
        effects_out = temp_root / "cv_fold_effects_multiplicity_corrected.tsv"
        gates_out = temp_root / "cv_gate_status_multiplicity_corrected.tsv"
        manifest_out = temp_root / "cv_multiplicity_correction_manifest.tsv"
        write_rows(effects_out, corrected_effects, list(corrected_effects[0]))
        write_rows(gates_out, corrected_gates, list(corrected_gates[0]))
        manifest = [
            {
                "source_id": source_id,
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for source_id, path in (
                ("upstream_cv_seal", upstream_seal_path),
                ("upstream_cv_effects", upstream_effects_path),
                ("upstream_cv_gate", upstream_gate_path),
            )
        ]
        write_rows(manifest_out, manifest, list(manifest[0]))
        outputs = {
            "cv_fold_effects_multiplicity_corrected": effects_out,
            "cv_gate_status_multiplicity_corrected": gates_out,
            "cv_multiplicity_correction_manifest": manifest_out,
        }
        correction_seal = {
            "status": "sealed_lineage_cv_multiplicity_correction",
            "candidate_id": CANDIDATE_ID,
            "upstream_candidate_id": UPSTREAM_ID,
            "reason": (
                "The first sealed CV gate used per-lineage p<0.05 without correcting "
                "the five-lineage primary family; correction occurred before target "
                "selection or experimental outcome access."
            ),
            "primary_family": "five continuous-all GSE244832 held-out lineage tests",
            "method": "Benjamini-Hochberg",
            "threshold": 0.05,
            "n_lineages_passing": sum(
                row["multiplicity_corrected_cv_gate_pass"] == "TRUE"
                for row in corrected_gates
            ),
            "experimental_outcomes_inspected": False,
            "experimental_targets_frozen": False,
            "output_sha256": {name: sha256_file(path) for name, path in outputs.items()},
        }
        seal_path = temp_root / "LINEAGE_DISEASE_CV_FDR_CORRECTION_SEALED.json"
        seal_path.write_text(
            json.dumps(correction_seal, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temp_root.rename(ROOT)
    except Exception:
        if temp_root.exists():
            shutil.rmtree(temp_root)
        raise
    print(
        "Plan 45 CV multiplicity correction sealed: "
        f"lineages_passing={correction_seal['n_lineages_passing']}/5; targets not frozen"
    )


if __name__ == "__main__":
    main()
