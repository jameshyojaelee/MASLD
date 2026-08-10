#!/usr/bin/env python3
"""Build the native-effect support matrix and mechanical Figure 5 verdict."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from public_functional_common import CANDIDATE_ROOT, atomic_write_text, read_tsv, require_sealed


DATASETS = ["GSE200418", "GSE207889", "GSE253380", "GSE106737"]
PRIMARY_SCHEMES = ["weighted", "equal", "leave_top"]


def read_frames(filename: str, required: bool = True) -> pd.DataFrame:
    frames = []
    for dataset in DATASETS:
        path = CANDIDATE_ROOT / "analyses" / dataset / filename
        if path.is_file():
            frame = pd.read_csv(path, sep="\t")
            if "dataset_id" not in frame.columns:
                frame["dataset_id"] = dataset
            frames.append(frame)
        elif required:
            raise FileNotFoundError(path)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def single(frame: pd.DataFrame, **filters: object) -> pd.Series:
    selected = frame
    for column, value in filters.items():
        selected = selected[selected[column] == value]
    if len(selected) != 1:
        raise RuntimeError(f"Expected one row for {filters}, found {len(selected)}")
    return selected.iloc[0]


def all_effects_direction(frame: pd.DataFrame, uid: str, dataset: str, contrast: str, expected: int, lineage: str | None = None) -> bool:
    selected = frame[(frame.program_uid == uid) & (frame.dataset_id == dataset) & (frame.contrast_id == contrast) & frame.scoring_scheme.isin(PRIMARY_SCHEMES)]
    if lineage is not None:
        selected = selected[selected.lineage == lineage]
    return len(selected) == 3 and bool(np.all(np.sign(selected.estimate.astype(float)) == expected))


def main() -> None:
    seal = require_sealed()
    programs = read_frames("program_effects.tsv")
    testability = read_frames("program_testability.tsv")
    sensitivities = read_frames("sensitivity.tsv")
    classes = read_frames("evidence_class_effects.tsv", required=False)
    class_sensitivities = read_frames("evidence_class_sensitivity.tsv", required=False)
    primary_programs = pd.DataFrame(read_tsv(CANDIDATE_ROOT / "frozen_inputs/external_test_programs.tsv"))
    primary_uids = set(primary_programs.program_uid)
    gates = pd.DataFrame(read_tsv(CANDIDATE_ROOT / "source_gate_status.tsv"))
    source_manifest_frames = []
    source_patterns = ["sources/*/source_manifest.tsv", "sources/*/source_amendment_manifest.tsv", "sources/*/annotation/annotation_source_manifest.tsv"]
    for pattern in source_patterns:
        for path in sorted(CANDIDATE_ROOT.glob(pattern)):
            frame = pd.read_csv(path, sep="\t")
            frame["manifest_path"] = str(path.relative_to(CANDIDATE_ROOT))
            source_manifest_frames.append(frame)
    if not source_manifest_frames:
        raise RuntimeError("No source provenance manifests found")
    source_manifest = pd.concat(source_manifest_frames, ignore_index=True, sort=False)

    support = programs[(programs.program_uid.isin(primary_uids)) & (programs.scoring_scheme == "weighted")].copy()
    support["object_type"] = "program"
    support["object_id"] = support.program_uid
    support["evidence_state"] = np.where(
        (support.primary_family_q.astype(float) < 0.05) & np.isfinite(support.primary_family_q.astype(float)),
        "nominal_family_supported_not_yet_cross_assay_adjudicated",
        "indeterminate",
    )
    support["source_dependence"] = np.where(support.dataset_id == "GSE253380", "model_specific_no_public_line_key", "source_independent_within_deposited_design")

    verdict_rows = []
    for uid in sorted(primary_uids):
        pcls = single(programs, program_uid=uid, dataset_id="GSE200418", contrast_id="PCLS_GFIPO_GFI_48h", scoring_scheme="weighted")
        pcls_gate = float(pcls.primary_family_q) < 0.05 and float(pcls.estimate) > 0 and int(pcls.n_positive_donors) >= 5
        score_sensitivity = all_effects_direction(programs, uid, "GSE200418", "PCLS_GFIPO_GFI_48h", 1)
        primary_all_schemes = programs[(programs.program_uid == uid) & (programs.dataset_id == "GSE200418") & (programs.contrast_id == "PCLS_GFIPO_GFI_48h") & programs.scoring_scheme.isin(PRIMARY_SCHEMES)]
        score_sensitivity = score_sensitivity and len(primary_all_schemes) == 3 and bool(np.all(primary_all_schemes.primary_family_q.astype(float) < 0.05))
        lodo = sensitivities[(sensitivities.dataset_id == "GSE200418") & (sensitivities.program_uid == uid) & (sensitivities.contrast_id == "PCLS_GFIPO_GFI_48h") & (sensitivities.sensitivity == "leave_one_donor_out")]
        culture = single(programs, program_uid=uid, dataset_id="GSE200418", contrast_id="PCLS_CTR_48h_24h", scoring_scheme="weighted")
        pcls_robust = score_sensitivity and len(lodo) == 7 and bool(np.all(lodo.estimate.astype(float) > 0)) and abs(float(pcls.estimate)) > abs(float(culture.estimate))

        schlo_gate = True
        for contrast in ["SCHLO_PA_CONTROL", "SCHLO_OA_CONTROL"]:
            row = single(programs, program_uid=uid, dataset_id="GSE207889", contrast_id=contrast, scoring_scheme="weighted", lineage="Hepatocytes")
            schlo_gate = schlo_gate and float(row.estimate) > 0 and int(row.n_direction_agree) == 2

        biopsy = single(programs, program_uid=uid, dataset_id="GSE106737", contrast_id="LSI_RESPONSE_DIFF", scoring_scheme="weighted")
        biopsy_schemes = all_effects_direction(programs, uid, "GSE106737", "LSI_RESPONSE_DIFF", -1)
        lopo = sensitivities[(sensitivities.dataset_id == "GSE106737") & (sensitivities.program_uid == uid) & (sensitivities.contrast_id == "LSI_RESPONSE_DIFF") & (sensitivities.sensitivity == "leave_one_participant_out")]
        biopsy_gate = float(biopsy.primary_family_q) < 0.05 and float(biopsy.estimate) < 0 and biopsy_schemes and len(lopo) == 20 and bool(np.all(lopo.estimate.astype(float) < 0))

        rygb_gate = all_effects_direction(programs, uid, "GSE106737", "RYGB_RESPONSE", -1)
        acmsd_gate = all(all_effects_direction(programs, uid, "GSE253380", contrast, -1) for contrast in ["ACMSD_CC_INTERACTION", "ACMSD_TT_INTERACTION"])
        orthogonal_gate = rygb_gate and acmsd_gate
        provenance_gate = (
            single(gates, dataset_id="GSE200418").inference_authorized == "true"
            and single(gates, dataset_id="GSE207889").status == "pass_with_source_limit"
            and single(gates, dataset_id="GSE106737").inference_authorized == "true"
            and single(gates, dataset_id="GSE253380").status == "pass_model_specific_no_line_key"
        )
        eligible = all([pcls_gate, pcls_robust, schlo_gate, biopsy_gate, orthogonal_gate, provenance_gate])
        verdict_rows.append(
            {
                "object_type": "program", "object_id": uid,
                "PCLS_INDUCTION": str(pcls_gate).lower(), "PCLS_ROBUSTNESS": str(pcls_robust).lower(),
                "SCHLO_LOCALIZATION": str(schlo_gate).lower(), "HUMAN_REVERSAL": str(biopsy_gate).lower(),
                "ORTHOGONAL_REVERSAL": str(orthogonal_gate).lower(), "PROVENANCE": str(provenance_gate).lower(),
                "main_figure_eligible": str(eligible).lower(),
                "release_terminal_state": "accepted_main" if eligible else "accepted_supplement",
                "reason": "all_six_preregistered_gates_pass" if eligible else "one_or_more_preregistered_gates_failed",
                "specification_sha256": seal["specification_sha256"],
            }
        )

    for evidence_class in ["disease_state_only", "genetic_only", "convergent"]:
        verdict_rows.append(
            {
                "object_type": "evidence_class", "object_id": evidence_class,
                "PCLS_INDUCTION": "evaluated", "PCLS_ROBUSTNESS": "matched_sensitivity_evaluated",
                "SCHLO_LOCALIZATION": "not_prespecified_for_class_inference", "HUMAN_REVERSAL": "evaluated",
                "ORTHOGONAL_REVERSAL": "CLCC1_skipped_source_unavailable;Myojin_reuse_only", "PROVENANCE": "true",
                "main_figure_eligible": "false", "release_terminal_state": "accepted_supplement",
                "reason": "evidence_class_does_not_have_complete_lineage_to_direct_phenotype_chain",
                "specification_sha256": seal["specification_sha256"],
            }
        )
    any_program = any(row["main_figure_eligible"] == "true" for row in verdict_rows if row["object_type"] == "program")
    verdict_rows.append(
        {
            "object_type": "workstream", "object_id": "public_functional_map",
            "PCLS_INDUCTION": "see_program_rows", "PCLS_ROBUSTNESS": "see_program_rows",
            "SCHLO_LOCALIZATION": "see_program_rows", "HUMAN_REVERSAL": "see_program_rows",
            "ORTHOGONAL_REVERSAL": "see_program_rows", "PROVENANCE": "see_program_rows",
            "main_figure_eligible": str(any_program).lower(),
            "release_terminal_state": "accepted_main" if any_program else "accepted_supplement",
            "reason": "at_least_one_frozen_object_passes_all_gates" if any_program else "no_frozen_object_passes_all_gates",
            "specification_sha256": seal["specification_sha256"],
        }
    )

    programs.to_csv(CANDIDATE_ROOT / "program_effects.tsv", sep="\t", index=False, na_rep="")
    testability.to_csv(CANDIDATE_ROOT / "program_testability.tsv", sep="\t", index=False, na_rep="")
    sensitivities.to_csv(CANDIDATE_ROOT / "sensitivity.tsv", sep="\t", index=False, na_rep="")
    classes.to_csv(CANDIDATE_ROOT / "evidence_class_effects.tsv", sep="\t", index=False, na_rep="")
    class_sensitivities.to_csv(CANDIDATE_ROOT / "evidence_class_sensitivity.tsv", sep="\t", index=False, na_rep="")
    support.to_csv(CANDIDATE_ROOT / "cross_assay_support.tsv", sep="\t", index=False, na_rep="")
    pd.DataFrame(verdict_rows).to_csv(CANDIDATE_ROOT / "figure5_verdict.tsv", sep="\t", index=False)
    source_manifest.to_csv(CANDIDATE_ROOT / "source_manifest.tsv", sep="\t", index=False, na_rep="")
    atomic_write_text(
        CANDIDATE_ROOT / "integration_status.json",
        json.dumps({"status": "complete", "main_figure_eligible": any_program, "release_terminal_state": "accepted_main" if any_program else "accepted_supplement", "n_program_effect_rows": len(programs), "n_class_effect_rows": len(classes)}, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps({"main_figure_eligible": any_program, "n_program_effect_rows": len(programs), "n_class_effect_rows": len(classes)}, indent=2))


if __name__ == "__main__":
    main()
