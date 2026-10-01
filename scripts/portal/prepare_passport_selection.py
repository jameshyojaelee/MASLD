#!/usr/bin/env python3
"""Write the exact, unsigned production PASS input-selection table.

The row roster and parent/terminal links are frozen here so the coordinator
does not have to hand-edit 28-column TSV rows.  This helper computes read-only
artifact/producer hashes but deliberately does not create the sibling
attestation; the coordinating agent must inspect the table and run
``prepare_passport_attestations.py selection`` separately.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import generate_evidence_passports as contract


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = contract.PROJECT_ROOT
MANUSCRIPT_RELEASE = "program-context-v2-candidate-2026-08-07"
PASSPORT_RELEASE = contract.PRODUCTION_ANALYSIS_RELEASE_ID
MM_ROOT = (
    PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates" / MANUSCRIPT_RELEASE
)
GEN_ROOT = MM_ROOT / "genetics_context"
# PASSPORT_GEN_RAW_ROOT points GEN_RAW at frozen classes rebuilt on new bulk and
# COLOC inputs (01_freeze_and_rederive.py under GEN_CONTEXT_ROOT_REL). The GEN
# terminal closure, READY, phenotype and provenance rows stay on GEN_ROOT: they
# carry no bulk or COLOC value.
GEN_RAW_ROOT = Path(os.environ.get("PASSPORT_GEN_RAW_ROOT", GEN_ROOT))
GEN_RAW_RELEASE = os.environ.get("PASSPORT_GEN_RAW_RELEASE", f"{MANUSCRIPT_RELEASE}:GEN")
STALE_WITH_NEW_GEN_RAW = {
    "GEN_CLOSURE": (
        "STALE: allowed_downstream_use 'all_joint_34_of_447_descriptive_interface' "
        "counts the v1 frozen classes, not the rebuilt GEN_RAW."
    ),
    "MYOJIN_CLASS": (
        "STALE: the evidence-class test used the v1 frozen classes "
        "(228 genetic_only, 463 disease_state_only, 15 convergent); not rerun."
    ),
}
HOTSPOT_ROOT = MM_ROOT / "hotspot"
PLAN13_ROOT = MM_ROOT / "spatial_context_semantic_v2_2026-08-08"
MYOJIN_ROOT = MM_ROOT / "myojin_hlf"
GEN_ID_ROOT = contract.DEFAULT_CANDIDATE_ROOT / "preflight/gen_identity_adjudication_v1"
DEFAULT_OUTPUT = contract.DEFAULT_CANDIDATE_ROOT / "passport_input_selection.tsv"


def _row(
    input_id: str,
    workstream_id: str,
    artifact_role: str,
    adapter_id: str,
    artifact_grain: str,
    artifact_path: Path,
    producer_script: Path,
    source_release_id: str,
    gate_verdict: str,
    terminal_gate_input_id: str,
    allowed_claim_wording: str,
    limitation: str,
    source_datasets: str,
    source_cohorts: str,
    source_publication: str,
    source_url: str,
    biological_unit: str,
    provenance_state: str,
    claim_scope: str,
    inclusion_destination: str,
    parent_input_ids: str,
) -> dict[str, Any]:
    artifact_path = artifact_path.resolve()
    producer_script = producer_script.resolve()
    if not artifact_path.is_file():
        raise contract.PassportContractError(
            "SELECTION_ARTIFACT_MISSING", f"{input_id}:{artifact_path}"
        )
    if not producer_script.is_file():
        raise contract.PassportContractError(
            "SELECTION_PRODUCER_MISSING", f"{input_id}:{producer_script}"
        )
    return {
        "input_id": input_id,
        "workstream_id": workstream_id,
        "artifact_role": artifact_role,
        "adapter_id": adapter_id,
        "artifact_grain": artifact_grain,
        "artifact_path": str(artifact_path),
        "artifact_sha256": contract.sha256_file(artifact_path),
        "artifact_bytes": artifact_path.stat().st_size,
        "producer_script": str(producer_script),
        "producer_sha256": contract.sha256_file(producer_script),
        "source_release_id": source_release_id,
        "gate_verdict": gate_verdict,
        "terminal_gate_input_id": terminal_gate_input_id,
        "allowed_claim_wording": allowed_claim_wording,
        "limitation": limitation,
        "source_datasets": source_datasets,
        "source_cohorts": source_cohorts,
        "source_publication": source_publication,
        "source_url": source_url,
        "biological_unit": biological_unit,
        "n_biological_units": "",
        "n_technical_units": "",
        "provenance_state": provenance_state,
        "claim_scope": claim_scope,
        "inclusion_destination": inclusion_destination,
        "selection_decision": "accepted",
        "parent_input_ids": parent_input_ids,
        "fixture_only": "false",
    }


def production_rows() -> list[dict[str, Any]]:
    gen_scripts = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2"
    hotspot_scripts = PROJECT_ROOT / "Analysis/SingleCell/scripts/hotspot_modules"
    plan13_scripts = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
    myojin_scripts = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/scripts/myojin_hlf"
    adjudicator = PROJECT_ROOT / "scripts/portal/adjudicate_gen_ensembl_identity.py"
    rows = [
        _row(
            "GENCODE_IDENTITY", "PASS01", "gencode_v49_identity",
            "gencode_v49_identity_v1", "reference",
            PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz",
            PROJECT_ROOT / "scripts/build_gencode_metadata.R", "GENCODE_v49_GRCh38p14",
            "include", "", "Stable GENCODE-v49 gene identity reference.",
            "Identity reference only; symbols are display aliases.", "GENCODE_v49", "",
            "GENCODE release 49", "https://www.gencodegenes.org/human/release_49.html",
            "annotated gene", "independent", "identity_only", "passport_only", "",
        ),
        _row(
            "GEN_CLOSURE", "PLAN30", "gen_terminal_closure", "gen_terminal_closure_v1",
            "gate", GEN_ROOT / "terminal_closure.tsv", gen_scripts / "09_close_genetics_terminal_gate.py",
            f"{MANUSCRIPT_RELEASE}:GEN", "include", "",
            "GEN coverage-limited terminal; positive evidence only.",
            "Context rescue and powered-negative genetic claims are prohibited.",
            "GEN_frozen_interface", "", "GEN terminal closure", "", "gate",
            "independent", "gate_only", "passport_only", "",
        ),
        _row(
            "GEN_READY", "PLAN30", "gen_terminal_ready", "gen_terminal_ready_v1", "gate",
            GEN_ROOT / "GEN_TERMINAL_CLOSURE_READY",
            gen_scripts / "10_validate_genetics_terminal_closure.py",
            f"{MANUSCRIPT_RELEASE}:GEN", "include", "", "Validated GEN terminal closure.",
            "Not a gene-level evidence result.", "GEN_frozen_interface", "",
            "GEN terminal closure", "", "gate", "independent", "gate_only",
            "passport_only", "GEN_CLOSURE",
        ),
        _row(
            "GEN_PHENOTYPES", "PLAN30", "gen_phenotype_registry", "source_provenance_v1",
            "provenance", GEN_ROOT / "phenotype_registry.tsv",
            gen_scripts / "02_build_phenotype_registry.py", f"{MANUSCRIPT_RELEASE}:GEN",
            "include", "GEN_READY", "Source-defined phenotype and trait provenance.",
            "Provenance only; creates no gene call.", "GEN_GWAS_registry", "source-specific",
            "GEN phenotype registry", "", "source-defined study", "reused_source",
            "source_provenance_only", "passport_only", "GEN_READY",
        ),
        _row(
            "GEN_PROVENANCE", "PLAN30", "gen_source_provenance", "source_provenance_v1",
            "provenance", GEN_ROOT / "source_provenance.tsv",
            gen_scripts / "02_build_phenotype_registry.py", f"{MANUSCRIPT_RELEASE}:GEN",
            "include", "GEN_READY", "Source and assay provenance for the frozen GEN interface.",
            "Provenance only; source-specific denominators must not be collapsed.",
            "canonical_human_bulk;GWAS_finemapping;Broadaway_liver_eQTL", "source-specific",
            "GEN source provenance", "", "source-defined study", "reused_source",
            "source_provenance_only", "passport_only", "GEN_READY",
        ),
        _row(
            "GEN_RAW", "PLAN30", "gen_frozen_classes_raw_identity_source",
            "gen_raw_identity_source_v1", "provenance", GEN_RAW_ROOT / "frozen_evidence_classes.tsv",
            gen_scripts / "01_freeze_and_rederive.py", GEN_RAW_RELEASE,
            "include", "GEN_READY", "Frozen source calls retained verbatim for identity adjudication.",
            "Mixed or duplicate Ensembl fields require outcome-blind v49 adjudication.",
            "canonical_human_bulk;GWAS_finemapping;Broadaway_liver_eQTL", "source-specific",
            "GEN frozen evidence classes", "", "gene/source-defined study", "reused_source",
            "source_provenance_only", "passport_only",
            "GEN_READY;GEN_PHENOTYPES;GEN_PROVENANCE;GENCODE_IDENTITY",
        ),
        _row(
            "GEN_ID_READY", "PASS02A", "gen_identity_adjudication_ready",
            "gen_identity_adjudication_ready_v1", "gate",
            GEN_ID_ROOT / "GEN_IDENTITY_ADJUDICATION_READY", adjudicator,
            PASSPORT_RELEASE, "include", "GEN_READY",
            "Outcome-blind GENCODE-v49 identity adjudication completed.",
            "Candidate-only terminal; canonical promotion is false.", "GEN_frozen_interface",
            "", "PASS GEN identity adjudication", "", "source row", "reused_source",
            "gate_only", "passport_only", "GEN_READY;GEN_RAW;GENCODE_IDENTITY",
        ),
        _row(
            "GEN_ID_AUDIT", "PASS02A", "gen_identity_adjudication_audit",
            "gen_identity_adjudication_audit_v1", "provenance",
            GEN_ID_ROOT / "gen_identity_adjudication_audit.tsv", adjudicator,
            PASSPORT_RELEASE, "include", "GEN_ID_READY",
            "Identity resolution and exclusion counts.",
            "Audit only; creates no evidence call.", "GEN_frozen_interface", "",
            "PASS GEN identity adjudication", "", "source row", "reused_source",
            "source_provenance_only", "passport_only", "GEN_ID_READY",
        ),
        _row(
            "GEN_ID_QUARANTINE", "PASS02A", "gen_identity_adjudication_quarantine",
            "gen_identity_adjudication_quarantine_v1", "provenance",
            GEN_ID_ROOT / "gen_identity_quarantine.tsv", adjudicator, PASSPORT_RELEASE,
            "include", "GEN_ID_READY",
            "Source rows excluded for identity ambiguity, absent-v49 ID, or call conflict.",
            "Excluded rows never create Gene Catalog entries or negative calls.", "GEN_frozen_interface",
            "", "PASS GEN identity adjudication", "", "source row", "reused_source",
            "source_provenance_only", "passport_only", "GEN_ID_READY",
        ),
        _row(
            "GEN_CLASSES", "PASS03", "gen_frozen_classes_ensembl_adjudicated",
            "gen_adjudicated_classes_v1", "gene",
            GEN_ID_ROOT / "gen_frozen_classes_ensembl_adjudicated.tsv", adjudicator,
            PASSPORT_RELEASE, "include", "GEN_ID_READY",
            "Genetically anchored and established-state-associated evidence in their frozen roles.",
            "No context rescue, causal/reactive binary, or powered genetic null; genetic denominators remain source-specific.",
            "canonical_human_bulk;GWAS_finemapping;Broadaway_liver_eQTL", "source-specific",
            "GEN frozen interface", "", "gene/source-defined study", "reused_source",
            "gene_evidence", "passport_only", "GEN_ID_READY;GEN_PHENOTYPES;GEN_PROVENANCE",
        ),
        _row(
            "HOTSPOT_READY", "PLAN20", "hotspot_semantic_ready", "hotspot_semantic_ready_v1",
            "gate", HOTSPOT_ROOT / "SEMANTIC_ADJUDICATION_READY",
            hotspot_scripts / "516_validate_hotspot_v2_semantics.py",
            f"{MANUSCRIPT_RELEASE}:Hotspot-v2", "include", "",
            "Validated Hotspot-v2 semantic adjudication.", "No tested-negative program is authorized.",
            "integrated_scRNA_atlas", "integrated_scRNA_donors", "Hotspot v2 candidate", "",
            "program", "reused_source", "gate_only", "passport_only", "",
        ),
        _row(
            "HOTSPOT_SEMANTICS", "PLAN20", "hotspot_program_semantics",
            "hotspot_program_semantics_v1", "program",
            HOTSPOT_ROOT / "program_registry_v2_semantic_adjudication.tsv",
            hotspot_scripts / "515_adjudicate_hotspot_v2_semantics.py",
            f"{MANUSCRIPT_RELEASE}:Hotspot-v2", "include", "HOTSPOT_READY",
            "Frozen program definitions and internal stage-association semantics.",
            "Program grain only; nonsignificance is indeterminate and membership cannot expand to gene calls.",
            "integrated_scRNA_atlas", "integrated_scRNA_donors", "Hotspot v2 candidate", "",
            "program", "reused_source", "program_context", "passport_only", "HOTSPOT_READY",
        ),
        _row(
            "HOTSPOT_MEMBERSHIP", "PLAN20", "hotspot_program_membership",
            "hotspot_program_membership_v1", "program_membership",
            HOTSPOT_ROOT / "program_membership_v2.tsv",
            hotspot_scripts / "513_freeze_hotspot_v2_registry.R",
            f"{MANUSCRIPT_RELEASE}:Hotspot-v2", "include", "HOTSPOT_READY",
            "Frozen program membership and weights.",
            "Membership only; unmapped source members remain in audit and no member-gene call is authorized.",
            "integrated_scRNA_atlas", "integrated_scRNA_donors", "Hotspot v2 candidate", "",
            "program member", "reused_source", "program_membership_only", "passport_only",
            "HOTSPOT_READY",
        ),
        _row(
            "PLAN13_FINAL_READY", "PLAN13", "plan13_final_ready", "plan13_final_ready_v1",
            "gate", PLAN13_ROOT / "final_integration/READY",
            plan13_scripts / "14_validate_final_integration.py",
            f"{MANUSCRIPT_RELEASE}:Plan13-semantic-v2", "include", "",
            "Complete real Plan13 integration terminal.",
            "Candidate source tables only; no PDF or canonical promotion.", "Plan13_integrated_context",
            "source-specific", "Plan13 final integration", "", "source-defined donor",
            "partially_dependent", "gate_only", "passport_only",
            "HOTSPOT_SEMANTICS;HOTSPOT_MEMBERSHIP",
        ),
        _row(
            "PLAN13_SEMANTIC_READY", "PLAN13", "plan13_semantic_v2_ready",
            "plan13_terminal_ready_v1", "gate", PLAN13_ROOT / "SEMANTIC_V2_READY",
            plan13_scripts / "17_seal_semantic_v2.py",
            f"{MANUSCRIPT_RELEASE}:Plan13-semantic-v2", "include", "",
            "Sealed semantic-v2 Plan13 terminal with zero tested-negative calls.",
            "Historical semantic labels are ineligible; canonical promotion is false.",
            "Plan13_integrated_context", "source-specific", "Plan13 semantic-v2 seal", "",
            "source-defined donor", "partially_dependent", "gate_only", "passport_only",
            "PLAN13_FINAL_READY",
        ),
        _row(
            "PLAN13_EFFECTS", "PLAN13", "plan13_integrated_program_context",
            "plan13_program_context_v1", "program",
            PLAN13_ROOT / "final_integration/integrated_program_effects.tsv",
            plan13_scripts / "13_build_final_integration.py",
            f"{MANUSCRIPT_RELEASE}:Plan13-semantic-v2", "include", "PLAN13_FINAL_READY",
            "Assay-native program context with source-specific biological units.",
            "Different effect units are not comparable; no member-gene expansion.",
            "GSE192741;GSE244832;GSE281367;GSE287826;Govaere2026_CosMx;Govaere2026_GeoMx;PXD051911;Vu_et_al_2025;Yakubovsky_2026",
            "source-specific", "Plan13 integrated source manifest", "", "source-defined donor",
            "partially_dependent", "program_context", "passport_only",
            "PLAN13_FINAL_READY;HOTSPOT_SEMANTICS;HOTSPOT_MEMBERSHIP",
        ),
        _row(
            "PLAN13_DATASETS", "PLAN13", "plan13_dataset_status", "dataset_status_v1",
            "dataset", PLAN13_ROOT / "final_integration/figure4_dataset_verdict.tsv",
            plan13_scripts / "13_build_final_integration.py",
            f"{MANUSCRIPT_RELEASE}:Plan13-semantic-v2", "include", "PLAN13_FINAL_READY",
            "Dataset-level include, nonconfirmatory, or skipped-source-gate verdict.",
            "GSE287826 skip is dataset-level and creates no gene pseudo-results.",
            "GSE192741;GSE244832;GSE281367;GSE287826;Govaere2026_CosMx;Govaere2026_GeoMx;PXD051911;Vu_et_al_2025;Yakubovsky_2026",
            "source-specific", "Plan13 integrated source manifest", "", "source-defined donor",
            "partially_dependent", "dataset_status_only", "passport_only", "PLAN13_FINAL_READY",
        ),
        _row(
            "MYOJIN_READY", "PLAN40", "myojin_phase_c_ready", "myojin_terminal_ready_v1",
            "gate", MYOJIN_ROOT / "PHASE_C_VALIDATED",
            myojin_scripts / "11_validate_phase_c.py", f"{MANUSCRIPT_RELEASE}:Myojin-HLF",
            "include", "", "Validated outcome-blind Myojin Phase C terminal.",
            "Biological positivity was not required for software acceptance.",
            "Myojin_HLF_palmitate;DepMap_24Q4", "HLF_cell_line", "Myojin HLF screen", "",
            "source-defined gene or program", "independent", "gate_only", "passport_only",
            "GEN_CLASSES;HOTSPOT_SEMANTICS;HOTSPOT_MEMBERSHIP",
        ),
        _row(
            "MYOJIN_ASSAY", "PLAN40", "myojin_figure5_verdict", "myojin_assay_verdict_v1",
            "assay", MYOJIN_ROOT / "fig5_verdict.tsv", myojin_scripts / "10_execute_one_pass.py",
            f"{MANUSCRIPT_RELEASE}:Myojin-HLF", "complete_nonconfirmatory", "MYOJIN_READY",
            "The prespecified Myojin assay did not qualify for a main figure.",
            "Assay-grain non-support only; no gene-level negative.",
            "Myojin_HLF_palmitate;DepMap_24Q4", "HLF_cell_line", "Myojin HLF screen", "",
            "assay", "independent", "assay_status_only", "supplement", "MYOJIN_READY",
        ),
        _row(
            "MYOJIN_CLASS", "PLAN40", "myojin_class_non_support",
            "myojin_class_nonsupport_v1", "assay", MYOJIN_ROOT / "class_effects.tsv",
            myojin_scripts / "10_execute_one_pass.py", f"{MANUSCRIPT_RELEASE}:Myojin-HLF",
            "complete_nonconfirmatory", "MYOJIN_READY",
            "Complete prespecified evidence-class test retained as nonconfirmatory.",
            "Evidence-class grain only; no tested-negative gene call.",
            "Myojin_HLF_palmitate;DepMap_24Q4", "HLF_cell_line", "Myojin HLF screen", "",
            "tested gene", "independent", "assay_status_only", "supplement", "MYOJIN_READY",
        ),
        _row(
            "MYOJIN_PROGRAM", "PLAN40", "myojin_program_non_support",
            "myojin_program_nonsupport_v1", "assay", MYOJIN_ROOT / "program_effects.tsv",
            myojin_scripts / "10_execute_one_pass.py", f"{MANUSCRIPT_RELEASE}:Myojin-HLF",
            "complete_nonconfirmatory", "MYOJIN_READY",
            "Complete prespecified program test retained as nonconfirmatory.",
            "Program grain only; no member-gene or tested-negative call.",
            "Myojin_HLF_palmitate;DepMap_24Q4", "HLF_cell_line", "Myojin HLF screen", "",
            "program", "independent", "assay_status_only", "supplement", "MYOJIN_READY",
        ),
    ]
    if len({row["input_id"] for row in rows}) != len(rows):
        raise contract.PassportContractError("SELECTION_INPUT_ID", "duplicate frozen row ID")
    if GEN_RAW_ROOT.resolve() != GEN_ROOT.resolve():
        # Reused rows whose content was computed on the v1 frozen classes.
        for row in rows:
            if row["input_id"] in STALE_WITH_NEW_GEN_RAW:
                row["limitation"] = f"{row['limitation']} {STALE_WITH_NEW_GEN_RAW[row['input_id']]}"
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    if output.exists():
        raise contract.PassportContractError(
            "SELECTION_EXISTS", f"refusing to overwrite existing selection: {output}"
        )
    if output != DEFAULT_OUTPUT.resolve():
        raise contract.PassportContractError(
            "SELECTION_OUTPUT_ROOT", f"production selection must be {DEFAULT_OUTPUT.resolve()}"
        )
    rows = production_rows()
    contract.atomic_write_tsv(output, rows, contract.SELECTION_REQUIRED_COLUMNS)
    print(
        f"PASS: unsigned production selection rows={len(rows)};path={output};"
        f"sha256={contract.sha256_file(output)};attestation_required=true"
    )


if __name__ == "__main__":
    try:
        main()
    except contract.PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
