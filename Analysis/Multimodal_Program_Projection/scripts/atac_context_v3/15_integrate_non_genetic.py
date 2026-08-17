#!/usr/bin/env python3
"""Candidate-only Figure 4F and Gene Catalog adapters after NON_GENETIC_READY."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
PROGRAM_RELEASE_ID = "program-context-v2-candidate-2026-08-07"
ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
SPATIAL_SOURCE = ROOT / (
    "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/"
    "figures/main/fig4_validation/data/fig4f_two_program_source_matrix.tsv"
)
MEMBERSHIP = ROOT / (
    "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/"
    "hotspot/program_membership_v2.tsv"
)


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    opener = __import__("gzip").open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, fields: tuple[str, ...] | list[str], rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_non_genetic_gate() -> dict[str, str]:
    gate = read(CANDIDATE / "NON_GENETIC_READY")
    if not gate:
        raise RuntimeError("NON_GENETIC_READY is empty")
    hashes = {}
    for row in gate:
        if row["release_id"] != RELEASE_ID or row["gate"] != "NON_GENETIC_READY":
            raise RuntimeError("NON_GENETIC_READY release or gate mismatch")
        if row["status"] != "READY":
            raise RuntimeError("NON_GENETIC_READY contains a non-ready artifact")
        artifact = (CANDIDATE / row["artifact"]).resolve()
        if CANDIDATE.resolve() not in artifact.parents or not artifact.is_file():
            raise RuntimeError(f"unsafe or missing sealed artifact: {artifact}")
        observed = sha256(artifact)
        if observed != row["sha256"]:
            raise RuntimeError(f"sealed artifact hash drift: {row['artifact']}")
        hashes[row["artifact"]] = observed
    return hashes


def main() -> None:
    validate_non_genetic_gate()
    out = CANDIDATE / "integration"
    if out.exists():
        raise RuntimeError(f"refusing to overwrite candidate integration: {out}")
    out.mkdir(parents=True)
    input_manifest = read(CANDIDATE / "input_manifest.tsv")
    frozen = {row["role"]: row["sha256"] for row in input_manifest}
    if sha256(SPATIAL_SOURCE) != frozen.get("spatial_fig4f_source"):
        raise RuntimeError("protected spatial Figure 4F source changed after input freeze")

    original_bytes = SPATIAL_SOURCE.read_bytes()
    original_copy = out / "fig4f_source_matrix_original.tsv"
    original_copy.write_bytes(original_bytes)
    if original_copy.read_bytes() != original_bytes:
        raise RuntimeError("protected Figure 4F copy is not byte-identical")
    with SPATIAL_SOURCE.open("r", encoding="utf-8", newline="") as handle:
        fields = next(csv.reader(handle, delimiter="\t"))

    program_rows = read(CANDIDATE / "programs/program_atac_results.tsv")
    coverage_rows = read(CANDIDATE / "programs/program_measurement_coverage.tsv")
    coverage_by_key = {
        (row["cohort"], row["program_uid"]): row for row in coverage_rows
    }
    if len(program_rows) != 234 or len(coverage_by_key) != 234:
        raise RuntimeError("ATAC program integration requires 234 unique cohort-program rows")
    appended = []
    state_rows = []
    for row in program_rows:
        coverage = coverage_by_key[(row["cohort"], row["program_uid"])]
        contrast = row["contrast_testable"] == "TRUE"
        measured = row["program_score_testable"] == "TRUE"
        effect = float(row["effect"]) if row["effect"] not in {"", "NA", "NaN"} else None
        equal = float(row["equal_weight_effect"]) if row["equal_weight_effect"] not in {"", "NA", "NaN"} else None
        leave = float(row["leave_top_gene_effect"]) if row["leave_top_gene_effect"] not in {"", "NA", "NaN"} else None
        payload = {
            "release_id": RELEASE_ID,
            "program_release_id": PROGRAM_RELEASE_ID,
            "registry_sha256": frozen["program_registry"],
            "dataset_id": row["cohort"],
            "assay_id": "snatac_promoter_accessibility_v3",
            "program_uid": row["program_uid"],
            "program_label": row["program_name"],
            "membership_sha256": row["membership_sha256"],
            "biological_unit": "donor",
            "biological_unit_resolution": "resolved",
            "n_biological": int(row["n_normal"]) + int(row["n_mash"]),
            "technical_unit": "donor_lineage_pseudobulk",
            "n_technical": int(row["n_normal"]) + int(row["n_mash"]),
            "source_dependence": (
                "source_dependent_combined_t_nk_label"
                if row["lineage_observability_state"] == "source_dependent_combined_t_nk_label"
                else "independent"
            ),
            "dataset_gate": "pass" if measured else "measurement_coverage_fail",
            "estimand": "MASH_minus_NORMAL" if contrast else "not_testable",
            "effect_unit": "donor_program_score_mean_difference",
            "estimate": row["effect"],
            "matched_null_sd": "",
            "pvalue": row["pvalue"],
            "qvalue": row["qvalue"],
            "multiplicity_family": "BH_complete_117_program_family_within_cohort",
            "n_null_draws": "",
            "n_genes_measured": row["n_measured_genes"],
            "retained_l1_weight": row["retained_l1_weight"],
            "equal_weight_sign_agree": str(
                effect is not None and equal is not None and effect * equal > 0
            ).upper(),
            "leave_top_weighted_gene_sign_agree": str(
                effect is not None and leave is not None and effect * leave > 0
            ).upper(),
            "within_source_result": row["within_source_state"],
            "evidence_state": row["cross_cohort_state"],
            "testability_reason": row["testability_reason"],
            "uncertainty_semantics": "limma_robust_empirical_bayes_standard_error",
            "source_release_id": RELEASE_ID,
            "source_row_sha256": "",
        }
        payload["source_row_sha256"] = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        appended.append(payload)
        state_rows.append({
            "release_id": RELEASE_ID,
            "cohort": row["cohort"],
            "program_uid": row["program_uid"],
            "membership_sha256": row["membership_sha256"],
            "lineage": row["lineage"],
            "n_program_genes": coverage["n_program_genes"],
            "n_promoter_measured_genes": row["n_measured_genes"],
            "retained_l1_weight": row["retained_l1_weight"],
            "promoter_measurement_state": (
                "observed" if int(row["n_measured_genes"]) > 0 else "unobserved"
            ),
            "program_score_state": "testable" if measured else "untestable",
            "contrast_state": "testable" if contrast else "untestable",
            "within_source_state": row["within_source_state"],
            "cross_cohort_state": row["cross_cohort_state"],
            "lineage_observability_state": row["lineage_observability_state"],
            "testability_reason": row["testability_reason"],
        })

    extended = out / "fig4f_source_matrix_with_atac_v3.tsv"
    extended.write_bytes(original_bytes)
    with extended.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writerows(appended)
    if extended.read_bytes()[: len(original_bytes)] != original_bytes:
        raise RuntimeError("existing spatial/protein/v2 rows were not preserved byte-for-byte")
    state_path = out / "fig4f_atac_v3_state_contract.tsv"
    write(state_path, tuple(state_rows[0]), state_rows)

    da = read(CANDIDATE / "da/da_peak_results.tsv.gz")
    promoter: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in da:
        for gene in filter(None, row["promoter_genes"].split(";")):
            key = (gene, row["lineage"], row["evidence_state"])
            if key not in promoter:
                promoter[key] = {
                    "release_id": RELEASE_ID,
                    "gene_symbol": gene,
                    "lineage": row["lineage"],
                    "promoter_da_state": row["evidence_state"],
                    "n_overlapping_tested_peaks": 1,
                    "source_dependence": str(
                        row["evidence_state"] == "source_dependent"
                    ).upper(),
                }
            else:
                promoter[key]["n_overlapping_tested_peaks"] = int(promoter[key]["n_overlapping_tested_peaks"]) + 1
    promoter_rows = sorted(
        promoter.values(),
        key=lambda row: (
            str(row["gene_symbol"]), str(row["lineage"]), str(row["promoter_da_state"])
        ),
    )
    write(
        out / "gene_catalog_promoter_da_adapter.tsv",
        (
            "release_id", "gene_symbol", "lineage", "promoter_da_state",
            "n_overlapping_tested_peaks", "source_dependence",
        ),
        promoter_rows,
    )

    membership = read(MEMBERSHIP)
    measured_symbols = [row for row in membership if row.get("mapped_symbol") and row.get("mapped_symbol_status") != "unmapped_or_ambiguous"]
    program_by_uid_cohort = {(row["program_uid"], row["cohort"]): row for row in program_rows}
    program_gene_rows = []
    for member in measured_symbols:
        for cohort in ("GSE244832", "GSE281367"):
            result = program_by_uid_cohort[(member["program_uid"], cohort)]
            program_gene_rows.append({
                "release_id": RELEASE_ID,
                "gene_symbol": member["mapped_symbol"],
                "program_uid": member["program_uid"],
                "cohort": cohort,
                "lineage": result["lineage"],
                "program_atac_promoter_coverage_state": (
                    "observed" if int(result["n_measured_genes"]) > 0 else "unobserved"
                ),
                "program_atac_score_state": (
                    "testable" if result["program_score_testable"] == "TRUE" else "untestable"
                ),
                "program_atac_contrast_state": result["within_source_state"],
                "program_atac_cross_cohort_state": result["cross_cohort_state"],
                "testability_reason": result["testability_reason"],
            })
    write(out / "gene_catalog_program_atac_adapter.tsv", tuple(program_gene_rows[0]), program_gene_rows)
    artifacts = (
        original_copy, extended, state_path, out / "gene_catalog_promoter_da_adapter.tsv",
        out / "gene_catalog_program_atac_adapter.tsv",
    )
    write(
        out / "integration_manifest.tsv",
        ("release_id", "role", "artifact", "sha256"),
        [
            {
                "release_id": RELEASE_ID,
                "role": "integration_artifact",
                "artifact": path.name,
                "sha256": sha256(path),
            }
            for path in artifacts
        ] + [
            {
                "release_id": RELEASE_ID,
                "role": "sealed_non_genetic_input_manifest",
                "artifact": "input_manifest.tsv",
                "sha256": sha256(CANDIDATE / "input_manifest.tsv"),
            },
            {
                "release_id": RELEASE_ID,
                "role": "post_gate_integration_producer",
                "artifact": str(Path(__file__).resolve().relative_to(ROOT)),
                "sha256": sha256(Path(__file__).resolve()),
            },
        ],
    )
    print(f"Appended {len(appended)} ATAC v3 program rows while preserving the protected source prefix")


if __name__ == "__main__":
    main()
