#!/usr/bin/env python3
"""Build the dependency-gated, explicitly provisional Plan 45 candidate audit.

This script does not freeze an experimental target list.  It records the exact
legacy candidates, current corrected-coloc availability, lead-variant allele
semantics, baseline lineage observability, and dated Cas13 coverage.  Final
selection is prohibited until the full provenance-corrected coloc portfolio is
complete and a new evidence/locus registry has been rebuilt from it.
"""

from __future__ import annotations

import csv
import datetime as dt
import math
from collections import defaultdict
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
    write_tsv,
)


PLAN43 = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "chronic-state-risk-bridge-2026-08-09"
)
LEGACY_LOCI = PLAN43 / "corrected_genetic_locus_registry.tsv"
LEGACY_CLASSES = PLAN43 / "integrity_correction_v2/corrected_evidence_classes_v2.tsv"
PHENOTYPES = PLAN43 / "frozen_inputs/phenotype_registry__phenotype_registry.tsv"
GWAS_REGISTRY = PROJECT_ROOT / "GWAS/finemapping/config/gwas_registry.tsv"
RERUN_ROOT = PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc_rerun"
EQTL_ROOT = PROJECT_ROOT / "data/broadaway_eqtl"
CELL_ORIGIN = PROJECT_ROOT / "RNA-seq/results/glp1ra/celltype_origin_axis.csv"
CAS13_FREEZE = (
    PROJECT_ROOT
    / "Cas13_Library_Design/scripts/guides/guide_overlap_fix_20260713"
    / "frozen_cas13_gene_library_20260731_v1/run_v1/library"
    / "frozen_cas13_gene_library_20260731_3504_targets.csv"
)

EXPECTED_RERUN_FILES = 50 * 22
DIRECT_STRATA = {
    "direct_masld_mash_diagnosis",
    "mri_pdff_or_histologic_steatosis",
}


def as_float(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) else float("nan")


def bool_text(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def read_csv_index(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row[key]: row
            for row in csv.DictReader(handle)
            if row.get(key, "").strip()
        }


def stream_selected_rows(
    path: Path, key_field: str, requested: set[str], delimiter: str = "\t"
) -> dict[str, dict[str, str]]:
    found: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter=delimiter):
            key = row.get(key_field, "")
            if key in requested:
                found[key] = row
                if len(found) == len(requested):
                    break
    return found


def exact_position_rows(
    path: Path,
    positions: set[int],
    position_field: str,
    delimiter: str = "\t",
) -> dict[int, dict[str, str]]:
    found: dict[int, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter=delimiter):
            try:
                position = int(row[position_field])
            except (KeyError, TypeError, ValueError):
                continue
            if position in positions:
                found[position] = row
                if len(found) == len(positions):
                    break
    return found


def parse_top_snp(value: str) -> tuple[str, int]:
    parts = str(value).replace("_", ":").split(":")
    if len(parts) < 2:
        raise ValueError(f"Unparseable SNP: {value}")
    return parts[0].removeprefix("chr"), int(parts[1])


def variant_orientation(
    gwas: dict[str, str] | None, eqtl: dict[str, str] | None
) -> dict[str, object]:
    blank = {
        "gwas_effect_allele": "",
        "gwas_other_allele": "",
        "gwas_beta": "",
        "gwas_p": "",
        "eqtl_effect_allele": "",
        "eqtl_other_allele": "",
        "eqtl_beta": "",
        "eqtl_p": "",
        "risk_allele": "",
        "risk_allele_effect_on_expression": "",
        "provisional_risk_expression_direction": "unresolved",
        "lead_orientation_status": "missing_source_variant",
    }
    if not gwas or not eqtl:
        return blank
    a1 = gwas.get("allele1", "").upper()
    a2 = gwas.get("allele2", "").upper()
    ea = eqtl.get("EA", "").upper()
    nea = eqtl.get("NEA", "").upper()
    g_beta = as_float(gwas.get("beta"))
    e_beta = as_float(eqtl.get("Beta"))
    if ea == a1 and nea == a2:
        e_on_a1 = e_beta
    elif ea == a2 and nea == a1:
        e_on_a1 = -e_beta
    else:
        blank["lead_orientation_status"] = "alleles_unresolved"
        return blank
    risk = a1 if g_beta > 0 else a2
    risk_effect = e_on_a1 if g_beta > 0 else -e_on_a1
    eqtl_p = as_float(eqtl.get("PVAL"))
    if not math.isfinite(eqtl_p) or eqtl_p >= 0.05:
        status = "lead_snp_eqtl_uninformative"
    else:
        status = "lead_snp_only_provisional"
    direction = "increase" if risk_effect > 0 else "decrease" if risk_effect < 0 else "zero"
    return {
        "gwas_effect_allele": a1,
        "gwas_other_allele": a2,
        "gwas_beta": gwas.get("beta", ""),
        "gwas_p": gwas.get("pval", ""),
        "eqtl_effect_allele": ea,
        "eqtl_other_allele": nea,
        "eqtl_beta": eqtl.get("Beta", ""),
        "eqtl_p": eqtl.get("PVAL", ""),
        "risk_allele": risk,
        "risk_allele_effect_on_expression": risk_effect,
        "provisional_risk_expression_direction": direction,
        "lead_orientation_status": status,
    }


def main() -> None:
    preliminary_path = CANDIDATE_ROOT / "PRELIMINARY.json"
    if preliminary_path.exists():
        raise RuntimeError(
            "Refusing to mutate a sealed Plan 45 snapshot. Set PLAN45_CANDIDATE_ID "
            "to a new dated identifier and rerun both preflight scripts there."
        )

    required = [
        LEGACY_LOCI,
        LEGACY_CLASSES,
        PHENOTYPES,
        GWAS_REGISTRY,
        CELL_ORIGIN,
        CAS13_FREEZE,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Plan 45 preflight inputs: " + "; ".join(missing))

    phenotypes = {row["study_name"]: row for row in read_tsv(PHENOTYPES)}
    classes = {row["gene_symbol"]: row for row in read_tsv(LEGACY_CLASSES)}
    registry = {row["study_name"]: row for row in read_tsv(GWAS_REGISTRY)}
    cell_origin = read_csv_index(CELL_ORIGIN, "gene")

    csv.field_size_limit(2**31 - 1)
    cas13 = read_csv_index(CAS13_FREEZE, "gene_symbol_human")

    provisional: list[dict[str, str]] = []
    for row in read_tsv(LEGACY_LOCI):
        study = row["tier12_driving_gwas"]
        meta = phenotypes.get(study, {})
        class_row = classes.get(row["gene_symbol"], {})
        if not bool_text(row["is_representative"]):
            continue
        if class_row.get("static_class") != "genetic_only":
            continue
        if meta.get("tier") != "1" or meta.get("phenotype_stratum") not in DIRECT_STRATA:
            continue
        provisional.append(row)
    provisional.sort(key=lambda row: row["gene_symbol"])
    if len(provisional) != 9:
        raise RuntimeError(f"Expected 9 legacy direct candidates, observed {len(provisional)}")

    by_sumstats: dict[Path, list[tuple[str, int]]] = defaultdict(list)
    by_eqtl: dict[Path, list[tuple[str, int]]] = defaultdict(list)
    parsed: dict[str, tuple[str, int]] = {}
    for row in provisional:
        gene = row["gene_symbol"]
        chrom, position = parse_top_snp(row["tier12_top_snp"])
        parsed[gene] = (chrom, position)
        study = row["tier12_driving_gwas"]
        sumstats = PROJECT_ROOT / "GWAS/finemapping" / registry[study]["sumstats_path"]
        by_sumstats[sumstats].append((gene, position))
        by_eqtl[EQTL_ROOT / f"chr{chrom}_marginal_summary_results.tsv"].append(
            (gene, position)
        )

    gwas_rows: dict[tuple[str, int], dict[str, str]] = {}
    for path, requests in by_sumstats.items():
        positions = {position for _, position in requests}
        rows = exact_position_rows(path, positions, "position")
        for gene, position in requests:
            if position in rows:
                gwas_rows[(gene, position)] = rows[position]

    eqtl_rows: dict[tuple[str, int], dict[str, str]] = {}
    for path, requests in by_eqtl.items():
        positions = {position for _, position in requests}
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    position = int(row["POS"])
                except (KeyError, TypeError, ValueError):
                    continue
                key = (row.get("GeneSymbol", ""), position)
                if key in {(gene, pos) for gene, pos in requests}:
                    eqtl_rows[key] = row

    expected_rerun_files = {
        RERUN_ROOT / study / f"susie_coloc_chr{chrom}.csv"
        for study in registry
        for chrom in range(1, 23)
    }
    rerun_files = set(RERUN_ROOT.glob("*/susie_coloc_chr*.csv"))
    missing_rerun_files = expected_rerun_files - rerun_files
    unexpected_rerun_files = rerun_files - expected_rerun_files
    empty_rerun_files = {path for path in rerun_files if path.stat().st_size == 0}
    valid_rerun_files = (rerun_files & expected_rerun_files) - empty_rerun_files
    outputs: list[dict[str, object]] = []
    for row in provisional:
        gene = row["gene_symbol"]
        study = row["tier12_driving_gwas"]
        chrom, position = parsed[gene]
        rerun_file = RERUN_ROOT / study / f"susie_coloc_chr{chrom}.csv"
        rerun_row: dict[str, str] = {}
        if rerun_file.is_file():
            rerun_row = stream_selected_rows(rerun_file, "gene", {gene}, delimiter=",").get(
                gene, {}
            )
        orientation = variant_orientation(
            gwas_rows.get((gene, position)), eqtl_rows.get((gene, position))
        )
        class_row = classes[gene]
        origin = cell_origin.get(gene, {})
        cas = cas13.get(gene, {})
        outputs.append(
            {
                "gene_symbol": gene,
                "ensembl_id": row["ensembl_id"],
                "legacy_study": study,
                "trait": row["tier12_driving_trait"],
                "phenotype_stratum": phenotypes[study]["phenotype_stratum"],
                "ancestry": phenotypes[study]["ancestry"],
                "legacy_locus_uid": row["tier12_locus_uid"],
                "legacy_top_snp": row["tier12_top_snp"],
                "legacy_susie_pp4": row["tier12_susie_pp4"],
                "legacy_abf_pp4": row["tier12_abf_pp4"],
                "corrected_rerun_file_status": (
                    "candidate_row_present"
                    if rerun_row
                    else "file_present_gene_absent"
                    if rerun_file.is_file()
                    else "pending"
                ),
                "corrected_rerun_susie_pp4": rerun_row.get("PP.H4.susie", ""),
                "corrected_rerun_abf_pp4": rerun_row.get("PP.H4.abf", ""),
                "corrected_rerun_method": rerun_row.get("method", ""),
                **orientation,
                "bulk_AveExpr": class_row.get("bulk_AveExpr", ""),
                "bulk_logFC": class_row.get("bulk_logFC", ""),
                "bulk_treat_fdr": class_row.get("bulk_treat_fdr", ""),
                "baseline_dominant_celltype": origin.get("dominant_celltype", ""),
                "baseline_hepatocyte_fraction": origin.get("hep_fraction", ""),
                "baseline_hepatocyte_cpm": origin.get("hep_cpm", ""),
                "baseline_expression_source_status": "legacy_donor_key_requires_refreeze",
                "dated_cas13_freeze_member": str(bool(cas)).lower(),
                "dated_cas13_mouse_gene": cas.get("gene_symbol_mouse", ""),
                "target_freeze_status": "prohibited_until_corrected_coloc_registry",
            }
        )

    fields = list(outputs[0])
    candidate_path = CANDIDATE_ROOT / "provisional_target_audit.tsv"
    write_tsv(candidate_path, outputs, fields)

    complete = (
        len(valid_rerun_files) == EXPECTED_RERUN_FILES
        and not missing_rerun_files
        and not unexpected_rerun_files
        and not empty_rerun_files
    )
    gate_rows = [
        {
            "gate_id": "corrected_coloc_portfolio",
            "status": "pass" if complete else "waiting",
            "observed": len(valid_rerun_files),
            "required": EXPECTED_RERUN_FILES,
            "detail": (
                "One exact, nonempty per-study/per-chromosome output is required "
                f"before registry rebuild; missing={len(missing_rerun_files)}, "
                f"unexpected={len(unexpected_rerun_files)}, empty={len(empty_rerun_files)}."
            ),
        },
        {
            "gate_id": "corrected_evidence_registry_v3",
            "status": "not_started" if not complete else "required_next",
            "observed": 0,
            "required": 1,
            "detail": "Must be rebuilt from corrected rerun; Plan 43's 326-locus registry is not final for experiment design.",
        },
        {
            "gate_id": "risk_orientation",
            "status": "provisional_only",
            "observed": sum(
                row["lead_orientation_status"] == "lead_snp_only_provisional"
                for row in outputs
            ),
            "required": "credible-set consensus for every frozen target",
            "detail": "Lead-SNP signs are diagnostic only and cannot freeze perturbation direction.",
        },
        {
            "gate_id": "target_selection",
            "status": "prohibited",
            "observed": 0,
            "required": "all upstream gates pass",
            "detail": "No target has been frozen or exposed to experimental outcomes.",
        },
    ]
    gate_path = CANDIDATE_ROOT / "source_gate_status.tsv"
    write_tsv(
        gate_path,
        gate_rows,
        ["gate_id", "status", "observed", "required", "detail"],
    )

    input_manifest = []
    for path in required:
        input_manifest.append(
            {
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "size_bytes": path.stat().st_size,
                "mtime_utc": dt.datetime.fromtimestamp(
                    path.stat().st_mtime, tz=dt.timezone.utc
                ).isoformat(),
                "sha256": sha256_file(path) if path.stat().st_size < 25_000_000 else "deferred_large_source",
            }
        )
    manifest_path = CANDIDATE_ROOT / "preflight_input_manifest.tsv"
    write_tsv(
        manifest_path,
        input_manifest,
        ["source_path", "size_bytes", "mtime_utc", "sha256"],
    )

    payload = {
        "status": "dependency_gated_preliminary",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "legacy_candidate_count": len(outputs),
        "corrected_rerun_outputs": len(valid_rerun_files),
        "corrected_rerun_expected": EXPECTED_RERUN_FILES,
        "target_list_frozen": False,
        "scientific_outcomes_inspected": False,
        "provisional_target_audit_sha256": sha256_file(candidate_path),
        "source_gate_status_sha256": sha256_file(gate_path),
        "preflight_input_manifest_sha256": sha256_file(manifest_path),
    }
    atomic_write_json(preliminary_path, payload)
    print(
        f"Plan 45 preflight: {len(outputs)} legacy candidates; "
        f"corrected coloc {len(valid_rerun_files)}/{EXPECTED_RERUN_FILES}; no targets frozen"
    )


if __name__ == "__main__":
    main()
