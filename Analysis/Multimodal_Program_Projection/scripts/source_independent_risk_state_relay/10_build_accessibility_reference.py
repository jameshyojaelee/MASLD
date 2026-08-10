#!/usr/bin/env python3
"""Seal the outcome-blind adult-liver chromatin-accessibility routing reference.

The release inventories existing hg38 peak sets and freezes how donor-collapsed
RNA lineages may map to them.  It deliberately does not inspect corrected
genetic loci, shared posteriors, or guideability and therefore cannot select an
experimental target.
"""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    sha256_file,
    write_tsv,
)


GSE244832_ROOT = PROJECT_ROOT / (
    "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2"
)
GSE281367_ROOT = PROJECT_ROOT / (
    "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/"
    "gse281367_hepatocyte_5fold_v2"
)

PRIMARY_PEAKS = {
    "Cholangiocytes": "Cholangiocytes_peaks.bed",
    "Circulating_NK_NKT": "Circulating_NK_NKT_peaks.bed",
    "Endothelial_cells": "Endothelial_cells_peaks.bed",
    "Fibroblasts": "Fibroblasts_peaks.bed",
    "Hepatocytes": "Hepatocytes_peaks.bed",
    "Macrophages": "Macrophages_peaks.bed",
    "Plasma_cells": "Plasma_cells_peaks.bed",
    "Resident_NK": "Resident_NK_peaks.bed",
    "T_cells": "T_cells_peaks.bed",
}

# Exact maps are eligible for a source-lineage gate.  Proxy and unavailable
# rows remain visible but cannot silently satisfy it.
LINEAGE_MAP = [
    ("B_cells", "Plasma_cells", "proxy", "B-cell peak set is not deposited"),
    ("Basophils", "", "unavailable", "no basophil/mast-cell ATAC reference"),
    ("Cholangiocytes", "Cholangiocytes", "exact", "matched label"),
    ("Circulating_NK_NKT", "Circulating_NK_NKT", "exact", "matched label"),
    ("Endothelial_cells", "Endothelial_cells", "exact", "matched label"),
    ("Fibroblasts", "Fibroblasts", "exact", "matched mesenchymal label"),
    ("Hepatocytes", "Hepatocytes", "exact", "matched label; independent hepatocyte replication required"),
    ("Macrophages", "Macrophages", "exact", "matched myeloid label"),
    ("Mono+mono_derived_cells", "Macrophages", "proxy", "monocyte-specific peak set is not deposited"),
    ("Neutrophils", "Macrophages", "proxy", "granulocyte peak set is not deposited"),
    ("Plasma_cells", "Plasma_cells", "exact", "matched label"),
    ("Resident_NK", "Resident_NK", "exact", "matched label"),
    ("T_cells", "T_cells", "exact", "matched label"),
    ("cDC1s", "Macrophages", "proxy", "dendritic-specific peak set is not deposited"),
    ("cDC2s", "Macrophages", "proxy", "dendritic-specific peak set is not deposited"),
    ("pDCs", "Macrophages", "proxy", "dendritic-specific peak set is not deposited"),
]


def audit_bed(path: Path, minimum_columns: int = 3) -> dict[str, int]:
    """Validate BED coordinates and return row/uniqueness counts."""
    row_count = 0
    unique_coordinates: set[tuple[str, int, int]] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < minimum_columns:
                raise RuntimeError(f"Malformed BED row {path}:{line_number}")
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise RuntimeError(f"Non-integer BED coordinate {path}:{line_number}") from exc
            if not chrom.startswith("chr") or start < 0 or end <= start:
                raise RuntimeError(f"Invalid BED interval {path}:{line_number}")
            row_count += 1
            unique_coordinates.add((chrom, start, end))
    if row_count == 0:
        raise RuntimeError(f"Empty BED source: {path}")
    return {
        "n_peak_rows": row_count,
        "n_unique_peak_coordinates": len(unique_coordinates),
        "n_duplicate_peak_rows": row_count - len(unique_coordinates),
    }


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite accessibility candidate: {CANDIDATE_ROOT}")

    pairing = PROJECT_ROOT / "data/GSE244832/metadata/donor_pairing.csv"
    label_summary = PROJECT_ROOT / (
        "Analysis/ATAC/Human_Multiome/results/label_transfer/label_transfer_summary.txt"
    )
    comparison = PROJECT_ROOT / (
        "Analysis/ATAC/Human_Multiome/results/label_transfer/peak_comparison.csv"
    )
    gse281367_peak = GSE281367_ROOT / "hepatocyte.idr_pooled_summit.narrowPeak"
    gse281367_contract = GSE281367_ROOT / "run_contract.json"
    gse281367_idr = GSE281367_ROOT / "peak_reproducibility_qc.json"
    gse281367_gate = GSE281367_ROOT / "gate_verdict.json"
    required = [
        pairing,
        label_summary,
        comparison,
        gse281367_peak,
        gse281367_contract,
        gse281367_idr,
        gse281367_gate,
        *(GSE244832_ROOT / name for name in PRIMARY_PEAKS.values()),
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing accessibility source: {path}")

    contract = json.loads(gse281367_contract.read_text(encoding="utf-8"))
    idr = json.loads(gse281367_idr.read_text(encoding="utf-8"))
    gate = json.loads(gse281367_gate.read_text(encoding="utf-8"))
    if contract.get("cohort") != "GSE281367" or contract.get("n_donors") != 12:
        raise RuntimeError("GSE281367 run contract drift")
    if idr.get("idr_threshold") != 0.05 or gate.get("comparison_ready") is not True:
        raise RuntimeError("GSE281367 reproducibility/comparison gate is not passed")
    if gate.get("saturation_authority") is not False:
        raise RuntimeError("GSE281367 reference unexpectedly claims saturation authority")

    manifest: list[dict[str, object]] = []
    for lineage, filename in sorted(PRIMARY_PEAKS.items()):
        path = GSE244832_ROOT / filename
        audit = audit_bed(path)
        manifest.append(
            {
                "source_id": f"GSE244832_{lineage}",
                "cohort": "GSE244832",
                "assay": "combinatorial_indexing_scATAC_MACS3",
                "reference_build": "GRCh38",
                "accessibility_lineage": lineage,
                "n_biological_donors": 18,
                "source_path": relative(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                **audit,
                "independent_replication_role": "primary_multilineage_reference",
                "known_boundary": (
                    "RNA-label-transfer overall concordance=0.484; binary accessibility "
                    "reference only, not a quantitative cell-type effect"
                ),
            }
        )

    audit = audit_bed(gse281367_peak, minimum_columns=10)
    manifest.append(
        {
            "source_id": "GSE281367_Hepatocytes_IDR",
            "cohort": "GSE281367",
            "assay": "standalone_snATAC_MACS3_IDR",
            "reference_build": "GRCh38",
            "accessibility_lineage": "Hepatocytes",
            "n_biological_donors": 12,
            "source_path": relative(gse281367_peak),
            "size_bytes": gse281367_peak.stat().st_size,
            "sha256": sha256_file(gse281367_peak),
            **audit,
            "independent_replication_role": "hepatocyte_replication_reference",
            "known_boundary": (
                "hepatocyte marker-argmax labels; RNA transfer unconfirmed; "
                "comparison-ready but saturation authority=false"
            ),
        }
    )

    map_rows = []
    for expression_lineage, accessibility_lineage, status, note in LINEAGE_MAP:
        map_rows.append(
            {
                "expression_lineage": expression_lineage,
                "accessibility_lineage": accessibility_lineage,
                "mapping_status": status,
                "eligible_for_source_lineage_gate": str(status == "exact").lower(),
                "hepatocyte_independent_replication_required": str(
                    expression_lineage == "Hepatocytes"
                ).lower(),
                "mapping_note": note,
            }
        )

    provenance = []
    for role, path in [
        ("GSE244832_authoritative_donor_pairing", pairing),
        ("GSE244832_label_transfer_quality", label_summary),
        ("GSE244832_peak_recall_audit", comparison),
        ("GSE281367_run_contract", gse281367_contract),
        ("GSE281367_IDR_reproducibility", gse281367_idr),
        ("GSE281367_comparison_gate", gse281367_gate),
    ]:
        provenance.append(
            {
                "role": role,
                "source_path": relative(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )

    CANDIDATE_ROOT.mkdir(parents=True)
    manifest_path = CANDIDATE_ROOT / "accessibility_source_manifest.tsv"
    map_path = CANDIDATE_ROOT / "lineage_accessibility_map.tsv"
    provenance_path = CANDIDATE_ROOT / "accessibility_provenance_manifest.tsv"
    gate_path = CANDIDATE_ROOT / "accessibility_gate_status.tsv"
    write_tsv(manifest_path, manifest, list(manifest[0]))
    write_tsv(map_path, map_rows, list(map_rows[0]))
    write_tsv(provenance_path, provenance, list(provenance[0]))
    write_tsv(
        gate_path,
        [
            {
                "gate": "reference_integrity",
                "status": "pass",
                "interpretation": "outcome_blind_accessibility_routing_reference",
                "experimental_targets_frozen": "false",
                "next_gate": (
                    "intersect pair-specific shared-95 variants after RSR-02, then "
                    "require exact lineage expression, accessibility, and guideability"
                ),
            }
        ],
        [
            "gate",
            "status",
            "interpretation",
            "experimental_targets_frozen",
            "next_gate",
        ],
    )
    seal = {
        "status": "accessibility_routing_reference_targets_not_selected",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "candidate_id": CANDIDATE_ROOT.name,
        "reference_build": "GRCh38",
        "n_peak_sources": len(manifest),
        "n_expression_lineages": len(map_rows),
        "n_exact_lineage_maps": sum(row["mapping_status"] == "exact" for row in map_rows),
        "gse244832_label_transfer_overall_concordance": 0.484,
        "hepatocyte_rule": (
            "shared-95 variant must overlap both GSE244832 and independent "
            "GSE281367 hepatocyte peak references"
        ),
        "nonhepatocyte_rule": (
            "exact GSE244832 lineage overlap is necessary but remains single-source"
        ),
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": False,
        "known_boundary": (
            "open chromatin is necessary but not sufficient for a functional element; "
            "binary peak overlap cannot establish causal cell of origin"
        ),
        "output_sha256": {
            manifest_path.name: sha256_file(manifest_path),
            map_path.name: sha256_file(map_path),
            provenance_path.name: sha256_file(provenance_path),
            gate_path.name: sha256_file(gate_path),
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "ACCESSIBILITY_REFERENCE_SEALED.json", seal)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
