#!/usr/bin/env python3
"""Record why the cross-assay coding-versus-lncRNA comparison is not yet valid."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    require(path.is_file() and not path.is_symlink(), f"missing input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(all(None not in row for row in rows), f"malformed row width: {path}")
    return rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--biotype-counts", type=Path, required=True)
    parser.add_argument("--biotype-review", type=Path, required=True)
    parser.add_argument("--eqtl-readiness", type=Path, required=True)
    parser.add_argument("--spatial-axis-audit", type=Path, required=True)
    parser.add_argument("--singlecell-gene-universe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists")

    counts = read_tsv(args.biotype_counts)
    require(
        {row["display_biotype"] for row in counts}
        == {"protein_coding", "lncRNA", "other"},
        "bulk biotype family drift",
    )
    review = json.loads(args.biotype_review.read_text(encoding="utf-8"))
    require(review.get("status") == "pass", "bulk biotype review did not pass")

    eqtl = {row["requirement"]: row for row in read_tsv(args.eqtl_readiness)}
    require(
        eqtl["complete_source_tested_gene_background"]["status"] == "missing",
        "eQTL tested-background status changed; reopen the readiness decision",
    )
    spatial = read_tsv(args.spatial_axis_audit)
    require(len(spatial) == 7, "accepted spatial assay family drift")
    sc_rows = read_tsv(args.singlecell_gene_universe)
    sc_ids = {row["gene_name"] for row in sc_rows}
    require(
        any(value.startswith("ENSG") for value in sc_ids)
        and any(not value.startswith("ENSG") for value in sc_ids),
        "single-cell identifier-mix premise changed",
    )

    audit = [
        {
            "assay": "five_cohort_bulk_rna",
            "native_scope": "versioned_GENCODE_v49_tested_gene_family",
            "readiness": "ready",
            "lncrna_applicability": "applicable",
            "reason": "complete versioned all-gene family and matching covariates are available",
        },
        {
            "assay": "deposited_liver_eqtl",
            "native_scope": "source_significant_eGenes_only",
            "readiness": "blocked",
            "lncrna_applicability": "applicable",
            "reason": "complete source-tested eGene background is not deposited; absence cannot be called unobservable",
        },
        {
            "assay": "donor_collapsed_single_cell_rna",
            "native_scope": "mixed_symbol_and_Ensembl_gene_axis",
            "readiness": "blocked",
            "lncrna_applicability": "applicable",
            "reason": "a frozen one-to-one versioned GENCODE adapter has not been produced for the full gene axis",
        },
        {
            "assay": "accepted_spatial_rna",
            "native_scope": "seven_assay_gene_axes",
            "readiness": "blocked",
            "lncrna_applicability": "applicable",
            "reason": "accepted integrated gene-context rows cover the frozen program union, not the complete versioned GENCODE universe",
        },
        {
            "assay": "promoter_accessible_chromatin",
            "native_scope": "locus_level_accessibility",
            "readiness": "blocked",
            "lncrna_applicability": "locus_context_only",
            "reason": "chromatin observes loci; a complete versioned promoter-locus capability table is not frozen",
        },
        {
            "assay": "liver_proteomics",
            "native_scope": "detected_protein_products",
            "readiness": "partial",
            "lncrna_applicability": "not_applicable",
            "reason": "protein coverage can describe coding genes, but cannot test lncRNA molecules and is not a matched biotype comparison",
        },
    ]

    raw = []
    for row in counts:
        raw.append(
            {
                "assay": "five_cohort_bulk_rna",
                "display_biotype": row["display_biotype"],
                "n_observable": int(row["n_tested"]),
                "denominator": int(row["n_tested"]),
                "raw_fraction": 1.0,
                "state": "observable_in_tested_bulk_family",
                "qualification": "the denominator is the assay-native filtered family, not all GENCODE genes",
            }
        )
    raw.append(
        {
            "assay": "liver_proteomics",
            "display_biotype": "lncRNA",
            "n_observable": "NA",
            "denominator": "NA",
            "raw_fraction": "NA",
            "state": "not_applicable",
            "qualification": "proteomics does not observe mature lncRNA molecules",
        }
    )

    args.output.mkdir(parents=True)
    write_tsv(
        args.output / "assay_source_audit.tsv",
        audit,
        ["assay", "native_scope", "readiness", "lncrna_applicability", "reason"],
    )
    write_tsv(
        args.output / "raw_observability.tsv",
        raw,
        [
            "assay",
            "display_biotype",
            "n_observable",
            "denominator",
            "raw_fraction",
            "state",
            "qualification",
        ],
    )
    write_tsv(
        args.output / "matching_gate_status.tsv",
        [
            {
                "gate": "abundance_adjusted_cross_assay_biotype_matching",
                "status": "skipped_incomplete_assay_matrix",
                "ready_assays": 1,
                "required_assays": 6,
                "interpretation": "no adjusted coding-versus-lncRNA observability claim or Figure_1 panel is authorized",
            }
        ],
        ["gate", "status", "ready_assays", "required_assays", "interpretation"],
    )
    input_paths = [
        args.biotype_counts,
        args.biotype_review,
        args.eqtl_readiness,
        args.spatial_axis_audit,
        args.singlecell_gene_universe,
    ]
    write_tsv(
        args.output / "source_manifest.tsv",
        [
            {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in input_paths
        ],
        ["path", "size_bytes", "sha256"],
    )
    artifacts = sorted(path for path in args.output.iterdir() if path.is_file())
    write_tsv(
        args.output / "output_manifest.tsv",
        [
            {
                "relative_path": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in artifacts
        ],
        ["relative_path", "size_bytes", "sha256"],
    )
    (args.output / "VALIDATED.json").write_text(
        json.dumps(
            {
                "status": "pass",
                "scientific_gate": "skipped_incomplete_assay_matrix",
                "inference_run": False,
                "figure_1_adjusted_observability_authorized": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
