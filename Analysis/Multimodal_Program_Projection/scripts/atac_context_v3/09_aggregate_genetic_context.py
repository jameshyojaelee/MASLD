#!/usr/bin/env python3
"""Validate lifted alleles and aggregate variant-consistent lineage mass."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import pysam

from atac_context_v3_lib import (
    LINEAGES,
    MIN_MAPPED_POSTERIOR,
    RELEASE_ID,
    STANDARD_CHROMS,
    ContractError,
    candidate_root,
    classify_genetic,
    default_candidate_root,
    write_tsv,
)
from genetics_context import normalize_against_reference


HG38_FASTA = Path(
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PeakIndex:
    def __init__(self, rows: list[dict[str, str]], cohort_field: str):
        values: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
        for row in rows:
            if row[cohort_field] != "TRUE":
                continue
            values[(row["lineage"], row["chrom"])].append(
                (int(row["start0"]), int(row["end"]))
            )
        self.starts: dict[tuple[str, str], list[int]] = {}
        self.prefix_end: dict[tuple[str, str], list[int]] = {}
        for key, intervals in values.items():
            intervals.sort()
            self.starts[key] = [value[0] for value in intervals]
            running = -1
            prefix = []
            for _, end in intervals:
                running = max(running, end)
                prefix.append(running)
            self.prefix_end[key] = prefix

    def overlaps(self, lineage: str, chrom: str, position_1based: int) -> bool:
        key = (lineage, chrom)
        starts = self.starts.get(key, [])
        upper = bisect.bisect_left(starts, position_1based)
        return upper > 0 and self.prefix_end[key][upper - 1] >= position_1based


def main() -> None:
    args = arguments()
    root = candidate_root(args.candidate_root)
    source = root / "genetics/liftover/variant_liftover_raw.tsv.gz"
    peak_path = root / "consensus_peak_manifest.tsv"
    if not source.is_file() or not peak_path.is_file():
        raise ContractError("liftover and consensus peaks must exist before genetics aggregation")
    out = root / "genetics/context"
    if out.exists():
        raise ContractError(f"refusing to overwrite genetics context: {out}")
    out.mkdir(parents=True)

    rows = read_rows(source)
    peaks = read_rows(peak_path)
    index244 = PeakIndex(peaks, "supported_gse244832")
    index281 = PeakIndex(peaks, "supported_gse281367")
    fasta = pysam.FastaFile(str(HG38_FASTA))
    audit = []
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    for number, row in enumerate(rows, start=1):
        posterior = float(row["SNP.PP.H4"])
        status = "mapped"
        normalized = None
        if row.get("unique_liftover", "").upper() != "TRUE":
            status = "unmapped_or_multimapped"
        elif row.get("hg38_chrom", "") not in STANDARD_CHROMS:
            status = "nonstandard_chromosome"
        else:
            normalized = normalize_against_reference(
                row["hg38_chrom"],
                int(row["hg38_position_1based"]),
                row["allele1"],
                row["allele2"],
                fasta.fetch,
            )
            if normalized is None:
                status = "hg38_allele_validation_failed"
        record = {
            "release_id": RELEASE_ID,
            "posterior_row_id": number,
            "gwas_name": row["gwas_name"],
            "gene": row["gene"],
            "ensembl": row["ensembl"],
            "signal_pair_index": row["signal_pair_index"],
            "hg19_variant_id": row["snp"],
            "allele1": row["allele1"],
            "allele2": row["allele2"],
            "snp_pp_h4": posterior,
            "mapping_status": status,
            "hg38_chrom": normalized and row["hg38_chrom"] or "",
            "hg38_position_1based": normalized and normalized.position_1based or "",
            "hg38_ref": normalized and normalized.ref or "",
            "hg38_alt": normalized and normalized.alt or "",
            "orientation": normalized and normalized.orientation or "",
        }
        audit.append(record)
        grouped[(row["gwas_name"], row["ensembl"], int(row["signal_pair_index"]))].append(
            {
                "posterior": posterior,
                "mapped": normalized is not None,
                "chrom": row.get("hg38_chrom", ""),
                "position_1based": normalized.position_1based if normalized else 0,
            }
        )
    fasta.close()

    summaries: dict[tuple[str, str, int], dict[str, str]] = {}
    for path in sorted((root / "genetics/replay_exports").glob("*/chr*/*/signal_pairs.tsv")):
        for row in read_rows(path):
            key = (row["gwas_name"], row["ensembl"], int(row["signal_pair_index"]))
            summaries[key] = row
    if set(grouped) - set(summaries):
        raise ContractError("variant posterior rows lack signal-pair summaries")
    plan_rows = read_rows(root / "genetics/replay_plan.tsv")
    pair_metadata = {
        (row["gwas_name"], row["ensembl"]): row
        for row in plan_rows
    }

    lineage_rows = []
    for key, variants in sorted(grouped.items()):
        original_mass = sum(float(row["posterior"]) for row in variants)
        if abs(original_mass - 1.0) > 1e-6:
            raise ContractError(f"SNP.PP.H4 does not sum to one for {key}: {original_mass}")
        summary = summaries[key]
        metadata = pair_metadata.get((key[0], key[1]))
        if metadata is None:
            raise ContractError(f"signal pair absent from frozen replay plan: {key}")
        pp_h4 = float(summary["PP.H4.abf"])
        for lineage in LINEAGES:
            mapped = sum(float(row["posterior"]) for row in variants if row["mapped"])
            mass244 = sum(
                float(row["posterior"])
                for row in variants
                if row["mapped"] and index244.overlaps(lineage, str(row["chrom"]), int(row["position_1based"]))
            )
            mass281 = sum(
                float(row["posterior"])
                for row in variants
                if row["mapped"] and index281.overlaps(lineage, str(row["chrom"]), int(row["position_1based"]))
            )
            shared = sum(
                float(row["posterior"])
                for row in variants
                if row["mapped"]
                and index244.overlaps(lineage, str(row["chrom"]), int(row["position_1based"]))
                and index281.overlaps(lineage, str(row["chrom"]), int(row["position_1based"]))
            )
            any_mass = sum(
                float(row["posterior"])
                for row in variants
                if row["mapped"]
                and (
                    index244.overlaps(lineage, str(row["chrom"]), int(row["position_1based"]))
                    or index281.overlaps(lineage, str(row["chrom"]), int(row["position_1based"]))
                )
            )
            lineage_rows.append({
                "release_id": RELEASE_ID,
                "gwas_name": key[0],
                "gene": summary["gene"],
                "ensembl": key[1],
                "trait": metadata["trait"],
                "trait_class": metadata["trait_class"],
                "ancestry": metadata["ancestry"],
                "signal_pair_index": key[2],
                "gwas_signal": summary["gwas_signal"],
                "eqtl_signal": summary["eqtl_signal"],
                "pp_h4": pp_h4,
                "lineage": lineage,
                "mapped_posterior_mass": mapped,
                "lost_posterior_mass": 1.0 - mapped,
                "gse244832_accessible_mass": mass244,
                "gse281367_accessible_mass": mass281,
                "shared_accessible_mass": shared,
                "any_accessible_mass": any_mass,
                "joint_any_accessible_mass": pp_h4 * any_mass,
                "joint_shared_accessible_mass": pp_h4 * shared,
                "evidence_state": classify_genetic(mapped, shared, mass244, mass281),
                "testability_reason": "" if mapped >= MIN_MAPPED_POSTERIOR else "mapped_posterior_mass_below_0.95",
                "lineage_label_note": (
                    "GSE281367_combined_T_NK_label"
                    if lineage == "t_nk"
                    else ("broad_macrophage_lineage" if lineage == "macrophage" else "")
                ),
            })

    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in lineage_rows:
        by_pair[(str(row["gwas_name"]), str(row["ensembl"]))].append(row)
    primary_index: dict[tuple[str, str], int] = {}
    for key, values in by_pair.items():
        candidates = {(int(row["signal_pair_index"]), float(row["pp_h4"]), str(row["gwas_signal"]), str(row["eqtl_signal"])) for row in values}
        selected = sorted(candidates, key=lambda value: (-value[1], value[2], value[3], value[0]))[0]
        primary_index[key] = selected[0]
    for row in lineage_rows:
        row["primary_signal_pair"] = str(
            int(row["signal_pair_index"]) == primary_index[(str(row["gwas_name"]), str(row["ensembl"]))]
        ).upper()

    write_tsv(out / "variant_liftover_audit.tsv", tuple(audit[0]), audit)
    write_tsv(out / "genetic_lineage_context_all_pairs.tsv", tuple(lineage_rows[0]), lineage_rows)
    primary = [row for row in lineage_rows if row["primary_signal_pair"] == "TRUE"]
    write_tsv(out / "genetic_lineage_context_primary_pairs.tsv", tuple(primary[0]), primary)
    print(f"Aggregated {len(grouped)} signal pairs across {len(LINEAGES)} lineages")


if __name__ == "__main__":
    main()
