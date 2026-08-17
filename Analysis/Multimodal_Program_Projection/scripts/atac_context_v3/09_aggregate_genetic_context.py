#!/usr/bin/env python3
"""Validate lifted alleles and aggregate variant-consistent lineage mass."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import os
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
    project_root,
    sha256_file,
    write_tsv,
)
from genetics_context import normalize_against_reference


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
    stale_pending = sorted((root / "genetics").glob(".context.pending.*"))
    if stale_pending:
        raise ContractError(f"unresolved prior genetics context stage: {stale_pending}")
    pending = root / f"genetics/.context.pending.{os.getpid()}"
    if pending.exists() or pending.is_symlink():
        raise ContractError(f"pending genetics context exists: {pending}")
    pending.mkdir(parents=True)

    rows = read_rows(source)
    peaks = read_rows(peak_path)
    if not rows or not peaks:
        raise ContractError("liftover posterior or consensus peak table is empty")
    frozen = read_rows(root / "input_manifest.tsv")
    fasta_rows = [row for row in frozen if row["role"] == "hg38_fasta"]
    if len(fasta_rows) != 1:
        raise ContractError("frozen hg38 FASTA is not unique")
    fasta_path = Path(fasta_rows[0]["relative_or_absolute_path"])
    if not fasta_path.is_absolute():
        fasta_path = project_root() / fasta_path
    if not fasta_path.is_file() or sha256_file(fasta_path) != fasta_rows[0]["sha256"]:
        raise ContractError("frozen hg38 FASTA hash mismatch")
    index244 = PeakIndex(peaks, "supported_gse244832")
    index281 = PeakIndex(peaks, "supported_gse281367")
    fasta = pysam.FastaFile(str(fasta_path))
    audit = []
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    posterior_keys: set[tuple[str, str, int, str]] = set()
    for number, row in enumerate(rows, start=1):
        posterior = float(row["SNP.PP.H4"])
        if not math.isfinite(posterior) or not 0 <= posterior <= 1:
            raise ContractError("variant posterior contains an invalid probability")
        posterior_key = (
            row["gwas_name"], row["ensembl"], int(row["signal_pair_index"]), row["snp"]
        )
        if posterior_key in posterior_keys:
            raise ContractError(f"duplicated variant posterior row: {posterior_key}")
        posterior_keys.add(posterior_key)
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
    for path in sorted((root / "genetics/replay_execution").glob(
        "batch_*/exports/*/chr*/*/signal_pairs.tsv"
    )):
        for row in read_rows(path):
            key = (row["gwas_name"], row["ensembl"], int(row["signal_pair_index"]))
            if key in summaries:
                raise ContractError(f"duplicated signal-pair summary: {key}")
            summaries[key] = row
    if set(grouped) != set(summaries):
        raise ContractError("signal-pair summaries and variant posterior families differ")
    plan_rows = read_rows(root / "genetics/prepared/replay_plan.tsv")
    if not plan_rows or {row["trait_class"] for row in plan_rows} != {
        "direct_MASLD", "liver_enzyme"
    }:
        raise ContractError("replay plan lacks both prespecified trait classes")
    plan_keys = [(row["gwas_name"], row["ensembl"]) for row in plan_rows]
    if len(plan_keys) != len(set(plan_keys)):
        raise ContractError("replay plan contains duplicated gene-study pairs")
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
                "gse244832_unique_accessible_mass": mass244 - shared,
                "gse281367_unique_accessible_mass": mass281 - shared,
                "shared_accessible_mass": shared,
                "any_accessible_mass": any_mass,
                "joint_gse244832_accessible_mass": pp_h4 * mass244,
                "joint_gse281367_accessible_mass": pp_h4 * mass281,
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

    planned_pairs = {(row["gwas_name"], row["ensembl"]) for row in plan_rows}
    observed_pairs = {(key[0], key[1]) for key in grouped}
    if observed_pairs != planned_pairs:
        raise ContractError("replay exports do not cover every planned gene-study pair")
    if not audit or not lineage_rows:
        raise ContractError("genetics aggregation produced no biological rows")
    write_tsv(pending / "variant_liftover_audit.tsv", tuple(audit[0]), audit)
    write_tsv(
        pending / "genetic_lineage_context_all_pairs.tsv",
        tuple(lineage_rows[0]),
        lineage_rows,
    )
    primary = [row for row in lineage_rows if row["primary_signal_pair"] == "TRUE"]
    write_tsv(
        pending / "genetic_lineage_context_primary_pairs.tsv",
        tuple(primary[0]),
        primary,
    )
    write_tsv(
        pending / "genetics_context_input_manifest.tsv",
        ("release_id", "role", "path", "sha256"),
        [{
            "release_id": RELEASE_ID,
            "role": "hg38_fasta",
            "path": str(fasta_path),
            "sha256": fasta_rows[0]["sha256"],
        }],
    )
    os.replace(pending, out)
    print(f"Aggregated {len(grouped)} signal pairs across {len(LINEAGES)} lineages")


if __name__ == "__main__":
    main()
