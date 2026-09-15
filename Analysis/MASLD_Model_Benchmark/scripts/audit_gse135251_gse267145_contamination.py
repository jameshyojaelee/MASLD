"""Audit GSE135251 against GSE267145 for participant and sample overlap.

Different stated provenance is not evidence.  GSE135251 is the European NAFLD
Registry and GSE267145 is a Vanderbilt bariatric cohort, but a claim of
disjointness has to rest on identifiers that were actually compared and on
expression profiles that were actually correlated.

Two layers:

1.  Identifier disjointness across every deposited accession namespace on both
    sides - GEO sample, SRA experiment, BioSample, BioProject - as exact set
    intersections.  Namespaces present on only one side are reported as such
    rather than counted as evidence.

2.  Molecular near-duplicate detection.  A shared participant deposited twice
    would show up as an anomalously similar pair of profiles.  Cross-cohort
    correlation is depressed by platform and batch, so no absolute threshold is
    meaningful; what is meaningful is whether the top cross-cohort pair stands
    apart from the cross-cohort bulk, and how it compares with the within-cohort
    tails, where every pair is known to be a different participant.  Reciprocal
    best hits are reported because a genuine duplicate is mutually top-ranked,
    while the arg-max of a noisy matrix is not.

The frozen disjointness statement in cross_cohort_expansion.toml covers
GSE260666 against GSE267145 and does not extend to GSE135251; this output file
exists because that gap had to be closed with its own evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import numpy as np


NAMESPACES = {
    "geo_sample": re.compile(r"GSM\d+"),
    "sra_experiment": re.compile(r"SRX\d+"),
    "sra_run": re.compile(r"SRR\d+"),
    "biosample": re.compile(r"SAMN\d+"),
    "bioproject": re.compile(r"PRJNA\d+"),
    "sra_study": re.compile(r"SRP\d+"),
}


class ContaminationAuditError(RuntimeError):
    """Raised when the audit cannot establish its own preconditions."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ContaminationAuditError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def read_rowname_offset_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read R `write.table(row.names=TRUE)` output, asserting the N+1 shape.

    The header carries N names and every data row carries N+1 fields because
    field 0 is an unnamed rowname.  A plain DictReader silently shifts every
    column by one and returns the run accession where the fibrosis stage
    belongs, so the invariant is asserted rather than trusted.
    """

    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    ragged = [index for index, row in enumerate(rows) if len(row) != len(header) + 1]
    if ragged:
        raise ContaminationAuditError(
            f"{path} violates the rowname N+1 invariant at data rows {ragged[:5]}"
        )
    return header, rows


def scan_namespaces(text: str) -> dict[str, set[str]]:
    return {name: set(pattern.findall(text)) for name, pattern in NAMESPACES.items()}


def merge_namespaces(parts: Sequence[Mapping[str, set[str]]]) -> dict[str, set[str]]:
    merged: dict[str, set[str]] = {name: set() for name in NAMESPACES}
    for part in parts:
        for name, values in part.items():
            merged[name] |= values
    return merged


def per_sample_rank(matrix: np.ndarray) -> np.ndarray:
    """Within-sample rank over the shared gene axis, then standardised.

    Ranking inside one sample removes library size and any per-sample scale, so
    the comparison is not driven by depth differences between two cohorts
    sequenced years apart.  Correlating the standardised ranks is Spearman.
    """

    from scipy.stats import rankdata

    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)):
        raise ContaminationAuditError("expression matrix is not finite")
    ranks = rankdata(values, method="average", axis=1)
    centred = ranks - ranks.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(centred, axis=1, keepdims=True)
    if np.any(norm <= 0):
        raise ContaminationAuditError("a sample has no rank variation")
    return centred / norm


def describe(values: np.ndarray) -> dict[str, float]:
    flat = np.asarray(values, dtype=np.float64).ravel()
    return {
        "n_pairs": int(flat.size),
        "min": float(np.min(flat)),
        "p50": float(np.quantile(flat, 0.50)),
        "p95": float(np.quantile(flat, 0.95)),
        "p99": float(np.quantile(flat, 0.99)),
        "p999": float(np.quantile(flat, 0.999)),
        "max": float(np.max(flat)),
        "mean": float(np.mean(flat)),
        "sd": float(np.std(flat, ddof=1)) if flat.size > 1 else 0.0,
    }


def upper_triangle(matrix: np.ndarray) -> np.ndarray:
    index = np.triu_indices(matrix.shape[0], k=1)
    return matrix[index]


def run_audit(
    *, benchmark_root: Path, source_root: Path, target_root: Path, join_root: Path,
    output: Path,
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise ContaminationAuditError(f"refusing to overwrite audit: {output}")
    for root in (source_root / "molecular", source_root / "raw", target_root):
        verify_frozen_tree(root)

    # ---- layer 1: identifiers -------------------------------------------------
    source_meta = source_root / "raw" / "GSE135251_metadata.tsv"
    header, rows = read_rowname_offset_tsv(source_meta)
    source_text = source_meta.read_text(encoding="utf-8")
    source_ids = scan_namespaces(source_text)

    _, source_participants = read_tsv(source_root / "molecular" / "participant_axis.tsv")
    _, target_participants = read_tsv(target_root / "participant_axis.tsv")
    join_text = (join_root / "participant_join.tsv").read_text(encoding="utf-8")
    with gzip.open(join_root / "raw" / "GSE267145_family.soft.gz", "rt", errors="replace") as handle:
        soft_text = handle.read()
    target_ids = merge_namespaces(
        [
            scan_namespaces(join_text),
            scan_namespaces(soft_text),
            scan_namespaces(
                "\n".join(
                    f"{row.get('rna_source_sample_accession','')} "
                    f"{row.get('h3k27ac_source_sample_accession','')}"
                    for row in target_participants
                )
            ),
        ]
    )

    identifier_rows: list[dict[str, Any]] = []
    overlaps: dict[str, list[str]] = {}
    for name in sorted(NAMESPACES):
        left, right = source_ids[name], target_ids[name]
        shared = sorted(left & right)
        if shared:
            overlaps[name] = shared
        identifier_rows.append(
            {
                "namespace": name,
                "gse135251_unique": len(left),
                "gse267145_unique": len(right),
                "shared": len(shared),
                "comparable": bool(left and right),
                "evidence_state": (
                    "compared" if left and right else "absent_on_one_side_not_evidence"
                ),
            }
        )

    # Participant identifiers live in different namespaces, so a raw string
    # comparison of them proves nothing on its own; it is recorded, not counted.
    source_pids = {row["participant_id"] for row in source_participants}
    target_pids = {row["participant_id"] for row in target_participants}
    participant_string_overlap = sorted(source_pids & target_pids)

    # ---- reconcile the deposited counts columns -------------------------------
    with (source_root / "raw" / "GSE135251_counts.tsv").open(encoding="utf-8") as handle:
        counts_header = handle.readline().rstrip("\n").split("\t")
    counts_columns = counts_header[1:]
    metadata_runs = {row[header.index("Run") + 1] for row in rows}
    metadata_gsms = {row[header.index("GEO_Accession (exp)") + 1] for row in rows}
    matched = [c for c in counts_columns if c in metadata_runs or c in metadata_gsms]
    unmatched = [c for c in counts_columns if c not in metadata_runs and c not in metadata_gsms]
    unmatched_ids = scan_namespaces("\n".join(unmatched))
    unmatched_in_target = {
        name: sorted(values & target_ids[name]) for name, values in unmatched_ids.items()
    }

    # ---- layer 2: molecular near-duplicates -----------------------------------
    _, source_axis = read_tsv(source_root / "molecular" / "rna_feature_axis.tsv")
    _, target_axis = read_tsv(target_root / "rna_feature_axis.tsv")
    source_genes = [row["stable_gene_id"] for row in source_axis]
    target_genes = [row["stable_gene_id"] for row in target_axis]
    shared_genes = sorted(set(source_genes) & set(target_genes))
    if len(shared_genes) < 1000:
        raise ContaminationAuditError("shared gene axis is too small to audit")
    source_lookup = {value: index for index, value in enumerate(source_genes)}
    target_lookup = {value: index for index, value in enumerate(target_genes)}

    source_values = np.load(
        source_root / "molecular" / "rna_values.npy", mmap_mode="r", allow_pickle=False
    )
    target_values = np.load(
        target_root / "rna_values.npy", mmap_mode="r", allow_pickle=False
    )
    source_matrix = np.asarray(source_values)[:, [source_lookup[g] for g in shared_genes]]
    target_matrix = np.asarray(target_values)[:, [target_lookup[g] for g in shared_genes]]
    source_ranked = per_sample_rank(source_matrix)
    target_ranked = per_sample_rank(target_matrix)

    cross = source_ranked @ target_ranked.T
    within_source = source_ranked @ source_ranked.T
    within_target = target_ranked @ target_ranked.T

    cross_stats = describe(cross)
    within_source_stats = describe(upper_triangle(within_source))
    within_target_stats = describe(upper_triangle(within_target))

    flat = np.sort(cross.ravel())[::-1]
    top_gap = float(flat[0] - flat[1])
    typical_gap = float(np.median(flat[:200] - flat[1:201]))
    max_z = float((cross_stats["max"] - cross_stats["mean"]) / cross_stats["sd"])

    source_best = np.argmax(cross, axis=1)
    target_best = np.argmax(cross, axis=0)
    reciprocal = [
        (int(i), int(source_best[i]))
        for i in range(cross.shape[0])
        if target_best[source_best[i]] == i
    ]
    source_ids_axis = [row["participant_id"] for row in source_participants]
    target_ids_axis = [row["participant_id"] for row in target_participants]
    reciprocal_rows = sorted(
        (
            {
                "gse135251_participant": source_ids_axis[i],
                "gse267145_participant": target_ids_axis[j],
                "spearman": float(cross[i, j]),
                "z_vs_cross_cohort": float(
                    (cross[i, j] - cross_stats["mean"]) / cross_stats["sd"]
                ),
            }
            for i, j in reciprocal
        ),
        key=lambda row: -row["spearman"],
    )

    # A duplicate would have to exceed what different participants reach inside a
    # single cohort and single protocol.  That tail is the reference.
    duplicate_reference = max(
        within_source_stats["max"], within_target_stats["max"]
    )
    exceeds_reference = [
        row for row in reciprocal_rows if row["spearman"] >= duplicate_reference
    ]

    output.mkdir(parents=True)
    with (output / "identifier_overlap.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "namespace",
                "gse135251_unique",
                "gse267145_unique",
                "shared",
                "comparable",
                "evidence_state",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(identifier_rows)
    with (output / "reciprocal_best_hits.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "gse135251_participant",
                "gse267145_participant",
                "spearman",
                "z_vs_cross_cohort",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(reciprocal_rows)

    total_shared = sum(len(v) for v in overlaps.values())
    receipt = {
        "schema_version": "masld-bench-gse135251-gse267145-contamination-audit-v1",
        "status": (
            "pass_no_overlap_detected"
            if total_shared == 0 and not exceeds_reference
            else "fail_overlap_detected"
        ),
        "source_series": "GSE135251",
        "target_series": "GSE267145",
        "prior_frozen_statement_covers": "GSE260666_vs_GSE267145",
        "prior_frozen_statement_extends_to_gse135251": False,
        "identifier_layer": {
            "namespaces_compared": [
                row["namespace"] for row in identifier_rows if row["comparable"]
            ],
            "namespaces_absent_on_one_side": [
                row["namespace"] for row in identifier_rows if not row["comparable"]
            ],
            "shared_identifiers_total": total_shared,
            "shared_by_namespace": {k: v for k, v in sorted(overlaps.items())},
            "participant_id_string_overlap": participant_string_overlap,
            "participant_id_namespaces_differ": True,
        },
        "counts_column_reconciliation": {
            "counts_columns_deposited": len(counts_columns),
            "columns_matched_to_metadata": len(matched),
            "columns_without_metadata": len(unmatched),
            "unmatched_columns": sorted(unmatched),
            "unmatched_columns_present_in_gse267145": {
                k: v for k, v in unmatched_in_target.items() if v
            },
            "unmatched_columns_in_frozen_source": False,
        },
        "molecular_layer": {
            "shared_genes": len(shared_genes),
            "gse135251_samples": int(source_ranked.shape[0]),
            "gse267145_samples": int(target_ranked.shape[0]),
            "similarity": "within_sample_rank_then_pearson_equals_spearman",
            "cross_cohort": cross_stats,
            "within_gse135251": within_source_stats,
            "within_gse267145": within_target_stats,
            "top_cross_cohort_gap": top_gap,
            "median_gap_in_top_200": typical_gap,
            "top_pair_z_vs_cross_cohort": max_z,
            "reciprocal_best_hits": len(reciprocal_rows),
            "duplicate_reference_within_cohort_max": duplicate_reference,
            "reciprocal_best_hits_exceeding_reference": exceeds_reference,
            "interpretation": (
                "no cross-cohort pair reaches the similarity that different "
                "participants already reach inside a single cohort, so there is "
                "no molecular near-duplicate signal"
                if not exceeds_reference
                else "at least one cross-cohort pair is as similar as a "
                "within-cohort pair and must be resolved before any transfer"
            ),
        },
        "provenance_asserted_not_counted_as_evidence": (
            "GSE135251 is the European NAFLD Registry and GSE267145 a Vanderbilt "
            "bariatric cohort; differing stated provenance motivated this audit "
            "but is not part of its evidence"
        ),
        "inputs_sha256": {
            "gse135251_source_artifacts": sha256_file(source_root / "ARTIFACTS.json"),
            "gse267145_molecular_artifacts": sha256_file(target_root / "ARTIFACTS.json"),
            "gse135251_metadata": sha256_file(source_meta),
        },
    }
    with (output / "audit.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--target-root", required=True, type=Path)
    parser.add_argument("--join-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run_audit(
        benchmark_root=arguments.benchmark_root,
        source_root=arguments.source_root,
        target_root=arguments.target_root,
        join_root=arguments.join_root,
        output=arguments.output,
    )
    print(json.dumps({k: v for k, v in receipt.items() if k != "molecular_layer"}, sort_keys=True))
    print(json.dumps(receipt["molecular_layer"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
