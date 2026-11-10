#!/usr/bin/env python3
"""Annotate accessibility links with exact-peak frozen motif matches.

Motif compatibility is neither unique TF attribution, TF activity, measured
binding, nor evidence of causal regulation. No cross-peak overlap borrowing.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import socket
import sqlite3


DEFAULT_MOTIF_ROOT = Path(__file__).resolve().parents[3] / "Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/chromvar"


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jaspar_metadata(root):
    records = rows(root / "jaspar_database.tsv")
    if len(records) != 1:
        raise ValueError("Require exactly one frozen JASPAR database record")
    record = records[0]
    database = Path(record["database_path"])
    if digest(database) != record["database_sha256"]:
        raise ValueError("Frozen JASPAR database SHA256 differs")
    connection = sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        matrices = connection.execute("SELECT ID, BASE_ID, VERSION, NAME FROM MATRIX").fetchall()
        annotations = {}
        for ident, tag, value in connection.execute("SELECT ID, TAG, VAL FROM MATRIX_ANNOTATION WHERE TAG IN ('family', 'class')"):
            annotations.setdefault(ident, {}).setdefault(tag, set()).add(value or "")
        metadata = {}
        for ident, base, version, name in matrices:
            if f"{base}.{version}" in metadata:
                raise ValueError("Duplicate JASPAR motif ID")
            annotation = annotations.get(ident, {})
            metadata[f"{base}.{version}"] = dict(motif_name=name,
                motif_family=";".join(sorted(annotation.get("family", []))),
                motif_class=";".join(sorted(annotation.get("class", []))))
    finally:
        connection.close()
    return metadata, record


def selected_motif_matches(directory, selected_peaks):
    """One streaming sparse MatrixMarket pass; retain requested rows only."""
    peak_rows = rows(directory / "peaks.tsv")
    motif_rows = rows(directory / "motifs.tsv")
    peaks = [r["peak_coordinate"] for r in peak_rows]
    motifs = [r["motif_id"] for r in motif_rows]
    if len(set(peaks)) != len(peaks) or len(set(motifs)) != len(motifs):
        raise ValueError("Duplicate frozen peak or motif IDs")
    selected = {i + 1: p for i, p in enumerate(peaks) if p in selected_peaks}
    matches = {p: set() for p in selected.values()}
    with (directory / "matches.mtx").open() as handle:
        header = handle.readline().strip().lower().split()
        if len(header) != 5 or header[:3] != ["%%matrixmarket", "matrix", "coordinate"] or header[4] != "general":
            raise ValueError("Expected general coordinate MatrixMarket motif matrix")
        if header[3] not in {"pattern", "integer", "real"}:
            raise ValueError("Unsupported motif matrix type")
        line = next(handle)
        while line.startswith("%"):
            line = next(handle)
        n_peaks, n_motifs, n_entries = map(int, line.split())
        if (n_peaks, n_motifs) != (len(peaks), len(motifs)):
            raise ValueError("Frozen motif matrix dimensions disagree with sidecars")
        count = 0
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            values = line.split()
            row, col = map(int, values[:2])
            if not (1 <= row <= n_peaks and 1 <= col <= n_motifs):
                raise ValueError("Motif matrix index outside dimensions")
            count += 1
            if row in selected and (len(values) == 2 or float(values[2]) > 0):
                matches[selected[row]].add(motifs[col - 1])
        if count != n_entries:
            raise ValueError("Motif matrix entry count disagrees with header")
    return matches


def annotate(network_roots, motif_root, output):
    if output.exists():
        raise FileExistsError(f"Refusing overwrite: {output}")
    for network in network_roots:
        if not (network / "completion.json").exists():
            raise ValueError(f"Network run incomplete: {network}")
        if json.loads((network / "completion.json").read_text()).get("status") != "completed" or any(network.glob("INVALIDATED*.json")):
            raise ValueError(f"Network run incomplete or invalidated: {network}")
    metadata, record = jaspar_metadata(motif_root)
    output.mkdir(parents=True, exist_ok=False)
    cache, coverage = {}, []
    fields = ["network_root", "cohort", "lineage", "lineage_state", "source_lineage", "outer_fold", "gene_id", "peak_id",
              "method", "rank", "score", "score_unit", "link_evidence_state", "motif_id", "motif_name", "motif_family", "motif_class",
              "evidence_state", "motif_not_unique_tf", "tf_activity_inferred", "tf_binding_inferred", "causal_regulation_inferred",
              "jaspar_package", "jaspar_database_sha256"]
    # Collect contexts first so each frozen matrix is scanned only once.
    contexts = {}
    for network in network_roots:
        lineage_meta = {(r["cohort"], r["lineage"]): r for r in rows(network / "comparisons.tsv")} if (network / "comparisons.tsv").exists() else {}
        for directory in sorted(network.iterdir()):
            path = directory / "links.tsv"
            if not directory.is_dir() or not path.exists() or path.stat().st_size < 2:
                continue
            cohort, lineage = directory.name.split("_", 1)
            links = rows(path)
            required = {"outer_fold", "gene_id", "peak_id", "method", "rank", "score", "score_unit", "evidence_state"}
            if links and not required.issubset(links[0]):
                raise ValueError(f"Missing link columns: {path}")
            entry = dict(network_root=str(network.resolve()), cohort=cohort, lineage=lineage,
                         lineage_state=lineage_meta.get((cohort, lineage), {}).get("lineage_state", "source_label"),
                         source_lineage=lineage_meta.get((cohort, lineage), {}).get("source_lineage", lineage), links=links)
            contexts.setdefault((cohort, lineage), []).append(entry)
    with gzip.open(output / "motif_peak_gene.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (cohort, lineage), entries in contexts.items():
            directory = motif_root / "motif_matches" / cohort / lineage
            available = all((directory / name).exists() for name in ["matches.mtx", "peaks.tsv", "motifs.tsv"])
            peaks = {r["peak_id"] for entry in entries for r in entry["links"]}
            if available:
                cache[(cohort, lineage)] = selected_motif_matches(directory, peaks)
            matches = cache.get((cohort, lineage), {})
            for entry in entries:
                n_present = n_compatible = n_rows = 0
                for link in entry["links"]:
                    motif_ids = matches.get(link["peak_id"], set())
                    n_present += int(link["peak_id"] in matches)
                    n_compatible += bool(motif_ids)
                    for motif in sorted(motif_ids):
                        if motif not in metadata:
                            raise ValueError(f"Matched motif missing from frozen JASPAR database: {motif}")
                        row = {key: entry[key] for key in ["network_root", "cohort", "lineage", "lineage_state", "source_lineage"]}
                        row.update({key: link[key] for key in required if key != "evidence_state"})
                        row.update(metadata[motif])
                        row.update(motif_id=motif, link_evidence_state=link["evidence_state"],
                                   evidence_state="motif_compatible_accessibility_link", motif_not_unique_tf=True,
                                   tf_activity_inferred=False, tf_binding_inferred=False, causal_regulation_inferred=False,
                                   jaspar_package=record["jaspar_package"], jaspar_database_sha256=record["database_sha256"])
                        writer.writerow(row)
                        n_rows += 1
                coverage.append({key: entry[key] for key in ["network_root", "cohort", "lineage", "lineage_state", "source_lineage"]} |
                    dict(status="annotated" if available else "untestable", reason="" if available else "native_peak_motifscan_not_run",
                         n_links=len(entry["links"]), n_exact_peak_links=n_present, n_links_with_motif=n_compatible,
                         n_links_absent_from_motif_universe=len(entry["links"]) - n_present if available else "",
                         n_motif_peak_gene_rows=n_rows, peak_matching="exact_ID_only_no_overlap_borrowing"))
            # Bound retained sparse context to one cohort/lineage at a time.
            cache.pop((cohort, lineage), None)
    coverage_fields = ["network_root", "cohort", "lineage", "lineage_state", "source_lineage", "status", "reason", "n_links",
                       "n_exact_peak_links", "n_links_with_motif", "n_links_absent_from_motif_universe", "n_motif_peak_gene_rows", "peak_matching"]
    with (output / "coverage.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=coverage_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(coverage)
    (output / "completion.json").write_text(json.dumps(dict(status="completed", scientific_adoption=False,
        interpretation="motif_compatibility_only_no_unique_TF_activity_binding_or_causality",
        motif_root=str(motif_root.resolve()), jaspar_source=record), indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", "--network-root", type=Path, action="append", required=True)
    parser.add_argument("--motif-root", type=Path, default=DEFAULT_MOTIF_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname():
        raise RuntimeError("Run motif-matrix annotation in a compute allocation")
    annotate(args.network, args.motif_root, args.output)


if __name__ == "__main__":
    main()
