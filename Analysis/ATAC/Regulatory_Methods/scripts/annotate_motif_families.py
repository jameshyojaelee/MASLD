#!/usr/bin/env python3
"""Stream exact-peak motif compatibility into TF-family/gene summaries.

Only original link records are sorted on disk. Motif/peak/gene Cartesian rows
are never created. Each family summary retains the source peak and motif IDs,
and its maximum link score retains the original method's continuous units.
Compatibility does not establish TF identity, activity, binding or causality.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import socket
import sqlite3
import tempfile

import annotate_motif_context as context


REQUIRED = {"outer_fold", "gene_id", "peak_id", "method", "rank", "score", "score_unit", "evidence_state"}
GROUP = ["outer_fold", "gene_id", "method", "endpoint", "input_representation", "score_unit"]
FIELDS = ["network_root", "source_scope", "source_interpretation", "cohort", "lineage", "lineage_state", "source_lineage",
          "outer_fold", "gene_id", "endpoint", "endpoint_recorded", "input_representation", "method", "motif_family",
          "family_annotation_state", "distinct_supporting_peak_count", "supporting_peak_ids", "max_peak_link_score",
          "score_unit", "motif_ids", "motif_names", "motif_classes", "link_evidence_states", "evidence_state",
          "motif_not_unique_tf", "tf_activity_inferred", "tf_binding_inferred", "causal_regulation_inferred",
          "jaspar_package", "jaspar_database_sha256"]
COVERAGE = ["network_root", "source_scope", "source_interpretation", "cohort", "lineage", "lineage_state", "source_lineage",
            "endpoint_recorded", "status", "reason", "n_link_records", "n_unique_selected_peak_ids", "n_exact_peak_link_records",
            "n_links_with_motif", "n_links_absent_from_motif_universe", "n_family_gene_rows", "peak_matching"]


def validated_completion(network: Path) -> dict:
    path = network / "completion.json"
    if not path.is_file() or any(network.glob("INVALIDATED*")):
        raise ValueError(f"Incomplete or invalidated network source: {network}")
    completion = json.loads(path.read_text())
    if completion.get("status") != "completed":
        raise ValueError(f"Network source is not completed: {network}")
    return completion


def source_contexts(networks: list[Path]) -> tuple[list[dict], list[dict]]:
    contexts, sources = [], []
    if len(set(network.resolve() for network in networks)) != len(networks):
        raise ValueError("Repeated network root would duplicate compatibility summaries")
    for network in networks:
        network = network.resolve()
        completion = validated_completion(network)
        scope = "technical_pilot" if completion.get("pilot") is True else "development_candidate" if completion.get("pilot") is False else "scope_not_declared"
        common = dict(network_root=str(network), source_scope=scope,
                      source_interpretation=completion.get("interpretation", "not_declared"))
        sources.append(dict(**common, completion=completion))
        comparison_path = network / "comparisons.tsv"
        lineage_meta = {(row["cohort"], row["lineage"]): row for row in context.rows(comparison_path)} if comparison_path.exists() else {}
        for directory in sorted(network.iterdir()):
            path = directory / "links.tsv"
            if not directory.is_dir() or not path.is_file() or path.stat().st_size < 2:
                continue
            if any(directory.glob("INVALIDATED*")):
                raise ValueError(f"Invalidated network context: {directory}")
            cohort, lineage = directory.name.split("_", 1)
            with path.open(newline="") as handle:
                header = csv.DictReader(handle, delimiter="\t").fieldnames or []
            if not REQUIRED.issubset(header):
                raise ValueError(f"Missing link columns: {path}")
            meta = lineage_meta.get((cohort, lineage), {})
            contexts.append(dict(**common, cohort=cohort, lineage=lineage,
                lineage_state=meta.get("lineage_state", "source_label"), source_lineage=meta.get("source_lineage", lineage),
                endpoint_recorded="endpoint" in header, path=path))
    return contexts, sources


def stage_links(entry: dict, connection: sqlite3.Connection | None) -> tuple[set[str], int]:
    """Bound sorting memory; this table contains unexpanded source links only."""
    if connection is not None:
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA cache_size=-32768")
        connection.execute("CREATE TABLE links (outer_fold TEXT, gene_id TEXT, method TEXT, endpoint TEXT, input_representation TEXT, score_unit TEXT, peak_id TEXT, score REAL, evidence TEXT)")
    before = entry["path"].stat()
    selected, batch, count = set(), [], 0
    with entry["path"].open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("cohort", entry["cohort"]) != entry["cohort"] or row.get("lineage", entry["lineage"]) != entry["lineage"]:
                raise ValueError(f"Link row changes cohort/lineage: {entry['path']}")
            if any(not row.get(key, "") for key in REQUIRED):
                raise ValueError(f"Blank required link identity/score field: {entry['path']}")
            score = float(row["score"])
            if not math.isfinite(score):
                raise ValueError(f"Nonfinite continuous link score: {entry['path']}")
            if entry["endpoint_recorded"] and not row.get("endpoint"):
                raise ValueError("Source declares endpoint column but leaves an endpoint blank")
            selected.add(row["peak_id"])
            if connection is not None:
                batch.append(tuple(row.get(key, "") for key in GROUP) + (row["peak_id"], score, row["evidence_state"]))
            count += 1
            if len(batch) >= 10000:
                connection.executemany("INSERT INTO links VALUES (?,?,?,?,?,?,?,?,?)", batch)
                batch.clear()
    if batch:
        connection.executemany("INSERT INTO links VALUES (?,?,?,?,?,?,?,?,?)", batch)
    if connection is not None:
        connection.execute("CREATE INDEX groups_order ON links (outer_fold,gene_id,method,endpoint,input_representation,score_unit)")
        connection.commit()
    after = entry["path"].stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source links changed during annotation: {entry['path']}")
    return selected, count


def families_for_motif(motif: str, metadata: dict) -> list[str]:
    if motif not in metadata:
        raise ValueError(f"Matched motif missing from frozen JASPAR database: {motif}")
    return sorted(set(filter(None, metadata[motif]["motif_family"].split(";")))) or [""]


def grouped_rows(connection: sqlite3.Connection, matches: dict, metadata: dict):
    cursor = connection.execute("SELECT outer_fold,gene_id,method,endpoint,input_representation,score_unit,peak_id,score,evidence FROM links ORDER BY outer_fold,gene_id,method,endpoint,input_representation,score_unit")
    present = compatible = 0
    for key, links in itertools.groupby(cursor, key=lambda row: tuple(row[:6])):
        families = {}
        for link in links:
            peak, score, evidence = link[6:]
            present += int(peak in matches)
            motif_ids = matches.get(peak, set())
            compatible += bool(motif_ids)
            for motif in motif_ids:
                for family in families_for_motif(motif, metadata):
                    summary = families.setdefault(family, dict(peaks=set(), motifs=set(), states=set(), maximum=-math.inf))
                    summary["peaks"].add(peak)
                    summary["motifs"].add(motif)
                    summary["states"].add(evidence)
                    summary["maximum"] = max(summary["maximum"], score)
        for family, summary in sorted(families.items()):
            row = dict(zip(GROUP, key))
            motif_ids = sorted(summary["motifs"])
            row.update(motif_family=family, family_annotation_state="recorded" if family else "unavailable",
                distinct_supporting_peak_count=len(summary["peaks"]), supporting_peak_ids=json.dumps(sorted(summary["peaks"]), separators=(",", ":")),
                max_peak_link_score=summary["maximum"], motif_ids=json.dumps(motif_ids, separators=(",", ":")),
                motif_names=json.dumps(sorted({metadata[motif]["motif_name"] for motif in motif_ids}), separators=(",", ":")),
                motif_classes=json.dumps(sorted({part for motif in motif_ids for part in metadata[motif]["motif_class"].split(";") if part}), separators=(",", ":")),
                link_evidence_states=json.dumps(sorted(summary["states"]), separators=(",", ":")))
            yield row
    # Final counts are returned as a terminal record, never a motif/gene row.
    yield {"_coverage": (present, compatible)}


def annotate(networks: list[Path], motif_root: Path, output: Path, compression_level: int = 9) -> None:
    if type(compression_level) is not int or not 0 <= compression_level <= 9:
        raise ValueError("compression_level must be an integer from 0 through 9")
    if output.exists():
        raise FileExistsError(f"Refusing overwrite: {output}")
    entries, sources = source_contexts(networks)
    destination = output.resolve()
    for source in [motif_root.resolve()] + [network.resolve() for network in networks]:
        if destination == source or source in destination.parents:
            raise ValueError(f"Output must not modify completed network or frozen motif inputs: {source}")
    metadata, record = context.jaspar_metadata(motif_root)
    output.mkdir(parents=True, exist_ok=False)
    entries.sort(key=lambda entry: (entry["cohort"], entry["lineage"], entry["network_root"]))
    with gzip.open(output / "motif_family_gene.tsv.gz", "wt", newline="", compresslevel=compression_level) as target, (output / "coverage.tsv").open("x", newline="") as coverage_target:
        writer = csv.DictWriter(target, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        coverage_writer = csv.DictWriter(coverage_target, fieldnames=COVERAGE, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        coverage_writer.writeheader()
        for _, same_context in itertools.groupby(entries, key=lambda entry: (entry["cohort"], entry["lineage"])):
            same_context = list(same_context)  # Only small context descriptors, never link rows.
            selected, staged = set(), []
            cohort, lineage = same_context[0]["cohort"], same_context[0]["lineage"]
            directory = motif_root / "motif_matches" / cohort / lineage
            whole_native_unscanned = cohort == "GSE296875"
            available = not whole_native_unscanned and all((directory / name).is_file() for name in ("matches.mtx", "peaks.tsv", "motifs.tsv"))
            with tempfile.TemporaryDirectory(prefix="family_sort_", dir=output) as temporary:
                try:
                    for index, entry in enumerate(same_context):
                        connection = sqlite3.connect(str(Path(temporary) / f"source_{index}.sqlite")) if available else None
                        peaks, n_links = stage_links(entry, connection)
                        selected.update(peaks)
                        staged.append((entry, connection, n_links, len(peaks)))
                    matches = context.selected_motif_matches(directory, selected) if available and selected else {}
                    for entry, connection, n_links, n_peaks in staged:
                        common = {key: entry[key] for key in COVERAGE[:7]}
                        present = compatible = n_rows = 0
                        if available:
                            for row in grouped_rows(connection, matches, metadata):
                                if "_coverage" in row:
                                    present, compatible = row["_coverage"]
                                    continue
                                writer.writerow(dict(**common, **row, endpoint_recorded=entry["endpoint_recorded"],
                                    evidence_state="motif_compatible_accessibility_link", motif_not_unique_tf=True,
                                    tf_activity_inferred=False, tf_binding_inferred=False, causal_regulation_inferred=False,
                                    jaspar_package=record["jaspar_package"], jaspar_database_sha256=record["database_sha256"]))
                                n_rows += 1
                        coverage_writer.writerow(dict(**common, endpoint_recorded=entry["endpoint_recorded"],
                            status="annotated" if available else "untestable",
                            reason="" if available else "whole_native_peak_universe_motifscan_not_run_125_peak_pilot_not_global" if whole_native_unscanned else "native_peak_motifscan_not_run",
                            n_link_records=n_links, n_unique_selected_peak_ids=n_peaks,
                            n_exact_peak_link_records=present if available else "", n_links_with_motif=compatible if available else "",
                            n_links_absent_from_motif_universe=n_links-present if available else "", n_family_gene_rows=n_rows,
                            peak_matching="exact_ID_only_no_overlap_borrowing"))
                finally:
                    for _, connection, _, _ in staged:
                        if connection is not None:
                            connection.close()
    for source in sources:
        validated_completion(Path(source["network_root"]))
    with (output / "completion.json").open("x") as handle:
        json.dump(dict(status="completed", scientific_adoption=False, sources=sources,
            motif_root=str(motif_root.resolve()), jaspar_source=record,
            interpretation="TF_family_motif_compatibility_only_no_TF_activity_binding_or_causality",
            exact_ID_only=True, continuous_score="maximum_source_peak_link_score_in_unchanged_method_units",
            endpoint_policy="preserve_recorded_endpoint_else_blank_without_inference",
            duplicate_peak_policy="distinct_peak_IDs_within_gene_fold_method_endpoint_family",
            memory_model="disk_sorted_unexpanded_links_one_context_motif_map_and_one_group_in_memory",
            GSE296875_whole_native_motif_universe="unscanned_untestable_125_peak_pilot_not_global"), handle, indent=2)
        handle.write("\n")


def self_test() -> None:
    """Compare against original Cartesian helper using tiny synthetic sources."""
    with tempfile.TemporaryDirectory(prefix="motif_family_fixture_") as temporary:
        root = Path(temporary)
        motif_root = root / "motifs"
        motif_root.mkdir()
        database = root / "jaspar.sqlite"
        with sqlite3.connect(database) as connection:
            connection.executescript("CREATE TABLE MATRIX (ID INTEGER, BASE_ID TEXT, VERSION INTEGER, NAME TEXT); CREATE TABLE MATRIX_ANNOTATION (ID INTEGER,TAG TEXT,VAL TEXT);")
            connection.executemany("INSERT INTO MATRIX VALUES (?,?,?,?)", [(1,"MA1",1,"TF1"),(2,"MA2",1,"TF2"),(3,"MA3",1,"TF3")])
            connection.executemany("INSERT INTO MATRIX_ANNOTATION VALUES (?,?,?)", [(1,"family","A"),(2,"family","A"),(3,"family","B")])
        def table(path, fields, records):
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader(); writer.writerows(records)
        table(motif_root / "jaspar_database.tsv", ["database_path","database_sha256","jaspar_package"], [dict(database_path=str(database),database_sha256=context.digest(database),jaspar_package="synthetic")])
        match_dir = motif_root / "motif_matches/GSE244832/hepatocyte"
        table(match_dir / "peaks.tsv", ["peak_coordinate"], [dict(peak_coordinate=p) for p in ("p1","p2","p3")])
        table(match_dir / "motifs.tsv", ["motif_id"], [dict(motif_id=p) for p in ("MA1.1","MA2.1","MA3.1")])
        (match_dir / "matches.mtx").write_text("%%MatrixMarket matrix coordinate pattern general\n3 3 4\n1 1\n1 2\n2 2\n2 3\n")
        links = [dict(outer_fold="0",gene_id="G",peak_id=p,method="m",rank="1",score=str(s),score_unit="continuous_native",evidence_state="candidate",endpoint=e)
                 for p,s,e in [("p1",.2,"promoter"),("p1",.7,"rna"),("p2",.1,"promoter"),("p1",.3,"promoter"),("p3",.8,"rna"),("p1_alias",.9,"promoter")]]
        def network(path, records, cohort="GSE244832"):
            path.mkdir()
            (path / "completion.json").write_text(json.dumps(dict(status="completed",pilot=True,interpretation="synthetic_technical_pilot")))
            table(path / f"{cohort}_hepatocyte/links.tsv", list(links[0]), records)
        source = root / "network"
        network(source, links)
        output = root / "families"
        annotate([source], motif_root, output)
        with gzip.open(output / "motif_family_gene.tsv.gz", "rt") as handle:
            actual = list(csv.DictReader(handle, delimiter="\t"))
        expected = {}
        for endpoint in ("promoter", "rna"):
            clone = root / endpoint
            network(clone, [row for row in links if row["endpoint"] == endpoint])
            cartesian = root / f"cartesian_{endpoint}"
            context.annotate([clone], motif_root, cartesian)
            with gzip.open(cartesian / "motif_peak_gene.tsv.gz", "rt") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    key = (endpoint,row["motif_family"])
                    record = expected.setdefault(key,dict(peaks=set(),motifs=set(),names=set(),maximum=-math.inf))
                    record["peaks"].add(row["peak_id"]); record["motifs"].add(row["motif_id"]); record["names"].add(row["motif_name"])
                    record["maximum"] = max(record["maximum"],float(row["score"]))
        assert len(actual) == len(expected) == 3
        for row in actual:
            reference = expected[(row["endpoint"],row["motif_family"])]
            assert set(json.loads(row["supporting_peak_ids"])) == reference["peaks"]
            assert int(row["distinct_supporting_peak_count"]) == len(reference["peaks"])
            assert set(json.loads(row["motif_ids"])) == reference["motifs"]
            assert set(json.loads(row["motif_names"])) == reference["names"]
            assert float(row["max_peak_link_score"]) == reference["maximum"]
            assert row["tf_activity_inferred"] == row["tf_binding_inferred"] == row["causal_regulation_inferred"] == "False"
        no_endpoint = root / "no_endpoint"
        no_endpoint.mkdir()
        (no_endpoint / "completion.json").write_text(json.dumps(dict(status="completed",pilot=False)))
        without_endpoint = [{key:value for key,value in row.items() if key != "endpoint"} for row in links]
        table(no_endpoint / "GSE244832_hepatocyte/links.tsv", list(without_endpoint[0]), without_endpoint)
        annotate([no_endpoint], motif_root, root / "no_endpoint_families")
        with gzip.open(root / "no_endpoint_families/motif_family_gene.tsv.gz", "rt") as handle:
            unrecorded = list(csv.DictReader(handle, delimiter="\t"))
        assert all(row["endpoint"] == "" and row["endpoint_recorded"] == "False" for row in unrecorded)
        assert all(row["source_scope"] == "development_candidate" for row in unrecorded)
        paired = root / "paired"
        network(paired, links, "GSE296875")
        shutil.copytree(match_dir, motif_root / "motif_matches/GSE296875/hepatocyte")
        annotate([paired], motif_root, root / "paired_families")
        paired_coverage = context.rows(root / "paired_families/coverage.tsv")
        assert paired_coverage[0]["status"] == "untestable" and int(paired_coverage[0]["n_family_gene_rows"]) == 0
        (source / "INVALIDATED_FIXTURE.json").write_text("{}")
        try:
            annotate([source], motif_root, root / "invalidated_output")
        except ValueError:
            pass
        else:
            raise AssertionError("Invalidated source accepted")
        assert not (root / "invalidated_output").exists()
        print("PASS: Cartesian equivalence; endpoint-preserving duplicate-peak handling; exact-ID-only; G296875 pilot exclusion; invalidated-source rejection")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, action="append")
    parser.add_argument("--motif-root", type=Path, default=context.DEFAULT_MOTIF_ROOT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--compression-level", type=int, default=9)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname():
        raise RuntimeError("Run motif annotation or fixtures in a compute allocation")
    if args.self_test:
        self_test()
    elif not args.network or args.output is None:
        parser.error("--network (repeatable) and --output are required")
    else:
        annotate(args.network, args.motif_root, args.output, args.compression_level)


if __name__ == "__main__":
    main()
