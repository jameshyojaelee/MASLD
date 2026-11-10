#!/usr/bin/env python3
"""Independent context annotation and bounded canonical streaming merge."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import resource
import shutil
import time

import annotate_motif_families as family


def context_key(entry):
    return tuple(str(entry[key]) for key in ("cohort", "lineage", "network_root"))


def scientific_completion(completion):
    return {key: value for key, value in completion.items() if key != "sources"}


def merge(networks, motif_root, parts, output, compression_level=9):
    if output.exists():
        raise FileExistsError(f"Refusing overwrite: {output}")
    entries, sources = family.source_contexts(networks)
    expected = {context_key(entry): entry for entry in entries}
    if len(expected) != len(entries):
        raise ValueError("Duplicate source context")
    indexed, completion_template = {}, None
    for part in parts:
        if any(part.glob("INVALIDATED*")):
            raise ValueError("Invalidated annotation part")
        completion = json.loads((part / "completion.json").read_text())
        if completion.get("status") != "completed" or completion.get("motif_root") != str(motif_root.resolve()):
            raise ValueError("Incomplete or different-motif annotation part")
        coverage = family.context.rows(part / "coverage.tsv")
        if len(coverage) != 1:
            raise ValueError("Annotation part must contain exactly one context")
        key = context_key(coverage[0])
        if key not in expected or key in indexed:
            raise ValueError("Unexpected or duplicate annotation context")
        source = next(source for source in sources if source["network_root"] == key[2])
        if completion.get("sources") != [source]:
            raise ValueError("Part source identity differs from completed network")
        scientific = scientific_completion(completion)
        if completion_template is None:
            completion_template = scientific
        elif scientific != completion_template:
            raise ValueError("Different scientific annotation metadata")
        indexed[key] = (part, coverage[0])
    if set(indexed) != set(expected):
        raise ValueError("Missing annotation contexts")
    # Validate the unchanged frozen JASPAR identity independently of part labels.
    _, jaspar = family.context.jaspar_metadata(motif_root)
    if completion_template["jaspar_source"] != jaspar:
        raise ValueError("Frozen motif database identity differs")
    output.mkdir(parents=True, exist_ok=False)
    table_hash = hashlib.sha256()
    part_receipts = []
    with gzip.open(output / "motif_family_gene.tsv.gz", "wb", compresslevel=compression_level) as target, (output / "coverage.tsv").open("x", newline="") as cover:
        writer = csv.DictWriter(cover, family.COVERAGE, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        expected_header = ("\t".join(family.FIELDS) + "\n").encode()
        target.write(expected_header)
        table_hash.update(expected_header)
        for key in sorted(expected):
            part, coverage = indexed[key]
            rows, last_byte = 0, b""
            with gzip.open(part / "motif_family_gene.tsv.gz", "rb") as source:
                if source.readline() != expected_header:
                    raise ValueError("Different annotation table schema")
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    target.write(chunk)
                    table_hash.update(chunk)
                    rows += chunk.count(b"\n")
                    last_byte = chunk[-1:]
            if (last_byte and last_byte != b"\n") or rows != int(coverage["n_family_gene_rows"]):
                raise ValueError("Incomplete annotation rows or changed coverage count")
            writer.writerow(coverage)
            part_receipts.append(dict(context=key, part=str(part.resolve()), rows=rows))
    for source in sources:
        family.validated_completion(Path(source["network_root"]))
    (output / "completion.json").write_text(json.dumps(dict(**completion_template, sources=sources), indent=2) + "\n")
    (output / "merge_receipt.json").write_text(json.dumps(dict(status="completed", canonical_context_order=True,
        decompressed_table_sha256=table_hash.hexdigest(), parts=part_receipts,
        compression_level=compression_level, scientific_adoption=False), indent=2) + "\n")


def run_context(job):
    network, motif_root, output, level = job
    started = time.monotonic()
    entries, _ = family.source_contexts([network])
    if len(entries) != 1:
        raise ValueError("Parallel annotation requires exactly one context per network root")
    path = entries[0]["path"]
    before = path.stat()
    source_hash = family.context.digest(path)
    family.annotate([network], motif_root, output, level)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Source links changed during independent annotation")
    receipt = dict(status="completed", network_root=str(network.resolve()), links_sha256=source_hash,
        links_bytes=before.st_size, links_mtime_ns=before.st_mtime_ns, worker_pid=os.getpid(),
        elapsed_seconds=time.monotonic()-started, max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (output / "execution.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def parallel(networks, motif_root, output, workers=4, compression_level=1):
    allocation = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    if not 1 <= workers <= min(4, allocation):
        raise ValueError("Workers must be1–4 and within the compute allocation")
    if not 0 <= compression_level <= 9:
        raise ValueError("Invalid gzip compression level")
    entries, _ = family.source_contexts(networks)
    if len(entries) != len(networks) or len({e["network_root"] for e in entries}) != len(networks):
        raise ValueError("Require exactly one context per network root")
    output.mkdir(parents=True, exist_ok=False)
    snapshots = output / "source_snapshot"
    snapshots.mkdir()
    hashes = {}
    for source in [Path(__file__), Path(family.__file__), Path(family.context.__file__)]:
        shutil.copyfile(source, snapshots / source.name)
        hashes[source.name] = family.context.digest(source)
    (output / "source_sha256.json").write_text(json.dumps(hashes, indent=2) + "\n")
    parts = output / "parts"
    parts.mkdir()
    jobs = [(Path(entry["network_root"]), motif_root, parts / f'{entry["cohort"]}_{entry["lineage"]}', compression_level)
            for entry in sorted(entries, key=context_key)]
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers) as executor:
        receipts = list(executor.map(run_context, jobs))
    merge(networks, motif_root, [job[2] for job in jobs], output / "merged", compression_level)
    (output / "execution.json").write_text(json.dumps(dict(status="completed", workers=workers,
        compression_level=compression_level, elapsed_seconds=time.monotonic()-started, contexts=receipts,
        scientific_adoption=False), indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, action="append", required=True)
    parser.add_argument("--motif-root", type=Path, default=family.context.DEFAULT_MOTIF_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--part", type=Path, action="append", help="Merge existing completed parts instead of launching workers")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--compression-level", type=int, default=1)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or "login" in family.socket.gethostname():
        raise RuntimeError("Compute allocation required")
    if args.part:
        merge(args.network, args.motif_root, args.part, args.output, args.compression_level)
    else:
        parallel(args.network, args.motif_root, args.output, args.workers, args.compression_level)


if __name__ == "__main__":
    main()
