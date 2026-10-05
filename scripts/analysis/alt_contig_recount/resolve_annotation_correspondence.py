#!/usr/bin/env python3
"""Resolve the counting annotation against GRCh38.p14, without symbol folding.

Every non-chromosome annotation receives a disposition. Unresolved records
prevent emission of the grouped counting annotation and any primary RNA fit.
Coordinates are zero-based half-open internally; the official GFF Target/Gap
alignment, rather than its placement bounding box, determines correspondence.
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import urllib.request

BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/"
STRUCTURE = BASE + "GCA_000001405.29_GRCh38.p14_assembly_structure/"
PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
ATTR = re.compile(r'(\w+) "([^"]*)"')
DISPOSITIONS = {"primary_gene", "equivalent", "distinct_locus", "alternative_only", "duplicated", "unresolved"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def union(intervals):
    out = []
    for start, end in sorted(set(intervals)):
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(end, out[-1][1]))
        else:
            out.append((start, end))
    return out


def length(intervals):
    return sum(end - start for start, end in union(intervals))


def overlap(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


@dataclass
class Gene:
    gene_id: str
    name: str
    biotype: str
    chrom: str
    strand: str
    start: int
    end: int
    hgnc: str = ""
    havana: str = ""
    exons: list = field(default_factory=list)


def read_genes(path):
    genes = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9 or c[2] not in {"gene", "exon"}:
                continue
            a = dict(ATTR.findall(c[8]))
            gid = a["gene_id"]
            if gid not in genes:
                genes[gid] = Gene(gid, a.get("gene_name", ""), a.get("gene_type", ""),
                                  c[0], c[6], int(c[3])-1, int(c[4]),
                                  a.get("hgnc_id", ""), a.get("havana_gene", ""))
            g = genes[gid]
            if (g.chrom, g.strand) != (c[0], c[6]):
                raise ValueError(f"Gene ID appears at different loci: {gid}")
            g.start, g.end = min(g.start, int(c[3])-1), max(g.end, int(c[4]))
            if c[2] == "exon":
                g.exons.append((int(c[3])-1, int(c[4])))
    for g in genes.values():
        g.exons = union(g.exons)
    return genes


def fetch(url, target, cache=None):
    if target.exists():
        raise ValueError(f"Refusing to overwrite {target}")
    cached = cache.get(url) if cache else None
    if cached:
        source = Path(cached["path"])
        if digest(source) != cached["sha256"]:
            raise ValueError(f"Assembly support cache hash mismatch: {source}")
        data = source.read_bytes()
    else:
        with urllib.request.urlopen(url, timeout=45) as response:
            data = response.read(5 * 1024 * 1024 + 1)
    if len(data) > 5 * 1024 * 1024:
        raise ValueError(f"Unexpectedly large assembly support file: {url}")
    target.write_bytes(data)
    return {"url": url, "path": str(target.resolve()), "sha256": digest(target), "bytes": len(data)}


def read_report(path):
    result = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.rstrip().split("\t")
            if len(c) < 10:
                continue
            row = {"name": c[0], "role": c[1], "genbank": c[4], "refseq": c[6], "unit": c[7], "ucsc": c[9]}
            for alias in (c[0], c[4], c[6], c[9]):
                if alias != "na":
                    result[alias] = row
    return result


def read_placements(path):
    with path.open() as handle:
        header = handle.readline().lstrip("#").rstrip().split("\t")
        return {r["alt_scaf_acc"]: r for r in csv.DictReader(handle, fieldnames=header, delimiter="\t")}


def alignment_blocks(path):
    """Return alt_start, alt_end, primary_start, primary_end, orientation."""
    blocks = []
    parent = None
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.rstrip().split("\t")
            if len(c) != 9 or c[2] != "match":
                continue
            a = dict(v.split("=", 1) for v in c[8].split(";") if "=" in v)
            target = a["Target"].split()
            strand_q = target[3] if len(target) == 4 else "+"
            p, q = int(c[3])-1, int(target[1])-1
            pend, qend = int(c[4]), int(target[2])
            pstep, qstep = (-1 if c[6] == "-" else 1), (-1 if strand_q == "-" else 1)
            if pstep == -1: p = pend
            if qstep == -1: q = qend
            gaps = re.findall(r"([MID])(\d+)", a.get("Gap", f"M{pend-int(c[3])+1}"))
            for op, number in gaps:
                n = int(number)
                if op == "M":
                    pn, qn = p+pstep*n, q+qstep*n
                    blocks.append((min(q, qn), max(q, qn), min(p, pn), max(p, pn), pstep*qstep))
                    p, q = pn, qn
                elif op == "D": p += pstep*n
                else: q += qstep*n
            if p != (pend if pstep == 1 else int(c[3])-1) or q != (qend if qstep == 1 else int(target[1])-1):
                raise ValueError(f"GFF Gap does not consume its coordinates: {path}")
            if parent is not None and parent != c[0]:
                raise ValueError(f"Mixed parent sequences: {path}")
            parent = c[0]
    return parent, blocks


def project(exons, blocks):
    mapped, covered, orientations = [], [], set()
    for s, e in exons:
        for qs, qe, ps, pe, ori in blocks:
            a, b = max(s, qs), min(e, qe)
            if a >= b: continue
            covered.append((a, b)); orientations.add(ori)
            mapped.append((ps+a-qs, ps+b-qs) if ori == 1 else (pe-(b-qs), pe-(a-qs)))
    return union(mapped), length(covered), orientations


def write_tsv(path, rows, fields):
    with path.open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-cache-manifest", type=Path)
    parser.add_argument("--prior-correspondence", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    sources = args.output / "sources"; sources.mkdir()
    manifest = []
    cache = {}
    if args.source_cache_manifest:
        with args.source_cache_manifest.open() as handle:
            cache = {r["url"]: r for r in csv.DictReader(handle, delimiter="\t")}
    report_path = sources / "assembly_report.txt"
    placement_path = sources / "all_alt_scaffold_placement.txt"
    manifest.append(fetch(BASE+"GCA_000001405.29_GRCh38.p14_assembly_report.txt", report_path, cache))
    manifest.append(fetch(STRUCTURE+"all_alt_scaffold_placement.txt", placement_path, cache))
    report, placements = read_report(report_path), read_placements(placement_path)
    print("Reading exact counting GTF", file=sys.stderr, flush=True)
    genes = read_genes(args.gtf)
    by_chrom = defaultdict(list)
    for g in genes.values():
        if g.chrom in PRIMARY: by_chrom[g.chrom].append(g)
    for items in by_chrom.values(): items.sort(key=lambda g: g.start)
    starts = {c: [g.start for g in items] for c, items in by_chrom.items()}
    mappings, mapping_errors = {}, {}
    for chrom in sorted({g.chrom for g in genes.values()} - PRIMARY):
        seq = report.get(chrom)
        if not seq or seq["unit"] == "Primary Assembly": continue
        placement = placements.get(seq["genbank"])
        if not placement or placement["parent_acc"] == "na":
            mapping_errors[chrom] = "assembly_has_no_coordinate_alignment"; continue
        unit = placement["alt_asm_name"]
        name = f'{seq["genbank"]}_{placement["parent_acc"]}.gff'
        try:
            manifest.append(fetch(STRUCTURE+f"{unit}/alt_scaffolds/alignments/{name}", sources/name, cache))
            parent, blocks = alignment_blocks(sources/name)
            parent_seq = report.get(parent)
            if not parent_seq or parent_seq["ucsc"] not in PRIMARY:
                mapping_errors[chrom] = "parent_is_not_a_primary_chromosome"
            else: mappings[chrom] = (parent_seq["ucsc"], blocks)
        except (OSError, ValueError) as exc:
            mapping_errors[chrom] = f"alignment_unavailable_or_invalid:{exc}"
    rows = []
    for gid, g in sorted(genes.items()):
        row = {"gene_id": gid, "gene_name": g.name, "gene_biotype": g.biotype, "chromosome": g.chrom,
               "strand": g.strand, "hgnc_id": g.hgnc, "exonic_bases": length(g.exons),
               "mapped_exonic_bases": "", "primary_gene_candidates": "", "counting_gene_id": gid,
               "unmapped_exonic_bases": "", "coordinate_alignment_state": "not_applicable",
               "alternate_biotype_differs": "false", "disposition": "unresolved", "reason": "", "evidence_source": ""}
        seq = report.get(g.chrom)
        if g.chrom in PRIMARY:
            row.update(disposition="primary_gene", reason="primary_chromosome_annotation", evidence_source="exact_counting_gtf")
        elif seq and seq["unit"] == "Primary Assembly":
            row.update(disposition="distinct_locus", reason="nonredundant_primary_assembly_scaffold", evidence_source="assembly_report.txt")
        elif g.chrom not in mappings:
            row["reason"] = mapping_errors.get(g.chrom, "sequence_not_in_assembly_report")
        elif not g.exons:
            row["reason"] = "annotation_has_no_exons"
        else:
            parent, blocks = mappings[g.chrom]
            mapped, covered, orientations = project(g.exons, blocks)
            row["mapped_exonic_bases"] = covered
            row["unmapped_exonic_bases"] = length(g.exons)-covered
            row["coordinate_alignment_state"] = "complete" if covered == length(g.exons) and len(orientations) == 1 else "partial_or_mixed"
            row["evidence_source"] = f'{seq["genbank"]}_{placements[seq["genbank"]]["parent_acc"]}.gff'
            if not covered or len(orientations) != 1:
                row["reason"] = "exonic_coordinate_correspondence_incomplete_or_mixed"
            else:
                ori = next(iter(orientations)); strand = g.strand if ori == 1 else ("-" if g.strand == "+" else "+")
                body_candidates = [p for p in by_chrom[parent][:bisect.bisect_left(starts[parent], mapped[-1][1])]
                                   if p.end > mapped[0][0] and p.strand == strand
                                   and overlap(mapped, [(p.start, p.end)])]
                candidates = [p for p in body_candidates if overlap(mapped, p.exons)]
                row["primary_gene_candidates"] = ";".join(p.gene_id for p in candidates)
                exact = [p for p in candidates if covered == length(g.exons) and mapped == p.exons and p.biotype == g.biotype]
                identified = [p for p in body_candidates if (g.hgnc and g.hgnc == p.hgnc)
                              or (g.havana and p.havana and g.havana.split(".")[0] == p.havana.split(".")[0])]
                selected = {p.gene_id: p for p in exact+identified}
                if len(selected) == 1:
                    p = next(iter(selected.values()))
                    row.update(disposition="equivalent", counting_gene_id=p.gene_id,
                               alternate_biotype_differs=str(p.biotype != g.biotype).lower(),
                               reason="mapped_unique_locus_and_authoritative_gene_identity" if p in identified else "exact_reciprocal_exon_correspondence")
                elif len(selected) > 1:
                    row.update(disposition="duplicated", reason="multiple_authoritative_primary_gene_correspondences")
                elif covered != length(g.exons):
                    row["reason"] = "partial_exon_map_without_authoritative_gene_identity"
                elif not candidates:
                    row.update(disposition="alternative_only", reason="fully_mapped_exons_have_no_same_strand_primary_gene_annotation")
                elif g.hgnc and all(p.hgnc and p.hgnc != g.hgnc for p in candidates):
                    row.update(disposition="distinct_locus", reason="mapped_locus_with_distinct_HGNC_identity")
                else:
                    row["reason"] = "coordinate_overlap_does_not_resolve_gene_identity"
        rows.append(row)
    # Two copies within one haplotype cannot both silently become one primary gene.
    groups = defaultdict(list)
    for r in rows:
        if r["disposition"] == "equivalent": groups[(r["chromosome"], r["counting_gene_id"])].append(r)
    for group in groups.values():
        if len(group) > 1:
            for r in group:
                r.update(disposition="duplicated", counting_gene_id=r["gene_id"], reason="multiple_gene_copies_within_one_haplotype")
    fields = list(rows[0])
    write_tsv(args.output/"annotation_correspondence.tsv", rows, fields)
    write_tsv(args.output/"source_manifest.tsv", manifest, ["url", "path", "sha256", "bytes"])
    counts = dict(Counter(r["disposition"] for r in rows))
    unresolved = [r for r in rows if r["disposition"] == "unresolved"]
    write_tsv(args.output/"unresolved_annotations.tsv", unresolved, fields)
    if args.prior_correspondence:
        with args.prior_correspondence.open() as handle:
            prior = {r["gene_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
        transitions = [{"gene_id": r["gene_id"], "previous_disposition": prior[r["gene_id"]]["disposition"],
                        "previous_reason": prior[r["gene_id"]]["reason"], "disposition": r["disposition"], "reason": r["reason"]}
                       for r in rows if r["gene_id"] in prior and
                       (r["disposition"], r["reason"]) != (prior[r["gene_id"]]["disposition"], prior[r["gene_id"]]["reason"])]
        write_tsv(args.output/"reason_transitions.tsv", transitions,
                  ["gene_id", "previous_disposition", "previous_reason", "disposition", "reason"])
    summary = {"status": "resolved" if not unresolved else "held_unresolved_annotations", "candidate_not_adopted": True,
               "assembly": "GCA_000001405.29_GRCh38.p14", "gtf": str(args.gtf.resolve()), "gtf_sha256": digest(args.gtf),
               "n_annotations": len(rows), "dispositions": counts, "n_unresolved": len(unresolved),
               "annotation_resolution_complete": not unresolved, "n_unresolved_annotations": len(unresolved),
               "counts_complete": False, "counts_sha256": None,
               "primary_refit_allowed": False, "symbol_only_merging_used": False}
    (args.output/"resolution_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    if not unresolved:
        group_ids = {r["gene_id"]: r["counting_gene_id"] for r in rows}
        with gzip.open(args.gtf, "rt") as src, gzip.open(args.output/"grouped_counting_annotation.gtf.gz", "wt") as dest:
            for line in src:
                if line.startswith("#"): dest.write(line); continue
                c = line.rstrip("\n").split("\t")
                if len(c) != 9 or c[2] != "exon": continue
                gid = dict(ATTR.findall(c[8]))["gene_id"]
                c[8] = re.sub(r'gene_id "[^"]+"', f'gene_id "{group_ids[gid]}"', c[8])
                dest.write("\t".join(c)+"\n")
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__": main()
