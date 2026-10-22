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
DISPOSITIONS = {"primary_gene", "equivalent", "distinct_locus", "alternative_only", "duplicated", "excluded_reference", "unresolved",
                "held_untestable"}


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


def verified_exon_segments(exons, blocks):
    """ALT coordinates supported by alignment matches; insertion bases excluded."""
    return union([(max(s, qs), min(e, qe)) for s, e in exons
                  for qs, qe, _, _, _ in blocks if max(s, qs) < min(e, qe)])


def write_tsv(path, rows, fields):
    with path.open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def read_allele_groups(root, genes, report):
    """Exact-version Ensembl groups, independently checked against the GTF."""
    def rows(name):
        with gzip.open(root/name, "rt") as handle:
            return [line.rstrip("\n").split("\t") for line in handle]
    meta = {r[2]: r[3] for r in rows("ensembl115_meta.txt.gz")}
    if meta.get("gencode.version") != "GENCODE 49" or meta.get("assembly.accession") != "GCA_000001405.29":
        raise ValueError("Allele groups must use the exact GENCODE49/GRCh38.p14 assembly")
    sequences = {r[0]: r[1] for r in rows("ensembl115_seq_region.txt.gz")}
    alleles = rows("alt_allele.txt.gz")
    allele_core_ids = {r[2] for r in alleles}
    core = {}
    for r in rows("gene.txt.gz"):
        if r[0] not in allele_core_ids: continue
        gid = f"{r[12]}.{r[13]}"
        core[r[0]] = gid
        if gid not in genes: continue
        g = genes[gid]; name = sequences[r[3]]
        seq = report.get(name) or report.get("chrM" if name == "MT" else "chr"+name)
        if not seq or seq != report.get(g.chrom) or (int(r[4])-1, int(r[5]), "+" if r[6] == "1" else "-") != (g.start, g.end, g.strand):
            raise ValueError(f"Ensembl/GTF exact-version locus mismatch: {gid}")
    attrs = defaultdict(set)
    for aid, flag in rows("alt_allele_attrib.txt.gz"): attrs[aid].add(flag)
    groups, memberships = defaultdict(list), defaultdict(list)
    for aid, group, core_id in alleles:
        gid = core.get(core_id)
        if gid:
            member = (gid, attrs[aid])
            groups[group].append(member); memberships[gid].append(group)
    return groups, memberships


def allele_identity_candidates(gene, body_candidates, genes, groups, memberships):
    """Membership alone never establishes equivalence or absence of a gene."""
    candidates = set()
    for group in memberships.get(gene.gene_id, []):
        members = groups[group]
        if any(gid not in genes for gid, _ in members):
            continue
        primary = [gid for gid, _ in members if genes[gid].chrom in PRIMARY]
        same_scaffold = [gid for gid, _ in members if genes[gid].chrom == gene.chrom]
        flags = set().union(*(flag for gid, flag in members if gid == gene.gene_id))
        if len(primary) != 1 or len(same_scaffold) != 1:
            continue
        if flags & {"IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY", "IN_BROKEN_ASSEMBLY_REGION", "IS_PAR"}:
            continue
        p = genes[primary[0]]
        if p not in body_candidates or (gene.hgnc and p.hgnc and gene.hgnc != p.hgnc):
            continue
        candidates.add(p.gene_id)
    return [genes[gid] for gid in sorted(candidates)]


def local_source_records(root, names):
    """Bind downloaded source bytes to their recorded official URL/hash."""
    records = {}
    for path in root.glob("*manifest.json"):
        value = json.loads(path.read_text())
        if isinstance(value, list):
            for row in value:
                if row.get("status", "downloaded") == "downloaded" and row.get("path") and row.get("sha256"):
                    records[Path(row["path"]).name] = row
    result = []
    for name in names:
        path = root/name; record = records.get(name)
        if not record or digest(path) != record["sha256"]:
            raise ValueError(f"Official source missing or altered: {path}")
        result.append({"url": record["url"], "path": str(path.resolve()), "sha256": record["sha256"], "bytes": path.stat().st_size})
    return result


def added_patch_component(gene, seq, placement, source_root):
    """Positive AGP evidence that a patch component replaces only primary N gaps."""
    if seq["unit"] != "PATCHES" or not placement or placement["parent_acc"] == "na":
        return False
    patch_agp = source_root/"alt.scaf.agp.gz"
    parent_agp = source_root/f"chr{placement['parent_name']}.agp.gz"
    if not patch_agp.exists() or not parent_agp.exists(): return False
    with gzip.open(patch_agp, "rt") as handle:
        components = [line.rstrip().split("\t") for line in handle if not line.startswith("#")]
    components = [r for r in components if len(r) == 9 and r[0] == seq["genbank"] and r[4] in {"A", "D", "F", "G", "O", "P", "W"}
                  and int(r[1])-1 <= gene.start and int(r[2]) >= gene.end]
    if len(components) != 1: return False
    start, end = int(placement["parent_start"])-1, int(placement["parent_stop"])
    with gzip.open(parent_agp, "rt") as handle:
        gaps = [line.rstrip().split("\t") for line in handle if not line.startswith("#")]
    return any(len(r) == 9 and r[4] == "N" and r[6] == "short_arm"
               and int(r[1])-1 <= start and int(r[2]) >= end for r in gaps)


def reference_exclusion(gene, report, intervals):
    """Exclude complete bodies/exon unions in official invalid reference regions.

    BED and gene coordinates are zero-based half-open. An intron-only overlap
    does not remove valid expression features. Assembly aliases are required.
    """
    seq = report.get(gene.chrom)
    if not seq or not gene.exons:
        return None
    for kind in ("contamination", "false_duplication"):
        matched = union([(start, end) for chrom, start, end, label in intervals
                         if label == kind and report.get(chrom) == seq])
        if (overlap([(gene.start, gene.end)], matched) == gene.end-gene.start
                and overlap(gene.exons, matched) == length(gene.exons)):
            return kind
    return None


def hold_unresolved(rows):
    """Keep unresolved annotations as their own counting targets and mark them not testable.

    A read shared between a held annotation and another gene then has two targets and is
    rejected; the held gene itself is reported but never enters the tested matrix.
    """
    held = [r for r in rows if r["disposition"] == "unresolved"]
    for r in held:
        r.update(disposition="held_untestable", counting_gene_id=r["gene_id"])
    return held


def write_grouped_annotation(gtf, output, rows, mappings):
    group_ids = {r["gene_id"]: r["counting_gene_id"] for r in rows}
    excluded = {r["gene_id"] for r in rows if r["disposition"] == "excluded_reference"}
    partial_equivalents = {r["gene_id"] for r in rows if r["disposition"] == "equivalent"
                          and r["coordinate_alignment_state"] != "complete"}
    with gzip.open(gtf, "rt") as src, gzip.open(output, "wt") as dest:
        for line in src:
            if line.startswith("#"): dest.write(line); continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9 or c[2] != "exon": continue
            gid = dict(ATTR.findall(c[8]))["gene_id"]
            if gid in excluded:
                continue
            if not group_ids[gid]:
                raise ValueError(f"Included annotation lacks a counting target: {gid}")
            c[8] = re.sub(r'gene_id "[^"]+"', f'gene_id "{group_ids[gid]}"', c[8])
            segments = [(int(c[3])-1, int(c[4]))]
            if gid in partial_equivalents:
                segments = verified_exon_segments(segments, mappings[c[0]][1])
            for start, end in segments:
                c[3], c[4] = str(start+1), str(end)
                dest.write("\t".join(c)+"\n")


def exact_transcript_correspondence(source_exons, target_exons, source_strand, target_strand, blocks):
    mapped, covered, orientations = project(source_exons, blocks)
    if len(orientations) != 1 or covered != length(source_exons):
        return False
    strand = source_strand if orientations == {1} else ("-" if source_strand == "+" else "+")
    return strand == target_strand and mapped == union(target_exons)


def resolve_exact_RNA_transcripts(gtf, metadata, genes, mappings, rows, allele_groups, memberships):
    """Unique curated RNA ID plus complete exact isoform coordinates, never xref alone."""
    accessions = defaultdict(set)
    with gzip.open(metadata, "rt") as handle:
        for line in handle:
            c = line.rstrip().split("\t")
            if len(c) >= 2 and re.fullmatch(r"N[MR]_\d+\.\d+", c[1]):
                accessions[c[0]].add(c[1])
    unresolved = {r["gene_id"]: r for r in rows if r["disposition"] == "unresolved"}
    wanted = set(unresolved)
    for r in unresolved.values(): wanted.update(filter(None, r["primary_gene_candidates"].split(";")))
    models, primary_accessions = {}, defaultdict(set)
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"): continue
            c = line.rstrip().split("\t")
            if len(c) != 9 or c[2] not in {"transcript", "exon"}: continue
            a = dict(ATTR.findall(c[8])); tid = a.get("transcript_id")
            if tid not in accessions: continue
            gid = a["gene_id"]
            if c[2] == "transcript" and c[0] in PRIMARY:
                for acc in accessions[tid]: primary_accessions[acc].add(gid)
            if gid not in wanted: continue
            model = models.setdefault((gid, tid), {"strand": c[6], "exons": []})
            if c[2] == "exon": model["exons"].append((int(c[3])-1, int(c[4])))
    by_gene = defaultdict(list)
    for (gid, tid), model in models.items():
        model["exons"] = union(model["exons"])
        by_gene[gid].append((tid, model))
    dispositions = {r["gene_id"]: r["disposition"] for r in rows}
    for gid, row in unresolved.items():
        g = genes[gid]
        if g.chrom not in mappings: continue
        candidates, proof = set(), []
        group_primary = set()
        flags = set()
        for group in memberships.get(gid, []):
            for member, attributes in allele_groups[group]:
                if member in genes and genes[member].chrom in PRIMARY: group_primary.add(member)
                if member == gid: flags.update(attributes)
        if flags & {"IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY", "IN_BROKEN_ASSEMBLY_REGION", "IS_PAR"}: continue
        for primary_id in filter(None, row["primary_gene_candidates"].split(";")):
            p = genes[primary_id]
            if dispositions[primary_id] != "primary_gene" or p.biotype != g.biotype: continue
            if g.hgnc and p.hgnc and g.hgnc != p.hgnc: continue
            if group_primary and primary_id not in group_primary: continue
            for tid, model in by_gene[gid]:
                for ptid, primary_model in by_gene[primary_id]:
                    shared = {acc for acc in accessions[tid] & accessions[ptid]
                              if primary_accessions[acc] == {primary_id}}
                    if shared and exact_transcript_correspondence(model["exons"], primary_model["exons"],
                        model["strand"], primary_model["strand"], mappings[g.chrom][1]):
                        candidates.add(primary_id)
                        proof.append(f"{tid}|{ptid}|{','.join(sorted(shared))}")
        if len(candidates) == 1:
            row.update(disposition="equivalent", counting_gene_id=next(iter(candidates)),
                       reason="unique_curated_RNA_accession_and_complete_exact_transcript_coordinate_correspondence",
                       transcript_identity_evidence=";".join(sorted(set(proof))))
        elif len(candidates) > 1:
            row.update(reason="conflicting_exact_transcript_gene_targets")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-cache-manifest", type=Path)
    parser.add_argument("--prior-correspondence", type=Path)
    parser.add_argument("--ensembl-source-root", type=Path)
    parser.add_argument("--added-patch-source-root", type=Path)
    parser.add_argument("--grc-exclusions", type=Path)
    parser.add_argument("--refseq-transcript-metadata", type=Path)
    parser.add_argument("--hold-unresolved-as-untestable", action="store_true",
                        help="Authorized 2026-10-07: count the resolved annotations and report the rest as not testable")
    parser.add_argument("--apply-grc-reference-exclusions", action="store_true",
                        help="Apply the authorized complete-body/exon contamination and false-duplication exclusions")
    args = parser.parse_args()
    if args.apply_grc_reference_exclusions and not args.grc_exclusions:
        parser.error("Reference exclusions require the official GRC BED")
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
    allele_groups, allele_memberships = ({}, {})
    if args.ensembl_source_root:
        manifest.extend(local_source_records(args.ensembl_source_root,
            ["ensembl115_meta.txt.gz", "ensembl115_seq_region.txt.gz", "gene.txt.gz", "alt_allele.txt.gz", "alt_allele_attrib.txt.gz", "homo_sapiens_core_115_38.sql.gz"]))
        allele_groups, allele_memberships = read_allele_groups(args.ensembl_source_root, genes, report)
    if args.added_patch_source_root:
        manifest.extend(local_source_records(args.added_patch_source_root,
            ["alt.scaf.agp.gz", "chr13.agp.gz", "chr14.agp.gz", "chr15.agp.gz", "chr22.agp.gz"]))
    if args.refseq_transcript_metadata:
        manifest.extend(local_source_records(args.refseq_transcript_metadata.parent, [args.refseq_transcript_metadata.name]))
    contamination = []
    reference_intervals = []
    if args.grc_exclusions:
        manifest.extend(local_source_records(args.grc_exclusions.parent, [args.grc_exclusions.name]))
        for line in args.grc_exclusions.read_text().splitlines():
            c = line.split("\t")
            if len(c) >= 4 and c[3] == "contamination":
                contamination.append((c[0], int(c[1]), int(c[2])))
            if len(c) >= 4 and c[3] in {"contamination", "false_duplication"}:
                reference_intervals.append((c[0], int(c[1]), int(c[2]), c[3]))
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
            # An official alignment can exist despite a missing placement row.
            matches = [url for url in cache if url.startswith(STRUCTURE)
                       and url.rsplit("/", 1)[-1].startswith(seq["genbank"]+"_") and url.endswith(".gff")]
            if len(matches) != 1:
                mapping_errors[chrom] = "assembly_has_no_coordinate_alignment"; continue
            alignment_url = matches[0]; name = alignment_url.rsplit("/", 1)[-1]
        else:
            unit = placement["alt_asm_name"]
            name = f'{seq["genbank"]}_{placement["parent_acc"]}.gff'
            alignment_url = STRUCTURE+f"{unit}/alt_scaffolds/alignments/{name}"
        try:
            manifest.append(fetch(alignment_url, sources/name, cache))
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
        row["ensembl_allele_group_ids"] = ";".join(allele_memberships.get(gid, []))
        row["transcript_identity_evidence"] = ""
        seq = report.get(g.chrom)
        excluded_reason = reference_exclusion(g, report, reference_intervals) if args.apply_grc_reference_exclusions else None
        if excluded_reason:
            row.update(disposition="excluded_reference", counting_gene_id="",
                       reason=f"GRC_{excluded_reason}_complete_gene_body_and_exons",
                       evidence_source="GRC_exclusions.bed")
        elif any(report.get(g.chrom) == report.get(chrom) and report.get(chrom)
               and g.start < end and g.end > start for chrom, start, end in contamination):
            row.update(reason="GRC_contamination_reference_policy_pending", evidence_source="GRC_exclusions.bed")
        elif g.chrom in PRIMARY:
            row.update(disposition="primary_gene", reason="primary_chromosome_annotation", evidence_source="exact_counting_gtf")
        elif seq and seq["unit"] == "Primary Assembly":
            row.update(disposition="distinct_locus", reason="nonredundant_primary_assembly_scaffold", evidence_source="assembly_report.txt")
        elif g.chrom not in mappings:
            if seq and args.added_patch_source_root and added_patch_component(g, seq, placements.get(seq["genbank"]), args.added_patch_source_root):
                row.update(disposition="alternative_only", reason="added_patch_component_wholly_replaces_primary_N_short_arm_gap",
                           evidence_source="official_patch_and_primary_AGP")
            else:
                row["reason"] = mapping_errors.get(g.chrom, "sequence_not_in_assembly_report")
        elif not g.exons:
            row["reason"] = "annotation_has_no_exons"
        else:
            parent, blocks = mappings[g.chrom]
            mapped, covered, orientations = project(g.exons, blocks)
            row["mapped_exonic_bases"] = covered
            row["unmapped_exonic_bases"] = length(g.exons)-covered
            row["coordinate_alignment_state"] = "complete" if covered == length(g.exons) and len(orientations) == 1 else "partial_or_mixed"
            row["evidence_source"] = next(r["path"] for r in manifest if Path(r["path"]).name.startswith(seq["genbank"]+"_") and r["path"].endswith(".gff"))
            mapped_body, body_covered, body_orientations = project([(g.start, g.end)], blocks)
            allele_body_candidates = []
            if body_covered and len(body_orientations) == 1:
                body_ori = next(iter(body_orientations))
                body_strand = g.strand if body_ori == 1 else ("-" if g.strand == "+" else "+")
                allele_body_candidates = [p for p in by_chrom[parent][:bisect.bisect_left(starts[parent], mapped_body[-1][1])]
                    if p.end > mapped_body[0][0] and p.strand == body_strand
                    and overlap(mapped_body, [(p.start, p.end)])]
            if not covered or len(orientations) != 1:
                identified = allele_identity_candidates(g, allele_body_candidates, genes, allele_groups, allele_memberships)
                if not covered and not orientations and len(identified) == 1:
                    p = identified[0]
                    row.update(disposition="equivalent", counting_gene_id=p.gene_id,
                        alternate_biotype_differs=str(p.biotype != g.biotype).lower(),
                        reason="mapped_gene_body_and_exact_Ensembl_allele_identity_no_rescuable_ALT_exons")
                else:
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
                coordinate_selected = {p.gene_id for p in exact+identified}
                allele_identified = allele_identity_candidates(g, allele_body_candidates, genes, allele_groups, allele_memberships)
                identified += allele_identified
                selected = {p.gene_id: p for p in exact+identified}
                if len(selected) == 1:
                    p = next(iter(selected.values()))
                    row.update(disposition="equivalent", counting_gene_id=p.gene_id,
                               alternate_biotype_differs=str(p.biotype != g.biotype).lower(),
                               reason="mapped_unique_locus_and_authoritative_gene_identity" if p in identified else "exact_reciprocal_exon_correspondence")
                elif len(selected) > 1:
                    if allele_identified and len(coordinate_selected) <= 1:
                        row["reason"] = "conflicting_allele_group_and_coordinate_or_gene_identity"
                    else:
                        row.update(disposition="duplicated", reason="multiple_authoritative_primary_gene_correspondences")
                elif covered != length(g.exons):
                    row["reason"] = "partial_exon_map_without_authoritative_gene_identity"
                elif not body_candidates:
                    row.update(disposition="alternative_only", reason="fully_mapped_exons_have_no_same_strand_primary_gene_body_annotation")
                elif not candidates:
                    row["reason"] = "mapped_exons_within_primary_gene_body_without_resolved_gene_identity"
                elif g.hgnc and all(p.hgnc and p.hgnc != g.hgnc for p in candidates):
                    row.update(disposition="distinct_locus", reason="mapped_locus_with_distinct_HGNC_identity")
                else:
                    row["reason"] = "coordinate_overlap_does_not_resolve_gene_identity"
        rows.append(row)
    if args.refseq_transcript_metadata:
        resolve_exact_RNA_transcripts(args.gtf, args.refseq_transcript_metadata, genes, mappings, rows,
                                     allele_groups, allele_memberships)
    excluded_targets = {r["gene_id"] for r in rows if r["disposition"] == "excluded_reference"}
    for r in rows:
        if r["disposition"] == "equivalent" and r["counting_gene_id"] in excluded_targets:
            r.update(disposition="unresolved", reason="equivalence_targets_excluded_reference_annotation")
    # Two copies within one haplotype cannot both silently become one primary gene.
    groups = defaultdict(list)
    for r in rows:
        if r["disposition"] == "equivalent": groups[(r["chromosome"], r["counting_gene_id"])].append(r)
    for group in groups.values():
        if len(group) > 1:
            for r in group:
                r.update(disposition="duplicated", counting_gene_id=r["gene_id"], reason="multiple_gene_copies_within_one_haplotype")
    fields = list(rows[0])
    held = hold_unresolved(rows) if args.hold_unresolved_as_untestable else []
    if held: write_tsv(args.output/"untestable_annotations.tsv", held, fields)
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
    status = "held_unresolved_annotations" if unresolved else ("resolved_with_held_untestable" if held else "resolved")
    summary = {"status": status, "candidate_not_adopted": True,
               "assembly": "GCA_000001405.29_GRCh38.p14", "gtf": str(args.gtf.resolve()), "gtf_sha256": digest(args.gtf),
               "n_annotations": len(rows), "dispositions": counts, "n_unresolved": len(unresolved),
               "annotation_resolution_complete": not unresolved, "n_unresolved_annotations": len(unresolved),
               "counts_complete": False, "counts_sha256": None,
               "primary_refit_allowed": False, "symbol_only_merging_used": False,
               "partial_equivalent_exon_policy": "coordinate_verified_segments_only",
               "reference_exclusions_applied": args.apply_grc_reference_exclusions,
               "excluded_reference_annotations": sorted(excluded_targets),
               "excluded_reference_count": len(excluded_targets),
               "excluded_reference_evidence_state": "not_applicable_reference_artifact",
               "valid_and_excluded_alignment_policy": "any_unassigned_alignment_rejects_complete_read_or_pair",
               "annotation_correspondence_sha256": digest(args.output/"annotation_correspondence.tsv"),
               "n_held_untestable_annotations": len(held),
               "held_untestable_policy": ("own_counting_target_excluded_from_tested_matrix;"
                                          "shared_reads_rejected_in_primary;"
                                          "single_candidate_credit_in_sensitivity") if held else None,
               "grouped_counting_annotation": None, "grouped_counting_annotation_sha256": None}
    if not unresolved:
        write_grouped_annotation(args.gtf, args.output/"grouped_counting_annotation.gtf.gz", rows, mappings)
        grouped_path = args.output/"grouped_counting_annotation.gtf.gz"
        summary.update(grouped_counting_annotation=str(grouped_path.resolve()),
                       grouped_counting_annotation_sha256=digest(grouped_path))
    (args.output/"resolution_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__": main()
