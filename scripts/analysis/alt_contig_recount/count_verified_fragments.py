#!/usr/bin/env python3
"""Count one molecule only when its complete STAR alignment set identifies a gene.

Input: featureCounts -M -R BAM assignments on the fully resolved grouped GTF,
then samtools sort -N. Paired libraries retain --countReadPairs -B; strandedness
and overlap decisions are those of the original featureCounts cohort command.
The featureCounts multimapper count table is deliberately not consumed.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from itertools import groupby
import json
from pathlib import Path


def molecule_assignment(records, paired):
    """Return (gene or None, mutually exclusive disposition) for one QNAME/RG."""
    if not records: return None, "empty_group"
    if any(r.get("supplementary", False) for r in records):
        return None, "supplementary_alignment"
    if any(r.get("chimeric", False) for r in records):
        return None, "chimeric_alignment"
    if any(r.get("unmapped", False) for r in records):
        return None, "unmapped_or_missing_mate"
    nh = {r.get("NH") for r in records}
    if len(nh) != 1 or None in nh or next(iter(nh)) < 1:
        return None, "missing_or_inconsistent_NH"
    expected = next(iter(nh))
    hits = defaultdict(list)
    for r in records:
        hi = r.get("HI", 1 if expected == 1 else None)
        if hi is None: return None, "missing_HI_for_multimapper"
        hits[hi].append(r)
    if len(hits) != expected:
        return None, "incomplete_NH_alignment_set"
    if set(hits) != set(range(1, expected+1)):
        return None, "invalid_HI_alignment_set"
    for hit in hits.values():
        if paired:
            if len(hit) != 2 or {r.get("mate") for r in hit} != {1, 2}:
                return None, "incomplete_or_duplicated_pair"
            if any(not r.get("paired") or r.get("mate_unmapped") for r in hit):
                return None, "unmapped_or_missing_mate"
            a, b = hit
            if a.get("reference_id") is None or b.get("reference_id") is None:
                return None, "missing_pair_coordinates"
            if a["reference_id"] != b["reference_id"]:
                return None, "cross_contig_pair"
            if any(r.get("next_reference_id") != mate["reference_id"] or
                   r.get("next_reference_start") != mate.get("reference_start")
                   for r, mate in ((a, b), (b, a))):
                return None, "nonreciprocal_mate_coordinates"
        elif len(hit) != 1:
            return None, "duplicated_single_end_hit"
    targets = set()
    for r in records:
        if r.get("XS") != "Assigned" or r.get("XN") != 1 or not r.get("XT"):
            return None, "unassigned_or_multigene_candidate"
        if "," in r["XT"]:
            return None, "unassigned_or_multigene_candidate"
        targets.add(r["XT"])
    if len(targets) != 1:
        return None, "cross_gene_alignment_set"
    return next(iter(targets)), "accepted_unique" if expected == 1 else "accepted_same_gene_multimapper"


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for b in iter(lambda:handle.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def main():
    import pysam
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotation-summary",type=Path,required=True)
    p.add_argument("--assignment-bam",type=Path,required=True)
    p.add_argument("--sample-id",required=True)
    p.add_argument("--layout",choices=["SINGLE","PAIRED"],required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    summary=json.loads(args.annotation_summary.read_text())
    if summary.get("annotation_resolution_complete") is not True or summary.get("n_unresolved_annotations") != 0:
        raise ValueError("Primary counts require full annotation correspondence resolution")
    with (args.annotation_summary.parent/"annotation_correspondence.tsv").open() as f:
        permitted={r["counting_gene_id"] for r in csv.DictReader(f,delimiter="\t")}
    args.output.mkdir(parents=True,exist_ok=False)
    counts=Counter(); dispositions=Counter(); previous=None
    with pysam.AlignmentFile(str(args.assignment_bam),"rb") as bam:
        if bam.header.to_dict().get("HD",{}).get("SO") != "queryname":
            raise ValueError("Assignments must be name-sorted with samtools sort -N")
        for name, group in groupby(bam.fetch(until_eof=True),key=lambda r:r.query_name):
            if previous is not None and name < previous:
                raise ValueError("QNAMEs are not lexically ordered; use samtools sort -N")
            previous=name; groups=defaultdict(list)
            for r in group:
                tags={k:v for k,v in r.tags}
                groups[tags.get("RG","")].append({"NH":tags.get("NH"),"HI":tags.get("HI"),
                    "XS":tags.get("XS"),"XN":tags.get("XN"),"XT":tags.get("XT"),
                    "unmapped":r.is_unmapped,"supplementary":r.is_supplementary,
                    "chimeric":r.has_tag("SA") or r.has_tag("ch"),
                    "paired":r.is_paired,"mate_unmapped":r.mate_is_unmapped,
                    "reference_id":r.reference_id,"reference_start":r.reference_start,
                    "next_reference_id":r.next_reference_id,"next_reference_start":r.next_reference_start,
                    "mate":1 if r.is_read1 else (2 if r.is_read2 else 0)})
            for records in groups.values():
                gene, status=molecule_assignment(records,args.layout=="PAIRED")
                dispositions[status]+=1
                if gene is not None:
                    if gene not in permitted: raise ValueError(f"Assigned gene absent from resolved GTF: {gene}")
                    counts[gene]+=1
    out=args.output/"gene_counts.tsv"
    with out.open("x") as f:
        writer=csv.writer(f,delimiter="\t",lineterminator="\n")
        writer.writerow(["gene_id",args.sample_id])
        for gene in sorted(permitted): writer.writerow([gene,counts[gene]])
    certificate={"sample_id":args.sample_id,"layout":args.layout,"candidate_not_adopted":True,
        "annotation_summary_sha256":digest(args.annotation_summary),"counts_tsv_sha256":digest(out),
        "n_molecules":sum(dispositions.values()),"n_assigned_molecules":sum(counts.values()),
        "dispositions":dict(dispositions),"maximum_contribution_per_molecule":1,"fractional_renormalization":False}
    (args.output/"count_completion.json").write_text(json.dumps(certificate,indent=2)+"\n")
    print(json.dumps(certificate,sort_keys=True))


if __name__=="__main__": main()
