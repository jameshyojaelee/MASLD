#!/usr/bin/env python3
"""Count one read/fragment only when its complete STAR alignment set identifies a gene.

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
import shutil


def molecule_assignment(records, paired, alias=None):
    """Return one read/pair assignment for a QNAME/RG, without UMI deduplication.

    alias maps a held annotation to its single candidate primary gene (sensitivity only).
    """
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
        targets.add(alias.get(r["XT"], r["XT"]) if alias else r["XT"])
    if len(targets) != 1:
        return None, "cross_gene_alignment_set"
    return next(iter(targets)), "accepted_unique" if expected == 1 else "accepted_same_gene_multimapper"


def held_alias(rows):
    """Held annotation -> counting gene of its only primary candidate, when exactly one exists."""
    held=[r for r in rows if r.get("disposition")=="held_untestable"]
    if not held: return {}
    target={r["gene_id"]:r["counting_gene_id"] for r in rows if r.get("disposition") not in ("held_untestable","excluded_reference")}
    alias={}
    for r in held:
        candidates=[c for c in r.get("primary_gene_candidates","").split(";") if c]
        if len(candidates)==1 and target.get(candidates[0]): alias[r["gene_id"]]=target[candidates[0]]
    return alias


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for b in iter(lambda:handle.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def file_stat(path):
    path=Path(path).resolve(strict=True)
    info=path.stat()
    return {"path":str(path),"size_bytes":info.st_size,"mtime_ns":str(info.st_mtime_ns)}


def validate_count_execution(summary_path, evidence_path, assignment_bam, sample, layout):
    """Bind assignments to the annotation and actual command without hashing BAMs."""
    summary=json.loads(Path(summary_path).read_text())
    grouped=Path(summary["grouped_counting_annotation"]).resolve(strict=True)
    correspondence=Path(summary_path).parent/"annotation_correspondence.tsv"
    grouped_sha=digest(grouped); correspondence_sha=digest(correspondence)
    if grouped_sha!=summary.get("grouped_counting_annotation_sha256") or correspondence_sha!=summary.get("annotation_correspondence_sha256"):
        raise ValueError("Grouped annotation or correspondence hash differs")
    evidence=json.loads(Path(evidence_path).read_text())
    if evidence.get("status")!="complete" or evidence.get("sample_id")!=sample or evidence.get("layout")!=layout or evidence.get("featurecounts_returncode")!=0 or evidence.get("samtools_sort_returncode")!=0:
        raise ValueError("Assignment execution is not complete for this library")
    if evidence.get("annotation_summary_sha256")!=digest(summary_path) or evidence.get("grouped_counting_annotation_sha256")!=grouped_sha or evidence.get("annotation_correspondence_sha256")!=correspondence_sha:
        raise ValueError("Executed annotation evidence differs")
    argv=evidence.get("featurecounts_argv",[])
    def option(name):
        if argv.count(name)!=1 or argv.index(name)+1>=len(argv):
            raise ValueError("Required featureCounts option absent or repeated: "+name)
        return argv[argv.index(name)+1]
    if not argv or Path(argv[0]).name!="featureCounts" or option("-a")!=str(grouped) or option("-s")!="2" or option("-R")!="BAM" or argv.count("-M")!=1:
        raise ValueError("Executed featureCounts annotation/strand/multimapper flags differ")
    paired_flags=("-p","--countReadPairs","-B")
    if any(argv.count(flag)!=(1 if layout=="PAIRED" else 0) for flag in paired_flags) or any(flag in argv for flag in ("--fraction","-O","-C","--primary","--donotsort")):
        raise ValueError("Executed fragment eligibility flags differ")
    value_options={"-T","-s","-R","-a","-o","--Rpath"}
    switches={"-M",*paired_flags}
    position=1
    while position<len(argv)-1:
        token=argv[position]
        if token in value_options:
            option(token)
            position+=2
        elif token in switches:
            position+=1
        else:
            raise ValueError("Unrecognized featureCounts eligibility option: "+token)
    before=evidence.get("source_bam_before",{})
    if not before.get("path") or before!=evidence.get("source_bam_after") or before!=file_stat(before["path"]) or argv[-1]!=before["path"]:
        raise ValueError("Source BAM path/stat changed during assignment")
    if Path(before["path"]).name.split(".",1)[0]!=sample:
        raise ValueError("Source BAM run ID differs from the counted sample")
    if evidence.get("assignment_bam")!=file_stat(assignment_bam):
        raise ValueError("Assignment BAM does not match executed output")
    sort=evidence.get("samtools_sort_argv",[])
    if len(sort)<3 or Path(sort[0]).name!="samtools" or sort[1]!="sort" or "-N" not in sort or sort.count("-o")!=1 or sort[sort.index("-o")+1]!=str(Path(assignment_bam).resolve()):
        raise ValueError("Assignment BAM lacks lexical name-sort execution")
    assigned=Path(sort[-1])
    if assigned.name!=Path(before["path"]).name+".featureCounts.bam" or assigned.parent.resolve()!=Path(option("--Rpath")).resolve():
        raise ValueError("Name-sort input does not identify featureCounts assignment output")
    versions=evidence.get("versions",{})
    if not all(isinstance(versions.get(name),str) and versions[name] for name in ("python","featureCounts","samtools")):
        raise ValueError("Assignment executable versions missing")
    return {"grouped_counting_annotation":str(grouped),"grouped_counting_annotation_sha256":grouped_sha,
            "annotation_correspondence_sha256":correspondence_sha,
            "source_bam":before,"assignment_execution_sha256":digest(evidence_path)}


def main():
    import pysam
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotation-summary",type=Path,required=True)
    p.add_argument("--assignment-bam",type=Path,required=True)
    p.add_argument("--execution-evidence",type=Path,required=True)
    p.add_argument("--sample-id",required=True)
    p.add_argument("--layout",choices=["SINGLE","PAIRED"],required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    reader_source_sha256=digest(Path(__file__))
    summary=json.loads(args.annotation_summary.read_text())
    if summary.get("annotation_resolution_complete") is not True or summary.get("n_unresolved_annotations") != 0:
        raise ValueError("Primary counts require full annotation correspondence resolution")
    execution=validate_count_execution(args.annotation_summary,args.execution_evidence,args.assignment_bam,args.sample_id,args.layout)
    with (args.annotation_summary.parent/"annotation_correspondence.tsv").open() as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    permitted={r["counting_gene_id"] for r in rows
               if r["counting_gene_id"] and r.get("disposition") != "excluded_reference"}
    held={r["gene_id"] for r in rows if r.get("disposition") == "held_untestable"}
    if len(held) != summary.get("n_held_untestable_annotations",0):
        raise ValueError("Held annotations differ from the annotation summary")
    alias=held_alias(rows); testable=permitted-held
    args.output.mkdir(parents=True,exist_ok=False)
    counts=Counter(); sensitivity=Counter(); dispositions=Counter(); previous=None
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
                # The sensitivity differs only for molecules that touch an aliased held annotation.
                if alias and any(r["XT"] in alias for r in records):
                    gene=molecule_assignment(records,args.layout=="PAIRED",alias)[0]
                if gene is not None: sensitivity[gene]+=1
    def write_counts(name,values,genes):
        path=args.output/name
        with path.open("x") as f:
            writer=csv.writer(f,delimiter="\t",lineterminator="\n")
            writer.writerow(["gene_id",args.sample_id])
            for gene in sorted(genes): writer.writerow([gene,values[gene]])
        return path
    out=write_counts("gene_counts.tsv",counts,testable)
    sensitivity_out=write_counts("gene_counts_single_candidate_credit.tsv",sensitivity,testable)
    held_out=write_counts("held_untestable_counts.tsv",counts,held)
    certificate={"sample_id":args.sample_id,"layout":args.layout,"candidate_not_adopted":True,
        "annotation_summary_sha256":digest(args.annotation_summary),"counts_tsv_sha256":digest(out),
        "n_molecules":sum(dispositions.values()),"n_assigned_molecules":sum(counts[g] for g in testable),
        "n_held_untestable_annotations":len(held),"n_held_untestable_molecules":sum(counts[g] for g in held),
        "n_single_candidate_credit_annotations":len(alias),
        "n_single_candidate_credit_molecules":sum(sensitivity[g] for g in testable),
        "single_candidate_credit_counts_tsv_sha256":digest(sensitivity_out),
        "held_untestable_counts_tsv_sha256":digest(held_out),
        "dispositions":dict(dispositions),"maximum_contribution_per_molecule":1,"fractional_renormalization":False}
    certificate["sequencing_unit"]="QNAME/RG read or read pair; not UMI-deduplicated biological molecules"
    copied_evidence=args.output/"assignment_execution.json"
    shutil.copyfile(args.execution_evidence,copied_evidence)
    certificate.update(execution)
    certificate["assignment_execution"]=str(copied_evidence.resolve())
    certificate["reader_versions"]={"pysam":pysam.__version__}
    certificate["count_reader_source_sha256"]=reader_source_sha256
    (args.output/"count_completion.json").write_text(json.dumps(certificate,indent=2)+"\n")
    print(json.dumps(certificate,sort_keys=True))


if __name__=="__main__": main()
