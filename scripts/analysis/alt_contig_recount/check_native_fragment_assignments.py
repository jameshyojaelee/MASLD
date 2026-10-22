#!/usr/bin/env python3
"""Check strict molecule accounting against installed featureCounts assignments."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import subprocess
import pysam
from count_verified_fragments import molecule_assignment


def write_bam(path, paired):
    header={"HD":{"VN":"1.6","SO":"unsorted"},"SQ":[{"SN":c,"LN":1000} for c in ["chr1","ALT1","chr2","excluded_reference"]]}
    with pysam.AlignmentFile(str(path),"wb",header=header) as bam:
        examples=[("unique",[(0,1,1)]),("same_gene",[(0,2,1),(1,2,2)]),
                  ("cross_gene",[(0,2,1),(2,2,2)]),("incomplete",[(0,2,1)]),
                  ("valid_and_excluded",[(0,2,1),(3,2,2)]),("excluded_only",[(3,1,1)])]
        for name,hits in examples:
            for chrom,nh,hi in hits:
                for mate in ([1,2] if paired else [0]):
                    r=pysam.AlignedSegment();r.query_name=name;r.query_sequence="A"*20
                    r.flag=(83 if mate==1 else 163) if paired else 16
                    if hi>1:r.flag|=256
                    r.reference_id=chrom;r.reference_start=110 if mate!=2 else 140
                    r.mapping_quality=255;r.cigarstring="20M"
                    r.query_qualities=pysam.qualitystring_to_array("I"*20)
                    if paired:
                        r.next_reference_id=chrom;r.next_reference_start=140 if mate==1 else 110
                        r.template_length=50 if mate==1 else -50
                    r.set_tags([("NH",nh),("HI",hi)])
                    bam.write(r)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--featurecounts",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    annotation=args.output/"fixture.gtf"
    annotation.write_text('chr1\tfixture\texon\t101\t200\t.\t+\t.\tgene_id "G1";\n'
                          'ALT1\tfixture\texon\t101\t200\t.\t+\t.\tgene_id "G1";\n'
                          'chr2\tfixture\texon\t101\t200\t.\t+\t.\tgene_id "G2";\n')
    expected={"unique":("G1","accepted_unique"),"same_gene":("G1","accepted_same_gene_multimapper"),
              "cross_gene":(None,"cross_gene_alignment_set"),"incomplete":(None,"incomplete_NH_alignment_set"),
              "valid_and_excluded":(None,"unassigned_or_multigene_candidate"),
              "excluded_only":(None,"unassigned_or_multigene_candidate")}
    result={"fixture_only":True,"scientific_outputs":False,"checks":[]}
    for paired in [False,True]:
        label="paired" if paired else "single";bam=args.output/f"{label}.bam";write_bam(bam,paired)
        cmd=[str(args.featurecounts),"-T","1","-s","2","-M","-R","BAM","-a",str(annotation),"-o",str(args.output/f"{label}.counts.tsv")]
        if paired:cmd.extend(["-p","--countReadPairs","-B"])
        cmd.append(str(bam))
        completed=subprocess.run(cmd,text=True,capture_output=True)
        (args.output/f"{label}.featurecounts.log").write_text(completed.stdout+completed.stderr)
        completed.check_returncode();groups=defaultdict(list)
        with pysam.AlignmentFile(str(bam)+".featureCounts.bam","rb") as assigned:
            for r in assigned.fetch(until_eof=True):
                tags=dict(r.tags)
                groups[r.query_name].append({"NH":tags.get("NH"),"HI":tags.get("HI"),"XS":tags.get("XS"),
                    "XN":tags.get("XN"),"XT":tags.get("XT"),"unmapped":r.is_unmapped,
                    "supplementary":r.is_supplementary,"mate":1 if r.is_read1 else (2 if r.is_read2 else 0)})
                groups[r.query_name][-1].update(chimeric=r.has_tag("SA") or r.has_tag("ch"),
                    paired=r.is_paired,mate_unmapped=r.mate_is_unmapped,
                    reference_id=r.reference_id,reference_start=r.reference_start,
                    next_reference_id=r.next_reference_id,next_reference_start=r.next_reference_start)
        observed={name:molecule_assignment(records,paired) for name,records in groups.items()}
        if observed!=expected:raise AssertionError(f"{label}: expected={expected},observed={observed}")
        result["checks"].append({"layout":label,"command":cmd,"n_fixture_molecules":len(observed),
            "accepted_molecules":sum(gene is not None for gene,status in observed.values()),
            "observed":observed,"status":"pass"})
    result["status"]="pass"
    (args.output/"native_fixture_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":main()
