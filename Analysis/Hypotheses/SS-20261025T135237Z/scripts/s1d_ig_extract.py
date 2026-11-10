#!/usr/bin/env python3
"""SS job S1d: IGH constant-region fragment counts for H38 from one STAR BAM.

The IGH block of scripts/analysis/hypotheses_2026q4/reads_extraction/extract_reads_one_library.py (lines 146-177,
364-371), unchanged, but reading only the IGH constant-exon span on chr14 and KI270846.1 through the BAM index
instead of the whole BAM. Rule, as in the original: every mapped alignment (multi-mappers included, flag 0x200
excluded) that overlaps a constant exon of ig_constant_exons.tsv adds its genes and HI tag to its fragment (read
name); a fragment is "all alignments in constant exons" when it has as many distinct HI as its NH. Output is the
original's ig_fragments.tsv, plus library_qc.tsv holding the STAR "Number of input reads" (H38 rate denominator).
Region reads give the same fragments as a full pass: an alignment outside the span cannot overlap a constant exon.

Usage: s1d_ig_extract.py --library S --bam B --log-final LOG --exons ig_constant_exons.tsv --out-dir DIR
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from collections import Counter

import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--library", required=True)
ap.add_argument("--bam", required=True)
ap.add_argument("--log-final", required=True)
ap.add_argument("--exons", required=True)
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()

t_start = time.time()
SCRIPT_SHA256 = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
if os.path.exists(args.out_dir):
    sys.exit(f"output exists; refusing: {args.out_dir}")
tmp_dir = f"{args.out_dir}.partial-{os.environ.get('SLURM_JOB_ID', os.getpid())}"
if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
os.makedirs(tmp_dir)

with open(args.exons) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    ig_rows = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]

af = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
tid_of = {n: i for i, n in enumerate(af.references)}
ig_by_tid, ig_range = {}, {}
for r_ in ig_rows:
    t_ = tid_of[r_["chrom"]]
    ig_by_tid.setdefault(t_, []).append((int(r_["start"]), int(r_["end"]), r_["gene"]))
for t_, lst in ig_by_tid.items():
    ig_range[t_] = (min(e[0] for e in lst), max(e[1] for e in lst))
ig_frag = {}  # qname -> [genes, his, nh]

n_rec = 0
for tid in sorted(ig_range):
    lo, hi = ig_range[tid]
    for r in af.fetch(af.references[tid], max(lo - 1, 0), hi):
        n_rec += 1
        flag = r.flag
        if flag & 4:
            continue
        start = r.reference_start
        end = r.reference_end
        rng = ig_range[tid]
        if start < rng[1] and end > rng[0] and not flag & 0x200:
            genes = {g for es, ee, g in ig_by_tid[tid] if es <= end and ee >= start + 1}
            if genes:
                rec = ig_frag.get(r.query_name)
                if rec is None:
                    rec = ig_frag[r.query_name] = [set(), set(), r.get_tag("NH") if r.has_tag("NH") else 1]
                rec[0].update(genes)
                rec[1].add(r.get_tag("HI") if r.has_tag("HI") else 1)

ig_count = Counter()
for qn, (genes, his, nh) in ig_frag.items():
    ig_count[(",".join(sorted(genes)), int(len(his) >= nh), "unique" if nh == 1 else "multi")] += 1
with open(os.path.join(tmp_dir, "ig_fragments.tsv"), "w") as fi:
    fi.write("library\tgene_set\tall_alignments_in_constant_exons\tmapping\tn_fragments\n")
    for (gs, comp, mp), n in sorted(ig_count.items()):
        fi.write(f"{args.library}\t{gs}\t{comp}\t{mp}\t{n}\n")

log = {}
with open(args.log_final) as fh:
    for line in fh:
        if "|" in line:
            k, v = line.split("|", 1)
            log[k.strip()] = v.strip()
qc = {"library": args.library, "bam_bytes": os.path.getsize(args.bam), "n_records_in_ig_span": n_rec,
      "n_ig_fragments": len(ig_frag), "star_input_reads": log.get("Number of input reads", "NA"),
      "star_unique_pct": log.get("Uniquely mapped reads %", "NA"), "seconds_total": round(time.time() - t_start, 1)}
with open(os.path.join(tmp_dir, "library_qc.tsv"), "w") as fq:
    fq.write("\t".join(qc) + "\n" + "\t".join(str(v) for v in qc.values()) + "\n")
with open(os.path.join(tmp_dir, "run_info.json"), "w") as fr:
    json.dump({"args": vars(args), "python": sys.version.split()[0], "pysam": pysam.__version__,
               "script_sha256": SCRIPT_SHA256, "exons_realpath": os.path.realpath(args.exons),
               "slurm_job": os.environ.get("SLURM_JOB_ID")}, fr, indent=1)
os.rename(tmp_dir, args.out_dir)
print(f"done {args.library}: {n_rec} records in IGH span, {len(ig_frag)} fragments, {qc['seconds_total']} s")
