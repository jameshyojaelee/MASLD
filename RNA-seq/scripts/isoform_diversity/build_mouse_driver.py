#!/usr/bin/env python3
"""
Build the mouse kallisto requant driver by harvesting exact FASTQ paths from the
existing mouse kallisto run_info.json 'call' records, restricted to the samples in
the unified mouse integration metadata (the QC set), with per-sample layout.

Output: mouse_driver.tsv  cols: dataset, sample_id, layout, fastq_r1, fastq_r2
"""
import csv, glob, json, os, re, sys, collections

PROJ = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
META = f"{PROJ}/RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv"
OUT  = f"{PROJ}/RNA-seq/scripts/isoform_diversity/mouse_driver.tsv"

# metadata: sample_id -> (dataset, layout)
meta = {}
with open(META) as fh:
    for r in csv.DictReader(fh):
        meta[r["sample_id"]] = (r["dataset"], r.get("layout", "").upper())

# harvest fastq paths from existing run_info.json 'call' fields
fq_re = re.compile(r'(/\S+?\.(?:fastq|fq)(?:\.gz)?)')
found = {}   # sample_id -> [fastqs]
for ri in glob.glob(f"{PROJ}/RNA-seq/Mouse/*/quant/kallisto/*/run_info.json") + \
          glob.glob(f"{PROJ}/RNA-seq/Mouse/*/*/quant/kallisto/*/run_info.json"):
    sid = os.path.basename(os.path.dirname(ri))
    if sid not in meta:
        continue
    try:
        call = json.load(open(ri)).get("call", "")
    except Exception:
        continue
    fqs = fq_re.findall(call)
    if fqs:
        found[sid] = fqs

rows, missing = [], []
for sid, (dataset, layout) in meta.items():
    fqs = found.get(sid)
    if not fqs:
        missing.append(sid); continue
    fqs = [f for f in fqs if os.path.exists(f)]
    if not fqs:
        missing.append(sid); continue
    # heuristic R1/R2: sort; if 2+ and looks paired use first two
    fqs_sorted = sorted(fqs)
    if layout == "PAIRED" and len(fqs_sorted) >= 2:
        r1, r2 = fqs_sorted[0], fqs_sorted[1]
    else:
        r1, r2 = fqs_sorted[0], ""
    rows.append((dataset, sid, layout or ("PAIRED" if r2 else "SINGLE"), r1, r2))

with open(OUT, "w") as out:
    out.write("dataset\tsample_id\tlayout\tfastq_r1\tfastq_r2\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

byds = collections.Counter(r[0] for r in rows)
lay  = collections.Counter(r[2] for r in rows)
sys.stderr.write(f"[mouse_driver] {len(rows)} samples ; datasets={dict(byds)} ; layout={dict(lay)}\n")
sys.stderr.write(f"[mouse_driver] missing fastq for {len(missing)} samples: {missing[:10]}\n")
sys.stderr.write(f"[mouse_driver] wrote {OUT}\n")
