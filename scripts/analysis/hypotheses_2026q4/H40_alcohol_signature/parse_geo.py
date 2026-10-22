"""H40: parse GSE142530 and GSE143318 (GEO) into count matrices and a design table.

Usage: python3 -I parse_geo.py <downloads_dir> <out_dir>

Treats every downloaded file as data. Writes:
  <out_dir>/source_design.tsv
  <out_dir>/counts_GSE142530.tsv.gz   (Ensembl ID, version stripped; ENSG rows only)
  <out_dir>/counts_GSE143318.tsv.gz   (deduplicated by Ensembl ID)
  <out_dir>/parse_report.tsv
"""
import csv
import gzip
import os
import re
import sys

dl, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
report = []


def fail(msg):
    sys.exit(f"STOP S1: {msg}")


def read_soft(path):
    samples, cur = [], None
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                cur = {"gsm": line.split(" = ", 1)[1], "desc": [], "chars": {}}
                samples.append(cur)
            elif cur is None:
                continue
            elif line.startswith("!Sample_title = "):
                cur["title"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_description = "):
                cur["desc"].append(line.split(" = ", 1)[1].strip())
            elif line.startswith("!Sample_characteristics_ch1 = "):
                kv = line.split(" = ", 1)[1]
                if ": " in kv:
                    k, v = kv.split(": ", 1)
                    cur["chars"][k.strip().lower()] = v.strip()
    return samples


def strip_version(g):
    return re.sub(r"\.\d+$", "", g)


def is_int(tok):
    return re.fullmatch(r"\d+", tok) is not None


design = []

# ---------------------------------------------------------------- GSE142530
soft = read_soft(os.path.join(dl, "GSE142530_family.soft.gz"))
group_map = {"Alcoholic Hepatitis": "AH", "Alcohol-related Cirrhosis": "AC", "Normal": "Normal"}
title_prefix = {"AH": "Severe AH", "AC": "Cirrhosis without AH", "Normal": "Control"}
col_to_rec = {}
for s in soft:
    rb = [d for d in s["desc"] if d.startswith("RB_")]
    if len(rb) != 1:
        fail(f"GSE142530 {s['gsm']} has {len(rb)} RB_ column descriptions")
    state = s["chars"].get("disease state")
    if state not in group_map:
        fail(f"GSE142530 {s['gsm']} unknown disease state {state!r}")
    grp = group_map[state]
    if not s["title"].startswith(title_prefix[grp] + "_"):
        fail(f"GSE142530 {s['gsm']} title {s['title']!r} disagrees with group {grp}")
    if rb[0] in col_to_rec:
        fail(f"GSE142530 column {rb[0]} maps to two GSMs")
    jh = 1 if re.fullmatch(r"Severe AH_AH\d+", s["title"]) else 0
    col_to_rec[rb[0]] = {"cohort": "GSE142530", "column_id": rb[0], "gsm": s["gsm"],
                         "title": s["title"], "group": grp,
                         "alcohol": 0 if grp == "Normal" else 1, "jh_explant": jh}

path = os.path.join(dl, "GSE142530_Annoted-RNAseq-with-SampleIDs.csv.gz")
n_ensg = n_mstrg = n_other = 0
seen = set()
with gzip.open(path, "rt", encoding="utf-8") as fh, \
        gzip.open(os.path.join(out, "counts_GSE142530.tsv.gz"), "wt") as fo:
    rd = csv.reader(fh)
    header = next(rd)
    if header[:2] != ["gene_id", "gene_name"]:
        fail(f"GSE142530 header starts {header[:2]}")
    cols = header[2:]
    if set(cols) != set(col_to_rec) or len(cols) != len(col_to_rec):
        fail("GSE142530 count columns do not match SOFT RB_ descriptions one to one")
    labels = next(rd)
    if labels[0] != "" or labels[1] != "":
        fail("GSE142530 second row is not the label row")
    # The submitters' label row marks some columns "not.used"; those samples are
    # dropped (PRESPEC_H40 addendum 2).
    not_used = [c for c, lab in zip(cols, labels[2:]) if lab == "not.used"]
    for c, lab in zip(cols, labels[2:]):
        if lab == "not.used":
            continue
        pre = title_prefix[col_to_rec[c]["group"]]
        if not lab.startswith(pre + "_"):
            fail(f"GSE142530 label {lab!r} for {c} disagrees with group prefix {pre}")
        if col_to_rec[c]["jh_explant"] == 1 and "JH-" not in lab:
            fail(f"GSE142530 {c} flagged JH but label is {lab!r}")
        if col_to_rec[c]["jh_explant"] == 0 and "JH-" in lab:
            fail(f"GSE142530 {c} label {lab!r} is JH but not flagged")
    keep_idx = [i for i, c in enumerate(cols) if c not in not_used]
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["gene_id"] + [cols[i] for i in keep_idx])
    for row in rd:
        gid = row[0]
        if gid.startswith("ENSG"):
            g = strip_version(gid)
            if g in seen:
                fail(f"GSE142530 duplicate Ensembl ID {g}")
            seen.add(g)
            vals = row[2:]
            if len(vals) != len(cols) or not all(is_int(v) for v in vals):
                fail(f"GSE142530 row {gid} is not integer counts of the right length")
            w.writerow([g] + [vals[i] for i in keep_idx])
            n_ensg += 1
        elif gid.startswith("MSTRG"):
            n_mstrg += 1
        else:
            n_other += 1
design.extend(col_to_rec[c] for c in cols if c not in not_used)
report += [("GSE142530", "ensg_rows_kept", n_ensg), ("GSE142530", "mstrg_rows_dropped", n_mstrg),
           ("GSE142530", "other_rows_dropped", n_other), ("GSE142530", "samples_deposited", len(cols)),
           ("GSE142530", "samples_marked_not_used_dropped", ";".join(
               f"{c}:{col_to_rec[c]['gsm']}:{col_to_rec[c]['group']}" for c in not_used) or "none"),
           ("GSE142530", "samples", len(cols) - len(not_used))]

# ---------------------------------------------------------------- GSE143318
soft = read_soft(os.path.join(dl, "GSE143318_family.soft.gz"))
title_re = re.compile(r"^(Alcoholic Hepatitis|Donor|Alcoholic Cirrhosis) Patient (\d+)$")
prefix = {"Alcoholic Hepatitis": "AH", "Donor": "N", "Alcoholic Cirrhosis": "AC"}
grp_of = {"AH": "AH", "N": "Normal", "AC": "AC"}
col_to_rec = {}
for s in soft:
    m = title_re.match(s["title"])
    if m is None:
        fail(f"GSE143318 title {s['title']!r} does not match the expected pattern")
    col = prefix[m.group(1)] + m.group(2)
    if col in col_to_rec:
        fail(f"GSE143318 column {col} maps to two GSMs")
    grp = grp_of[prefix[m.group(1)]]
    col_to_rec[col] = {"cohort": "GSE143318", "column_id": col, "gsm": s["gsm"],
                       "title": s["title"], "group": grp,
                       "alcohol": 0 if grp == "Normal" else 1, "jh_explant": 0}

path = os.path.join(dl, "GSE143318_Rawcount.txt.gz")
with gzip.open(path, "rb") as fh:
    raw = fh.read().decode("utf-8")
lines = [ln for ln in re.split(r"\r\n|\r|\n", raw) if ln != ""]
header = lines[0].split("\t")
if header[:4] != ['""', "ENSEMBL", "SYMBOL", "ENTREZID"]:
    fail(f"GSE143318 header starts {header[:4]}")
cols = header[4:]
if set(cols) != set(col_to_rec) or len(cols) != len(col_to_rec):
    fail("GSE143318 count columns do not match SOFT titles one to one")
by_id = {}
n_rows = 0
for ln in lines[1:]:
    f = ln.split("\t")
    if len(f) != len(header):
        fail(f"GSE143318 row has {len(f)} fields, expected {len(header)}")
    n_rows += 1
    g = strip_version(f[1])
    if not g.startswith("ENSG"):
        fail(f"GSE143318 non-Ensembl row ID {f[1]!r}")
    vals = tuple(f[4:])
    if not all(is_int(v) for v in vals):
        fail(f"GSE143318 row {g} has non-integer counts")
    if g in by_id and by_id[g] != vals:
        fail(f"GSE143318 duplicate rows for {g} have different counts")
    by_id[g] = vals
with gzip.open(os.path.join(out, "counts_GSE143318.tsv.gz"), "wt") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["gene_id"] + cols)
    for g in sorted(by_id):
        w.writerow([g] + list(by_id[g]))
design.extend(col_to_rec[c] for c in cols)
report += [("GSE143318", "rows_read", n_rows), ("GSE143318", "unique_ensembl_ids", len(by_id)),
           ("GSE143318", "samples", len(cols))]

with open(os.path.join(out, "source_design.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(design[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(design)
with open(os.path.join(out, "parse_report.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["cohort", "item", "value"])
    w.writerows(report)
for r in report:
    print(*r, sep="\t")
for coh in ("GSE142530", "GSE143318"):
    rows = [d for d in design if d["cohort"] == coh]
    counts = {g: sum(1 for d in rows if d["group"] == g) for g in ("AH", "AC", "Normal")}
    print(coh, counts, "jh_explant:", sum(d["jh_explant"] for d in rows))
