"""Write the Tukiainen 2017 chrX XCI classification (Suppl. Tables 1 and 13) to one TSV.
Usage: python -I 01_parse_tukiainen.py EXTRACTED_DIR/nature24265-s3 OUT_TSV
Refuses to overwrite OUT_TSV. Coordinates in the source are GRCh37 (Gencode v19); only Ensembl IDs are used downstream."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xlsx_min import read_sheets

base, out = sys.argv[1], sys.argv[2]
if os.path.exists(out):
    sys.exit("refusing to overwrite " + out)

t1 = list(read_sheets(os.path.join(base, "Suppl.Table.1.xlsx")).values())[0]
t13 = list(read_sheets(os.path.join(base, "Suppl.Table.13.xlsx")).values())[0]
h1 = t1[1]
i = {h: k for k, h in enumerate(h1)}
rows = {}
for r in t1[2:]:
    if len(r) < 8 or not r[1]:
        continue
    ensg = r[1].split(".")[0]
    rows[ensg] = dict(gene=r[0], ensg=ensg, chr37=r[2], start37=r[3], end37=r[4], biotype=r[5],
                      combined_status=r[i["Combined XCI status"]], region_t13="", reported_t13="")
h13 = t13[0]
j = {h: k for k, h in enumerate(h13)}
n13_only = 0
for r in t13[1:]:
    if len(r) < 7 or not r[0]:
        continue
    ensg = r[0].split(".")[0]
    if ensg not in rows:
        n13_only += 1
        rows[ensg] = dict(gene=r[1], ensg=ensg, chr37=r[2], start37=r[3], end37=r[4], biotype="",
                          combined_status="", region_t13="", reported_t13="")
    rows[ensg]["region_t13"] = r[j["Region"]]
    rows[ensg]["reported_t13"] = r[j["Reported XCI status"]]

cols = ["gene", "ensg", "chr37", "start37", "end37", "biotype", "combined_status", "region_t13", "reported_t13"]
with open(out, "x") as fh:
    fh.write("\t".join(cols) + "\n")
    for v in rows.values():
        fh.write("\t".join(v[c] for c in cols) + "\n")

print("genes", len(rows), "| only in Table 13:", n13_only)
print("Table 1 columns:", h1)
print("combined_status:", collections.Counter(v["combined_status"] for v in rows.values()))
print("region_t13:", collections.Counter(v["region_t13"] for v in rows.values()))
print("reported_t13:", collections.Counter(v["reported_t13"] for v in rows.values()))
print("chr37:", collections.Counter(v["chr37"] for v in rows.values()))
