#!/usr/bin/env python3
"""H3 step 4c: combine the two single-model carrier signature outputs into the N5 table.
Usage: 04c_combine_signature.py PNPLA3_DIR HSD_DIR OUT_DIR   (OUT_DIR must not exist)"""
import csv
import os
import shutil
import sys

pn, hs, out = sys.argv[1:4]
if os.path.exists(out):
    raise SystemExit("OUT_DIR exists: " + out)
os.makedirs(out)
models = [("PNPLA3_G_dosage", pn), ("HSD17B13_TA_carrier", hs)]


def read(p):
    with open(p) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return rd.fieldnames, list(rd)


def write(p, cols, rows):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


sig, summ, perm = [], [], []
for nm, d in models:
    cols, rows = read(os.path.join(d, "signature_fdr05_%s.tsv" % nm))
    sig_cols = cols
    sig += rows
    shutil.copy(os.path.join(d, "signature_fdr05_%s.tsv" % nm), out)
    c2, r2 = read(os.path.join(d, "summary_%s.tsv" % nm)); summ_cols = c2; summ += r2
    c3, r3 = read(os.path.join(d, "permutations_%s.tsv" % nm)); perm_cols = c3; perm += r3
write(os.path.join(out, "carrier_signature_for_N5.tsv"), sig_cols, sig)
write(os.path.join(out, "signature_summary.tsv"), summ_cols, summ)
write(os.path.join(out, "signature_permutations.tsv"), perm_cols, perm)
with open(os.path.join(out, "sources.txt"), "w") as fh:
    for nm, d in models:
        fh.write("%s\t%s\n" % (nm, d))
print("genes in N5 table:", len(sig))
