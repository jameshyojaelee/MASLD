#!/usr/bin/env python3
"""H30 detectability gate from observed depth (run before any outcome is read).

For each library and FOXO1 hotspot (R21L, S22W):
  n   = number of distinct fragments covering the hotspot (distinct start, mate start and strand);
        PCR duplicates are unmarked in the STAR BAMs and are not independent observations
  e   = background rate of the same plus-strand substitution at the other positions of the
        201 bp window around c.65 (codon 21-22 positions excluded), pooled over both strands,
        per read base
  k*  = smallest alt-fragment count with P(Binomial(n, e) >= k*) < alpha
  MDV = smallest true alt fraction v with P(Binomial(n, v) >= k*) >= power
The test is one-sided against the background; alpha and power are arguments, set in
PRESPEC_H30 before outcomes are read. The output holds depth and background only, never the
alt count at the hotspots, so the gate is blind to the outcome. Needs numpy and scipy
(both in the rnaseq env).

Usage: h30_power_gate.py LIBRARIES_DIR OUT_TSV [alpha=0.001] [power=0.8]
"""
import csv
import glob
import os
import sys

import numpy as np
from scipy.stats import binom

root, out_path = sys.argv[1:3]
alpha = float(sys.argv[3]) if len(sys.argv) > 3 else 0.001
power = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
if os.path.exists(out_path):
    sys.exit("refusing to overwrite " + out_path)
EXCLUDE = set(range(40666147, 40666153))  # codons 21 and 22 on the plus strand
ALT_COL = {"C": "C", "A": "A", "T": "T", "G": "G"}


def mdv(depth, k):
    if depth == 0:
        return float("nan")
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if binom.sf(k - 1, depth, mid) >= power:
            hi = mid
        else:
            lo = mid
    return hi


rows = []
for d in sorted(glob.glob(os.path.join(root, "*", "*"))):
    hp, pu = os.path.join(d, "foxo1_hotspots.tsv"), os.path.join(d, "foxo1_pileup.tsv")
    if not (os.path.exists(hp) and os.path.exists(pu)):
        continue
    with open(pu) as fh:
        pile = list(csv.DictReader(fh, delimiter="\t"))
    with open(hp) as fh:
        hots = list(csv.DictReader(fh, delimiter="\t"))
    for h in hots:
        if h["hotspot"] == "S22X":
            continue
        ref, alt = h["ref_plus"], h["alt_plus"]
        bg_alt = bg_depth = 0
        for p in pile:
            if int(p["pos"]) in EXCLUDE or p["ref_plus"] != ref:
                continue
            cnt = {b: int(p[f"{b}_p"]) + int(p[f"{b}_m"]) for b in "ACGT"}
            bg_alt += cnt[alt]
            bg_depth += sum(cnt.values())
        e = (bg_alt + 0.5) / (bg_depth + 1.0)
        n_reads = int(h["depth"])
        depth = int(h["n_distinct_fragments_all"])
        k = int(binom.isf(alpha, depth, e)) + 1 if depth else 0
        rows.append([os.path.basename(os.path.dirname(d)), os.path.basename(d), h["hotspot"], n_reads, depth,
                     bg_depth, bg_alt, f"{e:.3g}", k, f"{mdv(depth, k):.4g}" if k else "NA"])
with open(out_path, "w") as fo:
    fo.write("cohort\tlibrary\thotspot\tdepth_reads\tdepth_distinct_fragments\tbackground_depth\tbackground_alt\tbackground_rate\tk_star\tmin_detectable_fraction\n")
    for r in rows:
        fo.write("\t".join(map(str, r)) + "\n")
print(f"{len(rows)} library-hotspot rows; alpha={alpha}; power={power}")
