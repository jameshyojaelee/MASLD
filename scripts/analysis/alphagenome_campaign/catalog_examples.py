#!/usr/bin/env python3
"""Reproducible identity-selected queries and a synthetic phased-sequence request."""
import argparse
import json
from pathlib import Path
import sqlite3

from catalog_query import KINDS, query


def examples(db, out):
    out.mkdir(parents=True, exist_ok=False)
    connection = sqlite3.connect(db.resolve().as_uri()+"?mode=ro", uri=True)
    for kind in KINDS:
        # Choose the smallest stable identity, never the largest effect.
        identity = connection.execute("SELECT identifier FROM entity WHERE kind=? ORDER BY identifier LIMIT 1",(kind,)).fetchone()[0]
        if kind == "variant":
            identity = connection.execute("SELECT identifier FROM entity WHERE kind='variant' AND identifier LIKE 'GRCh38:%' ORDER BY identifier LIMIT 1").fetchone()[0]
        (out/(kind+".json")).write_text(json.dumps(query(db,kind,identity,limit=10),indent=2,allow_nan=False)+"\n")
    variants = []
    for (key,) in connection.execute("SELECT identifier FROM entity WHERE kind='variant' AND identifier LIKE 'GRCh38:%'"):
        p = key.split(":")
        if len(p)==5 and len(p[3])==len(p[4])==1:
            variants.append((p[1],int(p[2]),key))
    variants.sort()
    pair = next((a[2],b[2]) for a,b in zip(variants,variants[1:]) if a[0]==b[0] and 0 < b[1]-a[1] < 500 and a[1]>1024)
    (out/"synthetic_haplotype_request.json").write_text(json.dumps({
        "variants":list(pair), "phase":"same_haplotype", "phase_source":"synthetic_sequence_demonstration_not_observed_phase",
        "selection":"first_coordinate_sorted_adjacent_SNV_pair_within500bp; no_outcome_selection",
        "natural_haplotype_validation":False},indent=2)+"\n")
    connection.close()


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    examples(a.db,a.out)
