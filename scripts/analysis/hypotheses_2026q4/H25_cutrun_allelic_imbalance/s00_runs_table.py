#!/usr/bin/env python3
"""H25: join SRA runs (PRJNA1109796) to GSE267119 samples and histology metadata.

Grade = NAS = steatosis + ballooning + lobular inflammation; high grade = NAS >= 4.
Reads metadata only.
"""
import argparse
import csv
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runinfo", required=True)
    ap.add_argument("--geo-samples", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs = {r["SampleName"]: r for r in csv.DictReader(open(a.runinfo)) if r["BioProject"] == "PRJNA1109796"}
    rows = []
    for g in csv.DictReader(open(a.geo_samples), delimiter="\t"):
        ch = dict(c.split(": ", 1) for c in json.loads(g["characteristics_json"]))
        r = runs[g["accession"]]
        nas = int(ch["steatosis"]) + int(ch["balloning"]) + int(ch["lobular inflammation"])
        rows.append(dict(srr=r["Run"], gsm=g["accession"], participant=g["candidate_participant_id"], title=g["title"],
                         stage=ch["Stage"], sex=ch["Sex"], steatosis=ch["steatosis"], ballooning=ch["balloning"],
                         lobular_inflammation=ch["lobular inflammation"], fibrosis=ch["fibrosis nas_tidy"], nas=nas,
                         high_grade=int(nas >= 4), spots=r["spots"], size_mb=r["size_MB"], consent=r["Consent"]))
    assert len(rows) == 99 and len({r["participant"] for r in rows}) == 99
    rows.sort(key=lambda r: r["srr"])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(len(rows), "high_grade", sum(r["high_grade"] for r in rows))


if __name__ == "__main__":
    main()
