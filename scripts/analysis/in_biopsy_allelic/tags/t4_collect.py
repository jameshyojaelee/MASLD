#!/usr/bin/env python3
"""Model A pivot, T4 collect: one long table of allele depths at the final tag SNVs.

Reads the per cohort x chunk `bcftools query` tables written by t4_pileup.sbatch and
keeps, per library and tag SNV, the reads for the tag's REF and for the tag's ALT
(0 when the ALT base was not seen in that cohort). Rows with zero depth are dropped.
Fails if a chunk is missing or if the pileup REF differs from the tag REF. Tag SNVs
with no read in any library of a cohort are absent from that cohort's pileup and
are counted. Output stays in the restricted dir.
"""
import argparse
import re
from pathlib import Path

import pandas as pd

COHORTS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523",
           "GSE174478", "GSE193066", "GSE213621", "GSE240729", "PRJNA512027"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ad-dir", required=True)
    ap.add_argument("--tags", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tags = pd.read_csv(a.tags, sep="\t").drop_duplicates(["chrom", "pos"])[["chrom", "pos", "ref", "alt"]]
    tag_key = {(c, p): (r, x) for c, p, r, x in tags.itertuples(index=False)}
    rows, report = [], []
    for c in COHORTS:
        seen = set()
        for k in range(20):
            f = Path(a.ad_dir) / f"{c}.chunk{k}.ad.tsv.gz"
            if not f.exists():
                raise SystemExit(f"missing {f}")
            d = pd.read_csv(f, sep="\t", dtype=str)
            libs = [re.sub(r"\.bam:AD$", "", col.split("]", 1)[1]).rsplit("/", 1)[-1] for col in d.columns[4:]]
            for rec in d.itertuples(index=False):
                chrom, pos, ref, alts = rec[0], int(rec[1]), rec[2], rec[3].split(",")
                if (chrom, pos) not in tag_key:
                    continue
                # mpileup also writes an indel record at some SNV positions (a deletion has a longer REF,
                # an insertion a longer ALT with a one-base REF); skip it
                if len(ref) != 1 or any(len(x) != 1 for x in alts if x != "<*>"):
                    continue
                tref, talt = tag_key[(chrom, pos)]
                if ref != tref:
                    raise SystemExit(f"REF mismatch at {chrom}:{pos}: pileup {ref}, tag {tref}")
                seen.add((chrom, pos))
                j = 1 + alts.index(talt) if talt in alts else None
                for lib, ad in zip(libs, rec[4:]):
                    v = [int(x) if x != "." else 0 for x in ad.split(",")]
                    r_reads, a_reads = v[0], (v[j] if j is not None else 0)
                    other = sum(v[1:]) - a_reads
                    if r_reads + a_reads + other > 0:
                        rows.append((c, lib, chrom, pos, tref, talt, r_reads, a_reads, other))
        absent = len(tag_key) - len(seen)
        report.append({"cohort": c, "tag_snvs_in_pileup": len(seen), "tag_snvs_absent": absent})
    out = pd.DataFrame(rows, columns=["cohort", "run", "chrom", "pos", "ref", "alt",
                                      "ref_reads", "alt_reads", "other_reads"])
    if out.duplicated(["cohort", "run", "chrom", "pos"]).any():
        raise SystemExit("more than one record per library x tag SNV")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False, compression="gzip")
    rep = pd.DataFrame(report)
    rep["libraries"] = [out.loc[out["cohort"] == c, "run"].nunique() for c in rep["cohort"]]
    rep["rows"] = [int((out["cohort"] == c).sum()) for c in rep["cohort"]]
    print(rep.to_string(index=False))


if __name__ == "__main__":
    main()
