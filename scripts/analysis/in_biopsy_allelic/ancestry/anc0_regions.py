#!/usr/bin/env python3
"""Ancestry step 0: site lists and the 1kGP founder list for the panel extraction.

Inputs
  --targets   A1 identity targets (chrom, pos, "REF,ALT"; bgzipped, no header).
  --tags      tags_final.tsv (tag chrom/pos/ref/alt and GTEx lead_variant_id).
  --ped       1kGP 3,202-sample ped file.
Outputs (in --out)
  founders.txt            1kGP samples with no parent listed (FatherID = MotherID = 0).
  regions.<chrom>.tsv     sorted (chrom, pos) of identity sites, tag SNVs and tag leads.
  tag_regions.<chrom>.tsv sorted (chrom, pos) of tag SNVs only (anc1 pulls every panel record whose
                          span covers one, to count other alleles and deletions at the tag base).
  identity_sites.tsv      autosomal identity targets as chrom, pos, ref, alt, id (id = chr:pos:ref:alt
                          with the chromosome number, matching plink2 '@:#:$r:$a').
Rules: autosomes only (chr1-22). Reads no MASLD genotype, allele count or outcome.
"""
import argparse
from pathlib import Path

import pandas as pd

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def lead_sites(tags):
    """GTEx variant ids 'chr1_827221_T_C_b38' -> (chrom, pos)."""
    parts = tags["lead_variant_id"].str.split("_", expand=True)
    return pd.DataFrame({"chrom": parts[0], "pos": parts[1].astype(int)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--tags", required=True)
    ap.add_argument("--ped", required=True)
    ap.add_argument("--chroms", default=",".join(AUTOSOMES), help="comma list, e.g. chr20,chr21,chr22 for a smoke test")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    chroms = a.chroms.split(",")

    ped = pd.read_csv(a.ped, sep=r"\s+", dtype=str, keep_default_na=False)
    founders = ped[(ped["FatherID"] == "0") & (ped["MotherID"] == "0")]["SampleID"]
    founders.to_csv(out / "founders.txt", index=False, header=False)

    tg = pd.read_csv(a.targets, sep="\t", header=None, names=["chrom", "pos", "alleles"], dtype={"pos": int})
    tg[["ref", "alt"]] = tg["alleles"].str.split(",", expand=True)
    tg = tg[tg["chrom"].isin(chroms)].drop(columns="alleles")
    tg["id"] = tg["chrom"].str.replace("chr", "", regex=False) + ":" + tg["pos"].astype(str) + ":" + tg["ref"] + ":" + tg["alt"]
    tg.to_csv(out / "identity_sites.tsv", sep="\t", index=False)

    tags = pd.read_csv(a.tags, sep="\t", dtype={"pos": int}, keep_default_na=False)
    sites = pd.concat([tg[["chrom", "pos"]], tags[["chrom", "pos"]], lead_sites(tags)], ignore_index=True)
    sites = sites[sites["chrom"].isin(chroms)].drop_duplicates().sort_values(["chrom", "pos"])
    tag_sites = tags.loc[tags["chrom"].isin(chroms), ["chrom", "pos"]].drop_duplicates().sort_values(["chrom", "pos"])
    for c in chroms:
        sites[sites["chrom"] == c].to_csv(out / f"regions.{c}.tsv", sep="\t", index=False, header=False)
        tag_sites[tag_sites["chrom"] == c].to_csv(out / f"tag_regions.{c}.tsv", sep="\t", index=False, header=False)
    print(f"founders {len(founders)}; identity sites {len(tg)}; region positions {len(sites)} on {len(chroms)} chromosomes")


if __name__ == "__main__":
    main()
