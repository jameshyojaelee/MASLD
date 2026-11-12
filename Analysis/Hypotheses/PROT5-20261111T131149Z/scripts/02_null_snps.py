#!/usr/bin/env python3
"""PROT5 step 2: fixed null SNP set and extraction windows (no pQTL statistic is read).

Null SNPs: 100,000 1000G EUR autosomal non-palindromic SNVs drawn with seed 20261111, excluding
chr6:25-34 Mb, +/- 1 Mb of the 169 ALT leads and two control variants, +/- 1 Mb of each panel gene.
Windows: +/- 500 kb of the 169 leads and the control variants.
Usage: python -I 02_null_snps.py <project_root> <loci.tsv> <panels.tsv> <out_dir>
"""
import sys

import numpy as np
import pandas as pd

CONTROLS = [("ctrl_PNPLA3_rs738409", 22, 44324727, "G", "C", "G"), ("ctrl_GOT1_lead", 10, 101152720, "T", "G", "T")]


def main(root, loci_path, panels_path, out):
    loci = pd.read_csv(loci_path, sep="\t", dtype={"locus": str})[["locus", "chr", "pos", "a1", "a2", "alt_raising_allele"]]
    ctrl = pd.DataFrame(CONTROLS, columns=["locus", "chr", "pos", "a1", "a2", "alt_raising_allele"])
    leads = pd.concat([loci, ctrl], ignore_index=True)
    leads.to_csv(f"{out}/leads.tsv", sep="\t", index=False)
    w = leads.assign(s=leads["pos"] - 500_000, e=leads["pos"] + 500_000)[["chr", "s", "e"]]
    w.to_csv(f"{out}/windows.tsv", sep="\t", index=False, header=False)
    pan = pd.read_csv(panels_path, sep="\t", dtype={"chr": str})
    pan = pan[pan["panel"].isin(["hepatocyte", "negative_control"])]
    excl = [(int(r.chr), r.pos - 1_000_000, r.pos + 1_000_000) for r in leads.itertuples()]
    excl += [(int(r.chr), r.start - 1_000_000, r.end + 1_000_000) for r in pan.itertuples() if r.chr.isdigit()]
    excl.append((6, 25_000_000, 34_000_000))
    rows = []
    for c in range(1, 23):
        b = pd.read_csv(f"{root}/GWAS/finemapping/data/ld_ref/1kg_eur/chr{c}_eur.bim", sep="\t", header=None,
                        names=["chr", "id", "cm", "pos", "a1", "a2"], usecols=["chr", "pos", "a1", "a2"])
        b = b[b["a1"].isin(list("ACGT")) & b["a2"].isin(list("ACGT"))]
        pal = ((b["a1"] + b["a2"]).isin(["AT", "TA", "CG", "GC"]))
        b = b[~pal & ~b["pos"].duplicated(keep=False)]
        keep = np.ones(len(b), bool)
        for ec, s, e in excl:
            if ec == c:
                keep &= ~((b["pos"].values >= s) & (b["pos"].values <= e))
        rows.append(b[keep])
    allb = pd.concat(rows, ignore_index=True)
    rng = np.random.default_rng(20261111)
    idx = np.sort(rng.choice(len(allb), size=100_000, replace=False))
    null = allb.iloc[idx]
    null[["chr", "pos", "a1", "a2"]].to_csv(f"{out}/null_snps.tsv", sep="\t", index=False, header=False)
    print(f"eligible null SNPs {len(allb)}; drawn {len(null)}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
