#!/usr/bin/env python
"""GSET8: gene sets, size-matched random sets (seed 20261110) and AveExpr covariate file.

Reads the MAGMA genes.out of both traits only for the list of scored gene IDs (no statistics).
"""
import argparse
import os

import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PART = f"{ROOT}/Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient/final_partition_posthoc.tsv.gz"
SYM = f"{ROOT}/tools/pops/example/data/utils/gene_annot_jun10.txt"
CG0G = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv"
SEED = 20261110
N_RANDOM = 1000
TRAITS = ["pdff", "cirrhosis"]


def scored_genes(path):
    return set(pd.read_csv(path, sep=r"\s+", usecols=["GENE"])["GENE"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = ap.parse_args().out

    part = pd.read_csv(PART, sep="\t")
    hep = part[part["final"] == "within_hepatocyte"]
    counts = hep["dir"].value_counts().to_dict()
    assert counts == {"up": 689, "down": 466}, counts
    assert not hep["ensembl_id"].duplicated().any()

    sym = pd.read_csv(SYM, sep="\t")
    cg0 = pd.read_csv(CG0G, sep="\t", usecols=["g", "gene_name", "AveExpr"])
    # DIG1 Complex I rule (27_replicate_ci_pattern.R line 44); duplicated symbols keep the CG0 bulk ID
    ci = sym[sym["NAME"].str.match(r"^NDUF[ABSCV][0-9]") & (sym["NAME"] != "NDUFA4")].copy()
    dup = ci["NAME"].duplicated(keep=False)
    ci = ci[~dup | ci["ENSGID"].isin(set(cg0["g"]))]
    assert not ci["NAME"].duplicated().any() and len(ci) == 35, ci

    sets = {
        "hep_up": sorted(hep.loc[hep["dir"] == "up", "ensembl_id"]),
        "hep_down": sorted(hep.loc[hep["dir"] == "down", "ensembl_id"]),
        "complexI": sorted(ci["ENSGID"]),
    }
    scored = {t: scored_genes(f"{out}/{t}.genes.out") for t in TRAITS}
    both = sorted(scored["pdff"] & scored["cirrhosis"])

    rows = []
    for name, genes in sets.items():
        for g in genes:
            rows.append(dict(set=name, gene=g, scored_pdff=g in scored["pdff"],
                             scored_cirrhosis=g in scored["cirrhosis"]))
    memb = pd.DataFrame(rows).merge(sym[["ENSGID", "NAME"]].rename(columns={"ENSGID": "gene", "NAME": "symbol"}),
                                    on="gene", how="left")
    memb.to_csv(f"{out}/set_membership.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    universe = np.array(both)
    lines_primary = [name + " " + " ".join(genes) for name, genes in sets.items()]
    lines_random = []
    sizes = {}
    for name, genes in sets.items():
        n = len(set(genes) & set(both))
        sizes[name] = n
        for i in range(N_RANDOM):
            pick = rng.choice(universe, size=n, replace=False)
            lines_random.append(f"rnd_{name}_{i:04d} " + " ".join(sorted(pick)))
    with open(f"{out}/sets_primary.txt", "w") as f:
        f.write("\n".join(lines_primary) + "\n")
    with open(f"{out}/sets_with_random.txt", "w") as f:
        f.write("\n".join(lines_primary + lines_random) + "\n")

    cov = cg0.dropna(subset=["AveExpr"]).drop_duplicates("g")
    cov[["g", "AveExpr"]].rename(columns={"g": "GENE"}).to_csv(f"{out}/covar_aveexpr.txt", sep=" ", index=False)

    with open(f"{out}/set_sizes.tsv", "w") as f:
        f.write("set\tn_listed\tn_in_annotation_scored_pdff\tn_scored_cirrhosis\tn_scored_both_random_size\n")
        for name, genes in sets.items():
            f.write(f"{name}\t{len(genes)}\t{len(set(genes) & scored['pdff'])}\t"
                    f"{len(set(genes) & scored['cirrhosis'])}\t{sizes[name]}\n")
    print(open(f"{out}/set_sizes.tsv").read())
    print("universe scored in both traits:", len(both), "; genes with AveExpr:", len(cov))


if __name__ == "__main__":
    main()
