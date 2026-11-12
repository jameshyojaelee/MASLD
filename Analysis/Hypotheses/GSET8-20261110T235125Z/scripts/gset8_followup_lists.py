#!/usr/bin/env python
"""GSET8 follow-up gene lists: locus exclusion, top-gene removal from hep_down, hepatocyte-tested background."""
import argparse

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SYM = f"{ROOT}/tools/pops/example/data/utils/gene_annot_jun10.txt"
PART = f"{ROOT}/Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient/final_partition_posthoc.tsv.gz"
STOUF = f"{ROOT}/Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient/hep_corrected_stouffer.tsv.gz"
LOCI = {"PNPLA3": "ENSG00000100344", "TM6SF2": "ENSG00000213996", "HSD17B13": "ENSG00000170509",
        "APOE": "ENSG00000130203", "GPAM": "ENSG00000119927", "MTARC1": "ENSG00000186205"}
WIN = 1_000_000
HLA = ("6", 28_477_797, 33_448_354)  # extended MHC, GRCh37


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prim", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ann = pd.read_csv(SYM, sep="\t", dtype={"CHR": str})
    rows = []
    for name, eid in LOCI.items():
        r = ann[ann.ENSGID == eid].iloc[0]
        hit = ann[(ann.CHR == r.CHR) & (ann.END >= r.START - WIN) & (ann.START <= r.END + WIN)]
        rows += [dict(gene=g, reason=f"within 1 Mb of {name}") for g in hit.ENSGID]
    hit = ann[(ann.CHR == HLA[0]) & (ann.END >= HLA[1]) & (ann.START <= HLA[2])]
    rows += [dict(gene=g, reason="extended MHC chr6:28.48-33.45 Mb") for g in hit.ENSGID]
    ex = pd.DataFrame(rows).drop_duplicates("gene")
    ex.merge(ann[["ENSGID", "NAME"]], left_on="gene", right_on="ENSGID").drop(columns="ENSGID") \
        .to_csv(f"{a.out}/exclude_loci_genes_annotated.tsv", sep="\t", index=False)
    ex["gene"].to_csv(f"{a.out}/exclude_loci_genes.txt", index=False, header=False)

    memb = pd.read_csv(f"{a.prim}/set_membership.tsv", sep="\t")
    memb["excluded_by_loci"] = memb.gene.isin(set(ex.gene))
    memb.to_csv(f"{a.out}/set_membership_loci_flag.tsv", sep="\t", index=False)

    # top hep_down genes by Ghouse cirrhosis gene z
    g = pd.read_csv(f"{a.prim}/cirrhosis.genes.out", sep=r"\s+")
    hd = memb[memb.set == "hep_down"].merge(g, left_on="gene", right_on="GENE")
    hd = hd.sort_values("ZSTAT", ascending=False)
    keep = ["gene", "symbol", "CHR", "START", "NSNPS", "ZSTAT", "P", "excluded_by_loci"]
    hd[keep].head(20).to_csv(f"{a.out}/hep_down_top20_cirrhosis.tsv", sep="\t", index=False)
    order = list(hd.gene)
    lines = [f"hep_down_minus_top{k} " + " ".join(order[k:]) for k in (1, 5, 10, 20)]
    sig = set(hd.loc[hd.P < 2.5e-6, "gene"])
    lines.append("hep_down_minus_geneP2.5e-6 " + " ".join(x for x in order if x not in sig))
    with open(f"{a.out}/sets_drop_top.txt", "w") as f:
        f.write("\n".join(lines) + "\n")

    # hepatocyte-tested background: the stage genes CONF1 C3b tested in hepatocyte pseudobulk
    part = pd.read_csv(PART, sep="\t")
    st = pd.read_csv(STOUF, sep="\t")
    sym2id = dict(zip(part.symbol, part.ensembl_id))
    st["gene"] = st.symbol.map(sym2id)
    assert st.gene.notna().all()
    st["gene"].to_csv(f"{a.out}/include_heptested.txt", index=False, header=False)
    down = set(part.loc[part.dir == "down", "ensembl_id"])
    st.loc[st.gene.isin(down), "gene"].to_csv(f"{a.out}/include_heptested_down.txt", index=False, header=False)
    print("excluded loci genes:", len(ex), "; hep_down genes excluded:",
          int(memb[(memb.set == "hep_down")].excluded_by_loci.sum()),
          "; heptested:", len(st), "; heptested down:", int(st.gene.isin(down).sum()),
          "; hep_down gene P<2.5e-6:", len(sig))


if __name__ == "__main__":
    main()
