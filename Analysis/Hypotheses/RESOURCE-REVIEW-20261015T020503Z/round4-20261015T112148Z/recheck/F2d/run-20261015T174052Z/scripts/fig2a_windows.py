import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(__file__).resolve().parent.parent
tier = pd.read_csv(ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
main = set(tier[tier.placement == "main"].study_name)
cs = pd.read_csv(ROOT / "GWAS/finemapping/results/credible_sets.csv", dtype={"locus": str, "chromosome": str}, low_memory=False)
cs = cs[cs.study.isin(main)]
res = {}
for lab, sub in [("all_rows", cs), ("cs_rows", cs[cs.susie_cs.notna()])]:
    w = sub.groupby(["study", "chromosome", "locus"]).position.agg(["min", "max"]).reset_index()
    res[f"windows_{lab}"] = len(w)
    res[f"max_span_{lab}"] = int((w["max"] - w["min"]).max())
    res[f"ukbb_windows_{lab}"] = w[w.study.str.startswith("UKBB_")].groupby("study").size().to_dict()
    res[f"strata_with_cs_{lab}"] = int(w.study.nunique())
    c = pd.read_csv(ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
                    usecols=["gwas_name", "gene", "PP.H4.susie", "top_snp"], low_memory=False) if lab == "all_rows" else c
    pos = c[c.gwas_name.isin(main) & c.gene.notna() & (c["PP.H4.susie"] > 0.5)].copy()
    pos["chrom"] = pos.top_snp.str.split(":").str[0]
    pos["bp"] = pos.top_snp.str.split(":").str[1].astype(float)
    for pad in [0, 250_000]:
        m = pos.merge(w, left_on=["gwas_name", "chrom"], right_on=["study", "chromosome"], how="left")
        m["hit"] = (m.bp >= m["min"] - pad) & (m.bp <= m["max"] + pad)
        genes = m[m.hit].gene.nunique()
        res[f"genes_same_study_pad{pad}_{lab}"] = int(genes)
    res["n_genes"] = int(pos.gene.nunique()); res["n_pos_pairs"] = len(pos)
    nocs = main - set(w.study)
    res[f"main_without_cs_{lab}"] = sorted(nocs)
    res[f"genes_only_nocs_{lab}"] = int(pos.groupby("gene").gwas_name.apply(lambda s: set(s) <= nocs).sum())
json.dump(res, open(OUT / "fig2a_windows.json", "w"), indent=1)
print(json.dumps(res, indent=1))
