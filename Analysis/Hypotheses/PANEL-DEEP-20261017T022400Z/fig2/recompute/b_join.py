"""Fig 2B: join each gene's max-SuSiE-PP.H4 lead to SuSiE PIP from BOTH fine-mapping tables
(combined_finemapping.csv and combined_finemapping_coloc_targeted.csv)."""
import json, pandas as pd, numpy as np
U = "upstream"; O = "out_ABCEF"
bl = pd.read_csv(f"{O}/B_gene_leads.tsv", sep="\t", dtype={"c": str})
bl = bl.drop(columns=[c for c in bl.columns if c in ("study", "chromosome", "n_rec", "susie_pip_max", "susie_pip_min", "susie_pip_clean_max")])
fa = pd.read_csv(f"{U}/combined_finemapping_at_leads.tsv.gz", sep="\t", dtype={"chromosome": str}); fa["src"] = "combined"
fb = pd.read_csv(f"{U}/combined_finemapping_coloc_targeted_at_leads.tsv.gz", sep="\t", dtype={"chromosome": str}); fb["src"] = "coloc_targeted"
fm = pd.concat([fa, fb], ignore_index=True)
fm["p"] = fm.position.astype(int)
res = {}
for name, f in [("combined_only", fa.assign(p=fa.position.astype(int))), ("targeted_only", fb.assign(p=fb.position.astype(int))), ("both", fm)]:
    agg = f.groupby(["study", "chromosome", "p"]).agg(n_rec=("susie_pip", "size"), n_finite=("susie_pip", lambda s: s.notna().sum()),
                                                      pip_max=("susie_pip", "max"), pip_min=("susie_pip", "min"),
                                                      rec_max=("recommended_pip", "max")).reset_index()
    j = bl.merge(agg, left_on=["gwas_name", "c", "p"], right_on=["study", "chromosome", "p"], how="left")
    sh = j[j.pip_max.notna() & ((j.s4 > 0.05) | (j.pip_max >= 0.99))]
    shmin = j[j.pip_min.notna() & ((j.s4 > 0.05) | (j.pip_min >= 0.99))]
    res[name] = {"with_record": int(j.n_rec.notna().sum()), "with_pip": int(j.pip_max.notna().sum()),
                 "multi_record_genes": int((j.n_rec > 1).sum()),
                 "shown": len(sh), "corner_pip>=0.9_pp4<0.5": int(((sh.pip_max >= 0.9) & (sh.s4 < 0.5)).sum()),
                 "topright_pip>=0.9_pp4>0.5": int(((sh.pip_max >= 0.9) & (sh.s4 > 0.5)).sum()),
                 "pp4>0.5_pip<0.1": int(((sh.pip_max < 0.1) & (sh.s4 > 0.5)).sum()),
                 "shown_min": len(shmin), "corner_min": int(((shmin.pip_min >= 0.9) & (shmin.s4 < 0.5)).sum()),
                 "topright_min": int(((shmin.pip_min >= 0.9) & (shmin.s4 > 0.5)).sum())}
    if name == "both":
        j.to_csv(f"{O}/B_gene_leads_pip_both.tsv", sep="\t", index=False)
        sh.to_csv(f"{O}/B_shown_both.tsv", sep="\t", index=False)
        cm = sh[(sh.pip_max >= 0.9) & (sh.s4 < 0.5)]
        res[name]["corner_genes"] = sorted(cm.gene.tolist())
        diff = shmin[(shmin.pip_min >= 0.9) & (shmin.s4 < 0.5)]
        res[name]["corner_genes_min"] = sorted(diff.gene.tolist())
        res[name]["dropped_if_min"] = sorted(set(sh.gene) - set(shmin.gene))
        res[name]["topright_genes"] = sorted(sh[(sh.pip_max >= 0.9) & (sh.s4 > 0.5)].gene.tolist())
        for gname in ["RORA", "HSD17B13", "PNPLA3", "TM6SF2", "MTARC1", "GCKR", "THBS3", "POLL", "APOE", "PNPLA2", "SAMM50"]:
            r = j[j.gene == gname]
            if len(r):
                r = r.iloc[0]
                res[name][gname] = {"gwas": r.gwas_name, "lead": r.top_snp, "pp4": round(float(r.s4), 4),
                                    "pip_max": None if pd.isna(r.pip_max) else round(float(r.pip_max), 4),
                                    "pip_min": None if pd.isna(r.pip_min) else round(float(r.pip_min), 4), "n_rec": None if pd.isna(r.n_rec) else int(r.n_rec)}
            else:
                res[name][gname] = "no finite SuSiE PP.H4 in tier-1/2"
json.dump(res, open(f"{O}/B_join_summary.json", "w"), indent=1)
print(json.dumps(res, indent=1))
