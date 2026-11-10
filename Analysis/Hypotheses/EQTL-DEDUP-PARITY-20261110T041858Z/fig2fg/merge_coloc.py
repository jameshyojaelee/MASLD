# Replace r2 MVP EUR ALT/AST COLOC rows for the 90 refit genes with refit rows; drop MPRIP (untestable by prespec).
import glob, pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = f"{ROOT}/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z"
cols = ["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "top_snp"]
r2 = pd.read_csv(f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv", usecols=cols)
r2 = r2[r2.gwas_name.isin(["MVP_ALT_EUR", "MVP_AST_EUR"])]
cmp = pd.read_csv(f"{RUN}/results/eqtl_dedup_gene_comparison.tsv", sep="\t")
refit_ens = set(cmp.ensembl)
new = pd.concat([pd.read_csv(f, usecols=cols) for g in ["MVP_ALT_EUR", "MVP_AST_EUR"]
                 for f in glob.glob(f"{RUN}/coloc/{g}/susie_coloc_chr*.csv")])
new = new[new.ensembl.str.split(".").str[0].isin(refit_ens)]
r2e = r2.ensembl.str.split(".").str[0]
out = pd.concat([r2[~r2e.isin(refit_ens)], new])
print("r2 rows", len(r2), "replaced", int(r2e.isin(refit_ens).sum()), "refit rows", len(new), "out", len(out))
for g in ["CTSS", "IDI2", "MPRIP", "CDK5RAP3", "IRF5", "ABHD12"]:
    print(g, r2[r2.gene == g][["gwas_name", "PP.H4.abf", "PP.H4.susie"]].values.tolist(), "->",
          out[out.gene == g][["gwas_name", "PP.H4.abf", "PP.H4.susie"]].values.tolist())
out.to_csv("coloc_mvp_eur_alt_ast_refit.csv", index=False)
