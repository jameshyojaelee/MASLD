"""Figure S2J plotted rows: gene maximum of the meSuSiE EUR-specific component PIP and of the 11
shared components over ALT/AST/NAFLD loci, variants linked to genes by Broadaway liver cis-eQTL
pairs (canon_eqtl_pairs.tsv, either allele orientation). No randomness."""
import itertools
from pathlib import Path
import pandas as pd
R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(__file__).resolve().parent
ANC = ["EUR", "EAS", "AFR", "AMR"]
shared = ["pip_" + "_".join(c) for k in (2, 3, 4) for c in itertools.combinations(ANC, k)]
mv = pd.read_csv(R / "GWAS/finemapping/results/mesusie_mvp/mesusie_variant_summary_mvp.csv")
mv = mv[mv["trait_pair"].isin(["ALT", "AST", "NAFLD"])].copy()
mv["shared_max"] = mv[shared].max(axis=1)
mv["akey"] = [":".join(sorted([str(a), str(b)])) for a, b in zip(mv["a1"], mv["a2"])]
eq = pd.read_csv(OUT / "canon_eqtl_pairs.tsv", sep="\t", header=None,
                 names=["chr", "pos", "NEA", "EA", "gene", "ENSG"], dtype={"chr": int, "pos": int})
eq["akey"] = [":".join(sorted([str(a), str(b)])) for a, b in zip(eq["NEA"], eq["EA"])]
eq = eq[["chr", "pos", "akey", "gene"]].drop_duplicates()
mv["chr"] = mv["chr"].astype(int); mv["pos"] = mv["pos"].astype(int)
j = mv.merge(eq, on=["chr", "pos", "akey"], how="inner")
print("meSuSiE ALT/AST/NAFLD variant rows:", len(mv), "; linked gene-variant rows (canonical genes):", len(j))
g = j.groupby("gene").agg(eur=("pip_EUR", "max"), shared=("shared_max", "max"),
                          n_var=("variant_id", "nunique"), loci=("locus_id", "nunique"),
                          traits=("trait_pair", lambda s: ";".join(sorted(set(s))))).reset_index()
src = pd.read_csv(R / "figures/candidates/supp-redraw-20261020T0300Z/S2b/panels/FigS2J_crossancestry_pip_concentration_source.csv")
g = g.merge(src[["gene", "eur_pip", "joint_pip", "plotted"]], on="gene", how="left")
g["eur_4dp"] = g["eur"].round(4); g["shared_4dp"] = g["shared"].round(4)
g["match"] = (g["eur_4dp"] == g["eur_pip"]) & (g["shared_4dp"] == g["joint_pip"])
print(g.sort_values("shared", ascending=False).to_string(index=False))
print("canonical genes with no linked meSuSiE variant:",
      sorted(set(open(OUT / "canon_genes.txt").read().split()) - set(g["gene"])))
g.to_csv(OUT / "S2J_canon_gene_maxima.csv", index=False)
