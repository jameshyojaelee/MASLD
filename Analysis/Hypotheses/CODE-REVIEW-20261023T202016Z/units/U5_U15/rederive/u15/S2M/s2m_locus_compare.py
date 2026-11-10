"""Per-locus comparison of meSuSiE credible-set counts: locus summary n_cs vs CS present in the
per-variant table (in_cs members), ALT/AST/NAFLD loci. No randomness."""
from pathlib import Path
import pandas as pd
R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/mesusie_mvp")
OUT = Path(__file__).resolve().parent
loc = pd.read_csv(R / "mesusie_locus_summary_mvp.csv")
loc = loc[loc["trait_pair"].isin(["ALT", "AST", "NAFLD"])]
cs = pd.read_csv(OUT / "S2M_credible_sets.csv")
mine = cs.groupby("locus_id").size().rename("n_cs_variant_table")
m = loc.set_index("locus_id")[["n_cs", "max_pip", "n_variants", "converged"]].join(mine).fillna({"n_cs_variant_table": 0})
d = m[m["n_cs"] != m["n_cs_variant_table"]]
print("loci where counts differ:", len(d)); print(d.to_string())
ncols = [c for c in loc.columns if c.startswith("n_cs_") and not c.startswith("n_cs_variant")]
print(loc.set_index("locus_id").loc[d.index, [c for c in ncols if loc.set_index("locus_id").loc[d.index, c].notna().any()]].to_string())
v = pd.read_csv(R / "mesusie_variant_summary_mvp.csv")
v = v[v["locus_id"].isin(d.index)]
print("variant rows at these loci:", len(v), "; cs_index values:", v["cs_index"].value_counts().to_dict(),
      "; in_cs:", v["in_cs"].value_counts().to_dict())
print(v.sort_values("pip", ascending=False).groupby("locus_id").head(3)[["locus_id", "variant_id", "pip", "pip_EUR", "cs_index", "cs_category", "in_cs"]].to_string(index=False))
