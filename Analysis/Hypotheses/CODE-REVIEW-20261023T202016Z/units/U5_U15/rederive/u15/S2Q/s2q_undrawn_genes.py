"""COLOC maxima (35 strata) for the three Elison genes not drawn in Figure S2Q. No randomness."""
from pathlib import Path
import pandas as pd
U = Path(__file__).resolve().parents[1]
c = pd.read_csv(U / "common/coloc35.csv.gz")
s = c[c["gene"].isin(["CEBPA", "ONECUT1", "XBP1"])]
print(s.groupby(["gene", "tier"]).agg(susie=("PP.H4.susie", "max"), abf=("PP.H4.abf", "max"), n=("gwas_name", "nunique")).to_string())
