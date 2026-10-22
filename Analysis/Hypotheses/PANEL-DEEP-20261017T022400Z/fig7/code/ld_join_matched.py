"""LD-block join (7D) on matched rows: the held-out gene mean exists only for genes with training variants in
two or more LD blocks; compare left-in and held-out Pearson on that same row set."""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
FIX = EX / "variant_program_substrate_rebuild-20260830T160000Z"
full = pd.read_parquet(FIX / "tables/training_matrix.parquet", columns=["variant_key", "ENSG", "beta_alt"]).merge(
    pd.read_csv(FIX / "splits/variant_splits.tsv", sep="\t")[["variant_key", "block", "ld_fold"]], on="variant_key", how="inner")
x = full[(full.ld_fold != 4) & full.ENSG.notna()].copy()
gs = x.groupby("ENSG").beta_alt.agg(["sum", "count"]); gb = x.groupby(["ENSG", "block"]).beta_alt.agg(["sum", "count"]).rename(columns={"sum": "bs", "count": "bc"})
x = x.join(gs, on="ENSG").join(gb, on=["ENSG", "block"])
x["left_in"] = x["sum"] / x["count"]; x["held"] = (x["sum"] - x.bs) / (x["count"] - x.bc)
ok = np.isfinite(x.held)
m = x[ok]
res = {"rows_all": int(len(x)), "rows_matched": int(ok.sum()), "genes_all": int(x.ENSG.nunique()), "genes_matched": int(m.ENSG.nunique()),
       "left_in_all": float(np.corrcoef(x.left_in, x.beta_alt)[0, 1]),
       "left_in_matched": float(np.corrcoef(m.left_in, m.beta_alt)[0, 1]),
       "held_out_matched": float(np.corrcoef(m.held, m.beta_alt)[0, 1]),
       "own_block_mean_matched": float(np.corrcoef(m.bs / m.bc, m.beta_alt)[0, 1]),
       "frac_rows_singleton_variant_gene": float((x["count"] == 1).mean())}
json.dump(res, open(sys.argv[1] + "/ld_join_matched.json", "w"), indent=1); print(res)
