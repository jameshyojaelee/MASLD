#!/usr/bin/env python3
"""AN1 validation: compare share-only rebuild with the adopted r2 signal-pair table.

Exact-match criteria (all must hold for every validation row present in the run):
  1. every (study, gene, eQTL CS idx2) in the r2 pairs file is rebuilt, and no extra CS appears;
  2. |prop_eqtl(rebuilt) - prop_eqtl(r2)| <= 1e-9 for every CS (reported max);
  3. per-pair eligibility identical (pair eligible = prop_gwas >= 0.5 and prop_eqtl >= 0.5,
     with prop_gwas taken from the r2 file, where it is 1 for every pair);
  4. per-row n_cs_pairs_eligible (master) = n_gwas_cs(r2 pairs) x n_eqtl_cs_eligible(rebuilt);
  5. per-row testability identical: method == susie  <=>  >= 1 eligible pair;
  6. n_snps (rows merged) identical to the master n_snps.
Usage: an1_validate.py <run_out_dir> [<run_out_dir> ...]
"""
import sys, json
import numpy as np
import pandas as pd

R2 = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
outs = sys.argv[1:]
rows = pd.concat([pd.read_csv(f"{o}/row_shares.tsv", sep="\t", dtype={"ensembl": str}) for o in outs])
cs = pd.concat([pd.read_csv(f"{o}/eqtl_cs_shares.tsv", sep="\t", dtype={"ensembl": str}) for o in outs])
pairs = pd.read_csv(f"{R2}/susie_coloc_signal_pairs.tsv.gz", sep="\t", dtype={"ensembl": str})
pairs = pairs.rename(columns={"gwas": "gwas_name"})
key = ["gwas_name", "ensembl"]
pr = pairs.merge(rows[key], on=key)            # pairs of the rows in this run
res = {"n_rows": int(len(rows)), "n_pairs_r2": int(len(pr)),
       "row_status": rows.status.value_counts(dropna=False).to_dict()}

# 1-2: per eQTL CS
r2cs = pr.groupby(key + ["idx2"]).agg(prop_eqtl_r2=("prop_eqtl", "first"),
                                     nuniq=("prop_eqtl", "nunique")).reset_index()
assert (r2cs.nuniq == 1).all(), "r2 prop_eqtl not constant within an eQTL CS"
m = r2cs.merge(cs[key + ["idx2", "prop_eqtl", "eligible_cs"]], on=key + ["idx2"], how="outer", indicator=True)
res["cs_only_r2"] = int((m._merge == "left_only").sum())
res["cs_only_rebuilt"] = int((m._merge == "right_only").sum())
both = m[m._merge == "both"]
d = (both.prop_eqtl - both.prop_eqtl_r2).abs()
res["n_cs_compared"] = int(len(both))
res["max_abs_diff_prop_eqtl"] = float(d.max()) if len(d) else None
res["n_cs_diff_gt_1e-9"] = int((d > 1e-9).sum())
res["n_cs_eligibility_differs"] = int(((both.prop_eqtl_r2 >= 0.5) != (both.prop_eqtl >= 0.5)).sum())

# 3: per pair
pp = pr.merge(cs[key + ["idx2", "prop_eqtl"]].rename(columns={"prop_eqtl": "prop_eqtl_rb"}),
              on=key + ["idx2"], how="left")
pp["elig_rb"] = ~((pp.prop_gwas < 0.5) | (pp.prop_eqtl_rb < 0.5))
res["n_pairs_compared"] = int(pp.prop_eqtl_rb.notna().sum())
res["n_pairs_eligibility_differs"] = int((pp.elig_rb != pp.eligible).sum())
res["prop_gwas_min_r2"] = float(pr.prop_gwas.min())

# 4-6: per row
ng = pr.groupby(key).idx1.nunique().rename("n_gwas_cs_r2").reset_index()
rr = rows.merge(ng, on=key, how="left")
rr["pairs_elig_rb"] = rr.n_gwas_cs_r2 * rr.n_eqtl_cs_eligible
rr["testable_rb"] = rr.n_eqtl_cs_eligible > 0
rr["testable_r2"] = rr.method == "susie"
res["n_rows_n_pairs_eligible_differs"] = int((rr.pairs_elig_rb != rr.n_cs_pairs_eligible).sum())
res["n_rows_testability_differs"] = int((rr.testable_rb != rr.testable_r2).sum())
res["n_rows_n_snps_differs"] = int((rr.n_snps_rebuilt != rr.n_snps).sum())
res["n_rows_n_cs_pairs_differs"] = int((rr.n_gwas_cs_r2 * rr.n_eqtl_cs != rr.n_cs_pairs).sum())
res["n_rows_untestable_r2"] = int((~rr.testable_r2).sum())
res["n_rows_untestable_rebuilt"] = int((~rr.testable_rb).sum())
res["gwas_ids_subset_of_eqtl_false"] = int((rows.gwas_ids_subset_of_eqtl == False).sum())
res["ld_route"] = rows.ld_route.value_counts(dropna=False).to_dict()
res["EXACT_MATCH"] = bool(res["cs_only_r2"] == 0 and res["cs_only_rebuilt"] == 0 and
                          res["n_cs_diff_gt_1e-9"] == 0 and res["n_cs_eligibility_differs"] == 0 and
                          res["n_pairs_eligibility_differs"] == 0 and res["n_rows_n_pairs_eligible_differs"] == 0 and
                          res["n_rows_testability_differs"] == 0 and res["n_rows_n_snps_differs"] == 0 and
                          res["n_rows_n_cs_pairs_differs"] == 0 and len(rows) == rows.status.notna().sum())
print(json.dumps(res, indent=1, default=str))
bad = rr[(rr.pairs_elig_rb != rr.n_cs_pairs_eligible) | (rr.testable_rb != rr.testable_r2) |
         (rr.n_snps_rebuilt != rr.n_snps)]
if len(bad):
    print(bad[key + ["chr", "status", "n_snps", "n_snps_rebuilt", "n_cs_pairs", "n_cs_pairs_eligible",
                     "n_gwas_cs_r2", "n_eqtl_cs", "n_eqtl_cs_eligible", "method"]].head(40).to_string())
