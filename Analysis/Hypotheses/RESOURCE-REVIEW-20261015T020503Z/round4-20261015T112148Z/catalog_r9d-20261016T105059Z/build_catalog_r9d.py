#!/usr/bin/env python3
"""Gene Catalog r9d: the adopted r9 Catalog (user approval 2026-10-16, items H, I, R, S).
r9d = r9c candidate tables unchanged + three columns from H79 and H81 + release_state and
blocked_on updated for the approved items. Reads saved tables only; no random draws.
Usage: python build_catalog_r9d.py <empty output dir>"""
import sys, hashlib, shutil, subprocess
from pathlib import Path
import numpy as np, pandas as pd
np.random.seed(20260923)
OUT = Path(sys.argv[1])
R4 = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z")
R9C = R4 / "catalog_r9c_candidate-20261015T191932Z/build_out"
HYP = R4.parents[1]
H79 = HYP / "H79-20261015T125125Z/results/run-20261015T132031Z/D1_flagged_genes.tsv"
H79_RES = HYP / "H79-20261015T125125Z/RESULT_H79.md"
H81 = HYP / "H81-20261015T125233Z/run-20261015T141242Z/summary/D5_detection_limit_lookup_not_supported_genes.tsv.gz"
H81_RES = HYP / "H81-20261015T125233Z/run-20261015T141242Z/RESULT_H81.md"
RS = "r9d adopted 2026-10-16"
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
W = pd.read_parquet(R9C / "catalog_r9c_candidate_gene_table.parquet")
ids = W.ensembl_id.astype(str)
assert ids.is_unique and len(W) == 28233
# H79: ancestry-sensitive flag (D1_rule_M1); EUR-only flags (D1_rule_ME) are not used
f = pd.read_csv(H79, sep="\t", float_precision="round_trip")
f["eid"] = f.gene.astype(str).str.split(".").str[0]
m1 = f[f.D1_rule_M1.astype(str).str.upper() == "TRUE"].set_index("eid")
assert len(m1) == 64 and m1.M1_status.value_counts().to_dict() == {"lost": 40, "gained": 24}
has_rna = W["rna_log2fc"].notna()
W["rna_ancestry_sensitive_h79"] = pd.Series(np.where(ids.isin(m1.index), True, False), index=W.index).astype("boolean").where(has_rna)
W["rna_ancestry_change_h79"] = ids.map(m1.M1_status).fillna("").where(has_rna, "")
# H81: per-decile detection limit, genes without a supported bulk call
d = pd.read_csv(H81, sep="\t", float_precision="round_trip")
d["eid"] = d.gene.astype(str).str.split(".").str[0]
assert d.eid.is_unique
d = d.set_index("eid")
W["rna_expression_decile_h81"] = ids.map(d.decile).astype("Int64")
W["rna_detection_limit_log2fc_h81"] = ids.map(d.detection_limit_log2FC)
W["rna_below_detection_limit_h81"] = W.rna_expression_decile_h81.isin([1, 2]).astype("boolean").where(W.rna_expression_decile_h81.notna())
# dictionary
cd = pd.read_csv(R9C / "column_dictionary.tsv", sep="\t", keep_default_na=False)
cd["table"] = cd["table"].replace("catalog_r9c_candidate_gene_table", "catalog_r9d_gene_table")
clear = cd.blocked_on.str.contains("H79|corrected 50-study COLOC rerun|H51", regex=True)
cd.loc[clear, "caveat"] = (cd.loc[clear, "caveat"] + " ").str.strip() + np.where(cd.loc[clear, "blocked_on"].str.contains("COLOC"), " Block cleared: COLOC r2 adopted 2026-10-14.", " Block cleared: user approval 2026-10-16.")
cd.loc[clear, "caveat"] = cd.loc[clear, "caveat"].str.strip()
cd.loc[clear, "blocked_on"] = ""
cd.loc[cd["release_state"] == "candidate", "release_state"] = RS
new = [("rna_ancestry_sensitive_h79", "boolean", "TRUE for the 64 genes meeting the H79 ancestry-sensitive rule (M1: four ancestry PCs); NA without a bulk test", H79, H79_RES, "Call count is ancestry-sensitive; log2FC agreement over adopted calls Pearson 0.9972. EUR-only (ME) flags are not used."),
       ("rna_ancestry_change_h79", "lost | gained | empty", "Direction of the call change under M1 for flagged genes (40 lost, 24 gained)", H79, H79_RES, ""),
       ("rna_expression_decile_h81", "integer 1-10", "Expression decile used by H81, genes without a supported bulk call", H81, H81_RES, "Spike-in calibration on fixed 774 units; rates describe the recipe on these units"),
       ("rna_detection_limit_log2fc_h81", "log2 fold change", "Smallest uniform log2FC detected with 80% recall in the gene's decile (H81 D5 lookup)", H81, H81_RES, "Lookup only; no gene is relabelled"),
       ("rna_below_detection_limit_h81", "boolean", "TRUE in deciles 1-2, where no call does not exclude a 2-fold change (H81 Reading)", H81, H81_RES, "")]
rows = [{"table": "catalog_r9d_gene_table", "column": c, "source": "H79 (2026-10-16)" if "h79" in c else "H81 (2026-10-16)", "run_dir": str(sf.parents[2] if "h79" in c else sf.parents[2]),
         "result": str(rf), "source_file": str(sf), "unit": u, "allowed_values": "", "release_state": RS, "blocked_on": "", "caveat": cv} for c, u, desc, sf, rf, cv in new]
for r, (_, _, desc, _, _, _) in zip(rows, new): r["unit"] = r["unit"] + "; " + desc
cd = pd.concat([cd, pd.DataFrame(rows)], ignore_index=True)
bl = pd.read_csv(R9C / "blocked_on.tsv", sep="\t", keep_default_na=False)
keep = ~bl.blocked_on.str.contains("H79|corrected 50-study COLOC rerun|user approval of the H51", regex=True)
bl = bl[keep]
# write
W.to_parquet(OUT / "catalog_r9d_gene_table.parquet", index=False)
W.to_csv(OUT / "catalog_r9d_gene_table.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
cd.to_csv(OUT / "column_dictionary.tsv", sep="\t", index=False)
bl.to_csv(OUT / "blocked_on.tsv", sep="\t", index=False)
for f_ in R9C.glob("sidecar_*"): shutil.copy2(f_, OUT / f_.name)
for f_ in ["catalog_rules_r9.tsv", "universe_counts.tsv", "invariants.tsv"]: shutil.copy2(R9C / f_, OUT / f_)
inp = [R9C / "catalog_r9c_candidate_gene_table.parquet", R9C / "column_dictionary.tsv", R9C / "blocked_on.tsv", H79, H81]
pd.DataFrame({"file": [str(p) for p in inp], "sha256": [sha(p) for p in inp]}).to_csv(OUT / "inputs_sha256.tsv", sep="\t", index=False)
(OUT / "pip_freeze.txt").write_text(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout)
print("rows", len(W), "cols", W.shape[1], "flag", int(W.rna_ancestry_sensitive_h79.sum()), "below_DL", int(W.rna_below_detection_limit_h81.sum()), "blocked rows", len(bl))
