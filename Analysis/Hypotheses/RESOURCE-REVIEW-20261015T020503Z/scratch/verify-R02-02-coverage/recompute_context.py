"""Reproduce the adopted noncoding-context estimates, then recompute them after
annotating every member by position (no study-specific lookup gap).

Uses the producer's own functions from scripts/figures/fig2_multiplicity_correction.py
(read-only import) so the estimator, clustering and seeds are identical.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-02-coverage"
sys.path.insert(0, str(ROOT / "scripts/figures"))
import fig2_multiplicity_correction as fm  # noqa: E402

base = ROOT / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
arch = pd.read_csv(base / "credible_set_pip_architecture.tsv", sep="\t", dtype={"locus": str})
members = pd.read_csv(base / "credible_set_members_annotated.tsv", sep="\t",
                      low_memory=False, dtype={"locus": str})
arch = arch[arch.trait_scope.isin(fm.SCOPES)].copy()
members = members[members.trait_scope.isin(fm.SCOPES)].copy()
arch, members, excluded = fm.apply_gate(arch, members)
print("gated member rows", len(members), "instances", members.credible_set_uid.nunique())

sites = pd.read_csv(OUT / "site_context_reannotated.tsv", sep="\t")
sites["chromosome"] = sites["chromosome"].astype(str)
members["chromosome"] = members["chromosome"].astype(str)
members = members.merge(sites, on=["chromosome", "position"], how="left", validate="m:1")
assert members.pos_hg38.notna().all()

# 1. Agreement of position-based re-annotation with recorded values on resolved rows
res = members[members.promoter_context_status == "resolved"]
def tf(s):
    return s.astype(str).str.lower() == "true"
print("resolved rows", len(res))
print("promoter agreement", (tf(res.promoter_proximal) == res.promoter_new.astype(bool)).mean())
print("accessible agreement", (tf(res.lineage_accessible) == res.accessible_new.astype(bool)).mean())
rec_abc = tf(members.abc_enhancer_overlap)
print("recorded ABC True rows", rec_abc.sum(), "of which re-annotated True", members.loc[rec_abc, "abc_new"].sum())
print("ABC agreement on resolved rows", (tf(res.abc_enhancer_overlap) == res.abc_new.astype(bool)).mean())

# 2. Repaired members table: every row resolved by position
rep = members.copy()
rep["promoter_context_status"] = "resolved"
rep["lineage_accessibility_status"] = "resolved"
rep["promoter_proximal"] = rep.promoter_new.map({True: "true", False: "false"})
rep["lineage_accessible"] = rep.accessible_new.map({True: "true", False: "false"})
rep["abc_enhancer_overlap"] = rep.abc_new.map({True: "true", False: "false"})

def context_family(mem: pd.DataFrame, tag: str) -> pd.DataFrame:
    per_set = fm.assign_clusters(fm.parse_locus(fm.build_context_masses(mem), mem))
    rows = []
    for c in fm.CTX_COLS:
        r = fm.analyse(fm.cluster_table(per_set, c, "noncoding_pip_mass"), fm.fit_rng(c))
        rows.append({"version": tag, "test": c,
                     "n_direct": r["n_clusters_direct"], "n_enzyme": r["n_clusters_enzyme"],
                     "fraction_direct": r["fraction_direct"], "fraction_enzyme": r["fraction_enzyme"],
                     "difference": r["difference"], "diff_lo": r["diff_boot_lo"], "diff_hi": r["diff_boot_hi"],
                     "ratio": r["ratio"], "ratio_lo": r["ratio_boot_lo"], "ratio_hi": r["ratio_boot_hi"]})
    # instance-weighted (credible-set-weighted) display fractions
    inst = per_set.groupby("trait_scope")[fm.CTX_COLS + ["noncoding_pip_mass"]].sum()
    for c in fm.CTX_COLS:
        for sc in fm.SCOPES:
            rows.append({"version": tag, "test": c + "__instance_weighted__" + sc,
                         "fraction_direct" if sc == fm.SCOPES[0] else "fraction_enzyme":
                         inst.loc[sc, c] / inst.loc[sc, "noncoding_pip_mass"]})
    return pd.DataFrame(rows)

orig = context_family(members, "adopted_reproduced")
repd = context_family(rep, "position_repaired")
out = pd.concat([orig, repd])
out.to_csv(OUT / "context_family_adopted_vs_repaired.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
print(out[~out.test.str.contains("__")].round(6).to_string(index=False))
iw = out[out.test.str.contains("__")].copy()
iw["value"] = iw.fraction_direct.fillna(iw.fraction_enzyme)
print(iw[["version", "test", "value"]].round(4).to_string(index=False))

# 3. Instance-weighted unresolved and promoter fractions by ancestry (adopted table)
nc = fm.build_context_masses(members.assign(_a=members.ancestry))
anc = members.drop_duplicates("credible_set_uid").set_index("credible_set_uid").ancestry
nc["ancestry"] = nc.credible_set_uid.map(anc)
t = nc.groupby(["trait_scope", "ancestry"])[fm.CTX_COLS + ["noncoding_pip_mass"]].sum()
for c in fm.CTX_COLS:
    t[c] = t[c] / t["noncoding_pip_mass"]
t["n_instances"] = nc.groupby(["trait_scope", "ancestry"]).size()
print(t.round(4).to_string())
