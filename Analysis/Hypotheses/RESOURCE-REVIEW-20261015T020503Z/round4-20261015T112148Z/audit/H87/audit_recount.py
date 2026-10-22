#!/usr/bin/env python3
"""H87 audit: recount the numbers stated in PRESPEC_H87.md from inputs (no gene mapping, no outcome)."""
import csv
import json
from collections import Counter

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
MEM = f"{ROOT}/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv"
VEP = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/candidates/noncoding-dna-vep-2026-08-12/run_19782908/credible_set_variant_consequences_complete.tsv"
CAT = f"{ROOT}/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
CAT_B = f"{ROOT}/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/passport_gene_index.parquet"
GATE = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/R15-01"

d = pd.read_csv(MEM, sep="\t", dtype=str, keep_default_na=False)
d["pip"] = pd.to_numeric(d.normalized_susie_pip)
d["site"] = d.chromosome + ":" + d.position
print("member rows by gate", d.pip_sum_gate_passed.value_counts().to_dict())
print("gate=false rows with empty normalized_susie_pip", (d[d.pip_sum_gate_passed == "false"].normalized_susie_pip == "").sum())
c = d[d.consequence_category.isin(["protein_altering", "canonical_splice"])].copy()
g = c[c.pip_sum_gate_passed == "true"]
u = c[c.pip_sum_gate_passed == "false"]
for lab, s in [("gated", g), ("any", c)]:
    print(lab, "rows", len(s), "sites", s.site.nunique(), "variant_ids", s.variant_id.nunique(), "studies", s.study.nunique(),
          "sites_pip>=0.1", s[s.pip >= 0.1].site.nunique(), "studies_pip>=0.1", s[s.pip >= 0.1].study.nunique(),
          "sites_pip>=0.5", s[s.pip >= 0.5].site.nunique(), "studies_pip>=0.5", s[s.pip >= 0.5].study.nunique())
h = g[g.pip >= 0.5]
print(h.groupby("trait_scope").agg(sites=("site", "nunique"), variant_ids=("variant_id", "nunique"), studies=("study", "nunique")).to_string())
print("sites with >1 variant_id string (gated)", g.groupby("site").variant_id.nunique().loc[lambda x: x > 1].to_dict())
sp = pd.to_numeric(u.susie_pip)
rs = pd.to_numeric(u.raw_cs_pip_sum)
print("ungated coding rows", len(u), "sites", u.site.nunique(), "sites only ungated", len(set(u.site) - set(g.site)),
      "max susie_pip", sp.max(), "max susie_pip/raw_cs_pip_sum", (sp / rs).max(), "raw_cs_pip_sum range", rs.min(), rs.max())
print("canonical_splice rows", c[c.consequence_category == "canonical_splice"][["study", "variant_id", "pip_sum_gate_passed", "normalized_susie_pip"]].values.tolist())
indels = [v for v in g.variant_id.unique() if len(v.split(":")[2]) != 1 or len(v.split(":")[3]) != 1]
print("gated coding indels", indels)
for s_, lab in [("1:220970028", "MTARC1"), ("19:19379549", "TM6SF2"), ("22:44324727", "PNPLA3")]:
    x = g[g.site == s_].sort_values("pip", ascending=False).iloc[0]
    print("PC-1", lab, s_, x.study, x.normalized_susie_pip, "gated studies", g[g.site == s_].study.nunique(), "any studies", c[c.site == s_].study.nunique())
print("PNPLA3 2022_36402844_PDFF_EUR", g[(g.site == "22:44324727") & (g.study == "2022_36402844_PDFF_EUR")].normalized_susie_pip.tolist())
nc3 = d[(d.pip_sum_gate_passed == "true") & (d.consequence_category == "synonymous_or_utr") & (d.pip >= 0.5)]
print("NC-3 gated synonymous_or_utr PIP>=0.5 rows", len(nc3), "variant_ids", nc3.variant_id.nunique(), "sites", nc3.site.nunique())
v = pd.read_csv(VEP, sep="\t", dtype=str, keep_default_na=False)
print("gated coding variant_ids in repo VEP table", len(set(g.variant_id) & set(v.variant_id)), "of", g.variant_id.nunique())
print("repo VEP row 4:88231392:T:TA", v[v.variant_id == "4:88231392:T:TA"].iloc[0, :14].tolist())

a = pd.read_parquet(f"{CAT}/passport_gene_index.parquet")
b = pd.read_parquet(CAT_B)
ka, kb = a.ensembl_id.str.split(".").str[0], b.ensembl_id.str.split(".").str[0]
m = a.assign(k=ka).merge(b.assign(k=kb)[["k", "primary_evidence_class"]], on="k", suffixes=("", "_b"))
print("catalog rows", len(a), len(b), "same id set", set(ka) == set(kb), "class disagreements", int((m.primary_evidence_class != m.primary_evidence_class_b).sum()))
print("primary_evidence_class", a.primary_evidence_class.value_counts().to_dict())
for sym in ["PNPLA3", "TM6SF2", "MTARC1", "GCKR", "APOE", "SERPINA1", "HSD17B13", "MBOAT7", "TRIB1", "EPPK1", "GLP1R", "HKDC1", "THRB"]:
    print("catalog", sym, a[a.symbol == sym][["ensembl_id", "primary_evidence_class"]].values.tolist())
e = pd.read_parquet(f"{CAT}/passport_evidence_long.parquet")
ge = e[e.evidence_domain == "genetics"].copy()
ge["est"] = pd.to_numeric(ge.estimate)
print("genetics assay", ge.assay.value_counts().to_dict(), "call_state", ge.call_state.value_counts().to_dict())
print("PP.H4>0.5 by call_state", ge[ge.est > 0.5].call_state.value_counts().to_dict(), "min supported PP.H4", ge[ge.call_state == "supported"].est.min())
print("TM6SF2 genetics", ge[ge.symbol == "TM6SF2"][["call_state", "estimate", "phenotype"]].values.tolist())

rows = list(csv.DictReader(open(f"{GATE}/out/coding_gene_summary_full.tsv"), delimiter="\t"))
def f(x):
    try:
        y = float(x)
        return None if y != y else y
    except ValueError:
        return None
for t in [0, 0.1, 0.5]:
    s = [r for r in rows if f(r["max_pip"]) is not None and f(r["max_pip"]) >= t]
    inc = [r for r in s if r["primary_evidence_class"]]
    print("gate genes max_pip>=", t, len(s), "in catalog", len(inc), dict(Counter(r["primary_evidence_class"] for r in inc)))
print("gate genes with no PIP (ungated-only)", [r["symbol"] for r in rows if f(r["max_pip"]) is None])
print("gate HSD17B13 present", any(r["symbol"] == "HSD17B13" for r in rows))
j = json.load(open(f"{GATE}/tmp/vep/8_144942477.json"))
print("gate VEP 8:144942477", [(t["gene_id"], t.get("gene_symbol"), t["biotype"], t["consequence_terms"]) for r in j for t in r.get("transcript_consequences", [])])
