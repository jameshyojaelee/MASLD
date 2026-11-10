# DIG1 job 74: Gene Catalog r8 sidecar (user-approved 2026-10-25 after section 45). For every genetically anchored or
# concordant gene in the reviewed r8 catalog, does the anchoring colocalized SNP also carry liver-disease liability?
# Read before this script: job 73 (section 45). 414 of 426 anchored genes are anchored by ALT, AST or GGT.
# Rule (fixed here; independent of the job 73 benchmark ratio so that GGT anchors are handled the same way):
# 1. Anchor row = the catalog genetics row (dataset_id = anchoring GWAS). Lead SNP = top_snp of that GWAS x gene row in
#    COLOC r2 (hg19); ancestry from the same row.
# 2. z = beta/se for allele1, aligned to the anchor file's alleles, oriented to the allele that raises the anchor trait
#    (for NAFLD/NASH/PDFF anchors: the risk-raising allele).
# 3. Disease GWAS of the same ancestry: EUR MVP_NAFLD_EUR, MVP_Cirrhosis_EUR, FinnGen_NAFLD, Pazoki_PDFF;
#    AFR MVP_NAFLD_AFR, MVP_Cirrhosis_AFR; AMR MVP_NAFLD_AMR, MVP_Cirrhosis_AMR; other ancestries none.
#    Allele convention check: in each file PNPLA3 rs738409 G or TM6SF2 rs58542926 T must have z > 3, else stop.
# 4. Class: disease_liability_detected if any disease z >= 3.29; opposite_direction if any <= -3.29 (and none >= 3.29);
#    no_disease_association_detected if all available |z| < 1.96 and anchor |z| >= 5.45; weak_or_indeterminate otherwise;
#    not_classifiable if no same-ancestry disease GWAS, missing anchor file or SNP. The catalog has no negative-call
#    rule, so "no association detected" is not a powered null. The threshold 3.29 (two-sided p 1e-3) is per GWAS;
#    the union over up to 4 GWAS is not corrected.
# Usage: python -I 74_catalog_enzyme_locus_sidecar.py <new_out_dir>
import sys, os
import numpy as np, pandas as pd

out = sys.argv[1]
assert not os.path.exists(out)
os.makedirs(out)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CAT = f"{ROOT}/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
SS = f"{ROOT}/GWAS/finemapping/data/sumstats"
COLOC = f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
J73 = f"{ROOT}/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci"
DIS = {"EUR": ["MVP_NAFLD_EUR", "MVP_Cirrhosis_EUR", "FinnGen_NAFLD", "Pazoki_PDFF"],
       "AFR": ["MVP_NAFLD_AFR", "MVP_Cirrhosis_AFR"], "AMR": ["MVP_NAFLD_AMR", "MVP_Cirrhosis_AMR"]}
BENCH = [(22, 44324727, "G"), (19, 19379549, "T")]

g = pd.read_parquet(f"{CAT}/passport_gene_index.parquet")
e = pd.read_parquet(f"{CAT}/passport_evidence_long.parquet")
A = g[g.primary_evidence_class.isin(["genetically_anchored", "concordant"])][["ensembl_id", "symbol", "primary_evidence_class"]]
ge = e[(e.evidence_domain == "genetics")][["ensembl_id", "dataset_id", "phenotype", "estimate"]]
A = A.merge(ge, on="ensembl_id", how="left").rename(columns={"dataset_id": "anchor_gwas", "phenotype": "anchor_trait",
                                                             "estimate": "anchor_pp4"})
C = pd.read_csv(COLOC, usecols=["gwas_name", "ensembl", "PP.H4.abf", "PP.H4.susie", "top_snp", "ancestry"])
C["h4"] = C[["PP.H4.abf", "PP.H4.susie"]].max(axis=1, skipna=True)
C = C.sort_values("h4", ascending=False).drop_duplicates(["gwas_name", "ensembl"])
A = A.merge(C, left_on=["anchor_gwas", "ensembl_id"], right_on=["gwas_name", "ensembl"], how="left")
sp = A.top_snp.astype(str).str.replace("chr", "").str.split(":", expand=True)
A["chr"] = pd.to_numeric(sp[0], errors="coerce"); A["pos"] = pd.to_numeric(sp[1].str.split("_").str[0], errors="coerce")
print("anchored genes:", len(A), " with COLOC lead SNP:", A.chr.notna().sum())
print("anchor PP4 agrees with COLOC row (|diff|<1e-6):", int((abs(A.anchor_pp4 - A.h4) < 1e-6).sum()))

want = set(zip(A.chr.dropna().astype(int), A.pos.dropna().astype(int))) | {(c, p) for c, p, _ in BENCH}
def extract(t):
    f = f"{SS}/{t}_reformatted_hg19.tsv"
    if not os.path.exists(f): return None
    keep = []
    for ch in pd.read_csv(f, sep="\t", chunksize=2_000_000, usecols=["chromosome", "position", "allele1", "allele2", "beta", "se"],
                          dtype={"chromosome": str}):
        ch["chromosome"] = pd.to_numeric(ch.chromosome.str.replace("chr", ""), errors="coerce")
        ch = ch.dropna(subset=["chromosome"])
        keep.append(ch[[(int(a), int(b)) in want for a, b in zip(ch.chromosome, ch.position)]])
    d = pd.concat(keep).rename(columns={"chromosome": "chr", "position": "pos"})
    d["chr"] = d.chr.astype(int); d["z"] = d.beta / d.se
    d["allele1"] = d.allele1.str.upper(); d["allele2"] = d.allele2.str.upper()
    return d.drop_duplicates(["chr", "pos"], keep=False)

files = set(A.anchor_gwas.dropna()) | {t for v in DIS.values() for t in v}
Z = {t: extract(t) for t in sorted(files)}
chk = []
for t, d in Z.items():
    if d is None: chk.append((t, "missing_file", np.nan)); continue
    best = -np.inf
    for c, p, a in BENCH:
        r = d[(d.chr == c) & (d.pos == p)]
        if len(r): best = max(best, r.z.iloc[0] if r.allele1.iloc[0] == a else -r.z.iloc[0] if r.allele2.iloc[0] == a else -np.inf)
    chk.append((t, "ok" if best > 3 else "FAIL", best))
chk = pd.DataFrame(chk, columns=["gwas", "allele_check", "max_benchmark_z"]); chk.to_csv(f"{out}/allele_convention_check.tsv", sep="\t", index=False)
print(chk.to_string(index=False))
bad = chk[(chk.allele_check == "FAIL") & chk.gwas.isin([t for v in DIS.values() for t in v])]
assert bad.empty, f"disease file allele check failed: {bad.gwas.tolist()}"
fail_anchor = set(chk.gwas[chk.allele_check != "ok"])

rows = []
for r in A.itertuples():
    o = dict(ensembl_id=r.ensembl_id, symbol=r.symbol, primary_evidence_class=r.primary_evidence_class, anchor_gwas=r.anchor_gwas,
             anchor_trait=r.anchor_trait, anchor_pp4=r.anchor_pp4, ancestry=r.ancestry, lead_snp_hg19=r.top_snp)
    d = Z.get(r.anchor_gwas)
    reason = None
    if pd.isna(r.chr): reason = "no_coloc_lead_snp"
    elif d is None: reason = "anchor_sumstats_missing"
    elif r.anchor_gwas in fail_anchor: reason = "anchor_allele_convention_unverified"
    elif r.ancestry not in DIS: reason = "no_same_ancestry_disease_gwas"
    if reason is None:
        s = d[(d.chr == int(r.chr)) & (d.pos == int(r.pos))]
        if s.empty: reason = "lead_snp_absent_in_anchor_file"
    if reason is not None:
        o.update(sidecar_class="not_classifiable", reason=reason); rows.append(o); continue
    a1, a2, za = s.allele1.iloc[0], s.allele2.iloc[0], s.z.iloc[0]
    sg = 1.0 if za > 0 else -1.0
    o.update(raising_allele=a1 if sg > 0 else a2, z_anchor=abs(za))
    dz = {}
    for t in DIS[r.ancestry]:
        q = Z[t][(Z[t].chr == int(r.chr)) & (Z[t].pos == int(r.pos))]
        if q.empty: continue
        if q.allele1.iloc[0] == a1 and q.allele2.iloc[0] == a2: dz[t] = q.z.iloc[0] * sg
        elif q.allele1.iloc[0] == a2 and q.allele2.iloc[0] == a1: dz[t] = -q.z.iloc[0] * sg
    for t in DIS[r.ancestry]: o[f"z_{t}"] = dz.get(t, np.nan)
    v = np.array(list(dz.values()))
    if len(v) == 0: cls, why = "not_classifiable", "lead_snp_absent_in_disease_files"
    elif (v >= 3.29).any(): cls, why = "disease_liability_detected", ""
    elif (v <= -3.29).any(): cls, why = "opposite_direction", ""
    elif (np.abs(v) < 1.96).all() and abs(za) >= 5.45: cls, why = "no_disease_association_detected", ""
    else: cls, why = "weak_or_indeterminate", ""
    o.update(sidecar_class=cls, reason=why, n_disease_gwas=len(v)); rows.append(o)
S = pd.DataFrame(rows)
G73 = pd.read_csv(f"{J73}/coloc_genes_lead_z.tsv", sep="\t").merge(pd.read_csv(f"{J73}/loci.tsv", sep="\t")[["locus", "cls"]], on="locus")
S = S.merge(G73[["ensembl", "cls"]].rename(columns={"ensembl": "ensembl_id", "cls": "job73_mvp_class"}), on="ensembl_id", how="left")
WORD = {"disease_liability_detected": "The anchoring colocalized SNP is also associated with liver disease or liver fat in the same ancestry.",
        "no_disease_association_detected": "The anchoring SNP is associated with the liver-enzyme trait but not detectably with liver disease or liver fat; not a powered null.",
        "weak_or_indeterminate": "Disease association at the anchoring SNP is weak or below threshold; indeterminate.",
        "opposite_direction": "The allele raising the anchor trait lowers liver-disease risk.",
        "not_classifiable": "No same-ancestry disease GWAS or no usable anchor SNP; not assessed."}
S["allowed_wording"] = S.sidecar_class.map(WORD)
S["catalog_release_id"] = "gene-catalog-corrected-r8-20261009T142027Z"
S.to_csv(f"{out}/catalog_r8_anchor_disease_liability.tsv", sep="\t", index=False)
T = pd.crosstab(S.anchor_trait, S.sidecar_class, margins=True); T.to_csv(f"{out}/summary_by_trait.tsv", sep="\t")
print(T.to_string())
print(S[S.sidecar_class == "not_classifiable"].reason.value_counts().to_string())
print(pd.crosstab(S.job73_mvp_class.fillna("-"), S.sidecar_class).to_string())
pd.set_option("display.width", 250)
print(S[S.sidecar_class == "disease_liability_detected"].sort_values("z_anchor", ascending=False)
      [["symbol", "anchor_gwas", "z_anchor"] + [c for c in S.columns if c.startswith("z_") and c != "z_anchor"]].round(1).to_string(index=False))
os.system(f"/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/pip freeze > {out}/pip_freeze.txt")
