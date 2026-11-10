"""RV6: allele directions at three-way loci (GWAS x liver caQTL x liver eQTL). Follows PRESPEC_RV6.md.
Unit: GWAS locus (H53 locus id). Usage: python -I rv6_allele_direction.py <new_out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H53 = ROOT / "Analysis/Hypotheses/H53-20261015T112655Z"
CQ = ROOT / "data/external/currin_2025_caqtl"
EQ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl")
FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
DISC = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
OUT = Path(sys.argv[1])
PAL = {frozenset("AT"), frozenset("CG")}

tw = pd.read_csv(H53 / "results/three_way.tsv", sep="\t")
tw = tw[(tw["tw_PP.H4"] > 0.8) & (~tw.tw_dropped)]
lp = pd.read_csv(H53 / "results/lead_peak_pairs.tsv.gz", sep="\t",
                 usecols=["lead_key", "peak", "caqtl_lead_variant", "caqtl_lead_snp19", "locus", "gwas", "PP.H4", "lead_ld_r2", "arm"])
T = tw.merge(lp, on=["lead_key", "peak"], how="left", validate="many_to_one")
assert len(T) == 972 and T.caqtl_lead_variant.notna().all()
v = T.caqtl_lead_variant.str.split(":", expand=True)
T["chr"] = v[0].str.replace("chr", "").astype(int)
T["pos38"] = v[1].astype(int)
T["ref"], T["eff"] = v[2].str.upper(), v[3].str.upper()
T["pos19"] = T.caqtl_lead_snp19.astype(str).str.split(":").str[1].astype(int)
T["palindromic"] = [frozenset((a, b)) in PAL for a, b in zip(T.ref, T.eff)]
T["ensembl"] = T.ensembl.astype(str).str.replace(r"\..*$", "", regex=True)
tier = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t").set_index("study_name")
reg = pd.read_csv(FM / "config/gwas_registry.tsv", sep="\t").set_index("study_name")
T["trait"] = T.gwas.map(tier.trait)
T["tier"] = T.gwas.map(tier.tier)

# ---- caQTL lookups and the S23A allele check ----
s23 = pd.read_csv(H53 / "downloads_untrusted/currin_S23A.tsv", sep="\t", dtype=str)
s23 = s23[s23.caQTL_variant.str.contains(":", na=False)].copy()
s23["chr"] = s23.caQTL_variant.str.split(":").str[0]
s23["pos"] = s23.caQTL_variant.str.split(":").str[1]
s23["key"] = "chr" + s23.chr + ":" + s23.pos + "|" + s23.caPeak
s23 = s23.drop_duplicates("key")
need = set(T.caqtl_lead_variant + "|" + T.peak)
ca, chk = [], []
for c in sorted(set(T.chr) | set(pd.to_numeric(s23.chr, errors="coerce").dropna().astype(int))):
    f = list(CQ.glob(f"caQTL_hg38_chr{c}_*.txt.gz"))
    if len(f) != 1:
        continue
    peaks = set(T.peak[T.chr == c]) | set(s23.caPeak[s23.chr == str(c)])
    for ch in pd.read_csv(f[0], sep="\t", usecols=["variant", "peak", "beta", "varbeta"], chunksize=5_000_000):
        ch = ch[ch.peak.isin(peaks)]
        k = ch.variant + "|" + ch.peak
        ca.append(ch[k.isin(need)])
        k2 = ch.variant.str.rsplit(":", n=2).str[0] + "|" + ch.peak
        chk.append(ch[k2.isin(set(s23.key))].assign(key=k2[k2.isin(set(s23.key))]))
    print("caqtl chr", c, flush=True)
ca = pd.concat(ca).drop_duplicates(["variant", "peak"])
chk = pd.concat(chk).merge(s23[["key", "caQTL_EA", "caQTL_NEA", "caQTL_beta"]], on="key")
chk["caQTL_beta"] = pd.to_numeric(chk.caQTL_beta, errors="coerce")
chk["last"] = chk.variant.str.split(":").str[3]
chk["third"] = chk.variant.str.split(":").str[2]
chk["same_value"] = (chk.beta - chk.caQTL_beta).abs() < 1e-4
chk["neg_value"] = (chk.beta + chk.caQTL_beta).abs() < 1e-4
chk = chk[chk.same_value | chk.neg_value]
chk["agrees_last_allele"] = (chk.same_value & (chk.caQTL_EA == chk["last"])) | (chk.neg_value & (chk.caQTL_EA == chk.third))
chk.to_csv(OUT / "allele_check_S23A.tsv", sep="\t", index=False)
print("S23A allele check:", int(chk.agrees_last_allele.sum()), "of", len(chk), flush=True)
assert len(chk) >= 10 and chk.agrees_last_allele.mean() >= 0.95, "caQTL allele convention check failed"
T = T.merge(ca.rename(columns={"variant": "caqtl_lead_variant"}), on=["caqtl_lead_variant", "peak"], how="left")
T["z_ca"] = T.beta / np.sqrt(T.varbeta)

# ---- eQTL lookups (hg19) ----
eqr = []
for c, nd in T.groupby("chr"):
    keys = set(zip(nd.ensembl, nd.pos19))
    for ch in pd.read_csv(EQ / f"chr{c}_marginal_summary_results.tsv", sep="\t", chunksize=3_000_000,
                          usecols=["POS", "NEA", "EA", "Beta", "SE", "ENSG"]):
        ch = ch[ch.POS.isin(nd.pos19)]
        eqr.append(ch.loc[np.array([(e, p) in keys for e, p in zip(ch.ENSG, ch.POS)], dtype=bool)].assign(chr=c))
    print("eqtl chr", c, flush=True)
eqr = pd.concat(eqr).rename(columns={"ENSG": "ensembl", "POS": "pos19"})
eqr["z_eqtl_EA"] = eqr.Beta / eqr.SE
eqr = eqr.drop_duplicates(["chr", "ensembl", "pos19", "EA", "NEA"])
T = T.merge(eqr[["chr", "ensembl", "pos19", "EA", "NEA", "z_eqtl_EA"]], on=["chr", "ensembl", "pos19"], how="left")
ea, nea = T.EA.astype(str).str.upper(), T.NEA.astype(str).str.upper()
T["z_eqtl"] = np.where((ea == T.eff) & (nea == T.ref), T.z_eqtl_EA, np.where((ea == T.ref) & (nea == T.eff), -T.z_eqtl_EA, np.nan))

# ---- GWAS lookups (hg19) ----
gwr = []
for s, ps in T.groupby("gwas"):
    pos = set(zip(ps.chr, ps.pos19)); posset = {p for _, p in pos}
    for ch in pd.read_csv(FM / reg.loc[s, "sumstats_path"], sep="\t", chunksize=4_000_000, dtype={"chromosome": str}):
        ch = ch[ch.position.isin(posset)]
        ch = ch.loc[np.array([(int(c) if str(c).isdigit() else -1, p) in pos for c, p in zip(ch.chromosome, ch.position)], dtype=bool)]
        gwr.append(ch.assign(gwas=s, chr=pd.to_numeric(ch.chromosome)))
    print("gwas", s, flush=True)
gwr = pd.concat(gwr).rename(columns={"position": "pos19"})
gwr["z_a1"] = gwr.beta / gwr.se
gwr = gwr.drop_duplicates(["gwas", "chr", "pos19", "allele1", "allele2"])
T = T.merge(gwr[["gwas", "chr", "pos19", "allele1", "allele2", "z_a1"]], on=["gwas", "chr", "pos19"], how="left")
a1, a2 = T.allele1.astype(str).str.upper(), T.allele2.astype(str).str.upper()
T["z_gwas"] = np.where((a1 == T.eff) & (a2 == T.ref), T.z_a1, np.where((a1 == T.ref) & (a2 == T.eff), -T.z_a1, np.nan))
T.loc[T.palindromic, ["z_eqtl", "z_gwas"]] = np.nan

# ---- disease direction ----
d = pd.read_csv(DISC, usecols=["gene", "logFC", "padj"])
d["ensembl"] = d.gene.str.replace(r"\..*$", "", regex=True)
T = T.merge(d[["ensembl", "logFC", "padj"]], on="ensembl", how="left")

def label(df, zc, zg):
    s_ca = np.where(df.z_ca.abs() >= zc, np.sign(df.z_ca), np.nan)
    s_e = np.where(df.z_eqtl.abs() >= zc, np.sign(df.z_eqtl), np.nan)
    s_g = np.where(df.z_gwas.abs() >= zg, np.sign(df.z_gwas), np.nan)
    out = pd.DataFrame(index=df.index)
    out["ca_e_same"] = s_ca * s_e
    out["trait_allele_on_access"] = s_g * s_ca
    out["trait_allele_on_expr"] = s_g * s_e
    out["genetic_disease_concordant"] = np.where(np.isfinite(s_g * s_e) & (df.padj < 0.05), (s_g * s_e) == np.sign(df.logFC), np.nan)
    return out

T = T.join(label(T, 3, 2))
T["trait_raising_allele"] = np.where(T.z_gwas > 0, T.eff, np.where(T.z_gwas < 0, T.ref, ""))
T.to_csv(OUT / "three_way_alleles.tsv", sep="\t", index=False)

def tests(T, zc, zg, tag):
    L = T.sort_values(["tw_PP.H4", "PP.H4"], ascending=False).drop_duplicates("locus").copy()
    lab = label(L, zc, zg)
    t1 = lab.ca_e_same.dropna()
    k1 = int((t1 > 0).sum()); n1 = len(t1)
    p1 = stats.binomtest(k1, n1, 0.5).pvalue if n1 else np.nan
    liver = L.tier <= 2
    t2 = lab.genetic_disease_concordant[liver].dropna()
    k2 = int((t2 == 1).sum()); n2 = len(t2)
    p2 = stats.binomtest(k2, n2, 0.5).pvalue if n2 else np.nan
    gs = (lab.trait_allele_on_expr[liver] > 0)[t2.index]; ds = (L.logFC[liver] > 0)[t2.index]
    fisher = stats.fisher_exact(pd.crosstab(gs, ds).reindex(index=[False, True], columns=[False, True], fill_value=0).values) if n2 else (np.nan, np.nan)
    r = pd.DataFrame([dict(analysis=tag, test="T1_access_expression_same", n_loci=n1, k=k1, fraction=k1 / n1 if n1 else np.nan, p=p1),
                      dict(analysis=tag, test="T2_genetic_disease_concordant", n_loci=n2, k=k2, fraction=k2 / n2 if n2 else np.nan, p=p2,
                           fisher_or=fisher[0], fisher_p=fisher[1])])
    return L, r  # L carries the primary (z3/z2) labels from T

L, prim = tests(T, 3, 2, "primary_z3_z2")
# Holm over 2
o = prim.p.sort_values()
h = [min(1, 2 * o.iloc[0]), min(1, max(2 * o.iloc[0], o.iloc[1]))]
prim["holm_p"] = pd.Series(h, index=o.index).reindex(prim.index)
_, sens = tests(T, 2, 2, "sensitivity_z2_z2")
pd.concat([prim, sens]).to_csv(OUT / "tests.tsv", sep="\t", index=False)
L.to_csv(OUT / "locus_table.tsv", sep="\t", index=False)

with open(OUT / "summary.txt", "w") as fh:
    fh.write(f"Three-way rows {len(T)}; loci {T.locus.nunique()}; genes {T.gene.nunique()}\n")
    fh.write(f"Rows with caQTL z {T.z_ca.notna().sum()}, eQTL z {T.z_eqtl.notna().sum()}, GWAS z {T.z_gwas.notna().sum()}; palindromic {int(T.palindromic.sum())}\n")
    fh.write(f"S23A allele check: {int(chk.agrees_last_allele.sum())} of {len(chk)} agree with last-allele convention\n\n")
    fh.write(pd.concat([prim, sens]).to_string(index=False) + "\n\n")
    cols = ["locus", "gene", "gwas", "trait", "caqtl_lead_variant", "trait_raising_allele", "z_gwas", "z_ca", "z_eqtl", "logFC", "padj",
            "ca_e_same", "trait_allele_on_access", "trait_allele_on_expr", "genetic_disease_concordant"]
    fh.write("Locus table (one row per locus)\n" + L[cols].round(3).to_string(index=False) + "\n")
print("done")
