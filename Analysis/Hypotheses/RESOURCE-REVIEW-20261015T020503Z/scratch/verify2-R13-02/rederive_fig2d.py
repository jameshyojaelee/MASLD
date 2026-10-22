"""Independent re-derivation of Fig 2D quantities for claim R13-02.
Sentinel = position of the highest-PIP member in each credible set (independent of the
producer's decimal-restoration route). Single-linkage 1 Mb clustering across the gated
portfolio, as the caption describes."""
import numpy as np, pandas as pd
base = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
arch = pd.read_csv(base + "credible_set_pip_architecture.tsv", sep="\t", dtype={"locus": str})
mem = pd.read_csv(base + "credible_set_members_annotated.tsv", sep="\t", low_memory=False, dtype={"locus": str})
arch = arch[arch.pip_sum_gate_passed.astype(str).str.lower() == "true"].copy()
mem = mem[mem.pip_sum_gate_passed.astype(str).str.lower() == "true"].copy()
C = ["protein_altering_pip_mass", "canonical_splice_pip_mass", "synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]
print("rows after gate:", arch.trait_scope.value_counts().to_dict())
# Instance weighted
for s, a in arch.groupby("trait_scope"):
    tot = a[C].to_numpy().sum()
    f = {c: a[c].sum() / tot for c in C}
    print("INSTANCE", s, {k: round(v, 6) for k, v in f.items()}, "coding+splice=", round(f[C[0]] + f[C[1]], 6))
    sh = a.study.value_counts(normalize=True)
    print("  studies", a.study.nunique(), "traits", a.trait.nunique(), "largest study", sh.index[0], round(sh.iloc[0], 6))
# Sentinel: max-PIP member
top = mem.sort_values("normalized_susie_pip", ascending=False).drop_duplicates("credible_set_uid")
top = top.set_index("credible_set_uid")
arch["chrom"] = arch.credible_set_uid.map(top.chromosome).astype(str).str.replace("chr", "")
arch["pos"] = arch.credible_set_uid.map(top.position).astype("int64")
arch = arch.sort_values(["chrom", "pos"])
ids = []
for ch, b in arch.groupby("chrom", sort=False):
    new = np.concatenate([[True], np.diff(b.pos.to_numpy()) > 1_000_000])
    ids.append(pd.Series([f"{ch}:{i}" for i in np.cumsum(new)], index=b.index))
arch["region"] = pd.concat(ids)
print("total regions", arch.region.nunique())
res = {}
for s, a in arch.groupby("trait_scope"):
    pl = a.groupby("region")[C].mean()
    tot = pl.to_numpy().sum()
    f = {c: pl[c].sum() / tot for c in C}
    res[s] = f
    print("REGION", s, "n_regions", a.region.nunique(), {k: round(v, 6) for k, v in f.items()}, "coding+splice=", round(f[C[0]] + f[C[1]], 6))
d = res["tier1_direct_masld_pdff"]; e = res["tier2_liver_enzyme"]
r_reg = (d[C[0]] + d[C[1]]) / (e[C[0]] + e[C[1]])
print("REGION ratio", round(r_reg, 6))
shared = set(arch[arch.trait_scope == "tier1_direct_masld_pdff"].region) & set(arch[arch.trait_scope == "tier2_liver_enzyme"].region)
print("shared regions", len(shared))
# Noncoding context, producer definitions
m = mem.copy()
cons = m.consequence.astype(str)
utr = cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)", regex=True)
syn = cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)", regex=True)
nc = m[(m.consequence_category == "other_noncoding") | ((m.consequence_category == "synonymous_or_utr") & utr & ~syn)].copy()
fl = lambda c: nc[c].astype(str).str.lower() == "true"
p = nc.normalized_susie_pip.to_numpy()
pr = nc.promoter_context_status == "resolved"; lr = nc.lineage_accessibility_status == "resolved"
nc["nc"] = p
nc["promoter"] = p * (pr & fl("promoter_proximal"))
nc["abc"] = p * fl("abc_enhancer_overlap")
nc["access"] = p * (lr & fl("lineage_accessible"))
nc["resolved_no_ctx"] = p * (pr & lr & ~fl("promoter_proximal") & ~fl("abc_enhancer_overlap") & ~fl("lineage_accessible"))
nc["unresolved"] = p * (~pr | ~lr)
K = ["nc", "promoter", "abc", "access", "resolved_no_ctx", "unresolved"]
ps = nc.groupby(["credible_set_uid", "trait_scope"], as_index=False)[K].sum()
ps["region"] = ps.credible_set_uid.map(arch.set_index("credible_set_uid").region)
for s, a in ps.groupby("trait_scope"):
    iw = {k: a[k].sum() / a.nc.sum() for k in K[1:]}
    pl = a.groupby("region")[K].mean()
    rw = {k: pl[k].sum() / pl.nc.sum() for k in K[1:]}
    print("CTX", s, "n_sets", len(a), "\n  instance", {k: round(v, 6) for k, v in iw.items()}, "\n  region  ", {k: round(v, 6) for k, v in rw.items()},
          "\n  region unresolved+resolved_no_ctx =", round(rw["unresolved"] + rw["resolved_no_ctx"], 6))
