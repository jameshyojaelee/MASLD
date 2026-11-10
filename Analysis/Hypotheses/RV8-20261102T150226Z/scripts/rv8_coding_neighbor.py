"""RV8: flag colocalized genes anchored at, or in LD with, another gene's coding credible-set member. Follows PRESPEC_RV8.md.
Usage: PLINK=/path/to/plink python -I rv8_coding_neighbor.py <new_out_dir>
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
PAIRS = FM / "results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz"
MEM = ROOT / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv"
GTF = ROOT / "data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf"
OUT = Path(sys.argv[1])
PLINK = os.environ["PLINK"]
WIN = 1_000_000
CODING = ["protein_altering", "canonical_splice"]

tier = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t").set_index("study_name")
reg = pd.read_csv(FM / "config/gwas_registry.tsv", sep="\t").set_index("study_name")

# ---- colocalized pairs ----
p = pd.read_csv(PAIRS, sep="\t")
p = p[(p["PP.H4"] > 0.5) & p.eligible].copy()
p["tier"] = p.gwas.map(tier.tier)
p = p[p.tier <= 2].copy()
p["ancestry"] = p.gwas.map(reg.ancestry).str.lower()
p["chr"] = p.hit1.str.split(":").str[0].astype(int)
p["anchor_pos"] = p.hit1.str.split(":").str[1].astype(int)
p["ens"] = p.ensembl.str.replace(r"\..*$", "", regex=True)
print(f"pairs {len(p)}; genes {p.gene.nunique()}; unique gene-gwas {len(p[['gene','gwas']].drop_duplicates())}", flush=True)
assert p.gene.nunique() == 436

# ---- coding credible-set members, mapped to genes ----
m = pd.read_csv(MEM, sep="\t", low_memory=False, dtype={"locus": str})
m = m[m.pip_sum_gate_passed.astype(str).str.upper().eq("TRUE") & m.consequence_category.isin(CODING)].copy()
m["tier"] = m.study.map(tier.tier)
m = m[m.tier <= 2].copy()
m["chr"] = m.chromosome.astype(str).str.replace("chr", "").astype(int)

g = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 8], names=["chr", "feat", "start", "end", "attr"])
g = g[g.feat.isin(["CDS", "exon"]) & g.attr.str.contains('gene_type "protein_coding"')].copy()
g["chr"] = g.chr.str.replace("chr", "")
g = g[g.chr.isin([str(i) for i in range(1, 23)])].copy()
g["chr"] = g.chr.astype(int)
g["gene_name"] = g.attr.str.extract(r'gene_name "([^"]+)"')[0]
g["gene_id"] = g.attr.str.extract(r'gene_id "(ENSG\d+)')[0]
g = g.drop(columns="attr")

def genes_at(c, pos, cat):
    s = g[g.chr == c]
    if cat == "canonical_splice":
        e = s[s.feat == "exon"]
        hit = e[((e.start - pos).abs() <= 2) | ((e.end - pos).abs() <= 2)]
    else:
        hit = s[(s.feat == "CDS") & (s.start <= pos) & (s.end >= pos)]
        if hit.empty:
            hit = s[(s.feat == "exon") & (s.start <= pos) & (s.end >= pos)]
    return sorted(set(zip(hit.gene_id, hit.gene_name)))

u = m[["chr", "position", "consequence_category"]].drop_duplicates()
mp = {(r.chr, r.position): genes_at(r.chr, r.position, r.consequence_category) for r in u.itertuples()}
m["coding_genes"] = [";".join(n for _, n in mp[(c, q)]) for c, q in zip(m.chr, m.position)]
m["coding_gene_ids"] = [";".join(i for i, _ in mp[(c, q)]) for c, q in zip(m.chr, m.position)]
m.drop(columns=[c for c in m.columns if c.startswith(("promoter", "abc", "lineage"))]).to_csv(OUT / "coding_members_mapped.tsv", sep="\t", index=False)
print(f"coding members {len(m)} rows, {len(u)} positions, unmapped {sum(1 for v in mp.values() if not v)}", flush=True)
cm = (m.groupby(["chr", "position", "coding_genes", "coding_gene_ids"], as_index=False)
        .agg(max_pip=("normalized_susie_pip", "max"), consequence=("consequence", "first"),
             studies=("study", lambda x: ",".join(sorted(set(x))))))

# ---- LD, per ancestry and chromosome ----
need = []
for r in p.itertuples():
    near = cm[(cm.chr == r.chr) & ((cm.position - r.anchor_pos).abs() <= WIN)]
    for q in near.itertuples():
        need.append((r.Index, q.Index))
need = pd.DataFrame(need, columns=["pair", "cm"])
r2 = {}
for (anc, c), idx in p.loc[need.pair.unique()].groupby(["ancestry", "chr"]).groups.items():
    pos = set(p.loc[idx, "anchor_pos"]) | set(cm.loc[need[need.pair.isin(idx)].cm, "position"])
    tmp = OUT / f"tmp_{anc}_{c}"
    with open(f"{tmp}.range", "w") as fh:
        for x in sorted(pos):
            fh.write(f"{c}\t{x}\t{x}\tv{x}\n")
    bf = FM / f"data/ld_ref/1kg_{anc}/chr{c}_{anc}"
    run = subprocess.run([PLINK, "--bfile", str(bf), "--extract", "range", f"{tmp}.range", "--recode", "A", "--out", str(tmp)],
                         stdout=subprocess.DEVNULL)
    raw = pd.read_csv(f"{tmp}.raw", sep=r"\s+").iloc[:, 6:] if run.returncode == 0 else pd.DataFrame()
    bim = pd.read_csv(f"{bf}.bim", sep="\t", header=None, names=["c", "id", "cm", "pos", "a1", "a2"]).set_index("id")
    ids = [col.rsplit("_", 1)[0] for col in raw.columns]
    X = raw.fillna(raw.mean()).to_numpy(float)
    cols = {}
    for j, i in enumerate(ids):
        cols.setdefault(int(bim.loc[i, "pos"]) if np.ndim(bim.loc[i, "pos"]) == 0 else int(bim.loc[i, "pos"].iloc[0]), j)
    for i in idx:
        a = p.loc[i, "anchor_pos"]
        for k in need[need.pair == i].cm:
            b = cm.loc[k, "position"]
            if a in cols and b in cols:
                x, y = X[:, cols[a]], X[:, cols[b]]
                r2[(i, k)] = np.nan if x.std() == 0 or y.std() == 0 else np.corrcoef(x, y)[0, 1] ** 2
            else:
                r2[(i, k)] = 1.0 if a == b else np.nan
    for ext in (".range", ".raw", ".log", ".nosex"):
        Path(f"{tmp}{ext}").unlink(missing_ok=True)
    print("ld", anc, c, flush=True)

need["r2"] = [r2.get((a, b), np.nan) for a, b in zip(need.pair, need.cm)]
need = need.join(cm, on="cm").join(p[["ens", "gene", "anchor_pos"]], on="pair")
need["other_gene"] = [e not in ids.split(";") for e, ids in zip(need.ens, need.coding_gene_ids)]
need["same_variant"] = need.position == need.anchor_pos

def classify(rows, min_pip):
    r = rows[rows.max_pip >= min_pip]
    oth, own = r[r.other_gene & (r.coding_gene_ids != "")], r[~r.other_gene]
    if oth.same_variant.any():
        return "anchor_is_other_coding", oth[oth.same_variant].iloc[0]
    if (oth.r2 >= 0.8).any():
        return "not_separable", oth.sort_values("r2", ascending=False).iloc[0]
    if (oth.r2 >= 0.2).any():
        return "moderate_LD", oth.sort_values("r2", ascending=False).iloc[0]
    if own.same_variant.any() or (own.r2 >= 0.8).any():
        return "own_gene_coding", own.sort_values("r2", ascending=False).iloc[0]
    tag = "ld_missing" if (oth.r2.isna() & ~oth.same_variant).any() else "none"
    return tag, None

ORDER = ["anchor_is_other_coding", "not_separable", "moderate_LD", "own_gene_coding", "ld_missing", "none"]
rows = []
for i, r in p.iterrows():
    n = need[need.pair == i]
    for label, min_pip in (("primary", 0.1), ("any_pip", 0.0)):
        cls, best = classify(n, min_pip)
        rows.append(dict(pair=i, analysis=label, gene=r.gene, ensembl=r.ens, gwas=r.gwas, tier=r.tier, ancestry=r.ancestry,
                         anchor=r.hit1, pp4=r["PP.H4"], cls=cls,
                         coding_gene=None if best is None else best.coding_genes,
                         coding_variant=None if best is None else f"{r.chr}:{best.position}",
                         coding_pip=None if best is None else best.max_pip,
                         consequence=None if best is None else best.consequence,
                         r2=None if best is None else (1.0 if best.same_variant else best.r2),
                         coding_studies=None if best is None else best.studies))
pf = pd.DataFrame(rows)
pf.to_csv(OUT / "pair_flags.tsv", sep="\t", index=False)
pf["rank"] = pf.cls.map(ORDER.index)
gf = pf.sort_values(["rank", "pp4"], ascending=[True, False]).drop_duplicates(["analysis", "gene"]).copy()
gf["direct"] = gf.gene.isin(pf[pf.tier == 1].gene)
gf.drop(columns="rank").to_csv(OUT / "gene_flags.tsv", sep="\t", index=False)

M2 = {"SAMM50": "PNPLA3", "MAU2": "TM6SF2", "SPATA31H1": "GCKR", "SUPT7L": "GCKR", "IFI27L2": "SERPINA1",
      "SERPINA11": "SERPINA1", "CNOT1": "GOT2", "SPP1": "HSD17B13"}
with open(OUT / "summary.txt", "w") as fh:
    fh.write(f"Pairs {len(p)}; genes {p.gene.nunique()}; direct-trait (tier 1) genes {gf[(gf.analysis=='primary') & gf.direct].shape[0]}\n")
    fh.write(f"Coding members (tier 1-2, gate passed) {len(m)} rows at {len(u)} positions\n\n")
    for a in ("primary", "any_pip"):
        x = gf[gf.analysis == a]
        fh.write(f"Gene-level classes ({a})\n" + pd.crosstab(x.cls, x.direct, margins=True).reindex(ORDER + ["All"]).fillna(0).astype(int).to_string() + "\n\n")
    x = gf[(gf.analysis == "primary") & gf.cls.isin(ORDER[:4])]
    fh.write("Flagged genes (primary)\n" + x[["gene", "direct", "cls", "gwas", "anchor", "coding_gene", "coding_variant", "consequence", "coding_pip", "r2"]]
             .sort_values(["cls", "gene"]).round(3).to_string(index=False) + "\n\n")
    fh.write("M2 check (hand-found examples)\n")
    for gn, cg in M2.items():
        y = gf[(gf.analysis == "primary") & (gf.gene == gn)]
        fh.write(f"  {gn} (expected {cg}): " + ("not among the 436" if y.empty else f"{y.cls.iloc[0]} / {y.coding_gene.iloc[0]} / r2 {y.r2.iloc[0]}") + "\n")
print("done")
