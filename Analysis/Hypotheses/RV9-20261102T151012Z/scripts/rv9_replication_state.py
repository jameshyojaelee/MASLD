"""RV9: per-gene replication of the 436 SuSiE-colocalized genes in an independent same-trait European GWAS.
Follows PRESPEC_RV9.md. Usage: python -I rv9_replication_state.py <new_out_dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
R2 = FM / "results/susie_coloc_r2_20260924T053244Z"
OUT = Path(sys.argv[1])
SEED, NBOOT = 20261102, 10_000
PAIRS = [("MVP_ALT_EUR", "UKBB_ALT"), ("MVP_AST_EUR", "UKBB_AST"), ("MVP_NAFLD_EUR", "FinnGen_NAFLD")]
ORDER = ["replicated", "distinct variant", "neither", "untested"]

tier = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t").set_index("study_name")
a = pd.read_csv(R2 / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "gene", "ensembl", "PP.H3.susie", "PP.H4.susie"])
a["tier"] = a.gwas_name.map(tier.tier)
genes = a[(a.tier <= 2) & (a["PP.H4.susie"] > 0.5)].drop_duplicates("ensembl")[["ensembl", "gene"]]
assert len(genes) == 436, len(genes)
a = a.drop_duplicates(["gwas_name", "ensembl"]).set_index(["gwas_name", "ensembl"])

sp = pd.read_csv(R2 / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
sp = sp[(sp["PP.H4"] > 0.5) & sp.eligible].sort_values("PP.H4", ascending=False).drop_duplicates(["gwas", "ensembl"])
anchor = sp.set_index(["gwas", "ensembl"]).hit1

rows = []
for d1, d2 in PAIRS:
    for d, r in ((d1, d2), (d2, d1)):
        disc = a.loc[d]
        disc = disc[(disc["PP.H4.susie"] > 0.5) & disc.index.isin(genes.ensembl)]
        for e, x in disc.iterrows():
            if (r, e) in a.index and pd.notna(a.loc[(r, e), "PP.H4.susie"]):
                h3, h4 = a.loc[(r, e), "PP.H3.susie"], a.loc[(r, e), "PP.H4.susie"]
                st = "replicated" if h4 > 0.5 else "distinct variant" if h3 > 0.5 else "neither"
            else:
                h3 = h4 = np.nan
                st = "untested"
            rows.append(dict(discovery=d, replication=r, ensembl=e, gene=x.gene, disc_pp4=x["PP.H4.susie"],
                             rep_pp3=h3, rep_pp4=h4, state=st, anchor=anchor.get((d, e), np.nan)))
rows = pd.DataFrame(rows)
rows["band"] = pd.cut(rows.disc_pp4, [0.5, 0.8, 0.95, 1.0000001], right=False, labels=["0.5-0.8", "0.8-0.95", ">=0.95"])

# locus clusters: anchors on one chromosome within 500 kb, single linkage
an = rows.anchor.dropna().drop_duplicates()
ch = an.str.split(":").str[0]; po = an.str.split(":").str[1].astype(int)
cl, k = {}, 0
for c in sorted(ch.unique()):
    s = po[ch == c].sort_values()
    prev = None
    for v, p in zip(an[s.index], s):
        if prev is None or p - prev > 500_000:
            k += 1
        cl[v] = k; prev = p
rows["cluster"] = rows.anchor.map(cl)
miss = rows.cluster.isna()
rows.loc[miss, "cluster"] = [f"noanchor_{i}" for i in range(miss.sum())]
rows["dbgap_mvp"] = rows.discovery.str.startswith("MVP") | rows.replication.str.startswith("MVP")
rows.to_csv(OUT / "rows.tsv", sep="\t", index=False)

# gene states
rank = {s: i for i, s in enumerate(ORDER)}
best = rows.assign(r=rows.state.map(rank)).sort_values("r").drop_duplicates("ensembl")
gs = genes.merge(best[["ensembl", "state", "discovery", "replication", "disc_pp4", "rep_pp3", "rep_pp4"]], on="ensembl", how="left")
gs["state"] = gs.state.fillna("no comparison")
gs["direct"] = gs.ensembl.isin(a.reset_index().query("tier == 1 and `PP.H4.susie` > 0.5").ensembl)
gs.to_csv(OUT / "gene_states.tsv", sep="\t", index=False)

# bands with locus-cluster bootstrap
rng = np.random.default_rng(SEED)
t = rows[rows.state != "untested"].copy()
t["rep"] = (t.state == "replicated").astype(int)
t["cl"] = t.cluster.astype(str)
cls = t.cl.unique()
BANDS = ["0.5-0.8", "0.8-0.95", ">=0.95", "all"]
# per-cluster tested and replicated counts per band; a bootstrap draw weights clusters by how often they are picked
nmat = np.array([[len(t[(t.cl == c) & ((t.band == b) | (b == "all"))]) for b in BANDS] for c in cls])
kmat = np.array([[t[(t.cl == c) & ((t.band == b) | (b == "all"))].rep.sum() for b in BANDS] for c in cls])
w = rng.multinomial(len(cls), np.full(len(cls), 1 / len(cls)), size=NBOOT)
with np.errstate(invalid="ignore", divide="ignore"):
    boot = (w @ kmat) / (w @ nmat)
bands = []
for j, b in enumerate(BANDS):
    sub = t if b == "all" else t[t.band == b]
    lo, hi = np.nanpercentile(boot[:, j], [2.5, 97.5])
    bands.append(dict(band=b, n_tested_rows=len(sub), replicated=int(sub.rep.sum()), distinct=int((sub.state == "distinct variant").sum()),
                      neither=int((sub.state == "neither").sum()), fraction=sub.rep.mean(), ci_lo=lo, ci_hi=hi))
bands = pd.DataFrame(bands)
bands.to_csv(OUT / "bands.tsv", sep="\t", index=False)

with open(OUT / "summary.txt", "w") as fh:
    fh.write(f"Genes 436; discovery rows {len(rows)}; locus clusters {len(cls)} (tested rows)\n\nBy direction\n")
    fh.write(pd.crosstab([rows.discovery, rows.replication], rows.state).reindex(columns=ORDER, fill_value=0).to_string() + "\n\n")
    fh.write("Gene states (all 436; direct = tier 1)\n" + pd.crosstab(gs.state, gs.direct, margins=True).to_string() + "\n\n")
    fh.write("Replication fraction among tested rows, locus-cluster bootstrap 95% interval\n" + bands.round(3).to_string(index=False) + "\n\n")
    fh.write("Every row involves MVP summary statistics (dbGaP phs002453).\n")
print("done")
