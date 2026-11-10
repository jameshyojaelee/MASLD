#!/usr/bin/env python3
"""U12 independent re-derivation of Figure 5D and its text.

Printed: 767 gene-study pairs over 436 genes (22 direct genes / 28 pairs;
424 enzyme genes / 739 pairs; 10 in both); 12 of 22 and 167 of 424 genes with
<5% of posterior in peaks in every lineage; enzyme obs/expected 2.5-3.1;
94-212 of 424 genes above their expectation; panel percentages.

Steps, written from METHODS (not from the producer code):
 1. Tier-1/2 gene-study pairs with SuSiE PP.H4 > 0.5 from the adopted COLOC r2
    master and, for comparison, from the superseded r1 master.
 2. Primary signal pair per gene-study: highest PP.H4 among eligible r2 signal
    pairs; matched to the ATAC primary-pair table by (gwas, gene, hit1, hit2).
 3. Posterior mass in lineage peaks from the variant posterior table and the
    consensus peak manifest, per source and shared.
 4. Gene-level states, the 5% floor, locus-matched expectation.
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = R / "GWAS/finemapping"
R2 = FM / "results/susie_coloc_r2_20260924T053244Z"
R1M = FM / "results/susie_coloc/susie_coloc_all_gwas.csv"
CAND = R / "Analysis/Multimodal_Program_Projection/candidates"
ATAC_R1 = CAND / "atac-context-v3-candidate-2026-08-11-r1"
PRIM_R2 = CAND / "atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
PRIM_R1 = ATAC_R1 / "genetics/context/genetic_lineage_context_primary_pairs.tsv"
AUDIT = ATAC_R1 / "genetics/context/variant_liftover_audit.tsv"
PEAKS = ATAC_R1 / "consensus_peak_manifest.tsv"
PANEL = R / "figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5D.tsv"
LINEAGES = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]
log = {}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


for p in [R2 / "susie_coloc_all_gwas.csv", R1M, R2 / "susie_coloc_signal_pairs.tsv.gz", PRIM_R2, PRIM_R1, PEAKS, PANEL]:
    log[f"sha256:{p}"] = sha(p)

tiers = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t")
t12 = tiers[tiers.tier.isin([1, 2])].set_index("study_name")["tier_label"]
log["n_tier12_studies"] = int(len(t12))


def positive_pairs(master):
    d = pd.read_csv(master, usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "method"],
                    dtype={"gwas_name": str, "gene": str, "ensembl": str, "method": str})
    d = d[d.gwas_name.isin(t12.index) & d.ensembl.str.match(r"^ENSG\d{11}$", na=False)]
    d = d[pd.to_numeric(d["PP.H4.susie"], errors="coerce") > 0.5].copy()
    d["trait_class"] = d.gwas_name.map(t12)
    return d.drop_duplicates(["gwas_name", "ensembl"])


def census(d, tag):
    g = d.groupby("trait_class")
    out = {"pairs": int(len(d)), "genes": int(d.ensembl.nunique())}
    for k, v in g:
        out[f"pairs_{k}"] = int(len(v)); out[f"genes_{k}"] = int(v.ensembl.nunique())
    both = set(d[d.trait_class == "direct_MASLD"].ensembl) & set(d[d.trait_class == "liver_enzyme"].ensembl)
    out["genes_in_both_classes"] = len(both)
    log[f"census_{tag}"] = out
    return out


pos_r2 = positive_pairs(R2 / "susie_coloc_all_gwas.csv"); census(pos_r2, "coloc_r2_master")
pos_r1 = positive_pairs(R1M); census(pos_r1, "coloc_r1_master")

prim2 = pd.read_csv(PRIM_R2, sep="\t", dtype={"signal_pair_index": int})
prim1 = pd.read_csv(PRIM_R1, sep="\t", dtype={"signal_pair_index": int})
for tag, prim, pos in [("atac_r2_table", prim2, pos_r2), ("atac_r1_table", prim1, pos_r1)]:
    keys = prim[["gwas_name", "ensembl", "trait_class"]].drop_duplicates()
    census(keys, tag)
    a = set(map(tuple, keys[["gwas_name", "ensembl"]].values))
    b = set(map(tuple, pos[["gwas_name", "ensembl"]].values))
    log[f"{tag}_vs_master: only_table"] = len(a - b)
    log[f"{tag}_vs_master: only_master"] = len(b - a)

# Primary signal pair: highest PP.H4 among eligible r2 signal pairs.
sp = pd.read_csv(R2 / "susie_coloc_signal_pairs.tsv.gz", sep="\t", dtype={"eligible": str})
sp = sp[sp.eligible.str.upper() == "TRUE"]
best = sp.sort_values("PP.H4", ascending=False).drop_duplicates(["gwas", "ensembl"])
k2 = prim2[["gwas_name", "ensembl", "gwas_signal", "eqtl_signal", "pp_h4", "signal_pair_index"]].drop_duplicates()
m = k2.merge(best, left_on=["gwas_name", "ensembl"], right_on=["gwas", "ensembl"], how="left")
log["primary_pair_missing_in_r2_signal_pairs"] = int(m["PP.H4"].isna().sum())
log["primary_pair_hit_mismatch"] = int(((m.gwas_signal != m.hit1) | (m.eqtl_signal != m.hit2)).sum())
log["primary_pair_pp4_max_abs_diff"] = float((m.pp_h4 - m["PP.H4"]).abs().max())
log["primary_pair_pp4_le_0.5"] = int((m["PP.H4"] <= 0.5).sum())
mm = m.merge(pos_r2[["gwas_name", "ensembl", "PP.H4.susie"]], on=["gwas_name", "ensembl"], how="left")
log["primary_pp4_vs_master_gene_study_max_abs_diff"] = float((mm.pp_h4 - pd.to_numeric(mm["PP.H4.susie"])).abs().max())

# Posterior mass from variants.
keyset = set(zip(k2.gwas_name, k2.ensembl, k2.signal_pair_index.astype(str)))
parts = []
cols = ["gwas_name", "ensembl", "signal_pair_index", "snp_pp_h4", "mapping_status", "hg38_chrom", "hg38_position_1based"]
for ch in pd.read_csv(AUDIT, sep="\t", usecols=cols, chunksize=3_000_000,
                      dtype={"gwas_name": str, "ensembl": str, "signal_pair_index": str,
                             "mapping_status": str, "hg38_chrom": str, "hg38_position_1based": str}):
    sel = np.fromiter((k in keyset for k in zip(ch.gwas_name, ch.ensembl, ch.signal_pair_index)), bool, len(ch))
    if sel.any():
        parts.append(ch.loc[sel])
var = pd.concat(parts, ignore_index=True); del parts
var["signal_pair_index"] = var.signal_pair_index.astype(int)
var["mapped"] = var.mapping_status == "mapped"
var["pos"] = pd.to_numeric(var.hg38_position_1based, errors="coerce")
log["variant_rows_primary"] = int(len(var))

pk = pd.read_csv(PEAKS, sep="\t", usecols=["lineage", "chrom", "start0", "end", "supported_gse244832", "supported_gse281367"],
                 dtype={"supported_gse244832": str, "supported_gse281367": str})
pk["s244"] = pk.supported_gse244832.str.upper() == "TRUE"
pk["s281"] = pk.supported_gse281367.str.upper() == "TRUE"


def interval_sets(df):
    out = {}
    for chrom, g in df.groupby("chrom"):
        iv = g[["start0", "end"]].to_numpy(); iv = iv[np.argsort(iv[:, 0])]
        s, e = [iv[0, 0]], [iv[0, 1]]
        for a, b in iv[1:]:
            if a <= e[-1]:
                e[-1] = max(e[-1], b)
            else:
                s.append(a); e.append(b)
        out[chrom] = (np.array(s), np.array(e))
    return out


def inside(chrom, pos, sets):
    res = np.zeros(len(pos), bool)
    for c in np.unique(chrom):
        if c not in sets:
            continue
        s, e = sets[c]; idx = np.where(chrom == c)[0]; p = pos[idx]
        j = np.searchsorted(s, p, side="left") - 1  # last start0 < p (1-based p in (start0, end])
        ok = j >= 0
        res[idx[ok]] = p[ok] <= e[j[ok]]
    return res


v = var[var.mapped & var.pos.notna()].copy()
v["pos"] = v.pos.astype(np.int64)
ch = v.hg38_chrom.to_numpy(); po = v.pos.to_numpy()
rows = []
mapped_mass = var[var.mapped].groupby(["gwas_name", "ensembl", "signal_pair_index"]).snp_pp_h4.sum()
total_mass = var.groupby(["gwas_name", "ensembl", "signal_pair_index"]).snp_pp_h4.sum()
n_mapped = v.groupby(["gwas_name", "ensembl", "signal_pair_index"]).size()
for lin in LINEAGES:
    L = pk[pk.lineage == lin]
    sets = {"244": interval_sets(L[L.s244]), "281": interval_sets(L[L.s281]),
            "both": interval_sets(L[L.s244 & L.s281]), "any": interval_sets(L[L.s244 | L.s281])}
    for k, s in sets.items():
        v[f"in_{k}"] = inside(ch, po, s)
    for k in ["244", "281", "both", "any"]:
        v[f"m_{k}"] = v.snp_pp_h4 * v[f"in_{k}"]
    g = v.groupby(["gwas_name", "ensembl", "signal_pair_index"])
    agg = g[["m_244", "m_281", "m_both", "m_any"]].sum()
    agg["n_any_peak"] = g["in_any"].sum()
    agg = agg.reset_index(); agg["lineage"] = lin
    rows.append(agg)
mass = pd.concat(rows, ignore_index=True)
mass["mapped_mass"] = mass.set_index(["gwas_name", "ensembl", "signal_pair_index"]).index.map(mapped_mass)
mass["total_mass"] = mass.set_index(["gwas_name", "ensembl", "signal_pair_index"]).index.map(total_mass)
mass["n_mapped"] = mass.set_index(["gwas_name", "ensembl", "signal_pair_index"]).index.map(n_mapped)
mass["exp_any"] = mass.n_any_peak / mass.n_mapped


def state(r):
    if r.mapped_mass < 0.95:
        return "untestable"
    if r.m_both >= 0.5:
        return "replicated_accessible"
    if max(r.m_244, r.m_281) >= 0.5:
        return "source_dependent"
    if max(r.m_both, r.m_244, r.m_281) > 0:
        return "partial"
    return "indeterminate"


mass["state"] = mass.apply(state, axis=1)
cmp = mass.merge(prim2, on=["gwas_name", "ensembl", "signal_pair_index", "lineage"], how="outer", indicator=True)
log["mass_join"] = cmp["_merge"].value_counts().to_dict()
log["max_abs_diff_m244"] = float((cmp.m_244 - cmp.gse244832_accessible_mass).abs().max())
log["max_abs_diff_m281"] = float((cmp.m_281 - cmp.gse281367_accessible_mass).abs().max())
log["max_abs_diff_shared"] = float((cmp.m_both - cmp.shared_accessible_mass).abs().max())
log["max_abs_diff_any"] = float((cmp.m_any - cmp.any_accessible_mass).abs().max())
log["max_abs_diff_mapped"] = float((cmp.mapped_mass - cmp.mapped_posterior_mass).abs().max())
log["state_disagreements"] = int((cmp.state != cmp.evidence_state).sum())
log["untestable_rows"] = int((mass.state == "untestable").sum())
mass = mass.merge(prim2[["gwas_name", "ensembl", "trait_class", "gene"]].drop_duplicates(), on=["gwas_name", "ensembl"])
mass.to_csv(OUT / "fig5d_pair_lineage_mass_rederived.tsv.gz", sep="\t", index=False)

# Gene level.
mass["maxsrc"] = mass[["m_244", "m_281"]].max(axis=1)
gl = mass.groupby(["trait_class", "lineage", "ensembl"]).agg(
    ge50_any=("maxsrc", lambda x: bool((x >= 0.5).any())),
    ge50_both=("m_both", lambda x: bool((x >= 0.5).any())),
    max_src=("maxsrc", "max"),
    obs_mean=("m_any", "mean"), exp_mean=("exp_any", "mean"),
    obs_max=("m_any", "max"), exp_max=("exp_any", "max")).reset_index()
gl["band"] = np.select([gl.ge50_both, gl.ge50_any, gl.max_src >= 0.05], ["both", "one_source", "5_to_50"], "under5")
counts = gl.groupby(["trait_class", "lineage"]).agg(
    n_genes=("ensembl", "size"), ge50_any=("ge50_any", "sum"), ge50_both=("ge50_both", "sum"),
    n_5_to_50=("band", lambda b: int((b == "5_to_50").sum())), n_under5=("band", lambda b: int((b == "under5").sum()))
).reset_index()
counts["pct_any"] = 100 * counts.ge50_any / counts.n_genes
counts["pct_both"] = 100 * counts.ge50_both / counts.n_genes
panel = pd.read_csv(PANEL, sep="\t")
pw = panel.pivot_table(index=["trait_class", "lineage"], columns="measure", values="n_genes_meeting").reset_index()
pw = pw.rename(columns={c: f"{c}_panel" for c in ["ge50_any", "ge50_both", "ref50_any"]})
pn = panel.drop_duplicates(["trait_class", "lineage"])[["trait_class", "lineage", "n_under50", "n_5_to_50"]]
pn = pn.rename(columns={"n_under50": "n_under50_panel", "n_5_to_50": "n_5_to_50_panel"})
counts = counts.merge(pw, on=["trait_class", "lineage"]).merge(pn, on=["trait_class", "lineage"])
counts.to_csv(OUT / "fig5d_gene_level_counts_rederived.tsv", sep="\t", index=False)
log["panel_mismatch_any"] = int((counts.ge50_any != counts.ge50_any_panel).sum())
log["panel_mismatch_both"] = int((counts.ge50_both != counts.ge50_both_panel).sum())
log["panel_mismatch_5_to_50"] = int((counts.n_5_to_50 != counts.n_5_to_50_panel).sum()) if "n_5_to_50_panel" in counts else None

# Text: genes with < 5% in every lineage (max over pairs and lineages, per source).
gmax = mass.groupby(["trait_class", "ensembl"]).maxsrc.max().reset_index()
log["genes_under5_every_lineage"] = gmax.groupby("trait_class").apply(lambda d: f"{int((d.maxsrc < 0.05).sum())} of {len(d)}").to_dict()

# Locus-matched expectation: gene mean of pair obs/exp, then mean over genes; 2,000 gene bootstrap (seed fixed).
rng = np.random.default_rng(20261023)
bg = []
for (tc, lin), d in gl.groupby(["trait_class", "lineage"]):
    o, e = d.obs_mean.to_numpy(), d.exp_mean.to_numpy()
    idx = rng.integers(0, len(o), size=(2000, len(o)))
    rb = o[idx].mean(1) / e[idx].mean(1)
    bg.append({"trait_class": tc, "lineage": lin, "n_genes": len(o), "ratio": o.mean() / e.mean(),
               "ratio_lo": np.quantile(rb, 0.025), "ratio_hi": np.quantile(rb, 0.975),
               "n_obs_gt_exp": int((o > e).sum())})
bg = pd.DataFrame(bg); bg.to_csv(OUT / "fig5d_locus_expectation_rederived.tsv", sep="\t", index=False)

(OUT / "fig5d_log.json").write_text(json.dumps(log, indent=1, default=str))
print(json.dumps(log, indent=1, default=str)); print(counts.to_string()); print(bg.to_string())
