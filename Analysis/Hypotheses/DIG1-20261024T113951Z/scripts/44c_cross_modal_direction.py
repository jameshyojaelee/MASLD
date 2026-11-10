"""DIG1 job 44c (exploratory; combines 41c, 43a, 44a, 44b). Written before 44a and 44b outputs were read.
Genetic sign per colocalized gene (41c: SuSiE top pair, tier 1-2 liver traits, one gene per 1 Mb locus):
+1 = the allele that raises liver expression raises the trait (worse). "Agree" = a modality moves the
gene in the direction genetics links to more injury.
Primary family (Holm over 2), Fisher exact, two-sided:
  G1  genes with hepatocyte injury coupling BH < 0.05 (44a): genetic sign vs sign of the within-donor
      coupling excess. If fewer than 20 genes pass, genes with p < 0.05 are used instead.
  G2  genes in PXD051911 with disease-vs-control Welch p < 0.05 (44b): genetic sign vs sign of d. If
      fewer than 15 genes pass, all detected genes are used.
Descriptive: G2 in PXD065892 NAFLD and PXD027722 NASH; hepatocyte disease effect (44a, p < 0.05);
bulk composition-adjusted (43a, BH < 0.05). Convergence: per gene, the number of modalities that
agree and disagree (bulk adjusted BH < 0.05; hepatocyte coupling p < 0.05; hepatocyte disease p < 0.05;
protein p < 0.05 in any of the three disease contrasts, majority sign). Null for the count of genes
with >= 3 agreeing and 0 disagreeing modalities: 10,000 permutations of the genetic sign across genes
(seed 20261024).
Usage: python -I 44c_cross_modal_direction.py <new_out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
rng = np.random.default_rng(20261024)

G = pd.read_csv(RUN / "out_direction_posthoc_controls/gene_directions_susie_top_dedup.tsv", sep="\t")
G = G.drop_duplicates("gene")
C = pd.read_csv(RUN / "out_bulk_composition_adjusted/gene_disease_composition_adjusted.tsv.gz", sep="\t")[["ensembl", "beta_adj", "bh_adj"]]
cp = pd.read_csv(RUN / "out_coloc_genes_cells/hep_within_donor_summary.tsv", sep="\t").rename(columns={"mean_excess": "coup", "p_excess": "coup_p", "bh": "coup_bh"})
cd = pd.read_csv(RUN / "out_coloc_genes_cells/hep_disease_effect.tsv", sep="\t").rename(columns={"beta": "hdis", "p_exact": "hdis_p", "bh": "hdis_bh"})
pr = pd.read_csv(RUN / "out_coloc_genes_protein/protein_group_per_gene.tsv", sep="\t")
pr["key"] = pr.cohort + " " + pr.contrast
W = pr.pivot_table(index="gene", columns="key", values=["d", "p"])
W.columns = [f"{a}|{b}" for a, b in W.columns]; W = W.reset_index()
T = G.merge(C, on="ensembl", how="left").merge(cp[["gene", "coup", "n_excess_positive", "coup_p", "coup_bh"]], on="gene", how="left")
T = T.merge(cd[["gene", "hdis", "hdis_p", "hdis_bh"]], on="gene", how="left").merge(W, on="gene", how="left")


def fisher(gs, ms):
    gs, ms = np.sign(gs), np.sign(ms)
    a = int(((gs > 0) & (ms > 0)).sum()); b = int(((gs > 0) & (ms < 0)).sum()); c = int(((gs < 0) & (ms > 0)).sum()); d = int(((gs < 0) & (ms < 0)).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return dict(n=a + b + c + d, gen_up_mod_up=a, gen_up_mod_down=b, gen_down_mod_up=c, gen_down_mod_down=d,
                frac_agree=(a + d) / max(1, a + b + c + d), expected_agree=((a + b) * (a + c) + (c + d) * (b + d)) / max(1, (a + b + c + d)) ** 2,
                odds_ratio=orr, p_fisher=p)


rows = []
s1 = T[T.coup_bh < 0.05]; rule1 = "coupling BH < 0.05"
if len(s1) < 20:
    s1 = T[T.coup_p < 0.05]; rule1 = "coupling p < 0.05 (fallback)"
rows.append(dict(id="G1", test="genetic sign vs hepatocyte injury coupling", rule=rule1, **fisher(s1.gen_sign, s1.coup)))
k51 = "PXD051911 disease vs control"
s2 = T[T[f"p|{k51}"] < 0.05]; rule2 = "PXD051911 p < 0.05"
if len(s2) < 15:
    s2 = T[T[f"d|{k51}"].notna()]; rule2 = "PXD051911 all detected (fallback)"
rows.append(dict(id="G2", test="genetic sign vs PXD051911 protein disease direction", rule=rule2, **fisher(s2.gen_sign, s2[f"d|{k51}"])))
prim = pd.DataFrame(rows)
o = np.argsort(prim.p_fisher.to_numpy()); ph = np.maximum.accumulate(prim.p_fisher.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)

desc = []
for k in sorted(pr.key.unique()):
    s = T[T[f"p|{k}"] < 0.05]
    if len(s):
        desc.append(dict(modality=f"protein {k} p<0.05", **fisher(s.gen_sign, s[f"d|{k}"])))
    s = T[T[f"d|{k}"].notna()]
    desc.append(dict(modality=f"protein {k} all detected", **fisher(s.gen_sign, s[f"d|{k}"])))
s = T[T.hdis_p < 0.05]; desc.append(dict(modality="hepatocyte disease p<0.05", **fisher(s.gen_sign, s.hdis)))
s = T[T.hdis.notna()]; desc.append(dict(modality="hepatocyte disease all", **fisher(s.gen_sign, s.hdis)))
s = T[T.coup.notna()]; desc.append(dict(modality="hepatocyte coupling all", **fisher(s.gen_sign, s.coup)))
s = T[T.coup_p < 0.05]; desc.append(dict(modality="hepatocyte coupling p<0.05", **fisher(s.gen_sign, s.coup)))
s = T[T.bh_adj < 0.05]; desc.append(dict(modality="bulk composition-adjusted BH<0.05", **fisher(s.gen_sign, s.beta_adj)))
pd.DataFrame(desc).to_csv(OUT / "descriptive_by_modality.tsv", sep="\t", index=False)

# convergence
pk = [k for k in sorted(pr.key.unique()) if "AIH" not in k]
M = pd.DataFrame(index=T.index)
M["bulk"] = np.where(T.bh_adj < 0.05, np.sign(T.beta_adj), 0)
M["cell_coupling"] = np.where(T.coup_p < 0.05, np.sign(T.coup), 0)
M["cell_disease"] = np.where(T.hdis_p < 0.05, np.sign(T.hdis), 0)
ps = np.column_stack([np.where(T[f"p|{k}"] < 0.05, np.sign(T[f"d|{k}"]), 0) for k in pk])
M["protein"] = np.sign(np.nan_to_num(ps).sum(1))
Mv = M.fillna(0).to_numpy()


def conv(gs):
    ag = ((Mv * gs[:, None]) > 0).sum(1); dg = ((Mv * gs[:, None]) < 0).sum(1)
    return ag, dg


gs = T.gen_sign.to_numpy(float)
ag, dg = conv(gs)
T["n_agree"] = ag; T["n_disagree"] = dg
T = pd.concat([T, M.add_prefix("sign_")], axis=1)
T.to_csv(OUT / "cross_modal_table.tsv", sep="\t", index=False)
obs_full = int(((ag >= 3) & (dg == 0)).sum()); obs_rev = int(((dg >= 3) & (ag == 0)).sum())
NULL = np.array([[((a >= 3) & (d == 0)).sum(), ((d >= 3) & (a == 0)).sum()] for a, d in (conv(rng.permutation(gs)) for _ in range(10000))])
cv = pd.DataFrame([dict(statistic="genes with >=3 agree and 0 disagree", observed=obs_full, null_mean=NULL[:, 0].mean(), p_upper=(1 + (NULL[:, 0] >= obs_full).sum()) / 10001),
                   dict(statistic="genes with >=3 disagree and 0 agree", observed=obs_rev, null_mean=NULL[:, 1].mean(), p_upper=(1 + (NULL[:, 1] >= obs_rev).sum()) / 10001)])
cv.to_csv(OUT / "convergence_null.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(prim.round(4).to_string()); print(pd.DataFrame(desc).round(3).to_string()); print(cv.to_string())
cols = ["gene", "gwas", "PP4", "gen_sign", "logFC", "beta_adj", "coup", "coup_p", "hdis", "hdis_p"] + [f"d|{k}" for k in pk] + ["n_agree", "n_disagree"]
print(T[(T.n_agree >= 2) | (T.n_disagree >= 2)].sort_values(["n_agree", "n_disagree"], ascending=False)[cols].round(3).to_string())
