"""F4 vs F0 design inflation and per-dataset input differences for L03-04."""
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from scipy.stats import norm

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
NEW = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = ROOT + "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv"
SM = ROOT + "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_sample_manifest.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]

new = pd.read_csv(NEW, sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]
old = pd.concat([pd.read_csv(OLD.format(g=g), sep="\t", index_col=0) for g in DS]).loc[new.index, ct]

# per-dataset difference between inputs
for g in DS:
    idx = new.index[new.dataset == g]
    d = (new.loc[idx, ct] - old.loc[idx, ct]).abs()
    print(g, "n", len(idx), "max abs diff", f"{d.to_numpy().max():.4f}",
          "mean abs diff fibroblast", f"{d['Fibroblasts'].mean():.4f}", "hepatocyte", f"{d['Hepatocytes'].mean():.4f}")

def ilr(f):
    return np.log(np.maximum(f.to_numpy(float), 1e-6)) @ (-helmert(16).T)

sm = pd.read_csv(SM, sep="\t")
s = sm[sm.contrast == "F4_vs_F0"].set_index("sample_id")
print("F4_vs_F0 units", len(s), s.group.value_counts().to_dict())
y = (s.group == "F4").astype(float).to_numpy()
X0 = pd.get_dummies(s[["dataset", "inferred_sex"]], drop_first=True).astype(float)
X0.insert(0, "int", 1.0)
X0 = X0.to_numpy()

def rss(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return ((y - X @ b) ** 2).sum()

for name, f in [("rebuilt", new), ("earlier", old)]:
    Z = ilr(f.loc[s.index, ct])
    X1 = np.hstack([X0, Z])
    r0, r1 = rss(X0, y), rss(X1, y)
    infl = np.sqrt(r0 / r1)
    print(f"F4 vs F0 {name}: partial R2 of F4 on 15 ILR = {1 - r1 / r0:.3f}; SE inflation = {infl:.3f}; residual df = {len(y) - X1.shape[1]}")
    for gene, se in [("THRB", 0.0895533366472077), ("PPARA", 0.0899269049637686)]:
        print(f"   {gene}: unadjusted SE {se:.3f}, adjusted SE about {se * infl:.3f}; SE of adjusted minus unadjusted under no composition effect about {se * np.sqrt(infl**2 - 1):.3f}")

# main contrast: SE of the paired difference (adjusted minus unadjusted)
for infl, name in [(1.157, "rebuilt"), (1.177, "earlier")]:
    sd = 0.107 * np.sqrt(infl ** 2 - 1)
    print(f"disease vs control {name}: median-call SE of (adjusted - unadjusted) about {sd:.3f}; 80% power, two-sided 0.05: {2.80 * sd:.3f} log2; Bonferroni over 1854 calls: {(norm.ppf(1 - 0.025 / 1854) + norm.ppf(0.8)) * sd:.3f}")
