"""Dynamic-range composition coding: merge lineages (arm C map, then a fixed second-level map) until every
part has 5th percentile >= 1e-4 in every included cohort, under both fraction inputs. Inputs only."""
import pandas as pd, numpy as np
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
NEW = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = ROOT + "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
new = pd.read_csv(NEW, sep="\t").set_index("sample_id"); ds = new.dataset
ct = [c for c in new.columns if c != "dataset"]
old = pd.concat([pd.read_csv(OLD.format(g=g), sep="\t", index_col=0) for g in DS]).loc[new.index, ct]
new = new[ct]
LIN = {"hepatocytes": ["Hepatocytes"], "cholangiocytes": ["Cholangiocytes"], "endothelial": ["Endothelial cells"],
       "fibroblasts": ["Fibroblasts"], "mac_mono": ["Macrophages", "Mono+mono derived cells"],
       "other_myeloid": ["cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"],
       "T_NK": ["T cells", "Circulating NK/NKT", "Resident NK"], "B_plasma": ["B cells", "Plasma cells"]}
assert sorted(sum(LIN.values(), [])) == sorted(ct)
# arm C nearest-lineage map (RESULTS rare-type arms PRESPEC), then second level for merged groups that still fail
MAP1 = {"cholangiocytes": "hepatocytes", "other_myeloid": "mac_mono", "mac_mono": "other_myeloid", "B_plasma": "T_NK",
        "T_NK": "B_plasma", "fibroblasts": "endothelial", "endothelial": "fibroblasts", "hepatocytes": "cholangiocytes"}
MAP2 = {"stromal": "lymphoid", "lymphoid": "stromal", "myeloid": None, "hepatic": None}
def level2(name):
    if set(name.split("+")) & {"fibroblasts", "endothelial"} and not set(name.split("+")) & {"T_NK", "B_plasma"}: return "stromal"
    if set(name.split("+")) & {"T_NK", "B_plasma"} and not set(name.split("+")) & {"fibroblasts", "endothelial"}: return "lymphoid"
    if set(name.split("+")) & {"mac_mono", "other_myeloid"}: return "myeloid"
    return "hepatic"
def p05(groups, cohorts):
    out = {}
    for gname, members in groups.items():
        cols = sum((LIN[m] for m in members), [])
        out[gname] = min(fr.loc[ds == c, cols].sum(1).quantile(0.05) for fr in (new, old) for c in cohorts)
    return out
def run(cohorts):
    groups = {k: [k] for k in LIN}
    steps = []
    while True:
        q = p05(groups, cohorts)
        fail = {k: v for k, v in q.items() if v < 1e-4}
        if not fail:
            return groups, q, steps
        worst = min(fail, key=fail.get)
        if len(groups[worst]) == 1:
            tgt_lin = MAP1[worst]
            tgt = next(k for k, v in groups.items() if tgt_lin in v)
        else:
            lv = MAP2[level2(worst)]
            if lv is None:
                raise SystemExit(f"cannot merge {worst}")
            tgt = next(k for k in groups if k != worst and level2(k) == lv)
        steps.append(f"{worst} (min p05 {fail[worst]:.1e}) -> {tgt}")
        groups["+".join([tgt, worst])] = groups.pop(tgt) + groups.pop(worst)
for name, coh in [("DR-main", ["GSE130970", "GSE162694", "GSE213621"]), ("DR-sep", ["GSE126848", "GSE135251"]),
                  ("all five", DS)]:
    g, q, steps = run(coh)
    print(f"{name} cohorts {coh} units {int(ds.isin(coh).sum())}")
    for s in steps: print("   merge:", s)
    for k in g: print(f"   part {k}: {g[k]}  min p05 over cohorts x inputs {q[k]:.2e}")
