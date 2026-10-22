"""F4abcd R04-09: uncertainty of the Fig 4D discovery-versus-content lineage counts (80 match / 37 differ).

Content-lineage call re-implemented as in review job 1041 (which reproduced the published
plc-20260907T105813Z calls 117/117): donor-collapsed per-lineage pseudobulk, log2(CPM+1),
lineage profile = mean over donors, gene z across six lineages, program score = L1-weighted
mean z over member genes (>=5 genes), argmax lineage.

Two resampling schemes, both with donors as the unit:
  joint   : one draw of donors (with replacement) from the union of donors; each lineage uses
            the drawn donors that have that lineage. Keeps a donor's lineages together.
            Primary here. Seed 20260923, 1000 draws.
  within  : donors resampled independently within each lineage (review job 1041 scheme).
            Run with seed 20261015 and 200 draws to reproduce the review numbers, and with
            seed 20260923 and 1000 draws.
Stable call: same call as the full-data call in >= 95% of draws.
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
SCDIR = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
REG = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"
MEM = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
PUB = ROOT / "Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv"
L2L = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "cholangiocytes": "Cholangiocytes",
       "macrophages": "Macrophages", "tcells": "T_cells"}
LIN = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]

s2d = {}
for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103"):
    dp = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for _, r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s != "nan":
                s2d[s.strip()] = f"{ds}_{r['donor_id']}"
logcpm = {}
for ln in LIN:
    x = pd.read_csv(SCDIR / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
    x = x.groupby(level=0).sum()
    cols = pd.Series(x.columns).map(s2d).fillna(pd.Series(x.columns)).values
    x = x.T.groupby(cols).sum().T
    cpm = x.div(x.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
    logcpm[ln] = np.log2(cpm + 1).astype(np.float32)
genes = sorted(set.intersection(*[set(v.index) for v in logcpm.values()]))
logcpm = {ln: v.loc[genes] for ln, v in logcpm.items()}
reg = pd.read_csv(REG, sep="\t"); mem = pd.read_csv(MEM, sep="\t")
W = {}
for uid, g in mem.groupby("program_uid"):
    w = g.dropna(subset=["mapped_symbol"]).copy()
    w["wt"] = pd.to_numeric(w["original_l1_weight"], errors="coerce")
    w = w.dropna(subset=["wt"])
    W[uid] = w.groupby("mapped_symbol")["wt"].sum()
declared = pd.Series({r.program_uid: L2L[r.cell_type] for r in reg.itertuples()})


def content(prof):
    P = pd.DataFrame(prof).dropna(how="all")
    Z = P.sub(P.mean(axis=1), axis=0).div(P.std(axis=1).replace(0, np.nan), axis=0)
    out = {}
    for uid, wt in W.items():
        wt = wt[wt.index.isin(Z.index)]
        if len(wt) < 5:
            out[uid] = None; continue
        wt = wt / wt.sum()
        s = (Z.loc[wt.index].T * wt.values).T.sum(axis=0)
        out[uid] = s.idxmax()
    return pd.Series(out)


full = content({ln: v.mean(axis=1) for ln, v in logcpm.items()})
pub = pd.read_csv(PUB, sep="\t").set_index("program_uid").content_lineage
agree_pub = int((full.reindex(pub.index) == pub).sum())
match_full = int((full.reindex(declared.index) == declared).sum())
print(f"[full] content call agrees with published for {agree_pub}/117; matches discovery lineage {match_full}/117", flush=True)

union = sorted(set().union(*[set(v.columns) for v in logcpm.values()]))
print(f"[donors] union {len(union)}; per lineage", {ln: int(v.shape[1]) for ln, v in logcpm.items()}, flush=True)


def run(scheme, seed, nboot):
    rng = np.random.default_rng(seed)
    calls = []
    for b in range(nboot):
        prof = {}
        if scheme == "joint":
            drawn = pd.Series(np.array(union)[rng.integers(0, len(union), len(union))])
            for ln, v in logcpm.items():
                d = drawn[drawn.isin(v.columns)]
                prof[ln] = v[d.tolist()].mean(axis=1)
        else:
            for ln, v in logcpm.items():
                j = rng.integers(0, v.shape[1], v.shape[1])
                prof[ln] = v.iloc[:, j].mean(axis=1)
        calls.append(content(prof))
    C = pd.concat(calls, axis=1)
    stab = C.eq(full, axis=0).mean(axis=1)
    dec = declared.reindex(C.index)
    match_frac = C.eq(dec, axis=0).mean(axis=1)
    per_draw = C.eq(dec, axis=0).sum(axis=0).astype(int)
    B = pd.DataFrame(dict(discovery_lineage=dec, content_full=full.reindex(C.index), boot_same_call=stab,
                          boot_match_discovery=match_frac))
    B.index.name = "program_uid"
    B.to_csv(OUT / f"content_lineage_bootstrap_{scheme}_seed{seed}_B{nboot}.tsv", sep="\t")
    ms = int(((B.content_full == B.discovery_lineage) & (B.boot_same_call >= .95)).sum())
    ds_ = int(((B.content_full != B.discovery_lineage) & (B.boot_same_call >= .95)).sum())
    q = np.percentile(per_draw, [2.5, 50, 97.5])
    row = dict(scheme=scheme, seed=seed, n_draws=nboot, full_match=match_full, full_differ=117 - match_full,
               stable_match=ms, stable_differ=ds_, unstable=117 - ms - ds_,
               per_draw_match_median=float(q[1]), per_draw_match_q025=float(q[0]), per_draw_match_q975=float(q[2]),
               per_draw_match_min=int(per_draw.min()), per_draw_match_max=int(per_draw.max()),
               stromal_ecm_same_call=float(B.loc["hotspot_hepatocytes_f05c535ae5bbc0b9", "boot_same_call"]),
               ductular_same_call=float(B.loc["hotspot_hepatocytes_48f39dd4d817a10e", "boot_same_call"]))
    print("[boot]", row, flush=True)
    return row


rows = [run("within", 20261015, 200), run("joint", 20260923, 1000), run("within", 20260923, 1000)]
pd.DataFrame(rows).to_csv(OUT / "content_lineage_bootstrap_summary.tsv", sep="\t", index=False)
json.dump(dict(agree_pub=agree_pub, match_full=match_full, n_union=len(union)), open(OUT / "full_call_check.json", "w"))
