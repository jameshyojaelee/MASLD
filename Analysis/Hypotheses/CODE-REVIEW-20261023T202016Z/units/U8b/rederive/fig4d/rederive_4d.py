"""Independent re-derivation of the Figure 4D counts (80 match / 37 differ / 7 unstable).

Written from METHODS.md line 359 and the 4D legend, not from the producer:
  - unit: biological donor; runs are summed to donors with each dataset's donor_pairing.csv
  - per donor and lineage: log2(CPM + 1) of the summed pseudobulk counts
  - lineage profile: mean over donors that have the lineage
  - gene z across the six lineages (ddof = 0 here; argmax is invariant to the ddof constant)
  - program score per lineage: source-weight-weighted mean z over member genes present (>= 5 genes)
  - content lineage: argmax over the six lineages
Bootstrap: donors resampled with replacement. To allow a same-seed comparison the draw protocol
matches the one the 4D table reports (numpy default_rng(seed); joint scheme draws len(union) integers
over the sorted union of donor keys per draw; within scheme draws per lineage in a fixed lineage order).
Profiles are computed in float64 by count-weighted matrix products, not by the producer's code.

Usage: python -I rederive_4d.py OUTDIR
"""
import json
import sys
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
SC = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
HS = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
PUB = ROOT / "Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv"
BOOT = ROOT / "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/source_tables/R04-09_content_lineage_bootstrap_joint_seed20260923_B1000.tsv"
REC = ROOT / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/r4d_programs.tsv"
LINEAGES = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]
TO_LINEAGE = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "cholangiocytes": "Cholangiocytes",
              "macrophages": "Macrophages", "tcells": "T_cells"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---- run -> donor key
run2donor = {}
for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103"):
    p = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for d, runs in zip(p["donor_id"], p["rna_srrs"]):
        for r in str(runs).split(";"):
            r = r.strip()
            if r and r != "nan":
                run2donor[r] = f"{ds}_{d}"

# ---- membership (source weights; renormalized over genes present)
mem = pd.read_csv(HS / "program_membership_v2.tsv", sep="\t")
mem = mem[mem.mapped_symbol.notna()]
reg = pd.read_csv(HS / "program_registry_v2.tsv", sep="\t")
declared = {u: TO_LINEAGE[c] for u, c in zip(reg.program_uid, reg.cell_type)}
member_genes = set(mem.mapped_symbol)

# ---- donor x gene log2 CPM per lineage
logcpm, donors = {}, {}
universe = None
for ln in LINEAGES:
    x = pd.read_csv(SC / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
    n_dup_genes = int(x.index.duplicated().sum())
    x = x.groupby(level=0).sum()
    key = [run2donor.get(c, c) for c in x.columns]
    x = x.T.groupby(key).sum().T                       # donor-level counts
    lib = x.sum(axis=0)
    zero_lib = int((lib == 0).sum())
    x = x.loc[:, lib > 0]
    cpm = x.div(x.sum(axis=0), axis=1) * 1e6
    lc = np.log2(cpm + 1.0)
    universe = set(lc.index) if universe is None else universe & set(lc.index)
    logcpm[ln] = lc
    donors[ln] = list(lc.columns)
    print(f"[{ln}] donors={lc.shape[1]} genes={lc.shape[0]} dup_gene_rows={n_dup_genes} zero_lib_donors={zero_lib}", flush=True)

genes = sorted(universe & member_genes)
gidx = {g: i for i, g in enumerate(genes)}
X = {ln: logcpm[ln].loc[genes].to_numpy(dtype=np.float64) for ln in LINEAGES}   # genes x donors

# weight matrix: programs x genes, rows sum to 1 over present genes
uids = sorted(declared)
Wm = np.zeros((len(uids), len(genes)))
n_scored = {}
for i, u in enumerate(uids):
    g = mem[(mem.program_uid == u) & mem.mapped_symbol.isin(gidx)]
    w = g.groupby("mapped_symbol").source_weight.sum()
    n_scored[u] = len(w)
    if len(w) >= 5:
        Wm[i, [gidx[s] for s in w.index]] = w.values / w.values.sum()
assert all(v >= 5 for v in n_scored.values()), "program with < 5 scored genes"


def call(profile):
    """profile: genes x 6 lineages -> content lineage per program."""
    mu = profile.mean(axis=1, keepdims=True)
    sd = profile.std(axis=1, ddof=0, keepdims=True)
    Z = np.where(sd > 0, (profile - mu) / np.where(sd > 0, sd, 1.0), 0.0)
    S = Wm @ Z                                   # programs x 6
    return S.argmax(axis=1), S


full_prof = np.column_stack([X[ln].mean(axis=1) for ln in LINEAGES])
full_idx, full_S = call(full_prof)
full = pd.Series([LINEAGES[k] for k in full_idx], index=uids)
dec = pd.Series(declared).reindex(uids)

# ---- compare with the producer's published content calls
pub = pd.read_csv(PUB, sep="\t").set_index("program_uid")
agree = int((pub.content_lineage.reindex(uids) == full).sum())
zcols = [f"z_{ln}" for ln in LINEAGES]
# z ddof differs by a constant factor sqrt(6/5): compare rescaled scores
max_dz = float(np.nanmax(np.abs(pub.loc[uids, zcols].to_numpy() - full_S * np.sqrt(5 / 6))))
match = int((full == dec).sum())
hep_from_other = int(((full == "Hepatocytes") & (dec != "Hepatocytes")).sum())
tab = pd.crosstab(dec, full)
print(tab)
print(f"[full] agree with published {agree}/117; max |z diff| {max_dz:.2e}; match {match}; differ {117 - match}; "
      f"hepatocyte content from non-hepatocyte discovery {hep_from_other}", flush=True)

# ---- bootstrap
union = sorted(set().union(*[set(donors[ln]) for ln in LINEAGES]))
pos = {ln: np.array([union.index(d) for d in donors[ln]]) for ln in LINEAGES}
print(f"[donors] union {len(union)}", {ln: len(donors[ln]) for ln in LINEAGES}, flush=True)


def boot(scheme, seed, B):
    rng = np.random.default_rng(seed)
    same = np.zeros(len(uids)); per_draw_match = np.zeros(B, dtype=int)
    for b in range(B):
        cols = []
        if scheme == "joint":
            c = np.bincount(rng.integers(0, len(union), len(union)), minlength=len(union)).astype(float)
            for ln in LINEAGES:
                w = c[pos[ln]]
                cols.append(X[ln] @ w / w.sum())
        else:
            for ln in LINEAGES:
                n = X[ln].shape[1]
                w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
                cols.append(X[ln] @ w / w.sum())
        k, _ = call(np.column_stack(cols))
        lab = np.array(LINEAGES)[k]
        same += (lab == full.values)
        per_draw_match[b] = int((lab == dec.values).sum())
    frac = same / B
    stable = frac >= 0.95
    is_match = (full == dec).values
    q = np.percentile(per_draw_match, [2.5, 50, 97.5])
    row = dict(scheme=scheme, seed=seed, B=B, stable_match=int((stable & is_match).sum()),
               stable_differ=int((stable & ~is_match).sum()), unstable=int((~stable).sum()),
               per_draw_q025=float(q[0]), per_draw_median=float(q[1]), per_draw_q975=float(q[2]),
               per_draw_min=int(per_draw_match.min()), per_draw_max=int(per_draw_match.max()))
    print("[boot]", row, flush=True)
    return row, pd.Series(frac, index=uids, name=f"{scheme}_{seed}_B{B}")


if len(sys.argv) > 2 and sys.argv[2] == "--sweep":
    # seed sweep: how often each printed count changes across seeds (B = 1000 each)
    sweep = []
    for scheme in ("joint", "within"):
        for seed in range(1, int(sys.argv[3]) + 1):
            r, f = boot(scheme, seed, 1000)
            r["frac_0086633e"] = float(f["hotspot_hepatocytes_0086633e4bf81d56"])
            r["frac_2c616f34"] = float(f["hotspot_macrophages_2c616f342f705a8f"])
            sweep.append(r)
    S = pd.DataFrame(sweep)
    S.to_csv(OUT / "bootstrap_seed_sweep.tsv", sep="\t", index=False)
    print(S.groupby("scheme")[["stable_match", "stable_differ", "unstable"]].agg(lambda s: s.value_counts().to_dict()))
    sys.exit(0)

runs = [("joint", 20260923, 1000), ("joint", 20261023, 1000), ("joint", 777, 1000),
        ("within", 20260923, 1000), ("within", 20261023, 1000), ("within", 777, 1000)]
rows, fracs = [], []
for s in runs:
    r, f = boot(*s)
    rows.append(r); fracs.append(f)
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "bootstrap_summary_independent.tsv", sep="\t", index=False)

per = pd.concat([dec.rename("discovery"), full.rename("content_independent"),
                 pub.content_lineage.reindex(uids).rename("content_published"),
                 pd.Series(n_scored).reindex(uids).rename("n_scored_genes")] + fracs, axis=1)
prod = pd.read_csv(BOOT, sep="\t").set_index("program_uid")
per["producer_joint_20260923"] = prod.boot_same_call.reindex(uids)
per["abs_diff_same_seed"] = (per["joint_20260923_B1000"] - per["producer_joint_20260923"]).abs()
per.index.name = "program_uid"
per.to_csv(OUT / "fig4d_program_calls_independent.tsv", sep="\t")
tab.to_csv(OUT / "fig4d_crosstab_independent.tsv", sep="\t")
print("[same seed] programs whose same-call fraction differs from the producer:",
      int((per.abs_diff_same_seed > 0).sum()), "max", float(per.abs_diff_same_seed.max()), flush=True)
cols = [c for c in per.columns if c.endswith("_B1000")]
flip = per[(per[cols] < 0.95).any(axis=1) & (per[cols] >= 0.95).any(axis=1)]
print("[seed flips]\n", flip[["discovery", "content_independent"] + cols].to_string(), flush=True)
flip.to_csv(OUT / "fig4d_seed_flipping_programs.tsv", sep="\t")

json.dump(dict(agree_published=agree, max_abs_z_diff=max_dz, match=match, differ=117 - match,
               hepatocyte_content_from_other=hep_from_other, n_union=len(union),
               donors_per_lineage={ln: len(donors[ln]) for ln in LINEAGES}, n_genes_used=len(genes),
               inputs={str(p): sha(p) for p in [PUB, BOOT, HS / "program_membership_v2.tsv",
                                               HS / "program_registry_v2.tsv"]}),
          open(OUT / "full_call_independent.json", "w"), indent=1)
