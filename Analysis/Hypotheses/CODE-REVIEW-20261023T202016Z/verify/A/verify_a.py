#!/usr/bin/env python3
"""Verifier A: independent checks of findings U6-01 (Figure S3G control pool) and
U6-03 (Figure 1B gray-diamond p).

Matchers are written here from the stated rules, not copied from the producers.
S3G rule (recompute_genetics_state_interface.py): genetic genes in random order, each
takes the nearest unused pool gene on AveExpr, accepted only within 0.2 SD of AveExpr
(SD over the complete-case joint table); otherwise the genetic gene is unmatched.
Fig 1B rule (H57 NC-1, PRESPEC_H57.md): genetic genes in random order, each takes the
nearest unused non-genetic gene in U_eval on AveExpr, no caliper, ties to the lower side.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
BULK = R / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
BULK_REPRO = R / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
COLOC = R / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
PRINTED_PAIRS = R / "figures/candidates/fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_expression_matched_pairs.tsv"
N_ORDERS = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
MY_SEED = 20261023          # verifier seed for the 1,000 orders
PRODUCER_SEED = 20260817    # S3G producer seed, used only to try to rebuild the printed pairs
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_bulk(p):
    d = pd.read_csv(p)
    d["g"] = d["gene"].astype(str).str.split(".").str[0]
    n0 = len(d)
    d = d.drop_duplicates("g", keep="first")
    log(f"bulk {p.parent.parent.name}: rows {n0}, unique base IDs {len(d)}")
    return d


coloc = pd.read_csv(COLOC)
coloc["g"] = coloc["ensembl"].astype(str).str.split(".").str[0]
log("coloc rows", len(coloc), "unique base IDs", coloc["g"].nunique())
log("coloc states", coloc["susie_state_t12"].value_counts(dropna=False).to_dict())


def joint(bulk):
    j = bulk.merge(coloc[["g", "coloc_best_susie_pp4", "susie_state_t12"]], on="g", how="inner")
    j["susie"] = pd.to_numeric(j["coloc_best_susie_pp4"], errors="coerce")
    j["gen"] = (j["susie"] > 0.5).fillna(False)
    j["treat"] = j["treat_fdr"] < 0.05
    j["deg"] = (j["padj"] < 0.05) & (j["logFC"].abs() > 0.5)
    log(f"joint rows {len(j)}; genetic {int(j.gen.sum())}; states {j.susie_state_t12.value_counts().to_dict()}")
    return j


def match_greedy(case_ae, pool_ae, caliper, tie_left=True):
    """Return pool positions (into pool_ae sorted ascending) per case; -1 if unmatched.
    Nearest unused pool gene; within caliper or unmatched."""
    M = len(pool_ae)
    used = np.zeros(M, dtype=bool)
    out = np.full(len(case_ae), -1)
    for i, v in enumerate(case_ae):
        pos = int(np.searchsorted(pool_ae, v))
        l = pos - 1
        while l >= 0 and used[l]:
            l -= 1
        r = pos
        while r < M and used[r]:
            r += 1
        dl = v - pool_ae[l] if l >= 0 else np.inf
        dr = pool_ae[r] - v if r < M else np.inf
        c = l if (dl <= dr if tie_left else dl < dr) else r
        if min(dl, dr) > caliper:
            continue
        used[c] = True
        out[i] = c
    return out


def paired(case_call, ctrl_call):
    b = int((case_call & ~ctrl_call).sum())
    c = int((~case_call & ctrl_call).sum())
    both = int((case_call & ctrl_call).sum())
    nei = int((~case_call & ~ctrl_call).sum())
    p = stats.binomtest(b, b + c, 0.5).pvalue if b + c else np.nan
    return dict(n_pairs=len(case_call), gen_only=b, ctrl_only=c, neither=nei, both=both,
                OR=b / c if c else np.nan, mcnemar_exact_p=p)


def s3g_run(d, pool_mask, caliper, case_order_idx):
    """d: complete-case joint table. case_order_idx: row labels of genetic genes in visiting order."""
    pool = d[pool_mask].sort_values("AveExpr", kind="stable")
    pae = pool["AveExpr"].to_numpy()
    cae = d.loc[case_order_idx, "AveExpr"].to_numpy()
    pos = match_greedy(cae, pae, caliper)
    ok = pos >= 0
    case_rows = np.asarray(case_order_idx)[ok]
    ctrl_rows = pool.index.to_numpy()[pos[ok]]
    r = paired(d.loc[case_rows, "treat"].to_numpy(), d.loc[ctrl_rows, "treat"].to_numpy())
    r["ctrl_states"] = d.loc[ctrl_rows, "susie_state_t12"].value_counts().to_dict()
    return r, case_rows, ctrl_rows


# ------------------------------------------------------------------ S3G (U6-01)
res = {"inputs": {str(p): sha(p) for p in (BULK, BULK_REPRO, COLOC, PRINTED_PAIRS)}}
for label, bpath in (("adopted", BULK), ("reproduction", BULK_REPRO)):
    j = joint(load_bulk(bpath))
    d = j.dropna(subset=["AveExpr", "logFC", "padj", "treat_fdr"]).copy()
    log(f"[{label}] complete-case joint {len(d)}; AveExpr duplicates {int(d.AveExpr.duplicated().sum())}")
    cal_full = 0.2 * d["AveExpr"].std(ddof=1)
    ueval = d[d.susie_state_t12.isin(["positive", "tested_le_0.5"])]
    cal_ueval = 0.2 * ueval["AveExpr"].std(ddof=1)
    log(f"[{label}] caliper 0.2 SD: joint {cal_full:.5f}; U_eval {cal_ueval:.5f}")
    gen = d["gen"]
    pools = {
        "printed_pool_all_nongenetic": ~gen,
        "restricted_tested_le_0.5": (~gen) & (d.susie_state_t12 == "tested_le_0.5"),
    }
    cases = d[gen]

    # Try to rebuild the printed pairs: producer draws one integer per matched_ladder call
    # (official arm x2 defs, then fragment susie), so try each of the first four draws.
    rng0 = np.random.default_rng(PRODUCER_SEED)
    draws = [int(rng0.integers(1 << 31)) for _ in range(4)]
    printed = pd.read_csv(PRINTED_PAIRS, sep="\t")
    pset = set(zip(printed.genetic_gene, printed.matched_control_gene))
    for k, rs in enumerate(draws):
        order = cases.sample(frac=1.0, random_state=rs).index
        r, cr, kr = s3g_run(d, pools["printed_pool_all_nongenetic"], cal_full, order)
        mine = set(zip(d.loc[cr, "g"], d.loc[kr, "g"]))
        log(f"[{label}] producer draw {k}: {json.dumps({kk: vv for kk, vv in r.items()})}; "
            f"pairs identical to printed {len(mine & pset)}/{len(pset)}")
        if k == 2:
            res[f"{label}_printed_design_draw2"] = r
            # same visiting order, restricted pool, both calipers
            for cname, cal in (("joint_sd", cal_full), ("ueval_sd", cal_ueval)):
                rr, _, _ = s3g_run(d, pools["restricted_tested_le_0.5"], cal, order)
                res[f"{label}_restricted_draw2_{cname}"] = rr
                log(f"[{label}] restricted pool, same order, caliper {cname}: {json.dumps(rr)}")

    if label != "adopted":
        continue
    # 1,000 visiting orders, verifier seed, both pools (joint-SD caliper) and restricted with U_eval-SD caliper
    rng = np.random.default_rng(MY_SEED)
    rows = []
    for o in range(N_ORDERS):
        order = cases.index.to_numpy()[rng.permutation(len(cases))]
        for pname, cal, pm in (("printed_pool_all_nongenetic", cal_full, pools["printed_pool_all_nongenetic"]),
                               ("restricted_joint_sd", cal_full, pools["restricted_tested_le_0.5"]),
                               ("restricted_ueval_sd", cal_ueval, pools["restricted_tested_le_0.5"])):
            r, _, _ = s3g_run(d, pm, cal, order)
            r.pop("ctrl_states")
            rows.append(dict(order=o, pool=pname, **r))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "s3g_1000_orders.tsv", sep="\t", index=False)
    for pname, g in t.groupby("pool"):
        s = {"n_pairs": g.n_pairs.value_counts().to_dict(),
             "gen_only": g.gen_only.value_counts().to_dict(),
             "ctrl_only": g.ctrl_only.value_counts().to_dict(),
             "OR_min_median_max": [g.OR.min(), g.OR.median(), g.OR.max()],
             "p_min_median_max": [g.mcnemar_exact_p.min(), g.mcnemar_exact_p.median(), g.mcnemar_exact_p.max()],
             "frac_p_lt_0.05": float((g.mcnemar_exact_p < 0.05).mean())}
        res[f"s3g_1000_{pname}"] = s
        log(f"[1000 orders] {pname}: {json.dumps(s, default=float)}")

# ------------------------------------------------------------------ Fig 1B (U6-03)
j = joint(load_bulk(BULK)).sort_values("g").reset_index(drop=True)
ue = j[j.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)
gen = ue["gen"].to_numpy()
deg = ue["deg"].to_numpy()
ae = ue["AveExpr"].to_numpy()
log(f"U_eval N {len(ue)} G {gen.sum()} D {deg.sum()} shared {(gen & deg).sum()}")
bins = pd.qcut(ue["AveExpr"].rank(method="first"), 10, labels=False).to_numpy()
nk = np.array([(bins == k).sum() for k in range(10)])
Dk = np.array([deg[bins == k].sum() for k in range(10)])


def strat_exact(sel):
    Gk = np.array([sel[bins == k].sum() for k in range(10)])
    O = int((sel & deg).sum())
    E = float((Gk * Dk / nk).sum())
    pmf = np.array([1.0])
    for n_, g_, d_ in zip(nk, Gk, Dk):
        x = np.arange(0, min(g_, d_) + 1)
        pmf = np.convolve(pmf, hypergeom.pmf(x, n_, d_, g_))
    x = np.arange(len(pmf))
    p = float(min(1.0, pmf[np.abs(x - E) >= abs(O - E) - 1e-9].sum()))
    return O, E, p


O, E, p = strat_exact(gen)
log(f"genetic set: O {O} E_strat {E:.4f} exact p {p:.4f}")
res["fig1b_genetic"] = dict(O=O, E=E, p=p)
gi = np.where(gen)[0]
pool = np.where(~gen)[0]
pool = pool[np.argsort(ae[pool], kind="stable")]
pae = ae[pool]
rng = np.random.default_rng(MY_SEED + 1)
rows = []
for o in range(N_ORDERS):
    order = gi[rng.permutation(len(gi))]
    pos = match_greedy(ae[order], pae, np.inf)
    ctrl = pool[pos]
    sel = np.zeros(len(ue), dtype=bool)
    sel[ctrl] = True
    Oc, Ec, pc = strat_exact(sel)
    pr = paired(deg[order], deg[ctrl])
    rows.append(dict(order=o, O_ctrl=Oc, E_ctrl=Ec, p_strat_exact=pc, **pr))
t = pd.DataFrame(rows)
t.to_csv(OUT / "fig1b_nc_1000_orders.tsv", sep="\t", index=False)
s = {c: t[c].value_counts().to_dict() for c in ("O_ctrl", "gen_only", "ctrl_only", "both")}
s.update(E_ctrl=[t.E_ctrl.min(), t.E_ctrl.max()], p_strat=[t.p_strat_exact.min(), t.p_strat_exact.max()],
         p_mcnemar=[t.mcnemar_exact_p.min(), t.mcnemar_exact_p.median(), t.mcnemar_exact_p.max()])
# unpaired comparison of 26/426 vs 13/426 for reference
s["fisher_26_vs_13"] = float(stats.fisher_exact([[26, 400], [13, 413]])[1])
res["fig1b_nc_1000"] = s
log(f"[Fig1B NC 1000 orders] {json.dumps(s, default=float)}")

(OUT / "results.json").write_text(json.dumps(res, indent=2, default=float))
(OUT / "log.txt").write_text("\n".join(LOG) + "\n")
(OUT / "env.txt").write_text(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
                             f"scipy {stats.__name__} seeds my={MY_SEED},{MY_SEED + 1} producer={PRODUCER_SEED}\n")
