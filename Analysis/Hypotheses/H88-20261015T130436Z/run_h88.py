#!/usr/bin/env python3
# H88 run script. Implements PRESPEC_H88.md (sha256 3785b5e3a8f4002c672dc130423d5a29cf63aa6b2e8bb2e049804dcfb534b803).
# Question: within MASLD-cohort livers (HMSMA Visium, 32 primary donors), do cross-lineage program pairs that are
# coordinated in bulk co-vary across spots more than non-coordinated pairs of the same lineage pairs?
# Unit of inference: donor (array). Spots and program pairs are not replicates.
# Pipeline steps 1-6 (scoring) are copied line for line from the gate engine count_hmsma.py, which is copied unchanged
# into this run directory and hash-checked. Steps 7-8 and every statistic follow the prespecification.
# Usage: run_h88.py OUTDIR            (final run)
#        run_h88.py OUTDIR --debug [--bscale F] [--cache-in F] [--cache-out F]
#   --debug replaces the bulk coordinated labels with a seeded within-stratum random relabeling at step 4, so the
#   outcome code runs only on permuted labels; results go to DEBUG_RESULT_FAKE_LABELS.md.
import os, sys, re, json, time, hashlib, argparse, pickle, itertools
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:8.1f}s]", *a, flush=True)


ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = os.path.dirname(os.path.abspath(__file__))
PRESPEC_SHA = "3785b5e3a8f4002c672dc130423d5a29cf63aa6b2e8bb2e049804dcfb534b803"
ORIG_PRESPEC = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H88.md"
SP = f"{ROOT}/Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled"
T4R = ("RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/"
       "workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction")
R4 = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
GOUT = f"{R4}/gate/GAP-tissue/out"
META = {
    "coord": (f"{T4R}/coordinated_edges.tsv", 1114095, "af6a58989639aa382dbdb7a39d3c155fbb933773a5e1200b87eb3c17684716a3"),
    "edges": (f"{T4R}/program_interaction_edges.tsv.gz", 1869047, "645fce85f9d9af20656de53594857999485c3de9729627e6e1aa1e294b769776"),
    "comm": (f"{T4R}/program_communities.tsv", 4459, "f4b1b6e51b6e0c3edbdcdcff80654b34ea192957e918c8de4cd5571171bc0480"),
    "mem": ("Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv",
            1639956, "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"),
    "dmap": ("Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv", 1240, "be6bb8750d07b7b8020cf11fef4e8857a3c061654380b2a0a457190108f1dbfd"),
    "engine": (f"{R4}/gate/GAP-tissue/count_hmsma.py", 6100, "5fb3ad4b96f7724a4ceaf2513eb6732a5ea84b6cb43d9b1714877831b5314c21"),
    "gate_meas": (f"{GOUT}/program_measurability.tsv", 286511, "795e00a131d28dae7c834afc6166e3d66e38f21a67aed9259b0559413502406e"),
    "gate_summary": (f"{GOUT}/summary.json", 839, "ff2d3d2c4acff18644a79e12fad53fdb11536e1ed0eb2f2e3a8c121ae8419eb5"),
    "gate_pc": (f"{GOUT}/positive_controls_per_donor.tsv", 2760, "9923b8881da4db2d36fa4fa7a2e3b06c210577deaf4c6d425856730b92230400"),
    "gate_null": (f"{GOUT}/null_pair_scale_per_donor.tsv", 2884, "68f9b8004ea015da90dbe2d4cebc1a24bc1146e48ce807fc91547522da69fef2"),
    "count_py": (f"{R4}/write/H88/count_h88.py", 3296, "f8ad4401ee913f342bbba55e5a1af37f67c6f561e5294fbf769e8894aba5f188"),
    "count_out": (f"{R4}/write/H88/count_h88.out", 1910, "9f11f4227e9bcf1b389e80bd5bd0e09d1df6eec426163bd3ead861d2c69c54d0"),
    "joint_rule": ("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/JOINT_RULE_H52_H54.md", 2677,
                   "77255d390eea70bb9a1cc3632654b711e2b61dac58a2d8b13f01d6566f386b0b"),
}
PERIC = ["CYP2E1", "GLUL", "OAT", "CYP1A2", "CYP3A4"]
PERIP = ["CPS1", "HAL", "SDS", "ASS1", "PCK1"]
CHOL = ["KRT7", "KRT19", "EPCAM", "SOX9"]
FIBM = ["COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"]
DUCT = "hotspot_hepatocytes_48f39dd4d817a10e"
COLL = "hotspot_fibroblasts_98c7913eef5f9315"
# G1 targets written in the prespecification (P and P+ coordinated / non-coordinated per lineage-pair stratum)
EXP_P = {"cholangiocytes-fibroblasts": (11, 31), "cholangiocytes-hepatocytes": (26, 42), "cholangiocytes-macrophages": (19, 43),
         "cholangiocytes-tcells": (7, 16), "fibroblasts-hepatocytes": (20, 57), "fibroblasts-macrophages": (19, 49),
         "fibroblasts-tcells": (7, 13), "hepatocytes-macrophages": (34, 71), "hepatocytes-tcells": (9, 25), "macrophages-tcells": (8, 25)}
EXP_PP = {"cholangiocytes-fibroblasts": (8, 7), "cholangiocytes-hepatocytes": (23, 12), "cholangiocytes-macrophages": (14, 11),
          "cholangiocytes-tcells": (5, 5), "fibroblasts-hepatocytes": (19, 15), "fibroblasts-macrophages": (15, 7),
          "fibroblasts-tcells": (6, 0), "hepatocytes-macrophages": (29, 18), "hepatocytes-tcells": (7, 4), "macrophages-tcells": (8, 1)}
Z80 = 1.960 + 0.842   # one-sided alpha 0.025 (first Holm step), power 0.80
Z80_NP = 1.645 + 0.842

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--debug", action="store_true")
ap.add_argument("--bscale", type=float, default=1.0)
ap.add_argument("--cache-in")
ap.add_argument("--cache-out")
A = ap.parse_args()
if not A.debug and (A.bscale != 1.0 or A.cache_in or A.cache_out):
    sys.exit("REFUSE: the final run takes no debug options")
OUT = os.path.abspath(A.out)
os.makedirs(OUT, exist_ok=True)
if os.listdir(OUT):
    sys.exit(f"REFUSE: {OUT} is not empty; nothing is overwritten")
RESULT_NAME = "DEBUG_RESULT_FAKE_LABELS.md" if A.debug else "RESULT_H88.md"


def nB(n):
    return n if A.bscale == 1.0 else max(20, int(round(n * A.bscale)))


# Seeds: one SeedSequence(20260923).spawn(15); child k is used only for the component named in the prespec.
SS = np.random.SeedSequence(20260923).spawn(15)
RNG = [np.random.default_rng(s) for s in SS]
# 0 G2 | 1 PC2 | 2 NC1 | 3 NC2 | 4 T1 (+ e_i, per-stratum) | 5 T2 | 6 donor sign-flip | 7 F2a | 8 F2b | 9 F2c
# 10 donor bootstrap | 11 NP | 12 F2d | 13 neighborhood | 14 alternative measurability

KEY = {"mode": "DEBUG (fake labels at step 4)" if A.debug else "final", "bscale": A.bscale}
DIAG = []


def diag(name, scope, key, value):
    DIAG.append(dict(diagnostic=name, scope=scope, key=key, value=value))


def wtsv(df, name):
    p = f"{OUT}/{name}"
    if os.path.exists(p):
        sys.exit(f"REFUSE: {p} exists")
    df.to_csv(p, sep="\t", index=False)
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(x, d=4):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if abs(x) < 1e-3 and x != 0:
        return f"{x:.3e}"
    return f"{x:.{d}f}"


def flush_diag():
    p = f"{OUT}/diagnostics.tsv"
    if not os.path.exists(p):
        pd.DataFrame(DIAG, columns=["diagnostic", "scope", "key", "value"]).to_csv(p, sep="\t", index=False)


def stop(rule, text):
    KEY["stop_rule"] = rule
    flush_diag()
    lines = [f"# H88 result{' (DEBUG: fake labels, not a result)' if A.debug else ''}", "",
             f"Verdict: Untestable (stop rule {rule}).", "", text, "",
             "No class statistic was computed after this rule fired. The title clause stays as the joint rule "
             "decides, without spatial support.", "", "## Numbers reached", "", "```",
             json.dumps(KEY, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o)), "```", ""]
    with open(f"{OUT}/{RESULT_NAME}", "w") as f:
        f.write("\n".join(lines))
    log(f"STOP rule {rule}: {text}")
    sys.exit(0)


# ---------------------------------------------------------------- step 0: hashes
log("STEP 0: verify input hashes", "DEBUG MODE" if A.debug else "FINAL MODE")
fails, hrows = [], []
for path in [f"{RUN}/PRESPEC_H88.md", ORIG_PRESPEC]:
    h = sha256(path)
    hrows.append(dict(role="prespec", path=path, bytes=os.path.getsize(path), sha256=h, expected=PRESPEC_SHA, ok=h == PRESPEC_SHA))
    if h != PRESPEC_SHA:
        fails.append(f"prespec {path} sha256 {h}")
prespec_text = open(f"{RUN}/PRESPEC_H88.md").read()
h = sha256(f"{RUN}/count_hmsma.py")
hrows.append(dict(role="engine_copy_in_run_dir", path=f"{RUN}/count_hmsma.py", bytes=os.path.getsize(f"{RUN}/count_hmsma.py"),
                  sha256=h, expected=META["engine"][2], ok=h == META["engine"][2]))
if h != META["engine"][2]:
    fails.append("engine copy in run dir differs")
for k, (rel, nb, hx) in META.items():
    if hx not in prespec_text:
        fails.append(f"code error: hash for {k} not found in the prespec text")
    p = f"{ROOT}/{rel}"
    size, hh = os.path.getsize(p), sha256(p)
    hrows.append(dict(role=k, path=rel, bytes=size, sha256=hh, expected=hx, ok=(hh == hx and size == nb)))
    if hh != hx or size != nb:
        fails.append(f"{k}: bytes {size} (expected {nb}) sha256 {hh}")
ARR = {}
for m in re.finditer(r"(HRA_\d+\.h5ad)\s*\|?\s*([\d,]+)\s*\|?\s*`([0-9a-f]{64})`", prespec_text):
    ARR[m.group(1)] = (int(m.group(2).replace(",", "")), m.group(3))
if len(ARR) != 35:
    fails.append(f"code error: parsed {len(ARR)} arrays from the prespec, expected 35")
for name, (nb, hx) in sorted(ARR.items()):
    p = f"{SP}/{name}"
    size, hh = os.path.getsize(p), sha256(p)
    hrows.append(dict(role="visium_array", path=p, bytes=size, sha256=hh, expected=hx, ok=(hh == hx and size == nb)))
    if hh != hx or size != nb:
        fails.append(f"{name}: bytes {size} (expected {nb}) sha256 {hh}")
wtsv(pd.DataFrame(hrows), "input_hashes.tsv")
with open(f"{OUT}/script_sha256.txt", "w") as f:
    for p in [os.path.abspath(__file__), f"{RUN}/count_hmsma.py", f"{RUN}/PRESPEC_H88.md"]:
        f.write(f"{sha256(p)}  {p}\n")
KEY["input_hash_checks"] = f"{sum(r['ok'] for r in hrows)} of {len(hrows)} match"
log(KEY["input_hash_checks"])
if fails:
    stop(1, "Input hash mismatch: " + "; ".join(fails))
log("INPUT HASHES OK")

# ---------------------------------------------------------------- metadata (engine lines)
dmap = pd.read_csv(f"{ROOT}/{META['dmap'][0]}", sep="\t")
comm = pd.read_csv(f"{ROOT}/{META['comm'][0]}", sep="\t")
feats = set(comm.feature_id)
mem = pd.read_csv(f"{ROOT}/{META['mem'][0]}", sep="\t"); mem = mem[mem.program_uid.isin(feats)]
W = {p: g.groupby("mapped_symbol").original_l1_weight.sum() for p, g in mem.groupby("program_uid")}
coord = pd.read_csv(f"{ROOT}/{META['coord'][0]}", sep="\t")
need = set(mem.mapped_symbol) | set(PERIC) | set(PERIP) | set(CHOL) | set(FIBM)   # step 8 adds the PC3 markers
PRIM = sorted(dmap[dmap.group == "primary"].array_id)
CTRL = sorted(dmap[dmap.group == "inferred_control"].array_id)
ALLARR = sorted(dmap.array_id)
assert set(ALLARR) == {n.replace(".h5ad", "") for n in ARR}, "donor map and prespec array list differ"
SEX = dict(zip(dmap.array_id, dmap.sex))

e = pd.read_csv(f"{ROOT}/{META['edges'][0]}", sep="\t")
assert e.same_cell_type.dtype == bool and e.membership_entangled.dtype == bool
e["key"] = [tuple(sorted((a, b))) for a, b in zip(e.program_a, e.program_b)]
ck = {tuple(sorted((a, b))) for a, b in zip(coord.program_a, coord.program_b)}
e["coord"] = [k in ck for k in e.key]
e["lp"] = ["-".join(sorted([a, b])) for a, b in zip(e.cell_type_a, e.cell_type_b)]
LIN = {}
for a, ca, b, cb in zip(e.program_a, e.cell_type_a, e.program_b, e.cell_type_b):
    LIN[a] = ca
    LIN[b] = cb
assert all(LIN[p] == p.split("_")[1] for p in LIN)
gm = pd.read_csv(f"{ROOT}/{META['gate_meas'][0]}", sep="\t")
g44s = gm[gm.array_id.isin(PRIM)].groupby("program").measurable_080.all()
G44 = sorted(g44s[g44s].index)
lo_df = e[e.program_a.isin(G44) & e.program_b.isin(G44) & ~e.same_cell_type & (e.n_shared > 0)]
LOPAIRS = sorted(lo_df.key)
SH = {k: set(W[k[0]].index) & set(W[k[1]].index) for k in LOPAIRS}
log(f"gate 44-program set {len(G44)}; gene-sharing cross-lineage pairs for step 7: {len(LOPAIRS)}")


def resid(s, cov):  # engine line
    b, *_ = np.linalg.lstsq(cov, s, rcond=None); r = s - cov @ b; return r


# ---------------------------------------------------------------- step 1: score all 35 arrays (pipeline steps 1-8)
if A.cache_in:
    log("DEBUG: loading pass-1 cache", A.cache_in)
    with open(A.cache_in, "rb") as f:
        PASS1, meas_rows, lo_rows = pickle.load(f)
else:
    PASS1, meas_rows, lo_rows = {}, [], []
    for arr in ALLARR:
        a = ad.read_h5ad(f"{SP}/{arr}.h5ad")
        X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
        tot = np.asarray(X.sum(1)).ravel()
        sym = a.var.gene_symbol.astype(str).values
        keep = np.isin(sym, list(need))
        Xs = X[:, keep]; s2 = sym[keep]
        u, inv = np.unique(s2, return_inverse=True)
        M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(u)))
        C = (Xs @ M).toarray().astype(np.float32)
        det = (C > 0).mean(0)
        L = np.log1p(C / tot[:, None] * 1e4)
        Z = (L - L.mean(0)) / np.where(L.std(0) > 0, L.std(0), 1)
        idx = {g: i for i, g in enumerate(u)}
        cov = np.column_stack([np.ones(len(tot)), np.log1p(a.obs.total_counts.values), np.log1p(a.obs.n_genes.values)])
        raw, scores = {}, {}
        for p, w in W.items():
            g = [x for x in w.index if x in idx and det[idx[x]] >= 0.01]
            ret = w[g].sum() / w.sum() if len(g) else 0.0
            meas_rows.append(dict(array_id=arr, program=p, n_genes=len(w), n_measured=len(g), weight_retained=ret,
                                  measurable_080=ret >= 0.80, measurable_alt=(len(g) >= 8) and (ret >= 0.20), scored=len(g) >= 8))
            if len(g) >= 8:
                ww = w[g].values / w[g].sum()
                s = Z[:, [idx[x] for x in g]] @ ww
                raw[p] = s
                scores[p] = resid(s, cov)

        def mk(gl, thr=0.01):
            g = [x for x in gl if x in idx and det[idx[x]] >= thr]
            if not g:
                return None, []
            return resid(Z[:, [idx[x] for x in g]].mean(1), cov), g

        pc_s, gpc = mk(PERIC); pp_s, gpp = mk(PERIP)
        r1 = np.corrcoef(pc_s, pp_s)[0, 1]
        ch_s, gch = mk(CHOL, 0.005); fb_s, gfb = mk(FIBM, 0.005)
        r3 = np.corrcoef(ch_s, fb_s)[0, 1] if (ch_s is not None and fb_s is not None) else np.nan
        progs = sorted(scores)
        Sres = np.column_stack([scores[p] for p in progs])
        Sraw = np.column_stack([raw[p] for p in progs])
        with np.errstate(divide="ignore"):
            Zp = np.arctanh(np.corrcoef(Sres, rowvar=False))
        np.fill_diagonal(Zp, np.nan)
        # step 7: leave-shared-genes-out scores for the cross-lineage gene-sharing pairs among the 44 programs
        for k in LOPAIRS:
            sh = SH[k]
            outs, ns = [], []
            for p in k:
                w = W[p]
                g = [x for x in w.index if x not in sh and x in idx and det[idx[x]] >= 0.01]
                ns.append(len(g))
                if len(g) < 8:
                    outs.append(None)
                else:
                    ww = w[g].values / w[g].sum()
                    outs.append(resid(Z[:, [idx[x] for x in g]] @ ww, cov))
            zlo = np.arctanh(np.corrcoef(outs[0], outs[1])[0, 1]) if all(o is not None for o in outs) else np.nan
            lo_rows.append(dict(array_id=arr, program_a=k[0], program_b=k[1], n_shared_symbols=len(sh),
                                n_genes_left_a=ns[0], n_genes_left_b=ns[1], z_leaveout=zlo))
        PASS1[arr] = dict(n=len(tot), tot_median=float(np.median(tot)), n_genes_median=float(np.median(a.obs.n_genes.values)),
                          total_counts_median=float(np.median(a.obs.total_counts.values)), cov=cov,
                          row=a.obs.array_row.values.astype(int), col=a.obs.array_col.values.astype(int),
                          spatial=np.asarray(a.obsm["spatial"], dtype=float),
                          progs=progs, pidx={p: i for i, p in enumerate(progs)}, Sres=Sres, Sraw=Sraw, Zp=Zp,
                          pc_s=pc_s, pp_s=pp_s, r1=r1, gpc=gpc, gpp=gpp, ch_s=ch_s, fb_s=fb_s, r3=r3, gch=gch, gfb=gfb)
        log(arr, len(tot), "spots | scored programs", len(progs), "| PC1 r", round(r1, 4), "| PC3 markers", len(gch), len(gfb))
        del a, X, Xs, C, L, Z
    if A.cache_out:
        with open(A.cache_out, "wb") as f:
            pickle.dump((PASS1, meas_rows, lo_rows), f, protocol=4)

meas = pd.DataFrame(meas_rows)
wtsv(meas, "scored_genes.tsv")
mp = meas[meas.array_id.isin(PRIM)].groupby("program").measurable_080.all()
PROG44 = sorted(mp[mp].index)
NPROG = len(PROG44)
PI = {p: i for i, p in enumerate(PROG44)}
PAIRS946 = list(itertools.combinations(PROG44, 2))


def zpair(arr, a, b):
    d = PASS1[arr]
    if a not in d["pidx"] or b not in d["pidx"]:
        return np.nan
    return d["Zp"][d["pidx"][a], d["pidx"][b]]


Z946 = pd.DataFrame({f"{a}|{b}": [zpair(arr, a, b) for arr in ALLARR] for a, b in PAIRS946})
Z946.insert(0, "array_id", ALLARR)
p = wtsv(Z946, "z_pairs_35x946.tsv.gz")
with open(f"{OUT}/z_pairs_35x946.tsv.gz.sha256", "w") as f:
    f.write(f"{sha256(p)}  z_pairs_35x946.tsv.gz\n")
lo = pd.DataFrame(lo_rows)
wtsv(lo, "z_leaveout_pairs.tsv.gz")
log(f"STEP 1 done: 35 arrays scored; {NPROG} programs measurable in all 32 primary; {len(PAIRS946)} pairs; "
    f"{len(LOPAIRS)} leave-out pairs, {int(lo.z_leaveout.isna().sum())} missing donor values")

# ---------------------------------------------------------------- populations (metadata only)
x = e[e.program_a.isin(PROG44) & e.program_b.isin(PROG44) & ~e.same_cell_type].copy()


def make_pop(df):
    df = df.sort_values(["lp", "key"]).reset_index(drop=True)
    df["pa"] = [k[0] for k in df.key]
    df["pb"] = [k[1] for k in df.key]
    return df


POP_P = make_pop(x[(x.n_shared == 0) & ~x.membership_entangled])
POP_PP = make_pop(POP_P[POP_P.r_disjoint > 0])
POP_F2C = make_pop(x[~x.membership_entangled])


def strata_of(df):
    names = sorted(df.lp.unique())
    return names, [np.where(df.lp.values == s)[0] for s in names]


def zmat(df, arrays):
    return np.array([[zpair(arr, a, b) for arr in arrays] for a, b in zip(df.pa, df.pb)])   # pairs x donors


ZP = zmat(POP_P, PRIM)
assert not np.isnan(ZP).any(), "missing z in population P"
NCMASK = ~POP_P.coord.values          # non-coordinated pairs of P (used before step 4 without touching coordinated z)
ZNC = ZP[NCMASK]
NC_df = POP_P[NCMASK].reset_index(drop=True)
nc_names, nc_strata = strata_of(NC_df)
KS = [int(((POP_P.lp == s) & POP_P.coord).sum()) for s in nc_names]   # coordinated count per stratum (metadata)

# ---------------------------------------------------------------- step 2: PC1, stop rule 2, diagnostics D1-D5, D8
log("STEP 2: PC1 zonation reproduction")
gpc = pd.read_csv(f"{ROOT}/{META['gate_pc'][0]}", sep="\t").set_index("array_id")
pc1 = {arr: np.arctanh(PASS1[arr]["r1"]) for arr in ALLARR}
pc1_prim = np.array([pc1[a] for a in PRIM])
pc1_gate = np.array([np.arctanh(gpc.loc[a, "pc1_pericentral_vs_periportal_r"]) for a in ALLARR])
pc1_diff = np.abs(np.array([pc1[a] for a in ALLARR]) - pc1_gate)
KEY["PC1"] = dict(mean_z=float(pc1_prim.mean()), sd_z=float(pc1_prim.std(ddof=1)), n_neg=int((pc1_prim < 0).sum()),
                  n=len(pc1_prim), max_abs_diff_vs_gate_35_arrays=float(pc1_diff.max()))
for arr, d in zip(ALLARR, pc1_diff):
    diag("D2_pc1_reproduction", arr, "abs_diff_z_vs_gate", float(d))
log("PC1", KEY["PC1"])
if not (abs(round(pc1_prim.mean(), 4) - (-0.1705)) < 1e-9 and int((pc1_prim < 0).sum()) == 31):
    stop(2, f"PC1 did not reproduce the gate: mean z {pc1_prim.mean():.6f}, {int((pc1_prim < 0).sum())} of 32 negative "
            "(required -0.1705 to 4 decimals and 31 of 32). Pipeline discrepancy.")
log("PC1 REPRODUCED (stop rule 2 passed)")

# extra reproduction checks against the gate (diagnostics only)
mg = meas.merge(gm, on=["array_id", "program"], suffixes=("", "_gate"))
diag("repro_gate_measurability", "all", "rows_compared", len(mg))
diag("repro_gate_measurability", "all", "n_measured_mismatches", int((mg.n_measured != mg.n_measured_gate).sum()))
diag("repro_gate_measurability", "all", "max_abs_diff_weight_retained", float((mg.weight_retained - mg.weight_retained_gate).abs().max()))
diag("repro_gate_measurability", "all", "measurable_080_mismatches", int((mg.measurable_080 != mg.measurable_080_gate).sum()))
diag("repro_gate_measurability", "all", "set_44_equals_gate", bool(PROG44 == G44))
alle = pd.read_csv(f"{ROOT}/{META['edges'][0]}", sep="\t", usecols=["program_a", "program_b", "n_shared", "same_cell_type"])
ckf = set(map(frozenset, zip(coord.program_a, coord.program_b)))
pool = alle[(~alle.same_cell_type) & (alle.n_shared == 0)]
pool = pool[[frozenset(xx) not in ckf for xx in zip(pool.program_a, pool.program_b)]]
pool = pool[pool.program_a.isin(feats) & pool.program_b.isin(feats)]
nullpairs = pool.sample(n=min(300, len(pool)), random_state=20260923)
gnull = pd.read_csv(f"{ROOT}/{META['gate_null'][0]}", sep="\t").set_index("array_id")
nd, sdd, cnt = [], [], []
for arr in ALLARR:
    zs = [zpair(arr, a, b) for a, b in zip(nullpairs.program_a, nullpairs.program_b)]
    zs = [v for v in zs if not np.isnan(v)]
    nd.append(abs(np.mean(zs) - gnull.loc[arr, "null_mean_z"])); sdd.append(abs(np.std(zs) - gnull.loc[arr, "null_sd_z"]))
    cnt.append(len(zs) == gnull.loc[arr, "n_null_pairs"])
diag("repro_gate_null_pairs", "all", "max_abs_diff_null_mean_z", float(max(nd)))
diag("repro_gate_null_pairs", "all", "max_abs_diff_null_sd_z", float(max(sdd)))
diag("repro_gate_null_pairs", "all", "n_pairs_match_all_arrays", bool(all(cnt)))
an = [abs(np.tanh(zpair(arr, DUCT, COLL)) - gpc.loc[arr, "pc2_ductular_vs_tgfb_collagen_r"]) for arr in ALLARR
      if not np.isnan(gpc.loc[arr, "pc2_ductular_vs_tgfb_collagen_r"])]
diag("repro_gate_anchor_pair", "all", "n_arrays_compared", len(an))
diag("repro_gate_anchor_pair", "all", "max_abs_diff_r", float(max(an)) if an else np.nan)
KEY["gate_reproduction"] = {d["key"]: d["value"] for d in DIAG if d["diagnostic"].startswith("repro_")}
log("gate reproduction", KEY["gate_reproduction"])

# D1
for arr in ALLARR:
    d = PASS1[arr]
    ng = meas[(meas.array_id == arr) & meas.program.isin(PROG44)].n_measured
    for k, v in [("spots", d["n"]), ("median_library_size", d["tot_median"]), ("median_obs_total_counts", d["total_counts_median"]),
                 ("median_genes_detected", d["n_genes_median"]), ("scored_genes_min_44", int(ng.min())), ("scored_genes_median_44", float(ng.median()))]:
        diag("D1", arr, k, v)
# D3
zall_prim = Z946[Z946.array_id.isin(PRIM)].set_index("array_id").loc[PRIM].values
for j, arr in enumerate(PRIM):
    diag("D3", arr, "mean_z_946", float(np.nanmean(zall_prim[j]))); diag("D3", arr, "sd_z_946", float(np.nanstd(zall_prim[j], ddof=1)))
    diag("D3", arr, "mean_z_372_noncoord_P", float(ZNC[:, j].mean())); diag("D3", arr, "sd_z_372_noncoord_P", float(ZNC[:, j].std(ddof=1)))
KEY["D3"] = dict(mean_over_donors_of_mean_z_946=float(np.nanmean(zall_prim)), mean_over_donors_of_mean_z_372=float(ZNC.mean()),
                 median_sd_z_372=float(np.median(ZNC.std(0, ddof=1))))


# D4 and the F2a / F2b / neighborhood score transforms (label-free)
def std_cols(S):
    sd = S.std(0)
    return (S - S.mean(0)) / np.where(sd > 0, sd, 1)


def zP_from_scores(S44):   # S44: spots x 44 in PROG44 order -> z for population P pairs
    Ss = std_cols(S44)
    R = (Ss.T @ Ss) / Ss.shape[0]
    R = np.clip(R, -1, 1)
    return np.arctanh(R[[PI[a] for a in POP_P.pa], [PI[b] for b in POP_P.pb]])


D4 = []
for arr in PRIM:
    d = PASS1[arr]
    S44 = d["Sres"][:, [d["pidx"][p] for p in PROG44]]
    Ss = std_cols(S44)
    U, s, Vt = np.linalg.svd(Ss, full_matrices=False)
    t = U[:, 0] * s[0]
    share = s[0] ** 2 / (s ** 2).sum()
    zon = d["pc_s"] - d["pp_s"]
    r = np.corrcoef(t, zon)[0, 1]
    D4.append((share, r))
    diag("D4", arr, "pc1_variance_share_44", float(share)); diag("D4", arr, "corr_pc1_with_zonation_axis", float(r))
KEY["D4"] = dict(median_pc1_share=float(np.median([v[0] for v in D4])), median_abs_corr_zonation=float(np.median([abs(v[1]) for v in D4])))
# D5: between-pair versus within-pair (across donors) variance of z among the 372 non-coordinated pairs of P
pm = ZNC.mean(1); g = ZNC.mean()
ss_b = ZNC.shape[1] * ((pm - g) ** 2).sum(); ss_w = ((ZNC - pm[:, None]) ** 2).sum()
vw = ZNC.var(1, ddof=1).mean(); vb = pm.var(ddof=1) - vw / ZNC.shape[1]
ZNCc = ZNC - ZNC.mean(0, keepdims=True)   # donor means removed
pmc = ZNCc.mean(1); vwc = ZNCc.var(1, ddof=1).mean(); vbc = pmc.var(ddof=1) - vwc / ZNC.shape[1]
KEY["D5"] = dict(ss_share_between_pairs=float(ss_b / (ss_b + ss_w)), varcomp_share_between_pairs=float(vb / (vb + vw)),
                 varcomp_share_between_pairs_donor_means_removed=float(vbc / (vbc + vwc)),
                 sd_donor_averaged_pair_z=float(pm.std(ddof=1)))
for k, v in KEY["D5"].items():
    diag("D5", "372_noncoord_P", k, v)
# D8 tissue side: additive program + stratum fit of z over the 532 pairs of P (label-free; this is the F2d fit)
p_names, p_strata = strata_of(POP_P)
Xd = np.zeros((len(POP_P), len(p_names) + NPROG))
for si, idx in enumerate(p_strata):
    Xd[idx, si] = 1
for r_, (a, b) in enumerate(zip(POP_P.pa, POP_P.pb)):
    Xd[r_, len(p_names) + PI[a]] = 1; Xd[r_, len(p_names) + PI[b]] = 1
bF, *_ = np.linalg.lstsq(Xd, ZP, rcond=None)
ZP_F2d = ZP - Xd @ bF
r2 = 1 - (ZP_F2d ** 2).sum(0) / ((ZP - ZP.mean(0)) ** 2).sum(0)
for arr, v in zip(PRIM, r2):
    diag("D8_tissue", arr, "r2_program_plus_stratum_terms", float(v))
KEY["D8_tissue_median_r2"] = float(np.median(r2))
log("D3", KEY["D3"]); log("D4", KEY["D4"]); log("D5", KEY["D5"]); log("D8 tissue median R2", KEY["D8_tissue_median_r2"])

# ---------------------------------------------------------------- step 3: G1, G2, PC2, NC1, NC2, PC3
log("STEP 3: G1")
cntP = {s: (int(((POP_P.lp == s) & POP_P.coord).sum()), int(((POP_P.lp == s) & ~POP_P.coord).sum())) for s in sorted(POP_P.lp.unique())}
cntPP = {s: (int(((POP_PP.lp == s) & POP_PP.coord).sum()), int(((POP_PP.lp == s) & ~POP_PP.coord).sum())) for s in sorted(POP_PP.lp.unique())}
n_prim_read = sum(a in PASS1 for a in PRIM)
G1 = dict(primary_arrays_read=n_prim_read, programs_measurable_all32=NPROG,
          P=(int(POP_P.coord.sum()), int((~POP_P.coord).sum())), Pplus=(int(POP_PP.coord.sum()), int((~POP_PP.coord).sum())),
          F2c=(int(POP_F2C.coord.sum()), int((~POP_F2C.coord).sum())),
          strata_P_match=cntP == EXP_P, strata_Pplus_match=cntPP == EXP_PP,
          lineage_counts={l: sum(LIN[p] == l for p in PROG44) for l in sorted(set(LIN[p] for p in PROG44))})
KEY["G1"] = G1
log("G1", G1)
ok_g1 = (n_prim_read == 32 and NPROG == 44 and PROG44 == G44 and G1["P"] == (160, 372) and G1["Pplus"] == (134, 80)
         and cntP == EXP_P and cntPP == EXP_PP)
if not ok_g1:
    stop(3, f"G1 counts differ from the prespecification: {G1}; strata P {cntP}; strata P+ {cntPP}")
# D7 (metadata)
neg = POP_P[POP_P.coord & (POP_P.r_disjoint < 0)]
KEY["D7"] = dict(n_coord_P_negative_r_disjoint=len(neg), by_stratum=neg.lp.value_counts().to_dict())
log("D7", KEY["D7"])


def draw_subsets(rng, strata, ks, B, n):
    """B labelings; in each stratum (fixed sorted-name order) a uniformly random subset of size k is labeled True."""
    M = np.zeros((B, n), dtype=bool)
    for idx, k in zip(strata, ks):
        if k == 0:
            continue
        if k == len(idx):
            M[:, idx] = True
            continue
        sel = np.argsort(rng.random((B, len(idx))), axis=1)[:, :k]
        M[np.repeat(np.arange(B), k), idx[sel].ravel()] = True
    return M


def relabel(rng, strata, lab, B):
    return draw_subsets(rng, strata, [int(lab[i].sum()) for i in strata], B, len(lab))


def Dvec(M, zbar):   # D for each labeling (rows of M); valid when no z is missing (D is linear in donor-averaged z)
    Mf = M.astype(np.float64); n1 = Mf.sum(1); n0 = M.shape[1] - n1
    return Mf @ zbar / n1 - (1.0 - Mf) @ zbar / n0


def dmat(M, Z):      # per-donor d_i for each labeling: B x donors (no missing z)
    Mf = M.astype(np.float64); n1 = Mf.sum(1)[:, None]; n0 = M.shape[1] - n1
    return Mf @ Z / n1 - (1.0 - Mf) @ Z / n0


def dmat_nan(M, Z):  # as dmat; a missing z drops that pair from its class mean in that donor
    Mf = M.astype(np.float64); O = (~np.isnan(Z)).astype(np.float64); Z0 = np.nan_to_num(Z)
    return (Mf @ Z0) / (Mf @ O) - ((1.0 - Mf) @ Z0) / ((1.0 - Mf) @ O)


zbar_nc = ZNC.mean(1)
nnc = len(NC_df)
# G2 (stream 0)
B = nB(2000)
Mg = draw_subsets(RNG[0], nc_strata, KS, B, nnc)
Dg = Dvec(Mg, zbar_nc)
KEY["G2"] = dict(B=B, SD_G2=float(Dg.std(ddof=1)), MDE_G2=float(Z80 * Dg.std(ddof=1)), mean=float(Dg.mean()))
log("G2", KEY["G2"])
# PC2 (stream 1): planted class difference
Brep, Bperm = nB(500), nB(2000)
det03, det015 = [], []
for rep in range(Brep):
    pseudo = draw_subsets(RNG[1], nc_strata, KS, 1, nnc)[0]
    Mr = draw_subsets(RNG[1], nc_strata, KS, Bperm, nnc)
    for delta, store in [(0.03, det03), (0.015, det015)]:
        zb = zbar_nc + delta * pseudo
        Dobs = Dvec(pseudo[None, :], zb)[0]
        pval = (1 + (Dvec(Mr, zb) >= Dobs).sum()) / (Bperm + 1)
        store.append(pval < 0.025)
KEY["PC2"] = dict(replicates=Brep, relabelings=Bperm, detection_delta_030=float(np.mean(det03)), detection_delta_015=float(np.mean(det015)))
log("PC2", KEY["PC2"])
if np.mean(det03) < 0.50:
    stop(4, f"PC2 detection at delta 0.03 is {np.mean(det03):.3f} (< 0.50): underpowered. MDE_G2 {KEY['G2']['MDE_G2']:.4f}; "
            f"detection at 0.015 {np.mean(det015):.3f}.")
# NC1 (stream 2): pseudo-labels, T1 relabeling and donor sign-flip
Blab, Bperm, Bsign = nB(1000), nB(1000), nB(10000)
t1_rej, sf_rej = [], []
for lab in range(Blab):
    pseudo = draw_subsets(RNG[2], nc_strata, KS, 1, nnc)[0]
    Mr = draw_subsets(RNG[2], nc_strata, KS, Bperm, nnc)
    Dobs = Dvec(pseudo[None, :], zbar_nc)[0]
    t1_rej.append((1 + (Dvec(Mr, zbar_nc) >= Dobs).sum()) / (Bperm + 1) < 0.05)
    d_i = dmat(pseudo[None, :], ZNC)[0]
    sg = RNG[2].integers(0, 2, size=(Bsign, len(d_i))) * 2 - 1
    sf_rej.append((1 + ((sg * d_i).mean(1) >= d_i.mean()).sum()) / (Bsign + 1) < 0.05)
KEY["NC1"] = dict(labelings=Blab, relabelings=Bperm, sign_vectors=Bsign, T1_false_positive_rate=float(np.mean(t1_rej)),
                  donor_signflip_rejection_rate=float(np.mean(sf_rej)))
log("NC1", KEY["NC1"])
if np.mean(t1_rej) > 0.08:
    stop(5, f"NC1 label-permutation false-positive rate {np.mean(t1_rej):.3f} > 0.08: the code does not implement the stated null.")
# NC2 (stream 3): spot shuffle of the second program of every pair, independent permutation per pair
labP_true = POP_P.coord.values
pa_i = np.array([PI[a] for a in POP_P.pa]); pb_i = np.array([PI[b] for b in POP_P.pb])
nshuf = 20
nc2_D, nc2_maxd = [], []
STD44 = {}
for arr in PRIM:
    d = PASS1[arr]
    STD44[arr] = std_cols(d["Sres"][:, [d["pidx"][p] for p in PROG44]])
for k in range(nshuf):
    dk = []
    for arr in PRIM:
        Ss = STD44[arr]; n = Ss.shape[0]
        perms = np.column_stack([RNG[3].permutation(n) for _ in range(len(pa_i))])
        r = (Ss[:, pa_i] * Ss[perms, pb_i[None, :]]).mean(0)
        zz = np.arctanh(r)
        dk.append(zz[labP_true].mean() - zz[~labP_true].mean())
    dk = np.array(dk)
    nc2_D.append(dk.mean()); nc2_maxd.append(np.abs(dk).max())
KEY["NC2"] = dict(shuffles=nshuf, max_abs_D=float(np.max(np.abs(nc2_D))), max_abs_d_i=float(np.max(nc2_maxd)),
                  D_values=[float(v) for v in nc2_D])
log("NC2", {k: v for k, v in KEY["NC2"].items() if k != "D_values"})
if np.max(np.abs(nc2_D)) > 0.002 or np.max(nc2_maxd) > 0.01:
    stop(6, f"NC2 spot shuffle leaks structure: max |D| {np.max(np.abs(nc2_D)):.5f} (bound 0.002), max |d_i| {np.max(nc2_maxd):.5f} (bound 0.01).")
# PC3 (portal-tract co-location)
pc3 = {arr: (np.arctanh(PASS1[arr]["r3"]) if not np.isnan(PASS1[arr]["r3"]) else np.nan) for arr in ALLARR}
v3 = np.array([pc3[a] for a in PRIM])
n_pos3 = int(np.nansum(v3 > 0))
mean3 = float(np.nanmean(v3))
PC3_PASS = bool(mean3 > 0.05 and n_pos3 >= 24)
KEY["PC3"] = dict(mean_z=mean3, sd_z=float(np.nanstd(v3, ddof=1)), n_positive=n_pos3, n_with_value=int((~np.isnan(v3)).sum()),
                  n_negative=int(np.nansum(v3 < 0)), passed=PC3_PASS, rule="mean z > 0.05 and positive in >= 24 of 32")
log("PC3", KEY["PC3"])
KEY["stop_rules_1_to_6"] = "all passed"
log("ALL STOP RULES PASSED")

# ---------------------------------------------------------------- step 4: join class labels; T1, T2, NP, donor consistency
log("STEP 4: class labels joined")
if A.debug:
    fr = np.random.default_rng(777)
    fake = {}
    names_f, strata_f = strata_of(POP_F2C)
    labf = POP_F2C.coord.values.copy()
    for idx in strata_f:
        labf[idx] = fr.permutation(labf[idx])
    fake = dict(zip(POP_F2C.key, labf))
    for df in (POP_P, POP_PP, POP_F2C):
        df["coord"] = [fake[k] for k in df.key]
    log("DEBUG: coordinated labels replaced by a within-stratum random relabeling (seed 777); no real outcome is computed")
labP = POP_P.coord.values
labPP = POP_PP.coord.values
ZPP = zmat(POP_PP, PRIM)
ZPc = zmat(POP_P, CTRL)


def ttest(rng, df, lab, Z, B, nan=False):
    names, strata = strata_of(df)
    M = relabel(rng, strata, lab, B)
    f = dmat_nan if nan else dmat
    d_obs = f(lab[None, :], Z)[0]
    D = float(d_obs.mean())
    dnull = f(M, Z)
    Dn = dnull.mean(1)
    res = dict(n_coord=int(lab.sum()), n_noncoord=int((~lab).sum()), B=B, D=D, null_mean=float(Dn.mean()),
               D_excess=float(D - Dn.mean()), null_sd=float(Dn.std(ddof=1)),
               null_q025=float(np.quantile(Dn, 0.025)), null_q975=float(np.quantile(Dn, 0.975)),
               p=float((1 + (Dn >= D).sum()) / (B + 1)))
    return res, d_obs, dnull, M, Dn


def holm(ps, alpha=0.05):
    ps = np.asarray(ps, float); m = len(ps); order = np.argsort(ps); adj = np.empty(m); run = 0.0
    for rank, j in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[j])); adj[j] = run
    return adj, adj <= alpha


T1, d1, dnull1, M1, Dn1 = ttest(RNG[4], POP_P, labP, ZP, nB(10000))
T2, d2, dnull2, M2, Dn2 = ttest(RNG[5], POP_PP, labPP, ZPP, nB(10000))
adjF1, rejF1 = holm([T1["p"], T2["p"]])
T1.update(holm_p=float(adjF1[0]), reject_holm=bool(rejF1[0]), MDE=float(Z80 * T1["null_sd"]))
T2.update(holm_p=float(adjF1[1]), reject_holm=bool(rejF1[1]), MDE=float(Z80 * T2["null_sd"]))
log("T1", T1); log("T2", T2)
# donor consistency (stream 6) on e_i = d_i minus donor i's mean d_i over the T1 relabelings
m_i = dnull1.mean(0)
e_i = d1 - m_i
Bs = nB(200000)
sg = RNG[6].integers(0, 2, size=(Bs, len(e_i))) * 2 - 1
p_dc = float((1 + ((sg * e_i).mean(1) >= e_i.mean()).sum()) / (Bs + 1))
n_epos = int((e_i > 0).sum())
DC = dict(signflip_p=p_dc, n_e_positive=n_epos, n_donors=len(e_i), sign_vectors=Bs, holds=bool(p_dc < 0.05 and n_epos >= 22),
          mean_e=float(e_i.mean()))
log("donor consistency", DC)
# NP (stream 11): permute program identities within lineage; same draws serve P and P+
LINS = sorted(set(LIN[p] for p in PROG44))
LIDX = {l: np.array([PI[p] for p in PROG44 if LIN[p] == l]) for l in LINS}
Z44 = np.full((NPROG, NPROG, len(PRIM)), np.nan)
LOMISS = np.zeros((NPROG, NPROG), dtype=int)
shared_keys = set(LOPAIRS)
lo_map = {(r.program_a, r.program_b, r.array_id): r.z_leaveout for r in lo.itertuples()}
LINARR = np.array([LIN[p] for p in PROG44])
for a, b in PAIRS946:
    if LIN[a] == LIN[b]:
        continue
    if (a, b) in shared_keys:
        v = np.array([lo_map[(a, b, arr)] for arr in PRIM])
    else:
        v = np.array([zpair(arr, a, b) for arr in PRIM])
    Z44[PI[a], PI[b]] = v; Z44[PI[b], PI[a]] = v
    LOMISS[PI[a], PI[b]] = LOMISS[PI[b], PI[a]] = int(np.isnan(v).sum())


def np_test(rng, B):
    pops = {"P": (POP_P, labP), "Pplus": (POP_PP, labPP)}
    out = {k: [] for k in pops}
    miss = {k: [] for k in pops}
    ident = np.arange(NPROG)
    obs = {}
    for k, (df, lab) in pops.items():
        pa = np.array([PI[a] for a in df.pa]); pb = np.array([PI[b] for b in df.pb])
        Zm = Z44[ident[pa], ident[pb], :]
        obs[k] = float((np.nanmean(Zm[lab], 0) - np.nanmean(Zm[~lab], 0)).mean())
    for c0 in range(0, B, 500):
        b = min(500, B - c0)
        pi = np.tile(np.arange(NPROG), (b, 1))
        for l in LINS:
            idx = LIDX[l]
            pi[:, idx] = idx[np.argsort(rng.random((b, len(idx))), axis=1)]
        for k, (df, lab) in pops.items():
            pa = np.array([PI[a] for a in df.pa]); pb = np.array([PI[b] for b in df.pb])
            A_, B_ = pi[:, pa], pi[:, pb]
            assert not (LINARR[A_] == LINARR[B_]).any(), "a mapped pair is within one lineage"
            Zm = Z44[A_, B_, :]                      # b x pairs x donors
            nmiss = np.isnan(Zm).sum((1, 2))
            assert (nmiss == LOMISS[A_, B_].sum(1)).all(), "missing z outside the leave-out rule"
            m1 = np.nanmean(Zm[:, lab, :], axis=1); m0 = np.nanmean(Zm[:, ~lab, :], axis=1)
            out[k].append((m1 - m0).mean(1)); miss[k].append(nmiss)
    res = {}
    for k in pops:
        Dn = np.concatenate(out[k]); mm = np.concatenate(miss[k])
        res[k] = dict(B=B, D=obs[k], null_mean=float(Dn.mean()), null_sd=float(Dn.std(ddof=1)),
                      D_excess_vs_NP=float(obs[k] - Dn.mean()), p=float((1 + (Dn >= obs[k]).sum()) / (B + 1)),
                      MDE_NP=float(Z80_NP * Dn.std(ddof=1)), holds=bool((1 + (Dn >= obs[k]).sum()) / (B + 1) < 0.05),
                      missing_z_total=int(mm.sum()), draws_with_missing=int((mm > 0).sum()), max_missing_per_draw=int(mm.max()),
                      null_q025=float(np.quantile(Dn, 0.025)), null_q975=float(np.quantile(Dn, 0.975)))
        res[k]["_null"] = Dn
    return res


NP = np_test(RNG[11], nB(10000))
assert abs(NP["P"]["D"] - T1["D"]) < 1e-10 and abs(NP["Pplus"]["D"] - T2["D"]) < 1e-10, "NP identity does not reproduce T1/T2 D"
log("NP P", {k: v for k, v in NP["P"].items() if k != "_null"}); log("NP P+", {k: v for k, v in NP["Pplus"].items() if k != "_null"})

T1_rej, T2_rej = T1["reject_holm"], T2["reject_holm"]
NP_P, NP_PP = NP["P"]["holds"], NP["Pplus"]["holds"]
if not T1_rej or not DC["holds"]:
    VERDICT = "Not supported" if PC3_PASS else "Not supported (PC3 failed: not informative at spot scale)"
elif not NP_P:
    VERDICT = "Program alignment not excluded"
elif T2_rej and NP_PP:
    VERDICT = "Supported" if PC3_PASS else "Supported, portal control failed"
else:
    VERDICT = "Bulk-sign only"
log("VERDICT", VERDICT)

# per-stratum D_excess with label-null interval (from the T1 relabelings)
strat_rows = []
for s, idx in zip(*strata_of(POP_P)):
    lab_s = labP[idx]
    if lab_s.all() or (~lab_s).all():
        continue
    Zs = ZP[idx]
    dobs = dmat(lab_s[None, :], Zs)[0].mean()
    dn = dmat(M1[:, idx], Zs).mean(1)
    strat_rows.append(dict(output="per_stratum_T1", population="P", stratum=s, n_coord=int(lab_s.sum()), n_noncoord=int((~lab_s).sum()),
                           D=float(dobs), null_mean=float(dn.mean()), D_excess=float(dobs - dn.mean()),
                           excess_null_q025=float(np.quantile(dn - dn.mean(), 0.025)), excess_null_q975=float(np.quantile(dn - dn.mean(), 0.975)),
                           p_nominal=float((1 + (dn >= dobs).sum()) / (len(dn) + 1))))

# ---------------------------------------------------------------- step 5: F2a, F2b, F2c, F2d
log("STEP 5: F2")
ZF2a, ZF2b = [], []
for j_, arr in enumerate(PRIM):
    d = PASS1[arr]
    S44 = d["Sres"][:, [d["pidx"][p] for p in PROG44]]
    assert np.abs(zP_from_scores(S44) - ZP[:, j_]).max() < 1e-9, "transform path does not reproduce the step-6 z"
    Ss = std_cols(S44)
    U, s, Vt = np.linalg.svd(Ss, full_matrices=False)
    t = U[:, 0] * s[0]
    Xt = np.column_stack([np.ones(len(t)), t])
    bt, *_ = np.linalg.lstsq(Xt, Ss, rcond=None)
    ZF2a.append(zP_from_scores(Ss - Xt @ bt))
    cov2 = np.column_stack([d["cov"], d["pc_s"], d["pp_s"]])
    R44 = np.column_stack([resid(d["Sraw"][:, d["pidx"][p]], cov2) for p in PROG44])
    ZF2b.append(zP_from_scores(R44))
ZF2a = np.column_stack(ZF2a); ZF2b = np.column_stack(ZF2b)
ZF2c = np.array([[lo_map[(a, b, arr)] if (a, b) in shared_keys else zpair(arr, a, b) for arr in PRIM]
                 for a, b in zip(POP_F2C.pa, POP_F2C.pb)])
labF2c = POP_F2C.coord.values
F2 = {}
F2["F2a"], *_ = ttest(RNG[7], POP_P, labP, ZF2a, nB(10000))
F2["F2b"], *_ = ttest(RNG[8], POP_P, labP, ZF2b, nB(10000))
F2["F2c"], *_ = ttest(RNG[9], POP_F2C, labF2c, ZF2c, nB(10000), nan=True)
F2["F2d"], *_ = ttest(RNG[12], POP_P, labP, ZP_F2d, nB(10000))
F2["F2c"]["missing_z"] = int(np.isnan(ZF2c).sum())
adjF2, rejF2 = holm([F2[k]["p"] for k in ["F2a", "F2b", "F2c", "F2d"]])
for k, ad_, rj in zip(["F2a", "F2b", "F2c", "F2d"], adjF2, rejF2):
    F2[k].update(holm_p=float(ad_), reject_holm=bool(rj),
                 ratio_D_excess_to_T1=float(F2[k]["D_excess"] / T1["D_excess"]) if T1["D_excess"] != 0 else np.nan,
                 diff_D_excess_minus_T1=float(F2[k]["D_excess"] - T1["D_excess"]))
    log(k, F2[k])

# ---------------------------------------------------------------- step 6: descriptive outputs
log("STEP 6: descriptive")
desc = list(strat_rows)
# neighborhood scale (stream 13)
OFFS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
ZNB, nbcount, nbdist = [], [], []
for arr in PRIM:
    d = PASS1[arr]
    pos = {(r, c): i for i, (r, c) in enumerate(zip(d["row"], d["col"]))}
    rows, cols = list(range(d["n"])), list(range(d["n"]))
    for i, (r, c) in enumerate(zip(d["row"], d["col"])):
        for dr, dc in OFFS:
            j = pos.get((r + dr, c + dc))
            if j is not None:
                rows.append(i); cols.append(j)
    Aj = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(d["n"], d["n"]))
    deg = np.asarray(Aj.sum(1)).ravel()
    nbcount.append(deg.mean() - 1)
    rr, cc = np.array(rows[d["n"]:]), np.array(cols[d["n"]:])
    if len(rr):
        dist = np.sqrt(((d["spatial"][rr] - d["spatial"][cc]) ** 2).sum(1))
        nbdist.append((dist.min(), dist.max()))
    Aj = sp.diags(1.0 / deg) @ Aj
    S44 = d["Sres"][:, [d["pidx"][p] for p in PROG44]]
    ZNB.append(zP_from_scores(Aj @ S44))
ZNB = np.column_stack(ZNB)
NB, *_ = ttest(RNG[13], POP_P, labP, ZNB, nB(10000))
NB.update(mean_neighbors_per_spot=float(np.mean(nbcount)), neighbor_distance_min=float(min(v[0] for v in nbdist)),
          neighbor_distance_max=float(max(v[1] for v in nbdist)))
desc.append(dict(output="neighborhood_scale_T1", population="P", **NB))
log("neighborhood", NB)
# alternative measurability rule (stream 14)
ma = meas[meas.array_id.isin(PRIM)].groupby("program").measurable_alt.all()
ALT = sorted(ma[ma].index)
xa = e[e.program_a.isin(ALT) & e.program_b.isin(ALT) & ~e.same_cell_type]
POP_ALT = make_pop(xa[(xa.n_shared == 0) & ~xa.membership_entangled])
if A.debug:
    POP_ALT["coord"] = [fake.get(k, False) for k in POP_ALT.key]
ZALT = zmat(POP_ALT, PRIM)
ALTR, *_ = ttest(RNG[14], POP_ALT, POP_ALT.coord.values, ZALT, nB(10000), nan=True)
ALTR.update(n_programs=len(ALT), missing_z=int(np.isnan(ZALT).sum()),
            n_strata_both_classes=int(sum(POP_ALT.coord.values[i].any() and (~POP_ALT.coord.values[i]).any() for i in strata_of(POP_ALT)[1])))
desc.append(dict(output="alternative_measurability_T1", population="P_alt", **ALTR))
log("alternative measurability", ALTR)
# donor bootstrap (stream 10): interval for D_excess = mean of e_i over resampled donors
Bb = nB(10000)
bi = RNG[10].integers(0, len(e_i), size=(Bb, len(e_i)))
boot = e_i[bi].mean(1)
BOOT = dict(B=Bb, D_excess=float(e_i.mean()), q025=float(np.quantile(boot, 0.025)), q975=float(np.quantile(boot, 0.975)))
desc.append(dict(output="donor_bootstrap_D_excess", population="P", **BOOT))
# inferred controls and sex (no test)
dc_obs = dmat_nan(labP[None, :], ZPc)[0]
dc_null = dmat_nan(M1, ZPc).mean(0)
ec = dc_obs - dc_null
for k_, arr in enumerate(CTRL):
    desc.append(dict(output="inferred_control", population="P", stratum=arr, D=float(dc_obs[k_]), D_excess=float(ec[k_]),
                     pc1_z=float(pc1[arr]), pc3_z=float(pc3[arr])))
for sx in ["female", "male"]:
    msk = np.array([SEX[a] == sx for a in PRIM])
    desc.append(dict(output="by_sex", population="P", stratum=sx, n_donors=int(msk.sum()), D=float(d1[msk].mean()), D_excess=float(e_i[msk].mean()),
                     n_e_positive=int((e_i[msk] > 0).sum())))
# D8 label side (metadata)
deg_c = pd.concat([POP_P[POP_P.coord].pa, POP_P[POP_P.coord].pb]).value_counts().reindex(PROG44).fillna(0)
deg_t = pd.concat([POP_P.pa, POP_P.pb]).value_counts().reindex(PROG44)
KEY["D8_label"] = dict(sd_per_program_coord_fraction=float((deg_c / deg_t).std(ddof=1)), max_coord_degree=int(deg_c.max()),
                       programs_with_fraction_ge_half=int(((deg_c / deg_t) >= 0.5).sum()))

# ---------------------------------------------------------------- outputs
per_donor = pd.DataFrame(dict(array_id=PRIM, sex=[SEX[a] for a in PRIM], group="primary", n_spots=[PASS1[a]["n"] for a in PRIM],
                              d_T1=d1, m_T1_null_mean=m_i, e_T1=e_i, d_T2=d2, e_T2=d2 - dnull2.mean(0),
                              d_F2a=dmat(labP[None, :], ZF2a)[0], d_F2b=dmat(labP[None, :], ZF2b)[0], d_F2d=dmat(labP[None, :], ZP_F2d)[0],
                              pc1_z=[pc1[a] for a in PRIM], pc3_z=[pc3[a] for a in PRIM],
                              pc3_chol_markers=[",".join(PASS1[a]["gch"]) for a in PRIM], pc3_fib_markers=[",".join(PASS1[a]["gfb"]) for a in PRIM]))
ctrl_df = pd.DataFrame(dict(array_id=CTRL, sex=[SEX[a] for a in CTRL], group="inferred_control", n_spots=[PASS1[a]["n"] for a in CTRL],
                            d_T1=dc_obs, m_T1_null_mean=dc_null, e_T1=ec, pc1_z=[pc1[a] for a in CTRL], pc3_z=[pc3[a] for a in CTRL],
                            pc3_chol_markers=[",".join(PASS1[a]["gch"]) for a in CTRL], pc3_fib_markers=[",".join(PASS1[a]["gfb"]) for a in CTRL]))
wtsv(pd.concat([per_donor, ctrl_df], ignore_index=True), "per_donor.tsv")
pcr = [dict(control="PC1", **{k: v for k, v in KEY["PC1"].items()}, passed=True),
       dict(control="PC2", **KEY["PC2"], passed=bool(KEY["PC2"]["detection_delta_030"] >= 0.5)),
       dict(control="PC3", **{k: v for k, v in KEY["PC3"].items() if k != "rule"})]
wtsv(pd.DataFrame(pcr), "pc_results.tsv")
ncr = [dict(control="NC1", **KEY["NC1"], passed=bool(KEY["NC1"]["T1_false_positive_rate"] <= 0.08)),
       dict(control="NC2", **{k: v for k, v in KEY["NC2"].items() if k != "D_values"}, passed=True)]
wtsv(pd.DataFrame(ncr), "nc_results.tsv")
f1 = [dict(test="T1", population="P", family="F1 (Holm over 2, alpha 0.05, one-sided)", **T1),
      dict(test="T2", population="Pplus", family="F1 (Holm over 2, alpha 0.05, one-sided)", **T2),
      dict(test="NP", population="P", family="condition (nominal, p < 0.05)", **{k: v for k, v in NP["P"].items() if k != "_null"}),
      dict(test="NP", population="Pplus", family="condition (nominal, p < 0.05)", **{k: v for k, v in NP["Pplus"].items() if k != "_null"}),
      dict(test="donor_consistency", population="P", family="condition (nominal p < 0.05 and e_i > 0 in >= 22 of 32)", **DC),
      dict(test="G2", population="P noncoord", family="power", **KEY["G2"])]
wtsv(pd.DataFrame(f1), "f1_results.tsv")
wtsv(pd.DataFrame([dict(test=k, population="F2c" if k == "F2c" else "P", family="F2 (Holm over 4, alpha 0.05, one-sided)", **F2[k])
                   for k in ["F2a", "F2b", "F2c", "F2d"]]), "f2_results.tsv")
wtsv(pd.DataFrame(desc), "descriptive.tsv")
wtsv(pd.DataFrame(dict(draw=np.arange(len(Dn1)), T1_D_null=Dn1)), "null_T1.tsv.gz")
wtsv(pd.DataFrame(dict(draw=np.arange(len(Dn2)), T2_D_null=Dn2)), "null_T2.tsv.gz")
wtsv(pd.DataFrame(dict(draw=np.arange(len(NP["P"]["_null"])), NP_P_D_null=NP["P"]["_null"], NP_Pplus_D_null=NP["Pplus"]["_null"])), "null_NP.tsv.gz")
pd.DataFrame(DIAG + [dict(diagnostic="D6", scope="all", key=k, value=v) for k, v in
                     [("SD_G2", KEY["G2"]["SD_G2"]), ("MDE_G2", KEY["G2"]["MDE_G2"]), ("MDE_T1", T1["MDE"]), ("MDE_T2", T2["MDE"]),
                      ("MDE_NP_P", NP["P"]["MDE_NP"]), ("MDE_NP_Pplus", NP["Pplus"]["MDE_NP"]),
                      ("NP_missing_z_total_P", NP["P"]["missing_z_total"]), ("NP_missing_z_total_Pplus", NP["Pplus"]["missing_z_total"])]],
             columns=["diagnostic", "scope", "key", "value"]).to_csv(f"{OUT}/diagnostics.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- RESULT
x_ = fmt
rd = []
if VERDICT == "Supported" or VERDICT == "Supported, portal control failed":
    rd.append(f"In 32 HMSMA Visium livers, program pairs from different lineages that co-vary across people in bulk also co-vary across "
              f"55-micron spots within the same liver more than non-coordinated pairs of the same lineage pairs (excess {x_(T1['D_excess'])} "
              f"Fisher z; label-permutation p {x_(T1['p'])}; program-relabeling p {x_(NP['P']['p'])}; positive in {n_epos} of 32 livers), "
              f"including among pairs with positive bulk covariance (T2 {x_(T2['D_excess'])}, p {x_(T2['p'])}).")
    if F2["F2d"]["D_excess"] < 0.5 * T1["D_excess"]:
        rd.append("Most of the excess comes from programs that co-vary with many partners in tissue and also have many coordinated "
                  f"partners in bulk, not from which specific pairs are joined (F2d excess {x_(F2['F2d']['D_excess'])}).")
    if VERDICT != "Supported":
        rd.append("PC3 (portal-tract co-location) failed; this does not count as spatial evidence for fixed condition 1.")
elif VERDICT == "Bulk-sign only":
    if T2_rej and not NP_PP:
        rd.append(f"Program pairs with positive bulk covariance co-vary across spots, but coordinated pairs do not exceed other positively "
                  f"covarying pairs beyond what within-lineage program relabeling gives (T2 {x_(T2['D_excess'])}, p {x_(T2['p'])}; "
                  f"program-relabeling p {x_(NP['Pplus']['p'])}; an excess of {x_(NP['Pplus']['MDE_NP'])} or more would have been detected "
                  f"with 80% power). Not evidence for fixed condition 1.")
    else:
        rd.append(f"Program pairs with positive bulk covariance co-vary across spots, but coordinated pairs do not exceed other positively "
                  f"covarying pairs (T2 {x_(T2['D_excess'])}, p {x_(T2['p'])}; an excess of {x_(T2['MDE'])} or more would have been "
                  f"detected with 80% power). Not evidence for fixed condition 1.")
elif VERDICT == "Program alignment not excluded":
    rd.append(f"Bulk-coordinated program pairs co-varied across spots more than non-coordinated pairs of the same lineage pairs "
              f"(excess {x_(T1['D_excess'])}; label-permutation p {x_(T1['p'])}), but no more than when program identities were shuffled "
              f"within lineage (program-relabeling p {x_(NP['P']['p'])}; an excess of {x_(NP['P']['MDE_NP'])} or more would have been "
              f"detected with 80% power). Not evidence for fixed condition 1.")
else:
    if T1_rej and not DC["holds"]:
        rd.append(f"Bulk-coordinated program pairs exceeded non-coordinated pairs on the label-permutation test (excess {x_(T1['D_excess'])}; "
                  f"p {x_(T1['p'])}) but not liver by liver (donor sign-flip p {x_(DC['signflip_p'])}; positive in {n_epos} of 32; "
                  "required p < 0.05 and at least 22). Fixed condition 1 is not met.")
    elif PC3_PASS:
        rd.append(f"Bulk-coordinated program pairs did not co-vary across spots more than non-coordinated pairs (excess {x_(T1['D_excess'])}; "
                  f"label-permutation p {x_(T1['p'])}); an excess of {x_(KEY['G2']['MDE_G2'])} or more would have been detected with 80% "
                  "power. Fixed condition 1 is not met; under any outcome, the title may not claim coordination between cell types, and "
                  "the clause reads \"maps cross-sectional multicellular remodeling\".")
    if not PC3_PASS:
        rd.append(f"PC3 failed (mean z {x_(mean3)}, positive in {n_pos3} of 32): the result is not informative at spot scale, and is "
                  "not reported as absence of co-location. Fixed condition 1 is not met.")
if F2["F2b"]["D_excess"] < 0.5 * T1["D_excess"] and T1["D_excess"] > 0 and T1_rej:
    rd.append("F2b: the zonation-adjusted excess fell by more than half, so the co-location is largely shared lobular position.")

L_ = [f"# H88 result{' (DEBUG: fake labels, not a result)' if A.debug else ''}", "",
      f"Verdict (computed in code by the prespecified rule): **{VERDICT}**", "",
      f"Run: `{OUT}`; prespec sha256 {PRESPEC_SHA}; seed 20260923 (SeedSequence spawn 15).",
      "Unit: donor (32 primary HMSMA Visium arrays, 9 female and 23 male). Spots and program pairs are not replicates.", "",
      "## Reading", ""] + rd + ["",
      "## Gates and stop rules (all passed)", "",
      f"- Stop rule 1, input hashes: {KEY['input_hash_checks']}.",
      f"- Stop rule 2, PC1 zonation: mean z {x_(KEY['PC1']['mean_z'])} (gate -0.170490), {KEY['PC1']['n_neg']} of 32 negative; "
      f"max per-array |difference| from the gate {x_(KEY['PC1']['max_abs_diff_vs_gate_35_arrays'])}.",
      f"- Stop rule 3, G1: {G1['primary_arrays_read']} arrays, {G1['programs_measurable_all32']} programs, P {G1['P'][0]}/{G1['P'][1]}, "
      f"P+ {G1['Pplus'][0]}/{G1['Pplus'][1]}, F2c {G1['F2c'][0]}/{G1['F2c'][1]}; stratum tables match.",
      f"- G2: SD_G2 {x_(KEY['G2']['SD_G2'])}, MDE_G2 {x_(KEY['G2']['MDE_G2'])} Fisher z ({KEY['G2']['B']} pool draws).",
      f"- Stop rule 4, PC2 planted difference: detection {x_(KEY['PC2']['detection_delta_030'], 3)} at delta 0.03 (stop below 0.50; "
      f"predicted 0.70 to 1.00), {x_(KEY['PC2']['detection_delta_015'], 3)} at delta 0.015 (predicted 0.25 to 0.80).",
      f"- Stop rule 5, NC1: T1 false-positive rate {x_(KEY['NC1']['T1_false_positive_rate'], 3)} (stop above 0.08; predicted 0.03 to 0.07); "
      f"donor sign-flip on d_i rejected in {x_(KEY['NC1']['donor_signflip_rejection_rate'], 3)} of pseudo-labelings (predicted 0.20 to 0.50).",
      f"- Stop rule 6, NC2 spot shuffle: max |D| {x_(KEY['NC2']['max_abs_D'], 5)} (bound 0.002), max |d_i| {x_(KEY['NC2']['max_abs_d_i'], 5)} (bound 0.01).",
      f"- PC3 portal-tract co-location: mean z {x_(mean3)}, positive in {n_pos3} of 32 (values in {KEY['PC3']['n_with_value']}); "
      f"predicted mean > 0.05 and >= 24 positive: {'passed' if PC3_PASS else 'failed'}.", "",
      "## Primary family F1 (one-sided, Holm over 2 at alpha 0.05)", "",
      "| Test | Pairs (coord / non) | D | Label-null mean | D_excess | Label-null 95% range of D | p | Holm p | Rejects | MDE |",
      "|---|---|---:|---:|---:|---|---:|---:|---|---:|",
      f"| T1, P | {T1['n_coord']} / {T1['n_noncoord']} | {x_(T1['D'])} | {x_(T1['null_mean'])} | {x_(T1['D_excess'])} | "
      f"{x_(T1['null_q025'])} to {x_(T1['null_q975'])} | {x_(T1['p'])} | {x_(T1['holm_p'])} | {T1['reject_holm']} | {x_(T1['MDE'])} |",
      f"| T2, P+ | {T2['n_coord']} / {T2['n_noncoord']} | {x_(T2['D'])} | {x_(T2['null_mean'])} | {x_(T2['D_excess'])} | "
      f"{x_(T2['null_q025'])} to {x_(T2['null_q975'])} | {x_(T2['p'])} | {x_(T2['holm_p'])} | {T2['reject_holm']} | {x_(T2['MDE'])} |", "",
      f"Donor bootstrap 95% interval for the T1 D_excess (donor sampling only): {x_(BOOT['q025'])} to {x_(BOOT['q975'])}.", "",
      "## Conditions (nominal, not family members)", "",
      f"- NP on P: D minus program-relabeling null mean {x_(NP['P']['D_excess_vs_NP'])}; p {x_(NP['P']['p'])}; holds {NP_P}; "
      f"MDE_NP {x_(NP['P']['MDE_NP'])}; missing mapped-pair z {NP['P']['missing_z_total']} over {NP['P']['B']} draws.",
      f"- NP on P+: D minus program-relabeling null mean {x_(NP['Pplus']['D_excess_vs_NP'])}; p {x_(NP['Pplus']['p'])}; holds {NP_PP}; "
      f"MDE_NP {x_(NP['Pplus']['MDE_NP'])}; missing mapped-pair z {NP['Pplus']['missing_z_total']}.",
      f"- Donor consistency: mean e_i {x_(DC['mean_e'])}; sign-flip p {x_(DC['signflip_p'])}; e_i > 0 in {n_epos} of 32; holds {DC['holds']}.", "",
      "## Secondary family F2 (one-sided, Holm over 4 at alpha 0.05; qualifies wording only)", "",
      "| Test | Pairs (coord / non) | D | D_excess | Ratio to T1 D_excess | p | Holm p | Rejects |", "|---|---|---:|---:|---:|---:|---:|---|"] + [
      f"| {k} | {F2[k]['n_coord']} / {F2[k]['n_noncoord']} | {x_(F2[k]['D'])} | {x_(F2[k]['D_excess'])} | {x_(F2[k]['ratio_D_excess_to_T1'], 2)} | "
      f"{x_(F2[k]['p'])} | {x_(F2[k]['holm_p'])} | {F2[k]['reject_holm']} |" for k in ["F2a", "F2b", "F2c", "F2d"]] + ["",
      "Predictions: F2a, F2b, F2d D_excess falls by less than half if T1 rejects; F2c D_excess within 0.02 of T1's "
      f"(observed difference {x_(F2['F2c']['diff_D_excess_minus_T1'])}).", "",
      "## Diagnostics", "",
      f"- D3: mean z over the 946 pairs {x_(KEY['D3']['mean_over_donors_of_mean_z_946'])}; over the 372 non-coordinated pairs of P "
      f"{x_(KEY['D3']['mean_over_donors_of_mean_z_372'])}; median within-donor SD {x_(KEY['D3']['median_sd_z_372'])}.",
      f"- D4: spatial first component explains a median {x_(KEY['D4']['median_pc1_share'], 3)} of program-score variance; median |r| with "
      f"the zonation axis {x_(KEY['D4']['median_abs_corr_zonation'], 3)}.",
      f"- D5: share of pair-level variance between pairs (variance components) {x_(KEY['D5']['varcomp_share_between_pairs'], 3)}; "
      f"with donor means removed {x_(KEY['D5']['varcomp_share_between_pairs_donor_means_removed'], 3)}.",
      f"- D7: coordinated pairs in P with negative bulk r_disjoint: {KEY['D7']['n_coord_P_negative_r_disjoint']}.",
      f"- D8: per-program coordinated fraction SD {x_(KEY['D8_label']['sd_per_program_coord_fraction'], 3)}; median R2 of program plus "
      f"stratum terms over the 532 pairs {x_(KEY['D8_tissue_median_r2'], 3)}.",
      f"- Gate reproduction: {json.dumps(KEY['gate_reproduction'], default=str)}.", "",
      "## Descriptive (nominal p only, no test claimed)", "",
      f"- Neighborhood scale (spot plus up to 6 hex neighbors; mean {x_(NB['mean_neighbors_per_spot'], 2)} neighbors): D {x_(NB['D'])}, "
      f"D_excess {x_(NB['D_excess'])}, nominal p {x_(NB['p'])}.",
      f"- Alternative measurability rule ({ALTR['n_programs']} programs; P_alt {ALTR['n_coord']} / {ALTR['n_noncoord']}): "
      f"D_excess {x_(ALTR['D_excess'])}, nominal p {x_(ALTR['p'])}.",
      "- Per-stratum D_excess (label-null 95% range of the excess):"] + [
      f"  - {r['stratum']} ({r['n_coord']} / {r['n_noncoord']}): {x_(r['D_excess'])} ({x_(r['excess_null_q025'])} to {x_(r['excess_null_q975'])}), nominal p {x_(r['p_nominal'])}"
      for r in strat_rows] + [
      "- Inferred controls (not in D): " + "; ".join(f"{arr} d {x_(dc_obs[k_])}, e {x_(ec[k_])}" for k_, arr in enumerate(CTRL)) + ".",
      "- By sex: " + "; ".join(f"{r['stratum']} n {r['n_donors']}, mean e {x_(r['D_excess'])}, positive {r['n_e_positive']}"
                              for r in desc if r.get("output") == "by_sex") + ".", "",
      "## Limits", "",
      "Spot-scale co-location (55 micron spots hold several cells), not shown to be in different cells or in contact. No stage or "
      "disease contrast: HMSMA has no array-to-clinical key, and 4 of the 32 are unidentified controls. Only 44 of 113 programs are "
      "measurable in all 32 donors, so P covers 160 of the 1,776 cross-lineage coordinated edges. One cohort, one platform.", "",
      f"Wall time {time.time() - T0:.0f} s."]
with open(f"{OUT}/{RESULT_NAME}", "w") as f:
    f.write("\n".join(L_) + "\n")
with open(f"{OUT}/key_numbers.json", "w") as f:
    json.dump(dict(KEY, VERDICT=VERDICT, T1=T1, T2=T2, NP_P={k: v for k, v in NP["P"].items() if k != "_null"},
                   NP_Pplus={k: v for k, v in NP["Pplus"].items() if k != "_null"}, DC=DC, F2=F2, NB=NB, ALT=ALTR, BOOT=BOOT),
              f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else (bool(o) if isinstance(o, np.bool_) else str(o)))
log("DONE", VERDICT)
