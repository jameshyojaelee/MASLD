#!/usr/bin/env python3
"""H57: expression-stratified overlap of colocalized genes and canonical DEGs inside the
multi-signal-COLOC-evaluated universe (PRESPEC_H57.md, sha256 b9fb6782...).

Order follows the prespec: hashes (S1), set sizes (G1/S3), PC-A (S2), bins and G2 (S4),
identity and role-swap checks, PC-B (S5), PC-C, NC-1 (S6), permutation null mean (S7 part),
D3 formula check and D10 cluster counts, then the primary T1 result, then D1-D10.

--debug_null permutes DEG labels within 10 AveExpr bins of U_full before the outcome
section (used only for debugging the outcome code; never for the reported run).
"""
import argparse, hashlib, json, os, sys, time
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable
import statsmodels.formula.api as smf

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H57 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RR = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
PRESPEC_SHA = "b9fb6782bd05fdecabfc5e6ffd6c5da8630ea6f6ca9235426d127813fb45817b"
C = "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
BR = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
RDS = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds"
A = "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
TIER = "GWAS/finemapping/config/gwas_trait_tier.tsv"
CAT = "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/"
IFACE = "RNA-seq/results/genetics_state_interface/20261009T142027Z_coloc_r2_corrected_controls_r8/overlap_by_arm.tsv"
INPUTS = [
    (C, "6360984750864e7ab56593bce90e8ce665d781ec72a1ba52790d6ba53e83b8d6"),
    (B, "0491ed10024388c3497dccb9f99bf669683edd5903c9d21c5308850d7c10bde6"),
    (BR, "e30509ac9c46449cfbc916191332890af7fef06e3208b34a91d7cdf90f9658e4"),
    (RDS, "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b"),
    (A, "acdfa306f1ee45383758bc2522602fc2828c4060080b7ad2d96134fc4a399887"),
    (TIER, "865f108de0782a708469b7a4e42943075de79b5ae062dc452f4050034bc15fc7"),
    (CAT + "passport_gene_index.parquet", "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b"),
    (CAT + "passport_program_membership.parquet", "082c32b2df635379d8b1b45c502c8869bbb5fa97747f619076ebeb42919515ad"),
    (CAT + "passport_program_index.parquet", "b8fd7307818755ccf68a654f9aa725e34b78d6447846f11e8f0cde0a717cb4d7"),
    (IFACE, "4925d27cb01654bdd2b2852f7781902cdc185327f200f60536d9c805f42ad83f"),
    (RR + "scratch/verify-R01-01-reproduce/reproduce.py", "fb572178b511e7f335154aec881752cbbf87372df37a41037f43379bf687a98d"),
    (RR + "scratch/verify-R07-03-reproduce/strat_overlap.py", "399bcfc2f895fd099087f9ebb4c11554de8a08eeb2003354bde7b95e64e20cf8"),
    ("figures/main/fig1_atlas_overview/fig1b_complementary_maps.pdf", "05c3bfb2c34cd1c915cb56d63bb1f6e4f6cb781f991477b551a2c2de7b9b0421"),
    ("figures/main/fig3_bulk_transcriptomics/panels/fig3g_genetic_state_complementarity.pdf", "8989f5fe71fb993aa3b4ff970a6aaaa389c8b350e9583cb213ff1cfc0f7eb344"),
    (RR + "scratch/plan-H57/count_n.py", "4fdd6285fb0e9cb2443d69fffbc7c6360ed86f6f39f8ca605036bccca979b44a"),
    (RR + "scratch/plan-H57/count_n.out", "b0f9d057a92a239e869616453577083a53e00d04616e65c3bf8dc0c256648b28"),
    (RR + "scratch/plan-H57/mde.py", "5b32742a802dee8fe4f71fb9d198e8602516df0f47507d0de5a7c656c82efc3c"),
    (RR + "scratch/plan-H57/mde.out", "759821cc3de555c741f5ef23bd34ba57c87214bbe17cffa6e5993eb49bcc028d"),
]
TOL = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20260923)
ap.add_argument("--nperm", type=int, default=10000)
ap.add_argument("--npcb", type=int, default=1000)
ap.add_argument("--nnc", type=int, default=1000)
ap.add_argument("--debug_null", action="store_true")
args = ap.parse_args()
OUT = args.out
os.chdir(ROOT)
DIAG = []          # ordered diagnostics rows (section, key, value)
RES = {"debug_null": args.debug_null, "nperm": args.nperm, "npcb": args.npcb, "nnc": args.nnc}
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def diag(section, key, value):
    DIAG.append((section, key, value))
    log(f"DIAG {section} | {key} = {value}")


def write_diag():
    pd.DataFrame(DIAG, columns=["section", "key", "value"]).to_csv(f"{OUT}/diagnostics.tsv", sep="\t", index=False)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def finish(status, verdict, lines):
    RES["status"] = status
    RES["verdict"] = verdict
    write_diag()
    with open(f"{OUT}/results.json", "w") as f:
        json.dump(RES, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    hdr = [f"# RESULT_H57", "", f"Run directory: `{OUT}`",
           f"Prespec: `{H57}/PRESPEC_H57.md` (sha256 {PRESPEC_SHA})",
           f"Debug null mode: {args.debug_null}", "", f"**Status:** {status}", "", f"**Verdict:** {verdict}", ""]
    with open(f"{OUT}/RESULT_H57.md", "w") as f:
        f.write("\n".join(hdr + lines) + "\n")
    log("STATUS", status, "|", verdict)
    sys.exit(0)


def stop(rule, msg):
    finish("stopped_by_rule", f"{rule}: {msg}. No number from this run is used.",
           [f"Stop rule {rule} fired: {msg}", "", "See diagnostics.tsv for the values printed before the stop."])


# ---------------------------------------------------------------- 1. hashes (S1)
got = sha(f"{H57}/PRESPEC_H57.md")
diag("1_hash", "PRESPEC_H57.md", got)
if got != PRESPEC_SHA:
    stop("S1", f"prespec sha256 {got} != {PRESPEC_SHA}")
bad = []
for p, h in INPUTS:
    g = sha(p)
    diag("1_hash", p, "OK" if g == h else f"MISMATCH {g}")
    if g != h:
        bad.append(p)
if bad:
    stop("S1", "sha256 mismatch for " + "; ".join(bad))

# ---------------------------------------------------------------- load
def load_bulk(path):
    d = pd.read_csv(path)
    d["g"] = d["gene"].astype(str).str.split(".").str[0]
    return d

coloc = pd.read_csv(C)
coloc["g"] = coloc["ensembl"].astype(str).str.split(".").str[0]
diag("2_sizes", "coloc rows / unique g", f"{len(coloc)} / {coloc.g.nunique()}")
coloc = coloc.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
bulk = load_bulk(B)
bulk_r = load_bulk(BR)
diag("2_sizes", "bulk rows / unique g (adopted)", f"{len(bulk)} / {bulk.g.nunique()}")

# agreement of the two DEG tables
def is_deg(df):
    return (df["padj"] < 0.05) & (df["logFC"].abs() > 0.5)

mm = bulk.merge(bulk_r, on="gene", suffixes=("", "_r"))
diag("1_hash", "DEG tables: genes adopted / reproduction / merged", f"{len(bulk)} / {len(bulk_r)} / {len(mm)}")
diag("1_hash", "DEG tables: canonical calls adopted / reproduction", f"{int(is_deg(bulk).sum())} / {int(is_deg(bulk_r).sum())}")
dis = int((is_deg(mm) != ((mm.padj_r < 0.05) & (mm.logFC_r.abs() > 0.5))).sum())
diag("1_hash", "DEG tables: calls disagreeing", dis)
diag("1_hash", "DEG tables: max |AveExpr diff|", f"{np.abs(mm.AveExpr - mm.AveExpr_r).max():.3g}")


def build(bulk_df):
    jt = bulk_df.drop_duplicates("g").merge(
        coloc[["g", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_state_t12", "susie_driving_gwas"]],
        on="g", how="inner")
    jt["susie"] = pd.to_numeric(jt.coloc_best_susie_pp4, errors="coerce")
    jt["abf"] = pd.to_numeric(jt.coloc_best_abf_pp4, errors="coerce")
    jt["gen"] = (jt.susie >= 0.5).fillna(False).astype(bool)
    jt["deg"] = is_deg(jt).astype(bool)
    return jt

# primary frames: sorted by Ensembl ID (prespec "Bins")
full = build(bulk).sort_values("g").reset_index(drop=True)
ueval = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)
ueval7 = full[full.susie.notna()].reset_index(drop=True)

# ---------------------------------------------------------------- 2. set sizes (G1, S3)
def sizes(df):
    return len(df), int(df.gen.sum()), int(df.deg.sum()), int((df.gen & df.deg).sum())

for nm, df in [("U_full", full), ("U_eval", ueval), ("U_eval7", ueval7)]:
    diag("2_sizes", f"{nm} N/G/D/shared", "/".join(map(str, sizes(df))))
ct = pd.crosstab(full.susie_state_t12.fillna("NA"), full.deg)
diag("2_sizes", "U_full susie_state x DEG", json.dumps({k: {str(c): int(v) for c, v in r.items()} for k, r in ct.iterrows()}))
diag("2_sizes", "AveExpr NA / SE NA in U_full", f"{int(full.AveExpr.isna().sum())} / {int(full.SE.isna().sum())}")
if sizes(ueval) != (3769, 426, 220, 26):
    stop("S3", f"U_eval N/G/D/shared {sizes(ueval)} != (3769, 426, 220, 26)")
diag("5_gates", "G1", "pass (3769 / 426 / 220; shared 26)")

# ---------------------------------------------------------------- 3. PC-A reproduction (S2)
jr = build(bulk_r)  # verifier row order: merge order of the reproduction table
pca = {}
N, G, D, O = sizes(jr)
orr, p = stats.fisher_exact([[O, G - O], [D - O, N - G - D + O]])
pca["U_full"] = dict(N=N, G=G, D=D, O=O, E=G * D / N, OR=orr, p=p)
f7 = jr[jr.susie.notna()]
N7, G7, D7, O7 = sizes(f7)
orr7, p7 = stats.fisher_exact([[O7, G7 - O7], [D7 - O7, N7 - G7 - D7 + O7]])
pca["U_eval7"] = dict(N=N7, G=G7, D=D7, O=O7, E=G7 * D7 / N7, OR=orr7, p=p7)
strat_exp_full = {}
for k in (5, 10, 20, 50):
    b = pd.qcut(jr.AveExpr.rank(method="first"), k, labels=False)
    g = jr.groupby(b).agg(n=("gen", "size"), G=("gen", "sum"), D=("deg", "sum"))
    strat_exp_full[k] = float((g.G * g.D / g.n).sum())
iface = pd.read_csv(IFACE, sep="\t")
irow = iface[(iface.genetic_definition == "susie_pp4_0.5") & (iface.bulk_arm == "fragment_23370")].iloc[0]
checks = [
    ("U_full N/G/D/O", (N, G, D, O) == (14704, 426, 1155, 26)),
    ("U_full expected 33.4623", round(pca["U_full"]["E"], 4) == 33.4623),
    ("U_full OR 0.7570", round(orr, 4) == 0.7570),
    ("U_full Fisher p 0.2001", round(p, 4) == 0.2001),
    ("overlap_by_arm row matches", (int(irow.n_joint_testable), int(irow.n_genetic_joint), int(irow.n_deg_joint),
                                    int(irow.n_overlap_observed)) == (N, G, D, O)
     and abs(irow.n_overlap_expected - G * D / N) < 1e-9 and abs(irow.fisher_or - orr) < 1e-9 and abs(irow.fisher_p - p) < 1e-9),
    ("U_eval7 N/G/D/O", (N7, G7, D7, O7) == (3776, 426, 220, 26)),
    ("U_eval7 expected 24.8199", round(pca["U_eval7"]["E"], 4) == 24.8199),
    ("U_eval7 OR 1.0574", round(orr7, 4) == 1.0574),
    ("U_eval7 Fisher p 0.7429", round(p7, 4) == 0.7429),
] + [(f"U_full AveExpr {k}-bin expected {t}", abs(strat_exp_full[k] - t) <= 0.01)
     for k, t in zip((5, 10, 20, 50), (22.66, 22.32, 22.13, 22.28))]
for nm, ok in checks:
    diag("3_PC-A", nm, "pass" if ok else "FAIL")
diag("3_PC-A", "values", json.dumps({**{k: {kk: round(float(vv), 6) for kk, vv in v.items()} for k, v in pca.items()},
                                     "strat_full": {k: round(v, 4) for k, v in strat_exp_full.items()}}))
RES["PC_A"] = {"pass": all(ok for _, ok in checks), "checks": {nm: bool(ok) for nm, ok in checks}}
if not RES["PC_A"]["pass"]:
    stop("S2", "PC-A reproduction failed: " + ", ".join(nm for nm, ok in checks if not ok))

# ---------------------------------------------------------------- helpers
def make_bins(df, col, k):
    """df must already be sorted by Ensembl ID; ties in col broken by that order."""
    return pd.qcut(df[col].rank(method="first"), k, labels=False).to_numpy()


def bin_margins(bins, gen, deg):
    ks = np.unique(bins)
    n = np.array([(bins == k).sum() for k in ks])
    Gk = np.array([gen[bins == k].sum() for k in ks])
    Dk = np.array([deg[bins == k].sum() for k in ks])
    Ok = np.array([(gen & deg)[bins == k].sum() for k in ks])
    return n, Gk, Dk, Ok


def exact_null(n, Gk, Dk):
    """pmf of the total overlap: convolution of per-bin hypergeometrics (pop n_k, D_k successes, G_k draws)."""
    pmf = np.array([1.0])
    for nk, gk, dk in zip(n, Gk, Dk):
        x = np.arange(0, min(gk, dk) + 1)
        pmf = np.convolve(pmf, hypergeom.pmf(x, nk, dk, gk))
    return pmf


def exact_p(pmf, E, O):
    x = np.arange(len(pmf))
    return float(min(1.0, pmf[np.abs(x - E) >= abs(O - E) - TOL].sum()))


def mh(n, Gk, Dk, Ok):
    tabs = [np.array([[o, g - o], [d - o, nn - g - d + o]], dtype=float) for nn, g, d, o in zip(n, Gk, Dk, Ok)]
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint(alpha=0.05)
    cmh = st.test_null_odds(correction=False)
    bd = st.test_equal_odds(adjust=False)
    return dict(MH_OR=float(st.oddsratio_pooled), MH_lo=float(lo), MH_hi=float(hi),
                CMH_chi2=float(cmh.statistic), CMH_p=float(cmh.pvalue), BD_chi2=float(bd.statistic), BD_p=float(bd.pvalue))


def strat_test(bins, gen, deg, with_mh=True):
    n, Gk, Dk, Ok = bin_margins(bins, gen, deg)
    E = float((Gk * Dk / n).sum())
    O = int(Ok.sum())
    pmf = exact_null(n, Gk, Dk)
    r = dict(n_strata=len(n), O=O, E_strat=E, O_over_E=O / E if E > 0 else np.nan, p_exact=exact_p(pmf, E, O),
             null_mean=float((np.arange(len(pmf)) * pmf).sum()))
    if with_mh:
        r.update(mh(n, Gk, Dk, Ok))
    return r


def perm_null(bins, gen, deg, nperm, rng):
    idx = [np.where(bins == k)[0] for k in np.unique(bins)]
    sims = np.empty(nperm, dtype=int)
    for r in range(nperm):
        tot = 0
        for ii in idx:
            tot += int((rng.permutation(gen[ii]) & deg[ii]).sum())
        sims[r] = tot
    return sims


def perm_p(sims, E, O):
    return (int((np.abs(sims - E) >= abs(O - E) - TOL).sum()) + 1) / (len(sims) + 1)


def unstrat(gen, deg):
    N = len(gen); G = int(gen.sum()); D = int(deg.sum()); O = int((gen & deg).sum())
    orr, p = stats.fisher_exact([[O, G - O], [D - O, N - G - D + O]])
    return dict(N=N, G=G, D=D, O=O, E=G * D / N, O_over_E=O / (G * D / N), fisher_OR=float(orr), fisher_p=float(p))


ss = np.random.SeedSequence(args.seed).spawn(6)
RNG = {nm: np.random.default_rng(s) for nm, s in zip(
    ["T1_perm", "PCB_enrich", "PCB_deplete", "NC1_eval", "NC1_full", "D3_perm"], ss)}

gen_e = ueval.gen.to_numpy(); deg_e = ueval.deg.to_numpy()
gen_f = full.gen.to_numpy(); deg_f = full.deg.to_numpy()

# ---------------------------------------------------------------- 4/5. bins and G2
def bin_table(df, bins, label):
    n, Gk, Dk, _ = bin_margins(bins, df.gen.to_numpy(), df.deg.to_numpy())
    edges = [(float(df.AveExpr[bins == k].min()), float(df.AveExpr[bins == k].max())) for k in np.unique(bins)]
    t = pd.DataFrame(dict(universe=label, bin=np.unique(bins), AveExpr_min=[e[0] for e in edges],
                          AveExpr_max=[e[1] for e in edges], n_k=n, G_k=Gk, D_k=Dk, G_frac=Gk / n, D_frac=Dk / n))
    return t

bins_e10 = make_bins(ueval, "AveExpr", 10)
bins_f10 = make_bins(full, "AveExpr", 10)
bt = pd.concat([bin_table(ueval, bins_e10, "U_eval"), bin_table(full, bins_f10, "U_full")])
bt.to_csv(f"{OUT}/bins_aveexpr10.tsv", sep="\t", index=False)
for r in bt.itertuples():
    diag("4_bins", f"{r.universe} bin {r.bin}", f"AveExpr [{r.AveExpr_min:.3f}, {r.AveExpr_max:.3f}] n {r.n_k} G {r.G_k} D {r.D_k} "
         f"G/n {r.G_frac:.4f} D/n {r.D_frac:.4f}")

_, Gk, Dk, _ = bin_margins(bins_e10, gen_e, deg_e)
ok10 = int(((Gk >= 5) & (Dk >= 5)).sum())
K = None
if ok10 >= 7:
    K = 10
else:
    b5 = make_bins(ueval, "AveExpr", 5)
    _, G5, D5, _ = bin_margins(b5, gen_e, deg_e)
    ok5 = int(((G5 >= 5) & (D5 >= 5)).sum())
    diag("5_gates", "G2 5-bin bins meeting rule", ok5)
    K = 5 if ok5 >= 4 else None
diag("5_gates", "G2 10-bin bins meeting rule (G>=5 and D>=5)", ok10)
diag("5_gates", "G2 result / bins used", f"{'pass' if K else 'FAIL'} / {K}")
RES["G2"] = {"bins_meeting_rule_10": ok10, "K": K}

def untestable(rule, msg):
    d1 = unstrat(gen_e, deg_e)
    RES["D1"] = d1
    finish("completed", f"untestable ({rule}): {msg}. Only unstratified U_eval numbers are reported (D1).",
           [f"Stop rule {rule}: {msg}.", "",
            f"D1 U_eval unstratified: O {d1['O']}, expected {d1['E']:.4f}, Fisher OR {d1['fisher_OR']:.4f}, nominal p {d1['fisher_p']:.4f}.",
            f"Expected fractions for the Figure 1B markers: {100 * d1['D'] / d1['N']:.2f}% of the genetic row, {100 * d1['G'] / d1['N']:.2f}% of the DEG row."])

if K is None:
    untestable("S4", "G2 failed at 5 bins")
bins_e = bins_e10 if K == 10 else make_bins(ueval, "AveExpr", K)
n_e, G_e, D_e, _ = bin_margins(bins_e, gen_e, deg_e)
PMF_E = exact_null(n_e, G_e, D_e)
E_STRAT_NULLSPACE = float((G_e * D_e / n_e).sum())  # margins only; no overlap count used

# ---------------------------------------------------------------- 6/7. identity and role swap
one = np.zeros(len(ueval), dtype=int)
n1, G1, D1, _ = bin_margins(one, gen_e, deg_e)
diag("6_identity", "one-bin E_strat - G*D/N", f"{float((G1 * D1 / n1).sum()) - G1.sum() * D1.sum() / len(ueval):.3g}")
id_ok = abs(float((G1 * D1 / n1).sum()) - G1.sum() * D1.sum() / len(ueval)) < 1e-9
pmf_swap = exact_null(n_e, D_e, G_e)
m1 = (np.arange(len(PMF_E)) * PMF_E).sum(); m2 = (np.arange(len(pmf_swap)) * pmf_swap).sum()
L = max(len(PMF_E), len(pmf_swap))
pad = lambda v: np.pad(v, (0, L - len(v)))
diag("7_roleswap", "null mean diff / max pmf diff", f"{m1 - m2:.3g} / {np.abs(pad(PMF_E) - pad(pmf_swap)).max():.3g}")
diag("7_roleswap", "null mean - analytic E_strat (exact)", f"{m1 - E_STRAT_NULLSPACE:.3g}")
rs_ok = abs(m1 - m2) < 1e-9
RES["identity_ok"] = bool(id_ok); RES["roleswap_ok"] = bool(rs_ok)

# ---------------------------------------------------------------- 8. PC-B (S5)
x_all = np.arange(len(PMF_E))
P_LOOKUP = np.array([exact_p(PMF_E, E_STRAT_NULLSPACE, xx) for xx in x_all])


def pcb(s, rng):
    idx = [np.where(bins_e == k)[0] for k in np.unique(bins_e)]
    det = 0; failed = 0; totals = []
    for _ in range(args.npcb):
        gen = np.zeros(len(gen_e), dtype=bool)
        for ii in idx:
            gen[ii] = rng.permutation(gen_e[ii])
        ok = True
        for _step in range(abs(s)):
            if s > 0:  # genetic non-DEG <-> non-genetic DEG
                A_ = [ii[gen[ii] & ~deg_e[ii]] for ii in idx]; B_ = [ii[~gen[ii] & deg_e[ii]] for ii in idx]
            else:      # genetic DEG <-> non-genetic non-DEG
                A_ = [ii[gen[ii] & deg_e[ii]] for ii in idx]; B_ = [ii[~gen[ii] & ~deg_e[ii]] for ii in idx]
            w = np.array([len(a) * len(b) for a, b in zip(A_, B_)], dtype=float)
            if w.sum() == 0:
                ok = False; break
            k = rng.choice(len(idx), p=w / w.sum())
            a = rng.choice(A_[k]); b = rng.choice(B_[k])
            gen[a] = False; gen[b] = True
        if not ok:
            failed += 1; totals.append(-1); continue
        tot = int((gen & deg_e).sum()); totals.append(tot)
        det += int(P_LOOKUP[tot] < 0.05)
    return det / args.npcb, failed, totals

rate_up, fail_up, tot_up = pcb(+14, RNG["PCB_enrich"])
rate_dn, fail_dn, tot_dn = pcb(-12, RNG["PCB_deplete"])
pd.DataFrame(dict(rep=np.arange(args.npcb), total_enrich14=tot_up, total_deplete12=tot_dn)).to_csv(
    f"{OUT}/pcb_replicates.tsv", sep="\t", index=False)
diag("8_PC-B", "enrichment s=+14 detection rate / failed swaps", f"{rate_up:.3f} / {fail_up}")
diag("8_PC-B", "depletion s=-12 detection rate / failed swaps", f"{rate_dn:.3f} / {fail_dn}")
diag("8_PC-B", "prediction 0.70-0.90 each direction", f"enrich {'met' if 0.70 <= rate_up <= 0.90 else ('above' if rate_up > 0.90 else 'NOT met')}; "
     f"deplete {'met' if 0.70 <= rate_dn <= 0.90 else ('above' if rate_dn > 0.90 else 'NOT met')}")
RES["PC_B"] = dict(rate_enrich=rate_up, rate_deplete=rate_dn, failed_enrich=fail_up, failed_deplete=fail_dn,
                   pass_=bool(rate_up >= 0.5 and rate_dn >= 0.5))

# ---------------------------------------------------------------- 8. PC-C
pi = pd.read_parquet(CAT + "passport_program_index.parquet")
pm = pd.read_parquet(CAT + "passport_program_membership.parquet")
pm["g"] = pm.ensembl_id.astype(str).str.split(".").str[0]
fib = pi[(pi.cell_type == "fibroblasts") & pi.program_label.str.contains("ECM|collagen", regex=True)]
fset = set(pm.g[pm.program_uid.isin(fib.program_uid)])
prog_f = full.g.isin(fset).to_numpy()
diag("8_PC-C", "programs / union genes / in U_full", f"{len(fib)} / {len(fset)} / {int(prog_f.sum())}")
pcc = strat_test(bins_f10, prog_f, deg_f)
pcc_pass = bool(pcc["O_over_E"] > 1.5 and pcc["p_exact"] < 0.01)
diag("8_PC-C", "O / E_strat / O/E / exact p", f"{pcc['O']} / {pcc['E_strat']:.3f} / {pcc['O_over_E']:.3f} / {pcc['p_exact']:.3g}")
diag("8_PC-C", "pass (O/E > 1.5 and p < 0.01)", pcc_pass)
RES["PC_C"] = {**pcc, "pass_": pcc_pass, "labels": fib.program_label.tolist()}

# ---------------------------------------------------------------- 8. NC-1 (S6)
def nc1(df, bins, rng, label):
    gen = df.gen.to_numpy(); deg = df.deg.to_numpy(); ae = df.AveExpr.to_numpy()
    gi = np.where(gen)[0]
    pool = np.where(~gen)[0]
    order = np.argsort(ae[pool], kind="stable")  # df is Ensembl-sorted, so ties keep Ensembl order
    pool = pool[order]; pv = ae[pool]; M = len(pool)
    D = int(deg.sum()); Nn = len(df); G = len(gi)
    n_b, _, D_b, _ = bin_margins(bins, gen, deg)
    rows = []
    for _ in range(args.nnc):
        # disjoint-set "next unused" pointers to the right and left, with path compression
        right = np.arange(M + 1); left = np.arange(-1, M)  # left[i+1] holds pointer for position i
        def fr(i):
            root = i
            while right[root] != root: root = right[root]
            while right[i] != root: right[i], i = root, right[i]
            return root
        def fl(i):  # i in [-1, M-1]; returns largest unused <= i or -1
            root = i
            while root >= 0 and left[root + 1] != root: root = left[root + 1]
            j = i
            while j >= 0 and left[j + 1] != root: left[j + 1], j = root, left[j + 1]
            return root
        chosen = []; gaps = []
        for gidx in rng.permutation(gi):
            v = ae[gidx]
            pos = int(np.searchsorted(pv, v))
            r = fr(pos) if pos < M else M
            l = fl(pos - 1)
            dr = pv[r] - v if r < M else np.inf
            dl = v - pv[l] if l >= 0 else np.inf
            c = l if dl <= dr else r   # tie -> lower AveExpr side
            chosen.append(pool[c]); gaps.append(min(dl, dr))
            right[c] = c + 1; left[c + 1] = c - 1
        lst = np.zeros(Nn, dtype=bool); lst[np.array(chosen)] = True
        a1 = ae[lst]; a0 = ae[gi]
        smd = (a1.mean() - a0.mean()) / np.sqrt((a1.var(ddof=1) + a0.var(ddof=1)) / 2)
        O = int((lst & deg).sum())
        _, Gk, _, Ok = bin_margins(bins, lst, deg)
        E = float((Gk * D_b / n_b).sum())
        pmf = exact_null(n_b, Gk, D_b)
        rows.append(dict(max_gap=float(max(gaps)), smd=float(smd), O=O, unstrat_OE=O / (G * D / Nn),
                         strat_OE=O / E, p_exact=exact_p(pmf, E, O)))
    t = pd.DataFrame(rows); t.insert(0, "universe", label)
    return t

nc_e = nc1(ueval, bins_e, RNG["NC1_eval"], "U_eval")
nc_f = nc1(full, bins_f10, RNG["NC1_full"], "U_full")
pd.concat([nc_e, nc_f]).to_csv(f"{OUT}/nc1_lists.tsv", sep="\t", index=False)
NC = {}
for lab, t in [("U_eval", nc_e), ("U_full", nc_f)]:
    NC[lab] = dict(mean_strat_OE=float(t.strat_OE.mean()), mean_unstrat_OE=float(t.unstrat_OE.mean()),
                   fpr=float((t.p_exact < 0.05).mean()), max_gap_median=float(t.max_gap.median()),
                   max_gap_max=float(t.max_gap.max()), smd_mean=float(t.smd.mean()), smd_maxabs=float(t.smd.abs().max()))
    diag("8_NC-1", lab, json.dumps({k: round(v, 4) for k, v in NC[lab].items()}))
pe = NC["U_eval"]; pf = NC["U_full"]
diag("8_NC-1", "predictions U_eval", f"strat O/E in 0.95-1.05: {0.95 <= pe['mean_strat_OE'] <= 1.05}; fpr in 0.02-0.08: "
     f"{0.02 <= pe['fpr'] <= 0.08}; unstrat O/E < 1: {pe['mean_unstrat_OE'] < 1}")
diag("8_NC-1", "predictions U_full", f"strat O/E in 0.95-1.05: {0.95 <= pf['mean_strat_OE'] <= 1.05}; unstrat O/E (predicted ~0.66): {pf['mean_unstrat_OE']:.3f}")
S6 = bool(pe["fpr"] > 0.10 or not (0.90 <= pe["mean_strat_OE"] <= 1.10))
diag("8_NC-1", "S6 triggered", S6)
RES["NC_1"] = NC; RES["S6"] = S6

# ---------------------------------------------------------------- 9. T1 permutation null mean (S7 part)
sims = perm_null(bins_e, gen_e, deg_e, args.nperm, RNG["T1_perm"])
diag("9_T1perm", "permutation mean - analytic E_strat", f"{sims.mean() - E_STRAT_NULLSPACE:.4f}")

# ---------------------------------------------------------------- 10. D3 formula check, D10 cluster counts
ce = pd.read_csv(f"{OUT}/control_expression.tsv", sep="\t")
chk = bulk[["gene", "AveExpr"]].merge(ce, on="gene", how="inner")
d3_diff = float(np.abs(chk.ave_all_units - chk.AveExpr).max())
d3_ok = bool(len(chk) == 24115 and d3_diff < 1e-8)
diag("10_D3", "genes matched / max |all-unit mean - AveExpr|", f"{len(chk)} / {d3_diff:.3g}")
diag("10_D3", "formula check pass (< 1e-8)", d3_ok)
diag("10_D3", "Spearman control-only vs AveExpr", f"{stats.spearmanr(chk.ave_control, chk.AveExpr)[0]:.4f}")

tier = pd.read_csv(TIER, sep="\t")
t12 = set(tier.study_name[tier.tier.astype(str).isin(["1", "2"])])
allg = pd.read_csv(A, usecols=["gwas_name", "ensembl", "chr", "PP.H4.susie", "method", "top_snp", "eligibility_source"],
                   dtype={"chr": str, "top_snp": str})
allg["g"] = allg.ensembl.astype(str).str.split(".").str[0]
genes_gen = ueval[ueval.gen][["g", "susie", "susie_driving_gwas"]]
pos = genes_gen.merge(allg, left_on=["g", "susie_driving_gwas"], right_on=["g", "gwas_name"], how="left").drop_duplicates("g")
diag("10_D10", "genetic genes / with top_snp", f"{len(genes_gen)} / {int(pos.top_snp.notna().sum())}")
sp = pos.top_snp.str.split(":", expand=True); pos["c"] = sp[0]; pos["p"] = pd.to_numeric(sp[1])


def clusters(w):
    lab = pd.Series(index=pos.g.values, dtype=object)
    for ch, grp in pos.sort_values(["c", "p", "g"]).groupby("c"):
        cid = np.concatenate([[0], np.cumsum(np.diff(grp.p.to_numpy()) > w)])
        lab[grp.g.values] = [f"{ch}_{i}" for i in cid]
    return lab

CL = {w: clusters(w) for w in (250_000, 500_000, 1_000_000)}
for w, lab in CL.items():
    diag("10_D10", f"clusters among {len(pos)} genetic genes at {w // 1000} kb", lab.nunique())
RES["D3_formula_ok"] = d3_ok

write_diag()
log("Diagnostics written; entering outcome section.")

# ======================================================================== OUTCOME
if args.debug_null:
    rng_dbg = np.random.default_rng(999)
    newdeg = deg_f.copy()
    for k in np.unique(bins_f10):
        ii = np.where(bins_f10 == k)[0]
        newdeg[ii] = rng_dbg.permutation(deg_f[ii])
    full["deg"] = newdeg
    ueval = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)
    ueval7 = full[full.susie.notna()].reset_index(drop=True)
    gen_e = ueval.gen.to_numpy(); deg_e = ueval.deg.to_numpy(); deg_f = full.deg.to_numpy()
    n_e, G_e, D_e, _ = bin_margins(bins_e, gen_e, deg_e)
    E_STRAT_NULLSPACE = float((G_e * D_e / n_e).sum())
    sims = perm_null(bins_e, gen_e, deg_e, args.nperm, np.random.default_rng(998))
    log("DEBUG NULL: DEG labels permuted within U_full AveExpr deciles. Outcome numbers below are not real.")

if not RES["PC_B"]["pass_"]:
    untestable("S5", f"PC-B detection rate enrich {rate_up:.3f}, deplete {rate_dn:.3f} (need >= 0.50 each)")

# ---- primary T1
T1 = strat_test(bins_e, gen_e, deg_e)
p_perm = perm_p(sims, T1["E_strat"], T1["O"])
band = 3 * np.sqrt(T1["p_exact"] * (1 - T1["p_exact"]) / args.nperm) + 1 / (args.nperm + 1)
perm_mean_diff = float(sims.mean() - T1["E_strat"])
s7 = bool(abs(p_perm - T1["p_exact"]) > band or abs(perm_mean_diff) >= 0.1)
T1.update(p_perm=p_perm, perm_band=band, perm_mean_minus_E=perm_mean_diff, S7=s7, K=K)
RES["T1"] = T1
log("T1", json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in T1.items()}))
if s7:
    stop("S7", f"permutation p {p_perm:.5f} vs exact {T1['p_exact']:.5f} (band {band:.5f}); perm mean - E {perm_mean_diff:.4f}")
rbg_excl = bool(T1["MH_lo"] > 1 or T1["MH_hi"] < 1)
t1_sig = T1["p_exact"] < 0.05

# ---- descriptives
D = {}
D["D1"] = unstrat(gen_e, deg_e)
D["D1"]["exp_frac_genetic_row"] = D["D1"]["D"] / D["D1"]["N"]; D["D1"]["exp_frac_deg_row"] = D["D1"]["G"] / D["D1"]["N"]
D2 = {}
for k in (5, 20):
    D2[f"AveExpr_{k}bins"] = strat_test(make_bins(ueval, "AveExpr", k), gen_e, deg_e)
D2["SE_10bins"] = strat_test(make_bins(ueval, "SE", 10), gen_e, deg_e)
joint = make_bins(ueval, "AveExpr", 5) * 5 + make_bins(ueval, "SE", 5)
D2["AveExpr5_x_SE5"] = strat_test(joint, gen_e, deg_e)
lg = ueval.assign(gen_i=gen_e.astype(int), deg_i=deg_e.astype(int))
m = smf.logit("deg_i ~ gen_i + cr(AveExpr, df=5)", lg).fit(disp=0, maxiter=500)
ci = m.conf_int().loc["gen_i"]
D2["logit_spline5"] = dict(OR=float(np.exp(m.params["gen_i"])), lo=float(np.exp(ci[0])), hi=float(np.exp(ci[1])), p=float(m.pvalues["gen_i"]), converged=bool(m.mle_retvals["converged"]))
D["D2"] = D2
# D3
ue3 = ueval.merge(ce[["gene", "ave_control"]], on="gene", how="left")
assert ue3.g.tolist() == ueval.g.tolist() and ue3.ave_control.notna().all()
b3 = make_bins(ue3, "ave_control", 10)
D3 = strat_test(b3, gen_e, deg_e)
s3 = perm_null(b3, gen_e, deg_e, args.nperm, RNG["D3_perm"])
D3.update(p_perm=perm_p(s3, D3["E_strat"], D3["O"]), perm_mean_minus_E=float(s3.mean() - D3["E_strat"]), formula_ok=d3_ok)
D["D3"] = D3
# D4
b7 = make_bins(ueval7, "AveExpr", K)
D["D4"] = strat_test(b7, ueval7.gen.to_numpy(), ueval7.deg.to_numpy())
# D5
D["D5"] = {f"{k}bins": strat_test(make_bins(full, "AveExpr", k), gen_f, deg_f) for k in (5, 10, 20, 50)}
D["D5"]["unstratified"] = unstrat(gen_f, deg_f)
# D6
gi_ = pd.read_parquet(CAT + "passport_gene_index.parquet")
gi_["g"] = gi_.ensembl_id.astype(str).str.split(".").str[0]
cls = gi_[gi_.primary_evidence_class.isin(["established_state_associated", "concordant", "genetically_anchored"])]
d6 = cls.merge(full[["g", "susie_state_t12", "deg", "gen"]], on="g", how="left")
D["D6"] = {c: d6[d6.primary_evidence_class == c].susie_state_t12.fillna("not_in_U_full").value_counts().to_dict()
           for c in ["established_state_associated", "concordant", "genetically_anchored"]}
est = d6[d6.primary_evidence_class == "established_state_associated"]
D["D6"]["established_equals_deg_only_set"] = bool(set(est.g) == set(full.g[full.deg & ~full.gen]))
# D7
nores = full[(full.susie_state_t12 == "no_susie_result") & full.deg]
a12 = allg[allg.gwas_name.isin(t12) & allg.g.isin(nores.g)]
meth = a12.groupby("g").method.apply(lambda s: "+".join(sorted(set(s))))
elig = a12.groupby("g").eligibility_source.apply(lambda s: "+".join(sorted(set(s.astype(str)))))
tested_ng = ueval[ueval.deg & ~ueval.gen]
q = [0.1, 0.25, 0.5, 0.75, 0.9]
D["D7"] = dict(n_genes=len(nores), n_with_t12_rows=int(meth.shape[0]),
               method_combos=meth.value_counts().to_dict(), eligibility_combos=elig.value_counts().to_dict(),
               AveExpr_no_result=nores.AveExpr.quantile(q).round(3).to_dict(), AveExpr_no_result_mean=float(nores.AveExpr.mean()),
               AveExpr_tested_nongenetic_DEG=tested_ng.AveExpr.quantile(q).round(3).to_dict(),
               AveExpr_tested_nongenetic_DEG_mean=float(tested_ng.AveExpr.mean()), n_tested_nongenetic_DEG=len(tested_ng))
# D8
gabf = (full.abf >= 0.5).fillna(False).to_numpy()
D["D8"] = {"unstratified": unstrat(gabf, deg_f), "AveExpr_10bins": strat_test(bins_f10, gabf, deg_f)}
# D9
O_ = T1["O"]; Gt = int(gen_e.sum()); Dt = int(deg_e.sum())
D["D9"] = dict(genetic_row=Gt, genetic_both=O_, genetic_only=Gt - O_, deg_row=Dt, deg_both=O_, deg_only=Dt - O_,
               strat_marker_genetic_row=T1["E_strat"] / Gt, strat_marker_deg_row=T1["E_strat"] / Dt,
               unstrat_marker_genetic_row=Dt / len(ueval), unstrat_marker_deg_row=Gt / len(ueval),
               never_evaluated_degs=int(((full.susie_state_t12 == "no_susie_result") & full.deg).sum()),
               neither_U_eval=int((~gen_e & ~deg_e).sum()))
# D10
shared = set(ueval.g[gen_e & deg_e])
D10 = {}
for w, lab in CL.items():
    D10[f"clusters_{w // 1000}kb"] = int(lab.nunique())
    D10[f"clusters_with_shared_{w // 1000}kb"] = int(lab[[g for g in lab.index if g in shared]].nunique())
lab = CL[500_000]
pp = pos.assign(cl=lab[pos.g].values).sort_values(["cl", "susie", "g"], ascending=[True, False, True])
keep = set(pp.drop_duplicates("cl").g)
drop = set(pos.g) - keep
u10 = ueval[~ueval.g.isin(drop)].reset_index(drop=True)
D10["U_eval_after_pruning_N_G_D_O"] = list(sizes(u10))
D10["test_500kb"] = strat_test(make_bins(u10, "AveExpr", 10), u10.gen.to_numpy(), u10.deg.to_numpy())
D["D10"] = D10
RES["descriptives"] = D

# ---------------------------------------------------------------- verdict (prespec "Reading")
side = np.sign(T1["O"] - T1["E_strat"])
d3_same = bool(np.sign(D3["O"] - D3["E_strat"]) == side and side != 0)
d10_same = bool(np.sign(D10["test_500kb"]["O"] - D10["test_500kb"]["E_strat"]) == side and side != 0)
conds = {"PC-A pass": RES["PC_A"]["pass"], "PC-B pass": RES["PC_B"]["pass_"], "S6 not triggered": not S6,
         "PC-C pass": pcc_pass, "D3 formula check pass and same side": d3_ok and d3_same,
         "D10 (500 kb) same side": d10_same}
RES["reading_conditions"] = conds
mhs = f"MH OR {T1['MH_OR']:.3f} (RBG 95% CI {T1['MH_lo']:.3f} to {T1['MH_hi']:.3f})"
core = (f"O {T1['O']}, E_strat {T1['E_strat']:.3f} ({K} AveExpr bins), O/E_strat {T1['O_over_E']:.3f}, {mhs}, "
        f"exact two-sided T1 p {T1['p_exact']:.4g} (family size 1, alpha 0.05, no correction)")
if t1_sig and all(conds.values()):
    direction = "less often" if side < 0 else "more often"
    verdict = (f"Directional: among genes evaluated by multi-signal COLOC, colocalized genes are {direction} canonical DEGs "
               f"than genes of the same expression level. {core}.")
elif t1_sig:
    failing = [k for k, v in conds.items() if not v]
    verdict = f"T1 p < 0.05 but no directional wording allowed; failing checks: {', '.join(failing)}. {core}."
else:
    verdict = (f"Null: {core}. ORs between about 0.46 and 1.77 were not detectable at 80% power (G3); "
               f"PC-B stratified detection rates {rate_up:.2f} (OR up) and {rate_dn:.2f} (OR down).")
if t1_sig != rbg_excl:
    verdict += " T1 and the RBG CI disagree on excluding 1; T1 governs."
if S6:
    verdict += " Expression not fully controlled (S6)."

fmt = lambda r: (f"O {r['O']}, E_strat {r['E_strat']:.3f}, O/E {r['O_over_E']:.3f}, MH OR {r.get('MH_OR', np.nan):.3f} "
                 f"({r.get('MH_lo', np.nan):.3f} to {r.get('MH_hi', np.nan):.3f}), exact p {r['p_exact']:.4g}")
L_ = ["## Primary (T1, U_eval, 3,769 genes)", "",
      f"- {core}.",
      f"- Permutation check ({args.nperm} within-bin permutations, stream 0): p {p_perm:.4g}, band {band:.4g}, "
      f"perm mean - E_strat {perm_mean_diff:.4f}. S7 not triggered.",
      f"- Not tested, reported beside T1: CMH chi-square {T1['CMH_chi2']:.3f} (nominal p {T1['CMH_p']:.4g}); "
      f"Breslow-Day {T1['BD_chi2']:.3f} (nominal p {T1['BD_p']:.4g}).",
      f"- RBG CI excludes 1: {rbg_excl}.", "",
      "## Controls", "",
      f"- PC-A reproduction: pass ({len(checks)} checks).",
      f"- PC-B detection rates: enrichment s=+14 {rate_up:.3f}, depletion s=-12 {rate_dn:.3f} (predicted 0.70 to 0.90; S5 needs >= 0.50).",
      f"- PC-C fibroblast ECM/collagen programs in U_full ({int(prog_f.sum())} genes): {fmt(pcc)}; pass {pcc_pass} (predicted O/E > 1.5, p < 0.01).",
      f"- NC-1 U_eval: mean stratified O/E {pe['mean_strat_OE']:.3f}, false-positive rate {pe['fpr']:.3f}, mean unstratified O/E {pe['mean_unstrat_OE']:.3f}, "
      f"median max gap {pe['max_gap_median']:.4f}, max |SMD| {pe['smd_maxabs']:.4f}. S6 triggered: {S6}.",
      f"- NC-1 U_full: mean stratified O/E {pf['mean_strat_OE']:.3f}, false-positive rate {pf['fpr']:.3f}, mean unstratified O/E {pf['mean_unstrat_OE']:.3f}.",
      f"- Identity check {id_ok}; role-swap check {rs_ok}.", "",
      "## Reading conditions for directional wording", ""] + [f"- {k}: {v}" for k, v in conds.items()] + [
      "", "## Descriptives (nominal p values, outside the family, never called significant)", "",
      f"- D1 U_eval unstratified: O {D['D1']['O']}, expected {D['D1']['E']:.3f}, Fisher OR {D['D1']['fisher_OR']:.4f}, p {D['D1']['fisher_p']:.4f}; "
      f"Figure 1B unstratified markers {100 * D['D1']['exp_frac_genetic_row']:.2f}% (genetic row), {100 * D['D1']['exp_frac_deg_row']:.2f}% (DEG row).",
      f"- D2 AveExpr 5 bins: {fmt(D2['AveExpr_5bins'])}.", f"- D2 AveExpr 20 bins: {fmt(D2['AveExpr_20bins'])}.",
      f"- D2 SE deciles: {fmt(D2['SE_10bins'])}.", f"- D2 AveExpr 5 x SE 5: {fmt(D2['AveExpr5_x_SE5'])}.",
      f"- D2 logistic DEG ~ genetic + cr(AveExpr, df=5): OR {D2['logit_spline5']['OR']:.3f} ({D2['logit_spline5']['lo']:.3f} to {D2['logit_spline5']['hi']:.3f}), Wald p {D2['logit_spline5']['p']:.4g}, converged {D2['logit_spline5']['converged']}.",
      f"- D3 control-only stratifier (138 control units): {fmt(D3)}; permutation p {D3['p_perm']:.4g}; formula check max diff {d3_diff:.3g}.",
      f"- D4 U_eval7 (3,776): {fmt(D['D4'])}.",
      "- D5 U_full: " + "; ".join(f"{k} {fmt(v)}" for k, v in D["D5"].items() if k != "unstratified") +
      f"; unstratified expected {D['D5']['unstratified']['E']:.3f}.",
      f"- D6 Catalog r8 classes by susie_state_t12: {json.dumps(D['D6'])}.",
      f"- D7 {D['D7']['n_genes']} never-evaluated DEGs: method combos {json.dumps(D['D7']['method_combos'])}; eligibility {json.dumps(D['D7']['eligibility_combos'])}; "
      f"mean AveExpr {D['D7']['AveExpr_no_result_mean']:.3f} against {D['D7']['AveExpr_tested_nongenetic_DEG_mean']:.3f} for the {D['D7']['n_tested_nongenetic_DEG']} tested non-genetic DEGs.",
      f"- D8 ABF map U_full: unstratified O {D['D8']['unstratified']['O']} / E {D['D8']['unstratified']['E']:.3f}; 10 bins {fmt(D['D8']['AveExpr_10bins'])}.",
      f"- D9 Figure 1B display: {json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in D['D9'].items()})}.",
      f"- D10 clusters: {json.dumps({k: v for k, v in D10.items() if k.startswith('clusters')})}; one gene per 500 kb cluster "
      f"(U_eval N/G/D/O {D10['U_eval_after_pruning_N_G_D_O']}): {fmt(D10['test_500kb'])}.",
      "", "Files: diagnostics.tsv, results.json, bins_aveexpr10.tsv, pcb_replicates.tsv, nc1_lists.tsv, control_expression.tsv."]
finish("completed", verdict, L_)
