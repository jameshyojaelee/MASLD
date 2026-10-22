"""Independent refit of the 117-program primary stage model (score ~ stage_ordinal + dataset).

Checks claim R04-01: primary q for Ductular injury (hepatocytes 20), HC3 q,
GSE189600 drop, single-donor deletions, and a seeded within-dataset permutation null.
"""
import sys
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SCORES = ROOT + "/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
META = ROOT + "/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
PRIMARY = ROOT + "/figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv"
STAGE = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
TARGET = ("hepatocytes", 20)
N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
SEED = 20261015

from scipy import stats


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def load():
    s = pd.read_csv(SCORES, sep="\t", keep_default_na=False)
    m = pd.read_csv(META, sep="\t", keep_default_na=False)
    m = m[["sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"]].rename(
        columns={"sample": "donor", "dataset": "dataset_meta", "disease_stage_coarse": "stage_corr",
                 "exclude_stage_analysis": "excl_corr"})
    d = s.merge(m, on="donor", how="left", validate="many_to_one")
    assert d.stage_corr.notna().all()
    assert (d.dataset == d.dataset_meta).all()
    excl = d.excl_corr.astype(str).str.upper().eq("TRUE")
    d = d[d.stage_corr.isin(STAGE) & ~excl].copy()
    d["stage"] = d.stage_corr.map(STAGE).astype(float)
    return d[["cell_type", "module", "donor", "dataset", "stage", "score"]]


def fit_one(y, x, g):
    """OLS of y on x + dataset fixed effects via within-dataset centering (FWL)."""
    df = pd.DataFrame({"y": y, "x": x, "g": g})
    ng = df.groupby("g")["y"].transform("size").to_numpy()
    yc = (df.y - df.groupby("g").y.transform("mean")).to_numpy()
    xc = (df.x - df.groupby("g").x.transform("mean")).to_numpy()
    n = len(y)
    k = df.g.nunique()
    p = k + 1
    sxx = (xc ** 2).sum()
    beta = (xc * yc).sum() / sxx
    e = yc - beta * xc
    rdf = n - p
    se = np.sqrt((e ** 2).sum() / rdf / sxx)
    h = 1.0 / ng + xc ** 2 / sxx
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(h < 1 - 1e-10, e ** 2 / (1 - h) ** 2, 0.0)
    se3 = np.sqrt((xc ** 2 * w).sum()) / sxx
    t = beta / se
    t3 = beta / se3
    return dict(beta=beta, se=se, t=t, p=2 * stats.t.sf(abs(t), rdf), hc3_se=se3, t3=t3,
                p3=2 * stats.t.sf(abs(t3), rdf), rdf=rdf, n=n, k=k, maxlev=h.max())


def family(d):
    rows = []
    for (ct, mod), sub in d.groupby(["cell_type", "module"], sort=True):
        r = fit_one(sub.score.to_numpy(), sub.stage.to_numpy(), sub.dataset.to_numpy())
        r.update(cell_type=ct, module=mod)
        rows.append(r)
    f = pd.DataFrame(rows)
    f["q"] = bh(f.p)
    f["q3"] = bh(f.p3)
    return f


def target(f):
    return f[(f.cell_type == TARGET[0]) & (f.module == TARGET[1])].iloc[0]


def main():
    d = load()
    f = family(d)
    ref = pd.read_csv(PRIMARY, sep="\t")
    cmp = f.merge(ref, on=["cell_type", "module"], suffixes=("", "_ref"))
    print("programs fitted", len(f), "matched to primary_effects", len(cmp))
    print("max|beta-ref|", np.nanmax(np.abs(cmp.beta - cmp.beta_ref)),
          "max|q-ref|", np.nanmax(np.abs(cmp.q - cmp.qvalue)),
          "max|q3-ref|", np.nanmax(np.abs(cmp.q3 - cmp.hc3_qvalue)),
          "max|n-ref|", np.nanmax(np.abs(cmp.n - cmp.n_donors)),
          "max|maxlev-ref|", np.nanmax(np.abs(cmp.maxlev - cmp.max_leverage)))
    t = target(f)
    print("FULL 62: beta=%.6f se=%.6f p=%.6g q=%.6f hc3_p=%.6g hc3_q=%.6f n=%d k=%d" %
          (t.beta, t.se, t.p, t.q, t.p3, t.q3, t.n, t.k))
    print("n programs q<0.05:", int((f.q < 0.05).sum()), " hc3 q<0.05:", int((f.q3 < 0.05).sum()))

    td = d[(d.cell_type == TARGET[0]) & (d.module == TARGET[1])]
    print("target donors by dataset x stage:\n", pd.crosstab(td.dataset, td.stage))
    print(td[td.dataset == "GSE189600"][["donor", "dataset", "stage", "score"]].to_string())

    # (2) drop GSE189600
    f2 = family(d[d.dataset != "GSE189600"])
    t2 = target(f2)
    print("DROP GSE189600: beta=%.6f p=%.6g q=%.6f hc3_q=%.6f n=%d; programs q<0.05=%d" %
          (t2.beta, t2.p, t2.q, t2.q3, t2.n, int((f2.q < 0.05).sum())))

    # (3) single-donor deletions within the target program's donors (refit full 117 family)
    out = []
    for dn in td.donor.unique():
        fi = family(d[d.donor != dn])
        ti = target(fi)
        out.append(dict(donor=dn, dataset=td[td.donor == dn].dataset.iloc[0],
                        stage=td[td.donor == dn].stage.iloc[0], beta=ti.beta, p=ti.p, q=ti.q,
                        q3=ti.q3, n_sel=int((fi.q < 0.05).sum())))
    lo = pd.DataFrame(out).sort_values("q", ascending=False)
    print("single-donor deletions:", len(lo), " target q>=0.05 in", int((lo.q >= 0.05).sum()))
    print(lo.head(16).to_string(index=False))
    lo.to_csv("loo_donor_deletions.tsv", sep="\t", index=False)
    f.to_csv("refit_full_family.tsv", sep="\t", index=False)
    f2.to_csv("refit_drop_gse189600_family.tsv", sep="\t", index=False)

    if N_PERM > 0:
        permute(d, f, N_PERM)


def permute(d, f, n_perm):
    rng = np.random.default_rng(SEED)
    donors = np.sort(d.donor.unique())
    # per-permutation random key per donor; each program permutes its own donors within dataset by key order
    keys = rng.random((n_perm, len(donors)))
    key_idx = {dn: i for i, dn in enumerate(donors)}
    progs = list(d.groupby(["cell_type", "module"], sort=True))
    absT = np.zeros((n_perm, len(progs)))
    absT3 = np.zeros((n_perm, len(progs)))
    obs = np.zeros(len(progs))
    obs3 = np.zeros(len(progs))
    tgt = None
    for j, ((ct, mod), sub) in enumerate(progs):
        if (ct, mod) == TARGET:
            tgt = j
        sub = sub.sort_values(["dataset", "donor"]).reset_index(drop=True)
        g = sub.dataset.to_numpy()
        y = sub.score.to_numpy()
        x = sub.stage.to_numpy()
        n = len(y)
        ugs = np.unique(g)
        ng = pd.Series(g).map(pd.Series(g).value_counts()).to_numpy()
        yc = y - pd.Series(y).groupby(g).transform("mean").to_numpy()
        xc0 = x - pd.Series(x).groupby(g).transform("mean").to_numpy()
        sxx = (xc0 ** 2).sum()
        rdf = n - (len(ugs) + 1)
        # build permuted centered x: (n_perm, n)
        XC = np.empty((n_perm, n))
        cols = np.array([key_idx[dn] for dn in sub.donor])
        K = keys[:, cols]
        for gg in ugs:
            idx = np.where(g == gg)[0]
            order = np.argsort(K[:, idx], axis=1)  # rank order of donors in this dataset
            # assign original labels (in fixed order) to donors in key order
            perm_x = np.empty((n_perm, len(idx)))
            rows = np.arange(n_perm)[:, None]
            perm_x[rows, order] = xc0[idx][None, :]
            XC[:, idx] = perm_x
        allx = np.vstack([xc0[None, :], XC])
        beta = (allx * yc[None, :]).sum(1) / sxx
        e = yc[None, :] - beta[:, None] * allx
        se = np.sqrt((e ** 2).sum(1) / rdf / sxx)
        h = 1.0 / ng[None, :] + allx ** 2 / sxx
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(h < 1 - 1e-10, e ** 2 / (1 - h) ** 2, 0.0)
        se3 = np.sqrt((allx ** 2 * w).sum(1)) / sxx
        t = np.abs(beta / se)
        t3 = np.abs(beta / se3)
        obs[j], obs3[j] = t[0], t3[0]
        absT[:, j], absT3[:, j] = t[1:], t3[1:]
    pp = ((absT >= obs[None, :] - 1e-12).sum(0) + 1) / (n_perm + 1)
    pp3 = ((absT3 >= obs3[None, :] - 1e-12).sum(0) + 1) / (n_perm + 1)
    qp = bh(pp)
    qp3 = bh(pp3)
    mx = absT.max(1)
    mx3 = absT3.max(1)
    fw = ((mx >= obs[tgt] - 1e-12).sum() + 1) / (n_perm + 1)
    fw3 = ((mx3 >= obs3[tgt] - 1e-12).sum() + 1) / (n_perm + 1)
    print("PERM n=%d seed=%d floor=%.2g" % (n_perm, SEED, 1 / (n_perm + 1)))
    print("target |t| obs OLS=%.4f HC3=%.4f" % (obs[tgt], obs3[tgt]))
    print("target perm p OLS=%.6g BH117 q=%.6g maxT FWER=%.6g" % (pp[tgt], qp[tgt], fw))
    print("target perm p HC3=%.6g BH117 q=%.6g maxT FWER=%.6g" % (pp3[tgt], qp3[tgt], fw3))
    print("programs perm-BH q<0.05: OLS=%d HC3=%d" % ((qp < 0.05).sum(), (qp3 < 0.05).sum()))
    pd.DataFrame({"prog": [f"{c}_{m}" for (c, m), _ in progs], "obs_t": obs, "obs_t3": obs3,
                  "perm_p": pp, "perm_q": qp, "perm_p3": pp3, "perm_q3": qp3}).to_csv(
        "perm_null_family.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
