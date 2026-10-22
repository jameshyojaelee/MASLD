#!/usr/bin/env python3
"""NG07 (L01-03) feasibility gate. BLINDED: reports GGT h2 gate statistics,
per-trait tau* SEs, cross-trait jackknife correlation and SE of
Delta = tau*_GGT - tau*_ALT for SPEC_cholangiocyte and SPEC_hepatocyte.
It does NOT write GGT tau* point estimates or Delta point estimates.
ALT tau z (positive control, already known from H22) is reported.
Logic copied from H22 05_run_h2.py (pairing) and 06_summarize.py (tau*, jk).
Usage: gate.py <G> <H22> <alt_munged> <ggt_munged>
"""
import gzip, json, os, re, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

G, H22, ALT, GGT = sys.argv[1:5]
L = "/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin"
REF = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/ldsc_ref/sldsc"
BASE = f"{REF}/baselineLD_v2.2/baselineLD."
W = f"{REF}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FRQ = f"{REF}/1000G_Phase3_frq/1000G.EUR.QC."
MODELS = ["CT_cholangiocyte", "CT_hepatocyte"]


def jk_cov(a, b):
    B = len(a)
    return (B - 1) / B * np.sum((a - a.mean()) * (b - b.mean()))


def pair():
    ref = pd.concat([pd.read_csv(f"{BASE}{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)])
    w = pd.concat([pd.read_csv(f"{W}{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)])
    keep = set(ref["SNP"]) & set(w["SNP"])
    ss = {}
    for t, p in (("alt", ALT), ("ggt", GGT)):
        x = pd.read_csv(p, sep="\t").dropna(subset=["Z", "N"])
        ss[t] = x[x["SNP"].isin(keep)].drop_duplicates("SNP").set_index("SNP")
    snps = ss["alt"].index.intersection(ss["ggt"].index)
    rep = {"n_common_before_chisq": int(len(snps))}
    for it in range(20):
        ok = pd.Series(True, index=snps)
        for t in ss:
            x = ss[t].loc[snps]
            thr = max(0.001 * x["N"].max(), 80.0)
            ok &= ~((x["Z"] ** 2) >= thr)
            rep[f"{t}_chisq_max"] = float(thr)
        if ok.all():
            break
        snps = snps[ok.values]
    rep["n_pair_snps"] = int(len(snps))
    out = {}
    for t in ss:
        out[t] = os.path.join(G, "ldsc", f"{t}_pair.sumstats.gz")
        with gzip.open(out[t], "wt") as fh:
            ss[t].loc[snps].reset_index()[["SNP", "A1", "A2", "Z", "N"]].to_csv(fh, sep="\t", index=False)
    json.dump(rep, open(os.path.join(G, "ldsc", "pair_report.json"), "w"), indent=2)
    return out


def run(args):
    ss, model, out = args
    pre = [BASE] + ([] if model == "BASE" else [os.path.join(H22, "ldsc/models", model, f"{model}.")])
    cmd = [f"{L}/python", f"{L}/ldsc.py", "--h2", ss, "--ref-ld-chr", ",".join(pre), "--w-ld-chr", W,
           "--print-delete-vals", "--out", out, "--overlap-annot", "--frqfile-chr", FRQ, "--print-coefficients"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    assert r.returncode == 0 and "Analysis finished" in open(out + ".log").read(), out
    return out


def load(trait, model, base_m550, base_names, sd, m550, M_tot):
    pre = os.path.join(G, "ldsc/h2", trait, model)
    cols = ["UNION", "SPEC_" + model[3:]]
    names = base_names + cols
    M = np.concatenate([base_m550, [m550[c] for c in cols]])
    res = pd.read_csv(pre + ".results", sep="\t")
    coef = res["Coefficient"].to_numpy(float)
    se = res["Coefficient_std_error"].to_numpy(float)
    dele = np.loadtxt(pre + ".part_delete", ndmin=2)
    assert dele.shape[1] == len(names) == len(coef)
    h2 = float(M @ coef)
    h2_del = dele @ M
    i = names.index(cols[1])
    k = sd[cols[1]] * M_tot
    log = open(pre + ".log").read()
    n_snp = int(re.findall(r"(\d+) SNPs remain", log)[-1])
    return dict(tau_z=coef[i] / se[i], ts=coef[i] * k / h2, ts_del=dele[:, i] * k / h2_del, n_snp=n_snp, B=dele.shape[0])


def main():
    ps = pair()
    jobs = [(ps["ggt"], "BASE", os.path.join(G, "ldsc/h2/ggt/BASE")),
            (ps["alt"], "BASE", os.path.join(G, "ldsc/h2/alt/BASE"))]
    for m in MODELS:
        for t in ("alt", "ggt"):
            jobs.append((ps[t], m, os.path.join(G, "ldsc/h2", t, m)))
    with ProcessPoolExecutor(4) as ex:
        list(ex.map(run, jobs))
    bm = [[float(x) for x in open(f"{BASE}{c}.l2.M_5_50").readline().split()] for c in range(1, 23)]
    base_m550 = np.sum(bm, axis=0)
    hdr = pd.read_csv(f"{BASE}1.l2.ldscore.gz", sep="\t", nrows=0).columns[3:]
    base_names = [h[:-2] for h in hdr]
    M_tot = float(base_m550[0])
    mm = pd.concat([pd.read_csv(os.path.join(H22, f"ldsc/packed/M.{c}.tsv"), sep="\t") for c in range(1, 23)])
    m550 = mm.groupby("annotation", sort=False)["M_5_50"].sum()
    st = pd.concat([pd.read_csv(os.path.join(H22, f"ldsc/annot_stats/{c}.tsv"), sep="\t") for c in range(1, 23)])
    g = st.groupby("annotation", sort=False)[["n_common_in", "n_common_total"]].sum()
    p = g["n_common_in"] / g["n_common_total"]
    sd = np.sqrt(p * (1 - p))
    out = {}
    for t in ("ggt", "alt"):
        log = open(os.path.join(G, "ldsc/h2", t, "BASE.log")).read()
        h = re.search(r"Total Observed scale h2: ([-\d.e]+) \(([\d.e-]+)\)", log)
        it = re.search(r"Intercept: ([\d.]+) \(([\d.]+)\)", log)
        ra = re.search(r"Ratio: ([-\d.]+) \(([\d.]+)\)", log)
        mc = re.search(r"Mean Chi\^2: ([\d.]+)", log)
        out[f"{t}_base"] = dict(h2=float(h.group(1)), h2_se=float(h.group(2)),
                               h2_z=float(h.group(1)) / float(h.group(2)), intercept=float(it.group(1)),
                               ratio=float(ra.group(1)) if ra else None, mean_chisq=float(mc.group(1)))
    for m in MODELS:
        a = load("alt", m, base_m550, base_names, sd, m550, M_tot)
        b = load("ggt", m, base_m550, base_names, sd, m550, M_tot)
        assert a["n_snp"] == b["n_snp"] and a["B"] == b["B"]
        va, vb, cab = jk_cov(a["ts_del"], a["ts_del"]), jk_cov(b["ts_del"], b["ts_del"]), jk_cov(a["ts_del"], b["ts_del"])
        se_d = float(np.sqrt(va + vb - 2 * cab))
        out[m] = dict(alt_tau_z=float(a["tau_z"]), alt_tau_star=float(a["ts"]), alt_tau_star_se=float(np.sqrt(va)),
                      ggt_tau_star_se=float(np.sqrt(vb)), jk_corr=float(cab / np.sqrt(va * vb)),
                      se_delta=se_d, se_delta_if_independent=float(np.sqrt(va + vb)), n_snp=a["n_snp"])
    json.dump(out, open(os.path.join(G, "gate_blinded.json"), "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
