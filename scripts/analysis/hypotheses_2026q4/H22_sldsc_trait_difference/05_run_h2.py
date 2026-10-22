#!/usr/bin/env python3
"""H22 step 5: partitioned-heritability regressions (native ldsc.py --h2).

Phase 1 (stop rules): baselineLD-only h2 for the four traits, and the cell-type
model for hepatocyte peaks with ALT (positive control) and height (negative
control). If PDFF or cirrhosis h2 z < 4, or the ALT hepatocyte coefficient
z < 1.96, write ldsc/STOP.json and run nothing else.

Phase 2: all other models. PDFF and cirrhosis use pair-restricted summary
statistics with an identical SNP set, so the 200 jackknife blocks coincide and
their delete values can be paired. Null models are materialised on the fly in
ldsc/tmp_null/<MODEL>/ from the packed float32 LD scores and their files are
removed (file by file) after both traits have run.

Usage: 05_run_h2.py <OUT> <CONFIG_TSV> <N_WORKERS>
CONFIG_TSV has key<TAB>value rows: ldsc_python, ldsc_py, baseline, weights,
frq, eur_w_ld, pdff, cirr, alt, height.
"""
import gzip
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

PAIR = ("pdff", "cirr")


def cfg_read(path):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#"))


def finished(log):
    return os.path.exists(log) and "Analysis finished" in open(log).read()


def run_ldsc(cfg, sumstats, prefixes, out, overlap=True):
    log = out + ".log"
    if finished(log):
        return out
    cmd = [cfg["ldsc_python"], cfg["ldsc_py"], "--h2", sumstats,
           "--ref-ld-chr", ",".join(prefixes), "--w-ld-chr", cfg["weights"],
           "--print-delete-vals", "--out", out]
    if overlap:
        cmd += ["--overlap-annot", "--frqfile-chr", cfg["frq"], "--print-coefficients"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0 or not finished(log):
        raise RuntimeError(f"ldsc failed: {out}\n{r.stdout.decode()[-2000:]}")
    return out


def pair_sumstats(out, cfg):
    """Identical regression SNP set for PDFF and cirrhosis."""
    d = os.path.join(out, "ldsc/sumstats")
    os.makedirs(d, exist_ok=True)
    paths = {t: os.path.join(d, f"{t}_pair.sumstats.gz") for t in PAIR}
    if all(os.path.exists(p) for p in paths.values()):
        return paths
    ref = pd.concat([pd.read_csv(f"{cfg['baseline']}{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"])
                     for c in range(1, 23)])
    w = pd.concat([pd.read_csv(f"{cfg['weights']}{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"])
                   for c in range(1, 23)])
    keep = set(ref["SNP"]) & set(w["SNP"])
    ss = {}
    for t in PAIR:
        x = pd.read_csv(cfg[t], sep="\t").dropna(subset=["Z", "N"])
        ss[t] = x[x["SNP"].isin(keep)].drop_duplicates("SNP").set_index("SNP")
    common = ss["pdff"].index.intersection(ss["cirr"].index)
    rep = {"n_common_before_chisq": int(len(common))}
    # LDSC recomputes chisq_max = max(0.001 * N_max, 80) on the SNPs it receives,
    # so iterate until removing SNPs no longer changes either threshold.
    snps = common
    for it in range(20):
        ok = pd.Series(True, index=snps)
        for t in PAIR:
            x = ss[t].loc[snps]
            thr = max(0.001 * x["N"].max(), 80.0)
            ok &= ~((x["Z"] ** 2) >= thr)
            rep[f"{t}_chisq_max"] = float(thr)
        if ok.all():
            break
        snps = snps[ok.values]
    rep["chisq_iterations"] = it + 1
    for t in PAIR:
        rep[f"{t}_n_removed_chisq_total"] = int(len(common) - len(snps))
    rep["n_pair_snps"] = int(len(snps))
    for t in PAIR:
        x = ss[t].loc[snps].reset_index()[["SNP", "A1", "A2", "Z", "N"]]
        with gzip.open(paths[t], "wt") as fh:
            x.to_csv(fh, sep="\t", index=False)
    json.dump(rep, open(os.path.join(d, "pair_report.json"), "w"), indent=2)
    return paths


def materialise_null(out, model, cols):
    d = os.path.join(out, "ldsc/tmp_null", model)
    os.makedirs(d, exist_ok=True)
    meta = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t")
    names = list(meta.loc[meta["testable"], "annotation"])  # row order of packed l2.<c>.npy
    ii = [names.index(x) for x in cols]
    written = []
    for c in range(1, 23):
        arr = np.asarray(np.load(os.path.join(out, f"ldsc/packed/l2.{c}.npy"), mmap_mode="r")[ii, :]).T
        snps = pd.read_csv(os.path.join(out, f"ldsc/packed/snps.{c}.tsv.gz"), sep="\t")
        df = snps.copy()
        for j, n in enumerate(cols):
            df[n + "L2"] = arr[:, j]
        pre = os.path.join(d, f"{model}.{c}")
        with gzip.open(pre + ".l2.ldscore.gz", "wt", compresslevel=1) as fh:
            df.to_csv(fh, sep="\t", index=False, float_format="%.3f")
        mm = pd.read_csv(os.path.join(out, f"ldsc/packed/M.{c}.tsv"), sep="\t").set_index("annotation")
        for suf, col in ((".l2.M", "M"), (".l2.M_5_50", "M_5_50")):
            with open(pre + suf, "w") as fh:
                fh.write("\t".join(repr(float(mm.loc[x, col])) for x in cols) + "\n")
        written += [pre + s for s in (".l2.ldscore.gz", ".l2.M", ".l2.M_5_50")]
    return os.path.join(d, f"{model}."), written


def task_real(out, cfg, model, trait, ss):
    res = os.path.join(out, "ldsc/h2", trait)
    os.makedirs(res, exist_ok=True)
    prefixes = [cfg["baseline"]]
    if model != "BASE":
        prefixes.append(os.path.join(out, "ldsc/models", model, f"{model}."))
    return run_ldsc(cfg, ss, prefixes, os.path.join(res, model))


def task_null(out, cfg, model, cols, ss_pair):
    res_paths = [os.path.join(out, "ldsc/h2", t, model) for t in PAIR]
    if all(finished(p + ".log") for p in res_paths):
        return model
    prefix, written = materialise_null(out, model, cols)
    try:
        for t, rp in zip(PAIR, res_paths):
            os.makedirs(os.path.dirname(rp), exist_ok=True)
            run_ldsc(cfg, ss_pair[t], [cfg["baseline"], prefix], rp, overlap=False)
    finally:
        for f in written:
            if os.path.exists(f):
                os.remove(f)
    return model


def tau_z(results_path, category_suffix):
    df = pd.read_csv(results_path, sep="\t")
    row = df[df["Category"].str.startswith(category_suffix + "L2")]
    return float(row["Coefficient_z-score"].iloc[0])


def h2_z(log):
    txt = open(log).read()
    m = re.search(r"Total Observed scale h2: ([-\d.eE]+) \(([-\d.eE]+)\)", txt)
    return float(m.group(1)) / float(m.group(2)), float(m.group(1)), float(m.group(2))


def main(out, cfg_path, n_workers):
    cfg = cfg_read(cfg_path)
    n_workers = int(n_workers)
    ss_pair = pair_sumstats(out, cfg)
    ss = {"pdff": ss_pair["pdff"], "cirr": ss_pair["cirr"], "alt": cfg["alt"], "height": cfg["height"]}
    models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", dtype=str, keep_default_na=False)

    # ---- phase 1 ----
    p1 = [("BASE", t) for t in ss] + [("CT_hepatocyte", "alt"), ("CT_hepatocyte", "height")]
    with ProcessPoolExecutor(n_workers) as ex:
        futs = [ex.submit(task_real, out, cfg, m, t, ss[t]) for m, t in p1]
        for f in as_completed(futs):
            print("done", f.result(), flush=True)
    stop = {}
    for t in ss:
        z, h2, se = h2_z(os.path.join(out, "ldsc/h2", t, "BASE.log"))
        stop[f"h2_z_{t}"] = z
    stop["alt_hepatocyte_tau_z"] = tau_z(os.path.join(out, "ldsc/h2/alt/CT_hepatocyte.results"), "SPEC_hepatocyte")
    stop["height_hepatocyte_tau_z"] = tau_z(os.path.join(out, "ldsc/h2/height/CT_hepatocyte.results"), "SPEC_hepatocyte")
    stop["stop_h2"] = bool(stop["h2_z_cirr"] < 4 or stop["h2_z_pdff"] < 4)
    stop["stop_positive_control"] = bool(stop["alt_hepatocyte_tau_z"] < 1.96)
    json.dump(stop, open(os.path.join(out, "ldsc/phase1_checks.json"), "w"), indent=2)
    print(json.dumps(stop, indent=2), flush=True)
    if stop["stop_h2"] or stop["stop_positive_control"]:
        json.dump(stop, open(os.path.join(out, "ldsc/STOP.json"), "w"), indent=2)
        print("STOP rule met; contrast models not run")
        return

    # ---- phase 2 ----
    jobs = []
    for _, r in models.iterrows():
        if r["kind"] == "null":
            continue
        traits = ["pdff", "cirr", "alt", "height"] if r["kind"] in ("ct", "stea") else ["pdff", "cirr"]
        for t in traits:
            if (r["model"], t) in p1:
                continue
            jobs.append(("real", r["model"], t))
    nulls = models[models["kind"] == "null"]
    with ProcessPoolExecutor(n_workers) as ex:
        futs = [ex.submit(task_real, out, cfg, m, t, ss[t]) for _, m, t in jobs]
        futs += [ex.submit(task_null, out, cfg, r["model"], r["columns"].split(","), ss_pair)
                 for _, r in nulls.iterrows()]
        n_done = 0
        for f in as_completed(futs):
            f.result()
            n_done += 1
            if n_done % 25 == 0:
                print(f"{n_done}/{len(futs)} tasks done", flush=True)

    # ---- overlap check: cross-trait LDSC ----
    rg_out = os.path.join(out, "ldsc/rg/pdff_cirr")
    os.makedirs(os.path.dirname(rg_out), exist_ok=True)
    if not finished(rg_out + ".log"):
        subprocess.run([cfg["ldsc_python"], cfg["ldsc_py"], "--rg", f"{cfg['pdff']},{cfg['cirr']}",
                        "--ref-ld-chr", cfg["eur_w_ld"], "--w-ld-chr", cfg["eur_w_ld"], "--out", rg_out],
                       check=True, stdout=subprocess.DEVNULL)
    print("phase 2 done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
