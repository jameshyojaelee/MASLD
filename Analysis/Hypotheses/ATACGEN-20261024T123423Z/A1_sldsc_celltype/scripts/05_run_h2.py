#!/usr/bin/env python3
"""A1 step 5: partitioned-heritability regressions (native ldsc.py --h2), adapted from H22 05.

Phase 1: baselineLD-only h2 for every GWAS in the trait table. Testable = h2 z >= 4
(z = h2 / SE printed by LDSC). Untestable GWAS run nothing else.
Phase 2, testable GWAS only:
  every role       ONLY_<L>, ALL_<L> (v3, 5 lineages); CT296875_<L> (H22 files, 7 lineages)
  test, descriptive PD_<L>; CONT/SIGN/TERT/TOP10_<L> for hepatocyte and stellate
  test             null models NULLCONT/NULLTERT_<L>_<k>, materialised on the fly from
                   the packed float32 LD scores and deleted after all test GWAS ran.

Usage: 05_run_h2.py <OUT> <CONFIG_TSV> <TRAITS_TSV> <N_WORKERS>
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

H22_LINEAGES = ["b_cell", "cholangiocyte", "endothelial_cell", "fibroblast", "hepatocyte", "macrophage", "t_cell"]
NULL_KINDS = ("nullcont", "nulltert")


def cfg_read(path):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#"))


def finished(log):
    if not os.path.exists(log):
        return False
    txt = open(log).read()  # ldsc prints "Analysis finished" even after an error
    return "Analysis finished" in txt and "Traceback" not in txt


def run_ldsc(cfg, sumstats, prefixes, out, overlap=True):
    log = out + ".log"
    if finished(log):
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
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


def model_prefix(out, cfg, model):
    if model == "BASE":
        return [cfg["baseline"]]
    if model.startswith("CT296875_"):
        L = model[len("CT296875_"):]
        return [cfg["baseline"], os.path.join(cfg["h22"], "ldsc/models", f"CT_{L}", f"CT_{L}.")]
    return [cfg["baseline"], os.path.join(out, "ldsc/models", model, f"{model}.")]


def task_real(out, cfg, model, trait, ss):
    return run_ldsc(cfg, ss, model_prefix(out, cfg, model), os.path.join(out, "ldsc/h2", trait, model))


def materialise(out, model, cols):
    d = os.path.join(out, "ldsc/tmp_null", model)
    os.makedirs(d, exist_ok=True)
    ann = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t", keep_default_na=False)
    names = list(ann.loc[ann["testable"].astype(str).str.upper().eq("TRUE"), "annotation"])
    ii = [names.index(x) for x in cols]
    written = []
    for c in range(1, 23):
        arr = np.asarray(np.load(os.path.join(out, f"ldsc/packed/l2.{c}.npy"), mmap_mode="r")[ii, :]).T
        df = pd.read_csv(os.path.join(out, f"ldsc/packed/snps.{c}.tsv.gz"), sep="\t")
        for j, n in enumerate(cols):
            df[n + "L2"] = arr[:, j]
        pre = os.path.join(d, f"{model}.{c}")
        with gzip.open(pre + ".l2.ldscore.gz", "wt", compresslevel=1) as fh:
            df.to_csv(fh, sep="\t", index=False, float_format="%.4f")
        mm = pd.read_csv(os.path.join(out, f"ldsc/packed/M.{c}.tsv"), sep="\t").set_index("annotation")
        for suf, col in ((".l2.M", "M"), (".l2.M_5_50", "M_5_50")):
            with open(pre + suf, "w") as fh:
                fh.write("\t".join(repr(float(mm.loc[x, col])) for x in cols) + "\n")
        written += [pre + s for s in (".l2.ldscore.gz", ".l2.M", ".l2.M_5_50")]
    return os.path.join(d, f"{model}."), written


def task_null(out, cfg, model, cols, traits):
    res = [os.path.join(out, "ldsc/h2", t, model) for t, _ in traits]
    if all(finished(p + ".log") for p in res):
        return model
    prefix, written = materialise(out, model, cols)
    try:
        for (t, ss), rp in zip(traits, res):
            run_ldsc(cfg, ss, [cfg["baseline"], prefix], rp, overlap=False)
    finally:
        for f in written:
            if os.path.exists(f):
                os.remove(f)
    return model


def h2_z(log):
    txt = open(log).read()
    m = re.search(r"Total Observed scale h2: ([-\d.eE]+) \(([-\d.eE]+)\)", txt)
    return float(m.group(1)) / float(m.group(2)), float(m.group(1)), float(m.group(2))


def run_pool(n_workers, jobs):
    n_done = 0
    with ProcessPoolExecutor(n_workers) as ex:
        futs = [ex.submit(*j) for j in jobs]
        for f in as_completed(futs):
            f.result()
            n_done += 1
            if n_done % 50 == 0 or n_done == len(futs):
                print(f"{n_done}/{len(futs)} tasks done", flush=True)


def main(out, cfg_path, traits_path, n_workers):
    cfg = cfg_read(cfg_path)
    n_workers = int(n_workers)
    tr = pd.read_csv(traits_path, sep="\t")
    models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", dtype=str, keep_default_na=False)
    ss = dict(zip(tr["trait_id"], tr["sumstats"]))

    # ---- phase 1: testability ----
    run_pool(n_workers, [(task_real, out, cfg, "BASE", t, ss[t]) for t in tr["trait_id"]])
    rows = []
    for t in tr["trait_id"]:
        z, h2, se = h2_z(os.path.join(out, "ldsc/h2", t, "BASE.log"))
        rows.append(dict(trait_id=t, h2=h2, h2_se=se, h2_z=z, testable=z >= 4))
    base = pd.DataFrame(rows)
    base.to_csv(os.path.join(out, "ldsc/phase1_testability.tsv"), sep="\t", index=False)
    print(base.to_string(index=False), flush=True)
    tr = tr.merge(base[["trait_id", "testable"]], on="trait_id")
    ok = tr[tr["testable"]]

    # ---- phase 2 ----
    real = models[~models["kind"].isin(NULL_KINDS)]
    jobs = []
    for _, t in ok.iterrows():
        for _, m in real.iterrows():
            if m["kind"] in ("only", "all") or t["role"] in ("test", "descriptive"):
                jobs.append((task_real, out, cfg, m["model"], t["trait_id"], t["sumstats"]))
        for L in H22_LINEAGES:
            jobs.append((task_real, out, cfg, f"CT296875_{L}", t["trait_id"], t["sumstats"]))
    test_traits = [(r["trait_id"], r["sumstats"]) for _, r in ok[ok["role"] == "test"].iterrows()]
    for _, m in models[models["kind"].isin(NULL_KINDS)].iterrows():
        jobs.append((task_null, out, cfg, m["model"], m["columns"].split(","), test_traits))
    print(f"phase 2: {len(jobs)} tasks", flush=True)
    run_pool(n_workers, jobs)
    print("phase 2 done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:5])
