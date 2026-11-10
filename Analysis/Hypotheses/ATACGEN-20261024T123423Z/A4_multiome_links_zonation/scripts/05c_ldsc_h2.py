#!/usr/bin/env python3
"""A4 step 5c: S-LDSC regressions for the zonation peak sets (native ldsc.py --h2).

Adapted copy of H22 05_run_h2.py (run_ldsc, materialise_null); H22 files are not edited.
  BASE: baselineLD only, every trait (stop rule S3: h2 z >= 4).
  MAIN: baselineLD + ALL_HEP + PC + PP, every trait (--overlap-annot, --print-coefficients).
  SEC:  baselineLD + ALL_HEP + PCq + PPq, if testable.
  NULL_d: baselineLD + ALL_HEP + NULL_PC_d + NULL_PP_d, d = 1..20, for S3-eligible liver
          traits and height (coefficients only; materialised from packed LD scores).
Usage: 05c_ldsc_h2.py <A4> <CONFIG_TSV> <N_WORKERS>
CONFIG_TSV: key<TAB>value rows ldsc_python, ldsc_py, baseline, weights, frq, and one row per
trait sumstats (alt, ast, ggt, pdff, cirr, nafld, height).
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

TRAITS = ["alt", "ast", "ggt", "pdff", "cirr", "nafld", "height"]
LIVER = TRAITS[:-1]
H2_Z_MIN = 4.0


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


def materialise(a4, model, cols):
    d = os.path.join(a4, "work/ldsc/tmp_null", model)
    os.makedirs(d, exist_ok=True)
    meta = pd.read_csv(os.path.join(a4, "work/ldsc/annotations.tsv"), sep="\t")
    names = list(meta.loc[meta["testable"], "annotation"])
    ii = [names.index(x) for x in cols]
    written = []
    for c in range(1, 23):
        arr = np.asarray(np.load(os.path.join(a4, f"work/ldsc/packed/l2.{c}.npy"), mmap_mode="r")[ii, :]).T
        df = pd.read_csv(os.path.join(a4, f"work/ldsc/packed/snps.{c}.tsv.gz"), sep="\t")
        for j, n in enumerate(cols):
            df[n + "L2"] = arr[:, j]
        pre = os.path.join(d, f"{model}.{c}")
        with gzip.open(pre + ".l2.ldscore.gz", "wt", compresslevel=1) as fh:
            df.to_csv(fh, sep="\t", index=False, float_format="%.3f")
        mm = pd.read_csv(os.path.join(a4, f"work/ldsc/packed/M.{c}.tsv"), sep="\t").set_index("annotation")
        for suf, col in ((".l2.M", "M"), (".l2.M_5_50", "M_5_50")):
            with open(pre + suf, "w") as fh:
                fh.write("\t".join(repr(float(mm.loc[x, col])) for x in cols) + "\n")
        written += [pre + s for s in (".l2.ldscore.gz", ".l2.M", ".l2.M_5_50")]
    return os.path.join(d, f"{model}."), written


def task_real(a4, cfg, model, trait):
    res = os.path.join(a4, "work/ldsc/h2", trait)
    os.makedirs(res, exist_ok=True)
    prefixes = [cfg["baseline"]]
    if model != "BASE":
        prefixes.append(os.path.join(a4, "work/ldsc/models", model, f"{model}."))
    return run_ldsc(cfg, cfg[trait], prefixes, os.path.join(res, model))


def task_null(a4, cfg, model, cols, traits):
    outs = [os.path.join(a4, "work/ldsc/h2", t, model) for t in traits]
    if all(finished(o + ".log") for o in outs):
        return model
    prefix, written = materialise(a4, model, cols)
    try:
        for t, o in zip(traits, outs):
            os.makedirs(os.path.dirname(o), exist_ok=True)
            run_ldsc(cfg, cfg[t], [cfg["baseline"], prefix], o, overlap=False)
    finally:
        for f in written:
            if os.path.exists(f):
                os.remove(f)
    return model


def h2_z(log):
    m = re.search(r"Total Observed scale h2: ([-\d.eE]+) \(([-\d.eE]+)\)", open(log).read())
    return float(m.group(1)) / float(m.group(2)), float(m.group(1)), float(m.group(2))


def main(a4, cfg_path, n_workers):
    cfg = cfg_read(cfg_path)
    n_workers = int(n_workers)
    meta = pd.read_csv(os.path.join(a4, "work/ldsc/annotations.tsv"), sep="\t")
    sec = bool(meta.loc[meta["annotation"] == "PCq", "testable"].iloc[0])
    with ProcessPoolExecutor(n_workers) as ex:
        for f in as_completed([ex.submit(task_real, a4, cfg, "BASE", t) for t in TRAITS]):
            print("done", f.result(), flush=True)
    chk = {}
    for t in TRAITS:
        z, h2, se = h2_z(os.path.join(a4, "work/ldsc/h2", t, "BASE.log"))
        chk[t] = {"h2": h2, "se": se, "z": z, "eligible": bool(z >= H2_Z_MIN) if t in LIVER else False}
    json.dump(chk, open(os.path.join(a4, "work/ldsc/s3_checks.json"), "w"), indent=2)
    print(json.dumps(chk, indent=2), flush=True)
    null_traits = [t for t in LIVER if chk[t]["eligible"]] + ["height"]
    draws = sorted(meta["draw"].dropna().astype(int).unique())
    with ProcessPoolExecutor(n_workers) as ex:
        futs = [ex.submit(task_real, a4, cfg, "MAIN", t) for t in TRAITS]
        if sec:
            futs += [ex.submit(task_real, a4, cfg, "SEC", t) for t in TRAITS]
        futs += [ex.submit(task_null, a4, cfg, f"NULL_{d}", ["ALL_HEP", f"NULL_PC_{d}", f"NULL_PP_{d}"], null_traits)
                 for d in draws]
        n_done = 0
        for f in as_completed(futs):
            f.result()
            n_done += 1
            print(f"{n_done}/{len(futs)} tasks done", flush=True)
    print("h2 runs done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
