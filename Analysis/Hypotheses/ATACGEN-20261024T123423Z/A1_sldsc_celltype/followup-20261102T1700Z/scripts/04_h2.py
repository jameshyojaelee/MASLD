#!/usr/bin/env python3
"""Follow-up step 4: S-LDSC regressions (ldsc.py --h2 --overlap-annot), every trait x model,
plus the Q7 paired AST/ALT runs on the H22 CT_endothelial_cell model.
Usage: 04_h2.py <F> <CONFIG_TSV> <N_WORKERS>
CONFIG_TSV rows: ldsc_python, ldsc_py, baseline, weights, frq, h22, and one row per trait.
"""
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

TRAITS = ["ast_panukbb", "ast_mvp", "alt_panukbb", "alt_mvp", "ggt_panukbb", "pdff", "height_panukbb"]
MODELS = ["J7", "RM_E122", "RM_E100", "RM_E108", "RM_E066", "RMJ"]
PAIRS = ["pair_panukbb_ast", "pair_panukbb_alt", "pair_mvp_ast", "pair_mvp_alt"]


def cfg_read(path):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#"))


def finished(log):
    if not os.path.exists(log):
        return False
    txt = open(log).read()  # ldsc prints "Analysis finished" even after an error
    return "Analysis finished" in txt and "Traceback" not in txt


def run(cfg, sumstats, prefix, out):
    if finished(out + ".log"):
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [cfg["ldsc_python"], cfg["ldsc_py"], "--h2", sumstats, "--ref-ld-chr", f"{cfg['baseline']},{prefix}",
           "--w-ld-chr", cfg["weights"], "--overlap-annot", "--frqfile-chr", cfg["frq"],
           "--print-coefficients", "--print-delete-vals", "--out", out]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0 or not finished(out + ".log"):
        raise RuntimeError(f"ldsc failed: {out}\n{r.stdout.decode()[-2000:]}")
    return out


def main(f, cfg_path, n_workers):
    cfg = cfg_read(cfg_path)
    jobs = [(cfg[t], os.path.join(f, "ldsc/models", m, f"{m}."), os.path.join(f, "ldsc/h2", t, m))
            for t in TRAITS for m in MODELS]
    ct = os.path.join(cfg["h22"], "ldsc/models/CT_endothelial_cell/CT_endothelial_cell.")
    jobs += [(os.path.join(f, "sumstats", f"{p}.sumstats.gz"), ct, os.path.join(f, "ldsc/h2", p, "CT_endothelial_cell"))
             for p in PAIRS]
    with ProcessPoolExecutor(int(n_workers)) as ex:
        futs = [ex.submit(run, cfg, *j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            print(f"{i}/{len(futs)} {fu.result()}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
