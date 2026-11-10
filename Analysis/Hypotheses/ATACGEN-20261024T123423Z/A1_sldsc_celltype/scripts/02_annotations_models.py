#!/usr/bin/env python3
"""A1 step 2: annotation catalogue, seeded Q2 permutation nulls, model table.

Reads peaks/peaks.tsv (step 1). Writes
  peaks/annotations.tsv   one row per annotation column: name, kind (binary or
                          continuous), set, lineage, draw, n_lifted, testable
  peaks/nulls_<L>.npz     for L in hepatocyte, stellate: tested-peak row index
                          into peaks.tsv, permuted z (K x n), null RESP values,
                          null tertile label (-1 lost, 0 stable, 1 gain)
  ldsc/models.tsv         model, kind, lineage, draw, columns

Null rule (PRESPEC): K = 50 per lineage; within each accessibility x GC stratum
z_comb is permuted among that stratum's tested peaks; one generator
numpy.random.default_rng(20261024); lineages hepatocyte then stellate; draws
1..50. Null RESP is the standardized |z| of the permuted values (same multiset,
so the same mean and SD as the observed RESP); null gain/lost are tertiles of the
permuted z within stratum, ties broken by peak order as in step 1.

Usage: 02_annotations_models.py <OUT>
"""
import os
import sys

import numpy as np
import pandas as pd

SEED = 20261024
K = 50
MIN_PEAKS = 1000
LINEAGES = ["cholangiocyte", "hepatocyte", "macrophage", "stellate", "t_nk"]
DA_LINEAGES = ["hepatocyte", "stellate"]


def tertile_labels(z, strata):
    lab = np.zeros(len(z), dtype=np.int8)
    for s in np.unique(strata):
        i = np.where(strata == s)[0]
        n = len(i)
        r = np.empty(n, int)
        r[np.argsort(z[i], kind="stable")] = np.arange(1, n + 1)  # ties: first in peak order
        k = n // 3
        lab[i[r <= k]] = -1
        lab[i[r > n - k]] = 1
    return lab


def main(out):
    pk = pd.read_csv(os.path.join(out, "peaks/peaks.tsv"), sep="\t")
    for c in ("lifted", "only", "prom", "tested", "gain", "lost", "gain10", "lost10"):
        pk[c] = pk[c].astype(str).str.upper().eq("TRUE")
    lifted = pk["lifted"].values
    rows = []

    def add(name, kind, setname, lineage, mask, draw=""):
        n = int((mask & lifted).sum())
        rows.append(dict(annotation=name, kind=kind, set=setname, lineage=lineage, draw=draw,
                         n_peaks=int(mask.sum()), n_lifted=n, testable=n >= MIN_PEAKS))

    add("UNION", "binary", "UNION", "", np.ones(len(pk), bool))
    for L in LINEAGES:
        isL = pk["lineage"].values == L
        add(f"ALL_{L}", "binary", "ALL", L, isL)
        add(f"ONLY_{L}", "binary", "ONLY", L, isL & pk["only"].values)
        add(f"PROM_{L}", "binary", "PROM", L, isL & pk["only"].values & pk["prom"].values)
        add(f"DIST_{L}", "binary", "DIST", L, isL & pk["only"].values & ~pk["prom"].values)
    for L in DA_LINEAGES:
        t = (pk["lineage"].values == L) & pk["tested"].values
        add(f"TESTED_{L}", "binary", "TESTED", L, t)
        for nm in ("ACCESS", "GC", "RESP", "SIGN"):
            add(f"{nm}_{L}", "continuous", nm, L, t)
        for nm in ("GAIN", "LOST", "GAIN10", "LOST10"):
            add(f"{nm}_{L}", "binary", nm, L, t & pk[nm.lower()].values)

    rng = np.random.default_rng(SEED)
    for L in DA_LINEAGES:
        idx = np.where((pk["lineage"].values == L) & pk["tested"].values)[0]
        z = pk["z_comb"].values[idx]
        strata = pk["stratum"].values[idx].astype(int)
        absz = np.abs(z)
        mu, sd = absz.mean(), absz.std(ddof=1)
        # observed sets must be reproduced by the same tertile rule
        lab_obs = tertile_labels(z, strata)
        assert np.array_equal(lab_obs == 1, pk["gain"].values[idx]), f"{L}: gain set not reproduced"
        assert np.array_equal(lab_obs == -1, pk["lost"].values[idx]), f"{L}: lost set not reproduced"
        zp = np.empty((K, len(idx)), np.float64)
        for k in range(K):
            zk = z.copy()
            for s in np.unique(strata):
                i = np.where(strata == s)[0]
                zk[i] = z[i][rng.permutation(len(i))]
            zp[k] = zk
        resp = ((np.abs(zp) - mu) / sd).astype(np.float32)
        lab = np.stack([tertile_labels(zp[k], strata) for k in range(K)])
        np.savez(os.path.join(out, f"peaks/nulls_{L}.npz"), idx=idx, z=zp.astype(np.float32),
                 resp=resp, lab=lab)
        t = (pk["lineage"].values == L) & pk["tested"].values
        for k in range(1, K + 1):
            add(f"NULLRESP_{L}_{k}", "continuous", "NULLRESP", L, t, draw=str(k))
            g = np.zeros(len(pk), bool)
            g[idx[lab[k - 1] == 1]] = True
            add(f"NULLGAIN_{L}_{k}", "binary", "NULLGAIN", L, g, draw=str(k))
            g = np.zeros(len(pk), bool)
            g[idx[lab[k - 1] == -1]] = True
            add(f"NULLLOST_{L}_{k}", "binary", "NULLLOST", L, g, draw=str(k))
        print(f"{L}: {len(idx)} tested peaks, {len(np.unique(strata))} strata, {K} null draws")

    ann = pd.DataFrame(rows)
    ann.to_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t", index=False)
    ok = set(ann.loc[ann["testable"], "annotation"])

    mrows = []

    def model(name, kind, L, cols, draw=""):
        if all(c in ok for c in cols):
            mrows.append((name, kind, L, draw, ",".join(cols)))
        else:
            print(f"model {name} not built: untestable column(s) {[c for c in cols if c not in ok]}")

    for L in LINEAGES:
        model(f"ONLY_{L}", "only", L, ["UNION", f"ONLY_{L}"])
        model(f"ALL_{L}", "all", L, ["UNION", f"ALL_{L}"])
        model(f"PD_{L}", "pd", L, ["UNION", f"PROM_{L}", f"DIST_{L}"])
    for L in DA_LINEAGES:
        cov = [f"ALL_{L}", f"TESTED_{L}", f"ACCESS_{L}", f"GC_{L}"]
        model(f"CONT_{L}", "cont", L, cov + [f"RESP_{L}"])
        model(f"SIGN_{L}", "sign", L, cov + [f"SIGN_{L}"])
        model(f"TERT_{L}", "tert", L, cov + [f"GAIN_{L}", f"LOST_{L}"])
        model(f"TOP10_{L}", "top10", L, cov + [f"GAIN10_{L}", f"LOST10_{L}"])
        for k in range(1, K + 1):
            model(f"NULLCONT_{L}_{k}", "nullcont", L, cov + [f"NULLRESP_{L}_{k}"], str(k))
            model(f"NULLTERT_{L}_{k}", "nulltert", L, cov + [f"NULLGAIN_{L}_{k}", f"NULLLOST_{L}_{k}"], str(k))
    md = pd.DataFrame(mrows, columns=["model", "kind", "lineage", "draw", "columns"])
    os.makedirs(os.path.join(out, "ldsc"), exist_ok=True)
    md.to_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", index=False)
    print(ann[ann["draw"].astype(str) == ""][["annotation", "n_peaks", "n_lifted", "testable"]].to_string(index=False))
    print(md["kind"].value_counts().to_string())


if __name__ == "__main__":
    main(sys.argv[1])
