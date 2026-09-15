#!/usr/bin/env python3
"""De novo two-axis test in PXD051911 liver protein.

Per protein, the partial Spearman with NAS adjusting for fibrosis, and with
fibrosis adjusting for NAS; batch and sex are covariates in both directions.
Counts are tested against a batch-stratified permutation null of the axis under
test. Verifies the held-back prespecification digest before opening any matrix.

No RNA gene list is read, mapped or compared. The question is whether the
STRUCTURE recurs, not whether the molecules match.
"""
import json, sys, hashlib, pathlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as student_t

PRESPEC_SHA256 = "6e7ac857e0403bd183d86cfd51a665b6a7b5d93bbfcbd2d62d1a3578179e60ed"
BENCH = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark")
PRESPEC = BENCH / "config/evaluation/pxd051911_two_axis_protein_prespecification.json"
SRC = BENCH / "executions/pxd051911-activation-readiness-21109461/source"
MIN_OBSERVED = 50
N_PERM = int(__import__("os").environ.get("PXD_PERM", "1000"))
SEED = 20260829


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    v = p[ok]
    n = v.size
    if n == 0:
        return q
    o = np.argsort(v)
    ranked = v[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def residualize(Y, C):
    """Residualize columns of Y on design C (with intercept already included)."""
    beta, *_ = np.linalg.lstsq(C, Y, rcond=None)
    return Y - C @ beta


def partial_spearman(Rprot, axis_rank, C, n_obs, k_cov):
    """Rprot: n x p already rank-transformed and residualized on C.
    axis_rank: n-vector of ranks. Returns (rho, p) per column."""
    a = residualize(axis_rank.reshape(-1, 1), C).ravel()
    sa = a.std()
    if sa == 0:
        return np.full(Rprot.shape[1], np.nan), np.full(Rprot.shape[1], np.nan)
    a = (a - a.mean()) / sa
    sp = Rprot.std(axis=0)
    good = sp > 0
    rho = np.full(Rprot.shape[1], np.nan)
    Z = (Rprot[:, good] - Rprot[:, good].mean(axis=0)) / sp[good]
    rho[good] = (Z * a[:, None]).sum(axis=0) / len(a)
    rho = np.clip(rho, -0.999999, 0.999999)
    df = n_obs - 2 - k_cov
    tstat = rho * np.sqrt(df / (1 - rho ** 2))
    pv = 2 * student_t.sf(np.abs(tstat), df)
    return rho, pv


def main(out_dir):
    out = pathlib.Path(out_dir)
    if out.exists():
        sys.exit(f"Refusing to overwrite: {out}")
    out.mkdir(parents=True)

    got = sha256(PRESPEC)
    if got != PRESPEC_SHA256:
        sys.exit(f"PRESPEC DIGEST MISMATCH\n  expected {PRESPEC_SHA256}\n  got      {got}")
    spec = json.loads(PRESPEC.read_text())
    print(f"prespec digest verified: {got}")
    print(f"workstream: {spec['workstream_id']}  release_state: {spec['release_state']}")

    meta = pd.read_csv(SRC / "meta_data.txt", sep="\t", dtype=str)
    liv = meta[meta.liver_proteomics_filename.notna() & meta.liver_proteomics_filename.ne("NA")].copy()
    assert len(liv) == 58 and liv.patient_name.nunique() == 58, "liver arm is not 58 participants at 1:1"
    liv["batch"] = np.where(liv.liver_proteomics_filename.str.contains("ATLASLiver"), 1.0, 0.0)
    liv["fibrosis"] = liv.kleiner_fibrosis_grade.str.extract(r"F(\d)").astype(float)
    liv["NAS"] = pd.to_numeric(liv.nafld_activity_score, errors="coerce")
    liv["sex"] = (liv.gender == "Male").astype(float)
    assert liv.fibrosis.notna().all() and liv.NAS.notna().all()
    assert set(liv.fibrosis.unique()) <= {0.0, 1.0, 2.0, 3.0}, "unexpected fibrosis level in the liver arm"

    q = pd.read_csv(SRC / "liver_protein_quant.txt", sep="\t", low_memory=False)
    sample_cols = [c for c in q.columns if c in set(liv.liver_proteomics_filename)]
    assert len(sample_cols) == 58, f"matrix join is {len(sample_cols)}, expected 58"
    liv = liv.set_index("liver_proteomics_filename").loc[sample_cols].reset_index()
    X = q[sample_cols].to_numpy(dtype=float).T          # 58 x proteins
    gene = q["Genes"].astype(str).to_numpy()
    acc = q["ProteinAccessions"].astype(str).to_numpy()
    print(f"matrix {X.shape[0]} samples x {X.shape[1]} proteins; batch sizes "
          f"Qexn2={int((liv.batch==0).sum())} ATLASLiver={int((liv.batch==1).sum())}")

    obs = np.isfinite(X).sum(axis=0)
    testable = obs >= MIN_OBSERVED
    print(f"testable proteins (>= {MIN_OBSERVED}/58 observed): {int(testable.sum())} of {X.shape[1]}")

    rows, null_rows = [], []
    rng = np.random.default_rng(SEED)
    directions = [("activity_NAS_adj_fibrosis", "NAS", "fibrosis"),
                  ("fibrosis_adj_NAS", "fibrosis", "NAS")]

    for scope, mask_rows in [("primary_all_58", np.ones(58, bool)),
                             ("sensitivity_qexn2_only", (liv.batch.values == 0))]:
        sub = liv.loc[mask_rows]
        Xs = X[mask_rows]
        n = int(mask_rows.sum())
        use_batch = scope == "primary_all_58"
        for dname, axis, other in directions:
            # complete-case per protein within this scope
            keep = np.isfinite(Xs).sum(axis=0) >= min(MIN_OBSERVED, n - 5)
            idx = np.where(keep)[0]
            # build covariate design
            cols = [np.ones(n), sub[other].to_numpy(float), sub["sex"].to_numpy(float)]
            if use_batch:
                cols.append(sub["batch"].to_numpy(float))
            C = np.column_stack(cols)
            k_cov = C.shape[1] - 1
            # rank-transform proteins with complete data only; drop any with NaN in scope
            Xk = Xs[:, idx]
            full = np.isfinite(Xk).all(axis=0)
            idx = idx[full]
            Xk = Xk[:, full]
            if Xk.shape[1] == 0:
                continue
            Rprot = np.apply_along_axis(rankdata, 0, Xk)
            Rprot = residualize(Rprot, C)
            axis_rank = rankdata(sub[axis].to_numpy(float))
            rho, pv = partial_spearman(Rprot, axis_rank, C, n, k_cov)
            qv = bh(pv)
            for j, col in enumerate(idx):
                rows.append((scope, dname, gene[col], acc[col], n, k_cov,
                             rho[j], pv[j], qv[j]))
            n_disc = int(np.nansum(qv < 0.05))
            print(f"  {scope:24s} {dname:26s} family={len(idx):5d}  BH q<0.05 = {n_disc}")

            # batch-stratified permutation null of the discovery count
            strata = sub["batch"].to_numpy(float) if use_batch else np.zeros(n)
            groups = [np.where(strata == s)[0] for s in np.unique(strata)]
            avals = sub[axis].to_numpy(float)
            null_counts = np.empty(N_PERM, int)
            for b in range(N_PERM):
                perm = avals.copy()
                for g in groups:
                    perm[g] = rng.permutation(perm[g])
                r_, p_ = partial_spearman(Rprot, rankdata(perm), C, n, k_cov)
                null_counts[b] = int(np.nansum(bh(p_) < 0.05))
            pval = (1 + int((null_counts >= n_disc).sum())) / (1 + N_PERM)
            null_rows.append((scope, dname, len(idx), n, n_disc,
                              float(null_counts.mean()), float(np.percentile(null_counts, 95)),
                              int(null_counts.max()), pval, N_PERM))
            print(f"      null mean {null_counts.mean():.1f}  p95 {np.percentile(null_counts,95):.0f}"
                  f"  max {null_counts.max()}  permutation p = {pval:.4g}")

    pd.DataFrame(rows, columns=["scope", "direction", "gene", "protein_accession",
                                "n", "k_covariates", "partial_spearman", "p_value", "bh_q"]
                 ).to_csv(out / "two_axis_protein_per_protein.tsv.gz", sep="\t", index=False)
    nd = pd.DataFrame(null_rows, columns=["scope", "direction", "family_size", "n",
                                          "observed_bh_discoveries", "null_mean", "null_p95",
                                          "null_max", "permutation_p", "replicates"])
    nd.to_csv(out / "two_axis_protein_counts_vs_null.tsv", sep="\t", index=False)

    prim = nd[nd.scope == "primary_all_58"]
    act = prim[prim.direction == "activity_NAS_adj_fibrosis"].iloc[0]
    fib = prim[prim.direction == "fibrosis_adj_NAS"].iloc[0]
    a_sig, f_sig = act.permutation_p < 0.05, fib.permutation_p < 0.05
    if a_sig and f_sig:
        cell = "TWO_AXIS_REPLICATES_IN_PROTEIN"
    elif a_sig or f_sig:
        cell = "SINGLE_AXIS_ONLY:" + ("activity" if a_sig else "fibrosis")
    else:
        cell = "INDETERMINATE"
    (out / "OUTCOME_CELL.txt").write_text(
        cell + "\n"
        "A null here is INDETERMINATE, never a negative: at MDE partial rho ~ 0.48 against\n"
        "Track A effects estimated at n=180, it cannot distinguish structure absent in\n"
        "protein from structure present but underpowered.\n")
    print("\n=== OUTCOME CELL:", cell, "===")
    print(nd.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: run_pxd051911_two_axis_protein.py OUTPUT_DIR")
    main(sys.argv[1])
