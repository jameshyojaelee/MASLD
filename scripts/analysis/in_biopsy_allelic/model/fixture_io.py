"""Read the B-FIX shared fixture (model-a-bfix-*) into a Model A Problem, parameter vectors
and draw files. Checks the row order, the oriented count, s_t and the design columns."""
import csv
import json
from pathlib import Path

import numpy as np

from het_rule import het
from model_a import Problem, Rows

BASE_Z = ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10", "z_c"]
PARAM_FILES = ["theta1_start_like", "theta2_truth_v0", "theta3_truth_vpos", "theta4_perturbed"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _num(s):
    return np.nan if s in ("NA", "", "nan", "NaN") else float(s)


def load_problem(fdir):
    fdir = Path(fdir)
    glob = json.loads((fdir / "globals.json").read_text())
    cohorts = {r["cohort"]: r for r in read_tsv(fdir / "cohorts.tsv")}
    genes = {r["gene_id"]: r for r in read_tsv(fdir / "genes.tsv")}
    gene_names = sorted(glob["gene_order"])
    assert gene_names == glob["gene_order"], "gene order is not gene_id ascending"
    ind = read_tsv(fdir / "individuals.tsv")
    ids = [r["individual_id"] for r in ind]
    pos = {i: k for k, i in enumerate(ids)}
    ind_cohort = np.array([r["cohort"] for r in ind])
    zraw = np.array([[_num(r[c]) for c in BASE_Z] for r in ind])
    rows = read_tsv(fdir / "rows.tsv")
    keys = [(r["individual_id"], r["gene_id"]) for r in rows]
    assert keys == sorted(keys) and len(set(keys)) == len(keys), "rows not ordered by (individual, gene)"
    col = lambda name, f=int: np.array([f(r[name]) for r in rows])  # noqa: E731
    n, ref, alt, a = col("n"), col("ref_reads"), col("alt_reads"), col("a")
    orient, s_t = col("orientation"), col("s_t")
    assert np.all(n == ref + alt)
    assert np.all(a == np.where(orient == 1, alt, ref)), "oriented count"
    assert np.all(s_t == -orient), "s_t"
    row_cohort = [r["cohort"] for r in rows]
    row_ind = np.array([pos[r["individual_id"]] for r in rows])
    assert all(ind_cohort[i] == c for i, c in zip(row_ind, row_cohort))
    rule = np.array([[int(cohorts[c][k]) for k in ("min_dp", "min_minor_reads", "min_minor_pct")]
                     for c in row_cohort])
    assert het(ref, alt, rule[:, 0], rule[:, 1], rule[:, 2]).all(), "a row is not a het call"
    e_row = col("e_i", float)
    e_ind = np.array([float(r["e_i"]) for r in ind])
    assert np.allclose(e_row, e_ind[row_ind], rtol=0, atol=0)
    R = Rows(n, a, rule[:, 0], rule[:, 1], rule[:, 2], orient, col("p_ref", float),
             col("p_alt", float), col("p_x", float), e_row, glob["phi"])
    gidx = {g: k for k, g in enumerate(gene_names)}
    prob = Problem(ids, ind_cohort, np.array([float(r["S"]) for r in ind]), zraw, BASE_Z,
                   np.array([r["cluster_id"] for r in ind]), R, row_ind,
                   np.array([gidx[r["gene_id"]] for r in rows]), gene_names,
                   np.array([float(cohorts[c]["omega_c"]) for c in row_cohort]),
                   w_fib=[float(genes[g]["w_fib"]) for g in gene_names],
                   w_mye=[float(genes[g]["w_mye"]) for g in gene_names])
    assert prob.design_columns == glob["design_columns"], (prob.design_columns, glob["design_columns"])
    return prob, glob


def load_params(fdir, L):
    fdir = Path(fdir)
    return {name: L.from_dict(json.loads((fdir / "params" / f"{name}.json").read_text()))
            for name in PARAM_FILES}


def load_draws(fdir, prob):
    fdir = Path(fdir) / "draws"
    s_rows = read_tsv(fdir / "p1_Sstar.tsv")
    Sstar = np.array([[float(r[i]) for i in prob.ind_ids] for r in s_rows])
    assert [int(r["draw"]) for r in s_rows] == list(range(len(s_rows)))
    res = read_tsv(fdir / "p2_resample.tsv")
    B2 = max(int(r["draw"]) for r in res) + 1
    cl_codes = {c: k for k, c in enumerate(np.unique(prob.cluster))}
    ind_cl = np.array([cl_codes[c] for c in prob.cluster])
    mult = np.zeros((B2, len(prob.ind_ids)))
    strata_ok = True
    first_member = {}
    for i in np.argsort(prob.ind_ids):
        first_member.setdefault(prob.cluster[i], prob.cohort[i])
    for r in res:
        b, c = int(r["draw"]), r["cluster_id"]
        strata_ok &= first_member[c] == r["stratum"]
        counts = np.zeros(len(cl_codes))
        counts[cl_codes[c]] = 1
        mult[b] += counts[ind_cl]
    lam_rows = read_tsv(fdir / "lambda_draws.tsv")
    sources = sorted({r["source"] for r in lam_rows})
    lambdas = {s: np.full(max(int(r["draw"]) for r in lam_rows) + 1, np.nan) for s in sources}
    for r in lam_rows:
        lambdas[r["source"]][int(r["draw"])] = _num(r["lambda"])
    U = np.load(fdir / "p3_uniforms.npy")
    return dict(Sstar=Sstar, mult=mult, strata_match_rule=bool(strata_ok), lambdas=lambdas, U=U)
