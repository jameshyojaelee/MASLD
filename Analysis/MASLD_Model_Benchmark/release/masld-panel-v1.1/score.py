#!/usr/bin/env python3
"""masld-panel-v1.1 — 20-gene-pair liver fibrosis severity score (abundant protein-coding pairs, sample-independent NaN handling).

Reads EITHER RNA-seq counts OR qPCR Cq values for the panel genes and returns one number per sample.

WHAT CHANGED FROM v1 (N07):
  1. Abundance- and biotype-qualified panel fitted label-free on POOL_V2 (median CPM >= 5.0,
     frac_zero <= 0.05, protein_coding, non-chrY, excluding olfactory/testis families such as OR2I1).
  2. Direct Cq evaluation (is_cq=True): averages technical replicates per sample and evaluates
     Cq_b - Cq_a directly without 2^(-Cq) / nansum library-total coupling.
  3. Strictly sample-independent NaN / missing-pair handling:
     ok_ps = isfinite(lg[ra, s]) & isfinite(lg[rb, s]) is evaluated PER SAMPLE s (not across the whole
     batch via np.all). Missing pairs are centered at their frozen training pair_mean (m_p) via
     b0_raw + (P / used_s) * sum_{p in ok_s} coef_p * (diff_{p, s} - m_p).
     One patient's NaN or non-detect causes 0.0 drift on any co-batched patient.
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WEIGHTS = HERE / "weights/masld_panel_v1_1.npz"

MIN_PAIR_COVERAGE = 0.70


class ScoreError(RuntimeError):
    """Everything raises. Nothing warns and continues with a number."""


def load_panel(path=WEIGHTS):
    d = np.load(path, allow_pickle=True)
    coef = np.asarray(d["coef"], float)
    pair_mean = np.asarray(d["pair_mean"], float)
    intercept = float(d["intercept"])
    b0_raw = float(d["b0_raw"]) if "b0_raw" in d.files else float(intercept + np.dot(coef, pair_mean))
    return {
        "a": [str(x) for x in d["pair_gene_a"]],
        "b": [str(x) for x in d["pair_gene_b"]],
        "coef": coef,
        "intercept": intercept,
        "pair_mean": pair_mean,
        "b0_raw": b0_raw,
        "floor_cpm": float(d["floor_cpm"]),
        "formulation": str(d["formulation"]),
        "interval_calibration": str(d["interval_calibration"]),
    }


def normalise_ids(ids):
    out = []
    for g in ids:
        g = str(g).strip().upper()
        if g.endswith("_PAR_Y"):
            g = g[: -len("_PAR_Y")]
        if "." in g:
            head, _, tail = g.rpartition(".")
            if tail.isdigit():
                g = head
        out.append(g)
    return out


def read_table(path):
    """Rows = genes, columns = samples. First column is the gene id."""
    import csv

    p = Path(path)
    opener = __import__("gzip").open if p.suffix == ".gz" else open
    with opener(p, "rt") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        samples = header[1:]
        genes, rows = [], []
        for r in rd:
            if not r:
                continue
            if len(r) != len(header):
                raise ScoreError(
                    f"ragged row for gene {r[0]!r}: {len(r)} fields against a {len(header)}-field header."
                )
            genes.append(r[0])
            rows.append([float(x) for x in r[1:]])
    return normalise_ids(genes), samples, np.asarray(rows, float)


def score(values, gene_ids, panel, is_cq=False):
    """values: genes x samples. Returns (scores, report)."""
    values = np.asarray(values, dtype=float)
    gene_ids = normalise_ids(gene_ids)
    idx = {}
    for i, g in enumerate(gene_ids):
        idx.setdefault(g, []).append(i)
    needed = sorted(set(panel["a"]) | set(panel["b"]))
    present = [g for g in needed if g in idx]
    missing = [g for g in needed if g not in idx]
    need_pos = {g: i for i, g in enumerate(needed)}

    n = values.shape[1]
    if is_cq:
        # Direct Cq mode: mean replicates per gene per sample without 2^(-Cq) / nansum coupling
        cq_mat = np.full((len(needed), n), np.nan, dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for r, g in enumerate(needed):
                if g in idx:
                    cq_mat[r] = np.nanmean(values[idx[g]], axis=0)
        pair_diffs = np.full((len(panel["coef"]), n), np.nan, dtype=float)
        ok_ps = np.zeros((len(panel["coef"]), n), dtype=bool)
        for p_idx, (ga, gb) in enumerate(zip(panel["a"], panel["b"])):
            ra, rb = need_pos[ga], need_pos[gb]
            fin = np.isfinite(cq_mat[ra]) & np.isfinite(cq_mat[rb])
            ok_ps[p_idx] = fin
            pair_diffs[p_idx, fin] = cq_mat[rb, fin] - cq_mat[ra, fin]
    else:
        mat = np.full((len(needed), n), np.nan, dtype=float)
        for r, g in enumerate(needed):
            if g in idx:
                sub = values[idx[g]]
                has_fin = np.isfinite(sub).any(axis=0)
                if np.any(sub[np.isfinite(sub)] < 0):
                    raise ScoreError("negative values in a count matrix")
                mat[r, has_fin] = np.nansum(sub[:, has_fin], axis=0)

        # Sample-independent library total and CPM-of-panel floor
        tot = np.nansum(mat, axis=0)
        if np.any(tot <= 0):
            bad = np.flatnonzero(tot <= 0)[:5].tolist()
            raise ScoreError(f"panel total is zero for sample index {bad}: no panel gene detected")

        lg = np.full_like(mat, np.nan)
        fin_mat = np.isfinite(mat)
        for s in range(n):
            fs = fin_mat[:, s]
            comp_s = np.maximum(mat[fs, s] / tot[s] * 1e6, panel["floor_cpm"])
            lg[fs, s] = np.log2(comp_s)

        pair_diffs = np.full((len(panel["coef"]), n), np.nan, dtype=float)
        ok_ps = np.zeros((len(panel["coef"]), n), dtype=bool)
        for p_idx, (ga, gb) in enumerate(zip(panel["a"], panel["b"])):
            ra, rb = need_pos[ga], need_pos[gb]
            fin = np.isfinite(lg[ra]) & np.isfinite(lg[rb])
            ok_ps[p_idx] = fin
            pair_diffs[p_idx, fin] = lg[ra, fin] - lg[rb, fin]

    n_pairs = len(panel["coef"])
    coef = panel["coef"]
    pair_mean = panel["pair_mean"]
    b0_raw = panel["b0_raw"]
    scores = np.full(n, np.nan, dtype=float)
    used_per_sample = ok_ps.sum(axis=0)

    for s in range(n):
        used_s = int(used_per_sample[s])
        cov_s = used_s / n_pairs
        if cov_s < MIN_PAIR_COVERAGE:
            missing_s = [
                f"{panel['a'][p]}-{panel['b'][p]}" for p in range(n_pairs) if not ok_ps[p, s]
            ]
            raise ScoreError(
                f"sample index {s} pair coverage {cov_s:.3f} is below the {MIN_PAIR_COVERAGE} contract "
                f"({used_s} of {n_pairs} pairs evaluable). Missing pairs: {missing_s}. "
                f"REFUSING rather than returning a score from a partial panel."
            )
        mask_s = ok_ps[:, s]
        # Missing pairs are centered at pair_mean (m_p), so delta = 0; evaluable centered deltas
        # are scaled by (n_pairs / used_s) around b0_raw. When used_s == n_pairs this equals
        # intercept + sum_p coef_p * diff_{p,s} exactly.
        centered = pair_diffs[mask_s, s] - pair_mean[mask_s]
        scores[s] = b0_raw + float(np.dot(coef[mask_s], centered)) * (n_pairs / used_s)

    return scores, {
        "n_samples": int(n),
        "input_scale": "cq" if is_cq else "counts",
        "panel_genes_required": len(needed),
        "panel_genes_present": len(present),
        "panel_genes_missing": missing,
        "pairs_evaluable_min": int(used_per_sample.min()),
        "pairs_evaluable_max": int(used_per_sample.max()),
        "pairs_total": int(n_pairs),
        "pair_coverage_min": float(used_per_sample.min() / n_pairs),
        "pair_coverage_contract": MIN_PAIR_COVERAGE,
        "interval_calibration": panel["interval_calibration"],
        "score_is_an_ORDERING_not_a_stage": True,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--counts", help="TSV/TSV.GZ, rows=genes, cols=samples, RNA-seq counts")
    ap.add_argument("--cq", help="TSV/TSV.GZ, rows=genes, cols=samples, qPCR Cq values")
    ap.add_argument("--out", default="-")
    a = ap.parse_args()
    if bool(a.counts) == bool(a.cq):
        raise SystemExit("give exactly one of --counts or --cq")

    panel = load_panel()
    genes, samples, vals = read_table(a.counts or a.cq)
    try:
        s, rep = score(vals, genes, panel, is_cq=bool(a.cq))
    except ScoreError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        raise SystemExit(2)

    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    print("sample_id\tpanel_severity_score", file=fh)
    for sid, v in zip(samples, s):
        print(f"{sid}\t{v:.8f}", file=fh)
    if fh is not sys.stdout:
        fh.close()
    print(json.dumps(rep, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
