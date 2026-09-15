#!/usr/bin/env python3
"""masld-panel-v1 — a 20-gene-pair liver fibrosis severity score.

Reads EITHER RNA-seq counts OR qPCR Cq values for 20 genes and returns one number per sample.

WHAT THIS SHIPS, AND WHY IT IS SHORT
    a pair list, one coefficient per pair, one intercept.
No frozen quantile reference. No per-gene mean or standard deviation. No 26,629-length array.
That is the whole point: a dense z-scored panel silently requires three frozen vectors and a
transcriptome-scale input, and cannot be run on a qPCR machine. This can.

THE SCORE
    score = intercept + sum_p coef_p * (log2 x_{a_p} - log2 x_{b_p})
where x is the gene's share of the 20-gene panel total, in CPM-of-panel, floored at 1.0.

THE Cq IDENTITY
    With x = 2^(-Cq), log2 x_a - log2 x_b = Cq_b - Cq_a.
So, ABOVE THE DETECTION FLOOR, the score is exactly a weighted sum of delta-Cq values and needs no
normaliser, no housekeeping gene and no reference sample. The floor is the only thing that breaks
that identity, and it only bites on non-detects. This is the SeptiCyte RAPID form (an FDA-cleared
two-gene delta-Cq) generalised to twenty pairs.

WHAT IT IS AND IS NOT
    It ORDERS samples by fibrosis severity. It is NOT calibrated to the Kleiner scale, it emits no
    interval, and it is not a diagnostic. See MODEL_CARD.md section 4.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
WEIGHTS = HERE / "weights/masld_panel_v1.npz"

MIN_PAIR_COVERAGE = 0.70   # refuse below this; measured basis in MODEL_CARD section 5


class ScoreError(RuntimeError):
    """Everything raises. Nothing warns and continues with a number."""


def load_panel(path=WEIGHTS):
    d = np.load(path, allow_pickle=True)
    return {
        "a": [str(x) for x in d["pair_gene_a"]],
        "b": [str(x) for x in d["pair_gene_b"]],
        "coef": np.asarray(d["coef"], float),
        "intercept": float(d["intercept"]),
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
                    f"ragged row for gene {r[0]!r}: {len(r)} fields against a {len(header)}-field "
                    f"header. An unnamed pandas index column is the usual cause.")
            genes.append(r[0])
            rows.append([float(x) for x in r[1:]])
    return normalise_ids(genes), samples, np.asarray(rows, float)


def score(values, gene_ids, panel, is_cq=False):
    """values: genes x samples. Returns (scores, report)."""
    gene_ids = normalise_ids(gene_ids)
    idx = {}
    for i, g in enumerate(gene_ids):
        idx.setdefault(g, []).append(i)
    needed = sorted(set(panel["a"]) | set(panel["b"]))
    present = [g for g in needed if g in idx]
    missing = [g for g in needed if g not in idx]

    n = values.shape[1]
    if is_cq:
        # x proportional to 2^(-Cq); a non-detect is NaN or a Cq at/above the machine max
        mat = np.full((len(needed), n), np.nan)
        for r, g in enumerate(needed):
            if g in idx:
                v = values[idx[g]].mean(axis=0)      # mean of technical replicates
                mat[r] = np.power(2.0, -v)
    else:
        mat = np.full((len(needed), n), np.nan)
        for r, g in enumerate(needed):
            if g in idx:
                mat[r] = values[idx[g]].sum(axis=0)  # sum duplicate rows (PAR_Y etc.)
        if np.nanmin(mat) < 0:
            raise ScoreError("negative values in a count matrix")

    # panel composition, in CPM-of-panel, floored. The floor is the ONLY thing that breaks the
    # exact delta-Cq identity, and it only bites on non-detects.
    tot = np.nansum(mat, axis=0)
    if np.any(tot <= 0):
        bad = np.flatnonzero(tot <= 0)[:5].tolist()
        raise ScoreError(f"panel total is zero for sample index {bad}: no panel gene detected")
    comp = np.maximum(mat / tot * 1e6, panel["floor_cpm"])
    lg = np.log2(comp)

    scores = np.full(n, np.nan)
    used = 0
    acc = np.zeros(n)
    for p, (ga, gb, c) in enumerate(zip(panel["a"], panel["b"], panel["coef"])):
        ra, rb = needed.index(ga), needed.index(gb)
        if np.all(np.isfinite(lg[ra])) and np.all(np.isfinite(lg[rb])):
            acc += c * (lg[ra] - lg[rb])
            used += 1
    cov = used / len(panel["coef"])
    if cov < MIN_PAIR_COVERAGE:
        raise ScoreError(
            f"pair coverage {cov:.3f} is below the {MIN_PAIR_COVERAGE} contract "
            f"({used} of {len(panel['coef'])} pairs evaluable). Missing genes: {missing}. "
            f"REFUSING rather than returning a score from a partial panel.")
    # renormalise to the evaluable pairs so a dropped pair does not silently shrink the score
    scores = panel["intercept"] + acc * (len(panel["coef"]) / used)

    return scores, {
        "n_samples": int(n),
        "input_scale": "cq" if is_cq else "counts",
        "panel_genes_required": len(needed),
        "panel_genes_present": len(present),
        "panel_genes_missing": missing,
        "pairs_evaluable": int(used),
        "pairs_total": int(len(panel["coef"])),
        "pair_coverage": float(cov),
        "pair_coverage_contract": MIN_PAIR_COVERAGE,
        "interval_calibration": panel["interval_calibration"],
        "score_is_an_ORDERING_not_a_stage": True,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
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
