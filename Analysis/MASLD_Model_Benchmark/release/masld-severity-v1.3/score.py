#!/usr/bin/env python3
"""masld-pairwise-conformal-v1 (CANDIDATE, not released) - within-cohort pairwise ordering calls
for bulk liver RNA-seq, with a pool-calibrated abstention threshold.

Reads a genes-by-samples COUNTS matrix. For every pair of samples in one cohort (batch) it returns
either an ordering call ("i above j" in fibrosis severity) or an abstention. It does not return an
absolute stage and never compares samples across cohorts.

    python score.py --counts counts.tsv --out scores.tsv --pairs pairs.tsv [--batches batches.tsv]

Recipe (PRESPEC_FROZEN.md section 11 with AMENDMENT_A1.md A1.4, B1 of model-program-20261015T223347Z):
  1. ids: `_PAR_Y` and version suffixes stripped, duplicates summed (masld-severity-v1.2 recipe);
  2. the matched gene set must equal one of the three shipped planes EXACTLY (full 26,629 genes;
     GSE225740-type 16,108; GSE213621-type 26,531). Each plane's weights and threshold were fitted
     and sealed on that gene set only. Any other gene set is refused; --plane restricts a superset;
  3. log2 CPM over the plane's genes, ranks averaged within ties, mapped onto the frozen
     qn_reference re-interpolated to the plane's gene count, frozen per-gene mu and sd
     (masld-severity-v1.1/v1.2 `qn_apply`; one sample at a time, no batch statistics);
  4. refit score f = intercept + z . w (ridge refit on the 449-unit development pool);
  5. for i, j in one batch, d = f_i - f_j: call "i above j" if d > q, "j above i" if d < -q,
     abstain if |d| <= q, with q the sealed 0.90 threshold of the plane.

The coverage the 0.90 threshold gives is empirical and marginal over discordant within-cohort
pairs. It has no finite-sample guarantee and is not P(correct | call). Read MODEL_CARD.md.
"""
import argparse
import gzip
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WEIGHTS = os.path.join(HERE, "weights", "b1_pairwise_refit_weights.npz")

DENY = {"GSE193084", "GSE192959", "GSE193080", "GSE200460"}
VERSION_SUFFIX = re.compile(r"\.\d+$")
PAR_Y_SUFFIX = re.compile(r"_PAR_Y$", re.IGNORECASE)
LOG_LIKE_MAX = 25.0


class ScoreError(RuntimeError):
    """Every refusal raises."""


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def normalise_ids(raw):
    """As masld-severity-v1.2 score.py: `_PAR_Y` stripped, then the version suffix, uppercased."""
    out = np.empty(len(raw), dtype=object)
    for i, g in enumerate(raw):
        s = str(g).strip().upper()
        out[i] = VERSION_SUFFIX.sub("", PAR_Y_SUFFIX.sub("", s))
    return out.astype(str)


def check_denylist(*blobs):
    hits = sorted({a for a in DENY for b in blobs if b and a in str(b).upper()})
    if hits:
        raise ScoreError(f"REFUSED: input references {hits} (HCC SuperSeries GSE193084 and its "
                         f"siblings are never ingested).")


def read_counts(counts_path):
    """Return (counts genes-by-samples, gene ids, sample ids). Same TSV rules as v1.2."""
    with _open(counts_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows, gid = [], []
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            gid.append(p[0])
            rows.append(p[1:])
    if not rows:
        raise ScoreError(f"{counts_path} has a header and no data rows")
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise ScoreError(f"ragged data rows: field counts {sorted(widths)}")
    n_field = widths.pop()
    samples = np.array(header[1:] if len(header) == n_field + 1 else header)
    if len(samples) != n_field:
        raise ScoreError(f"header carries {len(header)} names for {n_field} data fields")
    mat = np.array(rows, dtype=np.float64)
    if not np.isfinite(mat).all():
        raise ScoreError("counts contain NaN or inf")
    if (mat < 0).any():
        raise ScoreError("counts contain negative values")
    if len(set(samples.tolist())) != len(samples):
        raise ScoreError("sample ids are not unique")
    return mat, np.array(gid), samples


def collapse_duplicates(mat, genes):
    uniq, inv = np.unique(genes, return_inverse=True)
    if len(uniq) == len(genes):
        return mat, genes, 0
    out = np.zeros((len(uniq), mat.shape[1]), dtype=np.float64)
    np.add.at(out, inv, mat)
    return out, uniq, len(genes) - len(uniq)


def qn_reference_at(m, ref):
    p = len(ref)
    if m == p:
        return ref
    return np.interp((np.arange(1, m + 1) - 0.5) / m, (np.arange(1, p + 1) - 0.5) / p, ref)


def qn_apply(X, ref):
    """Per sample: ranks averaged within ties, mapped onto the frozen reference (v1.2 qn_apply)."""
    n, m = X.shape
    r = qn_reference_at(m, ref)
    pos = np.arange(1, m + 1, dtype=np.float64)
    out = np.empty_like(X)
    for i in range(n):
        x = X[i]
        order = np.argsort(x, kind="stable")
        rk = np.empty(m, dtype=np.float64)
        rk[order] = np.arange(1, m + 1, dtype=np.float64)
        uniq, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
        if len(uniq) != m:
            sums = np.zeros(len(uniq))
            np.add.at(sums, inv, rk)
            rk = (sums / cnt)[inv]
        out[i] = np.interp(rk, pos, r)
    return out


def load_weights(path):
    W = np.load(path, allow_pickle=False)
    axis = np.asarray(W["gene_axis"]).astype(str)
    planes = {}
    for k in np.asarray(W["planes"]).astype(str):
        pos = np.asarray(W[f"{k}__axis_positions"], np.int64)
        planes[k] = dict(pos=pos, genes=axis[pos], w=np.asarray(W[f"{k}__w"], float),
                         intercept=float(W[f"{k}__intercept"]), lam=float(W[f"{k}__lambda"]),
                         q=float(W[f"{k}__q_0.9"]))
        if len(planes[k]["w"]) != len(pos):
            raise ScoreError(f"weights file inconsistent for plane {k}")
    return dict(axis=axis, mu=np.asarray(W["mu"], float), sd=np.asarray(W["sd"], float),
                ref=np.asarray(W["qn_reference"], float), planes=planes,
                seal_sha=str(W["thresholds_sealed_sha256"]))


def score(counts, genes, weights_path=DEFAULT_WEIGHTS, plane=None):
    """counts: genes-by-samples raw counts; genes: raw ids. Returns a result dict."""
    W = load_weights(weights_path)
    g = normalise_ids(genes)
    n_par_y = int(sum(1 for x in genes if PAR_Y_SUFFIX.search(str(x).strip().upper())))
    n_versioned = int(sum(1 for x in genes
                          if VERSION_SUFFIX.search(PAR_Y_SUFFIX.sub("", str(x).strip().upper()))))
    counts, g, n_collapsed = collapse_duplicates(counts, g)
    row = {gene: i for i, gene in enumerate(g)}
    axis_set = set(W["axis"].tolist())
    matched = {gene for gene in g if gene in axis_set}
    if not matched:
        raise ScoreError(f"ZERO GENE JOIN: no input id matched the axis (unversioned Ensembl gene "
                         f"ids, e.g. {W['axis'][:3].tolist()}); input starts {g[:3].tolist()}")
    if plane is None:
        hits = [k for k, P in W["planes"].items() if set(P["genes"].tolist()) == matched]
        if not hits:
            near = {k: dict(n_plane=len(P["genes"]),
                            n_plane_missing=len(set(P["genes"].tolist()) - matched),
                            n_extra=len(matched - set(P["genes"].tolist())))
                    for k, P in W["planes"].items()}
            raise ScoreError(
                f"MATCHED GENE SET ({len(matched)} axis genes) equals no shipped plane. Each plane's "
                f"weights and sealed threshold hold for that exact gene set only: {near}. If your "
                f"input is a superset of one plane, pass --plane <name> to restrict it.")
        plane = hits[0]
    if plane not in W["planes"]:
        raise ScoreError(f"unknown plane {plane}; shipped planes {sorted(W['planes'])}")
    P = W["planes"][plane]
    missing = [x for x in P["genes"] if x not in row]
    if missing:
        raise ScoreError(f"plane {plane} needs {len(P['genes'])} genes; {len(missing)} are absent "
                         f"(e.g. {missing[:3]})")
    full_lib = counts.sum(axis=0)
    if counts.shape[1] >= 2 and full_lib.min() > 0 and full_lib.std() / full_lib.mean() < 1e-6:
        raise ScoreError("every column sums to the same value: an already-normalized matrix, not counts")
    if float(counts.max()) < LOG_LIKE_MAX:
        raise ScoreError(f"largest value {counts.max():.4g} < {LOG_LIKE_MAX:g}: looks log-transformed")
    sub = counts[[row[x] for x in P["genes"]], :]
    lib = sub.sum(axis=0)
    if (lib <= 0).any():
        raise ScoreError(f"samples at column index {np.where(lib <= 0)[0].tolist()} have zero "
                         f"counts over the plane's genes")
    logcpm = np.log2(sub / lib[None, :] * 1e6 + 1.0).T            # samples x plane genes
    z = (qn_apply(logcpm, W["ref"]) - W["mu"][P["pos"]]) / W["sd"][P["pos"]]
    # one dot product per sample, so a sample's score does not depend on which other samples share
    # the call (a matrix-vector product may change summation order with the batch size)
    f = P["intercept"] + np.array([float(np.dot(zi, P["w"])) for zi in z])
    return dict(refit_score=f, plane=plane, q=P["q"], lam=P["lam"], n_plane_genes=int(len(P["genes"])),
                n_matched_axis_genes=int(len(matched)), n_input_genes=int(len(genes)),
                n_versioned_ids_stripped=n_versioned, n_par_y_ids_stripped=n_par_y,
                n_duplicate_rows_summed=int(n_collapsed),
                zero_fraction_over_plane=(sub == 0).mean(axis=0),
                thresholds_sealed_sha256=W["seal_sha"], finite_sample_guarantee=False,
                nominal=0.90)


def pair_calls(f, q):
    """All i<j pairs. Returns (i, j, d, call) with call in {+1: i above j, -1: j above i, 0: abstain}.
    The set holds 'i above j' iff d >= -q and 'j above i' iff d <= q; q > 0 here, so no set is
    empty and a single ordering is exactly |d| > q."""
    iu, ju = np.triu_indices(len(f), 1)
    d = f[iu] - f[ju]
    in_ij, in_ji = d >= -q, d <= q
    call = np.where(in_ij & ~in_ji, 1, np.where(in_ji & ~in_ij, -1, 0))
    return iu, ju, d, call


def main(argv=None):
    ap = argparse.ArgumentParser(description="Within-cohort pairwise fibrosis ordering calls "
                                             "(candidate). No labels are accepted. Read MODEL_CARD.md.")
    ap.add_argument("--counts", required=True, help="genes-by-samples raw counts TSV/TSV.GZ")
    ap.add_argument("--out", required=True, help="per-sample output TSV")
    ap.add_argument("--pairs", help="pair output TSV (within batch only)")
    ap.add_argument("--batches", help="TSV with columns sample_id, batch. Pairs are formed only "
                                      "within a batch. Default: all samples are one cohort.")
    ap.add_argument("--json", help="JSON report")
    ap.add_argument("--plane", help="restrict a superset input to one shipped plane")
    ap.add_argument("--weights", default=DEFAULT_WEIGHTS)
    ap.add_argument("--accession", help="GEO accession, checked against the denylist")
    a = ap.parse_args(argv)

    check_denylist(a.accession, a.counts)
    mat, genes, samples = read_counts(a.counts)
    check_denylist(*samples.tolist())
    r = score(mat, genes, a.weights, a.plane)
    f = r["refit_score"]
    if a.batches:
        bmap = {}
        with open(a.batches) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            si, bi = hdr.index("sample_id"), hdr.index("batch")
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                bmap[p[si]] = p[bi]
        miss = [s for s in samples if s not in bmap]
        if miss:
            raise ScoreError(f"{len(miss)} samples have no batch, e.g. {miss[:3]}")
        batch = np.array([bmap[s] for s in samples])
    else:
        batch = np.array(["all"] * len(samples))
        print("no --batches given: every sample is treated as one cohort. Pairs across cohorts "
              "are outside what this threshold was calibrated on.", file=sys.stderr)
    with open(a.out, "w") as fh:
        fh.write("sample_id\tbatch\trefit_score\tplane\n")
        for s, b, v in zip(samples, batch, f):
            fh.write(f"{s}\t{b}\t{v:.17g}\t{r['plane']}\n")
    n_pairs = n_calls = 0
    if a.pairs:
        lab = {1: "i_above_j", -1: "j_above_i", 0: "abstain"}
        with open(a.pairs, "w") as fh:
            fh.write("sample_i\tsample_j\tbatch\td\tcall\n")
            for b in sorted(set(batch.tolist())):
                idx = np.where(batch == b)[0]
                if len(idx) < 2:
                    continue
                iu, ju, d, call = pair_calls(f[idx], r["q"])
                n_pairs += len(d)
                n_calls += int((call != 0).sum())
                for x, y_, dd, c in zip(idx[iu], idx[ju], d, call):
                    fh.write(f"{samples[x]}\t{samples[y_]}\t{b}\t{dd:.17g}\t{lab[int(c)]}\n")
    rep = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()}
    rep.update(n_samples=int(len(samples)), n_batches=int(len(set(batch.tolist()))),
               n_pairs=n_pairs, n_pair_calls=n_calls)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=2)
    print(f"scored {len(samples)} samples on plane {r['plane']} ({r['n_plane_genes']} genes); "
          f"q(0.90) {r['q']:.6f}; {n_calls} of {n_pairs} within-batch pairs called", file=sys.stderr)
    print("coverage at 0.90 is empirical, marginal over discordant within-cohort pairs, with no "
          "finite-sample guarantee; it is not P(correct | call)", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except ScoreError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)
