#!/usr/bin/env python3
"""masld-liver-protein-severity-v1 -- score a liver PROTEOME for fibrosis severity.

Reads a proteins-by-samples ABUNDANCE matrix (raw intensities) and emits one severity
value per sample. No outcome labels, no reference batch, no second sample: every
constant the transform uses was frozen from the 58 training participants and ships in
weights/.

    python score.py --abundance my_intensities.tsv --out scores.tsv

Read MODEL_CARD.md before using any number this produces. In particular this instrument
has NO EXTERNAL VALIDATION -- one cohort, 58 participants, zero held-out cohorts -- and
its ORDERING is meaningful while its absolute level is not a Kleiner stage.
"""
import argparse
import gzip
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WEIGHTS = os.path.join(HERE, "weights", "masld_liver_protein_severity_v1_weights.npz")

# No proteomics denylist is known to this release. The field is empty rather than
# populated by guesswork; see MODEL_CARD.md deviation D2.
DENY = set()

ABSOLUTE_COVERAGE_FLOOR = 0.01
# "Is this already logged?" cannot be asked of intensities the way the RNA release asks it
# of counts. log2CPM is bounded by 19.9, so a small maximum gives counts away; protein
# intensities reach 3.087e8, whose log2 is 28.2, ABOVE any such bound. The discriminator
# that does work is DYNAMIC RANGE, measured on the training matrix:
#     raw intensities   max/median = 1.24e4
#     the same, logged  max/median = 1.93
# Anything under 20 is logged (or otherwise compressed) and is refused.
MIN_DYNAMIC_RANGE = 20.0
TRAIN_DYNAMIC_RANGE = 1.24e4


class ScoreError(RuntimeError):
    """Every refusal in this module raises. Nothing here warns and continues."""


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def normalise_ids(raw):
    """Uppercased and whitespace-stripped. Multi-accession entries keep their full form;
    the join also tries the primary token before the first ';', because a caller's search
    engine commonly reports one accession where the deposit reported a group."""
    return np.array([str(g).strip().upper() for g in raw])


def primary_token(x):
    return str(x).split(";")[0].strip().upper()


def check_denylist(*blobs):
    hits = sorted({a for a in DENY for b in blobs if b and a in str(b).upper()})
    if hits:
        raise ScoreError(f"REFUSED: input references {hits}.")


def read_matrix(path, ids_path=None, samples_path=None):
    """Returns (matrix proteins-by-samples, ids, sample names)."""
    if path.endswith(".npz"):
        z = np.load(path, allow_pickle=True)
        if "abundance" not in z.files:
            raise ScoreError(f"npz has no 'abundance' array; keys are {z.files}")
        mat = np.asarray(z["abundance"], dtype=np.float64)
        if "ids" in z.files:
            ids = np.asarray(z["ids"]).astype(str)
        elif ids_path:
            ids = np.array([l.strip() for l in _open(ids_path) if l.strip()])
        else:
            raise ScoreError("npz has no 'ids' array and --ids was not given")
        samples = (np.asarray(z["samples"]).astype(str) if "samples" in z.files
                   else np.array([f"S{i+1}" for i in range(mat.shape[1])]))
    else:
        rows, header = [], None
        with _open(path) as fh:
            for i, line in enumerate(fh):
                parts = line.rstrip("\n").split("\t")
                if i == 0:
                    header = parts
                    continue
                rows.append(parts)
        if header is None or not rows:
            raise ScoreError(f"{path} has a header and no data rows")
        widths = {len(r) for r in rows}
        if len(widths) != 1:
            raise ScoreError(f"ragged data rows: field counts {sorted(widths)}. "
                             f"A ragged file is usually an unnamed index column.")
        w = widths.pop()
        if w != len(header) and w != len(header) + 1:
            raise ScoreError(f"header has {len(header)} fields but data rows have {w}")
        off = 1 if w == len(header) + 1 else 0
        ids = np.array([r[0] for r in rows])
        samples = np.array(header[1 + off:]) if off == 0 else np.array(header)
        if len(samples) != w - 1:
            samples = np.array(header[-(w - 1):])
        vals = [[np.nan if (v.strip() in ("", "NA", "NaN", "nan", "NULL")) else float(v)
                 for v in r[1:]] for r in rows]
        mat = np.asarray(vals, dtype=np.float64)
    if ids_path and path.endswith(".npz") is False:
        ids = np.array([l.strip() for l in _open(ids_path) if l.strip()])
    if samples_path:
        samples = np.array([l.strip() for l in _open(samples_path) if l.strip()])
    if mat.ndim != 2:
        raise ScoreError(f"abundance must be 2-D; got shape {mat.shape}")
    if len(ids) != mat.shape[0]:
        raise ScoreError(f"{len(ids)} identifiers for {mat.shape[0]} matrix rows")
    if len(samples) != mat.shape[1]:
        raise ScoreError(f"abundance has {mat.shape[1]} columns and {len(samples)} sample ids")
    if np.isinf(mat).any():
        raise ScoreError("abundance contains inf")
    if np.nanmin(mat) < 0:
        raise ScoreError("abundance contains negative values; this is an intensity matrix, "
                         "not a log ratio")
    return mat, ids, samples


def collapse_duplicates(mat, ids):
    """Duplicate identifiers are SUMMED, matching the RNA release's convention. NaN is
    treated as absent, not as zero: a sum over all-NaN stays NaN."""
    uniq, inv = np.unique(ids, return_inverse=True)
    if len(uniq) == len(ids):
        return mat, ids, 0
    out = np.full((len(uniq), mat.shape[1]), np.nan)
    for j in range(len(uniq)):
        block = mat[inv == j]
        allnan = np.isnan(block).all(axis=0)
        s = np.nansum(block, axis=0)
        s[allnan] = np.nan
        out[j] = s
    return out, uniq, len(ids) - len(uniq)


def qn_apply(X, ref):
    """Quantile-normalise each row (sample) of a samples-by-features matrix onto a frozen
    reference distribution, interpolating the reference at (i-0.5)/m so the same
    distribution is used at any axis coverage."""
    ref = np.sort(np.asarray(ref, dtype=np.float64))
    n = len(ref)
    out = np.empty_like(X)
    for i in range(X.shape[0]):
        row = X[i]
        m = row.shape[0]
        order = np.argsort(row, kind="mergesort")
        pos = (np.arange(m) + 0.5) / m
        out[i, order] = np.interp(pos, (np.arange(n) + 0.5) / n, ref)
    return out


def score(abundance, ids, weights_path=DEFAULT_WEIGHTS, min_coverage=None,
          id_type="accession", quantile_normalize=False):
    """abundance: proteins-by-samples intensities. ids: raw identifiers."""
    W = np.load(weights_path, allow_pickle=True)
    axis = np.asarray(W["protein_axis"]).astype(str)
    axis_symbol = np.asarray(W["protein_axis_symbol"]).astype(str)
    w = np.asarray(W["w"], dtype=np.float64)
    med = np.asarray(W["impute_median"], dtype=np.float64)
    mu = np.asarray(W["mu"], dtype=np.float64)
    sd = np.asarray(W["sd"], dtype=np.float64)
    intercept = float(W["intercept"])
    measured_threshold = float(W["min_axis_coverage"])
    qn_ref = np.asarray(W["qn_reference"], dtype=np.float64) if "qn_reference" in W.files else None

    if not (len(axis) == len(axis_symbol) == len(w) == len(med) == len(mu) == len(sd)):
        raise ScoreError("weights file is inconsistent: axis, symbol, w, median, mu and sd "
                         "differ in length")
    if quantile_normalize and qn_ref is None:
        raise ScoreError("--quantile-normalize was asked for but the weights carry no reference")

    thr = measured_threshold if min_coverage is None else float(min_coverage)
    if thr < ABSOLUTE_COVERAGE_FLOOR:
        raise ScoreError(f"--min-coverage {thr} is below the absolute floor "
                         f"{ABSOLUTE_COVERAGE_FLOOR}. Below that the join is degenerate and "
                         f"no score is meaningful.")
    if id_type not in ("accession", "symbol"):
        raise ScoreError(f"--id-type must be accession or symbol, not {id_type!r}")

    g = normalise_ids(ids)
    abundance, g, n_collapsed = collapse_duplicates(abundance, g)

    key = axis if id_type == "accession" else axis_symbol
    pos = {k: i for i, k in enumerate(g)}
    matched, axis_hit, n_primary = [], [], 0
    for i, a in enumerate(key):
        a = a.strip().upper()
        if a in pos:
            matched.append(pos[a]); axis_hit.append(i)
        elif id_type == "accession":
            p = primary_token(a)
            if p in pos:
                matched.append(pos[p]); axis_hit.append(i); n_primary += 1
    matched = np.asarray(matched, dtype=int)
    axis_hit = np.asarray(axis_hit, dtype=int)
    n_matched = len(matched)
    frac = n_matched / len(axis)

    if n_matched == 0:
        raise ScoreError(
            f"ZERO PROTEIN JOIN: none of the {len(axis)} axis entries matched your {len(g)} "
            f"identifiers under --id-type {id_type}. The accession axis is UniProt "
            f"accessions as deposited, e.g. {axis[:3].tolist()}; the symbol axis is HGNC "
            f"symbols, e.g. {axis_symbol[:3].tolist()}. Your first three are "
            f"{g[:3].tolist()}. Convert your identifiers, or pass --id-type symbol.")
    if frac < thr:
        raise ScoreError(
            f"AXIS COVERAGE {frac:.4f} ({n_matched}/{len(axis)}) is below the threshold "
            f"{thr:.2f}. The threshold was MEASURED, not chosen: MODEL_CARD.md section 4 "
            f"reports leave-one-out Spearman as a function of coverage, and below this point "
            f"the score stops beating a single-protein baseline. Refusing rather than "
            f"returning a number.")

    if np.nanmin(abundance) < 0:
        raise ScoreError(
            f"the matrix contains negative values (minimum {np.nanmin(abundance):.4g}). These "
            f"are intensities, not log ratios or differences. Pass raw abundances.")
    mx, md = float(np.nanmax(abundance)), float(np.nanmedian(abundance))
    if not np.isfinite(md) or md <= 0:
        raise ScoreError(f"the median over the matched proteins is {md:.4g}; a matrix of "
                         f"intensities has a positive median.")
    dyn = mx / md
    if dyn < MIN_DYNAMIC_RANGE:
        raise ScoreError(
            f"dynamic range max/median is {dyn:.4g}, below {MIN_DYNAMIC_RANGE:g}. The training "
            f"intensities have {TRAIN_DYNAMIC_RANGE:.3g}; the SAME matrix after log2 has 1.93. "
            f"A matrix this compressed is already log-transformed or otherwise scaled. The "
            f"model applies its own log2, and passing logged values squares the transform. "
            f"Pass raw intensities.")

    sub = abundance[matched, :].T                                   # samples x matched
    allnan_samples = np.isnan(sub).all(axis=1)
    if allnan_samples.any():
        bad = [str(i) for i in np.where(allnan_samples)[0]]
        raise ScoreError(f"samples at column index {bad} are entirely missing over the "
                         f"matched proteins")
    n_imputed = int(np.isnan(sub).sum())
    sub = np.where(np.isnan(sub), np.broadcast_to(med[axis_hit], sub.shape), sub)
    logv = np.log2(np.clip(sub, 1e-6, None))
    if quantile_normalize:
        logv = qn_apply(logv, qn_ref)
    z = (logv - mu[axis_hit]) / sd[axis_hit]
    severity = intercept + z @ w[axis_hit]

    return dict(
        severity=severity, n_axis=len(axis), n_matched=n_matched, axis_coverage=frac,
        id_type=id_type, n_matched_via_primary_accession=n_primary,
        measured_threshold=measured_threshold, applied_threshold=thr,
        below_measured_threshold=bool(thr < measured_threshold),
        n_input_ids=len(ids), n_duplicate_rows_summed=n_collapsed,
        n_missing_values_imputed=n_imputed,
        quantile_normalized=bool(quantile_normalize),
        n_samples=int(severity.shape[0]))


def main(argv=None):
    ap = argparse.ArgumentParser(description="score a liver proteome for fibrosis severity")
    ap.add_argument("--abundance", required=True,
                    help="proteins-by-samples raw intensities, .tsv/.tsv.gz/.npz")
    ap.add_argument("--ids", help="one identifier per line, if not the first column")
    ap.add_argument("--samples", help="one sample id per line, if not the header")
    ap.add_argument("--out", required=True)
    ap.add_argument("--weights", default=DEFAULT_WEIGHTS)
    ap.add_argument("--id-type", default="accession", choices=("accession", "symbol"))
    ap.add_argument("--quantile-normalize", action="store_true",
                    help="map each sample onto the frozen training distribution first. Off by "
                         "default: the shipped recipe is the one that was validated.")
    ap.add_argument("--min-coverage", type=float, default=None,
                    help="override the measured axis-coverage threshold. Lowering it is "
                         "recorded in the output and is the caller's responsibility.")
    ap.add_argument("--json", help="write the join diagnostics here")
    a = ap.parse_args(argv)

    check_denylist(a.abundance, a.ids, a.samples)
    mat, ids, samples = read_matrix(a.abundance, a.ids, a.samples)
    r = score(mat, ids, a.weights, a.min_coverage, a.id_type, a.quantile_normalize)

    with open(a.out, "w") as fh:
        fh.write("sample_id\tseverity\taxis_coverage\n")
        for s, v in zip(samples, r["severity"]):
            fh.write(f"{s}\t{v:.8f}\t{r['axis_coverage']:.6f}\n")
    if r["below_measured_threshold"]:
        sys.stderr.write(
            f"WARNING: scored at coverage {r['axis_coverage']:.4f} with the threshold "
            f"overridden to {r['applied_threshold']:.2f}, below the measured "
            f"{r['measured_threshold']:.2f}. This is recorded in the output.\n")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                       for k, v in r.items()}, fh, indent=2)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ScoreError as e:
        sys.stderr.write(f"REFUSED: {e}\n")
        sys.exit(2)
