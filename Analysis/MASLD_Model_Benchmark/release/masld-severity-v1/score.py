#!/usr/bin/env python3
"""masld-severity-v1 — score bulk liver RNA-seq for latent fibrosis severity.

Reads a genes-by-samples COUNTS matrix and emits one latent severity value per
sample with an interval. Needs no outcome labels, no reference batch and no
second sample: one biopsy is enough, because every constant the transform uses
was frozen from the training pool and ships in weights/.

    python score.py --counts counts.tsv --out scores.tsv

Read MODEL_CARD.md before using any number this produces. In particular the
score is CROSS-SECTIONAL and carries an uncorrected per-cohort offset, so its
ordering is meaningful and its absolute level is not.
"""
import argparse
import gzip
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WEIGHTS = os.path.join(HERE, "weights", "masld_severity_v1_weights.npz")
DEFAULT_CALIB = os.path.join(HERE, "weights", "calibration.json")

# Cohorts that must never be ingested. GSE193084 is an HCC SuperSeries and the
# other three are its siblings. GSE193066, the reliability validation set, is a
# SubSeries of that SuperSeries, which is exactly why the SuperSeries and its
# other children are named here rather than left to judgement.
DENY = {"GSE193084", "GSE192959", "GSE193080", "GSE200460"}

VERSION_SUFFIX = re.compile(r"\.\d+$")
ABSOLUTE_COVERAGE_FLOOR = 0.01
MIN_MEDIAN_LIBRARY = 1e5
# log2CPM cannot exceed log2(1e6 + 1) = 19.93, so anything under this is logged data.
LOG_LIKE_MAX = 25.0


class ScoreError(RuntimeError):
    """Every refusal in this module raises. Nothing here warns and continues."""


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def normalise_ids(raw):
    """Version suffix stripped, uppercased. This is the training recipe and it is
    not optional: the axis is unversioned Ensembl stable gene ids."""
    return np.array([VERSION_SUFFIX.sub("", str(g).strip()).upper() for g in raw])


def check_denylist(*blobs):
    hits = sorted({a for a in DENY for b in blobs if b and a in str(b).upper()})
    if hits:
        raise ScoreError(
            f"REFUSED: input references {hits}. GSE193084 is an HCC SuperSeries and "
            f"{sorted(DENY - {'GSE193084'})} are its siblings. These are never ingested by this "
            f"model. If the accession appears only in a file path, rename the path.")
    return hits


def read_counts(counts_path, genes_path=None, samples_path=None, orientation="genes_x_samples"):
    """Return (counts genes-by-samples, gene ids, sample ids). Counts, not CPM, not logs."""
    if counts_path.endswith(".npz"):
        z = np.load(counts_path, allow_pickle=True)
        if "counts" not in z.files:
            raise ScoreError(f"npz has no 'counts' array; keys are {z.files}")
        mat = np.asarray(z["counts"], dtype=np.float64)
        if genes_path:
            genes = np.array([ln.strip() for ln in _open(genes_path) if ln.strip()])
        elif "genes" in z.files:
            genes = np.asarray(z["genes"]).astype(str)
        else:
            raise ScoreError("npz has no 'genes' array and --genes was not given")
        if samples_path:
            samples = np.array([ln.strip() for ln in _open(samples_path) if ln.strip()])
        elif "samples" in z.files:
            samples = np.asarray(z["samples"]).astype(str)
        else:
            samples = np.array([f"sample_{i+1}" for i in range(mat.shape[1])])
    else:
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
        # A header one name short of the data is an unnamed index column, which is
        # endemic. Resolve it here rather than letting every column shift by one.
        samples = np.array(header[1:] if len(header) == n_field + 1 else header)
        if len(samples) != n_field:
            raise ScoreError(
                f"header carries {len(header)} names for {n_field} data fields; cannot align "
                f"columns. Give a matrix whose header names either every column or every column "
                f"but the gene id.")
        mat = np.array(rows, dtype=np.float64)
        genes = np.array(gid)
        if genes_path:
            genes = np.array([ln.strip() for ln in _open(genes_path) if ln.strip()])
            if len(genes) != mat.shape[0]:
                raise ScoreError(f"--genes has {len(genes)} ids for {mat.shape[0]} matrix rows")

    if orientation == "samples_x_genes":
        mat = mat.T
    if mat.ndim != 2:
        raise ScoreError(f"counts must be 2-D; got shape {mat.shape}")
    if mat.shape[0] != len(genes):
        raise ScoreError(
            f"counts has {mat.shape[0]} rows and {len(genes)} gene ids. If your matrix is "
            f"samples-by-genes, pass --orientation samples_x_genes.")
    if mat.shape[1] != len(samples):
        raise ScoreError(f"counts has {mat.shape[1]} columns and {len(samples)} sample ids")
    if not np.isfinite(mat).all():
        raise ScoreError("counts contain NaN or inf")
    if (mat < 0).any():
        raise ScoreError("counts contain negative values; this is a counts matrix, not a log ratio")
    return mat, genes, samples


def collapse_duplicates(mat, genes):
    """Sum rows sharing a stable id. Not cosmetic: three pseudoautosomal genes
    (P2RY8, CSF2RA, ASMTL-AS1) appear as X and PAR_Y rows with identical counts,
    and leaving both in double-counts them and corrupts the library size."""
    uniq, inv = np.unique(genes, return_inverse=True)
    if len(uniq) == len(genes):
        return mat, genes, 0
    out = np.zeros((len(uniq), mat.shape[1]), dtype=np.float64)
    np.add.at(out, inv, mat)
    return out, uniq, len(genes) - len(uniq)


def qn_apply(X, ref):
    """Map each sample's ranks onto the frozen reference distribution.

    Needs ONE sample and the shipped reference vector: no batch, no second sample.
    With m != len(ref) genes present, the reference is evaluated at the m plotting
    positions (i-0.5)/m, so the same distribution is used at any axis coverage."""
    n, m = X.shape
    p = len(ref)
    if m == p:
        r = ref
    else:
        r = np.interp((np.arange(1, m + 1) - 0.5) / m,
                      (np.arange(1, p + 1) - 0.5) / p, ref)
    pos = np.arange(1, m + 1, dtype=np.float64)
    out = np.empty_like(X)
    for i in range(n):
        x = X[i]
        order = np.argsort(x, kind="stable")
        rk = np.empty(m, dtype=np.float64)
        rk[order] = np.arange(1, m + 1, dtype=np.float64)
        # average ranks within ties, so tied genes receive one common value
        uniq, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
        if len(uniq) != m:
            sums = np.zeros(len(uniq))
            np.add.at(sums, inv, rk)
            rk = (sums / cnt)[inv]
        out[i] = np.interp(rk, pos, r)
    return out


def score(counts, genes, weights_path=DEFAULT_WEIGHTS, min_coverage=None, calib_path=DEFAULT_CALIB):
    """counts: genes-by-samples. genes: raw identifiers. Returns a result dict."""
    W = np.load(weights_path, allow_pickle=True)
    axis = np.asarray(W["gene_axis"]).astype(str)
    w = np.asarray(W["w"], dtype=np.float64)
    mu = np.asarray(W["mu"], dtype=np.float64)
    sd = np.asarray(W["sd"], dtype=np.float64)
    intercept = float(W["intercept"])
    measured_threshold = float(W["min_axis_coverage"])
    # Batch-centred interval calibration. There is deliberately NO absolute-level
    # interval: it was measured uncalibrated on both held-out cohorts, and the score
    # has no between-cohort discrimination, so no absolute interval can be calibrated.
    centred_min_batch = int(W["centred_min_batch"])
    centred_grid = np.asarray(W["centred_grid_n"], dtype=np.int64)
    centred_hw80 = np.asarray(W["centred_hw_80"], dtype=np.float64)
    centred_hw90 = np.asarray(W["centred_hw_90"], dtype=np.float64)
    standardization = str(W["standardization"])
    qn_ref = np.asarray(W["qn_reference"], dtype=np.float64) if "qn_reference" in W.files else None
    excluded = set(np.asarray(W["excluded_by_detection_filter"]).astype(str).tolist()) \
        if "excluded_by_detection_filter" in W.files else set()
    if not (len(axis) == len(w) == len(mu) == len(sd)):
        raise ScoreError("weights file is inconsistent: axis, w, mu and sd differ in length")
    if standardization == "quantile" and qn_ref is None:
        raise ScoreError("weights declare the quantile transform but carry no qn_reference")

    thr = measured_threshold if min_coverage is None else float(min_coverage)
    if thr < ABSOLUTE_COVERAGE_FLOOR:
        raise ScoreError(
            f"--min-coverage {thr} is below the absolute floor {ABSOLUTE_COVERAGE_FLOOR}. "
            f"Below that the join is degenerate and no score is meaningful.")

    g = normalise_ids(genes)
    n_versioned = int(sum(1 for x in genes if VERSION_SUFFIX.search(str(x))))
    counts, g, n_collapsed = collapse_duplicates(counts, g)

    pos = {gene: i for i, gene in enumerate(g)}
    matched = np.array([pos[a] for a in axis if a in pos])
    axis_hit = np.array([i for i, a in enumerate(axis) if a in pos])
    n_matched = len(matched)
    frac = n_matched / len(axis)

    # The detection filter is ENFORCED here, not documented elsewhere. These genes fail
    # the training-pool detection rule; their frozen sd is near zero, and leaving them in
    # produced standardized values up to 1.5e6 and destroyed test-retest reliability.
    # They are dropped from the design AND from the library-size denominator.
    n_excluded_supplied = int(sum(1 for gene in g if gene in excluded)) if excluded else 0
    if excluded:
        leaked = [a for a in axis if a in excluded]
        if leaked:
            raise ScoreError(
                f"{len(leaked)} genes are on BOTH the release axis and the detection-filter "
                f"exclusion list, e.g. {leaked[:3]}. The weights file is internally inconsistent "
                f"and must not be used.")

    if n_matched == 0:
        raise ScoreError(
            f"ZERO GENE JOIN: none of the {len(axis)} axis genes matched your {len(g)} identifiers. "
            f"The axis is unversioned Ensembl stable gene ids, e.g. {axis[:3].tolist()}. "
            f"Your first three are {g[:3].tolist()}. Convert symbols or RefSeq ids to ENSG first.")
    if frac < thr:
        raise ScoreError(
            f"AXIS COVERAGE {frac:.4f} ({n_matched}/{len(axis)}) is below the threshold {thr:.2f}. "
            f"The threshold was measured, not chosen: see MODEL_CARD.md, which reports the transfer "
            f"Spearman as a function of coverage. Below it the score degrades past the point where "
            f"it beats an in-cohort baseline. Refusing rather than returning a number.")

    # ---- three checks that the input really is counts. Each raises; none warns.
    # A CPM matrix is per-sample proportional to counts and would in fact score
    # identically, but CPM cannot be told apart from TPM or FPKM, which carry a
    # length correction and would score WRONGLY. So anything already normalized is
    # refused rather than guessed at.
    full_lib = counts.sum(axis=0)
    if counts.shape[1] >= 2 and full_lib.min() > 0:
        cv = float(full_lib.std(ddof=0) / full_lib.mean())
        if cv < 1e-6:
            raise ScoreError(
                f"every column sums to the same value ({full_lib[0]:.6g}, CV {cv:.1e}). This is an "
                f"already-normalized matrix (CPM, TPM or similar), not counts. CPM would score "
                f"correctly but is indistinguishable here from TPM and FPKM, which carry a length "
                f"correction and would score WRONGLY. Pass raw or expected counts.")
    if counts.shape[1] == 1 and abs(float(full_lib[0]) - 1e6) < 1.0:
        raise ScoreError(
            f"the single column sums to {full_lib[0]:.6g}, i.e. 1e6. This is a normalized matrix, "
            f"not counts. Pass raw or expected counts.")
    mx = float(counts.max())
    if mx < LOG_LIKE_MAX:
        raise ScoreError(
            f"the largest value in the matrix is {mx:.4g}, below {LOG_LIKE_MAX:g}. log2CPM is bounded "
            f"by log2(1e6+1) = 19.9, so a matrix whose maximum is this small is almost certainly "
            f"already log-transformed. The model applies its own log2, and passing logged values "
            f"squares the transform. Pass raw or expected counts.")

    sub = counts[matched, :]
    lib = sub.sum(axis=0)
    if (lib <= 0).any():
        bad = [str(i) for i in np.where(lib <= 0)[0]]
        raise ScoreError(f"samples at column index {bad} have zero counts over the matched genes")
    med = float(np.median(lib))
    if med < MIN_MEDIAN_LIBRARY:
        raise ScoreError(
            f"median library size over the matched genes is {med:.3g}, below {MIN_MEDIAN_LIBRARY:.0g}. "
            f"Either this is not a counts matrix, or these libraries are too shallow to score.")

    logcpm = np.log2(sub / lib[None, :] * 1e6 + 1.0).T          # samples x matched
    if standardization == "quantile":
        logcpm = qn_apply(logcpm, qn_ref)
    z = (logcpm - mu[axis_hit]) / sd[axis_hit]
    latent = intercept + z @ w[axis_hit]

    # The batch-centred score and its interval. The estimand is the sample's deviation
    # from the MEDIAN of the batch it was scored with, in Kleiner stage units - not a
    # Kleiner stage. The batch median is an unknown constant the caller never learns.
    n_samples = latent.shape[0]
    centred = latent - np.median(latent)
    if n_samples >= centred_min_batch:
        j = int(np.searchsorted(centred_grid, n_samples, side="right") - 1)
        j = max(0, min(j, len(centred_grid) - 1))
        hw80, hw90 = float(centred_hw80[j]), float(centred_hw90[j])
        hw_grid_n = int(centred_grid[j])
    else:
        hw80 = hw90 = hw_grid_n = None

    calib = {}
    if os.path.exists(calib_path):
        with open(calib_path) as fh:
            calib = json.load(fh)
    return dict(
        latent_severity=latent, centred_severity=centred,
        n_samples=n_samples, centred_min_batch=centred_min_batch,
        interval_emitted=bool(hw80 is not None),
        centred_hw80=hw80, centred_hw90=hw90, hw_calibrated_at_n=hw_grid_n,
        standardization=standardization,
        n_genes_dropped_by_detection_filter=n_excluded_supplied,
        n_detection_filter_list=len(excluded),
        n_axis=len(axis), n_matched=n_matched, axis_coverage=frac,
        measured_threshold=measured_threshold, applied_threshold=thr,
        below_measured_threshold=bool(thr < measured_threshold),
        n_input_genes=len(genes), n_versioned_ids_stripped=n_versioned,
        n_duplicate_rows_summed=n_collapsed,
        median_library_size_over_matched_genes=med,
        calibration=calib)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Score bulk liver RNA-seq counts for latent MASLD fibrosis severity.",
        epilog="No outcome labels are required or accepted. Read MODEL_CARD.md.")
    ap.add_argument("--counts", required=True,
                    help="genes-by-samples counts matrix: TSV/TSV.GZ (first column gene ids, "
                         "header sample ids) or NPZ with arrays counts/genes[/samples]")
    ap.add_argument("--genes", help="optional text file of gene ids, one per row, overriding the "
                                    "ids in --counts")
    ap.add_argument("--samples", help="optional text file of sample ids, one per row")
    ap.add_argument("--orientation", default="genes_x_samples",
                    choices=["genes_x_samples", "samples_x_genes"])
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--json", help="optional JSON report with the join diagnostics")
    ap.add_argument("--weights", default=DEFAULT_WEIGHTS)
    ap.add_argument("--min-coverage", type=float, default=None,
                    help="override the measured axis-coverage threshold. Lowering it is recorded "
                         "in the output and is your responsibility.")
    ap.add_argument("--accession", help="GEO accession of the input, checked against the denylist")
    a = ap.parse_args(argv)

    check_denylist(a.accession, a.counts, a.genes)
    mat, genes, samples = read_counts(a.counts, a.genes, a.samples, a.orientation)
    check_denylist(*samples[:2000].tolist())
    r = score(mat, genes, a.weights, a.min_coverage)

    r_ = r
    with open(a.out, "w") as fh:
        if r_["interval_emitted"]:
            fh.write("sample_id\tlatent_severity\tseverity_vs_batch_median\t"
                     "vs_batch_median_lo80\tvs_batch_median_hi80\t"
                     "vs_batch_median_lo90\tvs_batch_median_hi90\taxis_coverage\n")
            c, h8, h9 = r_["centred_severity"], r_["centred_hw80"], r_["centred_hw90"]
            for i, s_ in enumerate(samples):
                fh.write(f"{s_}\t{r_['latent_severity'][i]:.6f}\t{c[i]:.6f}\t"
                         f"{c[i]-h8:.6f}\t{c[i]+h8:.6f}\t{c[i]-h9:.6f}\t{c[i]+h9:.6f}\t"
                         f"{r_['axis_coverage']:.6f}\n")
        else:
            # Refuse to emit an interval rather than emit one that is not calibrated
            # at this batch size. The score itself is still returned.
            fh.write("sample_id\tlatent_severity\tseverity_vs_batch_median\taxis_coverage\n")
            c = r_["centred_severity"]
            for i, s_ in enumerate(samples):
                fh.write(f"{s_}\t{r_['latent_severity'][i]:.6f}\t{c[i]:.6f}\t"
                         f"{r_['axis_coverage']:.6f}\n")

    rep = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()}
    rep["n_samples"] = int(len(samples))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=2)

    print(f"scored {len(samples)} samples", file=sys.stderr)
    print(f"axis coverage {r['axis_coverage']:.4f} ({r['n_matched']}/{r['n_axis']}); "
          f"threshold applied {r['applied_threshold']:.2f} (measured {r['measured_threshold']:.2f})",
          file=sys.stderr)
    if r["n_versioned_ids_stripped"]:
        print(f"stripped a version suffix from {r['n_versioned_ids_stripped']} identifiers",
              file=sys.stderr)
    if r["n_duplicate_rows_summed"]:
        print(f"summed {r['n_duplicate_rows_summed']} duplicate stable-id rows", file=sys.stderr)
    if r["n_genes_dropped_by_detection_filter"]:
        print(f"dropped {r['n_genes_dropped_by_detection_filter']} genes that fail the training-pool "
              f"detection filter (of {r['n_detection_filter_list']} on the exclusion list); they are "
              f"excluded from the design AND from the library-size denominator", file=sys.stderr)
    print(f"standardization: {r['standardization']}", file=sys.stderr)
    if r["below_measured_threshold"]:
        print("WARNING: you lowered the coverage threshold below the measured one", file=sys.stderr)
    if r["interval_emitted"]:
        print(f"interval: on severity_vs_batch_median (NOT on an absolute Kleiner stage), "
              f"half-widths {r['centred_hw80']:.3f} at 80% and {r['centred_hw90']:.3f} at 90%, "
              f"calibrated at batch size {r['hw_calibrated_at_n']}", file=sys.stderr)
    else:
        print(f"NO INTERVAL EMITTED: {r['n_samples']} samples is below the certified minimum "
              f"batch size of {r['centred_min_batch']}. The scores are still returned. An "
              f"interval at this batch size would not be calibrated, so none is given.",
              file=sys.stderr)
    print("SCOPE: cross-sectional ordering only. No interval on the ABSOLUTE level is shipped: "
          "it was measured uncalibrated on both held-out cohorts. See MODEL_CARD.md.",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ScoreError as e:
        print(f"score.py: {e}", file=sys.stderr)
        sys.exit(2)
