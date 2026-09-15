#!/usr/bin/env python3
"""masld-liver-chromatin-state-v1.2 — predict liver H3K27ac chromatin STATE, and (new) a per-region H3K27ac PROFILE, from bulk RNA counts.

Give it a genes-by-samples counts matrix; get ten chromatin state scores per sample, each with the
out-of-fold accuracy measured for that component and an interval built from that component's
out-of-fold residual spread.

    python score.py --counts counts.tsv --out states.tsv
    python score.py --counts counts.tsv --out states.tsv --profile profile.npz [--profile-form states_cis] [--transport marginal]

PROFILE (v1.2). --profile writes, per sample, the predicted residualised H3K27ac over the 96,460 training regions
(float64), the region keys, and the 4,847-region reliably-predictable mask. Read MODEL_CARD.md section 3b: the
in-cohort out-of-fold skill is +0.14 of the residual variance, it transfers to other people's chromatin with a
within-well permutation p of 0.001 at n = 39, and per-region values are NOT individually validated.

WHAT THIS IS. Ten axes of liver H3K27ac that bulk RNA can recover beyond histology. They were
defined by cross-covariance between paired RNA and H3K27ac in 99 participants (GSE267145), and
every one of them was shown, out of fold, to be predictable from RNA after conditioning on
steatosis, ballooning, lobular inflammation, fibrosis, sex and 14 RNA library descriptors.

WHAT THIS IS NOT. It is not a diagnostic, a stage, a subtype, or a measure of any individual
regulatory element. It is cross-sectional. Its scores are in units of the training cohort's
component standard deviation and are comparable WITHIN a batch you score together, not against an
external number. Read MODEL_CARD.md before using any value it produces.
"""
import argparse, gzip, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WEIGHTS = os.path.join(HERE, "weights", "chromatin_state_v1_2.npz")
VERSION_SUFFIX = re.compile(r"\.\d+$")
PAR_Y_SUFFIX = re.compile(r"_PAR_Y$", re.IGNORECASE)
MIN_COVERAGE = 0.60          # measured floor, see MODEL_CARD.md section 4
MIN_MEDIAN_LIBRARY = 1e5
LOG_LIKE_MAX = 25.0
DENY = {"GSE193084", "GSE192959", "GSE193080", "GSE200460"}


class StateError(RuntimeError):
    """Every refusal raises. Nothing here warns and returns a number anyway."""


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def normalise_ids(raw):
    out = np.empty(len(raw), dtype=object)
    for i, g in enumerate(raw):
        s = str(g).strip().upper()
        out[i] = VERSION_SUFFIX.sub("", PAR_Y_SUFFIX.sub("", s))
    return out.astype(str)


def read_counts(path):
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows, gid = [], []
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            gid.append(p[0]); rows.append(p[1:])
    if not rows:
        raise StateError(f"{path} has a header and no data rows")
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise StateError(f"ragged data rows: field counts {sorted(widths)}")
    nf = widths.pop()
    samples = np.array(header[1:] if len(header) == nf + 1 else header)
    if len(samples) != nf:
        raise StateError(f"header names {len(header)} columns for {nf} data fields")
    mat = np.array(rows, dtype=np.float64)
    if not np.isfinite(mat).all():
        raise StateError("counts contain NaN or inf")
    if (mat < 0).any():
        raise StateError("counts contain negative values; this is not a counts matrix")
    return mat, np.array(gid), samples


def collapse_duplicates(mat, genes):
    uniq, inv = np.unique(genes, return_inverse=True)
    if len(uniq) == len(genes):
        return mat, genes, 0
    out = np.zeros((len(uniq), mat.shape[1]))
    np.add.at(out, inv, mat)
    return out, uniq, len(genes) - len(uniq)


def score(counts, genes, weights_path=DEFAULT_WEIGHTS, min_coverage=None):
    W = np.load(weights_path, allow_pickle=True)
    axis = np.asarray(W["gene_axis"]).astype(str)
    g = normalise_ids(genes)
    counts, g, n_collapsed = collapse_duplicates(np.asarray(counts, float), g)
    pos = {x: i for i, x in enumerate(g)}
    hit = np.array([i for i, a in enumerate(axis) if a in pos])
    matched = np.array([pos[axis[i]] for i in hit])
    frac = len(hit) / len(axis)
    thr = MIN_COVERAGE if min_coverage is None else float(min_coverage)
    if len(hit) == 0:
        raise StateError(
            f"ZERO GENE JOIN: none of the {len(axis)} axis genes matched. The axis is unversioned "
            f"Ensembl stable ids, e.g. {axis[:3].tolist()}; yours start {g[:3].tolist()}.")
    if frac < thr:
        raise StateError(
            f"AXIS COVERAGE {frac:.4f} ({len(hit)}/{len(axis)}) is below {thr:.2f}. The PCA "
            f"representation is a fixed rotation of the full axis, so a partial gene set projects "
            f"onto it with an error this release has not measured below that floor. Refusing.")
    sub = counts[matched, :]
    lib = sub.sum(0)
    if (lib <= 0).any():
        raise StateError("a sample has zero counts over the matched genes")
    med = float(np.median(lib))
    if med < MIN_MEDIAN_LIBRARY:
        raise StateError(f"median library over matched genes {med:.3g} < {MIN_MEDIAN_LIBRARY:.0g}")
    mx = float(counts.max())
    if mx < LOG_LIKE_MAX:
        raise StateError(
            f"largest value {mx:.4g} < {LOG_LIKE_MAX:g}: this looks already log-transformed. "
            f"Pass raw or expected counts.")
    full = counts.sum(0)
    if counts.shape[1] >= 2 and full.min() > 0 and float(full.std() / full.mean()) < 1e-6:
        raise StateError("every column sums to the same value: an already-normalised matrix, not counts")
    x = np.log2(sub / lib[None, :] * 1e6 + 1.0).T                      # samples x matched
    # project onto the frozen PCA basis, restricted to the genes present
    mu = np.asarray(W["pca_mean"], float)[hit]
    Vp = np.asarray(W["pca_components"], float)[hit]
    zsd = np.asarray(W["pca_score_sd"], float)
    t, *_ = np.linalg.lstsq(Vp, (x - mu).T, rcond=None)                 # least squares, so partial coverage is handled
    Z = (t.T) / zsd
    smu, ssd = np.asarray(W["state_mu"], float), np.asarray(W["state_sd"], float)
    states = ((Z - smu) / ssd) @ np.asarray(W["state_W"], float) + np.asarray(W["state_b"], float)
    # v1.1: the fixed steatosis-associated chromatin head, same PCA scores, own ridge (see MODEL_CARD.md section 3)
    ste = ((Z - np.asarray(W["ste_mu"], float)) / np.asarray(W["ste_sd"], float)) @ np.asarray(W["ste_W"], float) + np.asarray(W["ste_b"], float)
    states = np.column_stack([states, ste[:, 0] if ste.ndim == 2 else ste])
    return dict(states=states, n_axis=len(axis), n_matched=len(hit), axis_coverage=frac,
                n_duplicate_rows_summed=n_collapsed, median_library=med,
                oof_residual_sd=np.asarray(W["oof_residual_sd"], float),
                n_train=int(W["n_train"]), cohort=str(W["cohort"]), n_pcs=int(W["n_pcs"]))


# ---------------------------------------------------------------------------------------------------------------
# v1.2: the H3K27ac PROFILE head. Everything above this line is v1.1's score.py, byte for byte (asserted by test_release.sh).
PROFILE_FORMS = ("states_cis", "rrr_cis", "rrr", "rrr_offset_cis", "rrr_offset_states_cis")
MIN_BATCH_MARGINAL = 10


def profile(counts, genes, weights_path=DEFAULT_WEIGHTS, form="states_cis", transport="marginal", min_coverage=None):
    """Per-sample predicted RESIDUALISED H3K27ac over the 96,460 training regions (float64, samples x regions; float32 storage
    changed donor ranks in 2 of 83,627 regions, so the release keeps double precision).

    form      states_cis  the primary head: per region, a ridge on [10 predicted states, <=10 nearest genes within 100 kb]
              rrr_cis     per region, a ridge on [rank-20 reduced-rank RNA prediction, the same cis genes]
              rrr         the reduced-rank prediction alone
              rrr_offset_cis         RRR held fixed + a per-region ridge on the cis genes fitted to RRR's residual
              rrr_offset_states_cis  RRR held fixed + a per-region ridge on [10 states, cis genes] fitted to RRR's residual
    transport marginal    every head input is z-scored across the samples you score together and mapped onto the
                          training marginal (needs a batch of >= MIN_BATCH_MARGINAL samples; the card numbers were
                          measured this way)
              raw         training mean/sd applied directly to your log-CPM (no batch needed; transfers worse)
    The output is a residual: the training cohort's four H3K27ac concentration descriptors were regressed out of the
    target before fitting, so the profile is 'H3K27ac beyond signal concentration', in log2-CPM units of GSE267145.
    """
    if form not in PROFILE_FORMS:
        raise StateError(f"unknown profile form {form!r}; one of {PROFILE_FORMS}")
    if transport not in ("marginal", "raw"):
        raise StateError(f"unknown transport {transport!r}; 'marginal' or 'raw'")
    W = np.load(weights_path, allow_pickle=True)
    axis = np.asarray(W["gene_axis"]).astype(str)
    g = normalise_ids(genes)
    counts, g, _ = collapse_duplicates(np.asarray(counts, float), g)
    pos = {x: i for i, x in enumerate(g)}
    hit = np.array([i for i, a in enumerate(axis) if a in pos])
    if len(hit) == 0:
        raise StateError("ZERO GENE JOIN for the profile head")
    matched = np.array([pos[axis[i]] for i in hit])
    frac = len(hit) / len(axis)
    thr = MIN_COVERAGE if min_coverage is None else float(min_coverage)
    if frac < thr:
        raise StateError(f"AXIS COVERAGE {frac:.4f} below {thr:.2f}; the profile head refuses")
    sub = counts[matched, :]
    lib = sub.sum(0)
    if (lib <= 0).any():
        raise StateError("a sample has zero counts over the matched genes")
    if float(counts.max()) < LOG_LIKE_MAX:
        raise StateError("input looks already log-transformed; pass counts")
    ns = counts.shape[1]
    if transport == "marginal" and ns < MIN_BATCH_MARGINAL:
        raise StateError(f"marginal transport z-scores every input across the batch you score together and needs "
                         f">= {MIN_BATCH_MARGINAL} samples; got {ns}. Use --transport raw, whose transfer numbers are "
                         f"lower (see MODEL_CARD.md section 3b).")
    x = np.log2(sub / lib[None, :] * 1e6 + 1.0).T                       # samples x matched
    xmean, xsd = np.asarray(W["prof_x_mean"], float), np.asarray(W["prof_x_sd"], float)
    Xt = np.tile(xmean, (ns, 1))
    present = np.zeros(len(axis), bool); present[hit] = True
    if transport == "marginal":
        mu_d, sd_d = x.mean(0), x.std(0)
        z = (x - mu_d) / np.where(sd_d < 1e-12, 1.0, sd_d)
        Xt[:, hit] = xmean[hit] + z * xsd[hit]
    else:
        Xt[:, hit] = x
    PCs = ((Xt - np.asarray(W["prof_pca_mean"], float)) @ np.asarray(W["prof_pca_components"], float)) / np.asarray(W["prof_pca_score_sd"], float)
    ST = ((PCs - np.asarray(W["prof_states_mu"], float)) / np.asarray(W["prof_states_sd"], float)) @ np.asarray(W["prof_states_W"], float) + np.asarray(W["prof_states_ym"], float)
    XZ = (Xt - xmean) / xsd
    if transport == "marginal":
        ST = np.asarray(W["prof_glob_mean"], float) + (ST - ST.mean(0)) / np.where(ST.std(0) < 1e-12, 1.0, ST.std(0)) * np.asarray(W["prof_glob_sd"], float)
        XZ = (XZ - XZ.mean(0)) / np.where(XZ.std(0) < 1e-12, 1.0, XZ.std(0))
        XZ[:, ~present] = 0.0

    def apply_head(prefix, GLOB):
        cg = np.asarray(W[f"{prefix}_coef_glob"], float)
        pred = np.tile(np.asarray(W[f"{prefix}_intercept"], float), (GLOB.shape[0], 1))
        if cg.shape[1] == 1:
            pred += GLOB * cg[:, 0][None, :]
        elif cg.shape[1] > 1:
            pred += GLOB @ cg.T
        idx_all, coef_all = np.asarray(W[f"{prefix}_cis_idx"]), np.asarray(W[f"{prefix}_coef_cis"], float)
        for k in range(idx_all.shape[1]):
            idx = idx_all[:, k]; ok = idx >= 0
            if ok.any():
                pred[:, ok] += XZ[:, idx[ok]] * coef_all[ok, k][None, :]
        return pred

    if form == "states_cis":
        pred = apply_head("prof_F1", ST)
    else:
        rrr = ((((Xt - np.asarray(W["prof_rrr_xm"], float)) / np.asarray(W["prof_rrr_xs"], float)) @ np.asarray(W["prof_rrr_VtDv"], float))
               @ np.asarray(W["prof_rrr_WVr"], float)) @ np.asarray(W["prof_rrr_Vr"], float).T + np.asarray(W["prof_rrr_ym"], float)
        if transport == "marginal":
            rrr = np.asarray(W["prof_rrr_feat_mean"], float) + (rrr - rrr.mean(0)) / np.where(rrr.std(0) < 1e-12, 1.0, rrr.std(0)) * np.asarray(W["prof_rrr_feat_sd"], float)
        if form == "rrr":
            pred = rrr
        elif form == "rrr_cis":
            pred = apply_head("prof_F3", rrr)
        elif form == "rrr_offset_cis":
            pred = rrr + apply_head("prof_Ocis", ST)
        else:
            pred = rrr + apply_head("prof_Ounion", ST)
    return dict(profile=np.asarray(pred, np.float64), region_key=np.asarray(W["prof_region_key"]).astype(str),
                reliable=np.asarray(W["prof_reliable"], bool), axis_coverage=frac, n_matched=len(hit), form=form, transport=transport)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="genes-by-samples counts TSV (first column gene ids)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--json", help="optional JSON diagnostics report")
    ap.add_argument("--weights", default=DEFAULT_WEIGHTS)
    ap.add_argument("--min-coverage", type=float, default=None)
    ap.add_argument("--accession", help="GEO accession, checked against the denylist")
    ap.add_argument("--profile", help="(v1.2) also write the per-region H3K27ac profile to this .npz")
    ap.add_argument("--profile-form", default="rrr_offset_cis", choices=list(PROFILE_FORMS))
    ap.add_argument("--transport", default="marginal", choices=["marginal", "raw"])
    a = ap.parse_args(argv)
    for blob in (a.accession, a.counts):
        hits = sorted({x for x in DENY if blob and x in str(blob).upper()})
        if hits:
            raise StateError(f"REFUSED: input references {hits}, which are never ingested by this model.")
    mat, genes, samples = read_counts(a.counts)
    hits = sorted({x for x in DENY for s in samples if x in str(s).upper()})
    if hits:
        raise StateError(f"REFUSED: sample ids reference {hits}")
    r = score(mat, genes, a.weights, a.min_coverage)
    K = r["states"].shape[1] - 1                     # ten state components + the steatosis-associated head
    with open(a.out, "w") as fh:
        fh.write("sample_id\t" + "\t".join(f"chromatin_state_k{k+1}" for k in range(K)) + "\tsteatosis_chromatin_axis\taxis_coverage\n")
        for i, s in enumerate(samples):
            fh.write(f"{s}\t" + "\t".join(f"{v:.6f}" for v in r['states'][i]) + f"\t{r['axis_coverage']:.6f}\n")
    rep = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()}
    if a.profile:
        pr = profile(mat, genes, a.weights, a.profile_form, a.transport, a.min_coverage)
        np.savez_compressed(a.profile, profile=pr["profile"], region_key=pr["region_key"], reliable=pr["reliable"], sample_id=samples,
                            form=pr["form"], transport=pr["transport"], axis_coverage=pr["axis_coverage"])
        rep["profile"] = {"file": a.profile, "form": pr["form"], "transport": pr["transport"], "n_regions": int(pr["profile"].shape[1]), "n_reliable": int(pr["reliable"].sum())}
        print(f"profile: {pr['profile'].shape[0]} samples x {pr['profile'].shape[1]} regions ({pr['form']}, {pr['transport']} transport); "
              f"{int(pr['reliable'].sum())} reliably predictable regions flagged", file=sys.stderr)
    rep["n_samples"] = int(len(samples))
    if a.json:
        json.dump(rep, open(a.json, "w"), indent=2)
    print(f"scored {len(samples)} samples for {K} chromatin state components + the steatosis-associated head", file=sys.stderr)
    print(f"axis coverage {r['axis_coverage']:.4f} ({r['n_matched']}/{r['n_axis']})", file=sys.stderr)
    if r["n_duplicate_rows_summed"]:
        print(f"summed {r['n_duplicate_rows_summed']} duplicate stable-id rows", file=sys.stderr)
    print(f"trained on {r['n_train']} participants ({r['cohort']}); {r['n_pcs']}-dim RNA representation",
          file=sys.stderr)
    print("SCOPE: within-batch ordering of chromatin state. Not a stage, not a subtype, not a "
          "diagnostic. Accuracy is DEVELOPMENT-GRADE (one cohort); see MODEL_CARD.md sections 3 and 6.",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except StateError as e:
        print(f"score.py: {e}", file=sys.stderr)
        sys.exit(2)
