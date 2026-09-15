#!/usr/bin/env python3
"""C1 endpoint 2, step 5: the FOUR-model recompute, with the Borzoi column present.

Successor to `c1e2_03_recompute.py`, which ran while the two Borzoi GPU jobs (21764806,
21764807) were still queued and therefore delivered endpoint 2 as a three-family comparison.
Those jobs COMPLETED on 2026-09-14 at 21:35 and 21:40 UTC, 3,102 variants scored and 0
unresolved, after the c1-endpoint2-gpu package had already written its report. This script is
the rerun that package's own RESULTS.md says is required when Borzoi lands.

Differences from `c1e2_03_recompute.py`, and only these:
  1. `--scoredir` (read-only, the c1-endpoint2-gpu package that holds the GPU outputs) is
     separated from `--outdir` (this package). Nothing is written into `--scoredir`.
  2. The Holm family of three is completed: the E2 transfer contrast enters at its measured raw
     p instead of as a placeholder row.
  3. The Borzoi column's provenance (archived src/79 panel vs this program's GPU run) is counted
     on every matched set, and the endpoint-1 primary contrast is additionally computed on the
     wave-1 archived-Borzoi-only subset so the enlargement is a one-component change.
  4. A prediction receipt is written from the computed numbers rather than transcribed.

Everything else -- loaders, statistic, weighting, block definition, seed, draw count, matched-set
rule -- is byte-identical to `c1e2_03_recompute.py`, so the three-family numbers this reproduces
are a check on that run rather than a new recipe.

Executes IMPLEMENTATION_SPEC.md section 5 (C1). Resampling unit: the 1-Mb block. 10,000 draws,
seed 20260914, every model recomputed inside the same draw so contrasts are paired.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import time

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
ATLAS = os.path.join(ROOT, "GWAS/finemapping/results/alphagenome_atlas")

CBP_ADULT = os.path.join(
    SF, "adult_liver_chrombpnet/hepatocyte_5fold_v2/currin_calibration/currin_scored_primary.tsv.gz")
CBP_281367 = os.path.join(
    SF, "adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/currin_calibration/"
        "currin_scored_primary.tsv.gz")
DIRFEAT = os.path.join(SF, "direction_features_caqtl.tsv")
P3A_CHUNKS = os.path.join(ATLAS, "p3a-benchmarks-20260909T190214Z/raw/atlas_caqtl")
TRACK0 = os.path.join(ATLAS, "run-20260909T153939Z/tables/prediction_records")
CALIB_JSON = os.path.join(
    SF, "adult_liver_chrombpnet/hepatocyte_5fold_v2/currin_calibration/"
        "currin_external_calibration.json")
E2_RECEIPT = os.path.join(
    ROOT, "GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z/"
          "tables/defined_contrast_receipt.json")

BOOT_SEED = 20260914
N_BOOT = 10000

AG_ATAC_TRACKS = ["UBERON:0001114 ATAC-seq", "UBERON:0001115 ATAC-seq", "UBERON:0002107 ATAC-seq"]
AG_DNASE_TRACKS = ["UBERON:0001114 DNase-seq", "UBERON:0001115 DNase-seq"]
T0_ATAC_COLS = [f"raw|primary_liver|{t}" for t in AG_ATAC_TRACKS]
T0_DNASE_COLS = [f"raw|primary_liver|{t}" for t in AG_DNASE_TRACKS]

NT_RANK = {"A": 0, "C": 1, "G": 2, "T": 3}

MODELS_E2 = ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367",
             "borzoi_atac", "borzoi_dnase", "alphagenome_atac_liver",
             "alphagenome_dnase_liver", "control_allele_identity", "control_position"]
MATCH_E2 = ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367", "borzoi_atac",
            "borzoi_dnase", "alphagenome_atac_liver", "alphagenome_dnase_liver"]

MODELS_E1_TIERA = ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367", "chrombpnet_hepg2",
                   "borzoi_atac", "borzoi_dnase", "alphagenome_atac_liver",
                   "alphagenome_dnase_liver", "control_allele_identity", "control_position"]
MODELS_E1_TIERA4 = ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367",
                    "borzoi_atac", "borzoi_dnase", "alphagenome_atac_liver",
                    "alphagenome_dnase_liver", "control_allele_identity", "control_position"]

REJ = "reject_at_fwer_0.05"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def sha256(path, max_bytes=400_000_000):
    if os.path.getsize(path) > max_bytes:
        return "not_hashed_too_large"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def block_of(chrom, pos):
    c = str(chrom)
    c = c if c.startswith("chr") else f"chr{c}"
    return f"{c}~{int(pos) // 1_000_000}"


def key_of(chrom, pos, ref, alt):
    c = str(chrom)
    c = c[3:] if c.startswith("chr") else c
    return f"{c}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


# --------------------------------------------------------------------------- statistics
def _pearson(a, b):
    a = a - a.mean()
    b = b - b.mean()
    da = np.sqrt((a * a).sum())
    db = np.sqrt((b * b).sum())
    if da <= 0 or db <= 0:
        return np.nan
    return float((a * b).sum() / (da * db))


def signed_spearman(score, label):
    if len(score) < 3:
        return np.nan
    return _pearson(stats.rankdata(score), stats.rankdata(label))


def weighted_spearman(score, label, w):
    if len(score) < 3:
        return np.nan
    rs = stats.rankdata(score)
    rl = stats.rankdata(label)
    sw = w.sum()
    if sw <= 0:
        return np.nan
    ms = (w * rs).sum() / sw
    ml = (w * rl).sum() / sw
    cov = (w * (rs - ms) * (rl - ml)).sum() / sw
    vs = (w * (rs - ms) ** 2).sum() / sw
    vl = (w * (rl - ml) ** 2).sum() / sw
    if vs <= 0 or vl <= 0:
        return np.nan
    return float(cov / np.sqrt(vs * vl))


def sign_concordance(score, label, w=None):
    ok = (score != 0) & (label != 0) & np.isfinite(score) & np.isfinite(label)
    if ok.sum() == 0:
        return np.nan
    agree = (np.sign(score[ok]) == np.sign(label[ok])).astype(float)
    if w is None:
        return float(agree.mean())
    ww = w[ok]
    return float((agree * ww).sum() / ww.sum())


class BlockBootstrap:
    def __init__(self, blocks, n_boot, seed):
        order = np.argsort(blocks, kind="mergesort")
        if not np.array_equal(order, np.arange(len(blocks))):
            raise ValueError("rows must be pre-sorted by block")
        uniq, starts, sizes = np.unique(blocks, return_index=True, return_counts=True)
        self.starts = starts.astype(np.int64)
        self.sizes = sizes.astype(np.int64)
        self.n_blocks = len(uniq)
        self.n_boot = n_boot
        self.rng = np.random.default_rng(seed)

    def draws(self):
        for _ in range(self.n_boot):
            pick = self.rng.integers(0, self.n_blocks, size=self.n_blocks)
            s = self.sizes[pick]
            st = self.starts[pick]
            total = int(s.sum())
            ends = np.cumsum(s)
            offs = np.arange(total) - np.repeat(ends - s, s)
            yield np.repeat(st, s) + offs


def ci(vals):
    v = vals[np.isfinite(vals)]
    if len(v) == 0:
        return (np.nan, np.nan)
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))


def two_sided_boot_p(d):
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return (np.nan, "no_valid_draws")
    p = 2.0 * min(float((d <= 0).mean()), float((d >= 0).mean()))
    p = min(1.0, p)
    floor = 1.0 / len(d)
    if p < floor:
        return (floor, f"floored_at_1/{len(d)}_draws")
    return (float(p), "")


# --------------------------------------------------------------------------- loaders
def load_cbp(path, tag):
    cols = ["variant_id_hg38", "chr_x", "pos_hg38", "ref", "alt", "label", "beta_alt",
            "q_value", "p_nominal", "peak_id", "reference_match", "heldout_fold", "logfc",
            "jsd", "active_allele_quantile", "model_peak_overlap"]
    d = pd.read_csv(path, sep="\t", usecols=cols)
    d["key"] = [key_of(c, p, r, a) for c, p, r, a in zip(d.chr_x, d.pos_hg38, d.ref, d.alt)]
    d = d.rename(columns={"logfc": tag})
    log(f"  {tag}: {len(d)} scored SNVs, {(d.label == 1).sum()} with beta_alt")
    return d


def load_atlas_chunks():
    import anndata as ad
    frames = []
    chunks = sorted(glob.glob(os.path.join(P3A_CHUNKS, "chunk_*")))
    log(f"  reading {len(chunks)} archived Atlas caQTL chunks")
    for d in chunks:
        per = {}
        for stem, tracks in (("ATAC", AG_ATAC_TRACKS), ("DNASE", AG_DNASE_TRACKS)):
            a = ad.read_h5ad(os.path.join(d, f"{stem}.h5ad"))
            names = a.var["name"].astype(str).values
            idx = [int(np.where(names == t)[0][0]) for t in tracks]
            X = np.asarray(a.X, dtype=float)[:, idx]
            Q = np.asarray(a.layers["quantiles"], dtype=float)[:, idx]
            v = pd.Series(a.obs["variant"].astype(str).values)
            parts = v.str.replace(">", ":", regex=False).str.split(":", expand=True)
            per[stem] = pd.DataFrame({
                "key": [key_of(c, p, r, a2) for c, p, r, a2 in
                        zip(parts[0], parts[1], parts[2], parts[3])],
                f"{stem}_raw": np.nanmean(X, axis=1),
                f"{stem}_q": np.nanmean(Q, axis=1)})
        m = per["ATAC"].merge(per["DNASE"], on="key", how="outer", validate="one_to_one")
        frames.append(m)
    out = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["key"])
    return out.rename(columns={"ATAC_raw": "alphagenome_atac_liver",
                               "DNASE_raw": "alphagenome_dnase_liver",
                               "ATAC_q": "alphagenome_atac_liver_quantile",
                               "DNASE_q": "alphagenome_dnase_liver_quantile"})


def load_track0_ase(uids):
    import pyarrow.parquet as pq
    want = set(uids)
    frames = {}
    for stem, cols in (("ATAC", T0_ATAC_COLS), ("DNASE", T0_DNASE_COLS)):
        f = pq.ParquetFile(os.path.join(TRACK0, f"{stem}.parquet"))
        keep = []
        for batch in f.iter_batches(batch_size=50_000, columns=["variant_uid"] + cols):
            b = batch.to_pandas()
            keep.append(b[b.variant_uid.isin(want)])
        d = pd.concat(keep, ignore_index=True).drop_duplicates(subset=["variant_uid"])
        d[f"alphagenome_{stem.lower()}_liver"] = d[cols].mean(axis=1)
        frames[stem] = d[["variant_uid", f"alphagenome_{stem.lower()}_liver"]]
        log(f"  Track0 {stem}: {len(d)} of {len(want)} ASE sites served")
    return frames["ATAC"].merge(frames["DNASE"], on="variant_uid", how="outer")


def load_new_chrombpnet(scored, fam, tag):
    """Assemble the new held-out-fold ChromBPNet scores for the ASE sites."""
    base = os.path.join(scored, "scoring", "chrombpnet", fam)
    assign = pd.read_csv(os.path.join(base, "fold_assignment.tsv"), sep="\t")
    frames = []
    for f in range(5):
        p = os.path.join(base, "scores", f"fold{f}.variant_scores.tsv")
        d = pd.read_csv(p, sep="\t")
        d["heldout_fold"] = f
        frames.append(d)
    s = pd.concat(frames, ignore_index=True)
    s = s.rename(columns={"variant_id": "uid"})
    m = assign.merge(s[["uid", "chr", "pos", "allele1", "allele2", "logfc", "jsd",
                        "active_allele_quantile", "heldout_fold"]],
                     on="uid", how="left", suffixes=("", "_scored"))
    bad_allele = int(((m.allele1.str.upper() != m.ref.str.upper())
                      | (m.allele2.str.upper() != m.alt.str.upper())).sum())
    bad_fold = int((m.heldout_fold != m.heldout_fold_scored).sum())
    log(f"  {tag}: {len(m)} rows, allele mismatches {bad_allele}, fold mismatches {bad_fold}, "
        f"missing logfc {int(m.logfc.isna().sum())}")
    if bad_allele or bad_fold:
        raise SystemExit(f"{tag}: allele or fold mismatch in the new ChromBPNet scores")
    m = m.rename(columns={"logfc": tag})
    return m[["uid", tag, "heldout_fold"]].rename(
        columns={"heldout_fold": f"{tag}_heldout_fold"}), dict(
        n=len(m), fold_counts=m.groupby("heldout_fold").size().to_dict())


def load_new_borzoi(scored):
    """Borzoi scores from the two GPU jobs; an empty frame if they are still absent."""
    p = os.path.join(scored, "scoring", "borzoi", "borzoi_scores.tsv")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        log("  new Borzoi scores ABSENT: the b6k GPU job has not landed")
        return pd.DataFrame(columns=["canon", "borzoi_atac_new", "borzoi_dnase_new",
                                     "borzoi_atac_delta", "orient_mult"])
    d = pd.read_csv(p, sep="\t")
    d["borzoi_atac_new"] = d.borzoi_atac_delta * d.orient_mult
    d["borzoi_dnase_new"] = d.borzoi_dnase_delta * d.orient_mult
    log(f"  new Borzoi rows {len(d)}; scored {int(d.borzoi_atac_delta.notna().sum())}; "
        f"unresolved {int((d.orient_mult == 0).sum())}; "
        f"orient_mult counts {d.orient_mult.value_counts().to_dict()}")
    return d


# --------------------------------------------------------------------------- endpoint runner
def run_endpoint(df, models, label_col, weight_col, name, out_rows, boot_store):
    df = df.sort_values("block_1mb", kind="mergesort").reset_index(drop=True)
    blocks = df["block_1mb"].values
    lab = df[label_col].values.astype(float)
    w = df[weight_col].values.astype(float) if weight_col else None
    n_blocks = int(pd.unique(blocks).size)
    log(f"[{name}] n={len(df)} variants, {n_blocks} 1-Mb blocks, {len(models)} models, "
        f"weighted={'yes' if w is not None else 'no'}")

    for m in models:
        s = df[m].values.astype(float)
        if w is not None:
            out_rows.append(dict(endpoint=name, model=m, metric="signed_spearman_weighted",
                                 value=weighted_spearman(s, lab, w),
                                 n_variants=len(df), n_blocks=n_blocks))
            out_rows.append(dict(endpoint=name, model=m, metric="sign_concordance_weighted",
                                 value=sign_concordance(s, lab, w),
                                 n_variants=len(df), n_blocks=n_blocks))
        out_rows.append(dict(endpoint=name, model=m, metric="signed_spearman_unweighted",
                             value=signed_spearman(s, lab),
                             n_variants=len(df), n_blocks=n_blocks))
        out_rows.append(dict(endpoint=name, model=m, metric="sign_concordance_unweighted",
                             value=sign_concordance(s, lab),
                             n_variants=len(df), n_blocks=n_blocks))

    kinds = (["weighted", "unweighted"] if w is not None else ["unweighted"])
    draws = {k: {m: np.full(N_BOOT, np.nan) for m in models} for k in kinds}
    scores = {m: df[m].values.astype(float) for m in models}
    bb = BlockBootstrap(blocks, N_BOOT, BOOT_SEED)
    t0 = time.time()
    for i, idx in enumerate(bb.draws()):
        y = lab[idx]
        ww = w[idx] if w is not None else None
        for m in models:
            s = scores[m][idx]
            draws["unweighted"][m][i] = signed_spearman(s, y)
            if ww is not None:
                draws["weighted"][m][i] = weighted_spearman(s, y, ww)
        if (i + 1) % 2000 == 0:
            log(f"  [{name}] bootstrap {i+1}/{N_BOOT} ({time.time()-t0:.0f}s)")
    for k in kinds:
        boot_store[f"{name}|{k}"] = draws[k]
        for m in models:
            lo, hi = ci(draws[k][m])
            out_rows.append(dict(endpoint=name, model=m,
                                 metric=f"signed_spearman_{k}_boot_ci_lo", value=lo,
                                 n_variants=len(df), n_blocks=n_blocks))
            out_rows.append(dict(endpoint=name, model=m,
                                 metric=f"signed_spearman_{k}_boot_ci_hi", value=hi,
                                 n_variants=len(df), n_blocks=n_blocks))
    return df


def contrasts(boot_store, name, kind, pairs, rows, plook):
    d0 = boot_store.get(f"{name}|{kind}")
    if d0 is None:
        return
    for a, b in pairs:
        if a not in d0 or b not in d0:
            continue
        d = d0[a] - d0[b]
        lo, hi = ci(d)
        p, note = two_sided_boot_p(d)
        pa = plook.get((name, a, kind), np.nan)
        pb = plook.get((name, b, kind), np.nan)
        rows.append(dict(endpoint=name, weighting=kind, contrast=f"{a} - {b}",
                         point=pa - pb, boot_mean=float(np.nanmean(d)),
                         ci_lo=lo, ci_hi=hi, excludes_zero=bool(lo > 0 or hi < 0),
                         boot_p_two_sided=p, note=note,
                         n_draws=int(np.isfinite(d).sum())))


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="this package, written")
    ap.add_argument("--scoredir", required=True,
                    help="c1-endpoint2-gpu package holding the GPU scores; read only")
    args = ap.parse_args()
    out = args.outdir
    scored = args.scoredir
    if os.path.abspath(out) == os.path.abspath(scored):
        raise SystemExit("--outdir must differ from --scoredir; --scoredir is read-only")
    tables = os.path.join(out, "tables")
    os.makedirs(tables, exist_ok=True)
    qc = {}

    log("=== loading ===")
    cbp_a = load_cbp(CBP_ADULT, "chrombpnet_adult_hep")
    cbp_b = load_cbp(CBP_281367, "chrombpnet_adult_hep_gse281367")
    dirfeat = pd.read_csv(DIRFEAT, sep="\t")
    dirfeat["key"] = dirfeat["canon"].astype(str)
    dirfeat = dirfeat.rename(columns={
        "cbp_logfc_labelframe": "chrombpnet_hepg2",
        "borzoi_atac_delta_labelframe": "borzoi_atac_archived",
        "borzoi_dnase_delta_labelframe": "borzoi_dnase_archived",
        "ag_atac_delta_labelframe": "ag_atac_src79",
        "ag_dnase_delta_labelframe": "ag_dnase_src79"})
    atlas = load_atlas_chunks()
    bz = load_new_borzoi(scored)
    bz_atac = dict(zip(bz.canon, bz.borzoi_atac_new))
    bz_dnase = dict(zip(bz.canon, bz.borzoi_dnase_new))

    # ---------------- the reproduction check, FIRST ---------------------------
    # The enlarged endpoint-1 Borzoi column mixes two runs. Whether they may be merged into one
    # column is decided here, before any endpoint statistic is read.
    rp = os.path.join(scored, "scoring", "borzoi", "borzoi_repro_scores.tsv")
    if os.path.exists(rp):
        r = pd.read_csv(rp, sep="\t")
        r["borzoi_atac_repro"] = r.borzoi_atac_delta * r.orient_mult
        r["borzoi_dnase_repro"] = r.borzoi_dnase_delta * r.orient_mult
        r = r.merge(dirfeat[["key", "borzoi_atac_archived", "borzoi_dnase_archived",
                             "borzoi_orient_mult"]],
                    left_on="canon", right_on="key", how="left")
        r = r[r.borzoi_atac_repro.notna() & r.borzoi_atac_archived.notna()]
        qc["borzoi_repro_run_n"] = int(len(r))
        if len(r) >= 10:
            for ch, new_c, old_c in (("atac", "borzoi_atac_repro", "borzoi_atac_archived"),
                                     ("dnase", "borzoi_dnase_repro", "borzoi_dnase_archived")):
                x = r[new_c].values.astype(float)
                y = r[old_c].values.astype(float)
                dd = np.abs(x - y)
                qc[f"borzoi_repro_run_pearson_{ch}"] = _pearson(x, y)
                qc[f"borzoi_repro_run_spearman_{ch}"] = signed_spearman(x, y)
                qc[f"borzoi_repro_run_max_abs_diff_{ch}"] = float(np.nanmax(dd))
                qc[f"borzoi_repro_run_median_abs_diff_{ch}"] = float(np.nanmedian(dd))
                qc[f"borzoi_repro_run_median_abs_archived_{ch}"] = float(np.nanmedian(np.abs(y)))
                qc[f"borzoi_repro_run_sign_agreement_{ch}"] = float(
                    (np.sign(x) == np.sign(y)).mean())
            qc["borzoi_repro_orient_mult_agreement"] = float(
                (r.orient_mult.astype(float) == r.borzoi_orient_mult.astype(float)).mean())
            r.to_csv(os.path.join(tables, "borzoi_reproduction_vs_archive.tsv"),
                     sep="\t", index=False)
        # merge rule, decided from the check and recorded either way
        tol = 1e-6
        ok = (qc.get("borzoi_repro_run_n", 0) >= 10
              and qc.get("borzoi_repro_run_max_abs_diff_atac", np.inf) < tol
              and qc.get("borzoi_repro_run_max_abs_diff_dnase", np.inf) < tol)
        qc["borzoi_merge_rule"] = (
            "MERGED_single_column" if ok else "SPLIT_separate_provenance_columns")
        qc["borzoi_merge_rule_tolerance"] = tol
        log(f"  Borzoi reproduction: n={qc.get('borzoi_repro_run_n')} "
            f"max|d| ATAC {qc.get('borzoi_repro_run_max_abs_diff_atac')} "
            f"DNase {qc.get('borzoi_repro_run_max_abs_diff_dnase')} -> "
            f"{qc['borzoi_merge_rule']}")
    else:
        qc["borzoi_repro_run_n"] = 0
        qc["borzoi_merge_rule"] = "SPLIT_separate_provenance_columns_no_reproduction_run"

    merged_ok = qc["borzoi_merge_rule"] == "MERGED_single_column"

    # ---------------- endpoint 1 ---------------------------------------------
    lab1 = cbp_a[cbp_a.label == 1].copy()
    lab1 = lab1[["key", "chr_x", "pos_hg38", "ref", "alt", "beta_alt", "q_value", "peak_id",
                 "heldout_fold", "chrombpnet_adult_hep", "model_peak_overlap"]]
    lab1 = lab1.merge(cbp_b[["key", "chrombpnet_adult_hep_gse281367"]], on="key", how="left")
    lab1 = lab1.merge(dirfeat[["key", "chrombpnet_hepg2", "borzoi_atac_archived",
                               "borzoi_dnase_archived", "ag_atac_src79", "ag_dnase_src79",
                               "effect", "feature_id"]], on="key", how="left")
    lab1 = lab1.merge(atlas, on="key", how="left")
    lab1["borzoi_atac_new"] = lab1.key.map(bz_atac)
    lab1["borzoi_dnase_new"] = lab1.key.map(bz_dnase)
    if merged_ok:
        lab1["borzoi_atac"] = lab1.borzoi_atac_archived.combine_first(lab1.borzoi_atac_new)
        lab1["borzoi_dnase"] = lab1.borzoi_dnase_archived.combine_first(lab1.borzoi_dnase_new)
    else:
        # the reproduction failed: the new run is its own provenance, the archived column stays
        # the one the wave-1 matched set used, and no row mixes the two.
        lab1["borzoi_atac"] = lab1.borzoi_atac_archived
        lab1["borzoi_dnase"] = lab1.borzoi_dnase_archived
    lab1["borzoi_provenance"] = np.where(
        lab1.borzoi_atac_archived.notna(), "archived_src79",
        np.where(lab1.borzoi_atac_new.notna() & merged_ok, "this_program_gpu", "none"))
    lab1["block_1mb"] = [block_of(c, p) for c, p in zip(lab1.chr_x, lab1.pos_hg38)]
    lab1["control_allele_identity"] = [
        NT_RANK.get(str(a).upper(), np.nan) - NT_RANK.get(str(r).upper(), np.nan)
        for r, a in zip(lab1.ref, lab1.alt)]
    lab1["control_position"] = (lab1.pos_hg38 % 1000) / 1000.0 - 0.5

    qc["endpoint1_leads"] = int(len(lab1))
    qc["endpoint1_borzoi_coverage_archived"] = int(lab1.borzoi_atac_archived.notna().sum())
    qc["endpoint1_borzoi_coverage_new_run"] = int(lab1.borzoi_atac_new.notna().sum())
    qc["endpoint1_borzoi_coverage_merged"] = int(lab1.borzoi_atac.notna().sum())
    ov = lab1[lab1.borzoi_atac_archived.notna() & lab1.borzoi_atac_new.notna()]
    qc["borzoi_incidental_overlap_n"] = int(len(ov))

    e1_tierA = [m for m in MODELS_E1_TIERA if int(lab1[m].notna().sum()) >= 10]
    e1_tierA4 = [m for m in MODELS_E1_TIERA4 if int(lab1[m].notna().sum()) >= 10]
    qc["endpoint1_models_absent_tierA4"] = [m for m in MODELS_E1_TIERA4 if m not in e1_tierA4]
    tierA = lab1.dropna(subset=e1_tierA + ["beta_alt"]).copy()
    tierA4 = lab1.dropna(subset=e1_tierA4 + ["beta_alt"]).copy()
    tierA_wave1 = lab1[lab1.borzoi_atac_archived.notna()].dropna(
        subset=e1_tierA + ["beta_alt"]).copy()
    tierA4_wave1 = lab1[lab1.borzoi_atac_archived.notna()].dropna(
        subset=e1_tierA4 + ["beta_alt"]).copy()
    for nm, d in (("tierA", tierA), ("tierA4", tierA4),
                  ("tierA_wave1", tierA_wave1), ("tierA4_wave1", tierA4_wave1)):
        log(f"  {nm} n={len(d)} blocks={d.block_1mb.nunique()} "
            f"provenance={d.borzoi_provenance.value_counts().to_dict()}")
        qc[f"endpoint1_{nm}_provenance"] = d.borzoi_provenance.value_counts().to_dict()

    # ---------------- endpoint 2 ---------------------------------------------
    ase = pd.read_csv(os.path.join(scored, "scoring", "ase_union_sites.tsv"), sep="\t")
    ase["block_1mb"] = [block_of(c, p) for c, p in zip(ase.chrom, ase.pos)]
    ase["precision_weight"] = 1.0 / ase.se_log2
    t0tab = load_track0_ase(ase.uid.tolist())
    ase = ase.merge(t0tab.rename(columns={"variant_uid": "uid"}), on="uid", how="left")

    cbp_new = {}
    for fam, tag in (("hepatocyte_5fold_v2", "chrombpnet_adult_hep"),
                     ("gse281367_hepatocyte_5fold_v2", "chrombpnet_adult_hep_gse281367")):
        tab, info = load_new_chrombpnet(scored, fam, tag)
        cbp_new[tag] = info
        ase = ase.merge(tab, on="uid", how="left")
    qc["endpoint2_chrombpnet_new"] = cbp_new

    # endpoint 2 Borzoi is single-provenance: every one of these sites was scored by this
    # program's GPU run, none by the archived panel.
    ase["borzoi_atac"] = ase.key.map(bz_atac)
    ase["borzoi_dnase"] = ase.key.map(bz_dnase)
    ase = ase.merge(dirfeat[["key", "chrombpnet_hepg2", "borzoi_atac_archived",
                            "borzoi_dnase_archived"]], on="key", how="left")
    ase = ase.merge(cbp_a[["key", "label", "beta_alt"]].rename(
        columns={"label": "in_currin_positive_set", "beta_alt": "currin_beta_alt"}),
        on="key", how="left")
    ase["control_allele_identity"] = [
        NT_RANK.get(str(a).upper(), np.nan) - NT_RANK.get(str(r).upper(), np.nan)
        for r, a in zip(ase.ref, ase.alt)]
    ase["control_position"] = (ase.pos % 1000) / 1000.0 - 0.5
    cov = {m: int(ase[m].notna().sum()) for m in
           ["alphagenome_atac_liver", "alphagenome_dnase_liver", "chrombpnet_adult_hep",
            "chrombpnet_adult_hep_gse281367", "borzoi_atac", "borzoi_dnase",
            "chrombpnet_hepg2", "borzoi_atac_archived"]}
    qc["endpoint2_coverage_after_gpu"] = cov
    log(f"  endpoint 2 coverage after GPU scoring: {cov}")
    ase.to_csv(os.path.join(tables, "endpoint2_ase_union_scores.tsv"), sep="\t", index=False)

    match_e2 = [m for m in MATCH_E2 if int(ase[m].notna().sum()) >= 10]
    models_e2 = [m for m in MODELS_E2 if m in match_e2 or m.startswith("control_")]
    dropped_e2 = [m for m in MATCH_E2 if m not in match_e2]
    qc["endpoint2_models_absent"] = dropped_e2
    if dropped_e2:
        log(f"  endpoint 2 models with no usable scores, dropped and NOT imputed: {dropped_e2}")
    ase_m = ase.dropna(subset=match_e2 + ["mean_log2_alt_over_ref", "se_log2"]).copy()
    log(f"  endpoint 2 matched n={len(ase_m)} blocks={ase_m.block_1mb.nunique()}")
    ase_all = ase.dropna(subset=["mean_log2_alt_over_ref", "se_log2"]).copy()
    # the three-family set the superseded run used, so the enlargement is one component
    e2_three = [m for m in match_e2 if not m.startswith("borzoi_")]
    ase_m3 = ase.dropna(subset=e2_three + ["mean_log2_alt_over_ref", "se_log2"]).copy()

    # ---------------- run -----------------------------------------------------
    rows, boot = [], {}
    tierA4 = run_endpoint(tierA4, e1_tierA4, "beta_alt", None,
                          "endpoint1_currin_tierA4", rows, boot)
    tierA = run_endpoint(tierA, e1_tierA, "beta_alt", None,
                         "endpoint1_currin_tierA", rows, boot)
    tierA4_wave1 = run_endpoint(tierA4_wave1, e1_tierA4, "beta_alt", None,
                                "endpoint1_currin_tierA4_wave1_subset", rows, boot)
    ase_m = run_endpoint(ase_m, models_e2, "mean_log2_alt_over_ref", "precision_weight",
                         "endpoint2_allelic_imbalance", rows, boot)
    ase_m3 = run_endpoint(ase_m3, [m for m in models_e2 if not m.startswith("borzoi_")],
                          "mean_log2_alt_over_ref", "precision_weight",
                          "endpoint2_three_family_495", rows, boot)
    ase_all = run_endpoint(ase_all, ["control_allele_identity", "control_position"],
                           "mean_log2_alt_over_ref", "precision_weight",
                           "endpoint2_all500_controls", rows, boot)

    tierA4.to_csv(os.path.join(tables, "endpoint1_matched_tierA4.tsv.gz"), sep="\t",
                  index=False, compression="gzip")
    tierA.to_csv(os.path.join(tables, "endpoint1_matched_tierA.tsv.gz"), sep="\t",
                 index=False, compression="gzip")
    ase_m.to_csv(os.path.join(tables, "endpoint2_matched.tsv"), sep="\t", index=False)

    metrics = pd.DataFrame(rows)
    metrics.to_csv(os.path.join(tables, "c1e2_model_metrics.tsv"), sep="\t", index=False)

    plook = {}
    for r in metrics.itertuples():
        if r.metric == "signed_spearman_weighted":
            plook[(r.endpoint, r.model, "weighted")] = r.value
        elif r.metric == "signed_spearman_unweighted":
            plook[(r.endpoint, r.model, "unweighted")] = r.value

    e1_pairs = [("chrombpnet_adult_hep", "alphagenome_atac_liver"),
                ("chrombpnet_adult_hep", "alphagenome_dnase_liver"),
                ("chrombpnet_adult_hep", "borzoi_atac"),
                ("chrombpnet_adult_hep", "borzoi_dnase"),
                ("chrombpnet_adult_hep", "chrombpnet_hepg2"),
                ("chrombpnet_adult_hep", "control_allele_identity"),
                ("borzoi_atac", "alphagenome_atac_liver"),
                ("borzoi_atac", "control_allele_identity"),
                ("alphagenome_atac_liver", "control_allele_identity"),
                ("chrombpnet_adult_hep_gse281367", "alphagenome_atac_liver")]
    e2_pairs = [("chrombpnet_adult_hep", "alphagenome_atac_liver"),
                ("chrombpnet_adult_hep", "alphagenome_dnase_liver"),
                ("chrombpnet_adult_hep", "borzoi_atac"),
                ("chrombpnet_adult_hep", "borzoi_dnase"),
                ("borzoi_atac", "alphagenome_atac_liver"),
                ("borzoi_atac", "borzoi_dnase"),
                ("chrombpnet_adult_hep", "control_allele_identity"),
                ("borzoi_atac", "control_allele_identity"),
                ("alphagenome_atac_liver", "control_allele_identity"),
                ("chrombpnet_adult_hep_gse281367", "alphagenome_atac_liver")]
    crows = []
    for ep in ("endpoint1_currin_tierA4", "endpoint1_currin_tierA4_wave1_subset"):
        contrasts(boot, ep, "unweighted",
                  [p for p in e1_pairs if "chrombpnet_hepg2" not in p], crows, plook)
    contrasts(boot, "endpoint1_currin_tierA", "unweighted", e1_pairs, crows, plook)
    for kind in ("weighted", "unweighted"):
        contrasts(boot, "endpoint2_allelic_imbalance", kind, e2_pairs, crows, plook)
        contrasts(boot, "endpoint2_three_family_495", kind,
                  [p for p in e2_pairs if "borzoi_atac" not in p and "borzoi_dnase" not in p],
                  crows, plook)
    con = pd.DataFrame(crows)
    con.to_csv(os.path.join(tables, "c1e2_contrasts.tsv"), sep="\t", index=False)
    log("contrasts\n" + con.to_string())

    # ---------------- Holm over the confirmatory family of three --------------
    def getp(ep, contrast, kind):
        s = con[(con.endpoint == ep) & (con.contrast == contrast) & (con.weighting == kind)]
        return float(s.boot_p_two_sided.iloc[0]) if len(s) else np.nan

    with open(E2_RECEIPT) as fh:
        e2rec = json.load(fh)
    p3 = float(e2rec["primary_gain_bootstrap_p"])
    qc["e2_transfer_raw_p_source"] = E2_RECEIPT
    qc["e2_transfer_raw_p"] = p3

    prim = "chrombpnet_adult_hep - alphagenome_atac_liver"
    p1 = getp("endpoint1_currin_tierA4", prim, "unweighted")
    p2 = getp("endpoint2_allelic_imbalance", prim, "weighted")

    def holm(members):
        fam = sorted(members, key=lambda t: t[1])
        m_family, res, prev = len(fam), [], 0.0
        for i, (k, p) in enumerate(fam):
            adj = max(prev, min(1.0, p * (m_family - i)))
            prev = adj
            res.append({"member": k, "raw_p": p, "holm_step_threshold": 0.05 / (m_family - i),
                        "holm_adjusted_p": adj, REJ: bool(adj <= 0.05)})
        return res

    members = [("C1_endpoint1_chrombpnet_vs_atlas_atac_tierA4", p1),
               ("C1_endpoint2_chrombpnet_vs_atlas_atac_weighted", p2),
               ("E2_transfer_contrast", p3)]
    hd = pd.DataFrame(holm(members))
    hd.to_csv(os.path.join(tables, "c1e2_holm_family.tsv"), sep="\t", index=False)
    log("Holm family of three\n" + hd.to_string())

    p2u = getp("endpoint2_allelic_imbalance", prim, "unweighted")
    sens = {
        "endpoint1_tierA4_raw_p": p1,
        "endpoint1_tierA_raw_p": getp("endpoint1_currin_tierA", prim, "unweighted"),
        "endpoint1_tierA4_wave1_subset_raw_p": getp(
            "endpoint1_currin_tierA4_wave1_subset", prim, "unweighted"),
        "endpoint2_weighted_raw_p": p2,
        "endpoint2_unweighted_raw_p": p2u,
        "endpoint2_three_family_495_weighted_raw_p": getp(
            "endpoint2_three_family_495", prim, "weighted"),
        "e2_transfer_raw_p": p3,
        "holm_with_weighted_endpoint2": holm(members),
        "holm_with_unweighted_endpoint2": holm(
            [members[0], ("C1_endpoint2_chrombpnet_vs_atlas_atac_unweighted", p2u), members[2]]),
    }
    with open(os.path.join(out, "holm_sensitivity.json"), "w") as fh:
        json.dump(sens, fh, indent=2, default=float)

    # ---------------- counts --------------------------------------------------
    counts = []
    for nm, d, ms in (("endpoint1_currin_tierA4", tierA4, e1_tierA4),
                      ("endpoint1_currin_tierA", tierA, e1_tierA),
                      ("endpoint1_currin_tierA4_wave1_subset", tierA4_wave1, e1_tierA4),
                      ("endpoint1_currin_tierA_wave1_subset", tierA_wave1, e1_tierA),
                      ("endpoint2_allelic_imbalance", ase_m, models_e2),
                      ("endpoint2_three_family_495", ase_m3,
                       [m for m in models_e2 if not m.startswith("borzoi_")])):
        counts.append(dict(matched_set=nm, n_variants=len(d),
                           n_blocks_1mb=int(d.block_1mb.nunique()), models=",".join(ms)))
    counts.append(dict(matched_set="endpoint1_label_universe_currin_leads_with_beta",
                       n_variants=int(len(lab1)),
                       n_blocks_1mb=int(lab1.block_1mb.nunique()), models="label_only"))
    counts.append(dict(matched_set="endpoint2_union_sites", n_variants=int(len(ase)),
                       n_blocks_1mb=int(ase.block_1mb.nunique()), models="label_only"))
    for m in ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367", "chrombpnet_hepg2",
              "borzoi_atac", "borzoi_dnase", "borzoi_atac_archived", "borzoi_atac_new",
              "alphagenome_atac_liver", "alphagenome_dnase_liver"]:
        counts.append(dict(matched_set=f"endpoint1_coverage:{m}",
                           n_variants=int(lab1[m].notna().sum()),
                           n_blocks_1mb=int(lab1.loc[lab1[m].notna(), "block_1mb"].nunique()),
                           models=m))
    for m, n in cov.items():
        counts.append(dict(matched_set=f"endpoint2_coverage:{m}", n_variants=n,
                           n_blocks_1mb=int(ase.loc[ase[m].notna(), "block_1mb"].nunique()),
                           models=m))
    pd.DataFrame(counts).to_csv(os.path.join(tables, "c1e2_matched_counts.tsv"),
                                sep="\t", index=False)

    # ---------------- loader reproduction against the deposit -----------------
    with open(CALIB_JSON) as fh:
        calib = json.load(fh)
    full = lab1.dropna(subset=["chrombpnet_adult_hep", "beta_alt"])
    qc["loader_reproduction_full_currin_spearman"] = signed_spearman(
        full.chrombpnet_adult_hep.values.astype(float), full.beta_alt.values.astype(float))
    qc["loader_reproduction_full_currin_sign_concordance"] = sign_concordance(
        full.chrombpnet_adult_hep.values.astype(float), full.beta_alt.values.astype(float))
    qc["deposited_calibration_keys"] = {
        k: v for k, v in calib.items()
        if isinstance(v, (int, float))
        and ("spearman" in k.lower() or "concord" in k.lower())}

    dep = cbp_a[["key", "chrombpnet_adult_hep", "heldout_fold"]].rename(
        columns={"chrombpnet_adult_hep": "cbp_deposited", "heldout_fold": "fold_deposited"})
    rep = ase.merge(dep, on="key", how="inner")
    rep = rep[rep.chrombpnet_adult_hep.notna() & rep.cbp_deposited.notna()]
    qc["chrombpnet_reproduction_n"] = int(len(rep))
    if len(rep) >= 10:
        diff = np.abs(rep.chrombpnet_adult_hep.values - rep.cbp_deposited.values)
        qc["chrombpnet_reproduction_max_abs_diff"] = float(np.nanmax(diff))
        qc["chrombpnet_reproduction_median_abs_diff"] = float(np.nanmedian(diff))
        qc["chrombpnet_reproduction_pearson"] = _pearson(
            rep.chrombpnet_adult_hep.values.astype(float),
            rep.cbp_deposited.values.astype(float))
        qc["chrombpnet_reproduction_fold_agreement"] = float(
            (rep.chrombpnet_adult_hep_heldout_fold == rep.fold_deposited).mean())
        rep[["key", "uid", "chrombpnet_adult_hep", "cbp_deposited",
             "chrombpnet_adult_hep_heldout_fold", "fold_deposited"]].to_csv(
            os.path.join(tables, "chrombpnet_reproduction_vs_deposit.tsv"),
            sep="\t", index=False)

    # ---------------- prediction receipt, computed not transcribed ------------
    def val(ep, model, metric):
        s = metrics[(metrics.endpoint == ep) & (metrics.model == model)
                    & (metrics.metric == metric)]
        return float(s.value.iloc[0]) if len(s) else np.nan

    def crow(ep, contrast, kind):
        s = con[(con.endpoint == ep) & (con.contrast == contrast) & (con.weighting == kind)]
        return s.iloc[0] if len(s) else None

    E2M = "endpoint2_allelic_imbalance"
    E1A4 = "endpoint1_currin_tierA4"
    preds = []
    n_e2 = int(len(ase_m))
    preds.append(dict(
        id="G1", prediction="four-model endpoint-2 matched set reaches >=450 of 500 sites",
        observed=f"{n_e2} sites, {int(ase_m.block_1mb.nunique())} blocks, "
                 f"families={','.join(match_e2)}",
        met=bool(n_e2 >= 450 and len(dropped_e2) == 0)))
    g2w = val(E2M, "control_allele_identity", "signed_spearman_weighted")
    g2u = val(E2M, "control_allele_identity", "signed_spearman_unweighted")
    g2lo = val(E2M, "control_allele_identity", "signed_spearman_weighted_boot_ci_lo")
    g2hi = val(E2M, "control_allele_identity", "signed_spearman_weighted_boot_ci_hi")
    preds.append(dict(
        id="G2",
        prediction="allele-identity control positive weighted Spearman, interval excludes 0, "
                   "unweighted smaller",
        observed=f"weighted {g2w:.4f} [{g2lo:.4f}, {g2hi:.4f}]; unweighted {g2u:.4f}",
        met=bool(g2w > 0 and g2lo > 0 and g2u < g2w)))
    c3w = crow(E2M, prim, "weighted")
    c3u = crow(E2M, prim, "unweighted")
    preds.append(dict(
        id="G3", prediction="endpoint-2 primary contrast keeps its negative sign",
        observed=f"weighted {c3w.point:.4f} [{c3w.ci_lo:.4f}, {c3w.ci_hi:.4f}]; "
                 f"unweighted {c3u.point:.4f} [{c3u.ci_lo:.4f}, {c3u.ci_hi:.4f}]",
        met=bool(c3w.point < 0 and c3u.point < 0)))
    g4e2 = val(E2M, "borzoi_atac", "signed_spearman_unweighted")
    g4e1 = val("endpoint1_currin_tierA", "borzoi_atac", "signed_spearman_unweighted")
    preds.append(dict(
        id="G4", prediction="Borzoi ATAC endpoint-2 Spearman below its endpoint-1 tier-A value",
        observed=f"endpoint 2 {g4e2:.4f} vs endpoint 1 tier A {g4e1:.4f}",
        met=bool(np.isfinite(g4e2) and np.isfinite(g4e1) and g4e2 < g4e1)))
    c5 = crow(E1A4, prim, "unweighted")
    preds.append(dict(
        id="G5", prediction="tier A4 reaches >=3,800 leads and the C1.1 contrast sign survives",
        observed=f"tier A4 {len(tierA4)} leads in {int(tierA4.block_1mb.nunique())} blocks; "
                 f"contrast {c5.point:.4f} [{c5.ci_lo:.4f}, {c5.ci_hi:.4f}]",
        met=bool(len(tierA4) >= 3800 and c5.point < 0 and c5.excludes_zero)))
    g6lo = val(E2M, "borzoi_atac", "signed_spearman_weighted_boot_ci_lo")
    g6hi = val(E2M, "borzoi_atac", "signed_spearman_weighted_boot_ci_hi")
    g6 = val(E2M, "borzoi_atac", "signed_spearman_weighted")
    preds.append(dict(
        id="G6", prediction="Borzoi ATAC endpoint-2 weighted interval excludes 0",
        observed=f"{g6:.4f} [{g6lo:.4f}, {g6hi:.4f}]",
        met=bool(np.isfinite(g6lo) and (g6lo > 0 or g6hi < 0))))
    preds.append(dict(
        id="C1.1",
        prediction="ChromBPNet beats Atlas ATAC at Currin leads by 0.00-0.10, interval excludes 0",
        observed=f"tier A4 {c5.point:.4f} [{c5.ci_lo:.4f}, {c5.ci_hi:.4f}]",
        met=bool(0.0 <= c5.point <= 0.10 and c5.excludes_zero)))
    seq_models = [m for m in match_e2]
    c12_first, c12_detail = True, []
    for m in seq_models:
        v2 = val(E2M, m, "signed_spearman_unweighted")
        v1 = val("endpoint1_currin_tierA4", m, "signed_spearman_unweighted")
        c12_detail.append(f"{m} {v2:.4f} vs {v1:.4f}")
        if not (np.isfinite(v2) and np.isfinite(v1) and v2 < v1):
            c12_first = False
    lows = [val(E2M, m, "signed_spearman_weighted_boot_ci_lo") for m in seq_models]
    his = [val(E2M, m, "signed_spearman_weighted_boot_ci_hi") for m in seq_models]
    c12_second = any((lo < 0 < hi) for lo, hi in zip(lows, his))
    preds.append(dict(
        id="C1.2",
        prediction="at imbalance sites every model below its caQTL value and >=1 interval "
                   "includes 0",
        observed="; ".join(c12_detail) + f"; lowest weighted CI lo {min(lows):.4f}",
        met=bool(c12_first and c12_second)))
    e1c = val(E1A4, "control_allele_identity", "signed_spearman_unweighted")
    preds.append(dict(
        id="C1.3", prediction="allele-identity control at chance on both endpoints |rho|<0.05",
        observed=f"endpoint 1 tier A4 {e1c:.4f}; endpoint 2 weighted {g2w:.4f}, "
                 f"unweighted {g2u:.4f}",
        met=bool(abs(e1c) < 0.05 and abs(g2w) < 0.05 and abs(g2u) < 0.05)))
    pd.DataFrame(preds).to_csv(os.path.join(tables, "c1e2_predictions.tsv"),
                               sep="\t", index=False)
    log("predictions\n" + pd.DataFrame(preds).to_string())

    with open(os.path.join(out, "qc_checks.json"), "w") as fh:
        json.dump(qc, fh, indent=2, default=float)
    log("QC " + json.dumps(qc, indent=1, default=float))

    man_paths = [CBP_ADULT, CBP_281367, DIRFEAT, CALIB_JSON, E2_RECEIPT,
                 os.path.join(TRACK0, "ATAC.parquet"), os.path.join(TRACK0, "DNASE.parquet"),
                 os.path.join(scored, "scoring", "ase_union_sites.tsv"),
                 os.path.join(scored, "scoring", "borzoi", "borzoi_scores.tsv"),
                 os.path.join(scored, "scoring", "borzoi", "borzoi_repro_scores.tsv"),
                 os.path.join(scored, "scoring", "borzoi", "borzoi_targets.tsv"),
                 os.path.join(scored, "scoring", "borzoi", "borzoi_repro_targets.tsv")]
    for fam in ("hepatocyte_5fold_v2", "gse281367_hepatocyte_5fold_v2"):
        for f in range(5):
            man_paths.append(os.path.join(scored, "scoring", "chrombpnet", fam, "scores",
                                          f"fold{f}.variant_scores.tsv"))
    man = [dict(path=p, bytes=os.path.getsize(p), sha256=sha256(p))
           for p in man_paths if os.path.exists(p)]
    man += [dict(path=p, bytes=-1, sha256="absent_not_produced")
            for p in man_paths if not os.path.exists(p)]
    pd.DataFrame(man).to_csv(os.path.join(out, "MANIFEST.tsv"), sep="\t", index=False)
    log("done")


if __name__ == "__main__":
    main()
