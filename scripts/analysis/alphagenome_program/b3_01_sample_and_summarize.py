#!/usr/bin/env python3
"""B3 step 1: fix the sample from the archive snapshot, then stream one summary row per region.

Implements section 4 (B3) of scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md under the
hashed prespecification in <out>/prespec/PRESPECIFICATION.md.

Read-only with respect to GWAS/finemapping/results/alphagenome_atlas/. Never concatenates a
substitution-level matrix across regions: each region is opened, reduced to a handful of scalars, and
closed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd

ARCHIVE = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/"
           "alphagenome_atlas/p4-saturation-20260909T191917Z/raw/saturation/U1_h3k27ac")
COVARIATES = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/"
              "alphagenome_atlas/run-20260909T153939Z/tables/region_covariates.tsv.gz")

SCORERS = ["ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "AVI_SCORE"]
REQUIRED = [s + ".h5ad" for s in SCORERS] + ["request.json"]

# primary_liver ontology terms, matching scripts/analysis/alphagenome_atlas/lib_atlas.track_class.
# CL:0000182 (hepatocyte) is deliberately NOT here: in this archive it is embryonic in-vitro
# differentiated cells for DNASE/CHIP_HISTONE and an unstaged FANTOM primary cell for CAGE.
LIVER_TISSUE = {"UBERON:0002107", "UBERON:0001114", "UBERON:0001115"}

SEED = 20260914
N_PER_QUINTILE = 600
N_QUINTILES = 5
STALE_SECONDS = 600  # a directory touched more recently than this is treated as still being written


# --------------------------------------------------------------------------- frame and sample
def scan_completed(snapshot_path: str, scan_start: float) -> pd.DataFrame:
    names = [ln.strip() for ln in open(snapshot_path) if ln.strip()]
    rows = []
    for name in names:
        d = os.path.join(ARCHIVE, name)
        try:
            ents = {e.name: e.stat() for e in os.scandir(d)}
        except OSError as exc:
            rows.append({"dirname": name, "complete": False, "why": f"scandir:{exc.__class__.__name__}"})
            continue
        missing = [r for r in REQUIRED if r not in ents or ents[r].st_size == 0]
        if missing:
            rows.append({"dirname": name, "complete": False, "why": "missing:" + ",".join(missing)})
            continue
        newest = max(ents[r].st_mtime for r in REQUIRED)
        if newest > scan_start - STALE_SECONDS:
            rows.append({"dirname": name, "complete": False, "why": "mtime_within_%ds" % STALE_SECONDS})
            continue
        rows.append({"dirname": name, "complete": True, "why": "",
                     "bytes": sum(ents[r].st_size for r in REQUIRED)})
    return pd.DataFrame(rows)


def build_sample(out_dir: str, snapshot_path: str) -> pd.DataFrame:
    scan_start = time.time()
    scan = scan_completed(snapshot_path, scan_start)
    scan["region_key"] = scan.dirname.str.replace("_", ":", n=1)
    n_snap, n_done = len(scan), int(scan.complete.sum())
    print(f"[frame] snapshot {n_snap} dirs; completed {n_done}; incomplete {n_snap - n_done}", flush=True)
    if n_snap - n_done:
        print(scan.loc[~scan.complete, "why"].value_counts().to_string(), flush=True)

    cov = pd.read_csv(COVARIATES, sep="\t")
    done = scan.loc[scan.complete, ["dirname", "region_key", "bytes"]]
    frame = done.merge(cov, on="region_key", how="inner", validate="one_to_one")
    if len(frame) != n_done:
        raise SystemExit(f"[frame] {n_done - len(frame)} completed regions have no covariate row")

    frame = frame.sort_values("region_key", kind="mergesort").reset_index(drop=True)
    edges_frame = np.quantile(frame.skill_shipped_mean, np.linspace(0, 1, N_QUINTILES + 1))
    edges_full = np.quantile(cov.skill_shipped_mean, np.linspace(0, 1, N_QUINTILES + 1))
    pd.DataFrame({"q": range(N_QUINTILES + 1), "edge_completed_frame": edges_frame,
                  "edge_full_universe": edges_full}).to_csv(
        os.path.join(out_dir, "tables", "quintile_edges.tsv"), sep="\t", index=False)

    frame["skill_quintile"] = pd.qcut(frame.skill_shipped_mean, N_QUINTILES, labels=[1, 2, 3, 4, 5]).astype(int)
    rng = np.random.default_rng(SEED)
    picks = []
    for q in range(1, N_QUINTILES + 1):
        pool = frame.index[frame.skill_quintile == q].to_numpy()
        if len(pool) < N_PER_QUINTILE:
            raise SystemExit(f"[sample] quintile {q} has only {len(pool)} completed regions")
        picks.append(rng.choice(pool, size=N_PER_QUINTILE, replace=False))
    sample = frame.loc[np.sort(np.concatenate(picks))].reset_index(drop=True)

    keep = ["region_key", "dirname", "chrom", "start", "end", "width", "mb_block", "gc", "promoter",
            "dist_nearest_tss", "log10_dist_tss", "signal_mean", "signal_sd", "zero_fraction",
            "skill_shipped_mean", "skill_quintile", "increment_local_mean", "reliable_shipped_form",
            "bytes"]
    sample[keep].to_csv(os.path.join(out_dir, "prespec", "sample_regions.tsv"), sep="\t", index=False)

    # also record the completed frame composition so the chromosome truncation is on file
    comp = (cov.assign(on_disk=cov.region_key.isin(set(scan.loc[scan.complete, "region_key"])))
            .groupby("chrom").agg(total=("region_key", "size"), completed=("on_disk", "sum")))
    comp["fraction"] = comp.completed / comp.total
    comp.to_csv(os.path.join(out_dir, "tables", "frame_by_chromosome.tsv"), sep="\t")

    meta = {"snapshot": os.path.basename(snapshot_path), "dirs_in_snapshot": n_snap,
            "completed_regions": n_done, "seed": SEED, "n_per_quintile": N_PER_QUINTILE,
            "sampled": len(sample), "stale_seconds": STALE_SECONDS,
            "scan_started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(scan_start)),
            "sample_blocks": int(sample.mb_block.nunique()),
            "sample_bytes_estimate": int(sample.bytes.sum())}
    with open(os.path.join(out_dir, "prespec", "sample_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    h = hashlib.sha256(open(os.path.join(out_dir, "prespec", "sample_regions.tsv"), "rb").read()).hexdigest()
    with open(os.path.join(out_dir, "prespec", "sample_regions.sha256"), "w") as fh:
        fh.write(f"{h}  sample_regions.tsv\n")
    print(f"[sample] {len(sample)} regions, {meta['sample_blocks']} 1-Mb blocks, "
          f"{meta['sample_bytes_estimate']/1e9:.1f} GB on disk; sha256 {h}", flush=True)
    print("[sample] region list hashed BEFORE the first h5ad open", flush=True)
    return sample


# --------------------------------------------------------------------------- per-region reduction
def _cat(f, path):
    """Read an h5ad categorical or plain string column as a numpy array of str."""
    node = f[path]
    if isinstance(node, h5py.Group):
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]])
        codes = node["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "")
        return out.astype(str)
    arr = node[:]
    return np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in arr])


def _groups_for(f, scorer):
    """Return {group_name: (column indices, life_stage_unannotated)} from this file's own var table."""
    varkeys = set(f["var"].keys())
    n = f["X"].shape[1]
    if scorer == "AVI_SCORE":
        return {"avi_all": (np.arange(n), False)}
    onto = _cat(f, "var/ontology_curie")
    bname = _cat(f, "var/biosample_name")
    has_stage = "biosample_life_stage" in varkeys
    stage = _cat(f, "var/biosample_life_stage") if has_stage else np.full(n, "", dtype=object)
    liver = np.isin(onto, list(LIVER_TISSUE))
    adult = liver & (stage == "adult") if has_stage else liver
    hepg2 = np.char.find(np.char.lower(bname.astype(str)), "hepg2") >= 0
    return {"adult_liver": (np.flatnonzero(adult), not has_stage),
            "hepg2": (np.flatnonzero(hepg2), not has_stage)}


def _stats(a: np.ndarray, pos_inv: np.ndarray, n_pos: int):
    """mean, 95th percentile, and top-5%-of-positions concentration of a per-substitution vector."""
    if a.size == 0 or not np.isfinite(a).any():
        return math.nan, math.nan, math.nan
    mean_abs = float(np.mean(a))
    p95 = float(np.percentile(a, 95))
    w = np.bincount(pos_inv, weights=a, minlength=n_pos)
    tot = w.sum()
    if not np.isfinite(tot) or tot <= 0:
        return mean_abs, p95, math.nan
    k = max(1, int(math.ceil(0.05 * n_pos)))
    top = np.partition(w, n_pos - k)[n_pos - k:]
    return mean_abs, p95, float(top.sum() / tot)


def summarize_region(task):
    region_key, dirname, width = task
    out = []
    for scorer in SCORERS:
        path = os.path.join(ARCHIVE, dirname, scorer + ".h5ad")
        try:
            with h5py.File(path, "r") as f:
                variants = _cat(f, "obs/variant")
                pos = np.asarray([int(v.split(":")[1]) for v in variants], dtype=np.int64)
                upos, pos_inv = np.unique(pos, return_inverse=True)
                n_pos, n_sub = len(upos), len(pos)
                names = _cat(f, "var/name")
                for gname, (cols, unstaged) in _groups_for(f, scorer).items():
                    base = {"region_key": region_key, "scorer": scorer, "group": gname,
                            "n_tracks": int(len(cols)), "n_sub": n_sub, "n_pos": n_pos,
                            "width": width, "sub_eq_3w": bool(n_sub == 3 * width),
                            "life_stage_unannotated": bool(unstaged),
                            "track_set": "|".join(sorted(names[cols])) if len(cols) else ""}
                    if len(cols) == 0:
                        base.update(dict.fromkeys(
                            ["mean_abs", "p95_abs", "conc_top5", "mean_abs_q", "p95_abs_q", "conc_top5_q"],
                            math.nan))
                        out.append(base)
                        continue
                    ci = np.sort(cols)
                    a = np.abs(f["X"][:, ci]).mean(axis=1)
                    aq = np.abs(f["layers/quantiles"][:, ci]).mean(axis=1)
                    m, p, c = _stats(a, pos_inv, n_pos)
                    mq, pq, cq = _stats(aq, pos_inv, n_pos)
                    base.update({"mean_abs": m, "p95_abs": p, "conc_top5": c,
                                 "mean_abs_q": mq, "p95_abs_q": pq, "conc_top5_q": cq})
                    out.append(base)
        except Exception as exc:  # noqa: BLE001 - a bad region is recorded, never silently dropped
            out.append({"region_key": region_key, "scorer": scorer, "group": "ERROR",
                        "n_tracks": 0, "n_sub": 0, "n_pos": 0, "width": width, "sub_eq_3w": False,
                        "life_stage_unannotated": False, "track_set": f"{exc.__class__.__name__}: {exc}",
                        "mean_abs": math.nan, "p95_abs": math.nan, "conc_top5": math.nan,
                        "mean_abs_q": math.nan, "p95_abs_q": math.nan, "conc_top5_q": math.nan})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    sample = build_sample(args.out, args.snapshot)
    tasks = list(zip(sample.region_key, sample.dirname, sample.width.astype(int)))

    t0 = time.time()
    rows, done = [], 0
    with Pool(args.workers) as pool:
        for res in pool.imap_unordered(summarize_region, tasks, chunksize=8):
            rows.extend(res)
            done += 1
            if done % 250 == 0:
                print(f"[read] {done}/{len(tasks)} regions, {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    dest = os.path.join(args.out, "tables", "region_summaries.tsv.gz")
    df.to_csv(dest, sep="\t", index=False)
    print(f"[read] wrote {dest}: {len(df)} rows in {time.time()-t0:.0f}s", flush=True)
    print("[read] errors:", int((df.group == "ERROR").sum()), flush=True)
    print(df.groupby(["scorer", "group"]).agg(n=("n_tracks", "size"),
                                              tracks_modal=("n_tracks", lambda s: s.mode().iat[0]),
                                              tracks_min=("n_tracks", "min"),
                                              tracks_max=("n_tracks", "max"),
                                              nan_mean=("mean_abs", lambda s: int(s.isna().sum()))
                                              ).to_string(), flush=True)
    hetero = df.groupby(["scorer", "group"]).track_set.nunique()
    print("[read] distinct track sets per scorer x group:\n" + hetero.to_string(), flush=True)


if __name__ == "__main__":
    sys.exit(main())
