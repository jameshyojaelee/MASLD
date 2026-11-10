#!/usr/bin/env python3
"""A3 step 3: per-motif p < 5e-5 thresholds, validation against motifmatchr, and
allelic PWM scores (ref, alt) for every in-peak variant x 879 JASPAR2024 motifs.

Score = max log2 PWM score over all windows that cover the variant, both strands,
with the PWMs exactly as motifmatchr builds them (exported by step 2).
"""
import os
import sys

import numpy as np
import pandas as pd

A3 = os.environ.get("A3_OUT", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption")
WORK = f"{A3}/work"
RES = f"{A3}/results"
P_CUT = 5e-5
EPS = 1e-3
BASES = "ACGT"
CHUNK = 20000


def load_pwms():
    pw = pd.read_csv(f"{WORK}/motif_pwm_log2.tsv", sep="\t")
    meta = pd.read_csv(f"{WORK}/motif_meta.tsv", sep="\t")
    ids = meta.motif_id.tolist()
    mats = {}
    for mid, g in pw.groupby("motif_id", sort=False):
        g = g.sort_values("pos")
        mats[mid] = g[list(BASES)].to_numpy(dtype=np.float64)  # L x 4
    return ids, [mats[i] for i in ids], meta


def threshold(pwm, bg):
    """Smallest score s with P(S >= s) <= P_CUT under iid bg (scores discretised to EPS);
    the maximum score when no score reaches P_CUT."""
    ints = np.round(pwm / EPS).astype(np.int64)  # L x 4
    lo = ints.min(axis=1)
    dist = np.array([1.0])
    base = 0
    for j in range(ints.shape[0]):
        off = ints[j] - lo[j]
        new = np.zeros(dist.size + off.max())
        for b in range(4):
            new[off[b]:off[b] + dist.size] += bg[b] * dist
        dist = new
        base += lo[j]
    tail = np.cumsum(dist[::-1])[::-1]  # P(S >= base + k)
    ok = np.nonzero(tail <= P_CUT)[0]
    if ok.size == 0:
        # p unattainable (short motifs): motifmatchr/MOODS then calls only windows at the
        # maximum score (checked on the validation sample, PRESPEC amendment 2026-10-25)
        return pwm.max(axis=1).sum() - 1e-4
    return (base + ok[0]) * EPS


def onehot(seqs):
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), -1)
    out = np.zeros(arr.shape + (4,), dtype=np.float32)
    for k, ch in enumerate(BASES.encode()):
        out[..., k] = arr == ch
    return out  # N x W x 4; N/other bases score 0 in every column (rare in peaks)


def rc_pwm(pwm):
    return pwm[::-1, ::-1]  # reverse positions; A,C,G,T -> T,G,C,A


def group_by_length(pwms):
    groups = {}
    for i, p in enumerate(pwms):
        groups.setdefault(p.shape[0], []).append(i)
    return groups


def score_center(X, pwms, groups, center):
    """Max score over windows covering column `center`, both strands. X: N x W x 4."""
    N = X.shape[0]
    out = np.full((N, len(pwms)), -np.inf, dtype=np.float32)
    for L, idx in groups.items():
        sub = X[:, center - L + 1:center + L, :]  # N x (2L-1) x 4
        win = np.lib.stride_tricks.sliding_window_view(sub, (L, 4), axis=(1, 2))  # N x L x 1 x L x 4
        win = win.reshape(N * L, 4 * L)
        Wf = np.stack([pwms[i].astype(np.float32).reshape(-1) for i in idx], axis=1)
        Wr = np.stack([rc_pwm(pwms[i]).astype(np.float32).reshape(-1) for i in idx], axis=1)
        sf = (win @ Wf).reshape(N, L, len(idx)).max(axis=1)
        sr = (win @ Wr).reshape(N, L, len(idx)).max(axis=1)
        out[:, idx] = np.maximum(sf, sr)
    return out


def score_full(X, pwms, groups):
    """Max score over every window of the sequence, both strands (validation only)."""
    N, W, _ = X.shape
    out = np.full((N, len(pwms)), -np.inf, dtype=np.float32)
    for L, idx in groups.items():
        win = np.lib.stride_tricks.sliding_window_view(X, (L, 4), axis=(1, 2))
        nw = W - L + 1
        win = win.reshape(N * nw, 4 * L)
        Wf = np.stack([pwms[i].astype(np.float32).reshape(-1) for i in idx], axis=1)
        Wr = np.stack([rc_pwm(pwms[i]).astype(np.float32).reshape(-1) for i in idx], axis=1)
        sf = (win @ Wf).reshape(N, nw, len(idx)).max(axis=1)
        sr = (win @ Wr).reshape(N, nw, len(idx)).max(axis=1)
        out[:, idx] = np.maximum(sf, sr)
    return out


def main():
    ids, pwms, meta = load_pwms()
    bgt = pd.read_csv(f"{WORK}/bg.tsv", sep="\t").set_index("base").freq
    bg = np.array([bgt[b] for b in BASES])
    thr = np.array([threshold(p, bg) for p in pwms])
    meta["threshold_p5e-5"] = thr
    meta["max_score"] = [p.max(axis=1).sum() for p in pwms]
    meta["min_score"] = [p.min(axis=1).sum() for p in pwms]
    meta.to_csv(f"{WORK}/motif_meta_thresholds.tsv", sep="\t", index=False)
    print(f"thresholds: {(thr < meta.max_score.to_numpy() - 1e-7).sum()} below max score of {len(thr)}", flush=True)
    groups = group_by_length(pwms)

    # validation against motifmatchr on reference sequences
    vs = pd.read_csv(f"{WORK}/validation_seqs.tsv", sep="\t")
    vm = pd.read_csv(f"{WORK}/validation_motifmatchr_matches.tsv", sep="\t")
    vsc = pd.read_csv(f"{WORK}/validation_motifmatchr_scores.tsv", sep="\t")
    S = score_full(onehot(vs.seq.tolist()), pwms, groups)
    mine = S >= thr[None, :]
    col = {m: j for j, m in enumerate(ids)}
    mm = np.zeros_like(mine)
    mm[vm.row.to_numpy() - 1, vm.motif_id.map(col).to_numpy()] = True
    agree = (mine == mm).mean()
    d = np.abs(S[vsc.row.to_numpy() - 1, vsc.motif_id.map(col).to_numpy()] - vsc.score.to_numpy())
    per_motif = pd.DataFrame({"motif_id": ids, "n_mine": mine.sum(0), "n_motifmatchr": mm.sum(0),
                              "n_disagree": (mine != mm).sum(0)})
    per_motif.to_csv(f"{RES}/validation_motifmatchr_per_motif.tsv", sep="\t", index=False)
    val = pd.DataFrame([{
        "n_sequences": len(vs), "n_motifs": len(ids), "call_agreement": agree,
        "n_calls_mine": int(mine.sum()), "n_calls_motifmatchr": int(mm.sum()),
        "n_disagree": int((mine != mm).sum()), "max_abs_score_diff_on_matches": float(d.max()),
        "median_abs_score_diff_on_matches": float(np.median(d)),
    }])
    val.to_csv(f"{RES}/validation_motifmatchr_agreement.tsv", sep="\t", index=False)
    print(val.to_string(), flush=True)
    if agree < 0.99:
        sys.exit("call agreement with motifmatchr below 0.99; stop and amend")

    v = pd.read_csv(f"{WORK}/variants_hg38.tsv.gz", sep="\t", usecols=["key", "in_any_peak", "seq", "ref", "alt"],
                    dtype={"key": np.int64})
    v = v[v.in_any_peak].reset_index(drop=True)
    W = len(v.seq.iloc[0])
    center = (W - 1) // 2
    assert (v.seq.str[center] == v.ref).all()
    n = len(v)
    s_ref = np.empty((n, len(ids)), dtype=np.float32)
    s_alt = np.empty((n, len(ids)), dtype=np.float32)
    for st in range(0, n, CHUNK):
        en = min(n, st + CHUNK)
        seqs = v.seq.iloc[st:en].tolist()
        X = onehot(seqs)
        s_ref[st:en] = score_center(X, pwms, groups, center)
        alt_idx = np.array([BASES.index(a) for a in v.alt.iloc[st:en]])
        X[:, center, :] = 0
        X[np.arange(en - st), center, alt_idx] = 1
        s_alt[st:en] = score_center(X, pwms, groups, center)
        print(f"scored {en}/{n}", flush=True)
    np.save(f"{WORK}/s_ref.npy", s_ref)
    np.save(f"{WORK}/s_alt.npy", s_alt)
    v[["key"]].to_csv(f"{WORK}/scored_keys.tsv.gz", sep="\t", index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
