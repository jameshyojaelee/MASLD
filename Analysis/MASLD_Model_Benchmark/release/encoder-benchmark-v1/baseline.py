#!/usr/bin/env python3
"""
The baselines a new encoder has to beat.

    python baseline.py --counts my_counts.npz --row-contract my_rows.tsv --out preds/

Two families, both computed on the caller's own data with no pretraining:

  hvg_pca_<head>     the 50-dimensional PCA of highly variable genes. NOT a
                     strawman. On the reference substrate it is one of only two
                     blocks with all seven held-out studies clean, which makes
                     it the most trustworthy row in the benchmark, and with a
                     two-layer head it is a dead tie with a liver-specific
                     encoder trained on 102 donors.

  library_shape_<head>   ten library-complexity descriptors and nothing else:
                     entropy, Gini, the top-N count shares, detection rate,
                     library size. No gene identity at all. This exists because
                     library-complexity descriptors explained apparent
                     biological signal in three separate analyses in this
                     project, and a benchmark that does not report its own
                     library-shape floor cannot tell you whether an encoder is
                     reading biology or sequencing depth.

Everything is fitted on the outer-training studies of each fold and applied to
the held-out study, so the outputs drop straight into evaluate.py.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from metrics import ROSTER  # noqa: E402

N_HVG = 2000
N_PC = 50
SEED = 20260824

DESCRIPTOR_NAMES = (
    "log10_library_size", "n_detected_genes", "frac_zero_genes",
    "frac_counts_top10", "frac_counts_top50", "frac_counts_top500",
    "shannon_entropy", "gini", "median_lognorm_detected", "iqr_lognorm_detected",
)


def library_shape_descriptors(X) -> tuple[np.ndarray, list[str]]:
    """Ten descriptors of the shape of a cell's count vector. No gene identity.

    Computed row by row so a sparse matrix never has to be densified. If the
    matrix carries no library size -- because it arrived already normalised --
    log10_library_size is reported as unavailable rather than silently faked.
    """
    import scipy.sparse as sp
    csr = sp.csr_matrix(X)
    n, g = csr.shape
    out = np.zeros((n, len(DESCRIPTOR_NAMES)))
    lib = np.asarray(csr.sum(axis=1)).ravel()
    normalised = bool(np.allclose(lib, lib[0], rtol=1e-3)) if n > 1 else False
    for i in range(n):
        s, e = csr.indptr[i], csr.indptr[i + 1]
        v = csr.data[s:e].astype(np.float64)
        v = v[v > 0]
        if v.size == 0:
            continue
        tot = v.sum()
        p = v / tot
        srt = -np.sort(-v)
        cum = np.cumsum(srt)
        sa = np.sort(v)
        k = sa.size
        gini = ((2 * np.arange(1, k + 1) - k - 1) * sa).sum() / (k * sa.sum())
        lg = np.log1p(v)
        out[i] = (
            np.log10(tot) if tot > 0 else 0.0,
            float(v.size),
            1.0 - v.size / g,
            cum[min(9, k - 1)] / tot,
            cum[min(49, k - 1)] / tot,
            cum[min(499, k - 1)] / tot,
            float(-(p * np.log(p)).sum()),
            float(gini),
            float(np.median(lg)),
            float(np.subtract(*np.percentile(lg, [75, 25]))),
        )
    names = list(DESCRIPTOR_NAMES)
    if normalised:
        # every library sums to the same total, so depth is not in this matrix
        out[:, 0] = 0.0
        names[0] = "log10_library_size__UNAVAILABLE_matrix_is_already_normalised"
    return out, names


def hvg_pca(train_X, apply_X, n_hvg=N_HVG, n_pc=N_PC):
    """Select HVGs and fit the PCA on TRAINING rows only, then project both."""
    import scipy.sparse as sp
    tr = sp.csr_matrix(train_X)
    lib = np.asarray(tr.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    trn = sp.diags(1e4 / lib) @ tr
    trn.data = np.log1p(trn.data)
    mean = np.asarray(trn.mean(axis=0)).ravel()
    sq = np.asarray(trn.multiply(trn).mean(axis=0)).ravel()
    var = np.maximum(sq - mean ** 2, 0.0)
    hvg = np.argsort(-var)[:n_hvg]
    A = np.asarray(trn[:, hvg].todense(), dtype=np.float64)
    mu = A.mean(axis=0)
    A -= mu
    # Eigendecompose the n_hvg x n_hvg covariance rather than SVD-ing the
    # n_cells x n_hvg matrix. Identical principal axes, and 2000^3 beats
    # 40000 x 2000^2 by two orders of magnitude.
    cov = (A.T @ A) / max(A.shape[0] - 1, 1)
    vals, vecs = np.linalg.eigh(cov)
    V = vecs[:, np.argsort(-vals)[:n_pc]]

    def project(M):
        m = sp.csr_matrix(M)
        l = np.asarray(m.sum(axis=1)).ravel()
        l[l == 0] = 1.0
        mn = sp.diags(1e4 / l) @ m
        mn.data = np.log1p(mn.data)
        return (np.asarray(mn[:, hvg].todense()) - mu) @ V

    return project(train_X), project(apply_X)


# ---------------------------------------------------------------- the head
# The benchmark's whole design is "several feature blocks through ONE shared
# head", so a baseline scored by a different head is not measuring the same
# thing. This is the head from the run that produced the reference results:
# scripts/fit_predict_common_cell_heads_study_50000.py, same architecture, same
# optimiser, same epoch count, same donor-class weights.
#
# ⛔ It is NOT bit-identical to the frozen shards. Those were trained on CUDA
# under torch.use_deterministic_algorithms(True); this runs wherever you are.
# MODEL_CARD.md reports the measured agreement. Same recipe, not same bits.
BATCH_SIZE = 1024
FIXED_EPOCHS = 30
LEARNING_RATE = 0.001
WEIGHT_DECAY = 0.01


def donor_class_weights(donors, labels) -> np.ndarray:
    """Equal mass per class, then equal mass per donor within a class."""
    counts, donors_by_class = {}, {}
    for d, y in zip(donors, labels):
        counts[(str(d), int(y))] = counts.get((str(d), int(y)), 0) + 1
        donors_by_class.setdefault(int(y), set()).add(str(d))
    k = len(donors_by_class)
    w = np.asarray([1.0 / (k * len(donors_by_class[int(y)]) * counts[(str(d), int(y))])
                    for d, y in zip(donors, labels)], dtype=np.float32)
    return w / w.sum() * len(w)


def build_head(head_id: str, width: int):
    from torch import nn
    if head_id == "linear":
        return nn.Linear(width, len(ROSTER))
    if head_id == "two_layer_mlp":
        hidden = min(256, width)
        return nn.Sequential(nn.Linear(width, hidden), nn.GELU(), nn.Dropout(0.1),
                             nn.Linear(hidden, len(ROSTER)))
    raise ValueError(f"unknown head {head_id!r}")


def fit_head(Xtr, ytr, Xte, head, donors_tr, seed=SEED):
    import torch
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(max(1, int(__import__("os").environ.get("SLURM_CPUS_PER_TASK", "1"))))
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    A = Xtr.astype(np.float32, copy=False)
    B = Xte.astype(np.float32, copy=False)
    mean = A.mean(axis=0, dtype=np.float64).astype(np.float32)
    scale = A.std(axis=0, dtype=np.float64).astype(np.float32)
    scale[scale == 0.0] = 1.0
    A, B = (A - mean) / scale, (B - mean) / scale
    w = donor_class_weights(donors_tr, ytr)

    model = build_head(head, A.shape[1]).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    x = torch.from_numpy(A).to(dev)
    y = torch.from_numpy(ytr.astype(np.int64, copy=False)).to(dev)
    ww = torch.from_numpy(w).to(dev)
    gen = torch.Generator().manual_seed(seed)
    for _ in range(FIXED_EPOCHS):
        model.train()
        order = torch.randperm(len(x), generator=gen)
        for st in range(0, len(x), BATCH_SIZE):
            b = order[st:st + BATCH_SIZE].to(dev)
            opt.zero_grad(set_to_none=True)
            loss = torch.nn.functional.cross_entropy(model(x[b]), y[b], reduction="none")
            ((loss * ww[b]).sum() / ww[b].sum()).backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        logits = model(torch.from_numpy(B).to(dev)).float().cpu().numpy()
    e = np.exp(logits - logits.max(axis=1, keepdims=True))
    return e / e.sum(axis=1, keepdims=True)


def read_contract(path):
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh]
    return rows


def run(counts_path, contract_path, out, blocks):
    import scipy.sparse as sp
    out = pathlib.Path(out); out.mkdir(parents=True, exist_ok=True)
    rows = read_contract(contract_path)
    z = np.load(counts_path, allow_pickle=False)
    if "format" in z.files:                      # scipy sparse npz
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    else:
        X = sp.csr_matrix(z[z.files[0]])
    if X.shape[0] != len(rows):
        raise ValueError(f"counts have {X.shape[0]} rows, row contract has {len(rows)}")
    y = np.asarray([ROSTER.index(r["true_class"]) for r in rows])
    fold = np.asarray([int(r["outer_fold"]) for r in rows])
    row_ids = np.asarray([r["row_id"] for r in rows])
    donor = np.asarray([r["donor_id"] for r in rows])

    feature_sets = {}
    _pca_cache: dict = {}
    if "library_shape" in blocks:
        D, names = library_shape_descriptors(X)
        feature_sets["library_shape"] = ("precomputed", D)
        print(f"  library-shape descriptors: {names}")
    preds = {}
    for block in blocks:
        for head in ("linear", "two_layer_mlp"):
            P = np.zeros((len(rows), len(ROSTER)))
            for f in sorted(set(fold.tolist())):
                te = np.flatnonzero(fold == f)
                tr = np.flatnonzero(fold != f)
                if block == "hvg_pca":
                    if f not in _pca_cache:
                        _pca_cache[f] = hvg_pca(X[tr], X[te])
                    Ftr, Fte = _pca_cache[f]
                else:
                    D = feature_sets["library_shape"][1]
                    Ftr, Fte = D[tr], D[te]
                P[te] = fit_head(Ftr, y[tr], Fte, head, donor[tr])
                print(f"    {block}::{head} fold {f}  n_train {len(tr)} n_test {len(te)}",
                      flush=True)
            preds[f"{block}::{head}"] = P.astype(np.float32)
    np.savez_compressed(out / "baseline_predictions.npz", row_id=row_ids, **preds)
    (out / "baseline_receipt.json").write_text(json.dumps({
        "blocks": sorted(preds), "n_hvg": N_HVG, "n_pc": N_PC, "seed": SEED,
        "descriptor_panel": list(DESCRIPTOR_NAMES),
        "fitted_on": "outer-training rows of each fold only",
        "head": {"architecture": "same as the reference run: linear, or "
                          "Linear->GELU->Dropout(0.1)->Linear with hidden=min(256,width)",
                 "optimizer": "AdamW", "lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY,
                 "epochs": FIXED_EPOCHS, "batch_size": BATCH_SIZE,
                 "sample_weights": "donor-class balanced",
                 "bitwise_identical_to_the_frozen_shards": False,
                 "why": ("the frozen shards were trained on CUDA under "
                         "torch.use_deterministic_algorithms(True). Same recipe, not "
                         "same bits. MODEL_CARD.md reports the measured agreement.")},
    }, indent=2) + "\n")
    print(f"wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="scipy sparse npz, cells x genes")
    ap.add_argument("--row-contract", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--blocks", default="hvg_pca,library_shape")
    a = ap.parse_args()
    run(a.counts, a.row_contract, a.out, [b for b in a.blocks.split(",") if b])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
