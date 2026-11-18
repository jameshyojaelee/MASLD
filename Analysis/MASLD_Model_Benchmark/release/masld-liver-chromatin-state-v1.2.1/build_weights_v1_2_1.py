#!/usr/bin/env python3
"""Build weights/chromatin_state_v1_2_1.npz with form-specific prof_reliable_<form> masks.

Preserves every v1.1 and v1.2 array bitwise (except adding prof_reliable_<form> across
all 5 profile forms and updating version='1.2.1').
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
V12_W = BENCH / "release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
STABLE_OOF = BENCH / "executions/chromatin-stable-rrr-forms-20260908T192024Z/out/stable_oof.npz"
CIS_OOF = BENCH / "executions/chromatin-cis-head-20260908T102817Z/out/cis_head_oof.npz"
EXACT_OOF = BENCH / "executions/chromatin-v12-exact-shipped-recipe-20260908T163444Z/out/fixed_config_oof.npz"
DISC_OUT = BENCH / "executions/chromatin-discrepancy-map-v12form-20260908T222353Z/out"
LANE0 = BENCH / "executions/chromatin-beyond-histology-20260907T161521Z"
sys.path.insert(0, str(LANE0))
import common as C0  # noqa: E402


def _conc_block(raw):
    lib = raw.sum(1)
    ncol = raw.shape[1]
    Pm = raw / lib[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(Pm > 0, Pm * np.log(Pm), 0.0), axis=1)
    sa = np.sort(raw, axis=1)
    gini = ((2 * np.arange(1, ncol + 1) - ncol - 1) * sa).sum(1) / (ncol * sa.sum(1))
    iqr = np.array([float(np.subtract(*np.percentile(np.log2(r[r > 0]), [75, 25]))) for r in raw])
    kt = int(round(0.05 * ncol))
    top5 = (-np.sort(-raw, axis=1))[:, :kt].sum(1) / lib
    return np.column_stack([gini, ent, iqr, top5])


def main():
    W12 = np.load(V12_W, allow_pickle=True)
    ZS = np.load(STABLE_OOF, allow_pickle=True)
    ZC = np.load(CIS_OOF, allow_pickle=True)
    ZF = np.load(EXACT_OOF, allow_pickle=True)

    Cres = ZS["Cres_oof"].astype(np.float64)
    fold = ZS["fold"]

    # Reload raw substrate exactly as run_disc.py to compute training-fold mean baseline
    C0.OUT = HERE / "weights"
    S = C0.load_substrate()
    assert np.array_equal(S["fold"], fold), "fold order mismatch"
    _Craw = C0.logcpm(S["h3_raw"])
    _CONC = _conc_block(S["h3_raw"])

    form_oof = {
        "states_cis": ZC["pred_nested"].astype(np.float64),
        "rrr": ZS["pred_rrr_wide_stable"].astype(np.float64),
        "rrr_cis": ZS["pred_rrr_cis_feature"].astype(np.float64),
        "rrr_offset_cis": ZS["pred_rrr_offset_cis_wide"].astype(np.float64),
        "rrr_offset_states_cis": ZS["pred_rrr_offset_union_wide"].astype(np.float64),
    }

    pfs_by_form = {k: [] for k in form_oof}
    for f in sorted(set(fold.tolist())):
        tr, te = np.where(fold != f)[0], np.where(fold == f)[0]
        Zi = np.column_stack([np.ones(len(_Craw)), _CONC])
        b, *_ = np.linalg.lstsq(Zi[tr], _Craw[tr], rcond=None)
        Cf = _Craw - Zi @ b
        assert np.allclose(Cf[te], Cres[te], atol=1e-6, rtol=0), f"fold {f}: Cres mismatch"
        den = ((Cf[te] - Cf[tr].mean(0)[None, :]) ** 2).sum(0)
        den_safe = np.where(den < 1e-12, np.nan, den)
        for form, PR in form_oof.items():
            num = ((Cres[te] - PR[te]) ** 2).sum(0)
            pfs_by_form[form].append(1.0 - num / den_safe)

    reliable_masks = {}
    counts = {}
    for form, folds_list in pfs_by_form.items():
        PFS = np.asarray(folds_list, dtype=np.float64)
        mask = np.all(np.isfinite(PFS) & (PFS > 0.2), axis=0)
        reliable_masks[form] = mask
        counts[form] = int(mask.sum())
        print(f"  form {form:22s}: {counts[form]:5d} reliable regions (all 5 folds skill > 0.2)", flush=True)

    assert np.array_equal(reliable_masks["states_cis"], np.asarray(W12["prof_reliable"], bool)), "states_cis mismatch vs v1.2 prof_reliable"
    assert np.array_equal(reliable_masks["states_cis"], np.asarray(ZC["reliable"], bool)), "states_cis mismatch vs cis_head_oof reliable"
    assert np.array_equal(reliable_masks["rrr_offset_cis"], np.asarray(W12["prof_reliable_offset_form"], bool)), "rrr_offset_cis mismatch vs v1.2 prof_reliable_offset_form"
    assert counts["states_cis"] == 4847
    assert counts["rrr_offset_cis"] == 4364

    # Diagnostic count for fixed-100kb exact_shipped_recipe OOF as well
    pfs_fx = []
    PR_fx = ZF["pred_fixed"].astype(np.float64)
    for f in sorted(set(fold.tolist())):
        tr, te = np.where(fold != f)[0], np.where(fold == f)[0]
        Zi = np.column_stack([np.ones(len(_Craw)), _CONC])
        b, *_ = np.linalg.lstsq(Zi[tr], _Craw[tr], rcond=None)
        Cf = _Craw - Zi @ b
        den = ((Cf[te] - Cf[tr].mean(0)[None, :]) ** 2).sum(0)
        den_safe = np.where(den < 1e-12, np.nan, den)
        num = ((Cres[te] - PR_fx[te]) ** 2).sum(0)
        pfs_fx.append(1.0 - num / den_safe)
    n_fx = int(np.all(np.isfinite(pfs_fx) & (np.asarray(pfs_fx) > 0.2), axis=0).sum())
    print(f"  diagnostic fixed_100kb states_cis: {n_fx} reliable regions", flush=True)

    disc_npz = DISC_OUT / "shipped_form_region_skills.npz"
    if disc_npz.exists():
        dz = np.load(disc_npz, allow_pickle=True)
        assert np.array_equal(reliable_masks["rrr_offset_cis"], np.asarray(dz["reliable_offset_form"], bool))
        assert np.array_equal(reliable_masks["states_cis"], np.asarray(dz["reliable_states_cis"], bool))

    payload = {k: W12[k] for k in W12.files}
    payload["version"] = np.array("1.2.1")
    for form, mask in reliable_masks.items():
        payload[f"prof_reliable_{form}"] = mask

    out_npz = HERE / "weights/chromatin_state_v1_2_1.npz"
    np.savez_compressed(out_npz, **payload)
    (HERE / "weights/reliable_masks_summary.json").write_text(
        json.dumps({"counts": counts, "diagnostic_fixed_100kb_states_cis_count": n_fx, "rule": "all 5 outer folds per-region OOF skill > 0.2"}, indent=2) + "\n"
    )
    print(f"Wrote {out_npz} cleanly.", flush=True)


if __name__ == "__main__":
    main()
