#!/usr/bin/env python3
"""AG2X Stage C scoring (PRESPEC_AG2X_stageC.md): AlphaGenome 2 SC_ATAC under cell conditioning.

Order, with gates (a fixed threshold that fails writes stageC/NOT_TESTABLE.json and stops before more calls):
  1. encode_cells, one request per donor set (3 calls) -> control (a).
  2. pilot: 10 tier 1 members, set1_otherhep (10 calls) -> SC_ATAC tracks returned and usable.
  3. otherhep scoring, sets 1-3 -> control (b) against the Stage B reference conditions on disk.
  4. ccl20pos scoring, sets 1-3.
Hard cap: STAGEC_CAP hosted calls (encode requests, first sends and every resend), read from the
persistent ledger stageC/raw/call_ledger.jsonl before every send; resumable (chunks written once, skipped).
Requests carry GRCh38 intervals, alleles and anonymous count matrices (obs index c0.., Ensembl var IDs).
The key is loaded by ag2x_lib.create_v2() and never printed.
Usage: c10_score_stageC.py [--stop-after controls|otherhep|all]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402
import c01_stageC_lib as S  # noqa: E402


def ledger(event: str, n: int, label: str):
    if S.spent() + n > S.STAGEC_CAP:
        raise SystemExit(f"Stage C hosted-call cap: {S.spent()} spent, next {n} ({label}), cap {S.STAGEC_CAP}; stopping")
    with open(S.LEDGER, "a") as fh:
        fh.write(json.dumps({"utc": L.utc(), "event": event, "n": n, "label": label}) + "\n")


def send(fn, n: int, label: str):
    """First send and every resend pass the cap check and are written to the ledger before they go out."""
    ledger("send", n, label)
    transient = 0
    while True:
        try:
            return fn()
        except BaseException as e:  # noqa: BLE001
            if L.is_quota_error(e):
                L.log(f"quota error on {label}: sleeping {L.QUOTA_SLEEP:.0f}s")
                time.sleep(L.QUOTA_SLEEP)
            elif L.is_transient_error(e) and transient < L.MAX_TRANSIENT:
                transient += 1
                wait = min(L.TRANSIENT_SLEEP * 2 ** (transient - 1), L.TRANSIENT_SLEEP_CAP)
                L.log(f"transient error {transient}/{L.MAX_TRANSIENT} on {label} ({e.code().name}): sleeping {wait:.0f}s")
                time.sleep(wait)
            else:
                raise
            ledger("resend", n, label)


def load_cells(s: str, state: str):
    """Anonymous AnnData for the request: integer counts, Ensembl var IDs, obs index c0..; nothing else."""
    import anndata as ad
    from scipy import sparse
    from alphagenome.models.v2 import dna_client, gene_counts_condition
    a = ad.read_h5ad(S.CELLS / f"{s}_{state}.h5ad")
    X = sparse.csr_matrix(a.X).astype(np.int32)
    assert 0 < X.shape[0] <= dna_client.MAX_CONDITIONING_CELLS
    assert (np.diff(X.indptr) > 0).all() and (X.data > 0).all()
    out = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])]),
                     var=pd.DataFrame(index=pd.Index(a.var_names.astype(str))))
    nbytes = gene_counts_condition.from_anndata(out).ByteSize()
    return out, nbytes


def encode_and_control_a(client, cells):
    import anndata as ad
    rows = []
    for s in S.SETS:
        p = S.RAW / "encode" / f"{s}.npz"
        p.parent.mkdir(parents=True, exist_ok=True)
        labels = np.concatenate([[st] * cells[(s, st)][0].n_obs for st in S.STATES])
        if not p.exists():
            A = ad.concat([cells[(s, st)][0] for st in S.STATES], join="outer", fill_value=0)
            A.obs_names = [f"c{i}" for i in range(A.n_obs)]
            assert A.n_obs <= 1000  # one EncodeCells request (MAX_ENCODE_CELLS_PER_REQUEST)
            E = send(lambda: client.encode_cells(A, progress_bar=False, max_workers=1), 1, f"encode {s}")
            np.savez(p.with_suffix(".tmp.npz"), E=E, labels=labels)
            os.replace(p.with_suffix(".tmp.npz"), p)
        z = np.load(p, allow_pickle=False)
        E, lab = z["E"], z["labels"]
        hep = np.isin(lab, ["ccl20pos", "otherhep"])
        y_hs = np.where(hep, "hep", "stellate")
        m = {st: E[lab == st].mean(0) for st in S.STATES}
        hm = E[hep].mean(0)
        rows.append({"set": s, "n_ccl20pos": int((lab == "ccl20pos").sum()), "n_otherhep": int((lab == "otherhep").sum()),
                     "n_stellate": int((lab == "stellate").sum()), "latent_dim": int(E.shape[1]),
                     "bal_acc_hep_vs_stellate": S.loo_nearest_centroid_bal_acc(E, y_hs),
                     "bal_acc_ccl20pos_vs_otherhep": S.loo_nearest_centroid_bal_acc(E[hep], lab[hep]),
                     "cosdist_mean_hep_vs_stellate": S.cosine_dist(hm, m["stellate"]),
                     "cosdist_mean_ccl20pos_vs_otherhep": S.cosine_dist(m["ccl20pos"], m["otherhep"])})
    t = pd.DataFrame(rows)
    t["pass"] = t.bal_acc_hep_vs_stellate >= S.CTRL_A_MIN_BAL_ACC
    t.to_csv(S.TAB / "control_a_encode_cells.tsv", sep="\t", index=False)
    L.log("control (a):\n" + t.round(4).to_string())
    return bool(t["pass"].all())


def score_cells(client, uids, adata, name, out_dir, chunk=50, workers=4):
    from alphagenome.models.v2 import dna_client, variant_scorers as vs
    out_dir.mkdir(parents=True, exist_ok=True)
    scorer = [vs.RECOMMENDED_VARIANT_SCORERS["SC_ATAC"]]
    for i in range(0, len(uids), chunk):
        part = uids[i:i + chunk]
        p = out_dir / f"chunk_{i // chunk:05d}.parquet"
        if p.exists():
            continue
        variants = [L.to_variant(u) for u in part]
        intervals = [v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB) for v in variants]
        t0, call_t = time.time(), L.utc()
        res = send(lambda: client.score_variants(intervals=intervals, variants=variants, variant_scorers=scorer,
                                                 cell_conditionings=[adata] * len(part), cell_condition_names=[name] * len(part),
                                                 progress_bar=False, max_workers=workers),
                   len(part), f"{out_dir.name} chunk {i // chunk}")
        el = time.time() - t0
        tidy = vs.tidy_scores(res, match_gene_strand=True, include_extended_metadata=True)
        tidy = (tidy if tidy is not None else pd.DataFrame()).copy()
        for col in tidy.columns:
            if tidy[col].dtype == object:
                tidy[col] = tidy[col].astype(str)
        tidy = tidy.assign(conditioning=name, predicted=True, imputed=tidy.get("imputation_mode", pd.Series(dtype=str)) != "observed",
                           candidate_not_adopted=True, client_version="v2", package_version=L.package_version(),
                           model_version="ALL_FOLDS (requested; no version string returned)", call_utc=call_t)
        tmp = p.with_suffix(".tmp")
        tidy.to_parquet(tmp, index=False)
        os.replace(tmp, p)
        with open(out_dir / "chunk_log.jsonl", "a") as fh:
            fh.write(json.dumps({"chunk": i // chunk, "n_variants": len(part), "elapsed_s": el, "call_utc": call_t,
                                 "rows": int(len(tidy)), "conditioning": name}) + "\n")
        L.log(f"{out_dir.name} chunk {i // chunk}: {len(part)} in {el:.1f}s ({len(part) / max(el, 1e-9):.2f}/s); spent {S.spent()}")


def not_testable(reason: str, **kw):
    json.dump({"utc": L.utc(), "stageC": "not testable", "reason": reason, **kw}, open(S.C / "NOT_TESTABLE.json", "w"), indent=1)
    L.log(f"NOT TESTABLE: {reason}")
    raise SystemExit(0)


def pilot_check(tidy: pd.DataFrame, name: str):
    sc = tidy[tidy.output_type == "SC_ATAC"]
    tracks = sc.groupby(["track_name", "imputation_mode"]).size().reset_index(name="rows")
    tracks["protocol"] = tracks.track_name.map(lambda n: S.protocol_of(n, name))
    tracks.to_csv(S.TAB / "pilot_sc_atac_tracks.tsv", sep="\t", index=False)
    usable = tracks[tracks.protocol.notna() & tracks.imputation_mode.isin(S.USABLE_MODES)]
    L.log("pilot tracks:\n" + tracks.to_string())
    return usable.protocol.nunique(), sorted(tidy.output_type.unique().tolist())


def main():
    stop = sys.argv[sys.argv.index("--stop-after") + 1] if "--stop-after" in sys.argv else "all"
    assert stop in ("controls", "otherhep", "all")
    for d in (S.RAW, S.TAB):
        d.mkdir(parents=True, exist_ok=True)
    if (S.C / "NOT_TESTABLE.json").exists():
        raise SystemExit("NOT_TESTABLE.json exists; Stage C stopped earlier")
    v = pd.read_csv(S.C / "inputs/variants_stageC.tsv", sep="\t")
    bad = L.assert_ref(v.uid.tolist())
    assert not bad, f"REF mismatch: {bad[:3]}"
    uids = sorted(v.uid)
    pilot = sorted(v.loc[v.pilot, "uid"])
    n_plan = 3 + len(pilot) + len(uids) * len(S.SETS) * len(S.SCORED_STATES)
    cells = {(s, st): load_cells(s, st) for s in S.SETS for st in S.STATES}
    plan = {"variants": len(uids), "members": int((v.stageC_role == "member").sum()),
            "background": int((v.stageC_role == "background").sum()), "planned_calls": n_plan, "cap": S.STAGEC_CAP,
            "spent_before_this_run": S.spent(), "request_bytes_cell_conditioning": {f"{s}_{st}": cells[(s, st)][1] for s, st in cells},
            "n_cells": {f"{s}_{st}": int(cells[(s, st)][0].n_obs) for s, st in cells}}
    json.dump(plan, open(S.RAW / "scoring_plan.json", "w"), indent=1)
    L.log(json.dumps(plan))
    assert n_plan <= S.STAGEC_CAP, "planned calls exceed the Stage C cap"

    client = L.create_v2()
    # 1. control (a)
    if not encode_and_control_a(client, cells):
        not_testable("control (a): encode_cells does not separate hepatocyte from stellate sets",
                     table="stageC/tables/control_a_encode_cells.tsv")
    # 2. pilot
    name = S.cond_name("set1", "otherhep")
    score_cells(client, pilot, cells[("set1", "otherhep")][0], name, S.RAW / "pilot_set1_otherhep", chunk=len(pilot))
    n_proto, otypes = pilot_check(S.read_dir(S.RAW / "pilot_set1_otherhep"), name)
    if n_proto == 0:
        not_testable("no usable SC_ATAC track (observed or imputed_supported) of the reference protocols under cell conditioning",
                     output_types=otypes, table="stageC/tables/pilot_sc_atac_tracks.tsv")
    if stop == "controls":
        L.log("stop after controls (a) and pilot")
        return
    # 3. otherhep, then control (b)
    for s in S.SETS:
        nm = S.cond_name(s, "otherhep")
        score_cells(client, uids, cells[(s, "otherhep")][0], nm, S.RAW / nm)
    ref_h, ref_s = S.reference(S.REF_HEP), S.reference(S.REF_STELLATE)
    rows = []
    for s in S.SETS:
        nm = S.cond_name(s, "otherhep")
        r = S.readout(S.read_dir(S.RAW / nm), nm)
        j = pd.concat([r, ref_h, ref_s], axis=1, keys=["cond", "ref_hep", "ref_stellate"]).reindex(uids).dropna()
        rh, rs = spearmanr(j.cond, j.ref_hep)[0], spearmanr(j.cond, j.ref_stellate)[0]
        rows.append({"set": s, "n_variants": len(j), "rho_vs_ref_hepatocyte": rh, "rho_vs_ref_stellate": rs,
                     "pass_rho": rh >= S.CTRL_B_MIN_RHO, "pass_hep_gt_stellate": rh > rs})
    t = pd.DataFrame(rows)
    t["pass"] = t.pass_rho & t.pass_hep_gt_stellate
    t.to_csv(S.TAB / "control_b_reference.tsv", sep="\t", index=False)
    L.log("control (b):\n" + t.round(4).to_string())
    # determinism: pilot vs main run, same conditioning
    a = S.read_dir(S.RAW / "pilot_set1_otherhep")
    b = S.read_dir(S.RAW / "set1_otherhep")
    k = ["variant_id", "track_name"]
    dm = a[k + ["raw_score"]].merge(b[k + ["raw_score"]], on=k, suffixes=("_pilot", "_main"))
    det = {"rows": int(len(dm)), "identical_share": float((dm.raw_score_pilot == dm.raw_score_main).mean()),
           "max_abs_diff": float((dm.raw_score_pilot - dm.raw_score_main).abs().max())}
    json.dump(det, open(S.TAB / "determinism_pilot_vs_main.json", "w"), indent=1)
    L.log(f"determinism {det}")
    if not bool(t["pass"].all()):
        not_testable("control (b): conditioned other-hepatocyte SC_ATAC does not reproduce the reference hepatocyte condition",
                     table="stageC/tables/control_b_reference.tsv")
    if stop == "otherhep":
        L.log("stop after otherhep and control (b)")
        return
    # 4. ccl20pos
    for s in S.SETS:
        nm = S.cond_name(s, "ccl20pos")
        score_cells(client, uids, cells[(s, "ccl20pos")][0], nm, S.RAW / nm)
    summ = {"utc": L.utc(), "package_version": L.package_version(), "hosted_calls_stageC": S.spent(), "cap": S.STAGEC_CAP}
    json.dump(summ, open(S.RAW / "scoring_summary.json", "w"), indent=1)
    L.log(json.dumps(summ))


if __name__ == "__main__":
    main()
