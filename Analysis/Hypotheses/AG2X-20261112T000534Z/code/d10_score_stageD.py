#!/usr/bin/env python3
"""Stage D hosted calls (PRESPEC_AG2X_stageD.md). Accessibility outputs only; no expression output is requested.

Order (so a cap stop drops the optional part last): for every variant (SORT1 control first)
  1. tracks: predict_variant, 1-Mb interval centred on the variant, one call per condition
     (bulk liver ATAC+DNASE, hepatocyte SC_ATAC, CD14 monocyte SC_ATAC); observed tracks kept, +/- 5 kb saved.
  2. ISM on the reference sequence: 64 bp (pos-31 .. pos+32), 3 substitutions per base, bulk liver
     (ATAC, DNASE scorers) and hepatocyte (SC_ATAC scorer); requests of 10 bp (30 substitutions).
  3. ISM on the ALT background (interval_variant = the variant): central 32 bp (pos-15 .. pos+16), hepatocyte.
Each request's result is written once and skipped on resume.
Ledger: stageD/raw/call_ledger.jsonl gets one line per send (first sends and resends) BEFORE the send; the
Stage D cap (6,000 substitution-or-prediction calls) is checked against the ledger before every send.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

D = L.PKG / "stageD"
RAW = D / "raw"
LEDGER = RAW / "call_ledger.jsonl"
CAP_D = 6000
PACKAGE_SPENT_BEFORE_D = 147_779  # with Stage B resends (RESULT_AG2X_stageB.md section 8)
PACKAGE_CAP = 170_000
HALF_KEEP = 5000
CONDS = {
    "bulk_liver": ("UBERON:0002107", ["ATAC", "DNASE"]),
    "hepatocyte": ("UBERON:0002107_CL:0000182", ["SC_ATAC"]),
    "monocyte_ctrl": ("UBERON:0000178_CL:0001054", ["SC_ATAC"]),
}
ISM_REF_CONDS = ["bulk_liver", "hepatocyte"]
ISM_ALT_CONDS = ["hepatocyte"]
ISM_REF_HALF = 32   # [pos-32, pos+32) 0-based = 64 bp; variant at offset 31
ISM_ALT_HALF = 16   # [pos-16, pos+16) 0-based = 32 bp; variant at offset 15
ISM_STEP = 10


def ledger_spent() -> int:
    if not LEDGER.exists():
        return 0
    return sum(json.loads(x)["n"] for x in open(LEDGER))


def ledger_add(kind: str, label: str, n: int):
    spent = ledger_spent()
    if spent + n > CAP_D:
        raise SystemExit(f"Stage D cap: {spent} used, next send {n}, cap {CAP_D}; stopping")
    if PACKAGE_SPENT_BEFORE_D + spent + n > PACKAGE_CAP:
        raise SystemExit("package cap would be exceeded; stopping")
    with open(LEDGER, "a") as fh:
        fh.write(json.dumps({"utc": L.utc(), "kind": kind, "label": label, "n": n}) + "\n")


def send(fn, label: str, n: int):
    """Ledger entry before every send; quota and transient errors resend (each resend is a new ledger entry)."""
    transient = 0
    kind = "first"
    while True:
        ledger_add(kind, label, n)
        try:
            return fn()
        except BaseException as e:  # noqa: BLE001
            if L.is_quota_error(e):
                L.log(f"quota error on {label}: sleeping {L.QUOTA_SLEEP:.0f}s")
                time.sleep(L.QUOTA_SLEEP)
                kind = "resend_quota"
                continue
            if L.is_transient_error(e) and transient < L.MAX_TRANSIENT:
                transient += 1
                wait = min(L.TRANSIENT_SLEEP * 2 ** (transient - 1), L.TRANSIENT_SLEEP_CAP)
                L.log(f"transient error {transient} on {label} ({e.code().name}): sleeping {wait:.0f}s")
                time.sleep(wait)
                kind = "resend_transient"
                continue
            raise


def plan(variants):
    rows = []
    for u in variants:
        for c in CONDS:
            rows.append(("tracks", u, c, None, 1))
    for u in variants:
        c0, p, r, a = L.parse_uid(u)
        for c in ISM_REF_CONDS:
            for s in range(p - ISM_REF_HALF, p + ISM_REF_HALF, ISM_STEP):
                e = min(s + ISM_STEP, p + ISM_REF_HALF)
                rows.append(("ism_ref", u, c, (s, e), 3 * (e - s)))
    for u in variants:
        c0, p, r, a = L.parse_uid(u)
        for c in ISM_ALT_CONDS:
            for s in range(p - ISM_ALT_HALF, p + ISM_ALT_HALF, ISM_STEP):
                e = min(s + ISM_STEP, p + ISM_ALT_HALF)
                rows.append(("ism_alt", u, c, (s, e), 3 * (e - s)))
    return pd.DataFrame(rows, columns=["step", "uid", "condition", "ism_span0", "n_calls"])


def safe(u):
    return u.replace(":", "_")


def do_tracks(client, u, cname):
    from alphagenome.models.v2 import dna_client
    from alphagenome.models.v2.dna_output import OutputType
    out = RAW / "tracks" / f"{safe(u)}__{cname}.npz"
    if out.exists():
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    cid, otypes = CONDS[cname]
    v = L.to_variant(u)
    req = [getattr(OutputType, o) for o in otypes]
    length = dna_client.SEQUENCE_LENGTH_1MB
    t0 = time.time()
    try:
        iv = v.reference_interval.resize(length)
        res = send(lambda: client.predict_variant(interval=iv, variant=v, requested_outputs=req, condition=cid),
                   f"tracks {u} {cname} 1MB", 1)
    except BaseException as e:  # noqa: BLE001
        if not L.is_message_size_error(e):
            raise
        L.log(f"message size error for {u} {cname} at 1 Mb; falling back to 131,072 bp (prespecified)")
        length = dna_client.SEQUENCE_LENGTH_100KB
        iv = v.reference_interval.resize(length)
        res = send(lambda: client.predict_variant(interval=iv, variant=v, requested_outputs=req, condition=cid),
                   f"tracks {u} {cname} 100KB", 1)
    c, p, r, a = L.parse_uid(u)
    save = {}
    meta = {"uid": u, "condition": cname, "condition_id": cid, "interval_start0": iv.start, "interval_end0": iv.end,
            "context_bp": length, "call_utc": L.utc(), "elapsed_s": time.time() - t0,
            "package_version": L.package_version(), "model_version": "ALL_FOLDS (requested; no version string returned)",
            "candidate_not_adopted": True, "outputs": {}}
    for o in otypes:
        tref = res.reference.get(getattr(OutputType, o))
        talt = res.alternate.get(getattr(OutputType, o))
        md = tref.metadata.reset_index(drop=True)
        obs = (md["imputation_mode"] == "observed").to_numpy() if "imputation_mode" in md else np.ones(len(md), bool)
        resl = int(tref.resolution)
        nbin = tref.values.shape[0]
        assert nbin * resl == iv.width, (nbin, resl, iv.width)
        vb = (p - 1 - iv.start) // resl
        lo, hi = max(0, vb - HALF_KEEP // resl), min(nbin, vb + HALF_KEEP // resl + 1)
        save[f"{o}_ref"] = tref.values[lo:hi][:, obs].astype(np.float32)
        save[f"{o}_alt"] = talt.values[lo:hi][:, obs].astype(np.float32)
        meta["outputs"][o] = {"resolution": resl, "slice_start0": iv.start + lo * resl, "variant_bin": vb - lo,
                              "tracks_observed": md.loc[obs, "name"].tolist() if "name" in md else [],
                              "tracks_all": md.get("name", pd.Series(dtype=str)).tolist(),
                              "imputation_modes_all": md.get("imputation_mode", pd.Series(dtype=str)).tolist()}
    tmp = out.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **save)
    tmp.rename(out)
    json.dump(meta, open(out.with_suffix(".json"), "w"), indent=1)
    L.log(f"tracks {u} {cname}: {time.time() - t0:.1f}s")


def scorer_of(adata, otypes):
    """score_ism_variants does not set uns['variant_scorer'] (score_variant does). Identify the output type of
    each AnnData from its track names: DNASE tracks are named '*DNase-seq*', bulk ATAC '*ATAC-seq-encode*'."""
    if len(otypes) == 1:
        return otypes[0]
    names = list(adata.var["name"])
    kinds = {"DNASE" if "DNase-seq" in n else "ATAC" if "ATAC-seq-encode" in n else "?" for n in names}
    if len(kinds) != 1 or "?" in kinds or kinds.pop() not in otypes:
        raise SystemExit(f"cannot identify scorer from track names {names}")
    return "DNASE" if "DNase-seq" in names[0] else "ATAC"


def label_ism_scorers(res, otypes, scorers):
    from alphagenome.models.v2 import variant_scorers as vs
    for per_variant in res:
        seen = []
        for ad in per_variant:
            o = scorer_of(ad, otypes)
            ad.uns["variant_scorer"] = vs.RECOMMENDED_VARIANT_SCORERS[o]
            seen.append(o)
        if sorted(seen) != sorted(otypes):
            raise SystemExit(f"scorer outputs per variant {seen}, expected {otypes}")


def do_ism(client, step, u, cname, span):
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers as vs
    s, e = span
    out = RAW / step / f"{safe(u)}__{cname}__{s}_{e}.parquet"
    if out.exists():
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    cid, otypes = CONDS[cname]
    scorers = [vs.RECOMMENDED_VARIANT_SCORERS[o] for o in otypes]
    v = L.to_variant(u)
    iv = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    ism_iv = genome.Interval(v.chromosome, s, e)
    iv_var = v if step == "ism_alt" else None
    t0 = time.time()
    call_t = L.utc()
    import pickle
    pk = out.with_suffix(".pkl")
    if pk.exists():  # response already received in an earlier run; parse it without a new send
        res = pickle.load(open(pk, "rb"))
    else:
        res = send(lambda: client.score_ism_variants(interval=iv, ism_interval=ism_iv, condition=cid,
                                                     variant_scorers=scorers, interval_variant=iv_var,
                                                     progress_bar=False, max_workers=1),
                   f"{step} {u} {cname} {s}-{e}", 3 * (e - s))
        with open(pk, "wb") as fh:  # raw response kept before parsing, so a parse error loses no call
            pickle.dump(res, fh)
    label_ism_scorers(res, otypes, scorers)
    tidy = vs.tidy_scores(res, match_gene_strand=True, include_extended_metadata=True)
    tidy = (tidy if tidy is not None else pd.DataFrame()).copy()
    for col in tidy.columns:
        if tidy[col].dtype == object:
            tidy[col] = tidy[col].astype(str)
    tidy["uid"] = u
    tidy["step"] = step
    tidy["condition_name"] = cname
    tidy["condition_id"] = cid
    tidy["ism_start0"] = s
    tidy["ism_end0"] = e
    tidy["package_version"] = L.package_version()
    tidy["model_version"] = "ALL_FOLDS (requested; no version string returned)"
    tidy["call_utc"] = call_t
    tidy["candidate_not_adopted"] = True
    tmp = out.with_suffix(".tmp")
    tidy.to_parquet(tmp, index=False)
    tmp.rename(out)
    L.log(f"{step} {u} {cname} {s}-{e}: {len(tidy)} rows, {time.time() - t0:.1f}s")


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(D / "inputs/stageD_variants.tsv", sep="\t")
    order = v.sort_values("role", kind="stable").uid.tolist()  # control_SORT1 sorts before selected
    bad = L.assert_ref(order)
    if bad:
        raise SystemExit(f"REF mismatch: {bad}")
    pl = plan(order)
    total = int(pl.n_calls.sum())
    pl.assign(ism_span0=pl.ism_span0.astype(str)).to_csv(RAW / "request_plan.tsv", sep="\t", index=False)
    L.log(json.dumps({"variants": order, "requests": len(pl), "planned_calls": total,
                      "by_step": pl.groupby("step").n_calls.sum().to_dict(), "ledger_spent": ledger_spent(), "cap": CAP_D}))
    if total > CAP_D:
        raise SystemExit("planned calls exceed the Stage D cap; stopping before any call")
    from alphagenome.models.v2 import variant_scorers as vs
    from alphagenome.models.v2.dna_output import OutputType
    for o in ["ATAC", "DNASE", "SC_ATAC"]:
        assert o in vs.RECOMMENDED_VARIANT_SCORERS and hasattr(OutputType, o)
    if "--plan-only" in sys.argv:
        L.log("plan only: no hosted call made")
        return
    client = L.create_v2()
    for _, r in pl.iterrows():
        if r.step == "tracks":
            do_tracks(client, r.uid, r.condition)
        else:
            do_ism(client, r.step, r.uid, r.condition, r.ism_span0)
    summ = {"utc": L.utc(), "package_version": L.package_version(), "stageD_calls_ledger": ledger_spent(),
            "package_total_with_resends": PACKAGE_SPENT_BEFORE_D + ledger_spent()}
    json.dump(summ, open(RAW / "scoring_summary.json", "w"), indent=1)
    L.log(json.dumps(summ))


if __name__ == "__main__":
    main()
