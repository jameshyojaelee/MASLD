#!/usr/bin/env python3
"""Stage 1 of the four-arm aspect decomposition.

Verifies the Stage 0 prespecification by digest BEFORE any matrix is opened,
then runs one instrument over four arms:

  Family H  harmonised pairwise directions, identical in form in every arm that
            deposits the labels. These are the only directions that are pooled
            across arms, and only at the level of the decision statistic.
  Family F  the full conditional decomposition, each component against all the
            others. Per arm only; the adjustment sets differ between arms
            because the deposited label sets differ, so pooling them would pool
            four different questions.
  Marginal  the channel-power check. An aspect whose nuisance-only direction is
            not itself UNIQUE cannot carry a negative.
  Controls  sex as a known-positive axis, an exact duplicate as a known-null,
            and a fully shuffled exposure as a known-null.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from aspect_axis_arms import (load_arm_a, load_arm_b, load_arm_c,  # noqa: E402
                              load_arm_d)
from aspect_axis_common import (BENCH, mde_spearman, rank_columns,  # noqa: E402
                                run_direction, sha256_file, stouffer,
                                tie_ceiling)
from masld_bench.artifacts import freeze_tree, verify_frozen_tree  # noqa: E402

COMPOSITE = "composite_nas"
FIBROSIS = "fibrosis"
METAB = "metabolism_steatosis"
INFLAM = "inflammation_lobular"
BALLOON = "ballooning"
NAMED = [METAB, INFLAM, FIBROSIS]


def design(arm, adjust, extra=None):
    cols = [np.ones(len(arm.participants))]
    names = []
    for a in adjust:
        cols.append(rankdata(arm.aspects[a].to_numpy(float)))
        names.append(a)
    for c in arm.nuisance.columns:
        cols.append(arm.nuisance[c].to_numpy(float))
        names.append(c)
    if extra is not None:
        for nm, v in extra.items():
            cols.append(np.asarray(v, float))
            names.append(nm)
    return np.column_stack(cols), names


def harmonised_plan(arm):
    have = set(arm.aspects.columns)
    plan = []
    def add(fam, x, z):
        if x in have and all(a in have for a in z):
            plan.append((fam, x, list(z)))
    add("H1", COMPOSITE, [FIBROSIS]); add("H1", FIBROSIS, [COMPOSITE])
    add("H2", METAB, [FIBROSIS]);     add("H2", FIBROSIS, [METAB])
    add("H3", INFLAM, [FIBROSIS]);    add("H3", FIBROSIS, [INFLAM])
    add("H4", METAB, [INFLAM]);       add("H4", INFLAM, [METAB])
    add("H5", BALLOON, [FIBROSIS])
    return plan


def full_plan(arm):
    comps = [a for a in arm.aspects.columns if a != COMPOSITE]
    if len(comps) < 2:
        return []
    return [("F", x, [a for a in comps if a != x]) for x in comps]


def evaluate_arm(arm, *, n_perm, seed, secondary_perm, verbose=True):
    rng = np.random.default_rng(seed)
    Rfeat = rank_columns(arm.features)
    n = len(arm.participants)
    family = Rfeat.shape[1]
    rows = []

    def record(fam, x, z, extra=None, tag=None, x_values=None, secondary=False):
        Cd, names = design(arm, z, extra)
        xr = rankdata(arm.aspects[x].to_numpy(float)) if x_values is None \
            else rankdata(np.asarray(x_values, float))
        raw = None
        builder = None
        if secondary:
            raw = (arm.aspects[x].to_numpy(float) if x_values is None
                   else np.asarray(x_values, float))
            def builder(xp, _z=z, _extra=extra):
                return design(arm, _z, _extra)[0]
        res = run_direction(Rfeat, xr, Cd, exposure=x, adjusted_for=names,
                            rng=rng, n_perm=n_perm,
                            x_raw_for_secondary=raw if secondary else None,
                            C_secondary_builder=builder if secondary else None)
        d = res.as_dict()
        d.update({"arm": arm.arm_id, "cohort": arm.cohort, "family": fam,
                  "tag": tag or f"{x}|{'+'.join(z) if z else 'nuisance_only'}"})
        rows.append(d)
        if verbose:
            print(f"  [{arm.arm_id}/{fam}] {d['tag']:<52} "
                  f"count={d['count_bh05']:<6} null_p95={d['null_count_p95']:<8.1f} "
                  f"p={d['perm_p']:.4f} resid_frac={d['residual_rank_variance_fraction']:.3f} "
                  f"-> {d['verdict']}", flush=True)
        return d

    for aspect in arm.aspects.columns:
        record("MARGINAL", aspect, [], tag=f"{aspect}|nuisance_only")
    for fam, x, z in harmonised_plan(arm):
        record(fam, x, z, secondary=(fam == "H1" and secondary_perm))
    for fam, x, z in full_plan(arm):
        record(fam, x, z)

    # --- controls -------------------------------------------------------
    if "sex_is_male" in arm.nuisance.columns:
        adj = [a for a in [COMPOSITE, FIBROSIS] if a in arm.aspects.columns]
        Cd = np.column_stack(
            [np.ones(n)] + [rankdata(arm.aspects[a].to_numpy(float)) for a in adj]
            + [arm.nuisance[c].to_numpy(float) for c in arm.nuisance.columns
               if c != "sex_is_male"])
        res = run_direction(Rfeat, rankdata(arm.nuisance.sex_is_male.to_numpy()),
                            Cd, exposure="sex_is_male", adjusted_for=adj,
                            rng=rng, n_perm=n_perm)
        d = res.as_dict()
        d.update({"arm": arm.arm_id, "cohort": arm.cohort,
                  "family": "CONTROL_POSITIVE",
                  "tag": f"sex_is_male|{'+'.join(adj)}"})
        rows.append(d)
        if verbose:
            print(f"  [{arm.arm_id}/POS] {d['tag']:<52} count={d['count_bh05']:<6} "
                  f"p={d['perm_p']:.4f} -> {d['verdict']}", flush=True)

    x_dup = arm.aspects[FIBROSIS].to_numpy(float)
    Cd = np.column_stack([np.ones(n), rankdata(x_dup)]
                         + [arm.nuisance[c].to_numpy(float)
                            for c in arm.nuisance.columns])
    res = run_direction(Rfeat, rankdata(x_dup), Cd, exposure=FIBROSIS,
                        adjusted_for=["fibrosis_exact_copy"], rng=rng,
                        n_perm=min(n_perm, 200))
    d = res.as_dict()
    d.update({"arm": arm.arm_id, "cohort": arm.cohort,
              "family": "CONTROL_NEGATIVE_N1",
              "tag": "fibrosis|fibrosis_exact_copy"})
    rows.append(d)
    if verbose:
        print(f"  [{arm.arm_id}/N1 ] {d['tag']:<52} -> {d['verdict']}", flush=True)

    shuffled = rng.permutation(arm.aspects[COMPOSITE].to_numpy(float))
    record("CONTROL_NEGATIVE_N2", COMPOSITE, [FIBROSIS], x_values=shuffled,
           tag="shuffled_composite|fibrosis")
    return rows


def testability(arm) -> dict:
    n = len(arm.participants)
    family = arm.features.shape[1]
    n_nuis = arm.nuisance.shape[1]
    out = {}
    comps = [a for a in arm.aspects.columns if a != COMPOSITE]
    for aspect in arm.aspects.columns:
        others = [a for a in comps if a != aspect]
        k = n_nuis + len(others)
        ceil = tie_ceiling(arm.aspects[aspect].to_numpy(float))
        lib = mde_spearman(n, k, 0.05)
        fam = mde_spearman(n, k, 0.05 / max(family, 1))
        call = ("untestable_ceiling_below_liberal_mde" if ceil < lib else
                "low_power_ceiling_below_family_mde" if ceil < fam else "testable")
        out[aspect] = {"tie_ceiling": round(float(ceil), 6),
                       "mde_liberal": round(float(lib), 4),
                       "mde_family": round(float(fam), 4),
                       "k_conditional": int(k), "call": call}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prespecification", required=True)
    ap.add_argument("--prespecification-sha256", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260901)
    ap.add_argument("--no-secondary-null", action="store_true")
    ap.add_argument("--arm-d-extra-lineages", default="Endothelial cells,Fibroblasts")
    args = ap.parse_args()

    spec_dir = pathlib.Path(args.prespecification)
    spec_file = spec_dir / "aspect_axis_prespecification.json"
    got = sha256_file(spec_file)
    if got != args.prespecification_sha256:
        sys.exit("PRESPECIFICATION DIGEST MISMATCH -- refusing to open any "
                 f"matrix.\n  expected {args.prespecification_sha256}\n"
                 f"  got      {got}")
    verify_frozen_tree(spec_dir)
    spec = json.loads(spec_file.read_text())
    print(f"prespecification verified before any matrix was opened: {got}")
    print(f"prespec_id: {spec['prespec_id']}")

    out = pathlib.Path(args.output)
    if out.exists():
        sys.exit(f"refusing to overwrite {out}")
    out.mkdir(parents=True)

    arms = {"A": load_arm_a(), "B": load_arm_b(), "C": load_arm_c(),
            "D": load_arm_d()}
    sealed = {(r["arm"], r["aspect"]): r for r in spec["ceiling_and_power"]}
    drift = []
    for arm in arms.values():
        for aspect, t in testability(arm).items():
            s = sealed.get((arm.arm_id, aspect))
            if s is None:
                drift.append(f"{arm.arm_id}/{aspect} absent from the sealed table")
            elif abs(s["tie_ceiling"] - t["tie_ceiling"]) > 1e-6:
                drift.append(f"{arm.arm_id}/{aspect} ceiling drifted "
                             f"{s['tie_ceiling']} -> {t['tie_ceiling']}")
    if drift:
        sys.exit("the substrate moved between Stage 0 and Stage 1:\n  "
                 + "\n  ".join(drift))
    print("every sealed ceiling re-derived identically from the substrate")

    rows = []
    for key in ["A", "B", "C", "D"]:
        arm = arms[key]
        print(f"\n=== arm {key} {arm.cohort} ({arm.assay}) n={len(arm.participants)} "
              f"family={arm.features.shape[1]} ===", flush=True)
        rows += evaluate_arm(arm, n_perm=args.n_perm, seed=args.seed + ord(key),
                             secondary_perm=not args.no_secondary_null)

    extra_lineages = [s.strip() for s in args.arm_d_extra_lineages.split(",")
                      if s.strip()]
    lineage_rows = []
    for lineage in extra_lineages:
        try:
            armx = load_arm_d(lineage=lineage, min_donors=20)
        except Exception as error:                       # noqa: BLE001
            print(f"arm D sensitivity lineage {lineage}: not admitted ({error})")
            lineage_rows.append({"arm": "D", "lineage": lineage,
                                 "admitted": False, "reason": str(error)})
            continue
        print(f"\n=== arm D sensitivity lineage {lineage} "
              f"n={len(armx.participants)} family={armx.features.shape[1]} ===",
              flush=True)
        sub = evaluate_arm(armx, n_perm=min(args.n_perm, 2000),
                           seed=args.seed + 900, secondary_perm=False)
        for r in sub:
            r["arm"] = f"D::{lineage}"
        lineage_rows += sub

    frame = pd.DataFrame(rows + lineage_rows)
    frame.to_csv(out / "directions.tsv", sep="\t", index=False)

    # ---- per-arm status and axis count --------------------------------
    def lookup(arm_id, family, tag):
        m = frame[(frame.arm == arm_id) & (frame.family == family)
                  & (frame.tag == tag)]
        return None if m.empty else m.iloc[0].to_dict()

    summary = {}
    for key in ["A", "B", "C", "D"]:
        arm = arms[key]
        t = testability(arm)
        pos = lookup(key, "CONTROL_POSITIVE",
                     f"sex_is_male|{'+'.join([a for a in [COMPOSITE, FIBROSIS] if a in arm.aspects.columns])}")
        pos_ok = None if pos is None else (pos["verdict"] == "UNIQUE")
        n1 = lookup(key, "CONTROL_NEGATIVE_N1", "fibrosis|fibrosis_exact_copy")
        n2 = lookup(key, "CONTROL_NEGATIVE_N2", "shuffled_composite|fibrosis")
        statuses, detail = {}, {}
        comps = [a for a in arm.aspects.columns if a != COMPOSITE]
        for aspect in NAMED + [BALLOON]:
            if aspect not in arm.aspects.columns:
                statuses[aspect] = "untestable"
                detail[aspect] = {"reason": "no label is deposited in this arm"}
                continue
            call = t[aspect]["call"]
            if call == "untestable_ceiling_below_liberal_mde":
                statuses[aspect] = "untestable"
                detail[aspect] = {"reason": call, **t[aspect]}
                continue
            marg = lookup(key, "MARGINAL", f"{aspect}|nuisance_only")
            cond = lookup(key, "F", f"{aspect}|{'+'.join(a for a in comps if a != aspect)}")
            channel_ok = marg is not None and marg["verdict"] == "UNIQUE"
            if cond is None:
                statuses[aspect] = "untestable"
                detail[aspect] = {"reason": "no full conditional direction exists "
                                            "in this arm (fewer than two components)",
                                  **t[aspect]}
                continue
            if cond["verdict"] == "UNIQUE":
                status = "supported"
            elif not channel_ok:
                status = "indeterminate"
            elif pos_ok is False:
                status = "indeterminate"
            elif call == "low_power_ceiling_below_family_mde":
                status = "indeterminate"
            elif cond["residual_rank_variance_fraction"] < spec[
                    "conditional_informativeness_floor"]:
                status = "indeterminate"
            elif cond["verdict"] == "EMPTY_NOT_INDEPENDENT":
                status = "not_independent"
            else:
                status = "indeterminate"
            statuses[aspect] = status
            detail[aspect] = {
                "ceiling_and_power": t[aspect],
                "marginal_channel": None if marg is None else {
                    "count": marg["count_bh05"], "perm_p": marg["perm_p"],
                    "verdict": marg["verdict"]},
                "conditional": {"tag": cond["tag"], "count": cond["count_bh05"],
                                "null_p95": cond["null_count_p95"],
                                "perm_p": cond["perm_p"],
                                "residual_rank_variance_fraction":
                                    cond["residual_rank_variance_fraction"],
                                "verdict": cond["verdict"]},
                "channel_power_known_positive_holds": bool(channel_ok),
            }
        n_axes = sum(1 for a in NAMED if statuses.get(a) == "supported")
        tested = [a for a in NAMED if statuses.get(a) not in (None, "untestable")]
        outcome = ("INDETERMINATE" if not tested else
                   {3: "THREE_AXES", 2: "TWO_AXES", 1: "ONE_AXIS"}.get(
                       n_axes, "NO_AXIS_RESOLVED"))
        if outcome == "NO_AXIS_RESOLVED" and all(
                statuses.get(a) == "indeterminate" for a in tested):
            outcome = "INDETERMINATE"
        h1 = {tag: lookup(key, "H1", tag) for tag in
              [f"{COMPOSITE}|{FIBROSIS}", f"{FIBROSIS}|{COMPOSITE}"]}
        two_axis = ("TWO_AXIS_ACTIVITY_FIBROSIS"
                    if all(v is not None and v["verdict"] == "UNIQUE"
                           for v in h1.values()) else "NOT_TWO_AXIS")
        summary[key] = {
            "cohort": arm.cohort, "assay": arm.assay,
            "n": len(arm.participants), "family_size": int(arm.features.shape[1]),
            "aspect_status": statuses,
            "named_aspect_axes_supported": n_axes,
            "named_aspects_tested": tested,
            "three_aspect_outcome": outcome,
            "activity_versus_fibrosis": two_axis,
            "controls": {
                "positive_sex": None if pos is None else {
                    "count": pos["count_bh05"], "perm_p": pos["perm_p"],
                    "verdict": pos["verdict"],
                    "reading": ("the instrument recovers a known separable axis "
                                "in this arm" if pos_ok else
                                "the instrument does NOT recover a known "
                                "separable axis here, so this arm cannot carry "
                                "a negative")},
                "negative_N1_duplicate": None if n1 is None else {
                    "verdict": n1["verdict"],
                    "passes": n1["verdict"] == "NOT_APPLICABLE_COLLINEAR"},
                "negative_N2_shuffled": None if n2 is None else {
                    "count": n2["count_bh05"], "perm_p": n2["perm_p"],
                    "verdict": n2["verdict"],
                    "passes": n2["verdict"] != "UNIQUE"},
            },
            "detail": detail,
            "notes": arm.notes,
        }

    # ---- harmonised meta ----------------------------------------------
    meta = {}
    for tag in [f"{COMPOSITE}|{FIBROSIS}", f"{FIBROSIS}|{COMPOSITE}",
                f"{METAB}|{FIBROSIS}", f"{FIBROSIS}|{METAB}",
                f"{INFLAM}|{FIBROSIS}", f"{FIBROSIS}|{INFLAM}",
                f"{METAB}|{INFLAM}", f"{INFLAM}|{METAB}",
                f"{BALLOON}|{FIBROSIS}"]:
        harmonised = {"H1", "H2", "H3", "H4", "H5"}
        got_rows = [r for _, r in frame.iterrows()
                    if r.tag == tag and r.arm in ("A", "B", "C", "D")
                    and r.family in harmonised]
        seen, deduped = set(), []
        for r in got_rows:
            if r.arm in seen:
                raise SystemExit(f"meta for {tag}: arm {r.arm} appears twice")
            seen.add(r.arm)
            deduped.append(r)
        got_rows = deduped
        if len(got_rows) < 2:
            meta[tag] = {"arms": [r.arm for r in got_rows],
                         "pooled": None,
                         "reason": "fewer than two arms deposit these labels"}
            continue
        ps = [max(float(r.perm_p), 1e-12) for r in got_rows]
        ws = [np.sqrt(float(r.n)) for r in got_rows]
        z, p, shares = stouffer(ps, ws, n_perm=args.n_perm)
        loo = {}
        for i, r in enumerate(got_rows):
            keep = [j for j in range(len(got_rows)) if j != i]
            if len(keep) >= 2:
                zl, pl, _ = stouffer([ps[j] for j in keep], [ws[j] for j in keep],
                                     n_perm=args.n_perm)
            elif keep:
                zl, pl = float(-np.log(ps[keep[0]])), ps[keep[0]]
            else:
                zl, pl = float("nan"), float("nan")
            loo[f"without_{r.arm}"] = {"combined_z": round(float(zl), 3),
                                       "combined_p": float(pl)}
        meta[tag] = {
            "arms": [r.arm for r in got_rows],
            "pooled_n": int(sum(int(r.n) for r in got_rows)),
            "per_arm": [{"arm": r.arm, "cohort": r.cohort, "n": int(r.n),
                         "family_size": int(r.family_size),
                         "count": int(r.count_bh05),
                         "null_p95": float(r.null_count_p95),
                         "perm_p": float(r.perm_p),
                         "verdict": r.verdict,
                         "z_share": round(float(s), 3)}
                        for r, s in zip(got_rows, shares)],
            "combined_z": round(float(z), 3), "combined_p": float(p),
            "leave_one_arm_out": loo,
            "concentration": {
                "largest_single_arm_z_share": round(float(max(shares) / z), 3)
                    if z > 0 else None,
                "reading": ("a share near 1 means the pooled statistic is one "
                            "arm's result wearing a pooled label"),
            },
            "instrument_is_identical_across_these_arms": True,
        }

    payload = {
        "prespec_id": spec["prespec_id"],
        "prespecification_sha256": got,
        "prespecification_verified_before_any_matrix_was_opened": True,
        "n_permutations": args.n_perm,
        "seed": args.seed,
        "per_arm": summary,
        "harmonised_meta": meta,
        "arms_are_four_instruments_and_are_never_pooled_at_the_feature_level": True,
        "arm_d_sensitivity_lineages": [
            {"lineage": r.get("lineage"), "admitted": r.get("admitted"),
             "reason": r.get("reason")}
            for r in lineage_rows if isinstance(r, dict) and "admitted" in r],
    }
    (out / "aspect_axis_result.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str),
        encoding="utf-8")

    print("\n=== per-arm outcome ===")
    for key, s in summary.items():
        print(f"{key} {s['cohort']:<12} n={s['n']:<4} "
              f"{s['three_aspect_outcome']:<16} axes={s['named_aspect_axes_supported']} "
              f"{s['activity_versus_fibrosis']:<28} "
              + " ".join(f"{a.split('_')[0]}={s['aspect_status'][a]}"
                         for a in NAMED))
    print("\n=== harmonised meta ===")
    for tag, m in meta.items():
        if m.get("pooled") is None and "combined_z" not in m:
            print(f"{tag:<44} {m.get('reason')}")
            continue
        print(f"{tag:<44} arms={'+'.join(m['arms'])} n={m['pooled_n']} "
              f"z={m['combined_z']} p={m['combined_p']:.2e} "
              f"max_arm_share={m['concentration']['largest_single_arm_z_share']}")

    freeze_tree(out, {"artifact_class": "aspect_axis_count_four_arm",
                      "stage": "stage_1_evaluation",
                      "prespec_id": spec["prespec_id"],
                      "prespecification_sha256": got})
    verify_frozen_tree(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
