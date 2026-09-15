#!/usr/bin/env python3
"""Defect quantification v2: the Track 0 direction columns pool track sets that are not the same panel.

`direction_variant_level.tsv.gz` takes its predicted allele direction from `liver_summaries.tsv.gz`, whose
liver median pools the `primary_liver` group (5 RNA tracks) with `hepatocyte` = CL:0000182 (2 tracks).
Every other liver column in the same dossier row is restricted to `q|primary_liver|`.

The wave-1 run called the two arms 'pooled' and 'adult'. That label was wrong. By this build's own
raw/scorer_metadata the 5 primary_liver RNA tracks are 3 adult, 1 child (UBERON:0002107 polyA plus RNA-seq)
and 1 embryonic (UBERON:0002107 total RNA-seq), and the 2 CL:0000182 tracks are embryonic. So the wave-1
increment is 'minus CL:0000182', not 'pooled versus adult'. v2 states it that way and adds the increment the
corrected label implies: adult-only versus all primary_liver.

Three arms, each differing from its comparator in ONE component, the track set:
  A  all 7 tracks   primary_liver 5 + CL:0000182 2        = the panel the deposit used
  B  5 tracks       primary_liver only                    = A minus CL:0000182
  C  3 tracks       primary_liver tracks whose life stage is adult
Increment 1 = B - A. Increment 2 = C - B. Arm A must first REPRODUCE the deposited `pred_allele1_liver_rna`
through this loader before either increment is believed.

Denominators, stated because the wave-1 deposit printed 877,512 without defining it:
  3,452,688  rows in direction_variant_level.tsv.gz
    877,512  panel rows: atlas_scored, a non-empty variant_uid, a signal with an ensembl gene, and an RNA
             parquet median available -- the rows this script can recompute
    741,737  of those that carry a deposited pred_allele1_liver_rna value
Sign-flip counts are reported on both the 877,512 and the 741,737 denominators.

Outputs into the response-layer directory:
  direction_track_panel_check.tsv.gz   per variant row: all three arms and their directions
  direction_track_panel_summary.json

Usage: b1_direction_track_panel_v2.py <response layer directory>
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import pathlib
import re
import statistics
import sys
from collections import defaultdict

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TRACK0 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables"
RECORDS = TRACK0 / "prediction_records"
TRACK_META = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/raw/scorer_metadata"
ALL_TRACKS = ["raw|primary_liver|UBERON:0001114 total RNA-seq",
              "raw|primary_liver|UBERON:0001115 total RNA-seq",
              "raw|primary_liver|UBERON:0002107 polyA plus RNA-seq",
              "raw|primary_liver|UBERON:0002107 total RNA-seq",
              "raw|primary_liver|UBERON:0001114 gtex Liver polyA plus RNA-seq",
              "raw|hepatocyte|CL:0000182 total RNA-seq",
              "raw|hepatocyte|CL:0000182 polyA plus RNA-seq"]
PRIMARY_LIVER = [c for c in ALL_TRACKS if c.startswith("raw|primary_liver|")]
CL0000182 = [c for c in ALL_TRACKS if c.startswith("raw|hepatocyte|")]


def life_stages():
    """track column -> life stage, from this build's own scorer metadata. Never inferred from the name."""
    meta = {}
    with open(TRACK_META / "RNA_SEQ.track_metadata.tsv") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            meta[r["name"]] = r["biosample_life_stage"]
    out = {}
    for c in ALL_TRACKS:
        name = c.split("|", 2)[2]
        if name not in meta:
            raise RuntimeError(f"no scorer metadata row for RNA track {name!r}")
        out[c] = meta[name]
    return out


STAGES = life_stages()
ADULT_PRIMARY_LIVER = [c for c in PRIMARY_LIVER if STAGES[c] == "adult"]
ARMS = {"all_7_tracks_primary_liver_plus_CL0000182": ALL_TRACKS,
        "primary_liver_5_tracks": PRIMARY_LIVER,
        "adult_primary_liver_3_tracks": ADULT_PRIMARY_LIVER}


def rows_of(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as h:
        yield from csv.DictReader(h, delimiter="\t")


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def _median(vals):
    v = [x for x in vals if x is not None]
    return statistics.median(v) if v else None


SOURCE_ID_RE = re.compile(r"^(?P<asm>[A-Za-z0-9]+):(?:chr)?(?P<chrom>[0-9]{1,2}|[XYM]|MT):(?P<pos>\d+):"
                          r"(?P<ref>[ACGTNacgtn]+):(?P<alt>[ACGTNacgtn]+)$")


def norm_src(s):
    m = SOURCE_ID_RE.match(str(s).strip())
    return None if m is None else f"{m['asm']}:chr{m['chrom']}:{int(m['pos'])}:{m['ref'].upper()}:{m['alt'].upper()}"


def main():
    outdir = pathlib.Path(sys.argv[1]).resolve()
    signals = {r["signal_uid"]: r for r in rows_of(TRACK0 / "eligible_signals.tsv")}

    need, drows, pairs = set(), [], set()
    denom = {"rows_in_direction_variant_level": 0, "dropped_not_atlas_scored": 0,
             "dropped_no_variant_uid": 0, "dropped_signal_has_no_ensembl_gene": 0,
             "eligible_before_parquet_join": 0}
    for r in rows_of(TRACK0 / "direction_variant_level.tsv.gz"):
        denom["rows_in_direction_variant_level"] += 1
        if r["atlas_scored"] != "True":
            denom["dropped_not_atlas_scored"] += 1
            continue
        if not r["variant_uid"]:
            denom["dropped_no_variant_uid"] += 1
            continue
        g = signals.get(r["signal_uid"], {}).get("ensembl", "")
        if not g:
            denom["dropped_signal_has_no_ensembl_gene"] += 1
            continue
        denom["eligible_before_parquet_join"] += 1
        need.add((r["variant_uid"], g))
        pairs.add((r["signal_uid"], r["variant_uid"]))
        drows.append({"signal_uid": r["signal_uid"], "variant_uid": r["variant_uid"], "gene_id": g,
                      "weight": r["weight"], "deposited_pred_allele1": r["pred_allele1_liver_rna"],
                      "eqtl_beta_allele1": r["eqtl_beta_allele1"],
                      "deposited_predicted_direction": r["predicted_direction"],
                      "measured_direction": r["measured_direction"],
                      "deposited_concordant": r["concordant"]})
    print(f"direction rows scored by the Atlas: {len(drows)}; distinct (variant, gene) pairs: {len(need)}",
          file=sys.stderr, flush=True)

    # allele_swap belongs to the (source id, hg38 uid) PAIR, not to the uid: 501,576 of 1,229,204 uids in the
    # crosswalk carry both 'none' and 'swapped' because different studies list the two alleles in either
    # order. Keying it on the uid alone inverts the predicted direction for whichever rows lost the draw.
    src_of = {}
    for r in rows_of(TRACK0 / "signal_variant_weights.tsv.gz"):
        k = (r["signal_uid"], r["variant_uid"])
        if r["variant_uid"] and k in pairs:
            src_of[k] = norm_src(r["source_variant_id"])
    swap_by_src = {}
    for r in rows_of(TRACK0 / "variant_crosswalk.tsv.gz"):
        k = norm_src(r["source_variant_id"])
        if k:
            swap_by_src[k] = r["allele_swap"]
    swap = {k: swap_by_src.get(v, "") for k, v in src_of.items()}
    print(f"allele_swap resolved for {sum(1 for v in swap.values() if v)} of {len(pairs)} (signal, variant) pairs",
          file=sys.stderr, flush=True)

    import pyarrow.parquet as pq
    arm_med = {a: {} for a in ARMS}
    cl_only = {}
    pf = pq.ParquetFile(RECORDS / "RNA_SEQ.parquet")
    for batch in pf.iter_batches(batch_size=200_000, columns=["variant_uid", "gene_id"] + ALL_TRACKS):
        d = batch.to_pydict()
        for i, uid in enumerate(d["variant_uid"]):
            k = (uid, d["gene_id"][i])
            if k not in need:
                continue
            for a, cols in ARMS.items():
                arm_med[a][k] = _median([_f(d[c][i]) for c in cols])
            cl_only[k] = _median([_f(d[c][i]) for c in CL0000182])
    base = "all_7_tracks_primary_liver_plus_CL0000182"
    print(f"medians recovered for {len(arm_med[base])} pairs", file=sys.stderr, flush=True)

    def sgn(v):
        return 0 if v is None or v == 0 else (1 if v > 0 else -1)

    out, repro_ok, repro_n, repro_maxdiff = [], 0, 0, 0.0
    for r in drows:
        k = (r["variant_uid"], r["gene_id"])
        sw = swap.get((r["signal_uid"], r["variant_uid"]), "")
        if arm_med[base].get(k) is None:
            continue

        def orient(v):
            return None if v is None else (-v if sw == "none" else v)
        vals = {a: orient(arm_med[a].get(k)) for a in ARMS}
        pa = vals[base]
        dep = _f(r["deposited_pred_allele1"])
        if dep is not None:
            repro_n += 1
            diff = abs(pa - dep)
            repro_maxdiff = max(repro_maxdiff, diff)
            if diff <= 1e-9 or (dep != 0 and diff / abs(dep) <= 1e-9):
                repro_ok += 1
        meas = _f(r["eqtl_beta_allele1"])
        row = {**r, "allele_swap": sw, "has_deposited_direction": str(dep is not None),
               "CL0000182_only_pred_allele1": ("" if orient(cl_only.get(k)) is None
                                               else repr(orient(cl_only.get(k))))}
        for a in ARMS:
            row[f"pred_allele1__{a}"] = ("" if vals[a] is None else repr(vals[a]))
            row[f"direction__{a}"] = sgn(vals[a])
            row[f"concordant__{a}"] = ("" if sgn(vals[a]) == 0 or sgn(meas) == 0
                                       else str(sgn(vals[a]) == sgn(meas)))
        row["sign_flips_increment1_minus_CL0000182"] = str(
            sgn(vals[base]) != sgn(vals["primary_liver_5_tracks"]))
        row["sign_flips_increment2_adult_only_vs_all_primary_liver"] = str(
            sgn(vals["primary_liver_5_tracks"]) != sgn(vals["adult_primary_liver_3_tracks"]))
        out.append(row)

    cols = list(out[0].keys())
    path = outdir / "direction_track_panel_check.tsv.gz"
    if path.exists():
        raise RuntimeError(f"refusing to overwrite: {path}")
    with gzip.open(path, "wt") as h:
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow(r)

    def conc(key, sel=None):
        s = [r for r in (sel if sel is not None else out) if r[key] in ("True", "False")]
        k = sum(r[key] == "True" for r in s)
        return {"n": len(s), "n_concordant": k, "fraction": (k / len(s)) if s else None}

    def wconc(key, sel=None):
        s = [r for r in (sel if sel is not None else out) if r[key] in ("True", "False")]
        tot = sum(float(r["weight"]) for r in s)
        pos = sum(float(r["weight"]) for r in s if r[key] == "True")
        return {"posterior_weighted_fraction": (pos / tot) if tot else None, "total_weight": tot}

    dep_rows = [r for r in out if r["has_deposited_direction"] == "True"]

    def flips(col):
        n_all = sum(r[col] == "True" for r in out)
        n_dep = sum(r[col] == "True" for r in dep_rows)
        return {"n_panel_rows": len(out), "flips_on_panel_rows": n_all,
                "share_of_panel_rows": (n_all / len(out)) if out else None,
                "n_rows_carrying_a_deposited_direction": len(dep_rows),
                "flips_on_rows_carrying_a_deposited_direction": n_dep,
                "share_of_rows_carrying_a_deposited_direction": (n_dep / len(dep_rows)) if dep_rows else None}

    summary = {
        "what_changed": "one component in each increment: the RNA track set behind the predicted direction",
        "arms": {a: {"tracks": cols, "n_tracks": len(cols),
                     "life_stages": {c: STAGES[c] for c in cols}} for a, cols in ARMS.items()},
        "CL0000182_tracks_dropped_in_increment_1": {c: STAGES[c] for c in CL0000182},
        "increment_1": "primary_liver_5_tracks minus all_7_tracks (i.e. drop CL:0000182)",
        "increment_2": "adult_primary_liver_3_tracks minus primary_liver_5_tracks (drop child + embryonic)",
        "label_correction": ("the wave-1 deposit called these arms 'pooled 7 tracks' and 'adult 5 tracks'. "
                             "The 5 primary_liver RNA tracks are 3 adult, 1 child and 1 embryonic, so the "
                             "wave-1 increment is 'minus CL:0000182', not 'pooled versus adult'."),
        "allele_swap_keyed_on": "(signal_uid, variant_uid) -> normalised source id -> crosswalk allele_swap",
        "denominators": {**denom, "panel_rows_recomputable": len(out),
                         "panel_rows_carrying_a_deposited_direction": len(dep_rows),
                         "what_877512_is": ("rows of direction_variant_level.tsv.gz that are atlas_scored, "
                                            "carry a variant_uid, belong to a signal with an ensembl gene, "
                                            "and have an RNA parquet median; NOT all direction rows")},
        "reproduction_of_deposited_pred_allele1": {
            "arm": base, "n_rows_with_a_deposited_value": repro_n, "n_reproduced_to_1e-9": repro_ok,
            "max_absolute_difference": repro_maxdiff,
            "why": ("neither increment is interpretable unless the deposited panel regenerates the deposited "
                    "number through this loader")},
        "n_rows": len(out),
        "sign_flips_increment1_minus_CL0000182": flips("sign_flips_increment1_minus_CL0000182"),
        "sign_flips_increment2_adult_only_vs_all_primary_liver": flips(
            "sign_flips_increment2_adult_only_vs_all_primary_liver"),
        "concordance_per_arm": {a: conc(f"concordant__{a}") for a in ARMS},
        "posterior_weighted_concordance_per_arm": {a: wconc(f"concordant__{a}") for a in ARMS},
        "concordance_per_arm_on_rows_with_a_deposited_direction": {
            a: conc(f"concordant__{a}", dep_rows) for a in ARMS},
        "scope": ("sign concordance over Atlas-scored credible-set variants; these are LD-correlated variants "
                  "within signals, not independent tests, and no p-value is computed here"),
    }
    json.dump(summary, (outdir / "direction_track_panel_summary.json").open("w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
