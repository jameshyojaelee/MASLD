#!/usr/bin/env python3
"""Defect quantification: the Track 0 direction columns use a liver panel that includes EMBRYONIC tracks.

`direction_variant_level.tsv.gz` takes its predicted allele direction from `liver_summaries.tsv.gz`, whose
liver median pools `primary_liver` (5 adult RNA tracks) with `hepatocyte` = CL:0000182 (2 tracks), and
CL:0000182 is embryonic in this Atlas build. Every other adult-liver column in the same dossier row is
restricted to `q|primary_liver|`. The two are therefore not on the same track panel, and the deposited
direction concordance is a pooled-panel number.

This script is a nested increment: one component changes, the track set. Both arms are computed here, with
the same loader, from the same parquet raw columns, and the pooled arm must first REPRODUCE the deposited
`pred_allele1_liver_rna` before the adult-only arm is believed.

Outputs into the response-layer directory:
  direction_track_panel_check.tsv   per variant row: pooled and adult-only prediction, both directions
  direction_track_panel_summary.json

Usage: b1_direction_track_panel.py <response layer directory>
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
POOLED = ["raw|primary_liver|UBERON:0001114 total RNA-seq",
          "raw|primary_liver|UBERON:0001115 total RNA-seq",
          "raw|primary_liver|UBERON:0002107 polyA plus RNA-seq",
          "raw|primary_liver|UBERON:0002107 total RNA-seq",
          "raw|primary_liver|UBERON:0001114 gtex Liver polyA plus RNA-seq",
          "raw|hepatocyte|CL:0000182 total RNA-seq",
          "raw|hepatocyte|CL:0000182 polyA plus RNA-seq"]
ADULT = [c for c in POOLED if c.startswith("raw|primary_liver|")]
EMBRYONIC = [c for c in POOLED if c.startswith("raw|hepatocyte|")]


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
    for r in rows_of(TRACK0 / "direction_variant_level.tsv.gz"):
        if r["atlas_scored"] != "True" or not r["variant_uid"]:
            continue
        g = signals.get(r["signal_uid"], {}).get("ensembl", "")
        if not g:
            continue
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
    pooled, adult, embryo = {}, {}, {}
    pf = pq.ParquetFile(RECORDS / "RNA_SEQ.parquet")
    for batch in pf.iter_batches(batch_size=200_000, columns=["variant_uid", "gene_id"] + POOLED):
        d = batch.to_pydict()
        for i, uid in enumerate(d["variant_uid"]):
            k = (uid, d["gene_id"][i])
            if k not in need:
                continue
            pooled[k] = _median([_f(d[c][i]) for c in POOLED])
            adult[k] = _median([_f(d[c][i]) for c in ADULT])
            embryo[k] = _median([_f(d[c][i]) for c in EMBRYONIC])
    print(f"medians recovered for {len(pooled)} pairs", file=sys.stderr, flush=True)

    out, repro_ok, repro_n, repro_maxdiff = [], 0, 0, 0.0
    for r in drows:
        k = (r["variant_uid"], r["gene_id"])
        sw = swap.get((r["signal_uid"], r["variant_uid"]), "")
        p, a, e = pooled.get(k), adult.get(k), embryo.get(k)
        if p is None:
            continue

        def orient(v):
            return None if v is None else (-v if sw == "none" else v)
        pa, aa, ea = orient(p), orient(a), orient(e)
        dep = _f(r["deposited_pred_allele1"])
        if dep is not None:
            repro_n += 1
            diff = abs(pa - dep)
            repro_maxdiff = max(repro_maxdiff, diff)
            if diff <= 1e-9 or (dep != 0 and diff / abs(dep) <= 1e-9):
                repro_ok += 1
        meas = _f(r["eqtl_beta_allele1"])
        def sgn(v):
            return 0 if v is None or v == 0 else (1 if v > 0 else -1)
        out.append({**r, "allele_swap": sw,
                    "pooled_pred_allele1": ("" if pa is None else repr(pa)),
                    "adult_only_pred_allele1": ("" if aa is None else repr(aa)),
                    "embryonic_only_pred_allele1": ("" if ea is None else repr(ea)),
                    "pooled_direction": sgn(pa), "adult_only_direction": sgn(aa),
                    "sign_flips_when_embryonic_tracks_removed": str(sgn(pa) != sgn(aa)),
                    "pooled_concordant": ("" if sgn(pa) == 0 or sgn(meas) == 0 else str(sgn(pa) == sgn(meas))),
                    "adult_only_concordant": ("" if sgn(aa) == 0 or sgn(meas) == 0 else str(sgn(aa) == sgn(meas)))})

    cols = list(out[0].keys())
    path = outdir / "direction_track_panel_check.tsv.gz"
    if path.exists():
        raise RuntimeError(f"refusing to overwrite: {path}")
    with gzip.open(path, "wt") as h:
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow(r)

    def conc(key):
        s = [r for r in out if r[key] in ("True", "False")]
        k = sum(r[key] == "True" for r in s)
        return {"n": len(s), "n_concordant": k, "fraction": (k / len(s)) if s else None}

    def wconc(key):
        s = [r for r in out if r[key] in ("True", "False")]
        tot = sum(float(r["weight"]) for r in s)
        pos = sum(float(r["weight"]) for r in s if r[key] == "True")
        return {"posterior_weighted_fraction": (pos / tot) if tot else None, "total_weight": tot}

    summary = {
        "what_changed": "one component: the RNA track set behind the predicted allele direction",
        "pooled_tracks": POOLED, "adult_only_tracks": ADULT, "embryonic_tracks_dropped": EMBRYONIC,
        "allele_swap_keyed_on": "(signal_uid, variant_uid) -> normalised source id -> crosswalk allele_swap",
        "reproduction_of_deposited_pred_allele1": {
            "n_rows_with_a_deposited_value": repro_n, "n_reproduced_to_1e-9": repro_ok,
            "max_absolute_difference": repro_maxdiff,
            "why": ("the adult-only arm is only interpretable if the pooled arm regenerates the deposited "
                    "number through this loader")},
        "n_rows": len(out),
        "n_sign_flips_when_embryonic_tracks_removed": sum(
            r["sign_flips_when_embryonic_tracks_removed"] == "True" for r in out),
        "share_sign_flips": (sum(r["sign_flips_when_embryonic_tracks_removed"] == "True" for r in out) / len(out)
                             if out else None),
        "concordance_pooled": conc("pooled_concordant"),
        "concordance_adult_only": conc("adult_only_concordant"),
        "posterior_weighted_concordance_pooled": wconc("pooled_concordant"),
        "posterior_weighted_concordance_adult_only": wconc("adult_only_concordant"),
        "scope": ("sign concordance over Atlas-scored credible-set variants; these are LD-correlated variants "
                  "within signals, not independent tests, and no p-value is computed here"),
    }
    json.dump(summary, (outdir / "direction_track_panel_summary.json").open("w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
