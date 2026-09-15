#!/usr/bin/env python3
"""Step 65 (P6d): a per-variant clinical dossier for the direct MASLD trait signals.

One row per (signal_uid, variant_uid) for the variants that carry the posterior mass at a direct MASLD
signal. Each row joins, for that exact variant: the Atlas adult-liver predictions, the chromatin we measured
in donors, the liver eQTL direction, which ancestries' credible sets contain it, whether its mechanism is
already known from experiments (the P2 anchor panel), and whether the Atlas annotation at its parent signal
is complete enough to be read at all.

The last column is the point of the lane. An uncovered signal still has a profile row, and that row's
top-variant columns are unweighted maxima over near-zero-weight SNVs while the posterior-weighted columns
correctly read ~1e-23. A reader who takes a populated column for a measured-quiet signal is reading an
artefact, so every row here travels with the parent signal's coverage state and its own served flag.

Nothing here is a new threshold, ranking, nomination or evidence class: every column is an export of an
adopted number or an existing prediction, plus the joins between them.

Prespecification: 65_clinical_dossier_prespec.json (written before any dossier row was read).
Outputs (tables/): clinical_variant_dossier.tsv, clinical_variant_dossier_summary.json
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import statistics
from collections import defaultdict

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
TRACK0 = la.track0_root() / "tables"
RECORDS = TRACK0 / "prediction_records"
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
PRESPEC = SCRIPT_DIR / "65_clinical_dossier_prespec.json"

WEIGHT_FLOOR = 0.01
DIRECT_UNIVERSES = ("A_direct", "B_direct")

# Adult liver only. The hepatocyte ontology CL:0000182 is EMBRYONIC in this Atlas build, so it is excluded
# here even though liver_summaries pools it; the two tables will therefore not print identical numbers.
ADULT_LIVER_PREFIX = "q|primary_liver|"
TRACK_PANEL = {
    "atac_liver_quantile": ("ATAC.parquet", lambda c: c.startswith(ADULT_LIVER_PREFIX) and c.endswith("ATAC-seq")),
    "dnase_liver_quantile": ("DNASE.parquet", lambda c: c.startswith(ADULT_LIVER_PREFIX) and c.endswith("DNase-seq")),
    "h3k27ac_liver_quantile": ("CHIP_HISTONE.parquet", lambda c: c.startswith(ADULT_LIVER_PREFIX) and c.endswith("H3K27ac")),
}
ATLAS_QUANTILE_COLUMNS = ["atac_liver_quantile", "dnase_liver_quantile", "h3k27ac_liver_quantile",
                          "rna_target_gene_quantile", "rna_strongest_abs_quantile",
                          "splice_site_usage_liver_quantile", "avi_quantile"]
MODEL_API_COLUMNS = ["model_api_rna_log2", "model_api_splice_log2", "model_api_atac_log2",
                     "model_api_dnase_log2", "model_api_h3k27ac_log2"]

# One implementation of both identities, shared with every other step.
normalise_source_id = la.normalise_source_id
row_key = la.posterior_key


def rescue_index(rescue_rows) -> dict:
    """Normalised source id -> scored step-63 indel row, under EVERY id of that variant.

    Step 63 scores one row per distinct variant under its heaviest id and lists all ids in
    `member_source_ids`; a row written before that grouping is indexed by its own id.
    """
    out = {}
    for r in rescue_rows:
        if r.get("arm") != "indel" or r.get("state") != "scored":
            continue
        for sid in str(r.get("member_source_ids") or r["source_variant_id"]).split(";"):
            out[normalise_source_id(sid)] = r
    return out


def rescue_for_row(index: dict, source_variant_id: str, ref: str, alt: str) -> tuple:
    """(rescue row, state) for a dossier row; an effect scored on other alleles is withheld, not attached."""
    r = index.get(normalise_source_id(source_variant_id))
    if not r:
        return {}, ""
    if (str(r["ref"]).upper(), str(r["alt"]).upper()) != (str(ref).upper(), str(alt).upper()):
        return {}, "orientation_differs"
    return r, "scored"


def variant_class(ref: str, alt: str) -> str:
    r, a = str(ref).upper(), str(alt).upper()
    if len(r) == len(a):
        return "snv" if len(r) == 1 else "mnv"
    return "insertion" if len(a) > len(r) else "deletion"


def orientation_source(dbsnp_uid: str, recovered_uid: str) -> str:
    """Which table named this variant's reference allele, so a reader can stratify on it."""
    if str(dbsnp_uid or "").strip():
        return "dbsnp"
    if str(recovered_uid or "").strip():
        return "reference_or_source_order"
    return ""


def alleles_of(variant_uid: str, source_variant_id: str, recovered_uid: str = "",
               dbsnp_uid: str = "") -> tuple:
    """Reference and alternate alleles, in hg38 reference orientation wherever that is known.

    Precedence: the Atlas-served uid, then dbSNP, then the reference-only repair, then the raw source id.
    dbSNP outranks the repair because it holds each variant's canonical REF/ALT; on the 98,947
    reference-ambiguous indels it resolves it flips the repair's call on 53.6%, and against the allele order
    of files with a fixed convention it is right on 99.9% (step 68b). A source id's allele order depends on
    the file it came from (UKBB lists ALT first, MVP and BBJ REF first, others mixed), so it is the last
    resort.
    """
    for uid in (str(variant_uid or "").strip(), str(dbsnp_uid or "").strip(),
                str(recovered_uid or "").strip()):
        if uid:
            parts = uid.split(":")
            if len(parts) != 4:
                raise la.ContractError(f"variant uid is not chrom:pos:ref:alt: {uid!r}")
            return parts[2], parts[3]
    parts = normalise_source_id(source_variant_id).split(":")
    return parts[3], parts[4]


def rank_and_cumulative(weights: dict) -> dict:
    """Rank and cumulative share over the FULL posterior of a signal, never over the exported subset.

    The weights are not renormalised: a signal whose posterior sums to 0.5 keeps that denominator, so the
    printed share never overstates how much of the signal a variant carries.
    """
    if not weights:
        raise la.ContractError("a signal with no posterior weights cannot be ranked")
    order = sorted(weights.items(), key=lambda kv: (-float(kv[1]), kv[0]))
    out, run = {}, 0.0
    for i, (uid, w) in enumerate(order, start=1):
        run += float(w)
        out[uid] = {"rank": i, "weight": float(w), "cumulative_share": run}
    return out


def select_dossier_variants(weights: dict, floor: float = WEIGHT_FLOOR) -> set:
    """Variants above the floor, plus the top-weight variant whatever its weight.

    A signal whose whole posterior sits on deferred indels has no variant above the floor; without the
    second rule it would simply vanish from the export, which is the opposite of the point.
    """
    if not weights:
        raise la.ContractError("a signal with no posterior weights has nothing to export")
    sel = {u for u, w in weights.items() if float(w) >= floor}
    top = sorted(weights.items(), key=lambda kv: (-float(kv[1]), kv[0]))[0][0]
    sel.add(top)
    return sel


def row_readable(coverage_state: str, served: bool) -> tuple:
    """A row is readable only when the variant has a prediction AND its parent signal's mass is covered."""
    if not served:
        return "False", "variant not served by the Atlas point-query API (indel or unmapped)"
    if coverage_state != "covered":
        return "False", f"parent signal coverage state is {coverage_state}; its signal-level columns rest on incomplete mass"
    return "True", ""


def guard_unserved_row(row: dict) -> None:
    """An unserved variant may never print a numeric Atlas quantile a reader would take for a null effect."""
    if str(row.get("atlas_served", "")) == "False":
        for c in ATLAS_QUANTILE_COLUMNS:
            if str(row.get(c, "")).strip() != "":
                raise la.ContractError(f"unserved variant carries an Atlas quantile in {c}: {row.get(c)!r}")


def direction_concordance(rows, min_rows: int = 30) -> dict:
    scored = [r for r in rows if str(r.get("direction_concordant", "")) in ("True", "False")]
    n = len(scored)
    if n < min_rows:
        return {"fraction": None, "n": n, "min_rows": min_rows,
                "reason": f"withheld: {n} rows carry both a measured and a predicted direction, below the prespecified minimum of {min_rows}"}
    k = sum(str(r["direction_concordant"]) == "True" for r in scored)
    return {"fraction": k / n, "n": n, "n_concordant": k, "min_rows": min_rows, "reason": "reported"}


def measured_peak_share(n_in_peak: int, n_with_peak_data: int, n_rows: int) -> dict:
    """The share is reported over the rows the peak table actually covers, with that denominator printed."""
    if n_with_peak_data <= 0:
        return {"share_over_covered": None, "denominator": 0, "n_rows": n_rows,
                "coverage_too_thin_to_extrapolate": True}
    return {"share_over_covered": n_in_peak / n_with_peak_data, "denominator": n_with_peak_data,
            "n_rows": n_rows, "coverage_too_thin_to_extrapolate": n_with_peak_data < 0.5 * n_rows}


def _f(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _median(values):
    vals = [v for v in values if v is not None]
    return statistics.median(vals) if vals else None


def _fmt(value):
    return "" if value is None else repr(float(value))


def load_signals() -> dict:
    return {r["signal_uid"]: r for r in la.read_tsv(TRACK0 / "eligible_signals.tsv")
            if r["universe"] in DIRECT_UNIVERSES}


def load_weights(signals: dict) -> tuple:
    """Full per-signal weight vectors for the direct signals, plus per-variant ancestry membership.

    Ancestry membership is taken at the same 0.01 floor as inclusion: a variant present in a credible set
    only at PIP 1e-6 is not a member a clinician would act on.
    """
    per_signal, meta, ancestries = defaultdict(dict), {}, defaultdict(set)
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            sig = signals.get(r["signal_uid"])
            if sig is None:
                continue
            w = float(r["weight"])
            key = row_key(r["variant_uid"], r["source_variant_id"])
            if key in per_signal[r["signal_uid"]]:
                raise la.ContractError(f"{r['signal_uid']}: two posterior entries share the identity {key}")
            per_signal[r["signal_uid"]][key] = w
            meta[key] = {"variant_uid": r["variant_uid"], "source_variant_id": r["source_variant_id"],
                         "mapping_status": r["mapping_status"], "exclusion_reason": r["exclusion_reason"]}
            if w >= WEIGHT_FLOOR and sig["ancestry"]:
                ancestries[key].add(sig["ancestry"])
    return per_signal, meta, ancestries


def load_atlas_channels(uids: set) -> dict:
    import pyarrow.parquet as pq

    out = defaultdict(dict)
    for col, (fname, keep) in TRACK_PANEL.items():
        pf = pq.ParquetFile(RECORDS / fname)
        cols = [c for c in pf.schema.names if keep(c)]
        if not cols:
            raise la.ContractError(f"no adult-liver track column matched for {col} in {fname}")
        t = pf.read(columns=["variant_uid"] + cols).to_pydict()
        for i, uid in enumerate(t["variant_uid"]):
            if uid in uids:
                out[uid][col] = _median([_f(t[c][i]) for c in cols])
        out["__tracks__"][col] = cols

    pf = pq.ParquetFile(RECORDS / "AVI_SCORE.parquet")
    t = pf.read(columns=["variant_uid", "raw|feature|AVI_SCORE", "q|feature|AVI_SCORE"]).to_pydict()
    for i, uid in enumerate(t["variant_uid"]):
        if uid in uids:
            out[uid]["avi_raw"] = _f(t["raw|feature|AVI_SCORE"][i])
            out[uid]["avi_quantile"] = _f(t["q|feature|AVI_SCORE"][i])
    return out


def load_gene_channels(uids: set) -> tuple:
    """Per-variant, per-gene adult-liver RNA and splice-site-usage quantiles."""
    import pyarrow.parquet as pq

    rna, splice, tracks = defaultdict(dict), defaultdict(dict), {}
    for fname, sink, key in (("RNA_SEQ.parquet", rna, "rna"), ("SPLICE_SITE_USAGE.parquet", splice, "splice")):
        pf = pq.ParquetFile(RECORDS / fname)
        cols = [c for c in pf.schema.names if c.startswith(ADULT_LIVER_PREFIX)]
        if not cols:
            raise la.ContractError(f"no adult-liver track column in {fname}")
        tracks[key] = cols
        for batch in pf.iter_batches(batch_size=200_000, columns=["variant_uid", "gene_id"] + cols):
            d = batch.to_pydict()
            for i, uid in enumerate(d["variant_uid"]):
                if uid not in uids:
                    continue
                m = _median([_f(d[c][i]) for c in cols])
                if m is not None:
                    sink[uid][d["gene_id"][i]] = m
    return rna, splice, tracks


def main() -> None:
    prespec = json.load(PRESPEC.open())
    signals = load_signals()
    per_signal, meta, ancestries = load_weights(signals)
    if len(per_signal) != len(signals):
        raise la.ContractError(f"weights cover {len(per_signal)} of {len(signals)} direct signals")

    ranks = {s: rank_and_cumulative(w) for s, w in per_signal.items()}
    selected = {s: select_dossier_variants(w) for s, w in per_signal.items()}
    keys = {k for sel in selected.values() for k in sel}
    uids = {meta[k]["variant_uid"] for k in keys if meta[k]["variant_uid"]}
    la.log(f"P6d: {len(signals)} direct signals, {sum(len(v) for v in selected.values())} rows, "
           f"{len(keys)} distinct posterior entries, {len(uids)} of them with an hg38 uid")

    served = {}
    with la.open_text(TRACK0 / "atlas_availability.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["variant_uid"] in uids:
                served[r["variant_uid"]] = r

    channels = load_atlas_channels(uids)
    rna, splice, gene_tracks = load_gene_channels(uids)

    peaks = defaultdict(list)
    with la.open_text(TRACK0 / "variant_peak_overlaps.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["signal_uid"] in signals and r["variant_uid"] in uids:
                peaks[(r["signal_uid"], r["variant_uid"])].append(r)

    directions = {}
    with la.open_text(TRACK0 / "direction_variant_level.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["signal_uid"] in signals and r["variant_uid"] in uids:
                directions[(r["signal_uid"], r["variant_uid"])] = r

    coverage = {r["signal_uid"]: r for r in la.read_tsv(la.PROJECT / prespec["_coverage_source"])} \
        if "_coverage_source" in prespec else {}
    if not coverage:
        cov_path = sorted((ROOT.parent).glob("p6b-ancestry-*/tables/ancestry_portability.tsv"))[-1]
        coverage = {r["signal_uid"]: r for r in la.read_tsv(cov_path)}
        la.log(f"P6d: coverage states from {cov_path}")

    catalog = {r["signal_uid"]: r for r in la.read_tsv(TRACK0 / "gene_catalog_atlas_sidecar.tsv")}

    anchors = {}
    apath = sorted((ROOT.parent).glob("p2-anchors-*/tables/anchor_panel.tsv"))
    if apath:
        for r in la.read_tsv(apath[-1]):
            if r.get("variant_uid"):
                anchors[r["variant_uid"]] = r

    dbsnp_uid = {}
    dpath = sorted((ROOT.parent).glob("p0-dbsnp-orientation-*/tables/dbsnp_indel_orientation.tsv"))
    dpath = [q for q in dpath if not (q.parent.parent / "SUPERSEDED.txt").exists()]
    if dpath:
        for r in la.read_tsv(dpath[-1]):
            if r["resolved_variant_uid"]:
                dbsnp_uid[r["source_variant_id"]] = r["resolved_variant_uid"]
        la.log(f"P6d: {len(dbsnp_uid)} dbSNP-resolved orientations from {dpath[-1]}")

    repair = {}
    ppath = sorted((ROOT.parent).glob("p0-indel-uid-repair-*/tables/deferred_indel_uid_repair.tsv"))
    ppath = [q for q in ppath if not (q.parent.parent / "SUPERSEDED.txt").exists()]
    if ppath:
        for r in la.read_tsv(ppath[-1]):
            if r["repair_state"] == "recovered":
                repair[r["source_variant_id"]] = r["recovered_variant_uid"]
        la.log(f"P6d: {len(repair)} recovered hg38 identities from {ppath[-1]}")

    rescue = {}
    rpath = sorted((ROOT.parent).glob("p6f-indel-rescue-*/tables/indel_rescue_effects.tsv"))
    rpath = [q for q in rpath if not (q.parent.parent / "SUPERSEDED.txt").exists()]
    if rpath:
        rescue = rescue_index(la.read_tsv(rpath[-1]))
        la.log(f"P6d: model-API rescue rows under {len(rescue)} source ids from {rpath[-1]}")

    rows = []
    for sig_uid, sel in selected.items():
        s = signals[sig_uid]
        cov = coverage.get(sig_uid)
        if cov is None:
            raise la.ContractError(f"no coverage state for {sig_uid}")
        cov_state = cov["atlas_coverage_state"]
        for key in sorted(sel, key=lambda k: (-per_signal[sig_uid][k], k)):
            rk = ranks[sig_uid][key]
            mt = meta[key]
            uid, src = mt["variant_uid"], mt["source_variant_id"]
            rec = repair.get(src, "")
            dbu = dbsnp_uid.get(src, "")
            ref, alt = alleles_of(uid, src, rec, dbu)
            av = served.get(uid) if uid else None
            is_served = av is not None and av.get("AVI_SCORE") == "1"
            ch = channels.get(uid, {}) if is_served else {}
            tgt_ens = s.get("ensembl") or ""
            gene_rna = rna.get(uid, {}) if is_served else {}
            gene_spl = splice.get(uid, {}) if is_served else {}
            strongest = max(gene_rna.items(), key=lambda kv: abs(kv[1])) if gene_rna else None
            pk = peaks.get((sig_uid, uid), []) if uid else []
            dr = directions.get((sig_uid, uid), {}) if uid else {}
            res, res_state = rescue_for_row(rescue, src, ref, alt)
            cat = catalog.get(sig_uid, {})
            anc = anchors.get(uid, {}) if uid else {}
            readable, why = row_readable(cov_state, is_served)
            row = {
                "signal_uid": sig_uid, "universe": s["universe"], "posterior_definition": s["posterior_definition"],
                "gwas_name": s["gwas_name"], "trait": s["trait"], "ancestry": s["ancestry"],
                "gene": s["gene"], "ensembl": tgt_ens, "analysis_block": s["analysis_block"],
                "variant_uid": uid, "hg38_uid_recovered": rec, "hg38_uid_dbsnp": dbu,
                "orientation_source": orientation_source(dbu, rec),
                "source_variant_id": src, "variant_key": key,
                "mapping_status": mt["mapping_status"], "exclusion_reason": mt["exclusion_reason"],
                "variant_class": variant_class(ref, alt),
                "weight": repr(rk["weight"]), "weight_rank": rk["rank"],
                "cumulative_share_at_rank": repr(rk["cumulative_share"]),
                "atlas_served": str(bool(is_served)),
                "atac_liver_quantile": _fmt(ch.get("atac_liver_quantile")),
                "dnase_liver_quantile": _fmt(ch.get("dnase_liver_quantile")),
                "h3k27ac_liver_quantile": _fmt(ch.get("h3k27ac_liver_quantile")),
                "rna_target_gene_quantile": _fmt(gene_rna.get(tgt_ens) if tgt_ens else None),
                "rna_strongest_gene": (strongest[0] if strongest else ""),
                "rna_strongest_abs_quantile": _fmt(abs(strongest[1]) if strongest else None),
                "n_genes_rna_scored": len(gene_rna),
                "splice_site_usage_liver_quantile": _fmt(gene_spl.get(tgt_ens) if tgt_ens else
                                                         (max(gene_spl.values(), key=abs) if gene_spl else None)),
                "avi_raw": _fmt(ch.get("avi_raw")), "avi_quantile": _fmt(ch.get("avi_quantile")),
                "measured_lineages": ";".join(sorted({p["lineage"] for p in pk})),
                "measured_peak": ";".join(sorted({p["peak"] for p in pk if p["peak"]})),
                "measured_da_evidence_state": ";".join(sorted({p["evidence_state"] for p in pk if p["evidence_state"]})),
                "measured_atac_logFC_gse244832": ";".join(p["logFC_gse244832"] for p in pk if p["logFC_gse244832"]),
                "measured_atac_logFC_gse281367": ";".join(p["logFC_gse281367"] for p in pk if p["logFC_gse281367"]),
                "gwas_beta_allele1": dr.get("gwas_beta_allele1", ""),
                "eqtl_beta_allele1": dr.get("eqtl_beta_allele1", ""),
                "pred_allele1_liver_rna": dr.get("pred_allele1_liver_rna", ""),
                "measured_direction": dr.get("measured_direction", ""),
                "predicted_direction": dr.get("predicted_direction", ""),
                "direction_concordant": dr.get("concordant", ""),
                "ancestry_credible_sets": ";".join(sorted(ancestries.get(key, set()))),
                "n_ancestries": len(ancestries.get(key, set())),
                "in_eur_credible_set": str("EUR" in ancestries.get(key, set())),
                "anchor_rsid": anc.get("rsid", ""), "anchor_class": anc.get("class", ""),
                "anchor_expected_channel": anc.get("expected_channel", ""), "anchor_verdict": anc.get("verdict", ""),
                "signal_atlas_coverage_state": cov_state,
                "signal_queried_share": cov["queried_share"],
                "signal_excluded_share_indel_deferred": cov["excluded_share_indel_deferred"],
                "row_readable": readable, "row_not_readable_because": why,
                "model_api_rna_log2": res.get("rna_log2", ""), "model_api_splice_log2": res.get("splice_log2", ""),
                "model_api_atac_log2": res.get("atac_log2", ""), "model_api_dnase_log2": res.get("dnase_log2", ""),
                "model_api_h3k27ac_log2": res.get("h3k27ac_log2", ""),
                "model_api_state": res_state,
                "next_experiment_rule_id": cat.get("next_experiment_rule_id", ""),
                "atlas_extension_state": cat.get("atlas_extension_state", ""),
                "candidate_mechanism": cat.get("candidate_mechanism", ""),
            }
            guard_unserved_row(row)
            rows.append(row)

    # G1: every direct signal present.
    if len({r["signal_uid"] for r in rows}) != len(signals):
        raise la.ContractError("a direct signal is missing from the dossier")
    # G4: the two prediction services never share a column.
    if set(ATLAS_QUANTILE_COLUMNS) & set(MODEL_API_COLUMNS):
        raise la.ContractError("an Atlas quantile column and a model-API column share a name")
    # G5: the exported weights come from the same vector the coverage table summed.
    for sig_uid, w in per_signal.items():
        tot = float(coverage[sig_uid]["total_mass"])
        if abs(sum(w.values()) - tot) > 1e-6:
            raise la.ContractError(f"{sig_uid}: weight vector sums to {sum(w.values())!r}, coverage table says {tot!r}")

    rows.sort(key=lambda r: (r["universe"], r["signal_uid"], int(r["weight_rank"])))
    la.write_tsv_once(TABLES / "clinical_variant_dossier.tsv", rows, list(rows[0].keys()))

    n_peak_data = sum(1 for r in rows if r["measured_lineages"] or r["measured_peak"])
    n_in_peak = sum(1 for r in rows if r["measured_peak"])
    top_rows = [r for r in rows if int(r["weight_rank"]) == 1]
    summary = {
        "prespec_sha256": la.sha256_file(PRESPEC),
        "n_signals": len(signals), "n_rows": len(rows),
        "n_distinct_posterior_entries": len(keys), "n_with_an_hg38_uid": len(uids),
        "weight_floor": WEIGHT_FLOOR,
        "adult_liver_track_panel": {**{k: v for k, v in channels["__tracks__"].items()},
                                    "rna": gene_tracks["rna"], "splice_site_usage": gene_tracks["splice"]},
        "why_this_does_not_match_liver_summaries": ("liver_summaries pools primary liver with the CL:0000182 "
                                                    "hepatocyte tracks, which are EMBRYONIC in this Atlas build; "
                                                    "the dossier uses adult liver only"),
        "rows_by_variant_class": {c: sum(r["variant_class"] == c for r in rows)
                                  for c in ("snv", "insertion", "deletion", "mnv")},
        "rows_by_coverage_state": {s: sum(r["signal_atlas_coverage_state"] == s for r in rows)
                                   for s in sorted({r["signal_atlas_coverage_state"] for r in rows})},
        "signals_by_coverage_state": {s: len({r["signal_uid"] for r in rows if r["signal_atlas_coverage_state"] == s})
                                      for s in sorted({r["signal_atlas_coverage_state"] for r in rows})},
        "n_rows_unserved": sum(r["atlas_served"] == "False" for r in rows),
        "rows_by_exclusion_reason": {x: sum(r["exclusion_reason"] == x for r in rows)
                                     for x in sorted({r["exclusion_reason"] for r in rows})},
        "n_rows_with_a_recovered_hg38_identity": sum(1 for r in rows if r["hg38_uid_recovered"]),
        "rows_by_orientation_source": {s: sum(r["orientation_source"] == s for r in rows)
                                       for s in sorted({r["orientation_source"] for r in rows})},
        "n_rows_where_dbsnp_flips_the_reference_only_call": sum(
            1 for r in rows if r["hg38_uid_dbsnp"] and r["hg38_uid_recovered"]
            and r["hg38_uid_dbsnp"] != r["hg38_uid_recovered"]),
        "recovered_identity_note": ("step 67 recovers the hg38 reference orientation of a deferred indel; it "
                                    "does NOT make the Atlas point-query API accept it, so atlas_served "
                                    "stays False and every Atlas quantile stays empty"),
        "n_rows_readable": sum(r["row_readable"] == "True" for r in rows),
        "R3_every_non_covered_signal_shows_an_unserved_row": all(
            any(r["atlas_served"] == "False" for r in rows if r["signal_uid"] == s)
            for s in {r["signal_uid"] for r in rows if r["signal_atlas_coverage_state"] != "covered"}),
        "R4_signals_whose_top_variant_the_atlas_cannot_serve": sorted(
            {r["signal_uid"] for r in top_rows if r["atlas_served"] == "False"}),
        "n_signals_top_variant_unserved": len({r["signal_uid"] for r in top_rows if r["atlas_served"] == "False"}),
        "n_rows_rescued_by_model_api": sum(1 for r in rows if r["model_api_splice_log2"] != ""),
        "n_rows_model_api_withheld_orientation_differs": sum(1 for r in rows
                                                             if r["model_api_state"] == "orientation_differs"),
        "n_rows_with_an_anchor": sum(1 for r in rows if r["anchor_rsid"]),
        "measured_peak_share": measured_peak_share(n_in_peak, n_peak_data, len(rows)),
        "direction_concordance": direction_concordance(rows),
        "ancestry_membership": {"n_rows_multi_ancestry": sum(int(r["n_ancestries"]) > 1 for r in rows),
                                "n_rows_eur_only": sum(r["ancestry_credible_sets"] == "EUR" for r in rows),
                                "membership_floor": WEIGHT_FLOOR},
        "claim_boundary": prespec["claim_boundary"],
    }
    json.dump(summary, (TABLES / "clinical_variant_dossier_summary.json").open("w"), indent=1, default=float)
    la.log(f"P6d: wrote {len(rows)} rows; readable {summary['n_rows_readable']}; "
           f"top-variant unserved in {summary['n_signals_top_variant_unserved']} signals")


if __name__ == "__main__":
    main()
