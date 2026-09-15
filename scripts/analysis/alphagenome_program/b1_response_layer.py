#!/usr/bin/env python3
"""B1/B2 (spec section 4): the locus-by-context-by-readout research layer over all 989 Atlas signals.

B1 extends the P6d clinical dossier's per-variant row and its readable-row rule from the 202 direct-MASLD
signals to the whole eligible universe (A_direct 29, B_direct 173, C_enzyme 787). One row per
(signal_uid, variant_uid) for every variant carrying at least 0.01 of its signal's posterior, plus each
signal's top-weight variant whatever its weight.

Every row travels with the parent signal's Atlas coverage state and its own served flag, because a variant
the Atlas never scored has empty prediction columns that must not be read as a null effect, and a signal
whose posterior mass was largely deferred has signal-level columns computed over the mass that survived.

B2 rebuilds the per-context MPRA differential-allelic-variant calls from the source reproduction and writes
the re-key that corrects the TableS3/TableS4 context swap in the p3a deposit. The p3a deposit is NOT edited.

Nothing here adopts a Resource number, gene class, figure or claim. Atlas quantiles are predictions, not
measurements, and no molecular explanation class in the signal summary is assigned from an Atlas score.

Usage: b1_response_layer.py <output_directory>
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import pathlib
import random
import re
import statistics
import subprocess
import sys
from collections import defaultdict

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = PROJECT / "GWAS/finemapping/results/alphagenome_atlas"
TRACK0 = ATLAS / "run-20260909T153939Z/tables"
RECORDS = TRACK0 / "prediction_records"

# Deposits are pinned by name, never globbed to "latest": another session is writing into this tree and an
# in-flight deposit directory can exist with an empty tables/.
COVERAGE_TSV = ATLAS / "p6b-ancestry-20260913T153316Z/tables/ancestry_portability.tsv"
COVERAGE_FLAGS_TSV = ATLAS / "p6b-ancestry-20260913T153316Z/tables/atlas_coverage_flags.tsv"
ASE_281367 = ATLAS / "p3-ase-20260909T192926Z/tables/allelic_sites.tsv"
ASE_244832 = ATLAS / "p3-ase-20260909T192926Z/tables/gse244832_allelic_sites.tsv"
DOSSIER_P6D = ATLAS / "p6d-dossier-20260914T173952Z/tables/clinical_variant_dossier.tsv"
RESCUE_TSV = ATLAS / "p6f-indel-rescue-20260914T135052Z/tables/indel_rescue_effects.tsv"
REPAIR_TSV = ATLAS / "p0-indel-uid-repair-20260913T193750Z/tables/deferred_indel_uid_repair.tsv"
# Both dbSNP orientation deposits that hold a table are stamped SUPERSEDED and the newest deposit directory
# was empty (a rerun in flight) when this ran; see DEFECTS.md entry B1-D3. They never conflict where both
# resolve a variant (0 of 41,665), but neither covers the other: the later deposit fixed the position match
# and lost 35,115 calls to an unresolvable 'both_in_dbsnp' state for want of a frequency tie-break, while the
# earlier one resolved those and wrongly called 26,068 variants absent. So both are read, the later one wins
# where it speaks, and the deposit each row's call came from is printed beside it.
DBSNP_TSV = ATLAS / "p0-dbsnp-orientation-20260914T174542Z/tables/dbsnp_indel_orientation.tsv"
DBSNP_TSV_EARLIER = ATLAS / "p0-dbsnp-orientation-20260914T121856Z/tables/dbsnp_indel_orientation.tsv"
ANCHORS_TSV = ATLAS / "p2-anchors-20260909T191328Z/tables/anchor_panel.tsv"
MPRA_SRC = PROJECT / "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction"
MPRA_P3A = ATLAS / "p3a-benchmarks-20260909T190214Z/tables/mpra_dav_calls_S2_S5.tsv"

SEED = 20260914
WEIGHT_FLOOR = 0.01
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

SOURCE_ID_RE = re.compile(r"^(?P<asm>[A-Za-z0-9]+):(?:chr)?(?P<chrom>[0-9]{1,2}|[XYM]|MT):(?P<pos>\d+):"
                          r"(?P<ref>[ACGTNacgtn]+):(?P<alt>[ACGTNacgtn]+)$")


class ContractError(RuntimeError):
    """Raised when a table violates a prespecified invariant."""


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def sha256_file(path: pathlib.Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as h:
        for b in iter(lambda: h.read(1 << 20), b""):
            d.update(b)
    return d.hexdigest()


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_tsv(path) -> list:
    with open_text(path) as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv_once(path, rows, columns) -> int:
    path = pathlib.Path(path)
    if path.exists():
        raise ContractError(f"refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(path, "wt") as h:
        w = csv.DictWriter(h, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})
            n += 1
    return n


# ---------------------------------------------------------------- identity helpers
# Same semantics as scripts/analysis/alphagenome_atlas/lib_atlas.py; reimplemented here so this package does
# not depend on a module another session is actively editing. Both files' sha256 go into MANIFEST.tsv and the
# equivalence is checked by reproducing the P6d dossier row for row (check R1).

def normalise_source_id(source_id: str) -> str:
    m = SOURCE_ID_RE.match(str(source_id).strip())
    if m is None:
        raise ContractError(f"source id is not assembly:chrom:pos:ref:alt: {source_id!r}")
    return f"{m['asm']}:chr{m['chrom']}:{int(m['pos'])}:{m['ref'].upper()}:{m['alt'].upper()}"


def posterior_key(variant_uid: str, source_variant_id: str) -> str:
    """NEVER key on variant_uid alone: a deferred variant has an empty uid and they would all collapse."""
    uid = str(variant_uid or "").strip()
    if uid:
        return uid
    src = str(source_variant_id or "").strip()
    if not src:
        raise ContractError("a posterior entry with neither an hg38 uid nor a source id cannot be identified")
    return normalise_source_id(src)


def variant_class(ref: str, alt: str) -> str:
    r, a = str(ref).upper(), str(alt).upper()
    if len(r) == len(a):
        return "snv" if len(r) == 1 else "mnv"
    return "insertion" if len(a) > len(r) else "deletion"


def orientation_source(dbsnp_uid: str, recovered_uid: str) -> str:
    if str(dbsnp_uid or "").strip():
        return "dbsnp"
    if str(recovered_uid or "").strip():
        return "reference_or_source_order"
    return ""


def alleles_of(variant_uid, source_variant_id, recovered_uid="", dbsnp_uid=""):
    for uid in (str(variant_uid or "").strip(), str(dbsnp_uid or "").strip(), str(recovered_uid or "").strip()):
        if uid:
            parts = uid.split(":")
            if len(parts) != 4:
                raise ContractError(f"variant uid is not chrom:pos:ref:alt: {uid!r}")
            return parts[2], parts[3]
    parts = normalise_source_id(source_variant_id).split(":")
    return parts[3], parts[4]


def block_1mb(variant_uid: str) -> str:
    if not variant_uid:
        return ""
    chrom, pos = variant_uid.split(":")[0], int(variant_uid.split(":")[1])
    return f"{chrom}~{pos // 1_000_000}"


def rank_and_cumulative(weights: dict) -> dict:
    if not weights:
        raise ContractError("a signal with no posterior weights cannot be ranked")
    order = sorted(weights.items(), key=lambda kv: (-float(kv[1]), kv[0]))
    out, run = {}, 0.0
    for i, (uid, w) in enumerate(order, start=1):
        run += float(w)
        out[uid] = {"rank": i, "weight": float(w), "cumulative_share": run}
    return out


def select_variants(weights: dict, floor: float = WEIGHT_FLOOR) -> set:
    if not weights:
        raise ContractError("a signal with no posterior weights has nothing to export")
    sel = {u for u, w in weights.items() if float(w) >= floor}
    sel.add(sorted(weights.items(), key=lambda kv: (-float(kv[1]), kv[0]))[0][0])
    return sel


def row_readable(coverage_state: str, served: bool):
    if not served:
        return "False", "variant not served by the Atlas point-query API (indel or unmapped)"
    if coverage_state != "covered":
        return "False", (f"parent signal coverage state is {coverage_state}; its signal-level columns rest "
                         "on incomplete mass")
    return "True", ""


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def _median(vals):
    v = [x for x in vals if x is not None]
    return statistics.median(v) if v else None


def _fmt(v):
    return "" if v is None else repr(float(v))


# ---------------------------------------------------------------- loaders
def load_signals() -> dict:
    return {r["signal_uid"]: r for r in read_tsv(TRACK0 / "eligible_signals.tsv")}


def load_weights(signals: dict):
    """Posterior weights, per-row metadata and credible-set ancestry membership.

    meta is keyed on (signal_uid, posterior key), NOT on the key alone. The same hg38 variant appears in
    several signals and its `source_variant_id` carries that study's own allele ORDER, so a dict keyed on the
    variant identity alone hands a row the allele order of whichever signal was read last. That id is the
    join key for the dbSNP orientation, the deferred-indel repair and the model-API rescue, so the wrong
    order is a wrong join, not a cosmetic difference.

    Ancestry membership is collected twice: over every eligible signal, and over the direct-MASLD signals
    only, because the second is the definition the P6d dossier used and the two must not be confused.
    """
    per_signal, meta = defaultdict(dict), {}
    anc_all, anc_direct = defaultdict(set), defaultdict(set)
    with open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            sig = signals.get(r["signal_uid"])
            if sig is None:
                continue
            w = float(r["weight"])
            key = posterior_key(r["variant_uid"], r["source_variant_id"])
            if key in per_signal[r["signal_uid"]]:
                raise ContractError(f"{r['signal_uid']}: two posterior entries share the identity {key}")
            per_signal[r["signal_uid"]][key] = w
            meta[(r["signal_uid"], key)] = {
                "variant_uid": r["variant_uid"], "source_variant_id": r["source_variant_id"],
                "mapping_status": r["mapping_status"], "exclusion_reason": r["exclusion_reason"]}
            if w >= WEIGHT_FLOOR and sig["ancestry"]:
                anc_all[key].add(sig["ancestry"])
                if sig["universe"] in ("A_direct", "B_direct"):
                    anc_direct[key].add(sig["ancestry"])
    return per_signal, meta, anc_all, anc_direct


def load_atlas_channels(uids: set) -> dict:
    import pyarrow.parquet as pq
    out = defaultdict(dict)
    tracks = {}
    for col, (fname, keep) in TRACK_PANEL.items():
        pf = pq.ParquetFile(RECORDS / fname)
        cols = [c for c in pf.schema.names if keep(c)]
        if not cols:
            raise ContractError(f"no adult-liver track column matched for {col} in {fname}")
        tracks[col] = cols
        for batch in pf.iter_batches(batch_size=100_000, columns=["variant_uid"] + cols):
            d = batch.to_pydict()
            for i, uid in enumerate(d["variant_uid"]):
                if uid in uids:
                    out[uid][col] = _median([_f(d[c][i]) for c in cols])
    pf = pq.ParquetFile(RECORDS / "AVI_SCORE.parquet")
    for batch in pf.iter_batches(batch_size=100_000,
                                 columns=["variant_uid", "raw|feature|AVI_SCORE", "q|feature|AVI_SCORE"]):
        d = batch.to_pydict()
        for i, uid in enumerate(d["variant_uid"]):
            if uid in uids:
                out[uid]["avi_raw"] = _f(d["raw|feature|AVI_SCORE"][i])
                out[uid]["avi_quantile"] = _f(d["q|feature|AVI_SCORE"][i])
    return out, tracks


def load_gene_channels(uids: set):
    import pyarrow.parquet as pq
    rna, splice, tracks = defaultdict(dict), defaultdict(dict), {}
    for fname, sink, key in (("RNA_SEQ.parquet", rna, "rna"), ("SPLICE_SITE_USAGE.parquet", splice, "splice")):
        pf = pq.ParquetFile(RECORDS / fname)
        cols = [c for c in pf.schema.names if c.startswith(ADULT_LIVER_PREFIX)]
        if not cols:
            raise ContractError(f"no adult-liver track column in {fname}")
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


def load_allelic_sites():
    """Measured allelic imbalance, keyed on the hg38 uid. The sign refers to ALT of that uid."""
    out = defaultdict(dict)
    for label, path in (("gse281367", ASE_281367), ("gse244832", ASE_244832)):
        for r in read_tsv(path):
            out[r["uid"]][label] = r
    return out


# ---------------------------------------------------------------- B2: MPRA context re-key
def mpra_rekey(outdir: pathlib.Path) -> dict:
    """Rebuild the per-context DAV calls from the source reproduction and prove the p3a sheet/context swap.

    Proof is by content, not by row count: each p3a sheet's element set and its log2 fold changes are
    compared against all four reproduced contexts.
    """
    contexts = [("HepG2", "control", "HepG2_control"), ("HepG2", "PAOA", "HepG2_PAOA"),
                ("LX2", "control", "LX2_control"), ("LX2", "TGFb", "LX2_TGFb")]
    official = {}
    for cl, ctx, name in contexts:
        rows = read_tsv(MPRA_SRC / f"{cl}.{ctx}.official_comparison.tsv.gz")
        official[name] = {r["element_id"]: r for r in rows}
    qc = {}
    for cl in ("HepG2", "LX2"):
        for e in json.load((MPRA_SRC / f"{cl}.source_reproduction_qc.json").open()):
            qc[f"{cl}_{'control' if e['context'] == 'control' else e['context']}"] = e

    p3a = defaultdict(list)
    for r in read_tsv(MPRA_P3A):
        p3a[r["table"]].append(r)

    out_rows, mapping = [], {}
    for sheet in sorted(p3a):
        sub = p3a[sheet]
        labelled = sub[0]["context"]
        keys = {r["coordination"] for r in sub}
        best, best_ov = None, -1
        for name, d in official.items():
            ov = len(keys & set(d))
            if ov > best_ov:
                best, best_ov = name, ov
        d = official[best]
        common = sorted(keys & set(d))
        maxdiff = max((abs(float(r["log2(Fold change)"]) - float(d[r["coordination"]]["official_log2FC"]))
                       for r in sub if r["coordination"] in d), default=None)
        mapping[sheet] = best
        out_rows.append({
            "sheet": sheet,
            "p3a_context_label": labelled,
            "p3a_n_rows": len(sub),
            "rebuilt_context": best,
            "rebuilt_n_official": qc[best]["n_official"],
            "element_overlap_with_rebuilt": best_ov,
            "element_overlap_share": repr(best_ov / len(d)) if d else "",
            "max_abs_log2FC_difference_on_shared_elements": ("" if maxdiff is None else repr(maxdiff)),
            "n_shared_elements_compared": len(common),
            "label_is_correct": str(labelled.lower().replace("_ctrl", "_control") == best.lower()),
            "effect_direction": qc[best]["effect_direction"],
        })
    write_tsv_once(outdir / "mpra_context_rekey.tsv", out_rows, list(out_rows[0].keys()))
    return {"mapping": mapping, "rows": out_rows,
            "qc_n_official": {k: v["n_official"] for k, v in qc.items()}}


# ---------------------------------------------------------------- main
def main() -> None:
    outdir = pathlib.Path(sys.argv[1]).resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise ContractError(f"refusing to write into a non-empty directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    random.seed(SEED)

    signals = load_signals()
    log(f"B1: {len(signals)} eligible signals")
    per_signal, meta, anc_all, anc_direct = load_weights(signals)
    if len(per_signal) != len(signals):
        raise ContractError(f"weights cover {len(per_signal)} of {len(signals)} signals")

    ranks = {s: rank_and_cumulative(w) for s, w in per_signal.items()}
    selected = {s: select_variants(w) for s, w in per_signal.items()}
    keys = {k for sel in selected.values() for k in sel}
    uids = {meta[(sig, k)]["variant_uid"] for sig, sel in selected.items() for k in sel
            if meta[(sig, k)]["variant_uid"]}
    log(f"B1: {sum(len(v) for v in selected.values())} rows, {len(keys)} distinct entries, {len(uids)} with a uid")

    served = {}
    with open_text(TRACK0 / "atlas_availability.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["variant_uid"] in uids:
                served[r["variant_uid"]] = r

    channels, track_panel = load_atlas_channels(uids)
    rna, splice, gene_tracks = load_gene_channels(uids)

    peaks = defaultdict(list)
    with open_text(TRACK0 / "variant_peak_overlaps.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["signal_uid"] in signals and r["variant_uid"] in uids:
                peaks[(r["signal_uid"], r["variant_uid"])].append(r)

    directions = {}
    with open_text(TRACK0 / "direction_variant_level.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["signal_uid"] in signals and r["variant_uid"] in uids:
                directions[(r["signal_uid"], r["variant_uid"])] = r

    coverage = {r["signal_uid"]: r for r in read_tsv(COVERAGE_TSV)}
    flags = {r["signal_uid"]: r for r in read_tsv(COVERAGE_FLAGS_TSV)}
    catalog = {r["signal_uid"]: r for r in read_tsv(TRACK0 / "gene_catalog_atlas_sidecar.tsv")}
    anchors = {r["variant_uid"]: r for r in read_tsv(ANCHORS_TSV) if r.get("variant_uid")}
    ase = load_allelic_sites()

    dbsnp_late = {r["source_variant_id"]: r["resolved_variant_uid"]
                  for r in read_tsv(DBSNP_TSV) if r["resolved_variant_uid"]}
    dbsnp_early = {r["source_variant_id"]: r["resolved_variant_uid"]
                   for r in read_tsv(DBSNP_TSV_EARLIER) if r["resolved_variant_uid"]}
    conflicts = [k for k in set(dbsnp_late) & set(dbsnp_early) if dbsnp_late[k] != dbsnp_early[k]]
    if conflicts:
        raise ContractError(f"the two dbSNP deposits disagree on {len(conflicts)} variants, e.g. {conflicts[0]}; "
                            "a union is only defensible while they never conflict")
    dbsnp_uid = {**dbsnp_early, **dbsnp_late}
    dbsnp_from = {k: ("dbsnp_20260914T174542Z" if k in dbsnp_late else "dbsnp_20260914T121856Z")
                  for k in dbsnp_uid}
    repair = {r["source_variant_id"]: r["recovered_variant_uid"]
              for r in read_tsv(REPAIR_TSV) if r["repair_state"] == "recovered"}
    rescue = {normalise_source_id(r["source_variant_id"]): r for r in read_tsv(RESCUE_TSV)
              if r.get("arm") == "indel" and r.get("state") == "scored"}
    log(f"B1: {len(dbsnp_uid)} dbSNP orientations ({len(dbsnp_late)} later deposit, {len(dbsnp_early)} earlier, {len(set(dbsnp_late) & set(dbsnp_early))} in both, 0 conflicts), {len(repair)} recovered identities, {len(rescue)} rescued indels")

    rows = []
    for sig_uid, sel in selected.items():
        s = signals[sig_uid]
        cov = coverage.get(sig_uid)
        if cov is None:
            raise ContractError(f"no coverage state for {sig_uid}")
        cov_state = cov["atlas_coverage_state"]
        for key in sorted(sel, key=lambda k: (-per_signal[sig_uid][k], k)):
            rk = ranks[sig_uid][key]
            mt = meta[(sig_uid, key)]
            uid, src = mt["variant_uid"], mt["source_variant_id"]
            rec, dbu = repair.get(src, ""), dbsnp_uid.get(src, "")
            ref, alt = alleles_of(uid, src, rec, dbu)
            av = served.get(uid) if uid else None
            is_served = av is not None and av.get("AVI_SCORE") == "1"
            ch = channels.get(uid, {}) if is_served else {}
            tgt = s.get("ensembl") or ""
            gene_rna = rna.get(uid, {}) if is_served else {}
            gene_spl = splice.get(uid, {}) if is_served else {}
            strongest = max(gene_rna.items(), key=lambda kv: abs(kv[1])) if gene_rna else None
            pk = peaks.get((sig_uid, uid), []) if uid else []
            dr = directions.get((sig_uid, uid), {}) if uid else {}
            res = rescue.get(normalise_source_id(src), {})
            cat = catalog.get(sig_uid, {})
            anc = anchors.get(uid, {}) if uid else {}
            a1 = ase.get(uid, {}).get("gse281367", {}) if uid else {}
            a2 = ase.get(uid, {}).get("gse244832", {}) if uid else {}
            readable, why = row_readable(cov_state, is_served)
            row = {
                "signal_uid": sig_uid, "universe": s["universe"], "trait_class": s["trait_class"],
                "posterior_definition": s["posterior_definition"], "gwas_name": s["gwas_name"],
                "trait": s["trait"], "ancestry": s["ancestry"], "gene": s["gene"], "ensembl": tgt,
                "analysis_block": s["analysis_block"],
                "variant_uid": uid, "block_1mb": block_1mb(uid),
                "hg38_uid_recovered": rec, "hg38_uid_dbsnp": dbu,
                "orientation_source": orientation_source(dbu, rec),
                "dbsnp_deposit": (dbsnp_from.get(src, "") if dbu else ""),
                "source_variant_id": src, "variant_key": key,
                "mapping_status": mt["mapping_status"], "exclusion_reason": mt["exclusion_reason"],
                "variant_class": variant_class(ref, alt), "hg38_ref": ref, "hg38_alt": alt,
                "weight": repr(rk["weight"]), "weight_rank": rk["rank"],
                "cumulative_share_at_rank": repr(rk["cumulative_share"]),
                "atlas_served": str(bool(is_served)),
                "atac_liver_quantile": _fmt(ch.get("atac_liver_quantile")),
                "dnase_liver_quantile": _fmt(ch.get("dnase_liver_quantile")),
                "h3k27ac_liver_quantile": _fmt(ch.get("h3k27ac_liver_quantile")),
                "rna_target_gene_quantile": _fmt(gene_rna.get(tgt) if tgt else None),
                "rna_strongest_gene": (strongest[0] if strongest else ""),
                "rna_strongest_abs_quantile": _fmt(abs(strongest[1]) if strongest else None),
                "n_genes_rna_scored": len(gene_rna),
                "splice_site_usage_liver_quantile": _fmt(gene_spl.get(tgt) if tgt else
                                                         (max(gene_spl.values(), key=abs) if gene_spl else None)),
                "avi_raw": _fmt(ch.get("avi_raw")), "avi_quantile": _fmt(ch.get("avi_quantile")),
                "measured_lineages": ";".join(sorted({p["lineage"] for p in pk})),
                "measured_peak": ";".join(sorted({p["peak"] for p in pk if p["peak"]})),
                "measured_da_evidence_state": ";".join(sorted({p["evidence_state"] for p in pk if p["evidence_state"]})),
                "measured_atac_logFC_gse244832": ";".join(p["logFC_gse244832"] for p in pk if p["logFC_gse244832"]),
                "measured_atac_logFC_gse281367": ";".join(p["logFC_gse281367"] for p in pk if p["logFC_gse281367"]),
                "ase_gse281367_n_het_donors": a1.get("n_het_donors", ""),
                "ase_gse281367_mean_log2_alt_over_ref": a1.get("mean_log2_alt_over_ref", ""),
                "ase_gse281367_se_log2": a1.get("se_log2", ""),
                "ase_gse244832_n_het_donors": a2.get("n_het_donors", ""),
                "ase_gse244832_mean_log2_alt_over_ref": a2.get("mean_log2_alt_over_ref", ""),
                "ase_gse244832_se_log2": a2.get("se_log2", ""),
                "ase_sign_refers_to": ("ALT of variant_uid (hg38 reference orientation)" if (a1 or a2) else ""),
                "gwas_beta_allele1": dr.get("gwas_beta_allele1", ""),
                "eqtl_beta_allele1": dr.get("eqtl_beta_allele1", ""),
                "pred_allele1_liver_rna": dr.get("pred_allele1_liver_rna", ""),
                "direction_sign_refers_to": ("allele1 of source_variant_id (GWAS effect allele after "
                                             "harmonisation)" if dr.get("eqtl_beta_allele1") else ""),
                "measured_direction": dr.get("measured_direction", ""),
                "predicted_direction": dr.get("predicted_direction", ""),
                "direction_concordant": dr.get("concordant", ""),
                "ancestry_credible_sets": ";".join(sorted(anc_all.get(key, set()))),
                "n_ancestries": len(anc_all.get(key, set())),
                "in_eur_credible_set": str("EUR" in anc_all.get(key, set())),
                "ancestry_credible_sets_direct_only": ";".join(sorted(anc_direct.get(key, set()))),
                "n_ancestries_direct_only": len(anc_direct.get(key, set())),
                "in_eur_credible_set_direct_only": str("EUR" in anc_direct.get(key, set())),
                "anchor_rsid": anc.get("rsid", ""), "anchor_class": anc.get("class", ""),
                "anchor_expected_channel": anc.get("expected_channel", ""), "anchor_verdict": anc.get("verdict", ""),
                "signal_atlas_coverage_state": cov_state,
                "signal_queried_share": cov["queried_share"],
                "signal_excluded_share_indel_deferred": cov["excluded_share_indel_deferred"],
                "signal_in_coverage_flag_table": str(sig_uid in flags),
                "row_readable": readable, "row_not_readable_because": why,
                "model_api_rna_log2": res.get("rna_log2", ""), "model_api_splice_log2": res.get("splice_log2", ""),
                "model_api_atac_log2": res.get("atac_log2", ""), "model_api_dnase_log2": res.get("dnase_log2", ""),
                "model_api_h3k27ac_log2": res.get("h3k27ac_log2", ""),
                "model_api_sign_refers_to": ("log2(ALT/REF) in hg38 reference orientation, hosted model API, "
                                             "not comparable to an Atlas quantile" if res else ""),
                "model_api_orientation_source": res.get("orientation_source", ""),
                "model_api_orientation_ambiguous": res.get("orientation_ambiguous", ""),
                "next_experiment_rule_id": cat.get("next_experiment_rule_id", ""),
                "atlas_extension_state": cat.get("atlas_extension_state", ""),
                "candidate_mechanism_atlas_predicted": cat.get("candidate_mechanism", ""),
            }
            if row["atlas_served"] == "False":
                for c in ATLAS_QUANTILE_COLUMNS:
                    if str(row.get(c, "")).strip():
                        raise ContractError(f"unserved variant carries an Atlas quantile in {c}")
            rows.append(row)

    if len({r["signal_uid"] for r in rows}) != len(signals):
        raise ContractError("a signal is missing from the response layer")
    if set(ATLAS_QUANTILE_COLUMNS) & set(MODEL_API_COLUMNS):
        raise ContractError("an Atlas quantile column and a model-API column share a name")
    for sig_uid, w in per_signal.items():
        tot = float(coverage[sig_uid]["total_mass"])
        if abs(sum(w.values()) - tot) > 1e-6:
            raise ContractError(f"{sig_uid}: weight vector sums to {sum(w.values())!r}, coverage says {tot!r}")

    rows.sort(key=lambda r: (r["universe"], r["signal_uid"], int(r["weight_rank"])))
    write_tsv_once(outdir / "response_layer_variants.tsv.gz", rows, list(rows[0].keys()))
    log(f"B1: wrote {len(rows)} variant rows")

    # ------------------------------------------------------------ R1: reproduce the P6d dossier
    check_r1 = reproduce_p6d(rows)

    # ------------------------------------------------------------ per-signal summary
    by_sig = defaultdict(list)
    for r in rows:
        by_sig[r["signal_uid"]].append(r)
    sig_rows = [signal_summary(signals[s], coverage[s], catalog.get(s, {}), by_sig[s], per_signal[s])
                for s in sorted(by_sig)]
    write_tsv_once(outdir / "response_layer_signals.tsv", sig_rows, list(sig_rows[0].keys()))

    # ------------------------------------------------------------ coverage gate
    gate = []
    for s in sorted(signals):
        cov, sr = coverage[s], by_sig[s]
        gate.append({
            "signal_uid": s, "universe": signals[s]["universe"], "trait_class": signals[s]["trait_class"],
            "posterior_definition": signals[s]["posterior_definition"],
            "analysis_block": signals[s]["analysis_block"],
            "atlas_coverage_state": cov["atlas_coverage_state"],
            "queried_share": cov["queried_share"],
            "excluded_share_indel_deferred": cov["excluded_share_indel_deferred"],
            "excluded_share_liftover_failed": cov["excluded_share_liftover_failed"],
            "total_mass": cov["total_mass"], "n_cs_variants": cov["n_cs_variants"],
            "in_coverage_flag_table": str(s in flags),
            "n_exported_rows": len(sr),
            "n_rows_served": sum(r["atlas_served"] == "True" for r in sr),
            "n_rows_readable": sum(r["row_readable"] == "True" for r in sr),
            "top_variant_served": [r["atlas_served"] for r in sr if int(r["weight_rank"]) == 1][0],
            "signal_gate": ("readable" if cov["atlas_coverage_state"] == "covered"
                            and any(r["row_readable"] == "True" for r in sr) else "not_readable"),
        })
    write_tsv_once(outdir / "coverage_gate.tsv", gate, list(gate[0].keys()))

    # ------------------------------------------------------------ B2
    mp = mpra_rekey(outdir)

    summary = build_summary(rows, sig_rows, gate, signals, mp, check_r1, track_panel, gene_tracks)
    json.dump(summary, (outdir / "response_layer_summary.json").open("w"), indent=1, default=float)
    log("B1: summary written")
    return summary


def reproduce_p6d(rows) -> dict:
    """R1: the direct-universe subset of this layer must equal the P6d dossier column for column."""
    if not DOSSIER_P6D.exists():
        return {"ran": False, "reason": "P6d dossier absent"}
    ref = {(r["signal_uid"], r["variant_key"]): r for r in read_tsv(DOSSIER_P6D)}
    mine = {(r["signal_uid"], r["variant_key"]): r for r in rows
            if r["universe"] in ("A_direct", "B_direct")}
    # P6d counted credible-set ancestries over the direct signals only; this layer also reports the
    # all-universe count. The equivalence check compares like with like.
    alias = {"ancestry_credible_sets": "ancestry_credible_sets_direct_only",
             "n_ancestries": "n_ancestries_direct_only",
             "in_eur_credible_set": "in_eur_credible_set_direct_only"}
    mycols = set(next(iter(mine.values())).keys())
    shared_cols = [c for c in next(iter(ref.values())).keys() if alias.get(c, c) in mycols]
    diffs, diffs_non_snv, missing = defaultdict(int), defaultdict(int), 0
    examples = {}
    for k in set(ref) | set(mine):
        if k not in ref or k not in mine:
            missing += 1
            continue
        for c in shared_cols:
            a, b = str(ref[k].get(c, "")), str(mine[k].get(alias.get(c, c), ""))
            if a != b:
                diffs[c] += 1
                if mine[k]["variant_class"] != "snv":
                    diffs_non_snv[c] += 1
                examples.setdefault(c, {"key": list(k), "p6d": a, "b1": b})
    # A source_variant_id difference matters only when the two spellings are different VARIANTS. Separate the
    # two: same normalised id means the same variant written differently; a different one is an allele-order
    # or position difference and would send the dbSNP, repair and rescue joins to another row.
    sid_same_variant = sid_diff_variant = 0
    for k in set(ref) & set(mine):
        a, b = ref[k].get("source_variant_id", ""), mine[k].get("source_variant_id", "")
        if a != b:
            if normalise_source_id(a) == normalise_source_id(b):
                sid_same_variant += 1
            else:
                sid_diff_variant += 1
    return {"ran": True, "n_p6d_rows": len(ref), "n_b1_direct_rows": len(mine),
            "source_variant_id_differs_same_normalised_variant": sid_same_variant,
            "source_variant_id_differs_DIFFERENT_variant": sid_diff_variant,
            "n_rows_present_in_only_one": missing,
            "n_shared_columns": len(shared_cols), "columns_compared": shared_cols,
            "n_columns_differing": len(diffs),
            "differing_columns": dict(diffs),
            "differing_columns_rows_that_are_not_snv": dict(diffs_non_snv),
            "examples": examples,
            "column_aliases_used": alias,
            "identical": len(diffs) == 0 and missing == 0,
            "note": ("any difference should be confined to the orientation columns on non-SNV rows: P6d "
                     "globbed the then-newest dbSNP orientation deposit, this layer pins "
                     f"{DBSNP_TSV.parent.parent.name}. A difference on an SNV row is a defect, not a pin.")}


MEASURED_DA_STATES = {"supported", "source_dependent"}


def signal_summary(sig, cov, cat, srows, weights) -> dict:
    """Which molecular explanation classes a reader could DISTINGUISH at this locus.

    The rule is evidence-channel presence, never an Atlas score value. A class is `measured` only when a
    donor-level measurement or an independent sequence annotation speaks to it; `served_only` means the sole
    channel is a hosted Atlas prediction, which cannot by itself separate one explanation from another; and
    `no_channel` means this deposit holds no observation of that kind at this locus at all.
    """
    n = len(srows)
    served = [r for r in srows if r["atlas_served"] == "True"]
    readable = [r for r in srows if r["row_readable"] == "True"]
    in_peak = [r for r in srows if r["measured_peak"]]
    da = [r for r in srows if set(r["measured_da_evidence_state"].split(";")) & MEASURED_DA_STATES]
    ase_rows = [r for r in srows if r["ase_gse281367_mean_log2_alt_over_ref"] or r["ase_gse244832_mean_log2_alt_over_ref"]]
    eqtl_rows = [r for r in srows if r["eqtl_beta_allele1"]]
    spl_served = [r for r in served if r["splice_site_usage_liver_quantile"]]
    rescued = [r for r in srows if r["model_api_splice_log2"]]
    cons_cov = _f(cat.get("consequence_coverage")) or 0.0
    coding_mass = _f(cat.get("coding_mass")) or 0.0

    coding = ("measured" if cons_cov >= 0.5 else ("partial" if cons_cov > 0 else "no_channel"))
    local_reg = ("measured" if (da or ase_rows) else ("served_only" if (in_peak or served) else "no_channel"))
    rna_proc = ("served_only" if spl_served or rescued else "no_channel")
    long_range = "no_channel"
    distinguishable = [c for c, v in (("coding", coding), ("local_regulatory", local_reg),
                                      ("rna_processing", rna_proc), ("candidate_long_range", long_range))
                       if v == "measured"]
    return {
        "signal_uid": sig["signal_uid"], "universe": sig["universe"], "trait_class": sig["trait_class"],
        "trait": sig["trait"], "gwas_name": sig["gwas_name"], "gene": sig["gene"], "ensembl": sig["ensembl"],
        "posterior_definition": sig["posterior_definition"], "analysis_block": sig["analysis_block"],
        "atlas_coverage_state": cov["atlas_coverage_state"], "queried_share": cov["queried_share"],
        "total_mass": cov["total_mass"],
        "n_exported_rows": n, "n_rows_served": len(served), "n_rows_readable": len(readable),
        "n_distinct_1mb_blocks": len({r["block_1mb"] for r in srows if r["block_1mb"]}),
        "top_weight": max(weights.values()) if weights else "",
        "n_rows_in_measured_peak": len(in_peak),
        "n_rows_with_da_evidence": len(da),
        "n_rows_with_measured_allelic_imbalance": len(ase_rows),
        "n_rows_with_eqtl_direction": len(eqtl_rows),
        "n_rows_with_model_api_rescue": len(rescued),
        "consequence_coverage": cat.get("consequence_coverage", ""), "coding_mass": cat.get("coding_mass", ""),
        "class_coding": coding, "class_local_regulatory": local_reg,
        "class_rna_processing": rna_proc, "class_candidate_long_range": long_range,
        "n_classes_distinguishable_from_measurement": len(distinguishable),
        "classes_distinguishable_from_measurement": ";".join(distinguishable) if distinguishable else "none",
        "locus_state": ("unresolved" if not distinguishable else ";".join(distinguishable)),
        "candidate_mechanism_atlas_predicted": cat.get("candidate_mechanism", ""),
        "next_experiment_rule_id": cat.get("next_experiment_rule_id", ""),
        "atlas_extension_state": "candidate_not_adopted",
    }


def build_summary(rows, sig_rows, gate, signals, mp, check_r1, track_panel, gene_tracks) -> dict:
    def share(k, tot):
        return (k / tot) if tot else None
    n = len(rows)
    by_u = defaultdict(list)
    for r in rows:
        by_u[r["universe"]].append(r)
    return {
        "what_this_is": ("one row per (signal, variant) for every variant above 0.01 posterior plus each "
                         "signal's top-weight variant, over all eligible Atlas signals; candidate layer only"),
        "seed": SEED, "weight_floor": WEIGHT_FLOOR,
        "n_signals": len(signals), "n_rows": n,
        "rows_by_universe": {u: len(v) for u, v in sorted(by_u.items())},
        "signals_by_universe": {u: len({r["signal_uid"] for r in v}) for u, v in sorted(by_u.items())},
        "rows_by_trait_class": {t: sum(r["trait_class"] == t for r in rows)
                                for t in sorted({r["trait_class"] for r in rows})},
        "rows_by_posterior_definition": {p: sum(r["posterior_definition"] == p for r in rows)
                                         for p in sorted({r["posterior_definition"] for r in rows})},
        "rows_by_variant_class": {c: sum(r["variant_class"] == c for r in rows)
                                  for c in sorted({r["variant_class"] for r in rows})},
        "rows_by_coverage_state": {s: sum(r["signal_atlas_coverage_state"] == s for r in rows)
                                   for s in sorted({r["signal_atlas_coverage_state"] for r in rows})},
        "signals_by_coverage_state": {s: len({r["signal_uid"] for r in rows if r["signal_atlas_coverage_state"] == s})
                                      for s in sorted({r["signal_atlas_coverage_state"] for r in rows})},
        "n_rows_served": sum(r["atlas_served"] == "True" for r in rows),
        "n_rows_readable": sum(r["row_readable"] == "True" for r in rows),
        "readable_share": share(sum(r["row_readable"] == "True" for r in rows), n),
        "readable_share_by_universe": {u: share(sum(r["row_readable"] == "True" for r in v), len(v))
                                       for u, v in sorted(by_u.items())},
        "n_rows_not_readable_because": {w: sum(r["row_not_readable_because"] == w for r in rows)
                                        for w in sorted({r["row_not_readable_because"] for r in rows})},
        "n_signals_top_variant_unserved": len({r["signal_uid"] for r in rows
                                               if int(r["weight_rank"]) == 1 and r["atlas_served"] == "False"}),
        "n_rows_in_measured_peak": sum(1 for r in rows if r["measured_peak"]),
        "n_rows_with_peak_table_coverage": sum(1 for r in rows if r["measured_lineages"] or r["measured_peak"]),
        "rows_by_measured_lineage": {l: sum(l in r["measured_lineages"].split(";") for r in rows)
                                     for l in ("hepatocyte", "stellate")},
        "n_rows_with_da_evidence_supported_or_source_dependent": sum(
            1 for r in rows if set(r["measured_da_evidence_state"].split(";")) & MEASURED_DA_STATES),
        "n_rows_with_measured_allelic_imbalance": sum(
            1 for r in rows if r["ase_gse281367_mean_log2_alt_over_ref"] or r["ase_gse244832_mean_log2_alt_over_ref"]),
        "n_rows_with_eqtl_direction": sum(1 for r in rows if r["eqtl_beta_allele1"]),
        "n_rows_with_model_api_rescue": sum(1 for r in rows if r["model_api_splice_log2"]),
        "n_rows_with_an_anchor": sum(1 for r in rows if r["anchor_rsid"]),
        "n_distinct_1mb_blocks": len({r["block_1mb"] for r in rows if r["block_1mb"]}),
        "n_distinct_analysis_blocks": len({r["analysis_block"] for r in rows if r["analysis_block"]}),
        "signal_class_counts": {c: {v: sum(s[c] == v for s in sig_rows)
                                    for v in sorted({s[c] for s in sig_rows})}
                                for c in ("class_coding", "class_local_regulatory",
                                          "class_rna_processing", "class_candidate_long_range")},
        "signals_by_n_classes_distinguishable": {k: sum(s["n_classes_distinguishable_from_measurement"] == k
                                                        for s in sig_rows) for k in (0, 1, 2, 3, 4)},
        "signals_unresolved": sum(s["locus_state"] == "unresolved" for s in sig_rows),
        "gate_counts": {g: sum(x["signal_gate"] == g for x in gate)
                        for g in sorted({x["signal_gate"] for x in gate})},
        "adult_liver_track_panel": {**track_panel, "rna": gene_tracks["rna"], "splice_site_usage": gene_tracks["splice"]},
        "why_not_liver_summaries": ("liver_summaries.tsv.gz pools the CL:0000182 hepatocyte tracks, which are "
                                    "EMBRYONIC in this Atlas build; this layer reads adult q|primary_liver| only"),
        "R1_reproduces_p6d_dossier": check_r1,
        "B2_mpra_rekey": mp,
        "claim_boundary": ("candidate research layer; no evidence class, gene membership, threshold, hero gene, "
                           "figure or printed Resource number changes"),
    }


INPUTS = [
    ATLAS / "p0-dbsnp-orientation-20260914T121856Z/tables/dbsnp_indel_orientation.tsv",
    TRACK0 / "eligible_signals.tsv", TRACK0 / "gene_catalog_atlas_sidecar.tsv",
    TRACK0 / "signal_variant_weights.tsv.gz", TRACK0 / "variant_peak_overlaps.tsv.gz",
    TRACK0 / "direction_variant_level.tsv.gz", TRACK0 / "atlas_availability.tsv.gz",
    TRACK0 / "posterior_summaries.tsv.gz",
    RECORDS / "ATAC.parquet", RECORDS / "DNASE.parquet", RECORDS / "CHIP_HISTONE.parquet",
    RECORDS / "RNA_SEQ.parquet", RECORDS / "SPLICE_SITE_USAGE.parquet", RECORDS / "AVI_SCORE.parquet",
    COVERAGE_TSV, COVERAGE_FLAGS_TSV, ASE_281367, ASE_244832, DOSSIER_P6D, RESCUE_TSV, REPAIR_TSV,
    DBSNP_TSV, ANCHORS_TSV, MPRA_P3A,
    MPRA_SRC / "HepG2.control.official_comparison.tsv.gz", MPRA_SRC / "HepG2.PAOA.official_comparison.tsv.gz",
    MPRA_SRC / "LX2.control.official_comparison.tsv.gz", MPRA_SRC / "LX2.TGFb.official_comparison.tsv.gz",
    MPRA_SRC / "HepG2.source_reproduction_qc.json", MPRA_SRC / "LX2.source_reproduction_qc.json",
    pathlib.Path(__file__).resolve(),
    PROJECT / "scripts/analysis/alphagenome_atlas/65_clinical_dossier.py",
    PROJECT / "scripts/analysis/alphagenome_atlas/lib_atlas.py",
    PROJECT / "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
]
SHA_SIZE_LIMIT = 400 * 1024 * 1024


def write_provenance(outdir: pathlib.Path) -> None:
    rows = []
    for p in INPUTS:
        st = p.stat()
        rows.append({"path": str(p), "exists": "True", "bytes": st.st_size,
                     "sha256": (sha256_file(p) if st.st_size <= SHA_SIZE_LIMIT else "not_hashed_over_400MB"),
                     "superseded_marker": str((p.parent.parent / "SUPERSEDED.txt").exists())})
    write_tsv_once(outdir / "MANIFEST.tsv", rows, list(rows[0].keys()))
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True)
    (outdir / "pip_freeze.txt").write_text(
        f"# python {sys.version}\n# executable {sys.executable}\n{freeze.stdout}")


if __name__ == "__main__":
    s = main()
    write_provenance(pathlib.Path(sys.argv[1]).resolve())
    print(json.dumps({k: s[k] for k in ("n_signals", "n_rows", "n_rows_readable", "readable_share")}, indent=1))
