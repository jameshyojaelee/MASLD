#!/usr/bin/env python3
"""Step 50 (P5A): are two variants at one genetic signal additive in silico?

The Atlas serves single variants only, so haplotypes are scored with the AlphaGenome model API
(`dna_client.predict_sequence`), which accepts an arbitrary sequence. For each eligible pair the four sequences
REF, V1, V2 and V1+V2 are built from the GRCh38 FASTA over a fixed 1-Mb window centred on the pair midpoint and
predicted with the same requested outputs and liver ontology terms.

Prespecified before any prediction (this docstring plus 50_haplotype_prespec.json):
  eligible pair = two mapped, queried SNVs of one signal, each carrying posterior weight >= MIN_WEIGHT, at most
  MAX_SEP bp apart, in one 1-Mb analysis block; per signal the top PAIRS_PER_SIGNAL pairs by weight product.
  summary per output = log2((sum of predicted signal in the readout window + EPS) / (same for REF + EPS));
  readout windows: RNA_SEQ = the target gene's span from GENCODE v49 (or the central 20 kb when the gene is not
  in the window), ATAC / DNASE / CHIP_HISTONE(H3K27ac) = +/- 1 kb around each variant, taken as the union.
  additivity residual = effect(V1+V2) - effect(V1) - effect(V2). A residual is called non-additive only against
  the matched random-pair null (same chromosome, same separation decile, both variants common in 1000G EUR).
Nothing is renormalised; the null regenerates the whole statistic.

Outputs (tables/): haplotype_pairs.tsv, haplotype_effects.tsv, haplotype_null_effects.tsv, haplotype_summary.json
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import sys
from collections import defaultdict

import numpy as np
import pysam
from alphagenome.data import genome
from alphagenome.models import dna_client

import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
GTF = la.PROJECT / "Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz"
TRACK0 = la.track0_root() / "tables"
LIVER_TERMS = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
OUTPUTS = ["RNA_SEQ", "ATAC", "DNASE", "CHIP_HISTONE"]
WINDOW = dna_client.SEQUENCE_LENGTH_1MB
MIN_WEIGHT = 0.05
MAX_SEP = 100_000
PAIRS_PER_SIGNAL = 3
FLANK = 1_000
EPS = 1e-6
SEED = 20260909
BH_Q = 0.10          # family = every eligible pair scored in this run (prespec amendment 2026-09-10)
N_NULL = 30          # bound from argv in main(); a module-level argv read breaks test import
N_NULL_DEFAULT = 30


def n_null_from_argv(argv: list) -> int:
    """Null-pair count from the command line. Importing the module must not depend on argv (the tests do)."""
    if len(argv) > 1:
        return int(argv[1])
    return N_NULL_DEFAULT


def load_pairs() -> list[dict]:
    signals = {s["signal_uid"]: s for s in la.read_tsv(TRACK0 / "eligible_signals.tsv")}
    w = defaultdict(dict)
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["mapping_status"] == "mapped" and r["in_query_set"] == "True" and float(r["weight"]) >= MIN_WEIGHT:
                w[r["signal_uid"]][r["variant_uid"]] = float(r["weight"])
    pairs = []
    for sig, d in w.items():
        s = signals[sig]
        items = sorted(d.items(), key=lambda kv: -kv[1])
        cand = []
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                u1, w1 = items[i]; u2, w2 = items[j]
                c1, p1, r1, a1 = u1.split(":"); c2, p2, r2, a2 = u2.split(":")
                if c1 != c2 or len(r1) != 1 or len(a1) != 1 or len(r2) != 1 or len(a2) != 1:
                    continue
                sep = abs(int(p1) - int(p2))
                if sep == 0 or sep > MAX_SEP:
                    continue
                cand.append({"signal_uid": sig, "universe": s["universe"], "gene": s["gene"], "ensembl": s["ensembl"], "analysis_block": s["analysis_block"],
                             "v1": u1, "v2": u2, "w1": w1, "w2": w2, "weight_product": w1 * w2, "separation_bp": sep, "chrom": c1})
        pairs.extend(sorted(cand, key=lambda r: -r["weight_product"])[:PAIRS_PER_SIGNAL])
    return pairs


def gene_spans() -> dict[str, tuple[str, int, int]]:
    out = {}
    with gzip.open(GTF, "rt") as h:
        for line in h:
            if line.startswith("#"):
                continue
            p = line.split("\t", 9)
            if p[2] != "gene":
                continue
            m = re.search(r'gene_id "([^"]+)"', p[8])
            if m:
                out[m.group(1).split(".")[0]] = (p[0], int(p[3]), int(p[4]))
    return out


def apply_variants(seq: str, start0: int, variants: list[tuple[int, str, str]]) -> str:
    """1-based positions; substitutions only (SNVs), so length is preserved."""
    s = list(seq)
    for pos, ref, alt in variants:
        i = pos - 1 - start0
        if s[i].upper() != ref.upper():
            raise la.ContractError(f"reference mismatch at {pos}: sequence has {s[i]}, expected {ref}")
        s[i] = alt
    return "".join(s)


def window_mask(length: int, start0: int, regions: list[tuple[int, int]]) -> np.ndarray:
    m = np.zeros(length, bool)
    for a, b in regions:
        lo, hi = max(0, a - start0), min(length, b - start0)
        if hi > lo:
            m[lo:hi] = True
    return m


def effect_threshold(null_v1: np.ndarray, null_v2: np.ndarray) -> float:
    """95th percentile of |single-variant effect| over the matched random-pair null, both arms pooled."""
    a = np.abs(np.concatenate([np.asarray(null_v1, float), np.asarray(null_v2, float)]))
    a = a[~np.isnan(a)]
    return float(np.quantile(a, 0.95)) if a.size else float("nan")


def effect_stratum(v1: np.ndarray, v2: np.ndarray, threshold: float) -> np.ndarray:
    """Pairs where BOTH variants carry a predicted effect above the null threshold."""
    a, b = np.abs(np.asarray(v1, float)), np.abs(np.asarray(v2, float))
    with np.errstate(invalid="ignore"):
        return (~np.isnan(a)) & (~np.isnan(b)) & (a > threshold) & (b > threshold)


def empirical_two_sided_p(residual: float, null: np.ndarray) -> float:
    """Two-sided empirical p for one residual against the matched random-pair null.

    (count of |null| >= |residual|, plus one) / (draws plus one). The plus-one keeps p away from zero: with a
    finite null the smallest honest statement is 1/(n+1), not "impossible".
    """
    nul = np.asarray(null, float)
    nul = nul[~np.isnan(nul)]
    if nul.size == 0 or residual != residual:
        return float("nan")
    return float((np.sum(np.abs(nul) >= abs(float(residual))) + 1) / (nul.size + 1))


def bh_reject(pvals: np.ndarray, q: float) -> np.ndarray:
    """Benjamini-Hochberg step-up over the whole family; returns the rejection mask.

    The family is every eligible pair scored in this run. Comparing each residual with the null's 95th
    percentile instead would flag about 5% of 1,587 pairs by construction.
    """
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    out = np.zeros(p.shape, bool)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx], kind="mergesort")]
    m = order.size
    thresh = q * (np.arange(1, m + 1) / m)
    passing = np.flatnonzero(p[order] <= thresh)
    if passing.size:
        out[order[: passing[-1] + 1]] = True
    return out


def fold_mask(mask: np.ndarray, n_rows: int) -> np.ndarray:
    """Fold a base-resolution readout mask onto a track's own row count.

    The model returns ATAC, DNase and RNA_SEQ at 1 bp but CHIP_HISTONE at 128 bp (8,192 rows over a 1-Mb
    window), so a base-resolution mask cannot index every track. A bin is kept when ANY base inside it is
    selected, which is the behaviour a readout window wants: a variant's +/-1 kb flank should keep the bins
    it touches rather than being lost to rounding.
    """
    mask = np.asarray(mask, bool)
    if n_rows == mask.shape[0]:
        return mask
    if n_rows <= 0 or mask.shape[0] % n_rows:
        raise la.ContractError(f"track rows {n_rows} do not divide the base mask of {mask.shape[0]}")
    return mask.reshape(n_rows, mask.shape[0] // n_rows).any(axis=1)


def summarise(output, mask_rna: np.ndarray, mask_local: np.ndarray) -> dict[str, float]:
    vals = {}
    for name, td, mask in (("rna", output.rna_seq, mask_rna), ("atac", output.atac, mask_local), ("dnase", output.dnase, mask_local), ("h3k27ac", output.chip_histone, mask_local)):
        if td is None:
            vals[name] = math.nan
            continue
        v = np.asarray(td.values, np.float64)
        keep = np.ones(v.shape[1], bool)
        if name == "h3k27ac" and td.metadata is not None and "name" in td.metadata:
            keep = td.metadata["name"].astype(str).str.contains("H3K27ac").to_numpy()
            if not keep.any():
                vals[name] = math.nan
                continue
        vals[name] = float(v[fold_mask(mask, v.shape[0])][:, keep].sum())
    return vals


MIN_MATCHED_NULL = 20


def stratum_is_testable(n_observed: int, n_null_matched: int, min_null: int = MIN_MATCHED_NULL) -> dict:
    """Can a stratum be tested against the null, given how many null pairs meet the same selection?

    A stratum defined by large single-variant effects must be compared with null pairs that also carry large
    single-variant effects. With only 1-3 of 300 null pairs clearing that bar, the comparison degenerates
    into active-versus-inactive and rejects most of the stratum by construction.
    """
    if n_null_matched < min_null:
        return {"testable": False, "n_null_matched": int(n_null_matched), "min_null": int(min_null),
                "reason": (f"only {n_null_matched} of the null pairs meet this stratum's selection "
                           f"(minimum {min_null}); there is no matched null, so no verdict is emitted. "
                           "Testing it would need enough null draws to yield matched active pairs.")}
    return {"testable": True, "n_null_matched": int(n_null_matched), "min_null": int(min_null), "reason": "matched null available"}


def score_channel(ch: str, rows: list, null_rows: list) -> dict:
    """Additivity verdict for one channel, on both prespecified strata.

    (a) every eligible pair; (b) pairs where BOTH variants carry a single-variant |effect| above the matched
    null's 95th percentile. Additivity is undefined when neither variant does anything, and stratum (b)
    conditions on components of the statistic, so it is reported beside (a) and never instead of it. BH runs
    within each stratum with its own family size.
    """
    obs = np.array([r[f"{ch}_residual"] for r in rows], float)
    nul = np.array([r[f"{ch}_residual"] for r in null_rows], float) if null_rows else np.array([])
    v1 = np.array([r.get(f"{ch}_v1_log2", math.nan) for r in rows], float)
    v2 = np.array([r.get(f"{ch}_v2_log2", math.nan) for r in rows], float)
    nv1 = np.array([r.get(f"{ch}_v1_log2", math.nan) for r in null_rows], float) if null_rows else np.array([])
    nv2 = np.array([r.get(f"{ch}_v2_log2", math.nan) for r in null_rows], float) if null_rows else np.array([])
    thr = effect_threshold(nv1, nv2) if nv1.size else math.nan
    keep = effect_stratum(v1, v2, thr) if thr == thr else np.zeros(obs.shape, bool)
    q = float(np.nanquantile(np.abs(nul), 0.95)) if nul.size else math.nan
    pvals = np.array([empirical_two_sided_p(v, nul) for v in obs], float)

    out = {"median_abs_residual": float(np.nanmedian(np.abs(obs))) if obs.size else math.nan,
           "null_median_abs_residual": float(np.nanmedian(np.abs(nul))) if nul.size else math.nan,
           "null_q95_abs_residual": q,
           "n_pairs_exceeding_null_q95_UNCORRECTED": int(np.nansum(np.abs(obs) > q)) if q == q else None,
           "median_abs_single_variant_effect": float(np.nanmedian(np.abs(np.concatenate([v1, v2])))) if obs.size else math.nan,
           "single_variant_effect_threshold_null_q95": thr,
           "bh_q": BH_Q}
    null_keep = effect_stratum(nv1, nv2, thr) if (thr == thr and nv1.size) else np.zeros(nul.shape, bool)
    out["null_pairs_meeting_active_stratum"] = int(null_keep.sum())
    for label, mask in (("all_pairs", np.ones(obs.shape, bool)), ("both_variants_active", keep)):
        gate = stratum_is_testable(int(mask.sum()), int(null_keep.sum()) if label == "both_variants_active" else int(nul.size))
        if not gate["testable"]:
            out[label] = {"family_size_pairs": int(mask.sum()), "testable": False, **gate}
            continue
        pm = np.where(mask, pvals, np.nan)
        rej = bh_reject(pm, BH_Q)
        for i, r in enumerate(rows):
            r[f"{ch}_p_empirical"] = pvals[i]
            r[f"{ch}_both_active"] = bool(keep[i])
            if label == "both_variants_active":
                r[f"{ch}_bh_reject_active"] = bool(rej[i])
            else:
                r[f"{ch}_bh_reject_all"] = bool(rej[i])
        sigs = sorted({rows[i]["signal_uid"] for i in np.flatnonzero(rej)})
        out[label] = {"family_size_pairs": int(np.sum(~np.isnan(pm))), "testable": True, **gate,
                      "n_pairs_bh": int(rej.sum()), "n_signals_bh": len(sigs), "signals_bh": sigs,
                      "min_p_empirical": float(np.nanmin(pm)) if np.any(~np.isnan(pm)) else math.nan}
    if int(keep.sum()) < 20:
        out["verdict"] = (f"untestable: only {int(keep.sum())} pairs have a predicted effect above the null threshold "
                          f"in this channel, so additivity is undefined for essentially every Resource variant pair here")
    return out


def rescore_from_tables() -> None:
    """Recompute the summary from deposited effects tables, without re-querying the model."""
    rows = la.read_tsv(TABLES / "haplotype_effects.tsv")
    null_rows = la.read_tsv(TABLES / "haplotype_null_effects.tsv")
    for r in rows + null_rows:
        for k, v in list(r.items()):
            if k.endswith(("_log2", "_residual")):
                r[k] = float(v) if v not in ("", "nan", None) else math.nan
    summary = {"n_pairs": len(rows), "n_null": len(null_rows), "rescored_from_tables": True}
    for ch in ("rna", "atac", "dnase", "h3k27ac"):
        summary[ch] = score_channel(ch, rows, null_rows)
    summary["claim_boundary"] = "in silico only; predicted non-additivity is not an epistasis claim"
    json.dump(summary, (TABLES / "haplotype_summary.json").open("w"), indent=1, default=float)
    la.log(f"P5A rescored: {len(rows)} pairs, {len(null_rows)} null")


def main() -> None:
    global N_NULL
    N_NULL = n_null_from_argv(sys.argv)
    prespec = {"min_weight": MIN_WEIGHT, "max_separation_bp": MAX_SEP, "pairs_per_signal": PAIRS_PER_SIGNAL, "window_bp": WINDOW, "flank_bp": FLANK,
               "outputs": OUTPUTS, "ontology_terms": LIVER_TERMS, "seed": SEED, "n_null_pairs": N_NULL,
               "statistic": "effect = log2((sum predicted signal in readout window + 1e-6)/(REF sum + 1e-6)); residual = joint - v1 - v2",
               "null": "matched random SNV pairs: same chromosome, same separation decile, both alleles taken from the FASTA reference and a fixed alternate rule",
               "written_before": "any prediction was made"}
    path = TABLES / "haplotype_prespec.json"
    if not path.exists():
        json.dump(prespec, path.open("w"), indent=1)
    pairs = load_pairs()
    if not pairs:
        la.log("no eligible variant pairs"); return
    la.write_tsv_once(TABLES / "haplotype_pairs.tsv", pairs, list(pairs[0].keys()))
    la.log(f"eligible pairs: {len(pairs)} over {len({p['signal_uid'] for p in pairs})} signals")
    fasta = pysam.FastaFile(FASTA)
    spans = gene_spans()
    key, _ = la.load_api_key()
    model = dna_client.create(key)
    outs = [getattr(dna_client.OutputType, o) for o in OUTPUTS]
    rng = np.random.default_rng(SEED)

    def score_pair(chrom: str, p1: int, r1: str, a1: str, p2: int, r2: str, a2: str, ensembl: str, tag: str) -> dict:
        mid = (p1 + p2) // 2
        start0 = max(0, mid - WINDOW // 2)
        end = start0 + WINDOW
        seq = fasta.fetch(chrom, start0, end).upper()
        if len(seq) != WINDOW:
            return {}
        interval = genome.Interval(chromosome=chrom, start=start0, end=end)
        span = spans.get(ensembl)
        rna_regions = [(span[1], span[2])] if span and span[0] == chrom and span[2] > start0 and span[1] < end else [(mid - 10_000, mid + 10_000)]
        mask_rna = window_mask(WINDOW, start0, rna_regions)
        mask_local = window_mask(WINDOW, start0, [(p1 - FLANK, p1 + FLANK), (p2 - FLANK, p2 + FLANK)])
        arms = {"ref": [], "v1": [(p1, r1, a1)], "v2": [(p2, r2, a2)], "joint": [(p1, r1, a1), (p2, r2, a2)]}
        summ = {}
        for arm, vs in arms.items():
            s = apply_variants(seq, start0, vs)
            o = aq.call_with_quota_retry(lambda: model.predict_sequence(sequence=s, requested_outputs=outs, ontology_terms=LIVER_TERMS, interval=interval), label=f"{tag}:{arm}")
            summ[arm] = summarise(o, mask_rna, mask_local)
        row = {"tag": tag, "chrom": chrom, "pos1": p1, "pos2": p2, "separation_bp": abs(p2 - p1), "ensembl": ensembl,
               "rna_readout": "gene_span" if span and span[0] == chrom else "central_20kb"}
        for ch in ("rna", "atac", "dnase", "h3k27ac"):
            ref = summ["ref"][ch]
            for arm in ("v1", "v2", "joint"):
                row[f"{ch}_{arm}_log2"] = math.log2((summ[arm][ch] + EPS) / (ref + EPS)) if ref == ref and summ[arm][ch] == summ[arm][ch] else math.nan
            row[f"{ch}_residual"] = row[f"{ch}_joint_log2"] - row[f"{ch}_v1_log2"] - row[f"{ch}_v2_log2"]
        return row

    rows = []
    for p in pairs:
        c, p1, r1, a1 = p["v1"].split(":"); _, p2, r2, a2 = p["v2"].split(":")
        r = score_pair(c, int(p1), r1, a1, int(p2), r2, a2, p["ensembl"], f"{p['signal_uid']}|{p['v1']}|{p['v2']}")
        if r:
            rows.append({**{k: p[k] for k in ("signal_uid", "universe", "gene", "analysis_block", "v1", "v2", "w1", "w2", "weight_product")}, **r})
            la.log(f"pair {p['signal_uid']} {p['v1']} {p['v2']}: rna residual {r['rna_residual']:+.4f}, atac {r['atac_residual']:+.4f}")
    la.write_tsv_once(TABLES / "haplotype_effects.tsv", rows, list(rows[0].keys()))

    # matched random-pair null: same chromosome and separation, reference alleles from the FASTA, alt = fixed transversion rule
    null_rows = []
    seps = [r["separation_bp"] for r in rows]
    for i in range(N_NULL):
        base = rows[int(rng.integers(0, len(rows)))]
        chrom = base["chrom"]; sep = seps[int(rng.integers(0, len(seps)))]
        length = fasta.get_reference_length(chrom)
        for _ in range(50):
            p1 = int(rng.integers(WINDOW // 2 + 1, length - WINDOW // 2 - sep - 1))
            p2 = p1 + sep
            r1 = fasta.fetch(chrom, p1 - 1, p1).upper(); r2 = fasta.fetch(chrom, p2 - 1, p2).upper()
            if r1 in "ACGT" and r2 in "ACGT":
                break
        else:
            continue
        alt = {"A": "T", "T": "A", "C": "G", "G": "C"}
        r = score_pair(chrom, p1, r1, alt[r1], p2, r2, alt[r2], base["ensembl"], f"null{i}")
        if r:
            null_rows.append(r)
    if null_rows:
        la.write_tsv_once(TABLES / "haplotype_null_effects.tsv", null_rows, list(null_rows[0].keys()))
    summary = {"n_pairs": len(rows), "n_null": len(null_rows), "prespec": prespec}
    for ch in ("rna", "atac", "dnase", "h3k27ac"):
        summary[ch] = score_channel(ch, rows, null_rows)
    summary["claim_boundary"] = "in silico only; predicted non-additivity is not an epistasis claim and does not change the recorded `indeterminate` evidence state for variant combinations"
    json.dump(summary, (TABLES / "haplotype_summary.json").open("w"), indent=1, default=float)
    la.log(f"P5A: {summary['n_pairs']} pairs, {summary['n_null']} null pairs; rna {summary['rna']}")


if __name__ == "__main__":
    main()
