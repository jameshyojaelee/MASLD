#!/usr/bin/env python3
"""F2/F3 step 2: score the GNMT haplotype and the 200-pair family against the corrected matched null.

Four work lists, all through one scorer and one readout definition:
  gnmt      6 unordered pairs of the four GNMT-region variants
  gnmtnull  100 chr6 draws matched to the GNMT primary pair's separation decile
  obs       the 200 highest weight-product pairs of the P5 pair list, rescored here
  null      1,000 draws matched on the observed pair's chromosome and separation decile

The recipe is the byte-equivalent copy in f2_recipe.py. The quota-retry wrapper is reused from
scripts/analysis/alphagenome_atlas/atlas_query.py (read-only import; nothing in that tree is modified).

Checkpointing: every scored pair is appended to raw/scored_rows.jsonl as one JSON object keyed by `tag`.
A restart reads that file and skips tags already present, so no call is paid for twice.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import math
import os
import pathlib
import sys
import time

import numpy as np
import pysam
from alphagenome.data import genome
from alphagenome.models import dna_client

sys.path.insert(0, "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas")
import atlas_query as aq          # noqa: E402  (read-only reuse of the quota-retry wrapper)
import lib_atlas as la            # noqa: E402  (read-only reuse of load_api_key)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import f2_recipe as R             # noqa: E402

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
P5_PAIRS = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_pairs.tsv"
CHECKPOINT = OUT / "raw" / "scored_rows.jsonl"

SEED = 20260914
N_OBS = 200
N_NULL_DEFAULT = 1000          # bound from argv in main(); see f2_00_amendment_03.json for the extension rule
N_GNMT_NULL = 100
MIN_SECONDS_PER_CALL = 2.1     # 28.6 calls/min ceiling, measured between call STARTS; job 21611396 shares this key
MAX_DRAW_TRIES = 200

GNMT_ENSEMBL = "ENSG00000124713"
# REF is the GRCh38 FASTA base; ALT is the dbSNP b157 GRCh38.p14 frequency-bearing alternate.
GNMT_VARIANTS = [
    {"rsid": "rs9471976", "chrom": "chr6", "pos": 42952811, "ref": "G", "alt": "T", "palindromic": False},
    {"rsid": "rs11752813", "chrom": "chr6", "pos": 42960279, "ref": "C", "alt": "G", "palindromic": True},
    {"rsid": "rs2296805", "chrom": "chr6", "pos": 42961020, "ref": "T", "alt": "G", "palindromic": False},
    {"rsid": "rs2296804", "chrom": "chr6", "pos": 42963523, "ref": "C", "alt": "G", "palindromic": True},
]
GNMT_PRIMARY = ("rs2296805", "rs2296804")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now(dt.timezone.utc).strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- inputs
def load_top_pairs(n: int) -> list[dict]:
    with open(P5_PAIRS) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda r: (-float(r["weight_product"]), r["signal_uid"], r["v1"], r["v2"]))
    return rows[:n]


def load_panel(chrom: str) -> np.ndarray:
    """Common-EUR panel for one chromosome: structured array of hg38 pos, a1, a2."""
    path = OUT / "tables" / "null_panel" / f"{chrom}.tsv.gz"
    pos, a1, a2 = [], [], []
    with gzip.open(path, "rt") as handle:
        handle.readline()
        for line in handle:
            p = line.rstrip("\n").split("\t")
            pos.append(int(p[1])); a1.append(p[2]); a2.append(p[3])
    order = np.argsort(np.asarray(pos))
    return np.asarray(pos)[order], np.asarray(a1)[order], np.asarray(a2)[order]


def decile_bins(seps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = np.quantile(seps, [i / 10 for i in range(11)])
    idx = np.clip(np.searchsorted(edges[1:-1], seps, side="right"), 0, 9)
    return edges, idx


# --------------------------------------------------------------------------- scorer
class Scorer:
    def __init__(self) -> None:
        self.fasta = pysam.FastaFile(R.FASTA)
        self.spans = R.gene_spans()
        key, src = la.load_api_key()
        log(f"API key from {src}")
        self.model = dna_client.create(key)
        self.outs = [getattr(dna_client.OutputType, o) for o in R.OUTPUTS]
        self.calls = 0
        self.t_last = 0.0

    def _predict(self, seq: str, interval, label: str):
        # the ceiling is on the request RATE, so the clock starts when a call starts, not when it returns;
        # timing from the previous call's return adds its latency on top and halves the achieved rate.
        wait = MIN_SECONDS_PER_CALL - (time.time() - self.t_last)
        if wait > 0:
            time.sleep(wait)
        self.t_last = time.time()
        o = aq.call_with_quota_retry(
            lambda: self.model.predict_sequence(sequence=seq, requested_outputs=self.outs,
                                                ontology_terms=R.LIVER_TERMS, interval=interval),
            label=label)
        self.calls += 1
        return o

    def score_pair(self, chrom: str, p1: int, r1: str, a1: str, p2: int, r2: str, a2: str,
                   ensembl: str, tag: str, extra_ref_call: bool = False) -> dict:
        """Identical to the P5 recipe's score_pair, with the throttle and an optional determinism probe."""
        mid = (p1 + p2) // 2
        start0 = max(0, mid - R.WINDOW // 2)
        end = start0 + R.WINDOW
        seq = self.fasta.fetch(chrom, start0, end).upper()
        if len(seq) != R.WINDOW:
            return {}
        interval = genome.Interval(chromosome=chrom, start=start0, end=end)
        span = self.spans.get(ensembl)
        rna_regions = ([(span[1], span[2])] if span and span[0] == chrom and span[2] > start0 and span[1] < end
                       else [(mid - 10_000, mid + 10_000)])
        mask_rna = R.window_mask(R.WINDOW, start0, rna_regions)
        mask_local = R.window_mask(R.WINDOW, start0, [(p1 - R.FLANK, p1 + R.FLANK), (p2 - R.FLANK, p2 + R.FLANK)])
        arms = {"ref": [], "v1": [(p1, r1, a1)], "v2": [(p2, r2, a2)], "joint": [(p1, r1, a1), (p2, r2, a2)]}
        summ = {}
        for arm, vs in arms.items():
            s = R.apply_variants(seq, start0, vs)
            o = self._predict(s, interval, f"{tag}:{arm}")
            summ[arm] = R.summarise(o, mask_rna, mask_local)
        row = {"tag": tag, "chrom": chrom, "pos1": p1, "pos2": p2, "ref1": r1, "alt1": a1, "ref2": r2, "alt2": a2,
               "separation_bp": abs(p2 - p1), "ensembl": ensembl,
               "rna_readout": "gene_span" if (span and span[0] == chrom and span[2] > start0 and span[1] < end) else "central_20kb"}
        for ch in R.CHANNELS:
            ref = summ["ref"][ch]
            row[f"{ch}_ref_sum"] = ref
            for arm in ("v1", "v2", "joint"):
                row[f"{ch}_{arm}_sum"] = summ[arm][ch]
                row[f"{ch}_{arm}_log2"] = (math.log2((summ[arm][ch] + R.EPS) / (ref + R.EPS))
                                           if ref == ref and summ[arm][ch] == summ[arm][ch] else math.nan)
            row[f"{ch}_residual"] = row[f"{ch}_joint_log2"] - row[f"{ch}_v1_log2"] - row[f"{ch}_v2_log2"]
        if extra_ref_call:
            o2 = self._predict(R.apply_variants(seq, start0, []), interval, f"{tag}:ref_repeat")
            rep = R.summarise(o2, mask_rna, mask_local)
            for ch in R.CHANNELS:
                row[f"{ch}_ref_sum_repeat"] = rep[ch]
        return row


# --------------------------------------------------------------------------- null draws
def draw_null_pairs(rng, base_rows: list[dict], edges: np.ndarray, bin_idx: np.ndarray,
                    panels: dict, fasta: pysam.FastaFile, n: int, prefix: str,
                    force_chrom: str | None = None, force_bin: int | None = None) -> list[dict]:
    """Matched draws: same chromosome and same separation decile bin as a sampled observed pair.

    Both variants are common (EUR MAF >= 0.05) 1000G SNVs lifted to hg38; REF is the GRCh38 base at the
    lifted position and must equal one of the two 1000G alleles, ALT is the other.
    """
    out = []
    tries_total = 0
    while len(out) < n:
        tries_total += 1
        if tries_total > 400 * n:
            raise RuntimeError(f"{prefix}: could not draw {n} matched pairs")
        k = int(rng.integers(0, len(base_rows)))
        base = base_rows[k]
        chrom = force_chrom or base["chrom"]
        b = force_bin if force_bin is not None else int(bin_idx[k])
        lo, hi = float(edges[b]), float(edges[b + 1])
        pos, al1, al2 = panels[chrom]
        i = int(rng.integers(0, pos.size))
        p1 = int(pos[i])
        j_lo = int(np.searchsorted(pos, p1 + lo, side="left"))
        j_hi = int(np.searchsorted(pos, p1 + hi, side="right"))
        k_lo = int(np.searchsorted(pos, p1 - hi, side="left"))
        k_hi = int(np.searchsorted(pos, p1 - lo, side="right"))
        cands = np.concatenate([np.arange(j_lo, j_hi), np.arange(k_lo, k_hi)])
        cands = cands[cands != i]
        if cands.size == 0:
            continue
        j = int(cands[int(rng.integers(0, cands.size))])
        q1, q2 = (i, j) if pos[i] <= pos[j] else (j, i)
        x1, x2 = int(pos[q1]), int(pos[q2])
        if x2 - x1 < lo or x2 - x1 > hi or x2 == x1:
            continue
        length = fasta.get_reference_length(chrom)
        mid = (x1 + x2) // 2
        if mid - R.WINDOW // 2 < 0 or mid - R.WINDOW // 2 + R.WINDOW > length:
            continue
        placed = []
        ok = True
        for q, x in ((q1, x1), (q2, x2)):
            base_ref = fasta.fetch(chrom, x - 1, x).upper()
            a, b2 = str(al1[q]), str(al2[q])
            if base_ref == a:
                placed.append((x, a, b2))
            elif base_ref == b2:
                placed.append((x, b2, a))
            else:
                ok = False
                break
        if not ok:
            continue
        out.append({"tag": f"{prefix}{len(out)}", "chrom": chrom, "ensembl": base["ensembl"],
                    "base_tag": base.get("tag", base.get("signal_uid", "")), "decile_bin": b,
                    "bin_lo": lo, "bin_hi": hi,
                    "pos1": placed[0][0], "ref1": placed[0][1], "alt1": placed[0][2],
                    "pos2": placed[1][0], "ref2": placed[1][1], "alt2": placed[1][2]})
    log(f"{prefix}: {len(out)} matched null pairs drawn in {tries_total} attempts")
    return out


# --------------------------------------------------------------------------- main
def main() -> None:
    n_null = int(sys.argv[1]) if len(sys.argv) > 1 else N_NULL_DEFAULT
    (OUT / "raw").mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    top = load_top_pairs(N_OBS)
    seps = np.array([int(r["separation_bp"]) for r in top], float)
    edges, bins = decile_bins(seps)
    log(f"observed family: {len(top)} pairs, decile edges {np.round(edges, 1).tolist()}")

    fasta = pysam.FastaFile(R.FASTA)
    chroms_needed = sorted({r["chrom"] for r in top} | {"chr6"})
    panels = {c: load_panel(c) for c in chroms_needed}
    log("panel sizes: " + ", ".join(f"{c}={panels[c][0].size}" for c in chroms_needed))

    # ---- work list, fixed before any call
    gnmt_pairs = []
    for i in range(len(GNMT_VARIANTS)):
        for j in range(i + 1, len(GNMT_VARIANTS)):
            v1, v2 = GNMT_VARIANTS[i], GNMT_VARIANTS[j]
            gnmt_pairs.append({"tag": f"gnmt|{v1['rsid']}|{v2['rsid']}", "chrom": "chr6", "ensembl": GNMT_ENSEMBL,
                               "pos1": v1["pos"], "ref1": v1["ref"], "alt1": v1["alt"],
                               "pos2": v2["pos"], "ref2": v2["ref"], "alt2": v2["alt"],
                               "rsid1": v1["rsid"], "rsid2": v2["rsid"],
                               "palindromic1": v1["palindromic"], "palindromic2": v2["palindromic"],
                               "is_primary": (v1["rsid"], v2["rsid"]) == GNMT_PRIMARY})
    gnmt_sep = abs(GNMT_VARIANTS[3]["pos"] - GNMT_VARIANTS[2]["pos"])
    gnmt_bin = int(np.clip(np.searchsorted(edges[1:-1], gnmt_sep, side="right"), 0, 9))
    log(f"GNMT primary separation {gnmt_sep} bp -> decile bin {gnmt_bin} [{edges[gnmt_bin]:.1f}, {edges[gnmt_bin+1]:.1f}]")

    obs_pairs = []
    for r in top:
        c, p1, r1, a1 = r["v1"].split(":")
        _, p2, r2, a2 = r["v2"].split(":")
        obs_pairs.append({"tag": f"obs|{r['signal_uid']}|{r['v1']}|{r['v2']}", "chrom": c, "ensembl": r["ensembl"],
                          "pos1": int(p1), "ref1": r1, "alt1": a1, "pos2": int(p2), "ref2": r2, "alt2": a2,
                          "signal_uid": r["signal_uid"], "universe": r["universe"], "gene": r["gene"],
                          "analysis_block": r["analysis_block"], "v1": r["v1"], "v2": r["v2"],
                          "w1": r["w1"], "w2": r["w2"], "weight_product": r["weight_product"]})

    gnmt_null = draw_null_pairs(rng, [{"chrom": "chr6", "ensembl": GNMT_ENSEMBL, "tag": "gnmt_primary"}],
                                edges, np.array([gnmt_bin]), panels, fasta, N_GNMT_NULL, "gnmtnull",
                                force_chrom="chr6", force_bin=gnmt_bin)
    # the draw stream is sequential and seeded, so the first 1,000 of a 2,500-draw request are the same
    # 1,000 pairs; an extension run therefore only pays for the pairs it adds (amendment 03).
    null_pairs = draw_null_pairs(rng, obs_pairs, edges, bins, panels, fasta, n_null, "null")

    for name, rows in (("gnmt_pairs", gnmt_pairs), ("obs_pairs", obs_pairs),
                       ("gnmt_null_pairs", gnmt_null), (f"null_pairs_n{n_null}", null_pairs)):
        path = OUT / "tables" / f"{name}.tsv"
        if not path.exists():
            cols = sorted({k for r in rows for k in r})
            with open(path, "w") as handle:
                w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
                w.writeheader()
                for r in rows:
                    w.writerow(r)

    done = set()
    if CHECKPOINT.exists():
        with open(CHECKPOINT) as handle:
            for line in handle:
                try:
                    done.add(json.loads(line)["tag"])
                except Exception:  # noqa: BLE001 - a truncated final line on a killed job
                    pass
    log(f"checkpoint: {len(done)} pairs already scored")

    work = gnmt_pairs + gnmt_null + obs_pairs + null_pairs
    todo = [p for p in work if p["tag"] not in done]
    log(f"work list: {len(work)} pairs, {len(todo)} to score, about {4 * len(todo)} calls")

    sc = Scorer()
    t0 = time.time()
    with open(CHECKPOINT, "a") as ck:
        for n, p in enumerate(todo, 1):
            row = sc.score_pair(p["chrom"], p["pos1"], p["ref1"], p["alt1"], p["pos2"], p["ref2"], p["alt2"],
                                p["ensembl"], p["tag"], extra_ref_call=bool(p.get("is_primary")))
            if not row:
                log(f"skipped (window off chromosome): {p['tag']}")
                continue
            ck.write(json.dumps({**{k: v for k, v in p.items() if k != "tag"}, **row}) + "\n")
            ck.flush()
            if n % 25 == 0 or n <= 10:
                rate = sc.calls / max(1e-9, (time.time() - t0) / 60.0)
                log(f"{n}/{len(todo)} {p['tag']}: rna {row['rna_residual']:+.5f} atac {row['atac_residual']:+.5f} "
                    f"({sc.calls} calls, {rate:.1f}/min)")
    log(f"DONE: {sc.calls} calls in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
