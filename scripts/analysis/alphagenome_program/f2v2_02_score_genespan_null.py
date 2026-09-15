#!/usr/bin/env python3
"""f2-haplotype-v2 step 2: score the 100 gene-span-matched chr6 null draws. 400 predict_sequence calls.

Same recipe (f2_recipe.py), same four arms (ref, v1, v2, joint), same +/-1 kb local readout, same 1,048,576 bp
window, same liver ontology terms as v1. The ONLY difference from v1's gnmtnull set is that each draw's
Ensembl id is the gene that hosts it, so score_pair's existing gene-span branch fires and the RNA readout is
a gene span instead of the central 20 kb.

Rate ceiling 2.15 s between call STARTS = 27.9 calls/min, at or below the 28/min the task allows; the API key
is shared with job 21611396 in another session.

Checkpointing: every scored draw is appended to raw/scored_genespan_null.jsonl keyed by tag. A restart skips
tags already present, so no call is paid for twice.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import os
import pathlib
import sys
import time

import pysam
from alphagenome.data import genome
from alphagenome.models import dna_client

sys.path.insert(0, "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas")
import atlas_query as aq          # noqa: E402  (read-only reuse of the quota-retry wrapper)
import lib_atlas as la            # noqa: E402  (read-only reuse of load_api_key)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import f2_recipe as R             # noqa: E402

OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
DRAWS = OUT / "tables" / "gnmt_genespan_null_pairs.tsv"
CHECKPOINT = OUT / "raw" / "scored_genespan_null.jsonl"
MIN_SECONDS_PER_CALL = 2.15


def log(msg: str) -> None:
    print(f"[{dt.datetime.now(dt.timezone.utc).strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


class Scorer:
    """Byte-for-byte the same scoring path as v1's f2_02_score.Scorer, with a 27.9/min ceiling."""

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
                   ensembl: str, tag: str) -> dict:
        mid = (p1 + p2) // 2
        start0 = max(0, mid - R.WINDOW // 2)
        end = start0 + R.WINDOW
        seq = self.fasta.fetch(chrom, start0, end).upper()
        if len(seq) != R.WINDOW:
            return {}
        interval = genome.Interval(chromosome=chrom, start=start0, end=end)
        span = self.spans.get(ensembl)
        use_span = bool(span and span[0] == chrom and span[2] > start0 and span[1] < end)
        rna_regions = [(span[1], span[2])] if use_span else [(mid - 10_000, mid + 10_000)]
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
               "rna_readout": "gene_span" if use_span else "central_20kb",
               "rna_readout_bp": int(mask_rna.sum()), "local_readout_bp": int(mask_local.sum())}
        import math
        for ch in R.CHANNELS:
            ref = summ["ref"][ch]
            row[f"{ch}_ref_sum"] = ref
            for arm in ("v1", "v2", "joint"):
                row[f"{ch}_{arm}_sum"] = summ[arm][ch]
                row[f"{ch}_{arm}_log2"] = (math.log2((summ[arm][ch] + R.EPS) / (ref + R.EPS))
                                           if ref == ref and summ[arm][ch] == summ[arm][ch] else math.nan)
            row[f"{ch}_residual"] = row[f"{ch}_joint_log2"] - row[f"{ch}_v1_log2"] - row[f"{ch}_v2_log2"]
        return row


def main() -> None:
    (OUT / "raw").mkdir(parents=True, exist_ok=True)
    with open(DRAWS) as handle:
        draws = list(csv.DictReader(handle, delimiter="\t"))
    log(f"{len(draws)} gene-span-matched draws to score")

    done = set()
    if CHECKPOINT.exists():
        with open(CHECKPOINT) as handle:
            for line in handle:
                try:
                    done.add(json.loads(line)["tag"])
                except Exception:  # noqa: BLE001 - a truncated final line on a killed job
                    pass
    log(f"checkpoint: {len(done)} draws already scored")
    todo = [d for d in draws if d["tag"] not in done]
    log(f"{len(todo)} to score, {4 * len(todo)} calls at <= 27.9/min "
        f"(about {4 * len(todo) / 27.9:.1f} min of rate-limited time)")

    sc = Scorer()
    t0 = time.time()
    with open(CHECKPOINT, "a") as ck:
        for n, d in enumerate(todo, 1):
            row = sc.score_pair("chr6", int(d["pos1"]), d["ref1"], d["alt1"], int(d["pos2"]), d["ref2"], d["alt2"],
                                d["ensembl"], d["tag"])
            if not row:
                log(f"skipped (window off chromosome): {d['tag']}")
                continue
            keep = {k: d[k] for k in ("gene_name", "gene_type", "span_start", "span_end", "span_width_bp",
                                      "decile_bin", "bin_lo", "bin_hi", "frac_pos1_in_span", "frac_pos2_in_span",
                                      "n_qualifying_pairs_in_gene", "base_tag")}
            ck.write(json.dumps({**keep, **row}) + "\n")
            ck.flush()
            if n % 10 == 0 or n <= 5:
                rate = sc.calls / max(1e-9, (time.time() - t0) / 60.0)
                log(f"{n}/{len(todo)} {d['tag']} {d['gene_name']} readout {row['rna_readout']} "
                    f"({row['rna_readout_bp']} bp): rna {row['rna_residual']:+.5f} atac {row['atac_residual']:+.5f} "
                    f"({sc.calls} calls, {rate:.1f}/min)")
    log(f"DONE: {sc.calls} calls in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
