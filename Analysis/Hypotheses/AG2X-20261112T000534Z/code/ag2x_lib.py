"""AG2X shared helpers: key loading, hosted-client creation, paced and resumable variant scoring.

Network access happens only in score_chunked() and the client factories. The API key is read by the
existing lib_atlas.load_api_key() (ALPHA_GENOME_API_KEY, then ~/.alphagenome_key) and is never
printed or written; only its length is logged. Requests carry GRCh38 coordinates, alleles and
condition IDs only.

Pacing follows scripts/analysis/alphagenome_atlas/atlas_query.py: fixed chunks written once and
skipped on resume, 65-s sleep on a per-minute quota error, doubling back-off on transient errors,
oversized responses never retried.
"""
from __future__ import annotations

import datetime as dt
import importlib.metadata
import json
import os
import pathlib
import sys
import time

REPO = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
PKG = REPO / "Analysis/Hypotheses/AG2X-20261112T000534Z"
FASTA = REPO / "Analysis/MASLD_Model_Benchmark/executions/sequence-fasta-21063808/GRCh38.p14.sequence_model.fa"
SEED = 20261111

sys.path.insert(0, str(REPO / "scripts/analysis/alphagenome_atlas"))
import lib_atlas  # noqa: E402  (key loader only; no network at import)

QUOTA_SLEEP = 65.0
TRANSIENT_SLEEP = 30.0
TRANSIENT_SLEEP_CAP = 300.0
MAX_TRANSIENT = 14


def log(msg: str) -> None:
    print(f"[{dt.datetime.now(dt.timezone.utc).strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def package_version() -> str:
    return importlib.metadata.version("alphagenome")


def _key() -> str:
    key, src = lib_atlas.load_api_key()
    log(f"key loaded (source type {src.split(':')[0]}, length {len(key)})")
    return key


def create_v2():
    from alphagenome.models.v2 import dna_client
    return dna_client.create(_key(), model_version=dna_client.ModelVersion.ALL_FOLDS)


def create_v1():
    from alphagenome.models.v1 import dna_client
    return dna_client.create(_key())


# ------------------------------------------------------------------ error classes
def _is_rpc(e):
    import grpc
    return isinstance(e, grpc.RpcError)


def is_message_size_error(e) -> bool:
    import grpc
    return _is_rpc(e) and e.code() == grpc.StatusCode.RESOURCE_EXHAUSTED and "message larger than max" in (e.details() or "")


def is_quota_error(e) -> bool:
    import grpc
    return _is_rpc(e) and e.code() == grpc.StatusCode.RESOURCE_EXHAUSTED and not is_message_size_error(e)


def is_transient_error(e) -> bool:
    import grpc
    return _is_rpc(e) and e.code() in (grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.DEADLINE_EXCEEDED, grpc.StatusCode.INTERNAL)


def call_with_retry(fn, label: str, counter: dict, n_calls: int = 0):
    """counter['attempted_calls'] counts first sends; counter['resent_calls'] adds n_calls per retry."""
    transient = 0
    while True:
        try:
            return fn()
        except BaseException as e:  # noqa: BLE001
            if is_quota_error(e):
                log(f"quota error on {label}: sleeping {QUOTA_SLEEP:.0f}s")
                counter["quota_errors"] = counter.get("quota_errors", 0) + 1
                time.sleep(QUOTA_SLEEP)
                counter["resent_calls"] = counter.get("resent_calls", 0) + n_calls
                continue
            if is_transient_error(e) and transient < MAX_TRANSIENT:
                transient += 1
                wait = min(TRANSIENT_SLEEP * 2 ** (transient - 1), TRANSIENT_SLEEP_CAP)
                log(f"transient error {transient}/{MAX_TRANSIENT} on {label} ({e.code().name}): sleeping {wait:.0f}s")
                time.sleep(wait)
                counter["resent_calls"] = counter.get("resent_calls", 0) + n_calls
                continue
            raise


# ------------------------------------------------------------------ variants
def parse_uid(uid: str):
    chrom, pos, ref, alt = uid.split(":")
    return chrom, int(pos), ref.upper(), alt.upper()


def assert_ref(uids, fasta_path=FASTA):
    """Every REF must match the GRCh38 FASTA; returns the list of mismatches (empty = pass)."""
    import pyfaidx
    fa = pyfaidx.Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)
    bad = []
    for u in uids:
        c, p, r, a = parse_uid(u)
        seq = fa[c][p - 1:p - 1 + len(r)]
        if seq != r:
            bad.append((u, seq))
    return bad


def to_variant(uid: str):
    from alphagenome.data import genome
    c, p, r, a = parse_uid(uid)
    return genome.Variant(chromosome=c, position=p, reference_bases=r, alternate_bases=a, name=uid)


# ------------------------------------------------------------------ scoring
def score_chunked(client, version: str, uids: list[str], condition: str | None, scorer_names: list[str],
                  out_dir: pathlib.Path, chunk: int = 50, max_workers: int = 4, counter: dict | None = None,
                  min_seconds_per_chunk: float = 0.0, keep=None, cap: int | None = None, used_before: int = 0):
    """Score uids in fixed chunks; each chunk's tidy table is written once (parquet) and skipped on resume.

    version 'v2': one condition per call (condition required). version 'v1': condition None, all tracks.
    keep: optional function(tidy_df) -> tidy_df applied before writing (used only to drop rows of other
    genes for very wide v1 responses; the v2 tidy tables are written whole).
    cap: hosted-call ceiling; before each chunk, used_before + attempted + resent + chunk size must not exceed it.
    Returns list of chunk parquet paths.
    """
    import pandas as pd
    counter = counter if counter is not None else {}
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if version == "v2":
        from alphagenome.models.v2 import dna_client, variant_scorers as vs
    else:
        from alphagenome.models.v1 import dna_client, variant_scorers as vs
    scorers = [vs.RECOMMENDED_VARIANT_SCORERS[s] for s in scorer_names]
    paths = []
    for i in range(0, len(uids), chunk):
        part = uids[i:i + chunk]
        p = out_dir / f"chunk_{i // chunk:05d}.parquet"
        if p.exists():
            paths.append(p)
            continue
        if cap is not None:
            spent = used_before + counter.get("attempted_calls", 0) + counter.get("resent_calls", 0)
            if spent + len(part) > cap:
                raise SystemExit(f"hosted-call cap: {spent} used, next chunk {len(part)}, cap {cap}; stopping")
        variants = [to_variant(u) for u in part]
        intervals = [v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB) for v in variants]
        t0 = time.time()
        call_t = utc()
        if version == "v2":
            fn = lambda: client.score_variants(intervals=intervals, variants=variants, conditions=[condition] * len(part),
                                               variant_scorers=scorers, progress_bar=False, max_workers=max_workers)
        else:
            fn = lambda: client.score_variants(intervals=intervals, variants=variants, variant_scorers=scorers,
                                               progress_bar=False, max_workers=max_workers)
        counter["attempted_calls"] = counter.get("attempted_calls", 0) + len(part)
        res = call_with_retry(fn, f"chunk {i // chunk}", counter, n_calls=len(part))
        elapsed = time.time() - t0
        tidy = vs.tidy_scores(res, match_gene_strand=True, include_extended_metadata=True)
        if tidy is None:
            tidy = pd.DataFrame()
        if keep is not None and len(tidy):
            tidy = keep(tidy)
        tidy = tidy.copy()
        for col in tidy.columns:
            if tidy[col].dtype == object:
                tidy[col] = tidy[col].astype(str)
        tidy["condition_id"] = condition if condition is not None else "v1_all_tracks"
        tidy["client_version"] = version
        tidy["package_version"] = package_version()
        tidy["model_version"] = "ALL_FOLDS (requested; no version string returned)" if version == "v2" else "v1 server default (no version string returned)"
        tidy["call_utc"] = call_t
        tmp = p.with_suffix(".tmp")
        tidy.to_parquet(tmp, index=False)
        os.replace(tmp, p)
        meta = {"chunk": i // chunk, "n_variants": len(part), "elapsed_s": elapsed, "call_utc": call_t,
                "condition": condition, "scorers": scorer_names, "rows": int(len(tidy)), "bytes_parquet": p.stat().st_size}
        with open(out_dir / "chunk_log.jsonl", "a") as fh:
            fh.write(json.dumps(meta) + "\n")
        counter["completed_calls"] = counter.get("completed_calls", 0) + len(part)
        log(f"{out_dir.name} chunk {i // chunk}: {len(part)} variants in {elapsed:.1f}s ({len(part) / max(elapsed, 1e-9):.2f}/s), {len(tidy)} rows")
        paths.append(p)
        if elapsed < min_seconds_per_chunk:
            time.sleep(min_seconds_per_chunk - elapsed)
    return paths
