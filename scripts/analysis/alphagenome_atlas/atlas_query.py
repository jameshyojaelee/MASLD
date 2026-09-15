"""Throttled, resumable Atlas retrieval (the only module that talks to the network)."""

from __future__ import annotations

import datetime as dt
import importlib.metadata
import pathlib
import time

import grpc
from alphagenome.atlas import atlas
from alphagenome.data import genome

import atlas_archive as aa
import lib_atlas as la

CHUNK = 250
MIN_SECONDS_PER_CHUNK = 12.0     # ~1,250 requests/min; the observed per-minute quota tripped near 2,000
QUOTA_SLEEP = 65.0
MAX_WORKERS = 6


def to_variant(uid: str) -> genome.Variant:
    chrom, pos, ref, alt = uid.split(":")
    return genome.Variant(chromosome=chrom, position=int(pos), reference_bases=ref, alternate_bases=alt, name=uid)


def create_client():
    key, src = la.load_api_key()
    la.log(f"API key from {src} (len={len(key)})")
    # create_client_with_large_messages rather than atlas.create: the stock client caps received frames at
    # 4 MB, which any interval query over the wide TF/histone scorers exceeds.
    return create_client_with_large_messages(key), importlib.metadata.version("alphagenome")


MAX_MESSAGE_BYTES = 512 * 1024 * 1024   # gRPC defaults to 4 MB; a 2-kb window over the wide TF/histone
                                        # scorers returns tens of MB, so the client must accept a large frame.


def is_message_size_error(error: BaseException) -> bool:
    """An oversized RESPONSE, which gRPC reports with the same status code as a quota breach.

    Retrying it can never succeed: the same request returns the same too-large response. It must propagate
    (or trigger a window split), never enter the quota sleep loop.
    """
    if not isinstance(error, grpc.RpcError):
        return False
    if error.code() != grpc.StatusCode.RESOURCE_EXHAUSTED:
        return False
    return "message larger than max" in (error.details() or "")


def is_quota_error(error: BaseException) -> bool:
    return (isinstance(error, grpc.RpcError)
            and error.code() == grpc.StatusCode.RESOURCE_EXHAUSTED
            and not is_message_size_error(error))


def filter_tracks_to_ontology(res: dict, terms) -> tuple[dict, dict]:
    """Keep only tracks whose ontology_curie is in `terms`, preserving obs, var and every layer.

    Retrieval cannot filter server-side: passing a large ontology_terms list to query_interval fails with
    UNAVAILABLE. Filtering at archive time gives the same result for a saturation lane, where the analysis
    reads only the liver and proxy panel anyway, and cuts a region from 46.7 MB to 10.3 MB. A scorer with no
    ontology column (AVI_SCORE) is kept whole. Track counts before and after are returned so the filter is
    auditable from the deposit.
    """
    keep_terms = set(terms)
    out, counts = {}, {}
    for scorer, a in res.items():
        if a is None:
            out[scorer] = a
            continue
        n_before = int(a.shape[1])
        if "ontology_curie" not in getattr(a, "var", {}):
            out[scorer] = a
            counts[scorer] = {"tracks_returned": n_before, "tracks_archived": n_before}
            continue
        mask = a.var["ontology_curie"].isin(keep_terms).to_numpy()
        out[scorer] = a[:, mask].copy()
        counts[scorer] = {"tracks_returned": n_before, "tracks_archived": int(mask.sum())}
    return out, counts


def query_interval_split(query, start: int, end: int, min_width: int = 64) -> list:
    """Call query(start, end); if the response is too large, halve the window and recurse.

    Returns the list of per-piece results in coordinate order. Refuses to split below min_width so a request
    that is oversized for some other reason fails loudly instead of recursing forever.
    """
    try:
        return [query(start, end)]
    except grpc.RpcError as error:
        if not is_message_size_error(error) or (end - start) <= min_width:
            raise
    mid = start + (end - start) // 2
    return query_interval_split(query, start, mid, min_width) + query_interval_split(query, mid, end, min_width)


def create_client_with_large_messages(api_key: str, max_message_bytes: int = MAX_MESSAGE_BYTES):
    """atlas.create, but with the gRPC receive limit raised.

    atlas.create builds its channel with only the service config, leaving grpc.max_receive_message_length at
    the 4 MB default. This rebuilds the same three public objects (channel, stub, AtlasClient) with the limit
    raised; everything else matches atlas.create exactly.
    """
    import importlib.resources
    from alphagenome.protos import atlas_service_pb2_grpc

    service_config = (importlib.resources.files("alphagenome") / "protos/grpc_service_config.json").read_text()
    channel = grpc.secure_channel(
        "dns:///gdmscience.googleapis.com:443",
        grpc.ssl_channel_credentials(),
        options=(("grpc.service_config", service_config),
                 ("grpc.max_receive_message_length", int(max_message_bytes))),
    )
    grpc.channel_ready_future(channel).result(None)
    stub = atlas_service_pb2_grpc.AtlasServiceStub(channel=channel)
    return atlas.AtlasClient(stub, metadata=[("x-goog-api-key", api_key)])


def is_transient_error(error: BaseException) -> bool:
    return isinstance(error, grpc.RpcError) and error.code() in (grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.DEADLINE_EXCEEDED, grpc.StatusCode.INTERNAL)


TRANSIENT_SLEEP = 30.0        # first backoff step
TRANSIENT_SLEEP_CAP = 300.0   # 5 minutes between attempts once backed off
MAX_TRANSIENT = 14            # >= 30 minutes of total outage survived (see transient_sleep_seconds)


def transient_sleep_seconds(attempt: int) -> float:
    """Backoff for the nth consecutive transient service error: 30 s doubling to a 5-minute cap.

    A real Atlas outage outlasted a flat 10 x 30 s budget and killed the saturation pilot; the identical
    request succeeded minutes later, so the retry budget was the problem, not the request shape. Doubling
    with a cap covers a materially longer outage without hammering a service that is already struggling.
    """
    return float(min(TRANSIENT_SLEEP * (2 ** max(0, attempt - 1)), TRANSIENT_SLEEP_CAP))


def call_with_quota_retry(fn, is_quota=is_quota_error, is_transient=is_transient_error, sleep=time.sleep, label: str = "request", max_transient: int = MAX_TRANSIENT):
    """Call fn(); on a per-minute-quota error sleep QUOTA_SLEEP and retry without limit; on a transient service error
    (UNAVAILABLE / DEADLINE_EXCEEDED / INTERNAL) sleep TRANSIENT_SLEEP and retry at most max_transient times; anything else propagates."""
    transient = 0
    while True:
        try:
            return fn()
        except BaseException as error:  # noqa: BLE001 - the predicates decide
            if is_quota(error):
                la.log(f"quota: sleeping {QUOTA_SLEEP}s before retrying {label}")
                sleep(QUOTA_SLEEP)
                continue
            if is_transient(error) and transient < max_transient:
                transient += 1
                wait = transient_sleep_seconds(transient)
                la.log(f"transient service error ({transient}/{max_transient}): sleeping {wait:.0f}s before retrying {label}")
                sleep(wait)
                continue
            raise


def query_archived(client, sdk: str, uids: list[str], scorers: list[str], out_dir: pathlib.Path, extra: dict | None = None) -> list[pathlib.Path]:
    """Query `uids` in fixed chunks; archive each chunk once; skip chunks already on disk."""
    out_dir = pathlib.Path(out_dir)
    done: list[pathlib.Path] = []
    for i in range(0, len(uids), CHUNK):
        chunk_dir = out_dir / f"chunk_{i // CHUNK:05d}"
        part = uids[i:i + CHUNK]
        if (chunk_dir / "request.json").exists():
            done.append(chunk_dir)
            continue
        if chunk_dir.exists():
            raise la.ContractError(f"partial chunk on disk, remove it by hand: {chunk_dir}")
        t1 = time.time()
        res = call_with_quota_retry(lambda: client.query_variants([to_variant(u) for u in part], requested_scorers=scorers, progress_bar=False, max_workers=MAX_WORKERS),
                                    label=f"chunk {i // CHUNK}")
        elapsed = time.time() - t1
        request = {"variants": part, "requested_scorers": scorers, "sdk_version": sdk, "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                   "elapsed_seconds": elapsed, "chunk_size": CHUNK, "max_workers": MAX_WORKERS, **(extra or {})}
        aa.archive_scores(res, chunk_dir, request)
        la.log(f"chunk {i // CHUNK}: {len(part)} variants in {elapsed:.1f}s")
        done.append(chunk_dir)
        if elapsed < MIN_SECONDS_PER_CHUNK:
            time.sleep(MIN_SECONDS_PER_CHUNK - elapsed)
    return done
