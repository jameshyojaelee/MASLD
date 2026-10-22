#!/usr/bin/env python3
"""H87 step 1: fetch Ensembl GRCh37 REST VEP and /map responses (A1, G4, G5, D2, D3 inputs).

Writes raw JSON responses with sha256 into --out and a fetch_status.json that step 2 reads.
Never raises on a network failure: it records the stop rule (S4 or S7) in fetch_status.json.
"""

from __future__ import annotations

import argparse
import csv
import json
import socket
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import h87_common as C  # noqa: E402

RETRIES = 3          # prespec A1: up to 3 retries per batch, 60 s apart
RETRY_SLEEP = 60
MAX_BATCH = 200
TIMEOUT = 300


def http(method: str, url: str, body: dict | None = None) -> tuple[int, bytes, dict]:
    data = None if body is None else json.dumps(body, separators=(",", ":")).encode()
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return resp.status, resp.read(), dict(resp.headers)
    except urllib.error.HTTPError as err:
        return err.code, err.read() if hasattr(err, "read") else b"", dict(err.headers or {})


def with_retries(method, url, body=None, sleep=RETRY_SLEEP, log=None):
    last = None
    for attempt in range(RETRIES + 1):
        try:
            status, raw, headers = http(method, url, body)
            if status == 200:
                return raw, attempt, None
            last = f"HTTP {status}: {raw[:200]!r}"
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError, OSError) as err:
            last = f"{type(err).__name__}: {err}"
        if log:
            log(f"  attempt {attempt + 1} failed for {url}: {last}")
        if attempt < RETRIES:
            time.sleep(sleep)
    return None, RETRIES + 1, last


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--map-sleep", type=float, default=0.12)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise SystemExit(f"refusing: {out} is not empty")
    (out / "map").mkdir()
    logf = open(out / "fetch.log", "x")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    status = {"stop": None, "host": socket.gethostname(), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    manifest = []

    def save(name: str, raw: bytes, endpoint: str, extra=""):
        path = out / name
        with open(path, "xb") as handle:
            handle.write(raw)
        manifest.append(dict(file=name, bytes=len(raw), sha256=C.sha256_bytes(raw), endpoint=endpoint, note=extra))

    def finish():
        status["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with open(out / "fetch_manifest.tsv", "x", newline="") as handle:
            w = csv.DictWriter(handle, ["file", "bytes", "sha256", "endpoint", "note"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(manifest)
        C.write_json(out / "fetch_status.json", status)
        log(f"FETCH_STATUS stop={status['stop']}")

    # Prespec and the two inputs this step reads.
    C.verify_prespec()
    chk = C.check_inputs({"members", "repo_vep"})
    for row in chk:
        log(f"input {row['role']} ok={row['ok']} sha256={row['observed_sha256']}")
    if not all(r["ok"] for r in chk):
        status["stop"] = "S1: input sha256 differs: " + ",".join(r["role"] for r in chk if not r["ok"])
        finish()
        return

    members = read_tsv(C.INPUTS["members"][0])
    repo = {r["variant_id"]: r for r in read_tsv(C.INPUTS["repo_vep"][0])}

    gated = [r for r in members if r["pip_sum_gate_passed"] == "true"]
    ungated = [r for r in members if r["pip_sum_gate_passed"] == "false"]
    prim_ids = sorted({r["variant_id"] for r in gated if r["consequence_category"] in C.CODING_CATS})
    nc3_ids = sorted({r["variant_id"] for r in gated if r["consequence_category"] == "synonymous_or_utr"
                      and float(r["normalized_susie_pip"]) >= 0.5})
    d6_ids = sorted({r["variant_id"] for r in ungated if r["consequence_category"] in C.CODING_CATS} - set(prim_ids))
    log(f"primary coding variant_ids {len(prim_ids)}; NC-3 variant_ids {len(nc3_ids)}; D6 ungated-only coding variant_ids {len(d6_ids)}")
    status.update(n_primary=len(prim_ids), n_nc3=len(nc3_ids), n_d6=len(d6_ids))

    def vcf(vid):
        r = repo[vid]
        return f"{r['chromosome']} {r['position']} {vid} {r['reference_allele']} {r['alternate_allele']} . . ."

    missing = {k: [v for v in ids if v not in repo] for k, ids in
               (("primary", prim_ids), ("nc3", nc3_ids), ("d6", d6_ids))}
    status["missing_from_repo_vep"] = missing
    if missing["primary"]:
        status["stop"] = f"S2: G3 failed, {len(missing['primary'])} primary variant_ids lack oriented alleles"
        finish()
        return

    # G4 ping and release info.
    raw, _, err = with_retries("GET", C.ENSEMBL37 + "/info/ping", sleep=60, log=log)
    ping_ok = raw is not None and json.loads(raw).get("ping") == 1
    status["g4_ping"] = ping_ok
    log(f"G4 /info/ping ok={ping_ok} err={err}")
    if not ping_ok:
        status["stop"] = f"S4: Ensembl GRCh37 REST unreachable from {status['host']}: {err}"
        finish()
        return
    save("ping.json", raw, "/info/ping")
    for ep, name in (("/info/software", "info_software.json"), ("/info/data", "info_data.json")):
        raw, _, err = with_retries("GET", C.ENSEMBL37 + ep, sleep=10, log=log)
        if raw is not None:
            save(name, raw, ep)
        else:
            log(f"WARNING {ep} failed: {err}")

    # VEP POST batches.
    def post_set(label, ids, lines=None, stop_on_fail=True):
        lines = lines or [vcf(v) for v in ids]
        recs = []
        for b, start in enumerate(range(0, len(lines), MAX_BATCH), 1):
            chunk = lines[start:start + MAX_BATCH]
            want = [ln.split()[2] for ln in chunk]
            body = {"variants": chunk, "canonical": 1}
            raw, attempts, err = with_retries("POST", C.ENSEMBL37 + "/vep/homo_sapiens/region", body, log=log)
            if raw is None:
                msg = f"S4: VEP batch {label}#{b} failed after {RETRIES} retries: {err}"
                log(msg)
                if stop_on_fail:
                    status["stop"] = msg
                recs.append(dict(label=label, batch=b, ok=False, error=err))
                return recs
            name = f"vep_{label}_batch{b}.json"
            save(name, raw, "/vep/homo_sapiens/region", f"canonical=1; n_submitted={len(chunk)}")
            with open(out / f"vep_{label}_batch{b}_request.json", "x") as handle:
                json.dump(body, handle, indent=1)
            got = [str(item.get("id", "")) for item in json.loads(raw)]
            ident = sorted(got) == sorted(want)
            recs.append(dict(label=label, batch=b, ok=True, file=name, n_submitted=len(chunk), n_returned=len(got),
                             identity_ok=ident, attempts=attempts,
                             missing=sorted(set(want) - set(got)), extra=sorted(set(got) - set(want))))
            log(f"VEP {label} batch {b}: submitted {len(chunk)} returned {len(got)} identity_ok={ident}")
            if not ident and stop_on_fail:
                status["stop"] = f"S7: VEP response identity mismatch in {label}#{b}"
                return recs
        return recs

    status["vep_batches"] = post_set("primary", prim_ids)
    if status["stop"]:
        finish()
        return
    status["vep_batches"] += post_set("nc3", [v for v in nc3_ids if v in repo], stop_on_fail=False)
    status["vep_batches"] += post_set("d6", [v for v in d6_ids if v in repo], stop_on_fail=False)
    indel_lines = ["4 88231392 4:88231392:T:TA_del TA T . . .", "4 88231392 4:88231392:T:TA_ins T TA . . ."]
    status["vep_batches"] += post_set("d2indel", None, lines=indel_lines, stop_on_fail=False)

    # D3: GRCh37 -> GRCh38 position of every gated coding site (and the deleted base of the one indel).
    regions = sorted({(repo[v]["chromosome"], int(repo[v]["position"])) for v in prim_ids},
                     key=lambda x: (x[0], x[1]))
    extra = []
    for v in prim_ids:
        r = repo[v]
        ref, alt = r["reference_allele"], r["alternate_allele"]
        if len(ref) != len(alt):
            pos = int(r["position"])
            if len(ref) > len(alt):
                extra += [(r["chromosome"], p) for p in range(pos + 1, pos + len(ref))]
            else:
                extra += [(r["chromosome"], pos), (r["chromosome"], pos + 1)]
    regions = sorted(set(regions) | set(extra), key=lambda x: (x[0], x[1]))
    map_fail = []
    for chrom, pos in regions:
        name = f"map/{chrom}_{pos}.json"
        got = None
        for server in (C.ENSEMBL37, C.ENSEMBL_MAIN):
            url = f"{server}/map/human/GRCh37/{chrom}:{pos}..{pos}:1/GRCh38"
            raw, _, err = with_retries("GET", url, sleep=5)
            if raw is not None:
                got = (raw, server)
                break
        if got is None:
            map_fail.append(dict(site=f"{chrom}:{pos}", error=err))
            log(f"MAP failed {chrom}:{pos}: {err}")
        else:
            save(name, got[0], f"/map GRCh37->GRCh38 via {got[1]}")
        time.sleep(args.map_sleep)
    status["map_regions"] = len(regions)
    status["map_failures"] = map_fail
    log(f"MAP regions {len(regions)} failures {len(map_fail)}")
    finish()


if __name__ == "__main__":
    main()
