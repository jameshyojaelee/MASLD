#!/usr/bin/env python3
"""Read-only local molecular-evidence queries. No individual disease-risk score.

The database is a candidate attachment to the Gene Catalog. Units and source
records remain separate; missing evidence never becomes a zero effect.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sqlite3

KINDS = ("variant", "gene", "region", "event", "program", "context")


def variant_id(value):
    match = re.fullmatch(r"(GRCh38|hg38):(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M):([1-9][0-9]*):([ACGT]+):([ACGT]+)", value)
    if not match:
        raise ValueError("Provide an explicit GRCh38:chromosome:1-based-position:REF:ALT identity; no rsID or ambiguous allele resolution is inferred")
    _, chrom, pos, ref, alt = match.groups()
    if ref == alt:
        raise ValueError("Reference and alternate alleles are identical; no variant contrast")
    return f"GRCh38:chr{chrom}:{int(pos)}:{ref}:{alt}"


def connect(path):
    conn = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only=ON")
    return conn


def query(path, kind, identity, assay=None, limit=100):
    if kind not in KINDS:
        raise ValueError("Unknown molecular entry point")
    if limit < 1 or limit > 1000:
        raise ValueError("Limit must be 1–1000; total counts remain visible")
    if kind == "variant":
        identity = variant_id(identity)
    with connect(path) as conn:
        entities = conn.execute("SELECT * FROM entity WHERE kind=? AND (identifier=? OR display_name=?) ORDER BY identifier",
                                (kind, identity, identity)).fetchall()
        answer = {"entry_point": kind, "query": identity, "release_state": "candidate_not_adopted",
                  "entities": [], "evidence": [], "relationships": [], "total_evidence": 0,
                  "status": "unsupported", "unsupported_reason": "identity_absent_from_bounded_specimen",
                  "inference": "stored_molecular_evidence_only_no_new_model_call",
                  "disease_risk_score": "not_supported"}
        if not entities:
            return answer
        keys = [r["entity_key"] for r in entities]
        answer["entities"] = [{**dict(r), "attributes": json.loads(r["attributes"])} for r in entities]
        if len(keys) > 1:
            answer["unsupported_reason"] = "ambiguous_identifier_use_exact_versioned_identity"
            return answer
        key = keys[0]
        clause = " AND e.assay=?" if assay else ""
        values = (key, assay) if assay else (key,)
        sql = "FROM evidence e JOIN evidence_entity ee USING(evidence_id) WHERE ee.entity_key=?" + clause
        answer["total_evidence"] = conn.execute("SELECT count(DISTINCT e.evidence_id) " + sql, values).fetchone()[0]
        records = conn.execute("SELECT DISTINCT e.* " + sql + " ORDER BY e.evidence_id LIMIT ?", (*values, limit)).fetchall()
        for row in records:
            record = dict(row)
            record["metadata"] = json.loads(record["metadata"])
            record["source"] = dict(conn.execute("SELECT * FROM source WHERE source_id=?", (record["source_id"],)).fetchone())
            record["objects"] = [dict(r) for r in conn.execute(
                "SELECT en.kind,en.identifier,en.display_name,ee.role FROM evidence_entity ee JOIN entity en USING(entity_key) WHERE ee.evidence_id=? ORDER BY en.kind,en.identifier",
                (record["evidence_id"],))]
            answer["evidence"].append(record)
        answer["relationships"] = [dict(r) for r in conn.execute(
            "SELECT r.*,a.kind AS from_kind,a.identifier AS from_id,b.kind AS to_kind,b.identifier AS to_id FROM relationship r JOIN entity a ON a.entity_key=r.from_key JOIN entity b ON b.entity_key=r.to_key WHERE r.from_key=? OR r.to_key=? ORDER BY r.link_type,a.identifier,b.identifier", (key, key))]
        answer["status"] = "stored_evidence_available" if records else "unsupported"
        answer["unsupported_reason"] = "" if records else ("assay_not_available_for_this_identity" if assay else "no_effect_measurement_or_prediction_in_specimen")
        answer["truncated"] = answer["total_evidence"] > len(records)
        answer["interpretation"] = "Keep source-native units and conflicting targets; relationships are not proof of regulation"
        return answer


def haplotype_query(path, variants, phase, phase_source):
    identities = [variant_id(v) for v in variants]
    if len(identities) < 2 or len(set(identities)) != len(identities):
        raise ValueError("Supply at least two distinct variants")
    loci = [v.split(":") for v in identities]
    if len({r[1] for r in loci}) != 1:
        raise ValueError("A local haplotype must use one chromosome")
    ordered = sorted((int(r[2]), int(r[2]) + len(r[3])) for r in loci)
    if any(a[1] > b[0] for a, b in zip(ordered, ordered[1:])):
        raise ValueError("Overlapping variant edits require an explicitly normalized combined allele")
    known = phase == "same_haplotype" and bool(phase_source and phase_source.strip())
    return {"entry_point": "haplotype", "variants": identities, "phase": phase,
            "phase_source": phase_source, "status": "unsupported",
            "unsupported_reason": "joint_effect_not_stored_in_specimen" if known else "unknown_phase_no_joint_sequence_constructed",
            "single_variant_evidence": [query(path, "variant", v) for v in identities],
            "joint_effect": None, "uncertainty": None, "disease_risk_score": "not_supported",
            "interpretation": "Single-variant effects cannot be added to establish a haplotype effect",
            "local_inference_route": "catalog_infer.py; declared phase and supported SNV/reference mapping required"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    sub = parser.add_subparsers(dest="command", required=True)
    q = sub.add_parser("query")
    q.add_argument("kind", choices=KINDS)
    q.add_argument("identity")
    q.add_argument("--assay")
    q.add_argument("--limit", type=int, default=100)
    h = sub.add_parser("haplotype")
    h.add_argument("variants", nargs="+")
    h.add_argument("--phase", choices=["unknown", "same_haplotype"], default="unknown")
    h.add_argument("--phase-source", default="")
    args = parser.parse_args()
    result = query(args.db, args.kind, args.identity, args.assay, args.limit) if args.command == "query" else haplotype_query(args.db, args.variants, args.phase, args.phase_source)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
