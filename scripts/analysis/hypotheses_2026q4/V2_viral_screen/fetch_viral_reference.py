#!/usr/bin/env python3
"""V2 target screen, step 1: assemble the target reference set from NCBI.

Writes, into an otherwise empty download directory:
  viral_refs.fasta          concatenated sequences, one record per accession
  viral_refs_accessions.tsv accession, requested slot, returned DEFINITION, length, sha256

The slot list is fixed in PRESPEC_V2.md. Accessions are not hard-coded by the
author where a genotype representative is needed: they are resolved by an NCBI
esearch query whose text is recorded, and whatever NCBI returns is written to
the table verbatim. Anchor RefSeq accessions are hard-coded because they are
the stable RefSeq records for these species.

Downloaded sequence is untrusted data. Nothing here executes it; it is only
read as text and written to FASTA.

Usage:
  python3 fetch_viral_reference.py --out-dir <empty download dir>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Anchor RefSeq records (species-level reference genomes) and controls.
ANCHORS = [
    ("HBV_refseq_anchor", "NC_003977.2"),
    ("HCV_refseq_anchor", "NC_004102.1"),
    ("HEV_refseq_anchor", "NC_001434.1"),
    ("CONTROL_phiX", "NC_001422.1"),
    ("CONTROL_HPV16", "NC_001526.4"),
]

# Genotype slots resolved by query. Length windows keep the hit on a complete
# genome rather than a partial gene record. `alt_org` is a second organism node
# tried when the first returns nothing (NCBI has moved target class E under
# Paslahepevirus, and some target class E genotype records sit under
# "human/swine target class E genotype Nx" organism names).
GENOTYPE_SLOTS = [
    # (slot, organism, alt organisms, genotype token, min length, max length)
    *[("HBV_genotype_%s" % g, "Hepatitis B virus", [], g, 3000, 3400)
      for g in list("ABCDEFGH")],
    *[("HCV_genotype_%s" % g, "Hepatitis C virus", [], g, 9000, 9800)
      for g in ["1", "2", "3", "4", "5", "6", "7"]],
    *[("HEV_genotype_%s" % g, "Paslahepevirus",
       ["human hepatitis E virus genotype %sa" % g,
        "swine hepatitis E virus genotype %sa" % g,
        "Hepatitis E virus"], g, 6600, 7400)
      for g in ["1", "2", "3", "4"]],
]

# Patent-derived sequence records are not usable as species references.
PATENT_TITLE = re.compile(r"(^[A-Z]{2}\s*\d{4,})|(-A/\d)|(\bpatent\b)", re.I)


def _get(url: str, retries: int = 4) -> str:
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:  # network flake
            last = exc
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("NCBI request failed: %s (%s)" % (url, last))


def esearch(term: str, retmax: int = 20) -> list[str]:
    url = "%s/esearch.fcgi?%s" % (EUTILS, urllib.parse.urlencode({
        "db": "nuccore", "term": term, "retmax": retmax, "retmode": "json",
        "sort": "relevance"}))
    data = json.loads(_get(url))
    return data.get("esearchresult", {}).get("idlist", [])


def esummary_accessions(uids: list[str]) -> list[tuple[str, str, int, str]]:
    if not uids:
        return []
    url = "%s/esummary.fcgi?%s" % (EUTILS, urllib.parse.urlencode({
        "db": "nuccore", "id": ",".join(uids), "retmode": "json"}))
    data = json.loads(_get(url)).get("result", {})
    out = []
    for uid in data.get("uids", []):
        rec = data.get(uid, {})
        out.append((rec.get("accessionversion", ""), rec.get("title", ""),
                    int(rec.get("slen") or 0), rec.get("organism", "")))
    return out


def efetch_fasta(accession: str) -> str:
    url = "%s/efetch.fcgi?%s" % (EUTILS, urllib.parse.urlencode({
        "db": "nuccore", "id": accession, "rettype": "fasta",
        "retmode": "text"}))
    text = _get(url)
    if not text.startswith(">"):
        raise RuntimeError("non-FASTA reply for %s" % accession)
    return text


def resolve_slot(slot: str, organism: str, alt_orgs: list[str], genotype: str,
                 lo: int, hi: int) -> tuple[str, str, int, str]:
    """Return (accession, title, length, query) for one genotype slot.

    Query templates are tried in a fixed order and the first hit that passes
    the acceptance test wins. The accepted record must be inside the length
    window, must carry the genotype as a standalone token in its title or its
    organism name (optionally followed by a subtype letter), and must not be a
    patent or UNVERIFIED record.
    """
    # Genotype token must stand alone or be followed by a subtype letter:
    # "genotype 3", "genotype 3a", "genotype 1a" all pass for genotype 3/1.
    pat = re.compile(r"genotype[s]?\s+%s(?![0-9])" % re.escape(genotype), re.I)
    orgs = [organism] + list(alt_orgs)
    templates = [
        '{org}[Organism] AND genotype {g}[Title] AND complete genome[Title] '
        'AND {lo}:{hi}[SLEN] NOT srcdb_patent[PROP]',
        '{org}[Organism] AND genotype {g}[Title] AND {lo}:{hi}[SLEN] '
        'NOT srcdb_patent[PROP]',
        '"{org}"[Organism] AND {lo}:{hi}[SLEN] NOT srcdb_patent[PROP]',
    ]
    tried = []
    for org in orgs:
        for tpl in templates:
            term = tpl.format(org=org, g=genotype, lo=lo, hi=hi)
            tried.append(term)
            hits = esummary_accessions(esearch(term))
            time.sleep(0.4)
            for acc, title, slen, org_name in hits:
                text = "%s %s" % (title, org_name)
                if (lo <= slen <= hi and pat.search(text)
                        and "UNVERIFIED" not in title.upper()
                        and not PATENT_TITLE.search(title)):
                    return acc, title, slen, term
    return "", "", 0, " || ".join(tried)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)
    fasta_path = os.path.join(args.out_dir, "viral_refs.fasta")
    table_path = os.path.join(args.out_dir, "viral_refs_accessions.tsv")

    rows = []
    records = []

    for slot, acc in ANCHORS:
        fa = efetch_fasta(acc)
        head, seq = fa.split("\n", 1)
        seq_clean = "".join(seq.split())
        rows.append({
            "slot": slot, "accession": acc, "definition": head[1:].strip(),
            "length": len(seq_clean), "resolution": "hard_coded_refseq_anchor",
            "query": "", "sha256": hashlib.sha256(seq_clean.encode()).hexdigest(),
        })
        records.append((acc, slot, seq_clean))
        time.sleep(0.4)

    for slot, organism, alt_orgs, genotype, lo, hi in GENOTYPE_SLOTS:
        acc, title, slen, term = resolve_slot(slot, organism, alt_orgs,
                                              genotype, lo, hi)
        if not acc:
            rows.append({"slot": slot, "accession": "", "definition": "",
                         "length": 0, "resolution": "unresolved",
                         "query": term, "sha256": ""})
            sys.stderr.write("UNRESOLVED %s\n" % slot)
            continue
        fa = efetch_fasta(acc)
        head, seq = fa.split("\n", 1)
        seq_clean = "".join(seq.split())
        rows.append({
            "slot": slot, "accession": acc, "definition": head[1:].strip(),
            "length": len(seq_clean), "resolution": "esearch_resolved",
            "query": term,
            "sha256": hashlib.sha256(seq_clean.encode()).hexdigest(),
        })
        records.append((acc, slot, seq_clean))
        time.sleep(0.4)

    # Drop duplicate accessions (an anchor may also win a genotype slot).
    seen = set()
    with open(fasta_path, "w") as out:
        for acc, slot, seq in records:
            if acc in seen:
                continue
            seen.add(acc)
            out.write(">%s %s\n" % (acc, slot))
            for i in range(0, len(seq), 70):
                out.write(seq[i:i + 70] + "\n")

    cols = ["slot", "accession", "definition", "length", "resolution",
            "query", "sha256"]
    with open(table_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    n_ok = sum(1 for r in rows if r["accession"])
    sys.stderr.write("wrote %d records (%d slots requested) to %s\n"
                     % (len(seen), len(rows), fasta_path))
    return 0 if n_ok == len(rows) else 2


if __name__ == "__main__":
    sys.exit(main())
