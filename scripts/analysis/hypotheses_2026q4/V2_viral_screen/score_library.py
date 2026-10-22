#!/usr/bin/env python3
"""V2 target screen, step 4: candidate extraction and per-library scoring.

Subcommand `candidates` reads the target-hit SAM, collapses alignments to unique
read sequences, and writes a FASTA of those unique sequences for the human
re-check plus a TSV that records, for every unique sequence, the target contig,
position, CIGAR and alignment score of its best target alignment.

Subcommand `score` reads that TSV, the human re-check SAM and the HISAT2
library summary, drops every candidate whose human alignment score is greater
than or equal to its target alignment score, and applies the positivity rule
fixed in PRESPEC_V2.md:

  called positive for a target when ALL of
    retained fragments           >= 20
    covered unmasked bases       >= 500 on the single best-hit contig
    distinct 200-bp windows hit  >= 3  on that contig
    fragments per million reads  >= 0.20
  indeterminate when at least one criterion is met but not all
  negative otherwise

Fragments are counted once per read name, so a mate pair counts once.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import sys

MIN_FRAGMENTS = 20
MIN_COVERED_BASES = 500
MIN_WINDOWS_200BP = 3
MIN_RPM = 0.20
WINDOW_BP = 200

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_span(cigar: str) -> int:
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def aln_score(fields: list[str]) -> int:
    for tag in fields[11:]:
        if tag.startswith("AS:i:"):
            return int(tag[5:])
    return -10 ** 9


def iter_sam(path: str):
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            yield f


def cmd_candidates(sam: str, out_fasta: str, out_tsv: str) -> int:
    # best target alignment per unique read sequence
    best: dict[str, dict] = {}
    for f in iter_sam(sam):
        flag = int(f[1])
        if flag & 0x4 or flag & 0x100 or flag & 0x800:
            continue
        seq = f[9].upper()
        if not seq or seq == "*":
            continue
        # store the sequence as written in the SAM record (reference strand)
        rec = {"qname": f[0], "flag": flag, "contig": f[2], "pos": int(f[3]),
               "mapq": int(f[4]), "cigar": f[5], "seq": seq,
               "viral_as": aln_score(f), "is_read2": bool(flag & 0x80)}
        prev = best.get(seq)
        if prev is None or rec["viral_as"] > prev["viral_as"]:
            best[seq] = rec

    cols = ["seq_id", "qname", "contig", "pos", "mapq", "cigar", "viral_as",
            "read_len", "is_read2", "n_alignments_collapsed"]
    counts = collections.Counter()
    for f in iter_sam(sam):
        if int(f[1]) & 0x4:
            continue
        counts[f[9].upper()] += 1

    with open(out_fasta, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("\t".join(cols) + "\n")
        for i, (seq, rec) in enumerate(sorted(best.items())):
            seq_id = "cand%06d" % i
            fa.write(">%s\n%s\n" % (seq_id, seq))
            tsv.write("\t".join(str(x) for x in [
                seq_id, rec["qname"], rec["contig"], rec["pos"], rec["mapq"],
                rec["cigar"], rec["viral_as"], len(seq),
                int(rec["is_read2"]), counts[seq]]) + "\n")
    sys.stderr.write("candidates: %d unique sequences\n" % len(best))
    return 0


def parse_summary(path: str) -> int:
    """Total reads or pairs processed, from the HISAT2 summary file."""
    total = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"^(\d+)\s+reads; of these:$", line)
            if m:
                total = int(m.group(1))
                break
    return total


def load_slots(acc_table: str) -> dict[str, str]:
    """contig accession -> target label (target class B/target class C/target class E/CONTROL-phiX/...)."""
    out: dict[str, str] = {}
    with open(acc_table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row["accession"]
            if not acc:
                continue
            slot = row["slot"]
            if slot.startswith("CONTROL_"):
                label = slot
            else:
                label = slot.split("_")[0]
            out.setdefault(acc, label)
    return out


def load_usable(mask_report: str) -> dict[str, int]:
    out: dict[str, int] = {}
    with open(mask_report) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["contig"]] = int(row["usable_length"])
    return out


def cmd_score(args) -> int:
    slots = load_slots(args.accession_table)
    usable = load_usable(args.mask_report)
    total_reads = parse_summary(args.library_summary)

    cands: dict[str, dict] = {}
    with open(args.candidate_table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["pos"] = int(row["pos"])
            row["viral_as"] = int(row["viral_as"])
            row["mapq"] = int(row["mapq"])
            row["n_alignments_collapsed"] = int(row["n_alignments_collapsed"])
            cands[row["seq_id"]] = row

    human_as: dict[str, int] = {}
    for f in iter_sam(args.human_sam):
        if int(f[1]) & 0x4:
            continue
        s = aln_score(f)
        if s > human_as.get(f[0], -10 ** 9):
            human_as[f[0]] = s

    # Drop any candidate that the human genome explains at least as well.
    retained = []
    dropped_human = 0
    dropped_mapq = 0
    for seq_id, rec in cands.items():
        if rec["mapq"] < args.min_mapq:
            dropped_mapq += 1
            continue
        h = human_as.get(seq_id)
        if h is not None and h >= rec["viral_as"]:
            dropped_human += 1
            continue
        retained.append(rec)

    # Per-target aggregation on the single best-hit contig.
    by_contig: dict[str, list[dict]] = collections.defaultdict(list)
    for rec in retained:
        by_contig[rec["contig"]].append(rec)

    virus_contigs: dict[str, list[str]] = collections.defaultdict(list)
    for contig in by_contig:
        virus_contigs[slots.get(contig, "UNASSIGNED")].append(contig)

    results = []
    for virus in sorted(set(list(slots.values()))):
        contigs = virus_contigs.get(virus, [])
        if not contigs:
            results.append({
                "virus": virus, "best_contig": "", "fragments": 0,
                "unique_sequences": 0, "covered_bases": 0, "windows_200bp": 0,
                "usable_length": 0, "rpm": 0.0, "call": "negative",
                "criteria_met": 0,
            })
            continue
        # best-hit contig = most distinct read names
        def nfrag(c: str) -> int:
            return len({r["qname"] for r in by_contig[c]})
        best_contig = max(contigs, key=nfrag)
        recs = by_contig[best_contig]
        frags = len({r["qname"] for r in recs})
        covered: set[int] = set()
        for r in recs:
            covered.update(range(r["pos"], r["pos"] + ref_span(r["cigar"])))
        windows = {p // WINDOW_BP for p in covered}
        rpm = (frags / (total_reads / 1e6)) if total_reads else 0.0
        met = sum([frags >= MIN_FRAGMENTS,
                   len(covered) >= MIN_COVERED_BASES,
                   len(windows) >= MIN_WINDOWS_200BP,
                   rpm >= MIN_RPM])
        if met == 4:
            call = "positive"
        elif met >= 1:
            call = "indeterminate"
        else:
            call = "negative"
        results.append({
            "virus": virus, "best_contig": best_contig, "fragments": frags,
            "unique_sequences": len(recs), "covered_bases": len(covered),
            "windows_200bp": len(windows),
            "usable_length": usable.get(best_contig, 0),
            "rpm": round(rpm, 4), "call": call, "criteria_met": met,
        })

    out = {
        "library": args.library, "cohort": args.cohort,
        "total_reads_or_pairs": total_reads,
        "candidate_unique_sequences": len(cands),
        "dropped_low_mapq": dropped_mapq,
        "dropped_human_homologous_read": dropped_human,
        "retained_candidates": len(retained),
        "rule": {"min_fragments": MIN_FRAGMENTS,
                 "min_covered_bases": MIN_COVERED_BASES,
                 "min_windows_200bp": MIN_WINDOWS_200BP,
                 "min_rpm": MIN_RPM, "min_mapq": args.min_mapq},
        "per_virus": results,
    }
    with open(args.out_json, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)

    cols = ["library", "cohort", "virus", "best_contig", "fragments",
            "unique_sequences", "covered_bases", "windows_200bp",
            "usable_length", "rpm", "criteria_met", "call",
            "total_reads_or_pairs"]
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results:
            r2 = dict(r)
            r2["library"] = args.library
            r2["cohort"] = args.cohort
            r2["total_reads_or_pairs"] = total_reads
            fh.write("\t".join(str(r2[c]) for c in cols) + "\n")
    sys.stderr.write("scored %s: %s\n" % (
        args.library,
        ", ".join("%s=%s" % (r["virus"], r["call"]) for r in results)))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("candidates")
    c.add_argument("--viral-sam", required=True)
    c.add_argument("--out-fasta", required=True)
    c.add_argument("--out-tsv", required=True)

    s = sub.add_parser("score")
    s.add_argument("--library", required=True)
    s.add_argument("--cohort", required=True)
    s.add_argument("--candidate-table", required=True)
    s.add_argument("--human-sam", required=True)
    s.add_argument("--library-summary", required=True)
    s.add_argument("--accession-table", required=True)
    s.add_argument("--mask-report", required=True)
    s.add_argument("--min-mapq", type=int, default=1)
    s.add_argument("--out-json", required=True)
    s.add_argument("--out-tsv", required=True)

    args = ap.parse_args()
    if args.cmd == "candidates":
        return cmd_candidates(args.viral_sam, args.out_fasta, args.out_tsv)
    return cmd_score(args)


if __name__ == "__main__":
    sys.exit(main())
