#!/usr/bin/env python3
"""V2 target screen: build the seeded synthetic positive control library.

The control is a host background plus a known number of target class B-derived reads, so
the prespecified positivity rule can be checked against a library whose true
answer is known. Nothing here is a biological claim; it exists only to show
that the rule fires when target sequence is present at a stated level.

  host background: the first --host-reads reads of an existing pilot FASTQ,
                   taken deterministically from the head of the file
  spike:           --spike-reads fragments drawn uniformly at random from the
                   unmasked positions of one target contig, random strand, with
                   a per-base substitution probability of --error-rate

The random number generator is seeded with --seed, so the control is
reproducible. Read names of spiked reads are prefixed SYNTH_SPIKE so they can
never be confused with real data.

Usage:
  make_positive_control.py --host-fastq <in.fastq.gz> --viral-fasta <refs.fasta>
      --contig <accession> --out-fastq <out.fastq.gz> --manifest <out.json>
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta_one(path: str, contig: str) -> str:
    seq: list[str] = []
    want = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if want:
                    break
                want = line[1:].split()[0] == contig
                continue
            if want:
                seq.append(line.strip())
    if not seq:
        raise SystemExit("contig %r not found in %s" % (contig, path))
    return "".join(seq).upper()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-fastq", required=True)
    ap.add_argument("--viral-fasta", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--out-fastq", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--host-reads", type=int, default=2_000_000)
    ap.add_argument("--spike-reads", type=int, default=60)
    ap.add_argument("--error-rate", type=float, default=0.002)
    ap.add_argument("--seed", type=int, default=20261007)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    genome = read_fasta_one(args.viral_fasta, args.contig)

    read_len = None
    n_host = 0
    with gzip.open(args.host_fastq, "rt") as fin, \
            gzip.open(args.out_fastq, "wt", compresslevel=6) as fout:
        while n_host < args.host_reads:
            name = fin.readline()
            if not name:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            if read_len is None:
                read_len = len(seq.strip())
            fout.write(name + seq + plus + qual)
            n_host += 1

        if read_len is None:
            raise SystemExit("host FASTQ is empty")

        spike_records = []
        for i in range(args.spike_reads):
            start = rng.randrange(0, max(1, len(genome) - read_len + 1))
            frag = genome[start:start + read_len]
            if len(frag) < read_len:          # circular reference wrap
                frag = (frag + genome)[:read_len]
            strand = "+" if rng.random() < 0.5 else "-"
            if strand == "-":
                frag = frag.translate(COMPLEMENT)[::-1]
            bases = list(frag)
            n_err = 0
            for j, b in enumerate(bases):
                if rng.random() < args.error_rate:
                    bases[j] = rng.choice([x for x in "ACGT" if x != b])
                    n_err += 1
            out_seq = "".join(bases)
            name = "SYNTH_SPIKE_%s_%06d" % (args.contig, i)
            fout.write("@%s\n%s\n+\n%s\n" % (name, out_seq, "I" * read_len))
            spike_records.append({"read": name, "contig": args.contig,
                                  "start": start, "strand": strand,
                                  "substitutions": n_err})

    manifest = {
        "host_fastq": args.host_fastq,
        "host_reads_written": n_host,
        "read_length": read_len,
        "viral_fasta": args.viral_fasta,
        "contig": args.contig,
        "contig_length": len(genome),
        "spike_reads": args.spike_reads,
        "error_rate": args.error_rate,
        "seed": args.seed,
        "expected_rpm": round(args.spike_reads / (n_host / 1e6), 4)
        if n_host else None,
        "out_fastq": args.out_fastq,
        "out_fastq_sha256": "",
        "spike_records": spike_records,
    }
    h = hashlib.sha256()
    with open(args.out_fastq, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    manifest["out_fastq_sha256"] = h.hexdigest()

    with open(args.manifest, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    sys.stderr.write("control: %d host reads + %d spiked reads, length %d, "
                     "expected RPM %.2f\n"
                     % (n_host, args.spike_reads, read_len,
                        manifest["expected_rpm"] or 0.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
