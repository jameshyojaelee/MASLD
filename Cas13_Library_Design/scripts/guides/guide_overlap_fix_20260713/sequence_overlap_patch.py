"""Exact-sequence, isoform-aware physical overlap filtering (isolated patch)."""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

GENE_RE = re.compile(r'gene_id "([^"]+)"')
TX_RE = re.compile(r'transcript_id "([^"]+)"')


def base(value: str) -> str:
    return value.rsplit(".", 1)[0] if "." in value else value


@dataclass
class Transcript:
    gene: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    sequence: str = ""

    def genomic_intervals(self, start0: int, length: int):
        """Map a 0-based transcript slice to inclusive genomic exon intervals."""
        stop0, cursor, result = start0 + length, 0, []
        for exon_start, exon_end in sorted(self.exons, reverse=self.strand == "-"):
            exon_length = exon_end - exon_start + 1
            left, right = max(start0, cursor), min(stop0, cursor + exon_length)
            if left < right:
                a, b = left - cursor, right - cursor - 1
                if self.strand == "+":
                    x, y = exon_start + a, exon_start + b
                else:
                    x, y = exon_end - b, exon_end - a
                result.append((self.chrom, min(x, y), max(x, y)))
            cursor += exon_length
            if cursor >= stop0:
                break
        covered = sum(y - x + 1 for _, x, y in result)
        return tuple(result) if covered == length else ()


class SequenceReference:
    def __init__(self, fasta: Path, gtf: Path, gene_ids: Iterable[str]):
        self.fasta, self.gtf = Path(fasta), Path(gtf)
        self.gene_ids = {base(str(x)) for x in gene_ids}
        self.transcripts: dict[str, Transcript] = {}
        self.by_gene: dict[str, list[str]] = defaultdict(list)
        self.candidate_rows = self.candidate_rows_located = self.candidate_rows_unlocated = 0
        self._models()
        self._sequences()

    def _models(self):
        with gzip.open(self.gtf, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                col = line.rstrip().split("\t")
                if len(col) != 9 or col[2] != "exon":
                    continue
                gm, tm = GENE_RE.search(col[8]), TX_RE.search(col[8])
                if not gm or not tm or base(gm.group(1)) not in self.gene_ids:
                    continue
                gene, tx = base(gm.group(1)), tm.group(1)
                if tx not in self.transcripts:
                    self.transcripts[tx] = Transcript(gene, col[0], col[6])
                    self.by_gene[gene].append(tx)
                self.transcripts[tx].exons.append((int(col[3]), int(col[4])))

    def _sequences(self):
        wanted, current, chunks = set(self.transcripts), None, []

        def save():
            if current in wanted:
                self.transcripts[current].sequence = "".join(chunks).upper()

        with self.fasta.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    save()
                    current, chunks = line[1:].split()[0], []
                elif current in wanted:
                    chunks.append(line.strip())
            save()

    def annotate_rows(self, rows: list[dict]) -> list[dict]:
        """Exact-search candidate target sequences in every transcript of the gene."""
        if not rows:
            return []
        gene = base(str(rows[0]["gene_id_base"]))
        targets = {str(r.get("target_seq") or "").upper() for r in rows} - {""}
        by_length = {n: {s for s in targets if len(s) == n}
                     for n in {len(s) for s in targets}}
        sites = defaultdict(list)
        for tx_id in self.by_gene.get(gene, []):
            tx, seq = self.transcripts[tx_id], self.transcripts[tx_id].sequence
            for length, wanted in by_length.items():
                for start in range(len(seq) - length + 1):
                    target = seq[start:start + length]
                    if target in wanted:
                        intervals = tx.genomic_intervals(start, length)
                        if intervals:
                            sites[target].append((tx_id, start, intervals))
        self.candidate_rows += len(rows)
        located = []
        for row in rows:
            found = sites.get(str(row.get("target_seq") or "").upper(), [])
            if not found:
                row["sequence_location_status"] = "not_found_in_gene_transcripts"
                self.candidate_rows_unlocated += 1
                continue
            row["sequence_location_status"] = "located"
            row["_actual_sites"] = found
            row["_actual_intervals"] = tuple(sorted({iv for _, _, ivs in found for iv in ivs}))
            # Preserve the production selector's constitutive-first objective, but
            # base its isoform-coverage tier on exact sequence occurrences rather
            # than the upstream tx_id_set annotation.
            row["n_isoforms_targeted"] = len({tx_id for tx_id, _, _ in found})
            row["position"] = min(pos for _, pos, _ in found) + 1
            self.candidate_rows_located += 1
            located.append(row)
        return located


def intervals_overlap(a, b):
    return a[0] == b[0] and max(a[1], b[1]) <= min(a[2], b[2])


def sequence_far_enough(row: dict, selected: list[dict], _min_spacing: int) -> bool:
    intervals = row.get("_actual_intervals", ())
    if not intervals:
        return False
    for prior in selected:
        old = prior.get("_actual_intervals", ())
        if not old or any(intervals_overlap(a, b) for a in intervals for b in old):
            return False
    return True


def install(gs, reference: SequenceReference):
    original = gs._select_for_gene
    count = 0

    def select(rows, n):
        nonlocal count
        result = original(reference.annotate_rows(rows), n)
        count += 1
        if count % 100 == 0:
            print(f"[sequence-overlap] {count:,} genes; located "
                  f"{reference.candidate_rows_located:,}/{reference.candidate_rows:,} candidates",
                  flush=True)
        return result

    gs._far_enough = sequence_far_enough
    gs._select_for_gene = select


def overlap_pairs(rows: list[dict], reference: SequenceReference):
    from itertools import combinations
    result = []
    for first, second in combinations(reference.annotate_rows(rows), 2):
        hits = [(a, b) for a in first.get("_actual_intervals", ())
                for b in second.get("_actual_intervals", ()) if intervals_overlap(a, b)]
        if hits:
            result.append({
                "gene_id_mouse": first.get("gene_id_base") or first.get("gene_id_mouse"),
                "gene_symbol_mouse": first.get("symbol") or first.get("gene_symbol_mouse"),
                "guide_id_1": first.get("guide_id", ""),
                "guide_id_2": second.get("guide_id", ""),
                "shared_genomic_nt": sum(min(a[2], b[2]) - max(a[1], b[1]) + 1
                                          for a, b in hits),
                "overlap_intervals": repr(hits),
            })
    return result
