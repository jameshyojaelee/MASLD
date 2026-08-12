#!/usr/bin/env python3
"""Build stable-ID-first GENCODE gene and lncRNA relationship tables.

Coordinates are interpreted as one-based, closed intervals, matching GTF.
The implementation uses only the Python standard library so the production
environment does not acquire a new dependency.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
from typing import Iterable, Iterator, Mapping, Sequence


GENCODE_V49_GTF_SHA256 = (
    "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
)
GENCODE_V49_GTF_BYTES = 97_570_076
ANNOTATION_RELEASE = "GENCODE v49"
CANONICAL_CHROMOSOMES = frozenset(
    [f"chr{number}" for number in range(1, 23)] + ["chrX", "chrY"]
)
ENSEMBL_GENE_ID = re.compile(r"^(ENSG[0-9]+)(?:\.([0-9]+))?$")

GENE_IDENTITY_FIELDS = (
    "annotation_release",
    "gene_id_versioned",
    "gene_id_base",
    "gene_version",
    "gene_name",
    "gene_type",
    "chromosome",
    "start_1based",
    "end_1based",
    "strand",
    "source",
    "level",
    "is_canonical_chromosome",
    "n_versions_for_base_id",
    "base_id_mapping_status",
    "n_gene_ids_for_symbol",
    "symbol_mapping_status",
)

LNCRNA_CLASS_FIELDS = (
    "annotation_release",
    "gene_id_versioned",
    "gene_id_base",
    "gene_name",
    "chromosome",
    "start_1based",
    "end_1based",
    "strand",
    "is_canonical_chromosome",
    "lncrna_genomic_class",
    "mapping_status",
    "main_text_eligible",
    "main_text_exclusion_reason",
    "n_antisense_exonic_pc_genes",
    "n_same_strand_exonic_pc_genes",
    "n_antisense_gene_body_pc_genes",
    "n_same_strand_gene_body_pc_genes",
    "n_opposite_strand_promoter_pc_genes",
    "classification_partner_gene_ids",
    "classification_partner_gene_names",
    "all_related_pc_gene_ids",
    "relation_evidence",
)

LNCRNA_EXCLUSION_FIELDS = (
    "annotation_release",
    "gene_id_versioned",
    "gene_id_base",
    "gene_name",
    "chromosome",
    "lncrna_genomic_class",
    "mapping_status",
    "exclusion_reason",
)

CLASS_ANTISENSE_EXONIC = "antisense_exonic_overlap"
CLASS_SAME_EXONIC = "same_strand_exonic_overlap"
CLASS_ANTISENSE_BODY = "antisense_gene_body_overlap"
CLASS_SAME_BODY = "same_strand_gene_body_overlap"
CLASS_OPPOSITE_PROMOTER = "opposite_strand_promoter_proximal"
CLASS_INTERGENIC = "intergenic"
CLASS_COMPLEX = "complex"

APPROVED_CLASSES = (
    CLASS_ANTISENSE_EXONIC,
    CLASS_SAME_EXONIC,
    CLASS_ANTISENSE_BODY,
    CLASS_SAME_BODY,
    CLASS_OPPOSITE_PROMOTER,
    CLASS_INTERGENIC,
    CLASS_COMPLEX,
)


class IdentityBuildError(RuntimeError):
    """Raised when a source or identity invariant is violated."""


@dataclass(frozen=True)
class Gene:
    """One GENCODE gene feature."""

    gene_id: str
    gene_id_base: str
    gene_version: str
    gene_name: str
    gene_type: str
    chromosome: str
    start: int
    end: int
    strand: str
    source: str
    level: str


@dataclass(frozen=True)
class RelationEvidence:
    """Protein-coding relationships used to classify one lncRNA."""

    antisense_exonic: frozenset[str] = frozenset()
    same_strand_exonic: frozenset[str] = frozenset()
    antisense_body: frozenset[str] = frozenset()
    same_strand_body: frozenset[str] = frozenset()
    opposite_promoter: frozenset[str] = frozenset()


class BinnedIntervalIndex:
    """Small dependency-free interval index for one-based closed intervals."""

    def __init__(self, bin_size: int = 100_000) -> None:
        if bin_size <= 0:
            raise ValueError("bin_size must be positive")
        self.bin_size = bin_size
        self._bins: dict[tuple[str, int], list[tuple[int, int, str]]] = defaultdict(
            list
        )

    def add(self, chromosome: str, start: int, end: int, value: str) -> None:
        if start < 1 or end < start:
            raise ValueError(f"invalid closed interval: {chromosome}:{start}-{end}")
        first = (start - 1) // self.bin_size
        last = (end - 1) // self.bin_size
        for bin_number in range(first, last + 1):
            self._bins[(chromosome, bin_number)].append((start, end, value))

    def query(self, chromosome: str, start: int, end: int) -> set[str]:
        if start < 1 or end < start:
            raise ValueError(f"invalid closed interval: {chromosome}:{start}-{end}")
        first = (start - 1) // self.bin_size
        last = (end - 1) // self.bin_size
        hits: set[str] = set()
        for bin_number in range(first, last + 1):
            for candidate_start, candidate_end, value in self._bins.get(
                (chromosome, bin_number), ()
            ):
                if candidate_start <= end and candidate_end >= start:
                    hits.add(value)
        return hits


def require(condition: bool, message: str) -> None:
    if not condition:
        raise IdentityBuildError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(value: str) -> dict[str, str]:
    """Parse the semicolon-delimited attribute column of a GTF record."""

    attributes: dict[str, str] = {}
    for item in value.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(maxsplit=1)
        require(len(parts) == 2, f"malformed GTF attribute: {item!r}")
        key, raw = parts
        attributes[key] = raw.strip().strip('"')
    return attributes


def iter_gtf(
    path: Path,
) -> Iterator[tuple[str, str, str, int, int, str, dict[str, str]]]:
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == 9, f"GTF line {line_number} does not have 9 fields")
            chromosome, source, feature, start, end, _, strand, _, raw_attributes = (
                fields
            )
            require(
                strand in {"+", "-", "."}, f"invalid strand on GTF line {line_number}"
            )
            yield (
                chromosome,
                source,
                feature,
                int(start),
                int(end),
                strand,
                parse_attributes(raw_attributes),
            )


def split_ensembl_gene_id(gene_id: str) -> tuple[str, str]:
    match = ENSEMBL_GENE_ID.fullmatch(gene_id)
    require(match is not None, f"non-Ensembl or malformed gene_id: {gene_id!r}")
    return match.group(1), match.group(2) or ""


def load_genes(path: Path) -> dict[str, Gene]:
    genes: dict[str, Gene] = {}
    for chromosome, source, feature, start, end, strand, attributes in iter_gtf(path):
        if feature != "gene":
            continue
        for required in ("gene_id", "gene_type"):
            require(required in attributes, f"gene feature lacks {required}")
        gene_id = attributes["gene_id"]
        require(gene_id not in genes, f"duplicate versioned gene_id: {gene_id}")
        gene_id_base, gene_version = split_ensembl_gene_id(gene_id)
        require(strand in {"+", "-"}, f"gene {gene_id} has unstranded annotation")
        genes[gene_id] = Gene(
            gene_id=gene_id,
            gene_id_base=gene_id_base,
            gene_version=gene_version,
            gene_name=attributes.get("gene_name", ""),
            gene_type=attributes["gene_type"],
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand,
            source=source,
            level=attributes.get("level", ""),
        )
    require(bool(genes), f"no gene features found in {path}")
    return genes


def chromosome_sort_key(chromosome: str) -> tuple[int, int | str]:
    if chromosome.startswith("chr"):
        suffix = chromosome[3:]
        if suffix.isdigit() and 1 <= int(suffix) <= 22:
            return (0, int(suffix))
        if suffix == "X":
            return (0, 23)
        if suffix == "Y":
            return (0, 24)
    return (1, chromosome)


def gene_sort_key(gene: Gene) -> tuple[tuple[int, int | str], int, int, str]:
    return (chromosome_sort_key(gene.chromosome), gene.start, gene.end, gene.gene_id)


def build_interval_inputs(
    gtf_path: Path,
    genes: Mapping[str, Gene],
) -> tuple[BinnedIntervalIndex, BinnedIntervalIndex, dict[str, list[tuple[int, int]]]]:
    pc_body_index = BinnedIntervalIndex()
    pc_exon_index = BinnedIntervalIndex()
    lncrna_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)

    for gene in genes.values():
        if gene.gene_type == "protein_coding":
            pc_body_index.add(gene.chromosome, gene.start, gene.end, gene.gene_id)

    for chromosome, _, feature, start, end, _, attributes in iter_gtf(gtf_path):
        if feature != "exon":
            continue
        gene_id = attributes.get("gene_id", "")
        require(gene_id in genes, f"exon references unknown gene_id: {gene_id!r}")
        gene = genes[gene_id]
        require(
            chromosome == gene.chromosome,
            f"exon chromosome differs from gene for {gene_id}",
        )
        if gene.gene_type == "protein_coding":
            pc_exon_index.add(chromosome, start, end, gene_id)
        elif gene.gene_type == "lncRNA":
            lncrna_exons[gene_id].append((start, end))

    return pc_body_index, pc_exon_index, lncrna_exons


def tss(gene: Gene) -> int:
    return gene.start if gene.strand == "+" else gene.end


def partition_by_strand(
    gene: Gene,
    candidate_ids: Iterable[str],
    genes: Mapping[str, Gene],
) -> tuple[set[str], set[str]]:
    opposite: set[str] = set()
    same: set[str] = set()
    for candidate_id in candidate_ids:
        candidate = genes[candidate_id]
        if candidate.strand == gene.strand:
            same.add(candidate_id)
        else:
            opposite.add(candidate_id)
    return opposite, same


def collect_relation_evidence(
    gene: Gene,
    exons: Sequence[tuple[int, int]],
    genes: Mapping[str, Gene],
    pc_body_index: BinnedIntervalIndex,
    pc_exon_index: BinnedIntervalIndex,
    promoter_window: int = 2_000,
) -> RelationEvidence:
    exon_hits: set[str] = set()
    for start, end in exons:
        exon_hits.update(pc_exon_index.query(gene.chromosome, start, end))
    opposite_exonic, same_exonic = partition_by_strand(gene, exon_hits, genes)

    body_hits = pc_body_index.query(gene.chromosome, gene.start, gene.end) - exon_hits
    opposite_body, same_body = partition_by_strand(gene, body_hits, genes)

    query_start = max(1, tss(gene) - promoter_window)
    query_end = tss(gene) + promoter_window
    promoter_hits = pc_body_index.query(gene.chromosome, query_start, query_end)
    opposite_promoter = {
        candidate_id
        for candidate_id in promoter_hits
        if genes[candidate_id].strand != gene.strand
        and abs(tss(genes[candidate_id]) - tss(gene)) <= promoter_window
    }

    return RelationEvidence(
        antisense_exonic=frozenset(opposite_exonic),
        same_strand_exonic=frozenset(same_exonic),
        antisense_body=frozenset(opposite_body),
        same_strand_body=frozenset(same_body),
        opposite_promoter=frozenset(opposite_promoter),
    )


def classify_relation(evidence: RelationEvidence) -> tuple[str, frozenset[str]]:
    """Classify using the approved priority and a conflict-at-tier rule.

    Opposite- and same-strand evidence at the same highest occupied overlap
    tier is incompatible and therefore classified as complex. Evidence at a
    lower tier never overrides a higher-tier relation.
    """

    if evidence.antisense_exonic and evidence.same_strand_exonic:
        return CLASS_COMPLEX, evidence.antisense_exonic | evidence.same_strand_exonic
    if evidence.antisense_exonic:
        return CLASS_ANTISENSE_EXONIC, evidence.antisense_exonic
    if evidence.same_strand_exonic:
        return CLASS_SAME_EXONIC, evidence.same_strand_exonic
    if evidence.antisense_body and evidence.same_strand_body:
        return CLASS_COMPLEX, evidence.antisense_body | evidence.same_strand_body
    if evidence.antisense_body:
        return CLASS_ANTISENSE_BODY, evidence.antisense_body
    if evidence.same_strand_body:
        return CLASS_SAME_BODY, evidence.same_strand_body
    if evidence.opposite_promoter:
        return CLASS_OPPOSITE_PROMOTER, evidence.opposite_promoter
    return CLASS_INTERGENIC, frozenset()


def relation_evidence_string(evidence: RelationEvidence) -> str:
    values = (
        (CLASS_ANTISENSE_EXONIC, evidence.antisense_exonic),
        (CLASS_SAME_EXONIC, evidence.same_strand_exonic),
        (CLASS_ANTISENSE_BODY, evidence.antisense_body),
        (CLASS_SAME_BODY, evidence.same_strand_body),
        (CLASS_OPPOSITE_PROMOTER, evidence.opposite_promoter),
    )
    return ";".join(f"{label}={len(ids)}" for label, ids in values)


def join_sorted(values: Iterable[str]) -> str:
    return ";".join(sorted(set(values)))


def build_gene_identity_rows(
    genes: Mapping[str, Gene],
    annotation_release: str = ANNOTATION_RELEASE,
) -> list[dict[str, str]]:
    base_counts = Counter(gene.gene_id_base for gene in genes.values())
    symbol_counts = Counter(gene.gene_name for gene in genes.values() if gene.gene_name)
    rows: list[dict[str, str]] = []
    for gene in sorted(genes.values(), key=gene_sort_key):
        n_for_symbol = symbol_counts.get(gene.gene_name, 0)
        rows.append(
            {
                "annotation_release": annotation_release,
                "gene_id_versioned": gene.gene_id,
                "gene_id_base": gene.gene_id_base,
                "gene_version": gene.gene_version,
                "gene_name": gene.gene_name,
                "gene_type": gene.gene_type,
                "chromosome": gene.chromosome,
                "start_1based": str(gene.start),
                "end_1based": str(gene.end),
                "strand": gene.strand,
                "source": gene.source,
                "level": gene.level,
                "is_canonical_chromosome": str(
                    gene.chromosome in CANONICAL_CHROMOSOMES
                ).lower(),
                "n_versions_for_base_id": str(base_counts[gene.gene_id_base]),
                "base_id_mapping_status": (
                    "unique" if base_counts[gene.gene_id_base] == 1 else "duplicated"
                ),
                "n_gene_ids_for_symbol": str(n_for_symbol),
                "symbol_mapping_status": (
                    "missing"
                    if not gene.gene_name
                    else "unique"
                    if n_for_symbol == 1
                    else "duplicated"
                ),
            }
        )
    return rows


def build_lncrna_rows(
    genes: Mapping[str, Gene],
    pc_body_index: BinnedIntervalIndex,
    pc_exon_index: BinnedIntervalIndex,
    lncrna_exons: Mapping[str, Sequence[tuple[int, int]]],
    annotation_release: str = ANNOTATION_RELEASE,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    exclusions: list[dict[str, str]] = []
    lncrnas = sorted(
        (gene for gene in genes.values() if gene.gene_type == "lncRNA"),
        key=gene_sort_key,
    )
    for gene in lncrnas:
        exons = lncrna_exons.get(gene.gene_id, ())
        require(bool(exons), f"lncRNA gene lacks exon features: {gene.gene_id}")
        evidence = collect_relation_evidence(
            gene,
            exons,
            genes,
            pc_body_index,
            pc_exon_index,
        )
        genomic_class, classification_partners = classify_relation(evidence)
        mapping_status = (
            "mapping_ambiguous"
            if genomic_class in {CLASS_SAME_EXONIC, CLASS_COMPLEX}
            else "mapping_unambiguous"
        )
        reasons: list[str] = []
        if gene.chromosome not in CANONICAL_CHROMOSOMES:
            reasons.append("alternative_locus")
        if genomic_class == CLASS_SAME_EXONIC:
            reasons.append("mapping_ambiguous_same_strand_exonic_overlap")
        elif genomic_class == CLASS_COMPLEX:
            reasons.append("mapping_ambiguous_complex")
        exclusion_reason = ";".join(reasons)
        all_related = (
            evidence.antisense_exonic
            | evidence.same_strand_exonic
            | evidence.antisense_body
            | evidence.same_strand_body
            | evidence.opposite_promoter
        )
        partner_names = [
            genes[partner].gene_name for partner in classification_partners
        ]
        row = {
            "annotation_release": annotation_release,
            "gene_id_versioned": gene.gene_id,
            "gene_id_base": gene.gene_id_base,
            "gene_name": gene.gene_name,
            "chromosome": gene.chromosome,
            "start_1based": str(gene.start),
            "end_1based": str(gene.end),
            "strand": gene.strand,
            "is_canonical_chromosome": str(
                gene.chromosome in CANONICAL_CHROMOSOMES
            ).lower(),
            "lncrna_genomic_class": genomic_class,
            "mapping_status": mapping_status,
            "main_text_eligible": str(not reasons).lower(),
            "main_text_exclusion_reason": exclusion_reason,
            "n_antisense_exonic_pc_genes": str(len(evidence.antisense_exonic)),
            "n_same_strand_exonic_pc_genes": str(len(evidence.same_strand_exonic)),
            "n_antisense_gene_body_pc_genes": str(len(evidence.antisense_body)),
            "n_same_strand_gene_body_pc_genes": str(len(evidence.same_strand_body)),
            "n_opposite_strand_promoter_pc_genes": str(len(evidence.opposite_promoter)),
            "classification_partner_gene_ids": join_sorted(classification_partners),
            "classification_partner_gene_names": join_sorted(partner_names),
            "all_related_pc_gene_ids": join_sorted(all_related),
            "relation_evidence": relation_evidence_string(evidence),
        }
        rows.append(row)
        if reasons:
            exclusions.append(
                {
                    "annotation_release": annotation_release,
                    "gene_id_versioned": gene.gene_id,
                    "gene_id_base": gene.gene_id_base,
                    "gene_name": gene.gene_name,
                    "chromosome": gene.chromosome,
                    "lncrna_genomic_class": genomic_class,
                    "mapping_status": mapping_status,
                    "exclusion_reason": exclusion_reason,
                }
            )
    return rows, exclusions


def write_tsv(
    path: Path, rows: Sequence[Mapping[str, str]], fields: Sequence[str]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            require(set(row) == set(fields), f"schema drift while writing {path.name}")
            writer.writerow(row)


def verify_source(
    gtf_path: Path,
    *,
    expected_sha256: str,
    expected_bytes: int,
) -> tuple[str, int]:
    require(gtf_path.is_file(), f"GTF does not exist: {gtf_path}")
    observed_bytes = gtf_path.stat().st_size
    require(
        observed_bytes == expected_bytes,
        f"GTF byte count mismatch: expected {expected_bytes}, observed {observed_bytes}",
    )
    observed_sha256 = sha256_file(gtf_path)
    require(
        observed_sha256 == expected_sha256,
        f"GTF SHA256 mismatch: expected {expected_sha256}, observed {observed_sha256}",
    )
    return observed_sha256, observed_bytes


def build_release(
    gtf_path: Path,
    output_dir: Path,
    *,
    annotation_release: str = ANNOTATION_RELEASE,
    fixture_mode: bool = False,
) -> dict[str, int | str]:
    """Build an immutable candidate directory and return its census."""

    require(not output_dir.exists(), f"output directory already exists: {output_dir}")
    if fixture_mode:
        require(
            os.environ.get("MASLD_NONCODING_FIXTURE_MODE") == "1",
            "fixture mode requires MASLD_NONCODING_FIXTURE_MODE=1",
        )
        observed_bytes = gtf_path.stat().st_size
        observed_sha256 = sha256_file(gtf_path)
    else:
        require(
            annotation_release == ANNOTATION_RELEASE,
            f"production annotation release must be exactly {ANNOTATION_RELEASE!r}",
        )
        observed_sha256, observed_bytes = verify_source(
            gtf_path,
            expected_sha256=GENCODE_V49_GTF_SHA256,
            expected_bytes=GENCODE_V49_GTF_BYTES,
        )

    genes = load_genes(gtf_path)
    pc_body_index, pc_exon_index, lncrna_exons = build_interval_inputs(gtf_path, genes)
    identity_rows = build_gene_identity_rows(genes, annotation_release)
    lncrna_rows, exclusion_rows = build_lncrna_rows(
        genes,
        pc_body_index,
        pc_exon_index,
        lncrna_exons,
        annotation_release,
    )

    output_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(output_dir / "gene_identity.tsv", identity_rows, GENE_IDENTITY_FIELDS)
    write_tsv(
        output_dir / "lncrna_genomic_class.tsv",
        lncrna_rows,
        LNCRNA_CLASS_FIELDS,
    )
    write_tsv(
        output_dir / "lncrna_mapping_exclusions.tsv",
        exclusion_rows,
        LNCRNA_EXCLUSION_FIELDS,
    )

    source_rows = [
        {
            "source_role": "gencode_annotation",
            "annotation_release": annotation_release,
            "source_path": str(gtf_path.resolve()),
            "size_bytes": str(observed_bytes),
            "sha256": observed_sha256,
            "fixture_mode": str(fixture_mode).lower(),
        }
    ]
    write_tsv(
        output_dir / "source_manifest.tsv",
        source_rows,
        (
            "source_role",
            "annotation_release",
            "source_path",
            "size_bytes",
            "sha256",
            "fixture_mode",
        ),
    )

    class_counts = Counter(row["lncrna_genomic_class"] for row in lncrna_rows)
    execution_manifest = {
        "annotation_release": annotation_release,
        "coordinate_system": "GTF one-based closed",
        "fixture_mode": fixture_mode,
        "python_version": platform.python_version(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256_file(Path(__file__).resolve()),
        "n_genes": len(identity_rows),
        "n_lncrna_genes": len(lncrna_rows),
        "n_mapping_exclusions": len(exclusion_rows),
        "lncrna_class_counts": dict(sorted(class_counts.items())),
    }
    execution_path = output_dir / "execution_manifest.json"
    with execution_path.open("x", encoding="utf-8") as handle:
        json.dump(execution_manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    checksum_rows = []
    for path in sorted(output_dir.iterdir(), key=lambda item: item.name):
        if path.name == "checksum_manifest.tsv":
            continue
        checksum_rows.append(
            {
                "relative_path": path.name,
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        )
    write_tsv(
        output_dir / "checksum_manifest.tsv",
        checksum_rows,
        ("relative_path", "size_bytes", "sha256"),
    )
    return {
        "n_genes": len(identity_rows),
        "n_lncrna_genes": len(lncrna_rows),
        "n_mapping_exclusions": len(exclusion_rows),
        "output_dir": str(output_dir.resolve()),
    }
