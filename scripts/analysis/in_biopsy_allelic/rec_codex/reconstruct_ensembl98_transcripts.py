"""Reconstruct all frozen GTF transcripts; qualify against native FASTA metadata.

Reference-only. No expression, sequencing acquisition, quantification or index fit.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
GENOME = REC / "ensembl98_primary_genome_21998416"
GENOME_SUMMARY_SHA = "9f93a8291cf676bbf1bd1bc552d2097b3e0f9850a05f86431e1d93df65392f6b"
GENOME_ACQUISITION_SHA = "72c4190cc3fbd521b9ae81326f50c7f3a517f2251a6d0844c9704ebad74c17a7"
GENOME_SHA = {"compressed": "197607411dc63c2a5638ab5f24ea23e160fabf46ceaf482a1b01daaf0ad1045a",
              "fasta": "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c",
              "faidx": "57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad"}
GTF = REC / "ensembl98_reference_21998355/Homo_sapiens.GRCh38.98.gtf.gz"
META = REC / "ensembl98_reference_census_21998377/fasta_all_record_metadata.tsv"
AXIS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv"
GUARDS = {GTF: "7966e412a6c5864f232ebcb40d67f9ab873064af090dac28a8cab013aaea18b7",
          META: "a87ab0997a4bec813e474e1e34b198c14b5dceae017592af910181722826af13",
          AXIS: "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd"}
MAX_SEQUENCE = 10_000_000
MAX_FASTA = 1_000_000_000
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def normalize(identifier):
    return re.sub(r"\.\d+$", "", re.sub(r"_PAR_Y$", "", identifier.strip().upper()))


def versioned(identifier, version):
    require(version.isdigit(), "Absent/non-numeric GTF version: " + identifier)
    suffix = re.search(r"\.(\d+)$", identifier)
    require(not suffix or suffix[1] == version, "Conflicting embedded GTF version")
    return identifier if suffix else identifier+"."+version


def parse_annotation():
    transcripts, exons, genes, errors, headers = {}, {}, {}, [], []
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            require(len(line) <= 1_000_000, "GTF line cap exceeded")
            if line.startswith("#"):
                headers.append(line.rstrip()); continue
            cols = line.rstrip().split("\t")
            require(len(cols) == 9, "Malformed GTF")
            if cols[2] not in ("gene", "transcript", "exon"):
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', cols[8]))
            gene = attrs["gene_id"]
            geom = (cols[0], int(cols[3]), int(cols[4]), cols[6])
            if cols[2] == "gene":
                require(gene not in genes, "Duplicate GTF gene ID")
                genes[gene] = geom
            elif cols[2] == "transcript":
                tx = attrs["transcript_id"]
                require(tx not in transcripts, "Duplicate native GTF transcript ID")
                transcripts[tx] = dict(raw_transcript_id=tx, transcript_id=versioned(tx, attrs.get("transcript_version", "")),
                                       raw_gene_id=gene, gene_id=versioned(gene, attrs.get("gene_version", "")),
                                       seqname=geom[0], start_1based=geom[1], end_inclusive=geom[2], strand=geom[3],
                                       gene_biotype=attrs.get("gene_biotype", ""), transcript_biotype=attrs.get("transcript_biotype", ""))
            else:
                tx = attrs["transcript_id"]
                exons.setdefault(tx, []).append((int(attrs.get("exon_number", "0")), geom,
                                                  gene, attrs.get("transcript_version", ""), attrs.get("gene_version", ""),
                                                  attrs.get("exon_id", ""), attrs.get("exon_version", "")))
    require(len(transcripts) == 227368, "Frozen complete GTF transcript count changed")
    for tx in sorted(set(exons)-set(transcripts)):
        errors.append(dict(transcript_id=tx, category="exon_without_transcript", detail=""))
    return transcripts, exons, genes, errors, headers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Reference reconstruction only on compute")
    require(sha(GENOME/"summary.json") == GENOME_SUMMARY_SHA, "Root-reviewed genome summary changed")
    require(sha(GENOME/"acquisition_receipt.json") == GENOME_ACQUISITION_SHA, "Root-reviewed acquisition receipt changed")
    genome_receipt = json.loads((GENOME/"summary.json").read_text())
    require(genome_receipt["coordinate_mismatches"] == 0, "Genome coordinate qualification not complete")
    require(genome_receipt["acquisition"]["source_receipt_sha256"] ==
            "1918a50008667cfb58f25db901b4e3c1ca62aadcb065e6c279c5f0cc81224927", "Genome source differs")
    require(json.loads((GENOME/"acquisition_receipt.json").read_text()) == genome_receipt["acquisition"], "Genome receipts disagree")
    compressed = GENOME/"Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    require(genome_receipt["acquisition"]["sha256"] == GENOME_SHA["compressed"] and
            genome_receipt["fasta_sha256"] == GENOME_SHA["fasta"] and genome_receipt["faidx_sha256"] == GENOME_SHA["faidx"],
            "Frozen genome hashes disagree with receipt")
    require(compressed.stat().st_size == 881211416 and sha(compressed) == GENOME_SHA["compressed"], "Verified compressed genome changed")
    fasta = GENOME/"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    require(fasta.stat().st_size == 3151425857 and sha(fasta) == GENOME_SHA["fasta"] and
            sha(Path(str(fasta)+".fai")) == GENOME_SHA["faidx"], "Verified genome/faidx changed")
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen reference metadata changed: " + str(path))
    require(importlib.metadata.version("pysam") == "0.23.3", "Installed pysam version changed")
    import pysam
    transcripts, exons, genes, errors, headers = parse_annotation()
    by_version = {r["transcript_id"]:r for r in transcripts.values()}
    require(len(by_version) == len(transcripts), "Duplicate versioned transcript IDs")
    existing, absent_source_ids = {}, []
    with META.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tx = row["transcript_id"]
            if tx in by_version:
                require(tx not in existing, "Duplicate matching native FASTA transcript ID")
                existing[tx] = (int(row["length"]), row["sequence_sha256"], row["fasta_gene_id"], row["sources"])
            else:
                absent_source_ids.append(dict(transcript_id=tx, fasta_gene_id=row["fasta_gene_id"], source=row["sources"],
                                              same_native_base_in_gtf=row["normalized_transcript_id"] in transcripts))
    args.output.mkdir(parents=True, exist_ok=False)
    out_fasta = args.output/"ensembl98_gtf_complete_transcripts.fa.partial"
    axis = args.output/"transcript_gene_axis.tsv"
    fields = list(next(iter(transcripts.values())))+["exon_count", "exon_geometry_json", "length", "sequence_sha256", "comparison"]
    statuses, represented_genes, written_bytes = Counter(), set(), 0
    with pysam.FastaFile(str(fasta)) as genome, out_fasta.open("xb") as fa, axis.open("x", newline="") as ax:
        writer = csv.DictWriter(ax, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for tx, row in sorted(transcripts.items()):
            coords = sorted(exons.get(tx, []))
            categories = []
            if not coords:
                categories.append("no_exons")
            if [x[0] for x in coords] != list(range(1, len(coords)+1)):
                categories.append("exon_number_missing_or_duplicate")
            if row["raw_gene_id"] not in genes:
                categories.append("gene_absent")
            elif not (genes[row["raw_gene_id"]][0] == row["seqname"] and
                      genes[row["raw_gene_id"]][3] == row["strand"] and
                      genes[row["raw_gene_id"]][1] <= row["start_1based"] <= row["end_inclusive"] <= genes[row["raw_gene_id"]][2]):
                categories.append("transcript_outside_gene_geometry")
            if row["strand"] not in ("+", "-") or row["seqname"] not in genome.references:
                categories.append("invalid_strand_or_missing_contig")
            exon_id_geometries = {}
            for _, geom, gene, tv, gv, exon_id, exon_version in coords:
                if exon_id in exon_id_geometries and exon_id_geometries[exon_id] != (geom, exon_version):
                    categories.append("within_transcript_exon_id_inconsistency")
                exon_id_geometries[exon_id] = (geom, exon_version)
                if (geom[0], geom[3], gene, versioned(tx, tv), versioned(gene, gv)) != (
                        row["seqname"], row["strand"], row["raw_gene_id"], row["transcript_id"], row["gene_id"]):
                    categories.append("exon_gene_strand_contig_or_version_conflict")
                if not (row["start_1based"] <= geom[1] <= geom[2] <= row["end_inclusive"]):
                    categories.append("exon_outside_transcript")
            if coords and (min(x[1][1] for x in coords), max(x[1][2] for x in coords)) != (row["start_1based"], row["end_inclusive"]):
                categories.append("exon_span_differs_from_transcript")
            for left, right in zip(coords, coords[1:]):
                if not (left[1][2] < right[1][1] if row["strand"] == "+" else left[1][1] > right[1][2]):
                    categories.append("exon_overlap_or_order_error")
            if categories:
                for category in sorted(set(categories)):
                    errors.append(dict(transcript_id=row["transcript_id"], category=category, detail=""))
                statuses["geometry_rejected"] += 1; continue
            length = sum(x[1][2]-x[1][1]+1 for x in coords)
            require(0 < length <= MAX_SEQUENCE, "Transcript sequence length cap exceeded")
            fragments = []
            for _, geom, _, _, _, _, _ in coords:
                require(geom[2] <= genome.get_reference_length(geom[0]), "Exon beyond genome contig")
                fragment = genome.fetch(geom[0], geom[1]-1, geom[2]).upper()
                require(len(fragment) == geom[2]-geom[1]+1 and re.fullmatch(r"[ACGTRYSWKMBDHVN]+", fragment), "Invalid fetched reference exon")
                fragments.append(fragment if row["strand"] == "+" else fragment.translate(COMPLEMENT)[::-1])
            sequence = "".join(fragments)
            require(len(sequence) == length, "Spliced sequence length differs")
            seq_sha = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            comparison = "no_existing_raw_version_id"
            if row["transcript_id"] in existing:
                old_length, old_sha, old_gene, _ = existing[row["transcript_id"]]
                failures = [label for ok, label in ((old_length == length, "length_mismatch"),
                            (old_sha == seq_sha, "sequence_hash_mismatch"), (old_gene == row["gene_id"], "native_gene_version_mismatch")) if not ok]
                comparison = ";".join(failures) if failures else "exact_existing_identity"
                for category in failures:
                    errors.append(dict(transcript_id=row["transcript_id"], category=category,
                                       detail=json.dumps(dict(existing_length=old_length, reconstructed_length=length,
                                                              existing_sha256=old_sha, reconstructed_sha256=seq_sha,
                                                              existing_gene_id=old_gene, reconstructed_gene_id=row["gene_id"]))))
            statuses[comparison] += 1
            encoded = (">"+row["transcript_id"]+" gene:"+row["gene_id"]+"\n"+
                       "\n".join(sequence[i:i+60] for i in range(0, length, 60))+"\n").encode("ascii")
            require(written_bytes+len(encoded) <= MAX_FASTA, "Transcript FASTA1GB cap exceeded; partial preserved")
            fa.write(encoded); written_bytes += len(encoded)
            writer.writerow(dict(**row, exon_count=len(coords), exon_geometry_json=json.dumps(coords),
                                 length=length, sequence_sha256=seq_sha, comparison=comparison))
            represented_genes.add(normalize(row["raw_gene_id"]))
    with AXIS.open(newline="") as handle:
        modeled = list(csv.DictReader(handle, delimiter="\t"))
    require(len(modeled) == 42163, "Frozen modeled gene count changed")
    missing = [r for r in modeled if normalize(r["stable_gene_id"]) not in represented_genes]
    for name, data, fieldnames in (("reconstruction_diagnostics.tsv", errors, ["transcript_id", "category", "detail"]),
            ("source_fasta_ids_outside_gtf_reference.tsv", absent_source_ids, ["transcript_id", "fasta_gene_id", "source", "same_native_base_in_gtf"]),
            ("missing_modeled_genes.tsv", missing, list(modeled[0]))):
        with (args.output/name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t"); writer.writeheader(); writer.writerows(data)
    qualified = not errors and not missing and sum(statuses.values()) == 227368
    final_name = "ensembl98_gtf_complete_transcripts.fa" if qualified else "ensembl98_gtf_complete_transcripts.unqualified.fa"
    out_fasta.rename(args.output/final_name)
    summary = dict(gtf_transcripts=len(transcripts), reconstruction_statuses=dict(statuses),
                   diagnostic_categories=dict(Counter(r["category"] for r in errors)), missing_modeled_genes=len(missing),
                   represented_modeled_genes=42163-len(missing), source_fasta_matching_raw_version_ids=len(existing),
                   source_fasta_ids_outside_gtf=len(absent_source_ids), fasta_bytes=written_bytes,
                   reference_identity_qualified=qualified, quantification_or_indexing_authorized=False,
                   private_pisces_reconstruction_asserted=False, receiving_values_read=False, gene_zero_fill=False,
                   gtf_headers=headers, genome_summary_sha256=GENOME_SUMMARY_SHA,
                   genome_acquisition_receipt_sha256=GENOME_ACQUISITION_SHA,
                   genome_asset_sha256=dict(compressed=genome_receipt["acquisition"]["sha256"],
                                             fasta=genome_receipt["fasta_sha256"], faidx=genome_receipt["faidx_sha256"]),
                   input_sha256={str(p):v for p,v in GUARDS.items()}, script_sha256=sha(Path(__file__)),
                   launcher_sha256=sha(Path(__file__).with_name("run_reconstruct_ensembl98_transcripts.sbatch")),
                   environment=dict(python=sys.version, pysam=pysam.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"]),
                   output_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:summary[k] for k in ("gtf_transcripts", "reconstruction_statuses", "diagnostic_categories", "represented_modeled_genes", "reference_identity_qualified")}, sort_keys=True))
    if not qualified:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
