"""Acquire exactly three public Ensembl98 reference assets and census their axes.

No sequencing/count acquisition, index building or quantifier configuration.
Run on compute; reference sequences are streamed with a 10-MB record cap.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from urllib.parse import urlsplit
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
META = REC / "ensembl98_reference_metadata_20260930T215320Z"
AXIS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv"
GUARDS = {"receipt.json": "019ebbe68ea09f9ba675d1ce03cb5415b25d0d7853950bdb6d29d4467426b614",
          "detail_receipt.json": "3158098e520b70f573f092791317975c3a1724e02fc2415274f2bcfe736b9d49"}
AXIS_SHA = "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd"
REUSE = REC / "ensembl98_reference_21998355"
REUSE_RECEIPT_SHA = "946f870eec333b818250fb5736abc32c5ff14bcf3e1dcde41b3faa1c2944cf8f"
REUSE_ASSET_SHA = {
    "cdna": "6082e92b92fd5ea5e0c31216496904e8df5c93955c25cb08c47ab815463c1e1f",
    "ncrna": "4feb9291eae60d1c9185cf2dcbace34497e4b89bbe261e11cf6512c9f7387450",
    "gtf": "7966e412a6c5864f232ebcb40d67f9ab873064af090dac28a8cab013aaea18b7"}
MAX_FILE = 128_000_000
MAX_TOTAL = 500_000_000
MAX_SEQUENCE = 10_000_000
MAX_LINE = 1_000_000
MAX_TRANSCRIPTS = 1_000_000
MAX_GENES = 500_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def normalize(identifier):
    return re.sub(r"\.\d+$", "", re.sub(r"_PAR_Y$", "", identifier.strip().upper()))


def write_tsv(path, rows, fields):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def source_plan():
    for name, expected in GUARDS.items():
        require(sha(META / name) == expected, "Frozen metadata receipt changed")
    listings = json.loads((META / "receipt.json").read_text())
    details = json.loads((META / "detail_receipt.json").read_text())
    for row in listings:
        require(sha(META / (row["kind"] + "_listing.html")) == row["sha256"], "Archived listing changed")
    for row in details:
        if not row["name"].endswith("_HEAD"):
            require(sha(META / (row["name"] + ".txt")) == row["sha256"], "Archived source receipt changed")
    plan = []
    for kind in ("cdna", "ncrna", "gtf"):
        heads = [r for r in details if r["name"] == kind + "_HEAD"]
        require(len(heads) == 1, "Missing/duplicate reference HEAD")
        row = heads[0]
        parsed = urlsplit(row["url"])
        require(parsed.scheme == "https" and parsed.hostname == "ftp.ensembl.org" and
                parsed.path.startswith("/pub/release-98/"), "Unexpected reference URL")
        require(row["status"] == 200, "Source HEAD status not successful")
        size = int(row["content_length"])
        require(0 < size <= MAX_FILE, "Reference exceeds individual byte cap")
        name = Path(parsed.path).name
        checksum_rows = [line.split() for line in (META / (kind + "_CHECKSUMS.txt")).read_text().splitlines()]
        matching = [r for r in checksum_rows if len(r) == 3 and r[2] == name]
        require(len(matching) == 1, "Missing/ambiguous published BSD checksum")
        plan.append(dict(kind=kind, url=row["url"], filename=name, bytes=size,
                         bsd_sum=int(matching[0][0]), bsd_blocks=int(matching[0][1])))
    require(sum(p["bytes"] for p in plan) <= MAX_TOTAL, "Reference bundle exceeds total cap")
    require(sha(AXIS) == AXIS_SHA, "Modeled RNA metadata axis changed")
    return plan


def acquire(row, output):
    destination = output / row["filename"]
    partial = output / (row["filename"] + ".partial")
    require(not destination.exists() and not partial.exists(), "Refusing reference overwrite")
    h = hashlib.sha256()
    nbytes = 0
    with urlopen(row["url"], timeout=60) as response, partial.open("xb") as handle:
        require(response.geturl() == row["url"], "Unexpected reference redirect")
        require(response.status == 200 and response.headers.get("Content-Length") is not None,
                "Reference response lacks successful fixed length")
        require(int(response.headers["Content-Length"]) == row["bytes"], "Reference length changed since HEAD")
        headers = {k:response.headers.get(k) for k in ("Content-Length", "ETag", "Last-Modified", "Content-Type")}
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            nbytes += len(block)
            require(nbytes <= row["bytes"] and nbytes <= MAX_FILE, "Download byte cap exceeded")
            h.update(block); handle.write(block)
    require(nbytes == row["bytes"], "Truncated reference download")
    checksum = subprocess.run(["sum", "-r", str(partial)], check=True, capture_output=True, text=True)
    tokens = checksum.stdout.split()
    require(len(tokens) >= 2 and (int(tokens[0]), int(tokens[1])) == (row["bsd_sum"], row["bsd_blocks"]),
            "Published BSD Unix sum/1024-byte blocks disagree")
    with partial.open("rb") as handle:
        require(handle.read(2) == b"\x1f\x8b", "Reference is not gzip")
    partial.rename(destination)
    return dict(**row, sha256=h.hexdigest(), response_headers=headers, unix_sum=checksum.stdout.strip())


def gzip_lines(path):
    with gzip.open(path, "rb") as handle:
        while True:
            line = handle.readline(MAX_LINE + 1)
            if not line:
                break
            require(len(line) <= MAX_LINE, "Reference line exceeds cap")
            yield line.decode("utf-8").rstrip("\r\n")


def parse_gtf(path):
    transcripts, genes = {}, {}
    for line in gzip_lines(path):
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        require(len(columns) == 9, "Malformed GTF")
        if columns[2] not in ("gene", "transcript"):
            continue
        attrs = dict(re.findall(r'(\w+) "([^\"]*)"', columns[8]))
        require("gene_id" in attrs, "Missing GTF gene ID")
        # Preserve native PAR/locus IDs in annotation geometry; normalization
        # below is only for comparison to the released modeled-gene axis.
        gene = re.sub(r"\.\d+$", "", attrs["gene_id"].strip().upper())
        base = dict(gene_id=gene, raw_gene_id=attrs["gene_id"], seqname=columns[0], start_1based=int(columns[3]),
                    end_inclusive=int(columns[4]), strand=columns[6])
        require(base["start_1based"] > 0 and base["end_inclusive"] >= base["start_1based"], "Invalid GTF coordinates")
        if columns[2] == "gene":
            row = dict(**base, gene_version=attrs.get("gene_version", ""),
                       gene_name=attrs.get("gene_name", ""), gene_biotype=attrs.get("gene_biotype", ""))
            if row not in genes.setdefault(gene, []):
                genes[gene].append(row)
        else:
            require("transcript_id" in attrs, "Missing GTF transcript ID")
            transcript = re.sub(r"\.\d+$", "", attrs["transcript_id"].strip().upper())
            row = dict(**base, transcript_id=transcript, raw_transcript_id=attrs["transcript_id"],
                       transcript_version=attrs.get("transcript_version", ""), gene_version=attrs.get("gene_version", ""))
            if row not in transcripts.setdefault(transcript, []):
                transcripts[transcript].append(row)
        require(len(transcripts) <= MAX_TRANSCRIPTS and len(genes) <= MAX_GENES, "Annotation metadata cap exceeded")
    require(transcripts and genes, "Empty GTF axis")
    require(all(r["gene_id"] in genes for rows in transcripts.values() for r in rows), "GTF transcript gene absent")
    return transcripts, genes


def fasta_records(path):
    header, h, length = None, None, 0
    for line in gzip_lines(path):
        if line.startswith(">"):
            if header is not None:
                require(length > 0, "Empty reference transcript")
                yield header, length, h.hexdigest()
            header, h, length = line[1:], hashlib.sha256(), 0
        elif line:
            require(header is not None, "FASTA sequence precedes header")
            sequence = line.upper()
            require(re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence), "Unexpected transcript alphabet")
            length += len(sequence)
            require(length <= MAX_SEQUENCE, "Transcript exceeds 10-MB sequence cap")
            h.update(sequence.encode("ascii"))
    require(header is not None and length > 0, "Empty FASTA")
    yield header, length, h.hexdigest()


def census(output, plan, asset_root=None):
    paths = {r["kind"]:(asset_root or output)/r["filename"] for r in plan}
    annotation_headers = []
    for line in gzip_lines(paths["gtf"]):
        if not line.startswith("#"):
            break
        annotation_headers.append(line)
    (output/"gtf_header.txt").write_text("\n".join(annotation_headers)+"\n")
    gtf, genes = parse_gtf(paths["gtf"])
    union, sets, duplicate_rows, records, diagnostics = {}, {}, [], [], []
    ambiguous = set()
    for kind in ("cdna", "ncrna"):
        ids = set()
        for header, length, sequence_sha in fasta_records(paths[kind]):
            transcript = header.split()[0]
            stable = normalize(transcript)
            native_base = re.sub(r"\.\d+$", "", transcript.strip().upper())
            matches = re.findall(r"(?:^|\s)gene:([^\s]+)", header)
            require(len(matches) == 1, "Missing/ambiguous FASTA header gene")
            gene = normalize(matches[0])
            geometry = gtf.get(native_base, [])
            gene_base = re.sub(r"\.\d+$", "", matches[0].strip().upper())
            tx_version = re.search(r"\.(\d+)$", transcript)
            gene_version = re.search(r"\.(\d+)$", matches[0])
            categories = []
            if not geometry:
                categories.append("transcript_absent_from_gtf")
            elif any(r["gene_id"] != gene_base for r in geometry):
                categories.append("native_gene_conflict")
            if geometry and any(normalize(r["gene_id"]) != gene for r in geometry):
                categories.append("normalized_gene_conflict")
            if geometry and tx_version and any(r["transcript_version"] != tx_version[1] for r in geometry):
                categories.append("transcript_version_discrepancy")
            if geometry and gene_version and any(r["gene_version"] != gene_version[1] for r in geometry):
                categories.append("gene_version_discrepancy")
            mapping_confirmed = bool(geometry) and not any("conflict" in c for c in categories)
            row = dict(transcript_id=transcript, normalized_transcript_id=stable,
                       gene_id=gene, fasta_gene_id=matches[0], length=length, sequence_sha256=sequence_sha,
                       sources=kind, fasta_header=header, mapping_confirmed=mapping_confirmed,
                       diagnostic_categories=";".join(categories), gtf_geometry_json=json.dumps(geometry, sort_keys=True))
            records.append(dict(row))
            if categories:
                diagnostics.append(dict(row))
            if transcript in union:
                previous = union[transcript]
                identical = all(previous[k] == row[k] for k in ("gene_id", "fasta_gene_id", "length", "sequence_sha256"))
                if not identical:
                    ambiguous.add(transcript)
                duplicate_rows.append(dict(transcript_id=transcript, first_source=previous["sources"], duplicate_source=kind,
                                           identity_agrees=identical))
                if kind not in previous["sources"].split(";"):
                    previous["sources"] += ";" + kind
            else:
                union[transcript] = row
            ids.add(transcript)
            require(len(union) <= MAX_TRANSCRIPTS, "Transcript metadata cap exceeded")
        sets[kind] = ids
    # Header coverage is descriptive only. Conflicting duplicate IDs do not
    # receive an arbitrary union representative or a modeled-gene admission.
    union = {t:r for t,r in union.items() if t not in ambiguous}
    cdna_genes = {r["gene_id"] for r in records if r["sources"] == "cdna"}
    union_genes = {r["gene_id"] for r in records}
    confirmed_genes = {r["gene_id"] for r in records if r["mapping_confirmed"] and r["transcript_id"] not in ambiguous}
    affected_genes = {r["gene_id"] for r in diagnostics}
    affected_genes.update(normalize(g["gene_id"]) for r in diagnostics
                          for g in json.loads(r["gtf_geometry_json"]))
    normalized_gtf_genes = {normalize(g) for g in genes}
    with AXIS.open(newline="") as handle:
        modeled = list(csv.DictReader(handle, delimiter="\t"))
    require(len(modeled) == 42163, "Modeled axis size changed")
    coverage = []
    for source in modeled:
        target = normalize(source["stable_gene_id"])
        native = normalize(source["matrix_gene_id"])
        coverage.append(dict(**source, normalized_modeled_gene_id=target, normalized_source_gene_id=native,
                             modeled_in_gtf=target in normalized_gtf_genes, modeled_in_cdna=target in cdna_genes,
                             modeled_in_union=target in union_genes, source_in_gtf=native in normalized_gtf_genes,
                             source_in_cdna=native in cdna_genes, source_in_union=native in union_genes,
                             modeled_has_confirmed_transcript=target in confirmed_genes,
                             modeled_has_mapping_or_version_diagnostic=target in affected_genes))
    fields = list(records[0])
    write_tsv(output/"fasta_all_record_metadata.tsv", records, fields)
    write_tsv(output/"fasta_gtf_mapping_diagnostics.tsv", diagnostics, fields)
    write_tsv(output/"transcript_union_axis.tsv", sorted(union.values(), key=lambda r:r["transcript_id"]), list(next(iter(union.values()))))
    transcript_rows = [r for rows in gtf.values() for r in rows]
    gene_rows = [r for rows in genes.values() for r in rows]
    write_tsv(output/"gtf_transcript_gene_coordinates.tsv", sorted(transcript_rows, key=lambda r:r["transcript_id"]), list(transcript_rows[0]))
    write_tsv(output/"gtf_gene_coordinates.tsv", sorted(gene_rows, key=lambda r:r["gene_id"]), list(gene_rows[0]))
    write_tsv(output/"modeled_gene_reference_coverage.tsv", coverage, list(coverage[0]))
    write_tsv(output/"missing_modeled_union_genes.tsv", [r for r in coverage if not r["modeled_in_union"]], list(coverage[0]))
    write_tsv(output/"duplicate_transcript_ids.tsv", duplicate_rows, ["transcript_id", "first_source", "duplicate_source", "identity_agrees"])
    write_tsv(output/"modeled_genes_affected_by_diagnostics.tsv", [r for r in coverage if r["modeled_has_mapping_or_version_diagnostic"]], list(coverage[0]))
    return dict(gtf_genes=len(genes), gtf_normalized_genes=len(normalized_gtf_genes),
                gtf_header=annotation_headers,
                gtf_transcripts=len(gtf), gtf_gene_coordinate_rows=len(gene_rows),
                gtf_transcript_coordinate_rows=len(transcript_rows), cdna_transcripts=len(sets["cdna"]),
                ncrna_transcripts=len(sets["ncrna"]), union_transcripts=len(union), duplicate_records=len(duplicate_rows),
                cdna_genes=len(cdna_genes), union_genes=len(union_genes), modeled_genes=len(coverage),
                ambiguous_duplicate_transcript_ids=sorted(ambiguous),
                diagnostic_record_categories=dict(Counter(c for r in diagnostics for c in r["diagnostic_categories"].split(";"))),
                modeled_genes_with_confirmed_transcript=sum(r["modeled_has_confirmed_transcript"] for r in coverage),
                modeled_genes_affected_by_diagnostics=sum(r["modeled_has_mapping_or_version_diagnostic"] for r in coverage),
                modeled_cdna_covered=sum(r["modeled_in_cdna"] for r in coverage),
                modeled_union_covered=sum(r["modeled_in_union"] for r in coverage),
                source_union_covered=sum(r["source_in_union"] for r in coverage),
                missing_modeled_union_genes=[r["normalized_modeled_gene_id"] for r in coverage if not r["modeled_in_union"]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--census-only", action="store_true", help="Reuse only frozen job21998355 assets; no network")
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Acquire/census reference only in compute job")
    plan = source_plan()
    args.output.mkdir(parents=True, exist_ok=False)
    archived = args.output/"source_metadata"; archived.mkdir()
    for path in sorted(META.iterdir()):
        require(path.is_file(), "Unexpected metadata subdirectory")
        shutil.copyfile(path, archived/path.name)
    if args.census_only:
        receipt_path = REUSE / "acquisition_receipt.json"
        require(sha(receipt_path) == REUSE_RECEIPT_SHA, "Frozen acquisition receipt changed")
        receipts = json.loads(receipt_path.read_text())
        require(len(receipts) == 3 and {r["kind"] for r in receipts} == set(REUSE_ASSET_SHA), "Unexpected reuse receipt roster")
        for row in plan:
            prior = next(r for r in receipts if r["kind"] == row["kind"])
            require(all(prior[k] == row[k] for k in row), "Acquired asset metadata differs from frozen source plan")
            require(prior["sha256"] == REUSE_ASSET_SHA[row["kind"]], "Acquired SHA differs from frozen protocol")
            asset = REUSE / row["filename"]
            require(asset.stat().st_size == row["bytes"] and sha(asset) == prior["sha256"], "Reuse asset size/SHA mismatch")
        shutil.copyfile(receipt_path, args.output/"acquisition_receipt.json")
    else:
        receipts = []
        for row in plan:
            receipts.append(acquire(row, args.output))
            (args.output/"acquisition_receipt.json").write_text(json.dumps(receipts, indent=2)+"\n")
    stats = census(args.output, plan, REUSE if args.census_only else None)
    launcher = "run_ensembl98_reference_census.sbatch" if args.census_only else "run_acquire_ensembl98_reference.sbatch"
    report = dict(schema_version="ensembl98-reference-census-v2", acquisition=receipts, coverage=stats,
                  census_only=args.census_only, asset_root=str(REUSE if args.census_only else args.output),
                  reused_acquisition_receipt_sha256=REUSE_RECEIPT_SHA if args.census_only else None,
                  modeled_gene_quantification_admitted=False,
                  input_axis=dict(path=str(AXIS), sha256=AXIS_SHA), source_metadata=str(META),
                  script_sha256=sha(Path(__file__)), launcher_sha256=sha(Path(__file__).with_name(launcher)),
                  quantification_flags_selected=False, index_built=False, sequencing_or_count_values_acquired=False,
                  exact_private_pisces_reconstruction_asserted=False, zero_missing_rows_manufactured=False,
                  normalization="upper-case; strip terminal _PAR_Y, then terminal dot-version, matching released scorer",
                  duplicate_rule="deduplicate identical raw transcript IDs only if gene identity, sequence SHA256 and length agree; distinct raw IDs and GTF locus geometries preserved",
                  limits=dict(compressed_total=MAX_TOTAL, compressed_file=MAX_FILE, sequence_record=MAX_SEQUENCE,
                              text_line=MAX_LINE, transcript_records=MAX_TRANSCRIPTS, gene_records=MAX_GENES),
                  environment=dict(python=sys.version, platform=platform.platform(), slurm_job_id=os.environ["SLURM_JOB_ID"]),
                  output_sha256={str(p.relative_to(args.output)):sha(p) for p in sorted(args.output.rglob("*")) if p.is_file()})
    (args.output/"summary.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:v for k,v in stats.items() if not isinstance(v,list)}, sort_keys=True))


if __name__ == "__main__":
    main()
