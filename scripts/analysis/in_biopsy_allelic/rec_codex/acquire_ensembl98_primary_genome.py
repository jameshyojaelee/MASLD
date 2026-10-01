"""Acquire one frozen public genome and faidx it; no quantification index."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
META = REC / "ensembl98_primary_genome_metadata_20260930T222808Z"
RECEIPT_SHA = "1918a50008667cfb58f25db901b4e3c1ca62aadcb065e6c279c5f0cc81224927"
URL = "https://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
SIZE = 881211416
MAX_COMPRESSED = 1_000_000_000
MAX_UNCOMPRESSED = 3_300_000_000
SAMTOOLS = Path("/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools")
H3_AXIS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/h3k27ac_feature_axis.tsv"
H3_SHA = "cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986"
GTF_GENES = REC / "ensembl98_reference_census_21998377/gtf_gene_coordinates.tsv"
GTF_SHA = "594e4b78de02adc274d91d28954fbceffadba2481ca3b79e02e6da46ecc96e41"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Genome acquisition/index only on compute")
    require(sha(META/"receipt.json") == RECEIPT_SHA, "Frozen metadata receipt changed")
    receipt = json.loads((META/"receipt.json").read_text())
    for entry in receipt[:2]:
        name = entry["url"].rsplit("/", 1)[1]
        require(entry["status"] == 200 and sha(META/(name+".txt")) == entry["sha256"], "Frozen source metadata changed")
    require(len(receipt) == 3 and receipt[2]["url"] == URL and receipt[2]["status"] == 200 and
            receipt[2]["content_length"] == SIZE and SIZE <= MAX_COMPRESSED, "Unexpected genome HEAD")
    require("GCA_000001405.28" in (META/"README.txt").read_text(), "Source assembly accession differs")
    checksum = [r.split() for r in (META/"CHECKSUMS.txt").read_text().splitlines()
                if r.split() and r.split()[-1] == URL.rsplit("/", 1)[1]]
    require(len(checksum) == 1 and len(checksum[0]) == 3, "Ambiguous published genome checksum")
    expected_sum = (int(checksum[0][0]), int(checksum[0][1]))
    require(sha(H3_AXIS) == H3_SHA and sha(GTF_GENES) == GTF_SHA, "Fixed coordinate metadata changed")
    version = subprocess.run([str(SAMTOOLS), "--version"], check=True, capture_output=True, text=True).stdout
    require(version.splitlines()[0] == "samtools 1.22.1", "Installed samtools version changed")
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(META, args.output/"source_metadata")
    compressed = args.output/URL.rsplit("/", 1)[1]
    partial = compressed.with_name(compressed.name+".partial")
    digest, nbytes = hashlib.sha256(), 0
    with urlopen(URL, timeout=60) as response, partial.open("xb") as handle:
        require(response.status == 200 and response.geturl() == URL and
                int(response.headers.get("Content-Length", -1)) == SIZE, "Unexpected genome download response")
        response_headers = dict(response.headers)
        for block in iter(lambda: response.read(1 << 20), b""):
            nbytes += len(block)
            require(nbytes <= SIZE and nbytes <= MAX_COMPRESSED, "Compressed genome cap exceeded; partial preserved")
            digest.update(block)
            handle.write(block)
    require(nbytes == SIZE, "Truncated genome download; partial preserved")
    sum_result = subprocess.run(["sum", "-r", str(partial)], check=True, capture_output=True, text=True).stdout
    require(tuple(map(int, sum_result.split()[:2])) == expected_sum, "Published BSD checksum differs; partial preserved")
    require(sha(partial) == digest.hexdigest(), "Downloaded file differs from streaming SHA")
    partial.rename(compressed)
    acquisition = dict(url=URL, bytes=nbytes, sha256=digest.hexdigest(), published_bsd_sum=expected_sum,
                       unix_sum=sum_result.strip(), response_headers=response_headers,
                       source_receipt_sha256=RECEIPT_SHA)
    (args.output/"acquisition_receipt.json").write_text(json.dumps(acquisition, indent=2)+"\n")
    fasta = args.output/compressed.name.removesuffix(".gz")
    fasta_partial = fasta.with_name(fasta.name+".partial")
    decompressed_sha, uncompressed = hashlib.sha256(), 0
    with gzip.open(compressed, "rb") as source, fasta_partial.open("xb") as dest:
        for block in iter(lambda: source.read(1 << 20), b""):
            require(uncompressed+len(block) <= MAX_UNCOMPRESSED,
                    "Uncompressed genome exceeds3.3GB; partial preserved, no faidx")
            uncompressed += len(block)
            decompressed_sha.update(block)
            dest.write(block)
    require(uncompressed > 0, "Empty genome")
    fasta_partial.rename(fasta)
    subprocess.run([str(SAMTOOLS), "faidx", str(fasta)], check=True)
    fai = Path(str(fasta)+".fai")
    contigs = {}
    for line in fai.read_text().splitlines():
        cols = line.split("\t")
        require(len(cols) == 5 and cols[0] not in contigs, "Malformed/duplicate genome index contig")
        contigs[cols[0]] = int(cols[1])
    h3_counts, mismatches = {}, []
    with H3_AXIS.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["opaque_source_feature_key"]
            match = re.fullmatch(r"chr(\d+|X|Y):(\d+)-(\d+)", key)
            require(match is not None, "Fixed H3 chromosome convention changed")
            name, start, end = match[1], int(match[2]), int(match[3])
            h3_counts[name] = h3_counts.get(name, 0)+1
            if name not in contigs or not 1 <= start <= end <= contigs[name]:
                mismatches.append(dict(source="fixed_h3", identifier=key, seqname=name, start=start, end=end))
    gtf_rows = 0
    with GTF_GENES.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gtf_rows += 1
            name, start, end = row["seqname"], int(row["start_1based"]), int(row["end_inclusive"])
            if name not in contigs or not 1 <= start <= end <= contigs[name]:
                mismatches.append(dict(source="ensembl98_gtf_gene", identifier=row["gene_id"], seqname=name, start=start, end=end))
    (args.output/"coordinate_mismatches.json").write_text(json.dumps(mismatches, indent=2)+"\n")
    summary = dict(acquisition=acquisition, fasta_bytes=uncompressed, fasta_sha256=decompressed_sha.hexdigest(),
                   faidx_sha256=sha(fai), contig_lengths=contigs, h3_region_counts=h3_counts,
                   h3_regions=sum(h3_counts.values()), gtf_gene_rows=gtf_rows, coordinate_mismatches=len(mismatches),
                   samtools_version=version, samtools_binary_sha256=sha(SAMTOOLS), script_sha256=sha(Path(__file__)),
                   launcher_sha256=sha(Path(__file__).with_name("run_acquire_ensembl98_primary_genome.sbatch")),
                   input_axis_sha256=H3_SHA, gtf_gene_metadata_sha256=GTF_SHA, python=sys.version,
                   slurm_job_id=os.environ["SLURM_JOB_ID"], caps=dict(compressed=MAX_COMPRESSED, uncompressed=MAX_UNCOMPRESSED),
                   indexes_built=["samtools_faidx_only"], transcript_reconstruction=False,
                   sequencing_or_count_acquisition=False, private_pisces_reconstruction_asserted=False,
                   quantification_flags_selected=False)
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:summary[k] for k in ("fasta_bytes", "fasta_sha256", "h3_regions", "coordinate_mismatches")}, sort_keys=True))


if __name__ == "__main__":
    main()
