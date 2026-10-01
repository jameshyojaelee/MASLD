"""Build a GSE128072 paired raw-run roster without acquiring sequencing files."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import sys
from urllib.parse import urlsplit
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[4]
SOFT = ROOT / "GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/source_metadata/GSE128072_family.soft.gz"
SOFT_SHA = "92d47ccb362fec5913e899d765982f413c1b39327ab92b80a40af77536030def"
FIELDS = ("run_accession", "experiment_accession", "sample_accession", "library_strategy",
          "library_layout", "read_count", "base_count", "fastq_ftp", "fastq_bytes",
          "fastq_md5", "sra_bytes")
ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA526249"
           "&result=read_run&fields=" + ",".join(FIELDS) + "&format=tsv&download=false")
ENA_CAP = 2_000_000
SOFT_TEXT_CAP = 16_000_000
EXPECTED_DONORS = {"B1", "B7", "B8", "B9", "B10", "B12", "B15", "B19", "B21",
                   "B22", "B24", "B26", "B36", "B38", "B41", "B46", "B47"}
TITLE = re.compile(r"^(B\d+)(-H3K27ac_ChIP-Seq|-Input_ChIP-Seq|_RNA-Seq) \[([^\[\]]+)\]$")
TOKEN = re.compile(r"^(FGC\d+)_s_(\d+)(_[0-9]+)?_([ACGTN]+-[ACGTN]+)$")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_tsv(path, rows, fields):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def parse_soft(data):
    require(digest(data) == SOFT_SHA, "Frozen SOFT checksum changed")
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as handle:
        text = handle.read(SOFT_TEXT_CAP + 1)
    require(len(text) <= SOFT_TEXT_CAP, "SOFT metadata exceeds decompression cap")
    samples = []
    all_gsm = set()
    metadata_only_rna = []
    for block in text.decode("utf-8").split("^SAMPLE = ")[1:]:
        gsm = block.splitlines()[0].strip()
        require(gsm not in all_gsm, "Duplicate GEO sample")
        all_gsm.add(gsm)
        metadata = defaultdict(list)
        for line in block.splitlines()[1:]:
            if line.startswith(("!Sample_characteristics_ch1 = ", "!Sample_title = ", "!Sample_relation = ")):
                key, value = line.split(" = ", 1)
                metadata[key].append(value)
        characteristics = defaultdict(list)
        for value in metadata["!Sample_characteristics_ch1"]:
            key, separator, val = value.partition(": ")
            if separator:
                characteristics[key].append(val)
        if characteristics["cohort"] != ["Penn Cohort 1"]:
            continue
        require(len(metadata["!Sample_title"]) == 1, "Ambiguous source title")
        title = metadata["!Sample_title"][0]
        # H3K4me3 is excluded by source assay title, before intersecting donors.
        match = TITLE.fullmatch(title)
        if not match:
            if gsm == "GSM3663220" and title == "B50_living tissue donor_RNA-Seq":
                require(characteristics["individual"] == ["B50"], "Metadata-only donor disagreement")
                require(not any(v.startswith("SRA:") for v in metadata["!Sample_relation"]),
                        "Metadata-only RNA record now has a sequencing relation")
                metadata_only_rna.append(dict(donor="B50", GSM=gsm, title=title,
                                              exclusion="RNA-only donor; no published SRX relation"))
                continue
            require("-H3K4me3_ChIP-Seq [" in title, "Unexpected Penn Cohort 1 assay title")
            continue
        donor, assay_title, token = match.groups()
        native = TOKEN.fullmatch(token)
        require(native is not None, "Unexpected native title token format")
        require(characteristics["individual"] == [donor], "Title/individual donor disagreement")
        assay = {"_RNA-Seq": "RNA", "-H3K27ac_ChIP-Seq": "H3K27ac",
                 "-Input_ChIP-Seq": "Input"}[assay_title]
        relations = [re.fullmatch(r"SRA: https://www\.ncbi\.nlm\.nih\.gov/sra\?term=(SRX\d+)", v)
                     for v in metadata["!Sample_relation"] if v.startswith("SRA:")]
        require(len(relations) == 1 and relations[0] is not None, "Missing/ambiguous SRX relation")
        biosamples = [re.fullmatch(r"BioSample: https://www\.ncbi\.nlm\.nih\.gov/biosample/(SAMN\d+)", v)
                      for v in metadata["!Sample_relation"] if v.startswith("BioSample:")]
        require(len(biosamples) == 1 and biosamples[0] is not None, "Missing/ambiguous BioSample relation")
        samples.append(dict(donor=donor, assay=assay, GSM=gsm,
                            experiment_accession=relations[0].group(1), title=title,
                            source_biosample=biosamples[0].group(1), native_library_token=token,
                            native_FGC_title_token=native.group(1), native_lane_title_token=native.group(2),
                            native_extra_title_token=native.group(3) or "", native_barcode_title_token=native.group(4)))
    rna = {r["donor"] for r in samples if r["assay"] == "RNA"} | {r["donor"] for r in metadata_only_rna}
    chip = {r["donor"] for r in samples if r["assay"] == "H3K27ac"}
    require(not chip.intersection(r["donor"] for r in metadata_only_rna),
            "Paired donor lacks a published RNA run relation")
    donors = rna & chip
    selected = [r for r in samples if r["donor"] in donors]
    require(donors == EXPECTED_DONORS, "Derived paired donor roster changed")
    require(len(selected) == 468, "Derived selected GSM count changed")
    require(len({r["experiment_accession"] for r in selected}) == len(selected), "Duplicate selected SRX")
    return selected, dict(rna_only=sorted(rna-chip), h3k27ac_only=sorted(chip-rna),
                         paired_donors=sorted(donors, key=lambda x: int(x[1:])),
                         metadata_only_rna=metadata_only_rna)


def fetch_ena():
    # This is the sole network operation. No FASTQ/SRA URL is opened.
    with urlopen(ENA_URL, timeout=60) as response:
        final_url = response.geturl()
        parsed = urlsplit(final_url)
        require(parsed.hostname == "www.ebi.ac.uk" and parsed.path == "/ena/portal/api/filereport",
                "Unexpected ENA metadata redirect")
        raw = response.read(ENA_CAP + 1)
        headers = {k: response.headers.get(k) for k in ("Content-Type", "Content-Length", "ETag", "Last-Modified")}
    require(len(raw) <= ENA_CAP, "ENA response exceeds metadata cap")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
    require(tuple(reader.fieldnames or ()) == FIELDS, "ENA metadata header changed")
    rows = list(reader)
    require(all(set(r) == set(FIELDS) and None not in r.values() for r in rows), "Malformed ENA metadata row")
    return raw, rows, dict(url=ENA_URL, final_url=final_url, bytes=len(raw), sha256=digest(raw),
                          response_headers=headers, response_cap_bytes=ENA_CAP)


def build_roster(selected, ena):
    by_experiment = defaultdict(list)
    for row in ena:
        by_experiment[row["experiment_accession"]].append(row)
    runs, files = [], []
    seen_runs, seen_files = set(), set()
    for source in sorted(selected, key=lambda r: (int(r["donor"][1:]), r["assay"], r["GSM"])):
        matches = by_experiment[source["experiment_accession"]]
        require(len(matches) == 1, "Missing/duplicate selected SRX-to-run join")
        row = matches[0]
        run = row["run_accession"]
        require(re.fullmatch(r"SRR\d+", run) and run not in seen_runs, "Invalid/duplicate selected SRR")
        seen_runs.add(run)
        require(re.fullmatch(r"SAMN\d+", row["sample_accession"]), "Missing native BioSample ID")
        require(row["sample_accession"] == source["source_biosample"], "GEO/ENA BioSample disagreement")
        layout = "PAIRED" if source["assay"] == "RNA" else "SINGLE"
        strategy = "RNA-Seq" if source["assay"] == "RNA" else "ChIP-Seq"
        require(row["library_layout"] == layout and row["library_strategy"] == strategy, "Assay/layout disagreement")
        for field in ("read_count", "base_count"):
            require(row[field].isdigit() and int(row[field]) > 0, "Missing/invalid " + field)
        urls = row["fastq_ftp"].split(";")
        sizes = row["fastq_bytes"].split(";")
        checksums = row["fastq_md5"].split(";")
        nfiles = 2 if layout == "PAIRED" else 1
        require(len(urls) == len(sizes) == len(checksums) == nfiles, "Incomplete FASTQ file fields")
        for number, (url, size, md5) in enumerate(zip(urls, sizes, checksums), 1):
            require(url.startswith("ftp.sra.ebi.ac.uk/") and url.endswith(".fastq.gz") and url not in seen_files,
                    "Invalid/duplicate published FASTQ URL")
            require(size.isdigit() and int(size) > 0 and re.fullmatch(r"[0-9a-fA-F]{32}", md5), "Invalid file size/checksum")
            seen_files.add(url)
            files.append(dict(**source, run_accession=run, file_number=number,
                              fastq_url=url, fastq_bytes=int(size), fastq_md5=md5.lower()))
        runs.append(dict(**source, **{k: v for k, v in row.items() if k != "experiment_accession"},
                         fastq_file_count=nfiles, total_fastq_bytes=sum(map(int, sizes))))
    require(len(runs) == 468 and len(files) == 604, "Derived run/file totals changed")
    require(sum(r["fastq_bytes"] for r in files) == 187_002_556_838, "Derived compressed-byte total changed")
    require(Counter(r["assay"] for r in runs) == {"RNA": 136, "H3K27ac": 202, "Input": 130}, "Derived assay totals changed")
    return runs, files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Run this metadata producer in its compute allocation")
    selected, selection = parse_soft(SOFT.read_bytes())
    args.output.mkdir(parents=True, exist_ok=False)
    raw, ena, receipt = fetch_ena()
    (args.output / "ena_read_run_metadata.tsv").write_bytes(raw)
    (args.output / "ena_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    runs, files = build_roster(selected, ena)
    source_fields = list(selected[0])
    write_tsv(args.output / "source_samples.tsv", selected, source_fields)
    write_tsv(args.output / "raw_runs.tsv", runs, list(runs[0]))
    write_tsv(args.output / "raw_files.tsv", files, list(files[0]))
    totals = []
    for donor in selection["paired_donors"]:
        for assay in ("RNA", "H3K27ac", "Input"):
            group = [r for r in runs if r["donor"] == donor and r["assay"] == assay]
            require(group, "Missing donor assay")
            totals.append(dict(donor=donor, assay=assay, n_runs=len(group),
                               n_title_tokens=len({r["native_library_token"] for r in group}),
                               n_fastq_files=sum(r["fastq_file_count"] for r in group),
                               fastq_bytes=sum(r["total_fastq_bytes"] for r in group),
                               read_count=sum(int(r["read_count"]) for r in group),
                               base_count=sum(int(r["base_count"]) for r in group)))
    write_tsv(args.output / "donor_assay_totals.tsv", totals, list(totals[0]))
    summary = dict(schema_version="gse128072-paired-raw-metadata-v1", selection=selection,
                   n_donors=len(selection["paired_donors"]), n_runs=len(runs), n_fastq_files=len(files),
                   total_fastq_bytes=sum(r["fastq_bytes"] for r in files),
                   assay_run_counts=dict(Counter(r["assay"] for r in runs)),
                   source_soft=dict(path=str(SOFT), sha256=SOFT_SHA), ena_receipt=receipt,
                   script_sha256=digest(Path(__file__).read_bytes()),
                   launcher_sha256=digest(Path(__file__).with_name("run_paired_h3k27ac_raw_roster.sbatch").read_bytes()),
                   native_library_token_semantics="Verbatim bracket token from source title; not an established unique biological library or Input-to-H3K27ac match",
                   read_unit="RNA read_count denotes paired spots; ChIP/Input single spots; sequencing runs are not biological donors",
                   raw_sequencing_acquired=False, protected_data_read=False,
                   environment=dict(python=sys.version, platform=platform.platform(), slurm_job_id=os.environ["SLURM_JOB_ID"]),
                   output_sha256={p.name: digest(p.read_bytes()) for p in sorted(args.output.iterdir()) if p.is_file()})
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("n_donors", "n_runs", "n_fastq_files", "total_fastq_bytes", "raw_sequencing_acquired")}, sort_keys=True))


if __name__ == "__main__":
    main()
