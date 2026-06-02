#!/usr/bin/env python3
import csv
import gzip
import io
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


DATASETS = {
    "GSE156918": {
        "runinfo_term": "SRP278840",
        "series_prefix": "GSE156nnn",
        "controls": {"CHOW"},
        "mcd_tokens": {"MCD"},
        "diet_field": "diet",
        "genotype_field": "genotype",
        "exclude_if_treatment_contains": None,
    },
    "GSE205974": {
        "runinfo_term": "PRJNA848688",
        "series_prefix": "GSE205nnn",
        "controls": {"MCS"},
        "mcd_tokens": {"MCD"},
        "diet_field": "treatment",
        "genotype_field": "genotype",
        "exclude_if_treatment_contains": {"HGPP"},
    },
}


def strip_quotes(value):
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def fetch_series_matrix(gse, series_prefix):
    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{series_prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"
    )
    with urllib.request.urlopen(url) as handle:
        data = handle.read()
    return gzip.decompress(data).decode("utf-8", errors="replace")


def parse_series_matrix(text):
    samples = []
    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = line.split("\t")
        key = parts[0][len("!Sample_") :]
        values = [strip_quotes(v) for v in parts[1:]]
        if key == "geo_accession":
            samples = [
                {"gsm": gsm, "characteristics": [], "relations": []}
                for gsm in values
            ]
            continue
        if not samples:
            continue
        for idx, value in enumerate(values):
            if key.startswith("characteristics"):
                samples[idx]["characteristics"].append(value)
            elif key == "relation":
                samples[idx]["relations"].append(value)
            else:
                samples[idx][key] = value
    return samples


def extract_characteristic(characteristics, label):
    label_lower = f"{label.lower()}:"
    for item in characteristics:
        if item.lower().startswith(label_lower):
            return item.split(":", 1)[1].strip()
    return None


def extract_biosample(relations):
    for entry in relations:
        match = re.search(r"SAMN\d+", entry)
        if match:
            return match.group(0)
    return None


def extract_srx(relations):
    for entry in relations:
        match = re.search(r"SRX\d+", entry)
        if match:
            return match.group(0)
    return None


def fetch_runinfo(term):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    query = urllib.parse.quote(term)
    url = f"{base}esearch.fcgi?db=sra&term={query}"
    root = ET.fromstring(urllib.request.urlopen(url, timeout=30).read())
    count_text = root.findtext(".//Count", default="0")
    count = int(count_text)
    ids = [elem.text for elem in root.findall(".//Id")]
    if count > len(ids):
        url = f"{base}esearch.fcgi?db=sra&term={query}&retmax={count}"
        root = ET.fromstring(urllib.request.urlopen(url, timeout=60).read())
        ids = [elem.text for elem in root.findall(".//Id")]
    if not ids:
        raise RuntimeError(f"No SRA IDs found for term {term}")
    url2 = (
        f"{base}efetch.fcgi?db=sra&id={','.join(ids)}"
        "&rettype=runinfo&retmode=text"
    )
    raw = urllib.request.urlopen(url2, timeout=60).read().decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(raw))
    return list(reader)


def diet_from_field(value, controls, mcd_tokens, exclude_tokens):
    if value is None:
        return None
    upper = value.upper()
    if exclude_tokens and any(token in upper for token in exclude_tokens):
        return None
    if any(token in upper for token in mcd_tokens):
        return "MCD"
    if any(token in upper for token in controls):
        return "Control"
    return None


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    root = Path(__file__).resolve().parents[1]
    for gse_id, cfg in DATASETS.items():
        dataset_root = root / gse_id
        if not dataset_root.exists():
            print(f"Skipping {gse_id}: {dataset_root} not found")
            continue

        print(f"Processing {gse_id}...")
        matrix = fetch_series_matrix(gse_id, cfg["series_prefix"])
        samples = parse_series_matrix(matrix)
        runinfo = fetch_runinfo(cfg["runinfo_term"])

        run_by_biosample = {}
        for row in runinfo:
            biosample = row.get("BioSample")
            if biosample:
                run_by_biosample.setdefault(biosample, []).append(row)

        selected = []
        for sample in samples:
            characteristics = sample.get("characteristics", [])
            relations = sample.get("relations", [])
            diet_field = extract_characteristic(characteristics, cfg["diet_field"])
            diet = diet_from_field(
                diet_field,
                controls=cfg["controls"],
                mcd_tokens=cfg["mcd_tokens"],
                exclude_tokens=cfg["exclude_if_treatment_contains"],
            )
            if diet is None:
                continue
            genotype = extract_characteristic(characteristics, cfg["genotype_field"])
            biosample = extract_biosample(relations)
            srx = extract_srx(relations)
            runs = []
            if biosample:
                runs = run_by_biosample.get(biosample, [])
            if not runs and srx:
                runs = [row for row in runinfo if row.get("Experiment") == srx]
            if not runs:
                label = biosample or srx or "unknown"
                print(f"Warning: no SRA run for {sample.get('gsm')} ({label})")
                continue
            run_ids = sorted({row["Run"] for row in runs})
            layout = runs[0].get("LibraryLayout", "")
            if layout == "PAIRED":
                fastq_r1 = f"fastq/{sample.get('gsm')}_1.fastq.gz"
                fastq_r2 = f"fastq/{sample.get('gsm')}_2.fastq.gz"
            else:
                fastq_r1 = f"fastq/{sample.get('gsm')}.fastq.gz"
                fastq_r2 = ""

            selected.append(
                {
                    "sample_id": sample.get("gsm"),
                    "title": sample.get("title"),
                    "diet": diet,
                    "diet_detail": diet_field,
                    "genotype": genotype,
                    "biosample": biosample,
                    "runs": ";".join(run_ids),
                    "layout": layout,
                    "fastq_r1": fastq_r1,
                    "fastq_r2": fastq_r2,
                }
            )

        selected.sort(key=lambda row: (row["diet"], row["sample_id"]))
        metadata_dir = dataset_root / "metadata"
        write_tsv(
            metadata_dir / "samples.tsv",
            selected,
            [
                "sample_id",
                "diet",
                "diet_detail",
                "genotype",
                "biosample",
                "runs",
                "layout",
                "fastq_r1",
                "fastq_r2",
                "title",
            ],
        )

        # Save a small manifest for downloads.
        runs_path = metadata_dir / "sra_runs.txt"
        run_ids_all = sorted(
            {run for row in selected for run in row["runs"].split(";") if run}
        )
        with runs_path.open("w") as handle:
            for run in run_ids_all:
                handle.write(f"{run}\n")

        # Save filtered runinfo for record-keeping.
        runinfo_rows = [row for row in runinfo if row.get("Run") in set(run_ids_all)]
        if runinfo_rows:
            runinfo_path = metadata_dir / "runinfo_selected.csv"
            with runinfo_path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=runinfo_rows[0].keys())
                writer.writeheader()
                writer.writerows(runinfo_rows)

        print(f"Wrote {len(selected)} samples to {metadata_dir / 'samples.tsv'}")


if __name__ == "__main__":
    main()
