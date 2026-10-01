#!/usr/bin/env python3
"""Check the public paired liver RNA/H3K27ac lead without opening protected effects."""
import argparse
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def soft_samples(path):
    samples = []
    sample = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                if sample:
                    samples.append(sample)
                sample = {"gsm": line.strip().split(" = ", 1)[1], "features": []}
            elif sample is not None and line.startswith("!Sample_characteristics_ch1 = "):
                sample["features"].append(line.strip().split(" = ", 1)[1])
            elif sample is not None and line.startswith("!Sample_library_strategy = "):
                sample["assay"] = line.strip().split(" = ", 1)[1]
        if sample:
            samples.append(sample)
    return samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft", type=Path, required=True)
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--h3", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    metadata = soft_samples(args.soft)
    per = defaultdict(set)
    for sample in metadata:
        fields = "|".join(sample["features"])
        match = re.search(r"individual: ([^|]+)", fields)
        if not match:
            continue
        donor = match.group(1)
        if "H3K27ac" in fields:
            per["h3"].add(donor)
        if sample.get("assay") == "RNA-Seq":
            per["rna"].add(donor)
    with gzip.open(args.rna, "rt") as handle:
        rna = handle.readline().rstrip("\n").split("\t")
    with gzip.open(args.h3, "rt") as handle:
        h3 = [x.removesuffix("-H3K27ac") for x in handle.readline().rstrip("\n").split("\t")[6:]]
        first = handle.readline().rstrip("\n").split("\t")
    if not first[1].startswith("chr"):
        raise ValueError("H3K27ac interval header/assembly assumptions changed")
    if len(rna) != len(set(rna)) or len(h3) != len(set(h3)):
        raise ValueError("Repeated donor column")
    paired = sorted(set(rna) & set(h3) & per["rna"] & per["h3"])
    if len(paired) < 10:
        raise ValueError("Insufficient paired public donors for the prespecified batch transformation")
    report = {
        "source": "GSE128072 Penn Cohort 1, normal deceased donor liver",
        "sample_records": len(metadata), "rna_columns": len(rna), "h3_columns": len(h3),
        "paired_distinct_donors": len(paired), "paired_donors": paired,
        "assay": "RNA TPM and H3K27ac ChIP feature counts",
        "h3_source_assembly": "hg19 per GEO; requires interval liftover to GSE267145 GRCh38 before region evaluation",
        "rna_unit_mismatch": "TPM is not the GSE267145 deposited fractional RNA estimate or source count; direct released-counts scoring is not valid",
        "histology": "not supplied for these paired donors; external within-histology claim unavailable",
        "external_status": "candidate normal-liver H3K27ac transport; no performance measured",
        "source_shas": {str(p): sha256(p) for p in (args.soft, args.rna, args.h3)},
    }
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("sample_records", "rna_columns", "h3_columns", "paired_distinct_donors", "external_status")}))


if __name__ == "__main__":
    main()
