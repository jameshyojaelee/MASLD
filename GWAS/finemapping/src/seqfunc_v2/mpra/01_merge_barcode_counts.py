#!/usr/bin/env python3
"""Merge the eight Hu MPRA count files for one cell line in lockstep.

The GEO files share an oligo/barcode row order. This script asserts that
contract on every row and emits a wide barcode-by-sample table without a lossy
pre-aggregation, preserving the substrate required by MPRAnalyze.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from contextlib import ExitStack
from pathlib import Path


def sample_info(path: Path) -> dict[str, str | int]:
    name = path.name
    if "G2_" in name:
        cell = "HepG2"
        condition = "PAOA" if "_PAOA_" in name else "control"
    elif "Lx2_" in name:
        cell = "LX2"
        condition = "TGFb" if "_TGFb_" in name else "control"
    else:
        raise ValueError(f"Unrecognized MPRA filename: {name}")
    match = re.search(r"[Rr]ep(\d+)", name)
    if not match:
        raise ValueError(f"Missing replicate in {name}")
    rep = int(match.group(1))
    return {"sample_id": f"{cell}_{condition}_r{rep}", "cell_line": cell, "condition": condition, "replicate": rep}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-line", choices=["HepG2", "LX2"], required=True)
    parser.add_argument("--threads", type=int, default=16, help="Recorded resource contract; decompression is multi-stream")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(os.environ.get("MASLD_PROJECT_ROOT", Path(__file__).resolve().parents[5])),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    src = root / "GWAS/finemapping/data/seqfunc_external/hu2025_mpra/raw_counts_v2"
    out = root / "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/counts"
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(src.glob("*DNA_RNA_counts.txt.gz"))
    selected = []
    for path in files:
        info = sample_info(path)
        if info["cell_line"] == args.cell_line:
            selected.append((path, info))
    selected.sort(key=lambda x: (str(x[1]["condition"] != "control"), int(x[1]["replicate"])))
    if len(selected) != 8:
        raise SystemExit(f"Expected 8 {args.cell_line} files, found {len(selected)}")
    if sorted((x[1]["condition"], x[1]["replicate"]) for x in selected) != sorted(
        [("control", i) for i in range(1, 5)] + [(("PAOA" if args.cell_line == "HepG2" else "TGFb"), i) for i in range(1, 5)]
    ):
        raise SystemExit("Cell-line condition/replicate manifest is incomplete")

    manifest_fields = ["sample_id", "cell_line", "condition", "replicate", "path", "size", "sha256"]
    manifest = []
    for path, info in selected:
        manifest.append({**info, "path": str(path), "size": path.stat().st_size, "sha256": sha256(path)})
    manifest_path = out / f"{args.cell_line}.sample_manifest.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)

    output = out / f"{args.cell_line}.barcode_counts.tsv.gz"
    tmp = output.with_suffix(output.suffix + ".tmp")
    n_rows = 0
    n_mut = 0
    with ExitStack() as stack:
        handles = [stack.enter_context(gzip.open(path, "rt")) for path, _ in selected]
        with gzip.open(tmp, "wt", compresslevel=1, newline="") as dst:
            header = ["construct", "element_id", "allele", "barcode"]
            for _, info in selected:
                header.extend([f"dna__{info['sample_id']}", f"rna__{info['sample_id']}"])
            dst.write("\t".join(header) + "\n")
            while True:
                lines = [handle.readline() for handle in handles]
                if not any(lines):
                    break
                if not all(lines):
                    raise SystemExit(f"Input row-count mismatch after {n_rows} rows")
                parsed = [line.rstrip("\n").split("\t") for line in lines]
                if any(len(row) != 4 for row in parsed):
                    raise SystemExit(f"Malformed count row after {n_rows} rows")
                key = parsed[0][:2]
                if any(row[:2] != key for row in parsed[1:]):
                    raise SystemExit(f"Oligo/barcode row-order mismatch after {n_rows} rows: {key}")
                construct, barcode = key
                is_mut = construct.endswith("_Mut")
                element = construct[:-4] if is_mut else construct
                allele = "alt" if is_mut else "ref"
                values = [construct, element, allele, barcode]
                for row in parsed:
                    int(row[2]); int(row[3])
                    values.extend(row[2:4])
                dst.write("\t".join(values) + "\n")
                n_rows += 1
                n_mut += int(is_mut)
    os.replace(tmp, output)
    contract = {
        "cell_line": args.cell_line,
        "n_samples": 8,
        "n_rows": n_rows,
        "n_ref_barcode_rows": n_rows - n_mut,
        "n_alt_barcode_rows": n_mut,
        "row_order_asserted_every_row": True,
        "aggregation": "none",
        "output": str(output),
        "output_sha256": sha256(output),
        "threads_requested": args.threads,
    }
    contract_path = out / f"{args.cell_line}.merge_contract.json"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2))


if __name__ == "__main__":
    main()
