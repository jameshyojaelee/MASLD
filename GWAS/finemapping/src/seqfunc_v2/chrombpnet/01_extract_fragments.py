#!/usr/bin/env python3
"""Parallel extraction into pooled and donor-intact pseudoreplicates."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def load_barcodes(path: Path) -> set[str]:
    with path.open() as handle:
        values = {line.strip() for line in handle if line.strip()}
    if not values:
        raise RuntimeError(f"empty barcode allowlist: {path}")
    return values


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_donor(donor: dict[str, str], part_dir: str) -> dict[str, object]:
    """One donor per process: independent source file and compression streams."""
    donor_id = donor["donor_id"]
    allow = load_barcodes(Path(donor["barcode_file"]))
    part = Path(part_dir)
    frag_path = part / f"{donor_id}.fragments.tsv.gz"
    tag_path = part / f"{donor_id}.tn5.tagAlign.gz"
    seen_barcodes = set()
    n_total = n_kept = multiplicity = 0
    with gzip.open(donor["fragment_file"], "rt") as source, \
         gzip.open(frag_path, "wt", compresslevel=6) as frag_out, \
         gzip.open(tag_path, "wt", compresslevel=6) as tag_out:
        for line in source:
            n_total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise RuntimeError(f"malformed fragment row in {donor_id}: {line[:100]}")
            chrom, start_s, end_s, barcode, count_s = fields[:5]
            if barcode not in allow or chrom not in AUTOSOMES:
                continue
            start, end, count = int(start_s), int(end_s), int(count_s)
            if end <= start + 9:
                continue
            seen_barcodes.add(barcode)
            n_kept += 1
            multiplicity += count
            frag_out.write(f"{chrom}\t{start}\t{end}\t{barcode}\t{count}\n")
            left, right = start + 4, end - 5
            # Unique-fragment Tn5 ends; PCR/optical multiplicity is QC only.
            tag_out.write(f"{chrom}\t{left}\t{left + 1}\t{barcode}\t1000\t+\n")
            tag_out.write(f"{chrom}\t{right}\t{right + 1}\t{barcode}\t1000\t-\n")
    if not seen_barcodes or n_kept == 0:
        raise RuntimeError(f"no hepatocyte fragments matched for {donor_id}")
    return {
        "donor_id": donor_id,
        "condition": donor["condition"],
        "pseudoreplicate": donor["pseudoreplicate"],
        "allowlist_barcodes": len(allow),
        "barcodes_with_autosomal_fragments": len(seen_barcodes),
        "source_fragment_rows": n_total,
        "retained_fragment_rows": n_kept,
        "retained_fragment_multiplicity": multiplicity,
        "fragment_part": str(frag_path),
        "tag_part": str(tag_path),
    }


def concatenate_gzip_members(paths: list[Path], dest: Path) -> None:
    """Binary concatenation is a standards-compliant multi-member gzip stream."""
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with tmp.open("wb") as target:
        for path in paths:
            with path.open("rb") as source:
                shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
    if tmp.stat().st_size == 0:
        raise RuntimeError(f"empty concatenated product: {tmp}")
    os.replace(tmp, dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "16")))
    args = ap.parse_args()
    out = args.out.resolve()
    manifest = out / "donor_manifest.tsv"
    if not manifest.is_file():
        raise SystemExit("run 00_preflight.py first")
    with manifest.open() as handle:
        donors = list(csv.DictReader(handle, delimiter="\t"))
    if len(donors) != 18:
        raise SystemExit("expected 18 donors")

    frag_dir = out / "fragments"
    part_dir = frag_dir / "donor_parts"
    part_dir.mkdir(parents=True, exist_ok=True)
    # Stale partials cannot be mistaken for completed donor products.
    for path in part_dir.glob("D*.*.gz"):
        path.unlink()
    results = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(donors))) as pool:
        futures = {pool.submit(extract_donor, donor, str(part_dir)): donor["donor_id"] for donor in donors}
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as exc:
                raise SystemExit(f"donor extraction failed for {futures[future]}: {exc}") from exc
    results.sort(key=lambda x: str(x["donor_id"]))
    for result in results:
        result["barcode_fragment_coverage"] = int(result["barcodes_with_autosomal_fragments"]) / int(result["allowlist_barcodes"])
        result["donor_extraction_qc_pass"] = (
            result["barcode_fragment_coverage"] >= 0.75
            and int(result["retained_fragment_rows"]) >= 1_000_000
        )
    failed = [str(r["donor_id"]) for r in results if not r["donor_extraction_qc_pass"]]
    if failed:
        raise SystemExit(f"donors failed extraction QC (coverage>=0.75 and fragments>=1e6): {failed}")

    targets = {"pooled": results}
    for rep in (1, 2):
        targets[f"pseudorep{rep}"] = [r for r in results if int(r["pseudoreplicate"]) == rep]
    products = []
    summaries = {}
    for target, selected in targets.items():
        frag = frag_dir / f"hepatocyte.{target}.fragments.tsv.gz"
        tag = frag_dir / f"hepatocyte.{target}.tn5.tagAlign.gz"
        concatenate_gzip_members([Path(str(r["fragment_part"])) for r in selected], frag)
        concatenate_gzip_members([Path(str(r["tag_part"])) for r in selected], tag)
        summaries[target] = {
            "n_donors": len(selected),
            "fragment_rows": sum(int(r["retained_fragment_rows"]) for r in selected),
            "tn5_tags": 2 * sum(int(r["retained_fragment_rows"]) for r in selected),
            "fragment_multiplicity": sum(int(r["retained_fragment_multiplicity"]) for r in selected),
        }
        for kind, path in (("fragments", frag), ("tn5_tagalign", tag)):
            products.append({"sample": target, "kind": kind, "path": str(path), "bytes": path.stat().st_size, "md5": md5(path)})
    ratio = summaries["pseudorep1"]["fragment_rows"] / summaries["pseudorep2"]["fragment_rows"]
    if not 0.5 <= ratio <= 2.0:
        raise SystemExit(f"pseudoreplicate fragment imbalance is excessive: ratio={ratio:.3f}")

    donor_fields = [x for x in results[0] if x not in {"fragment_part", "tag_part"}]
    with (out / "fragment_qc_by_donor.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=donor_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    with (out / "fragment_products.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(products[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(products)
    summary = {
        "parallel_workers": min(args.workers, len(donors)),
        "gzip_layout": "concatenated donor-level gzip members in deterministic donor-id order",
        "pools": summaries,
        "pseudorep1_to_pseudorep2_fragment_ratio": ratio,
    }
    with (out / "fragment_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    shutil.rmtree(part_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
