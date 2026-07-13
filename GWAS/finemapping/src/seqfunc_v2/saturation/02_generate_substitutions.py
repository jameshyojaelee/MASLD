#!/usr/bin/env python3
"""Generate all three substitutions in 501-bp anchor windows.

Each unique genomic substitution is scored once with the exact model that held
its chromosome out. Overlapping anchor windows are restored through an explicit
target-to-genomic map; no target context is silently discarded.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
from collections import defaultdict
from pathlib import Path

from common import CHROM_TO_FOLD, OUT, atomic_json, sha256

BASES = "ACGT"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", type=Path, default=OUT / "anchors/anchor_manifest.tsv")
    ap.add_argument("--fasta", type=Path, default=Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"))
    ap.add_argument("--out-dir", type=Path, default=OUT / "sequences")
    ap.add_argument("--window", type=int, default=501)
    ap.add_argument("--shard-size", type=int, default=25000)
    args = ap.parse_args()
    if args.window != 501 or args.window % 2 != 1:
        raise SystemExit("approved saturation contract requires exactly 501 bp")
    gate = OUT / "gates/upstream_gate.json"
    if not gate.is_file():
        raise SystemExit(f"missing PASS upstream gate: {gate}")
    from common import read_json
    if not read_json(gate).get("saturation_submission_allowed", False):
        raise SystemExit("upstream gate BLOCKED")
    with (OUT / "gates/model_manifest.locked.tsv").open() as handle:
        locked = list(csv.DictReader(handle, delimiter="\t"))
    locked_fastas = {Path(x["reference_fasta"]).resolve() for x in locked}
    if locked_fastas != {args.fasta.resolve()}:
        raise SystemExit(f"sequence FASTA differs from locked model FASTA: {args.fasta} versus {locked_fastas}")
    import pysam
    fa = pysam.FastaFile(str(args.fasta)); half = args.window // 2
    with args.anchors.open() as handle:
        anchors = list(csv.DictReader(handle, delimiter="\t"))
    if not 1 <= len(anchors) <= 500:
        raise SystemExit(f"anchor count outside approved 1..500: {len(anchors)}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "GRCh38.primary.chrom.sizes.tsv").open("w") as handle:
        for chrom in sorted(CHROM_TO_FOLD, key=lambda x: int(x[3:])):
            handle.write(f"{chrom}\t{fa.get_reference_length(chrom)}\n")
    locked_sizes = {Path(x["chrom_sizes_path"]).resolve() for x in locked}
    if len(locked_sizes) != 1 or (args.out_dir / "GRCh38.primary.chrom.sizes.tsv").read_text() != next(iter(locked_sizes)).read_text():
        raise SystemExit("FASTA-derived chromosome sizes differ from the locked producer chromosome sizes")
    target_map, genomic, rejected = [], {}, []
    for row in anchors:
        rank = int(row["anchor_rank"]); chrom = row["chr_hg38"]
        pos = int(row["pos_hg38"]); ref = row["ref_hg38"].upper()
        if chrom not in CHROM_TO_FOLD:
            rejected.append({"anchor_rank": rank, "reason": "no_frozen_fold", "detail": chrom}); continue
        start = pos - 1 - half; end = start + args.window
        if start < 0 or end > fa.get_reference_length(chrom):
            rejected.append({"anchor_rank": rank, "reason": "window_outside_contig", "detail": f"{start}:{end}"}); continue
        seq = fa.fetch(chrom, start, end).upper()
        if len(seq) != args.window or seq[half] != ref or any(x not in BASES for x in seq):
            rejected.append({"anchor_rank": rank, "reason": "sequence_reference_or_non_acgt", "detail": seq[half] if seq else ""}); continue
        target_id = f"sat{rank:04d}_{chrom}_{pos}_{ref}_{row['alt_hg38']}"
        for offset, base in enumerate(seq):
            mut_pos = start + offset + 1
            for alt in BASES:
                if alt == base: continue
                gid = f"{chrom}:{mut_pos}:{base}:{alt}"
                tid = f"{target_id}__{gid}"
                key = (chrom, mut_pos, base, alt, CHROM_TO_FOLD[chrom])
                if gid in genomic and genomic[gid] != key:
                    raise RuntimeError(f"inconsistent genomic identity: {gid}")
                genomic[gid] = key
                target_map.append({
                    "target_perturbation_id": tid, "genomic_variant_id": gid,
                    "target_id": target_id, "anchor_rank": rank,
                    "anchor_variant_id_hg38": row["variant_id_hg38"],
                    "chrom": chrom, "pos_hg38": mut_pos, "ref": base, "alt": alt,
                    "offset_from_anchor": mut_pos - pos,
                    "is_observed_anchor_alt": str(mut_pos == pos and alt == row["alt_hg38"]).upper(),
                    "model_fold_id": CHROM_TO_FOLD[chrom], "uniform_max_pip": row["uniform_max_pip"],
                })
    if rejected:
        with (args.out_dir / "anchor_rejections.tsv").open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(rejected[0]), delimiter="\t"); w.writeheader(); w.writerows(rejected)
        raise SystemExit(f"hard sequence gate failed for {len(rejected)} anchors")
    expected = len(anchors) * args.window * 3
    if len(target_map) != expected:
        raise SystemExit(f"substitution count mismatch: {len(target_map)} != {expected}")

    map_path = args.out_dir / "target_variant_map.tsv.gz"
    tmp = map_path.with_suffix(".tmp.gz")
    with gzip.open(tmp, "wt", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(target_map[0]), delimiter="\t"); w.writeheader(); w.writerows(target_map)
    os.replace(tmp, map_path)

    ordered = sorted(genomic.items(), key=lambda z: (int(z[1][0][3:]), z[1][1], z[1][2], z[1][3]))
    by_fold = defaultdict(list)
    for gid, key in ordered: by_fold[key[4]].append((gid, key))
    shard_rows = []
    for fold in sorted(by_fold):
        values = by_fold[fold]
        for start in range(0, len(values), args.shard_size):
            shard_id = f"{fold}_shard{start // args.shard_size:03d}"
            tab = args.out_dir / f"{shard_id}.variants.tsv"
            cbp = args.out_dir / f"{shard_id}.chrombpnet.tsv"
            with tab.open("w", newline="") as ht, cbp.open("w") as hc:
                w = csv.writer(ht, delimiter="\t"); w.writerow(["chrom","pos_hg38","ref","alt","variant_id","model_fold_id"])
                for gid, (chrom, pos, ref, alt, observed_fold) in values[start:start + args.shard_size]:
                    if observed_fold != fold: raise RuntimeError("fold routing invariant failed")
                    w.writerow([chrom,pos,ref,alt,gid,fold]); hc.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gid}\n")
            shard_rows.append({"array_index": len(shard_rows), "shard_id": shard_id, "fold_id": fold,
                               "test_chromosomes": ";".join(sorted({x[1][0] for x in values[start:start + args.shard_size]})),
                               "n_variants": len(values[start:start + args.shard_size]),
                               "variant_table": str(tab.resolve()), "scorer_list": str(cbp.resolve())})
    with (args.out_dir / "scoring_shards.tsv").open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(shard_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(shard_rows)
    atomic_json({
        "status": "PASS", "anchors": len(anchors), "sequence_window": args.window,
        "target_substitution_rows": len(target_map), "maximum_allowed": 751500,
        "unique_genomic_substitutions": len(genomic), "scoring_shards": len(shard_rows),
        "overlapping_genomic_substitutions_scored_once": True,
        "routing": "exact frozen chromosome-heldout fold",
        "anchor_manifest_sha256": sha256(args.anchors), "fasta": str(args.fasta),
        "chrom_sizes": str((args.out_dir / "GRCh38.primary.chrom.sizes.tsv").resolve()),
        "formal_fdr_allowed": False, "apply_only_firewall": True,
    }, args.out_dir / "sequence_contract.json")
    print(f"PASS: {len(target_map)} target substitutions, {len(genomic)} unique, {len(shard_rows)} shards")


if __name__ == "__main__":
    main()
