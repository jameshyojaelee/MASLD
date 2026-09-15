#!/usr/bin/env python3
"""D1 task 4: the mapping-bias control the existing allelic counts lack (WASP-style swap and realign).

The two cohorts were aligned with different aligners, so the control is run per cohort with that
cohort's aligner and settings:
  GSE281367  cellranger-atac, BWA-MEM, cellranger-arc GRCh38-2024-A index
  GSE244832  bowtie2 2.5.4 --very-sensitive, the GRCh38 index named in the BAM @PG line

Nested comparator. The swapped read is compared with the SAME read realigned unswapped through the
same command, not with its original BAM record: the originals were aligned as pairs and this
realignment is single-end, so comparing against the BAM would confound the base swap with the
change in alignment mode. The two arms differ in exactly one base.

A read fails if the swapped and unswapped realignments differ in (contig, position, strand), or if
they fall on opposite sides of the MAPQ 30 threshold the pileup uses (-q 30), or if either is
unmapped. Primary site filter: drop the site if ANY evaluated read fails. Secondary: drop the site
if more than 1 percent of evaluated reads fail. Both are reported.

Writes into <out>/tables:
  mapping_bias_sites.tsv      per site: reads evaluated, reads failing, both filter verdicts
  mapping_bias_effect.tsv     site-level imbalance recomputed with and without the filter
  mapping_bias.json           headline counts
"""
from __future__ import annotations

import json
import pathlib
import random
import subprocess
import sys

import numpy as np
import pandas as pd
import pysam

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
P3 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z"
TARGETS = P3 / "tables/atac_targets_in_peaks.tsv"
BWA_INDEX = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
BT2_INDEX = "/gpfs/commons/home/jameslee/reference_genome/bowtie2/GRCh38/GRCh38"
CAP = 300
SEED = 20260914
MAPQ_CLASS = 30
ROTATE = {"A": "C", "C": "G", "G": "T", "T": "A", "N": "A"}
COHORT_DONORS = {"GSE281367": [f"Z{i:02d}" for i in range(1, 13)],
                 "GSE244832": [f"D{i:02d}" for i in range(1, 19)]}


def load_targets() -> dict:
    tgt = {}
    with open(TARGETS) as h:
        for line in h:
            c, p, ref, alt = line.rstrip("\n").split("\t")
            tgt[(c, int(p))] = {"ref": ref, "alt": alt, "uid": f"{c}:{p}:{ref}:{alt}"}
    return tgt


def collect_reads(cohort: str, subdir: pathlib.Path, targets: dict) -> dict:
    """Reservoir-sample up to CAP reads per site across the cohort's donors.

    A read is a candidate only if it carries the reference or the target alternate base at the site
    (a third allele, a deletion or a skip there has no swap defined)."""
    by_chrom: dict[str, list[int]] = {}
    for (c, p) in targets:
        by_chrom.setdefault(c, []).append(p)
    for c in by_chrom:
        by_chrom[c].sort()
    arrays = {c: np.asarray(v) for c, v in by_chrom.items()}

    res: dict[str, list] = {}
    seen: dict[str, int] = {}
    rngs: dict[str, random.Random] = {}
    for donor in COHORT_DONORS[cohort]:
        bam = subdir / f"{donor}.bam"
        if not bam.exists():
            continue
        with pysam.AlignmentFile(str(bam), "rb") as fh:
            for read in fh.fetch(until_eof=False):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate \
                        or read.is_qcfail or read.mapping_quality < MAPQ_CLASS:
                    continue
                chrom = read.reference_name
                arr = arrays.get(chrom)
                if arr is None:
                    continue
                lo = np.searchsorted(arr, read.reference_start + 1, "left")
                hi = np.searchsorted(arr, read.reference_end, "right")
                if hi <= lo:
                    continue
                pairs = None
                seq = read.query_sequence
                if not seq:
                    continue
                qual = pysam.qualities_to_qualitystring(read.query_qualities) or ("I" * len(seq))
                for pos in arr[lo:hi]:
                    t = targets[(chrom, int(pos))]
                    uid = t["uid"]
                    if pairs is None:
                        pairs = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
                    q = pairs.get(int(pos) - 1)
                    if q is None:
                        continue
                    base = seq[q].upper()
                    if base not in (t["ref"], t["alt"]):
                        continue
                    swapped = t["alt"] if base == t["ref"] else t["ref"]
                    # matched non-allelic perturbation (amendment 01): one base changed at an in-read
                    # offset at least 10 bases from the target, so the control arm differs from the
                    # allele arm in one component only - whether the changed base is the allele.
                    qc = q + 15 if q + 15 < len(seq) else q - 15
                    if 0 <= qc < len(seq):
                        ctrl = seq[:qc] + ROTATE.get(seq[qc].upper(), "A") + seq[qc + 1:]
                    else:
                        ctrl = None
                    rec = (seq[:q] + swapped + seq[q + 1:], seq, qual, ctrl)
                    n = seen.get(uid, 0) + 1
                    seen[uid] = n
                    bucket = res.setdefault(uid, [])
                    if len(bucket) < CAP:
                        bucket.append(rec)
                    else:
                        rng = rngs.setdefault(uid, random.Random(f"{SEED}:{uid}"))
                        j = rng.randrange(n)
                        if j < CAP:
                            bucket[j] = rec
    return res


def write_fastq(res: dict, path_swap: pathlib.Path, path_orig: pathlib.Path,
                path_ctrl: pathlib.Path) -> int:
    n = 0
    with path_swap.open("w") as fs, path_orig.open("w") as fo, path_ctrl.open("w") as fc:
        for uid, bucket in sorted(res.items()):
            for i, (swap, orig, qual, ctrl) in enumerate(bucket):
                name = f"{uid}|{i}"
                fs.write(f"@{name}\n{swap}\n+\n{qual}\n")
                fo.write(f"@{name}\n{orig}\n+\n{qual}\n")
                if ctrl is not None:
                    fc.write(f"@{name}\n{ctrl}\n+\n{qual}\n")
                n += 1
    return n


def realign(cohort: str, fq: pathlib.Path, sam: pathlib.Path, threads: int) -> None:
    if cohort == "GSE281367":
        cmd = ["bwa", "mem", "-t", str(threads), BWA_INDEX, str(fq)]
    else:
        cmd = ["bowtie2", "--very-sensitive", "-p", str(threads), "-x", BT2_INDEX, "-U", str(fq)]
    with sam.open("w") as out, open(str(sam) + ".log", "w") as err:
        subprocess.run(cmd, stdout=out, stderr=err, check=True)


def parse_sam(sam: pathlib.Path) -> dict:
    out = {}
    with pysam.AlignmentFile(str(sam), "r") as fh:
        for r in fh:
            if r.is_secondary or r.is_supplementary:
                continue
            if r.query_name in out:
                continue
            out[r.query_name] = (None if r.is_unmapped else r.reference_name,
                                 -1 if r.is_unmapped else r.reference_start,
                                 bool(r.is_reverse), int(r.mapping_quality), bool(r.is_unmapped))
    return out


def main() -> None:
    out = pathlib.Path(sys.argv[1]).resolve()
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    tables, work = out / "tables", out / "raw" / "bias"
    work.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    subdir = out / "raw/subset_bam"

    rows = []
    readlen: dict[str, float] = {}
    for cohort in COHORT_DONORS:
        res = collect_reads(cohort, subdir, targets)
        fq_s, fq_o = work / f"{cohort}.swap.fq", work / f"{cohort}.orig.fq"
        fq_c = work / f"{cohort}.ctrl.fq"
        n = write_fastq(res, fq_s, fq_o, fq_c)
        print(f"[d1_04] {cohort}: {len(res)} sites, {n} reads (cap {CAP}/site)", flush=True)
        sam_s, sam_o = work / f"{cohort}.swap.sam", work / f"{cohort}.orig.sam"
        sam_c = work / f"{cohort}.ctrl.sam"
        for fq, sam in ((fq_o, sam_o), (fq_s, sam_s), (fq_c, sam_c)):
            realign(cohort, fq, sam, threads)
        A, B, C = parse_sam(sam_o), parse_sam(sam_s), parse_sam(sam_c)

        def failed(a, b) -> bool:
            if b is None:
                return True
            same_place = (a[0], a[1], a[2]) == (b[0], b[1], b[2])
            same_class = (a[3] >= MAPQ_CLASS) == (b[3] >= MAPQ_CLASS)
            return bool(a[4] or b[4] or not same_place or not same_class)

        per = {}
        for name, a in A.items():
            uid = name.split("|")[0]
            d = per.setdefault(uid, [0, 0, 0, 0])
            d[0] += 1
            d[1] += int(failed(a, B.get(name)))
            c = C.get(name)
            if c is not None:
                d[2] += 1
                d[3] += int(failed(a, c))
        readlen[cohort] = float(np.median([len(r[1]) for b in res.values() for r in b])) if res else float("nan")
        for uid, (n_eval, n_fail, n_ctrl, n_ctrl_fail) in per.items():
            rows.append({"cohort": cohort, "uid": uid, "n_reads_evaluated": n_eval,
                         "n_reads_failing": n_fail, "fail_fraction": n_fail / n_eval,
                         "n_control_reads": n_ctrl, "n_control_failing": n_ctrl_fail,
                         "control_fail_fraction": (n_ctrl_fail / n_ctrl) if n_ctrl else float("nan"),
                         "pass_primary_any_read": n_fail == 0,
                         "pass_secondary_1pct": (n_fail / n_eval) <= 0.01,
                         "control_pass_primary_any_read": n_ctrl_fail == 0,
                         "control_pass_secondary_1pct": ((n_ctrl_fail / n_ctrl) <= 0.01) if n_ctrl else False})
        for f in (fq_s, fq_o, fq_c, sam_s, sam_o, sam_c):
            f.unlink(missing_ok=True)
            pathlib.Path(str(f) + ".log").unlink(missing_ok=True)
    bias = pd.DataFrame(rows)
    bias.to_csv(tables / "mapping_bias_sites.tsv", sep="\t", index=False)

    # recompute the site-level imbalance with and without the filter, through the same loader
    sites = pd.read_csv(tables / "lineage_site_summary.tsv", sep="\t")
    eff = []
    for cohort in COHORT_DONORS:
        b = bias[bias["cohort"] == cohort].set_index("uid")
        for lineage in ("ALL_READS", "Hepatocyte", "ALL_LABELLED"):
            s = sites[(sites["scope"] == cohort) & (sites["lineage"] == lineage)].copy()
            if s.empty:
                continue
            s["n_reads_evaluated"] = s["uid"].map(b["n_reads_evaluated"])
            s["pass_primary"] = s["uid"].map(b["pass_primary_any_read"])
            s["pass_secondary"] = s["uid"].map(b["pass_secondary_1pct"])
            for label, mask in (("unfiltered", pd.Series(True, index=s.index)),
                                ("wasp_primary", s["pass_primary"].fillna(False)),
                                ("wasp_secondary", s["pass_secondary"].fillna(False))):
                g = s[mask]
                eff.append({"cohort": cohort, "lineage": lineage, "filter": label,
                            "n_sites": int(len(g)),
                            "mean_log2_alt_over_ref": float(g["mean_log2_alt_over_ref"].mean()) if len(g) else float("nan"),
                            "mean_abs_log2": float(g["mean_log2_alt_over_ref"].abs().mean()) if len(g) else float("nan"),
                            "frac_sites_alt_higher": float((g["mean_log2_alt_over_ref"] > 0).mean()) if len(g) else float("nan"),
                            "n_blocks": int(g["block"].nunique()) if len(g) else 0})
    effdf = pd.DataFrame(eff)
    effdf.to_csv(tables / "mapping_bias_effect.tsv", sep="\t", index=False)

    summary = {"cap_reads_per_site": CAP, "seed": SEED, "mapq_class_boundary": MAPQ_CLASS,
               "aligners": {"GSE281367": "bwa mem (cellranger-arc GRCh38-2024-A)",
                            "GSE244832": "bowtie2 --very-sensitive (cohort GRCh38 index)"},
               "control_arm": "amendment 01: one base changed at an in-read offset >= 10 from the "
                              "target (q+15, else q-15), A->C C->G G->T T->A; realigned and compared "
                              "by the same rule. Its failure rate is the floor any single-base change "
                              "imposes at this read length and aligner.",
               "per_cohort": []}
    for cohort in COHORT_DONORS:
        b = bias[bias["cohort"] == cohort]
        summary["per_cohort"].append({
            "cohort": cohort, "sites_evaluated": int(len(b)),
            "reads_evaluated": int(b["n_reads_evaluated"].sum()),
            "reads_failing": int(b["n_reads_failing"].sum()),
            "read_fail_rate": float(b["n_reads_failing"].sum() / max(1, b["n_reads_evaluated"].sum())),
            "sites_removed_primary": int((~b["pass_primary_any_read"]).sum()),
            "sites_removed_secondary": int((~b["pass_secondary_1pct"]).sum()),
            "median_read_length": readlen.get(cohort, float("nan")),
            "control_reads_evaluated": int(b["n_control_reads"].sum()),
            "control_read_fail_rate": float(b["n_control_failing"].sum() / max(1, b["n_control_reads"].sum())),
            "control_sites_removed_primary": int((~b["control_pass_primary_any_read"]).sum()),
            "control_sites_removed_secondary": int((~b["control_pass_secondary_1pct"]).sum()),
            "allele_minus_control_read_fail_rate": float(
                b["n_reads_failing"].sum() / max(1, b["n_reads_evaluated"].sum())
                - b["n_control_failing"].sum() / max(1, b["n_control_reads"].sum())),
        })
    json.dump(summary, (tables / "mapping_bias.json").open("w"), indent=1, default=float)
    print(json.dumps(summary, indent=1))
    print(effdf.to_string(index=False))


if __name__ == "__main__":
    main()
