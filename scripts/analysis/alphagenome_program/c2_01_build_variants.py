#!/usr/bin/env python3
"""C2 step 1: build the Currin caQTL peak-lead label table and its 4,096-bp windows.

Label orientation: `lead_variant_ID` is chr:pos:REF:ALT and the source FastQTL `beta` is the
ALT-dosage slope, so beta_alt = beta unchanged.  The file's EA column annotates the
accessibility-increasing allele, not the beta allele; re-orienting to EA would flip every row with
beta < 0.  The relation is verified here on the whole file before anything is written.

Outcome-blind downstream: this script writes the label, but the embedding extractor never reads it.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LEAD = ROOT / (
    "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/"
    "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz"
)
FOLDS = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv"
FASTA = Path(
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)
WINDOW = 4096
HALF = WINDOW // 2
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
AMBIGUOUS = ({"A", "T"}, {"C", "G"})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_leads() -> list[dict[str, str]]:
    with gzip.open(LEAD, "rt", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def verify_orientation(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = {"ea_is_alt_beta_positive": 0, "ea_is_ref_beta_negative": 0, "inconsistent": 0, "beta_zero": 0}
    for row in rows:
        chrom, pos, ref, alt = row["lead_variant_ID"].split(":", 3)
        beta = float(row["beta"])
        ea, nea = row["EA"].upper(), row["NEA"].upper()
        if {ea, nea} != {ref.upper(), alt.upper()}:
            counts["inconsistent"] += 1
        elif beta == 0:
            counts["beta_zero"] += 1
        elif beta > 0 and ea == alt.upper():
            counts["ea_is_alt_beta_positive"] += 1
        elif beta < 0 and ea == ref.upper():
            counts["ea_is_ref_beta_negative"] += 1
        else:
            counts["inconsistent"] += 1
    if counts["inconsistent"] or counts["beta_zero"]:
        raise SystemExit(f"Currin EA/beta orientation contract differs: {counts}")
    return counts


def fold_map() -> dict[str, int]:
    mapping: dict[str, int] = {}
    with FOLDS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for chrom in row["test_chromosomes"].split(","):
                if chrom in mapping:
                    raise SystemExit(f"chromosome assigned twice: {chrom}")
                mapping[chrom] = int(row["fold"])
    if set(mapping) != set(AUTOSOMES):
        raise SystemExit("fold manifest does not cover chr1-chr22 exactly once")
    return mapping


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    out = arguments.output
    out.mkdir(parents=True, exist_ok=True)

    import pysam

    rows = read_leads()
    orientation = verify_orientation(rows)
    census: dict[str, int] = {"source_lead_rows": len(rows)}

    frame = pd.DataFrame(
        {
            "peak_id": [row["peak_ID"] for row in rows],
            "lead_variant_id": [row["lead_variant_ID"] for row in rows],
            "beta_source": [float(row["beta"]) for row in rows],
            "q_val": [float(row["q_val"]) for row in rows],
            "p_nominal": [float(row["p_val_nominal"]) for row in rows],
            "maf": [float(row["all_MAF"]) for row in rows],
            "imputation_r2": [float(row["lead_imputation_r2"]) for row in rows],
            "ea": [row["EA"].upper() for row in rows],
            "peak_start_hg38": [int(row["peak_start_hg38"]) for row in rows],
            "peak_stop_hg38": [int(row["peak_stop_hg38"]) for row in rows],
        }
    )
    parts = frame["lead_variant_id"].str.split(":", n=3, expand=True)
    frame["chr"] = parts[0]
    frame["pos_hg38"] = parts[1].astype(int)
    frame["ref"] = parts[2].str.upper()
    frame["alt"] = parts[3].str.upper()
    frame["beta_alt"] = frame["beta_source"]

    frame = frame.loc[frame["q_val"] < 0.05].copy()
    census["after_q_lt_0.05"] = len(frame)
    frame = frame.loc[(frame["ref"].str.len() == 1) & (frame["alt"].str.len() == 1)].copy()
    census["after_snv_only"] = len(frame)
    frame = frame.loc[frame["chr"].isin(AUTOSOMES)].copy()
    census["after_autosomes"] = len(frame)

    genome = pysam.FastaFile(str(FASTA))
    lengths = dict(zip(genome.references, genome.lengths))
    keep_ref, keep_window = [], []
    for row in frame.itertuples(index=False):
        pos0 = int(row.pos_hg38) - 1
        base = genome.fetch(row.chr, pos0, pos0 + 1).upper()
        keep_ref.append(base == row.ref)
        start, end = pos0 - HALF, pos0 + HALF
        if start < 0 or end > lengths[row.chr] or base != row.ref:
            keep_window.append(False)
            continue
        seq = genome.fetch(row.chr, start, end).upper()
        keep_window.append(len(seq) == WINDOW and not (set(seq) - set("ACGT")))
    frame["reference_match"] = keep_ref
    census["reference_mismatch_dropped"] = int((~frame["reference_match"]).sum())
    frame["window_ok"] = keep_window
    frame = frame.loc[frame["reference_match"] & frame["window_ok"]].copy()
    census["after_reference_and_window"] = len(frame)

    before = len(frame)
    frame = frame.sort_values(["lead_variant_id", "q_val"]).drop_duplicates("lead_variant_id", keep="first")
    census["duplicate_variant_rows_removed"] = before - len(frame)
    census["unique_variants"] = len(frame)

    frame["strand_ambiguous"] = [
        {r, a} in AMBIGUOUS for r, a in zip(frame["ref"], frame["alt"])
    ]
    census["strand_ambiguous"] = int(frame["strand_ambiguous"].sum())

    mapping = fold_map()
    frame["heldout_fold"] = frame["chr"].map(mapping)
    frame["block_1mb"] = frame["chr"] + ":" + (frame["pos_hg38"] // 1_000_000).astype(str)
    frame["variant_index0"] = HALF
    frame["window_start0"] = frame["pos_hg38"] - 1 - HALF
    frame["window_end0"] = frame["pos_hg38"] - 1 + HALF
    frame = frame.sort_values(["chr", "pos_hg38"]).reset_index(drop=True)
    frame["row_index"] = np.arange(len(frame))

    census["blocks_1mb"] = int(frame["block_1mb"].nunique())
    census["peaks"] = int(frame["peak_id"].nunique())
    census["beta_alt_positive"] = int((frame["beta_alt"] > 0).sum())
    census["beta_alt_negative"] = int((frame["beta_alt"] < 0).sum())
    per_fold = {
        str(fold): {
            "variants": int((frame["heldout_fold"] == fold).sum()),
            "blocks": int(frame.loc[frame["heldout_fold"] == fold, "block_1mb"].nunique()),
        }
        for fold in range(5)
    }

    frame.to_csv(out / "currin_lead_labels.tsv.gz", sep="\t", index=False, compression="gzip")
    # Outcome-free manifest: the embedding extractor reads only this file plus the FASTA, so no
    # measured allele effect is visible while features are produced.
    manifest_columns = [
        "row_index",
        "lead_variant_id",
        "chr",
        "pos_hg38",
        "ref",
        "alt",
        "window_start0",
        "window_end0",
        "variant_index0",
    ]
    frame[manifest_columns].to_csv(
        out / "window_manifest.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    # Token-encoded REF windows so the GPU step needs only numpy + torch: no FASTA, no pandas,
    # no label column.  Token ids are the frozen probes' character map A7 C8 G9 T10.
    token_ids = {"A": 7, "C": 8, "G": 9, "T": 10}
    tokens = np.empty((len(frame), WINDOW), dtype=np.uint8)
    table = np.zeros(256, dtype=np.uint8)
    for base, value in token_ids.items():
        table[ord(base)] = value
    for position, row in enumerate(frame.itertuples(index=False)):
        seq = genome.fetch(row.chr, int(row.window_start0), int(row.window_end0)).upper()
        tokens[position] = table[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    if int(tokens.min()) < 7 or int(tokens.max()) > 10:
        raise SystemExit("token encoding produced a non-ACGT base")
    centre = tokens[:, HALF]
    expected = np.asarray([token_ids[base] for base in frame["ref"]], dtype=np.uint8)
    if not np.array_equal(centre, expected):
        raise SystemExit("token window centre does not carry the reference allele")
    np.save(out / "window_tokens.npy", tokens)
    np.savez_compressed(
        out / "window_meta.npz",
        row_index=frame["row_index"].to_numpy(np.int64),
        variant_id=frame["lead_variant_id"].to_numpy(dtype=str),
        ref_token=expected,
        alt_token=np.asarray([token_ids[base] for base in frame["alt"]], dtype=np.uint8),
        variant_index0=frame["variant_index0"].to_numpy(np.int64),
    )

    contract = {
        "schema_version": "agp-c2-currin-lead-label-v1",
        "source": str(LEAD.relative_to(ROOT)),
        "source_sha256": sha256(LEAD),
        "fold_manifest": str(FOLDS.relative_to(ROOT)),
        "fold_manifest_sha256": sha256(FOLDS),
        "fasta": str(FASTA),
        "genome_build": "GRCh38 (cellranger-arc GRCh38-2024-A)",
        "label": "beta_alt = source FastQTL beta, ALT-dosage slope; positive means ALT increases accessibility",
        "ea_column_is_not_the_beta_allele": True,
        "orientation_check_counts": orientation,
        "window_bp": WINDOW,
        "variant_index0": HALF,
        "pool_start0": 1280,
        "pool_end0": 2816,
        "census": census,
        "per_fold": per_fold,
        "resampling_unit": "1-Mb block chr:floor(pos/1e6)",
        "inferential_unit": "peak-lead variant",
        "outcome_read_by_extractor": False,
    }
    (out / "label_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
