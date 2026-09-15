#!/usr/bin/env python3
"""E2 step 1: windows and labels for the reporter-to-endogenous transfer test.

Three deposits are written:

* `window_tokens.npy` / `window_meta.npz` - 4,096-bp token windows for all 4,359 reconstructed
  GSE281364 paired SNVs, variant at 0-based index 2048, the published reporter geometry
  (`gse281364-hyenadna-window-ablation-20260905T192600Z`, arm `w4096_proportional`).  These two
  files carry no label, so the GPU extractor never opens an outcome.
* `reporter_labels.tsv.gz` - the per-element MPRA allele effect
  `d = log2((RNA+0.5)/(DNA+0.5))_alt - log2((RNA+0.5)/(DNA+0.5))_ref`, averaged over the four
  experimental replicates within a context, per context and in the prespecified combinations.
* `endogenous_labels.tsv.gz` - the Currin nominal caQTL beta from the E1 bridge, nearest-peak
  (primary) and min-p-peak (outcome-selected sensitivity), already oriented to the MPRA ALT allele
  by E1.  The multiplier is NOT re-applied.

Fold assignment is the ChromBPNet `hepatocyte_5fold_v2` chromosome grouping, so every held-out fold
is a set of whole chromosomes and no 1-Mb block straddles a fold boundary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
ELEMENTS = BENCH / "executions/gse281364-outcome-blind-splits-21068885/split/elements.tsv"
OUTCOMES = BENCH / "executions/gse281364-validated-reconstruction-21066470/validated/replicate_outcomes.tsv.gz"
GROUPMAP = BENCH / "executions/gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv"
BRIDGE = ROOT / "GWAS/finemapping/results/alphagenome_program/e1-bridge-20260914T194145Z/bridge_variants.tsv"
FOLDS = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv"
FASTA = Path(
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)

WINDOW = 4096
HALF = WINDOW // 2
PSEUDOCOUNT = 0.5
CONTEXTS = ("HepG2_control", "HepG2_PAOA", "LX2_control", "LX2_TGFb")
HEPG2 = ("HepG2_control", "HepG2_PAOA")
REPLICATES = (1, 2, 3, 4)
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
AMBIGUOUS = ({"A", "T"}, {"C", "G"})
TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10}


class BuildError(RuntimeError):
    """Raised when an input contract differs from what this package requires."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fold_map() -> dict[str, int]:
    mapping: dict[str, int] = {}
    with FOLDS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for chrom in row["test_chromosomes"].split(","):
                if chrom in mapping:
                    raise BuildError(f"chromosome assigned twice: {chrom}")
                mapping[chrom] = int(row["fold"])
    if set(mapping) != set(AUTOSOMES):
        raise BuildError("fold manifest does not cover chr1-chr22 exactly once")
    return mapping


def reporter_labels(elements: set[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    """Per-element, per-context allele effect from the validated replicate table."""
    counts = {"rows_read": 0, "rows_kept": 0, "below_qc_rows": 0}
    pairs: dict[tuple[str, str, int], dict[str, tuple[int, int]]] = defaultdict(dict)
    with gzip.open(OUTCOMES, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts["rows_read"] += 1
            if row["element_id"] not in elements or row["context_id"] not in CONTEXTS:
                continue
            if row["assay_state"] != "observed":
                counts["below_qc_rows"] += 1
                continue
            if (
                row["biological_unit"] != "experimental_replicate"
                or row["pairing"] != "same_sample_different_aliquot"
                or row["donor_id"] != "not_applicable"
                or row["allele"] not in {"ref", "alt"}
            ):
                raise BuildError("reporter outcome topology differs")
            key = (row["element_id"], row["context_id"], int(row["experimental_replicate"]))
            if row["allele"] in pairs[key]:
                raise BuildError("duplicate observed reporter allele")
            pairs[key][row["allele"]] = (int(row["DNA"]), int(row["RNA"]))
            counts["rows_kept"] += 1

    records = []
    for element in sorted(elements):
        record: dict[str, object] = {"element_id": element}
        for context in CONTEXTS:
            values = []
            for replicate in REPLICATES:
                pair = pairs.get((element, context, replicate), {})
                if set(pair) != {"ref", "alt"}:
                    values = []
                    break
                dna_r, rna_r = pair["ref"]
                dna_a, rna_a = pair["alt"]
                values.append(
                    np.log2((rna_a + PSEUDOCOUNT) / (dna_a + PSEUDOCOUNT))
                    - np.log2((rna_r + PSEUDOCOUNT) / (dna_r + PSEUDOCOUNT))
                )
            record[f"d_{context}"] = float(np.mean(values)) if values else np.nan
            record[f"complete_{context}"] = bool(values)
        records.append(record)
    frame = pd.DataFrame.from_records(records)
    frame["d_hepg2_mean"] = frame[[f"d_{c}" for c in HEPG2]].mean(axis=1, skipna=False)
    frame["d_all_context_mean"] = frame[[f"d_{c}" for c in CONTEXTS]].mean(axis=1, skipna=False)
    return frame, counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    out = arguments.output
    out.mkdir(parents=True, exist_ok=True)

    import pysam

    elements = pd.read_csv(ELEMENTS, sep="\t")
    if len(elements) != 4359:
        raise BuildError(f"element census differs: {len(elements)}")
    elements = elements.rename(columns={"genomic_ref": "mpra_ref", "genomic_alt": "mpra_alt"})
    if not ((elements.mpra_ref.str.len() == 1) & (elements.mpra_alt.str.len() == 1)).all():
        raise BuildError("non-SNV element present")
    if not elements.contig.isin(AUTOSOMES).all():
        raise BuildError("non-autosomal element present")

    group = pd.read_csv(GROUPMAP, sep="\t")
    elements = elements.merge(
        group[["source_locus_group_id", "borzoi_long_range_group_id"]].rename(
            columns={
                "source_locus_group_id": "outer_locus_sequence_group_id",
                "borzoi_long_range_group_id": "block_lr239",
            }
        ),
        on="outer_locus_sequence_group_id",
        how="left",
        validate="many_to_one",
    )
    if elements.block_lr239.isna().any():
        raise BuildError("long-range block join left gaps")

    mapping = fold_map()
    elements["chrom_fold"] = elements.contig.map(mapping)
    elements["block_1mb"] = elements.contig + "~" + (elements.variant_pos1 // 1_000_000).astype(str)
    elements["strand_ambiguous"] = [
        {r, a} in AMBIGUOUS for r, a in zip(elements.mpra_ref, elements.mpra_alt)
    ]
    elements = elements.sort_values(["contig", "variant_pos1"]).reset_index(drop=True)
    elements["row_index"] = np.arange(len(elements))

    straddle = elements.groupby("block_1mb").chrom_fold.nunique()
    straddle_outer = elements.groupby("block_1mb").outer_fold.nunique()
    if int((straddle > 1).sum()) != 0:
        raise BuildError("a 1-Mb block straddles a chromosome fold")

    # ---------------------------------------------------------------- windows
    genome = pysam.FastaFile(str(FASTA))
    lengths = dict(zip(genome.references, genome.lengths))
    table = np.zeros(256, dtype=np.uint8)
    for base, value in TOKEN_IDS.items():
        table[ord(base)] = value
    # The frozen HyenaDNA token contract is A7 C8 G9 T10 with no N token, so an element whose
    # 4,096-bp window is not pure ACGT cannot be encoded under the published recipe and is dropped.
    # The dropped ids are recorded in the contract, never silently removed.
    sequences: list[np.ndarray] = []
    unencodable: list[dict[str, str]] = []
    keep = []
    for row in elements.itertuples(index=False):
        start, end = int(row.variant_pos1) - 1 - HALF, int(row.variant_pos1) - 1 + HALF
        if start < 0 or end > lengths[row.contig]:
            unencodable.append({"element_id": row.element_id, "reason": "window_off_contig"})
            keep.append(False)
            continue
        sequence = genome.fetch(row.contig, start, end).upper()
        if len(sequence) != WINDOW or set(sequence) - set("ACGT"):
            unencodable.append({"element_id": row.element_id, "reason": "window_not_pure_ACGT"})
            keep.append(False)
            continue
        sequences.append(table[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)])
        keep.append(True)
    elements = elements.loc[np.asarray(keep)].reset_index(drop=True)
    elements["row_index"] = np.arange(len(elements))
    tokens = np.vstack(sequences).astype(np.uint8)
    if tokens.shape != (len(elements), WINDOW):
        raise BuildError("token matrix does not match the kept element census")
    reference_token = np.asarray([TOKEN_IDS[b] for b in elements.mpra_ref], dtype=np.uint8)
    alternate_token = np.asarray([TOKEN_IDS[b] for b in elements.mpra_alt], dtype=np.uint8)
    if not np.array_equal(tokens[:, HALF], reference_token):
        raise BuildError("window centre does not carry the MPRA reference allele")
    np.save(out / "window_tokens.npy", tokens)
    np.savez_compressed(
        out / "window_meta.npz",
        row_index=elements.row_index.to_numpy(np.int64),
        variant_id=(
            elements.contig
            + ":"
            + elements.variant_pos1.astype(str)
            + ":"
            + elements.mpra_ref
            + ":"
            + elements.mpra_alt
        ).to_numpy(dtype=str),
        ref_token=reference_token,
        alt_token=alternate_token,
        variant_index0=np.full(len(elements), HALF, dtype=np.int64),
    )

    # ---------------------------------------------------------------- reporter label
    reporter, reporter_counts = reporter_labels(set(elements.element_id))
    elements = elements.merge(reporter, on="element_id", how="left", validate="one_to_one")

    # ---------------------------------------------------------------- endogenous label
    bridge = pd.read_csv(BRIDGE, sep="\t")
    endogenous = []
    for label, tag in (
        ("currin_nominal_caqtl_nearest_peak", "nearest"),
        ("currin_nominal_caqtl_minp_peak", "minp"),
    ):
        part = bridge.loc[bridge.label == label].copy()
        if part.element_id.duplicated().any():
            raise BuildError(f"{label} is not one row per element")
        if not (part.match_orientation == "direct").all() or not (
            part.orientation_multiplier == 1.0
        ).all():
            raise BuildError(f"{label} orientation is not uniformly direct")
        part[f"beta_{tag}"] = part.label_value_oriented_to_mpra_alt.astype(float)
        part[f"se_{tag}"] = part.label_se.astype(float)
        part[f"p_{tag}"] = part.label_p.astype(float)
        part[f"peak_{tag}"] = part.label_feature.astype(str)
        part[f"distance_{tag}"] = (
            part.label_extra.str.extract(r"distance_from_peakCenter=(-?\d+)")[0].astype(float)
        )
        if part[f"distance_{tag}"].isna().any():
            raise BuildError(f"{label} carries a row without distance_from_peakCenter")
        endogenous.append(
            part[
                [
                    "element_id",
                    f"beta_{tag}",
                    f"se_{tag}",
                    f"p_{tag}",
                    f"peak_{tag}",
                    f"distance_{tag}",
                ]
            ]
        )
    labels = endogenous[0].merge(endogenous[1], on="element_id", how="outer", validate="one_to_one")
    elements = elements.merge(labels, on="element_id", how="left", validate="one_to_one")
    elements["abs_distance_nearest"] = elements.distance_nearest.abs()
    elements["has_endogenous"] = elements.beta_nearest.notna()

    keep = [
        "row_index",
        "element_id",
        "contig",
        "variant_pos1",
        "mpra_ref",
        "mpra_alt",
        "strand_ambiguous",
        "outer_locus_sequence_group_id",
        "outer_fold",
        "chrom_fold",
        "block_1mb",
        "block_lr239",
        "has_endogenous",
        "beta_nearest",
        "se_nearest",
        "p_nearest",
        "peak_nearest",
        "distance_nearest",
        "abs_distance_nearest",
        "beta_minp",
        "se_minp",
        "p_minp",
        "peak_minp",
        "distance_minp",
        *[f"d_{c}" for c in CONTEXTS],
        *[f"complete_{c}" for c in CONTEXTS],
        "d_hepg2_mean",
        "d_all_context_mean",
    ]
    elements[keep].to_csv(out / "e2_units.tsv.gz", sep="\t", index=False, compression="gzip")

    strata = {}
    for threshold in (500, 1000, 5000, 10000, 100000, None):
        mask = elements.has_endogenous & (
            elements.abs_distance_nearest <= (threshold if threshold is not None else np.inf)
        )
        subset = elements.loc[mask]
        per_fold = {
            str(fold): {
                "variants": int((subset.chrom_fold == fold).sum()),
                "blocks": int(subset.loc[subset.chrom_fold == fold, "block_1mb"].nunique()),
            }
            for fold in range(5)
        }
        strata["unlimited" if threshold is None else str(threshold)] = {
            "variants": int(mask.sum()),
            "blocks_1mb": int(subset.block_1mb.nunique()),
            "blocks_lr239": int(subset.block_lr239.nunique()),
            "per_fold": per_fold,
            "min_fold_blocks": min(v["blocks"] for v in per_fold.values()),
            "min_fold_variants": min(v["variants"] for v in per_fold.values()),
        }

    contract = {
        "schema_version": "agp-e2-transfer-inputs-v1",
        "elements": int(len(elements)),
        "elements_dropped_unencodable": unencodable,
        "elements_with_endogenous_label": int(elements.has_endogenous.sum()),
        "reporter_complete_by_context": {
            context: int(elements[f"complete_{context}"].sum()) for context in CONTEXTS
        },
        "reporter_complete_hepg2_both": int(elements.d_hepg2_mean.notna().sum()),
        "reporter_complete_all_four": int(elements.d_all_context_mean.notna().sum()),
        "reporter_row_counts": reporter_counts,
        "blocks_1mb_total": int(elements.block_1mb.nunique()),
        "blocks_1mb_straddling_chrom_fold": int((straddle > 1).sum()),
        "blocks_1mb_straddling_mpra_outer_fold": int((straddle_outer > 1).sum()),
        "per_fold_elements": {
            str(fold): int((elements.chrom_fold == fold).sum()) for fold in range(5)
        },
        "distance_strata": strata,
        "window_bp": WINDOW,
        "variant_index0": HALF,
        "pool_start0": 1280,
        "pool_end0": 2816,
        "genome_build": "GRCh38 (cellranger-arc GRCh38-2024-A)",
        "endogenous_label": (
            "Currin nominal caQTL beta, already oriented to the MPRA genomic_alt allele by E1; "
            "positive = the MPRA ALT allele raises normalised ATAC signal; multiplier not re-applied"
        ),
        "reporter_label": (
            "log2((RNA+0.5)/(DNA+0.5)) alt minus ref, mean over four experimental replicates "
            "within a context; positive = the MPRA ALT allele raises reporter activity"
        ),
        "fold_definition": "ChromBPNet hepatocyte_5fold_v2 chromosome groups",
        "inputs_sha256": {
            str(path): sha256(path)
            for path in (ELEMENTS, OUTCOMES, GROUPMAP, BRIDGE, FOLDS)
        },
        "outcome_read_by_extractor": False,
    }
    (out / "input_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
