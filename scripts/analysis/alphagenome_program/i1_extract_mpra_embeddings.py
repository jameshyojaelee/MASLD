#!/usr/bin/env python3
"""Adaptation arm 2, stage 1: frozen AlphaGenome trunk embeddings for the GSE281364 reporter set.

Prespecified in
`GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/PRESPECIFICATION.md`
(sha256 f44a912f73a1a30cddcfb60aacd5ddd0edb2e49838a0b8c6011c28be3ad267fc) and its addendum 01. Nothing
here reads an MPRA outcome; this stage cannot see a label.

Route to the trunk, from the runtime probe's ADDENDUM section 2: `create()` keeps the trunk private, so the
Orbax restore of `create()` is reproduced here and `create_model(metadata)`'s third return value,
`trunk_apply_fn(params, state, dna_sequence, organism_index) -> Embeddings`, is called directly.
`Embeddings.get_sequence_embeddings(128)` gives the (1, S/128, 3072) representation. The 1-bp
representation is not used.

Two prespecified input lengths (2,048 primary, 16,384 secondary) and two prespecified poolings
(pool-region mean primary, variant-centre bin secondary), all four written to the same npz per length.

Outputs of these LOCAL weights may be used as features. Hosted Atlas/API outputs may not, and the two are
never pooled.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import pathlib
import time

import numpy as np

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BENCH = PROJ / "Analysis/MASLD_Model_Benchmark"
CHECKPOINT = BENCH / "executions/alphagenome-weights-20260915T103442Z/checkpoints"
FASTA_PATH = (
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)
FIXTURE = BENCH / "executions/gse281364-dna-lm-common-fixture-21069076/fixture"
MANIFEST_TSV = FIXTURE / "sequence_manifest.tsv"
ALLELE_FA = FIXTURE / "common_4096.alleles.fa.gz"
ROW_UNIVERSE = BENCH / "executions/model-check-220-21088696/contract/row_universe.tsv"

ALLELE_ORDER = ("REF", "ALT", "REF_RC", "ALT_RC")
FIXTURE_WIDTH = 4096
FIXTURE_VARIANT_INDEX = 2048  # forward; RC is 2047
POOL_START0_4096 = 1280
POOL_END0_4096 = 2816
BIN_BP = 128
EMBED_WIDTH = 3072
ALPHABET = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}", flush=True)


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


class ExtractionError(RuntimeError):
    pass


# --------------------------------------------------------------------------------------------------
# substrate
# --------------------------------------------------------------------------------------------------
def read_manifest() -> list[dict[str, str]]:
    with MANIFEST_TSV.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1033:
        raise ExtractionError(f"manifest carries {len(rows)} rows, expected 1033")
    for r in rows:
        if int(r["input_end0"]) - int(r["input_start0"]) != FIXTURE_WIDTH:
            raise ExtractionError(f"{r['element_id']}: fixture window is not {FIXTURE_WIDTH} bp")
        if int(r["forward_variant_index0"]) != FIXTURE_VARIANT_INDEX:
            raise ExtractionError(f"{r['element_id']}: forward variant index is not {FIXTURE_VARIANT_INDEX}")
        if int(r["variant_pos0"]) - int(r["input_start0"]) != FIXTURE_VARIANT_INDEX:
            raise ExtractionError(f"{r['element_id']}: variant offset differs from the fixture index")
        if len(r["ref"]) != 1 or len(r["alt"]) != 1:
            raise ExtractionError(f"{r['element_id']}: not a single-base substitution")
    return rows


def read_fixture_fasta() -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with gzip.open(ALLELE_FA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    if len(seqs) != 4132:
        raise ExtractionError(f"fixture fasta carries {len(seqs)} records, expected 4132")
    return seqs


def read_row_universe_split() -> dict[str, tuple[int, str]]:
    """Folds and 1-Mb blocks come from the row authority, NOT from the sequence fixture.

    The fixture's own outer_fold column disagrees with the published comparator's fold for 830 of the
    1,033 elements; it is the obsolete 6-kb assignment the campaign receipts disclaim.
    """
    out: dict[str, tuple[int, str]] = {}
    with ROW_UNIVERSE.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            value = (int(r["outer_fold"].removeprefix("fold-")), r["long_range_block_id"])
            if out.setdefault(r["element_id"], value) != value:
                raise ExtractionError(f"{r['element_id']}: row universe fold/block is not unique")
    if len(out) != 1033:
        raise ExtractionError(f"row universe carries {len(out)} elements, expected 1033")
    return out


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def validate_alphabet(seq: str) -> set[str]:
    return set(seq) - ALPHABET


def build_windows(row, fixture_seqs, fasta, length: int):
    """The four allele windows at `length`, plus the variant index and pool bin range for each.

    At 2,048 bp the window is the central 2,048 bp of the comparator's own 4,096-bp fixture record, so the
    input is a strict subsequence of the comparator's input. At 16,384 bp it is native hg38 re-extracted
    around the same variant with the same single-base substitution applied at the centre.
    """
    prefix = row["fasta_record_prefix"]
    if length == FIXTURE_WIDTH // 2:
        lo = (FIXTURE_WIDTH - length) // 2
        hi = lo + length
        fwd_ref = fixture_seqs[f"{prefix}|REF"][lo:hi]
        fwd_alt = fixture_seqs[f"{prefix}|ALT"][lo:hi]
        rc_ref = fixture_seqs[f"{prefix}|REF_RC"][lo:hi]
        rc_alt = fixture_seqs[f"{prefix}|ALT_RC"][lo:hi]
        # The central slice of the reverse complement must equal the reverse complement of the central
        # slice; assert it rather than assume the fixture's RC record is the one we think it is.
        if rc_ref != revcomp(fwd_ref) or rc_alt != revcomp(fwd_alt):
            raise ExtractionError(f"{row['element_id']}: sliced RC record is not the RC of the sliced forward")
        var_fwd = FIXTURE_VARIANT_INDEX - lo
    else:
        contig = row["contig"]
        var_pos0 = int(row["variant_pos0"])
        start0 = var_pos0 - length // 2
        end0 = start0 + length
        if start0 < 0:
            raise ExtractionError(f"{row['element_id']}: {length} bp window runs off the contig start")
        fwd_ref = fasta.fetch(contig, start0, end0).upper()
        if len(fwd_ref) != length:
            raise ExtractionError(f"{row['element_id']}: fetched {len(fwd_ref)} bp, expected {length}")
        var_fwd = var_pos0 - start0
        if fwd_ref[var_fwd] != row["ref"]:
            raise ExtractionError(
                f"{row['element_id']}: reference base {fwd_ref[var_fwd]} at the centre, manifest says {row['ref']}"
            )
        fwd_alt = fwd_ref[:var_fwd] + row["alt"] + fwd_ref[var_fwd + 1 :]
        rc_ref = revcomp(fwd_ref)
        rc_alt = revcomp(fwd_alt)
    var_rc = length - 1 - var_fwd
    # The pool region is the comparator's, expressed as an offset from the variant: +/- 768 bp.
    half = (POOL_END0_4096 - POOL_START0_4096) // 2
    centre = length // 2
    pool_lo, pool_hi = centre - half, centre + half
    if pool_lo % BIN_BP or pool_hi % BIN_BP:
        raise ExtractionError(f"pool region {pool_lo}:{pool_hi} does not tile {BIN_BP}-bp bins at {length}")
    bin_lo, bin_hi = pool_lo // BIN_BP, pool_hi // BIN_BP
    records = {
        "REF": (fwd_ref, var_fwd),
        "ALT": (fwd_alt, var_fwd),
        "REF_RC": (rc_ref, var_rc),
        "ALT_RC": (rc_alt, var_rc),
    }
    return records, (bin_lo, bin_hi)


# --------------------------------------------------------------------------------------------------
# model
# --------------------------------------------------------------------------------------------------
def load_trunk(checkpoint_path=CHECKPOINT):
    """Restore the local checkpoint and return (trunk_apply_fn, params, state, device, seconds).

    Reproduces dna_model.create()'s restore so Orbax strict validation against the initialised parameter
    tree stays on, then calls create_model directly because create() keeps the trunk private.
    """
    import jax
    import jax.numpy as jnp
    import orbax.checkpoint as ocp
    from alphagenome.models import dna_model as dna_model_api
    from alphagenome_research.model import dna_model
    from alphagenome_research.model.metadata import metadata as metadata_lib

    device = jax.local_devices()[0]
    if device.platform != "gpu":
        raise ExtractionError(f"expected a GPU device, jax offers {jax.local_devices()}")

    # Metadata for BOTH organisms at the packaged default, exactly as create() does when no metadata
    # override is supplied. Dropping mouse would change the tree the restore is validated against.
    metadata = {
        organism: metadata_lib.load(organism)
        for organism in (dna_model_api.Organism.HOMO_SAPIENS, dna_model_api.Organism.MUS_MUSCULUS)
    }
    settings = dna_model.ModelSettings()
    init_fn, _apply_fn, trunk_apply_fn, _heads, _junctions = dna_model.create_model(
        metadata,
        num_splice_sites=settings.num_splice_sites,
        splice_site_threshold=settings.splice_site_threshold,
    )
    target_shapes = jax.eval_shape(
        init_fn,
        jax.random.PRNGKey(0),
        jax.ShapeDtypeStruct((1, 2048, 4), dtype=jnp.float32),
        jax.ShapeDtypeStruct((1,), dtype=jnp.int32),
    )
    t0 = time.time()
    params, state = ocp.StandardCheckpointer().restore(
        str(checkpoint_path), target=target_shapes, strict=True
    )
    params = jax.device_put(params, device)
    state = jax.device_put(state, device)
    load_s = time.time() - t0
    log(f"checkpoint restored with strict validation in {load_s:.1f}s onto {device}")
    return jax.jit(trunk_apply_fn), params, state, device, load_s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--length", type=int, required=True, choices=[2048, 16384])
    ap.add_argument("--out", required=True, help="npz path for the pooled embeddings")
    ap.add_argument("--receipt", required=True, help="json path for the extraction receipt")
    ap.add_argument("--determinism-elements", type=int, default=16)
    ap.add_argument("--card", required=True, help="GPU tag from nvidia-smi; must be l40s")
    args = ap.parse_args()

    import jax
    import jax.numpy as jnp
    import pysam
    from alphagenome_research.model import one_hot_encoder

    length = args.length
    rows = read_manifest()
    fixture_seqs = read_fixture_fasta()
    split = read_row_universe_split()
    if set(split) != {r["element_id"] for r in rows}:
        raise ExtractionError("row universe and sequence fixture cover different elements")
    fasta = pysam.FastaFile(FASTA_PATH)
    encoder = one_hot_encoder.DNAOneHotEncoder()

    trunk, params, state, device, load_s = load_trunk()
    gpu_name = str(jax.local_devices()[0].device_kind)
    card = args.card
    if "L40S" not in gpu_name.upper() or card != "l40s":
        raise ExtractionError(
            f"this arm scores on L40S only; jax reports {gpu_name!r} and the job passed --card {card!r}"
        )
    log(f"card {card}, jax device_kind {gpu_name!r}")

    organism_index = jax.device_put(np.zeros((1,), dtype=np.int32), device)

    n = len(rows)
    element_ids = [r["element_id"] for r in rows]
    pooled_region = np.full((n, 4, EMBED_WIDTH), np.nan, dtype=np.float32)
    pooled_centre = np.full((n, 4, EMBED_WIDTH), np.nan, dtype=np.float32)
    rejects: list[dict] = []
    bins_used = None
    centre_bins: dict[str, list[int]] = {a: [] for a in ALLELE_ORDER}
    first_call_s = None
    warm_times: list[float] = []
    determinism: list[dict] = []
    embedding_shape = None

    for i, row in enumerate(rows):
        records, (bin_lo, bin_hi) = build_windows(row, fixture_seqs, fasta, length)
        if bins_used is None:
            bins_used = (bin_lo, bin_hi)
        elif bins_used != (bin_lo, bin_hi):
            raise ExtractionError("pool bin range is not constant across elements")
        bad_alleles = {}
        for allele, (seq, _v) in records.items():
            extra = validate_alphabet(seq)
            if extra:
                bad_alleles[allele] = sorted(extra)
        if bad_alleles:
            rejects.append(
                {
                    "element_id": row["element_id"],
                    "contig": row["contig"],
                    "variant_pos0": int(row["variant_pos0"]),
                    "input_length_bp": length,
                    "alleles_with_non_acgt": bad_alleles,
                    "action": "rejected_not_replaced",
                }
            )
            log(f"REJECT {row['element_id']}: non-ACGT bytes {bad_alleles}")
            continue
        for j, allele in enumerate(ALLELE_ORDER):
            seq, var_idx = records[allele]
            if len(seq) != length:
                raise ExtractionError(f"{row['element_id']}/{allele}: window is {len(seq)} bp")
            x = jax.device_put(np.asarray(encoder.encode(seq))[np.newaxis], device)
            t0 = time.time()
            emb = trunk(params, state, x, organism_index)
            e128 = np.asarray(
                jax.device_get(emb.get_sequence_embeddings(128).astype(jnp.float32)),
                dtype=np.float32,
            )
            dt = time.time() - t0
            if first_call_s is None:
                first_call_s = dt
            else:
                warm_times.append(dt)
            if embedding_shape is None:
                embedding_shape = list(e128.shape)
                log(f"embeddings_128bp shape at {length} bp: {embedding_shape}")
            if e128.shape != (1, length // BIN_BP, EMBED_WIDTH):
                raise ExtractionError(f"unexpected embedding shape {e128.shape} at {length} bp")
            if not np.isfinite(e128).all():
                raise ExtractionError(f"{row['element_id']}/{allele}: non-finite embedding")
            pooled_region[i, j] = e128[0, bin_lo:bin_hi, :].mean(axis=0)
            cb = var_idx // BIN_BP
            pooled_centre[i, j] = e128[0, cb, :]
            if i == 0:
                centre_bins[allele].append(cb)
            del emb, e128
        if i < args.determinism_elements:
            seq, var_idx = records["REF"]
            x = jax.device_put(np.asarray(encoder.encode(seq))[np.newaxis], device)
            emb2 = trunk(params, state, x, organism_index)
            e2 = np.asarray(
                jax.device_get(emb2.get_sequence_embeddings(128).astype(jnp.float32)),
                dtype=np.float32,
            )
            repeat = e2[0, bins_used[0] : bins_used[1], :].mean(axis=0)
            determinism.append(
                {
                    "element_id": row["element_id"],
                    "max_abs_diff_pooled": float(np.max(np.abs(repeat - pooled_region[i, 0]))),
                    "bit_identical": bool(np.array_equal(repeat, pooled_region[i, 0])),
                }
            )
            del emb2, e2
        if (i + 1) % 100 == 0:
            log(f"{i + 1}/{n} elements")

    keep = ~np.isnan(pooled_region[:, 0, 0])
    n_keep = int(keep.sum())
    log(f"{n_keep} of {n} elements extracted at {length} bp; {len(rejects)} rejected")
    if n_keep + len(rejects) != n:
        raise ExtractionError("kept plus rejected does not equal the element census")
    if not np.isfinite(pooled_region[keep]).all() or not np.isfinite(pooled_centre[keep]).all():
        raise ExtractionError("non-finite pooled embedding survived")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        element_ids=np.asarray(element_ids),
        allele_order=np.asarray(ALLELE_ORDER),
        extracted=keep,
        pooled_region=pooled_region,
        pooled_centre_bin=pooled_centre,
        outer_fold=np.asarray([split[e][0] for e in element_ids], dtype=np.int64),
        long_range_block_id=np.asarray([split[e][1] for e in element_ids]),
        input_length_bp=np.asarray(length),
        pool_bin_lo=np.asarray(bins_used[0]),
        pool_bin_hi=np.asarray(bins_used[1]),
        gpu_card=np.asarray(card),
    )
    log(f"wrote {out}")

    receipt = {
        "stage": "ad_arm2_reporter_frozen_probe_embeddings",
        "prespecification_sha256": "f44a912f73a1a30cddcfb60aacd5ddd0edb2e49838a0b8c6011c28be3ad267fc",
        "addendum_01_sha256": "260914ce33d1aed30dd68fb856df64bf051979f84d34d1e07004824a2b5284bb",
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_length_bp": length,
        "window_content": (
            "central_2048_of_comparator_4096_fixture_record"
            if length == 2048
            else "native_hg38_16384_recentred_on_variant_with_snv_applied"
        ),
        "native_flanks": True,
        "padding_used": False,
        "representation": "Embeddings.get_sequence_embeddings(128)",
        "embedding_shape": embedding_shape,
        "embedding_width": EMBED_WIDTH,
        "bfloat16_upcast_before_arithmetic": True,
        "pooling_primary": f"mean_over_128bp_bins_{bins_used[0]}_to_{bins_used[1]}_exclusive",
        "pooling_primary_bp": [bins_used[0] * BIN_BP, bins_used[1] * BIN_BP],
        "pooling_primary_matches_comparator_pool_region": True,
        "pooling_secondary": "single_128bp_bin_containing_the_variant",
        "centre_bin_first_element": {k: v for k, v in centre_bins.items()},
        "elements_total": n,
        "elements_extracted": n_keep,
        "windows_alphabet_validated": n * 4,
        "windows_rejected": len(rejects),
        "rejected": rejects,
        "gpu_card": card,
        "gpu_kind_reported": str(gpu_name),
        "checkpoint": str(CHECKPOINT),
        "checkpoint_orbax_strict_validation": True,
        "checkpoint_restore_seconds": round(load_s, 2),
        "first_trunk_call_seconds": None if first_call_s is None else round(first_call_s, 3),
        "warm_trunk_call_seconds_median": (
            None if not warm_times else float(np.median(np.asarray(warm_times)))
        ),
        "warm_trunk_call_seconds_total": float(sum(warm_times)) if warm_times else None,
        "trunk_calls": n_keep * 4 + len(determinism),
        "determinism_repeat_subset": determinism,
        "determinism_all_bit_identical": (
            None if not determinism else all(d["bit_identical"] for d in determinism)
        ),
        "split_source": str(ROW_UNIVERSE),
        "split_source_sha256": sha256_file(ROW_UNIVERSE),
        "sequence_source_fasta": FASTA_PATH,
        "sequence_source_fixture": str(ALLELE_FA),
        "sequence_source_fixture_sha256": sha256_file(ALLELE_FA),
        "manifest_sha256": sha256_file(MANIFEST_TSV),
        "outcomes_read": False,
        "terms": "noncommercial_local_alphagenome_weights_derivatives_inherit",
    }
    rp = pathlib.Path(args.receipt)
    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text(json.dumps(receipt, indent=1, sort_keys=True, default=str) + "\n")
    log(f"wrote {rp}")


if __name__ == "__main__":
    main()
