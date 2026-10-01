#!/usr/bin/env python3
"""Score REF and ALT locally with the archived Atlas readout: 501-bp centre mask, DIFF_LOG2_SUM.

The archived column `alphagenome_atac_liver` is the mean over three named liver ATAC tracks of the
hosted Atlas `ATAC` scorer's raw value, and that scorer is
`CenterMaskScorer(ATAC, width=501, DIFF_LOG2_SUM)`. This script reproduces that quantity from the
local weights.

Route. `model.score_variant` cannot run in this offline build: its `_predict_variant` unconditionally
calls `jnp.asarray(splice_junction_masks.splice_sites)`, which is None whenever the model is built
without a splice-site extractor, and the packaged extractor's feather files sit behind a Google
Storage URL that compute nodes cannot reach. Measured, not assumed: the first pilot attempt failed
exactly there, at all five lengths, with `ValueError: None is not a valid value for jnp.array`.
So the two `predict_sequence` calls the runtime probe already validated are used, and the Atlas
readout is applied with the SAME packaged functions `score_variant` would have applied,
`create_center_mask` and `_apply_aggregation`. For a matched-reference substitution the two routes
see identical inputs; `i1_common.center_mask_and_aggregation` records why.

One process scores one input length, because the predict call jits on the sequence shape.

Guards, per ADAPTATION_ARM_AMENDMENT_01.md:
  - section 3: the card is queried once and written into every row; nothing merges rows from two
    cards.
  - section 4: the packaged one-hot encoder writes a zero vector for any byte that is not A/C/G/T and
    warns about nothing, so every extracted window is alphabet-validated on the interval the model
    will see; a window with a non-ACGT byte, a chromosome-edge overhang, or a reference mismatch is
    REJECTED, counted, written to the rejected table, and not silently replaced.

The score table is appended and flushed every `--flush-every` variants, so a wall-time kill leaves a
usable partial table and `--resume` skips the keys already on disk for this card.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import i1_common as C  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--flush-every", type=int, default=100)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--require-card", default="l40s")
    args = ap.parse_args()

    import jax.numpy as jnp
    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rej_path = out_path.with_name(out_path.stem + "_rejected.tsv")
    meta_path = out_path.with_name(out_path.stem + "_run.json")

    card = C.gpu_card()
    C.log(f"card {card}")
    if args.require_card and card.get("card_tag") != args.require_card:
        raise C.ContractError(
            f"scores are not reproducible across GPU architectures; this job requires "
            f"{args.require_card!r} and the node offers {card.get('card')!r}"
        )

    variants = pd.read_csv(args.variants, sep="\t")
    if args.limit:
        variants = variants.head(args.limit)

    done: set[str] = set()
    if args.resume and out_path.exists():
        prior = pd.read_csv(out_path, sep="\t")
        prior = prior[prior.card_tag.astype(str) == card["card_tag"]]
        prior = prior[prior.length_bp.astype(int) == args.length]
        done = set(prior.key.astype(str))
        C.log(f"resume: {len(done)} rows already on disk for this card and length")

    todo = variants[~variants.key.astype(str).isin(done)].reset_index(drop=True)
    C.log(f"{len(todo)} of {len(variants)} variants to score at {args.length} bp")

    # Assert the scorer settings the gate's quantity depends on before spending any GPU time.
    atac_scorer, dnase_scorer = C.atac_dnase_scorers()
    make_mask, aggregate, AggType, get_res = C.center_mask_and_aggregation()
    if get_res(atac_scorer.requested_output) != 1 or get_res(dnase_scorer.requested_output) != 1:
        raise C.ContractError("ATAC/DNASE are expected at 1-bp resolution")

    model, device, load_s = C.load_model()
    outs = [dna_client.OutputType.ATAC, dna_client.OutputType.DNASE]
    fasta = pysam.FastaFile(C.FASTA_PATH)

    rows: list[dict] = []
    rejected: list[dict] = []
    first_call_s = None
    warm: list[float] = []
    header_needed = not (out_path.exists() and done)
    rej_header_needed = not rej_path.exists()
    track_report = None
    n_rejected_total = 0

    def flush() -> None:
        nonlocal rows, rejected, header_needed, rej_header_needed
        if rows:
            pd.DataFrame(rows).to_csv(
                out_path,
                sep="\t",
                index=False,
                mode="w" if header_needed else "a",
                header=header_needed,
            )
            header_needed = False
            rows = []
        if rejected:
            pd.DataFrame(rejected).to_csv(
                rej_path,
                sep="\t",
                index=False,
                mode="w" if rej_header_needed else "a",
                header=rej_header_needed,
            )
            rej_header_needed = False
            rejected = []

    t_start = time.time()
    for n, r in enumerate(todo.itertuples(index=False), start=1):
        chrom = str(r.chr)
        pos1 = int(r.pos_hg38)
        ref, alt = str(r.ref).upper(), str(r.alt).upper()
        start0, end = C.window_bounds(pos1, args.length)
        check = C.validate_window(fasta, chrom, start0, end, pos1, ref)
        if not check["acgt_ok"] or not check["reference_match"]:
            n_rejected_total += 1
            rejected.append(
                {
                    "key": r.key,
                    "chr": chrom,
                    "pos_hg38": pos1,
                    "ref": ref,
                    "alt": alt,
                    "length_bp": args.length,
                    "reason": (
                        "out_of_bounds"
                        if check["out_of_bounds"]
                        else ("non_acgt" if check["n_non_acgt"] else "reference_mismatch")
                    ),
                    **{k: v for k, v in check.items() if k != "non_acgt_characters"},
                    "non_acgt_characters": str(check["non_acgt_characters"]),
                }
            )
            continue

        ref_seq, alt_seq = C.extract_ref_alt(fasta, chrom, start0, end, pos1, ref, alt)
        interval = genome.Interval(chromosome=chrom, start=start0, end=end)
        variant = genome.Variant(
            chromosome=chrom, position=pos1, reference_bases=ref, alternate_bases=alt
        )

        t0 = time.time()
        pred = {}
        for arm, seq in (("ref", ref_seq), ("alt", alt_seq)):
            pred[arm] = model.predict_sequence(
                sequence=seq,
                requested_outputs=outs,
                ontology_terms=C.LIVER_TERMS,
                interval=interval,
            )
        dt = time.time() - t0
        if first_call_s is None:
            first_call_s = dt
            C.log(f"first variant (includes compile) {dt:.2f}s")
        else:
            warm.append(dt)

        row = {
            "key": r.key,
            "chr": chrom,
            "pos_hg38": pos1,
            "ref": ref,
            "alt": alt,
            "block_1mb": r.block_1mb,
            "heldout_fold": r.heldout_fold,
            "beta_alt": r.beta_alt,
            "length_bp": args.length,
            "window_start0": start0,
            "window_end": end,
            "mask_width_bp": C.CENTER_MASK_WIDTH,
            "aggregation": "DIFF_LOG2_SUM",
            "card": card.get("card"),
            "card_tag": card.get("card_tag"),
            "driver_version": card.get("driver_version"),
            "compute_capability": card.get("compute_capability"),
            "seconds": round(dt, 4),
            "n_non_acgt": check["n_non_acgt"],
        }

        report = {}
        for chan, attr, scorer, want in (
            ("atac", "atac", atac_scorer, C.AG_ATAC_TRACKS),
            ("dnase", "dnase", dnase_scorer, C.AG_DNASE_TRACKS),
        ):
            td_ref = getattr(pred["ref"], attr)
            td_alt = getattr(pred["alt"], attr)
            if td_ref is None or td_ref.values is None or td_ref.values.size == 0:
                raise C.ContractError(f"{chan} returned no track values")
            if int(td_ref.resolution) != 1:
                raise C.ContractError(f"{chan} resolution is {td_ref.resolution}, expected 1")
            names = td_ref.metadata["name"].astype(str).values
            idx = []
            for t in want:
                hit = np.where(names == t)[0]
                if hit.size != 1:
                    raise C.ContractError(
                        f"track {t!r} appears {hit.size} times among the {len(names)} returned "
                        f"{chan} tracks; the pinned liver ontology must select it exactly once"
                    )
                idx.append(int(hit[0]))

            mask = make_mask(
                interval.as_unstranded(), variant, width=C.CENTER_MASK_WIDTH, resolution=1
            )
            if int(mask.sum()) != C.CENTER_MASK_WIDTH:
                raise C.ContractError(
                    f"centre mask covers {int(mask.sum())} positions, expected {C.CENTER_MASK_WIDTH}"
                )
            scores = np.asarray(
                aggregate(
                    jnp.asarray(np.asarray(td_ref.values, dtype=np.float32)),
                    jnp.asarray(np.asarray(td_alt.values, dtype=np.float32)),
                    jnp.asarray(mask),
                    aggregation_type=scorer.aggregation_type,
                ),
                dtype=np.float64,
            )
            row[f"local_{chan}_liver"] = float(np.mean(scores[idx]))
            for t, j in zip(want, idx):
                row[f"local_{chan}__{t.split()[0]}"] = float(scores[j])
            report[chan] = {
                "tracks_returned": int(len(names)),
                "track_names_returned": [str(x) for x in names],
                "liver_indices_used": idx,
                "mask_positions": int(mask.sum()),
            }
            if chan == "atac":
                row["mask_start0_offset"] = int(np.flatnonzero(mask.ravel())[0])
                row["mask_end0_offset"] = int(np.flatnonzero(mask.ravel())[-1]) + 1

        if track_report is None:
            track_report = report
            C.log(f"track selection {report}")

        rows.append(row)
        del pred

        if n % args.flush_every == 0:
            flush()
            rate = (time.time() - t_start) / n
            C.log(
                f"{n}/{len(todo)} done, {rate:.3f}s per variant, "
                f"eta {(len(todo) - n) * rate / 60:.1f} min, rejected {n_rejected_total}"
            )
    flush()

    meta = {
        "length_bp": args.length,
        "variants_requested": int(len(variants)),
        "variants_attempted": int(len(todo)),
        "variants_scored": int(len(todo) - n_rejected_total),
        "rejected_windows": int(n_rejected_total),
        "first_variant_seconds": first_call_s,
        "warm_call_seconds_median": float(np.median(warm)) if warm else None,
        "warm_call_seconds_mean": float(np.mean(warm)) if warm else None,
        "total_seconds": round(time.time() - t_start, 1),
        "checkpoint_load_seconds": round(load_s, 2),
        "checkpoint": str(C.CHECKPOINT),
        "fasta": C.FASTA_PATH,
        "ontology_terms": C.LIVER_TERMS,
        "requested_outputs": ["ATAC", "DNASE"],
        "atac_liver_track_names": C.AG_ATAC_TRACKS,
        "dnase_liver_track_names": C.AG_DNASE_TRACKS,
        "track_selection": track_report,
        "readout": (
            "packaged create_center_mask(width=501, resolution=1) with packaged _apply_aggregation "
            "DIFF_LOG2_SUM, applied to two predict_sequence calls; identical to the hosted Atlas "
            "ATAC/DNASE scorer settings taken from RECOMMENDED_VARIANT_SCORERS"
        ),
        "route_note": (
            "model.score_variant is unusable offline: _predict_variant calls "
            "jnp.asarray(splice_sites) and splice_sites is None without a splice-site extractor, "
            "whose feather files are behind a Google Storage URL"
        ),
        "memory_peak": C.memory_stats(device),
        **card,
    }
    C.write_json(meta_path, meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
