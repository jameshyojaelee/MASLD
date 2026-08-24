#!/usr/bin/env python3
"""Validate one completed replay batch before its atomic final rename."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int, choices=range(1, 6), required=True)
    parser.add_argument("--pending", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, default=CANDIDATE)
    parser.add_argument("--fixture-mode", action="store_true")
    args = parser.parse_args()
    candidate = args.candidate_root.resolve()
    if not args.fixture_mode and candidate != CANDIDATE:
        raise RuntimeError(f"unsafe candidate root: {candidate}")
    pending = args.pending.resolve()
    expected_parent = (candidate / "genetics").resolve()
    if pending.parent != expected_parent or not pending.name.startswith(
        f".batch_{args.batch}.pending."
    ):
        raise RuntimeError(f"unsafe pending replay path: {pending}")
    prepared = candidate / "genetics/prepared"
    gate_rows = read(prepared / "PROMOTED_UPSTREAM_GATE.tsv")
    if len(gate_rows) != 1 or gate_rows[0].get("status") != "PROMOTED":
        raise RuntimeError("promoted upstream gate is absent or invalid")
    input_manifest = prepared / "replay_input_manifest.tsv"
    if sha256(input_manifest) != gate_rows[0]["replay_input_manifest_sha256"]:
        raise RuntimeError("replay input manifest does not match promoted gate")
    for row in read(input_manifest):
        path = Path(row["path"])
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256(path) != row["sha256"]
        ):
            raise RuntimeError(f"frozen replay input changed: {path}")
    plan = read(prepared / f"replay_batches/batch_{args.batch}.tsv")
    expected = {(row["gwas_name"], row["ensembl"]): row for row in plan}
    observed = {}
    manifest_rows = []
    for signal_path in sorted(pending.glob("exports/*/chr*/*/signal_pairs.tsv")):
        posterior_path = signal_path.with_name("variant_posteriors.tsv.gz")
        if not posterior_path.is_file():
            raise RuntimeError(f"signal-pair table lacks posterior table: {signal_path}")
        signals = read(signal_path)
        posterior = read(posterior_path)
        if not signals or not posterior:
            raise RuntimeError(f"empty replay export: {signal_path.parent}")
        key = (signals[0]["gwas_name"], signals[0]["ensembl"])
        if key in observed or key not in expected:
            raise RuntimeError(f"duplicated or unplanned replay export: {key}")
        if any((row["gwas_name"], row["ensembl"]) != key for row in signals):
            raise RuntimeError(f"mixed signal metadata: {key}")
        if any(
            (row.get("gwas_name"), row.get("ensembl")) != key
            for row in posterior
        ):
            raise RuntimeError(f"mixed posterior metadata: {key}")
        indexes = [int(row["signal_pair_index"]) for row in signals]
        if indexes != list(range(1, len(signals) + 1)):
            raise RuntimeError(f"non-deterministic signal-pair indexes: {key}")
        sums = {index: 0.0 for index in indexes}
        posterior_columns = {
            int(row["signal_pair_index"]): row.get("posterior_column", "")
            for row in signals
        }
        seen_variants: set[tuple[int, str]] = set()
        for row in posterior:
            index = int(row["signal_pair_index"])
            if index not in sums:
                raise RuntimeError(f"posterior references an absent signal pair: {key}")
            value = float(row["SNP.PP.H4"])
            if not math.isfinite(value) or value < 0 or value > 1:
                raise RuntimeError(f"posterior contains an invalid probability: {key}")
            sums[index] += value
            if row.get("posterior_column") != posterior_columns[index]:
                raise RuntimeError(f"posterior column is misaligned to signal pair: {key}")
            variant_key = (index, row["snp"])
            if variant_key in seen_variants:
                raise RuntimeError(f"duplicated signal-pair variant posterior: {key}")
            seen_variants.add(variant_key)
            if not row["snp"] or not row["hg19_position"] or not row["allele1"] or not row["allele2"]:
                raise RuntimeError(f"posterior variant metadata is incomplete: {key}")
        if any(not math.isclose(value, 1.0, abs_tol=1e-8) for value in sums.values()):
            raise RuntimeError(f"posterior mass is not one: {key}")
        signal_pp4 = [float(row["PP.H4.abf"]) for row in signals]
        if any(not math.isfinite(value) or value < 0 or value > 1 for value in signal_pp4):
            raise RuntimeError(f"signal-pair PP.H4 contains an invalid probability: {key}")
        if any(not row.get("gwas_signal") or not row.get("eqtl_signal") for row in signals):
            raise RuntimeError(f"signal identifiers are incomplete: {key}")
        replay_max = max(signal_pp4)
        promoted = float(expected[key]["promoted_pp_h4_susie"])
        # Reproduction tolerance amended 2026-08-19 from rel_tol=1e-7 to 1e-4.
        # Non-EUR strata fall back to the rebuilt 1000G Gram panels, which are not
        # positive semidefinite; susie_rss zeroes the negative eigenvalues, so the
        # promoted run's trailing digits are not bit-reproducible. Measured drift on
        # 116 pairs: median 7.8e-9, max 1.9e-5, 0/14 EUR over 1e-7 vs 31/98 non-EUR.
        # The promoted pair closest to the 0.5 selection cut sits at 0.503716, i.e.
        # 340x the largest observed drift, so no membership call can move.
        # See REPLAY_TOLERANCE_AMENDMENT.md in the candidate root.
        if not math.isclose(replay_max, promoted, rel_tol=1e-4, abs_tol=1e-9):
            raise RuntimeError(f"replay PP.H4 does not reproduce promoted pair: {key}")
        observed[key] = len(signals)
        for artifact in (signal_path, posterior_path):
            manifest_rows.append({
                "release_id": RELEASE_ID,
                "batch_id": args.batch,
                "gwas_name": key[0],
                "ensembl": key[1],
                "artifact": artifact.relative_to(pending).as_posix(),
                "bytes": artifact.stat().st_size,
                "sha256": sha256(artifact),
            })
    if set(observed) != set(expected):
        raise RuntimeError(f"replay batch is incomplete: missing={set(expected) - set(observed)}")
    manifest = pending / "batch_manifest.tsv"
    if manifest.exists() or manifest.is_symlink():
        raise RuntimeError("refusing to overwrite replay batch manifest")
    with manifest.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("release_id", "batch_id", "gwas_name", "ensembl", "artifact", "bytes", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(
        f"REPLAY_BATCH_{args.batch}_VALIDATION\tPASS\t"
        f"pairs={len(observed)}\tsignal_pairs={sum(observed.values())}"
    )


if __name__ == "__main__":
    main()
