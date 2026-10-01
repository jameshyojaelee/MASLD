#!/usr/bin/env python3
"""Build the AD.1 gate variant set: the 3,845 matched Currin leads with an archived Atlas score.

Inputs, all read-only:
  - `c1-endpoint2-v2-.../tables/endpoint1_matched_tierA4.tsv.gz` — the 3,845-lead tier A4 matched
    set and its archived hosted-Atlas columns `alphagenome_atac_liver` / `alphagenome_dnase_liver`.
  - `c2-endogenous-head-.../inputs/currin_lead_labels.tsv.gz` — `beta_alt`, `heldout_fold`,
    `block_1mb` for all 32,322 leads, so the gate set carries the same fold and block assignment
    the zero-shot endpoint will be scored on.
  - the archived Atlas chunks, re-read here rather than trusted through tier A4, because a relayed
    number is a claim until it is found in its producing file.

Outputs three things into the run's `inputs/`:
  - `gate_variants.tsv`   one row per lead: key, chrom, pos, ref, alt, block, fold, beta_alt,
                          archived ATAC and DNase, per-track archived values.
  - `pilot_variants.tsv`  a block-disjoint pilot subset for the length sweep and the cost pilot.
  - `gate_set_census.json` counts and the archive reproduction check.

Nothing here touches the GPU.
"""

from __future__ import annotations

import argparse
import glob
import os
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import i1_common as C  # noqa: E402

P3A_CHUNKS = (
    C.PROJ
    / "GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/raw/atlas_caqtl"
)


def load_atlas_chunks() -> pd.DataFrame:
    """Re-read the archived Atlas chunks with the per-track values kept.

    Same selection as `c1_matched_endogenous.load_atlas_chunks`: the raw `.X` column of each named
    liver track, averaged. The per-track columns are kept here so the gate can report which track
    drives any disagreement rather than only the 3-track mean.
    """
    import anndata as ad

    frames = []
    chunks = sorted(glob.glob(os.path.join(P3A_CHUNKS, "chunk_*")))
    if not chunks:
        raise C.ContractError(f"no archived Atlas chunks under {P3A_CHUNKS}")
    C.log(f"reading {len(chunks)} archived Atlas caQTL chunks")
    for d in chunks:
        per = {}
        for stem, tracks in (("ATAC", C.AG_ATAC_TRACKS), ("DNASE", C.AG_DNASE_TRACKS)):
            a = ad.read_h5ad(os.path.join(d, f"{stem}.h5ad"))
            names = a.var["name"].astype(str).values
            idx = [int(np.where(names == t)[0][0]) for t in tracks]
            X = np.asarray(a.X, dtype=float)[:, idx]
            v = pd.Series(a.obs["variant"].astype(str).values)
            parts = v.str.replace(">", ":", regex=False).str.split(":", expand=True)
            frame = pd.DataFrame(
                {
                    "key": [
                        C.key_of(c, p, r, a2)
                        for c, p, r, a2 in zip(parts[0], parts[1], parts[2], parts[3])
                    ],
                    f"{stem}_raw": np.nanmean(X, axis=1),
                }
            )
            for j, t in enumerate(tracks):
                frame[f"archived_{stem.lower()}__{t.split()[0]}"] = X[:, j]
            per[stem] = frame
        frames.append(per["ATAC"].merge(per["DNASE"], on="key", how="outer", validate="one_to_one"))
    out = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["key"])
    out = out.rename(
        columns={
            "ATAC_raw": "archived_atac_liver",
            "DNASE_raw": "archived_dnase_liver",
        }
    )
    C.log(f"archived Atlas caQTL scores: {len(out)} unique variants")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="run directory")
    ap.add_argument("--pilot-blocks", type=int, default=60)
    ap.add_argument("--sweep-n", type=int, default=64)
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    (out / "inputs").mkdir(parents=True, exist_ok=True)

    tier = pd.read_csv(C.TIER_A4, sep="\t")
    C.log(f"tier A4: {len(tier)} rows, {tier.block_1mb.nunique()} blocks")

    labels = pd.read_csv(C.C2_LABELS, sep="\t")
    labels["key"] = [
        C.key_of(c, p, r, a)
        for c, p, r, a in zip(labels.chr, labels.pos_hg38, labels.ref, labels.alt)
    ]
    if labels.key.duplicated().any():
        raise C.ContractError("duplicate keys in the Currin lead label table")

    atlas = load_atlas_chunks()

    census: dict = {
        "tier_a4_rows": int(len(tier)),
        "tier_a4_blocks": int(tier.block_1mb.nunique()),
        "currin_label_rows": int(len(labels)),
        "archived_atlas_unique_variants": int(len(atlas)),
    }

    keep = [
        "key",
        "chr",
        "pos_hg38",
        "ref",
        "alt",
        "beta_alt",
        "heldout_fold",
        "block_1mb",
        "row_index",
        "strand_ambiguous",
        "reference_match",
        "peak_id",
    ]
    gate = labels[keep].merge(
        tier[["key", "alphagenome_atac_liver", "alphagenome_dnase_liver"]],
        on="key",
        how="inner",
        validate="one_to_one",
    )
    census["gate_rows_after_tierA4_join"] = int(len(gate))
    gate = gate.merge(atlas, on="key", how="left", validate="one_to_one")

    # The archive, re-read here, must reproduce the tier A4 column it was the source of.
    for chan in ("atac", "dnase"):
        a = gate[f"alphagenome_{chan}_liver"].to_numpy(float)
        b = gate[f"archived_{chan}_liver"].to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        census[f"archive_reread_vs_tierA4_{chan}_n"] = int(ok.sum())
        census[f"archive_reread_vs_tierA4_{chan}_max_abs_diff"] = (
            float(np.max(np.abs(a[ok] - b[ok]))) if ok.any() else None
        )
        census[f"archive_reread_vs_tierA4_{chan}_missing"] = int((~np.isfinite(b)).sum())

    if census["archive_reread_vs_tierA4_atac_max_abs_diff"] not in (None,) and (
        census["archive_reread_vs_tierA4_atac_max_abs_diff"] > 1e-9
    ):
        raise C.ContractError(
            "re-reading the archive does not reproduce the tier A4 ATAC column; "
            f"max abs diff {census['archive_reread_vs_tierA4_atac_max_abs_diff']}"
        )

    gate = gate.sort_values(["chr", "pos_hg38"]).reset_index(drop=True)
    census["gate_rows"] = int(len(gate))
    census["gate_blocks"] = int(gate.block_1mb.nunique())
    census["gate_per_fold"] = {
        int(f): int(n) for f, n in gate.heldout_fold.value_counts().sort_index().items()
    }
    census["gate_blocks_per_fold"] = {
        int(f): int(g.block_1mb.nunique()) for f, g in gate.groupby("heldout_fold")
    }
    census["gate_snv_only"] = bool(
        (gate.ref.str.len() == 1).all() and (gate.alt.str.len() == 1).all()
    )
    census["gate_reference_match_all"] = bool(gate.reference_match.astype(str).eq("True").all())
    census["gate_strand_ambiguous"] = int(gate.strand_ambiguous.astype(str).eq("True").sum())

    gate.to_csv(out / "inputs/gate_variants.tsv", sep="\t", index=False)
    C.log(f"gate set: {len(gate)} leads in {gate.block_1mb.nunique()} blocks")

    # Pilot subset: whole 1-Mb blocks, evenly spread across the sorted block list, so the pilot
    # correlation is itself block-resampleable and is not one chromosome's worth of leads.
    blocks = np.array(sorted(gate.block_1mb.unique()))
    take = blocks[np.linspace(0, blocks.size - 1, num=min(args.pilot_blocks, blocks.size)).astype(int)]
    pilot = gate[gate.block_1mb.isin(set(take.tolist()))].reset_index(drop=True)
    pilot.to_csv(out / "inputs/pilot_variants.tsv", sep="\t", index=False)
    census["pilot_rows"] = int(len(pilot))
    census["pilot_blocks"] = int(pilot.block_1mb.nunique())

    # Length sweep subset: the first sweep_n rows of the pilot, one variant per block where possible.
    sweep = pilot.drop_duplicates(subset=["block_1mb"]).head(args.sweep_n).reset_index(drop=True)
    sweep.to_csv(out / "inputs/sweep_variants.tsv", sep="\t", index=False)
    census["sweep_rows"] = int(len(sweep))
    census["sweep_blocks"] = int(sweep.block_1mb.nunique())

    census["inputs"] = {
        str(p): C.sha256_file(p)
        for p in (C.TIER_A4, C.C2_LABELS, C.C2_OOF, C.C2_METRICS)
    }
    census["archived_chunk_count"] = len(sorted(glob.glob(os.path.join(P3A_CHUNKS, "chunk_*"))))
    C.write_json(out / "inputs/gate_set_census.json", census)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
