#!/usr/bin/env python3
"""E3 step 1: measured reporter allele effects, long-range blocks, and the element selection.

No model-API call happens here. The selection is written to disk with its seed before any score is
requested, which is the condition the prespecification states.

Outputs (tables/): mpra_allele_effects.tsv, long_range_blocks.tsv, selected_elements.tsv,
selection_receipt.json
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BENCH = PROJ / "Analysis/MASLD_Model_Benchmark/executions"
VALID = BENCH / "gse281364-validated-reconstruction-21066470/validated"
FIXTURE = BENCH / "gse281364-borzoi-native-fixture-21083008/fixture"

SEED = 20260914
N_RANDOM = 200
N_TOP = 100
PRIMARY_CONTEXT = "HepG2_control"
LR_GAP = 524_288  # the Borzoi native fixture's long-range linkage distance


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def long_range_blocks(elements: pd.DataFrame, gap: int = LR_GAP) -> pd.Series:
    """Single-linkage clusters of variant positions on a contig, gap threshold `gap`.

    This rule was verified to reproduce all 239 `borzoi_long_range_group_id` values of the fixture
    manifest exactly on its 1,033 elements before it was used to extend blocks to all 4,359.
    """
    out = {}
    for contig, g in elements.groupby("contig"):
        g = g.sort_values(["variant_pos0", "element_id"])
        cid, last = 0, None
        for pos, eid in zip(g["variant_pos0"].to_numpy(), g["element_id"].to_numpy()):
            if last is not None and pos - last >= gap:
                cid += 1
            out[eid] = f"{contig}_lr{cid:04d}"
            last = pos
    return pd.Series(out, name="long_range_block")


def verify_block_rule(elements: pd.DataFrame) -> dict:
    mf = pd.read_csv(FIXTURE / "manifest.tsv", sep="\t", dtype={"element_id": str})
    sub = elements[elements["element_id"].isin(set(mf["element_id"]))]
    mine = long_range_blocks(sub)
    df = pd.DataFrame({"mine": mine, "truth": mf.set_index("element_id")["borzoi_long_range_group_id"]})
    df = df.dropna()
    exact = bool((df.groupby("truth")["mine"].nunique() == 1).all()
                 and (df.groupby("mine")["truth"].nunique() == 1).all())
    return {"n_fixture_elements": int(df.shape[0]), "n_blocks_reproduced": int(df["mine"].nunique()),
            "n_blocks_in_fixture": int(df["truth"].nunique()), "partition_identical": exact}


def measured_effects() -> pd.DataFrame:
    ro = pd.read_csv(VALID / "replicate_outcomes.tsv.gz", sep="\t")
    ro = ro[ro["assay_state"] == "observed"].copy()
    ro["activity"] = np.log2((ro["RNA"] + 0.5) / (ro["DNA"] + 0.5))
    wide = ro.pivot_table(index=["element_id", "context_id", "experimental_replicate"],
                          columns="allele", values="activity", aggfunc="first")
    if not {"ref", "alt"} <= set(wide.columns):
        raise SystemExit("replicate_outcomes lacks both alleles")
    wide["d"] = wide["alt"] - wide["ref"]
    per = wide.reset_index()
    eff = per.groupby(["element_id", "context_id"])["d"].agg(["mean", "std", "count"]).reset_index()
    eff = eff.rename(columns={"mean": "allele_effect", "std": "allele_effect_sd_over_reps",
                              "count": "n_replicates"})
    # also carry per-replicate values so a later reliability check needs no re-derivation
    reps = per.pivot_table(index=["element_id", "context_id"], columns="experimental_replicate",
                           values="d")
    reps.columns = [f"d_rep{int(c)}" for c in reps.columns]
    return eff.merge(reps.reset_index(), on=["element_id", "context_id"], how="left")


def main() -> None:
    out_root = pathlib.Path(sys.argv[1])
    tables = out_root / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    el = pd.read_csv(VALID / "elements.tsv", sep="\t", dtype=str)
    snv = el[el["pair_state"] == "paired_snv"].copy()
    for c in ("interval_start_as_published", "interval_end_as_published", "variant_pos0",
              "variant_pos1", "reference_match_start0", "reference_match_end0",
              "oligo_difference_index0"):
        snv[c] = snv[c].astype(int)
    if snv.shape[0] != 4359:
        raise SystemExit(f"expected 4359 paired_snv elements, got {snv.shape[0]}")
    if set(snv["reference_match_orientation"]) != {"forward"}:
        raise SystemExit("not every paired_snv element is forward-oriented")

    block_check = verify_block_rule(snv)
    if not block_check["partition_identical"]:
        raise SystemExit(f"long-range block rule does not reproduce the fixture: {block_check}")
    blocks = long_range_blocks(snv)
    snv["long_range_block"] = snv["element_id"].map(blocks)
    snv["mb_block"] = snv["contig"] + ":" + (snv["variant_pos0"] // 1_000_000).astype(str)

    eff = measured_effects()
    eff = eff[eff["element_id"].isin(set(snv["element_id"]))]
    eff.to_csv(tables / "mpra_allele_effects.tsv", sep="\t", index=False)

    snv[["element_id", "contig", "variant_pos0", "variant_pos1", "genomic_ref", "genomic_alt",
         "interval_start_as_published", "interval_end_as_published", "long_range_block",
         "mb_block"]].to_csv(tables / "long_range_blocks.tsv", sep="\t", index=False)

    ids = np.array(sorted(snv["element_id"]))
    rng = np.random.default_rng(SEED)
    random_ids = set(rng.choice(ids, size=N_RANDOM, replace=False).tolist())

    prim = eff[eff["context_id"] == PRIMARY_CONTEXT].copy()
    prim["abs_effect"] = prim["allele_effect"].abs()
    prim = prim.sort_values(["abs_effect", "element_id"], ascending=[False, True])
    top_ids = set(prim.head(N_TOP)["element_id"].tolist())

    rows = []
    for eid in sorted(set(random_ids) | set(top_ids)):
        r = snv[snv["element_id"] == eid].iloc[0]
        rows.append({
            "element_id": eid,
            "in_random_stratum": eid in random_ids,
            "in_topeffect_stratum": eid in top_ids,
            "contig": r["contig"],
            "oligo_start0": int(r["interval_start_as_published"]),
            "oligo_end0": int(r["interval_end_as_published"]),
            "variant_pos0": int(r["variant_pos0"]),
            "genomic_ref": r["genomic_ref"],
            "genomic_alt": r["genomic_alt"],
            "long_range_block": r["long_range_block"],
            "mb_block": r["mb_block"],
        })
    sel = pd.DataFrame(rows)
    sel.to_csv(tables / "selected_elements.tsv", sep="\t", index=False)

    receipt = {
        "written_before_any_api_score": True,
        "seed": SEED,
        "n_universe_paired_snv": int(snv.shape[0]),
        "n_random_stratum": int(sel["in_random_stratum"].sum()),
        "n_topeffect_stratum": int(sel["in_topeffect_stratum"].sum()),
        "n_overlap": int((sel["in_random_stratum"] & sel["in_topeffect_stratum"]).sum()),
        "n_unique_elements_to_score": int(sel.shape[0]),
        "api_calls_planned": int(sel.shape[0]) * 4,
        "primary_context_for_topeffect": PRIMARY_CONTEXT,
        "long_range_block_rule": {"gap_bp": LR_GAP, **block_check,
                                  "n_blocks_all_4359": int(snv["long_range_block"].nunique()),
                                  "n_blocks_selected": int(sel["long_range_block"].nunique())},
        "inputs": {str(p.relative_to(PROJ)): sha256_file(p) for p in
                   [VALID / "elements.tsv", FIXTURE / "manifest.tsv",
                    FIXTURE / "long_range_components.tsv", FIXTURE / "source_group_map.tsv"]},
        "input_no_sha_too_large": [str((VALID / "replicate_outcomes.tsv.gz").relative_to(PROJ))],
    }
    receipt["inputs"][str((VALID / "replicate_outcomes.tsv.gz").relative_to(PROJ))] = sha256_file(
        VALID / "replicate_outcomes.tsv.gz")
    receipt.pop("input_no_sha_too_large")
    json.dump(receipt, (tables / "selection_receipt.json").open("w"), indent=1)
    print(json.dumps(receipt, indent=1))


if __name__ == "__main__":
    main()
