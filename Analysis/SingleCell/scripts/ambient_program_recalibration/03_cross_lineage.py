#!/usr/bin/env python3
"""Biological-donor-collapsed same-atlas localization of the two hero programs."""

from __future__ import annotations

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

from common import (
    MEMBERSHIP,
    REGISTRY,
    SEED,
    bh_adjust,
    build_sample_to_donor,
    linear_fit,
    load_donor_metadata,
    refuse_existing,
    require,
    sha256,
)


CANDIDATE = Path(os.environ["CAND_ROOT"])
RESULTS = CANDIDATE / "results"
ATLAS = (
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
    / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad"
)
HEROES = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": {
        "short_name": "ECM/IGFBP7",
        "target_lineage": "Fibroblasts",
    },
    "hotspot_hepatocytes_48f39dd4d817a10e": {
        "short_name": "ductular-injury/BICC1",
        "target_lineage": "Cholangiocytes",
    },
}
LINEAGES = (
    "Hepatocytes",
    "Fibroblasts",
    "Cholangiocytes",
    "Endothelial cells",
    "Macrophages",
    "T cells",
)
MIN_CELLS = (50, 20)
ANNOTATION_FILTERS = ("all_annotated_cells", "cell_type_conf_ge_0.90")


def load_programs() -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    registry = pd.read_csv(REGISTRY, sep="\t")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    hero_registry = registry[registry["program_uid"].isin(HEROES)].copy()
    require(len(hero_registry) == 2, "hero program registry mismatch")
    programs: dict[str, dict[str, float]] = {}
    for uid in HEROES:
        rows = membership[membership["program_uid"] == uid]
        require(len(rows) >= 8, f"too few frozen genes for {uid}")
        weights = rows.groupby("source_gene")["original_l1_weight"].sum().astype(float)
        require(np.isclose(weights.sum(), 1.0, atol=1e-10), f"weight sum drift: {uid}")
        programs[uid] = weights.to_dict()
    return hero_registry, programs


def aggregate_donor_lineage(
    programs: dict[str, dict[str, float]], chunk_size: int = 5000
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], list[str]]:
    sample_to_donor = build_sample_to_donor()
    with h5py.File(ATLAS, "r") as handle:
        obs = read_elem(handle["obs"]).reset_index(drop=True)
        var = read_elem(handle["var"])
        require(
            {"sample", "dataset", "cell_type", "cell_type_conf"} <= set(obs),
            "global atlas metadata lacks cross-lineage fields",
        )
        obs["sample"] = obs["sample"].astype(str)
        obs["donor"] = obs["sample"].map(sample_to_donor).fillna(obs["sample"])
        obs["dataset"] = obs["dataset"].astype(str)
        obs["cell_type"] = obs["cell_type"].astype(str)
        relevant = obs["cell_type"].isin(LINEAGES).to_numpy()
        genes = sorted(set().union(*[set(weights) for weights in programs.values()]) & set(var.index))
        positions = var.index.get_indexer(genes)
        require(len(genes) >= 8 and np.all(positions >= 0), "hero gene axis is insufficient")

        group_key = obs["donor"] + "||" + obs["cell_type"]
        relevant_keys = group_key[relevant]
        levels = pd.Index(sorted(relevant_keys.unique()))
        level_map = pd.Series(np.arange(len(levels)), index=levels)
        codes = group_key.map(level_map).fillna(-1).astype(int).to_numpy()
        confidence = pd.to_numeric(obs["cell_type_conf"], errors="coerce").fillna(0).to_numpy()
        accumulators = {
            "all_annotated_cells": (
                np.zeros((len(levels), len(genes)), dtype=np.float64),
                np.zeros(len(levels), dtype=np.int64),
            ),
            "cell_type_conf_ge_0.90": (
                np.zeros((len(levels), len(genes)), dtype=np.float64),
                np.zeros(len(levels), dtype=np.int64),
            ),
        }
        matrix = handle["layers/counts"]
        matrix_shape = tuple(int(value) for value in matrix.attrs["shape"])
        require(matrix_shape[0] == len(obs), "global atlas cell axis drift")
        for start in range(0, matrix_shape[0], chunk_size):
            end = min(start + chunk_size, matrix_shape[0])
            indptr = np.asarray(matrix["indptr"][start : end + 1], dtype=np.int64)
            data_start, data_end = int(indptr[0]), int(indptr[-1])
            block = sp.csr_matrix(
                (
                    np.asarray(matrix["data"][data_start:data_end]),
                    np.asarray(matrix["indices"][data_start:data_end], dtype=np.int32),
                    indptr - data_start,
                ),
                shape=(end - start, matrix_shape[1]),
            )
            totals = np.asarray(block.sum(axis=1)).ravel()
            selected = block[:, positions].astype(np.float64).tocsr()
            scale = np.divide(
                1e4,
                totals,
                out=np.zeros_like(totals, dtype=np.float64),
                where=totals > 0,
            )
            selected = sp.diags(scale) @ selected
            selected.data = np.log1p(selected.data)
            block_codes = codes[start:end]
            for annotation_filter, keep_confidence in (
                ("all_annotated_cells", np.ones(end - start, dtype=bool)),
                ("cell_type_conf_ge_0.90", confidence[start:end] >= 0.90),
            ):
                keep = (block_codes >= 0) & keep_confidence
                if not keep.any():
                    continue
                kept_codes = block_codes[keep]
                indicator = sp.csr_matrix(
                    (
                        np.ones(len(kept_codes)),
                        (kept_codes, np.arange(len(kept_codes))),
                    ),
                    shape=(len(levels), len(kept_codes)),
                )
                sums, counts = accumulators[annotation_filter]
                sums += (indicator @ selected[keep, :]).toarray()
                counts += np.bincount(kept_codes, minlength=len(levels))

    donor_dataset = obs.drop_duplicates("donor").set_index("donor")["dataset"]
    require(
        obs.groupby("donor")["dataset"].nunique().max() == 1,
        "biological donor spans multiple datasets",
    )
    outputs: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for annotation_filter, (sums, counts) in accumulators.items():
        eligible = counts > 0
        meta = pd.DataFrame({"group_id": levels[eligible], "n_cells": counts[eligible]})
        meta[["donor", "lineage"]] = meta["group_id"].str.split(
            "||", n=1, expand=True, regex=False
        )
        meta["dataset"] = meta["donor"].map(donor_dataset)
        require(not meta["dataset"].isna().any(), "donor dataset join failed")
        means = sums[eligible, :] / counts[eligible, None]
        outputs[annotation_filter] = (pd.DataFrame(means, columns=genes), meta)
    return outputs, genes


def score_programs(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    programs: dict[str, dict[str, float]],
    minimum_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    keep = metadata["n_cells"] >= minimum_cells
    expression = expression.loc[keep].reset_index(drop=True)
    metadata = metadata.loc[keep].reset_index(drop=True)
    score = pd.DataFrame(index=expression.index)
    coverage_rows = []
    z = expression.copy()
    for dataset, index in metadata.groupby("dataset").groups.items():
        index = list(index)
        means = expression.loc[index].mean(axis=0)
        standard_deviations = expression.loc[index].std(axis=0, ddof=1).replace(0, np.nan)
        z.loc[index] = (expression.loc[index] - means) / standard_deviations
    for uid, weights in programs.items():
        genes = [gene for gene in weights if gene in z and np.isfinite(z[gene]).any()]
        retained_weight = float(sum(weights[gene] for gene in genes))
        testable = len(genes) >= 8 and retained_weight >= 0.80
        coverage_rows.append(
            {
                "program_uid": uid,
                "minimum_cells": minimum_cells,
                "n_frozen_genes": len(weights),
                "n_observed_genes": len(genes),
                "retained_l1_weight": retained_weight,
                "testable": testable,
            }
        )
        if not testable:
            score[uid] = np.nan
            continue
        normalized_weights = np.array([weights[gene] for gene in genes], dtype=float)
        normalized_weights /= normalized_weights.sum()
        score[uid] = np.asarray(z[genes].fillna(0.0)) @ normalized_weights
    return pd.concat([metadata, score], axis=1), pd.DataFrame(coverage_rows)


def paired_fit(
    scores: pd.DataFrame, uid: str, comparison_lineage: str
) -> tuple[dict[str, object], pd.DataFrame]:
    subset = scores[scores["lineage"].isin(["Hepatocytes", comparison_lineage])]
    wide = subset.pivot_table(
        index=["dataset", "donor"], columns="lineage", values=uid, aggfunc="first"
    ).reset_index()
    if "Hepatocytes" not in wide or comparison_lineage not in wide:
        raise RuntimeError("paired lineage columns unavailable")
    wide = wide.dropna(subset=["Hepatocytes", comparison_lineage]).copy()
    wide["difference"] = wide[comparison_lineage] - wide["Hepatocytes"]
    counts = wide["dataset"].value_counts()
    eligible_datasets = sorted(counts[counts >= 3].index)
    wide = wide[wide["dataset"].isin(eligible_datasets)].copy()
    require(len(wide) >= 6, "fewer than six paired donors")
    require(len(eligible_datasets) >= 2, "fewer than two paired datasets")
    design = pd.get_dummies(wide["dataset"], drop_first=False, dtype=float)
    contrast = np.repeat(1.0 / design.shape[1], design.shape[1])
    fit = linear_fit(wide["difference"].to_numpy(float), design.to_numpy(float), contrast)
    fit["n_paired_donors"] = int(wide["donor"].nunique())
    fit["n_datasets"] = int(len(eligible_datasets))
    fit["eligible_datasets"] = ";".join(eligible_datasets)
    return fit, wide


def main() -> None:
    outputs = {
        "scores": RESULTS / "hero_lineage_scores.tsv.gz",
        "contrasts": RESULTS / "hero_lineage_contrasts.tsv",
        "lodo": RESULTS / "hero_lineage_lodo.tsv",
        "promotion": RESULTS / "hero_lineage_promotion.tsv",
        "census": RESULTS / "hero_lineage_census.tsv",
        "coverage": RESULTS / "hero_lineage_weight_coverage.tsv",
        "provenance": RESULTS / "hero_lineage_provenance.json",
    }
    for path in outputs.values():
        refuse_existing(path)
    hero_registry, programs = load_programs()
    aggregates, genes = aggregate_donor_lineage(programs)
    donor_metadata = load_donor_metadata()[
        ["donor", "dataset", "disease_stage_coarse", "exclude"]
    ]

    all_score_tables = []
    coverage_tables = []
    contrast_rows = []
    lodo_rows = []
    census_rows = []
    for annotation_filter in ANNOTATION_FILTERS:
        expression, metadata = aggregates[annotation_filter]
        for minimum_cells in MIN_CELLS:
            scores, coverage = score_programs(
                expression, metadata, programs, minimum_cells
            )
            scores = scores.merge(
                donor_metadata,
                on=["donor", "dataset"],
                how="left",
                validate="many_to_one",
            )
            require(not scores["disease_stage_coarse"].isna().all(), "diagnosis join failed")
            scores.insert(0, "minimum_cells", minimum_cells)
            scores.insert(0, "annotation_filter", annotation_filter)
            all_score_tables.append(scores)
            coverage.insert(0, "annotation_filter", annotation_filter)
            coverage_tables.append(coverage)
            for lineage, group in scores.groupby("lineage"):
                census_rows.append(
                    {
                        "annotation_filter": annotation_filter,
                        "minimum_cells": minimum_cells,
                        "lineage": lineage,
                        "n_donor_lineage_profiles": int(len(group)),
                        "n_biological_donors": int(group["donor"].nunique()),
                        "n_datasets": int(group["dataset"].nunique()),
                    }
                )
            for uid in HEROES:
                for comparison_lineage in LINEAGES[1:]:
                    row = {
                        "program_uid": uid,
                        "program_name": HEROES[uid]["short_name"],
                        "comparison": f"{comparison_lineage}-minus-Hepatocytes",
                        "comparison_lineage": comparison_lineage,
                        "annotation_filter": annotation_filter,
                        "minimum_cells": minimum_cells,
                    }
                    try:
                        fit, paired = paired_fit(scores, uid, comparison_lineage)
                        row.update(fit)
                        row["failure_reason"] = ""
                        for held_out in sorted(paired["dataset"].unique()):
                            lodo = row.copy()
                            lodo["held_out_dataset"] = held_out
                            try:
                                lodo_fit, _ = paired_fit(
                                    scores[scores["dataset"] != held_out], uid, comparison_lineage
                                )
                                lodo.update(lodo_fit)
                                lodo["failure_reason"] = ""
                            except RuntimeError as error:
                                lodo["beta"] = np.nan
                                lodo["failure_reason"] = str(error)
                            lodo_rows.append(lodo)
                    except RuntimeError as error:
                        for key in (
                            "beta", "se", "pvalue", "hc3_se", "hc3_pvalue", "ci_low",
                            "ci_high", "hc3_ci_low", "hc3_ci_high", "n_paired_donors",
                            "n_datasets", "eligible_datasets",
                        ):
                            row[key] = np.nan
                        row["failure_reason"] = str(error)
                    contrast_rows.append(row)

    contrasts = pd.DataFrame(contrast_rows)
    for _, index in contrasts.groupby(["annotation_filter", "minimum_cells"]).groups.items():
        index = list(index)
        require(len(index) == 10, "cross-lineage multiplicity family is not ten")
        contrasts.loc[index, "qvalue"] = bh_adjust(contrasts.loc[index, "pvalue"])
        contrasts.loc[index, "hc3_qvalue"] = bh_adjust(
            contrasts.loc[index, "hc3_pvalue"]
        )

    all_scores = pd.concat(all_score_tables, ignore_index=True)
    lodo = pd.DataFrame(lodo_rows)
    promotion_rows = []
    for uid, definition in HEROES.items():
        target = definition["target_lineage"]
        primary = contrasts[
            (contrasts["program_uid"] == uid)
            & (contrasts["comparison_lineage"] == target)
            & (contrasts["annotation_filter"] == "all_annotated_cells")
            & (contrasts["minimum_cells"] == 50)
        ]
        require(len(primary) == 1, "primary hero contrast cardinality drift")
        primary = primary.iloc[0]
        top_lineage_by_threshold = {}
        for minimum_cells in MIN_CELLS:
            subset = all_scores[
                (all_scores["annotation_filter"] == "all_annotated_cells")
                & (all_scores["minimum_cells"] == minimum_cells)
            ]
            medians = subset.groupby("lineage")[uid].median()
            top_lineage_by_threshold[minimum_cells] = (
                str(medians.idxmax()) if len(medians.dropna()) else ""
            )
        high_confidence = contrasts[
            (contrasts["program_uid"] == uid)
            & (contrasts["comparison_lineage"] == target)
            & (contrasts["annotation_filter"] == "cell_type_conf_ge_0.90")
            & (contrasts["minimum_cells"] == 50)
        ]
        require(len(high_confidence) == 1, "high-confidence hero contrast drift")
        hero_lodo = lodo[
            (lodo["program_uid"] == uid)
            & (lodo["comparison_lineage"] == target)
            & (lodo["annotation_filter"] == "all_annotated_cells")
            & (lodo["minimum_cells"] == 50)
            & lodo["beta"].notna()
        ]
        primary_beta = pd.to_numeric(pd.Series([primary["beta"]]), errors="coerce").iloc[0]
        primary_qvalue = pd.to_numeric(
            pd.Series([primary["qvalue"]]), errors="coerce"
        ).iloc[0]
        primary_n_donors = pd.to_numeric(
            pd.Series([primary["n_paired_donors"]]), errors="coerce"
        ).iloc[0]
        primary_n_datasets = pd.to_numeric(
            pd.Series([primary["n_datasets"]]), errors="coerce"
        ).iloc[0]
        high_confidence_beta = pd.to_numeric(
            pd.Series([high_confidence.iloc[0]["beta"]]), errors="coerce"
        ).iloc[0]
        criteria = {
            "positive_primary_estimate": bool(np.isfinite(primary_beta) and primary_beta > 0),
            "primary_bh_q_lt_0.05": bool(
                np.isfinite(primary_qvalue) and primary_qvalue < 0.05
            ),
            "at_least_20_paired_donors": bool(
                np.isfinite(primary_n_donors) and primary_n_donors >= 20
            ),
            "at_least_3_datasets": bool(
                np.isfinite(primary_n_datasets) and primary_n_datasets >= 3
            ),
            "at_least_3_estimable_lodo": bool(len(hero_lodo) >= 3),
            "all_estimable_lodo_positive": bool(
                len(hero_lodo) >= 3 and (hero_lodo["beta"].astype(float) > 0).all()
            ),
            "target_top_lineage_min50": top_lineage_by_threshold[50] == target,
            "target_top_lineage_min20": top_lineage_by_threshold[20] == target,
            "high_confidence_direction_positive": bool(
                np.isfinite(high_confidence_beta) and high_confidence_beta > 0
            ),
        }
        promotion_rows.append(
            {
                "program_uid": uid,
                "program_name": definition["short_name"],
                "target_lineage": target,
                "top_lineage_min50": top_lineage_by_threshold[50],
                "top_lineage_min20": top_lineage_by_threshold[20],
                "primary_beta": primary["beta"],
                "primary_qvalue": primary["qvalue"],
                "primary_hc3_qvalue": primary["hc3_qvalue"],
                "n_paired_donors": primary["n_paired_donors"],
                "n_datasets": primary["n_datasets"],
                "n_estimable_lodo": int(len(hero_lodo)),
                **criteria,
                "program_pass": bool(all(criteria.values())),
            }
        )
    promotion = pd.DataFrame(promotion_rows)
    both_pass = bool(promotion["program_pass"].all())
    promotion["main_figure_family_eligible"] = both_pass
    promotion["destination"] = "Figure_4E" if both_pass else "Figure_S4"

    all_scores.to_csv(outputs["scores"], sep="\t", index=False, compression="gzip")
    contrasts.to_csv(outputs["contrasts"], sep="\t", index=False)
    lodo.to_csv(outputs["lodo"], sep="\t", index=False)
    promotion.to_csv(outputs["promotion"], sep="\t", index=False)
    pd.DataFrame(census_rows).to_csv(outputs["census"], sep="\t", index=False)
    pd.concat(coverage_tables, ignore_index=True).to_csv(
        outputs["coverage"], sep="\t", index=False
    )
    provenance = {
        "seed": SEED,
        "atlas": str(ATLAS.resolve()),
        "atlas_sha256": sha256(ATLAS),
        "registry_sha256": sha256(REGISTRY),
        "membership_sha256": sha256(MEMBERSHIP),
        "programs": list(HEROES),
        "lineages": list(LINEAGES),
        "minimum_cells": list(MIN_CELLS),
        "annotation_filters": list(ANNOTATION_FILTERS),
        "score_definition": "mean per-cell log1p(CP10K), gene-z within dataset retaining lineage contrasts, frozen L1-weighted mean",
        "inferential_unit": "biological donor",
        "multiplicity_family_size": 10,
        "main_figure_family_eligible": both_pass,
        "claim_boundary": "same-atlas descriptive localization; not lineage origin or cell autonomy",
    }
    outputs["provenance"].write_text(json.dumps(provenance, indent=2) + "\n")
    print(promotion.to_string(index=False))


if __name__ == "__main__":
    np.random.seed(SEED)
    main()
