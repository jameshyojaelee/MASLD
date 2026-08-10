#!/usr/bin/env python3
"""Donor-collapse the integrated atlas for the frozen M8/M20 identity audit."""

from __future__ import annotations

import csv
import math
from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy.stats import ttest_rel

from bridge_common import CANDIDATE_ROOT, PROJECT_ROOT, require_validated_seal, write_tsv


ATLAS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad"
M8 = "hotspot_hepatocytes_f05c535ae5bbc0b9"
M20 = "hotspot_hepatocytes_48f39dd4d817a10e"
AMBIENT = {"ALB", "APOA1", "APOA2", "APOC1", "FGA", "FGB", "FGG", "TTR", "HP", "ORM1", "ORM2"}


def bh_adjust(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment over finite values, preserving row order."""
    output = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.notna() & np.isfinite(values)
    if not finite.any():
        return output
    observed = values.loc[finite].astype(float)
    order = np.argsort(observed.to_numpy())
    ranked = observed.to_numpy()[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    original = np.empty_like(adjusted)
    original[order] = adjusted
    output.loc[observed.index] = original
    return output


def load_programs() -> tuple[dict[str, dict[str, float]], dict[str, dict[str, str]]]:
    membership_path = CANDIDATE_ROOT / "frozen_inputs/program_membership__program_membership_v2.tsv"
    registry_path = CANDIDATE_ROOT / "frozen_inputs/program_registry__program_registry_v2.tsv"
    programs: dict[str, dict[str, float]] = {M8: {}, M20: {}}
    with membership_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["program_uid"] in programs and row["mapped_symbol"]:
                programs[row["program_uid"]][row["mapped_symbol"]] = programs[row["program_uid"]].get(row["mapped_symbol"], 0.0) + float(row["original_l1_weight"])
    registry = {}
    with registry_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["program_uid"] in programs:
                registry[row["program_uid"]] = row
    return programs, registry


def aggregate_selected_genes(
    programs: dict[str, dict[str, float]], chunk_size: int = 5_000
) -> tuple[dict[str, tuple[pd.DataFrame, pd.DataFrame]], list[str]]:
    """Stream the 1.14-billion-entry CSR layer without loading the full atlas.

    The source H5AD stores both X and counts as large CSR matrices.  Loading the
    complete AnnData would needlessly materialize both.  This routine reads the
    count layer in row chunks, extracts only M8/M20 genes, normalizes within
    cell, and accumulates donor-by-cell-type means.
    """
    with h5py.File(ATLAS, "r") as handle:
        obs = read_elem(handle["obs"]).reset_index(drop=True)
        var = read_elem(handle["var"])
        genes = sorted(set().union(*[set(program) for program in programs.values()]) & set(var.index))
        positions = var.index.get_indexer(genes)
        if len(genes) < 8 or np.any(positions < 0):
            raise RuntimeError(f"Insufficient or unresolved M8/M20 genes: {len(genes)}")

        keys = obs["sample"].astype(str) + "||" + obs["cell_type"].astype(str)
        codes, levels = pd.factorize(keys, sort=True)
        n_groups = len(levels)
        confidence = np.asarray(pd.to_numeric(obs["cell_type_conf"], errors="coerce").fillna(0) >= 0.90)
        accumulators = {
            "primary": (np.zeros((n_groups, len(genes))), np.zeros(n_groups, dtype=np.int64)),
            "high_annotation_confidence_doublet_proxy": (
                np.zeros((n_groups, len(genes))), np.zeros(n_groups, dtype=np.int64)
            ),
        }
        matrix = handle["layers/counts"]
        matrix_shape = tuple(int(value) for value in matrix.attrs["shape"])
        for row_start in range(0, matrix_shape[0], chunk_size):
            row_end = min(row_start + chunk_size, matrix_shape[0])
            indptr = np.asarray(matrix["indptr"][row_start:row_end + 1], dtype=np.int64)
            data_start, data_end = int(indptr[0]), int(indptr[-1])
            block = sp.csr_matrix(
                (
                    np.asarray(matrix["data"][data_start:data_end]),
                    np.asarray(matrix["indices"][data_start:data_end], dtype=np.int32),
                    indptr - data_start,
                ),
                shape=(row_end - row_start, matrix_shape[1]),
            )
            totals = np.asarray(block.sum(axis=1)).ravel()
            selected = block[:, positions].astype(np.float64).tocsr()
            scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
            selected = sp.diags(scale) @ selected
            selected.data = np.log1p(selected.data)
            block_codes = codes[row_start:row_end]
            block_confidence = confidence[row_start:row_end]
            for sensitivity, keep in (
                ("primary", np.ones(row_end - row_start, dtype=bool)),
                ("high_annotation_confidence_doublet_proxy", block_confidence),
            ):
                group_sums, group_counts = accumulators[sensitivity]
                kept_codes = block_codes[keep]
                if not len(kept_codes):
                    continue
                indicator = sp.csr_matrix(
                    (np.ones(len(kept_codes)), (kept_codes, np.arange(len(kept_codes)))),
                    shape=(n_groups, len(kept_codes)),
                )
                group_sums += (indicator @ selected[keep, :]).toarray()
                group_counts += np.bincount(kept_codes, minlength=n_groups)

    donor_dataset = obs.drop_duplicates("sample").set_index("sample")["dataset"].astype(str)
    outputs: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for sensitivity, (group_sums, group_counts) in accumulators.items():
        eligible = group_counts > 0
        meta = pd.DataFrame({"group_id": levels[eligible], "n_cells": group_counts[eligible]})
        # ``||`` is a literal delimiter.  pandas treats strings as regular
        # expressions by default, where ``||`` matches the empty string.
        meta[["sample", "cell_type"]] = meta["group_id"].str.split("||", n=1, expand=True, regex=False)
        meta["dataset"] = meta["sample"].map(donor_dataset)
        mean = group_sums[eligible, :] / group_counts[eligible, None]
        outputs[sensitivity] = (pd.DataFrame(mean, columns=genes), meta)
    return outputs, genes


def program_scores(mean: pd.DataFrame, programs: dict[str, dict[str, float]], sensitivity: str) -> pd.DataFrame:
    z = (mean - mean.mean(axis=0)) / mean.std(axis=0, ddof=1).replace(0, np.nan)
    output = pd.DataFrame(index=mean.index)
    for uid, weights in programs.items():
        current = {gene: weight for gene, weight in weights.items() if gene in z.columns and (sensitivity != "ambient_gene_removed" or gene not in AMBIENT)}
        finite = [gene for gene in current if np.isfinite(z[gene]).any()]
        mass = sum(current[gene] for gene in finite)
        if len(finite) < 8 or mass < 0.20:
            output[uid] = np.nan
        else:
            w = np.array([current[gene] for gene in finite]); w = w / w.sum()
            output[uid] = np.asarray(z[finite].fillna(0)) @ w
    return output


def contrasts(scores: pd.DataFrame, meta: pd.DataFrame, sensitivity: str) -> list[dict[str, object]]:
    rows = []
    lineages = sorted(set(meta["cell_type"]) - {"Hepatocytes"})
    for uid in scores.columns:
        table = pd.concat([meta.reset_index(drop=True), scores[[uid]].reset_index(drop=True)], axis=1)
        for lineage in lineages:
            wide = table[table["cell_type"].isin(["Hepatocytes", lineage])].pivot_table(index="sample", columns="cell_type", values=uid)
            if not {"Hepatocytes", lineage}.issubset(wide.columns):
                wide = pd.DataFrame(columns=["Hepatocytes", lineage])
            else:
                wide = wide[["Hepatocytes", lineage]].dropna()
            if len(wide) >= 3:
                effect = float(np.mean(wide[lineage] - wide["Hepatocytes"]))
                test = ttest_rel(wide[lineage], wide["Hepatocytes"])
                p = float(test.pvalue)
            else: effect = p = math.nan
            rows.append({"program_uid": uid, "sensitivity": sensitivity, "lineage_vs_hepatocyte": lineage,
                         "estimate": effect, "p": p, "n_paired_donors": len(wide), "direction": "other_higher" if effect > 0 else "hepatocyte_higher"})
    return rows


def main() -> None:
    require_validated_seal()
    programs, registry = load_programs()
    aggregates, genes = aggregate_selected_genes(programs)
    score_rows = []; contrast_rows = []
    for sensitivity in ("primary", "high_annotation_confidence_doublet_proxy"):
        mean, meta = aggregates[sensitivity]
        scores = program_scores(mean, programs, sensitivity)
        for i in range(len(meta)):
            for uid in scores:
                score_rows.append({"program_uid": uid, "sensitivity": sensitivity, **meta.iloc[i].to_dict(), "program_score": scores.iloc[i][uid]})
        contrast_rows.extend(contrasts(scores, meta, sensitivity))
        if sensitivity == "primary":
            ambient_scores = program_scores(mean, programs, "ambient_gene_removed")
            contrast_rows.extend(contrasts(ambient_scores, meta, "ambient_gene_removed"))

    primary_scores = pd.DataFrame(score_rows)
    primary_scores = primary_scores[primary_scores.sensitivity == "primary"]
    loo_rows = []
    for dataset in sorted(primary_scores.dataset.unique()):
        subset_scores = primary_scores[primary_scores.dataset != dataset]
        for uid in (M8, M20):
            for lineage in sorted(set(subset_scores.cell_type) - {"Hepatocytes"}):
                wide = subset_scores[(subset_scores.program_uid == uid) & subset_scores.cell_type.isin(["Hepatocytes", lineage])].pivot_table(index="sample", columns="cell_type", values="program_score")
                if not {"Hepatocytes", lineage}.issubset(wide.columns):
                    wide = pd.DataFrame(columns=["Hepatocytes", lineage])
                else:
                    wide = wide[["Hepatocytes", lineage]].dropna()
                effect = float(np.mean(wide[lineage] - wide["Hepatocytes"])) if len(wide) else math.nan
                loo_rows.append({"program_uid": uid, "held_out_dataset": dataset, "lineage_vs_hepatocyte": lineage,
                                 "estimate": effect, "n_paired_donors": len(wide), "direction": "other_higher" if effect > 0 else "hepatocyte_higher"})

    summary = []
    contrast_frame = pd.DataFrame(contrast_rows)
    contrast_frame["q_within_sensitivity"] = contrast_frame.groupby("sensitivity", group_keys=False)["p"].apply(bh_adjust)
    contrast_rows = contrast_frame.to_dict(orient="records")
    score_frame = pd.DataFrame(score_rows)
    for uid in (M8, M20):
        medians = score_frame[(score_frame.program_uid == uid) & (score_frame.sensitivity == "primary")].groupby("cell_type").program_score.median()
        top = str(medians.idxmax())
        primary = contrast_frame[(contrast_frame.program_uid == uid) & (contrast_frame.sensitivity == "primary")]
        ambient = contrast_frame[(contrast_frame.program_uid == uid) & (contrast_frame.sensitivity == "ambient_gene_removed")]
        conf = contrast_frame[(contrast_frame.program_uid == uid) & (contrast_frame.sensitivity == "high_annotation_confidence_doublet_proxy")]
        if uid == M8:
            proposed = "hepatocyte-cluster stromal ECM program"
            identity_pass = top in {"Fibroblasts", "Endothelial cells"} or "External encapsulating structure" in registry[uid]["module_top_pathway"]
        else:
            ductular = primary[primary.lineage_vs_hepatocyte == "Cholangiocytes"]
            identity_pass = top == "Cholangiocytes" and len(ductular) == 1 and float(ductular.iloc[0].estimate) > 0
            proposed = "ductular injury program" if identity_pass else "BICC1-associated injury program"
        summary.append({"program_uid": uid, "legacy_module": registry[uid]["module"], "original_label": registry[uid]["module_name"],
            "membership_top_pathway": registry[uid]["module_top_pathway"], "top_donor_collapsed_lineage": top,
            "proposed_label": proposed, "identity_support_pass": identity_pass,
            "ambient_sensitivity_available": True, "annotation_confidence_sensitivity_available": True,
            "authoritative_doublet_field_available": False,
            "doublet_limit": "high_cell_type_confidence_is_an_annotation_ambiguity_proxy_not_a_deposited_doublet_call",
            "n_donors": int(score_frame[(score_frame.program_uid == uid) & (score_frame.sensitivity == "primary")]["sample"].nunique())})

    write_tsv(CANDIDATE_ROOT / "lineage/program_identity_scores.tsv", score_rows, list(score_rows[0]))
    write_tsv(CANDIDATE_ROOT / "lineage/program_identity_contrasts.tsv", contrast_rows, list(contrast_rows[0]))
    write_tsv(CANDIDATE_ROOT / "lineage/program_identity_dataset_loo.tsv", loo_rows, list(loo_rows[0]))
    write_tsv(CANDIDATE_ROOT / "program_identity_audit.tsv", summary, list(summary[0]))
    print("PROGRAM_IDENTITY_AUDIT_COMPLETE")


if __name__ == "__main__": main()
