"""Consensus aggregation — 2-of-N rule across model predictions per arm.

Used by every arm's Consensus Aggregator subagent. The aggregator pulls
Peer-Reviewer-validated outputs from Model Runners (which are BLINDED to each
other's outputs) and produces consensus predictions.

Independence of runners is preserved by the team architecture; this module
just enforces the consensus rule mathematically.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Sequence

import numpy as np
import pandas as pd

from output_schema import ConsensusResult, ModelRunOutput


def aggregate_downstream_consensus(
    *,
    arm: str,
    runs: Sequence[ModelRunOutput],
    min_models: int = 2,
    direction_must_agree: bool = True,
) -> list[ConsensusResult]:
    """For D1 / D5: per-target consensus on top-N downstream genes.

    A downstream gene 'passes' if ≥ min_models predict it in top-100 for the same
    target with consistent direction (if direction_must_agree).
    """
    # Index: target_gene → downstream_gene → [model -> (rank, direction, magnitude)]
    idx: dict[tuple[str, str], dict[str, tuple[int, str, float]]] = defaultdict(dict)
    for run in runs:
        for p in run.predictions:
            key = (p.target_gene, p.downstream_gene)
            idx[key][run.model] = (p.abs_rank, p.direction, p.logFC_predicted)

    results = []
    for (target, downstream), per_model in idx.items():
        n_models = len(per_model)
        if n_models < min_models:
            continue
        if direction_must_agree:
            directions = [d for _, d, _ in per_model.values()]
            if len(set(directions)) > 1:
                continue
            consensus_dir = directions[0]
        else:
            consensus_dir = "either"
        magnitudes = [m for _, _, m in per_model.values()]
        consensus_mag = float(np.mean(magnitudes))
        results.append(
            ConsensusResult(
                arm=arm,
                gene=f"{target}->{downstream}",
                n_models_agree=n_models,
                n_models_total=len(runs),
                consensus_direction=consensus_dir,
                consensus_magnitude=consensus_mag,
                pass_consensus=n_models >= min_models,
                contributing_models=list(per_model.keys()),
            )
        )
    return results


def aggregate_reversal_consensus(
    *,
    runs: Sequence[ModelRunOutput],
    top_pct: float = 0.05,
    min_models: int = 2,
) -> pd.DataFrame:
    """For D2: gene must rank top-5% in ≥ min_models.

    Returns DataFrame: gene, n_models_top, contributing_models, mean_rank, mean_score, pass.
    """
    if not runs:
        return pd.DataFrame()
    # Build per-model ranks
    per_model_ranks: dict[str, dict[str, int]] = {}
    per_model_scores: dict[str, dict[str, float]] = {}
    for run in runs:
        m = run.model
        per_model_ranks[m] = {}
        per_model_scores[m] = {}
        for p in run.predictions:
            per_model_ranks[m][p.gene] = p.reversal_rank
            per_model_scores[m][p.gene] = p.reversal_score

    # Identify universe of all genes
    all_genes = set()
    for d in per_model_ranks.values():
        all_genes.update(d.keys())

    rows = []
    for gene in all_genes:
        ranks = []
        scores = []
        contributing = []
        for m, d in per_model_ranks.items():
            r = d.get(gene)
            if r is None:
                continue
            total_in_model = len(d)
            if r <= max(1, int(total_in_model * top_pct)):
                ranks.append(r)
                scores.append(per_model_scores[m].get(gene, np.nan))
                contributing.append(m)
        rows.append(
            dict(
                gene=gene,
                n_models_top=len(ranks),
                contributing_models=";".join(contributing) if contributing else "",
                mean_rank=float(np.mean(ranks)) if ranks else np.nan,
                mean_score=float(np.nanmean(scores)) if scores else np.nan,
                pass_consensus=len(ranks) >= min_models,
            )
        )
    if not rows:
        return pd.DataFrame(
            columns=["gene", "n_models_top", "contributing_models", "mean_rank", "mean_score", "pass_consensus"]
        )
    return pd.DataFrame(rows).sort_values("mean_rank")


def aggregate_synergy_consensus(
    *,
    runs: Sequence[ModelRunOutput],
    sigma_threshold: float = 2.0,
    min_models: int = 2,
) -> pd.DataFrame:
    """For D3: synergy magnitude > sigma_threshold above additive in ≥ min_models."""
    if not runs:
        return pd.DataFrame()
    rows: list[dict] = []
    # Index by sorted gene tuple to handle gene order variation
    pair_acc: dict[tuple[str, ...], list[tuple[str, float, str]]] = defaultdict(list)
    for run in runs:
        for p in run.predictions:
            key = tuple(sorted([g for g in (p.gene1, p.gene2, p.gene3, p.gene4) if g]))
            if p.sigma_above_additive >= sigma_threshold:
                pair_acc[key].append((run.model, p.synergy_magnitude, p.synergy_class))
    for key, hits in pair_acc.items():
        models = [h[0] for h in hits]
        mags = [h[1] for h in hits]
        classes = [h[2] for h in hits]
        # Direction (synergistic vs antagonistic) must agree
        if len(set(classes)) > 1:
            continue
        rows.append(
            dict(
                genes="+".join(key),
                n_genes=len(key),
                synergy_class=classes[0],
                mean_magnitude=float(np.mean(mags)),
                n_models=len(models),
                contributing_models=";".join(models),
                pass_consensus=len(models) >= min_models,
            )
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "genes", "n_genes", "synergy_class", "mean_magnitude",
                "n_models", "contributing_models", "pass_consensus",
            ]
        )
    return pd.DataFrame(rows).sort_values("mean_magnitude", ascending=False)


def aggregate_circuit_consensus(
    *,
    runs: Sequence[ModelRunOutput],
    min_methods: int = 2,
) -> pd.DataFrame:
    """For D4: receiver response replicated in ≥ min_methods propagation pipelines."""
    if not runs:
        return pd.DataFrame()
    circ: dict[tuple[str, str], dict[str, list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for run in runs:
        method = run.model
        for p in run.predictions:
            key = (p.sender_gene, p.receiver_cell_type)
            circ[key][method].extend(p.receiver_response_genes)
    rows = []
    for (sender, receiver), per_method in circ.items():
        methods = list(per_method.keys())
        if len(methods) < min_methods:
            continue
        # Intersection across methods
        gene_lists = [set(genes) for genes in per_method.values()]
        consensus_genes = set.intersection(*gene_lists) if gene_lists else set()
        if not consensus_genes:
            continue
        rows.append(
            dict(
                sender_ligand=sender,
                receiver_cell_type=receiver,
                n_methods=len(methods),
                contributing_methods=";".join(methods),
                consensus_response_genes=";".join(sorted(consensus_genes)),
                n_response_genes=len(consensus_genes),
                pass_consensus=len(methods) >= min_methods,
            )
        )
    return pd.DataFrame(rows)
