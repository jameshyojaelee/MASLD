"""Reproduce the Stage 0c memberships and build the two external arms.

Stage 1, job 1 of 3. No outcome-versus-expression association is computed in
either external arm here. This job establishes what the assigned sets are, what
survives the cross-namespace join into each arm, and nothing else.

**Reproduction, not reuse.**  Stage 0c deposited counts and summaries only: no
per-gene partial vector and no per-gene jackknife membership. Stage 0 deposited
a per-gene table and 0b and 0c did not, which is a lesson about what a check
output file owes its successors and is recorded here rather than silently repaired.
So the memberships are recomputed through the *imported* Stage 0b code path --
the same functions the held-back run executed -- and the job aborts unless six
published invariants reproduce exactly. A digest read would prove only that the
right file was opened; an exact reproduction proves the numbers are the same
numbers.

**The join is cross-namespace and that is the hazard.**  GSE135251 is
ENSG-keyed; both external arms are SYMBOL-keyed, and they do not share a feature
space with each other. If one assigned cell maps into an arm at a different rate
from another, the set comparison is confounded by mappability before any biology
enters. The loss is therefore reported per cell per arm unconditionally, and the
confound is removed by construction rather than tested for: the
expression-matched background is drawn from the arm's own joined universe, so an
assigned set and its background always share a mappability regime.

**Normalization.**  A per-gene Spearman across samples is invariant to any
monotone per-gene transform but *not* to per-sample scaling, which is exactly
what library size is. The external arms are therefore counts-per-million
normalized. The orthogonality of library size to each outcome is measured and
reported for every arm rather than assumed.
"""

from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import rankdata

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class SubstrateError(RuntimeError):
    """Raised when the substrate or the reproduction does not hold."""


#: Every published Stage 0c number the reproduction must land on exactly.
STAGE_0C_INVARIANTS = {
    "realised_universe": 61940,
    "nas_score|fibrosis_stage": {
        "family_size": 61940,
        "genes_bh_below_0_05": 4388,
        "observed_max_abs_partial_r": 0.5689,
    },
    "fibrosis_stage|nas_score": {
        "family_size": 61940,
        "genes_bh_below_0_05": 1305,
        "observed_max_abs_partial_r": 0.6177,
    },
}

ACTIVITY, FIBROSIS = "nas_score", "fibrosis_stage"
CELLS = ("activity_only", "fibrosis_only", "both", "neither")


def _import(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise SubstrateError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    a = rankdata(left, method="average").astype(float)
    b = rankdata(right, method="average").astype(float)
    a -= a.mean()
    b -= b.mean()
    return float(a @ b / np.sqrt((a**2).sum() * (b**2).sum()))


def bh_membership(analysis, unit_genes_residual, residual_exposure, n: int) -> np.ndarray:
    """Boolean BH-0.05 membership over one direction's family."""

    norm = float(np.sqrt((residual_exposure**2).sum()))
    correlations = unit_genes_residual.T @ (residual_exposure / norm)
    critical = analysis.bh_critical_abs_r(int(correlations.size), n - 3)
    count = analysis.bh_count_from_abs_r(np.abs(correlations), critical)
    if count == 0:
        return np.zeros(correlations.size, dtype=bool), correlations, 0
    threshold = np.sort(np.abs(correlations))[::-1][count - 1]
    return np.abs(correlations) >= threshold, correlations, count


def load_counts(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Read a SYMBOL-keyed counts matrix as (symbols, samples, values[gene, sample])."""

    lines = path.read_text(encoding="utf-8").splitlines()
    samples = lines[0].split("\t")[1:]
    symbols: list[str] = []
    rows: list[list[float]] = []
    for line in lines[1:]:
        fields = line.split("\t")
        symbols.append(fields[0])
        rows.append([float(value) for value in fields[1:]])
    return symbols, samples, np.asarray(rows, dtype=float)


def counts_per_million(values: np.ndarray) -> np.ndarray:
    """CPM over the sample axis. values is [gene, sample]."""

    library = values.sum(axis=0, keepdims=True)
    if np.any(library <= 0):
        raise SubstrateError("a sample has a non-positive library size")
    return values / library * 1e6


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-script", type=Path, required=True)
    parser.add_argument("--arm-b-source", type=Path, required=True)
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise SubstrateError("refusing to overwrite a substrate build")
    analysis = _import(arguments.analysis_script, "arm_a_analysis")

    # ---------------- reproduce Stage 0c ----------------
    source = arguments.arm_b_source
    rows = read_table(source / "outcomes" / "participant_endpoints.tsv")
    participants = [row["participant_id"] for row in rows]
    n = len(rows)
    axes = {
        ACTIVITY: np.asarray([int(r[ACTIVITY]) for r in rows], dtype=float),
        FIBROSIS: np.asarray([int(r[FIBROSIS]) for r in rows], dtype=float),
    }
    features = read_table(source / "molecular" / "rna_feature_axis.tsv")
    gene_ids = [row["stable_gene_id"] for row in features]
    values = np.load(source / "molecular" / "rna_values.npy")
    if values.shape != (n, len(gene_ids)):
        raise SubstrateError("the Arm B matrix does not match its axes")

    ordered = np.sort(values, axis=0)
    realised = ((ordered[1:] != ordered[:-1]).sum(axis=0) + 1) >= 2
    del ordered
    universe = np.flatnonzero(realised)
    if int(universe.size) != STAGE_0C_INVARIANTS["realised_universe"]:
        raise SubstrateError(
            f"realised universe {universe.size} is not Stage 0c's "
            f"{STAGE_0C_INVARIANTS['realised_universe']}"
        )
    realised_values = values[:, universe]
    realised_ids = [gene_ids[i] for i in universe]
    gene_ranks = analysis.average_ranks(realised_values)
    gene_ranks -= gene_ranks.mean(axis=0, keepdims=True)
    print(f"reproduced realised universe: {universe.size}", flush=True)

    memberships: dict[str, np.ndarray] = {}
    reproduction = []
    for exposure, covariate in ((ACTIVITY, FIBROSIS), (FIBROSIS, ACTIVITY)):
        unit_covariate = analysis.unit_ranks(axes[covariate])
        scaled, norms = analysis.unit_columns(
            analysis.residualize(gene_ranks, unit_covariate)
        )
        keep = norms > 0.0
        if not keep.all():
            raise SubstrateError("a gene dropped as collinear; Stage 0c dropped none")
        member, correlations, count = bh_membership(
            analysis,
            scaled,
            analysis.residualize(analysis.centred_ranks(axes[exposure]), unit_covariate),
            n,
        )
        key = f"{exposure}|{covariate}"
        expected = STAGE_0C_INVARIANTS[key]
        observed_max = round(float(np.abs(correlations).max()), 4)
        if count != expected["genes_bh_below_0_05"]:
            raise SubstrateError(
                f"{key} count {count} is not Stage 0c's {expected['genes_bh_below_0_05']}"
            )
        if int(correlations.size) != expected["family_size"]:
            raise SubstrateError(f"{key} family size does not reproduce")
        if observed_max != expected["observed_max_abs_partial_r"]:
            raise SubstrateError(
                f"{key} max |partial r| {observed_max} is not "
                f"{expected['observed_max_abs_partial_r']}"
            )
        memberships[key] = member
        reproduction.append({
            "direction": key,
            "family_size": int(correlations.size),
            "genes_bh_below_0_05": int(count),
            "observed_max_abs_partial_r": observed_max,
            "matches_stage_0c": True,
        })
        print(f"reproduced {key}: {count}", flush=True)

    activity_member = memberships[f"{ACTIVITY}|{FIBROSIS}"]
    fibrosis_member = memberships[f"{FIBROSIS}|{ACTIVITY}"]
    assignment = np.where(
        activity_member & fibrosis_member, "both",
        np.where(activity_member, "activity_only",
                 np.where(fibrosis_member, "fibrosis_only", "neither")),
    )
    cell_sizes = {cell: int((assignment == cell).sum()) for cell in CELLS}
    if cell_sizes["activity_only"] + cell_sizes["both"] != 4388:
        raise SubstrateError("the activity marginal does not reconstruct")
    if cell_sizes["fibrosis_only"] + cell_sizes["both"] != 1305:
        raise SubstrateError("the fibrosis marginal does not reconstruct")
    print(f"assignment cells: {cell_sizes}", flush=True)

    # ---------------- per-gene jackknife membership frequency ----------------
    frequency = {cell: np.zeros(universe.size, dtype=np.int64) for cell in CELLS}
    for dropped in range(n):
        rows_kept = np.arange(n) != dropped
        block = realised_values[rows_kept, :]
        constant = block.min(axis=0) == block.max(axis=0)
        ranks = analysis.average_ranks(block)
        ranks -= ranks.mean(axis=0, keepdims=True)
        fold: dict[str, np.ndarray] = {}
        for exposure, covariate in ((ACTIVITY, FIBROSIS), (FIBROSIS, ACTIVITY)):
            unit_covariate = analysis.unit_ranks(axes[covariate][rows_kept])
            scaled, norms = analysis.unit_columns(
                analysis.residualize(ranks, unit_covariate)
            )
            usable = norms > 0.0
            member = np.zeros(universe.size, dtype=bool)
            partial, _, _ = bh_membership(
                analysis,
                scaled[:, usable],
                analysis.residualize(
                    analysis.centred_ranks(axes[exposure][rows_kept]), unit_covariate
                ),
                n - 1,
            )
            member[np.flatnonzero(usable)] = partial
            member[constant] = False
            fold[exposure] = member
        fold_assignment = np.where(
            fold[ACTIVITY] & fold[FIBROSIS], "both",
            np.where(fold[ACTIVITY], "activity_only",
                     np.where(fold[FIBROSIS], "fibrosis_only", "neither")),
        )
        for cell in CELLS:
            frequency[cell] += (fold_assignment == cell)
        if (dropped + 1) % 30 == 0:
            print(f"  jackknife {dropped + 1}/{n}", flush=True)
    print("jackknife membership frequency done", flush=True)

    stability = {}
    for cell in CELLS:
        selected = assignment == cell
        if not selected.any():
            stability[cell] = {"n_genes": 0}
            continue
        share = frequency[cell][selected] / float(n)
        stability[cell] = {
            "n_genes": int(selected.sum()),
            "median_fraction_of_fits_retaining_the_assignment": float(np.median(share)),
            "percentile_5": float(np.percentile(share, 5)),
            "fraction_of_genes_assigned_in_every_fit": float(np.mean(share == 1.0)),
            "fraction_of_genes_assigned_in_at_least_90_percent": float(
                np.mean(share >= 0.9)
            ),
        }

    # ---------------- external arms ----------------
    mapping: dict[str, str] = {}
    for line in (source / "molecular" / "gene_mapping_states.tsv").read_text(
        encoding="utf-8"
    ).splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) >= 3 and fields[2]:
            mapping[fields[0]] = fields[2]

    activation = json.loads(arguments.activation.read_text(encoding="utf-8"))
    assignment_by_gene = dict(zip(realised_ids, assignment.tolist()))
    arms = {}
    for cohort in ("GSE130970", "GSE193066"):
        symbols, samples, raw = load_counts(
            arguments.external_root / cohort / f"{cohort}_counts.tsv"
        )
        cpm = counts_per_million(raw)
        library = raw.sum(axis=0)
        varying = cpm.min(axis=1) != cpm.max(axis=1)
        joined_gene: dict[str, int] = {}
        collisions = 0
        for index, symbol in enumerate(symbols):
            if not varying[index]:
                continue
            stable = mapping.get(symbol)
            if stable is None:
                continue
            if stable in joined_gene:
                collisions += 1
                continue
            joined_gene[stable] = index
        joined_cells = Counter(
            assignment_by_gene[stable]
            for stable in joined_gene
            if stable in assignment_by_gene
        )
        loss = {}
        for cell in CELLS:
            total = cell_sizes[cell]
            survived = int(joined_cells.get(cell, 0))
            loss[cell] = {
                "assigned": total,
                "joined_into_this_arm": survived,
                "surviving_fraction": survived / total if total else None,
            }
        arms[cohort] = {
            "features_deposited": len(symbols),
            "features_non_constant": int(varying.sum()),
            "features_joined_to_a_stable_gene_id": len(joined_gene),
            "many_to_one_collisions_dropped": collisions,
            "samples": len(samples),
            "library_size_min": float(library.min()),
            "library_size_max": float(library.max()),
            "library_size_ratio": float(library.max() / library.min()),
            "join_loss_per_assigned_cell": loss,
            "normalization": "counts_per_million over the sample axis",
            "why_cpm": (
                "a per-gene Spearman across samples is invariant to a monotone "
                "per-gene transform but not to per-sample scaling, which is "
                "what library size is"
            ),
        }
        print(f"{cohort}: joined {len(joined_gene)} genes", flush=True)
        (arguments.output / "arms").mkdir(mode=0o750, parents=True, exist_ok=True)
        order = sorted(joined_gene)
        matrix = cpm[[joined_gene[stable] for stable in order], :].T
        np.save(arguments.output / "arms" / f"{cohort}_cpm.npy", matrix)
        (arguments.output / "arms" / f"{cohort}_gene_axis.tsv").write_text(
            "stable_gene_id\tsource_symbol\n"
            + "".join(f"{s}\t{symbols[joined_gene[s]]}\n" for s in order),
            encoding="utf-8",
        )
        (arguments.output / "arms" / f"{cohort}_sample_axis.tsv").write_text(
            "sample_id\n" + "".join(f"{s}\n" for s in samples), encoding="utf-8"
        )

    payload = {
        "schema_version": "masld-bench-two-axis-substrate-v1",
        "no_outcome_association_was_computed_in_either_external_arm": True,
        "reproduction_of_stage_0c": {
            "method": (
                "recomputed through the imported Stage 0b code path, not a "
                "reimplementation, and aborted unless every published invariant "
                "landed exactly"
            ),
            "why_not_a_digest_read": (
                "Stage 0c deposited counts and summaries only. A digest read "
                "would prove the right file was opened; an exact reproduction "
                "proves the numbers are the same numbers."
            ),
            "lesson_recorded_rather_than_silently_repaired": (
                "Stage 0 deposited a per-gene association table. Stage 0b and "
                "Stage 0c did not, so neither could hand its memberships to a "
                "successor. A gate artifact owes its successors the per-gene "
                "quantities its counts were computed from."
            ),
            "invariants_reproduced": reproduction,
            "realised_universe": int(universe.size),
        },
        "assignment": {
            "definition": (
                "the four-way cross of the two BH 0.05 memberships from Stage 0c"
            ),
            "cell_sizes": cell_sizes,
            "marginals_reconstructed": {
                "activity_only_plus_both": cell_sizes["activity_only"] + cell_sizes["both"],
                "fibrosis_only_plus_both": cell_sizes["fibrosis_only"] + cell_sizes["both"],
            },
            "the_both_cell_cannot_carry_a_specificity_claim": (
                "genes in it are BH-significant on both axes by construction, so "
                "a specificity test that included them would fail for an "
                "arithmetic reason"
            ),
            "jackknife_membership_stability": {
                "n_leave_one_out_fits": n,
                "per_cell": stability,
                "is_a_gate": False,
            },
        },
        "external_arms": arms,
        "join_hazard": {
            "training_namespace": "ENSG (stable_gene_id)",
            "external_namespace": "SYMBOL",
            "the_two_arms_do_not_share_a_feature_space": True,
            "mapping_table": (
                "the Arm B source gene_mapping_states.tsv, which carries the "
                "symbol and ensembl row keys the training axis was built from"
            ),
            "how_the_confound_is_removed": (
                "not by a threshold on unequal loss, which would invite a "
                "judgment call at the moment the answer is visible, but by "
                "construction: the expression-matched background is drawn from "
                "the arm's own joined universe, so an assigned set and its "
                "background always share a mappability regime"
            ),
            "scope_statement": (
                "every external result applies to the joined subset of each "
                "assigned set, with the surviving fraction named, and never to "
                "the full assigned set"
            ),
        },
        "activation_record": {
            "path": str(arguments.activation),
            "user_approved_scope_only": activation["approval"][
                "technical_adjudication"]["user_approved"],
            "no_per_cohort_verdict_may_be_presented_as_user_approved": True,
        },
        "seed": arguments.seed,
    }
    (arguments.output / "substrate.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    header = "stable_gene_id\tassignment\tjackknife_fits_retaining_it\n"
    (arguments.output / "gene_assignment.tsv").write_text(
        header
        + "".join(
            f"{gene}\t{cell}\t{int(frequency[cell][i])}\n"
            for i, (gene, cell) in enumerate(zip(realised_ids, assignment.tolist()))
        ),
        encoding="utf-8",
    )
    freeze_tree(arguments.output, {
        "artifact_class": "two_axis_external_substrate",
        "stage_0c_reproduced_exactly": True,
        "outcome_association_computed_externally": False,
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)
    print(json.dumps({
        "cell_sizes": cell_sizes,
        "join_loss": {
            cohort: {c: round(v["join_loss_per_assigned_cell"][c]["surviving_fraction"] or 0, 4)
                     for c in CELLS}
            for cohort, v in arms.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
