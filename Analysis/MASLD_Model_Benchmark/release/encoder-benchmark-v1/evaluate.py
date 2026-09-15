#!/usr/bin/env python3
"""
Exposure-resolved encoder benchmark.

    python evaluate.py --predictions reference/predictions.npz \
                       --row-contract reference/row_contract.tsv \
                       --exposure reference/exposure_declaration.json \
                       --reference scvi_liver_latent::two_layer_mlp \
                       --out results/

What comes back is never a bare score. For each model it is the paired delta
against the reference you name, restricted to that model's CLEAN held-out
studies, with a 95 percent interval from a donor cluster bootstrap in which
every model is resampled with the SAME donor multiplicities so the deltas are
paired on the same draws.

Models with no clean held-out study get no number. They are written to
not_comparable.tsv, which has no numeric column, and they cannot enter
leaderboard.tsv -- see exposure.py, where that refusal is enforced by type.

The donor is the unit. 50,000 cells come from 102 donors and cells are not
independent.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from exposure import (  # noqa: E402
    ExposureLedger, build_result, dump_results, partition,
    write_leaderboard, write_not_comparable,
)
from metrics import (  # noqa: E402
    BOOTSTRAP_RESAMPLES, BOOTSTRAP_SEED, ROSTER, build_multiplicities,
    endpoints_from_multiplicities, sufficient_statistics,
)

REQUIRED_CONTRACT_COLUMNS = ("row_id", "donor_id", "study_id", "outer_fold", "true_class")


class ContractError(ValueError):
    """Raised when the caller's inputs do not meet the documented contract."""


def read_row_contract(path: pathlib.Path) -> dict[str, np.ndarray]:
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in REQUIRED_CONTRACT_COLUMNS if c not in header]
        if missing:
            raise ContractError(
                f"{path}: row contract is missing {missing}. Required columns are "
                f"{list(REQUIRED_CONTRACT_COLUMNS)}."
            )
        idx = {c: header.index(c) for c in REQUIRED_CONTRACT_COLUMNS}
        cols: dict[str, list] = {c: [] for c in REQUIRED_CONTRACT_COLUMNS}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                raise ContractError(
                    f"{path}: ragged row, {len(parts)} fields against a {len(header)}-field "
                    f"header. An unnamed index column is the usual cause."
                )
            for c, i in idx.items():
                cols[c].append(parts[i])
    out = {c: np.asarray(v) for c, v in cols.items()}
    unknown = sorted(set(out["true_class"]) - set(ROSTER))
    if unknown:
        raise ContractError(f"{path}: true_class values outside the frozen roster: {unknown}")
    return out


def load_predictions(path: pathlib.Path) -> dict[str, np.ndarray]:
    z = np.load(path, allow_pickle=False)
    order = [k for k in z.files if k != "row_id"]
    if "row_id" not in z.files:
        raise ContractError(f"{path}: predictions npz must carry a row_id array")
    return {"__row_id__": z["row_id"].astype(str), **{k: z[k].astype(np.float64) for k in order}}


def align(pred: dict[str, np.ndarray], contract: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Join predictions to the row contract on row_id. A zero or partial join raises."""
    pred_ids = pred["__row_id__"]
    want = contract["row_id"]
    pos = {r: i for i, r in enumerate(pred_ids)}
    missing = [r for r in want if r not in pos]
    if len(missing) == len(want):
        raise ContractError(
            "row_id join produced ZERO matches. The predictions and the row contract "
            "do not describe the same cells."
        )
    if missing:
        raise ContractError(
            f"row_id join is partial: {len(missing)} of {len(want)} contract rows have no "
            f"prediction, e.g. {missing[:3]}. A partial join silently changes the metric."
        )
    take = np.asarray([pos[r] for r in want], dtype=np.int64)
    return {k: (v[take] if k != "__row_id__" else v[take]) for k, v in pred.items()}


def per_study_draws(
    probs: np.ndarray, truth: np.ndarray, donor_index: np.ndarray, study: np.ndarray,
    donors: int, multiplicities: np.ndarray, donor_study: np.ndarray,
    stats: dict | None = None,
) -> tuple[dict[str, float], dict[str, np.ndarray], dict]:
    """Donor sufficient statistics do not depend on the bootstrap seed, so a
    caller sweeping seeds may pass them back in rather than recomputing them."""
    if stats is None:
        predicted = probs.argmax(axis=1)
        stats = sufficient_statistics(truth, predicted, probs, donor_index, donors=donors)
    points: dict[str, float] = {}
    draws: dict[str, np.ndarray] = {}
    for s in sorted(set(study.tolist())):
        sel = np.flatnonzero(donor_study == s)
        roster = stats["present"][sel].any(axis=0)
        identity = np.ones((1, multiplicities.shape[1]), dtype=np.int16)
        points[s] = float(
            endpoints_from_multiplicities(stats, identity, sel, roster)["macro_f1"][0]
        )
        draws[s] = endpoints_from_multiplicities(stats, multiplicities, sel, roster)["macro_f1"]
    return points, draws, stats


def run(predictions, row_contract, exposure_path, reference_id, out, seed,
        stats_cache=None, head_policy="match"):
    """head_policy 'match' compares a model only against a reference carrying the
    SAME head. The whole design of this benchmark is several feature blocks
    through ONE shared head, so a linear-head model ranked against an MLP-head
    reference is measuring the head difference and calling it an encoder
    difference -- which is the exact confusion the reference result is about.
    Models with a different head are SKIPPED and listed, not silently ranked.
    'any' lifts the constraint for callers who know what they are asking for."""
    out = pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    contract = read_row_contract(pathlib.Path(row_contract))
    pred = align(load_predictions(pathlib.Path(predictions)), contract)
    declaration = json.loads(pathlib.Path(exposure_path).read_text(encoding="utf-8"))

    model_ids = [k for k in pred if k != "__row_id__"]
    if reference_id not in model_ids:
        raise ContractError(
            f"reference {reference_id!r} is not among the supplied models: {sorted(model_ids)}"
        )

    donor_labels = sorted(set(contract["donor_id"].tolist()))
    donor_pos = {d: i for i, d in enumerate(donor_labels)}
    donor_index = np.asarray([donor_pos[d] for d in contract["donor_id"]], dtype=np.int64)
    truth = np.asarray([ROSTER.index(c) for c in contract["true_class"]], dtype=np.int64)
    study = contract["study_id"]
    donor_study = np.empty(len(donor_labels), dtype=object)
    for d, s in zip(contract["donor_id"], study):
        donor_study[donor_pos[d]] = s
    donor_study = np.asarray(donor_study.tolist())

    # class-coverage present matrix for the multiplicity builder
    present = np.zeros((len(donor_labels), len(ROSTER)), dtype=bool)
    present[donor_index, truth] = True
    multiplicities = build_multiplicities(
        donor_study, present, n_resamples=BOOTSTRAP_RESAMPLES, seed=seed
    )

    cache = {} if stats_cache is None else stats_cache
    points, draws = {}, {}
    for m in model_ids:
        p, d, st = per_study_draws(
            pred[m], truth, donor_index, study, len(donor_labels), multiplicities,
            donor_study, stats=cache.get(m),
        )
        cache[m] = st
        points[m], draws[m] = p, d
        print(f"  scored {m}", flush=True)

    ref_head = declaration.get(reference_id, {}).get("head_id", "")
    results, skipped = [], []
    for m in model_ids:
        if m == reference_id:
            continue
        m_head = declaration[m].get("head_id", "")
        if head_policy == "match" and m_head != ref_head:
            skipped.append({"model_id": m, "head_id": m_head,
                            "reference_head_id": ref_head,
                            "reason": "head_differs_from_the_reference"})
            continue
        per_study = declaration[m]["per_study_exposure"]
        ledger = ExposureLedger.from_per_study(m, per_study)
        head = declaration[m].get("head_id", "")
        if not ledger.has_clean_study:
            results.append(build_result(
                model_id=m, reference_id=reference_id, head_id=head,
                metric="donor_class_balanced_macro_f1", ledger=ledger))
            continue
        clean = list(ledger.clean)
        mp = float(np.mean([points[m][s] for s in clean]))
        rp = float(np.mean([points[reference_id][s] for s in clean]))
        md = np.mean(np.stack([draws[m][s] for s in clean], axis=1), axis=1)
        rd = np.mean(np.stack([draws[reference_id][s] for s in clean], axis=1), axis=1)
        delta_draws = md - rd
        results.append(build_result(
            model_id=m, reference_id=reference_id, head_id=head,
            metric="donor_class_balanced_macro_f1", ledger=ledger,
            delta=mp - rp,
            lower=float(np.quantile(delta_draws, 0.025, method="linear")),
            median=float(np.quantile(delta_draws, 0.5, method="linear")),
            upper=float(np.quantile(delta_draws, 0.975, method="linear")),
            probability_model_better=float(np.mean(delta_draws > 0.0)),
            n_resamples=BOOTSTRAP_RESAMPLES))

    comparable, withheld = partition(results)
    write_leaderboard(out / "leaderboard.tsv", comparable)
    write_not_comparable(out / "not_comparable.tsv", withheld)
    dump_results(out / "results.json", results)
    with open(out / "skipped_for_head.tsv", "w", encoding="utf-8") as fh:
        fh.write("model_id\thead_id\treference_head_id\treason\n")
        for s_ in skipped:
            fh.write(f"{s_['model_id']}\t{s_['head_id']}\t{s_['reference_head_id']}\t"
                     f"{s_['reason']}\n")
    print(f"\n{len(comparable)} comparable, {len(withheld)} withheld for exposure, "
          f"{len(skipped)} skipped for a differing head")
    for w in withheld:
        print(f"  WITHHELD {w.model_id}: {w.reason}")
    print(f"wrote {out}")
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--row-contract", required=True)
    ap.add_argument("--exposure", required=True)
    ap.add_argument("--reference", required=True,
                    help="model_id of the task's own native baseline")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    ap.add_argument("--head-policy", choices=("match", "any"), default="match",
                    help="'match' (default) compares only same-head pairs")
    a = ap.parse_args()
    run(a.predictions, a.row_contract, a.exposure, a.reference, a.out, a.seed,
        head_policy=a.head_policy)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
