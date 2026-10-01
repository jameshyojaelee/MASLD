#!/usr/bin/env python3
"""Render an identity-selected molecular query with all stored predictions.

# KEY MESSAGE: A variant query separates its measured regional effect,
# prediction recipes, and missing target or model support.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from catalog_query import query


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Generate candidate figures on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    check = json.loads(args.check.read_text())
    if check["status"] != "pass" or not check["all_stored_effects_match_source"]:
        raise ValueError("Independent Catalog source reconstruction must pass")
    recorded = json.loads(args.example.read_text())
    result = query(args.db, "variant", recorded["query"], assay="caQTL", limit=100)
    if result != recorded:
        raise ValueError("The independently checked query changed")
    evidence = result["evidence"]
    measured = [r for r in evidence if r["kind"] == "measured"]
    predicted = [r for r in evidence if r["kind"] == "predicted"]
    if len(measured) != 1 or len(predicted) != 12 or result["total_evidence"] != 13:
        raise ValueError("One measurement and all twelve fixed recipes required")
    measurement = measured[0]
    if measurement["se"] is None or measurement["se"] <= 0:
        raise ValueError("Measured source SE unavailable")
    if any(r["se"] is not None or r["units"] != measurement["units"] for r in predicted):
        raise ValueError("Prediction uncertainty or effect units differ")
    target = {o["identifier"] for r in evidence for o in r["objects"] if o["role"] == "measured_target_region"}
    if len(target) != 1:
        raise ValueError("Query does not identify one measured target")
    for r in predicted:
        md = r["metadata"]
        if md["heldout_fold"] != 0 or md["training_folds"] != [1, 2, 3, 4]:
            raise ValueError("Stored development folds differ")
    pools = {"variant": 0, "symmetric": 1, "target": 2}
    predicted.sort(key=lambda r: (r["metadata"]["sequence_length_bp"], pools[r["metadata"]["pooling"]], r["metadata"]["head"]))
    rows = [{"label": "Measured association", "effect": measurement["effect"],
             "se": measurement["se"], "kind": "measured", "evidence_id": measurement["evidence_id"], "model": ""}]
    for r in predicted:
        md = r["metadata"]
        label = f'{md["sequence_length_bp"]:,} bp | {md["pooling"]} | ' + ("ridge" if md["head"] == "ridge" else "shared head")
        rows.append({"label": label, "effect": r["effect"], "se": None,
                     "kind": "predicted", "evidence_id": r["evidence_id"], "model": r["model_version"]})
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6,
        "axes.titlesize": 6, "axes.labelsize": 6, "xtick.labelsize": 6,
        "ytick.labelsize": 6, "legend.fontsize": 6, "font.weight": "normal",
        "axes.titleweight": "normal", "pdf.fonttype": 42, "axes.linewidth": .5,
        "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(5.3, 3.85))
    center, se = measurement["effect"], measurement["se"]
    ax.axvspan(center-se, center+se, color="#9E9E9E", alpha=.18, lw=0)
    ax.axvline(0, color="#9E9E9E", lw=.6, ls="--")
    for i, row in enumerate(rows):
        color = "#9E9E9E" if i == 0 else "#1565C0" if i <= 6 else "#C9265E"
        marker = "D" if i == 0 else "o" if i <= 6 else "s"
        ax.plot(row["effect"], i, marker, color=color, ms=3)
        if row["se"] is not None:
            ax.errorbar(row["effect"], i, xerr=row["se"], color=color, lw=.8, capsize=2)
        ax.text(1.025, i, f'{row["effect"]:.4f}', transform=ax.get_yaxis_transform(), va="center")
    ax.set_yticks(range(len(rows)), [r["label"] for r in rows])
    ax.set_ylim(len(rows)-.4, -.7)
    ax.set_xlabel("Accessibility slope per ALT dosage\n(source normalized units)")
    ax.set_title("Variant query: chr1:906982 C>T (GRCh38)", loc="left", pad=9)
    ax.text(1.025, 1.025, "Effect", transform=ax.transAxes, va="bottom")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.tick_params(length=2, pad=2)
    fig.text(.025, .12, "Measured target: Currin peak33; no gene target established by this query.")
    fig.text(.025, .075, "Gray band: measured slope ± 1 SE. Prediction uncertainty is unavailable.")
    fig.text(.025, .03, "Identity-selected example; all 12 stored recipes; candidate, not adopted.")
    fig.subplots_adjust(left=.37, right=.90, top=.89, bottom=.29)
    fig.savefig(args.out / "Figure6_measured_predicted_variant_query.pdf")
    plt.close(fig)
    with (args.out / "plotted_values.tsv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    (args.out / "query.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    (args.out / "sources.json").write_text(json.dumps({
        "selection": "same_first_deposited_matched_Currin_row_as_checked_Catalog_example; no_effect_size_or_prediction_accuracy_selection",
        "checks": "all thirteen plotted values match the independently checked read-only query",
        "source_hashes": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.example, args.check)},
        "database": str(args.db), "measured_target": sorted(target),
        "statistical_test": "none; one query illustration is not a model accuracy comparison"
    }, indent=2)+"\n")
    (args.out / "CAPTIONS.md").write_text("""# Candidate Figure 6: a measured and predicted regional effect

The figure uses the same first deposited Currin variant with all twelve matched
predictions as the independently checked Catalog extension. This identity-based
selection did not use effect size, accuracy or biological prominence. It is a
query illustration, not evidence that this variant represents the population.

The gray point and band show a marginal accessibility association and its
source standard error (±1 SE). Source-normalized slope per ALT dosage is not a
log2 accessibility change or a causal effect. The study describes 138 donors;
effective per-variant biological n is unavailable. Significance-selected lead
ascertainment and association-estimation uncertainty remain.

Every fixed recipe is displayed, ordered by sequence length, pooling and head
rather than prediction error. Blue circles use 2,048 bp and magenta squares
16,384 bp. All predictions hold out chromosome fold 0 and fit on folds 1–4;
foundation-model exposure remains unresolved. Individual prediction intervals
were not estimated, and measurement SE is not assigned to predictions. Stored
heads were not saved or deployed. This is not the adaptation comparison.

The associated molecular object is Currin peak33. The query does not establish
a gene target, inherited MASLD mechanism, treatment response or personal risk.
The full six-entry interface remains available through catalog_query.py. Old
Figure 6 specimens remain unchanged. This candidate is not adopted.
""")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("db", "example", "check", "out"):
        p.add_argument("--"+name, type=Path, required=True)
    main(p.parse_args())
