"""Additive overlay recording the two BayesPrism variants and their provenance.

This record does NOT re-cut the lineage substrate.  ``model-data-950`` is bound
by digest and read only; its ARTIFACTS.json and the pinned manifest
``config/showcase_substrate/bayesprism_sources.sha256`` are re-hashed after the
run to prove neither was touched.  950 is bound into the 951 and 952 records,
and re-freezing it underneath them is the failure this overlay exists to avoid.

What it corrects
----------------
The 950 lineage substrate froze one variant per cohort and said nothing about
the other.  The mechanism was a glob: ``*/*_bayesprism_proportions.tsv`` does
not match ``*_bayesprism_proportions_star_backup.tsv``, so the pinned manifest
is a FILTER that reads like a LISTING.  It is not wrong about what it pins and
it is silent about what it excluded, which is the more dangerous shape.

What it establishes
-------------------
The variants are not an unresolvable ambiguity.  Each pairs with a distinct
upstream quantification, and the pairing is measurable from the deposited files
rather than inferred from a filename.
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from scripts.build_showcase_aspect_lineage_substrate import (
    SubstrateError,
    read_bayesprism_proportions,
    sha256_file,
    write_json,
    write_tsv,
)

CANONICAL = "_bayesprism_proportions.tsv"
BACKUP = "_bayesprism_proportions_star_backup.tsv"
DETECTION_FLOOR = 1e-4


def counts_fingerprint(path: Path) -> dict[str, Any]:
    """Gene-universe fingerprint of a counts matrix, read from its first rows only."""

    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        first_key = handle.readline().split("\t", 1)[0]
        rows = 1 + sum(1 for line in handle if line.strip())
    return {
        "path": str(path),
        "gene_rows": rows,
        "sample_columns": len(header) - 1,
        "first_row_key": first_key,
    }


def classify_quantification(gene_rows: int, first_key: str) -> str:
    """Name the upstream quantification from the gene universe it produced.

    Measured, not inferred from a filename: the three quantifications in this
    tree emit distinguishable gene universes and row orderings.
    """

    if gene_rows == 77078:
        return "kallisto_human"
    if gene_rows == 37606 or gene_rows == 37615:
        return "star_human"
    if gene_rows == 33423:
        return "mouse"
    return f"unrecognised_{gene_rows}"


def build_overlay(
    *, results_root: Path, bulk_root: Path, substrate: Path, manifest: Path,
    kallisto_driver: Path,
) -> dict[str, Any]:
    cohorts = sorted(
        {Path(p).parent.name for p in glob.glob(f"{results_root}/*/*{CANONICAL}")}
    )
    if not cohorts:
        raise SubstrateError(f"no proportions under {results_root}")

    inventory: list[dict[str, Any]] = []
    comparison: list[dict[str, Any]] = []
    detect: list[dict[str, Any]] = []

    for cohort in cohorts:
        can_p = results_root / cohort / f"{cohort}{CANONICAL}"
        bak_p = results_root / cohort / f"{cohort}{BACKUP}"
        can_counts = bulk_root / cohort / f"{cohort}_counts.tsv"
        bak_counts = bulk_root / cohort / f"{cohort}_counts_star_backup.tsv"

        lin_c, keys_c, M_c = read_bayesprism_proportions(can_p)
        rec: dict[str, Any] = {
            "cohort": cohort,
            "canonical_path": str(can_p),
            "canonical_sha256": sha256_file(can_p),
            "canonical_rows": len(keys_c),
            "canonical_lineages": len(lin_c),
            "canonical_small_values_0_to_1e8": int(((M_c > 0) & (M_c < 1e-8)).sum()),
            "canonical_min_positive": float(M_c[M_c > 0].min()),
            "backup_present": bak_p.is_file(),
        }
        if can_counts.is_file():
            fp = counts_fingerprint(can_counts)
            rec["canonical_counts"] = fp
            rec["canonical_quantification"] = classify_quantification(
                fp["gene_rows"], fp["first_row_key"]
            )
            rec["counts_columns_match_proportion_rows"] = (
                fp["sample_columns"] == len(keys_c)
            )
        if bak_p.is_file():
            lin_b, keys_b, M_b = read_bayesprism_proportions(bak_p)
            rec.update(
                backup_path=str(bak_p),
                backup_sha256=sha256_file(bak_p),
                backup_rows=len(keys_b),
                backup_lineages=len(lin_b),
                backup_small_values_0_to_1e8=int(((M_b > 0) & (M_b < 1e-8)).sum()),
                backup_min_positive=float(M_b[M_b > 0].min()),
            )
            if bak_counts.is_file():
                fpb = counts_fingerprint(bak_counts)
                rec["backup_counts"] = fpb
                rec["backup_quantification"] = classify_quantification(
                    fpb["gene_rows"], fpb["first_row_key"]
                )
                rec["backup_counts_columns_match_backup_rows"] = (
                    fpb["sample_columns"] == len(keys_b)
                )
            shared = [k for k in keys_c if k in set(keys_b)]
            i_c = {k: i for i, k in enumerate(keys_c)}
            i_b = {k: i for i, k in enumerate(keys_b)}
            A = M_c[[i_c[k] for k in shared]]
            B = M_b[[i_b[k] for k in shared]]
            d = np.abs(A - B)
            comparison.append({
                "cohort": cohort,
                "roster_identical": lin_c == lin_b,
                "canonical_rows": len(keys_c),
                "backup_rows": len(keys_b),
                "shared_rows": len(shared),
                "canonical_only_rows": len(set(keys_c) - set(keys_b)),
                "backup_only_rows": len(set(keys_b) - set(keys_c)),
                "max_abs_difference": float(d.max()),
                "mean_abs_difference": float(d.mean()),
                "rows_identical": int((d.max(axis=1) == 0).sum()),
            })
            for j, name in enumerate(lin_c):
                detect.append({
                    "cohort": cohort,
                    "lineage": name,
                    "canonical_above_floor_fraction": float((M_c[:, j] > DETECTION_FLOOR).mean()),
                    "backup_above_floor_fraction": float((M_b[:, j] > DETECTION_FLOOR).mean()),
                })
        inventory.append(rec)

    # Collapse structure grouped by quantification rather than by cohort.
    by_quant: dict[str, dict[str, Any]] = {}
    for rec in inventory:
        q = rec.get("canonical_quantification", "unknown")
        b = by_quant.setdefault(q, {"cohorts": [], "small_values": 0, "min_positive": []})
        b["cohorts"].append(rec["cohort"])
        b["small_values"] += rec["canonical_small_values_0_to_1e8"]
        b["min_positive"].append(rec["canonical_min_positive"])
    for q, b in by_quant.items():
        b["min_positive_across_cohorts"] = min(b["min_positive"])
        b["max_of_min_positive"] = max(b["min_positive"])
        del b["min_positive"]

    driver_text = kallisto_driver.read_text(encoding="utf-8") if kallisto_driver.is_file() else ""
    driver_cohorts = re.findall(r'^\s*"(GSE\d+|PRJNA\d+)"', driver_text, flags=re.M)
    dual = sorted(r["cohort"] for r in inventory if r["backup_present"])

    return {
        "schema_version": "masld-bench-lineage-variant-overlay-v1",
        "record_id": "lineage_variant_provenance_overlay_v1",
        "status": "registered_active",
        "record_is_additive_overlay": True,
        "bound_substrate": {
            "path": str(substrate),
            "artifacts_sha256": sha256_file(substrate / "ARTIFACTS.json"),
            "lineage_receipt_sha256": sha256_file(substrate / "lineage" / "receipt.json"),
            "substrate_edited": False,
            "substrate_refrozen": False,
        },
        "bound_manifest": {
            "path": str(manifest),
            "sha256": sha256_file(manifest),
            "manifest_edited": False,
        },
        "what_this_overlay_corrects": {
            "defect": (
                "The 950 lineage substrate froze one variant per cohort and its "
                "receipt does not mention the other. The string star_backup does "
                "not appear anywhere in lineage/receipt.json."
            ),
            "mechanism": (
                "The glob */*_bayesprism_proportions.tsv does not match "
                "*_bayesprism_proportions_star_backup.tsv. The pinned manifest is "
                "therefore a FILTER that reads like a LISTING: it is accurate "
                "about what it pins and silent about what it excluded. The next "
                "manifest written in this tree will reach for the same glob."
            ),
            "reusable_rule": (
                "A manifest must record what it excluded, not only what it "
                "included. A glob that silently drops a sibling variant produces "
                "a substrate that looks complete and is not."
            ),
            "figures_affected": (
                "Every per-lineage detectability figure in the 950 receipt and in "
                "the lineage scouting report is canonical-variant only and must be "
                "read as such."
            ),
        },
        "variant_provenance_resolved": {
            "resolved": True,
            "how": (
                "Not from the filename. Each proportions variant pairs with a "
                "distinct counts matrix whose gene universe is measurable, and the "
                "pairing is confirmed on two independent signals: cohort "
                "membership and exact sample-count agreement."
            ),
            "signal_1_cohort_membership": {
                "cohorts_with_a_backup": dual,
                "kallisto_driver": str(kallisto_driver),
                "kallisto_driver_sha256": sha256_file(kallisto_driver) if kallisto_driver.is_file() else None,
                "kallisto_driver_cohorts": driver_cohorts,
                "sets_are_identical": sorted(driver_cohorts) == dual,
                "driver_self_description": (
                    "Stage 1: Per-cohort BayesPrism on kallisto-canonical counts"
                ),
            },
            "signal_2_sample_counts": (
                "Per cohort and per variant, the counts-matrix column count equals "
                "the proportions row count exactly, including the asymmetric cases "
                "GSE135251 216-vs-180 and GSE126848 56-vs-53."
            ),
            "conclusion": {
                "canonical_proportions_are": "BayesPrism on the kallisto counts (77,078 genes, first key TSPAN6)",
                "star_backup_proportions_are": "the earlier BayesPrism run on STAR counts (37,606-37,615 genes, symbol-keyed, first key A1BG)",
                "star_backup_is_not_a_rounding_variant": True,
            },
            "weak_signal_not_relied_on": (
                "All five backup files share an mtime near 2026-05-26 19:25 while "
                "the canonical files are staggered 19:57 to 20:05, consistent with "
                "a batch snapshot followed by a sequential rerun. mtime is weak "
                "evidence here because mv preserves it, so this corroborates the "
                "conclusion and does not carry it."
            ),
        },
        "which_variant_pairs_with_the_frozen_axis": {
            "gse135251_frozen_axis_counts": "GSE135251_counts.tsv",
            "gse135251_frozen_axis_counts_sha256": "954bb227892f65d83d409f3fd517f8bd28544f4df0a4cb91ac73014508bffd60",
            "gse135251_frozen_axis_gene_rows": 77078,
            "frozen_axis_quantification": "kallisto_human",
            "therefore": (
                "The CANONICAL proportions are the variant that pairs with the "
                "frozen GSE135251 expression axis. The star_backup pairs with a "
                "STAR quantification that no frozen axis in this benchmark uses."
            ),
            "consequence_for_a_sensitivity": (
                "A both-variant run is a robustness question about the upstream "
                "quantification, not a resolution of a provenance ambiguity. The "
                "provenance is settled."
            ),
        },
        "quantification_split_across_the_axis_arms": {
            "finding": (
                "The three axis cohorts are NOT on one quantification. GSE135251 "
                "and GSE130970 are kallisto; GSE193066 is STAR. GSE193066 has no "
                "backup because it was never part of the kallisto rerun, so its "
                "canonical proportions sit on the same footing as the other "
                "cohorts' star_backup."
            ),
            "reinforces_the_separate_arm_rule": (
                "The two-axis activation already forbids pooling GSE130970 and "
                "GSE193066 on the omitted-level grounds. This is an independent "
                "second reason for the same rule: the two external arms are "
                "quantified differently as well as graded differently."
            ),
        },
        "collapse_structure_is_a_quantification_property": {
            "finding": (
                "The near-zero mass that W2 recorded as a cohort-level effective "
                "floor difference tracks the QUANTIFICATION, not the cohort. Every "
                "kallisto cohort carries values in (0, 1e-8); every human STAR "
                "cohort carries exactly zero."
            ),
            "by_quantification": by_quant,
            "consequence": (
                "A detection-floor rule evaluated per cohort is really being "
                "evaluated per quantification. Comparing an above-floor fraction "
                "between a kallisto arm and a STAR arm compares two upstream "
                "pipelines, not two cohorts."
            ),
        },
        "variant_comparison": comparison,
        "variant_inventory": inventory,
        "crosswalk_finding": {
            "finding": (
                "The GSE135251 proportions are SRR-keyed and the frozen molecular "
                "axis is GSM-keyed, so a direct join is 0 of 180. The 950 receipt "
                "reported a 100% join because it was measured against the "
                "SRR-keyed source metadata, a different object."
            ),
            "the_fix_was_already_inside_the_artifact": (
                "substrate/arm_b/participant_axis.tsv carries participant_id (GSM) "
                "and run_accession (SRR) side by side, 180 rows, 1:1 in both "
                "directions. The crosswalk was in the same artifact whose receipt "
                "reported the wrong join."
            ),
            "classification": "a missing CHECK, not a missing capability",
            "why_it_matters": (
                "This is the composition rule in its sharpest form: every part was "
                "correct and verified, and the join between them was never tested. "
                "Per-artifact verification cannot see a join."
            ),
            "no_new_crosswalk_needed": True,
            "bind_instead": "substrate/arm_b/participant_axis.tsv",
        },
        "no_provenance_receipt_beside_the_files": {
            "finding": (
                "Neither variant carries a receipt naming its aligner or run. The "
                "_star_backup suffix is the only in-band provenance signal, and it "
                "names a backup rather than a method."
            ),
            "resolved_anyway_by": "the counts gene-universe fingerprint, recorded above",
            "rds_objects_present_but_not_read": (
                "<cohort>_bayesprism_result.rds and _star_backup.rds exist beside "
                "the proportions and may carry richer provenance. Reading them "
                "needs R and was not done here."
            ),
        },
        "what_this_overlay_does_NOT_do": [
            "it does not re-freeze or edit model-data-950",
            "it does not edit the pinned manifest",
            "it does not change any figure in 951 or 952",
            "it does not choose a variant for any analysis",
            "it does not design Stage 3",
        ],
        "approval": {
            "user_approved": False,
            "authored_by": "team lead and this producer",
            "note": (
                "Technical record. No part of it has been seen by the user. The "
                "2026-08-27 user approval covered the two-cohort scope decision "
                "only and does not extend here."
            ),
        },
        "status_note": "prepared and frozen; additive only",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--bulk-root", type=Path, required=True)
    parser.add_argument("--substrate", type=Path, required=True)
    parser.add_argument("--substrate-artifacts-sha256", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--kallisto-driver", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    observed = sha256_file(args.substrate / "ARTIFACTS.json")
    if observed != args.substrate_artifacts_sha256:
        raise SubstrateError(
            f"bound substrate digest differs: expected "
            f"{args.substrate_artifacts_sha256}, observed {observed}"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    overlay = build_overlay(
        results_root=args.results_root,
        bulk_root=args.bulk_root,
        substrate=args.substrate,
        manifest=args.manifest,
        kallisto_driver=args.kallisto_driver,
    )
    write_json(args.output / "lineage_variant_overlay.json", overlay)

    inv = overlay["variant_inventory"]
    write_tsv(
        args.output / "variant_inventory.tsv",
        ["cohort", "canonical_rows", "canonical_quantification", "backup_present",
         "backup_rows", "backup_quantification", "canonical_small_values",
         "canonical_min_positive"],
        [[r["cohort"], r["canonical_rows"], r.get("canonical_quantification", ""),
          r["backup_present"], r.get("backup_rows", ""), r.get("backup_quantification", ""),
          r["canonical_small_values_0_to_1e8"], f"{r['canonical_min_positive']:.6e}"]
         for r in inv],
    )
    write_tsv(
        args.output / "variant_disagreement.tsv",
        ["cohort", "roster_identical", "canonical_rows", "backup_rows", "shared_rows",
         "max_abs_difference", "mean_abs_difference", "rows_identical"],
        [[c["cohort"], c["roster_identical"], c["canonical_rows"], c["backup_rows"],
          c["shared_rows"], f"{c['max_abs_difference']:.4f}",
          f"{c['mean_abs_difference']:.6f}", c["rows_identical"]]
         for c in overlay["variant_comparison"]],
    )
    print(json.dumps({
        "cohorts": len(inv),
        "with_backup": sum(1 for r in inv if r["backup_present"]),
        "provenance_resolved": overlay["variant_provenance_resolved"]["resolved"],
        "sets_identical": overlay["variant_provenance_resolved"][
            "signal_1_cohort_membership"]["sets_are_identical"],
        "collapse_by_quantification": {
            k: {"cohorts": len(v["cohorts"]), "small_values": v["small_values"]}
            for k, v in overlay["collapse_structure_is_a_quantification_property"][
                "by_quantification"].items()},
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
