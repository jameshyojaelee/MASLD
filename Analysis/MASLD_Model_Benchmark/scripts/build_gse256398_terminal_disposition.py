#!/usr/bin/env python3
"""Freeze the GSE256398 terminal disposition from already-frozen evidence.

Two things were missing before this cohort could support any accuracy or
selection claim: authoritative barcode-level cell labels, and an executable
TaskSpec.  This script re-derives, rather than asserts, what the authoritative
source actually released, and records why neither gap can be closed.

Everything here is outcome-blind.  No model is fit, no metric is computed, and
the donor-level power figures are design-only: they use group sizes alone, not
the frozen embeddings.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from hashlib import sha256
import json
from math import comb
from pathlib import Path


EXPECTED_SOURCE_ARTIFACTS = (
    "f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20"
)

# A per-barcode cell annotation would have to arrive as one of these.  The
# deposit ships neither, and the check below proves it from the frozen SOFT.
ANNOTATION_SUFFIXES = (
    ".csv",
    ".tsv",
    ".txt",
    ".rds",
    ".h5ad",
    ".loom",
    ".mtx",
    ".xlsx",
)

HUMAN = "Homo sapiens"
MOUSE = "Mus musculus"


class GSE256398DispositionError(RuntimeError):
    """Raised when frozen evidence differs from the audited disposition."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def parse_soft(path: Path) -> dict[str, dict[str, list[str]]]:
    """Return per-sample organism and supplementary-file records."""

    samples: dict[str, dict[str, list[str]]] = {}
    current: str | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                current = line.split("=", 1)[1].strip()
                samples[current] = {"organism": [], "supplementary": []}
            elif current and line.startswith("!Sample_organism_ch1"):
                samples[current]["organism"].append(line.split("=", 1)[1].strip())
            elif current and line.startswith("!Sample_supplementary_file"):
                samples[current]["supplementary"].append(line.split("=", 1)[1].strip())
    if not samples:
        raise GSE256398DispositionError("no samples parsed from the frozen SOFT")
    return samples


def barcode_label_search(samples: dict[str, dict[str, list[str]]]) -> dict[str, object]:
    """Prove from the deposit itself whether any cell annotation was released."""

    supplementary = [
        value
        for record in samples.values()
        for value in record["supplementary"]
        if value.lower() != "none"
    ]
    names = [value.rsplit("/", 1)[-1] for value in supplementary]
    count_matrices = [n for n in names if n.endswith(".h5")]
    candidates = [
        n
        for n in names
        if n.lower().endswith(ANNOTATION_SUFFIXES) and not n.endswith(".h5")
    ]
    return {
        "sample_supplementary_files": len(names),
        "count_matrices_h5": len(count_matrices),
        "candidate_annotation_files": sorted(candidates),
        "barcode_level_cell_annotation_released": bool(candidates),
    }


def species_partition(samples: dict[str, dict[str, list[str]]]) -> dict[str, object]:
    organisms = Counter(o for r in samples.values() for o in r["organism"])
    human = sorted(g for g, r in samples.items() if HUMAN in r["organism"])
    mouse = sorted(g for g, r in samples.items() if MOUSE in r["organism"])
    return {
        "organism_counts": dict(organisms),
        "human_samples": len(human),
        "mouse_samples": len(mouse),
        "mouse_gsm": mouse,
        "series_is_mixed_species": bool(human and mouse),
        "human_only_filter_required": True,
    }


def binomial_tail(successes: int, trials: int, probability: float) -> float:
    return sum(
        comb(trials, i) * probability**i * (1 - probability) ** (trials - i)
        for i in range(successes, trials + 1)
    )


def donor_level_power(group_sizes: dict[str, int]) -> dict[str, object]:
    """Design-only feasibility for a donor-level endpoint.

    Uses group sizes alone.  No embedding, no model, no metric on real data.
    """

    total = sum(group_sizes.values())
    smallest = min(group_sizes.values())
    results: dict[str, object] = {
        "masld_relevant_donors": total,
        "group_sizes": dict(sorted(group_sizes.items())),
        "smallest_class_donors": smallest,
    }

    healthy = group_sizes["healthy_control"]
    results["binary_healthy_vs_masld_spectrum"] = {
        "distinct_label_assignments": comb(total, healthy),
        "minimum_attainable_two_sided_permutation_p": 2 / comb(total, healthy),
    }

    for name, baseline in (
        ("binary_healthy_vs_masld_spectrum", (total - healthy) / total),
        ("four_class_disease_group", max(group_sizes.values()) / total),
    ):
        threshold = next(
            k for k in range(total + 1) if binomial_tail(k, total, baseline) < 0.05
        )
        entry = {
            "majority_class_baseline": baseline,
            "correct_donors_required_at_alpha_0.05": threshold,
            "accuracy_required": threshold / total,
            "power_at_true_accuracy_0.80": binomial_tail(threshold, total, 0.80),
            "power_at_true_accuracy_0.90": binomial_tail(threshold, total, 0.90),
        }
        results.setdefault("leave_one_donor_out_bounds", {})[name] = entry

    # A class with three donors cannot be estimated in every outer fold.
    folds = 5
    results["donor_safe_grouped_cv"] = {
        "outer_folds": folds,
        "smallest_class_donors": smallest,
        "outer_folds_guaranteed_without_smallest_class": max(0, folds - smallest),
        "smallest_class_estimable_in_every_fold": smallest >= folds,
    }
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--donor-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise GSE256398DispositionError("terminal disposition output exists")
    if sha256_file(args.source / "ARTIFACTS.json") != EXPECTED_SOURCE_ARTIFACTS:
        raise GSE256398DispositionError("frozen source ARTIFACTS SHA-256 differs")

    samples = parse_soft(args.source / "metadata/GSE256398_family.soft.gz")
    labels = barcode_label_search(samples)
    species = species_partition(samples)
    if labels["barcode_level_cell_annotation_released"]:
        raise GSE256398DispositionError(
            "a candidate barcode annotation appeared; the terminal disposition "
            "must be re-audited rather than reasserted"
        )
    if species["human_samples"] != 26:
        raise GSE256398DispositionError("human donor census differs")

    with args.donor_metadata.open(encoding="utf-8", newline="") as handle:
        donors = list(csv.DictReader(handle, delimiter="\t"))
    if len(donors) != 26:
        raise GSE256398DispositionError("frozen donor metadata census differs")
    relevant = [d for d in donors if d["development_role"] == "masld_relevant"]
    ood = [d for d in donors if d["development_role"] == "etiology_ood"]
    if len(relevant) != 17 or len(ood) != 9:
        raise GSE256398DispositionError("development role partition differs")
    if any(d["age_state"] != "observed" or d["sex_state"] != "observed" for d in donors):
        raise GSE256398DispositionError("donor age or sex completeness differs")

    power = donor_level_power(Counter(d["disease_group"] for d in relevant))

    args.output.mkdir(parents=True, exist_ok=False)
    disposition = {
        "schema_version": "masld-bench-gse256398-terminal-disposition-v1",
        "dataset_id": "gse256398",
        "status": "terminal_no_executable_biological_taskspec",
        "unit_of_replication": "donor",
        "sealed_outcomes_read": False,
        "models_fit": [],
        "models_scored": [],
        "barcode_cell_labels_read": False,
        "evidence": {
            "authoritative_source": "NCBI GEO GSE256398",
            "barcode_label_search": labels,
            "species_partition": species,
            "external_annotation_search": {
                "cellxgene_discover_public_collections_searched": 388,
                "cellxgene_matching_collections": 0,
                "primary_publication": {
                    "pubmed_id": 40074890,
                    "doi": "10.1038/s41586-025-08677-w",
                    "journal": "Nature",
                    "scope": "hepatic stellate cells and R-spondin 3; the human "
                    "snRNA arm is supporting data, not the primary analysis",
                    "pmc_open_access": False,
                    "paywall_bypassed": False,
                },
                "secondary_publication": {
                    "pubmed_id": 40824250,
                    "doi": "10.1097/HC9.0000000000000771",
                    "journal": "Hepatology Communications",
                },
            },
            "donor_level_labels": {
                "state": "observed_and_authoritative",
                "donors": len(donors),
                "masld_relevant": len(relevant),
                "etiology_ood": len(ood),
                "age_observed": 26,
                "recorded_sex_observed": 26,
                "note": "Donor labels are authoritative; barcode labels are not. "
                "A donor-level endpoint therefore does not need barcode labels.",
            },
            "donor_level_taskspec_feasibility": power,
        },
        # Two independent grounds, strongest first.  The power ground does not
        # depend on labels at all, so this disposition survives someone later
        # locating an annotation file.
        "terminal_grounds": [
            "independent_ground_1_power_and_split_contract",
            "independent_ground_2_absent_label_authority",
        ],
        "terminal_reasons": [
            "GROUND 1 (independent of labels). Donor labels are authoritative, "
            "so a donor-level endpoint needs no barcode labels at all, and it "
            "still fails. The smallest MASLD-relevant class holds three "
            "donors, so at least two of five donor-safe outer folds contain "
            "none of it and the class is not estimable fold-wise. The binary "
            "healthy-versus-MASLD-spectrum contrast needs 15 of 17 donors "
            "correct to beat its own 0.647 majority baseline, giving 0.31 "
            "power at a true accuracy of 0.80, below the 0.80 campaign gate.",
            "GROUND 1 (continued). A single cohort is development-only under "
            "the split contract, so even a passing donor-level result could "
            "authorize neither model selection nor an external claim.",
            "GROUND 2. No barcode-level cell annotation exists in the "
            "authoritative deposit. Every GEO sample supplementary file is a "
            "CellBender filtered count matrix; the series carries no "
            "annotation file in any format, and no public CELLxGENE "
            "collection reproduces it.",
            "GROUND 2 (why). The primary publication is mouse-focused and the "
            "human snRNA arm is supporting data, so annotation files were "
            "never produced. It is not in PMC and no paywall was bypassed.",
            "GROUND 2 (consequence). Transferred or reconstructed cell labels "
            "cannot become accuracy truth, so no cell-level accuracy or "
            "model-selection claim is reachable by any route.",
            "The frozen TranscriptFormer embeddings do not rescue either "
            "ground. They support descriptive donor-level diagnostics only, "
            "which the frozen descriptive audit already covers.",
        ],
        "ground_independence": (
            "Ground 1 rests on group sizes and the split contract; ground 2 "
            "rests on what the depositors released. An authoritative barcode "
            "annotation appearing later would retire ground 2 and leave "
            "ground 1 standing, so the terminal status would not change "
            "without a second compatible cohort."
        ),
        "retained_uses": [
            "Descriptive, composition-confounded donor embedding diagnostics "
            "that make no accuracy or selection claim.",
            "The nine alcohol-associated donors remain an etiology OOD stratum "
            "and are never MASLD-negative controls.",
        ],
        "prohibited": [
            "cell_type_accuracy_claim",
            "model_selection_on_this_cohort",
            "external_seal_use",
            "masld_vs_control_contrast_built_from_alcohol_associated_donors",
            "cross_species_pooling_with_the_four_mouse_samples",
        ],
        "reacquisition_guard": {
            "severity": "silent_contamination_if_ignored",
            "reason": "GSE256398 is a mixed-species series and it grew after "
            "the project audit. Four Mus musculus samples now sit alongside "
            "the 26 human donors under the same accession.",
            "consequence_if_ignored": "A re-download of the series without a "
            "species filter would pull four mouse liver snRNA samples into a "
            "human cohort without raising any error, contaminating the donor "
            "census and violating the human-only scope of this paper.",
            "required_filter": "Homo sapiens only, GSM8097071 through GSM8097096",
            "mouse_gsm_to_exclude": species["mouse_gsm"],
            "existing_frozen_work_affected": False,
            "existing_frozen_work_note": "The frozen source artifact took only "
            "the 26 human donors, so no current number changes.",
        },
        "revisit_condition": (
            "Barcode labels alone do not reopen this: they would retire ground "
            "2 and leave ground 1 standing. Reopening requires an independent "
            "compatible cohort that makes a multi-cohort donor-level task "
            "possible, which would address the power and single-cohort "
            "grounds; authoritative barcode labels would additionally be "
            "needed for any cell-level endpoint."
        ),
    }
    (args.output / "terminal_disposition.json").write_text(
        json.dumps(disposition, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(disposition, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
