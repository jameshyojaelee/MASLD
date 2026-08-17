"""Compatibility audit for the pinned HLiCA external-reference subset."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


class ExternalCompatibilityError(ValueError):
    """Raised when external counts cannot enter the locked comparison."""


FROZEN_BROAD_LABELS = {
    "B cells", "Cholangiocytes", "Circulating NK/NKT", "Endothelial cells",
    "Fibroblasts", "Hepatocytes", "Macrophages", "Mono+mono derived cells",
    "Neutrophils", "Plasma cells", "Resident NK", "T cells", "cDC1s",
    "cDC2s", "pDCs",
}


def read_label_map(path: str | Path) -> dict[str, dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"author_cell_type", "frozen_broad_label", "include_reference", "rationale"}
    if not rows or required - rows[0].keys():
        raise ExternalCompatibilityError("external label map is empty or missing columns")
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        label = row["author_cell_type"]
        if label in result:
            raise ExternalCompatibilityError(f"duplicate external author label: {label}")
        if row["include_reference"] not in {"True", "False"}:
            raise ExternalCompatibilityError(f"invalid include_reference value for {label}")
        if row["include_reference"] == "True" and row["frozen_broad_label"] not in FROZEN_BROAD_LABELS:
            raise ExternalCompatibilityError(f"non-frozen broad label for {label}")
        if row["include_reference"] == "False" and row["frozen_broad_label"] != "Unknown":
            raise ExternalCompatibilityError(f"excluded label {label} must map to Unknown")
        result[label] = row
    return result


def gene_mapping_diagnostics(
    model_genes: Iterable[str], external_ids: Iterable[str], external_names: Iterable[str]
) -> dict[str, Any]:
    model_genes = list(map(str, model_genes))
    ids = list(map(str, external_ids))
    names = list(map(str, external_names))
    if len(ids) != len(names):
        raise ExternalCompatibilityError("external gene IDs and names differ in length")
    if len(model_genes) != len(set(model_genes)):
        raise ExternalCompatibilityError("model gene vocabulary contains duplicates")
    by_id: dict[str, list[int]] = {}
    by_name: dict[str, list[int]] = {}
    for index, (gene_id, gene_name) in enumerate(zip(ids, names)):
        by_id.setdefault(gene_id, []).append(index)
        by_name.setdefault(gene_name, []).append(index)
    indices: list[int] = []
    methods = Counter()
    missing: list[str] = []
    ambiguous: list[str] = []
    for gene in model_genes:
        if len(by_id.get(gene, [])) == 1:
            indices.append(by_id[gene][0])
            methods["external_var_id"] += 1
        elif len(by_name.get(gene, [])) == 1:
            indices.append(by_name[gene][0])
            methods["external_feature_name"] += 1
        elif gene not in by_id and gene not in by_name:
            missing.append(gene)
        else:
            ambiguous.append(gene)
    duplicate_targets = len(indices) - len(set(indices))
    return {
        "model_genes": len(model_genes),
        "mapping_methods": dict(sorted(methods.items())),
        "unique_external_indices": len(set(indices)),
        "resolved_indices": indices,
        "missing_genes": missing,
        "ambiguous_genes": ambiguous,
        "duplicate_target_count": duplicate_targets,
        "eligible": not missing and not ambiguous and duplicate_targets == 0,
    }


def resolve_gene_indices(
    model_genes: Iterable[str], external_ids: Iterable[str], external_names: Iterable[str]
) -> tuple[list[int], dict[str, Any]]:
    diagnostics = gene_mapping_diagnostics(model_genes, external_ids, external_names)
    if not diagnostics["eligible"]:
        raise ExternalCompatibilityError(
            "external gene mapping failed: "
            f"missing_count={len(diagnostics['missing_genes'])}, "
            f"ambiguous_count={len(diagnostics['ambiguous_genes'])}, "
            f"duplicate_target_count={diagnostics['duplicate_target_count']}, "
            f"missing={diagnostics['missing_genes'][:20]}, "
            f"ambiguous={diagnostics['ambiguous_genes'][:20]}"
        )
    indices = diagnostics.pop("resolved_indices")
    return indices, diagnostics


def audit_external_compatibility(
    external_h5ad: str | Path,
    prepared_h5ad: str | Path,
    policy: dict[str, Any],
    label_map_path: str | Path,
) -> dict[str, Any]:
    import anndata as ad
    import numpy as np

    labels = read_label_map(label_map_path)
    external = ad.read_h5ad(external_h5ad, backed="r")
    prepared = ad.read_h5ad(prepared_h5ad, backed="r")
    try:
        if external.raw is None:
            raise ExternalCompatibilityError("external H5AD has no raw count matrix")
        required_obs = {"STUDY", "donor_id", "author_cell_type", "suspension_type", "assay"}
        if required_obs - set(external.obs.columns):
            raise ExternalCompatibilityError("external H5AD lacks required donor/study/label fields")
        allowed_studies = set(policy["external_arm"]["included_studies"])
        overlap_aliases = set(policy["external_arm"]["current_reference_alias_overlap"].values())
        study_mask = external.obs["STUDY"].astype(str).isin(allowed_studies).to_numpy()
        donor_values = external.obs["donor_id"].astype(str)
        overlap_mask = donor_values.isin(overlap_aliases).to_numpy()
        target_mask = study_mask & ~overlap_mask
        target = external.obs.loc[target_mask].copy()
        observed_studies = set(target["STUDY"].astype(str))
        if observed_studies != allowed_studies:
            raise ExternalCompatibilityError("an allowed external study has no retained cells")
        donors = sorted(set(target["donor_id"].astype(str)))
        expected = int(policy["external_arm"]["new_unique_donors"])
        if len(donors) != expected:
            raise ExternalCompatibilityError(
                f"external unique donor count changed: {len(donors)} != {expected}"
            )
        observed_labels = set(target["author_cell_type"].astype(str))
        missing_labels = sorted(observed_labels - labels.keys())
        if missing_labels:
            raise ExternalCompatibilityError(f"unmapped external author labels: {missing_labels}")
        included_labels = {
            name for name, row in labels.items() if row["include_reference"] == "True"
        }
        included_mask = target["author_cell_type"].astype(str).isin(included_labels)
        retained = target.loc[included_mask].copy()
        retained["frozen_broad_label"] = [
            labels[value]["frozen_broad_label"]
            for value in retained["author_cell_type"].astype(str)
        ]
        feature_names = external.raw.var["feature_name"].astype(str)
        gene_summary = gene_mapping_diagnostics(
            prepared.var_names.astype(str), external.raw.var_names.astype(str), feature_names
        )
        indices = gene_summary.pop("resolved_indices")
        group_fields = ["STUDY", "suspension_type", "assay", "frozen_broad_label"]
        groups = (
            retained.groupby(group_fields, observed=True)
            .size().rename("cells").reset_index().sort_values(group_fields)
        )
        donor_cells = retained.groupby("donor_id", observed=True).size().sort_index()
        return {
            "schema_version": "masld-external-compatibility-v1",
            "source_cells_in_allowed_studies": int(study_mask.sum()),
            "overlap_cells_removed": int((study_mask & overlap_mask).sum()),
            "candidate_cells_before_label_filter": int(target_mask.sum()),
            "ambiguous_label_cells_removed": int((~included_mask).sum()),
            "candidate_cells": int(len(retained)),
            "candidate_donors": len(donors),
            "candidate_donor_ids": donors,
            "cells_by_donor": {str(key): int(value) for key, value in donor_cells.items()},
            "cells_by_study": {
                str(key): int(value)
                for key, value in retained.groupby("STUDY", observed=True).size().items()
            },
            "cells_by_frozen_broad_label": {
                str(key): int(value)
                for key, value in retained.groupby("frozen_broad_label", observed=True).size().items()
            },
            "cells_by_study_preparation_assay_label": groups.to_dict(orient="records"),
            "gene_mapping": gene_summary,
            "eligible_for_identical_4000_gene_comparison": bool(gene_summary["eligible"]),
            "external_gene_indices_sha256": (
                __import__("hashlib").sha256(
                    np.asarray(indices, dtype=np.int64).tobytes()
                ).hexdigest()
                if gene_summary["eligible"] else None
            ),
        }
    finally:
        external.file.close()
        prepared.file.close()
