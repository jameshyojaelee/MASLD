#!/usr/bin/env python3
"""Read identities and recount fixed ATAC peaks without changing named inputs.

The heavy ``pilot`` and ``count`` commands require a Slurm allocation. Counts
use the same SnapATAC2 fragment strategy as ATAC Context v3: one deduplicated
fragment record per overlapping peak, irrespective of fragment multiplicity.
``audit`` reads HDF5 observation metadata only, including the RNA run/library
identity problem; it never reads a gene-expression matrix.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite

ROOT = Path(__file__).resolve().parents[4]
V3 = ROOT / "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3"
sys.path.insert(0, str(V3))
from atac_context_v3_lib import LABEL_MAP  # noqa: E402

COHORTS = {
    "GSE244832": {
        "h5ad": "Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_relabeled.h5ad",
        "metadata": "data/GSE244832/metadata/donor_pairing.csv",
        "per_donor": "Analysis/ATAC/Human_Multiome/results/snapatac2/per_donor",
    },
    "GSE281367": {
        "h5ad": "Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad",
        "metadata": "data/GSE281367/metadata/donor_pairing.csv",
        "per_donor": "Analysis/ATAC/Human_External/snapatac2_fast/per_donor",
    },
}


def configuration(path: Path) -> dict:
    with path.open() as handle:
        if path.suffix.lower() == ".json":
            config = json.load(handle)
        else:
            import yaml
            config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("configuration must be a mapping")
    return config


def project_path(config: dict, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else Path(config.get("project_root", ROOT)) / path


def cohort_spec(config: dict, cohort: str) -> dict[str, Path]:
    spec = dict(COHORTS[cohort])
    overrides = config.get("cohorts", {}).get(cohort, {})
    if isinstance(overrides, dict):
        spec.update({key: value for key, value in overrides.items() if key in spec})
    return {key: project_path(config, str(value)) for key, value in spec.items()}


def obs_column(group: h5py.Group, key: str) -> np.ndarray:
    node = group[key]
    if isinstance(node, h5py.Dataset):
        return node.asstr()[:] if node.dtype.kind in "OSU" else node[:]
    categories = node["categories"].asstr()[:]
    codes = node["codes"][:]
    return np.asarray([categories[index] if index >= 0 else "" for index in codes])


def observations(path: Path, columns: tuple[str, ...]) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        index_key = obs.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        result = pd.DataFrame({key: obs_column(obs, key) for key in columns})
        result.insert(0, "raw_barcode", obs_column(obs, str(index_key)))
    return result


def donor_metadata(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if frame["donor_id"].duplicated().any():
        raise ValueError("duplicate donor identity in metadata")
    frame["condition"] = frame["condition"].str.upper()
    if not frame["condition"].isin(["NORMAL", "MASL", "MASH"]).all():
        raise ValueError("unrecognized source disease condition")
    return frame


def cell_roster(config: dict, cohort: str) -> pd.DataFrame:
    spec = cohort_spec(config, cohort)
    obs = observations(spec["h5ad"], ("donor_id", "cell_type"))
    meta = donor_metadata(spec["metadata"])
    if obs.duplicated(["donor_id", "raw_barcode"]).any():
        raise ValueError(f"duplicate donor/barcode identities in {cohort}")
    obs = obs.merge(meta[["donor_id", "condition"]], on="donor_id", how="left", validate="many_to_one", sort=False)
    if obs["condition"].isna().any():
        raise ValueError("cells lack authoritative donor condition")
    obs["cohort"] = cohort
    obs["lineage"] = obs["cell_type"].map(LABEL_MAP[cohort]).fillna("excluded")
    obs["cell_id"] = cohort + ":" + obs["donor_id"] + ":" + obs["raw_barcode"]
    sizes = obs.groupby(["donor_id", "lineage"])["cell_id"].transform("size")
    obs["donor_lineage_cells"] = sizes
    obs["da_eligible"] = (sizes >= 20) & obs["condition"].isin(["NORMAL", "MASH"]) & (obs["lineage"] != "excluded")
    network_min = int(config.get("network", {}).get("min_cells_per_donor", config.get("network_min_donor_cells", 100)))
    obs["network_eligible"] = (sizes >= network_min) & (obs["lineage"] != "excluded")
    return obs


def pilot_indices(donors: np.ndarray, n_cells: int, seed: int, min_per_donor: int = 100) -> np.ndarray:
    """Balanced, deterministic sampling without consulting disease labels."""
    if n_cells <= 0:
        raise ValueError("pilot cell limit must be positive")
    rng = np.random.default_rng(seed)
    eligible = [donor for donor in sorted(set(donors)) if np.count_nonzero(donors == donor) >= min_per_donor]
    max_donors = n_cells // min_per_donor
    if not eligible or max_donors < 1:
        raise ValueError("pilot cannot retain minimum cells in any donor")
    if len(eligible) > max_donors:
        eligible = sorted(rng.choice(eligible, max_donors, replace=False).tolist())
    pools = [rng.permutation(np.flatnonzero(donors == donor)).tolist() for donor in eligible]
    target = min(n_cells, sum(map(len, pools)))
    selected = []
    while len(selected) < target:
        for pool in pools:
            if pool and len(selected) < target:
                selected.append(pool.pop())
    return np.asarray(sorted(selected), dtype=np.int64)


def aggregate_counts(matrix: sp.spmatrix, donor_codes: np.ndarray, n_donors: int) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    matrix = matrix.tocsr(copy=True)
    if matrix.shape[0] != len(donor_codes):
        raise ValueError("matrix and donor roster differ")
    if np.any(matrix.data < 0) or np.any(matrix.data != np.rint(matrix.data)):
        raise ValueError("fragment counts must be nonnegative integers")
    matrix.sum_duplicates()
    membership = sp.csr_matrix((np.ones(len(donor_codes), dtype=np.int64), (donor_codes, np.arange(len(donor_codes)))), shape=(n_donors, len(donor_codes)))
    raw = (membership @ matrix.astype(np.int64)).tocsr()
    binary = matrix.copy().astype(np.int64)
    binary.eliminate_zeros()
    binary.data[:] = 1
    detected = (membership @ binary).tocsr()
    return raw, detected


def balanced_cell_mask(roster: pd.DataFrame, seed: int) -> tuple[np.ndarray, int]:
    """Equal nuclei per primary eligible donor; never thin peak entries."""
    eligible = roster["da_eligible"].astype(bool).to_numpy()
    sizes = roster.loc[eligible].groupby("donor_id").size()
    if sizes.empty:
        return np.zeros(len(roster), dtype=bool), 0
    target = int(sizes.min())
    if target < 20:
        raise ValueError("balanced primary donor has fewer than20nuclei")
    rng = np.random.default_rng(seed)
    selected = np.zeros(len(roster), dtype=bool)
    values = roster["donor_id"].to_numpy()
    for donor in sorted(sizes.index):
        indices = np.flatnonzero(eligible & (values == donor))
        selected[rng.choice(indices, target, replace=False)] = True
    return selected, target


def audit_rna(config: dict) -> tuple[pd.DataFrame, dict]:
    path = project_path(config, str(config.get("rna_reference", "Analysis/ATAC/Human_Multiome/results/label_transfer/rna_reference_GSE244832.h5ad")))
    obs = observations(path, ("sample", "condition"))
    meta = donor_metadata(cohort_spec(config, "GSE244832")["metadata"])
    lookup = {}
    for row in meta.to_dict("records"):
        for run in row["rna_srrs"].split(";"):
            if run in lookup:
                raise ValueError("RNA run maps to multiple source libraries")
            lookup[run] = (row["donor_id"], row["rna_gsm"], row["condition"])
    rows = []
    for row in obs.itertuples(index=False):
        sample = str(row.sample)
        identity = lookup.get(sample)
        prefix = sample + "_"
        barcode = row.raw_barcode[len(prefix):] if row.raw_barcode.startswith(prefix) else ""
        rows.append({"source_cell_id": row.raw_barcode, "rna_run": sample, "raw_barcode": barcode,
                     "donor_id": identity[0] if identity else "", "rna_library_gsm": identity[1] if identity else "",
                     "source_condition": row.condition, "donor_condition": identity[2] if identity else "",
                     "identity_resolved": bool(identity and barcode)})
    table = pd.DataFrame(rows)
    valid = table["identity_resolved"]
    duplicate = pd.Series(False, index=table.index)
    duplicate.loc[valid] = table.loc[valid].duplicated(["rna_library_gsm", "raw_barcode"], keep=False)
    table["duplicate_library_barcode"] = duplicate
    unresolved = int((~valid).sum())
    summary = {"source": str(path), "rows": len(table), "unresolved_identity_rows": unresolved,
               "duplicate_library_barcode_rows": int(duplicate.sum()),
               "unique_resolved_library_barcodes": int(table.loc[valid, ["rna_library_gsm", "raw_barcode"]].drop_duplicates().shape[0]),
               "source_conditions": sorted(obs["condition"].unique().tolist()),
               "expression_matrix_read": False,
               "rna_model_eligible": bool(not unresolved and not duplicate.any()),
               "next_action": "compare count vectors within duplicated library/barcode groups on compute; do not sum or drop automatically" if duplicate.any() else "verify raw-count and biological-library eligibility before modelling"}
    return table, summary


def write_json(path: Path, value: dict) -> None:
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def check_output_root(config: dict, output: Path) -> None:
    root = Path(config.get("project_root", ROOT)).resolve()
    protected = [root / "docs", root / "figures/main", project_path(config, str(config.get("v3_candidate_root", "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1")))]
    for cohort in COHORTS:
        spec = cohort_spec(config, cohort)
        protected.extend([spec["h5ad"].parent, spec["per_donor"], spec["metadata"].parent])
    destination = output.resolve()
    for path in protected:
        source = path.resolve()
        if destination == source or source in destination.parents:
            raise ValueError(f"output must not modify a protected source: {source}")


def verify_fragment_records(config: dict, output: Path, cohort: str) -> None:
    """Independent indexed-fragment check of one nonzero peak per donor."""
    if not os.environ.get("SLURM_JOB_ID") or not os.environ.get("SLURM_JOB_NODELIST"):
        raise RuntimeError("direct fragment checks require a compute allocation")
    import pysam
    cells = pd.read_csv(output / "cells.tsv", sep="\t", dtype=str)
    donors = pd.read_csv(output / "donors.tsv", sep="\t", dtype=str)
    peaks = pd.read_csv(output / "peaks.tsv", sep="\t")
    with gzip.open(output / "donor_counts.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsr()
    if matrix.shape != (len(donors), len(peaks)):
        raise ValueError("direct check input axes differ")
    rng = np.random.default_rng(int(config.get("seed", 42)))
    checks = []
    override = config.get("cohorts", {}).get(cohort, {})
    template = override.get("fragment_template") if isinstance(override, dict) else None
    for index, donor in enumerate(donors["donor_id"]):
        nonzero = matrix.getrow(index).indices
        if not len(nonzero):
            continue
        peak_index = int(rng.choice(nonzero))
        peak = peaks.iloc[peak_index]
        barcodes = set(cells.loc[cells["donor_id"] == donor, "raw_barcode"])
        if template:
            source = project_path(config, template.format(donor=donor))
        elif cohort == "GSE244832":
            source = project_path(config, f"Analysis/ATAC/Human_Multiome/results/fragments/{donor}_fragments.tsv.gz")
        else:
            source = project_path(config, f"Analysis/ATAC/Human_External/cellranger/{donor}/outs/fragments.tsv.gz")
        records = 0
        with pysam.TabixFile(str(source)) as tabix:
            for record in tabix.fetch(str(peak["chrom"]), int(peak["start0"]), int(peak["end"])):
                fields = record.split("\t")
                if int(fields[1]) < int(peak["end"]) and int(fields[2]) > int(peak["start0"]) and fields[3] in barcodes:
                    records += 1
        expected = int(matrix[index, peak_index])
        checks.append({"donor_id": donor, "peak_coordinate": str(peak["peak_coordinate"]), "matrix_count": expected, "direct_fragment_records": records, "matches": expected == records})
    result = {"cohort": cohort, "counting_unit": "deduplicated_fragment_record", "checks": checks,
              "all_match": bool(checks and all(row["matches"] for row in checks))}
    write_json(output / "direct_fragment_recount.json", result)
    if not result["all_match"]:
        raise ValueError("independent fragment record recount failed; diagnostic retained")


def audit_rna_vectors(config: dict, output: Path) -> None:
    """Compare thirty duplicated library/barcode groups without deduplicating."""
    if not os.environ.get("SLURM_JOB_ID") or not os.environ.get("SLURM_JOB_NODELIST"):
        raise RuntimeError("RNA stored-vector comparison requires a compute allocation")
    identity = pd.read_csv(output / "GSE244832.rna_identity.tsv.gz", sep="\t", dtype=str)
    identity["source_row"] = np.arange(len(identity))
    duplicate = identity["duplicate_library_barcode"].str.lower().eq("true")
    candidates = identity.loc[duplicate, ["rna_library_gsm", "raw_barcode"]].drop_duplicates()
    rng = np.random.default_rng(int(config.get("seed", 42)))
    chosen = candidates.iloc[np.sort(rng.choice(len(candidates), min(30, len(candidates)), replace=False))]
    selected = identity.merge(chosen, on=["rna_library_gsm", "raw_barcode"], how="inner", validate="many_to_one")
    del identity
    path = project_path(config, str(config.get("rna_reference", "Analysis/ATAC/Human_Multiome/results/label_transfer/rna_reference_GSE244832.h5ad")))
    checks = []
    with h5py.File(path, "r") as handle:
        matrix = handle["X"]
        encoding = matrix.attrs.get("encoding-type", "dense")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        header = {"X_encoding": str(encoding), "X_shape": list(map(int, matrix.attrs["shape"] if isinstance(matrix, h5py.Group) else matrix.shape)),
                  "X_dtype": str(matrix["data"].dtype if isinstance(matrix, h5py.Group) else matrix.dtype),
                  "layers_available": list(handle["layers"].keys()) if "layers" in handle else [], "raw_available": "raw" in handle,
                  "normalization_documented": False}
        if isinstance(matrix, h5py.Group) and encoding != "csr_matrix":
            raise ValueError("bounded RNA comparison requires row-sparse CSR or dense source")
        for (library, barcode), rows in selected.groupby(["rna_library_gsm", "raw_barcode"], sort=True):
            vectors = []
            for index in rows["source_row"].astype(int):
                if isinstance(matrix, h5py.Group):
                    first, last = map(int, matrix["indptr"][index:index + 2])
                    vector = sp.csr_matrix((matrix["data"][first:last], matrix["indices"][first:last], np.array([0, last - first])), shape=(1, header["X_shape"][1]))
                    vector.sum_duplicates()
                    vector.sort_indices()
                else:
                    vector = sp.csr_matrix(matrix[index:index + 1, :])
                vectors.append(vector)
            equal = all((vectors[0] != other).nnz == 0 for other in vectors[1:])
            values = np.concatenate([vector.data for vector in vectors])
            checks.append({"rna_library_gsm": library, "raw_barcode": barcode, "stored_rows": len(rows), "rna_runs": sorted(rows["rna_run"].unique()),
                           "nonzero_features_per_stored_row": [int(vector.nnz) for vector in vectors],
                           "all_stored_vectors_identical": equal,
                           "all_sampled_values_nonnegative_integers": bool(np.all(np.isfinite(values)) and np.all(values >= 0) and np.all(values == np.rint(values)))})
    all_equal = bool(checks and all(check["all_stored_vectors_identical"] for check in checks))
    write_json(output / "rna_duplicate_vector_check.json", {"source": str(path), "seed": int(config.get("seed", 42)), "matrix_header": header,
               "groups_checked": len(checks), "all_groups_identical": all_equal, "checks": checks,
               "interpretation": "sampled duplicated library/barcode rows repeat identical stored expression vectors" if all_equal else "at least one sampled group differs; biological-library reconstruction remains unresolved",
               "raw_count_status": "not established by stored-value integrality alone", "raw_counts_verified": False, "rows_deleted_or_summed": False,
               "rna_model_eligible": False, "independence_of_all_other_groups_verified": False})


def audit(config: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=False)
    summaries = {}
    for cohort in COHORTS:
        roster = cell_roster(config, cohort)
        roster.to_csv(output / f"{cohort}.cells.tsv.gz", sep="\t", index=False)
        donor_rows = roster.groupby(["donor_id", "condition", "lineage"], sort=True).agg(n_cells=("cell_id", "size"), da_eligible=("da_eligible", "first"), network_eligible=("network_eligible", "first")).reset_index()
        donor_rows.to_csv(output / f"{cohort}.donor_lineages.tsv", sep="\t", index=False)
        summaries[cohort] = {"cells": len(roster), "donors": int(roster["donor_id"].nunique()), "identity_duplicates": 0}
    rna, rna_summary = audit_rna(config)
    rna.to_csv(output / "GSE244832.rna_identity.tsv.gz", sep="\t", index=False)
    write_json(output / "diagnostics.json", {"cohorts": summaries, "rna": rna_summary, "matrix_counts_read": False})


def count(config: dict, output: Path, cohort: str, lineage: str, pilot: bool) -> None:
    if not os.environ.get("SLURM_JOB_ID") or not os.environ.get("SLURM_JOB_NODELIST"):
        raise RuntimeError("fragment reconstruction requires a compute-node Slurm allocation")
    import snapatac2 as snap
    import snap_helpers
    spec = cohort_spec(config, cohort)
    snap_helpers.COHORT_SPEC[cohort] = spec
    candidate = project_path(config, str(config.get("v3_candidate_root", "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1")))
    peak_bed = candidate / "peaks/consensus" / f"{lineage}.bed"
    roster = cell_roster(config, cohort)
    roster = roster.loc[roster["lineage"] == lineage].copy()
    if roster.empty:
        raise ValueError("no cells in requested lineage")
    output.mkdir(parents=True, exist_ok=False)
    chunk_cells = int(config.get("count_chunk_cells", 500))
    if chunk_cells <= 0:
        raise ValueError("chunk size must be positive")
    with tempfile.TemporaryDirectory(prefix="fragment_recount_", dir=output) as temporary:
        temporary = Path(temporary)
        dataset, sources = snap_helpers.build_dataset(cohort, temporary / "source.h5ads")
        peak_matrix = None
        selected_dataset = None
        try:
            source_lineages = np.asarray(dataset.obs["lineage_v3"]).astype(str)
            source_donors = np.asarray(dataset.obs["donor_id"]).astype(str)
            selected_cells = np.flatnonzero(source_lineages == lineage)
            if pilot:
                selected_cells = selected_cells[pilot_indices(source_donors[selected_cells], int(config.get("pilot_cells", 1000)), int(config.get("seed", 42)), int(config.get("network", {}).get("min_cells_per_donor", 100)))]
            print(f"Selecting {len(selected_cells)} {cohort}/{lineage} nuclei from fragment-backed donors", flush=True)
            # ``out`` is essential: subset must never modify backing inputs.
            # Peak recount reads fragment obsm, not the historical 6M-column
            # tile X. Dropping X columns avoids copying tile metadata per donor.
            selected_dataset, _ = dataset.subset(selected_cells.tolist(), var_indices=[], out=str(temporary / "selected"))
            peak_matrix = snap.pp.make_peak_matrix(selected_dataset, peak_file=peak_bed, counting_strategy="fragment", chunk_size=chunk_cells, file=temporary / "peaks.h5ad", inplace=False)
            print(f"Recounted {peak_matrix.n_obs} nuclei × {peak_matrix.n_vars} fixed peaks", flush=True)
            obs = peak_matrix.obs
            donor_values = np.asarray(obs["donor_id"]).astype(str)
            lineage_values = np.asarray(obs["lineage_v3"]).astype(str)
            names = np.asarray(peak_matrix.obs_names).astype(str)
            indices = np.flatnonzero(lineage_values == lineage)
            ids = pd.MultiIndex.from_arrays([donor_values[indices], names[indices]])
            roster = roster.set_index(["donor_id", "raw_barcode"]).loc[ids].reset_index()
            donors = sorted(roster["donor_id"].unique())
            donor_lookup = {donor: index for index, donor in enumerate(donors)}
            donor_codes = roster["donor_id"].map(donor_lookup).to_numpy(dtype=np.int64)
            raw_sum = sp.csr_matrix((len(donors), peak_matrix.n_vars), dtype=np.int64)
            binary_sum = raw_sum.copy()
            balanced_raw = raw_sum.copy()
            balanced_binary = raw_sum.copy()
            balanced_mask, balanced_n = balanced_cell_mask(roster, int(config.get("seed", 42))) if not pilot else (np.zeros(len(roster), dtype=bool), 0)
            parts = []
            chunks = output / "cell_chunks"
            chunks.mkdir()
            for first in range(0, len(indices), chunk_cells):
                selected = indices[first:first + chunk_cells]
                # AnnDataSet has already restricted rows to this lineage.
                # Contiguous slices avoid HDF5 point-selection overhead.
                row_selection = slice(int(selected[0]), int(selected[-1]) + 1) if len(selected) == 1 or np.all(np.diff(selected) == 1) else selected
                matrix = sp.csr_matrix(peak_matrix.X[row_selection, :])
                raw, binary = aggregate_counts(matrix, donor_codes[first:first + len(selected)], len(donors))
                matrix = matrix.astype(np.int64)
                raw_sum += raw
                binary_sum += binary
                balanced_chunk = balanced_mask[first:first + len(selected)]
                if balanced_chunk.any():
                    br, bb = aggregate_counts(matrix[balanced_chunk], donor_codes[first:first + len(selected)][balanced_chunk], len(donors))
                    balanced_raw += br
                    balanced_binary += bb
                path = chunks / f"chunk_{len(parts):06d}.npz"
                sp.save_npz(path, matrix)
                parts.append({"path": str(path.relative_to(output)), "row_start": first, "rows": matrix.shape[0], "nnz": matrix.nnz})
                if len(parts) % 20 == 0:
                    print(f"Exported {min(first + chunk_cells, len(indices))}/{len(indices)} nuclei", flush=True)
            peaks = pd.DataFrame({"column_index_1based": np.arange(1, peak_matrix.n_vars + 1), "peak_coordinate": np.asarray(peak_matrix.var_names).astype(str)})
            coordinates = peaks["peak_coordinate"].str.extract(r"^(?P<chrom>[^:]+):(?P<start0>[0-9]+)-(?P<end>[0-9]+)$")
            if coordinates.isna().any().any():
                raise ValueError("invalid consensus peak coordinates")
            coordinates[["start0", "end"]] = coordinates[["start0", "end"]].astype(np.int64)
            if (coordinates["start0"] >= coordinates["end"]).any():
                raise ValueError("invalid half-open consensus peak interval")
            peaks["peak_id"] = peaks["peak_coordinate"]
            peaks = pd.concat([peaks, coordinates], axis=1)
            reference_path = candidate / "counts" / cohort
            reference_peaks = pd.read_csv(reference_path / f"{lineage}.peaks.tsv", sep="\t")
            if not peaks["peak_coordinate"].equals(reference_peaks["peak_coordinate"]):
                raise ValueError("per-cell recount peak order differs from frozen v3")
            conservation = "pilot_subset_only"
            if not pilot:
                reference_donors = pd.read_csv(reference_path / f"{lineage}.donors.tsv", sep="\t")
                if donors != reference_donors["donor_id"].tolist():
                    raise ValueError("donor order differs from frozen v3")
                with gzip.open(reference_path / f"{lineage}.mtx.gz", "rb") as handle:
                    reference = mmread(handle).tocsr()
                if reference.shape != raw_sum.shape or (reference != raw_sum).nnz:
                    raise ValueError("per-cell sums differ from frozen v3 donor fragment counts")
                if roster.groupby("donor_id").size().tolist() != reference_donors["n_cells"].tolist():
                    raise ValueError("cell numbers differ from frozen v3")
                conservation = "exact_all_donor_peak_cells"
            # Only combine chunks when estimated peak CSR memory stays bounded.
            nnz = sum(part["nnz"] for part in parts)
            estimated_bytes = 3 * (nnz * 12 + (len(indices) + 1) * 8)
            limit = int(float(config.get("max_csr_assembly_gib", 8)) * 1024**3)
            if estimated_bytes <= limit:
                matrix = sp.vstack([sp.load_npz(output / part["path"]) for part in parts], format="csr")
                sp.save_npz(output / "cell_counts.npz", matrix)
                del matrix
            for name, matrix in (("donor_counts", raw_sum), ("donor_binary_counts", binary_sum)):
                with gzip.open(output / f"{name}.mtx.gz", "wb") as handle:
                    mmwrite(handle, matrix, field="integer")
            roster.to_csv(output / "cells.tsv", sep="\t", index=False)
            peaks.to_csv(output / "peaks.tsv", sep="\t", index=False)
            donor_table = roster.groupby(["donor_id", "condition"], sort=True).agg(n_cells=("cell_id", "size"), da_eligible=("da_eligible", "first"), network_eligible=("network_eligible", "first")).reset_index()
            donor_table["total_counted_fragments"] = np.asarray(raw_sum.sum(axis=1)).ravel()
            donor_table["total_detected_cell_peak_pairs"] = np.asarray(binary_sum.sum(axis=1)).ravel()
            donor_table["contrast_eligible"] = donor_table["n_cells"] >= 20
            donor_table["cohort"] = cohort
            donor_table["lineage"] = lineage
            donor_table["exclusion_reason"] = np.where(donor_table["contrast_eligible"], "", "fewer_than_20_cells")
            donor_table.to_csv(output / "donors.tsv", sep="\t", index=False)
            if not pilot:
                for name, matrix in (("fragment_counts", raw_sum), ("binary_counts", binary_sum)):
                    export = output.parents[2] / name / cohort
                    export.mkdir(parents=True, exist_ok=True)
                    paths = [export / f"{lineage}.{suffix}" for suffix in ("mtx.gz", "donors.tsv", "peaks.tsv")]
                    if any(path.exists() for path in paths):
                        raise FileExistsError(f"refusing to overwrite donor count export: {export}/{lineage}")
                    with gzip.open(paths[0], "wb") as handle:
                        mmwrite(handle, matrix, field="integer")
                    donor_table.to_csv(paths[1], sep="\t", index=False)
                    peaks.to_csv(paths[2], sep="\t", index=False)
                balanced_roster = roster.loc[balanced_mask]
                balanced_roster[["cell_id", "donor_id", "condition"]].to_csv(output / "balanced_cells.tsv", sep="\t", index=False)
                keep_donors = donor_table["donor_id"].isin(balanced_roster["donor_id"].unique()).to_numpy()
                balanced_table = donor_table.loc[keep_donors].copy()
                balanced_table["n_cells"] = balanced_n
                balanced_table["network_eligible"] = balanced_n >= int(config.get("network", {}).get("min_cells_per_donor", 100))
                balanced_table["total_counted_fragments"] = np.asarray(balanced_raw[keep_donors].sum(axis=1)).ravel()
                balanced_table["total_detected_cell_peak_pairs"] = np.asarray(balanced_binary[keep_donors].sum(axis=1)).ravel()
                for name, matrix in (("balanced_counts", balanced_raw), ("balanced_binary_counts", balanced_binary)):
                    export = output.parents[2] / name / cohort
                    export.mkdir(parents=True, exist_ok=True)
                    paths = [export / f"{lineage}.{suffix}" for suffix in ("mtx.gz", "donors.tsv", "peaks.tsv")]
                    if any(path.exists() for path in paths):
                        raise FileExistsError(f"refusing balanced donor count overwrite: {export}/{lineage}")
                    with gzip.open(paths[0], "wb") as handle:
                        mmwrite(handle, matrix[keep_donors], field="integer")
                    balanced_table.to_csv(paths[1], sep="\t", index=False)
                    peaks.to_csv(paths[2], sep="\t", index=False)
            write_json(output / "diagnostics.json", {"cohort": cohort, "lineage": lineage, "pilot": pilot, "seed": int(config.get("seed", 42)), "cells": len(indices), "peaks": peak_matrix.n_vars, "donors": donors, "counting_unit": "deduplicated_fragment_record", "multiplicity_used": False, "v3_conservation": conservation, "balanced_sensitivity_cells_per_eligible_donor": balanced_n, "balanced_sensitivity_total_cells": int(balanced_mask.sum()), "cell_count_format": "csr_npz" if (output / "cell_counts.npz").exists() else "chunked_csr_npz", "chunks": parts, "estimated_csr_assembly_bytes": estimated_bytes, "peak_bed_sha256": hashlib.sha256(peak_bed.read_bytes()).hexdigest(), "snapatac2_version": snap.__version__})
        finally:
            if peak_matrix is not None:
                peak_matrix.close()
            if selected_dataset is not None:
                selected_dataset.close()
            snap_helpers.close_dataset(dataset, sources)


def main() -> None:
    from threadpoolctl import threadpool_limits
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["audit", "pilot", "count", "verify", "rna-vectors"])
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--cohort", choices=tuple(COHORTS))
    parser.add_argument("--lineage", choices=["hepatocyte", "stellate", "macrophage"])
    parser.add_argument("--input-mode", choices=["pilot", "count"], default="pilot", help="which existing counts to verify")
    args = parser.parse_args()
    config = configuration(args.config)
    check_output_root(config, args.output_root)
    # scipy.io MatrixMarket uses its own native thread pool; OMP alone does
    # not constrain it. The controller also constrains BLAS operations.
    with threadpool_limits(limits=int(os.environ.get("OMP_NUM_THREADS", "2"))):
        if args.command == "audit":
            audit(config, args.output_root / "audit")
        elif args.command == "rna-vectors":
            audit_rna_vectors(config, args.output_root / "audit")
        else:
            if not args.cohort or not args.lineage:
                parser.error("pilot/count require --cohort and --lineage")
            if args.command == "verify":
                verify_fragment_records(config, args.output_root / args.input_mode / args.cohort / args.lineage, args.cohort)
            else:
                count(config, args.output_root / args.command / args.cohort / args.lineage, args.cohort, args.lineage, args.command == "pilot")


if __name__ == "__main__":
    main()
