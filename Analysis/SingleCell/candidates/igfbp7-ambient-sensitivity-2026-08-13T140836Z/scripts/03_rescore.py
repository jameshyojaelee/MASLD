#!/usr/bin/env python
"""Step 3/3 of the IGFBP7 ambient-RNA SENSITIVITY ANALYSIS.

SCOPE (binding): SENSITIVITY ANALYSIS ONLY. The 117 frozen Hotspot programs,
their memberships and their L1 weights are read-only here. Nothing is
rediscovered, refit, reweighted, renamed or re-selected. Every program is
re-SCORED with its EXACT frozen member-gene set on ambient-corrected counts and
compared against its stored score. `program_registry_v2.tsv` and
`program_membership_v2.tsv` are opened read-only and never written.

Method (mirrors the already-approved in-repo sensitivity recipe in
scripts/manuscript/noncoding_resource/rescore_programs_without_lncrna.py):
  * rebuild the frozen hepatocyte Hotspot run (same cells, same genes, same
    latent space, same k=30 kNN graph, same DANB model);
  * confirm the rebuild reproduces the STORED per-cell program scores at
    r >= 0.995 before any conclusion is drawn;
  * score each frozen program twice on the rebuilt graph, once on raw counts and
    once on decontX-corrected counts, and add only the correction delta to the
    exact stored score;
  * collapse cells to BIOLOGICAL DONORS by the pooled-cell donor mean (the
    frozen primary score definition), never by sequencing run;
  * refit the frozen primary model score ~ stage_ordinal + factor(dataset) with
    homoskedastic OLS SEs and HC3 sandwich SEs, and recompute BH over the
    complete 117-program family with the corrected hepatocyte p-values
    substituted in place of the stored ones.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import hotspot
import numpy as np
import pandas as pd
import scipy.sparse as sp
from hotspot import modules as hotspot_modules
from scipy import stats

SEED = 20260813
np.random.seed(SEED)

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
WORK = CAND / "work"
RES = CAND / "results"
RES.mkdir(parents=True, exist_ok=True)

HS_ROOT = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
PROGRAM_ROOT = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-candidate-2026-08-07/hotspot"
)
REGISTRY = PROGRAM_ROOT / "program_registry_v2.tsv"
MEMBERSHIP = PROGRAM_ROOT / "program_membership_v2.tsv"
READY = PROGRAM_ROOT / "READY"
SCORE_501 = ROOT / "Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py"
META = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"

PAIRINGS = {
    "GSE244832": ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}
PRIMARY_STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
CELL_TYPE = "hepatocytes"
N_NEIGHBORS = 30
MINIMUM_BASELINE_REPRODUCTION_R = 0.995

# ---------------------------------------------------------------- depth mode --
# The DANB null centers each cell against an expected value mu = t_gene * t_cell /
# total, where t_cell is the cell's UMI depth. The first run passed the atlas's
# RAW n_counts as t_cell when scoring BOTH the raw and the decontaminated matrix,
# so a cell that lost 20% of its UMIs was centered against its pre-correction
# depth and its residual went negative for every gene. Because decontX removes a
# stage-correlated amount (donor-level frac_removed ~ stage + dataset:
# beta = +0.0177, p = 0.025, n = 62), that mechanical shift is itself
# stage-correlated, and it depressed all 30 programs including zero-ambient ones.
# See AUDIT_FINDINGS.md F1.
#
#   raw_depth  reproduces the original run byte-for-byte. Default, so nothing
#              already written can change under a rerun.
#   dec_depth  scales each cell's depth by its measured retained fraction
#              umi_decont / umi_raw before scoring the decontaminated matrix.
#
# The retained FRACTION is applied to the atlas depth rather than substituting
# umi_decont outright: umi_decont is summed over each dataset's own exported gene
# set, which is not the frozen atlas gene universe, so the absolute totals are not
# interchangeable while the per-cell ratio is.
DEPTH_MODE = os.environ.get("RESCORE_DEPTH_MODE", "raw_depth")
if DEPTH_MODE not in ("raw_depth", "dec_depth"):
    raise RuntimeError(f"RESCORE_DEPTH_MODE must be raw_depth or dec_depth, got {DEPTH_MODE!r}")
OUT_SUFFIX = "" if DEPTH_MODE == "raw_depth" else "__dec_depth"


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError(msg)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_501():
    spec = importlib.util.spec_from_file_location("hotspot_run_501", SCORE_501)
    require(spec is not None and spec.loader is not None, "cannot load 501 producer")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def donor_map() -> dict[str, str]:
    out: dict[str, str] = {}
    for dataset, path in PAIRINGS.items():
        table = pd.read_csv(path, dtype=str)
        for row in table.itertuples(index=False):
            donor = f"{dataset}_{row.donor_id}"
            for sample in str(row.rna_srrs).split(";"):
                if sample.strip():
                    out[sample.strip()] = donor
    return out


def donor_metadata(mapping: dict[str, str]) -> pd.DataFrame:
    table = pd.read_csv(META, sep="\t", dtype=str)
    required = {"sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"}
    require(required <= set(table), "donor metadata schema drift")
    table = table[list(required)].copy()
    table["donor"] = table["sample"].map(mapping).fillna(table["sample"])
    table["exclude"] = (
        table["exclude_stage_analysis"].str.lower().isin({"true", "t", "1", "yes"})
    )
    rows = []
    for donor, group in table.groupby("donor", sort=False):
        datasets = sorted(set(group["dataset"].dropna()) - {""})
        stages = sorted(set(group["disease_stage_coarse"].dropna()) - {""})
        require(len(datasets) <= 1 and len(stages) <= 1, f"metadata conflict: {donor}")
        rows.append(
            {
                "donor": donor,
                "dataset": datasets[0] if datasets else None,
                "stage": stages[0] if stages else None,
                "exclude": bool(group["exclude"].any()),
            }
        )
    return pd.DataFrame(rows)


def stage_fit(scores: pd.DataFrame, metadata: pd.DataFrame) -> dict:
    """Frozen primary model: score ~ stage_ordinal + factor(dataset).

    Returns OLS beta/SE/p and the HC3 sandwich SE/p on the stage coefficient,
    with a two-sided residual-df t test, exactly as the frozen analysis
    specification defines them.
    """
    frame = scores.merge(metadata, on="donor", how="inner", validate="one_to_one")
    frame = frame[(~frame["exclude"]) & frame["stage"].isin(PRIMARY_STAGES)].copy()
    frame["stage_ordinal"] = frame["stage"].map(PRIMARY_STAGES).astype(float)
    dummies = pd.get_dummies(frame["dataset"], drop_first=True, dtype=float)
    X = np.column_stack(
        [np.ones(len(frame)), frame["stage_ordinal"].to_numpy(float), dummies.to_numpy(float)]
    )
    y = frame["score"].to_numpy(float)
    rank = np.linalg.matrix_rank(X)
    require(rank == X.shape[1], "rank-deficient stage model")
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    n, k = X.shape
    df = n - k
    sigma2 = float(resid @ resid) / df
    se_ols = float(np.sqrt(sigma2 * XtX_inv[1, 1]))
    H = X @ XtX_inv @ X.T
    h = np.clip(np.diag(H), 0.0, 1 - 1e-12)
    omega = (resid / (1.0 - h)) ** 2
    meat = X.T @ (omega[:, None] * X)
    V_hc3 = XtX_inv @ meat @ XtX_inv
    se_hc3 = float(np.sqrt(V_hc3[1, 1]))
    t_ols = float(beta[1]) / se_ols
    t_hc3 = float(beta[1]) / se_hc3
    return {
        "beta": float(beta[1]),
        "se": se_ols,
        "pvalue": float(2 * stats.t.sf(abs(t_ols), df)),
        "hc3_se": se_hc3,
        "hc3_pvalue": float(2 * stats.t.sf(abs(t_hc3), df)),
        "max_leverage": float(h.max()),
        "residual_df": int(df),
        "n_donors": int(n),
        "n_datasets": int(frame["dataset"].nunique()),
    }


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ps = p[order]
    q = np.minimum.accumulate((ps * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def read_corrected(dataset: str, n_score_genes: int) -> tuple[np.ndarray, sp.csc_matrix, sp.csc_matrix]:
    dsdir = WORK / dataset
    dims = json.loads((dsdir / "dec_dims.json").read_text(encoding="utf-8"))
    if isinstance(dims, list):
        dims = dims[0]
    ncells = int(dims["ncells"])
    require(int(dims["ngenes_score"]) == n_score_genes, "score gene axis drift")
    cells = pd.read_csv(dsdir / "dec_cells.csv")["cell_id"].astype(str).to_numpy()
    require(len(cells) == ncells, "corrected cell census drift")

    def load(prefix: str, nnz: int) -> sp.csc_matrix:
        i = np.fromfile(dsdir / f"{prefix}_indices.bin", dtype="<i4")
        p = np.fromfile(dsdir / f"{prefix}_indptr.bin", dtype="<i4")
        x = np.fromfile(dsdir / f"{prefix}_data.bin", dtype="<f8")
        require(len(i) == nnz and len(x) == nnz and len(p) == ncells + 1,
                f"{prefix} binary shape drift for {dataset}")
        return sp.csc_matrix((x, i, p), shape=(n_score_genes, ncells))

    dec = load("dec", int(dims["dec_nnz"]))
    raw = load("raw", int(dims["raw_nnz"]))
    require(np.isclose(dec.data.sum(), float(dims["dec_sum"]), rtol=1e-9),
            f"decontaminated checksum drift for {dataset}")
    require(np.isclose(raw.data.sum(), float(dims["raw_sum"]), rtol=1e-9),
            f"raw checksum drift for {dataset}")
    return cells, dec, raw


def retained_fraction(cell_index: pd.Index) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Per-cell umi_decont / umi_raw, aligned to the frozen atlas cell order.

    Returns (retained, uncorrected_positions, uncorrected_datasets).

    A dataset where decontX failed carries raw counts through, so every one of its
    cells has retained == 1 exactly and its depth is unchanged. Those cells are
    also identified here because the original provenance record reported
    `corrected_coverage: 1.0` and `n_uncorrected_cells: 0`, which measured census
    completeness rather than correction completeness.
    """
    retained = np.ones(len(cell_index), dtype=np.float64)
    seen = np.zeros(len(cell_index), dtype=bool)
    uncorrected_datasets: list[str] = []
    uncorrected_mask = np.zeros(len(cell_index), dtype=bool)

    for path in sorted(RES.glob("02_contamination_per_cell__*.tsv.gz")):
        dataset = path.name.split("__")[1].replace(".tsv.gz", "")
        d = pd.read_csv(path, sep="\t")
        pos = cell_index.get_indexer(d["cell_id"].astype(str))
        keep = pos >= 0
        require(keep.all(), f"{dataset}: contamination rows absent from the frozen universe")
        r = (d["umi_decont"].to_numpy(float) / d["umi_raw"].to_numpy(float))
        # decontX failure, a zero-depth cell, or any non-finite ratio leaves the
        # depth untouched rather than silently scaling it to nonsense.
        bad = ~np.isfinite(r) | (r <= 0)
        r[bad] = 1.0
        retained[pos] = r
        seen[pos] = True

        log = RES / f"02_decontx_run_log__{dataset}.tsv"
        if log.exists():
            flag = str(pd.read_csv(log, sep="\t").loc[0, "corrected"]).strip().upper()
            if flag in ("FALSE", "0", "NO"):
                uncorrected_datasets.append(dataset)
                uncorrected_mask[pos] = True
                require(
                    np.allclose(r, 1.0),
                    f"{dataset} is flagged uncorrected but its retained fraction is not 1",
                )

    require(seen.all(), "some frozen atlas cells have no contamination record")
    return retained, np.flatnonzero(uncorrected_mask), uncorrected_datasets


def main() -> None:
    # ---- frozen contract ------------------------------------------------
    ready = pd.read_csv(READY, sep="\t")
    require(len(ready) == 1, "READY seal malformed")
    require(str(ready.loc[0, "status"]) == "ready_for_external_testing", "READY status drift")
    require(sha256(REGISTRY) == str(ready.loc[0, "registry_sha256"]), "registry hash drift")
    require(sha256(MEMBERSHIP) == str(ready.loc[0, "membership_table_sha256"]), "membership hash drift")
    registry = pd.read_csv(REGISTRY, sep="\t")
    require(len(registry) == 117, "expected 117 frozen programs")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    print("[rescore] frozen registry + membership hashes verified against READY", flush=True)

    # BH reproduction check on the stored family before anything else.
    stored_q = registry["primary_qvalue"].to_numpy(float)
    require(
        np.allclose(bh(registry["primary_pvalue"].to_numpy(float)), stored_q, atol=1e-12),
        "cannot reproduce the stored BH family over 117 programs",
    )
    print("[rescore] stored BH family over 117 programs reproduced exactly", flush=True)

    # ---- rebuild the frozen hepatocyte Hotspot run ----------------------
    run501 = load_501()
    run_metadata = json.loads(
        (HS_ROOT / CELL_TYPE / "run_metadata.json").read_text(encoding="utf-8")
    )
    adata = run501.load_atlas(CELL_TYPE, smoke=False)
    excluded = run_metadata.get("exclude_datasets") or []
    if excluded:
        adata = adata[~adata.obs["dataset"].isin(set(excluded))].copy()
    adata = run501.strip_confounders(adata)
    adata = run501.filter_detected(adata, min_frac=0.01)
    require(adata.n_obs == int(run_metadata["n_cells"]), "frozen cell census drift")
    require(adata.n_vars == int(run_metadata["n_genes_kept"]), "frozen gene census drift")
    run501.ensure_raw_layer(adata)
    latent = run501.resolve_latent(adata)
    hs = hotspot.Hotspot(
        adata, layer_key="counts", model="danb",
        latent_obsm_key=latent, umi_counts_obs_key="n_counts",
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
    print(f"[rescore] rebuilt frozen run: {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

    cell_index = pd.Index(adata.obs_names.astype(str))

    # ---- assemble the ambient-corrected count matrix --------------------
    score_genes = (WORK / "score_genes.txt").read_text(encoding="utf-8").split()
    manifest = json.loads((WORK / "manifest.json").read_text(encoding="utf-8"))
    gene_pos = {g: i for i, g in enumerate(score_genes)}

    dec_blocks, raw_blocks, cell_blocks = [], [], []
    for dataset in sorted(manifest["datasets"]):
        if not (WORK / dataset / "dec_dims.json").exists():
            print(f"[rescore] {dataset}: no corrected counts written; skipping", flush=True)
            continue
        cells, dec, raw = read_corrected(dataset, len(score_genes))
        dec_blocks.append(dec)
        raw_blocks.append(raw)
        cell_blocks.append(cells)
        print(f"[rescore] {dataset}: loaded {len(cells):,} corrected cells", flush=True)

    all_cells = np.concatenate(cell_blocks)
    require(len(set(all_cells)) == len(all_cells), "duplicate corrected cells across datasets")
    pos = cell_index.get_indexer(all_cells)
    require((pos >= 0).all(), "corrected cells absent from the frozen universe")
    coverage = float(len(all_cells) / adata.n_obs)
    print(f"[rescore] corrected-count coverage of the frozen universe: {coverage:.6f}", flush=True)
    require(
        len(all_cells) == adata.n_obs,
        "corrected counts do not cover the frozen cell universe exactly; the "
        "frozen cell and donor census must be preserved",
    )
    # reorder the concatenated blocks into frozen atlas cell order
    inverse = np.empty(adata.n_obs, dtype=np.int64)
    inverse[pos] = np.arange(len(all_cells))
    # Column-permute as CSC, then hold as CSR: every downstream access slices
    # GENE ROWS, which is fast on CSR and pathological on CSC.
    dec_full = sp.hstack(dec_blocks, format="csc")[:, inverse].tocsr()
    raw_full = sp.hstack(raw_blocks, format="csc")[:, inverse].tocsr()
    del dec_blocks, raw_blocks
    covered = np.ones(adata.n_obs, dtype=bool)

    # ---- per-cell DANB depth -------------------------------------------
    retained, uncovered, uncorrected_datasets = retained_fraction(cell_index)
    umi_raw_depth = hs.umi_counts.values.astype(np.float64)
    if DEPTH_MODE == "dec_depth":
        umi_dec_depth = umi_raw_depth * retained
        require((umi_dec_depth > 0).all(), "a decontaminated depth is non-positive")
    else:
        umi_dec_depth = umi_raw_depth
    correction_coverage = float(1.0 - len(uncovered) / adata.n_obs)
    print(
        f"[rescore] depth mode {DEPTH_MODE}; retained fraction "
        f"median {np.median(retained):.4f} p05 {np.quantile(retained, 0.05):.4f}; "
        f"uncorrected cells {len(uncovered):,} ({100*(1-correction_coverage):.2f}%) "
        f"from {uncorrected_datasets or 'none'}",
        flush=True,
    )

    # sanity: the transported raw counts must equal the atlas counts
    check_genes = [g for g in ["IGFBP7", "BICC1", "ALB", "DCN"] if g in gene_pos and g in adata.var_names]
    raw_check = {}
    for g in check_genes:
        atlas_col = np.asarray(
            adata[:, g].layers["counts"].todense()
            if sp.issparse(adata[:, g].layers["counts"])
            else adata[:, g].layers["counts"]
        ).ravel()
        transported = np.asarray(raw_full[gene_pos[g], :].todense()).ravel()
        agree = float(np.abs(atlas_col[covered] - transported[covered]).max())
        raw_check[g] = agree
        print(f"[rescore] transport check {g}: max abs diff on covered cells = {agree:g}", flush=True)
        require(agree < 1e-6, f"transported raw counts disagree with the atlas for {g}")

    # ---- score every frozen hepatocyte program twice --------------------
    stored = pd.read_parquet(HS_ROOT / CELL_TYPE / "cell_scores.parquet")
    stored["module"] = stored["module"].astype(int)
    mapping = donor_map()
    metadata = donor_metadata(mapping)
    donors = (
        adata.obs["sample"].astype(str).map(mapping).fillna(adata.obs["sample"].astype(str))
    ).to_numpy()

    reg_hep = registry[registry["cell_type"] == CELL_TYPE].set_index("module")
    mem_hep = membership[membership["cell_type"] == CELL_TYPE]

    summary_rows = []
    donor_frames = []
    modules = sorted(int(m) for m in reg_hep.index)
    for module in modules:
        genes = mem_hep.loc[mem_hep["module"].astype(int) == module, "source_gene"].astype(str).tolist()
        require(
            len(genes) == int(reg_hep.loc[module, "n_source_genes"]),
            f"frozen membership size drift for module {module}",
        )
        require(set(genes) <= set(adata.var_names.astype(str)), f"missing member genes: {module}")
        require(set(genes) <= set(gene_pos), f"member genes absent from corrected export: {module}")

        gi = np.array([gene_pos[g] for g in genes])
        counts_raw = np.asarray(raw_full[gi, :].todense(), dtype=np.float64)
        counts_dec = np.asarray(dec_full[gi, :].todense(), dtype=np.float64)
        if len(uncovered):
            counts_dec[:, uncovered] = counts_raw[:, uncovered]

        # The raw arm always uses the atlas depth, so the reproduction gate below
        # still validates the rebuild against the stored scores. Only the
        # decontaminated arm's depth changes, and only under dec_depth.
        score_raw = hotspot_modules.compute_scores(
            counts_raw, hs.model, umi_raw_depth, hs.neighbors.values, hs.weights.values)
        score_dec = hotspot_modules.compute_scores(
            counts_dec, hs.model, umi_dec_depth, hs.neighbors.values, hs.weights.values)

        stored_module = (
            stored[stored["module"] == module].set_index("cell_id")["score"].reindex(cell_index)
        )
        require(not stored_module.isna().any(), f"stored score join failed: {module}")
        stored_vec = stored_module.to_numpy(float)
        repro_r = float(np.corrcoef(score_raw, stored_vec)[0, 1])
        repro_max = float(np.max(np.abs(score_raw - stored_vec)))
        print(
            f"REPRODUCTION\t{CELL_TYPE}\t{module}\tr={repro_r:.12g}\tmax_abs_error={repro_max:.12g}",
            flush=True,
        )
        require(
            repro_r >= MINIMUM_BASELINE_REPRODUCTION_R,
            f"frozen Hotspot score not reproduced for module {module}",
        )
        corrected_vec = stored_vec + (score_dec - score_raw)

        cell_frame = pd.DataFrame(
            {"donor": donors, "stored": stored_vec, "corrected": corrected_vec}
        )
        donor = cell_frame.groupby("donor", as_index=False)[["stored", "corrected"]].mean()
        donor = donor.sort_values("donor").reset_index(drop=True)

        fit_stored = stage_fit(donor[["donor", "stored"]].rename(columns={"stored": "score"}), metadata)
        fit_corr = stage_fit(donor[["donor", "corrected"]].rename(columns={"corrected": "score"}), metadata)

        reg_beta = float(reg_hep.loc[module, "primary_beta"])
        reg_se = float(reg_hep.loc[module, "primary_se"])
        reg_hc3 = float(reg_hep.loc[module, "primary_hc3_se"])
        require(
            np.isclose(fit_stored["beta"], reg_beta, atol=1e-8, rtol=1e-8),
            f"stage beta drift vs frozen registry for module {module}",
        )

        # donor-level agreement between stored and ambient-corrected scores
        pear = float(np.corrcoef(donor["stored"], donor["corrected"])[0, 1])
        spear = float(stats.spearmanr(donor["stored"], donor["corrected"]).statistic)
        cell_pear = float(np.corrcoef(stored_vec, corrected_vec)[0, 1])

        summary_rows.append(
            {
                "program_uid": reg_hep.loc[module, "program_uid"],
                "cell_type": CELL_TYPE,
                "module": module,
                "module_name": reg_hep.loc[module, "module_name"],
                "n_member_genes": len(genes),
                "robust_display": reg_hep.loc[module, "robust_display"],
                "primary_selected": reg_hep.loc[module, "primary_selected"],
                "score_reproduction_r": repro_r,
                "score_reproduction_max_abs_error": repro_max,
                "cell_pearson_stored_vs_corrected": cell_pear,
                "donor_pearson_stored_vs_corrected": pear,
                "donor_spearman_stored_vs_corrected": spear,
                "stored_registry_beta": reg_beta,
                "stored_registry_se": reg_se,
                "stored_registry_hc3_se": reg_hc3,
                "stored_registry_pvalue": float(reg_hep.loc[module, "primary_pvalue"]),
                "stored_registry_hc3_pvalue": float(reg_hep.loc[module, "primary_hc3_pvalue"]),
                "stored_registry_qvalue": float(reg_hep.loc[module, "primary_qvalue"]),
                "refit_stored_beta": fit_stored["beta"],
                "refit_stored_se": fit_stored["se"],
                "refit_stored_hc3_se": fit_stored["hc3_se"],
                "refit_stored_pvalue": fit_stored["pvalue"],
                "refit_stored_hc3_pvalue": fit_stored["hc3_pvalue"],
                "corrected_beta": fit_corr["beta"],
                "corrected_se": fit_corr["se"],
                "corrected_pvalue": fit_corr["pvalue"],
                "corrected_hc3_se": fit_corr["hc3_se"],
                "corrected_hc3_pvalue": fit_corr["hc3_pvalue"],
                "corrected_max_leverage": fit_corr["max_leverage"],
                "corrected_residual_df": fit_corr["residual_df"],
                "direction_preserved": bool(np.sign(fit_corr["beta"]) == np.sign(fit_stored["beta"])),
                "n_donors": fit_corr["n_donors"],
                "n_datasets": fit_corr["n_datasets"],
                "mean_stored_score": float(donor["stored"].mean()),
                "mean_corrected_score": float(donor["corrected"].mean()),
            }
        )
        d = donor.copy()
        d.insert(0, "module", module)
        d.insert(0, "cell_type", CELL_TYPE)
        donor_frames.append(d)
        print(
            f"[rescore] module {module:>2} {reg_hep.loc[module, 'module_name'][:34]:<36} "
            f"beta {fit_stored['beta']:+.4f} -> {fit_corr['beta']:+.4f}  "
            f"donor r={pear:.4f} rho={spear:.4f}",
            flush=True,
        )

    summary = pd.DataFrame(summary_rows)

    # ---- BH over the complete frozen 117-program family -----------------
    # The frozen multiplicity rule is BH over all 117 programs. The corrected
    # hepatocyte p-values are substituted for their stored counterparts; every
    # non-hepatocyte program keeps its stored p-value.
    for label, col in (("", "pvalue"), ("hc3_", "hc3_pvalue")):
        stored_p = registry[f"primary_{col}"].to_numpy(float).copy()
        key = registry["cell_type"].astype(str) + "::" + registry["module"].astype(str)
        pos = {k: i for i, k in enumerate(key)}
        p_new = stored_p.copy()
        for row in summary.itertuples(index=False):
            p_new[pos[f"{CELL_TYPE}::{row.module}"]] = getattr(row, f"corrected_{col}")
        q_new = bh(p_new)
        q_map = {k: q_new[i] for k, i in pos.items()}
        summary[f"corrected_{label}qvalue"] = [
            q_map[f"{CELL_TYPE}::{m}"] for m in summary["module"]
        ]

    # Outputs are suffixed under dec_depth so a rerun can never overwrite the
    # raw_depth results it is meant to be compared against.
    for path in (RES / f"03_program_stage_effects{OUT_SUFFIX}.tsv",
                 RES / f"03_donor_scores_stored_vs_corrected{OUT_SUFFIX}.tsv",
                 RES / f"03_provenance{OUT_SUFFIX}.json"):
        require(
            not path.exists() or DEPTH_MODE == "raw_depth",
            f"refusing to overwrite an existing output: {path}",
        )
    summary.to_csv(RES / f"03_program_stage_effects{OUT_SUFFIX}.tsv", sep="\t", index=False)
    pd.concat(donor_frames, ignore_index=True).to_csv(
        RES / f"03_donor_scores_stored_vs_corrected{OUT_SUFFIX}.tsv", sep="\t", index=False
    )

    provenance = {
        "seed": SEED,
        "cell_type": CELL_TYPE,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_neighbors": N_NEIGHBORS,
        "latent_key": latent,
        "depth_mode": DEPTH_MODE,
        "retained_fraction_median": float(np.median(retained)),
        "retained_fraction_p05": float(np.quantile(retained, 0.05)),
        # census_coverage is what the previous `corrected_coverage` actually
        # measured: every atlas cell is present. correction_coverage is the
        # fraction whose counts decontX actually modified. The two differ because
        # a dataset where decontX failed carries raw counts through.
        "census_coverage": coverage,
        "correction_coverage": correction_coverage,
        "n_uncorrected_cells": int(len(uncovered)),
        "uncorrected_datasets": uncorrected_datasets,
        "registry_sha256": sha256(REGISTRY),
        "membership_sha256": sha256(MEMBERSHIP),
        "raw_transport_check_max_abs_diff": raw_check,
        "bh_family_size": 117,
        "note": "SENSITIVITY ANALYSIS ONLY - frozen programs re-scored, never refit",
    }
    (RES / "03_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print("[rescore] DONE", flush=True)


if __name__ == "__main__":
    main()
